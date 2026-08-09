#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TInterpreter.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TView.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

// These functions are intentionally global: ROOT's TButton invokes its action
// through the interpreter while the animation loop observes these flags.
bool gSimulationPaused = false;
bool gExitRequested = false;
TButton* gStopButton = nullptr;

void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
void ExitSimulation() {
    gExitRequested = true;
    if (gApplication) gApplication->Terminate(0);
}

// Relativistic two-body electrodynamics in SI units.  Mutual fields are
// evaluated at retarded time; radiation reaction uses a local, order-reduced
// Landau-Lifshitz approximation.
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double epsilon0 = 8.8541878128e-12;
constexpr double c = 299792458.0;
constexpr double mu0 = 4.0 * pi * 1.0e-7;
constexpr double eCharge = 1.602176634e-19;
constexpr double electronMass = 9.1093837139e-31;
constexpr double positronMass = electronMass;
constexpr double bohrMagneton = 9.2740100657e-24;
constexpr double nuclearMagneton = 5.0507837461e-27;
constexpr double coulomb = 1.0 / (4.0 * pi * epsilon0);
constexpr double bohrRadius = 5.29177210903e-11;
constexpr double nuclearCutoff = 1.0e-14; // point-particle theory has failed below this scale
constexpr double hbar = 1.054571817e-34;

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    double squaredNorm() const { return x*x + y*y + z*z; }
    double norm() const { return std::sqrt(squaredNorm()); }
};

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
Vec3 unit(const Vec3& v) { return v / v.norm(); }

struct ElectromagneticField { Vec3 electric, magnetic; };

struct State {
    Vec3 electronPosition, positronPosition;
    Vec3 electronVelocity, positronVelocity;
    Vec3 electronAcceleration, positronAcceleration;
    Vec3 electronDipole, positronDipole; // classical magnetic moments, in J/T
    double time = 0;
    double radiatedEnergy = 0;
};

struct Frame {
    Vec3 electron, positron, electronDipole, positronDipole;
    double time, radius, radiatedEnergy, mechanicalEnergy, schottEnergy;
    double electronMechanicalEnergy, positronMechanicalEnergy;
};

enum class Phenomenon { DirectCollision, Scattering, ParaPositronium, OrthoPositronium };

struct InitialConditions {
    double relativeEnergy;
    double orbitalAngularMomentum;
    double predictedClosestApproach;
    double dipoleAlignment;
    Phenomenon phenomenon;
    std::uint64_t seed;
};

struct SimulationResult {
    std::vector<Frame> frames;
    InitialConditions initial;
};

struct PairGeometry {
    Vec3 electronMinusPositron;
    double distance;
    double inverseDistance;
    double inverseDistanceCubed;
    double inverseDistanceFourth;
};

PairGeometry pairGeometry(const State& s) {
    const Vec3 electronMinusPositron = s.electronPosition - s.positronPosition;
    const double distanceSquared = electronMinusPositron.squaredNorm();
    const double distance = std::sqrt(distanceSquared);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;
    return {electronMinusPositron, distance, inverseDistance,
            inverseDistanceSquared * inverseDistance,
            inverseDistanceSquared * inverseDistanceSquared};
}

double separation(const State& s) { return pairGeometry(s).distance; }

double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

struct ChargeKinematics { Vec3 position, velocity, acceleration; };

Vec3 lerp(const Vec3& first, const Vec3& second, double fraction) {
    return first + (second - first) * fraction;
}

ChargeKinematics interpolatedCharge(const State& older, const State& newer, bool electron, double time) {
    const double span = newer.time - older.time;
    const double fraction = span > 0.0 ? std::clamp((time - older.time) / span, 0.0, 1.0) : 1.0;
    const Vec3 oldPosition = electron ? older.electronPosition : older.positronPosition;
    const Vec3 newPosition = electron ? newer.electronPosition : newer.positronPosition;
    const Vec3 oldVelocity = electron ? older.electronVelocity : older.positronVelocity;
    const Vec3 newVelocity = electron ? newer.electronVelocity : newer.positronVelocity;
    const Vec3 oldAcceleration = electron ? older.electronAcceleration : older.positronAcceleration;
    const Vec3 newAcceleration = electron ? newer.electronAcceleration : newer.positronAcceleration;
    return {lerp(oldPosition, newPosition, fraction), lerp(oldVelocity, newVelocity, fraction),
            lerp(oldAcceleration, newAcceleration, fraction)};
}

// Mutual, retarded Lienard-Wiechert field of a moving point charge.
ElectromagneticField lienardWiechertField(const Vec3& observationPosition,
                                          const State& olderSourceState,
                                          const State& newerSourceState,
                                          bool sourceIsElectron,
                                          double sourceCharge) {
    ChargeKinematics source = interpolatedCharge(olderSourceState, newerSourceState, sourceIsElectron,
                                                  newerSourceState.time);
    double retardedTime = newerSourceState.time - (observationPosition - source.position).norm() / c;
    for (int iteration = 0; iteration < 6; ++iteration) {
        source = interpolatedCharge(olderSourceState, newerSourceState, sourceIsElectron, retardedTime);
        const double refinedTime = newerSourceState.time - (observationPosition - source.position).norm() / c;
        // The fixed-point iteration converges very quickly here (v/c is small).
        // Avoid doing all six rounds once the retarded time is already stable.
        if (std::abs(refinedTime - retardedTime) <= 1.0e-24) {
            retardedTime = refinedTime;
            break;
        }
        retardedTime = refinedTime;
    }
    const Vec3 displacement = observationPosition - source.position;
    const double distance = displacement.norm();
    const Vec3 direction = displacement / distance;
    const Vec3 beta = source.velocity / c;
    const double betaSquared = beta.squaredNorm();
    const double kappa = std::max(1.0e-8, 1.0 - dot(direction, beta));
    const Vec3 velocityField = (direction - beta) * ((1.0 - betaSquared) /
                              (kappa*kappa*kappa * distance*distance));
    const Vec3 accelerationField = cross(direction, cross(direction - beta, source.acceleration)) /
                                   (c*c * kappa*kappa*kappa * distance);
    const Vec3 electric = (velocityField + accelerationField) * (coulomb * sourceCharge);
    return {electric, cross(direction, electric) / c};
}

// Static magnetic field of a point dipole.  'at' is the field point.
Vec3 dipoleField(const Vec3& sourceToTarget, const PairGeometry& geometry, const Vec3& sourceDipole) {
    const Vec3 n = sourceToTarget * geometry.inverseDistance;
    return (n * (3.0 * sourceDipole.x*n.x + 3.0 * sourceDipole.y*n.y + 3.0 * sourceDipole.z*n.z) - sourceDipole)
           * (mu0 / (4.0 * pi) * geometry.inverseDistanceCubed);
}

// F = grad(m_target . B_source), with r directed from source to target.
Vec3 dipoleForce(const Vec3& sourceToTarget, const PairGeometry& geometry,
                 const Vec3& targetDipole, const Vec3& sourceDipole) {
    const Vec3 n = sourceToTarget * geometry.inverseDistance;
    const double mnTarget = targetDipole.x*n.x + targetDipole.y*n.y + targetDipole.z*n.z;
    const double mnSource = sourceDipole.x*n.x + sourceDipole.y*n.y + sourceDipole.z*n.z;
    const double dots = targetDipole.x*sourceDipole.x + targetDipole.y*sourceDipole.y + targetDipole.z*sourceDipole.z;
    return (sourceDipole * mnTarget + targetDipole * mnSource + n * (dots - 5.0*mnTarget*mnSource))
           * (3.0 * mu0 / (4.0 * pi) * geometry.inverseDistanceFourth);
}

void precessDipole(Vec3& dipole, const Vec3& field, double gyromagneticRatio, double dt) {
    const double fieldMagnitude = field.norm();
    if (fieldMagnitude == 0.0) return;
    const Vec3 axis = field / fieldMagnitude;
    const double angle = gyromagneticRatio * fieldMagnitude * dt;
    // Rodrigues rotation preserves the magnitude of the magnetic moment.
    const Vec3 rotated = dipole * std::cos(angle) + cross(axis, dipole) * std::sin(angle)
                       + axis * ((axis.x*dipole.x + axis.y*dipole.y + axis.z*dipole.z) * (1.0 - std::cos(angle)));
    dipole = rotated * (dipole.norm() / rotated.norm());
}

double gamma(const Vec3& velocity) {
    return 1.0 / std::sqrt(std::max(1.0e-15, 1.0 - velocity.squaredNorm() / (c*c)));
}

Vec3 momentum(const Vec3& velocity, double mass) { return velocity * (gamma(velocity) * mass); }

Vec3 velocityFromMomentum(const Vec3& momentum, double mass) {
    const double gammaFromMomentum = std::sqrt(1.0 + momentum.squaredNorm() / (mass*mass*c*c));
    return momentum / (gammaFromMomentum * mass);
}

Vec3 relativisticAcceleration(const Vec3& velocity, const Vec3& force, double mass) {
    const double velocityForce = dot(velocity, force);
    return (force - velocity * (velocityForce / (c*c))) / (gamma(velocity) * mass);
}

Vec3 lorentzForce(double charge, const Vec3& velocity, const ElectromagneticField& field) {
    return (field.electric + cross(velocity, field.magnetic)) * charge;
}

double lienardPower(double charge, const Vec3& velocity, const Vec3& acceleration) {
    const double gammaValue = gamma(velocity);
    const double transverseTerm = cross(velocity, acceleration).squaredNorm() / (c*c);
    return charge*charge * std::pow(gammaValue, 6) * (acceleration.squaredNorm() - transverseTerm) /
           (6.0 * pi * epsilon0 * c*c*c);
}

struct MutualForces { Vec3 electron, positron; };

MutualForces mutualForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 forceOnElectron = geometry.electronMinusPositron
                               * (-coulomb * eCharge * eCharge * geometry.inverseDistanceCubed);
    return {forceOnElectron, forceOnElectron * -1.0};
}

MutualForces landauLifshitzForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 relativeVelocity = s.electronVelocity - s.positronVelocity;
    const double radialVelocity = dot(geometry.electronMinusPositron, relativeVelocity);
    const double inverseDistanceFifth = geometry.inverseDistanceFourth * geometry.inverseDistance;
    const double strength = coulomb * eCharge * eCharge;
    // Total time derivative of F_e = -k e^2 r/r^3.  In the non-relativistic
    // reduction of order, F_RR = tau dF_ext/dt; it is derived from the
    // Abraham-Lorentz equation rather than fitted to the radiated power.
    const Vec3 electronForceDerivative = relativeVelocity * (-strength * geometry.inverseDistanceCubed)
                                       + geometry.electronMinusPositron
                                       * (3.0 * strength * radialVelocity * inverseDistanceFifth);
    const double tau = eCharge * eCharge / (6.0 * pi * epsilon0 * electronMass * c*c*c);
    const Vec3 electronReaction = electronForceDerivative * tau;
    return {electronReaction, electronReaction * -1.0};
}

// Relativistic predictor-corrector update. At positronium's initial v/c of
// about 0.005, the instantaneous Coulomb interaction is the controlled leading
// term; radiation reaction is the Landau-Lifshitz reduction of order.
void advance(State& s, double dt) {
    MutualForces forces = mutualForces(s);
    const MutualForces reactions = landauLifshitzForces(s);
    const Vec3 electronAcceleration = relativisticAcceleration(s.electronVelocity, forces.electron, electronMass);
    const Vec3 positronAcceleration = relativisticAcceleration(s.positronVelocity, forces.positron, positronMass);
    const double electronPower = lienardPower(-eCharge, s.electronVelocity, electronAcceleration);
    const double positronPower = lienardPower(eCharge, s.positronVelocity, positronAcceleration);
    Vec3 electronMomentum = momentum(s.electronVelocity, electronMass) + (forces.electron + reactions.electron) * (0.5 * dt);
    Vec3 positronMomentum = momentum(s.positronVelocity, positronMass) + (forces.positron + reactions.positron) * (0.5 * dt);
    State trial = s;
    trial.time += dt;
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    trial.electronPosition += trial.electronVelocity * dt;
    trial.positronPosition += trial.positronVelocity * dt;
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, forces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, forces.positron, positronMass);

    const MutualForces trialForces = mutualForces(trial);
    const MutualForces trialReactions = landauLifshitzForces(trial);
    const Vec3 trialElectronAcceleration = relativisticAcceleration(trial.electronVelocity, trialForces.electron,
                                                                    electronMass);
    const Vec3 trialPositronAcceleration = relativisticAcceleration(trial.positronVelocity, trialForces.positron,
                                                                    positronMass);
    const double trialElectronPower = lienardPower(-eCharge, trial.electronVelocity, trialElectronAcceleration);
    const double trialPositronPower = lienardPower(eCharge, trial.positronVelocity, trialPositronAcceleration);
    electronMomentum += (trialForces.electron + trialReactions.electron) * (0.5 * dt);
    positronMomentum += (trialForces.positron + trialReactions.positron) * (0.5 * dt);
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, trialForces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, trialForces.positron, positronMass);
    PairGeometry geometry = pairGeometry(trial);

    const Vec3 fieldAtElectron = dipoleField(geometry.electronMinusPositron, geometry, trial.positronDipole);
    const Vec3 fieldAtPositron = dipoleField(geometry.electronMinusPositron * -1.0, geometry, trial.electronDipole);
    precessDipole(trial.electronDipole, fieldAtElectron, -1.76085963023e11, dt);
    precessDipole(trial.positronDipole, fieldAtPositron, 1.76085963023e11, dt);
    trial.radiatedEnergy += 0.5 * (electronPower + positronPower
                                 + trialElectronPower + trialPositronPower) * dt;
    s = trial;
}

Frame makeFrame(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const double electronKinetic = (gamma(s.electronVelocity) - 1.0) * electronMass * c*c;
    const double positronKinetic = (gamma(s.positronVelocity) - 1.0) * positronMass * c*c;
    const double coulombPotential = -coulomb * eCharge*eCharge * geometry.inverseDistance;
    const MutualForces forces = mutualForces(s);
    const Vec3 electronAcceleration = relativisticAcceleration(s.electronVelocity, forces.electron, electronMass);
    const Vec3 positronAcceleration = relativisticAcceleration(s.positronVelocity, forces.positron, positronMass);
    const double tau = eCharge * eCharge / (6.0 * pi * epsilon0 * electronMass * c*c*c);
    // Near-field (Schott) energy required by the instantaneous radiation
    // balance: d(E_mech + E_rad + E_Schott)/dt = 0 to LL order.
    const double schottEnergy = -electronMass * tau
        * (dot(electronAcceleration, s.electronVelocity) + dot(positronAcceleration, s.positronVelocity));
    return {s.electronPosition, s.positronPosition, s.electronDipole, s.positronDipole,
            s.time, geometry.distance, s.radiatedEnergy,
            electronKinetic + positronKinetic + coulombPotential,
            schottEnergy,
            electronKinetic + 0.5 * coulombPotential,
            positronKinetic + 0.5 * coulombPotential};
}

const char* phenomenonName(Phenomenon phenomenon) {
    switch (phenomenon) {
        case Phenomenon::DirectCollision: return "Direct collision";
        case Phenomenon::Scattering: return "Scattering";
        case Phenomenon::ParaPositronium: return "Para-positronium";
        case Phenomenon::OrthoPositronium: return "Ortho-positronium";
    }
    return "Unknown";
}

SimulationResult simulate(std::uint64_t seed) {
    const double reducedMass = electronMass * positronMass / (electronMass + positronMass);
    const double circularSpeed = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * bohrRadius));
    const double escapeSpeed = std::sqrt(2.0) * circularSpeed;
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> unitRandom(0.0, 1.0);
    std::uniform_real_distribution<double> signedRandom(-1.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, 2.0*pi);
    std::uniform_int_distribution<int> scenario(0, 3);
    const auto randomDirection = [&]() {
        const double z = signedRandom(random);
        const double phi = azimuth(random);
        const double radial = std::sqrt(1.0 - z*z);
        return Vec3{radial*std::cos(phi), radial*std::sin(phi), z};
    };

    State s;
    s.electronPosition = {bohrRadius * positronMass / (electronMass + positronMass), 0, 0};
    s.positronPosition = {-bohrRadius * electronMass / (electronMass + positronMass), 0, 0};

    // Stratified random sampling gives every branch of the decision tree a
    // useful chance while all actual initial values remain random.
    const int sampledScenario = scenario(random);
    double radialSpeed = 0.0;
    double tangentialSpeed = 0.0;
    if (sampledScenario == 0) {
        radialSpeed = -escapeSpeed * (0.35 + 0.40 * unitRandom(random));
        tangentialSpeed = circularSpeed * (0.001 + 0.004 * unitRandom(random));
    } else if (sampledScenario == 1) {
        const double speed = escapeSpeed * (1.05 + 0.35 * unitRandom(random));
        const double tangentialFraction = 0.35 + 0.35 * unitRandom(random);
        tangentialSpeed = speed * tangentialFraction;
        radialSpeed = -speed * std::sqrt(1.0 - tangentialFraction*tangentialFraction);
    } else {
        radialSpeed = circularSpeed * (-0.12 + 0.24 * unitRandom(random));
        tangentialSpeed = circularSpeed * (0.72 + 0.25 * unitRandom(random));
    }
    const Vec3 relativeVelocity{radialSpeed, tangentialSpeed, 0.0};
    s.electronVelocity = relativeVelocity * (positronMass / (electronMass + positronMass));
    s.positronVelocity = relativeVelocity * (-electronMass / (electronMass + positronMass));

    s.electronDipole = randomDirection() * bohrMagneton;
    do {
        s.positronDipole = randomDirection() * bohrMagneton;
    } while ((sampledScenario == 2 && dot(s.electronDipole, s.positronDipole)
                                  / (bohrMagneton*bohrMagneton) < 0.5)
          || (sampledScenario == 3 && dot(s.electronDipole, s.positronDipole)
                                  / (bohrMagneton*bohrMagneton) >= 0.5));

    const double relativeEnergy = 0.5 * reducedMass * relativeVelocity.squaredNorm()
                                - coulomb * eCharge*eCharge / bohrRadius;
    const double orbitalAngularMomentum = reducedMass
                                        * cross(Vec3{bohrRadius, 0, 0}, relativeVelocity).norm();
    const double specificEnergy = relativeEnergy / reducedMass;
    const double specificAngularMomentum = orbitalAngularMomentum / reducedMass;
    const double attractionParameter = coulomb * eCharge*eCharge / reducedMass;
    const double eccentricity = std::sqrt(std::max(0.0, 1.0 + 2.0 * specificEnergy
        * specificAngularMomentum*specificAngularMomentum
        / (attractionParameter*attractionParameter)));
    const double predictedClosestApproach = specificAngularMomentum == 0.0 ? 0.0
        : specificAngularMomentum*specificAngularMomentum
          / (attractionParameter * (1.0 + eccentricity));
    const double dipoleAlignment = dot(s.electronDipole, s.positronDipole)
                                 / (bohrMagneton*bohrMagneton);
    Phenomenon phenomenon;
    if (radialSpeed < 0.0 && predictedClosestApproach < nuclearCutoff) {
        phenomenon = Phenomenon::DirectCollision;
    } else if (relativeEnergy >= 0.0) {
        phenomenon = Phenomenon::Scattering;
    } else if (dipoleAlignment >= 0.5) {
        phenomenon = Phenomenon::ParaPositronium;
    } else {
        phenomenon = Phenomenon::OrthoPositronium;
    }

    constexpr int frameCount = 1200;
    const double displayedLifetime = phenomenon == Phenomenon::DirectCollision ? 4.0e-16
                                   : phenomenon == Phenomenon::Scattering ? 2.0e-15
                                   : 1.50e-11;
    std::vector<Frame> frames;
    frames.reserve(frameCount);
    double nextFrame = 0.0;
    const double frameInterval = displayedLifetime / (frameCount - 1);

    while (frames.size() < frameCount && separation(s) > nuclearCutoff) {
        if (s.time >= nextFrame) {
            frames.push_back(makeFrame(s));
            nextFrame += frameInterval;
        }
        // Resolve each instantaneous orbit well enough to keep numerical
        // energy drift below the physical radiation loss.
        const double r = separation(s);
        const double omega = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * r*r*r));
        advance(s, std::min(2.0e-18, 2.0 * pi / (160.0 * omega)));
    }
    if (frames.empty()) throw std::runtime_error("No simulation frames were produced");
    return {std::move(frames), {relativeEnergy, orbitalAngularMomentum,
            predictedClosestApproach, dipoleAlignment, phenomenon, seed}};
}

std::string labelFor(const Frame& f) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << "t = " << f.time * 1.0e12 << " ps"
        << "     E_{e+p} = " << f.mechanicalEnergy / eCharge << " eV"
        << "     E_{rad} = " << f.radiatedEnergy / eCharge << " eV";
    return out.str();
}

std::string deltaLabelFor(const Frame& current, const Frame& initial) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << "#Deltat = " << (current.time - initial.time) * 1.0e12 << " ps"
        << "     #DeltaE_{e+p} = " << (current.mechanicalEnergy - initial.mechanicalEnergy) / eCharge << " eV"
        << "     #DeltaE_{rad} = " << (current.radiatedEnergy - initial.radiatedEnergy) / eCharge << " eV";
    return out.str();
}

std::string formatTableValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string spinLabel(const Frame& frame) {
    const char* electronArrow = frame.electronDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    const char* positronArrow = frame.positronDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    return std::string(electronArrow) + " " + positronArrow;
}

void setDipoleArrow(TPolyLine3D& shaft, TPolyLine3D& leftHead, TPolyLine3D& rightHead,
                    const Vec3& centre, const Vec3& dipole) {
    const Vec3 direction = unit(dipole);
    const Vec3 tip = centre + direction * 0.31;
    const Vec3 helper = std::abs(direction.z) < 0.8 ? Vec3{0, 0, 1} : Vec3{0, 1, 0};
    const Vec3 side = unit(cross(direction, helper));
    const Vec3 base = tip - direction * 0.085;
    shaft.SetPoint(0, centre.x, centre.y, centre.z);
    shaft.SetPoint(1, tip.x, tip.y, tip.z);
    leftHead.SetPoint(0, tip.x, tip.y, tip.z);
    leftHead.SetPoint(1, base.x + side.x*0.055, base.y + side.y*0.055, base.z + side.z*0.055);
    rightHead.SetPoint(0, tip.x, tip.y, tip.z);
    rightHead.SetPoint(1, base.x - side.x*0.055, base.y - side.y*0.055, base.z - side.z*0.055);
}
} // namespace

int main(int argc, char** argv) {
    std::random_device seedSource;
    std::uint64_t seed = (static_cast<std::uint64_t>(seedSource()) << 32) ^ seedSource();
    bool diagnose = false;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--diagnose") {
            diagnose = true;
        } else if (argument == "--seed" && i + 1 < argc) {
            seed = std::stoull(argv[++i]);
        }
    }
    SimulationResult simulation = simulate(seed);
    const std::vector<Frame>& frames = simulation.frames;
    const InitialConditions& initialConditions = simulation.initial;
    std::cout << phenomenonName(initialConditions.phenomenon) << " simulated for " << frames.back().time << " s; "
              << frames.size() << " animation frames.\n";
    if (diagnose) {
        const double initialTotal = frames.front().mechanicalEnergy + frames.front().radiatedEnergy
                                  + frames.front().schottEnergy;
        const double finalTotal = frames.back().mechanicalEnergy + frames.back().radiatedEnergy
                                + frames.back().schottEnergy;
        const double energyDrift = finalTotal - initialTotal;
        const double energyScale = std::max(std::abs(frames.back().radiatedEnergy), std::abs(initialTotal));
        const double relativeEnergyDrift = std::abs(energyDrift) / energyScale;
        const bool expectedMotion = initialConditions.phenomenon == Phenomenon::Scattering
            ? frames.back().radius > frames.front().radius
            : frames.back().radius < frames.front().radius;
        // At a direct collision the point-particle model terminates in its
        // singular regime, so only pre-contact motion is a meaningful test.
        const bool energyBalanced = initialConditions.phenomenon == Phenomenon::DirectCollision
                                  || relativeEnergyDrift < 0.05;
        const auto radiusBounds = std::minmax_element(frames.begin(), frames.end(),
            [](const Frame& a, const Frame& b) { return a.radius < b.radius; });
        std::cout << std::setprecision(8)
                  << "seed:           " << initialConditions.seed << "\n"
                  << "phenomenon:     " << phenomenonName(initialConditions.phenomenon) << "\n"
                  << "relative E:     " << initialConditions.relativeEnergy / eCharge << " eV\n"
                  << "orbital L:      " << initialConditions.orbitalAngularMomentum / hbar << " hbar\n"
                  << "predicted rmin: " << initialConditions.predictedClosestApproach * 1.0e12 << " pm\n"
                  << "initial radius: " << frames.front().radius * 1.0e12 << " pm\n"
                  << "final radius:   " << frames.back().radius * 1.0e12 << " pm\n"
                  << "radius range:   " << radiusBounds.first->radius * 1.0e12 << " .. "
                  << radiusBounds.second->radius * 1.0e12 << " pm\n"
                  << "radiated:       " << frames.back().radiatedEnergy / eCharge << " eV\n"
                  << "Schott energy:  " << frames.back().schottEnergy / eCharge << " eV\n"
                  << "energy drift:   " << energyDrift / eCharge << " eV ("
                  << relativeEnergyDrift * 100.0 << "%)\n"
                  << "diagnostic:     " << (expectedMotion && energyBalanced ? "PASS" : "FAIL") << '\n';
        return expectedMotion && energyBalanced ? 0 : 1;
    }

    TApplication app("classical_hydrogen", &argc, argv);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kBlack);
    gStyle->SetPadColor(kBlack);

    TCanvas canvas("atom", "Classical model of positronium", 1100, 820);
    canvas.SetFillColor(kBlack);
    canvas.SetSupportGL(kTRUE);
    TPad scene("scene", "Simulation", 0.0, 0.0, 1.0, 0.86);
    TPad controls("controls", "Controls", 0.0, 0.86, 1.0, 1.0);
    scene.SetFillColor(kBlack);
    controls.SetFillColor(kBlack);
    scene.Draw();
    controls.Draw();

    scene.cd();
    // ROOT rotates this 3D view when the user drags the mouse in the scene.
    // A non-zero Z span preserves perspective for the initially planar orbit.
    TView* view = TView::CreateView(1);
    const double scale = 1.0 / bohrRadius;
    double viewSpan = 1.10;
    for (const Frame& frame : frames) {
        viewSpan = std::max({viewSpan, frame.electron.norm() * scale * 1.10,
                            frame.positron.norm() * scale * 1.10});
    }
    view->SetRange(-viewSpan, -viewSpan, -0.45*viewSpan,
                   viewSpan, viewSpan, 0.45*viewSpan);
    gPad->SetTheta(70);
    gPad->SetPhi(25);

    std::vector<double> x(frames.size()), y(frames.size()), z(frames.size());
    for (size_t i = 0; i < frames.size(); ++i) {
        x[i] = frames[i].electron.x * scale;
        y[i] = frames[i].electron.y * scale;
        z[i] = frames[i].electron.z * scale;
    }
    TPolyLine3D path(static_cast<int>(frames.size()), x.data(), y.data(), z.data());
    path.SetLineColor(kOrange + 7);
    path.SetLineWidth(2);
    path.Draw();

    TPolyMarker3D electron(1), positron(1);
    electron.SetMarkerStyle(20); electron.SetMarkerSize(1.8); electron.SetMarkerColor(kAzure + 1);
    positron.SetMarkerStyle(20); positron.SetMarkerSize(2.5); positron.SetMarkerColor(kRed + 1);
    electron.Draw("same");
    positron.Draw("same");

    TPolyLine3D electronDipoleShaft(2), electronDipoleLeft(2), electronDipoleRight(2);
    TPolyLine3D positronDipoleShaft(2), positronDipoleLeft(2), positronDipoleRight(2);
    for (TPolyLine3D* arrow : {&electronDipoleShaft, &electronDipoleLeft, &electronDipoleRight}) {
        arrow->SetLineColor(kAzure + 1); arrow->SetLineWidth(3); arrow->Draw("same");
    }
    for (TPolyLine3D* arrow : {&positronDipoleShaft, &positronDipoleLeft, &positronDipoleRight}) {
        arrow->SetLineColor(kRed + 1); arrow->SetLineWidth(3); arrow->Draw("same");
    }

    const Frame& initialFrame = frames.front();
    TLatex phenomenonLabel;
    phenomenonLabel.SetNDC(); phenomenonLabel.SetTextColor(kOrange + 7);
    phenomenonLabel.SetTextSize(0.034); phenomenonLabel.SetTextFont(62);
    phenomenonLabel.SetText(0.03, 0.185,
        (std::string("Phenomenon: ") + phenomenonName(initialConditions.phenomenon)).c_str());
    phenomenonLabel.Draw();
    std::ostringstream conditionsText;
    conditionsText << std::fixed << std::setprecision(3)
                   << "E_{rel} = " << initialConditions.relativeEnergy / eCharge << " eV"
                   << "     L_{orb} = " << initialConditions.orbitalAngularMomentum / hbar << " #hbar"
                   << "     r_{min} = " << initialConditions.predictedClosestApproach * 1.0e12 << " pm";
    TLatex conditionsLabel;
    conditionsLabel.SetNDC(); conditionsLabel.SetTextColor(kWhite);
    conditionsLabel.SetTextSize(0.024); conditionsLabel.SetTextFont(42);
    conditionsLabel.SetText(0.03, 0.150, conditionsText.str().c_str());
    conditionsLabel.Draw();
    constexpr double bottomHeaderY = 0.125;
    constexpr double bottomInitialY = 0.090;
    constexpr double bottomCurrentY = 0.055;
    constexpr double bottomDeltaY = 0.020;
    const std::array<double, 7> bottomXs = {0.02, 0.13, 0.23, 0.33, 0.45, 0.57, 0.69};
    const std::array<const char*, 8> bottomHeaders = {
        "stan", "t [ps]", "Spin (e^{-},e^{+})", "r [pm]", "E_{e} [eV]", "E_{p} [eV]", "E_{sum} [eV]", "E_{rad} [eV]"
    };
    std::array<TLatex, 7> bottomHeaderLabels;
    for (size_t i = 0; i < bottomHeaderLabels.size(); ++i) {
        bottomHeaderLabels[i].SetNDC(); bottomHeaderLabels[i].SetTextColor(kWhite);
        bottomHeaderLabels[i].SetTextSize(0.024); bottomHeaderLabels[i].SetTextFont(62);
        bottomHeaderLabels[i].SetText(bottomXs[i], bottomHeaderY, bottomHeaders[i]);
        bottomHeaderLabels[i].Draw();
    }
    TLatex bottomRadHeader;
    bottomRadHeader.SetNDC(); bottomRadHeader.SetTextColor(kWhite);
    bottomRadHeader.SetTextSize(0.024); bottomRadHeader.SetTextFont(62);
    bottomRadHeader.SetText(0.84, bottomHeaderY, bottomHeaders.back());
    bottomRadHeader.Draw();
    std::array<TLatex, 7> initialBottomRow;
    std::array<TLatex, 7> currentBottomRow;
    std::array<TLatex, 7> deltaBottomRow;
    TLatex initialBottomRad;
    TLatex currentBottomRad;
    TLatex deltaBottomRad;
    const auto initializeBottomRow = [&](std::array<TLatex, 7>& row, int color, double textSize) {
        for (TLatex& cell : row) {
            cell.SetNDC(); cell.SetTextColor(color);
            cell.SetTextSize(textSize); cell.SetTextFont(62);
            cell.SetText(0, 0, "");
            cell.Draw();
        }
    };
    initializeBottomRow(initialBottomRow, kCyan + 1, 0.024);
    initializeBottomRow(currentBottomRow, kYellow + 1, 0.024);
    initializeBottomRow(deltaBottomRow, kGreen + 2, 0.024);
    initialBottomRad.SetNDC(); initialBottomRad.SetTextColor(kCyan + 1);
    initialBottomRad.SetTextSize(0.024); initialBottomRad.SetTextFont(62);
    currentBottomRad.SetNDC(); currentBottomRad.SetTextColor(kYellow + 1);
    currentBottomRad.SetTextSize(0.024); currentBottomRad.SetTextFont(62);
    deltaBottomRad.SetNDC(); deltaBottomRad.SetTextColor(kGreen + 2);
    deltaBottomRad.SetTextSize(0.024); deltaBottomRad.SetTextFont(62);

    const auto drawBottomRow = [&](std::array<TLatex, 7>& row, double y, const std::array<std::string, 7>& values) {
        for (size_t i = 0; i < row.size(); ++i) {
            row[i].SetText(bottomXs[i], y, values[i].c_str());
        }
    };

    drawBottomRow(initialBottomRow, bottomInitialY, std::array<std::string, 7>{
        "initial",
        formatTableValue(initialFrame.time * 1.0e12),
        spinLabel(initialFrame),
        formatTableValue(initialFrame.radius * 1.0e12),
        formatTableValue(initialFrame.electronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.positronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.mechanicalEnergy / eCharge)
    });
    initialBottomRad.SetText(0.84, bottomInitialY, formatTableValue(initialFrame.radiatedEnergy / eCharge).c_str());
    initialBottomRad.Draw();
    drawBottomRow(currentBottomRow, bottomCurrentY, std::array<std::string, 7>{"current", "0.00", spinLabel(initialFrame), "0.00", "0.00", "0.00", "0.00"});
    currentBottomRad.SetText(0.84, bottomCurrentY, "0.00");
    currentBottomRad.Draw();
    drawBottomRow(deltaBottomRow, bottomDeltaY, std::array<std::string, 7>{"delta", "0.00", "--", "0.00", "0.00", "0.00", "0.00"});
    deltaBottomRad.SetText(0.84, bottomDeltaY, "0.00");
    deltaBottomRad.Draw();

    controls.cd();
    gInterpreter->Declare("void ToggleSimulation(); void ExitSimulation();");
    TButton stopButton("STOP", "ToggleSimulation();", 0.73, 0.25, 0.85, 0.78);
    stopButton.SetFillColor(kOrange + 7);
    stopButton.SetTextFont(62);
    stopButton.Draw();
    gStopButton = &stopButton;
    TButton exitButton("EXIT", "ExitSimulation();", 0.87, 0.25, 0.97, 0.78);
    exitButton.SetFillColor(kRed + 1);
    exitButton.SetTextFont(62);
    exitButton.Draw();
    const auto syncStopButton = [&]() {
        if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
        controls.Modified();
    };
    syncStopButton();
    controls.Modified();
    controls.Update();

    const auto updateBottomRow = [&](const Frame& f) {
        const std::array<std::string, 7> currentValues = {
            "current",
            formatTableValue(f.time * 1.0e12),
            spinLabel(f),
            formatTableValue(f.radius * 1.0e12),
            formatTableValue(f.electronMechanicalEnergy / eCharge),
            formatTableValue(f.positronMechanicalEnergy / eCharge),
            formatTableValue(f.mechanicalEnergy / eCharge)
        };
        drawBottomRow(currentBottomRow, bottomCurrentY, currentValues);
        currentBottomRad.SetText(0.84, bottomCurrentY, formatTableValue(f.radiatedEnergy / eCharge).c_str());

        const std::array<std::string, 7> deltaValues = {
            "delta",
            formatTableValue((f.time - initialFrame.time) * 1.0e12),
            "--",
            formatTableValue((f.radius - initialFrame.radius) * 1.0e12),
            formatTableValue((f.electronMechanicalEnergy - initialFrame.electronMechanicalEnergy) / eCharge),
            formatTableValue((f.positronMechanicalEnergy - initialFrame.positronMechanicalEnergy) / eCharge),
            formatTableValue((f.mechanicalEnergy - initialFrame.mechanicalEnergy) / eCharge)
        };
        drawBottomRow(deltaBottomRow, bottomDeltaY, deltaValues);
        deltaBottomRad.SetText(0.84, bottomDeltaY,
                               formatTableValue((f.radiatedEnergy - initialFrame.radiatedEnergy) / eCharge).c_str());
    };

    updateBottomRow(initialFrame);

    constexpr size_t renderStride = 3;
    constexpr unsigned int renderDelayMilliseconds = 16;
    for (size_t i = 0; i < frames.size(); ++i) {
        while (gSimulationPaused && !gExitRequested) {
            syncStopButton();
            controls.Update();
            gSystem->ProcessEvents();
            syncStopButton();
            gSystem->Sleep(20);
        }
        if (gExitRequested) return 0;
        if (i % renderStride != 0 && i + 1 != frames.size()) continue;
        syncStopButton();
        const Frame& f = frames[i];
        electron.SetPoint(0, f.electron.x * scale, f.electron.y * scale, f.electron.z * scale);
        positron.SetPoint(0, f.positron.x * scale, f.positron.y * scale, f.positron.z * scale);
        const Vec3 electronPosition = f.electron * scale;
        const Vec3 positronPosition = f.positron * scale;
        setDipoleArrow(electronDipoleShaft, electronDipoleLeft, electronDipoleRight,
                       electronPosition, f.electronDipole);
        setDipoleArrow(positronDipoleShaft, positronDipoleLeft, positronDipoleRight,
                       positronPosition, f.positronDipole);
        updateBottomRow(f);
        scene.Modified();
        canvas.Modified();
        canvas.Update();
        gSystem->ProcessEvents();
        syncStopButton();
        gSystem->Sleep(renderDelayMilliseconds);
    }
    controls.cd();
    canvas.Modified();
    canvas.Update();
    app.Run();
    return 0;
}
