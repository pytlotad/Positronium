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
    Vec3 electronExternalForce, positronExternalForce;
    Vec3 electronDipole, positronDipole; // classical magnetic moments, in J/T
    double time = 0;
    double radiatedEnergy = 0;
};

struct Frame {
    Vec3 electron, positron, electronDipole, positronDipole;
    double time, radius, radiatedEnergy, mechanicalEnergy;
    double electronMechanicalEnergy, positronMechanicalEnergy;
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
        retardedTime = newerSourceState.time - (observationPosition - source.position).norm() / c;
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

// In the weak-self-field regime this is the local, order-reduced
// Landau-Lifshitz force.  The derivative of the external force is evaluated
// from consecutive integration steps, avoiding the unstable third derivative
// in the Lorentz-Abraham-Dirac equation.
Vec3 landauLifshitzForce(const Vec3& externalForce, const Vec3& previousExternalForce,
                         double mass, double dt, bool hasHistory) {
    if (!hasHistory) return {};
    const double tau = eCharge*eCharge / (6.0 * pi * epsilon0 * mass * c*c*c);
    return (externalForce - previousExternalForce) * (tau / dt);
}

struct MutualForces { Vec3 electron, positron; };

MutualForces mutualForces(const State& history, const State& s) {
    const ElectromagneticField positronField = lienardWiechertField(s.electronPosition, history, s, false, eCharge);
    const ElectromagneticField electronField = lienardWiechertField(s.positronPosition, history, s, true, -eCharge);
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 dipoleForceOnElectron = dipoleForce(geometry.electronMinusPositron, geometry,
                                                    s.electronDipole, s.positronDipole);
    return {lorentzForce(-eCharge, s.electronVelocity, positronField) + dipoleForceOnElectron,
            lorentzForce(eCharge, s.positronVelocity, electronField) - dipoleForceOnElectron};
}

// Relativistic momentum update driven by mutually retarded fields.
void advance(State& s, const State& history, double dt, bool hasHistory) {
    MutualForces forces = mutualForces(history, s);
    const Vec3 electronReaction = landauLifshitzForce(forces.electron, s.electronExternalForce,
                                                      electronMass, dt, hasHistory);
    const Vec3 positronReaction = landauLifshitzForce(forces.positron, s.positronExternalForce,
                                                      positronMass, dt, hasHistory);
    Vec3 electronMomentum = momentum(s.electronVelocity, electronMass) + (forces.electron + electronReaction) * (0.5 * dt);
    Vec3 positronMomentum = momentum(s.positronVelocity, positronMass) + (forces.positron + positronReaction) * (0.5 * dt);
    State trial = s;
    trial.time += dt;
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    trial.electronPosition += trial.electronVelocity * dt;
    trial.positronPosition += trial.positronVelocity * dt;
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, forces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, forces.positron, positronMass);

    const MutualForces trialForces = mutualForces(s, trial);
    const Vec3 trialElectronReaction = landauLifshitzForce(trialForces.electron, forces.electron,
                                                           electronMass, dt, true);
    const Vec3 trialPositronReaction = landauLifshitzForce(trialForces.positron, forces.positron,
                                                           positronMass, dt, true);
    electronMomentum += (trialForces.electron + trialElectronReaction) * (0.5 * dt);
    positronMomentum += (trialForces.positron + trialPositronReaction) * (0.5 * dt);
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, trialForces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, trialForces.positron, positronMass);
    trial.electronExternalForce = trialForces.electron;
    trial.positronExternalForce = trialForces.positron;

    PairGeometry geometry = pairGeometry(trial);

    const Vec3 fieldAtElectron = dipoleField(geometry.electronMinusPositron, geometry, trial.positronDipole);
    const Vec3 fieldAtPositron = dipoleField(geometry.electronMinusPositron * -1.0, geometry, trial.electronDipole);
    // The real dipole precession is extremely slow on an orbital timescale.
    // This factor advances only the displayed spin dynamics so its changing
    // 3D direction is observable; the mechanical dipole force above remains
    // at its physical, unamplified value.
    constexpr double visibleSpinTimeScale = 3.0e4;
    precessDipole(trial.electronDipole, fieldAtElectron, -1.76085963023e11, dt * visibleSpinTimeScale);
    precessDipole(trial.positronDipole, fieldAtPositron, 1.76085963023e11, dt * visibleSpinTimeScale);
    trial.radiatedEnergy += (lienardPower(-eCharge, trial.electronVelocity, trial.electronAcceleration)
                           + lienardPower(eCharge, trial.positronVelocity, trial.positronAcceleration)) * dt;
    s = trial;
}

Frame makeFrame(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const double electronKinetic = (gamma(s.electronVelocity) - 1.0) * electronMass * c*c;
    const double positronKinetic = (gamma(s.positronVelocity) - 1.0) * positronMass * c*c;
    const double coulombPotential = -coulomb * eCharge*eCharge * geometry.inverseDistance;
    const double dipolePotential = -dot(s.electronDipole,
                                        dipoleField(geometry.electronMinusPositron, geometry, s.positronDipole));
    return {s.electronPosition, s.positronPosition, s.electronDipole, s.positronDipole,
            s.time, geometry.distance, s.radiatedEnergy,
            electronKinetic + positronKinetic + coulombPotential + dipolePotential,
            electronKinetic + 0.5 * (coulombPotential + dipolePotential),
            positronKinetic + 0.5 * (coulombPotential + dipolePotential)};
}

std::vector<Frame> simulate() {
    const double reducedMass = electronMass * positronMass / (electronMass + positronMass);
    const double circularSpeed = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * bohrRadius));
    State s;
    // Centre of mass is at rest; both bodies therefore move, with the positron's
    // small reflex orbit included.
    s.electronPosition = {bohrRadius * positronMass / (electronMass + positronMass), 0, 0};
    s.positronPosition = {-bohrRadius * electronMass / (electronMass + positronMass), 0, 0};
    s.electronVelocity = {0, circularSpeed * positronMass / (electronMass + positronMass), 0};
    s.positronVelocity = {0, -circularSpeed * electronMass / (electronMass + positronMass), 0};
    // Fixed seed means the initial, randomly distributed directions are
    // repeatable between runs, which makes the simulation comparable.
    std::mt19937_64 random(0x484944524f47454eULL);
    std::uniform_real_distribution<double> distributed(-1.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, 2.0*pi);
    const auto randomDirection = [&]() {
        const double z = distributed(random);
        const double phi = azimuth(random);
        const double radial = std::sqrt(1.0 - z*z);
        return Vec3{radial*std::cos(phi), radial*std::sin(phi), z};
    };
    s.electronDipole = randomDirection() * bohrMagneton;
    s.positronDipole = randomDirection() * bohrMagneton;
    State history = s;
    bool hasHistory = false;

    constexpr int frameCount = 5000;
    constexpr double displayedLifetime = 1.50e-10;
    std::vector<Frame> frames;
    frames.reserve(frameCount);
    double nextFrame = 0.0;
    const double frameInterval = displayedLifetime / (frameCount - 1);

    while (frames.size() < frameCount && separation(s) > nuclearCutoff) {
        if (s.time >= nextFrame) {
            frames.push_back(makeFrame(s));
            nextFrame += frameInterval;
        }
        // At least 80 steps per instantaneous orbit.  This shortens with the
        // radius and keeps the Coulomb integration stable through the plunge.
        const double r = separation(s);
        const double omega = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * r*r*r));
        const State current = s;
        advance(s, history, std::min(2.0e-18, 2.0 * pi / (80.0 * omega)), hasHistory);
        history = current;
        hasHistory = true;
    }
    if (frames.empty()) throw std::runtime_error("No simulation frames were produced");
    return frames;
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
    std::vector<Frame> frames = simulate();
    std::cout << "Classical radiative collapse simulated for " << frames.back().time << " s; "
              << frames.size() << " animation frames.\n";

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
    view->SetRange(-1.10, -1.10, -0.50, 1.10, 1.10, 0.50);
    gPad->SetTheta(70);
    gPad->SetPhi(25);

    const double scale = 1.0 / bohrRadius;
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
    constexpr double bottomHeaderY = 0.125;
    constexpr double bottomInitialY = 0.090;
    constexpr double bottomCurrentY = 0.055;
    constexpr double bottomDeltaY = 0.020;
    const std::array<double, 6> bottomXs = {0.03, 0.16, 0.29, 0.41, 0.53, 0.66};
    const std::array<const char*, 7> bottomHeaders = {
        "stan", "t [ps]", "r [pm]", "E_{e} [eV]", "E_{p} [eV]", "E_{sum} [eV]", "E_{rad} [eV]"
    };
    std::array<TLatex, 6> bottomHeaderLabels;
    for (size_t i = 0; i < bottomHeaderLabels.size(); ++i) {
        bottomHeaderLabels[i].SetNDC(); bottomHeaderLabels[i].SetTextColor(kWhite);
        bottomHeaderLabels[i].SetTextSize(0.024); bottomHeaderLabels[i].SetTextFont(62);
        bottomHeaderLabels[i].SetText(bottomXs[i], bottomHeaderY, bottomHeaders[i]);
        bottomHeaderLabels[i].Draw();
    }
    TLatex bottomRadHeader;
    bottomRadHeader.SetNDC(); bottomRadHeader.SetTextColor(kWhite);
    bottomRadHeader.SetTextSize(0.024); bottomRadHeader.SetTextFont(62);
    bottomRadHeader.SetText(0.80, bottomHeaderY, bottomHeaders.back());
    bottomRadHeader.Draw();
    std::array<TLatex, 6> initialBottomRow;
    std::array<TLatex, 6> currentBottomRow;
    std::array<TLatex, 6> deltaBottomRow;
    TLatex initialBottomRad;
    TLatex currentBottomRad;
    TLatex deltaBottomRad;
    const auto initializeBottomRow = [&](std::array<TLatex, 6>& row, int color, double textSize) {
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

    const auto drawBottomRow = [&](std::array<TLatex, 6>& row, double y, const std::array<std::string, 6>& values) {
        for (size_t i = 0; i < row.size(); ++i) {
            row[i].SetText(bottomXs[i], y, values[i].c_str());
        }
    };

    drawBottomRow(initialBottomRow, bottomInitialY, std::array<std::string, 6>{
        "initial",
        formatTableValue(initialFrame.time * 1.0e12),
        formatTableValue(initialFrame.radius * 1.0e12),
        formatTableValue(initialFrame.electronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.positronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.mechanicalEnergy / eCharge)
    });
    initialBottomRad.SetText(0.80, bottomInitialY, formatTableValue(initialFrame.radiatedEnergy / eCharge).c_str());
    initialBottomRad.Draw();
    drawBottomRow(currentBottomRow, bottomCurrentY, std::array<std::string, 6>{"current", "0.00", "0.00", "0.00", "0.00", "0.00"});
    currentBottomRad.SetText(0.80, bottomCurrentY, "0.00");
    currentBottomRad.Draw();
    drawBottomRow(deltaBottomRow, bottomDeltaY, std::array<std::string, 6>{"delta", "0.00", "0.00", "0.00", "0.00", "0.00"});
    deltaBottomRad.SetText(0.80, bottomDeltaY, "0.00");
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
        const std::array<std::string, 6> currentValues = {
            "current",
            formatTableValue(f.time * 1.0e12),
            formatTableValue(f.radius * 1.0e12),
            formatTableValue(f.electronMechanicalEnergy / eCharge),
            formatTableValue(f.positronMechanicalEnergy / eCharge),
            formatTableValue(f.mechanicalEnergy / eCharge)
        };
        drawBottomRow(currentBottomRow, bottomCurrentY, currentValues);
        currentBottomRad.SetText(0.80, bottomCurrentY, formatTableValue(f.radiatedEnergy / eCharge).c_str());

        const std::array<std::string, 6> deltaValues = {
            "delta",
            formatTableValue((f.time - initialFrame.time) * 1.0e12),
            formatTableValue((f.radius - initialFrame.radius) * 1.0e12),
            formatTableValue((f.electronMechanicalEnergy - initialFrame.electronMechanicalEnergy) / eCharge),
            formatTableValue((f.positronMechanicalEnergy - initialFrame.positronMechanicalEnergy) / eCharge),
            formatTableValue((f.mechanicalEnergy - initialFrame.mechanicalEnergy) / eCharge)
        };
        drawBottomRow(deltaBottomRow, bottomDeltaY, deltaValues);
        deltaBottomRad.SetText(0.80, bottomDeltaY,
                               formatTableValue((f.radiatedEnergy - initialFrame.radiatedEnergy) / eCharge).c_str());
    };

    updateBottomRow(initialFrame);

    constexpr size_t renderStride = 12;
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
