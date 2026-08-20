// CREM attribution notice
//
// This program may be used by anyone. Publications based on results produced
// with this code, including modified versions of the code or CREM model,
// should identify the Positronium repository at
// https://github.com/pytlotad/Positronium and credit the CREM model author,
// Tadeusz Slawomir Pytlos (tadeusz.slawomir.pytlos@gmail.com).

#include <TApplication.h>
#include <TBox.h>
#include <THStack.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TInterpreter.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TView.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <exception>
#include <future>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "modules/vector3.hpp"
#include "modules/state.hpp"
#include "modules/physical_constants.hpp"
#include "modules/particle_species.hpp"
#include "modules/two_body_kinematics.hpp"
#include "modules/root_export.hpp"
#include "modules/statistics_archive.hpp"

// These functions are intentionally global: ROOT's TButton invokes its action
// through the interpreter while the animation loop observes these flags.
bool gSimulationPaused = false;
bool gExitRequested = false;
TButton* gStopButton = nullptr;
TCanvas* gVisualCanvas = nullptr;
int gVisualPhenomenon = 0;
bool gVisualExitSaveAttempted = false;

void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
void ExitSimulation() {
    gExitRequested = true;
    if(gVisualCanvas&&gVisualPhenomenon>=1&&gVisualPhenomenon<=4) {
        gVisualExitSaveAttempted=true;
        gVisualCanvas->Modified();
        gVisualCanvas->Update();
        const root_export::ExportResult screenshot=
            root_export::saveVisualScreenshot(
                *gVisualCanvas,gVisualPhenomenon);
        if(screenshot) {
            std::cout<<"Saved visual screenshot: "
                     <<screenshot.path.string()<<'\n';
        } else {
            std::cerr<<"Warning: could not save visual screenshot "
                     <<screenshot.path.string()<<": "
                     <<screenshot.error<<'\n';
        }
    }
    if (gApplication) gApplication->Terminate(0);
}

// Approximate relativistic two-body electrodynamics in SI units. Mutual
// Charge-charge forces use mutually retarded Lienard-Wiechert fields;
// dipole couplings remain low-velocity approximations and radiation reaction
// uses a local, order-reduced electric-dipole model.
namespace {
using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::StateHistory;
using positronium::objects::DipoleTensor;
using namespace positronium::parameters;
namespace two_body = positronium::kinematics;
constexpr double c=speedOfLight;
constexpr double eCharge=elementaryCharge;
constexpr double coulomb=coulombConstant;

// The two ROLES the dynamics integrates, as opposed to the species that fill
// them.  Every force law below refers to "first" and "second"; which particles
// those are is decided here and nowhere else.  Naming them electron/positron
// was accurate only for the default pair and actively misleading for any
// other, since the state fields are role slots, not species.
//
// These are variables, not constants, because --pair chooses the pair at
// startup.  Dropping constexpr costs nothing measurable: on a 30-event beam
// run it moved the wall clock by 0.07%, inside the 0.6% spread between
// repeats, because the inner loop is dominated by the retarded-history search
// and the field solver rather than by folding these scalars.  They are
// constant-initialized to the default pair, so a run that never passes --pair
// is bit-identical to the old build.
ParticlePair activePair=defaultPair;
ParticleSpecies firstSpecies=activePair.first;
ParticleSpecies secondSpecies=activePair.second;
double firstMass=firstSpecies.mass;
double secondMass=secondSpecies.mass;
double firstCharge=firstSpecies.charge;
double secondCharge=secondSpecies.charge;
double firstGFactor=firstSpecies.gFactor;
double secondGFactor=secondSpecies.gFactor;

// Gyromagnetic ratio relating the intrinsic moment to the intrinsic spin,
// mu = gamma S with gamma = g q / (2 m).  The g belongs here and was missing
// from four call sites that wrote q/(2m) instead.
//
// The model carries |mu| = (g/2) * magneton, so the correct ratio returns
// S = hbar/2 exactly, as a spin-1/2 particle must.  Dropping g returned
// 1.00116 hbar, i.e. the reported intrinsic angular momentum was inflated by
// g = 2.0023, and the dipole's response to a radiation-reaction torque,
// d(mu) = gamma tau dt, was weaker by the same factor.
double gyromagneticRatio(double charge,double mass,double gFactor) {
    return gFactor*charge/(2.0*mass);
}
double firstGyromagneticRatioOf() {
    return gyromagneticRatio(firstCharge,firstMass,firstGFactor);
}
double secondGyromagneticRatioOf() {
    return gyromagneticRatio(secondCharge,secondMass,secondGFactor);
}
// Magnitudes differ between roles for every pair except a particle and its
// own antiparticle: the proton carries 1.41e-26 J/T against the electron's
// 9.28e-24, four hundred times smaller.  Code that reached for one role's
// moment while dressing the other was correct only by that accident.
double firstMagneticMoment=magneticMoment(firstSpecies);
double secondMagneticMoment=magneticMoment(secondSpecies);

// Pair-level quantities.  These three replace the e^2 that used to be written
// out wherever two charges met, and they are NOT interchangeable:
//
//   pairChargeProduct   q1 q2, SIGNED.  Negative for an attracting pair, which
//                       is the only case with bound states.  Belongs in the
//                       Coulomb potential and force, where the sign is the
//                       physics.
//   pairCoulombStrength k|q1 q2|, a magnitude.  Belongs in the Kepler
//                       relations (period, semi-major axis, attraction
//                       parameter), which presuppose attraction.
//   pairDipoleCharge    the effective charge of the pair's electric dipole.
//                       With the centre of mass at rest, m1 r1 + m2 r2 = 0, so
//                       d = q1 r1 + q2 r2 = [(q1 m2 - q2 m1)/(m1+m2)] r.  For
//                       e+e- that collapses to -e and |d| = e|r|, which is why
//                       the old code could write e and be right; for a pair
//                       with unequal masses it does not.
double pairChargeProduct=firstCharge*secondCharge;
double pairCoulombStrength=coulomb*magnitude(pairChargeProduct);
double pairDipoleCharge=
    (firstCharge*secondMass-secondCharge*firstMass)/(firstMass+secondMass);
double pairReducedMass=firstMass*secondMass/(firstMass+secondMass);

// Point every role constant at a different pair.  Called once, from argument
// parsing, before anything integrates.
//
// The two conditions below used to be static_asserts on the hard-coded pair.
// As runtime checks they are worth strictly more: they now guard a value the
// USER supplies rather than one a developer edited into the header, and they
// are the reason --pair cannot be used to ask for a system this model has no
// business integrating.  A repelling pair has no bound states at all, so
// every experiment built around capture and inspiral would be meaningless
// rather than merely inaccurate.
void applyPair(const ParticlePair& pair) {
    if(!isAttracting(pair))
        throw std::invalid_argument(std::string("pair ")+pair.first.name+"+"
            +pair.second.name+" repels; it has no bound states");
    if(magnitude(magnitude(chargeProduct(pair))-eCharge*eCharge)
        >1.0e-12*eCharge*eCharge)
        throw std::invalid_argument(std::string("pair ")+pair.first.name+"+"
            +pair.second.name+" does not carry opposite unit charges");
    activePair=pair;
    firstSpecies=pair.first;
    secondSpecies=pair.second;
    firstMass=firstSpecies.mass;
    secondMass=secondSpecies.mass;
    firstCharge=firstSpecies.charge;
    secondCharge=secondSpecies.charge;
    firstGFactor=firstSpecies.gFactor;
    secondGFactor=secondSpecies.gFactor;
    firstMagneticMoment=magneticMoment(firstSpecies);
    secondMagneticMoment=magneticMoment(secondSpecies);
    pairChargeProduct=firstCharge*secondCharge;
    pairCoulombStrength=coulomb*magnitude(pairChargeProduct);
    pairDipoleCharge=
        (firstCharge*secondMass-secondCharge*firstMass)/(firstMass+secondMass);
    pairReducedMass=firstMass*secondMass/(firstMass+secondMass);
    // Derived in particle_species.hpp from BOTH moments and the lighter mass,
    // so it has to follow the pair rather than stay at the default's 68.47 fm.
    magneticRegularizationRadius=dipoleRegularizationRadius(pair);
    // Same reason: the boundary is a fraction of the pair's own orbit, and a
    // heavy pair starts inside hydrogen's 529 fm.
    collisionBoundaryRadius=collisionBoundaryOf(pair);
}

// Parse "first,second" as the command line spells it.
void applyPairFromOption(const std::string& value) {
    const std::size_t comma=value.find(',');
    if(comma==std::string::npos)
        throw std::invalid_argument(
            "--pair needs two species separated by a comma, e.g. proton,electron");
    const ParticleSpecies* first=speciesByName(value.substr(0,comma));
    const ParticleSpecies* second=speciesByName(value.substr(comma+1));
    if(!first||!second)
        throw std::invalid_argument("unknown species in --pair; known species: "
            +selectableSpeciesList());
    applyPair({*first,*second});
}
// Uniform external magnetic field, zero unless the run asks for one.  Zero is
// the historical behaviour and the default: the model's own list of excluded
// effects named external fields, and every committed result was produced
// without one.
//
// It enters in exactly three places, which are the three the pair can feel it
// through: the instantaneous force sum, the retarded force sum, and the local
// field each particle sees, the last of which carries it into Thomas-BMT
// precession for both roles.
Vec3 gExternalMagneticField;

// Earth's field is about this, and it is the value the startup question
// offers.  Worth knowing before reading the output: at 50 uT the cyclotron
// rate eB/m is 8.8e6 rad/s against an orbital rate near 3e15 rad/s, so the
// orbit itself is untouched at the ninth decimal.  What the field does reach
// is the dipoles, which precess at the same 8.8e6 rad/s -- about 3e-4 rad over
// a 35 ps collapse.  The effect is real, small, and mostly magnetic.
inline constexpr double earthScaleMagneticField=50.0e-6;

double dot(const Vec3& a,const Vec3& b);
double gamma(const Vec3& velocity);

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

// Classical zero-point field of stochastic electrodynamics: a random,
// homogeneous, isotropic radiation field with spectral energy density
//
//     rho(omega) = hbar omega^3 / (2 pi^2 c^3).
//
// It is NOT a drag.  Radiation reaction is the dissipative half and is already
// in the model; this is the fluctuating half, and it pumps energy INTO the
// orbit.  Together they are a Langevin pair whose stationary state is the SED
// ground state, which is where the L = hbar initial condition this project
// already uses comes from.  The omega^3 spectrum is the unique one invariant
// under boosts, which is exactly why it cannot act as a drag: a drag would
// single out a rest frame.
//
// Represented as equal-energy plane waves.  Sampling omega from the inverse
// CDF of omega^3 over the band, omega = (w1^4 + x(w2^4 - w1^4))^(1/4), makes
// every mode carry the same energy, so one amplitude serves all of them:
// u = eps0 <E^2> and <cos^2> = 1/2 give E = sqrt(2u/(N eps0)).
//
// The band cannot be the whole spectrum.  It is centred on the PAIR's own
// orbital frequency, because that is where the secular exchange happens, and
// its upper edge is what the integrator has to resolve -- the step is set by
// the trajectory error probe, so a wider band buys accuracy in the field at
// the price of a much shorter step.
struct ZeroPointField {
    // Directions and polarizations are fixed; the FREQUENCIES are not.  Each
    // mode carries a dimensionless factor f_n against the pair's osculating
    // orbital frequency, so the whole band rides up with the orbit as it
    // shrinks and stays in resonance instead of being left behind.
    std::vector<Vec3> direction, polarization, magneticDirection;
    std::vector<double> phase, frequencyFactor;
    // amplitude = coefficient * omega_orb^2, because the band energy density
    // integral of hbar omega^3 over [f_lo, f_hi]*omega_orb scales as
    // omega_orb^4 and the amplitude is its square root.  The field therefore
    // strengthens as the orbit tightens, which is what the omega^3 spectrum
    // says it should do.
    double amplitudeCoefficient=0.0;
    bool active() const {
        return amplitudeCoefficient>0.0&&!direction.empty();
    }
    void sample(const Vec3& position,double orbitalFrequency,
                double accumulatedPhase,Vec3& electric,Vec3& magnetic) const {
        electric=Vec3{};
        magnetic=Vec3{};
        if(!(orbitalFrequency>0.0)) return;
        const double amplitude=amplitudeCoefficient
            *orbitalFrequency*orbitalFrequency;
        for(std::size_t n=0;n<direction.size();++n) {
            const double waveNumber=frequencyFactor[n]*orbitalFrequency/c;
            // The temporal phase is f_n times the INTEGRATED orbital phase,
            // so its time derivative is exactly the mode's instantaneous
            // frequency and nothing jumps when the band moves.  The spatial
            // term is kept for completeness but is tiny: k*r is about 3e-3 rad
            // across the pair, which is the long-wavelength limit.
            const double wave=std::cos(
                dot(direction[n],position)*waveNumber
                -frequencyFactor[n]*accumulatedPhase+phase[n]);
            electric+=polarization[n]*(amplitude*wave);
            magnetic+=magneticDirection[n]*(amplitude*wave/c);
        }
    }
};
ZeroPointField gZeroPointField;

// scale multiplies the field AMPLITUDE, so 1 is the full zero-point level and
// the absorbed power scales as scale^2.  Anything other than 1 is a numerical
// experiment, not the physical field.
ZeroPointField makeZeroPointField(double lowFactor,double highFactor,
                                  int modeCount,double scale,
                                  std::uint64_t seed) {
    ZeroPointField field;
    if(!(scale>0.0)||modeCount<=0||!(highFactor>lowFactor)) return field;
    const double lowFourth=std::pow(lowFactor,4.0);
    const double highFourth=std::pow(highFactor,4.0);
    // u(omega_orb) = hbar (f_hi^4 - f_lo^4) omega_orb^4 / (8 pi^2 c^3), and
    // u = eps0 <E^2> with <cos^2> = 1/2 over N equal-energy modes.
    field.amplitudeCoefficient=scale*std::sqrt(
        2.0*hbar*(highFourth-lowFourth)
        /(8.0*pi*pi*c*c*c*static_cast<double>(modeCount)*epsilon0));
    std::mt19937_64 random(seed^0x5851f42d4c957f2dULL);
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    for(int mode=0;mode<modeCount;++mode) {
        const double factor=std::pow(lowFourth
            +uniform(random)*(highFourth-lowFourth),0.25);
        const double cosine=2.0*uniform(random)-1.0;
        const double azimuth=2.0*pi*uniform(random);
        const double transverse=std::sqrt(std::max(0.0,1.0-cosine*cosine));
        const Vec3 propagation{transverse*std::cos(azimuth),
                               transverse*std::sin(azimuth),cosine};
        const Vec3 seedAxis=std::abs(propagation.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
        Vec3 first=cross(propagation,seedAxis);
        first=first/first.norm();
        const Vec3 second=cross(propagation,first);
        const double polarizationAngle=2.0*pi*uniform(random);
        const Vec3 polarization=first*std::cos(polarizationAngle)
                               +second*std::sin(polarizationAngle);
        field.direction.push_back(propagation);
        field.polarization.push_back(polarization);
        field.magneticDirection.push_back(cross(propagation,polarization));
        field.phase.push_back(2.0*pi*uniform(random));
        field.frequencyFactor.push_back(factor);
    }
    return field;
}

DipoleTensor lorentzBoostDipole(const DipoleTensor& dipole,
                                const Vec3& frameVelocity) {
    const double speedSquared=frameVelocity.squaredNorm();
    if(speedSquared==0.0) return dipole;
    const double boostGamma=gamma(frameVelocity);
    const double longitudinalFactor=boostGamma*boostGamma
        /(boostGamma+1.0)/(c*c);
    return {
        (dipole.electric+cross(frameVelocity,dipole.magnetic)/(c*c))
            *boostGamma
            -frameVelocity*(longitudinalFactor
                *dot(frameVelocity,dipole.electric)),
        (dipole.magnetic-cross(frameVelocity,dipole.electric))*boostGamma
            -frameVelocity*(longitudinalFactor
                *dot(frameVelocity,dipole.magnetic))};
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Lorentz invariants of the polarization-magnetization tensor.  Only the
// covariance tests consume them, so they are not compiled into the production
// binary at all.
double dipoleFirstInvariant(const DipoleTensor& dipole) {
    return dipole.magnetic.squaredNorm()
        -c*c*dipole.electric.squaredNorm();
}

double dipoleSecondInvariant(const DipoleTensor& dipole) {
    return c*dot(dipole.electric,dipole.magnetic);
}
#endif
Vec3 unit(const Vec3& v) { return v / v.norm(); }
struct FourVector { double time=0.0; Vec3 space; }; // x^0=ct convention
struct ElectromagneticField { Vec3 electric, magnetic; };
double minkowskiDot(const FourVector& first,const FourVector& second) {
    return first.time*second.time-dot(first.space,second.space);
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
#include "modules/maxwell_validation_backend.hpp"
#endif

bool isFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

bool isFinite(const State& state) {
    return isFinite(state.firstPosition) && isFinite(state.secondPosition)
        && isFinite(state.firstVelocity) && isFinite(state.secondVelocity)
        && two_body::subluminal(state.firstVelocity)
        && two_body::subluminal(state.secondVelocity)
        && isFinite(state.firstAcceleration) && isFinite(state.secondAcceleration)
        && isFinite(state.firstDipole) && isFinite(state.secondDipole)
        &&isFinite(state.firstElectricDipole)
        &&isFinite(state.secondElectricDipole)
        &&isFinite(state.firstProperDipole)
        &&isFinite(state.secondProperDipole)
        && isFinite(state.radiatedMomentum)
        && isFinite(state.radiatedAngularMomentum)
        && isFinite(state.boundFieldMomentum)
        && isFinite(state.boundFieldAngularMomentum)
        && isFinite(state.reactionMomentumMismatch)
        && isFinite(state.reactionAngularMomentumMismatch)
        && std::isfinite(state.time) && std::isfinite(state.radiatedEnergy)
        && std::isfinite(state.dipoleConstraintEnergy)
        && std::isfinite(state.zeroPointPhase)
        && std::isfinite(state.boundFieldEnergy)
        && std::isfinite(state.reactionEnergyMismatch);
}

Vec3 interpolateVector(const Vec3& before, const Vec3& after, double fraction) {
    return before + (after - before) * fraction;
}

Vec3 interpolateDipole(const Vec3& before, const Vec3& after, double fraction) {
    Vec3 result = interpolateVector(before, after, fraction);
    const double targetNorm = before.norm() + (after.norm() - before.norm()) * fraction;
    const double resultNorm = result.norm();
    if (resultNorm > 0.0) result = result * (targetNorm / resultNorm);
    return result;
}

State interpolateState(const State& before, const State& after, double fraction) {
    State result;
    result.firstPosition = interpolateVector(
        before.firstPosition, after.firstPosition, fraction);
    result.secondPosition = interpolateVector(
        before.secondPosition, after.secondPosition, fraction);
    result.firstVelocity = interpolateVector(
        before.firstVelocity, after.firstVelocity, fraction);
    result.secondVelocity = interpolateVector(
        before.secondVelocity, after.secondVelocity, fraction);
    result.firstAcceleration = interpolateVector(
        before.firstAcceleration, after.firstAcceleration, fraction);
    result.secondAcceleration = interpolateVector(
        before.secondAcceleration, after.secondAcceleration, fraction);
    result.firstDipole = interpolateDipole(
        before.firstDipole, after.firstDipole, fraction);
    result.secondDipole = interpolateDipole(
        before.secondDipole, after.secondDipole, fraction);
    result.firstElectricDipole=interpolateVector(
        before.firstElectricDipole,after.firstElectricDipole,fraction);
    result.secondElectricDipole=interpolateVector(
        before.secondElectricDipole,after.secondElectricDipole,fraction);
    result.firstProperDipole=interpolateDipole(
        before.firstProperDipole,after.firstProperDipole,fraction);
    result.secondProperDipole=interpolateDipole(
        before.secondProperDipole,after.secondProperDipole,fraction);
    result.time = before.time + (after.time - before.time) * fraction;
    result.radiatedEnergy = before.radiatedEnergy
        + (after.radiatedEnergy - before.radiatedEnergy) * fraction;
    result.dipoleConstraintEnergy=before.dipoleConstraintEnergy
        +(after.dipoleConstraintEnergy-before.dipoleConstraintEnergy)*fraction;
    result.radiatedMomentum = interpolateVector(
        before.radiatedMomentum, after.radiatedMomentum, fraction);
    result.radiatedAngularMomentum = interpolateVector(
        before.radiatedAngularMomentum, after.radiatedAngularMomentum, fraction);
    result.boundFieldEnergy=before.boundFieldEnergy
        +(after.boundFieldEnergy-before.boundFieldEnergy)*fraction;
    result.boundFieldMomentum=interpolateVector(
        before.boundFieldMomentum,after.boundFieldMomentum,fraction);
    result.boundFieldAngularMomentum=interpolateVector(
        before.boundFieldAngularMomentum,after.boundFieldAngularMomentum,fraction);
    result.reactionEnergyMismatch=before.reactionEnergyMismatch
        +(after.reactionEnergyMismatch-before.reactionEnergyMismatch)*fraction;
    result.reactionMomentumMismatch=interpolateVector(
        before.reactionMomentumMismatch,after.reactionMomentumMismatch,fraction);
    result.reactionAngularMomentumMismatch=interpolateVector(
        before.reactionAngularMomentumMismatch,
        after.reactionAngularMomentumMismatch,fraction);
    return result;
}

double separationCrossingFraction(const State& before,const State& after,
                                  double targetSeparation) {
    const Vec3 initial=before.firstPosition-before.secondPosition;
    const Vec3 change=(after.firstPosition-after.secondPosition)-initial;
    const double quadratic=dot(change,change);
    const double linear=2.0*dot(initial,change);
    const double constant=dot(initial,initial)-targetSeparation*targetSeparation;
    if(quadratic<=1.0e-300) {
        if(std::abs(linear)<=1.0e-300) return 1.0;
        return std::clamp(-constant/linear,0.0,1.0);
    }
    const double discriminant=std::max(0.0,linear*linear-4.0*quadratic*constant);
    const double root1=(-linear-std::sqrt(discriminant))/(2.0*quadratic);
    const double root2=(-linear+std::sqrt(discriminant))/(2.0*quadratic);
    if(root1>=0.0&&root1<=1.0) return root1;
    if(root2>=0.0&&root2<=1.0) return root2;
    return std::clamp(root1,0.0,1.0);
}

struct Frame {
    Vec3 first, second, firstDipole, secondDipole;
    Vec3 noetherMomentum, noetherAngularMomentum;
    Vec3 radiatedMomentum, radiatedAngularMomentum;
    double canonicalMomentumScale;
    double time, radius, radiatedEnergy, mechanicalEnergy, schottEnergy;
    double firstMechanicalEnergy, secondMechanicalEnergy;
    double boundFieldEnergy=0.0,reactionEnergyMismatch=0.0;
    Vec3 boundFieldMomentum,boundFieldAngularMomentum;
    Vec3 reactionMomentumMismatch,reactionAngularMomentumMismatch;
};

enum class Phenomenon { DirectCollision, Scattering, ParaPositronium, OrthoPositronium };
enum class SimulationOutcome { ReachedCutoff, ObservationLimit, NumericalFailure };
enum class VisualStyle { Unselected, Line, Dot };

struct InitialConditions {
    double relativeEnergy;
    double orbitalAngularMomentum;
    double predictedClosestApproach;
    double dipoleAlignment;
    double timeToCutoff;
    Phenomenon phenomenon;
    std::uint64_t seed;
};

struct SimulationResult {
    std::vector<Frame> frames;
    InitialConditions initial;
    SimulationOutcome outcome;
    double minimumSeparation;
    double elapsedTime;
    double finalRadiatedEnergy;
    double maximumBeta;
};

struct SimulationOptions {
    bool collectFrames = true;
    int frameCount = 1200;
    double observationTime = 0.0; // zero selects the phenomenon's visual window
    // Positive values override the visual point-particle cutoff. Statistical
    // CREM lifetime studies use the agreed charge-cloud boundary 0.01*a0.
    double terminalSeparation = 0.0;
    std::function<void(const Frame&)> frameReady;
    // Called after every accepted integration step.  Unlike frameReady this
    // is independent of physical-time sampling and is therefore suitable for
    // keeping an interactive visualization responsive.
    std::function<void(const State&)> stepReady;
    std::function<bool()> stopRequested;
    // Set false by callers that never read radiatedEnergy / radiatedMomentum /
    // radiatedAngularMomentum / boundField* from the result.  It switches off
    // the far-zone Poynting quadrature, which is the dominant per-step cost
    // and does not influence the trajectory (see
    // ClassicalTrajectoryEngine::Accuracy::computeOutwardFlux).
    bool radiatedEnergyBookkeeping = true;
};

struct PairGeometry {
    Vec3 firstMinusSecond;
    double distance;
    double inverseDistance;
    double inverseDistanceCubed;
    double inverseDistanceFourth;
};

PairGeometry pairGeometry(const State& s) {
    const Vec3 firstMinusSecond = s.firstPosition - s.secondPosition;
    const double distanceSquared = firstMinusSecond.squaredNorm();
    const double distance = std::sqrt(distanceSquared);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;
    return {firstMinusSecond, distance, inverseDistance,
            inverseDistanceSquared * inverseDistance,
            inverseDistanceSquared * inverseDistanceSquared};
}

double separation(const State& s) { return pairGeometry(s).distance; }

double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

struct ChargeKinematics { Vec3 position, velocity, acceleration; };

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
Vec3 lerp(const Vec3& first,const Vec3& second,double fraction) {
    return first+(second-first)*fraction;
}
#endif

ChargeKinematics interpolatedCharge(const State& older, const State& newer, bool first, double time) {
    const double span = newer.time - older.time;
    if(!(span>0.0)) {
        return {first?newer.firstPosition:newer.secondPosition,
                first?newer.firstVelocity:newer.secondVelocity,
                first?newer.firstAcceleration:newer.secondAcceleration};
    }
    const double fraction=std::clamp((time-older.time)/span,0.0,1.0);
    const Vec3 oldPosition = first ? older.firstPosition : older.secondPosition;
    const Vec3 newPosition = first ? newer.firstPosition : newer.secondPosition;
    const Vec3 oldVelocity = first ? older.firstVelocity : older.secondVelocity;
    const Vec3 newVelocity = first ? newer.firstVelocity : newer.secondVelocity;
    const double s2=fraction*fraction;
    const double s3=s2*fraction;
    const double h00=2.0*s3-3.0*s2+1.0;
    const double h10=s3-2.0*s2+fraction;
    const double h01=-2.0*s3+3.0*s2;
    const double h11=s3-s2;
    const Vec3 position=oldPosition*h00+oldVelocity*(span*h10)
        +newPosition*h01+newVelocity*(span*h11);

    const double dh00=6.0*s2-6.0*fraction;
    const double dh10=3.0*s2-4.0*fraction+1.0;
    const double dh01=-dh00;
    const double dh11=3.0*s2-2.0*fraction;
    const Vec3 velocity=(oldPosition*dh00+oldVelocity*(span*dh10)
        +newPosition*dh01+newVelocity*(span*dh11))/span;

    const double d2h00=12.0*fraction-6.0;
    const double d2h10=6.0*fraction-4.0;
    const double d2h01=-d2h00;
    const double d2h11=6.0*fraction-2.0;
    const Vec3 acceleration=(oldPosition*d2h00+oldVelocity*(span*d2h10)
        +newPosition*d2h01+newVelocity*(span*d2h11))/(span*span);
    return {position,velocity,acceleration};
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
ChargeKinematics linearlyInterpolatedCharge(const State& older,
    const State& newer,bool first,double time) {
    const double span=newer.time-older.time;
    const double fraction=span>0.0
        ?std::clamp((time-older.time)/span,0.0,1.0):1.0;
    return {lerp(first?older.firstPosition:older.secondPosition,
                 first?newer.firstPosition:newer.secondPosition,fraction),
            lerp(first?older.firstVelocity:older.secondVelocity,
                 first?newer.firstVelocity:newer.secondVelocity,fraction),
            lerp(first?older.firstAcceleration:older.secondAcceleration,
                 first?newer.firstAcceleration:newer.secondAcceleration,
                 fraction)};
}
#endif

ChargeKinematics historicalCharge(const StateHistory& history,
                                   const State& present, bool first,
                                   double time) {
    const State& earliest = history.empty() ? present : history.front();
    if (time <= earliest.time) {
        const double delta = time - earliest.time;
        const Vec3 position = first ? earliest.firstPosition : earliest.secondPosition;
        const Vec3 velocity = first ? earliest.firstVelocity : earliest.secondVelocity;
        const Vec3 acceleration = first ? earliest.firstAcceleration
                                           : earliest.secondAcceleration;
        return {position + velocity * delta + acceleration * (0.5 * delta * delta),
                velocity + acceleration * delta, acceleration};
    }

    // Binary search, not a linear scan.  This lookup runs inside the Newton
    // iteration of every retarded-field evaluation, and the history reaches
    // thousands of entries when a trajectory turns around at short range: the
    // retention window shrinks like r while the step shrinks like r^{3/2}, so
    // the node count grows as the pair approaches.  historicalState() already
    // searched this way; this was the one place that did not.
    const auto newer = std::lower_bound(history.begin(), history.end(), time,
        [](const State& sample, double requested) {
            return sample.time < requested;
        });
    if (newer == history.end()) {
        const State& latest = history.back();
        if (present.time > latest.time) {
            return interpolatedCharge(latest, present, first, time);
        }
        return interpolatedCharge(latest, latest, first, time);
    }
    if (newer == history.begin()) {
        return interpolatedCharge(*newer, *newer, first, time);
    }
    return interpolatedCharge(*std::prev(newer), *newer, first, time);
}

// Mutual, retarded Lienard-Wiechert field of a moving point charge.
ElectromagneticField lienardWiechertField(const Vec3& observationPosition,
                                          double observationTime,
                                          const StateHistory& history,
                                          const State& presentState,
                                          bool sourceIsFirst,
                                          double sourceCharge,
                                          double regularizationRadius=0.0) {
    ChargeKinematics source = historicalCharge(
        history, presentState, sourceIsFirst, observationTime);
    double retardedTime = observationTime
                        - (observationPosition - source.position).norm() / c;
    for (int iteration = 0; iteration < 16; ++iteration) {
        source = historicalCharge(
            history, presentState, sourceIsFirst, retardedTime);
        const Vec3 retardedDisplacement=observationPosition-source.position;
        const double retardedDistance=retardedDisplacement.norm();
        const Vec3 retardedDirection=retardedDisplacement/retardedDistance;
        const double lightConeResidual=retardedTime+retardedDistance/c
                                      -observationTime;
        const double lightConeDerivative=std::max(1.0e-8,
            1.0-dot(retardedDirection,source.velocity/c));
        const double refinedTime=retardedTime
                               -lightConeResidual/lightConeDerivative;
        if (std::abs(refinedTime-retardedTime)
            <=1.0e-30+1.0e-14*std::abs(retardedTime)) {
            retardedTime = refinedTime;
            break;
        }
        retardedTime = refinedTime;
    }
    const Vec3 displacement = observationPosition - source.position;
    const double distance = displacement.norm();
    if(distance<=std::numeric_limits<double>::min()) return {};
    const Vec3 direction = displacement / distance;
    // Trial stages may cross the declared boundary before the enclosing event
    // locator clips the trajectory.  Never evaluate the singular point-charge
    // formula inside a domain whose result is discarded by the model.
    const double fieldDistance=std::max(distance,nuclearCutoff);
    const Vec3 beta = source.velocity / c;
    const double betaSquared = beta.squaredNorm();
    const double kappa = std::max(1.0e-8, 1.0 - dot(direction, beta));
    const Vec3 velocityField = (direction - beta) * ((1.0 - betaSquared) /
                              (kappa*kappa*kappa * fieldDistance*fieldDistance));
    const Vec3 accelerationField = cross(direction, cross(direction - beta, source.acceleration)) /
                                   (c*c * kappa*kappa*kappa * fieldDistance);
    double formFactor=1.0;
    if(regularizationRadius>0.0) {
        const double u=distance/(std::sqrt(2.0)*regularizationRadius);
        formFactor=std::erf(u)-2.0*u*std::exp(-u*u)/std::sqrt(pi);
    }
    const Vec3 electric = (velocityField + accelerationField)
                        * (coulomb * sourceCharge*formFactor);
    return {electric, cross(direction, electric) / c};
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
Vec3 regularizedDipoleVectorPotential(const Vec3& observationPosition,
                                      const Vec3& sourcePosition,
                                      const Vec3& dipole,double radius) {
    const Vec3 displacement=observationPosition-sourcePosition;
    const double distance=displacement.norm();
    if(distance<=std::numeric_limits<double>::min()||dipole.squaredNorm()==0.0)
        return {};
    const double u=distance/(std::sqrt(2.0)*radius);
    const double formFactor=std::erf(u)-2.0*u*std::exp(-u*u)/std::sqrt(pi);
    return cross(dipole,displacement)
        *(mu0/(4.0*pi)*formFactor/(distance*distance*distance));
}

void initializeRetardedPairFields(MaxwellBlock& block,const State& state,
                                  const StateHistory& history,
                                  int projectionIterations=400) {
    const RelativisticChargeCloud first{firstCharge,chargeCloudRestRadius};
    const RelativisticChargeCloud second{secondCharge,chargeCloudRestRadius};
    block.clearSources();
    block.depositCloud(first,state.firstPosition,state.firstVelocity);
    block.depositCloud(second,state.secondPosition,state.secondVelocity);
    block.depositCovariantDipole(state.firstPosition,state.firstVelocity,
                                 state.firstDipole);
    block.depositCovariantDipole(state.secondPosition,state.secondVelocity,
                                 state.secondDipole);
    block.finalizeBoundInstantaneous();
    block.clearFields();
    std::vector<Vec3> dipoleVectorPotential(block.cells().size());
    for(int k=0;k<block.cellsPerAxis();++k)
        for(int j=0;j<block.cellsPerAxis();++j)
            for(int i=0;i<block.cellsPerAxis();++i) {
                const Vec3 position=block.cellPosition(i,j,k);
                const ElectromagneticField fromFirst=lienardWiechertField(
                    position,state.time,history,state,true,firstCharge,
                    chargeCloudRestRadius);
                const ElectromagneticField fromSecond=lienardWiechertField(
                    position,state.time,history,state,false,secondCharge,
                    chargeCloudRestRadius);
                const Vec3 dipolePotential=
                    regularizedDipoleVectorPotential(position,state.firstPosition,
                        state.firstDipole,chargeCloudRestRadius)
                   +regularizedDipoleVectorPotential(position,state.secondPosition,
                        state.secondDipole,chargeCloudRestRadius);
                dipoleVectorPotential[(static_cast<std::size_t>(k)
                    *block.cellsPerAxis()+j)*block.cellsPerAxis()+i]=dipolePotential;
                block.replaceFieldAt(position,
                    fromFirst.electric+fromSecond.electric,
                    fromFirst.magnetic+fromSecond.magnetic);
            }
    block.addMagneticCurl(dipoleVectorPotential);
    // Preserve the retarded transverse content while matching the discrete
    // extended sources and eliminating grid-level magnetic monopoles.
    block.projectElectricGaussConstraint(projectionIterations);
    block.projectMagneticDivergenceConstraint(projectionIterations);
}
#endif

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Rodrigues precession used only by the grid-coupled pusher in the field
// validation build; production dipole transport goes through the covariant
// BMT integrator in advanceCovariantBmt().
void precessDipole(Vec3& dipole,const Vec3& field,
                   double gyromagneticRatio,double dt) {
    const double dipoleMagnitude=dipole.norm();
    if(dipoleMagnitude==0.0) return;
    const double fieldMagnitude = field.norm();
    if (fieldMagnitude == 0.0) return;
    const Vec3 axis = field / fieldMagnitude;
    // d(mu)/dt = gamma mu x B. rotatedAround uses B x mu for a positive
    // angle, hence the minus sign.
    const double angle = -gyromagneticRatio * fieldMagnitude * dt;
    // Rodrigues rotation preserves the magnitude of the magnetic moment.
    const Vec3 rotated = dipole * std::cos(angle) + cross(axis, dipole) * std::sin(angle)
                       + axis * ((axis.x*dipole.x + axis.y*dipole.y + axis.z*dipole.z) * (1.0 - std::cos(angle)));
    dipole=rotated*(dipoleMagnitude/rotated.norm());
}
#endif

double gamma(const Vec3& velocity) {
    return two_body::strictLorentzFactor(velocity);
}

Vec3 momentum(const Vec3& velocity, double mass) { return velocity * (gamma(velocity) * mass); }

double kineticEnergy(const Vec3& velocity,double mass) {
    return two_body::kineticEnergyFromVelocity(velocity,mass);
}

Vec3 velocityFromMomentum(const Vec3& momentum, double mass) {
    const double gammaFromMomentum = std::sqrt(1.0 + momentum.squaredNorm() / (mass*mass*c*c));
    return momentum / (gammaFromMomentum * mass);
}

void synchronizeCovariantDipoles(State& state) {
    if(state.firstProperDipole.squaredNorm()==0.0
        &&state.firstDipole.squaredNorm()>0.0)
        state.firstProperDipole=state.firstDipole;
    if(state.secondProperDipole.squaredNorm()==0.0
        &&state.secondDipole.squaredNorm()>0.0)
        state.secondProperDipole=state.secondDipole;
    const DipoleTensor firstLab=lorentzBoostDipole(
        {{},state.firstProperDipole},state.firstVelocity);
    const DipoleTensor secondLab=lorentzBoostDipole(
        {{},state.secondProperDipole},state.secondVelocity);
    state.firstElectricDipole=firstLab.electric;
    state.firstDipole=firstLab.magnetic;
    state.secondElectricDipole=secondLab.electric;
    state.secondDipole=secondLab.magnetic;
}
Vec3 thomasBmtEffectiveField(const Vec3& velocity,
                             const ElectromagneticField& field,
                             double gFactor);

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Relativistic Boris rotation written in momentum variables.  Electric
// half-kicks surround a magnetic rotation evaluated with the intermediate
// Lorentz factor; the map is time reversible and cannot produce |v|>=c.
Vec3 relativisticBorisPush(const Vec3& momentumBefore, double charge,
                           double mass, const Vec3& electric,
                           const Vec3& magnetic, double dt) {
    const Vec3 momentumMinus=momentumBefore+electric*(0.5*charge*dt);
    const double gammaMinus=std::sqrt(1.0
        +momentumMinus.squaredNorm()/(mass*mass*c*c));
    const Vec3 t=magnetic*(charge*dt/(2.0*gammaMinus*mass));
    const Vec3 s=t*(2.0/(1.0+t.squaredNorm()));
    const Vec3 momentumPrime=momentumMinus+cross(momentumMinus,t);
    const Vec3 momentumPlus=momentumMinus+cross(momentumPrime,s);
    return momentumPlus+electric*(0.5*charge*dt);
}

// Production-order particle/Yee coupling. The gathered staggered fields push
// both particles first; their complete old-to-new trajectories then deposit
// rho^{n+1} and J^{n+1/2}. Maxwell therefore receives a source satisfying its
// own discrete continuity equation and needs no Poisson repair.
void pushStateWithYeeField(State& state,YeeMaxwellBlock& field,double dt) {
    const State before=state;
    const auto [firstElectric,firstMagnetic]=
        field.interpolateField(before.firstPosition);
    const auto [secondElectric,secondMagnetic]=
        field.interpolateField(before.secondPosition);
    const Vec3 firstMomentum=relativisticBorisPush(
        momentum(before.firstVelocity,firstMass),firstCharge,firstMass,
        firstElectric,firstMagnetic,dt);
    const Vec3 secondMomentum=relativisticBorisPush(
        momentum(before.secondVelocity,secondMass),secondCharge,secondMass,
        secondElectric,secondMagnetic,dt);
    state.firstVelocity=velocityFromMomentum(firstMomentum,firstMass);
    state.secondVelocity=velocityFromMomentum(secondMomentum,secondMass);
    state.firstPosition=before.firstPosition+state.firstVelocity*dt;
    state.secondPosition=before.secondPosition+state.secondVelocity*dt;
    state.firstAcceleration=(state.firstVelocity-before.firstVelocity)/dt;
    state.secondAcceleration=(state.secondVelocity-before.secondVelocity)/dt;
    field.clearSources();
    field.depositGaussianEsirkepov(firstCharge,chargeCloudRestRadius,
        before.firstPosition,before.firstVelocity,state.firstPosition,
        state.firstVelocity,dt);
    field.depositGaussianEsirkepov(secondCharge,chargeCloudRestRadius,
        before.secondPosition,before.secondVelocity,state.secondPosition,
        state.secondVelocity,dt);
    field.depositCovariantDipoleYee(state.firstPosition,
        state.firstVelocity,state.firstDipole,chargeCloudRestRadius);
    field.depositCovariantDipoleYee(state.secondPosition,
        state.secondVelocity,state.secondDipole,chargeCloudRestRadius);
    field.finalizeBoundSources(dt);
    field.advance(dt);
    state.time+=dt;
}

void pushStateWithGridField(State& state, const MaxwellBlock& field,
                            double dt, bool reciprocalDipoles=false,
                            const State* samplingState=nullptr,
                            const ElectromagneticField& firstSelfField={},
                            const ElectromagneticField& secondSelfField={},
                            const DynamicSelfFieldCalibration* firstCalibration=nullptr,
                            const DynamicSelfFieldCalibration* secondCalibration=nullptr) {
    const State& sample=samplingState?*samplingState:state;
    auto [firstElectric,firstMagnetic]=
        field.interpolateField(sample.firstPosition);
    auto [secondElectric,secondMagnetic]=
        field.interpolateField(sample.secondPosition);
    const ElectromagneticField currentFirstSelf=firstCalibration
        ?firstCalibration->field(field,sample.firstPosition,
                                    sample.firstVelocity,firstCharge):firstSelfField;
    const ElectromagneticField currentSecondSelf=secondCalibration
        ?secondCalibration->field(field,sample.secondPosition,
                                    sample.secondVelocity,secondCharge):secondSelfField;
    firstElectric=firstElectric-currentFirstSelf.electric;
    firstMagnetic=firstMagnetic-currentFirstSelf.magnetic;
    secondElectric=secondElectric-currentSecondSelf.electric;
    secondMagnetic=secondMagnetic-currentSecondSelf.magnetic;
    const GridDipoleInteraction firstDipoleInteraction=reciprocalDipoles
        ?field.dipoleInteraction(sample.firstPosition,sample.firstDipole)
        :GridDipoleInteraction{};
    const GridDipoleInteraction secondDipoleInteraction=reciprocalDipoles
        ?field.dipoleInteraction(sample.secondPosition,sample.secondDipole)
        :GridDipoleInteraction{};
    if(!reciprocalDipoles) {
        precessDipole(state.firstDipole,
            thomasBmtEffectiveField(state.firstVelocity,
                {firstElectric,firstMagnetic},firstGFactor),
            firstCharge/firstMass,0.5*dt);
        precessDipole(state.secondDipole,
            thomasBmtEffectiveField(state.secondVelocity,
                {secondElectric,secondMagnetic},secondGFactor),
            secondCharge/secondMass,0.5*dt);
    }
    const Vec3 firstMomentum=relativisticBorisPush(
        momentum(state.firstVelocity,firstMass),firstCharge,firstMass,
        firstElectric,firstMagnetic,dt)+firstDipoleInteraction.force*dt;
    const Vec3 secondMomentum=relativisticBorisPush(
        momentum(state.secondVelocity,secondMass),secondCharge,secondMass,
        secondElectric,secondMagnetic,dt)+secondDipoleInteraction.force*dt;
    state.firstVelocity=velocityFromMomentum(firstMomentum,firstMass);
    state.secondVelocity=velocityFromMomentum(secondMomentum,secondMass);
    state.firstPosition+=state.firstVelocity*dt;
    state.secondPosition+=state.secondVelocity*dt;
    if(reciprocalDipoles) {
        const double firstGyromagneticRatio=firstGyromagneticRatioOf();
        const double secondGyromagneticRatio=secondGyromagneticRatioOf();
        const double firstNorm=state.firstDipole.norm();
        const double secondNorm=state.secondDipole.norm();
        state.firstDipole+=firstDipoleInteraction.torque
                            *(firstGyromagneticRatio*dt);
        state.secondDipole+=secondDipoleInteraction.torque
                            *(secondGyromagneticRatio*dt);
        if(state.firstDipole.norm()>0.0)
            state.firstDipole=state.firstDipole*(firstNorm/state.firstDipole.norm());
        if(state.secondDipole.norm()>0.0)
            state.secondDipole=state.secondDipole*(secondNorm/state.secondDipole.norm());
        state.dipoleConstraintEnergy+=(firstDipoleInteraction.fieldPower
            +secondDipoleInteraction.fieldPower
            -dot(firstDipoleInteraction.force,state.firstVelocity)
            -dot(secondDipoleInteraction.force,state.secondVelocity))*dt;
    } else {
        precessDipole(state.firstDipole,
            thomasBmtEffectiveField(state.firstVelocity,
                {firstElectric,firstMagnetic},firstGFactor),
            firstCharge/firstMass,0.5*dt);
        precessDipole(state.secondDipole,
            thomasBmtEffectiveField(state.secondVelocity,
                {secondElectric,secondMagnetic},secondGFactor),
            secondCharge/secondMass,0.5*dt);
    }
    state.time+=dt;
}
#endif

#include "modules/electrodynamics.hpp"

// Selected once from --radiation-reaction and read by every trajectory
// constructed below (visual, beam and interaction experiments alike).
// Defaults to individualLandauLifshitz (radiation ON): the electric-dipole
// charge self-force is the only channel that removes orbital energy (see
// individualLandauLifshitzSelfForces / coherentElectricDipoleReaction), and
// estimateCremCollapse now measures the classical inspiral mechanically
// rather than assuming it, so it needs that channel switched on to observe
// anything.
// [[maybe_unused]] because the validation executable's main() never reaches
// the trajectory constructors that read this, so GCC sees no use in that
// build.
// Order of the trajectory composition: 2 is the bare symmetric step, 4 the
// Yoshida composition of three of them.  Two is the default and every result
// in this repository was produced with it.
[[maybe_unused]] int gIntegratorOrder = 2;
[[maybe_unused]] ChargeRadiationReactionModel gRadiationReactionModel =
    ChargeRadiationReactionModel::individualLandauLifshitz;

#include "modules/crem_engine.hpp"


#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Used only by the annihilation-generator unit test in the validation build:
// the production panels now plot the exact isotropic reference instead of
// fitting a sampled one.
struct LegendreFitSummary {
    double anisotropy = std::numeric_limits<double>::quiet_NaN();
    double standardError = std::numeric_limits<double>::quiet_NaN();
    std::size_t count = 0;
    bool atBoundary = false;
};

LegendreFitSummary fitSecondLegendreAnisotropy(
    const std::vector<double>& cosines) {
    LegendreFitSummary result;
    std::vector<double> secondLegendre;
    secondLegendre.reserve(cosines.size());
    for (double cosine : cosines) {
        if (!std::isfinite(cosine)) continue;
        secondLegendre.push_back(0.5*(3.0*cosine*cosine - 1.0));
    }
    result.count = secondLegendre.size();
    if (result.count == 0) return result;

    // The event-level likelihood is proportional to
    // product_i [1 + a2 P2(cos(theta_i))], with the physical range
    // -1 <= a2 <= 2.  Its derivative is monotone, so bisection gives the
    // constrained unbinned MLE without treating the two back-to-back photons
    // as independent observations.
    constexpr double lower = -1.0;
    constexpr double upper = 2.0;
    constexpr double guard = 64.0*std::numeric_limits<double>::epsilon();
    const auto derivative = [&](double anisotropy) {
        double value = 0.0;
        for (double legendre : secondLegendre) {
            value += legendre/(1.0 + anisotropy*legendre);
        }
        return value;
    };
    const double safeLower = lower + guard;
    const double safeUpper = upper - guard;
    if (derivative(safeLower) <= 0.0) {
        result.anisotropy = lower;
        result.atBoundary = true;
    } else if (derivative(safeUpper) >= 0.0) {
        result.anisotropy = upper;
        result.atBoundary = true;
    } else {
        double left = safeLower;
        double right = safeUpper;
        for (int iteration = 0; iteration < 100; ++iteration) {
            const double middle = 0.5*(left + right);
            if (derivative(middle) > 0.0) left = middle;
            else right = middle;
        }
        result.anisotropy = 0.5*(left + right);
        double information = 0.0;
        for (double legendre : secondLegendre) {
            const double scaled = legendre
                /(1.0 + result.anisotropy*legendre);
            information += scaled*scaled;
        }
        if (result.count >= 2 && information > 0.0) {
            result.standardError = 1.0/std::sqrt(information);
        }
    }
    return result;
}
#endif

#ifndef POSITRONIUM_VALIDATION_EXECUTABLE
#include "modules/crem_trajectory.hpp"

std::string formatTableValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string cutoffTimeLabel(double timeToCutoff) {
    if (std::isinf(timeToCutoff)) return "not reached";
    const double picoseconds = timeToCutoff * 1.0e12;
    std::ostringstream out;
    if (picoseconds < 1.0e-3) out << std::scientific << std::setprecision(2);
    else out << std::fixed << std::setprecision(3);
    out << picoseconds << " ps";
    return out.str();
}

std::string spinLabel(const Frame& frame) {
    const char* firstArrow = frame.firstDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    const char* secondArrow = frame.secondDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    return std::string(firstArrow) + " " + secondArrow;
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

std::uint64_t splitMix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// The default ceiling of 30 suits the trajectory experiments, whose sample
// sizes are limited by integration cost.  The annihilation generator is cheap
// enough to run millions of events, so its kinematics histograms pass a larger
// ceiling instead of throwing that resolution away.
int histogramBins(size_t count, int maximumBins = 30) {
    return std::clamp(static_cast<int>(std::lround(2.0 * std::sqrt(count))),
                      6, maximumBins);
}

// ---------------------------------------------------------------------------
// Plot colour convention.
//
// Every line, marker set and legend entry is coloured by WHERE ITS NUMBERS COME
// FROM, so one glance at any panel says what kind of statement it makes:
//
//   CREM simulation     blue         trajectories integrated by this program
//   sampled input       orange       quantities drawn from a random distribution
//   experimental        bluish green published measurements
//   theory / analytic   vermillion   closed-form reference curves
//
// The hues are the Okabe-Ito colour-vision-deficiency-safe palette rather than
// plain red/green/blue.  Plain red and green are precisely the pair that the
// roughly 8% of men with deuteranopia cannot separate, and here they would have
// carried the two most important comparisons.  Line and marker STYLE repeats
// the same information, so the panels survive greyscale printing too.
//
// Fill colour is reserved for classification (experiment 5) and uses pale tints
// that cannot be mistaken for the saturated provenance outlines.
namespace plot_style {

int crem()         { return TColor::GetColor("#0072B2"); }
int sampled()      { return TColor::GetColor("#E69F00"); }
int experimental() { return TColor::GetColor("#009E73"); }
int theory()       { return TColor::GetColor("#D55E00"); }

// Pale classification fills, in the order of InteractionOutcome.
int classificationFill(std::size_t slot) {
    static const std::array<const char*,6> tints{
        "#E8A0A0",  // Collision
        "#8FC3E8",  // Scattering
        "#8FD1B0",  // Para-Positronium
        "#C9A0DC",  // Ortho-Positronium
        "#BEBEBE",  // Unresolved
        "#EBC17A"}; // NumericalFailure
    return TColor::GetColor(tints[std::min(slot, tints.size() - 1)]);
}

// Legend text naming only the provenances a given panel actually contains, so
// the key stays short and never claims something the panel does not show.
std::string key(bool hasCrem, bool hasSampled, bool hasExperimental,
                bool hasTheory) {
    std::string result = "colours: ";
    bool first = true;
    const auto add = [&](bool present, const char* text) {
        if (!present) return;
        if (!first) result += ", ";
        result += text;
        first = false;
    };
    add(hasCrem, "blue CREM");
    add(hasSampled, "orange sampled");
    add(hasExperimental, "green measured");
    add(hasTheory, "vermillion theory");
    return result;
}

} // namespace plot_style

void styleHistogram(TH1D& histogram, int color) {
    histogram.SetDirectory(nullptr);
    histogram.SetLineColor(color);
    histogram.SetFillColorAlpha(color, 0.45);
    histogram.SetLineWidth(2);
}

struct GaussianFitSummary {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    std::size_t count = 0;
};

GaussianFitSummary gaussianMaximumLikelihood(const std::vector<double>& values) {
    GaussianFitSummary result;
    result.mean = 0.0;
    double secondMoment = 0.0;
    for (double value : values) {
        if (!std::isfinite(value)) continue;
        ++result.count;
        const double delta = value - result.mean;
        result.mean += delta/static_cast<double>(result.count);
        const double updatedDelta = value - result.mean;
        secondMoment += delta*updatedDelta;
    }
    if (result.count == 0) {
        result.mean = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    result.sigma = std::sqrt(std::max(0.0, secondMoment)
        / static_cast<double>(result.count));
    return result;
}

// One trajectory's contribution to a right-censored survival sample.
// observed=true means the collapse was actually seen at `time`; observed=false
// means the run stopped at `time` with the pair still bound, so all that is
// known is T > time.
struct SurvivalObservation { double time=0.0; bool observed=false; };

struct KaplanMeierPoint {
    double time=0.0;
    double survival=1.0;
    double standardError=0.0; // Greenwood
    int atRisk=0;
    int events=0;
};

struct KaplanMeierEstimate {
    std::vector<KaplanMeierPoint> curve;
    double medianSurvival=std::numeric_limits<double>::quiet_NaN();
    bool medianReached=false;
    // Restricted mean survival time: the area under the curve out to
    // `horizon`.  With censored observations beyond the last collapse the
    // unrestricted mean is not identifiable, so RMST is the honest summary.
    double restrictedMean=std::numeric_limits<double>::quiet_NaN();
    double restrictedMeanError=std::numeric_limits<double>::quiet_NaN();
    double horizon=std::numeric_limits<double>::quiet_NaN();
    double survivalAtHorizon=std::numeric_limits<double>::quiet_NaN();
    double largestCensoredTime=std::numeric_limits<double>::quiet_NaN();
    int eventCount=0;
    int censoredCount=0;
};

// Product-limit estimator with Greenwood standard errors.  This is what the
// bound-decay experiments need: a trajectory stopped by the wall-clock budget
// has NOT told us nothing, it has told us the collapse time exceeds the
// simulated time it reached, and averaging only the completed runs throws that
// away -- badly, because the budget preferentially stops the widest orbits,
// which are exactly the slowest to collapse.
KaplanMeierEstimate kaplanMeier(std::vector<SurvivalObservation> sample) {
    KaplanMeierEstimate result;
    sample.erase(std::remove_if(sample.begin(),sample.end(),
        [](const SurvivalObservation& o) {
            return !std::isfinite(o.time)||o.time<0.0;
        }),sample.end());
    if(sample.empty()) return result;
    // Ties: events are ordered before censorings at the same time, the
    // standard convention -- a run censored at exactly t was still at risk
    // when the collapse at t happened.
    std::sort(sample.begin(),sample.end(),
        [](const SurvivalObservation& a,const SurvivalObservation& b) {
            if(a.time!=b.time) return a.time<b.time;
            return a.observed&&!b.observed;
        });
    for(const SurvivalObservation& o : sample) {
        if(o.observed) ++result.eventCount;
        else {
            ++result.censoredCount;
            result.largestCensoredTime=std::isfinite(result.largestCensoredTime)
                ?std::max(result.largestCensoredTime,o.time):o.time;
        }
    }
    double survival=1.0;
    double greenwoodSum=0.0;
    std::size_t index=0;
    const int total=static_cast<int>(sample.size());
    result.curve.push_back({0.0,1.0,0.0,total,0});
    while(index<sample.size()) {
        const double time=sample[index].time;
        const int atRisk=total-static_cast<int>(index);
        int events=0;
        std::size_t next=index;
        while(next<sample.size()&&sample[next].time==time) {
            if(sample[next].observed) ++events;
            ++next;
        }
        if(events>0) {
            survival*=1.0-static_cast<double>(events)/static_cast<double>(atRisk);
            if(atRisk>events) {
                greenwoodSum+=static_cast<double>(events)
                    /(static_cast<double>(atRisk)
                      *static_cast<double>(atRisk-events));
            } else {
                greenwoodSum=std::numeric_limits<double>::infinity();
            }
            // If the last risk-set member fails, Greenwood's sum is infinite
            // while S(t) is exactly zero.  The literal product 0*sqrt(inf) is
            // NaN; its limiting pointwise standard error is zero.
            const double standardError=survival==0.0
                ?0.0:survival*std::sqrt(greenwoodSum);
            result.curve.push_back(
                {time,survival,standardError,atRisk,events});
        }
        index=next;
    }
    // An all-censored sample is still a valid (degenerate) KM estimate: no
    // observed event lowers S(t), so S(t)=1 through the largest censoring time.
    // Continue through the common integration below to report a finite horizon,
    // RMST=horizon and Greenwood error 0.  The median correctly remains
    // unreached; nothing here extrapolates a lifetime beyond the observed data.
    // Horizon: the last time the estimator actually observed something, so
    // RMST is not extrapolated past the data.
    result.horizon=std::max(result.curve.back().time,
        std::isfinite(result.largestCensoredTime)?result.largestCensoredTime
                                                 :result.curve.back().time);
    result.survivalAtHorizon=result.curve.back().survival;
    for(const KaplanMeierPoint& point : result.curve) {
        if(!result.medianReached&&point.survival<=0.5) {
            result.medianSurvival=point.time;
            result.medianReached=true;
        }
    }
    // RMST = integral of the step function, plus its standard error from the
    // usual sum of squared tail areas weighted by the Greenwood increments.
    double area=0.0;
    for(std::size_t i=0;i+1<result.curve.size();++i) {
        area+=result.curve[i].survival
             *(result.curve[i+1].time-result.curve[i].time);
    }
    area+=result.curve.back().survival
         *(result.horizon-result.curve.back().time);
    result.restrictedMean=area;
    double varianceSum=0.0;
    for(std::size_t i=1;i<result.curve.size();++i) {
        if(result.curve[i].events==0) continue;
        double tailArea=0.0;
        for(std::size_t j=i;j+1<result.curve.size();++j) {
            tailArea+=result.curve[j].survival
                     *(result.curve[j+1].time-result.curve[j].time);
        }
        tailArea+=result.curve.back().survival
                 *(result.horizon-result.curve.back().time);
        const int atRisk=result.curve[i].atRisk;
        const int events=result.curve[i].events;
        if(atRisk>events) {
            varianceSum+=tailArea*tailArea*static_cast<double>(events)
                /(static_cast<double>(atRisk)
                  *static_cast<double>(atRisk-events));
        }
    }
    result.restrictedMeanError=std::sqrt(std::max(0.0,varianceSum));
    return result;
}

std::string compactNumber(double value, int precision = 5) {
    std::ostringstream output;
    output << std::setprecision(precision) << value;
    return output.str();
}

TPaveText* drawAnalysisBox(std::vector<std::unique_ptr<TPaveText>>& storage,
                           double x1, double y1, double x2, double y2,
                           const std::vector<std::string>& lines,
                           double textSize = 0.027) {
    auto box = std::make_unique<TPaveText>(x1, y1, x2, y2, "NDC");
    box->SetFillColorAlpha(kWhite, 0.84);
    box->SetTextAlign(12);
    box->SetTextFont(42);
    box->SetTextSize(textSize);
    for (const std::string& line : lines) box->AddText(line.c_str());
    box->Draw();
    TPaveText* result = box.get();
    storage.push_back(std::move(box));
    return result;
}

std::unique_ptr<TF1> gaussianMleOverlay(const std::string& name,
                                        const GaussianFitSummary& fit,
                                        const TH1D& histogram, int color) {
    if (fit.count < 2 || !(fit.sigma > 0.0) || !std::isfinite(fit.sigma)) {
        return nullptr;
    }
    const double lower = histogram.GetXaxis()->GetXmin();
    const double upper = histogram.GetXaxis()->GetXmax();
    auto function = std::make_unique<TF1>(name.c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)", lower, upper);
    const double amplitude = static_cast<double>(fit.count)
        * histogram.GetBinWidth(1)/(std::sqrt(2.0*pi)*fit.sigma);
    function->SetParameters(amplitude, fit.mean, fit.sigma);
    function->SetLineColor(color);
    function->SetLineWidth(2);
    function->SetLineStyle(2);
    function->Draw("SAME");
    return function;
}

void reportExports(const std::vector<root_export::ExportResult>& results) {
    for (const root_export::ExportResult& result : results) {
        if (result) {
            std::cout << "Saved plot: " << result.path.string() << '\n';
        } else {
            std::cerr << "Warning: could not save " << result.path.string()
                      << ": " << result.error << '\n';
        }
    }
}

bool reportArchiveOperation(const statistics_archive::OperationResult& result,
                            const char* description) {
    if (result) {
        std::cout << "Saved " << description << ": "
                  << result.path.string() << '\n';
        return true;
    }
    std::cerr << "Warning: could not save " << description << ' '
              << result.path.string() << ": " << result.error << '\n';
    return false;
}

#include "modules/crem_collapse.hpp"

// runCount counts full CREM trajectories, which cost seconds each and set the
// collapse-time panel's statistics.  decayEventCount counts events of the
// independent ideal-vacuum annihilation generator, which costs microseconds and
// sets the photon-kinematics panels.  These used to be the same number, so the
// cheap photon histograms inherited the expensive trajectory budget and were
// capped at 100 entries for no physical reason.
// The photon panels are ANALYTIC reference curves, not Monte Carlo output.
// The annihilation generator is a quantum prescription that is deliberately
// independent of the classical model, and sampling it here only reproduced the
// very distribution it draws from: the 2-gamma line is a compile-time constant,
// the polar distribution is isotropic by construction, and the 3-gamma spectrum
// is the same Ore-Powell density that the reference curve already plots.  The
// sampler is exercised where a self-consistency check belongs, in
// positronium_validation.  What is drawn here is the exact reference.
int showBoundDecayStatistics(std::uint64_t seed, int selectedPhenomenon,
                             int runCount, double wallClockBudgetSeconds) {
    // Experiments 1 and 2 are positronium experiments, not pair-general ones,
    // and --pair must not be able to pretend otherwise.
    //
    // The CREM inspiral they measure IS pair-general and correct for any pair:
    // run for p+e- it returns a mean collapse time of 17.8 ps against the
    // textbook classical hydrogen value of about 16 ps, with an initial
    // orbital period of 0.164 fs against the Bohr orbit's 0.152 fs.  What is
    // not pair-general is everything it is reported ALONGSIDE.  The
    // annihilation-time spectrum, the tau_exp scale reference and the
    // para/ortho labelling all come from the MEASURED positronium lifetimes
    // (Al-Ramadhan & Gidley, PRL 72 (1994); Vallery et al., PRL 90 (2003)).
    //
    // A proton and an electron do not annihilate at all -- they form hydrogen,
    // which is stable -- so for that pair the whole overlay is positronium's
    // data pasted onto someone else's dynamics.  Nor do the numbers transfer
    // to the pairs that DO annihilate: true muonium and protonium annihilate
    // through channels this model does not carry, so the guard is positronium
    // specifically and not "a particle with its own antiparticle".
    if(!isPositronium(activePair)) {
        std::cerr
            << "Experiments 1 and 2 measure a positronium CREM collapse "
               "against the measured positronium annihilation lifetimes, so "
               "they only mean anything for electron+positron.\n"
            << "The selected pair is " << firstSpecies.name << " + "
            << secondSpecies.name << ", which does not annihilate into those "
               "channels; the classical inspiral would be computed correctly "
               "but reported beside data that does not describe it.\n"
            << "Experiments 3, 4 and 5 carry no annihilation reference and "
               "run for any attracting pair.\n";
        return 2;
    }
    const bool isPara = selectedPhenomenon == 1;
    const double timeScale = isPara ? 1.0e12 : 1.0e9;
    const char* timeUnit = isPara ? "ps" : "ns";

    const std::vector<CremCollapseEstimate> collapseEstimates=
        runCremCollapseExperiment(seed,selectedPhenomenon,runCount,
                                  wallClockBudgetSeconds);
    std::vector<double> decayTimes;
    std::vector<double> calibrationPowers;
    // Right-censored sample for the product-limit estimator: every trajectory
    // contributes, a completed one as an observed collapse and a stopped one
    // as a lower bound at the simulated time it reached.
    std::vector<SurvivalObservation> survivalSample;
    survivalSample.reserve(static_cast<size_t>(runCount));
    // External-reference comparisons, collected only from completed runs so
    // that measurement and reference span the same stretch of orbit.
    std::vector<double> measuredCollapse, analyticCollapse; // plot units
    std::vector<double> larmorRatios;
    std::vector<double> dipoleCouplingsGHz;
    // Period at both ends of the inspiral and the revolutions between them,
    // collected only from trajectories that ran to the boundary so the three
    // numbers describe the same complete collapses as the lifetime above.
    std::vector<double> initialPeriods, finalPeriods, revolutionCounts;
    int reachedCutoffCount = 0;
    int observationLimitCount = 0;
    int noSecularLossCount = 0;
    int calibrationFailureCount = 0;
    double censoredSimulatedTimeMax = 0.0;
    double censoredSimulatedTimeSum = 0.0;
    decayTimes.reserve(static_cast<size_t>(runCount));
    calibrationPowers.reserve(static_cast<size_t>(runCount));

    for (int index = 0; index < runCount; ++index) {
        const CremCollapseEstimate& estimate =
            collapseEstimates[static_cast<size_t>(index)];
        if(std::isfinite(estimate.lifetimeSeconds))
            decayTimes.push_back(estimate.lifetimeSeconds*timeScale);
        if(estimate.calibrationOutcome==SimulationOutcome::ReachedCutoff) {
            if(std::isfinite(estimate.initialPeriodSeconds))
                initialPeriods.push_back(estimate.initialPeriodSeconds);
            if(std::isfinite(estimate.finalPeriodSeconds))
                finalPeriods.push_back(estimate.finalPeriodSeconds);
            if(std::isfinite(estimate.revolutions))
                revolutionCounts.push_back(estimate.revolutions);
            if(std::isfinite(estimate.lifetimeSeconds)
               &&std::isfinite(estimate.analyticCollapseSeconds)
               &&estimate.analyticCollapseSeconds>0.0) {
                measuredCollapse.push_back(estimate.lifetimeSeconds*timeScale);
                analyticCollapse.push_back(
                    estimate.analyticCollapseSeconds*timeScale);
            }
        }
        // Independent of completion: both are properties of the prepared
        // orbit and of orbits actually resolved, so a censored run still
        // contributes a valid measurement of them.
        if(std::isfinite(estimate.larmorPowerRatio)
           && estimate.larmorPowerRatio > 0.0) {
            larmorRatios.push_back(estimate.larmorPowerRatio);
        }
        if(std::isfinite(estimate.dipoleCouplingHz)) {
            dipoleCouplingsGHz.push_back(estimate.dipoleCouplingHz*1.0e-9);
        }
        if(std::isfinite(estimate.meanRadiatedPowerWatts)
           && estimate.meanRadiatedPowerWatts > 0.0) {
            calibrationPowers.push_back(estimate.meanRadiatedPowerWatts);
        }
        switch(estimate.calibrationOutcome) {
            case SimulationOutcome::ReachedCutoff:
                ++reachedCutoffCount;
                if(std::isfinite(estimate.lifetimeSeconds)) {
                    survivalSample.push_back(
                        {estimate.lifetimeSeconds*timeScale,true});
                }
                break;
            case SimulationOutcome::ObservationLimit:
                // Both stopped states are right-censored: the pair was still
                // bound when observation ended, so the collapse time is known
                // only to exceed the simulated time reached.  A non-decaying
                // trajectory is the limiting case of that.
                survivalSample.push_back(
                    {estimate.calibrationSeconds*timeScale,false});
                if(estimate.secularLossAbsent) { ++noSecularLossCount; break; }
                ++observationLimitCount;
                censoredSimulatedTimeMax=std::max(
                    censoredSimulatedTimeMax,estimate.calibrationSeconds);
                censoredSimulatedTimeSum+=estimate.calibrationSeconds;
                break;
            case SimulationOutcome::NumericalFailure:
                // Excluded outright: a non-finite state carries no bound on
                // the collapse time in either direction.
                ++calibrationFailureCount; break;
        }
    }
    const KaplanMeierEstimate survival=kaplanMeier(survivalSample);
    const int usableTrajectories=survival.eventCount+survival.censoredCount;
    const double completionPercent=usableTrajectories>0
        ?100.0*static_cast<double>(survival.eventCount)
             /static_cast<double>(usableTrajectories)
        :0.0;


    const GaussianFitSummary initialPeriodMoments=
        gaussianMaximumLikelihood(initialPeriods);
    const GaussianFitSummary finalPeriodMoments=
        gaussianMaximumLikelihood(finalPeriods);
    const GaussianFitSummary revolutionMoments=
        gaussianMaximumLikelihood(revolutionCounts);
    const GaussianFitSummary collapseMoments=gaussianMaximumLikelihood(decayTimes);
    const double estimatedLifetime=collapseMoments.mean/timeScale;
    // Standard error of the sample mean.  The previous mean/sqrt(N) is the MLE
    // error of an *exponential* sample, but these collapse times are not
    // exponentially distributed: every trajectory starts at exactly r=a0 and
    // only the initial velocity is randomized, so the sample is narrow and
    // unimodal.  sigma/sqrt(N) is the estimator that matches the actual data.
    const double estimatedError=collapseMoments.count>1
        ?collapseMoments.sigma/timeScale
            /std::sqrt(static_cast<double>(collapseMoments.count))
        :std::numeric_limits<double>::quiet_NaN();
    const double relativeSpread=collapseMoments.mean>0.0
        ?collapseMoments.sigma/collapseMoments.mean
        :std::numeric_limits<double>::quiet_NaN();
    const statistics_archive::ScientificValue& lifetimeReference =
        statistics_archive::scientificValue(isPara
            ? "para_lifetime_from_rate" : "ortho_lifetime_from_rate");
    // The primary measured quantity is the decay rate; the lifetime above is
    // its derived reciprocal.  Quote the rate too, since that is what was
    // actually measured.
    const statistics_archive::ScientificValue& rateReference =
        statistics_archive::scientificValue(isPara
            ? "para_decay_rate_measurement" : "ortho_decay_rate_measurement");
    const double experimentalLifetime = lifetimeReference.value;
    const double experimentalLifetimeError = lifetimeReference.totalUncertainty;
    std::cout << (isPara ? "Para-positronium" : "Ortho-positronium")
              << " study: " << runCount << " CREM trajectories ("
              << reachedCutoffCount << " reached the collision boundary, "
              << observationLimitCount << " censored by the per-event "
                 "wall-clock budget, " << noSecularLossCount
              << " not decaying, " << calibrationFailureCount
              << " numerical failures).\n"
              << "Photon panels are exact reference curves, not samples: the "
                 "annihilation generator is a quantum prescription independent\n"
                 "of the classical model, and its self-consistency is checked "
                 "in positronium_validation.\n"
              << "Collapse time over all "
              << (survival.eventCount+survival.censoredCount)
              << " usable trajectories (" << survival.eventCount
              << " observed collapses, " << survival.censoredCount
              << " right-censored\nat the simulated time they reached).  "
                 "Completion fraction: " << completionPercent << "%.\n";
    if(survival.medianReached) {
        std::cout << "  Kaplan-Meier median    " << survival.medianSurvival
                  << ' ' << timeUnit << '\n';
    } else {
        std::cout << "  Kaplan-Meier median    not reached: S(t) is still "
                  << survival.survivalAtHorizon << " at "
                  << survival.horizon << ' ' << timeUnit << '\n';
    }
    if(std::isfinite(survival.restrictedMean)) {
        std::cout << "  Kaplan-Meier RMST      " << survival.restrictedMean;
        if(survival.eventCount>0) {
            std::cout << " +/- " << survival.restrictedMeanError;
        }
        std::cout << ' ' << timeUnit;
        if(survival.eventCount==0)
            std::cout << "  (all observations censored; uncertainty not estimable)";
        std::cout << "  (area under S(t) out to " << survival.horizon
                  << ' ' << timeUnit << ")\n";
    }
    if(collapseMoments.count>0) {
        std::cout << "  mean of completed runs " << estimatedLifetime * timeScale
                  << " +/- " << estimatedError * timeScale << ' ' << timeUnit
                  << " (sigma/mean = " << relativeSpread << ")\n";
    } else {
        std::cout << "  mean of completed runs unavailable: no collapse was "
                     "observed\n";
    }
    if(completionPercent < 90.0 && std::isfinite(survival.restrictedMean)) {
        if(survival.eventCount==0) {
            std::cout << "  WARNING: every usable trajectory is censored.  "
                         "Kaplan-Meier therefore has S(t)=1 through "
                      << survival.horizon << ' ' << timeUnit
                      << " and RMST=" << survival.restrictedMean << ' '
                      << timeUnit << ".\n  This finite restricted estimate does "
                         "not identify the median or the unrestricted mean; "
                         "raise --crem-wallclock-budget-s until collapses are "
                         "observed.\n";
        } else {
            std::cout << "  WARNING: with " << completionPercent
                      << "% completion neither figure is an unbiased estimate.  "
                         "Their biases have OPPOSITE sign, so the\n  true value is "
                         "expected between them -- a plausibility range, not a "
                         "proven bound:\n    lower end "
                      << estimatedLifetime * timeScale << ' ' << timeUnit
                      << " -- the completed-run mean.  Collapse time scales as a^3 "
                         "while the orbits that must be\n      integrated scale as "
                         "a^(3/2), so the budget preferentially stops the widest "
                         "orbits, which are the\n      slowest to collapse; "
                         "averaging the survivors is biased low.\n    upper end "
                      << survival.restrictedMean << ' ' << timeUnit
                      << " -- the Kaplan-Meier RMST.  Once censoring begins the "
                         "estimator observes almost no further\n      collapses, "
                         "so S(t) is held up and the area under it is biased high.\n"
                         "  Censoring is INFORMATIVE here: the censoring time and "
                         "the collapse time are driven by the same\n  semi-major "
                         "axis, which breaks the independence Kaplan-Meier assumes."
                         "  Raise --crem-wallclock-budget-s\n  until the completion "
                         "fraction approaches 100%; the two ends then converge onto "
                         "the true value.\n";
        }
    }
    std::cout << "Model: full CREM mechanical integration under the active "
                 "--radiation-reaction model, run until the pair's periapsis\n"
                 "reaches 0.1*a0 or the per-event wall-clock budget is spent "
                 "(then censored, not extrapolated).  The last stretch down to "
                 "the\n0.01*a0 boundary is truncated deliberately: with t ~ a^3 "
                 "it is 0.1% of the collapse time, below the 3% per-jump\n"
                 "tolerance, and the one-period measurement window stops being "
                 "valid there.  External lifetime is comparison only.\n"
              << "Caution: the CREM collapse time is a classical inspiral time."
                 " Para and ortho differ here only through the initial dipole\n"
                 "alignment, whose coupling is ~1e-5 of the Coulomb potential,"
                 " so both channels yield the same collapse distribution while\n"
                 "their measured annihilation lifetimes differ by ~1000x.  The"
                 " comparison is a scale reference, not a prediction.\n";
    if(revolutionMoments.count>0) {
        std::cout<<"\nQuasi-closed orbit (the "<<revolutionMoments.count
                 <<" trajectories that reached the boundary):\n"
                 <<"  T of the sampled initial orbit  "
                 <<initialPeriodMoments.mean*1.0e15<<" +/- "
                 <<initialPeriodMoments.sigma*1.0e15<<" fs\n"
                 <<"  T of the last resolved orbit    "
                 <<finalPeriodMoments.mean*1.0e15<<" +/- "
                 <<finalPeriodMoments.sigma*1.0e15<<" fs\n"
                 <<"  revolutions between them        "
                 <<revolutionMoments.mean
                 <<" +/- "<<revolutionMoments.sigma<<'\n'
                 <<"The orbit is not closed: it is an inspiral, so T shrinks as"
                   " a^(3/2) while the pair sinks and no single period\n"
                   "describes the whole run.  T is the Kepler period "
                   "2*pi*sqrt(mu*a^3/(k e^2)) of the osculating orbit.  Every "
                   "trajectory starts at\nRADIUS a_0 but with a sub-circular "
                   "tangential speed, so its semi-major axis is already below "
                   "a_0 and T is below the\n107.5 as of a circular orbit at "
                   "a_0.  The run ends when the PERIAPSIS reaches 0.1*a_0.  "
                   "The revolution count is\naccumulated by the orbit-averaged "
                   "integrator, resolved and skipped orbits alike.\n";
    }
    if(observationLimitCount>0) {
        std::cout<<"Note: "<<observationLimitCount<<" of "<<runCount
                 <<" trajectories did not reach the boundary within the "
                    "per-event wall-clock budget and were excluded from the "
                    "mean above; raise --crem-wallclock-budget-s to reduce "
                    "censoring.\n"
                 <<"  Simulated time covered by censored trajectories: mean "
                 <<(censoredSimulatedTimeSum/observationLimitCount)*1.0e12
                 <<" ps, max "<<censoredSimulatedTimeMax*1.0e12<<" ps"
                    " (out of the classical estimate of a few ps to collapse).\n";
    }
    if(noSecularLossCount>0) {
        std::cout<<"Note: "<<noSecularLossCount<<" of "<<runCount
                 <<" trajectories showed no secular energy loss over the first "
                    "measured orbit, so they are not decaying at all and no\n"
                    "  collapse time exists for them.  This is the expected "
                    "result with --radiation-reaction disabled: the charge "
                    "self-force is the only\n  channel that removes orbital "
                    "energy.  Raising --crem-wallclock-budget-s cannot change "
                    "it.\n";
    }
    if(collapseMoments.count==0) {
        if(noSecularLossCount>=observationLimitCount) {
            std::cerr<<"No CREM trajectory decayed: the active "
                        "--radiation-reaction model removes no orbital energy, "
                        "so there is no\ninspiral to measure and no plots were "
                        "produced. Select individual, coherent or automatic.\n";
            return 2;
        }
        std::cerr<<"No CREM trajectory reached the collision boundary within "
                    "the wall-clock budget; no plots were produced. Raise "
                    "--crem-wallclock-budget-s and try again.\n";
        return 2;
    }
    // Statistical rendering is a pure batch job.  Select the virtual ROOT
    // backend before constructing any canvas; creating TApplication first can
    // permanently bind ROOT 6.40 to an X11 painter in this process.
    gROOT->SetBatch(kTRUE);
    root_export::preparePdfExporter();
    gStyle->SetOptStat(1110);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    const std::string stateName = isPara ? "Para-positronium" : "Ortho-positronium";
    const std::string canvasTitle = stateName + " decay observables (ideal vacuum)";
    TCanvas canvas("decay_statistics", canvasTitle.c_str(), 1280, 900);
    canvas.SetFillColor(kWhite);
    TPad distributionsPage("decay_distributions_page", "Decay distributions",
                           0.0, 0.0, 1.0, 1.0);
    TPad diagnosticsPage("decay_diagnostics_page", "Numerical closure diagnostics",
                         0.0, 0.0, 1.0, 1.0);
    for (TPad* page : {&distributionsPage, &diagnosticsPage}) {
        page->SetFillColor(kWhite);
        page->SetFillStyle(1001);
    }
    canvas.cd();
    distributionsPage.Draw();
    diagnosticsPage.Draw();
    // Four distribution pads, down from six: the two exact photon-kinematics
    // pads carried no CREM output at all and the sixth was a text card.
    distributionsPage.cd();
    distributionsPage.Divide(2, 2, 0.006, 0.006);
    diagnosticsPage.cd();
    diagnosticsPage.Divide(1, 2, 0.006, 0.006);
    std::vector<std::unique_ptr<TPaveText>> analysisBoxes;
    std::vector<std::unique_ptr<TF1>> analysisFunctions;

    std::vector<double> sortedCollapseTimes=decayTimes;
    std::sort(sortedCollapseTimes.begin(),sortedCollapseTimes.end());
    // Kaplan-Meier step curve, drawn as an explicit staircase so the flat
    // stretches between collapses are visible rather than interpolated away.
    // The old curve plotted (N-i)/N over the completed runs only, which drew
    // a line reaching zero survival even when 90% of the sample never
    // collapsed within the budget.
    std::vector<double> survivalStepTime, survivalStepValue;
    std::vector<double> survivalBandTime, survivalBandValue, survivalBandError;
    for(const KaplanMeierPoint& point : survival.curve) {
        if(!survivalStepTime.empty()) {
            survivalStepTime.push_back(point.time);
            survivalStepValue.push_back(survivalStepValue.back());
        }
        survivalStepTime.push_back(point.time);
        survivalStepValue.push_back(point.survival);
        if(point.events>0&&std::isfinite(point.standardError)) {
            survivalBandTime.push_back(point.time);
            survivalBandValue.push_back(point.survival);
            survivalBandError.push_back(point.standardError);
        }
    }
    // Extend the last step out to the horizon so the curve does not stop in
    // mid-air at the final observed collapse.
    if(!survivalStepTime.empty()&&std::isfinite(survival.horizon)
       &&survival.horizon>survivalStepTime.back()) {
        survivalStepTime.push_back(survival.horizon);
        survivalStepValue.push_back(survivalStepValue.back());
    }
    // Censoring marks: where trajectories left the risk set without collapsing.
    std::vector<double> censorTime, censorValue;
    {
        std::vector<double> censoredTimes;
        for(const SurvivalObservation& observation : survivalSample)
            if(!observation.observed) censoredTimes.push_back(observation.time);
        std::sort(censoredTimes.begin(),censoredTimes.end());
        for(double time : censoredTimes) {
            double atTime=1.0;
            for(const KaplanMeierPoint& point : survival.curve)
                if(point.time<=time) atTime=point.survival;
            censorTime.push_back(time);
            censorValue.push_back(atTime);
        }
    }
    // Scale the axis to the simulated data, not to the experimental lifetime.
    // For o-Ps the two differ by ~8 orders of magnitude, and forcing tau_exp
    // into the range collapsed the whole CREM sample onto a single pixel at
    // the left edge.  The experimental curve is still drawn: staying pinned at
    // survival=1 across this window is exactly the honest visual statement
    // that the classical inspiral is far faster than the measured decay.
    // The axis must cover the censored observations too, otherwise the very
    // trajectories the estimator now accounts for fall off the right edge.
    const double survivalHorizon=std::isfinite(survival.horizon)
        ?std::max(survival.horizon,sortedCollapseTimes.back())
        :sortedCollapseTimes.back();
    const double lifetimeLower=std::max(
        0.5*sortedCollapseTimes.front(),
        1.0e-6*std::max(survivalHorizon,1.0e-300));
    const double lifetimeUpper=std::max(1.15*survivalHorizon,
        4.0*lifetimeLower);
    // Both are already in plot units (ps for p-Ps, ns for o-Ps).
    const double experimentalRatio=collapseMoments.mean>0.0
        ?experimentalLifetime/collapseMoments.mean
        :std::numeric_limits<double>::quiet_NaN();
    TGraph lifetimeSurvival(static_cast<int>(survivalStepTime.size()),
        survivalStepTime.data(),survivalStepValue.data());
    std::ostringstream lifetimeTitle;
    lifetimeTitle<<"CREM collapse survival (Kaplan-Meier, censoring-aware);t ["
                 <<timeUnit<<"];Survival fraction";
    lifetimeSurvival.SetTitle(lifetimeTitle.str().c_str());
    lifetimeSurvival.SetLineColor(plot_style::crem());
    lifetimeSurvival.SetLineWidth(2);
    lifetimeSurvival.SetMarkerColor(plot_style::crem());
    lifetimeSurvival.SetMarkerStyle(1);
    // Greenwood pointwise standard errors at the observed collapses.
    TGraphErrors survivalBand(static_cast<int>(survivalBandTime.size()),
        survivalBandTime.data(),survivalBandValue.data(),
        nullptr,survivalBandError.data());
    survivalBand.SetLineColor(plot_style::crem());
    survivalBand.SetMarkerColor(plot_style::crem());
    survivalBand.SetMarkerStyle(20);
    survivalBand.SetMarkerSize(0.6);
    // Vertical ticks where a trajectory was censored: the information the old
    // completed-runs-only curve discarded.
    TGraph censorMarks(static_cast<int>(censorTime.size()),
        censorTime.data(),censorValue.data());
    censorMarks.SetMarkerColor(plot_style::crem());
    censorMarks.SetMarkerStyle(2); // upright cross, the survival-analysis tick
    censorMarks.SetMarkerSize(0.9);
    TF1 experimentalCurve("experimental_lifetime_distribution", "exp(-x/[0])",
                          lifetimeLower, lifetimeUpper);
    experimentalCurve.SetParameter(0,experimentalLifetime);
    experimentalCurve.SetLineColor(plot_style::experimental());
    experimentalCurve.SetLineWidth(3);
    experimentalCurve.SetLineStyle(3);
    // The descriptive exp() through the completed-run mean was removed: the
    // sample is not exponential (the panel said so itself), and the mean it
    // was drawn through is the biased one.  Drawing it next to a
    // censoring-aware estimate invited exactly the misreading this panel now
    // warns about.

    distributionsPage.cd(1);
    gPad->SetGrid();
    // A log axis spanning less than a decade draws no labelled major tick at
    // all -- the panel came out with a bare "t [ps]" and no numbers on it.
    // Use the log axis only when the sample really is decades wide.
    const bool survivalUsesLogAxis =
        lifetimeLower > 0.0 && lifetimeUpper/lifetimeLower > 20.0;
    const double survivalAxisLower =
        survivalUsesLogAxis ? lifetimeLower : 0.0;
    if(survivalUsesLogAxis) gPad->SetLogx();
    lifetimeSurvival.Draw("AL");
    lifetimeSurvival.GetXaxis()->SetLimits(survivalAxisLower,lifetimeUpper);
    lifetimeSurvival.SetMinimum(0.0);
    lifetimeSurvival.SetMaximum(1.08);
    survivalBand.Draw("P");
    if(censorMarks.GetN()>0) censorMarks.Draw("P");
    experimentalCurve.Draw("SAME");
    // Which corner is free depends on the shape the curve actually takes, so
    // pick it from the data instead of hard-coding one.  With high completion
    // the staircase descends from the top-left and the upper right is empty;
    // under heavy censoring it stays pinned near 1 across the whole width and
    // only the lower left is empty.  A fixed corner covered the curve in one
    // case or the other.
    const bool boxUpperRight = completionPercent >= 50.0;
    const double boxX1 = boxUpperRight ? 0.53 : 0.13;
    const double boxY1 = boxUpperRight ? 0.48 : 0.13;
    const double boxX2 = boxUpperRight ? 0.97 : 0.63;
    const double boxY2 = boxUpperRight ? 0.90 : 0.56;
    drawAnalysisBox(analysisBoxes, boxX1, boxY1, boxX2, boxY2, {
        "N = " + std::to_string(runCount) + ":  "
            + std::to_string(survival.eventCount) + " collapsed,  "
            + std::to_string(survival.censoredCount) + " censored  ("
            + compactNumber(completionPercent, 3) + "% complete)",
        "Kaplan-Meier; ticks = censored, bars = Greenwood",
        survival.medianReached
            ? "KM median = " + compactNumber(survival.medianSurvival)
                + " " + timeUnit
            : "KM median not reached; S = "
                + compactNumber(survival.survivalAtHorizon, 3) + " at "
                + compactNumber(survivalHorizon) + " " + timeUnit,
        "KM RMST = " + compactNumber(survival.restrictedMean)
            + " #pm " + compactNumber(survival.restrictedMeanError)
            + " " + timeUnit,
        "completed-run mean = " + compactNumber(collapseMoments.mean)
            + " " + timeUnit,
        completionPercent < 90.0
            ? "Informative censoring: the two biases have opposite"
            : "High completion: the two figures agree.",
        completionPercent < 90.0
            ? "sign. Raise --crem-wallclock-budget-s to close the gap."
            : "",
        plot_style::key(true, false, true, false),
        "#tau_{exp} = " + compactNumber(experimentalLifetime) + " " + timeUnit
            + " (scale reference only, not a prediction)",
        isPara ? "Al-Ramadhan & Gidley, PRL 72 (1994)"
               : "Vallery et al., PRL 90 (2003)"
    }, 0.019);

    // Annihilation-time spectrum built from the MEASURED decay rate, drawn
    // analytically.  The published rate is real data; turning it into a
    // spectrum needs no Monte Carlo, because N(t) = exp(-t/tau)/tau is exactly
    // what an exponential decay law with that rate implies.  Sampling it would
    // only have added noise to a curve we already know in closed form.
    const double annihilationUpper = 6.0*experimentalLifetime;
    TF1 annihilationSpectrum("annihilation_time_spectrum",
                             "exp(-x/[0])/[0]", 0.0, annihilationUpper);
    annihilationSpectrum.SetParameter(0, experimentalLifetime);
    annihilationSpectrum.SetParName(0, "tau");
    annihilationSpectrum.SetLineColor(plot_style::experimental());
    annihilationSpectrum.SetLineWidth(3);
    annihilationSpectrum.SetNpx(600);
    std::ostringstream annihilationTitle;
    annihilationTitle << "Annihilation-time spectrum from the measured rate;t ["
                      << timeUnit << "];(1/N) dN/dt [" << timeUnit << "^{-1}]";
    annihilationSpectrum.SetTitle(annihilationTitle.str().c_str());
    distributionsPage.cd(2);
    gPad->SetGrid();
    annihilationSpectrum.Draw("L");
    // The classical collapse time is orders of magnitude shorter, so the marker
    // sits hard against the left edge.  That is the honest picture and the box
    // gives the ratio in numbers.
    std::unique_ptr<TLine> collapseMarker;
    if (collapseMoments.mean > 0.0
        && collapseMoments.mean < annihilationUpper) {
        collapseMarker = std::make_unique<TLine>(
            collapseMoments.mean, 0.0,
            collapseMoments.mean, 1.0/experimentalLifetime);
        collapseMarker->SetLineColor(plot_style::crem());
        collapseMarker->SetLineWidth(2);
        collapseMarker->SetLineStyle(2);
        collapseMarker->Draw();
    }
    drawAnalysisBox(analysisBoxes, 0.40, 0.50, 0.95, 0.91, {
        plot_style::key(true, false, true, false),
        "Measured data, drawn analytically - no Monte Carlo.",
        "Decay law N(t) #propto exp(-t/#tau) with the published",
        "rate; the spectrum follows in closed form.",
        "#lambda = " + compactNumber(rateReference.value, 6) + " #pm "
            + compactNumber(rateReference.totalUncertainty, 2) + " "
            + rateReference.unit,
        "#tau_{exp} = " + compactNumber(experimentalLifetime) + " #pm "
            + compactNumber(experimentalLifetimeError) + " " + timeUnit,
        isPara ? "Al-Ramadhan & Gidley, PRL 72 (1994)"
               : "Vallery et al., PRL 90 (2003)",
        "blue dashed: CREM collapse #LTt#GT = "
            + compactNumber(collapseMoments.mean, 4) + " " + timeUnit,
        "#tau_{exp}/#LTt#GT #approx " + compactNumber(experimentalRatio, 3)
            + ": the classical inspiral is far faster.",
        "The two are different processes and are not",
        "expected to agree; this panel sets the scale."
    }, 0.0185);

    // --- Pad 3: measured collapse time against closed-form electrodynamics ---
    // The engine is compared with the textbook classical inspiral of a
    // radiating electric dipole, evaluated between the SAME two osculating
    // orbits each trajectory actually covered.  Nothing from CREM enters the
    // reference, so a departure from the diagonal is a statement about the
    // engine, not about positronium.
    distributionsPage.cd(3);
    gPad->SetGrid();
    double collapseAxisUpper = 1.0;
    for (std::size_t i = 0; i < measuredCollapse.size(); ++i) {
        collapseAxisUpper = std::max({collapseAxisUpper,
            measuredCollapse[i], analyticCollapse[i]});
    }
    collapseAxisUpper *= 1.15;
    TGraph collapseVersusTheory(static_cast<int>(measuredCollapse.size()),
        analyticCollapse.data(), measuredCollapse.data());
    std::ostringstream collapseCompareTitle;
    collapseCompareTitle << "CREM collapse time vs closed-form classical "
                            "inspiral;t_{classical} [" << timeUnit
                         << "];t_{CREM} [" << timeUnit << ']';
    collapseVersusTheory.SetTitle(collapseCompareTitle.str().c_str());
    collapseVersusTheory.SetMarkerColor(plot_style::crem());
    collapseVersusTheory.SetMarkerStyle(20);
    collapseVersusTheory.SetMarkerSize(0.7);
    std::unique_ptr<TF1> collapseDiagonal;
    std::vector<double> collapseRatios;
    for (std::size_t i = 0; i < measuredCollapse.size(); ++i) {
        if (analyticCollapse[i] > 0.0)
            collapseRatios.push_back(measuredCollapse[i]/analyticCollapse[i]);
    }
    const GaussianFitSummary collapseRatioMoments =
        gaussianMaximumLikelihood(collapseRatios);
    if (collapseVersusTheory.GetN() > 0) {
        collapseVersusTheory.Draw("AP");
        collapseVersusTheory.GetXaxis()->SetLimits(0.0, collapseAxisUpper);
        collapseVersusTheory.SetMinimum(0.0);
        collapseVersusTheory.SetMaximum(collapseAxisUpper);
        collapseDiagonal = std::make_unique<TF1>("collapse_diagonal", "x",
            0.0, collapseAxisUpper);
        collapseDiagonal->SetLineColor(plot_style::theory());
        collapseDiagonal->SetLineWidth(3);
        collapseDiagonal->SetLineStyle(2);
        collapseDiagonal->Draw("SAME");
    }
    // Lower right: the engine radiates less than the reference, so every point
    // sits ABOVE the diagonal and the region below it is empty.
    drawAnalysisBox(analysisBoxes, 0.44, 0.13, 0.96, 0.46, {
        plot_style::key(true, false, false, true),
        "Reference: da/dt = -C/a^{2}, C = 8ke^{4}/(6#pi#varepsilon_{0}c^{3}m^{2}),",
        "orbit-averaged with the DIPOLE factor (1+e^{2}/2)/(1-e^{2})^{5/2}",
        "evaluated at each trajectory's own a and e. Zero free parameters.",
        "completed trajectories: N = "
            + std::to_string(measuredCollapse.size()),
        "#LTt_{CREM}/t_{classical}#GT = "
            + compactNumber(collapseRatioMoments.mean, 4)
            + " #pm " + compactNumber(collapseRatioMoments.sigma, 3),
        "dashed: exact agreement. Eccentricity is held fixed in",
        "the reference; radiation actually circularizes the orbit."
    }, 0.019);

    // --- Pad 4: radiation sector against the Larmor rate ---
    // Per checkpoint the engine measures the orbital energy lost over one
    // resolved orbit; the Larmor power of the coherent electric dipole for
    // that same osculating orbit is the reference.  A ratio of 1 means the
    // engine reproduces coherent dipole radiation.  A ratio near 1/2 would
    // mean the two charges are radiating incoherently instead: their dipole
    // contributions add as 4 coherently against 2 incoherently.
    distributionsPage.cd(4);
    gPad->SetGrid();
    double larmorLower = std::numeric_limits<double>::infinity();
    double larmorUpper = 0.0;
    for (double value : larmorRatios) {
        larmorLower = std::min(larmorLower, value);
        larmorUpper = std::max(larmorUpper, value);
    }
    if (!(larmorLower < larmorUpper)) { larmorLower = 0.0; larmorUpper = 1.5; }
    // Always keep both reference marks on the axis.  Left to the data alone
    // the range collapsed onto the model's own 0.3% spread and the coherent
    // mark at 1 fell off the plot, hiding the very comparison the panel is for.
    larmorLower = std::min(larmorLower, 0.45);
    larmorUpper = std::max(larmorUpper, 1.05);
    const double larmorPadding = 0.08*(larmorUpper - larmorLower) + 1.0e-6;
    TH1D larmorRatioHistogram("crem_larmor_ratio",
        "Radiated power: CREM measured vs Larmor;"
        "P_{CREM} / P_{Larmor};Trajectories",
        histogramBins(larmorRatios.size()),
        larmorLower - larmorPadding, larmorUpper + larmorPadding);
    styleHistogram(larmorRatioHistogram, plot_style::crem());
    larmorRatioHistogram.SetStats(false);
    for (double value : larmorRatios) larmorRatioHistogram.Fill(value);
    larmorRatioHistogram.Draw("HIST");
    const GaussianFitSummary larmorMoments =
        gaussianMaximumLikelihood(larmorRatios);
    std::unique_ptr<TLine> coherentMarker, incoherentMarker;
    if (!larmorRatios.empty()) {
        const double markerTop = 1.05*larmorRatioHistogram.GetMaximum();
        if (larmorLower - larmorPadding <= 1.0
            && 1.0 <= larmorUpper + larmorPadding) {
            coherentMarker = std::make_unique<TLine>(1.0, 0.0, 1.0, markerTop);
            coherentMarker->SetLineColor(plot_style::theory());
            coherentMarker->SetLineWidth(3);
            coherentMarker->Draw();
        }
        if (larmorLower - larmorPadding <= 0.5
            && 0.5 <= larmorUpper + larmorPadding) {
            incoherentMarker = std::make_unique<TLine>(0.5, 0.0, 0.5, markerTop);
            incoherentMarker->SetLineColor(plot_style::theory());
            incoherentMarker->SetLineWidth(2);
            incoherentMarker->SetLineStyle(3);
            incoherentMarker->Draw();
        }
    }
    // Upper right: with the axis now forced to span both reference marks the
    // model's narrow peak sits on the left, leaving this corner clear.
    drawAnalysisBox(analysisBoxes, 0.48, 0.56, 0.96, 0.91, {
        plot_style::key(true, false, false, true),
        "Reference: P = |d''|^{2}/(6#pi#varepsilon_{0}c^{3}), d = e r,",
        "orbit-averaged at each checkpoint's own a and e.",
        "checkpoint samples: N = " + std::to_string(larmorRatios.size()),
        "#LTP_{CREM}/P_{Larmor}#GT = " + compactNumber(larmorMoments.mean, 4)
            + " #pm " + compactNumber(larmorMoments.sigma, 3),
        "solid line at 1: coherent electric dipole",
        "dotted line at 0.5: two charges radiating incoherently"
    }, 0.019);

    // The diagnostics page used to hold four closure histograms of the photon
    // generator.  Those measured the generator's own arithmetic, not the
    // classical model, and for the 2-gamma channel they were identically zero
    // by construction; they now live in positronium_validation.  What belongs
    // here is the calibration behaviour of the CREM trajectories that actually
    // produce the collapse time on the facing page.
    double powerLower = std::numeric_limits<double>::infinity();
    double powerUpper = 0.0;
    for (double value : calibrationPowers) {
        powerLower = std::min(powerLower, value);
        powerUpper = std::max(powerUpper, value);
    }
    if (!(powerLower < powerUpper)) {
        powerLower = 0.0;
        powerUpper = 1.0;
    }
    const double powerPadding = 0.05*(powerUpper - powerLower);
    TH1D calibrationPowerHistogram("crem_calibration_power",
        "Orbit-averaged radiated power of the CREM collapse;"
        "#LTP#GT [W];Trajectories",
        histogramBins(calibrationPowers.size()),
        powerLower - powerPadding, powerUpper + powerPadding);
    styleHistogram(calibrationPowerHistogram, plot_style::crem());
    calibrationPowerHistogram.SetStats(false);
    for (double value : calibrationPowers) calibrationPowerHistogram.Fill(value);
    const GaussianFitSummary powerMoments =
        gaussianMaximumLikelihood(calibrationPowers);

    diagnosticsPage.cd(1);
    gPad->SetGrid();
    calibrationPowerHistogram.Draw("HIST");
    drawAnalysisBox(analysisBoxes, 0.50, 0.58, 0.95, 0.91, {
        "Trajectories with a finite power: N = "
            + std::to_string(calibrationPowers.size()),
        "#LTP#GT = " + compactNumber(powerMoments.mean) + " W",
        "#sigma(P) = " + compactNumber(powerMoments.sigma) + " W",
        "Radiated energy / elapsed time, averaged over each",
        "trajectory's full mechanical run to the boundary",
        "(not extrapolated)."
    }, 0.021);

    // --- Diagnostics pad 2: magnetic sector against the measured splitting ---
    // The classical dipole-dipole coupling of the prepared pair is the only
    // channel by which this model distinguishes para from ortho.  The measured
    // o-Ps/p-Ps hyperfine splitting is what that channel would have to
    // reproduce.  The gap is the point: the real splitting is dominated by
    // virtual annihilation and the Fermi contact term, neither of which a
    // classical dipole model contains.
    diagnosticsPage.cd(2);
    gPad->SetGrid();
    double couplingLower = std::numeric_limits<double>::infinity();
    double couplingUpper = -std::numeric_limits<double>::infinity();
    for (double value : dipoleCouplingsGHz) {
        couplingLower = std::min(couplingLower, std::abs(value));
        couplingUpper = std::max(couplingUpper, std::abs(value));
    }
    if (!(couplingLower < couplingUpper)) { couplingLower = 0.0; couplingUpper = 1.0; }
    const double couplingPadding = 0.08*(couplingUpper - couplingLower) + 1.0e-9;
    TH1D dipoleCouplingHistogram("crem_dipole_coupling",
        "Classical dipole-dipole coupling vs measured hyperfine splitting;"
        "|U_{dd}| / h [GHz];Trajectories",
        histogramBins(dipoleCouplingsGHz.size()),
        std::max(0.0, couplingLower - couplingPadding),
        couplingUpper + couplingPadding);
    styleHistogram(dipoleCouplingHistogram, plot_style::crem());
    dipoleCouplingHistogram.SetStats(false);
    for (double value : dipoleCouplingsGHz)
        dipoleCouplingHistogram.Fill(std::abs(value));
    dipoleCouplingHistogram.Draw("HIST");
    std::vector<double> couplingMagnitudes;
    for (double value : dipoleCouplingsGHz)
        couplingMagnitudes.push_back(std::abs(value));
    const GaussianFitSummary couplingMoments =
        gaussianMaximumLikelihood(couplingMagnitudes);
    // Measured o-Ps/p-Ps splitting, 203.394 GHz.  Far off this axis, so it is
    // quoted rather than drawn: forcing it into range would flatten the
    // model's own distribution onto a single bin.
    constexpr double hyperfineSplittingGHz = 203.3941;
    const double couplingFraction = couplingMoments.mean > 0.0
        ? 100.0*couplingMoments.mean/hyperfineSplittingGHz : 0.0;
    drawAnalysisBox(analysisBoxes, 0.42, 0.55, 0.95, 0.91, {
        plot_style::key(true, false, true, false),
        "trajectories: N = " + std::to_string(couplingMagnitudes.size()),
        "#LT|U_{dd}|/h#GT = " + compactNumber(couplingMoments.mean, 4) + " GHz",
        "measured o-Ps/p-Ps splitting = 203.3941 GHz",
        "classical dipolar term covers "
            + compactNumber(couplingFraction, 3) + "% of it",
        "The remainder is virtual annihilation and the Fermi",
        "contact term. A classical point-dipole model has",
        "neither, so this gap is structural, not a fit residual."
    }, 0.020);

    canvas.cd();
    distributionsPage.Pop();
    canvas.Modified();
    canvas.Update();
    std::vector<root_export::NamedPad> plotsToSave{
        {distributionsPage.GetPad(1), 1, 1, "crem_collapse_time"},
        {distributionsPage.GetPad(2), 1, 2, "annihilation_time"},
        {distributionsPage.GetPad(3), 1, 3, "collapse_time_vs_theory"},
        {distributionsPage.GetPad(4), 1, 4, "radiated_power_vs_larmor"},
        {diagnosticsPage.GetPad(1), 2, 1, "diagnostic_calibration_power"},
        {diagnosticsPage.GetPad(2), 2, 2, "dipole_coupling_vs_hyperfine"}
    };
    reportExports(root_export::saveStatisticalPlots(
        selectedPhenomenon, plotsToSave));

    bool persistenceOk = reportArchiveOperation(
        statistics_archive::writeScientificReferencesText(),
        "scientific-reference catalogue");
    return persistenceOk ? 0 : 3;
}

enum class BeamOutcome {
    Escaped, ShortRange, Captured, Unresolved, NumericalFailure
};

double relativeConservativeParticleEnergy(const State& state);

struct BeamConfiguration {
    double centreOfMassKineticEnergy;
    double impactParameterMaximum;
    double matchingRadius;
    double maximumFlightTime;
    double maximumAsymptoticBeta;
    double analysisThetaMinimum;
    double coulombLength;
    double analyticCutoffImpactParameter;
    int angleBins;
    bool shortRangeFocus;
};

struct EndpointDiagnostics {
    Vec3 noetherMomentum;
    Vec3 noetherAngularMomentum;
    Vec3 radiatedMomentum;
    Vec3 radiatedAngularMomentum;
    Vec3 boundFieldMomentum;
    Vec3 boundFieldAngularMomentum;
    double canonicalMomentumScale = 0.0;
    double mechanicalEnergy = 0.0;
    double radiatedEnergy = 0.0;
    double boundFieldEnergy = 0.0;
    double schottEnergy = 0.0;
    double time = 0.0;
    double radius = 0.0;
    // Independent measurement of the mismatch between the individual
    // Landau-Lifshitz reaction and the coherent far-field flux.  Unlike the
    // bound-field reservoir above this is NOT constructed to close, so it is
    // the only genuine test of the radiation sector in this program.
    double reactionEnergyMismatch = 0.0;
    Vec3 reactionMomentumMismatch;
    Vec3 reactionAngularMomentumMismatch;
};

struct BeamEvent {
    BeamOutcome outcome = BeamOutcome::NumericalFailure;
    std::uint64_t eventSeed = 0;
    double impactParameter = 0.0;
    double impactAzimuth = std::numeric_limits<double>::quiet_NaN();
    double scatteringAngle = std::numeric_limits<double>::quiet_NaN();
    double outgoingEnergy = std::numeric_limits<double>::quiet_NaN();
    double radiatedEnergy = std::numeric_limits<double>::quiet_NaN();
    State initialState;
    State terminalState;
    bool initialStateValid = false;
    bool terminalStateValid = false;
    std::uint64_t integrationSteps = 0;
    double minimumSeparation = std::numeric_limits<double>::quiet_NaN();
    bool minimumSeparationValid = false;
    EndpointDiagnostics initialDiagnostics;
    EndpointDiagnostics finalDiagnostics;
    bool diagnosticsValid = false;

    BeamEvent() = default;
    BeamEvent(BeamOutcome outcomeValue, std::uint64_t seedValue,
              double impactParameterValue,
              double scatteringAngleValue = std::numeric_limits<double>::quiet_NaN(),
              double outgoingEnergyValue = std::numeric_limits<double>::quiet_NaN(),
              double radiatedEnergyValue = std::numeric_limits<double>::quiet_NaN())
        : outcome(outcomeValue), eventSeed(seedValue),
          impactParameter(impactParameterValue),
          scatteringAngle(scatteringAngleValue), outgoingEnergy(outgoingEnergyValue),
          radiatedEnergy(radiatedEnergyValue) {}
};

EndpointDiagnostics endpointDiagnostics(const State& state) {
    const Frame frame = makeFrame(state);
    return {frame.noetherMomentum, frame.noetherAngularMomentum,
            frame.radiatedMomentum, frame.radiatedAngularMomentum,
            frame.boundFieldMomentum,frame.boundFieldAngularMomentum,
            frame.canonicalMomentumScale, frame.mechanicalEnergy,
            frame.radiatedEnergy,frame.boundFieldEnergy,frame.schottEnergy,
            frame.time, frame.radius,
            frame.reactionEnergyMismatch,frame.reactionMomentumMismatch,
            frame.reactionAngularMomentumMismatch};
}

BeamConfiguration makeBeamConfiguration(int selectedPhenomenon,
                                        double energyEv,
                                        double thetaMinimumDegrees,
                                        int angleBins,
                                        double impactParameterMaximumPm,
                                        double matchingRadiusPm) {
    if (!(energyEv > 0.0) || !std::isfinite(energyEv)) {
        throw std::invalid_argument("beam energy must be finite and positive");
    }
    if (!(thetaMinimumDegrees > 0.0 && thetaMinimumDegrees < 180.0)
        || !std::isfinite(thetaMinimumDegrees)) {
        throw std::invalid_argument("theta-min must lie between 0 and 180 degrees");
    }
    if (angleBins < 4 || angleBins > 90) {
        throw std::invalid_argument("angle-bins must lie between 4 and 90");
    }
    if (!std::isfinite(impactParameterMaximumPm) || impactParameterMaximumPm < 0.0) {
        throw std::invalid_argument("bmax must be finite and non-negative (zero selects auto)");
    }
    if (!std::isfinite(matchingRadiusPm) || matchingRadiusPm < 0.0) {
        throw std::invalid_argument(
            "matching radius must be finite and non-negative (zero selects auto)");
    }

    const bool shortRangeFocus = selectedPhenomenon == 3;
    const double energy = energyEv * eCharge;
    const double coulombStrength = pairCoulombStrength;
    const double coulombLength = coulombStrength / (2.0 * energy);
    const double cutoffImpactParameter = std::sqrt(
        nuclearCutoff * nuclearCutoff + coulombStrength * nuclearCutoff / energy);
    const double requestedTheta = thetaMinimumDegrees * pi / 180.0;
    const double defaultImpactMaximum = shortRangeFocus
        ? 12.0 * cutoffImpactParameter
        : 1.25 * coulombLength / std::tan(0.5 * requestedTheta);
    const double impactMaximum = impactParameterMaximumPm > 0.0
        ? impactParameterMaximumPm * 1.0e-12 : defaultImpactMaximum;
    if (!(impactMaximum > 0.0) || !std::isfinite(impactMaximum)) {
        throw std::invalid_argument("bmax must be finite and positive");
    }
    const double samplingThetaMinimum = 2.0 * std::atan(coulombLength / impactMaximum);
    const double analysisThetaMinimum = std::max(requestedTheta, samplingThetaMinimum);
    const double automaticMatchingRadius = std::max(
        100.0 * impactMaximum, 1000.0 * coulombStrength / energy);
    const double matchingRadius = matchingRadiusPm > 0.0
        ? matchingRadiusPm * 1.0e-12 : automaticMatchingRadius;
    if (!(matchingRadius > impactMaximum) || !std::isfinite(matchingRadius)) {
        throw std::invalid_argument("matching radius must be finite and larger than bmax");
    }

    const two_body::CentreOfMomentumKinematics asymptotic =
        two_body::centreOfMomentumKinematics(energy,firstMass,secondMass);
    if (!asymptotic.valid() || !(asymptotic.relativeSpeed > 0.0)) {
        throw std::invalid_argument(
            "beam energy cannot be represented for the selected pair");
    }
    return {energy, impactMaximum, matchingRadius,
            8.0 * matchingRadius / asymptotic.relativeSpeed,
            std::max(asymptotic.firstSpeed,asymptotic.secondSpeed)/c,
            analysisThetaMinimum, coulombLength, cutoffImpactParameter,
            angleBins, shortRangeFocus};
}

BeamEvent simulateBeamEvent(
    std::uint64_t seed,const BeamConfiguration& configuration,
    ClassicalTrajectoryEngine::Accuracy accuracy) {
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double impactParameter = configuration.impactParameterMaximum
        * std::sqrt(uniform(random));
    const double azimuth = 2.0 * pi * uniform(random);
    const Vec3 impactDirection{0.0, std::cos(azimuth), std::sin(azimuth)};
    const Vec3 beamDirection{1.0, 0.0, 0.0};
    const double longitudinalDistance = std::sqrt(
        configuration.matchingRadius * configuration.matchingRadius
        - impactParameter * impactParameter);
    const Vec3 relativePosition = beamDirection * (-longitudinalDistance)
                                + impactDirection * impactParameter;
    const Vec3 radialDirection = relativePosition / configuration.matchingRadius;
    const double sineAsymptote = impactParameter / configuration.matchingRadius;
    const Vec3 tangentDirection = beamDirection * sineAsymptote
        + impactDirection * (longitudinalDistance / configuration.matchingRadius);

    const double coulombStrength = pairCoulombStrength;
    const two_body::IncomingTwoBodyKinematics incoming =
        two_body::incomingTwoBodyKinematics(
            configuration.centreOfMassKineticEnergy,
            coulombStrength/configuration.matchingRadius,
            impactParameter,configuration.matchingRadius,
            radialDirection,tangentDirection,
            firstMass,secondMass);
    if (!incoming.valid()) {
        BeamEvent result{BeamOutcome::NumericalFailure, seed, impactParameter};
        result.impactAzimuth = azimuth;
        return result;
    }

    const auto randomDirection = [&]() {
        const double cosine = 2.0 * uniform(random) - 1.0;
        const double phi = 2.0 * pi * uniform(random);
        const double transverse = std::sqrt(std::max(0.0, 1.0 - cosine*cosine));
        return Vec3{transverse * std::cos(phi), transverse * std::sin(phi), cosine};
    };
    State state;
    // The beam experiments are initialized directly in the pair's COM frame.
    // Equal and opposite momenta, rather than a mass-weighted split of one
    // relative velocity, keep the total momentum exactly zero for unequal
    // masses.  Energy weights put the free-particle centre of energy at the
    // origin; they reduce to the familiar mass weights at low velocity.
    const double finiteEnergy = incoming.finiteRadius.firstEnergy
                              + incoming.finiteRadius.secondEnergy;
    state.firstPosition =
        relativePosition * (incoming.finiteRadius.secondEnergy/finiteEnergy);
    state.secondPosition =
        relativePosition * (-incoming.finiteRadius.firstEnergy/finiteEnergy);
    state.firstVelocity = incoming.firstVelocity;
    state.secondVelocity = incoming.secondVelocity;
    state.firstDipole = randomDirection() * firstMagneticMoment;
    state.secondDipole = randomDirection() * secondMagneticMoment;
    ClassicalTrajectoryEngine trajectory(state,accuracy);
    const State initialState = state;
    const EndpointDiagnostics initialDiagnostics = endpointDiagnostics(state);
    std::uint64_t integrationSteps = 0;
    double minimumSeparation = separation(state);
    const auto makeResult = [&](BeamOutcome outcome,
                                double scatteringAngle,
                                double outgoingEnergy,
                                double emittedEnergy,
                                const State* terminalState) {
        BeamEvent result{outcome, seed, impactParameter, scatteringAngle,
                         outgoingEnergy, emittedEnergy};
        result.impactAzimuth = azimuth;
        result.initialState = initialState;
        result.initialStateValid = true;
        result.integrationSteps = integrationSteps;
        result.minimumSeparation = minimumSeparation;
        result.minimumSeparationValid = std::isfinite(minimumSeparation);
        result.initialDiagnostics = initialDiagnostics;
        if (terminalState && isFinite(*terminalState)) {
            result.terminalState = *terminalState;
            result.terminalStateValid = true;
        }
        return result;
    };

    const double reducedMass = firstMass * secondMass
                                 / (firstMass + secondMass);
    bool passedClosestApproach = false;
    while (state.time < configuration.maximumFlightTime) {
        const double radius = separation(state);
        const double relativeSpeed =
            (state.firstVelocity - state.secondVelocity).norm();
        const double omega = std::sqrt(pairCoulombStrength
                                      / (reducedMass * radius*radius*radius));
        const double transitStep = 0.02 * radius/std::max(relativeSpeed, 1.0);
        // Far from the interaction region the transit and orbital scales are
        // much longer than 2 as. A global attosecond cap forced tens of
        // thousands of needless steps before the particles reached one
        // another. The adaptive trajectory engine still subdivides this
        // ceiling whenever its local-error estimate requires it.
        const double dt = std::min({1.0e-15, 2.0*pi/(320.0*omega), transitStep,
            configuration.maximumFlightTime - state.time});
        if (!(dt > 0.0) || !std::isfinite(dt)) {
            return makeResult(BeamOutcome::NumericalFailure,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                state.radiatedEnergy, &state);
        }
        const State beforeStep = state;
        ++integrationSteps;
        if(!trajectory.advance(state,dt)) {
            return makeResult(BeamOutcome::NumericalFailure,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                beforeStep.radiatedEnergy, &beforeStep);
        }
        const Vec3 currentRelativePosition =
            state.firstPosition - state.secondPosition;
        const Vec3 currentRelativeVelocity =
            state.firstVelocity - state.secondVelocity;
        const double currentRadius = currentRelativePosition.norm();
        if (!(currentRadius > 0.0) || !std::isfinite(currentRadius)) {
            return makeResult(BeamOutcome::NumericalFailure,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                beforeStep.radiatedEnergy, &beforeStep);
        }
        minimumSeparation = std::min(minimumSeparation, currentRadius);
        if (currentRadius <= nuclearCutoff) {
            const double crossingFraction=separationCrossingFraction(
                beforeStep,state,nuclearCutoff);
            State cutoffState=interpolateState(
                beforeStep,state,crossingFraction);
            minimumSeparation=std::min(minimumSeparation,nuclearCutoff);
            return makeResult(BeamOutcome::ShortRange,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                cutoffState.radiatedEnergy,&cutoffState);
        }
        if (dot(currentRelativePosition, currentRelativeVelocity) > 0.0) {
            passedClosestApproach = true;
        }
        // Once an outgoing branch has negative conservative particle energy,
        // it is a captured/bound trajectory and cannot reach the matching
        // sphere without a later external energy source. Do not integrate its
        // many tiny bound-orbit periods all the way to the wall-clock flight
        // gate; report it as a distinct captured channel.
        if(passedClosestApproach&&relativeConservativeParticleEnergy(state)<0.0) {
            return makeResult(BeamOutcome::Captured,
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(),
                state.radiatedEnergy,&state);
        }
        if (passedClosestApproach && currentRadius >= configuration.matchingRadius) {
            const double crossingFraction=separationCrossingFraction(
                beforeStep,state,configuration.matchingRadius);
            const State diagnosticEndpoint=interpolateState(
                beforeStep,state,crossingFraction);
            const Vec3 endpointRelativeVelocity=
                diagnosticEndpoint.firstVelocity
                -diagnosticEndpoint.secondVelocity;
            const double endpointRelativeSpeed=endpointRelativeVelocity.norm();
            if (!(endpointRelativeSpeed > 0.0)
                || !std::isfinite(endpointRelativeSpeed)) {
                return makeResult(BeamOutcome::NumericalFailure,
                    std::numeric_limits<double>::quiet_NaN(),
                    std::numeric_limits<double>::quiet_NaN(),
                    diagnosticEndpoint.radiatedEnergy,&diagnosticEndpoint);
            }
            const double angle=std::acos(std::clamp(
                dot(beamDirection,endpointRelativeVelocity)
                    /endpointRelativeSpeed,-1.0,1.0));
            // At the finite matching sphere the kinetic energy still contains
            // the Coulomb acceleration.  Evaluate every terminal observable at
            // the same interpolated crossing event; using the post-step state
            // mixed two different radii in one reported result.
            const double outgoingEnergy=
                relativeConservativeParticleEnergy(diagnosticEndpoint);
            BeamEvent result = makeResult(BeamOutcome::Escaped, angle,
                outgoingEnergy,diagnosticEndpoint.radiatedEnergy,
                &diagnosticEndpoint);
            result.finalDiagnostics = endpointDiagnostics(diagnosticEndpoint);
            result.diagnosticsValid = true;
            return result;
        }
    }
    return makeResult(BeamOutcome::Unresolved,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(), state.radiatedEnergy, &state);
}

// ---------------------------------------------------------------------------
// Experiment 5 "Interactions": an first and a second are fired at each
// other with a Gaussian centre-of-mass energy and a Gaussian impact parameter
// centred on a head-on collision, both dipoles randomly oriented in space.
// Every trajectory is classified by what actually happens to it.
// ---------------------------------------------------------------------------

enum class InteractionOutcome {
    Collision, Scattering, ParaPositronium, OrthoPositronium,
    Unresolved, NumericalFailure
};

const char* interactionOutcomeName(InteractionOutcome outcome) {
    switch (outcome) {
        case InteractionOutcome::Collision:        return "Collision";
        case InteractionOutcome::Scattering:       return "Scattering";
        case InteractionOutcome::ParaPositronium:  return "Para-Positronium";
        case InteractionOutcome::OrthoPositronium: return "Ortho-Positronium";
        case InteractionOutcome::Unresolved:       return "Unresolved";
        case InteractionOutcome::NumericalFailure: return "NumericalFailure";
    }
    return "Unknown";
}

struct InteractionConfiguration {
    double meanKineticEnergy = 0.0;      // per-particle laboratory energy [J]
    double kineticEnergySigma = 0.0;     // [J]
    double minimumKineticEnergy = 0.0;   // truncation of the Gaussian [J]
    double impactParameterSigma = 0.0;   // [m], distribution centred on b=0
    double matchingRadius = 0.0;         // [m], initial and escape separation
    double maximumFlightTime = 0.0;      // [s]
    // The bound phase is measured in orbits, not in absolute time: the orbital
    // period scales as a^{3/2} and the adaptive step scales with it, so a fixed
    // time window costs a bounded number of steps for a loose capture and tens
    // of thousands for a tight one.  Orbits keep the cost per event flat.
    double boundObservationOrbits = 0.0;
    double boundObservationTimeCap = 0.0; // [s] guard for very loose captures
    // Wall-clock budget per trajectory.  The per-event cost is extremely
    // heavy-tailed: a trajectory that turns around inside the classical
    // dipole barrier near 0.2 pm hits a stiff region where single adaptive
    // steps cost about a second, so a handful of events would otherwise run
    // for hours.  Exceeding the budget censors the event as Unresolved, the
    // same censoring convention the beam experiments use for their flight gate.
    double eventWallClockBudgetSeconds = 0.0;
    // Separation at which the trajectory is declared to have collided.  This is
    // the model's own resolution limit, chargeCloudRestRadius = 0.01*a0 =
    // 529 fm, not the point-particle boundary nuclearCutoff = 10 fm.
    //
    // The point-particle boundary sits BELOW the classical dipole barrier at
    // r* = 193 fm, where the dipole and Coulomb energies cross, so no
    // trajectory could ever reach it: they turned around at 170-450 fm and
    // ground away there, needing steps a thousand times shorter than normal
    // until the wall-clock budget censored them.  The program was spending its
    // whole budget resolving a regime the model cannot describe anyway, since
    // 193 fm is half the reduced Compton wavelength and the charge sector
    // already assumes a source of radius 529 fm.  Terminating at the model's
    // own limit ends those trajectories promptly with a meaningful label.
    double collisionRadius = 0.0;
};

struct InteractionEvent {
    InteractionOutcome outcome = InteractionOutcome::NumericalFailure;
    std::uint64_t eventSeed = 0;
    double kineticEnergyEv = std::numeric_limits<double>::quiet_NaN();
    double impactParameter = std::numeric_limits<double>::quiet_NaN();
    double minimumSeparation = std::numeric_limits<double>::quiet_NaN();
    // Mean and spread of cos(mu_e, mu_p) sampled over the bound phase.
    double dipoleAlignment = std::numeric_limits<double>::quiet_NaN();
    double dipoleAlignmentSpread = std::numeric_limits<double>::quiet_NaN();
    double radiatedEnergyEv = std::numeric_limits<double>::quiet_NaN();
    double finalRelativeEnergyEv = std::numeric_limits<double>::quiet_NaN();
    double scatteringAngleDegrees = std::numeric_limits<double>::quiet_NaN();
    double elapsedTime = std::numeric_limits<double>::quiet_NaN();
    // Extrapolated CREM collapse time of a captured pair from the secular
    // model dE/dt = P, with E = -k e^2/(2a) and P proportional to a^-4.  It is
    // a classical inspiral time, not a quantum annihilation lifetime.
    //
    // This is NOT the same measurement experiments 1 and 2 report.  Those now
    // integrate the trajectory mechanically to the boundary; this one still
    // extrapolates from a measured power, and it only exists where the
    // measurement window closed a full orbit (see boundObservedOrbits), which
    // is the precondition for "orbit-averaged P" to mean anything.  The two
    // numbers describe different orbits by different methods and should not
    // be pooled.
    double collapseTimeSeconds = std::numeric_limits<double>::quiet_NaN();
    double boundRadiatedPowerWatts = std::numeric_limits<double>::quiet_NaN();
    // Kepler period of the captured quasi-closed orbit, taken as the longer of
    // the capture-time and window-close periods -- the same quantity the
    // orbit gate below divides by.
    double boundOrbitalPeriodSeconds = std::numeric_limits<double>::quiet_NaN();
    // Length of the bound observation window in units of the capture-time
    // Kepler period.  Below 1 no orbit average exists and collapseTimeSeconds
    // stays NaN.
    double boundObservedOrbits = std::numeric_limits<double>::quiet_NaN();
    EndpointDiagnostics initialDiagnostics;
    EndpointDiagnostics finalDiagnostics;
    bool diagnosticsValid = false;
};

// Kinetic plus Coulomb energy of the pair, deliberately excluding the
// short-range dipole-dipole term.  conservativeParticleEnergy() cannot be used
// to decide capture here: at the sub-picometre turning points the classical
// dipole interaction reaches 1e15 eV and swamps the Coulomb binding, so it
// reports "bound" with an energy that makes the Kepler relations meaningless
// (a = -k/2E collapses to 1e-25 m).  The secular collapse model is a Coulomb
// model, and experiments 1 and 2 calibrate it on exactly this energy.
double coulombPairEnergy(const State& state) {
    const PairGeometry geometry = pairGeometry(state);
    const two_body::FreeTwoBodyKinematics freePair =
        two_body::freeTwoBodyKinematics(
            state.firstVelocity,firstMass,
            state.secondVelocity,secondMass);
    if (!freePair.valid()) return std::numeric_limits<double>::quiet_NaN();
    return freePair.centreOfMomentumKineticEnergy
         - pairCoulombStrength*geometry.inverseDistance;
}

// Conservative model energy with the free translational kinetic energy of the
// pair removed.  This is used by the COM-initialized beam experiments 3/4.
// Experiment 5 deliberately uses coulombPairEnergy() instead, because its
// capture and secular-collapse branches are both Coulomb/Kepler models and
// must not let the short-range dipole term redefine their binding energy.
double relativeConservativeParticleEnergy(const State& state) {
    const two_body::FreeTwoBodyKinematics freePair =
        two_body::freeTwoBodyKinematics(
            state.firstVelocity,firstMass,
            state.secondVelocity,secondMass);
    if (!freePair.valid()) return std::numeric_limits<double>::quiet_NaN();
    return conservativeParticleEnergy(state)
         - freePair.laboratoryKineticEnergy
         + freePair.centreOfMomentumKineticEnergy;
}

double dipoleAlignmentOf(const State& state) {
    const double firstNorm = state.firstDipole.norm();
    const double secondNorm = state.secondDipole.norm();
    if (!(firstNorm > 0.0) || !(secondNorm > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::clamp(dot(state.firstDipole, state.secondDipole)
                      / (firstNorm*secondNorm), -1.0, 1.0);
}

InteractionConfiguration makeInteractionConfiguration(
    double meanEnergyEv, double energySigmaEv, double impactSigmaPm) {
    if (!(meanEnergyEv > 0.0) || !std::isfinite(meanEnergyEv)) {
        throw std::invalid_argument(
            "interaction mean energy must be finite and positive");
    }
    if (!(energySigmaEv >= 0.0) || !std::isfinite(energySigmaEv)) {
        throw std::invalid_argument(
            "interaction energy sigma must be finite and non-negative");
    }
    if (!(impactSigmaPm >= 0.0) || !std::isfinite(impactSigmaPm)) {
        throw std::invalid_argument(
            "interaction impact-parameter sigma must be finite and non-negative"
            " (zero selects the automatic Coulomb length)");
    }
    InteractionConfiguration configuration;
    configuration.meanKineticEnergy = meanEnergyEv*eCharge;
    configuration.kineticEnergySigma = energySigmaEv*eCharge;
    // Truncate well below the mean but strictly positive: a particle with
    // non-positive laboratory kinetic energy has no incoming beam asymptote.
    // Never put the cutoff above the requested mean.  In particular, a
    // zero-width beam below 1 meV must use its actual speed when the COM flight
    // window is sized, rather than a fictitious faster 1 meV particle.
    configuration.minimumKineticEnergy=std::min(
        configuration.meanKineticEnergy,
        std::max(0.02*configuration.meanKineticEnergy,1.0e-3*eCharge));
    // The scale that decides capture versus fly-past is the Coulomb length
    // l_C = k|q1 q2|/(2 K_CM).  meanKineticEnergy is a PER-PARTICLE lab
    // quantity, so first convert two representative opposing beams to their
    // invariant K_CM; using the lab mean directly is wrong even for equal
    // masses (where K_CM is approximately twice as large), and has no fixed
    // correction for unequal masses.
    const two_body::HeadOnLabKinematics representativeBeam =
        two_body::headOnLabKinematics(
            configuration.meanKineticEnergy,firstMass,
            configuration.meanKineticEnergy,secondMass,Vec3{1.0,0.0,0.0});
    if (!representativeBeam.valid()
        || !(representativeBeam.pair.centreOfMomentumKineticEnergy > 0.0)) {
        throw std::invalid_argument(
            "interaction energy cannot be represented for the selected pair");
    }
    const double coulombLength = pairCoulombStrength
        /(2.0*representativeBeam.pair.centreOfMomentumKineticEnergy);
    configuration.impactParameterSigma = impactSigmaPm > 0.0
        ? impactSigmaPm*1.0e-12 : coulombLength;
    // Start far enough out that the Coulomb potential is a small correction to
    // the sampled kinetic energy, but close enough that the approach phase does
    // not dominate the run time.
    //
    // Expressed in the PAIR's Bohr radius, not hydrogen's.  50*a0 was 25 pair
    // radii for positronium and nothing in particular for anyone else: for
    // mu+mu- it would have started the pair 5000 atom-widths out and spent the
    // whole run on the approach.
    //
    // Note that the criterion in the comment above is really about the Coulomb
    // length l_C = k|q1 q2|/(2K), which is computed a few lines up, and NOT
    // about the Bohr radius at all -- the Bohr radius is here only because the
    // experiment was first written for positronium.  At the default 0.6 eV,
    // l_C is 1200 pm against a start separation of 2646 pm, so "small
    // correction" holds by a factor of 2.2 and no more.  Retuning that is a
    // separate decision about experiment 5's defaults, not part of making the
    // constant pair-relative, so it is left alone and written down instead.
    // Keep the half-normal beam well inside the matching sphere for every
    // species.  A fixed 25 pair radii is only picometres for true muonium and
    // protonium; clipping a 60 pm beam at R/2 then piles almost every event on
    // one artificial impact parameter.  Five sigma fit below R/2 here.
    configuration.matchingRadius = std::max(
        25.0*pairBohrRadius(activePair),
        10.0*configuration.impactParameterSigma);
    const double slowestEnergy=configuration.minimumKineticEnergy;
    const two_body::HeadOnLabKinematics slowestLabPair=
        two_body::headOnLabKinematics(
            slowestEnergy,firstMass,slowestEnergy,secondMass,
            Vec3{1.0,0.0,0.0});
    const two_body::CentreOfMomentumKinematics slowestComPair=
        two_body::centreOfMomentumKinematics(
            slowestLabPair.pair.centreOfMomentumKineticEnergy,
            firstMass,secondMass);
    const double slowestClosingSpeed=slowestComPair.relativeSpeed;
    if (!slowestLabPair.valid() || !slowestComPair.valid()
        || !(slowestClosingSpeed > 0.0)
        || !std::isfinite(slowestClosingSpeed)) {
        throw std::invalid_argument(
            "interaction energy cannot be represented for the selected pair");
    }
    configuration.maximumFlightTime = 6.0*configuration.matchingRadius
                                    / slowestClosingSpeed;
    // Long enough to average the dipole orientation over several bound orbits.
    configuration.boundObservationOrbits = 8.0;
    configuration.boundObservationTimeCap = 2.0e-16;
    configuration.eventWallClockBudgetSeconds = 20.0;
    configuration.collisionRadius = collisionBoundaryRadius;
    return configuration;
}

InteractionEvent simulateInteractionEvent(
    std::uint64_t seed, const InteractionConfiguration& configuration,
    ClassicalTrajectoryEngine::Accuracy accuracy) {
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> energyGaussian(
        configuration.meanKineticEnergy, configuration.kineticEnergySigma);
    std::normal_distribution<double> impactGaussian(
        0.0, configuration.impactParameterSigma);

    InteractionEvent result;
    result.eventSeed = seed;
    // Each particle's kinetic energy is drawn on its own, so the two beams
    // carry independent spreads and the pair generally arrives with a NET
    // momentum.  That is the physical content of the change: for a two-body
    // system whose centre of mass is at rest the two energies cannot be
    // independent, because zero total momentum forces |p1| = |p2| and fixes
    // the second energy once the first and the two masses are known.  Drawing
    // them separately therefore necessarily sets the centre of mass moving.
    //
    // configuration.meanKineticEnergy is consequently a PER-PARTICLE lab
    // energy now, not the centre-of-mass energy it used to be.  The relative
    // energy the collision actually happens at is derived below and is roughly
    // twice as large for an equal-mass pair, so this moves the experiment's
    // operating point; both numbers are printed at startup.
    const auto sampleKinetic = [&]() {
        if (!(configuration.kineticEnergySigma > 0.0))
            return configuration.meanKineticEnergy;
        for (int attempt = 0; attempt < 1000; ++attempt) {
            const double sampled = energyGaussian(random);
            if (sampled >= configuration.minimumKineticEnergy) return sampled;
        }
        return configuration.minimumKineticEnergy;
    };
    const double firstKineticEnergy = sampleKinetic();
    const double secondKineticEnergy = sampleKinetic();
    // Gaussian around a head-on collision: the signed sample is folded, which
    // is the half-normal radial profile of a beam centred on b=0.  Reject the
    // negligible tail outside the matching sphere instead of clipping it;
    // clipping creates a delta-like pile-up at exactly R/2 for a narrow atomic
    // scale or a user-selected broad beam.
    double impactParameter = std::numeric_limits<double>::quiet_NaN();
    for (int attempt=0;attempt<1000;++attempt) {
        const double sampledImpact=std::abs(impactGaussian(random));
        if (sampledImpact <= 0.5*configuration.matchingRadius) {
            impactParameter=sampledImpact;
            break;
        }
    }
    if (!std::isfinite(impactParameter)) return result;
    result.impactParameter = impactParameter;

    const double azimuth = 2.0*pi*uniform(random);
    const Vec3 impactDirection{0.0, std::cos(azimuth), std::sin(azimuth)};
    const Vec3 beamDirection{1.0, 0.0, 0.0};
    const double longitudinalDistance = std::sqrt(
        configuration.matchingRadius*configuration.matchingRadius
        - impactParameter*impactParameter);
    const Vec3 relativePosition = beamDirection*(-longitudinalDistance)
                                + impactDirection*impactParameter;
    const Vec3 radialDirection = relativePosition/configuration.matchingRadius;
    const Vec3 tangentDirection = beamDirection
            *(impactParameter/configuration.matchingRadius)
        + impactDirection*(longitudinalDistance/configuration.matchingRadius);

    // Head-on lab beams: the two momenta point at each other along the beam
    // axis and do NOT cancel, because the energies were drawn independently.
    // Construct their exact four-momenta first; the invariant collision energy
    // and COM velocity then follow without subtracting nearly equal rest-energy
    // squares for a light-heavy pair.
    const two_body::HeadOnLabKinematics labIncoming =
        two_body::headOnLabKinematics(
            firstKineticEnergy,firstMass,
            secondKineticEnergy,secondMass,beamDirection);
    if (!labIncoming.valid()) return result;
    // The energy the collision actually happens at, which is what every
    // threshold and axis in this experiment means by "collision energy".
    const double kineticEnergy =
        labIncoming.pair.centreOfMomentumKineticEnergy;
    if (!(kineticEnergy > 0.0)) return result;
    result.kineticEnergyEv = kineticEnergy/eCharge;

    // Transform the independently sampled lab beams to their COM invariant,
    // then define and integrate the matching sphere in that COM frame.  A
    // lab-t=const sphere is not simultaneous in COM, while the Coulomb matching
    // potential, Kepler clock and capture energy below are all internal COM
    // quantities.  Keeping one frame for the entire trajectory avoids mixing
    // those quantities with a boosted laboratory geometry.
    const double coulombStrength = pairCoulombStrength;
    const two_body::IncomingTwoBodyKinematics incoming =
        two_body::incomingTwoBodyKinematics(
            kineticEnergy,coulombStrength/configuration.matchingRadius,
            impactParameter,configuration.matchingRadius,
            radialDirection,tangentDirection,firstMass,secondMass);
    if (!incoming.valid()) return result;

    const auto randomDirection = [&]() {
        const double cosine = 2.0*uniform(random) - 1.0;
        const double phi = 2.0*pi*uniform(random);
        const double transverse = std::sqrt(std::max(0.0, 1.0 - cosine*cosine));
        return Vec3{transverse*std::cos(phi), transverse*std::sin(phi), cosine};
    };
    State state;
    // Energy weights place the COM centre of energy at the origin while
    // preserving the requested relative position.  They reduce to the usual
    // mass weights in the nonrelativistic limit.
    const double finiteComEnergy=incoming.finiteRadius.firstEnergy
                                +incoming.finiteRadius.secondEnergy;
    state.firstPosition = relativePosition
        *(incoming.finiteRadius.secondEnergy/finiteComEnergy);
    state.secondPosition = relativePosition
        *(-incoming.finiteRadius.firstEnergy/finiteComEnergy);
    state.firstVelocity = incoming.firstVelocity;
    state.secondVelocity = incoming.secondVelocity;
    state.firstDipole = randomDirection()*firstMagneticMoment;
    state.secondDipole = randomDirection()*secondMagneticMoment;

    ClassicalTrajectoryEngine trajectory(state, accuracy);
    result.initialDiagnostics = endpointDiagnostics(state);
    double minimumSeparation = separation(state);
    bool passedClosestApproach = false;
    bool bound = false;
    int captureAttempts = 0;
    // Five per cent of THIS event's collision energy.  It used to be five per
    // cent of configuration.meanKineticEnergy, which was the centre-of-mass
    // energy; that field is now a per-particle lab energy, so reading it here
    // would have quietly halved the margin for an equal-mass pair.
    const double captureMargin = 0.05*kineticEnergy;
    double boundStartTime = 0.0;
    // Charge-sector radiated energy, i.e. the total outward flux minus the
    // magnetic-dipole (M1) part.  state.radiatedEnergy accumulates
    // outwardFlux.energy, which includes M1; state.dipoleConstraintEnergy
    // accumulates exactly -integral(P_M1 dt), so the sum is the E1+E2
    // radiation that actually recoils the orbit.  This is the same split
    // chargeMismatchRates() applies in the energy-balance diagnostic.
    double boundStartChargeRadiated = 0.0;
    double boundStartOrbitalPeriod = 0.0;
    double boundObservationTime = configuration.boundObservationTimeCap;
    // Welford accumulation of cos(mu_e, mu_p) over the bound phase.
    std::size_t alignmentCount = 0;
    double alignmentMean = 0.0;
    double alignmentSecondMoment = 0.0;

    const double reducedMass = firstMass*secondMass
                                 / (firstMass + secondMass);
    const auto wallClockStart = std::chrono::steady_clock::now();
    long stepCounter = 0;
    const auto finish = [&](InteractionOutcome outcome, const State& endpoint) {
        result.outcome = outcome;
        result.minimumSeparation = minimumSeparation;
        result.radiatedEnergyEv = endpoint.radiatedEnergy/eCharge;
        result.finalRelativeEnergyEv =
            coulombPairEnergy(endpoint)/eCharge;
        result.elapsedTime = endpoint.time;
        if (alignmentCount > 0) {
            result.dipoleAlignment = alignmentMean;
            result.dipoleAlignmentSpread = alignmentCount > 1
                ? std::sqrt(std::max(0.0, alignmentSecondMoment)
                            / static_cast<double>(alignmentCount))
                : 0.0;
        }
        if (isFinite(endpoint)) {
            result.finalDiagnostics = endpointDiagnostics(endpoint);
            result.diagnosticsValid = true;
        }
        return result;
    };

    while (state.time < configuration.maximumFlightTime) {
        const double radius = separation(state);
        const double relativeSpeed =
            (state.firstVelocity - state.secondVelocity).norm();
        const double omega = std::sqrt(coulombStrength
                                       / (reducedMass*radius*radius*radius));
        const double transitStep = 0.02*radius/std::max(relativeSpeed, 1.0);
        double dt = std::min({1.0e-15, 2.0*pi/(320.0*omega), transitStep,
                              configuration.maximumFlightTime - state.time});
        // A step clamped to the end of a window can fall below the floating
        // point resolution of state.time, after which state.time += dt is a
        // no-op: the clock freezes and the event spins until the wall-clock
        // budget censors it as Unresolved.  Measured on captured pairs with
        // E_Coulomb = -17 eV, so those were successful captures being thrown
        // away by a rounding artefact rather than trajectories needing work.
        const double clockResolution = std::max(std::abs(state.time),
            configuration.boundObservationTimeCap)
            * std::numeric_limits<double>::epsilon() * 8.0;
        if (bound) {
            const double remaining =
                boundStartTime + boundObservationTime - state.time;
            if (remaining <= clockResolution) {
                boundObservationTime = state.time - boundStartTime;
            } else {
                dt = std::min(dt, remaining);
            }
        }
        if (configuration.maximumFlightTime - state.time <= clockResolution) {
            return finish(InteractionOutcome::Unresolved, state);
        }
        if (!(dt > 0.0) || !std::isfinite(dt)) return finish(
            InteractionOutcome::NumericalFailure, state);
        // Checked every 8 steps: inside the dipole barrier a single adaptive
        // step costs about a second, so a coarser interval would let an event
        // overshoot its budget by a minute before the gate noticed.
        if ((++stepCounter & 0x7) == 0
            && configuration.eventWallClockBudgetSeconds > 0.0) {
            const double spent = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wallClockStart).count();
            if (spent > configuration.eventWallClockBudgetSeconds) {
                return finish(InteractionOutcome::Unresolved, state);
            }
        }
        const State beforeStep = state;
        if (!trajectory.advance(state, dt)) return finish(
            InteractionOutcome::NumericalFailure, beforeStep);
        // Belt and braces: any path that fails to advance the clock would
        // otherwise loop until the budget runs out.
        if (!(state.time > beforeStep.time)) {
            if (bound) boundObservationTime = state.time - boundStartTime;
            else return finish(InteractionOutcome::Unresolved, state);
        }

        const Vec3 relative = state.firstPosition - state.secondPosition;
        const double currentRadius = relative.norm();
        if (!(currentRadius > 0.0) || !std::isfinite(currentRadius)) {
            return finish(InteractionOutcome::NumericalFailure, beforeStep);
        }
        minimumSeparation = std::min(minimumSeparation, currentRadius);

        // A trajectory that reaches the point-particle boundary has collided,
        // whether or not it had already been captured.
        if (currentRadius <= configuration.collisionRadius) {
            const double crossingFraction = separationCrossingFraction(
                beforeStep, state, configuration.collisionRadius);
            const State cutoffState = interpolateState(
                beforeStep, state, crossingFraction);
            minimumSeparation = std::min(minimumSeparation,
                                          configuration.collisionRadius);
            return finish(InteractionOutcome::Collision, cutoffState);
        }

        if (dot(relative, state.firstVelocity - state.secondVelocity) > 0.0) {
            passedClosestApproach = true;
        }
        // Require the pair to be bound by a finite margin, not merely to have
        // crossed zero.  Near the turning point the Coulomb energy hovers
        // around zero, and a bare "< 0" test lets a marginal pair enter and
        // leave the bound branch hundreds of times, each cycle paying for a
        // full observation window.  The retry limit is a second backstop.
        if (!bound && passedClosestApproach && captureAttempts < 3) {
            const double boundEnergy = coulombPairEnergy(state);
            if (boundEnergy < -captureMargin) {
                ++captureAttempts;
                bound = true;
                boundStartTime = state.time;
                boundStartChargeRadiated =
                    state.radiatedEnergy + state.dipoleConstraintEnergy;
                // Kepler period of the captured orbit, from its semi-major
                // axis a = -k/(2E).  Averaging over a fixed number of orbits
                // costs a fixed number of adaptive steps at any binding energy.
                const double semiMajorAxis = -coulombStrength/(2.0*boundEnergy);
                const double orbitalPeriod = 2.0*pi*std::sqrt(
                    reducedMass*semiMajorAxis*semiMajorAxis*semiMajorAxis
                    / coulombStrength);
                boundStartOrbitalPeriod = std::isfinite(orbitalPeriod)
                        && orbitalPeriod > 0.0 ? orbitalPeriod : 0.0;
                boundObservationTime = boundStartOrbitalPeriod > 0.0
                    ? std::min(configuration.boundObservationOrbits
                                   *orbitalPeriod,
                               configuration.boundObservationTimeCap)
                    : configuration.boundObservationTimeCap;
            }
        }
        if (bound) {
            const double alignment = dipoleAlignmentOf(state);
            if (std::isfinite(alignment)) {
                ++alignmentCount;
                const double delta = alignment - alignmentMean;
                alignmentMean += delta/static_cast<double>(alignmentCount);
                alignmentSecondMoment += delta*(alignment - alignmentMean);
            }
            if (state.time >= boundStartTime + boundObservationTime) {
                // Orbit-averaged radiated power measured over the observation
                // window closes the secular Coulomb inspiral to 0.01*a0.
                const double observed = state.time - boundStartTime;
                const double radiated =
                    state.radiatedEnergy + state.dipoleConstraintEnergy
                    - boundStartChargeRadiated;
                const double finalEnergy = coulombPairEnergy(state);
                if (!(finalEnergy < 0.0)) {
                    // Marginal capture that drifted back above threshold: the
                    // pair was never really bound, so let it run on and be
                    // classified by whatever it does next.
                    bound = false;
                    alignmentCount = 0;
                    alignmentMean = 0.0;
                    alignmentSecondMoment = 0.0;
                    continue;
                }
                if (observed > 0.0 && radiated > 0.0) {
                    const double power = radiated/observed;
                    const double semiMajorAxis =
                        -coulombStrength/(2.0*finalEnergy);
                    result.boundRadiatedPowerWatts = power;
                    // The secular model dE/dt = P(a) with P proportional to
                    // a^-4 takes the power ORBIT-AVERAGED, and it extrapolates
                    // from the orbit the pair is on at the END of the window,
                    // so that orbit's period is what the window has to cover.
                    // Capture is only detected after closest approach, so a
                    // shorter window samples the outbound periapsis leg.  For
                    // these near-parabolic captures (e ~ 0.999) the periapsis
                    // power exceeds the orbit average by orders of magnitude,
                    // so a sub-orbit window does not measure the quantity the
                    // formula consumes.  boundObservationTimeCap truncates the
                    // window to well under one period for loose captures, so
                    // the gate is not hypothetical.
                    const double finalPeriod = 2.0*pi*std::sqrt(
                        reducedMass*semiMajorAxis*semiMajorAxis*semiMajorAxis
                        / coulombStrength);
                    const double gatePeriod = std::max(
                        boundStartOrbitalPeriod,
                        std::isfinite(finalPeriod) ? finalPeriod : 0.0);
                    result.boundOrbitalPeriodSeconds = gatePeriod > 0.0
                        ? gatePeriod : std::numeric_limits<double>::quiet_NaN();
                    result.boundObservedOrbits = gatePeriod > 0.0
                        ? observed/gatePeriod : 0.0;
                    if (semiMajorAxis > collisionBoundaryRadius
                        && result.boundObservedOrbits >= 1.0) {
                        const double ratio =
                            collisionBoundaryRadius/semiMajorAxis;
                        const double collapse = (-finalEnergy)/(3.0*power)
                            *(1.0 - ratio*ratio*ratio);
                        if (collapse > 0.0 && std::isfinite(collapse)) {
                            result.collapseTimeSeconds = collapse;
                        }
                    }
                }
                // Parallel magnetic moments correspond to antiparallel spins,
                // i.e. the singlet.  The 0.5 threshold on an initially isotropic
                // relative orientation also reproduces the 1:3 para:ortho
                // statistical weight, which is a useful consistency check.
                const bool parallelMoments = alignmentMean >= 0.5;
                return finish(parallelMoments
                    ? InteractionOutcome::ParaPositronium
                    : InteractionOutcome::OrthoPositronium, state);
            }
            continue;
        }
        if (passedClosestApproach
            && currentRadius >= configuration.matchingRadius) {
            const double crossingFraction=separationCrossingFraction(
                beforeStep,state,configuration.matchingRadius);
            const State escapeState=interpolateState(
                beforeStep,state,crossingFraction);
            const Vec3 outgoing=escapeState.firstVelocity
                               -escapeState.secondVelocity;
            if (!(outgoing.norm() > 0.0) || !std::isfinite(outgoing.norm())) {
                return finish(InteractionOutcome::NumericalFailure,escapeState);
            }
            result.scatteringAngleDegrees = std::acos(std::clamp(
                dot(beamDirection, outgoing)/outgoing.norm(), -1.0, 1.0))
                * 180.0/pi;
            return finish(InteractionOutcome::Scattering,escapeState);
        }
    }
    return finish(InteractionOutcome::Unresolved, state);
}

std::vector<InteractionEvent> runInteractionExperiment(
    std::uint64_t masterSeed, const InteractionConfiguration& configuration,
    int runCount) {
    std::vector<InteractionEvent> events(static_cast<size_t>(runCount));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::mutex outputMutex;
    const int workerCount = std::min(runCount,
        static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
    std::cout << "Running " << runCount << " interaction trajectories on "
              << workerCount << " worker" << (workerCount == 1 ? "" : "s")
              << ".\n";
    const auto worker = [&]() {
        while (true) {
            const int index = nextIndex.fetch_add(1);
            if (index >= runCount) break;
            const std::uint64_t eventSeed = splitMix64(
                masterSeed + 0x5bf03635ULL
                + static_cast<std::uint64_t>(index));
            // Same accuracy ladder as the beam experiment.  A tighter setting
            // is unaffordable here: the final plunge toward the collision
            // boundary is stiff, and every extra subdivision level doubles the
            // work in exactly the region where the step is already smallest.
            InteractionEvent event = simulateInteractionEvent(
                eventSeed, configuration, {.relativeTolerance=1.0e-3,
                                           .maximumDepth=8,
                                           .reactionModel=gRadiationReactionModel});
            if (event.outcome == InteractionOutcome::NumericalFailure) {
                event = simulateInteractionEvent(
                    eventSeed, configuration, {.relativeTolerance=1.0e-5,
                                               .maximumDepth=12,
                                               .reactionModel=gRadiationReactionModel});
            }
            events[static_cast<size_t>(index)] = std::move(event);
            const int done = completed.fetch_add(1) + 1;
            if (done % 10 == 0 || done == runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout << "Interaction trajectories: " << done << "/"
                          << runCount << '\n';
            }
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    for (int index = 0; index < workerCount; ++index) workers.emplace_back(worker);
    for (std::thread& thread : workers) thread.join();
    return events;
}

std::vector<BeamEvent> runBeamExperiment(std::uint64_t masterSeed,
                                         const BeamConfiguration& configuration,
                                         int runCount) {
    std::vector<BeamEvent> events(static_cast<size_t>(runCount));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::atomic<int> retried{0};
    std::atomic<int> recovered{0};
    std::mutex outputMutex;
    const int workerCount = std::min(runCount,
        static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    const auto worker = [&]() {
        while (true) {
            const int index = nextIndex.fetch_add(1);
            if (index >= runCount) break;
            const std::uint64_t eventSeed=splitMix64(
                masterSeed+static_cast<std::uint64_t>(index));
            BeamEvent event=simulateBeamEvent(eventSeed,configuration,
                {.relativeTolerance=1.0e-3,.maximumDepth=8,
                 .reactionModel=gRadiationReactionModel});
            if(event.outcome==BeamOutcome::NumericalFailure) {
                retried.fetch_add(1);
                event=simulateBeamEvent(eventSeed,configuration,
                    {.relativeTolerance=1.0e-5,.maximumDepth=12,
                     .reactionModel=gRadiationReactionModel});
                if(event.outcome!=BeamOutcome::NumericalFailure)
                    recovered.fetch_add(1);
            }
            events[static_cast<size_t>(index)]=std::move(event);
            const int done = completed.fetch_add(1) + 1;
            if (done % 10 == 0 || done == runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout << "Beam trajectories: " << done << "/" << runCount << '\n';
            }
        }
    };
    for (int index = 0; index < workerCount; ++index) workers.emplace_back(worker);
    for (std::thread& workerThread : workers) workerThread.join();
    if(retried.load()>0) {
        const char* retryNoun=retried.load()==1?" trajectory":" trajectories";
        const char* recoveredNoun=recovered.load()==1?" trajectory":" trajectories";
        std::cout<<"Adaptive recovery: retried "<<retried.load()<<retryNoun
                 <<", recovered "<<recovered.load()<<recoveredNoun<<".\n";
    }
    return events;
}

double binomialCrossSectionError(double sampledArea, int count, int total) {
    const double probability = static_cast<double>(count) / total;
    return sampledArea * std::sqrt(probability * (1.0 - probability) / total);
}

double sampleQuantile(std::vector<double> values, double probability) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double position = std::clamp(probability, 0.0, 1.0)
                          * static_cast<double>(values.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(position));
    const size_t upper = static_cast<size_t>(std::ceil(position));
    const double fraction = position - static_cast<double>(lower);
    return values[lower] + fraction * (values[upper] - values[lower]);
}

std::string crossSectionEstimate(double sampledArea, int count, int total,
                                 double unit) {
    std::ostringstream text;
    if (count == 0) {
        const double upper95 = sampledArea
            * (1.0 - std::pow(0.05, 1.0/static_cast<double>(total)));
        text << "0; < " << upper95/unit << " (95% CL)";
    } else if (count == total) {
        const double lower95 = sampledArea
            * std::pow(0.05, 1.0/static_cast<double>(total));
        text << sampledArea/unit << "; > " << lower95/unit << " (95% CL)";
    } else {
        text << sampledArea * count/(total*unit) << " +/- "
             << binomialCrossSectionError(sampledArea, count, total)/unit;
    }
    return text.str();
}

int showBeamStatistics(std::uint64_t seed, int selectedPhenomenon, int runCount,
                       double beamEnergyEv, double thetaMinimumDegrees,
                       int angleBins, double impactMaximumPm,
                       double matchingRadiusPm) {
    const BeamConfiguration configuration = makeBeamConfiguration(
        selectedPhenomenon, beamEnergyEv, thetaMinimumDegrees, angleBins,
        impactMaximumPm, matchingRadiusPm);
    const int workerCount = std::min(runCount,
        static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
    const double matchingPotentialRatio = pairCoulombStrength
        / (configuration.matchingRadius * configuration.centreOfMassKineticEnergy);
    const bool matchingRegionWarning = configuration.matchingRadius
            < 20.0*configuration.impactParameterMaximum
        || matchingPotentialRatio > 0.01;
    std::cout << "Running a " << firstSpecies.name << '+' << secondSpecies.name
              << " beam experiment with " << runCount
              << (runCount == 1 ? " trajectory on " : " trajectories on ")
              << workerCount << " worker"
              << (workerCount == 1 ? "" : "s") << ".\n"
              << "K_CM = " << beamEnergyEv << " eV, b_max = "
              << configuration.impactParameterMaximum * 1.0e12
              << " pm, R_match = " << configuration.matchingRadius * 1.0e12
              << " pm, theta acceptance >= "
              << configuration.analysisThetaMinimum * 180.0/pi
              << " deg, max beta_infinity = "
              << configuration.maximumAsymptoticBeta << ".\n";
    if (configuration.maximumAsymptoticBeta > 0.1) {
        std::cout << "Warning: beta > 0.1; the Darwin/instantaneous interaction model is "
                     "outside its intended low-velocity regime.\n";
    }
    if (matchingRegionWarning) {
        std::cout << "Warning: the matching sphere is close to the interaction region; "
                     "increase --matching-radius-pm and check convergence.\n";
    }
    const std::vector<BeamEvent> events = runBeamExperiment(seed, configuration, runCount);

    int escaped = 0;
    int fiducial = 0;
    int shortRange = 0;
    int captured = 0;
    int unresolved = 0;
    int failed = 0;
    std::vector<double> fiducialEnergyLossesEv;
    for (const BeamEvent& event : events) {
        switch (event.outcome) {
            case BeamOutcome::Escaped:
                ++escaped;
                if (event.scatteringAngle >= configuration.analysisThetaMinimum) {
                    ++fiducial;
                    fiducialEnergyLossesEv.push_back(
                        (configuration.centreOfMassKineticEnergy
                         - event.outgoingEnergy) / eCharge);
                }
                break;
            case BeamOutcome::ShortRange: ++shortRange; break;
            case BeamOutcome::Captured: ++captured; break;
            case BeamOutcome::Unresolved: ++unresolved; break;
            case BeamOutcome::NumericalFailure: ++failed; break;
        }
    }
    if (failed > 0) {
        std::cerr << "Outcomes: escaped=" << escaped << ", short-range=" << shortRange
                  << ", captured=" << captured
                  << ", unresolved=" << unresolved << ", failed=" << failed << '\n'
                  << "Cross-section report is INVALID because at least one trajectory "
                     "became numerically non-finite. No plots were produced.\n";
        bool persistenceOk = reportArchiveOperation(
            statistics_archive::writeScientificReferencesText(),
            "scientific-reference catalogue");
        return persistenceOk ? 2 : 3;
    }
    constexpr double diagnosticFloor = 64.0 * std::numeric_limits<double>::epsilon();
    const double energyFloor = diagnosticFloor * 2.0*firstMass*c*c;
    const double momentumFloor = diagnosticFloor * 2.0*firstMass*c;
    std::vector<double> relativeEnergyClosures;
    std::vector<double> absoluteEnergyClosures;
    std::vector<double> relativeMomentumClosures;
    std::vector<double> relativeAngularMomentumClosures;
    std::vector<double> logMomentumClosures;
    std::vector<double> logAngularMomentumClosures;
    // Independent physical residuals.  The three "closure" vectors above are
    // exact identities of the discrete bookkeeping (boundField* is defined as
    // the residual that closes them), so they only measure roundoff and the
    // endpoint interpolation.  These two measure the model itself.
    std::vector<double> reactionMismatchFractions;
    std::vector<double> boundReservoirFractions;
    // Same reservoir, normalized by the ORBIT's own mechanical energy instead
    // of by the radiated energy.  The radiated-energy form is unusable as a
    // validity signal because its denominator can be arbitrarily small; this
    // one measures the reservoir against a scale that does not vanish.
    std::vector<double> boundReservoirVsOrbit;
    double maximumSchottFraction = 0.0;
    for (const BeamEvent& event : events) {
        if (event.outcome != BeamOutcome::Escaped || !event.diagnosticsValid) continue;
        const EndpointDiagnostics& initial = event.initialDiagnostics;
        const EndpointDiagnostics& final = event.finalDiagnostics;
        const double initialBalance=initial.mechanicalEnergy
            +initial.radiatedEnergy+initial.boundFieldEnergy;
        const double finalBalance=final.mechanicalEnergy
            +final.radiatedEnergy+final.boundFieldEnergy;
        const double energyScale = std::max({std::abs(initial.mechanicalEnergy),
                                              std::abs(initialBalance), energyFloor});
        const double energyClosure = (finalBalance - initialBalance)/energyScale;
        const double momentumScale = std::max(
            initial.canonicalMomentumScale, momentumFloor);
        const double momentumClosure=((final.noetherMomentum
            +final.radiatedMomentum+final.boundFieldMomentum)
            -(initial.noetherMomentum+initial.radiatedMomentum
              +initial.boundFieldMomentum)).norm()/momentumScale;
        const double angularMomentumScale = std::max({
            initial.noetherAngularMomentum.norm(),
            0.5*event.impactParameter*initial.canonicalMomentumScale, hbar});
        const double angularMomentumClosure =
            ((final.noetherAngularMomentum+final.radiatedAngularMomentum
                +final.boundFieldAngularMomentum)
             - (initial.noetherAngularMomentum
                +initial.radiatedAngularMomentum
                +initial.boundFieldAngularMomentum)).norm()/angularMomentumScale;
        if (!std::isfinite(energyClosure) || !std::isfinite(momentumClosure)
            || !std::isfinite(angularMomentumClosure)) {
            continue;
        }
        relativeEnergyClosures.push_back(energyClosure);
        absoluteEnergyClosures.push_back(std::abs(energyClosure));
        relativeMomentumClosures.push_back(momentumClosure);
        relativeAngularMomentumClosures.push_back(angularMomentumClosure);
        logMomentumClosures.push_back(std::log10(std::max(momentumClosure,
                                                          diagnosticFloor)));
        logAngularMomentumClosures.push_back(std::log10(std::max(
            angularMomentumClosure, diagnosticFloor)));
        maximumSchottFraction = std::max({maximumSchottFraction,
            std::abs(initial.schottEnergy)/energyScale,
            std::abs(final.schottEnergy)/energyScale});
        // Scale both physical residuals by the energy the run claims to have
        // radiated over the same interval: the question they answer is "how
        // large is the unmodelled piece compared to the reported emission".
        const double radiatedChange = final.radiatedEnergy
                                    - initial.radiatedEnergy;
        const double radiationScale = std::max(std::abs(radiatedChange),
                                                energyFloor);
        const double reactionFraction = std::abs(final.reactionEnergyMismatch
            - initial.reactionEnergyMismatch)/radiationScale;
        const double reservoirFraction = std::abs(final.boundFieldEnergy
            - initial.boundFieldEnergy)/radiationScale;
        if (std::isfinite(reactionFraction)) {
            reactionMismatchFractions.push_back(reactionFraction);
        }
        if (std::isfinite(reservoirFraction)) {
            boundReservoirFractions.push_back(reservoirFraction);
        }
        const double orbitScale = std::max(std::abs(initial.mechanicalEnergy),
                                           1.0e-300);
        const double reservoirVsOrbit = std::abs(final.boundFieldEnergy
            - initial.boundFieldEnergy)/orbitScale;
        if (std::isfinite(reservoirVsOrbit)) {
            boundReservoirVsOrbit.push_back(reservoirVsOrbit);
        }
    }
    const double sampledArea = pi * configuration.impactParameterMaximum
                             * configuration.impactParameterMaximum;
    constexpr double barn = 1.0e-28;
    const double analyticFiducial = pi * configuration.coulombLength
        * configuration.coulombLength
        / std::pow(std::tan(0.5 * configuration.analysisThetaMinimum), 2);
    const double analyticCutoff = pi * configuration.analyticCutoffImpactParameter
        * configuration.analyticCutoffImpactParameter;
    const double fiducialProbability = static_cast<double>(fiducial)/runCount;
    const double rutherfordNormalization = analyticFiducial > 0.0
        ? sampledArea*fiducialProbability/analyticFiducial : 0.0;
    const double rutherfordNormalizationError = analyticFiducial > 0.0
        ? sampledArea*std::sqrt(fiducialProbability
            * (1.0 - fiducialProbability)/runCount)/analyticFiducial : 0.0;
    constexpr double normal95 = 1.95996398454005;
    const double wilsonDenominator = 1.0 + normal95*normal95/runCount;
    const double wilsonCenter = (fiducialProbability
        + normal95*normal95/(2.0*runCount))/wilsonDenominator;
    const double wilsonHalfWidth = normal95/wilsonDenominator
        * std::sqrt(fiducialProbability*(1.0-fiducialProbability)/runCount
            + normal95*normal95/(4.0*runCount*runCount));
    const double normalizationConversion = analyticFiducial > 0.0
        ? sampledArea/analyticFiducial : 0.0;
    const double rutherfordNormalizationLower95 = normalizationConversion
        * std::max(0.0, wilsonCenter - wilsonHalfWidth);
    const double rutherfordNormalizationUpper95 = normalizationConversion
        * std::min(1.0, wilsonCenter + wilsonHalfWidth);
    std::cout << "Outcomes: escaped=" << escaped << ", short-range=" << shortRange
              << ", captured=" << captured
              << ", unresolved=" << unresolved << ", failed=" << failed << '\n'
              << "Model sigma(theta >= acceptance) = "
              << crossSectionEstimate(sampledArea, fiducial, runCount, barn)
              << " barn; Rutherford reference = "
              << analyticFiducial/barn << " barn\n"
              << "Rutherford-shape normalization fit: C_R="
              << rutherfordNormalization;
    if (fiducial > 0 && fiducial < runCount) {
        std::cout << " +/- " << rutherfordNormalizationError;
    }
    std::cout << "; Wilson 95% CI [" << rutherfordNormalizationLower95
              << ", " << rutherfordNormalizationUpper95 << "]\n";
    if (configuration.shortRangeFocus) {
        std::cout << "Model sigma(reach cutoff) = "
                  << crossSectionEstimate(sampledArea, shortRange, runCount, barn)
                  << " barn; pure-Coulomb cutoff reference = "
                  << analyticCutoff/barn << " barn\n"
                  << "Important: cutoff reach is not a QED annihilation cross section.\n";
    }
    if(captured>0) {
        std::cout<<"Captured="<<captured
                 <<": negative-energy outgoing states are a separate bound channel "
                   "and are excluded from escaped/cutoff cross sections.\n";
    }
    std::cout << "Elastic data are a classical low-velocity model, not Bhabha scattering.\n";
    if (unresolved > 0) {
        std::cout << "Warning: " << unresolved
                  << " trajectories did not produce a resolved outgoing endpoint within "
                     "the flight-time/accuracy gate; reported channel cross sections "
                     "are censored finite-gate values.\n";
    }
    if (runCount < 1000) {
        std::cout << "Statistical note: " << runCount
                  << (runCount == 1 ? " trajectory is" : " trajectories are")
                  << " suitable for a preview; use at least 1000 "
                     "for a smoother differential cross section.\n";
    }
    if (!relativeEnergyClosures.empty()) {
        std::cout << "Bookkeeping identity residual on "
                  << relativeEnergyClosures.size()
                  << " escaped trajectories at R_match (roundoff/interpolation "
                     "only, NOT a conservation test): median |deltaE|="
                  << sampleQuantile(absoluteEnergyClosures, 0.5)
                  << ", median delta(P_N+P_rad)="
                  << sampleQuantile(relativeMomentumClosures, 0.5)
                  << ", median delta(J_N+J_rad)="
                  << sampleQuantile(relativeAngularMomentumClosures, 0.5)
                  << ".\n"
                  << "Independent physical residuals (relative to the reported "
                     "radiated energy): median |dE_reaction-vs-flux|/E_rad="
                  << sampleQuantile(reactionMismatchFractions, 0.5)
                  << ", median |E_bound|/E_rad="
                  << sampleQuantile(boundReservoirFractions, 0.5) << ".\n";
    }
    // Self-diagnosis of the regime, driven by measurement rather than by a
    // hard-coded experiment number.  The reservoir is the energy the balance
    // could not attribute to mechanics or radiation and had to invent; once it
    // rivals the orbit's own energy, the bookkeeping has stopped describing
    // the trajectory and no energy quantity from the run may be read as a
    // measurement.
    //
    // This is what the physical-completeness audit found for the short-range
    // channel, where the reservoir reaches 1.4 to 9.4 times the orbit energy
    // while the bound and wide-scattering channels stay near 1e-4.  Stating it
    // as a threshold on the measured value rather than as a ban on experiment
    // 3 means it also fires for any pair or any settings that stray into the
    // same regime, and stays silent when they do not.
    if (!boundReservoirVsOrbit.empty()) {
        // The median is the wrong statistic here and reporting it alone would
        // hide the problem: a beam ensemble is dominated by wide fly-bys whose
        // reservoir is negligible, while the close encounters this channel
        // exists to study sit in the tail.  For experiment 3 the median is
        // 0.057 against 6.7e-08 for experiment 4 -- six orders apart already,
        // but far below the level the audit found on close trajectories.  The
        // tail and the exceedance fraction are what decide.
        const double medianReservoir =
            sampleQuantile(boundReservoirVsOrbit, 0.5);
        const double tailReservoir =
            sampleQuantile(boundReservoirVsOrbit, 0.95);
        const std::size_t overOrbit = static_cast<std::size_t>(
            std::count_if(boundReservoirVsOrbit.begin(),
                          boundReservoirVsOrbit.end(),
                          [](double value){ return value > 1.0; }));
        const double overFraction = static_cast<double>(overOrbit)
            /static_cast<double>(boundReservoirVsOrbit.size());
        std::cout << "Reservoir against orbit energy |E_bound|/|E_rel|: median="
                  << medianReservoir << ", 95th=" << tailReservoir
                  << ", exceeding the orbit energy in " << overOrbit << " of "
                  << boundReservoirVsOrbit.size() << " trajectories.\n";
        // The bound is on the tail, not on the median and not on exceeding
        // the orbit energy outright.  Measured: experiment 3 gives a 95th
        // percentile of 0.35 against 1.8e-07 for experiment 4, so 0.1 sits six
        // orders clear of both, and it means one trajectory in twenty has a
        // tenth of its orbit energy unattributed.
        //
        // The 1.4-9.4 figures in the completeness audit come from a single
        // diagnostic trajectory, not from a beam ensemble; the ensemble is
        // dominated by wide fly-bys and comes out milder.
        if (tailReservoir > 0.1 || overFraction > 0.02) {
            std::cout << "WARNING: on part of this ensemble the reconstructed "
                         "bound-field reservoir reaches a substantial "
                         "fraction of the orbit's own energy, so those "
                         "trajectories are outside the regime where the "
                         "model's energy bookkeeping describes them.\n"
                      << "         Angular distributions rest on scattering "
                         "geometry and are not invalidated by this, but NO "
                         "energy quantity from this run may be read as a "
                         "measurement.\n";
        }
    }
    gROOT->SetBatch(kTRUE);
    root_export::preparePdfExporter();
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    const std::string canvasTitle = configuration.shortRangeFocus
        ? "e+e- short-range beam channel (model boundary)"
        : "e+e- elastic beam scattering (classical model)";
    TCanvas canvas("beam_statistics", canvasTitle.c_str(), 1280, 900);
    canvas.SetFillColor(kWhite);
    TPad distributionsPage("beam_distributions_page", "Beam distributions",
                           0.0, 0.0, 1.0, 1.0);
    TPad diagnosticsPage("beam_diagnostics_page", "Numerical closure diagnostics",
                         0.0, 0.0, 1.0, 1.0);
    for (TPad* page : {&distributionsPage, &diagnosticsPage}) {
        page->SetFillColor(kWhite);
        page->SetFillStyle(1001);
    }
    canvas.cd();
    distributionsPage.Draw();
    diagnosticsPage.Draw();
    distributionsPage.cd();
    distributionsPage.Divide(2, 2, 0.006, 0.006);
    diagnosticsPage.cd();
    diagnosticsPage.Divide(2, 2, 0.006, 0.006);
    std::vector<std::unique_ptr<TPaveText>> beamAnalysisBoxes;
    std::vector<std::unique_ptr<TF1>> beamAnalysisFunctions;

    const double thetaLowDegrees = configuration.analysisThetaMinimum * 180.0/pi;
    std::vector<double> angularEdges(static_cast<size_t>(configuration.angleBins + 1));
    const double acceptedReferenceImpact = configuration.coulombLength
        / std::tan(0.5 * configuration.analysisThetaMinimum);
    for (int edge = 0; edge < configuration.angleBins; ++edge) {
        const double remainingArea = 1.0
            - static_cast<double>(edge) / configuration.angleBins;
        const double referenceImpact = acceptedReferenceImpact
            * std::sqrt(remainingArea);
        angularEdges[static_cast<size_t>(edge)] = 2.0
            * std::atan(configuration.coulombLength / referenceImpact) * 180.0/pi;
    }
    angularEdges.back() = 180.0;
    TH1D angularCrossSection("angular_cross_section",
        "Differential elastic cross section;#theta_{CM} [deg];d#sigma/d#Omega [barn/sr]",
        configuration.angleBins, angularEdges.data());
    angularCrossSection.SetDirectory(nullptr);
    angularCrossSection.SetLineColor(plot_style::crem());
    angularCrossSection.SetMarkerColor(plot_style::crem());
    angularCrossSection.SetMarkerStyle(20);
    for (const BeamEvent& event : events) {
        if (event.outcome == BeamOutcome::Escaped
            && event.scatteringAngle >= configuration.analysisThetaMinimum) {
            angularCrossSection.Fill(event.scatteringAngle * 180.0/pi);
        }
    }
    for (int bin = 1; bin <= angularCrossSection.GetNbinsX(); ++bin) {
        const double thetaLow = angularCrossSection.GetXaxis()->GetBinLowEdge(bin) * pi/180.0;
        const double thetaHigh = angularCrossSection.GetXaxis()->GetBinUpEdge(bin) * pi/180.0;
        const double solidAngle = 2.0*pi*(std::cos(thetaLow) - std::cos(thetaHigh));
        const double count = angularCrossSection.GetBinContent(bin);
        const double scale = sampledArea / (runCount * solidAngle * barn);
        angularCrossSection.SetBinContent(bin, count * scale);
        angularCrossSection.SetBinError(bin,
            std::sqrt(count * (1.0 - count/runCount)) * scale);
    }

    TH1D rutherfordHistogram("rutherford_bin_average",
        "Rutherford bin average", configuration.angleBins, angularEdges.data());
    rutherfordHistogram.SetDirectory(nullptr);
    rutherfordHistogram.SetLineColor(plot_style::theory());
    rutherfordHistogram.SetLineWidth(3);
    rutherfordHistogram.SetLineStyle(2);
    rutherfordHistogram.SetFillStyle(0);
    double referenceMinimum = std::numeric_limits<double>::infinity();
    for (int bin = 1; bin <= angularCrossSection.GetNbinsX(); ++bin) {
        const double thetaLow = angularCrossSection.GetXaxis()->GetBinLowEdge(bin) * pi/180.0;
        const double thetaHigh = angularCrossSection.GetXaxis()->GetBinUpEdge(bin) * pi/180.0;
        const double solidAngle = 2.0*pi*(std::cos(thetaLow) - std::cos(thetaHigh));
        const double crossSectionInBin = pi * configuration.coulombLength
            * configuration.coulombLength
            * (1.0/std::pow(std::tan(0.5*thetaLow), 2)
             - 1.0/std::pow(std::tan(0.5*thetaHigh), 2));
        const double binAverage = crossSectionInBin/(solidAngle*barn);
        rutherfordHistogram.SetBinContent(bin, binAverage);
        referenceMinimum = std::min(referenceMinimum, binAverage);
    }
    TH1D fittedRutherfordHistogram(rutherfordHistogram);
    fittedRutherfordHistogram.SetName("rutherford_normalization_fit");
    fittedRutherfordHistogram.SetTitle("Fitted Rutherford-shape normalization");
    fittedRutherfordHistogram.Scale(rutherfordNormalization);
    fittedRutherfordHistogram.SetLineColor(plot_style::crem());
    fittedRutherfordHistogram.SetLineStyle(7);
    fittedRutherfordHistogram.SetLineWidth(3);

    constexpr size_t cumulativePointCount = 100;
    std::array<double, cumulativePointCount> cumulativeAngles{};
    std::array<double, cumulativePointCount> cumulativeModel{};
    std::array<double, cumulativePointCount> cumulativeErrors{};
    std::array<double, cumulativePointCount> cumulativeReference{};
    std::array<double, cumulativePointCount> cumulativeFitted{};
    std::array<double, cumulativePointCount> zeroAngleErrors{};
    const double cumulativeUpper = configuration.analysisThetaMinimum
        + 0.98 * (pi - configuration.analysisThetaMinimum);
    for (size_t point = 0; point < cumulativePointCount; ++point) {
        const double fraction = static_cast<double>(point)
                              / (cumulativePointCount - 1);
        const double theta = configuration.analysisThetaMinimum
            + fraction * (cumulativeUpper - configuration.analysisThetaMinimum);
        int accepted = 0;
        for (const BeamEvent& event : events) {
            if (event.outcome == BeamOutcome::Escaped
                && event.scatteringAngle >= theta) {
                ++accepted;
            }
        }
        cumulativeAngles[point] = theta * 180.0/pi;
        cumulativeModel[point] = sampledArea * accepted / (runCount*barn);
        cumulativeErrors[point] = binomialCrossSectionError(
            sampledArea, accepted, runCount) / barn;
        cumulativeReference[point] = pi * configuration.coulombLength
            * configuration.coulombLength
            / std::pow(std::tan(0.5*theta), 2) / barn;
        cumulativeFitted[point] = rutherfordNormalization
            * cumulativeReference[point];
    }
    TGraphErrors cumulativeGraph(static_cast<int>(cumulativePointCount),
        cumulativeAngles.data(), cumulativeModel.data(), zeroAngleErrors.data(),
        cumulativeErrors.data());
    cumulativeGraph.SetTitle(
        "Cumulative elastic cross section (correlated thresholds);minimum #theta_{CM} [deg];#sigma(#theta #geq #theta_{min}) [barn]");
    cumulativeGraph.SetLineColor(plot_style::crem());
    cumulativeGraph.SetMarkerColor(plot_style::crem());
    cumulativeGraph.SetMarkerStyle(20);
    cumulativeGraph.SetMarkerSize(0.55);
    TGraph cumulativeRutherford(static_cast<int>(cumulativePointCount),
        cumulativeAngles.data(), cumulativeReference.data());
    cumulativeRutherford.SetLineColor(plot_style::theory());
    cumulativeRutherford.SetLineWidth(3);
    cumulativeRutherford.SetLineStyle(2);
    TGraph cumulativeFit(static_cast<int>(cumulativePointCount),
        cumulativeAngles.data(), cumulativeFitted.data());
    cumulativeFit.SetLineColor(plot_style::crem());
    cumulativeFit.SetLineStyle(7);
    cumulativeFit.SetLineWidth(3);

    std::unique_ptr<TH1D> energyLossSpectrum;
    if (configuration.shortRangeFocus && !fiducialEnergyLossesEv.empty()) {
        const double incidentEnergyEv = configuration.centreOfMassKineticEnergy/eCharge;
        const auto [observedMinimum, observedMaximum] = std::minmax_element(
            fiducialEnergyLossesEv.begin(), fiducialEnergyLossesEv.end());
        double energyLower = std::min(-0.02*incidentEnergyEv, *observedMinimum);
        double energyUpper = std::max(1.02*incidentEnergyEv, *observedMaximum);
        const double energyPadding = std::max(1.0e-12,
            0.03*(energyUpper - energyLower));
        energyLower -= energyPadding;
        energyUpper += energyPadding;
        energyLossSpectrum = std::make_unique<TH1D>("energy_loss_cross_section",
            "Fiducial energy-loss cross section;#DeltaE=K_{CM}-E_{out} [eV];d#sigma(#theta#geq#theta_{min})/d#DeltaE [barn/eV]",
            histogramBins(fiducialEnergyLossesEv.size()),
            energyLower, energyUpper);
        energyLossSpectrum->SetDirectory(nullptr);
        energyLossSpectrum->SetLineColor(plot_style::crem());
        energyLossSpectrum->SetFillColorAlpha(plot_style::crem(), 0.35);
        energyLossSpectrum->Sumw2();
        const double energyWeight = sampledArea
            / (runCount * barn * energyLossSpectrum->GetBinWidth(1));
        for (double energyLoss : fiducialEnergyLossesEv) {
            energyLossSpectrum->Fill(energyLoss, energyWeight);
        }
    }

    distributionsPage.cd(1);
    gPad->SetGrid();
    gPad->SetLogy();
    angularCrossSection.SetMinimum(std::max(1.0e-12,
        referenceMinimum * 0.2));
    angularCrossSection.SetMaximum(std::max({rutherfordHistogram.GetMaximum(),
        fittedRutherfordHistogram.GetMaximum(), angularCrossSection.GetMaximum()}) * 1.5);
    angularCrossSection.Draw("E1");
    rutherfordHistogram.Draw("HIST SAME");
    fittedRutherfordHistogram.Draw("HIST SAME");
    TLegend angularLegend(0.55, 0.69, 0.89, 0.89);
    angularLegend.SetFillColorAlpha(kWhite, 0.88);
    angularLegend.AddEntry(&angularCrossSection, "trajectory model", "lep");
    angularLegend.AddEntry(&rutherfordHistogram, "Rutherford (bin average)", "l");
    angularLegend.AddEntry(&fittedRutherfordHistogram,
                           "C_{R} #times Rutherford fit", "l");
    angularLegend.Draw();
    std::string normalizationEstimate = "C_{R} = "
        + compactNumber(rutherfordNormalization);
    if (fiducial > 0 && fiducial < runCount) {
        normalizationEstimate += " #pm "
            + compactNumber(rutherfordNormalizationError);
    }
    drawAnalysisBox(beamAnalysisBoxes, 0.12, 0.10, 0.57, 0.34, {
        plot_style::key(true, false, false, true),
        "MC trajectories: N = " + std::to_string(runCount),
        "Fit: C_{R} d#sigma_{R}/d#Omega (binomial MLE)",
        normalizationEstimate,
        "95% Wilson: [" + compactNumber(rutherfordNormalizationLower95)
            + ", " + compactNumber(rutherfordNormalizationUpper95) + "]",
        "Reference: pure Coulomb C_{R}=1",
        configuration.shortRangeFocus
            ? "C_{R} applies only to escaped events; cutoff excluded"
            : "C_{R} applies to the escaped elastic channel",
        "Experimental:",
        "no matching free e^{+}e^{-} dataset loaded"
    }, 0.020);

    distributionsPage.cd(2);
    gPad->SetGrid();
    gPad->SetLogy();
    cumulativeGraph.SetMinimum(std::max(1.0e-12, cumulativeReference.back() * 0.2));
    cumulativeGraph.SetMaximum(std::max(cumulativeReference.front(),
        *std::max_element(cumulativeModel.begin(), cumulativeModel.end())) * 1.5);
    cumulativeGraph.Draw("ALP");
    cumulativeRutherford.Draw("L SAME");
    cumulativeFit.Draw("L SAME");
    TLegend cumulativeLegend(0.53, 0.68, 0.89, 0.89);
    cumulativeLegend.SetFillColorAlpha(kWhite, 0.88);
    cumulativeLegend.AddEntry(&cumulativeGraph, "trajectory model", "lep");
    cumulativeLegend.AddEntry(&cumulativeRutherford, "Rutherford", "l");
    cumulativeLegend.AddEntry(&cumulativeFit,
                              "projection of panel-1 fit", "l");
    cumulativeLegend.Draw();
    drawAnalysisBox(beamAnalysisBoxes, 0.12, 0.10, 0.58, 0.32, {
        plot_style::key(true, false, false, true),
        "MC trajectories: N = " + std::to_string(runCount),
        "Model: C_{R} #pi l_{C}^{2} cot^{2}(#theta/2)",
        "C_{R} from differential-panel fit = "
            + compactNumber(rutherfordNormalization),
        "No independent fit: thresholds are correlated",
        "Experimental:",
        "no matching cuts/covariance dataset loaded"
    }, 0.022);
    distributionsPage.cd(3);
    std::unique_ptr<TPaveText> thirdPanelInformation;
    if (energyLossSpectrum) {
        gPad->SetGrid();
        energyLossSpectrum->Draw("HIST E1");
        const GaussianFitSummary energyLossSummary =
            gaussianMaximumLikelihood(fiducialEnergyLossesEv);
        drawAnalysisBox(beamAnalysisBoxes, 0.45, 0.55, 0.94, 0.90, {
            "Accepted trajectories: N = "
                + std::to_string(fiducialEnergyLossesEv.size()),
            "Fit: none - escaped channel; cutoff-censored",
            "No justified universal Gaussian/Landau law",
            "mean = " + compactNumber(energyLossSummary.mean) + " eV",
            "median / p95 = "
                + compactNumber(sampleQuantile(fiducialEnergyLossesEv, 0.5))
                + " / " + compactNumber(sampleQuantile(fiducialEnergyLossesEv, 0.95))
                + " eV",
            "Experimental:",
            "no matching detector-response dataset"
        }, 0.022);
    } else if (!configuration.shortRangeFocus) {
        thirdPanelInformation = std::make_unique<TPaveText>(
            0.10, 0.18, 0.90, 0.82, "NDC");
        thirdPanelInformation->SetFillColorAlpha(kWhite, 0.94);
        thirdPanelInformation->SetTextAlign(12);
        thirdPanelInformation->SetTextFont(42);
        thirdPanelInformation->AddText("Fixed-energy two-body elastic channel");
        thirdPanelInformation->AddText("Ideal conservative result: E_{out}=K_{CM}");
        thirdPanelInformation->AddText("No detector response is modeled");
        thirdPanelInformation->AddText(
            "A narrow E_{out} histogram would display numerical drift");
        thirdPanelInformation->AddText("and is intentionally not drawn");
        thirdPanelInformation->AddText("Experimental: no matching response dataset");
        thirdPanelInformation->Draw();
    } else {
        thirdPanelInformation = std::make_unique<TPaveText>(
            0.12, 0.25, 0.88, 0.75, "NDC");
        thirdPanelInformation->SetFillColorAlpha(kWhite, 0.94);
        thirdPanelInformation->SetTextAlign(22);
        thirdPanelInformation->AddText(
            "No escaped event passed the angular acceptance");
        thirdPanelInformation->AddText(
            "No artificial energy-loss spectrum is drawn");
        thirdPanelInformation->AddText("Experimental: no matching dataset loaded");
        thirdPanelInformation->Draw();
    }
    distributionsPage.cd(4);
    TPaveText summary(0.07, 0.10, 0.93, 0.90, "NDC");
    summary.SetFillColorAlpha(kWhite, 0.94);
    summary.SetTextAlign(12);
    summary.SetTextFont(42);
    std::ostringstream setupLine;
    setupLine << "K_{CM} = " << beamEnergyEv << " eV, N = " << runCount
              << ", max #beta_{#infty} = "
              << configuration.maximumAsymptoticBeta;
    summary.AddText(setupLine.str().c_str());
    std::ostringstream acceptanceLine;
    acceptanceLine << "#theta_{CM} >= " << thetaLowDegrees << " deg";
    summary.AddText(acceptanceLine.str().c_str());
    std::ostringstream elasticLine;
    elasticLine << "#sigma_{fid} = "
                << crossSectionEstimate(sampledArea, fiducial, runCount, barn)
                << " barn";
    summary.AddText(elasticLine.str().c_str());
    std::ostringstream referenceLine;
    if (configuration.shortRangeFocus) {
        std::ostringstream cutoffLine;
        cutoffLine << "#sigma_{cutoff} = "
                   << crossSectionEstimate(sampledArea, shortRange, runCount, barn)
                   << " barn";
        summary.AddText(cutoffLine.str().c_str());
        referenceLine << "pure-Coulomb #sigma_{cutoff} = "
                      << analyticCutoff/barn << " barn";
    } else {
        referenceLine << "Rutherford #sigma_{fid} = "
                      << analyticFiducial/barn << " barn";
    }
    summary.AddText(referenceLine.str().c_str());
    std::ostringstream outcomeLine;
    outcomeLine << "escaped/cutoff/captured/gated/failed = " << escaped << '/'
                << shortRange << '/' << captured << '/' << unresolved << '/' << failed;
    summary.AddText(outcomeLine.str().c_str());
    summary.AddText("Classical low-v trajectory model");
    if (configuration.shortRangeFocus) {
        summary.AddText("#sigma_{cutoff} is not an annihilation cross section");
    }
    if (unresolved > 0) {
        summary.AddText("Finite flight-time/accuracy gate: channels are censored");
    }
    if (configuration.maximumAsymptoticBeta > 0.1) {
        summary.AddText("WARNING: outside the intended low-velocity regime");
    }
    if (matchingRegionWarning) {
        summary.AddText("WARNING: matching sphere may be too close");
    }
    summary.AddText("Experimental comparison:");
    summary.AddText("No matching free-pair dataset loaded");
    summary.AddText("Rutherford curves are analytic references, not data");
    summary.AddText("Relativistic experiment requires QED/Bhabha scattering");
    summary.Draw();

    double energyClosureLimit = 1.0e-15;
    for (double value : relativeEnergyClosures) {
        energyClosureLimit = std::max(energyClosureLimit, 1.10*std::abs(value));
    }
    double momentumLogUpper = 1.0;
    for (double value : logMomentumClosures) {
        momentumLogUpper = std::max(momentumLogUpper, value + 0.5);
    }
    double angularMomentumLogUpper = 1.0;
    for (double value : logAngularMomentumClosures) {
        angularMomentumLogUpper = std::max(angularMomentumLogUpper, value + 0.5);
    }
    TH1D energyBalanceHistogram("beam_energy_balance_closure",
        "Energy bookkeeping IDENTITY residual (not a conservation test)"
        ";#Delta(E_{N}+E_{rad}+E_{S})/E_{scale};Escaped trajectories",
        histogramBins(relativeEnergyClosures.size()),
        -energyClosureLimit, energyClosureLimit);
    TH1D momentumBalanceHistogram("beam_momentum_balance_closure",
        "Momentum bookkeeping IDENTITY residual (not a conservation test)"
        ";log_{10}(|#Delta(P_{N}+P_{rad})|/P_{scale});Escaped trajectories",
        histogramBins(logMomentumClosures.size()), -14.5, momentumLogUpper);
    TH1D angularMomentumBalanceHistogram("beam_angular_momentum_balance_closure",
        "Angular momentum bookkeeping IDENTITY residual (not a conservation test)"
        ";log_{10}(|#Delta(J_{N}+J_{rad})|/J_{scale});Escaped trajectories",
        histogramBins(logAngularMomentumClosures.size()), -14.5,
        angularMomentumLogUpper);
    styleHistogram(energyBalanceHistogram, plot_style::crem());
    styleHistogram(momentumBalanceHistogram, plot_style::crem());
    styleHistogram(angularMomentumBalanceHistogram, plot_style::crem());
    energyBalanceHistogram.SetStats(false);
    momentumBalanceHistogram.SetStats(false);
    angularMomentumBalanceHistogram.SetStats(false);
    for (double value : relativeEnergyClosures) energyBalanceHistogram.Fill(value);
    for (double value : logMomentumClosures) momentumBalanceHistogram.Fill(value);
    for (double value : logAngularMomentumClosures) {
        angularMomentumBalanceHistogram.Fill(value);
    }

    TPaveText diagnosticSummary(0.07, 0.08, 0.93, 0.92, "NDC");
    diagnosticSummary.SetFillColorAlpha(kWhite, 0.95);
    diagnosticSummary.SetTextAlign(12);
    diagnosticSummary.SetTextFont(42);
    std::ostringstream diagnosticCountLine;
    diagnosticCountLine << "same-R_{match} escaped endpoints: "
                        << relativeEnergyClosures.size() << '/' << runCount;
    diagnosticSummary.AddText(diagnosticCountLine.str().c_str());
    if (!relativeEnergyClosures.empty()) {
        std::ostringstream energyDiagnosticLine;
        energyDiagnosticLine << "median / p95 |#deltaE| = "
            << sampleQuantile(absoluteEnergyClosures, 0.5) << " / "
            << sampleQuantile(absoluteEnergyClosures, 0.95);
        diagnosticSummary.AddText(energyDiagnosticLine.str().c_str());
        std::ostringstream momentumDiagnosticLine;
        momentumDiagnosticLine << "median / p95 #delta(P_{N}+P_{rad}) = "
            << sampleQuantile(relativeMomentumClosures, 0.5) << " / "
            << sampleQuantile(relativeMomentumClosures, 0.95);
        diagnosticSummary.AddText(momentumDiagnosticLine.str().c_str());
        std::ostringstream angularDiagnosticLine;
        angularDiagnosticLine << "median / p95 #delta(J_{N}+J_{rad}) = "
            << sampleQuantile(relativeAngularMomentumClosures, 0.5) << " / "
            << sampleQuantile(relativeAngularMomentumClosures, 0.95);
        diagnosticSummary.AddText(angularDiagnosticLine.str().c_str());
        std::ostringstream schottDiagnosticLine;
        schottDiagnosticLine << "max |E_{S}|/E_{scale} = "
                             << maximumSchottFraction;
        diagnosticSummary.AddText(schottDiagnosticLine.str().c_str());
        diagnosticSummary.AddText(
            "The three residuals above are ALGEBRAIC IDENTITIES: E_{bound},");
        diagnosticSummary.AddText(
            "P_{bound}, J_{bound} are defined as the residual that closes them.");
        diagnosticSummary.AddText("Independent physical residuals:");
        std::ostringstream reactionDiagnosticLine;
        reactionDiagnosticLine << "  median / p95 |#deltaE_{LL-vs-flux}|/E_{rad} = "
            << sampleQuantile(reactionMismatchFractions, 0.5) << " / "
            << sampleQuantile(reactionMismatchFractions, 0.95);
        diagnosticSummary.AddText(reactionDiagnosticLine.str().c_str());
        std::ostringstream reservoirDiagnosticLine;
        reservoirDiagnosticLine << "  median / p95 |E_{bound}|/E_{rad} = "
            << sampleQuantile(boundReservoirFractions, 0.5) << " / "
            << sampleQuantile(boundReservoirFractions, 0.95);
        diagnosticSummary.AddText(reservoirDiagnosticLine.str().c_str());
        diagnosticSummary.AddText(
            "A reservoir comparable to E_{rad} means the reported radiated");
        diagnosticSummary.AddText(
            "energy is not resolved better than that fraction.");
    } else {
        diagnosticSummary.AddText("No escaped endpoint pair: closure histograms are empty");
    }
    diagnosticSummary.AddText("Cutoff, captured and unresolved outcomes are excluded");
    diagnosticSummary.AddText("P_{rad}, J_{rad}: integrated Maxwell flux on the control wavefront");
    diagnosticSummary.AddText("P_{N}, J_{N}: approximate near-field/particle Noether sector");
    if (configuration.maximumAsymptoticBeta > 0.1) {
        diagnosticSummary.AddText("WARNING: outside the intended low-velocity regime");
    }
    if (matchingRegionWarning) {
        diagnosticSummary.AddText("WARNING: matching sphere may be too close");
    }

    diagnosticsPage.cd(1);
    gPad->SetGrid();
    energyBalanceHistogram.Draw("HIST");
    // The plotted variable is an algebraic identity of the discrete update:
    // boundField{Energy,Momentum,AngularMomentum} is defined in
    // integrateElectrodynamicStep as exactly the residual that closes it.  The
    // box therefore states that plainly and quotes the independent physical
    // residual next to it, so the panel cannot be read as a conservation test.
    const auto drawBeamResidualFit = [&](TH1D& histogram,
                                         const std::vector<double>& values,
                                         const std::string& functionName,
                                         const std::string& independentLabel,
                                         const std::vector<double>& independent) {
        const GaussianFitSummary fit = gaussianMaximumLikelihood(values);
        std::unique_ptr<TF1> curve = gaussianMleOverlay(
            functionName, fit, histogram, plot_style::crem());
        if (curve) beamAnalysisFunctions.push_back(std::move(curve));
        std::vector<std::string> lines;
        if (fit.count == 0) {
            lines = {"Entries: N = 0",
                     "Fit: unavailable - no escaped endpoint pair",
                     "Experimental:",
                     "n/a - numerical bookkeeping diagnostic"};
        } else {
            lines = {
                "Entries: N = " + std::to_string(fit.count),
                fit.sigma > 0.0 ? "Fit: Gaussian MLE in plotted variable"
                                : "Fit: degenerate Gaussian (#sigma=0)",
                "#mu = " + compactNumber(fit.mean),
                "#sigma = " + compactNumber(fit.sigma),
                "IDENTITY by construction: E_{bound} is defined",
                "as the residual that closes this sum;",
                "this measures roundoff + endpoint interpolation.",
                "Independent physical residual:",
                independent.empty()
                    ? independentLabel + " = n/a"
                    : independentLabel + " median/p95 = "
                        + compactNumber(sampleQuantile(independent, 0.5))
                        + " / " + compactNumber(sampleQuantile(independent, 0.95))
            };
        }
        drawAnalysisBox(beamAnalysisBoxes, 0.44, 0.55, 0.96, 0.91,
                        lines, 0.019);
    };
    drawBeamResidualFit(energyBalanceHistogram, relativeEnergyClosures,
                        "beam_energy_balance_gaussian",
                        "|#deltaE_{LL-vs-flux}|/E_{rad}",
                        reactionMismatchFractions);
    diagnosticsPage.cd(2);
    gPad->SetGrid();
    momentumBalanceHistogram.Draw("HIST");
    drawBeamResidualFit(momentumBalanceHistogram, logMomentumClosures,
                        "beam_momentum_balance_gaussian",
                        "|E_{bound}|/E_{rad}", boundReservoirFractions);
    diagnosticsPage.cd(3);
    gPad->SetGrid();
    angularMomentumBalanceHistogram.Draw("HIST");
    drawBeamResidualFit(angularMomentumBalanceHistogram,
                        logAngularMomentumClosures,
                        "beam_angular_momentum_balance_gaussian",
                        "|E_{bound}|/E_{rad}", boundReservoirFractions);
    diagnosticsPage.cd(4);
    diagnosticSummary.Draw();

    canvas.cd();
    distributionsPage.Pop();
    canvas.Modified();
    canvas.Update();
    std::vector<root_export::NamedPad> plotsToSave{
        {distributionsPage.GetPad(1), 1, 1, "differential_cross_section"},
        {distributionsPage.GetPad(2), 1, 2, "cumulative_cross_section"},
        {diagnosticsPage.GetPad(1), 2, 1, "diagnostic_energy_balance"},
        {diagnosticsPage.GetPad(2), 2, 2, "diagnostic_momentum_balance"},
        {diagnosticsPage.GetPad(3), 2, 3, "diagnostic_angular_momentum_balance"}
    };
    if (energyLossSpectrum) {
        plotsToSave.push_back(
            {distributionsPage.GetPad(3), 1, 3, "energy_loss_cross_section"});
    }
    reportExports(root_export::saveStatisticalPlots(
        selectedPhenomenon, plotsToSave));

    bool persistenceOk = reportArchiveOperation(
        statistics_archive::writeScientificReferencesText(),
        "scientific-reference catalogue");
    if (failed != 0) return 2;
    return persistenceOk ? 0 : 3;
}

int showInteractionStatistics(std::uint64_t seed, int runCount,
                              double meanEnergyEv, double energySigmaEv,
                              double impactSigmaPm) {
    const InteractionConfiguration configuration = makeInteractionConfiguration(
        meanEnergyEv, energySigmaEv, impactSigmaPm);
    std::cout << "Interaction experiment: EACH particle's lab kinetic energy "
                 "~ N(" << meanEnergyEv << ", "
              << energySigmaEv << ") eV truncated at "
              << configuration.minimumKineticEnergy/eCharge
              << " eV, drawn independently, so the pair carries a net momentum "
                 "and the centre of mass moves; K_CM below is the invariant "
                 "derived from the two, and each trajectory is integrated in "
                 "that COM frame.  b ~ |N(0, "
              << configuration.impactParameterSigma*1.0e12
              << " pm)|, truncated by rejection at R_match/2; both dipoles "
                 "isotropic.\n"
              << "Start/escape separation " << configuration.matchingRadius*1.0e12
              << " pm; collision boundary "
              << configuration.collisionRadius*1.0e12
              << " pm; bound states observed for "
              << configuration.boundObservationOrbits << " orbits (capped at "
              << configuration.boundObservationTimeCap*1.0e15 << " fs).\n";
    const std::vector<InteractionEvent> events =
        runInteractionExperiment(seed, configuration, runCount);

        // One colour per classification, used consistently in every panel so a
    // class can be recognised across the whole page without a second legend.
    // Classification lives in the FILL, in pale tints; the outline carries the
    // provenance colour, so a panel says both "where these numbers came from"
    // and "which class each part is" without the two competing.
    std::array<int,6> outcomeColours{};
    for (std::size_t slot = 0; slot < outcomeColours.size(); ++slot) {
        outcomeColours[slot] = plot_style::classificationFill(slot);
    }
    constexpr std::array<InteractionOutcome,6> outcomeOrder{
        InteractionOutcome::Collision, InteractionOutcome::Scattering,
        InteractionOutcome::ParaPositronium,
        InteractionOutcome::OrthoPositronium,
        InteractionOutcome::Unresolved, InteractionOutcome::NumericalFailure};
    std::array<int,6> outcomeCounts{};
    std::vector<double> allEnergies, allImpacts;
    std::vector<double> boundEnergies, boundImpacts, boundAlignments;
    std::vector<double> boundAlignmentSpreads;
    std::vector<double> identityResiduals, reactionFractions;
    std::vector<double> collapseTimesPs;
    // Bound events whose observation window never closed a full Kepler orbit,
    // so no orbit-averaged power -- and hence no secular collapse time --
    // could be formed for them.
    int shortWindowBoundEvents = 0;
    // Quasi-closed orbit of every captured pair, whether or not it qualified
    // for a collapse time: the period itself is measurable from the capture
    // energy even when the window is too short to average power over.
    std::vector<double> boundPeriodsFs, boundObservedOrbitCounts;
    // Per-class means of the sampled beam parameters.  The contrast between
    // classes is the point: capture selects small impact parameters, so a
    // single overall mean would hide the mechanism.
    std::array<double,6> energySums{}, impactSums{};
    std::array<int,6> classSamples{};
    std::array<std::vector<double>,6> classEnergies, classImpacts;
    for (const InteractionEvent& event : events) {
        for (std::size_t slot = 0; slot < outcomeOrder.size(); ++slot) {
            if (event.outcome != outcomeOrder[slot]) continue;
            ++outcomeCounts[slot];
            if (std::isfinite(event.kineticEnergyEv)
                && std::isfinite(event.impactParameter)) {
                energySums[slot] += event.kineticEnergyEv;
                impactSums[slot] += event.impactParameter*1.0e12;
                ++classSamples[slot];
                classEnergies[slot].push_back(event.kineticEnergyEv);
                classImpacts[slot].push_back(event.impactParameter*1.0e12);
            }
        }
        if (std::isfinite(event.kineticEnergyEv)) {
            allEnergies.push_back(event.kineticEnergyEv);
        }
        if (std::isfinite(event.impactParameter)) {
            allImpacts.push_back(event.impactParameter*1.0e12);
        }
        const bool isBound =
            event.outcome == InteractionOutcome::ParaPositronium
         || event.outcome == InteractionOutcome::OrthoPositronium;
        if (isBound) {
            boundEnergies.push_back(event.kineticEnergyEv);
            boundImpacts.push_back(event.impactParameter*1.0e12);
            boundAlignments.push_back(event.dipoleAlignment);
            if (std::isfinite(event.dipoleAlignmentSpread)) {
                boundAlignmentSpreads.push_back(event.dipoleAlignmentSpread);
            }
            if (std::isfinite(event.collapseTimeSeconds)) {
                collapseTimesPs.push_back(event.collapseTimeSeconds*1.0e12);
            } else if (std::isfinite(event.boundObservedOrbits)
                       && event.boundObservedOrbits < 1.0) {
                ++shortWindowBoundEvents;
            }
            if (std::isfinite(event.boundOrbitalPeriodSeconds)) {
                boundPeriodsFs.push_back(
                    event.boundOrbitalPeriodSeconds*1.0e15);
            }
            if (std::isfinite(event.boundObservedOrbits)) {
                boundObservedOrbitCounts.push_back(event.boundObservedOrbits);
            }
        }
        if (!event.diagnosticsValid) continue;
        const EndpointDiagnostics& initial = event.initialDiagnostics;
        const EndpointDiagnostics& final = event.finalDiagnostics;
        const double initialBalance = initial.mechanicalEnergy
            + initial.radiatedEnergy + initial.boundFieldEnergy;
        const double finalBalance = final.mechanicalEnergy
            + final.radiatedEnergy + final.boundFieldEnergy;
        const double energyScale = std::max({std::abs(initial.mechanicalEnergy),
            std::abs(initialBalance), 64.0*std::numeric_limits<double>::epsilon()
                *2.0*firstMass*c*c});
        const double identity = (finalBalance - initialBalance)/energyScale;
        const double radiated = std::abs(final.radiatedEnergy
                                          - initial.radiatedEnergy);
        const double reaction = std::abs(final.reactionEnergyMismatch
            - initial.reactionEnergyMismatch)
            / std::max(radiated, 1.0e-300);
        if (std::isfinite(identity)) identityResiduals.push_back(identity);
        if (std::isfinite(reaction) && radiated > 0.0) {
            reactionFractions.push_back(reaction);
        }
    }

    // The summary the experiment exists to produce.
    std::cout << "\nInteraction outcome summary (" << runCount
              << " trajectories):\n";
    for (std::size_t slot = 0; slot < outcomeOrder.size(); ++slot) {
        if (outcomeCounts[slot] == 0
            && (outcomeOrder[slot] == InteractionOutcome::Unresolved
             || outcomeOrder[slot] == InteractionOutcome::NumericalFailure)) {
            continue;
        }
        std::ostringstream line;
        line << "  " << std::left << std::setw(18)
             << interactionOutcomeName(outcomeOrder[slot]) << std::right
             << std::setw(6) << outcomeCounts[slot] << "  ("
             << std::fixed << std::setprecision(1)
             << 100.0*outcomeCounts[slot]/runCount << "%)";
        std::cout << line.str() << '\n';
    }
    const int boundTotal = outcomeCounts[2] + outcomeCounts[3];
    if (boundTotal > 0) {
        std::cout << "  bound total       " << std::setw(6) << boundTotal
                  << "; para:ortho = " << outcomeCounts[2] << ":"
                  << outcomeCounts[3]
                  << " (isotropic expectation 1:3)\n";
    }
    std::cout << "Para = parallel magnetic moments (antiparallel spins, S=0); "
                 "ortho = antiparallel moments (S=1).\n";

    const auto mean = [](const std::vector<double>& values) {
        if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
        return std::accumulate(values.begin(), values.end(), 0.0)
             / static_cast<double>(values.size());
    };
    std::cout << "\nSampled beam parameters:\n"
              << "  all trajectories   <K_CM> = " << mean(allEnergies)
              << " eV,  <b> = " << mean(allImpacts) << " pm\n";
    for (std::size_t slot = 0; slot < 4; ++slot) {
        if (classSamples[slot] == 0) continue;
        std::ostringstream line;
        line << "  " << std::left << std::setw(18)
             << interactionOutcomeName(outcomeOrder[slot]) << std::right
             << "<K_CM> = " << energySums[slot]/classSamples[slot]
             << " eV,  <b> = " << impactSums[slot]/classSamples[slot] << " pm";
        std::cout << line.str() << '\n';
    }

    if (!boundPeriodsFs.empty()) {
        const auto periodExtremes = std::minmax_element(
            boundPeriodsFs.begin(), boundPeriodsFs.end());
        std::cout << "\nQuasi-closed orbit of the captured pairs ("
                  << boundPeriodsFs.size() << " of " << boundTotal
                  << " bound events):\n"
                  << "  T  minimum " << *periodExtremes.first << " fs\n"
                  << "  T  mean    " << mean(boundPeriodsFs) << " fs\n"
                  << "  T  maximum " << *periodExtremes.second << " fs\n";
        if (!boundObservedOrbitCounts.empty()) {
            const auto orbitExtremes = std::minmax_element(
                boundObservedOrbitCounts.begin(),
                boundObservedOrbitCounts.end());
            std::cout << "  revolutions observed: minimum "
                      << *orbitExtremes.first << ", mean "
                      << mean(boundObservedOrbitCounts) << ", maximum "
                      << *orbitExtremes.second << '\n';
        }
        std::cout << "T is the Kepler period of the osculating capture orbit, "
                     "the longer of its capture-time and window-close values.\n"
                     "The revolution count is what the bound observation window "
                     "actually covered, NOT a count to collapse: these pairs "
                     "are\nfollowed only long enough to classify them, so the "
                     "window is bounded by boundObservationTimeCap.\n";
    }
    if (!collapseTimesPs.empty()) {
        const auto extremes = std::minmax_element(
            collapseTimesPs.begin(), collapseTimesPs.end());
        std::cout << "\nBound-state CREM collapse time (" << collapseTimesPs.size()
                  << " of " << boundTotal << " bound events):\n"
                  << "  minimum " << *extremes.first << " ps\n"
                  << "  mean    " << mean(collapseTimesPs) << " ps\n"
                  << "  maximum " << *extremes.second << " ps\n"
                  << "This is the extrapolated classical inspiral time to "
                     "0.01*a0 from the orbit-averaged radiated power (secular\n"
                     "dE/dt=P model; experiments 1 and 2 now measure this "
                     "mechanically instead, so the two are NOT comparable).  It "
                     "is not a quantum\nannihilation lifetime.  Only bound "
                     "events whose observation window closed a full Kepler "
                     "orbit are included:\nbelow that the window samples the "
                     "post-periapsis leg, whose instantaneous power is not an "
                     "orbit average.\n";
        if (shortWindowBoundEvents > 0) {
            std::cout << "  " << shortWindowBoundEvents << " of " << boundTotal
                      << " bound events were excluded by that gate (window "
                         "shorter than one orbit).\n";
        }
    } else if (boundTotal > 0) {
        std::cout << "\nNo bound event yielded a finite collapse time: "
                  << shortWindowBoundEvents << " of " << boundTotal
                  << " had an observation window shorter than one Kepler "
                     "orbit,\nso no orbit-averaged power exists for them.  "
                     "These captures are loose enough that one period exceeds "
                     "the bound observation cap.\n";
    }
    if (outcomeCounts[4] > 0 || outcomeCounts[5] > 0) {
        std::cerr << "Warning: " << outcomeCounts[4] << " unresolved and "
                  << outcomeCounts[5] << " failed trajectories are excluded "
                     "from the physical classes.\n";
    }
    if (!identityResiduals.empty()) {
        std::cout << "Bookkeeping identity residual (not a conservation test): "
                     "median " << sampleQuantile(identityResiduals, 0.5)
                  << "; independent |dE_LL-vs-flux|/E_rad median "
                  << sampleQuantile(reactionFractions, 0.5) << ".\n";
    }

    gROOT->SetBatch(kTRUE);
    root_export::preparePdfExporter();
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    TCanvas canvas("interaction_statistics",
                   "e+e- interaction classification", 1280, 900);
    canvas.SetFillColor(kWhite);
    TPad distributionsPage("interaction_distributions_page",
                           "Interaction distributions", 0.0, 0.0, 1.0, 1.0);
    TPad diagnosticsPage("interaction_diagnostics_page",
                         "Numerical diagnostics", 0.0, 0.0, 1.0, 1.0);
    for (TPad* page : {&distributionsPage, &diagnosticsPage}) {
        page->SetFillColor(kWhite);
        page->SetFillStyle(1001);
    }
    canvas.cd();
    distributionsPage.Draw();
    diagnosticsPage.Draw();
    distributionsPage.cd();
    distributionsPage.Divide(2, 3, 0.006, 0.006);
    diagnosticsPage.cd();
    diagnosticsPage.Divide(1, 1, 0.006, 0.006);
    std::vector<std::unique_ptr<TPaveText>> analysisBoxes;
    std::vector<std::unique_ptr<TF1>> analysisFunctions;

    // All six classifications, including the censored ones: with a large
    // Unresolved fraction a four-bin chart would imply the physical classes
    // share 100% of the sample when they do not.
    TH1D outcomeHistogram("interaction_outcome_summary",
        "Interaction outcome classification;;Trajectories", 6, 0.0, 6.0);
    for (int slot = 0; slot < 6; ++slot) {
        outcomeHistogram.GetXaxis()->SetBinLabel(slot + 1,
            interactionOutcomeName(outcomeOrder[static_cast<std::size_t>(slot)]));
        outcomeHistogram.SetBinContent(slot + 1,
            outcomeCounts[static_cast<std::size_t>(slot)]);
        outcomeHistogram.SetBinError(slot + 1,
            std::sqrt(static_cast<double>(
                outcomeCounts[static_cast<std::size_t>(slot)])));
    }
    outcomeHistogram.SetStats(false);
    outcomeHistogram.SetLineColor(plot_style::crem());
    outcomeHistogram.SetLineWidth(2);
    outcomeHistogram.GetXaxis()->SetLabelSize(0.035);
    outcomeHistogram.SetMinimum(0.0);
    // Headroom so the analysis box does not sit on top of the tallest bar.
    outcomeHistogram.SetMaximum(1.75*std::max(1,
        *std::max_element(outcomeCounts.begin(), outcomeCounts.end())));
    distributionsPage.cd(1);
    gPad->SetGrid();
    gPad->SetBottomMargin(0.16);
    outcomeHistogram.Draw("HIST");
    // ROOT fills a whole TH1 with a single colour, and a per-class overlay
    // histogram would draw its empty bins as a coloured line along the axis.
    // Plain boxes give one colour per class with no such artefact.
    std::vector<std::unique_ptr<TBox>> outcomeBars;
    for (std::size_t slot = 0; slot < 6; ++slot) {
        const double height = outcomeCounts[slot];
        if (!(height > 0.0)) continue;
        auto bar = std::make_unique<TBox>(
            slot + 0.02, 0.0, slot + 0.98, height);
        bar->SetFillColorAlpha(outcomeColours[slot], 0.95);
        bar->SetLineColor(plot_style::crem());
        bar->SetLineWidth(2);
        bar->Draw("l");
        outcomeBars.push_back(std::move(bar));
    }
    std::vector<std::string> summaryLines{
        "Trajectories: N = " + std::to_string(runCount)};
    for (std::size_t slot = 0; slot < 4; ++slot) {
        summaryLines.push_back(std::string(
            interactionOutcomeName(outcomeOrder[slot])) + ": "
            + std::to_string(outcomeCounts[slot]) + " ("
            + compactNumber(100.0*outcomeCounts[slot]/runCount, 3) + "%)");
    }
    if (outcomeCounts[4] + outcomeCounts[5] > 0) {
        summaryLines.push_back("gated/failed: "
            + std::to_string(outcomeCounts[4] + outcomeCounts[5]));
    }
    summaryLines.push_back("Collision: r reached "
        + compactNumber(configuration.collisionRadius*1.0e12, 3)
        + " pm (model resolution limit 0.01a_{0})");
    summaryLines.push_back("Para: parallel #mu (S=0); Ortho: antiparallel #mu");
    summaryLines.push_back(plot_style::key(true, false, false, false));
    summaryLines.push_back("outline = provenance, fill = class:");
    summaryLines.push_back("Collision, Scattering, Para, Ortho,");
    summaryLines.push_back("Unresolved, NumericalFailure (left to right)");
    summaryLines.push_back("all: #LTK_{CM}#GT = "
        + compactNumber(mean(allEnergies), 4) + " eV, #LTb#GT = "
        + compactNumber(mean(allImpacts), 4) + " pm");
    if (!collapseTimesPs.empty()) {
        const auto extremes = std::minmax_element(
            collapseTimesPs.begin(), collapseTimesPs.end());
        summaryLines.push_back("bound CREM collapse time [ps]:");
        summaryLines.push_back("  min / mean / max = "
            + compactNumber(*extremes.first, 3) + " / "
            + compactNumber(mean(collapseTimesPs), 3) + " / "
            + compactNumber(*extremes.second, 3));
        summaryLines.push_back("  classical inspiral, not annihilation");
    }
    drawAnalysisBox(analysisBoxes, 0.42, 0.42, 0.96, 0.91, summaryLines, 0.019);

    // Built from the observed collision energies rather than from
    // meanKineticEnergy +/- 4 sigma.  That formula described the sampled
    // quantity only while ONE centre-of-mass energy was drawn per event; now
    // two per-particle lab energies are drawn and the collision energy is the
    // invariant derived from them, whose mean is roughly twice the
    // per-particle mean for an equal-mass pair and whose spread is not
    // available in closed form here.
    const double observedEnergyLow = allEnergies.empty() ? 0.0
        : *std::min_element(allEnergies.begin(), allEnergies.end());
    const double observedEnergyHigh = allEnergies.empty() ? 1.0
        : *std::max_element(allEnergies.begin(), allEnergies.end());
    const double energySpan = std::max(observedEnergyHigh - observedEnergyLow,
                                       1.0e-12);
    const double energyLower = std::max(0.0,
        observedEnergyLow - 0.05*energySpan);
    const double energyUpper = observedEnergyHigh + 0.05*energySpan;
    TH1D energyHistogram("interaction_collision_energy",
        "Collision energy from independent per-particle draws"
        ";K_{CM} [eV];Trajectories",
        histogramBins(allEnergies.size()), energyLower, energyUpper);
    energyHistogram.SetStats(false);
    energyHistogram.SetLineColor(plot_style::sampled());
    energyHistogram.SetLineWidth(3);
    energyHistogram.SetFillStyle(0);
    for (double value : allEnergies) energyHistogram.Fill(value);
    distributionsPage.cd(2);
    gPad->SetGrid();
    energyHistogram.Draw("HIST");
    // Stacked rather than overlaid: overlaid components draw their empty bins
    // as a coloured line along the axis, and stacking also makes the classes
    // visibly sum to the grey total.
    THStack energyStack("interaction_energy_stack", "");
    std::vector<std::unique_ptr<TH1D>> energyByClass;
    for (std::size_t slot = 0; slot < 6; ++slot) {
        if (classEnergies[slot].empty()) continue;
        auto part = std::make_unique<TH1D>(
            ("interaction_energy_class_" + std::to_string(slot)).c_str(), "",
            histogramBins(allEnergies.size()), energyLower, energyUpper);
        part->SetDirectory(nullptr);
        part->SetStats(false);
        part->SetLineColor(plot_style::crem());
        part->SetLineWidth(1);
        part->SetFillColorAlpha(outcomeColours[slot], 0.95);
        for (double value : classEnergies[slot]) part->Fill(value);
        energyStack.Add(part.get());
        energyByClass.push_back(std::move(part));
    }
    if (!energyByClass.empty()) energyStack.Draw("HIST SAME");
    energyHistogram.Draw("HIST SAME");
    const GaussianFitSummary energyMoments =
        gaussianMaximumLikelihood(allEnergies);
    drawAnalysisBox(analysisBoxes, 0.52, 0.60, 0.95, 0.91, {
        "Sampled: N = " + std::to_string(allEnergies.size()),
        "per-particle lab input #mu / #sigma = " + compactNumber(meanEnergyEv, 4)
            + " / " + compactNumber(energySigmaEv, 4) + " eV",
        "sample K_{CM} #mu / #sigma = " + compactNumber(energyMoments.mean, 4)
            + " / " + compactNumber(energyMoments.sigma, 4) + " eV",
        plot_style::key(false, true, false, false),
        "orange outline = sampled K_{CM}, fill = class",
        "bound fraction = " + compactNumber(
            100.0*boundTotal/std::max(runCount, 1), 3) + "%",
        "#LTK_{CM}#GT bound = " + compactNumber(mean(boundEnergies), 4)
            + " eV vs all = " + compactNumber(mean(allEnergies), 4) + " eV"
    }, 0.022);

    const double impactUpper = std::max(1.0e-3,
        4.0*configuration.impactParameterSigma*1.0e12);
    TH1D impactHistogram("interaction_impact_parameter",
        "Sampled impact parameter;b [pm];Trajectories",
        histogramBins(allImpacts.size()), 0.0, impactUpper);
    impactHistogram.SetStats(false);
    impactHistogram.SetLineColor(plot_style::sampled());
    impactHistogram.SetLineWidth(3);
    impactHistogram.SetFillStyle(0);
    for (double value : allImpacts) impactHistogram.Fill(value);
    distributionsPage.cd(3);
    gPad->SetGrid();
    impactHistogram.Draw("HIST");
    THStack impactStack("interaction_impact_stack", "");
    std::vector<std::unique_ptr<TH1D>> impactByClass;
    for (std::size_t slot = 0; slot < 6; ++slot) {
        if (classImpacts[slot].empty()) continue;
        auto part = std::make_unique<TH1D>(
            ("interaction_impact_class_" + std::to_string(slot)).c_str(), "",
            histogramBins(allImpacts.size()), 0.0, impactUpper);
        part->SetDirectory(nullptr);
        part->SetStats(false);
        part->SetLineColor(plot_style::crem());
        part->SetLineWidth(1);
        part->SetFillColorAlpha(outcomeColours[slot], 0.95);
        for (double value : classImpacts[slot]) part->Fill(value);
        impactStack.Add(part.get());
        impactByClass.push_back(std::move(part));
    }
    if (!impactByClass.empty()) impactStack.Draw("HIST SAME");
    impactHistogram.Draw("HIST SAME");
    drawAnalysisBox(analysisBoxes, 0.48, 0.58, 0.95, 0.91, {
        "Sampled: N = " + std::to_string(allImpacts.size()),
        "b #approx |N(0, " + compactNumber(
            configuration.impactParameterSigma*1.0e12, 4) + " pm)|"
            + (impactSigmaPm > 0.0 ? "" : " (auto = l_{C})"),
        "half-normal, centred on a head-on collision",
        "#LTb#GT = " + compactNumber(mean(allImpacts), 4)
            + " pm, median = " + compactNumber(
                sampleQuantile(allImpacts, 0.5), 4) + " pm",
        plot_style::key(false, true, false, false),
        "orange outline = sampled input, fill = class",
        "#LTb#GT bound = " + compactNumber(mean(boundImpacts), 4) + " pm"
    }, 0.022);

    TH1D alignmentHistogram("interaction_dipole_alignment",
        "Dipole alignment of bound states;"
        "#LTcos(#mu_{e},#mu_{p})#GT over the bound phase;Events",
        histogramBins(boundAlignments.size()), -1.0, 1.0);
    alignmentHistogram.SetStats(false);
    alignmentHistogram.SetLineColor(plot_style::crem());
    alignmentHistogram.SetLineWidth(2);
    alignmentHistogram.SetFillStyle(0);
    for (double value : boundAlignments) alignmentHistogram.Fill(value);
    distributionsPage.cd(4);
    gPad->SetGrid();
    alignmentHistogram.Draw("HIST");
    // The para/ortho threshold is exactly the classification boundary, so the
    // two coloured parts of this histogram are the two bound classes.
    TH1D paraAlignmentPart(alignmentHistogram);
    TH1D orthoAlignmentPart(alignmentHistogram);
    paraAlignmentPart.SetName("interaction_alignment_para");
    orthoAlignmentPart.SetName("interaction_alignment_ortho");
    for (int bin = 1; bin <= alignmentHistogram.GetNbinsX(); ++bin) {
        const bool paraSide =
            alignmentHistogram.GetXaxis()->GetBinCenter(bin) >= 0.5;
        if (paraSide) orthoAlignmentPart.SetBinContent(bin, 0.0);
        else paraAlignmentPart.SetBinContent(bin, 0.0);
    }
    // The copies inherit the hollow fill style of the outline histogram, so
    // the solid style has to be restored or the class colours never show.
    paraAlignmentPart.SetLineColor(plot_style::crem());
    paraAlignmentPart.SetFillColorAlpha(outcomeColours[2], 0.95);
    paraAlignmentPart.SetFillStyle(1001);
    paraAlignmentPart.SetLineWidth(2);
    orthoAlignmentPart.SetLineColor(plot_style::crem());
    orthoAlignmentPart.SetFillColorAlpha(outcomeColours[3], 0.95);
    orthoAlignmentPart.SetFillStyle(1001);
    orthoAlignmentPart.SetLineWidth(2);
    THStack alignmentStack("interaction_alignment_stack", "");
    alignmentStack.Add(&orthoAlignmentPart);
    alignmentStack.Add(&paraAlignmentPart);
    alignmentStack.Draw("HIST SAME");
    alignmentHistogram.Draw("HIST SAME");
    std::unique_ptr<TLine> thresholdLine;
    if (!boundAlignments.empty()) {
        const double lineHeight = std::max(1.0,
            1.05*alignmentHistogram.GetMaximum());
        thresholdLine = std::make_unique<TLine>(0.5, 0.0, 0.5, lineHeight);
        thresholdLine->SetLineColor(plot_style::theory());
        thresholdLine->SetLineWidth(3);
        thresholdLine->SetLineStyle(2);
        thresholdLine->Draw();
    }
    const double medianSpread = boundAlignmentSpreads.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : sampleQuantile(boundAlignmentSpreads, 0.5);
    drawAnalysisBox(analysisBoxes, 0.13, 0.55, 0.62, 0.91, {
        "Bound events: N = " + std::to_string(boundAlignments.size()),
        plot_style::key(true, false, false, true),
        "vermillion dashed: para/ortho threshold at +0.5",
        "para (#geq0.5) : ortho = " + std::to_string(outcomeCounts[2]) + " : "
            + std::to_string(outcomeCounts[3]),
        "isotropic expectation 1 : 3",
        "median drift within an event = " + compactNumber(medianSpread, 3),
        "Orientation is set at capture: this classical",
        "model has no spin-relaxation channel, so the",
        "label reflects the initial random orientation."
    }, 0.019);

    // Annihilation-time spectrum for each bound classification, drawn from the
    // MEASURED decay rate of that state.  Experiment 5 decides which of the two
    // a captured pair is, so each class gets the spectrum that actually applies
    // to it.  The two lifetimes differ by a factor of about 1135, which is why
    // they cannot share one axis and get a panel each.
    std::vector<std::unique_ptr<TF1>> boundSpectra;
    std::vector<std::unique_ptr<TLine>> boundCollapseMarkers;
    const std::array<std::size_t,2> boundSlots{2, 3};
    const std::array<const char*,2> boundValueIds{
        "para_lifetime_from_rate", "ortho_lifetime_from_rate"};
    const std::array<const char*,2> boundRateIds{
        "para_decay_rate_measurement", "ortho_decay_rate_measurement"};
    const std::array<const char*,2> boundUnits{"ps", "ns"};
    const std::array<const char*,2> boundSources{
        "Al-Ramadhan & Gidley, PRL 72 (1994)",
        "Vallery et al., PRL 90 (2003)"};
    const std::array<double,2> boundToPlotUnit{1.0, 1.0e-3};
    for (std::size_t index = 0; index < boundSlots.size(); ++index) {
        const std::size_t slot = boundSlots[index];
        const statistics_archive::ScientificValue& lifetime =
            statistics_archive::scientificValue(boundValueIds[index]);
        const statistics_archive::ScientificValue& rate =
            statistics_archive::scientificValue(boundRateIds[index]);
        const double upper = 6.0*lifetime.value;
        auto spectrum = std::make_unique<TF1>(
            ("interaction_annihilation_" + std::to_string(index)).c_str(),
            "exp(-x/[0])/[0]", 0.0, upper);
        spectrum->SetParameter(0, lifetime.value);
        spectrum->SetLineColor(plot_style::experimental());
        spectrum->SetLineWidth(3);
        spectrum->SetNpx(600);
        std::ostringstream title;
        title << interactionOutcomeName(outcomeOrder[slot])
              << " annihilation-time spectrum from the measured rate;t ["
              << boundUnits[index] << "];(1/N) dN/dt ["
              << boundUnits[index] << "^{-1}]";
        spectrum->SetTitle(title.str().c_str());
        distributionsPage.cd(static_cast<int>(5 + index));
        gPad->SetGrid();
        spectrum->Draw("L");
        // Mean CREM collapse time of the events actually classified into this
        // class, converted into the panel's time unit.
        double classCollapseMean = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> classCollapse;
        for (const InteractionEvent& event : events) {
            if (event.outcome != outcomeOrder[slot]) continue;
            if (std::isfinite(event.collapseTimeSeconds)) {
                classCollapse.push_back(event.collapseTimeSeconds*1.0e12
                                        *boundToPlotUnit[index]);
            }
        }
        if (!classCollapse.empty()) {
            classCollapseMean = std::accumulate(classCollapse.begin(),
                classCollapse.end(), 0.0)/classCollapse.size();
            if (classCollapseMean > 0.0 && classCollapseMean < upper) {
                auto marker = std::make_unique<TLine>(classCollapseMean, 0.0,
                    classCollapseMean, 1.0/lifetime.value);
                marker->SetLineColor(plot_style::crem());
                marker->SetLineWidth(2);
                marker->SetLineStyle(2);
                marker->Draw();
                boundCollapseMarkers.push_back(std::move(marker));
            }
        }
        drawAnalysisBox(analysisBoxes, 0.38, 0.50, 0.95, 0.91, {
            plot_style::key(true, false, true, false),
            "Measured data, drawn analytically - no Monte Carlo.",
            std::string("Applies to the ")
                + interactionOutcomeName(outcomeOrder[slot]) + " class:",
            "events classified here = "
                + std::to_string(outcomeCounts[slot]),
            "#lambda = " + compactNumber(rate.value, 6) + " #pm "
                + compactNumber(rate.totalUncertainty, 2) + " " + rate.unit,
            "#tau_{exp} = " + compactNumber(lifetime.value) + " #pm "
                + compactNumber(lifetime.totalUncertainty) + " "
                + boundUnits[index],
            boundSources[index],
            classCollapse.empty()
                ? std::string("no CREM collapse time in this class")
                : "blue dashed: CREM collapse #LTt#GT = "
                    + compactNumber(classCollapseMean, 3) + " "
                    + boundUnits[index],
            "Classical inspiral and quantum annihilation are",
            "different processes; this sets the scale."
        }, 0.0185);
        boundSpectra.push_back(std::move(spectrum));
    }

    TPaveText diagnosticSummary(0.07, 0.20, 0.93, 0.80, "NDC");
    diagnosticSummary.SetFillColorAlpha(kWhite, 0.95);
    diagnosticSummary.SetTextAlign(12);
    diagnosticSummary.SetTextFont(42);
    diagnosticSummary.AddText("Numerical diagnostics, interaction experiment");
    std::ostringstream diagnosticCount;
    diagnosticCount << "endpoints with diagnostics: "
                    << identityResiduals.size() << '/' << runCount;
    diagnosticSummary.AddText(diagnosticCount.str().c_str());
    if (!identityResiduals.empty()) {
        std::ostringstream identityLine;
        identityLine << "median bookkeeping IDENTITY residual = "
                     << sampleQuantile(identityResiduals, 0.5)
                     << " (roundoff only, NOT a conservation test)";
        diagnosticSummary.AddText(identityLine.str().c_str());
        std::ostringstream reactionLine;
        reactionLine << "median |#deltaE_{LL-vs-flux}|/E_{rad} = "
                     << sampleQuantile(reactionFractions, 0.5)
                     << "  (independent physical residual)";
        diagnosticSummary.AddText(reactionLine.str().c_str());
    }
    diagnosticSummary.AddText(
        "E_{bound} is defined as the residual that closes the balance, so the");
    diagnosticSummary.AddText(
        "identity above cannot be quoted as evidence of energy conservation.");
    diagnosticSummary.AddText("Classical low-velocity model; not QED.");
    diagnosticsPage.cd(1);
    diagnosticSummary.Draw();

    canvas.cd();
    distributionsPage.Pop();
    canvas.Modified();
    canvas.Update();
    reportExports(root_export::saveStatisticalPlots(5, {
        {distributionsPage.GetPad(1), 1, 1, "outcome_summary"},
        {distributionsPage.GetPad(2), 1, 2, "collision_energy"},
        {distributionsPage.GetPad(3), 1, 3, "impact_parameter"},
        {distributionsPage.GetPad(4), 1, 4, "dipole_alignment"},
        {distributionsPage.GetPad(5), 1, 5, "annihilation_time_para"},
        {distributionsPage.GetPad(6), 1, 6, "annihilation_time_ortho"},
        {diagnosticsPage.GetPad(1), 2, 1, "diagnostic_summary"}
    }));
    bool persistenceOk = reportArchiveOperation(
        statistics_archive::writeScientificReferencesText(),
        "scientific-reference catalogue");
    if (outcomeCounts[5] != 0) return 2;
    return persistenceOk ? 0 : 3;
}

int showStatisticalAnalysis(std::uint64_t seed, int selectedPhenomenon,
                            int runCount,
                            double beamEnergyEv,
                            double thetaMinimumDegrees, int angleBins,
                            double impactMaximumPm, double matchingRadiusPm,
                            double interactionEnergyEv,
                            double interactionEnergySigmaEv,
                            double interactionImpactSigmaPm,
                            double cremWallClockBudgetSeconds) {
    if (selectedPhenomenon == 5) {
        return showInteractionStatistics(seed, runCount, interactionEnergyEv,
                                         interactionEnergySigmaEv,
                                         interactionImpactSigmaPm);
    }
    if (selectedPhenomenon == 1 || selectedPhenomenon == 2) {
        return showBoundDecayStatistics(seed, selectedPhenomenon, runCount,
                                        cremWallClockBudgetSeconds);
    }
    return showBeamStatistics(seed, selectedPhenomenon, runCount,
                              beamEnergyEv, thetaMinimumDegrees,
                              angleBins, impactMaximumPm, matchingRadiusPm);
}
#endif

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
#include "modules/maxwell_validation.hpp"
#endif
} // namespace

int main(int argc, char** argv) {
    std::cout
        << "CREM attribution: publications based on this program or modified "
           "versions should cite https://github.com/pytlotad/Positronium and "
           "credit Tadeusz Slawomir Pytlos "
           "(tadeusz.slawomir.pytlos@gmail.com).\n";
#ifdef POSITRONIUM_VALIDATION_EXECUTABLE
    // The suite honours --pair too.  Sweeping pairs used to mean editing
    // defaultPair in the header and rebuilding for each one, which is three
    // minutes per pair and a working tree that has to be put back afterwards.
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument != "--pair") continue;
            if (i + 1 >= argc)
                throw std::invalid_argument("missing value for --pair");
            applyPairFromOption(argv[++i]);
        }
    } catch (const std::exception& error) {
        std::cerr << "Invalid command-line option: " << error.what() << '\n';
        return 1;
    }
    std::cout << "Pair: " << firstSpecies.name << " + " << secondSpecies.name
              << '\n';
    return runMaxwellSelfTest();
#else
    std::random_device seedSource;
    std::uint64_t seed = (static_cast<std::uint64_t>(seedSource()) << 32) ^ seedSource();
    bool diagnose = false;
    int selectedMode = 0;
    VisualStyle visualStyle = VisualStyle::Unselected;
    int selectedPhenomenon = 0;
    // One sample size for every statistical experiment.  The per-experiment
    // preview overrides that used to sit below made the default depend on which
    // channel was selected, which is a poor property for a number that appears
    // on every plot as "N =".
    int statisticalRuns = 1000;
    // Experiments 3 and 4 only.  Unrelated to interactionEnergyEv below: the
    // beam channel is a scattering measurement whose reported cross sections
    // scale as 1/K_CM^2, so this value fixes the axis range of every committed
    // beam plot and must not be retuned along with the experiment-5 energy.
    double beamEnergyEv = 20.0;
    double thetaMinimumDegrees = 5.0;
    int angleBins = 10;
    double impactParameterMaximumPm = 0.0;
    double matchingRadiusPm = 0.0;
    // Experiment 5.  A mean of 0.6 eV sits well below the 6.8 eV Ps binding
    // energy, so the pair has to shed correspondingly less energy to bind and
    // the captured fraction rises.  The impact-parameter width follows the
    // Coulomb length automatically, which scales as 1/K_CM.
    double interactionEnergyEv = 0.6;
    double interactionEnergySigmaEv = 0.4;
    // Fixed in absolute terms, NOT tied to the Coulomb length: l_C scales as
    // 1/K_CM, so an auto width would widen faster than the capture threshold
    // and lowering the energy would produce fewer bound states, not more.
    // Zero still selects the l_C rule for anyone who wants that scaling.
    double interactionImpactSigmaPm = 60.0;
    // Experiments 1/2 now integrate the full mechanical trajectory to the
    // collision boundary instead of extrapolating; a bound orbit near a0
    // needs a huge number of cheap orbits before the radiative loss becomes
    // visible, so each event is censored once it spends this long on the
    // wall clock rather than left to run indefinitely.
    double cremWallClockBudgetSeconds = 20.0;
    // Negative means "not stated on the command line", which is what lets the
    // startup question below stay silent for a fully specified batch run
    // instead of blocking it on stdin.
    double externalFieldMicroTesla = -1.0;
    // Amplitude of the stochastic-electrodynamics zero-point field, in units
    // of the physical level.  Off by default: it is an experiment, not part of
    // the model every committed result was produced with.
    double zeroPointScale = 0.0;
    // Band edges in units of the pair's orbital angular frequency.
    double zeroPointBandLow = 0.3;
    double zeroPointBandHigh = 3.0;
    // Mode count is a convergence knob, not physics: the real spectrum is a
    // continuum and 64 discrete waves are a sampling of it.
    int zeroPointModes = 64;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            const auto requireValue = [&](const std::string& option) {
                if (i + 1 >= argc) {
                    throw std::invalid_argument("missing value for " + option);
                }
                return std::string(argv[++i]);
            };
            if (argument == "--diagnose") {
                diagnose = true;
            } else if (argument == "--maxwell-test") {
                throw std::invalid_argument(
                    "Maxwell validation moved to ./positronium_validation");
            } else if (argument == "--no-gui") {
                // Retained for command-line compatibility. Statistical mode
                // is always batch-only and never opens a window.
            } else if (argument == "--pair") {
                applyPairFromOption(requireValue(argument));
            } else if (argument == "--seed") {
                seed = std::stoull(requireValue(argument));
            } else if (argument == "--phenomenon") {
                selectedPhenomenon = std::stoi(requireValue(argument));
            } else if (argument == "--runs") {
                statisticalRuns = std::stoi(requireValue(argument));
            } else if (argument == "--decay-events") {
                (void)requireValue(argument);
                throw std::invalid_argument(
                    "--decay-events was removed: the photon panels are exact "
                    "reference curves and no longer sample the generator");
            } else if (argument == "--stat-window-ps") {
                (void)requireValue(argument);
                throw std::invalid_argument(
                    "--stat-window-ps was removed: the CREM calibration window is "
                    "fixed by the orbit-averaged collapse estimator");
            } else if (argument == "--beam-energy-ev") {
                beamEnergyEv = std::stod(requireValue(argument));
            } else if (argument == "--theta-min-deg") {
                thetaMinimumDegrees = std::stod(requireValue(argument));
            } else if (argument == "--angle-bins") {
                angleBins = std::stoi(requireValue(argument));
            } else if (argument == "--bmax-pm") {
                impactParameterMaximumPm = std::stod(requireValue(argument));
            } else if (argument == "--matching-radius-pm") {
                matchingRadiusPm = std::stod(requireValue(argument));
            } else if (argument == "--interaction-energy-ev") {
                interactionEnergyEv = std::stod(requireValue(argument));
            } else if (argument == "--interaction-energy-sigma-ev") {
                interactionEnergySigmaEv = std::stod(requireValue(argument));
            } else if (argument == "--interaction-bsigma-pm") {
                interactionImpactSigmaPm = std::stod(requireValue(argument));
            } else if (argument == "--radiation-reaction") {
                const std::string model = requireValue(argument);
                if (model == "disabled") {
                    gRadiationReactionModel = ChargeRadiationReactionModel::disabled;
                } else if (model == "coherent") {
                    gRadiationReactionModel =
                        ChargeRadiationReactionModel::coherentElectricDipole;
                } else if (model == "individual") {
                    gRadiationReactionModel =
                        ChargeRadiationReactionModel::individualLandauLifshitz;
                } else if (model == "automatic") {
                    gRadiationReactionModel = ChargeRadiationReactionModel::automatic;
                } else {
                    throw std::invalid_argument("--radiation-reaction must be "
                        "disabled, coherent, individual or automatic");
                }
            } else if (argument == "--crem-wallclock-budget-s") {
                cremWallClockBudgetSeconds = std::stod(requireValue(argument));
                if (!(cremWallClockBudgetSeconds > 0.0)
                    || !std::isfinite(cremWallClockBudgetSeconds)) {
                    throw std::invalid_argument(
                        "--crem-wallclock-budget-s must be finite and positive");
                }
            } else if (argument == "--mode") {
                const std::string mode = requireValue(argument);
                if (mode == "visual" || mode == "1") selectedMode = 1;
                else if (mode == "statistical" || mode == "statistics" || mode == "2") {
                    selectedMode = 2;
                } else if (mode == "maxwell" || mode == "3") {
                    throw std::invalid_argument(
                        "Maxwell validation moved to ./positronium_validation");
                } else {
                    throw std::invalid_argument("mode must be visual or statistical");
                }
            } else if (argument == "--zpf-modes") {
                zeroPointModes = std::stoi(requireValue(argument));
                if (zeroPointModes <= 0)
                    throw std::invalid_argument("--zpf-modes must be positive");
            } else if (argument == "--integrator-order") {
                gIntegratorOrder = std::stoi(requireValue(argument));
                if (gIntegratorOrder != 2 && gIntegratorOrder != 4)
                    throw std::invalid_argument(
                        "--integrator-order must be 2 or 4");
            } else if (argument == "--zpf-band") {
                const std::string value = requireValue(argument);
                const std::size_t comma = value.find(',');
                if (comma == std::string::npos)
                    throw std::invalid_argument(
                        "--zpf-band needs lo,hi in units of the orbital frequency");
                zeroPointBandLow = std::stod(value.substr(0, comma));
                zeroPointBandHigh = std::stod(value.substr(comma + 1));
                if (!(zeroPointBandHigh > zeroPointBandLow)
                    || !(zeroPointBandLow > 0.0))
                    throw std::invalid_argument("--zpf-band needs 0 < lo < hi");
            } else if (argument == "--zpf") {
                zeroPointScale = std::stod(requireValue(argument));
                if (!std::isfinite(zeroPointScale) || zeroPointScale < 0.0) {
                    throw std::invalid_argument(
                        "--zpf scale must be finite and non-negative");
                }
            } else if (argument == "--external-field") {
                externalFieldMicroTesla = std::stod(requireValue(argument));
                if (!std::isfinite(externalFieldMicroTesla)
                    || externalFieldMicroTesla < 0.0) {
                    throw std::invalid_argument(
                        "external field must be finite and non-negative [uT]");
                }
            } else if (argument == "--visual-style") {
                const std::string style = requireValue(argument);
                if (style == "line" || style == "1") visualStyle = VisualStyle::Line;
                else if (style == "dot" || style == "2") visualStyle = VisualStyle::Dot;
                else throw std::invalid_argument("visual style must be line or dot");
            } else if (argument.rfind("--", 0) == 0) {
                throw std::invalid_argument("unknown option " + argument);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "Invalid command-line option: " << error.what() << '\n';
        return 1;
    }
    // Settle the mode first: --diagnose implies the visual engine and so is a
    // fully specified run whenever it also names a phenomenon.  Leaving this
    // after the question below made every --diagnose run stop on it.
    if (diagnose) selectedMode = 1;
    // Asked before every other question, and only when some other question is
    // going to be asked anyway.  A run that states both --mode and
    // --phenomenon is a batch run and must never stop here waiting for input.
    if (externalFieldMicroTesla < 0.0
        && (selectedMode == 0 || selectedPhenomenon == 0)) {
        std::cout << "Move the pair in an external magnetic field of 50 uT?\n"
                  << "1 -> No  (as before: no external field)\n"
                  << "2 -> Yes (uniform 50 uT, orientation drawn at random)\n"
                  << "Selection [1-2]: " << std::flush;
        int selectedExternalField = 0;
        if (!(std::cin >> selectedExternalField)
            || selectedExternalField < 1 || selectedExternalField > 2) {
            std::cerr << "Invalid selection.\n";
            return 1;
        }
        externalFieldMicroTesla =
            selectedExternalField == 2 ? earthScaleMagneticField*1.0e6 : 0.0;
    }
    if (externalFieldMicroTesla > 0.0) {
        // Isotropic direction, drawn from the run seed so the orientation is
        // part of what --seed reproduces rather than a hidden extra input.
        std::mt19937_64 fieldRandom(seed ^ 0x9e3779b97f4a7c15ULL);
        std::uniform_real_distribution<double> uniformField(0.0, 1.0);
        const double cosine = 2.0*uniformField(fieldRandom) - 1.0;
        const double azimuth = 2.0*pi*uniformField(fieldRandom);
        const double transverse = std::sqrt(std::max(0.0, 1.0 - cosine*cosine));
        gExternalMagneticField = Vec3{transverse*std::cos(azimuth),
                                      transverse*std::sin(azimuth), cosine}
                               * (externalFieldMicroTesla*1.0e-6);
        std::cout << "External magnetic field: " << externalFieldMicroTesla
                  << " uT along (" << gExternalMagneticField.x/gExternalMagneticField.norm()
                  << ", " << gExternalMagneticField.y/gExternalMagneticField.norm()
                  << ", " << gExternalMagneticField.z/gExternalMagneticField.norm()
                  << ").  At this strength the orbit is untouched and the "
                     "visible channel is dipole precession.\n";
    } else {
        std::cout << "External magnetic field: none.\n";
    }
    if (zeroPointScale > 0.0) {
        // Band centred on the pair's own orbital angular frequency, which for
        // a Kepler orbit at the pair Bohr radius is sqrt(k|q1 q2|/(mu a^3)).
        const double orbitRadius = pairBohrRadius(activePair);
        const double orbitalFrequency = std::sqrt(pairCoulombStrength
            /(pairReducedMass*orbitRadius*orbitRadius*orbitRadius));
        // 0.3 to 3 times that.  The upper edge is set by what the trajectory
        // integrator can resolve, not by physics: its step follows the error
        // probe, so a wider band costs a proportionally shorter step.  The
        // secular exchange lives near resonance anyway.
        // Factors, not frequencies: the band is now defined against the
        // pair's OSCULATING orbital frequency and is re-evaluated at every
        // force call, so it stays in resonance for the whole inspiral.  With
        // a fixed band the orbit outran it -- the measured period falls from
        // 0.327 fs to 0.0031 fs over a collapse, a factor of 105 -- which is
        // why the fixed-band lifetime depended so strongly on where the upper
        // edge was cut.
        gZeroPointField = makeZeroPointField(zeroPointBandLow,
            zeroPointBandHigh, zeroPointModes, zeroPointScale, seed);
        // Root-mean-square field, not the per-mode amplitude: N equal-energy
        // modes with <cos^2> = 1/2 give E_rms = amplitude*sqrt(N/2).
        const double initialAmplitude = gZeroPointField.amplitudeCoefficient
            *orbitalFrequency*orbitalFrequency
            *std::sqrt(zeroPointModes/2.0);
        std::cout << "Zero-point field: scale " << zeroPointScale << ", "
                  << zeroPointModes << " modes over ["
                  << zeroPointBandLow << ", " << zeroPointBandHigh
                  << "] x the osculating orbital frequency, which starts at "
                  << orbitalFrequency << " rad/s and RISES as the orbit "
                     "tightens, so the band follows it.\n"
                  << "  Starting E_rms " << initialAmplitude
                  << " V/m against a binding field of "
                  << pairCoulombStrength/(eCharge*orbitRadius*orbitRadius)
                  << " V/m; it grows as the square of the orbital frequency.\n"
                  << "  This is the FLUCTUATING half of stochastic "
                  "electrodynamics; it feeds energy in, and the radiation "
                     "reaction is the dissipative half already present.\n";
    }
    {
        const char* reactionModelName =
            gRadiationReactionModel == ChargeRadiationReactionModel::disabled
                ? "disabled (no orbital energy loss; the classical inspiral "
                  "does not run)"
            : gRadiationReactionModel
                    == ChargeRadiationReactionModel::coherentElectricDipole
                ? "coherent (Abraham-Lorentz on the pair's electric dipole)"
            : gRadiationReactionModel
                    == ChargeRadiationReactionModel::individualLandauLifshitz
                ? "individual (reduced-order Landau-Lifshitz per particle)"
                : "automatic (blended individual/coherent)";
        std::cout << "Charge radiation reaction: " << reactionModelName << ".\n";
        std::cout << "Trajectory composition order: " << gIntegratorOrder
                  << (gIntegratorOrder==4
                      ? " (Yoshida composition of three symmetric steps)"
                      : " (bare symmetric step)") << ".\n";
        std::cout << "Pair: " << firstSpecies.name << " + "
                  << secondSpecies.name << " (reduced mass "
                  << pairReducedMass/electronMass << " m_e, Bohr radius "
                  << pairBohrRadius(activePair)/bohrRadius << " a0, binding "
                  << pairBindingEnergy(activePair)/eCharge << " eV).\n";
    }
    if (selectedMode == 0) {
        std::cout << "Choose simulation mode:\n"
                  << "1 -> Visual simulation\n"
                  << "2 -> Statistical analysis\n"
                  << "Selection [1-2]: " << std::flush;
        if (!(std::cin >> selectedMode) || selectedMode < 1 || selectedMode > 2) {
            std::cerr << "Invalid selection. Enter 1 or 2.\n";
            return 1;
        }
    }
    if (selectedMode == 1 && visualStyle == VisualStyle::Unselected && !diagnose) {
        int selectedVisualStyle = 0;
        std::cout << "Choose visual simulation style:\n"
                  << "1 -> Line\n"
                  << "2 -> Dot\n"
                  << "Selection [1-2]: " << std::flush;
        if (!(std::cin >> selectedVisualStyle)
            || selectedVisualStyle < 1 || selectedVisualStyle > 2) {
            std::cerr << "Invalid selection. Enter 1 or 2.\n";
            return 1;
        }
        visualStyle = selectedVisualStyle == 1 ? VisualStyle::Line : VisualStyle::Dot;
    }
    if (visualStyle == VisualStyle::Unselected) visualStyle = VisualStyle::Line;
    // Experiment 5 exists only in Statistical: Visual integrates one prepared
    // trajectory, whereas 5 classifies an ensemble of collision outcomes.
    const int maximumPhenomenon = selectedMode == 2 ? 5 : 4;
    if (selectedPhenomenon == 5 && selectedMode != 2) {
        std::cerr << "Experiment 5 (Interactions) exists only in statistical "
                     "mode; use --mode statistical.\n";
        return 1;
    }
    if (selectedPhenomenon < 1 || selectedPhenomenon > maximumPhenomenon) {
        if (selectedMode == 2) {
            std::cout << "Choose statistical experiment:\n"
                      << "1 -> Para-positronium CREM collapse + 2 gamma kinematics\n"
                      << "2 -> Ortho-positronium CREM collapse + 3 gamma kinematics\n"
                      << "3 -> e+e- beam: short-range/cutoff channel\n"
                      << "4 -> e+e- beam: elastic scattering\n"
                      << "5 -> Interactions: classify collision outcomes\n";
        } else {
            std::cout << "Choose phenomenon to simulate:\n"
                      << "1 -> Para-positronium\n"
                      << "2 -> Ortho-positronium\n"
                      << "3 -> Direct collision\n"
                      << "4 -> Scattering\n";
        }
        std::cout << "Selection [1-" << maximumPhenomenon << "]: " << std::flush;
        if (!(std::cin >> selectedPhenomenon) || selectedPhenomenon < 1
            || selectedPhenomenon > maximumPhenomenon) {
            std::cerr << "Invalid selection. Enter a number from 1 to "
                      << maximumPhenomenon << ".\n";
            return 1;
        }
    }
    if (selectedMode == 2) {
        // Phenomena 1/2 now mechanically integrate every trajectory to the
        // collision boundary (orbit-averaged, but still two full-orbit
        // measurements per checkpoint) instead of a closed-form
        // extrapolation, so the old 10000 ceiling could mean hours of
        // wall-clock time; 1000 keeps a full run in the tens-of-minutes
        // range while still giving smooth statistics.
        const int maximumStatisticalRuns=selectedPhenomenon<=2?1000
            :(selectedPhenomenon==5?100000:100000);
        if (statisticalRuns < 1 || statisticalRuns > maximumStatisticalRuns) {
            std::cerr << "The number of CREM trajectories/beam trials must be from 1 to "
                      <<maximumStatisticalRuns<<" for this experiment.\n";
            return 1;
        }
        try {
            return showStatisticalAnalysis(seed, selectedPhenomenon,
                                           statisticalRuns,
                                           beamEnergyEv,
                                           thetaMinimumDegrees,
                                           angleBins, impactParameterMaximumPm,
                                           matchingRadiusPm,
                                           interactionEnergyEv,
                                           interactionEnergySigmaEv,
                                           interactionImpactSigmaPm,
                                           cremWallClockBudgetSeconds);
        } catch (const std::exception& error) {
            std::cerr << "Invalid statistical experiment: " << error.what() << '\n';
            return 1;
        }
    }
    // Diagnostics need only the completed trajectory. Interactive Visual uses
    // a one-step probe to obtain deterministic initial conditions, then runs
    // the real integrator only after the ROOT canvas has been painted.
    SimulationOptions initialOptions;
    if (!diagnose) {
        initialOptions.frameCount = 2;
        initialOptions.observationTime = 1.0e-24;
    }
    SimulationResult simulation = simulate(seed, selectedPhenomenon,
                                           initialOptions);
    const std::vector<Frame>& frames = simulation.frames;
    const InitialConditions initialConditions = simulation.initial;
    if (simulation.outcome == SimulationOutcome::NumericalFailure) {
        std::cerr << "Simulation stopped because the numerical state became invalid.\n";
        return 2;
    }
    if (diagnose) {
        std::cout << phenomenonName(initialConditions.phenomenon) << " simulated for "
                  << simulation.elapsedTime << " s; "
                  << frames.size() << " animation frames.\n";
        const double initialTotal=frames.front().mechanicalEnergy
            +frames.front().radiatedEnergy+frames.front().boundFieldEnergy;
        const double finalTotal=frames.back().mechanicalEnergy
            +frames.back().radiatedEnergy+frames.back().boundFieldEnergy;
        const double energyDrift = finalTotal - initialTotal;
        const double energyScale = std::max(std::abs(frames.back().radiatedEnergy), std::abs(initialTotal));
        const double relativeEnergyDrift = std::abs(energyDrift) / energyScale;
        const Vec3 momentumDrift =
            (frames.back().noetherMomentum+frames.back().radiatedMomentum
                +frames.back().boundFieldMomentum)
          -(frames.front().noetherMomentum+frames.front().radiatedMomentum
                +frames.front().boundFieldMomentum);
        const Vec3 angularMomentumDrift =
            (frames.back().noetherAngularMomentum
                +frames.back().radiatedAngularMomentum
                +frames.back().boundFieldAngularMomentum)
          - (frames.front().noetherAngularMomentum
             +frames.front().radiatedAngularMomentum
             +frames.front().boundFieldAngularMomentum);
        const double relativeMomentumDrift = momentumDrift.norm()
            / std::max(frames.front().canonicalMomentumScale, firstMass * c * 1.0e-15);
        const double angularMomentumScale =
            std::max(hbar, frames.front().noetherAngularMomentum.norm());
        const double relativeAngularMomentumDrift =
            angularMomentumDrift.norm() / angularMomentumScale;
        // "Ends closer than it started" is the right statement for a plunging
        // collision and its mirror image for scattering, but NOT for a bound
        // orbit.  It only held while the sampler drew exclusively sub-circular
        // tangential speeds, which put every bound pair at apoapsis on the
        // first frame so the radius could only decrease.  Now that the band
        // straddles the circular orbit, a pair can start near periapsis and
        // legitimately move outward.  The physical statement for a bound orbit
        // is that it stays inside its own apoapsis, a(1+e), which the slow
        // inspiral can only shrink.
        const bool boundPhenomenon =
            initialConditions.phenomenon == Phenomenon::ParaPositronium
            || initialConditions.phenomenon == Phenomenon::OrthoPositronium;
        bool expectedMotion;
        if (boundPhenomenon) {
            const double reducedMass =
                firstMass * secondMass / (firstMass + secondMass);
            const double specificEnergy =
                initialConditions.relativeEnergy / reducedMass;
            const double specificAngularMomentum =
                initialConditions.orbitalAngularMomentum / reducedMass;
            const double attractionParameter =
                pairCoulombStrength / reducedMass;
            const double eccentricity = std::sqrt(std::max(0.0, 1.0
                + 2.0 * specificEnergy * specificAngularMomentum
                    * specificAngularMomentum
                    / (attractionParameter * attractionParameter)));
            const double semiMajorAxis =
                -attractionParameter / (2.0 * specificEnergy);
            const double apoapsis = semiMajorAxis * (1.0 + eccentricity);
            const auto radiusExtremes = std::minmax_element(
                frames.begin(), frames.end(),
                [](const Frame& a, const Frame& b) { return a.radius < b.radius; });
            // 1% margin absorbs the finite frame sampling and the Darwin and
            // retardation corrections the Kepler apoapsis does not contain.
            expectedMotion = specificEnergy < 0.0
                && radiusExtremes.second->radius <= 1.01 * apoapsis;
        } else {
            expectedMotion =
                initialConditions.phenomenon == Phenomenon::Scattering
                ? frames.back().radius > frames.front().radius
                : frames.back().radius < frames.front().radius;
        }
        double maximumRelativeDipoleNormDrift = 0.0;
        for (const Frame& frame : frames) {
            maximumRelativeDipoleNormDrift = std::max({
                maximumRelativeDipoleNormDrift,
                std::abs(frame.firstDipole.norm() / firstMagneticMoment - 1.0),
                std::abs(frame.secondDipole.norm() / secondMagneticMoment - 1.0)});
        }
        const bool dipoleNormsValid = std::isfinite(maximumRelativeDipoleNormDrift)
                                  && maximumRelativeDipoleNormDrift < 1.0e-12;
        // At a direct collision the point-particle model terminates in its
        // singular regime, so only pre-contact motion is a meaningful test.
        const bool directCollision=
            initialConditions.phenomenon==Phenomenon::DirectCollision;
        const bool trajectoryValid = std::isfinite(frames.back().radius)
                                  && expectedMotion
                                  &&(directCollision||dipoleNormsValid);
        const auto radiusBounds = std::minmax_element(frames.begin(), frames.end(),
            [](const Frame& a, const Frame& b) { return a.radius < b.radius; });
        std::cout << std::setprecision(8)
                  << "seed:           " << initialConditions.seed << "\n"
                  << "phenomenon:     " << phenomenonName(initialConditions.phenomenon) << "\n"
                  << "relative E:     " << initialConditions.relativeEnergy / eCharge << " eV\n"
                  << "orbital L:      " << initialConditions.orbitalAngularMomentum / hbar << " hbar\n"
                  << "predicted rmin: " << initialConditions.predictedClosestApproach * 1.0e12 << " pm\n"
                  << "time to cutoff: ";
        if (std::isinf(initialConditions.timeToCutoff)) std::cout << "not reached\n";
        else std::cout << initialConditions.timeToCutoff * 1.0e12 << " ps\n";
        std::cout
                  << "initial radius: " << frames.front().radius * 1.0e12 << " pm\n"
                  << "final radius:   " << frames.back().radius * 1.0e12 << " pm\n"
                  << "radius range:   " << radiusBounds.first->radius * 1.0e12 << " .. "
                  << radiusBounds.second->radius * 1.0e12 << " pm\n"
                  << "radiated:       " << frames.back().radiatedEnergy / eCharge << " eV\n"
                  << "Schott energy:  " << frames.back().schottEnergy / eCharge << " eV\n"
                  << "bound field E:  " << frames.back().boundFieldEnergy/eCharge
                  << " eV (residual reservoir; |E_bound|/E_rad = "
                  << std::abs(frames.back().boundFieldEnergy)
                     /std::max(std::abs(frames.back().radiatedEnergy),1.0e-300)
                  << ")\n"
                  << "LL/coherent dE: " << frames.back().reactionEnergyMismatch/eCharge
                  << " eV (|dE|/E_rad = "
                  << std::abs(frames.back().reactionEnergyMismatch)
                     /std::max(std::abs(frames.back().radiatedEnergy),1.0e-300)
                  << ")\n"
                  << "identity resid: " << energyDrift / eCharge << " eV ("
                  << relativeEnergyDrift * 100.0 << "%) -- E_bound is DEFINED as\n"
                  << "                the residual that closes this sum, so a"
                     " zero here is\n"
                  << "                roundoff only, NOT a conservation test."
                     "  The two\n"
                  << "                ratios above are the honest residuals.\n"
                  << "field |P_rad|:  " << frames.back().radiatedMomentum.norm()
                  << " kg m/s\n"
                  << "field |J_rad|:  " << frames.back().radiatedAngularMomentum.norm()/hbar
                  << " hbar\n"
                  << "identity |dP|:  " << momentumDrift.norm() << " kg m/s ("
                  << relativeMomentumDrift * 100.0 << "% of characteristic p)\n"
                  << "identity |dJ|:  " << angularMomentumDrift.norm() / hbar << " hbar ("
                  << relativeAngularMomentumDrift * 100.0
                  << "% of characteristic J)\n"
                  << "                (both identities, as for the energy above)\n"
                  << "max |mu| drift: " << maximumRelativeDipoleNormDrift * 100.0
                  << "%\n"
                  << "trajectory:     " << (trajectoryValid ? "PASS" : "FAIL") << '\n';
        return trajectoryValid ? 0 : 1;
    }

    TApplication app("classical_hydrogen", &argc, argv);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kBlack);
    gStyle->SetPadColor(kBlack);

    TCanvas canvas("atom", "Classical model of positronium", 1100, 820);
    gVisualCanvas=&canvas;
    gVisualPhenomenon=selectedPhenomenon;
    gVisualExitSaveAttempted=false;
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
    // Positions in units of the PAIR's Bohr radius, so the camera frames the
    // orbit for every pair.  In hydrogen's a0 a muonium orbit is 1e-3 across
    // and the view floor below would have collapsed it to an invisible dot in
    // the middle of an empty box.  For e+e- this halves the numbers the span
    // is built from, so the default pair's camera sits twice as close -- the
    // orbit fills the view instead of occupying its middle half, which is the
    // same framing every other pair now gets.
    const double scale = 1.0 / pairBohrRadius(activePair);
    double viewSpan = 1.10;
    for (const Frame& frame : frames) {
        viewSpan = std::max({viewSpan, frame.first.norm() * scale * 1.10,
                            frame.second.norm() * scale * 1.10});
    }
    view->SetRange(-viewSpan, -viewSpan, -0.45*viewSpan,
                   viewSpan, viewSpan, 0.45*viewSpan);
    gPad->SetTheta(70);
    gPad->SetPhi(25);

    // Keep an owned snapshot: `simulation` is replaced by the complete run
    // after the canvas has been created.
    const Frame initialFrame = frames.front();
    const double initialFirstX = initialFrame.first.x * scale;
    const double initialFirstY = initialFrame.first.y * scale;
    const double initialFirstZ = initialFrame.first.z * scale;
    const double initialSecondX = initialFrame.second.x * scale;
    const double initialSecondY = initialFrame.second.y * scale;
    const double initialSecondZ = initialFrame.second.z * scale;
    TPolyLine3D firstPath(1, &initialFirstX,
                             &initialFirstY, &initialFirstZ);
    firstPath.SetLineColor(kBlue + 2);
    firstPath.SetLineWidth(2);
    TPolyLine3D secondPath(1, &initialSecondX,
                             &initialSecondY, &initialSecondZ);
    secondPath.SetLineColor(kOrange + 7);
    secondPath.SetLineWidth(2);

    TPolyMarker3D firstDots(1, 7);
    firstDots.SetMarkerColor(kBlue + 2);
    firstDots.SetMarkerSize(0.7);
    TPolyMarker3D secondDots(1, 7);
    secondDots.SetMarkerColor(kOrange + 7);
    secondDots.SetMarkerSize(0.7);
    firstDots.SetPoint(0, initialFirstX, initialFirstY, initialFirstZ);
    secondDots.SetPoint(0, initialSecondX, initialSecondY, initialSecondZ);

    if (visualStyle == VisualStyle::Line) {
        firstPath.Draw();
        secondPath.Draw("same");
    } else {
        firstDots.Draw();
        secondDots.Draw("same");
    }

    TPolyMarker3D first(1), second(1);
    first.SetMarkerStyle(20); first.SetMarkerSize(2.2); first.SetMarkerColor(kAzure + 1);
    second.SetMarkerStyle(20); second.SetMarkerSize(2.2); second.SetMarkerColor(kRed + 1);
    if (visualStyle == VisualStyle::Line) {
        first.Draw("same");
        second.Draw("same");
    }

    TPolyLine3D firstDipoleShaft(2), firstDipoleLeft(2), firstDipoleRight(2);
    TPolyLine3D secondDipoleShaft(2), secondDipoleLeft(2), secondDipoleRight(2);
    for (TPolyLine3D* arrow : {&firstDipoleShaft, &firstDipoleLeft, &firstDipoleRight}) {
        arrow->SetLineColor(kAzure + 1); arrow->SetLineWidth(3);
        if (visualStyle == VisualStyle::Line) arrow->Draw("same");
    }
    for (TPolyLine3D* arrow : {&secondDipoleShaft, &secondDipoleLeft, &secondDipoleRight}) {
        arrow->SetLineColor(kRed + 1); arrow->SetLineWidth(3);
        if (visualStyle == VisualStyle::Line) arrow->Draw("same");
    }

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
                   << "     r_{min} = " << initialConditions.predictedClosestApproach * 1.0e12 << " pm"
                   << "     Cutoff time = ";
    conditionsText << cutoffTimeLabel(initialConditions.timeToCutoff);
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
        formatTableValue(initialFrame.firstMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.secondMechanicalEnergy / eCharge),
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
            formatTableValue(f.firstMechanicalEnergy / eCharge),
            formatTableValue(f.secondMechanicalEnergy / eCharge),
            formatTableValue(f.mechanicalEnergy / eCharge)
        };
        drawBottomRow(currentBottomRow, bottomCurrentY, currentValues);
        currentBottomRad.SetText(0.84, bottomCurrentY, formatTableValue(f.radiatedEnergy / eCharge).c_str());

        const std::array<std::string, 7> deltaValues = {
            "delta",
            formatTableValue((f.time - initialFrame.time) * 1.0e12),
            "--",
            formatTableValue((f.radius - initialFrame.radius) * 1.0e12),
            formatTableValue((f.firstMechanicalEnergy - initialFrame.firstMechanicalEnergy) / eCharge),
            formatTableValue((f.secondMechanicalEnergy - initialFrame.secondMechanicalEnergy) / eCharge),
            formatTableValue((f.mechanicalEnergy - initialFrame.mechanicalEnergy) / eCharge)
        };
        drawBottomRow(deltaBottomRow, bottomDeltaY, deltaValues);
        deltaBottomRad.SetText(0.84, bottomDeltaY,
                               formatTableValue((f.radiatedEnergy - initialFrame.radiatedEnergy) / eCharge).c_str());
    };

    updateBottomRow(initialFrame);
    scene.Modified();
    canvas.Modified();
    canvas.Update();
    gSystem->ProcessEvents();

    std::size_t renderedFrame = 1;
    bool receivedLiveInitialFrame = false;
    double currentViewSpan = viewSpan;
    SimulationOptions visualOptions;
    // Stored frame sampling remains independent of the adaptive physical
    // step.  The live view below is driven by accepted solver steps, because
    // consecutive physical-time samples can be separated by millions of
    // small integration steps.
    visualOptions.frameCount = 240;
    visualOptions.stopRequested=[](){return gExitRequested;};
    using VisualClock = std::chrono::steady_clock;
    auto lastEventPump = VisualClock::now();
    auto lastRepaint = lastEventPump - std::chrono::milliseconds(40);
    const auto serviceVisualControls = [&]() {
        const auto now = VisualClock::now();
        if (now - lastEventPump >= std::chrono::milliseconds(8)) {
            gSystem->ProcessEvents();
            lastEventPump = now;
        }
        while (gSimulationPaused && !gExitRequested) {
            syncStopButton();
            controls.Update();
            gSystem->ProcessEvents();
            syncStopButton();
            gSystem->Sleep(20);
        }
    };
    const auto renderVisualFrame = [&](const Frame& f) {
        syncStopButton();
        const double firstX = f.first.x * scale;
        const double firstY = f.first.y * scale;
        const double firstZ = f.first.z * scale;
        const double secondX = f.second.x * scale;
        const double secondY = f.second.y * scale;
        const double secondZ = f.second.z * scale;
        // The canvas is seeded from a short deterministic probe so it can be
        // shown immediately.  Replace that seed with the full run's actual
        // first frame instead of retaining it as a detached history point.
        const std::size_t pointIndex = receivedLiveInitialFrame ? renderedFrame : 0;
        firstPath.SetPoint(static_cast<int>(pointIndex),
                              firstX, firstY, firstZ);
        secondPath.SetPoint(static_cast<int>(pointIndex),
                              secondX, secondY, secondZ);
        firstDots.SetPoint(static_cast<int>(pointIndex),
                              firstX, firstY, firstZ);
        secondDots.SetPoint(static_cast<int>(pointIndex),
                              secondX, secondY, secondZ);
        if (receivedLiveInitialFrame) {
            ++renderedFrame;
        } else {
            receivedLiveInitialFrame = true;
        }
        if (visualStyle == VisualStyle::Line) {
            first.SetPoint(0, firstX, firstY, firstZ);
            second.SetPoint(0, secondX, secondY, secondZ);
            const Vec3 firstPosition = f.first * scale;
            const Vec3 secondPosition = f.second * scale;
            setDipoleArrow(firstDipoleShaft, firstDipoleLeft, firstDipoleRight,
                           firstPosition, f.firstDipole);
            setDipoleArrow(secondDipoleShaft, secondDipoleLeft, secondDipoleRight,
                           secondPosition, f.secondDipole);
        }
        const double requiredSpan = 1.10 * std::max(
            f.first.norm(), f.second.norm()) * scale;
        if (requiredSpan > currentViewSpan) {
            currentViewSpan = requiredSpan;
            view->SetRange(-currentViewSpan, -currentViewSpan,
                           -0.45*currentViewSpan,
                            currentViewSpan,  currentViewSpan,
                            0.45*currentViewSpan);
        }
        updateBottomRow(f);
        scene.Modified();
        canvas.Modified();
        canvas.Update();
        gSystem->ProcessEvents();
        syncStopButton();
    };
    visualOptions.stepReady = [&](const State& state) {
        serviceVisualControls();
        if (gExitRequested) return;
        const auto now = VisualClock::now();
        if (now - lastRepaint < std::chrono::milliseconds(33)) return;
        lastRepaint = now;
        renderVisualFrame(makeFrame(state));
    };
    // A short trajectory may finish faster than one screen refresh.  Its
    // sampled frames are paced just enough to remain perceptible, while the
    // step callback above supplies progress during long integrations.
    visualOptions.frameReady = [&](const Frame& frame) {
        serviceVisualControls();
        if (gExitRequested) return;
        const auto now=VisualClock::now();
        if(now-lastRepaint<std::chrono::milliseconds(33)) return;
        renderVisualFrame(frame);
        lastRepaint=VisualClock::now();
    };
    const auto saveVisualCanvas=[&]() {
        canvas.cd();
        canvas.Modified();
        canvas.Update();
        const root_export::ExportResult screenshot =
            root_export::saveVisualScreenshot(canvas,selectedPhenomenon);
        if(screenshot) {
            std::cout<<"Saved visual screenshot: "
                     <<screenshot.path.string()<<'\n';
        } else {
            std::cerr<<"Warning: could not save visual screenshot "
                     <<screenshot.path.string()<<": "
                     <<screenshot.error<<'\n';
        }
        return screenshot.succeeded();
    };
    simulation = simulate(seed, selectedPhenomenon, visualOptions);
    if (gExitRequested) {
        if(!gVisualExitSaveAttempted) saveVisualCanvas();
        return 0;
    }
    if (simulation.outcome == SimulationOutcome::NumericalFailure) {
        std::cerr << "Simulation stopped because the numerical state became invalid.\n";
        return 2;
    }
    if (!frames.empty()) {
        // Ensure that even a very short run contributes its terminal point.
        renderVisualFrame(frames.back());
    }
    scene.Modified();
    canvas.Modified();
    canvas.Update();
    gSystem->ProcessEvents();
    std::cout << phenomenonName(initialConditions.phenomenon) << " simulated for "
              << simulation.elapsedTime << " s; "
              << frames.size() << " animation frames.\n";
    controls.cd();
    canvas.Modified();
    canvas.Update();
    saveVisualCanvas();
    app.Run();
    if(gExitRequested&&!gVisualExitSaveAttempted) saveVisualCanvas();
    gVisualCanvas=nullptr;
    return 0;
#endif
}
