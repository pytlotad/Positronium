#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
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

#include "bound_decay.hpp"
#include "objects/vector3.hpp"
#include "objects/state.hpp"
#include "parameters/physical_constants.hpp"
#include "root_export.hpp"
#include "statistics_archive.hpp"

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
constexpr double c=speedOfLight;
constexpr double eCharge=elementaryCharge;
constexpr double coulomb=coulombConstant;
double dot(const Vec3& a,const Vec3& b);
double gamma(const Vec3& velocity);

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

[[maybe_unused]] DipoleTensor lorentzBoostDipole(const DipoleTensor& dipole,
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

[[maybe_unused]] double dipoleFirstInvariant(const DipoleTensor& dipole) {
    return dipole.magnetic.squaredNorm()
        -c*c*dipole.electric.squaredNorm();
}

[[maybe_unused]] double dipoleSecondInvariant(const DipoleTensor& dipole) {
    return c*dot(dipole.electric,dipole.magnetic);
}
Vec3 unit(const Vec3& v) { return v / v.norm(); }
struct FourVector { double time=0.0; Vec3 space; }; // x^0=ct convention
struct ElectromagneticField { Vec3 electric, magnetic; };
double minkowskiDot(const FourVector& first,const FourVector& second) {
    return first.time*second.time-dot(first.space,second.space);
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
#include "fields/maxwell_validation_backend.hpp"
#endif

bool isFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

bool isFinite(const State& state) {
    return isFinite(state.electronPosition) && isFinite(state.positronPosition)
        && isFinite(state.electronVelocity) && isFinite(state.positronVelocity)
        && isFinite(state.electronAcceleration) && isFinite(state.positronAcceleration)
        && isFinite(state.electronDipole) && isFinite(state.positronDipole)
        &&isFinite(state.electronElectricDipole)
        &&isFinite(state.positronElectricDipole)
        &&isFinite(state.electronProperDipole)
        &&isFinite(state.positronProperDipole)
        && isFinite(state.radiatedMomentum)
        && isFinite(state.radiatedAngularMomentum)
        && isFinite(state.boundFieldMomentum)
        && isFinite(state.boundFieldAngularMomentum)
        && isFinite(state.reactionMomentumMismatch)
        && isFinite(state.reactionAngularMomentumMismatch)
        && std::isfinite(state.time) && std::isfinite(state.radiatedEnergy)
        && std::isfinite(state.dipoleConstraintEnergy)
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
    result.electronPosition = interpolateVector(
        before.electronPosition, after.electronPosition, fraction);
    result.positronPosition = interpolateVector(
        before.positronPosition, after.positronPosition, fraction);
    result.electronVelocity = interpolateVector(
        before.electronVelocity, after.electronVelocity, fraction);
    result.positronVelocity = interpolateVector(
        before.positronVelocity, after.positronVelocity, fraction);
    result.electronAcceleration = interpolateVector(
        before.electronAcceleration, after.electronAcceleration, fraction);
    result.positronAcceleration = interpolateVector(
        before.positronAcceleration, after.positronAcceleration, fraction);
    result.electronDipole = interpolateDipole(
        before.electronDipole, after.electronDipole, fraction);
    result.positronDipole = interpolateDipole(
        before.positronDipole, after.positronDipole, fraction);
    result.electronElectricDipole=interpolateVector(
        before.electronElectricDipole,after.electronElectricDipole,fraction);
    result.positronElectricDipole=interpolateVector(
        before.positronElectricDipole,after.positronElectricDipole,fraction);
    result.electronProperDipole=interpolateDipole(
        before.electronProperDipole,after.electronProperDipole,fraction);
    result.positronProperDipole=interpolateDipole(
        before.positronProperDipole,after.positronProperDipole,fraction);
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
    const Vec3 initial=before.electronPosition-before.positronPosition;
    const Vec3 change=(after.electronPosition-after.positronPosition)-initial;
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
    Vec3 electron, positron, electronDipole, positronDipole;
    Vec3 noetherMomentum, noetherAngularMomentum;
    Vec3 radiatedMomentum, radiatedAngularMomentum;
    double canonicalMomentumScale;
    double time, radius, radiatedEnergy, mechanicalEnergy, schottEnergy;
    double electronMechanicalEnergy, positronMechanicalEnergy;
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
    std::function<void(const Frame&)> frameReady;
    // Called after every accepted integration step.  Unlike frameReady this
    // is independent of physical-time sampling and is therefore suitable for
    // keeping an interactive visualization responsive.
    std::function<void(const State&)> stepReady;
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

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
Vec3 lerp(const Vec3& first,const Vec3& second,double fraction) {
    return first+(second-first)*fraction;
}
#endif

ChargeKinematics interpolatedCharge(const State& older, const State& newer, bool electron, double time) {
    const double span = newer.time - older.time;
    if(!(span>0.0)) {
        return {electron?newer.electronPosition:newer.positronPosition,
                electron?newer.electronVelocity:newer.positronVelocity,
                electron?newer.electronAcceleration:newer.positronAcceleration};
    }
    const double fraction=std::clamp((time-older.time)/span,0.0,1.0);
    const Vec3 oldPosition = electron ? older.electronPosition : older.positronPosition;
    const Vec3 newPosition = electron ? newer.electronPosition : newer.positronPosition;
    const Vec3 oldVelocity = electron ? older.electronVelocity : older.positronVelocity;
    const Vec3 newVelocity = electron ? newer.electronVelocity : newer.positronVelocity;
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
    const State& newer,bool electron,double time) {
    const double span=newer.time-older.time;
    const double fraction=span>0.0
        ?std::clamp((time-older.time)/span,0.0,1.0):1.0;
    return {lerp(electron?older.electronPosition:older.positronPosition,
                 electron?newer.electronPosition:newer.positronPosition,fraction),
            lerp(electron?older.electronVelocity:older.positronVelocity,
                 electron?newer.electronVelocity:newer.positronVelocity,fraction),
            lerp(electron?older.electronAcceleration:older.positronAcceleration,
                 electron?newer.electronAcceleration:newer.positronAcceleration,
                 fraction)};
}
#endif

ChargeKinematics historicalCharge(const StateHistory& history,
                                   const State& present, bool electron,
                                   double time) {
    const State& earliest = history.empty() ? present : history.front();
    if (time <= earliest.time) {
        const double delta = time - earliest.time;
        const Vec3 position = electron ? earliest.electronPosition : earliest.positronPosition;
        const Vec3 velocity = electron ? earliest.electronVelocity : earliest.positronVelocity;
        const Vec3 acceleration = electron ? earliest.electronAcceleration
                                           : earliest.positronAcceleration;
        return {position + velocity * delta + acceleration * (0.5 * delta * delta),
                velocity + acceleration * delta, acceleration};
    }

    const State* older = &earliest;
    for (const State& state : history) {
        if (state.time >= time) return interpolatedCharge(*older, state, electron, time);
        older = &state;
    }
    if (present.time > older->time) {
        return interpolatedCharge(*older, present, electron, time);
    }
    return interpolatedCharge(*older, *older, electron, time);
}

// Mutual, retarded Lienard-Wiechert field of a moving point charge.
ElectromagneticField lienardWiechertField(const Vec3& observationPosition,
                                          double observationTime,
                                          const StateHistory& history,
                                          const State& presentState,
                                          bool sourceIsElectron,
                                          double sourceCharge,
                                          double regularizationRadius=0.0) {
    ChargeKinematics source = historicalCharge(
        history, presentState, sourceIsElectron, observationTime);
    double retardedTime = observationTime
                        - (observationPosition - source.position).norm() / c;
    for (int iteration = 0; iteration < 16; ++iteration) {
        source = historicalCharge(
            history, presentState, sourceIsElectron, retardedTime);
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
    const RelativisticChargeCloud electron{-eCharge,chargeCloudRestRadius};
    const RelativisticChargeCloud positron{eCharge,chargeCloudRestRadius};
    block.clearSources();
    block.depositCloud(electron,state.electronPosition,state.electronVelocity);
    block.depositCloud(positron,state.positronPosition,state.positronVelocity);
    block.depositCovariantDipole(state.electronPosition,state.electronVelocity,
                                 state.electronDipole);
    block.depositCovariantDipole(state.positronPosition,state.positronVelocity,
                                 state.positronDipole);
    block.finalizeBoundInstantaneous();
    block.clearFields();
    std::vector<Vec3> dipoleVectorPotential(block.cells().size());
    for(int k=0;k<block.cellsPerAxis();++k)
        for(int j=0;j<block.cellsPerAxis();++j)
            for(int i=0;i<block.cellsPerAxis();++i) {
                const Vec3 position=block.cellPosition(i,j,k);
                const ElectromagneticField fromElectron=lienardWiechertField(
                    position,state.time,history,state,true,-eCharge,
                    chargeCloudRestRadius);
                const ElectromagneticField fromPositron=lienardWiechertField(
                    position,state.time,history,state,false,eCharge,
                    chargeCloudRestRadius);
                const Vec3 dipolePotential=
                    regularizedDipoleVectorPotential(position,state.electronPosition,
                        state.electronDipole,chargeCloudRestRadius)
                   +regularizedDipoleVectorPotential(position,state.positronPosition,
                        state.positronDipole,chargeCloudRestRadius);
                dipoleVectorPotential[(static_cast<std::size_t>(k)
                    *block.cellsPerAxis()+j)*block.cellsPerAxis()+i]=dipolePotential;
                block.replaceFieldAt(position,
                    fromElectron.electric+fromPositron.electric,
                    fromElectron.magnetic+fromPositron.magnetic);
            }
    block.addMagneticCurl(dipoleVectorPotential);
    // Preserve the retarded transverse content while matching the discrete
    // extended sources and eliminating grid-level magnetic monopoles.
    block.projectElectricGaussConstraint(projectionIterations);
    block.projectMagneticDivergenceConstraint(projectionIterations);
}
#endif

[[maybe_unused]] void precessDipole(Vec3& dipole,const Vec3& field,
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

double gamma(const Vec3& velocity) {
    return 1.0 / std::sqrt(std::max(1.0e-15, 1.0 - velocity.squaredNorm() / (c*c)));
}

Vec3 momentum(const Vec3& velocity, double mass) { return velocity * (gamma(velocity) * mass); }

Vec3 velocityFromMomentum(const Vec3& momentum, double mass) {
    const double gammaFromMomentum = std::sqrt(1.0 + momentum.squaredNorm() / (mass*mass*c*c));
    return momentum / (gammaFromMomentum * mass);
}

void synchronizeCovariantDipoles(State& state) {
    if(state.electronProperDipole.squaredNorm()==0.0
        &&state.electronDipole.squaredNorm()>0.0)
        state.electronProperDipole=state.electronDipole;
    if(state.positronProperDipole.squaredNorm()==0.0
        &&state.positronDipole.squaredNorm()>0.0)
        state.positronProperDipole=state.positronDipole;
    const DipoleTensor electronLab=lorentzBoostDipole(
        {{},state.electronProperDipole},state.electronVelocity);
    const DipoleTensor positronLab=lorentzBoostDipole(
        {{},state.positronProperDipole},state.positronVelocity);
    state.electronElectricDipole=electronLab.electric;
    state.electronDipole=electronLab.magnetic;
    state.positronElectricDipole=positronLab.electric;
    state.positronDipole=positronLab.magnetic;
}
Vec3 thomasBmtEffectiveField(const Vec3& velocity,
                             const ElectromagneticField& field);

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
    const auto [electronElectric,electronMagnetic]=
        field.interpolateField(before.electronPosition);
    const auto [positronElectric,positronMagnetic]=
        field.interpolateField(before.positronPosition);
    const Vec3 electronMomentum=relativisticBorisPush(
        momentum(before.electronVelocity,electronMass),-eCharge,electronMass,
        electronElectric,electronMagnetic,dt);
    const Vec3 positronMomentum=relativisticBorisPush(
        momentum(before.positronVelocity,positronMass),eCharge,positronMass,
        positronElectric,positronMagnetic,dt);
    state.electronVelocity=velocityFromMomentum(electronMomentum,electronMass);
    state.positronVelocity=velocityFromMomentum(positronMomentum,positronMass);
    state.electronPosition=before.electronPosition+state.electronVelocity*dt;
    state.positronPosition=before.positronPosition+state.positronVelocity*dt;
    state.electronAcceleration=(state.electronVelocity-before.electronVelocity)/dt;
    state.positronAcceleration=(state.positronVelocity-before.positronVelocity)/dt;
    field.clearSources();
    field.depositGaussianEsirkepov(-eCharge,chargeCloudRestRadius,
        before.electronPosition,before.electronVelocity,state.electronPosition,
        state.electronVelocity,dt);
    field.depositGaussianEsirkepov(eCharge,chargeCloudRestRadius,
        before.positronPosition,before.positronVelocity,state.positronPosition,
        state.positronVelocity,dt);
    field.depositCovariantDipoleYee(state.electronPosition,
        state.electronVelocity,state.electronDipole,chargeCloudRestRadius);
    field.depositCovariantDipoleYee(state.positronPosition,
        state.positronVelocity,state.positronDipole,chargeCloudRestRadius);
    field.finalizeBoundSources(dt);
    field.advance(dt);
    state.time+=dt;
}

void pushStateWithGridField(State& state, const MaxwellBlock& field,
                            double dt, bool reciprocalDipoles=false,
                            const State* samplingState=nullptr,
                            const ElectromagneticField& electronSelfField={},
                            const ElectromagneticField& positronSelfField={},
                            const DynamicSelfFieldCalibration* electronCalibration=nullptr,
                            const DynamicSelfFieldCalibration* positronCalibration=nullptr) {
    const State& sample=samplingState?*samplingState:state;
    auto [electronElectric,electronMagnetic]=
        field.interpolateField(sample.electronPosition);
    auto [positronElectric,positronMagnetic]=
        field.interpolateField(sample.positronPosition);
    const ElectromagneticField currentElectronSelf=electronCalibration
        ?electronCalibration->field(field,sample.electronPosition,
                                    sample.electronVelocity,-eCharge):electronSelfField;
    const ElectromagneticField currentPositronSelf=positronCalibration
        ?positronCalibration->field(field,sample.positronPosition,
                                    sample.positronVelocity,eCharge):positronSelfField;
    electronElectric=electronElectric-currentElectronSelf.electric;
    electronMagnetic=electronMagnetic-currentElectronSelf.magnetic;
    positronElectric=positronElectric-currentPositronSelf.electric;
    positronMagnetic=positronMagnetic-currentPositronSelf.magnetic;
    const GridDipoleInteraction electronDipoleInteraction=reciprocalDipoles
        ?field.dipoleInteraction(sample.electronPosition,sample.electronDipole)
        :GridDipoleInteraction{};
    const GridDipoleInteraction positronDipoleInteraction=reciprocalDipoles
        ?field.dipoleInteraction(sample.positronPosition,sample.positronDipole)
        :GridDipoleInteraction{};
    if(!reciprocalDipoles) {
        precessDipole(state.electronDipole,
            thomasBmtEffectiveField(state.electronVelocity,
                {electronElectric,electronMagnetic}),-eCharge/electronMass,0.5*dt);
        precessDipole(state.positronDipole,
            thomasBmtEffectiveField(state.positronVelocity,
                {positronElectric,positronMagnetic}),eCharge/positronMass,0.5*dt);
    }
    const Vec3 electronMomentum=relativisticBorisPush(
        momentum(state.electronVelocity,electronMass),-eCharge,electronMass,
        electronElectric,electronMagnetic,dt)+electronDipoleInteraction.force*dt;
    const Vec3 positronMomentum=relativisticBorisPush(
        momentum(state.positronVelocity,positronMass),eCharge,positronMass,
        positronElectric,positronMagnetic,dt)+positronDipoleInteraction.force*dt;
    state.electronVelocity=velocityFromMomentum(electronMomentum,electronMass);
    state.positronVelocity=velocityFromMomentum(positronMomentum,positronMass);
    state.electronPosition+=state.electronVelocity*dt;
    state.positronPosition+=state.positronVelocity*dt;
    if(reciprocalDipoles) {
        constexpr double electronGyromagneticRatio=-eCharge/(2.0*electronMass);
        constexpr double positronGyromagneticRatio=eCharge/(2.0*positronMass);
        const double electronNorm=state.electronDipole.norm();
        const double positronNorm=state.positronDipole.norm();
        state.electronDipole+=electronDipoleInteraction.torque
                            *(electronGyromagneticRatio*dt);
        state.positronDipole+=positronDipoleInteraction.torque
                            *(positronGyromagneticRatio*dt);
        if(state.electronDipole.norm()>0.0)
            state.electronDipole=state.electronDipole*(electronNorm/state.electronDipole.norm());
        if(state.positronDipole.norm()>0.0)
            state.positronDipole=state.positronDipole*(positronNorm/state.positronDipole.norm());
        state.dipoleConstraintEnergy+=(electronDipoleInteraction.fieldPower
            +positronDipoleInteraction.fieldPower
            -dot(electronDipoleInteraction.force,state.electronVelocity)
            -dot(positronDipoleInteraction.force,state.positronVelocity))*dt;
    } else {
        precessDipole(state.electronDipole,
            thomasBmtEffectiveField(state.electronVelocity,
                {electronElectric,electronMagnetic}),-eCharge/electronMass,0.5*dt);
        precessDipole(state.positronDipole,
            thomasBmtEffectiveField(state.positronVelocity,
                {positronElectric,positronMagnetic}),eCharge/positronMass,0.5*dt);
    }
    state.time+=dt;
}
#endif

#include "interactions/electrodynamics.hpp"

StateHistory causalInitialHistory(const State& initial,double spanFactor=8.0,
                                  int intervalCount=64) {
    State endpoint=initial;
    synchronizeCovariantDipoles(endpoint);
    const double lightCrossingTime=separation(endpoint)/c;
    const double historySpan=std::max(1.0e-24,spanFactor*lightCrossingTime);
    intervalCount=std::max(intervalCount,2);
    const auto buildHistory=[&](const Vec3& electronAcceleration,
                                const Vec3& positronAcceleration) {
        StateHistory result;
        for(int index=intervalCount;index>=0;--index) {
            const double offset=-historySpan*index/intervalCount;
            State sample=endpoint;
            sample.time=initial.time+offset;
            sample.electronPosition=initial.electronPosition
                +initial.electronVelocity*offset
                +electronAcceleration*(0.5*offset*offset);
            sample.positronPosition=initial.positronPosition
                +initial.positronVelocity*offset
                +positronAcceleration*(0.5*offset*offset);
            sample.electronVelocity=initial.electronVelocity
                +electronAcceleration*offset;
            sample.positronVelocity=initial.positronVelocity
                +positronAcceleration*offset;
            sample.electronAcceleration=electronAcceleration;
            sample.positronAcceleration=positronAcceleration;
            // Bookkeeping starts at the requested t=0 event, not in the
            // hidden causal preparation interval.
            sample.radiatedEnergy=initial.radiatedEnergy;
            sample.radiatedMomentum=initial.radiatedMomentum;
            sample.radiatedAngularMomentum=initial.radiatedAngularMomentum;
            synchronizeCovariantDipoles(sample);
            result.push_back(sample);
        }
        return result;
    };
    const MutualForces seedForces=allExternalForces(endpoint);
    Vec3 electronAcceleration=relativisticAcceleration(
        endpoint.electronVelocity,seedForces.electron,electronMass);
    Vec3 positronAcceleration=relativisticAcceleration(
        endpoint.positronVelocity,seedForces.positron,positronMass);
    StateHistory history=buildHistory(electronAcceleration,positronAcceleration);
    // Two inexpensive Picard updates make the hidden past consistent with
    // the same retarded interaction used at t=0 instead of freezing the
    // instantaneous Coulomb acceleration into the whole preparation span.
    for(int iteration=0;iteration<2;++iteration) {
        const MutualForces retarded=retardedExternalForces(endpoint,history);
        electronAcceleration=relativisticAcceleration(
            endpoint.electronVelocity,retarded.electron,electronMass);
        positronAcceleration=relativisticAcceleration(
            endpoint.positronVelocity,retarded.positron,positronMass);
        history=buildHistory(electronAcceleration,positronAcceleration);
    }
    return history;
}

class ClassicalTrajectoryEngine {
public:
    struct Accuracy {
        double relativeTolerance = 1.0e-7;
        int maximumDepth = 14;
        ChargeRadiationReactionModel reactionModel =
            ChargeRadiationReactionModel::individualLandauLifshitz;
    };
    explicit ClassicalTrajectoryEngine(const State& initial)
        :history_(causalInitialHistory(initial)) {}
    ClassicalTrajectoryEngine(const State& initial,Accuracy accuracy)
        :history_(causalInitialHistory(initial)),accuracy_(accuracy) {}
    ClassicalTrajectoryEngine(StateHistory history,Accuracy accuracy)
        :history_(std::move(history)),accuracy_(accuracy) {}
    bool advance(State& state,double dt) {
        if(!(dt>0.0)||!std::isfinite(dt)) return false;
        return advanceAdaptive(state,dt,0);
    }
    const StateHistory& history() const { return history_; }
private:
    static double normalizedStepError(const State& coarse,const State& fine) {
        const double lengthScale=std::max(separation(fine),nuclearCutoff);
        const double speedScale=std::max(
            (fine.electronVelocity-fine.positronVelocity).norm(),1.0e-6*c);
        return std::max({
            (coarse.electronPosition-fine.electronPosition).norm()/lengthScale,
            (coarse.positronPosition-fine.positronPosition).norm()/lengthScale,
            (coarse.electronVelocity-fine.electronVelocity).norm()/speedScale,
            (coarse.positronVelocity-fine.positronVelocity).norm()/speedScale});
    }

    bool advanceAdaptive(State& state,double dt,int depth) {
        const State start=state;
        const StateHistory startingHistory=history_;

        State coarse=start;
        integrateElectrodynamicStep(coarse,dt,startingHistory,false,
            accuracy_.reactionModel);

        State fine=start;
        StateHistory fineHistory=startingHistory;
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,false,
            accuracy_.reactionModel);
        if(isFinite(fine)) appendStateHistory(fineHistory,fine);
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,false,
            accuracy_.reactionModel);
        if(!isFinite(coarse)||!isFinite(fine)) return false;

        const double error=normalizedStepError(coarse,fine);
        if(error<=accuracy_.relativeTolerance||depth>=accuracy_.maximumDepth) {
            if(!std::isfinite(error)) return false;
            // Re-evaluate only the accepted two-half-step path with complete
            // far-zone bookkeeping. Trial paths above affect the local error
            // estimate but never become part of the physical history.
            state=start;
            history_=startingHistory;
            integrateElectrodynamicStep(state,0.5*dt,history_,true,
                accuracy_.reactionModel);
            if(!isFinite(state)) return false;
            appendStateHistory(history_,state);
            integrateElectrodynamicStep(state,0.5*dt,history_,true,
                accuracy_.reactionModel);
            if(!isFinite(state)) return false;
            appendStateHistory(history_,state);
            return true;
        }

        state=start;
        history_=startingHistory;
        return advanceAdaptive(state,0.5*dt,depth+1)
            &&advanceAdaptive(state,0.5*dt,depth+1);
    }
    StateHistory history_;
    Accuracy accuracy_;
};

#ifndef POSITRONIUM_VALIDATION_EXECUTABLE
Frame makeFrame(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const double electronKinetic = (gamma(s.electronVelocity) - 1.0) * electronMass * c*c;
    const double positronKinetic = (gamma(s.positronVelocity) - 1.0) * positronMass * c*c;
    const double coulombPotential = -coulomb * eCharge*eCharge * geometry.inverseDistance;
    const double dipolePotential = -dot(s.electronDipole,
        regularizedDipoleField(geometry.electronMinusPositron, s.positronDipole));
    const double darwinEnergy = darwinInteractionEnergy(s);
    const MutualForces forces = allExternalForces(s);
    const Vec3 electronAcceleration = relativisticAcceleration(s.electronVelocity, forces.electron, electronMass);
    const Vec3 positronAcceleration = relativisticAcceleration(s.positronVelocity, forces.positron, positronMass);
    // Low-velocity diagnostic proxy for the sum of the two individual Schott
    // near-field energies.  Cross terms belong to the mutual retarded field
    // and must not be inserted again through a collective dipole expression.
    const double schottEnergy = -eCharge * eCharge
        * (dot(electronAcceleration, s.electronVelocity)
         + dot(positronAcceleration, s.positronVelocity))
        / (6.0 * pi * epsilon0 * c*c*c);
    // Noether energy of the conservative approximate action.  The q-mu term
    // is homogeneous of degree one in both velocities and therefore cancels
    // in the Legendre transform.
    const double conservativeNoetherEnergy=conservativeParticleEnergy(s);
    return {s.electronPosition,s.positronPosition,
            s.electronProperDipole.squaredNorm()>0.0
                ?s.electronProperDipole:s.electronDipole,
            s.positronProperDipole.squaredNorm()>0.0
                ?s.positronProperDipole:s.positronDipole,
            noetherMomentum(s), noetherAngularMomentum(s),
            s.radiatedMomentum, s.radiatedAngularMomentum,
            canonicalMomentumScale(s),
            s.time, geometry.distance, s.radiatedEnergy,
            conservativeNoetherEnergy,
            schottEnergy,
            electronKinetic + 0.5 * (coulombPotential + dipolePotential
                + darwinEnergy + s.dipoleConstraintEnergy),
            positronKinetic + 0.5 * (coulombPotential + dipolePotential
                + darwinEnergy + s.dipoleConstraintEnergy),
            s.boundFieldEnergy,s.reactionEnergyMismatch,
            s.boundFieldMomentum,s.boundFieldAngularMomentum,
            s.reactionMomentumMismatch,s.reactionAngularMomentumMismatch};
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

SimulationResult simulate(std::uint64_t seed, int selectedPhenomenon,
                          SimulationOptions options = {}) {
    const double reducedMass = electronMass * positronMass / (electronMass + positronMass);
    const double circularSpeed = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * bohrRadius));
    const double escapeSpeed = std::sqrt(2.0) * circularSpeed;
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> unitRandom(0.0, 1.0);
    std::uniform_real_distribution<double> signedRandom(-1.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, 2.0*pi);
    const auto randomDirection = [&]() {
        const double z = signedRandom(random);
        const double phi = azimuth(random);
        const double radial = std::sqrt(1.0 - z*z);
        return Vec3{radial*std::cos(phi), radial*std::sin(phi), z};
    };

    State s;
    s.electronPosition = {bohrRadius * positronMass / (electronMass + positronMass), 0, 0};
    s.positronPosition = {-bohrRadius * electronMass / (electronMass + positronMass), 0, 0};

    // The selected branch chooses a physically useful sampling range while
    // all actual initial values inside that range remain random.
    // Menu order is para, ortho, direct collision, scattering. Internally the
    // samplers retain the order direct, scattering, para, ortho.
    const std::array<int, 5> scenarioForMenuChoice = {0, 2, 3, 0, 1};
    const int sampledScenario = scenarioForMenuChoice[selectedPhenomenon];
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

    const double defaultObservationTime = phenomenon == Phenomenon::DirectCollision ? 4.0e-16
                                        : phenomenon == Phenomenon::Scattering ? 2.0e-15
                                        : 1.50e-11;
    const double observationTime = options.observationTime > 0.0
                                 ? options.observationTime : defaultObservationTime;
    const int frameCount = std::max(2, options.frameCount);
    std::vector<Frame> frames;
    if (options.collectFrames) frames.reserve(frameCount + 1);
    const double frameInterval = observationTime / (frameCount - 1);
    double nextFrame = frameInterval;
    double minimumSeparation = separation(s);
    double maximumBeta = std::max(s.electronVelocity.norm(), s.positronVelocity.norm()) / c;
    double elapsedTime = s.time;
    double finalRadiatedEnergy = s.radiatedEnergy;
    SimulationOutcome outcome = SimulationOutcome::NumericalFailure;
    ClassicalTrajectoryEngine trajectory(s);
    if (options.collectFrames) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }

    while (s.time < observationTime && separation(s) > nuclearCutoff) {
        // Resolve each instantaneous orbit well enough to keep numerical
        // energy drift below the physical radiation loss.
        const State beforeStep = s;
        const double r = separation(beforeStep);
        const double omega = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * r*r*r));
        const double dt = std::min({2.0e-18, 2.0 * pi / (320.0 * omega),
                                    observationTime - s.time});
        if (!(dt > 0.0) || !std::isfinite(dt)) break;
        if(!trajectory.advance(s,dt)) break;
        const double currentSeparation = separation(s);
        if (!(currentSeparation > 0.0) || !std::isfinite(currentSeparation)) break;
        if (options.stepReady) options.stepReady(s);

        const bool reachedCutoff = currentSeparation <= nuclearCutoff;
        const double crossingFraction = reachedCutoff
            ? separationCrossingFraction(beforeStep,s,nuclearCutoff) : 1.0;
        const double validEndTime = beforeStep.time + crossingFraction * dt;
        if (options.collectFrames) {
            while (nextFrame <= validEndTime
                   && static_cast<int>(frames.size()) < frameCount) {
                const double sampleFraction = std::clamp(
                    (nextFrame - beforeStep.time) / dt, 0.0, crossingFraction);
                const State sampledState = interpolateState(beforeStep, s, sampleFraction);
                if (separation(sampledState) > nuclearCutoff) {
                    frames.push_back(makeFrame(sampledState));
                    if (options.frameReady) options.frameReady(frames.back());
                }
                nextFrame += frameInterval;
            }
        }

        if (reachedCutoff) {
            // Locate the terminal event linearly inside the last, already
            // very small adaptive step.  This avoids reporting a numerical
            // overshoot below the boundary of the point-particle model.
            elapsedTime = beforeStep.time + crossingFraction * dt;
            finalRadiatedEnergy = beforeStep.radiatedEnergy + crossingFraction
                * (s.radiatedEnergy - beforeStep.radiatedEnergy);
            const Vec3 eventElectronVelocity = beforeStep.electronVelocity
                + (s.electronVelocity - beforeStep.electronVelocity) * crossingFraction;
            const Vec3 eventPositronVelocity = beforeStep.positronVelocity
                + (s.positronVelocity - beforeStep.positronVelocity) * crossingFraction;
            minimumSeparation = std::min(minimumSeparation, nuclearCutoff);
            maximumBeta = std::max(maximumBeta,
                std::max(eventElectronVelocity.norm(), eventPositronVelocity.norm()) / c);
            outcome = SimulationOutcome::ReachedCutoff;
            break;
        }

        minimumSeparation = std::min(minimumSeparation, currentSeparation);
        maximumBeta = std::max(maximumBeta,
            std::max(s.electronVelocity.norm(), s.positronVelocity.norm()) / c);
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;
    }
    if (outcome != SimulationOutcome::ReachedCutoff
        && isFinite(s) && s.time >= observationTime) {
        outcome = SimulationOutcome::ObservationLimit;
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;
    }
    if (options.collectFrames && outcome == SimulationOutcome::ObservationLimit
        && (frames.empty() || frames.back().time < s.time)) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }
    const double timeToCutoff = outcome == SimulationOutcome::ReachedCutoff
                              ? elapsedTime : std::numeric_limits<double>::infinity();
    return {std::move(frames), {relativeEnergy, orbitalAngularMomentum,
            predictedClosestApproach, dipoleAlignment, timeToCutoff, phenomenon, seed},
            outcome, minimumSeparation, elapsedTime, finalRadiatedEnergy, maximumBeta};
}

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

std::uint64_t splitMix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

int histogramBins(size_t count) {
    return std::clamp(static_cast<int>(std::lround(2.0 * std::sqrt(count))), 6, 30);
}

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

double orePowellSpectrum(double* coordinate, double* parameters) {
    const double endpoint = parameters[1];
    if (!(endpoint > 0.0)) return 0.0;
    const double x = coordinate[0]/endpoint;
    if (!(x > 0.0) || x > 1.0) return 0.0;
    if (x >= 1.0 - 1.0e-12) return parameters[0];
    if (x < 1.0e-4) {
        return parameters[0] * (5.0*x/6.0 + x*x/6.0);
    }
    const double logarithm = std::log1p(-x);
    const double oneMinus = 1.0 - x;
    const double twoMinus = 2.0 - x;
    const double shape = x*oneMinus/(twoMinus*twoMinus)
        - 2.0*oneMinus*oneMinus*logarithm
            /(twoMinus*twoMinus*twoMinus)
        + twoMinus/x + 2.0*oneMinus*logarithm/(x*x);
    return parameters[0]*shape;
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

int showBoundDecayStatistics(std::uint64_t seed, int selectedPhenomenon,
                             int runCount) {
    const bool isPara = selectedPhenomenon == 1;
    const bound_decay::PositroniumState state = isPara
        ? bound_decay::PositroniumState::Para
        : bound_decay::PositroniumState::Ortho;
    const statistics_archive::ScientificValue& modelLifetimeReference =
        statistics_archive::scientificValue(isPara
            ? "pdg_para_lifetime_input" : "pdg_ortho_lifetime_input");
    const double referenceLifetime = modelLifetimeReference.value
        * (isPara ? 1.0e-12 : 1.0e-9);
    const double timeScale = isPara ? 1.0e12 : 1.0e9;
    const char* timeUnit = isPara ? "ps" : "ns";

    std::mt19937_64 random(seed);
    std::vector<bound_decay::DecayEvent> events;
    events.reserve(static_cast<size_t>(runCount));
    std::vector<double> decayTimes;
    std::vector<double> photonEnergies;
    std::vector<double> paraCosines;
    std::vector<double> paraIndependentCosines;
    std::vector<double> orthoMaximumEnergies;
    std::vector<double> orthoMiddleEnergies;
    std::vector<double> orthoLeadingAngles;
    std::vector<double> energyClosureEpsilon;
    std::vector<double> momentumClosureEpsilon;
    std::vector<double> photonMassShellEpsilon;
    std::vector<double> parentInvariantEpsilon;
    decayTimes.reserve(static_cast<size_t>(runCount));
    photonEnergies.reserve(static_cast<size_t>(runCount * (isPara ? 2 : 3)));

    double lifetimeSum = 0.0;
    for (int index = 0; index < runCount; ++index) {
        events.push_back(bound_decay::generateDecay(state, random));
        const bound_decay::DecayEvent& event = events.back();
        const double displayedTime = event.decayTimeSeconds * timeScale;
        decayTimes.push_back(displayedTime);
        lifetimeSum += event.decayTimeSeconds;
        double photonEnergySum = 0.0;
        bound_decay::Vec3 photonMomentumSum;
        for (size_t photon = 0; photon < event.photonCount; ++photon) {
            const bound_decay::Photon& currentPhoton = event.photons[photon];
            photonEnergies.push_back(
                bound_decay::energyKeV(currentPhoton.energyJoules));
            photonEnergySum += currentPhoton.energyJoules;
            const bound_decay::Vec3 photonMomentum = currentPhoton.momentum();
            photonMomentumSum = photonMomentumSum + photonMomentum;
            const double energySquared = currentPhoton.energyJoules
                                       * currentPhoton.energyJoules;
            const double shellResidual = (energySquared
                - bound_decay::speedOfLight*bound_decay::speedOfLight
                    * bound_decay::dot(photonMomentum, photonMomentum))
                / (energySquared * std::numeric_limits<double>::epsilon());
            photonMassShellEpsilon.push_back(shellResidual);
        }
        const double parentEnergy = bound_decay::positroniumRestEnergyJoules;
        const double machineEpsilon = std::numeric_limits<double>::epsilon();
        energyClosureEpsilon.push_back(
            (photonEnergySum - parentEnergy)/(parentEnergy*machineEpsilon));
        momentumClosureEpsilon.push_back(bound_decay::speedOfLight
            * bound_decay::norm(photonMomentumSum)/(parentEnergy*machineEpsilon));
        parentInvariantEpsilon.push_back(
            (photonEnergySum*photonEnergySum
             - bound_decay::speedOfLight*bound_decay::speedOfLight
                 * bound_decay::dot(photonMomentumSum, photonMomentumSum)
             - parentEnergy*parentEnergy)
            / (parentEnergy*parentEnergy*machineEpsilon));
        if (isPara) {
            paraIndependentCosines.push_back(event.paraPhotonCosPolar);
            paraCosines.push_back(event.paraPhotonCosPolar);
            paraCosines.push_back(-event.paraPhotonCosPolar);
        } else {
            orthoMaximumEnergies.push_back(
                bound_decay::energyKeV(event.orthoMaximumEnergyJoules));
            orthoMiddleEnergies.push_back(
                bound_decay::energyKeV(event.orthoMiddleEnergyJoules));
            orthoLeadingAngles.push_back(
                event.orthoLeadingPairAngleRadians * 180.0 / pi);
        }
    }

    const double estimatedLifetime = lifetimeSum / runCount;
    const double estimatedError = estimatedLifetime / std::sqrt(static_cast<double>(runCount));
    const statistics_archive::ScientificValue& lifetimeReference =
        statistics_archive::scientificValue(isPara
            ? "para_lifetime_from_rate" : "ortho_lifetime_from_rate");
    const statistics_archive::ScientificValue& photonEnergyReference =
        statistics_archive::scientificValue("pdg_two_photon_energy");
    const statistics_archive::ScientificValue& anisotropyReference =
        statistics_archive::scientificValue("ideal_para_anisotropy");
    const double experimentalLifetime = lifetimeReference.value;
    const double experimentalLifetimeError = lifetimeReference.totalUncertainty;
    std::cout << (isPara ? "Para-positronium" : "Ortho-positronium")
              << " ideal-vacuum decay sample: " << runCount << " events\n"
              << "Reference mean lifetime: " << referenceLifetime * timeScale << ' '
              << timeUnit << "; sample mean: " << estimatedLifetime * timeScale
              << " +/- " << estimatedError * timeScale << ' ' << timeUnit << "\n"
              << "Model: phenomenological lifetime + leading-order "
              << (isPara ? "2-gamma" : "Ore-Powell 3-gamma")
              << " kinematics; no detector or material effects.\n";
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
    distributionsPage.cd();
    distributionsPage.Divide(2, 2, 0.006, 0.006);
    diagnosticsPage.cd();
    diagnosticsPage.Divide(2, 2, 0.006, 0.006);
    std::vector<std::unique_ptr<TPaveText>> analysisBoxes;
    std::vector<std::unique_ptr<TF1>> analysisFunctions;

    const double maximumObservedTime = *std::max_element(decayTimes.begin(), decayTimes.end());
    const double lifetimeUpper = std::max(maximumObservedTime * 1.05,
                                           referenceLifetime * timeScale * 5.0);
    std::ostringstream lifetimeTitle;
    lifetimeTitle << "Annihilation-time spectrum;t [" << timeUnit << "];Events";
    TH1D lifetimeHistogram("decay_time_histogram", lifetimeTitle.str().c_str(),
        histogramBins(decayTimes.size()), 0.0, lifetimeUpper);
    styleHistogram(lifetimeHistogram, kOrange + 7);
    lifetimeHistogram.SetStats(false);
    for (double value : decayTimes) lifetimeHistogram.Fill(value);
    TF1 inputLifetime("input_lifetime_distribution", "[0]*exp(-x/[1])",
                      0.0, lifetimeUpper);
    inputLifetime.SetParameters(
        runCount * lifetimeHistogram.GetBinWidth(1)/(referenceLifetime*timeScale),
        referenceLifetime*timeScale);
    inputLifetime.SetLineColor(kRed + 1);
    inputLifetime.SetLineWidth(2);
    TF1 fittedLifetime("fitted_lifetime_distribution", "[0]*exp(-x/[1])",
                       0.0, lifetimeUpper);
    const double displayedFittedLifetime = estimatedLifetime*timeScale;
    fittedLifetime.SetParameters(
        runCount*lifetimeHistogram.GetBinWidth(1)/displayedFittedLifetime,
        displayedFittedLifetime);
    fittedLifetime.SetLineColor(kBlue + 1);
    fittedLifetime.SetLineWidth(2);
    fittedLifetime.SetLineStyle(2);
    lifetimeHistogram.SetMaximum(1.15 * std::max(
        lifetimeHistogram.GetMaximum(),
        std::max(inputLifetime.Eval(0.0), fittedLifetime.Eval(0.0))));

    const double photonEndpoint = 0.5 * bound_decay::energyKeV(
        bound_decay::positroniumRestEnergyJoules);
    TH1D photonHistogram("decay_photon_energy_histogram",
        isPara ? "Ideal 2#gamma energy line;E_{#gamma} [keV];Photons"
               : "Ore-Powell inclusive 3#gamma spectrum;E_{#gamma} [keV];Photons",
        isPara ? 40 : histogramBins(photonEnergies.size()),
        isPara ? photonEndpoint - 0.02 : 0.0,
        isPara ? photonEndpoint + 0.02 : photonEndpoint);
    styleHistogram(photonHistogram, kAzure + 1);
    photonHistogram.SetStats(false);
    for (double value : photonEnergies) photonHistogram.Fill(value);
    const GaussianFitSummary photonMoments = gaussianMaximumLikelihood(photonEnergies);
    std::unique_ptr<TF1> orePowellTemplate;
    if (!isPara) {
        orePowellTemplate = std::make_unique<TF1>(
            "ore_powell_spectrum_template", orePowellSpectrum,
            std::max(1.0e-6, photonEndpoint*1.0e-6), photonEndpoint, 2);
        const double expectedAmplitude = 3.0*runCount*photonHistogram.GetBinWidth(1)
            * 2.0/((pi*pi - 9.0)*photonEndpoint);
        orePowellTemplate->SetParameters(expectedAmplitude, photonEndpoint);
        orePowellTemplate->SetParName(0, "normalization");
        orePowellTemplate->SetParName(1, "endpoint_keV");
        orePowellTemplate->FixParameter(0, expectedAmplitude);
        orePowellTemplate->FixParameter(1, photonEndpoint);
        orePowellTemplate->SetLineColor(kBlue + 1);
        orePowellTemplate->SetLineWidth(2);
        orePowellTemplate->SetLineStyle(2);
    }

    distributionsPage.cd(1);
    gPad->SetGrid();
    lifetimeHistogram.Draw("HIST");
    inputLifetime.Draw("SAME");
    fittedLifetime.Draw("SAME");
    drawAnalysisBox(analysisBoxes, 0.43, 0.47, 0.94, 0.91, {
        "MC events: N = " + std::to_string(runCount),
        "Model input: #tau = "
            + compactNumber(referenceLifetime*timeScale) + " " + timeUnit,
        "Fit: N_{0} exp(-t/#tau), unbinned MLE",
        "#tau_{fit} = " + compactNumber(displayedFittedLifetime)
            + " #pm " + compactNumber(estimatedError*timeScale) + " " + timeUnit,
        "uncertainty: asymptotic SE; N_{0} derived = "
            + compactNumber(fittedLifetime.GetParameter(0)),
        "red: model input; blue dashed: MC fit",
        "Experimental:",
        "#tau_{exp} = " + compactNumber(experimentalLifetime)
            + " #pm " + compactNumber(experimentalLifetimeError) + " " + timeUnit,
        isPara ? "Al-Ramadhan & Gidley, PRL 72 (1994)"
               : "Vallery et al., PRL 90 (2003)"
    }, 0.021);

    distributionsPage.cd(2);
    gPad->SetGrid();
    photonHistogram.Draw("HIST");
    if (orePowellTemplate) orePowellTemplate->Draw("SAME");
    std::unique_ptr<TLine> simulatedPhotonLine;
    std::unique_ptr<TLine> experimentalPhotonLine;
    if (isPara) {
        const double lineHeight = std::max(1.0, 1.05*photonHistogram.GetMaximum());
        simulatedPhotonLine = std::make_unique<TLine>(
            photonMoments.mean, 0.0, photonMoments.mean, lineHeight);
        simulatedPhotonLine->SetLineColor(kBlue + 1);
        simulatedPhotonLine->SetLineWidth(2);
        simulatedPhotonLine->SetLineStyle(2);
        simulatedPhotonLine->Draw();
        experimentalPhotonLine = std::make_unique<TLine>(
            photonEnergyReference.value, 0.0,
            photonEnergyReference.value, lineHeight);
        experimentalPhotonLine->SetLineColor(kRed + 1);
        experimentalPhotonLine->SetLineWidth(2);
        experimentalPhotonLine->SetLineStyle(3);
        experimentalPhotonLine->Draw();
        drawAnalysisBox(analysisBoxes, 0.42, 0.57, 0.94, 0.91, {
            "MC events: N = " + std::to_string(runCount)
                + "; photons = " + std::to_string(photonEnergies.size()),
            "Fit: monoenergetic #delta(E-#mu)",
            "#mu_{fit} = " + compactNumber(photonMoments.mean, 9) + " keV",
            "#sigma_{truth} = " + compactNumber(photonMoments.sigma) + " keV",
            "blue: MC line; red: accepted-energy line",
            "Experimental / accepted mass:",
            "#mu_{exp} = " + compactNumber(photonEnergyReference.value, 10)
                + " #pm " + compactNumber(
                    photonEnergyReference.totalUncertainty, 2)
                + " keV",
            "PDG 2025 positronium mass"
        }, 0.023);
    } else {
        drawAnalysisBox(analysisBoxes, 0.08, 0.51, 0.60, 0.91, {
            "MC events: N = " + std::to_string(runCount)
                + "; photons = " + std::to_string(photonEnergies.size()),
            "Reference: LO Ore-Powell F(E/E_{max})",
            "A fixed by 3N = "
                + compactNumber(orePowellTemplate->GetParameter(0)),
            "E_{max} fixed = " + compactNumber(photonEndpoint, 9) + " keV",
            "free fit parameters: 0; LO truth-level shape",
            "Experimental / accepted endpoint:",
            "E_{max} = " + compactNumber(photonEnergyReference.value, 10)
                + " #pm " + compactNumber(
                    photonEnergyReference.totalUncertainty, 2)
                + " keV (PDG 2025)",
            "Continuum agrees with QED: Chang et al. (1985)"
        }, 0.021);
    }

    std::unique_ptr<TH1D> angularHistogram;
    std::unique_ptr<TH2D> dalitzHistogram;
    std::unique_ptr<TF1> paraAngularFit;
    LegendreFitSummary anisotropyFit;
    TPaveText channelInformation(0.10, 0.16, 0.90, 0.84, "NDC");
    if (isPara) {
        angularHistogram = std::make_unique<TH1D>("para_polar_histogram",
            "Photon polar angle (unpolarized p-Ps);cos(#theta_{#gamma});Photons",
            histogramBins(paraCosines.size()), -1.0, 1.0);
        styleHistogram(*angularHistogram, kGreen + 2);
        angularHistogram->SetStats(false);
        for (double value : paraCosines) angularHistogram->Fill(value);
        paraAngularFit = std::make_unique<TF1>(
            "para_angular_legendre_fit",
            "[0]*(1+[1]*0.5*(3*x*x-1))", -1.0, 1.0);
        anisotropyFit = fitSecondLegendreAnisotropy(paraIndependentCosines);
        paraAngularFit->SetParameters(
            static_cast<double>(paraCosines.size())/angularHistogram->GetNbinsX(),
            anisotropyFit.anisotropy);
        paraAngularFit->SetParName(0, "normalization");
        paraAngularFit->SetParName(1, "a2");
        paraAngularFit->SetLineColor(kBlue + 1);
        paraAngularFit->SetLineWidth(2);
        paraAngularFit->SetLineStyle(2);
        distributionsPage.cd(3);
        gPad->SetGrid();
        angularHistogram->Draw("HIST");
        paraAngularFit->Draw("SAME");
        std::vector<std::string> angularAnalysis{
            "MC events: N = " + std::to_string(runCount)
                + "; photons = " + std::to_string(paraCosines.size()),
            "Fit: C[1+a_{2}P_{2}(cos#theta)]",
            "constrained unbinned event-level MLE",
            "C = " + compactNumber(paraAngularFit->GetParameter(0)),
            "a_{2} = " + compactNumber(anisotropyFit.anisotropy),
            anisotropyFit.atBoundary
                ? "boundary solution; uncertainty unavailable"
                : "SE(a_{2}) = " + compactNumber(anisotropyFit.standardError),
            "Experimental:",
            "no apparatus-independent a_{2}; ideal unpolarized a_{2}="
                + compactNumber(anisotropyReference.value)
        };
        drawAnalysisBox(analysisBoxes, 0.43, 0.55, 0.94, 0.91,
                        angularAnalysis, 0.022);
        distributionsPage.cd(4);
        channelInformation.SetFillColorAlpha(kWhite, 0.92);
        channelInformation.SetTextAlign(12);
        channelInformation.SetTextFont(42);
        channelInformation.AddText("Ideal vacuum, truth level");
        channelInformation.AddText(
            ("MC events: N = " + std::to_string(runCount)).c_str());
        channelInformation.AddText("Dominant channel: p-Ps #rightarrow 2#gamma");
        channelInformation.AddText("Photon multiplicity: 2");
        channelInformation.AddText("Opening angle in CM: 180 degrees");
        std::ostringstream energyText;
        energyText << "E_{#gamma} = " << photonEndpoint << " keV";
        channelInformation.AddText(energyText.str().c_str());
        channelInformation.AddText("No detector broadening or material effects");
        channelInformation.Draw();
    } else {
        const int dalitzBins = std::clamp(static_cast<int>(std::lround(
            2.0*std::sqrt(std::sqrt(static_cast<double>(events.size()))))),
            6, 20);
        dalitzHistogram = std::make_unique<TH2D>("ortho_dalitz_histogram",
            "Three-photon Dalitz distribution;E_{max} [keV];E_{mid} [keV]",
            dalitzBins, (2.0/3.0)*photonEndpoint, photonEndpoint,
            dalitzBins, 0.5*photonEndpoint, photonEndpoint);
        dalitzHistogram->SetDirectory(nullptr);
        dalitzHistogram->SetStats(false);
        for (size_t index = 0; index < orthoMaximumEnergies.size(); ++index) {
            dalitzHistogram->Fill(orthoMaximumEnergies[index], orthoMiddleEnergies[index]);
        }
        angularHistogram = std::make_unique<TH1D>("ortho_leading_angle_histogram",
            "Angle between two leading photons;#theta_{12} [deg];Events",
            histogramBins(orthoLeadingAngles.size()), 120.0, 180.0);
        styleHistogram(*angularHistogram, kMagenta + 1);
        angularHistogram->SetStats(false);
        for (double value : orthoLeadingAngles) angularHistogram->Fill(value);
        distributionsPage.cd(3);
        gPad->SetRightMargin(0.14);
        dalitzHistogram->Draw("COLZ");
        drawAnalysisBox(analysisBoxes, 0.08, 0.64, 0.56, 0.91, {
            "MC events: N = " + std::to_string(runCount),
            "Generator/reference: Ore-Powell 2D template",
            "P=A #Sigma_{ij}(1-cos#theta_{ij})^{2}",
            "free shape parameters: 0",
            "Experimental:",
            "no matching acceptance-corrected dataset loaded"
        }, 0.023);
        distributionsPage.cd(4);
        gPad->SetGrid();
        angularHistogram->Draw("HIST");
        drawAnalysisBox(analysisBoxes, 0.08, 0.61, 0.58, 0.91, {
            "MC events: N = " + std::to_string(runCount),
            "Generator/reference: Ore-Powell angular projection",
            "free shape parameters: 0",
            "range: 120 deg #leq #theta_{12} #leq 180 deg",
            "Experimental:",
            "no matching detector-independent parameter loaded"
        }, 0.024);
    }

    TH1D energyClosureHistogram("decay_energy_closure",
        "Energy-conservation closure;(#SigmaE_{#gamma}-E_{Ps})/(E_{Ps}#epsilon);Events",
        129, -64.5, 64.5);
    TH1D momentumClosureHistogram("decay_momentum_closure",
        "Momentum-conservation closure;c|#Sigmap_{#gamma}|/(E_{Ps}#epsilon);Events",
        64, 0.0, 64.0);
    TH1D photonShellHistogram("decay_photon_mass_shell",
        "Photon mass-shell closure;(E_{#gamma}^{2}-c^{2}p_{#gamma}^{2})/(E_{#gamma}^{2}#epsilon);Photons",
        129, -64.5, 64.5);
    TH1D parentInvariantHistogram("decay_parent_invariant",
        "Parent invariant-mass-squared closure;(E_{sum}^{2}-c^{2}|p_{sum}|^{2}-E_{Ps}^{2})/(E_{Ps}^{2}#epsilon);Events",
        129, -64.5, 64.5);
    styleHistogram(energyClosureHistogram, kOrange + 7);
    styleHistogram(momentumClosureHistogram, kAzure + 1);
    styleHistogram(photonShellHistogram, kGreen + 2);
    styleHistogram(parentInvariantHistogram, kMagenta + 1);
    for (TH1D* histogram : {&energyClosureHistogram, &momentumClosureHistogram,
                            &photonShellHistogram, &parentInvariantHistogram}) {
        histogram->SetStats(false);
    }
    for (double value : energyClosureEpsilon) energyClosureHistogram.Fill(value);
    for (double value : momentumClosureEpsilon) momentumClosureHistogram.Fill(value);
    for (double value : photonMassShellEpsilon) photonShellHistogram.Fill(value);
    for (double value : parentInvariantEpsilon) parentInvariantHistogram.Fill(value);

    const auto maximumAbsolute = [](const std::vector<double>& values) {
        double maximum = 0.0;
        for (double value : values) maximum = std::max(maximum, std::abs(value));
        return maximum;
    };
    const double maximumClosure = std::max({
        maximumAbsolute(energyClosureEpsilon),
        maximumAbsolute(momentumClosureEpsilon),
        maximumAbsolute(photonMassShellEpsilon),
        maximumAbsolute(parentInvariantEpsilon)});
    const auto drawClosureAnalysis = [&](TH1D& histogram,
                                         const std::vector<double>& values,
                                         const std::string& functionName,
                                         bool includeGlobalNotes) {
        const GaussianFitSummary fit = gaussianMaximumLikelihood(values);
        std::unique_ptr<TF1> curve = gaussianMleOverlay(
            functionName, fit, histogram, kBlue + 1);
        if (curve) analysisFunctions.push_back(std::move(curve));
        std::vector<std::string> lines{
            "Entries: N = " + std::to_string(fit.count),
            fit.sigma > 0.0 ? "Fit: Gaussian MLE in plotted residual"
                            : "Fit: degenerate Gaussian (#sigma=0)",
            "#mu = " + compactNumber(fit.mean),
            "#sigma = " + compactNumber(fit.sigma),
            "Experimental:",
            "n/a - numerical closure, expected center = 0"
        };
        if (includeGlobalNotes) {
            lines.push_back("max all closures = "
                + compactNumber(maximumClosure) + " #epsilon");
            lines.push_back("J not tested: photon helicities absent");
        }
        drawAnalysisBox(analysisBoxes, includeGlobalNotes ? 0.58 : 0.52,
                        includeGlobalNotes ? 0.49 : 0.61,
                        0.97, 0.91, lines, includeGlobalNotes ? 0.020 : 0.024);
    };

    diagnosticsPage.cd(1);
    gPad->SetGrid();
    energyClosureHistogram.Draw("HIST");
    drawClosureAnalysis(energyClosureHistogram, energyClosureEpsilon,
                        "decay_energy_closure_gaussian", false);
    diagnosticsPage.cd(2);
    gPad->SetGrid();
    momentumClosureHistogram.Draw("HIST");
    drawClosureAnalysis(momentumClosureHistogram, momentumClosureEpsilon,
                        "decay_momentum_closure_gaussian", false);
    diagnosticsPage.cd(3);
    gPad->SetGrid();
    photonShellHistogram.Draw("HIST");
    drawClosureAnalysis(photonShellHistogram, photonMassShellEpsilon,
                        "decay_photon_shell_gaussian", false);
    diagnosticsPage.cd(4);
    gPad->SetGrid();
    parentInvariantHistogram.Draw("HIST");
    drawClosureAnalysis(parentInvariantHistogram, parentInvariantEpsilon,
                        "decay_parent_invariant_gaussian", true);

    canvas.cd();
    distributionsPage.Pop();
    canvas.Modified();
    canvas.Update();
    std::vector<root_export::NamedPad> plotsToSave{
        {distributionsPage.GetPad(1), 1, 1, "annihilation_time"},
        {distributionsPage.GetPad(2), 1, 2, "photon_energy"},
        {distributionsPage.GetPad(3), 1, 3, isPara ? "photon_polar_angle"
                                              : "three_photon_dalitz"},
        {diagnosticsPage.GetPad(1), 2, 1, "diagnostic_energy_closure"},
        {diagnosticsPage.GetPad(2), 2, 2, "diagnostic_momentum_closure"},
        {diagnosticsPage.GetPad(3), 2, 3, "diagnostic_photon_mass_shell"},
        {diagnosticsPage.GetPad(4), 2, 4, "diagnostic_parent_invariant_mass"}
    };
    if (!isPara) {
        plotsToSave.push_back(
            {distributionsPage.GetPad(4), 1, 4, "leading_photon_angle"});
    }
    reportExports(root_export::saveStatisticalPlots(
        selectedPhenomenon, plotsToSave));

    bool persistenceOk = reportArchiveOperation(
        statistics_archive::writeScientificReferencesText(),
        "scientific-reference catalogue");
    return persistenceOk ? 0 : 3;
}

enum class BeamOutcome { Escaped, ShortRange, Unresolved, NumericalFailure };

struct BeamConfiguration {
    double centreOfMassKineticEnergy;
    double impactParameterMaximum;
    double matchingRadius;
    double maximumFlightTime;
    double asymptoticRelativeSpeed;
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
            frame.time, frame.radius};
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
    const double coulombStrength = coulomb * eCharge * eCharge;
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

    const double gammaInfinity = 1.0 + energy / (2.0 * electronMass * c*c);
    const double particleSpeed = c * std::sqrt(1.0 - 1.0 / (gammaInfinity * gammaInfinity));
    const double relativeSpeed = 2.0 * particleSpeed;
    return {energy, impactMaximum, matchingRadius,
            8.0 * matchingRadius / relativeSpeed, relativeSpeed,
            analysisThetaMinimum, coulombLength, cutoffImpactParameter,
            angleBins, shortRangeFocus};
}

BeamEvent simulateBeamEvent(std::uint64_t seed, const BeamConfiguration& configuration) {
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

    const double coulombStrength = coulomb * eCharge * eCharge;
    const double finiteKineticEnergy = configuration.centreOfMassKineticEnergy
        + coulombStrength / configuration.matchingRadius;
    const double finiteGamma = 1.0
        + finiteKineticEnergy / (2.0 * electronMass * c*c);
    const double finiteParticleSpeed = c
        * std::sqrt(1.0 - 1.0 / (finiteGamma * finiteGamma));
    const double finiteRelativeSpeed = 2.0 * finiteParticleSpeed;
    const double asymptoticGamma = 1.0 + configuration.centreOfMassKineticEnergy
        / (2.0 * electronMass*c*c);
    const double asymptoticParticleSpeed = 0.5 * configuration.asymptoticRelativeSpeed;
    const double asymptoticMomentum = asymptoticGamma * electronMass
        * asymptoticParticleSpeed;
    const double finiteMomentum = finiteGamma * electronMass * finiteParticleSpeed;
    const double tangentialFraction = impactParameter * asymptoticMomentum
        / (configuration.matchingRadius * finiteMomentum);
    if (!(tangentialFraction >= 0.0 && tangentialFraction < 1.0)) {
        BeamEvent result{BeamOutcome::NumericalFailure, seed, impactParameter};
        result.impactAzimuth = azimuth;
        return result;
    }
    const double tangentialSpeed = finiteRelativeSpeed * tangentialFraction;
    const double radialSpeedSquared = finiteRelativeSpeed * finiteRelativeSpeed
                                    - tangentialSpeed * tangentialSpeed;
    if (!(radialSpeedSquared > 0.0)) {
        BeamEvent result{BeamOutcome::NumericalFailure, seed, impactParameter};
        result.impactAzimuth = azimuth;
        return result;
    }
    const Vec3 relativeVelocity = radialDirection * (-std::sqrt(radialSpeedSquared))
                                + tangentDirection * tangentialSpeed;

    const auto randomDirection = [&]() {
        const double cosine = 2.0 * uniform(random) - 1.0;
        const double phi = 2.0 * pi * uniform(random);
        const double transverse = std::sqrt(std::max(0.0, 1.0 - cosine*cosine));
        return Vec3{transverse * std::cos(phi), transverse * std::sin(phi), cosine};
    };
    State state;
    state.electronPosition = relativePosition * 0.5;
    state.positronPosition = relativePosition * -0.5;
    state.electronVelocity = relativeVelocity * 0.5;
    state.positronVelocity = relativeVelocity * -0.5;
    state.electronDipole = randomDirection() * bohrMagneton;
    state.positronDipole = randomDirection() * bohrMagneton;
    ClassicalTrajectoryEngine trajectory(state);
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

    constexpr double reducedMass = electronMass * positronMass
                                 / (electronMass + positronMass);
    bool passedClosestApproach = false;
    while (state.time < configuration.maximumFlightTime) {
        const double radius = separation(state);
        const double relativeSpeed =
            (state.electronVelocity - state.positronVelocity).norm();
        const double omega = std::sqrt(coulomb * eCharge*eCharge
                                      / (reducedMass * radius*radius*radius));
        const double transitStep = 0.02 * radius/std::max(relativeSpeed, 1.0);
        const double dt = std::min({2.0e-18, 2.0*pi/(320.0*omega), transitStep,
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
            state.electronPosition - state.positronPosition;
        const Vec3 currentRelativeVelocity =
            state.electronVelocity - state.positronVelocity;
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
        if (passedClosestApproach && currentRadius >= configuration.matchingRadius) {
            const double angle = std::acos(std::clamp(
                dot(beamDirection, currentRelativeVelocity)
                    / currentRelativeVelocity.norm(), -1.0, 1.0));
            // At the finite matching sphere the kinetic energy still contains
            // the Coulomb acceleration.  The conservative mechanical energy
            // is the asymptotic outgoing kinetic energy estimator.
            const double outgoingEnergy = makeFrame(state).mechanicalEnergy;
            const Vec3 previousRelativePosition =
                beforeStep.electronPosition - beforeStep.positronPosition;
            const Vec3 stepDisplacement =
                currentRelativePosition - previousRelativePosition;
            const double quadraticA = dot(stepDisplacement, stepDisplacement);
            const double quadraticB = 2.0 * dot(previousRelativePosition,
                                                 stepDisplacement);
            const double quadraticC = dot(previousRelativePosition,
                                           previousRelativePosition)
                                    - configuration.matchingRadius
                                      * configuration.matchingRadius;
            double crossingFraction = 1.0;
            const double discriminant = quadraticB*quadraticB
                                      - 4.0*quadraticA*quadraticC;
            if (quadraticA > 0.0 && discriminant >= 0.0) {
                crossingFraction = std::clamp(
                    (-quadraticB + std::sqrt(discriminant))/(2.0*quadraticA),
                    0.0, 1.0);
            }
            const State diagnosticEndpoint = interpolateState(
                beforeStep, state, crossingFraction);
            BeamEvent result = makeResult(BeamOutcome::Escaped, angle,
                outgoingEnergy, state.radiatedEnergy, &diagnosticEndpoint);
            result.finalDiagnostics = endpointDiagnostics(diagnosticEndpoint);
            result.diagnosticsValid = true;
            return result;
        }
    }
    return makeResult(BeamOutcome::Unresolved,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(), state.radiatedEnergy, &state);
}

std::vector<BeamEvent> runBeamExperiment(std::uint64_t masterSeed,
                                         const BeamConfiguration& configuration,
                                         int runCount) {
    std::vector<BeamEvent> events(static_cast<size_t>(runCount));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::mutex outputMutex;
    const int workerCount = std::min(runCount,
        static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    const auto worker = [&]() {
        while (true) {
            const int index = nextIndex.fetch_add(1);
            if (index >= runCount) break;
            events[static_cast<size_t>(index)] = simulateBeamEvent(
                splitMix64(masterSeed + static_cast<std::uint64_t>(index)), configuration);
            const int done = completed.fetch_add(1) + 1;
            if (done % 10 == 0 || done == runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout << "Beam trajectories: " << done << "/" << runCount << '\n';
            }
        }
    };
    for (int index = 0; index < workerCount; ++index) workers.emplace_back(worker);
    for (std::thread& workerThread : workers) workerThread.join();
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
    const double matchingPotentialRatio = coulomb*eCharge*eCharge
        / (configuration.matchingRadius * configuration.centreOfMassKineticEnergy);
    const bool matchingRegionWarning = configuration.matchingRadius
            < 20.0*configuration.impactParameterMaximum
        || matchingPotentialRatio > 0.01;
    std::cout << "Running an e+e- beam experiment with " << runCount
              << (runCount == 1 ? " trajectory on " : " trajectories on ")
              << workerCount << " worker"
              << (workerCount == 1 ? "" : "s") << ".\n"
              << "K_CM = " << beamEnergyEv << " eV, b_max = "
              << configuration.impactParameterMaximum * 1.0e12
              << " pm, R_match = " << configuration.matchingRadius * 1.0e12
              << " pm, theta acceptance >= "
              << configuration.analysisThetaMinimum * 180.0/pi << " deg, beta_infinity = "
              << configuration.asymptoticRelativeSpeed/(2.0*c) << ".\n";
    if (configuration.asymptoticRelativeSpeed/(2.0*c) > 0.1) {
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
            case BeamOutcome::Unresolved: ++unresolved; break;
            case BeamOutcome::NumericalFailure: ++failed; break;
        }
    }
    if (failed > 0) {
        std::cerr << "Outcomes: escaped=" << escaped << ", short-range=" << shortRange
                  << ", unresolved=" << unresolved << ", failed=" << failed << '\n'
                  << "Cross-section report is INVALID because at least one trajectory "
                     "became numerically non-finite. No plots were produced.\n";
        bool persistenceOk = reportArchiveOperation(
            statistics_archive::writeScientificReferencesText(),
            "scientific-reference catalogue");
        return persistenceOk ? 2 : 3;
    }
    constexpr double diagnosticFloor = 64.0 * std::numeric_limits<double>::epsilon();
    const double energyFloor = diagnosticFloor * 2.0*electronMass*c*c;
    const double momentumFloor = diagnosticFloor * 2.0*electronMass*c;
    std::vector<double> relativeEnergyClosures;
    std::vector<double> absoluteEnergyClosures;
    std::vector<double> relativeMomentumClosures;
    std::vector<double> relativeAngularMomentumClosures;
    std::vector<double> logMomentumClosures;
    std::vector<double> logAngularMomentumClosures;
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
    std::cout << "Elastic data are a classical low-velocity model, not Bhabha scattering.\n";
    if (unresolved > 0) {
        std::cout << "Warning: " << unresolved
                  << " trajectories did not leave the interaction region before the "
                     "flight-time gate; reported channel cross sections are finite-gate values.\n";
    }
    if (runCount < 1000) {
        std::cout << "Statistical note: " << runCount
                  << (runCount == 1 ? " trajectory is" : " trajectories are")
                  << " suitable for a preview; use at least 1000 "
                     "for a smoother differential cross section.\n";
    }
    if (!relativeEnergyClosures.empty()) {
        std::cout << "Numerical closure on " << relativeEnergyClosures.size()
                  << " escaped trajectories at R_match: median |deltaE|="
                  << sampleQuantile(absoluteEnergyClosures, 0.5)
                  << ", median delta(P_N+P_rad)="
                  << sampleQuantile(relativeMomentumClosures, 0.5)
                  << ", median delta(J_N+J_rad)="
                  << sampleQuantile(relativeAngularMomentumClosures, 0.5)
                  << ".\n";
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
    angularCrossSection.SetLineColor(kAzure + 1);
    angularCrossSection.SetMarkerColor(kAzure + 1);
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
    rutherfordHistogram.SetLineColor(kRed + 1);
    rutherfordHistogram.SetLineWidth(2);
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
    fittedRutherfordHistogram.SetLineColor(kBlue + 2);
    fittedRutherfordHistogram.SetLineStyle(2);
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
    cumulativeGraph.SetLineColor(kGreen + 3);
    cumulativeGraph.SetMarkerColor(kGreen + 3);
    cumulativeGraph.SetMarkerStyle(20);
    cumulativeGraph.SetMarkerSize(0.55);
    TGraph cumulativeRutherford(static_cast<int>(cumulativePointCount),
        cumulativeAngles.data(), cumulativeReference.data());
    cumulativeRutherford.SetLineColor(kRed + 1);
    cumulativeRutherford.SetLineWidth(2);
    TGraph cumulativeFit(static_cast<int>(cumulativePointCount),
        cumulativeAngles.data(), cumulativeFitted.data());
    cumulativeFit.SetLineColor(kBlue + 2);
    cumulativeFit.SetLineStyle(2);
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
        energyLossSpectrum->SetLineColor(kMagenta + 1);
        energyLossSpectrum->SetFillColorAlpha(kMagenta + 1, 0.45);
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
              << ", #beta_{#infty} = "
              << configuration.asymptoticRelativeSpeed/(2.0*c);
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
    outcomeLine << "escaped/cutoff/time-limit/failed = " << escaped << '/'
                << shortRange << '/' << unresolved << '/' << failed;
    summary.AddText(outcomeLine.str().c_str());
    summary.AddText("Classical low-v trajectory model");
    if (configuration.shortRangeFocus) {
        summary.AddText("#sigma_{cutoff} is not an annihilation cross section");
    }
    if (unresolved > 0) {
        summary.AddText("Finite flight-time gate: channel values are incomplete");
    }
    if (configuration.asymptoticRelativeSpeed/(2.0*c) > 0.1) {
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
        "Energy-balance closure;#Delta(E_{N}+E_{rad}+E_{S})/E_{scale};Escaped trajectories",
        histogramBins(relativeEnergyClosures.size()),
        -energyClosureLimit, energyClosureLimit);
    TH1D momentumBalanceHistogram("beam_momentum_balance_closure",
        "Particle-field momentum closure;log_{10}(|#Delta(P_{N}+P_{rad})|/P_{scale});Escaped trajectories",
        histogramBins(logMomentumClosures.size()), -14.5, momentumLogUpper);
    TH1D angularMomentumBalanceHistogram("beam_angular_momentum_balance_closure",
        "Particle-field angular momentum closure;log_{10}(|#Delta(J_{N}+J_{rad})|/J_{scale});Escaped trajectories",
        histogramBins(logAngularMomentumClosures.size()), -14.5,
        angularMomentumLogUpper);
    styleHistogram(energyBalanceHistogram, kOrange + 7);
    styleHistogram(momentumBalanceHistogram, kAzure + 1);
    styleHistogram(angularMomentumBalanceHistogram, kMagenta + 1);
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
    } else {
        diagnosticSummary.AddText("No escaped endpoint pair: closure histograms are empty");
    }
    diagnosticSummary.AddText("Cutoff and time-gate outcomes are excluded");
    diagnosticSummary.AddText("P_{rad}, J_{rad}: integrated Maxwell flux on the control wavefront");
    diagnosticSummary.AddText("P_{N}, J_{N}: approximate near-field/particle Noether sector");
    if (configuration.asymptoticRelativeSpeed/(2.0*c) > 0.1) {
        diagnosticSummary.AddText("WARNING: outside the intended low-velocity regime");
    }
    if (matchingRegionWarning) {
        diagnosticSummary.AddText("WARNING: matching sphere may be too close");
    }

    diagnosticsPage.cd(1);
    gPad->SetGrid();
    energyBalanceHistogram.Draw("HIST");
    const auto drawBeamResidualFit = [&](TH1D& histogram,
                                         const std::vector<double>& values,
                                         const std::string& functionName) {
        const GaussianFitSummary fit = gaussianMaximumLikelihood(values);
        std::unique_ptr<TF1> curve = gaussianMleOverlay(
            functionName, fit, histogram, kBlue + 1);
        if (curve) beamAnalysisFunctions.push_back(std::move(curve));
        std::vector<std::string> lines;
        if (fit.count == 0) {
            lines = {"Entries: N = 0",
                     "Fit: unavailable - no escaped endpoint pair",
                     "Experimental:",
                     "n/a - numerical conservation diagnostic"};
        } else {
            lines = {
                "Entries: N = " + std::to_string(fit.count),
                fit.sigma > 0.0 ? "Fit: Gaussian MLE in plotted variable"
                                : "Fit: degenerate Gaussian (#sigma=0)",
                "#mu = " + compactNumber(fit.mean),
                "#sigma = " + compactNumber(fit.sigma),
                "Experimental:",
                "n/a - numerical conservation diagnostic"
            };
        }
        drawAnalysisBox(beamAnalysisBoxes, 0.51, 0.62, 0.95, 0.91,
                        lines, 0.023);
    };
    drawBeamResidualFit(energyBalanceHistogram, relativeEnergyClosures,
                        "beam_energy_balance_gaussian");
    diagnosticsPage.cd(2);
    gPad->SetGrid();
    momentumBalanceHistogram.Draw("HIST");
    drawBeamResidualFit(momentumBalanceHistogram, logMomentumClosures,
                        "beam_momentum_balance_gaussian");
    diagnosticsPage.cd(3);
    gPad->SetGrid();
    angularMomentumBalanceHistogram.Draw("HIST");
    drawBeamResidualFit(angularMomentumBalanceHistogram,
                        logAngularMomentumClosures,
                        "beam_angular_momentum_balance_gaussian");
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

int showStatisticalAnalysis(std::uint64_t seed, int selectedPhenomenon,
                            int runCount,
                            double beamEnergyEv,
                            double thetaMinimumDegrees, int angleBins,
                            double impactMaximumPm, double matchingRadiusPm) {
    if (selectedPhenomenon == 1 || selectedPhenomenon == 2) {
        return showBoundDecayStatistics(seed, selectedPhenomenon, runCount);
    }
    return showBeamStatistics(seed, selectedPhenomenon, runCount,
                              beamEnergyEv, thetaMinimumDegrees,
                              angleBins, impactMaximumPm, matchingRadiusPm);
}
#endif

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
#include "tests/maxwell_validation.hpp"
#endif
} // namespace

int main(int argc, char** argv) {
#ifdef POSITRONIUM_VALIDATION_EXECUTABLE
    (void)argc;
    (void)argv;
    return runMaxwellSelfTest();
#else
    std::random_device seedSource;
    std::uint64_t seed = (static_cast<std::uint64_t>(seedSource()) << 32) ^ seedSource();
    bool diagnose = false;
    int selectedMode = 0;
    VisualStyle visualStyle = VisualStyle::Unselected;
    int selectedPhenomenon = 0;
    int statisticalRuns = 10000;
    double beamEnergyEv = 20.0;
    double thetaMinimumDegrees = 5.0;
    int angleBins = 10;
    double impactParameterMaximumPm = 0.0;
    double matchingRadiusPm = 0.0;
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
            } else if (argument == "--seed") {
                seed = std::stoull(requireValue(argument));
            } else if (argument == "--phenomenon") {
                selectedPhenomenon = std::stoi(requireValue(argument));
            } else if (argument == "--runs") {
                statisticalRuns = std::stoi(requireValue(argument));
            } else if (argument == "--stat-window-ps") {
                (void)requireValue(argument);
                throw std::invalid_argument(
                    "--stat-window-ps was removed: bound-state statistics now sample "
                    "a phenomenological ideal-vacuum annihilation-time distribution");
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
    if (diagnose) selectedMode = 1;
    if (selectedMode == 0) {
        std::cout << "Choose simulation mode:\n"
                  << "1 -> Visual simulation\n"
                  << "2 -> Statistical analysis (" << statisticalRuns
                  << " event" << (statisticalRuns == 1 ? "" : "s/trials") << ")\n"
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
    if (selectedPhenomenon < 1 || selectedPhenomenon > 4) {
        if (selectedMode == 2) {
            std::cout << "Choose statistical experiment:\n"
                      << "1 -> Para-positronium decay in vacuum (2 gamma)\n"
                      << "2 -> Ortho-positronium decay in vacuum (3 gamma)\n"
                      << "3 -> e+e- beam: short-range/cutoff channel\n"
                      << "4 -> e+e- beam: elastic scattering\n";
        } else {
            std::cout << "Choose phenomenon to simulate:\n"
                      << "1 -> Para-positronium\n"
                      << "2 -> Ortho-positronium\n"
                      << "3 -> Direct collision\n"
                      << "4 -> Scattering\n";
        }
        std::cout << "Selection [1-4]: " << std::flush;
        if (!(std::cin >> selectedPhenomenon) || selectedPhenomenon < 1 || selectedPhenomenon > 4) {
            std::cerr << "Invalid selection. Enter a number from 1 to 4.\n";
            return 1;
        }
    }
    if (selectedMode == 2) {
        if (statisticalRuns < 1 || statisticalRuns > 100000) {
            std::cerr << "The number of statistical events/trials must be from 1 to 100000.\n";
            return 1;
        }
        try {
            return showStatisticalAnalysis(seed, selectedPhenomenon,
                                           statisticalRuns, beamEnergyEv,
                                           thetaMinimumDegrees,
                                           angleBins, impactParameterMaximumPm,
                                           matchingRadiusPm);
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
            / std::max(frames.front().canonicalMomentumScale, electronMass * c * 1.0e-15);
        const double angularMomentumScale =
            std::max(hbar, frames.front().noetherAngularMomentum.norm());
        const double relativeAngularMomentumDrift =
            angularMomentumDrift.norm() / angularMomentumScale;
        const bool expectedMotion = initialConditions.phenomenon == Phenomenon::Scattering
            ? frames.back().radius > frames.front().radius
            : frames.back().radius < frames.front().radius;
        double maximumRelativeDipoleNormDrift = 0.0;
        for (const Frame& frame : frames) {
            maximumRelativeDipoleNormDrift = std::max({
                maximumRelativeDipoleNormDrift,
                std::abs(frame.electronDipole.norm() / bohrMagneton - 1.0),
                std::abs(frame.positronDipole.norm() / bohrMagneton - 1.0)});
        }
        const bool dipoleNormsValid = std::isfinite(maximumRelativeDipoleNormDrift)
                                  && maximumRelativeDipoleNormDrift < 1.0e-12;
        // At a direct collision the point-particle model terminates in its
        // singular regime, so only pre-contact motion is a meaningful test.
        const bool trajectoryValid = std::isfinite(frames.back().radius)
                                  && expectedMotion && dipoleNormsValid;
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
                  << " eV\n"
                  << "LL/coherent dE: " << frames.back().reactionEnergyMismatch/eCharge
                  << " eV\n"
                  << "energy drift:   " << energyDrift / eCharge << " eV ("
                  << relativeEnergyDrift * 100.0 << "%)\n"
                  << "energy balance: particles + bound/interference field"
                     " + outward radiation\n"
                  << "field |P_rad|:  " << frames.back().radiatedMomentum.norm()
                  << " kg m/s\n"
                  << "field |J_rad|:  " << frames.back().radiatedAngularMomentum.norm()/hbar
                  << " hbar\n"
                  << "total |dP|:     " << momentumDrift.norm() << " kg m/s ("
                  << relativeMomentumDrift * 100.0 << "% of characteristic p)\n"
                  << "total |dJ|:     " << angularMomentumDrift.norm() / hbar << " hbar ("
                  << relativeAngularMomentumDrift * 100.0
                  << "% of characteristic J)\n"
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

    // Keep an owned snapshot: `simulation` is replaced by the complete run
    // after the canvas has been created.
    const Frame initialFrame = frames.front();
    const double initialElectronX = initialFrame.electron.x * scale;
    const double initialElectronY = initialFrame.electron.y * scale;
    const double initialElectronZ = initialFrame.electron.z * scale;
    const double initialPositronX = initialFrame.positron.x * scale;
    const double initialPositronY = initialFrame.positron.y * scale;
    const double initialPositronZ = initialFrame.positron.z * scale;
    TPolyLine3D electronPath(1, &initialElectronX,
                             &initialElectronY, &initialElectronZ);
    electronPath.SetLineColor(kBlue + 2);
    electronPath.SetLineWidth(2);
    TPolyLine3D positronPath(1, &initialPositronX,
                             &initialPositronY, &initialPositronZ);
    positronPath.SetLineColor(kOrange + 7);
    positronPath.SetLineWidth(2);

    TPolyMarker3D electronDots(1, 7);
    electronDots.SetMarkerColor(kBlue + 2);
    electronDots.SetMarkerSize(0.7);
    TPolyMarker3D positronDots(1, 7);
    positronDots.SetMarkerColor(kOrange + 7);
    positronDots.SetMarkerSize(0.7);
    electronDots.SetPoint(0, initialElectronX, initialElectronY, initialElectronZ);
    positronDots.SetPoint(0, initialPositronX, initialPositronY, initialPositronZ);

    if (visualStyle == VisualStyle::Line) {
        electronPath.Draw();
        positronPath.Draw("same");
    } else {
        electronDots.Draw();
        positronDots.Draw("same");
    }

    TPolyMarker3D electron(1), positron(1);
    electron.SetMarkerStyle(20); electron.SetMarkerSize(2.2); electron.SetMarkerColor(kAzure + 1);
    positron.SetMarkerStyle(20); positron.SetMarkerSize(2.2); positron.SetMarkerColor(kRed + 1);
    if (visualStyle == VisualStyle::Line) {
        electron.Draw("same");
        positron.Draw("same");
    }

    TPolyLine3D electronDipoleShaft(2), electronDipoleLeft(2), electronDipoleRight(2);
    TPolyLine3D positronDipoleShaft(2), positronDipoleLeft(2), positronDipoleRight(2);
    for (TPolyLine3D* arrow : {&electronDipoleShaft, &electronDipoleLeft, &electronDipoleRight}) {
        arrow->SetLineColor(kAzure + 1); arrow->SetLineWidth(3);
        if (visualStyle == VisualStyle::Line) arrow->Draw("same");
    }
    for (TPolyLine3D* arrow : {&positronDipoleShaft, &positronDipoleLeft, &positronDipoleRight}) {
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
    visualOptions.frameCount = 480;
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
        const double electronX = f.electron.x * scale;
        const double electronY = f.electron.y * scale;
        const double electronZ = f.electron.z * scale;
        const double positronX = f.positron.x * scale;
        const double positronY = f.positron.y * scale;
        const double positronZ = f.positron.z * scale;
        // The canvas is seeded from a short deterministic probe so it can be
        // shown immediately.  Replace that seed with the full run's actual
        // first frame instead of retaining it as a detached history point.
        const std::size_t pointIndex = receivedLiveInitialFrame ? renderedFrame : 0;
        electronPath.SetPoint(static_cast<int>(pointIndex),
                              electronX, electronY, electronZ);
        positronPath.SetPoint(static_cast<int>(pointIndex),
                              positronX, positronY, positronZ);
        electronDots.SetPoint(static_cast<int>(pointIndex),
                              electronX, electronY, electronZ);
        positronDots.SetPoint(static_cast<int>(pointIndex),
                              positronX, positronY, positronZ);
        if (receivedLiveInitialFrame) {
            ++renderedFrame;
        } else {
            receivedLiveInitialFrame = true;
        }
        if (visualStyle == VisualStyle::Line) {
            electron.SetPoint(0, electronX, electronY, electronZ);
            positron.SetPoint(0, positronX, positronY, positronZ);
            const Vec3 electronPosition = f.electron * scale;
            const Vec3 positronPosition = f.positron * scale;
            setDipoleArrow(electronDipoleShaft, electronDipoleLeft, electronDipoleRight,
                           electronPosition, f.electronDipole);
            setDipoleArrow(positronDipoleShaft, positronDipoleLeft, positronDipoleRight,
                           positronPosition, f.positronDipole);
        }
        const double requiredSpan = 1.10 * std::max(
            f.electron.norm(), f.positron.norm()) * scale;
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
        renderVisualFrame(frame);
        lastRepaint = VisualClock::now();
        gSystem->Sleep(10);
    };
    simulation = simulate(seed, selectedPhenomenon, visualOptions);
    if (gExitRequested) return 0;
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
    const root_export::ExportResult screenshot =
        root_export::saveVisualScreenshot(canvas, selectedPhenomenon);
    if (screenshot) {
        std::cout << "Saved visual screenshot: "
                  << screenshot.path.string() << '\n';
    } else {
        std::cerr << "Warning: could not save visual screenshot "
                  << screenshot.path.string() << ": "
                  << screenshot.error << '\n';
    }
    app.Run();
    return 0;
#endif
}
