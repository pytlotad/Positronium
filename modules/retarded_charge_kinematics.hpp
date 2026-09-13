#pragma once

// Retarded point-charge kinematics: reconstructing where/how fast/how
// accelerated a charge WAS at an earlier time from the state history, and
// solving the light-cone equation that turns "earlier" into "retarded" for
// the mutual Lienard-Wiechert field.  This is the smallest self-contained
// slice of the engine that lienardWiechertField and its two Newton-loop
// siblings in electrodynamics.hpp (retardedElectricDipoleField,
// retardedMagneticDipoleField) all build on.
//
// Self-contained and order-independent.  It names what it needs through a
// using-directive on positronium::parameters and using-declarations for the
// object types, rather than reopening namespace positronium: the header is
// still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.

#include "pair_geometry.hpp"
#include "physical_constants.hpp"
#include "relativistic_field_types.hpp"
#include "state.hpp"
#include "vector3.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

// CREM_CAUSALITY audit counters.  A retarded field is causal only if every
// source sample it reads lies at or before the present, and only if the root
// the light-cone solve lands on is the RETARDED one (t_ret <= t_obs) rather
// than the advanced one, which the same quadratic also admits.  These count
// violations rather than asserting, so one run reports how often and how badly
// instead of stopping at the first.
//
// Atomic, and every test computed from local values.  The first version of
// this audit used a plain bool to mark "inside the converged read" and plain
// counters, which the collapse experiment's worker threads raced: it reported
// hundreds of acausal converged reads that a single-threaded run of the same
// seeds showed to be zero.  The lesson is in the instrument, not the code
// under test, so the instrument is now safe to run either way.
struct CausalityAudit {
    std::atomic<unsigned long long> historyCalls{0}, futureSamples{0};
    std::atomic<unsigned long long> fieldCalls{0}, advancedRoots{0};
    std::atomic<unsigned long long> unconverged{0}, futureAtConvergedRead{0};
    std::atomic<unsigned long long> observationAheadOfPresent{0};
    std::atomic<double> worstFutureSeconds{0.0};
    std::atomic<double> worstAdvancedSeconds{0.0};
    std::atomic<double> worstConvergedFutureSeconds{0.0};
    std::atomic<double> worstLightConeResidual{0.0};
    bool enabled=false;

    // Atomics are not copy-assignable, so clearing needs a member rather than
    // an assignment from a fresh instance.  The validation suite arms, runs
    // and asserts, so it has to be able to start from zero.
    void reset() {
        historyCalls=0; futureSamples=0; fieldCalls=0; advancedRoots=0;
        unconverged=0; futureAtConvergedRead=0; observationAheadOfPresent=0;
        worstFutureSeconds=0.0; worstAdvancedSeconds=0.0;
        worstConvergedFutureSeconds=0.0; worstLightConeResidual=0.0;
    }
};
inline CausalityAudit gCausalityAudit;

// max() on an atomic double, for the "worst seen" fields above.
inline void recordWorst(std::atomic<double>& worst,double candidate) {
    double seen=worst.load(std::memory_order_relaxed);
    while(candidate>seen
          &&!worst.compare_exchange_weak(seen,candidate,
                                         std::memory_order_relaxed)) {}
}

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::StateHistory;
using positronium::objects::cross;
using positronium::objects::dot;
using namespace positronium::parameters;

struct ChargeKinematics { Vec3 position, velocity, acceleration; };

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
inline Vec3 lerp(const Vec3& first,const Vec3& second,double fraction) {
    return first+(second-first)*fraction;
}
#endif

// clampToSegment=false continues the segment's cubic past its ends instead of
// freezing it there; see RetardedSegmentPin below for the one use.
inline ChargeKinematics interpolatedCharge(const State& older, const State& newer,
    bool first, double time, bool clampToSegment=true) {
    const double span = newer.time - older.time;
    if(!(span>0.0)) {
        return {first?newer.firstPosition:newer.secondPosition,
                first?newer.firstVelocity:newer.secondVelocity,
                first?newer.firstAcceleration:newer.secondAcceleration};
    }
    const double rawFraction=(time-older.time)/span;
    const double fraction=clampToSegment?std::clamp(rawFraction,0.0,1.0):rawFraction;
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
inline ChargeKinematics linearlyInterpolatedCharge(const State& older,
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

// One history segment, pinned for a group of retarded reads that must see the
// SAME polynomial.
//
// interpolatedCharge() is C1: velocity is continuous at a history node, the
// cubic's acceleration is not.  On a tilted para orbit at 1.1 r* (floor
// 0.25 r*) the jump was 1-6% of |a| at every one of 119 nodes, and the
// stored node accelerations sat halfway between the left and right limits.
// The Lienard-Wiechert and dipole radiation fields follow that acceleration.
// covariantDipoleGradientForce differences six probes 1e-4 r apart; when ONE
// probe's retarded time crossed a node, its coupling stepped by 5e-4 while the
// other five moved by 1e-7, the central difference jumped by 33%
// (3.68e-2 -> 2.21e-2) within 1e-28 s, and the adaptive integrator failed on
// it at every depth.  The gradient of a field is only meaningful on one
// smooth worldline, so while a pin is set every read of that source's
// charge inside [lower, upper] comes from the pinned segment, continued past
// its ends (the window extends by a few probe offsets over light speed,
// far below the node spacing).
//
// Making the whole history C2 instead (quintic Hermite on the stored
// accelerations) was measured and rejected: it removed the jump but moved
// trajectory-convergence (tol 1e-7 residual 3.9e-7 -> 3.3e-6) and the
// long-horizon radiative balance (0.059 -> 0.174) out of their bands, and
// the tilted trajectories came out no better than with this pin.  Adding
// the radiation-reaction kick to the stored acceleration changed neither
// (0.17363 -> 0.17365), so that is not why; the cause was not found.
struct RetardedSegmentPin {
    const StateHistory* history=nullptr;
    bool first=false;
    std::size_t newerIndex=0;   // == history->size(): the segment ends at present
    double lower=0.0, upper=0.0;
};
inline thread_local RetardedSegmentPin gRetardedSegmentPin;

struct RetardedSegmentPinGuard {
    RetardedSegmentPin saved;
    explicit RetardedSegmentPinGuard(const RetardedSegmentPin& pin)
        : saved(gRetardedSegmentPin) { gRetardedSegmentPin=pin; }
    ~RetardedSegmentPinGuard() { gRetardedSegmentPin=saved; }
    RetardedSegmentPinGuard(const RetardedSegmentPinGuard&)=delete;
    RetardedSegmentPinGuard& operator=(const RetardedSegmentPinGuard&)=delete;
};

// The pinned segment's end points when `time` falls inside an active pin for
// this history and source.
inline bool pinnedSegmentEnds(const StateHistory& history,const State& present,
    bool first,double time,const State*& older,const State*& newer) {
    const RetardedSegmentPin& pin=gRetardedSegmentPin;
    if(pin.history!=&history||pin.first!=first
       ||!(time>=pin.lower&&time<=pin.upper)) return false;
    older=&history[pin.newerIndex-1];
    newer=pin.newerIndex<history.size()?&history[pin.newerIndex]:&present;
    return true;
}

inline ChargeKinematics historicalCharge(const StateHistory& history,
                                   const State& present, bool first,
                                   double time) {
    if(gCausalityAudit.enabled) {
        gCausalityAudit.historyCalls.fetch_add(1,std::memory_order_relaxed);
        // Sampling a source later than the present state is reading the
        // future.  A transient overshoot of an intermediate Newton iterate is
        // discarded by the next one; what must never happen is the CONVERGED
        // read doing it, and that is checked at the call site below.
        const double ahead=time-present.time;
        if(ahead>0.0) {
            gCausalityAudit.futureSamples.fetch_add(1,std::memory_order_relaxed);
            recordWorst(gCausalityAudit.worstFutureSeconds,ahead);
        }
    }
    {
        const State* older=nullptr;
        const State* newer=nullptr;
        if(pinnedSegmentEnds(history,present,first,time,older,newer))
            return interpolatedCharge(*older,*newer,first,time,false);
    }
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

// Charge kinematics PLUS the third derivative, all taken from the one cubic
// Hermite segment historicalCharge() would use at `time`.
//
// Why it exists.  The two-charge limit of a point dipole places two auxiliary
// poles whose charge grows as 1/separation and which solve their retarded
// times separately, so whenever the central retarded time sits within
// separation/c of a node the poles read accelerations from DIFFERENT
// segments, and the node's jump stops cancelling.  Measured on a tilted para
// orbit at 1.1 r*: the retarded time landed 5.9e-30 s from node 106, whose
// acceleration jumps by 0.35%; the dipole field went wrong by a factor up to 9
// inside a window whose width scaled with the pole separation, and the
// adaptive integrator failed on it at every depth.
//
// A cubic is its own third-order Taylor series, so expanding (x, v, a, j)
// from this segment reproduces historicalCharge() EXACTLY anywhere inside the
// segment and continues it smoothly across its ends.
struct ChargeKinematicsWithJerk { Vec3 position, velocity, acceleration, jerk; };

inline ChargeKinematicsWithJerk interpolatedChargeWithJerk(
    const State& older,const State& newer,bool first,double time,
    bool clampToSegment=true) {
    const ChargeKinematics base=
        interpolatedCharge(older,newer,first,time,clampToSegment);
    const double span=newer.time-older.time;
    if(!(span>0.0)) return {base.position,base.velocity,base.acceleration,{}};
    const Vec3 oldPosition=first?older.firstPosition:older.secondPosition;
    const Vec3 newPosition=first?newer.firstPosition:newer.secondPosition;
    const Vec3 oldVelocity=first?older.firstVelocity:older.secondVelocity;
    const Vec3 newVelocity=first?newer.firstVelocity:newer.secondVelocity;
    // Third derivatives of the Hermite basis in the fraction s:
    // h00'''=12, h10'''=6, h01'''=-12, h11'''=6; d/dt = (1/span) d/ds.
    const Vec3 jerk=(oldPosition*12.0+oldVelocity*(6.0*span)
        -newPosition*12.0+newVelocity*(6.0*span))/(span*span*span);
    return {base.position,base.velocity,base.acceleration,jerk};
}

// Same segment selection as historicalCharge(), branch for branch, pin
// included; the extrapolation and zero-span branches have no cubic and
// report zero jerk.
inline ChargeKinematicsWithJerk historicalChargeWithJerk(
    const StateHistory& history,const State& present,bool first,double time) {
    {
        const State* older=nullptr;
        const State* newer=nullptr;
        if(pinnedSegmentEnds(history,present,first,time,older,newer))
            return interpolatedChargeWithJerk(*older,*newer,first,time,false);
    }
    const State& earliest=history.empty()?present:history.front();
    if(time<=earliest.time) {
        const ChargeKinematics base=historicalCharge(history,present,first,time);
        return {base.position,base.velocity,base.acceleration,{}};
    }
    const auto newer=std::lower_bound(history.begin(),history.end(),time,
        [](const State& sample,double requested) {
            return sample.time<requested;
        });
    if(newer==history.end()) {
        const State& latest=history.back();
        if(present.time>latest.time)
            return interpolatedChargeWithJerk(latest,present,first,time);
        return interpolatedChargeWithJerk(latest,latest,first,time);
    }
    if(newer==history.begin())
        return interpolatedChargeWithJerk(*newer,*newer,first,time);
    return interpolatedChargeWithJerk(*std::prev(newer),*newer,first,time);
}

// Pin for a group of reads centred on the retarded time of `observation`.
// Returns an inactive pin (history=nullptr) when that retarded time has no
// two-node segment or the window would reach beyond a quarter of it.
inline RetardedSegmentPin retardedSegmentPinAt(const StateHistory& history,
    const State& present,bool sourceIsFirst,const Vec3& observation,
    double observationTime,double halfWidth) {
    RetardedSegmentPin pin;
    if(history.empty()) return pin;
    double retardedTime=observationTime
        -(observation-historicalCharge(history,present,sourceIsFirst,
                                       observationTime).position).norm()/c;
    for(int iteration=0;iteration<32;++iteration) {
        const ChargeKinematics source=
            historicalCharge(history,present,sourceIsFirst,retardedTime);
        const Vec3 sight=observation-source.position;
        const double distance=sight.norm();
        const Vec3 direction=distance>0.0?sight/distance:Vec3{};
        const double residual=retardedTime+distance/c-observationTime;
        retardedTime-=residual/std::max(1.0e-8,
            1.0-dot(direction,source.velocity/c));
        if(std::abs(residual)<=1.0e-15*std::abs(observationTime)) break;
    }
    if(!(retardedTime>history.front().time)) return pin;
    const auto newer=std::lower_bound(history.begin(),history.end(),
        retardedTime,[](const State& sample,double requested) {
            return sample.time<requested;
        });
    const std::size_t newerIndex=
        static_cast<std::size_t>(newer-history.begin());
    const State& olderState=history[newerIndex-1];
    const double newerTime=newerIndex<history.size()
        ?history[newerIndex].time:present.time;
    if(newerIndex==history.size()&&!(present.time>olderState.time)) return pin;
    if(!(halfWidth<0.25*(newerTime-olderState.time))) return pin;
    pin.history=&history;
    pin.first=sourceIsFirst;
    pin.newerIndex=newerIndex;
    pin.lower=olderState.time-halfWidth;
    pin.upper=newerTime+halfWidth;
    return pin;
}

// Mutual, retarded Lienard-Wiechert field of a moving point charge.
inline ElectromagneticField lienardWiechertField(const Vec3& observationPosition,
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
        // A coincident observation point leaves the direction undefined, and
        // the bare division below then returns NaN, which propagates through
        // every later iterate and out of the loop (audit point 3.1).  At
        // coincidence the light cone is flat: the residual is the time
        // difference alone, and a zero direction makes lightConeDerivative
        // exactly 1 below, which is the correct Newton step for that case.
        // Every retardedDistance above the smallest normal double takes the
        // same arithmetic as before.
        const Vec3 retardedDirection=
            retardedDistance>std::numeric_limits<double>::min()
                ? retardedDisplacement/retardedDistance : Vec3{};
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
    // source was last evaluated at the PREVIOUS iterate, one Newton step
    // behind the retardedTime this loop just converged to (the other two
    // Newton loops on this same light-cone equation, in
    // retardedElectricDipoleField and retardedMagneticDipoleField, both
    // re-fetch their source at the converged time for exactly this reason).
    // Quadratic convergence keeps the resulting position/velocity/
    // acceleration error tiny, but it is real and free to remove.
    source = historicalCharge(history, presentState, sourceIsFirst, retardedTime);
    const Vec3 displacement = observationPosition - source.position;
    const double distance = displacement.norm();
    if(gCausalityAudit.enabled) {
        gCausalityAudit.fieldCalls.fetch_add(1,std::memory_order_relaxed);
        // The field must not be asked for later than the last committed state.
        if(observationTime>presentState.time)
            gCausalityAudit.observationAheadOfPresent.fetch_add(
                1,std::memory_order_relaxed);
        // The advanced root of the same light-cone equation sits after the
        // observation time.  Landing on it is an acausal solution, not a
        // convergence failure, and is counted separately from one.
        const double advanced=retardedTime-observationTime;
        if(advanced>0.0) {
            gCausalityAudit.advancedRoots.fetch_add(1,std::memory_order_relaxed);
            recordWorst(gCausalityAudit.worstAdvancedSeconds,advanced);
        }
        // The converged read itself, which is the one that enters the field.
        const double convergedAhead=retardedTime-presentState.time;
        if(convergedAhead>0.0) {
            gCausalityAudit.futureAtConvergedRead.fetch_add(
                1,std::memory_order_relaxed);
            recordWorst(gCausalityAudit.worstConvergedFutureSeconds,
                        convergedAhead);
        }
        // Light-cone closure: the separation must equal c times the delay.
        const double residual=std::abs(
            (observationTime-retardedTime)*c-distance);
        const double scale=std::max(distance,std::numeric_limits<double>::min());
        recordWorst(gCausalityAudit.worstLightConeResidual,residual/scale);
        if(residual/scale>1.0e-6)
            gCausalityAudit.unconverged.fetch_add(1,std::memory_order_relaxed);
    }
    if(distance<=std::numeric_limits<double>::min()) return {};
    const Vec3 direction = displacement / distance;
    // Trial stages may cross the declared boundary before the enclosing event
    // locator clips the trajectory.  Never evaluate the singular point-charge
    // formula inside a domain whose result is discarded by the model.
    // nuclearCutoff is the pre-existing numerical safety net, kept for every
    // pair.
    //
    // separationFloor() (e+e- only) softens the field the SAME way the
    // energies and the instantaneous forces are softened: the whole field is
    // scaled by (R/R_c)^3 with R_c = sqrt(R^2 + floor^2), so a static source
    // gives exactly -grad of the Plummer potential -k/R_c that
    // conservativeParticleEnergy and coulombForces use.  It used to be a hard
    // clamp, max(R, floor): constant magnitude k/floor^2 inside the floor and
    // unsoftened outside, while the energy was Plummer.  Measured at the
    // default floor on a para passage (audit section 64b): inside one floor
    // this force did +0.51 k/r0 of work while the Coulomb energy moved by
    // -0.001, and that gain, with the dipole sector's, ejected every bound
    // pair that crossed the core.
    const double fieldDistance=std::max(distance,nuclearCutoff);
    const Vec3 beta = source.velocity / c;
    const double betaSquared = beta.squaredNorm();
    const double kappa = std::max(1.0e-8, 1.0 - dot(direction, beta));
    double plummerScale=1.0;
    if(const double floor=separationFloor(); floor>0.0) {
        // The distance in the source's instantaneous rest frame,
        // gamma kappa R = u.(x - x_ret)/c, is a Lorentz scalar, so scaling the
        // field tensor by a function of it keeps the field covariant; for a
        // source at rest it is R itself, so the static limit is still the
        // Plummer gradient.
        const double restDistance=fieldDistance*kappa
            /std::sqrt(std::max(1.0-betaSquared,1.0e-300));
        const double ratio=restDistance
            /std::sqrt(restDistance*restDistance+floor*floor);
        plummerScale=ratio*ratio*ratio;
    }
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
                        * (coulomb * sourceCharge*formFactor*plummerScale);
    return {electric, cross(direction, electric) / c};
}
