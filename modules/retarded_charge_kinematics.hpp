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
#include <cmath>
#include <limits>

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

inline ChargeKinematics interpolatedCharge(const State& older, const State& newer, bool first, double time) {
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

inline ChargeKinematics historicalCharge(const StateHistory& history,
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
    if(distance<=std::numeric_limits<double>::min()) return {};
    const Vec3 direction = displacement / distance;
    // Trial stages may cross the declared boundary before the enclosing event
    // locator clips the trajectory.  Never evaluate the singular point-charge
    // formula inside a domain whose result is discarded by the model.
    // nuclearCutoff is the pre-existing numerical safety net (kept for every
    // pair); separationFloor() is the physically motivated one (e+e- only,
    // and larger, so it dominates there). Whichever is bigger wins.
    const double fieldDistance=std::max({distance,nuclearCutoff,separationFloor()});
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
