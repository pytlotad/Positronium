#pragma once

// Pair separation geometry and its regularization: the TRUE geometry used
// for termination/diagnostics, and the floored, Plummer-softened variant
// every conservative force/energy law reads instead so nothing diverges on
// a close encounter.  This is the other prerequisite (besides dot()/cross())
// that modules/retarded_charge_kinematics.hpp -- and essentially everything
// downstream of it -- builds on.
//
// Self-contained and order-independent.  It names what it needs through a
// using-directive on positronium::parameters and using-declarations for the
// object types, rather than reopening namespace positronium: the header is
// still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.

#include "pair_configuration.hpp"
#include "particle_species.hpp"
#include "physical_constants.hpp"
#include "state.hpp"
#include "vector3.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::dot;
using namespace positronium::parameters;

// Rescales a source-to-target separation to an effective distance
// sqrt(r^2+floor^2), keeping its direction -- Plummer-style gravitational
// softening, tied to a physically motivated radius (the Compton barrier for
// e+e-) instead of an arbitrary softening length. Force/field laws fed this
// vector see r_eff -> floor (bounded) as r -> 0 and r_eff -> r (unchanged)
// for r >> floor, with EVERY derivative continuous in between: unlike a hard
// clamp to floor (this function's previous form, r_eff = max(r,floor)),
// there is no kink in the force at r = floor. That kink was measured
// directly to stall the adaptive integrator's step-doubling error probe
// (repeated failed subdivision, maximumDepth exhausted) exactly where an
// inspiralling orbit first crossed the old hard floor -- a structural
// property of a piecewise-constant regularization, not a tolerance or
// step-size problem, so no amount of retrying fixed it.
//
// The correction this makes to the TRUE force is negligible except near
// floor: at r = 10*floor it is already only 0.5%, and at production
// (Bohr-radius) separations, ~500x floor, it is a few parts in 10^7 --
// several orders below every other tolerance already accepted in this
// engine. It is not confined to r < floor the way the hard clamp was; it is
// evaluated identically for every r, which is what removes the kink (a
// switch between "identity" and "constant" at one point is exactly where a
// derivative discontinuity has to live).
//
// True coincidence has no direction to preserve; any fixed one works, since
// every consumer only reads the returned MAGNITUDE from it in that case, and
// that magnitude is exactly floor there (sqrt(0+floor^2)), so the general
// formula already gives the right answer -- the special case below exists
// only to avoid dividing trueSeparation by a zero trueDistance.
inline Vec3 clampedSeparationVector(const Vec3& trueSeparation, double floor) {
    if (!(floor > 0.0)) return trueSeparation;
    const double trueDistance = trueSeparation.norm();
    const double effectiveDistance =
        std::sqrt(trueDistance*trueDistance + floor*floor);
    if (!(trueDistance > std::numeric_limits<double>::min())) {
        return Vec3{effectiveDistance, 0.0, 0.0};
    }
    return trueSeparation * (effectiveDistance / trueDistance);
}

// The floor itself: comptonBarrierRadius for e+e-, since it is derived from
// the electron's mass and g-factor specifically (see physical_constants.hpp)
// and would not mean the same thing for any other pair. No regularization
// for other pairs yet -- this is deliberately narrower than "every pair"
// until an analogous barrier is derived for them.
inline double separationFloor() {
    return isPositronium(activePair) ? comptonBarrierRadius : 0.0;
}

struct PairGeometry {
    Vec3 firstMinusSecond;
    double distance;
    double inverseDistance;
    double inverseDistanceCubed;
    double inverseDistanceFourth;
};

// TRUE geometry -- used for termination checks, minimumSeparation reporting,
// step-error normalization, and every other diagnostic that must say where
// the pair actually is, not where the force laws pretend it is.  Force laws
// use clampedPairGeometry() below instead; nothing here changes for them.
inline PairGeometry pairGeometry(const State& s) {
    const Vec3 firstMinusSecond = s.firstPosition - s.secondPosition;
    const double distanceSquared = firstMinusSecond.squaredNorm();
    const double distance = std::sqrt(distanceSquared);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;
    return {firstMinusSecond, distance, inverseDistance,
            inverseDistanceSquared * inverseDistance,
            inverseDistanceSquared * inverseDistanceSquared};
}

// Same geometry, but with the separation floored at separationFloor() before
// any inverse power is taken.  Every CONSERVATIVE force/energy law
// (Coulomb, Darwin) reads this instead of pairGeometry(), so none of them
// diverge on a close encounter; termination logic and diagnostics keep
// reading the true pairGeometry() above so they still report what actually
// happened.
inline PairGeometry clampedPairGeometry(const State& s) {
    const Vec3 firstMinusSecond = clampedSeparationVector(
        s.firstPosition - s.secondPosition, separationFloor());
    const double distanceSquared = firstMinusSecond.squaredNorm();
    const double distance = std::sqrt(distanceSquared);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;
    return {firstMinusSecond, distance, inverseDistance,
            inverseDistanceSquared * inverseDistance,
            inverseDistanceSquared * inverseDistanceSquared};
}

inline double separation(const State& s) { return pairGeometry(s).distance; }
