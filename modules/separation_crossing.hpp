#pragma once

// Solves for the interpolation fraction in [0,1] at which the separation
// between two sampled states first crosses a target distance -- used to
// locate a boundary crossing (a collision, a capture threshold) between two
// already-integrated steps without re-running the integrator at finer
// resolution.
//
// Self-contained and order-independent.  It names what it needs through
// using-declarations rather than reopening namespace positronium: the header
// is still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.

#include "state.hpp"
#include "vector3.hpp"

#include <algorithm>
#include <cmath>

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::dot;

inline double separationCrossingFraction(const State& before,const State& after,
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
