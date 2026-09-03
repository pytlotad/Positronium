#pragma once

// Basic relativistic field/vector primitives: a unit vector, the
// four-vector and electromagnetic-field types, and the Minkowski inner
// product.  ElectromagneticField in particular is the return type
// modules/retarded_charge_kinematics.hpp's lienardWiechertField and its
// electrodynamics.hpp siblings all share.
//
// Self-contained and order-independent: it pulls in its own Vec3 and dot()
// and names them through using-declarations rather than reopening
// namespace positronium.  That distinction matters, because this header is
// still textually included inside positronium.cpp's anonymous namespace:
// reopening a named namespace there would create {anonymous}::positronium
// and hide the real ::positronium from every later lookup.

#include "vector3.hpp"

#include <cmath>
#include <limits>

using positronium::objects::Vec3;
using positronium::objects::dot;

// The zero-vector guard replaces a bare v/v.norm(), which returned NaN for a
// null argument (audit point 3.2).  No production call site passes a zero
// vector, so the returned values are unchanged wherever unit() was defined.
//
// The test is > 0 rather than > DBL_MIN because Vec3::norm() no longer
// underflows: it now reports the true length of a vector whose components are
// too small to square, so this can divide by it.  Every component is at most
// the norm, so the ratios cannot overflow even when the norm is subnormal.
inline Vec3 unit(const Vec3& v) {
    const double n = v.norm();
    return n > 0.0 ? v / n : Vec3{};
}

struct FourVector { double time=0.0; Vec3 space; }; // x^0=ct convention
struct ElectromagneticField { Vec3 electric, magnetic; };

inline double minkowskiDot(const FourVector& first,const FourVector& second) {
    return first.time*second.time-dot(first.space,second.space);
}
