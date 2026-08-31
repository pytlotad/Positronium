#pragma once

// Basic relativistic field/vector primitives: a unit vector, the
// four-vector and electromagnetic-field types, and the Minkowski inner
// product.  ElectromagneticField in particular is the return type
// modules/retarded_charge_kinematics.hpp's lienardWiechertField and its
// electrodynamics.hpp siblings all share.
//
// Extracted verbatim from positronium.cpp (Stage 0 of splitting engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied, so it depends on that namespace
// already having in scope: Vec3 and dot().  Not yet a standalone,
// order-independent header.

Vec3 unit(const Vec3& v) { return v / v.norm(); }
struct FourVector { double time=0.0; Vec3 space; }; // x^0=ct convention
struct ElectromagneticField { Vec3 electric, magnetic; };
double minkowskiDot(const FourVector& first,const FourVector& second) {
    return first.time*second.time-dot(first.space,second.space);
}
