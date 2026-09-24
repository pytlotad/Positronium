#pragma once

// Forward-mode automatic differentiation, one independent variable.
//
// This exists for one job: the material derivative D/Ds of the retarded
// dipole field along a target's worldline (hiddenMomentumRateForce).  That
// field is built as the limit of two point charges whose separation is 1e-5
// of the source distance, so the two pole fields are about 1e5 times the
// dipole field they differ by.  Differencing THAT numerically divides a
// cancellation floor of ~1e-10 of the field by a step of ~1/2000 of the
// field's own time scale and lands at 2e-05 (audit section 161).
//
// A dual number carries the derivative alongside the value through the same
// arithmetic, so the subtraction happens on derivatives that were each formed
// exactly rather than on a difference quotient.  The floor then stays at the
// pole cancellation itself instead of being multiplied by 1/h.
//
// WHY NOT A NUMERIC ALTERNATIVE.  A higher-order stencil does not help: the
// error is round-off amplified by 1/h, not truncation, so Richardson
// extrapolation makes it worse.  Complex-step differentiation would remove
// the cancellation too, but the field path clamps with std::max and branches
// on comparisons, and complex-step is only valid through analytic operations;
// duals carry the selected branch's derivative correctly and are exact for
// the piecewise-smooth expression the code actually evaluates.
//
// WHAT IT IS NOT.  Non-smooth points are not handled specially: at a std::max
// switch the derivative is the selected branch's one-sided derivative, which
// is what the production field itself has there.  This is deliberate -- the
// clamps exist to bound the field and the code downstream already treats
// those geometries as regularized, not as differentiable.

#include <cmath>

namespace positronium::dual {

// value and its derivative with respect to the one independent variable.
struct Dual {
    double value=0.0, derivative=0.0;
    constexpr Dual()=default;
    constexpr Dual(double v):value(v),derivative(0.0) {}
    constexpr Dual(double v,double d):value(v),derivative(d) {}
};

constexpr Dual operator+(Dual a,Dual b) {
    return {a.value+b.value,a.derivative+b.derivative};
}
constexpr Dual operator-(Dual a,Dual b) {
    return {a.value-b.value,a.derivative-b.derivative};
}
constexpr Dual operator-(Dual a) { return {-a.value,-a.derivative}; }
constexpr Dual operator*(Dual a,Dual b) {
    return {a.value*b.value,a.derivative*b.value+a.value*b.derivative};
}
constexpr Dual operator/(Dual a,Dual b) {
    const double inverse=1.0/b.value;
    return {a.value*inverse,
            (a.derivative-a.value*inverse*b.derivative)*inverse};
}
constexpr Dual& operator+=(Dual& a,Dual b) { a=a+b; return a; }
constexpr Dual& operator-=(Dual& a,Dual b) { a=a-b; return a; }
constexpr Dual& operator*=(Dual& a,Dual b) { a=a*b; return a; }

inline Dual sqrtDual(Dual a) {
    const double root=std::sqrt(a.value);
    // d(sqrt x) = dx/(2 sqrt x).  At zero the derivative is unbounded; the
    // callers reach this only through squared norms they have already
    // rejected as degenerate, so returning zero there is a guard and not a
    // silent approximation.
    return {root,root>0.0?a.derivative/(2.0*root):0.0};
}

// max against a CONSTANT bound, carrying the selected branch's derivative.
// The bound is a clamp the field path applies (nuclearCutoff, a softening
// length); it does not move with the independent variable, so the derivative
// is either the argument's or exactly zero.
inline Dual maxDual(Dual a,double bound) {
    return a.value>=bound?a:Dual{bound,0.0};
}

// max between two quantities that BOTH move with the variable.  Ties pick
// the first, matching std::max's own tie-breaking, so a value comparison
// against the production path cannot disagree on the boundary.
inline Dual maxDual(Dual a,Dual b) { return a.value>=b.value?a:b; }

struct DualVec3 {
    Dual x,y,z;
    constexpr DualVec3()=default;
    constexpr DualVec3(Dual a,Dual b,Dual d):x(a),y(b),z(d) {}
};

constexpr DualVec3 operator+(const DualVec3& a,const DualVec3& b) {
    return {a.x+b.x,a.y+b.y,a.z+b.z};
}
constexpr DualVec3 operator-(const DualVec3& a,const DualVec3& b) {
    return {a.x-b.x,a.y-b.y,a.z-b.z};
}
constexpr DualVec3 operator*(const DualVec3& a,Dual s) {
    return {a.x*s,a.y*s,a.z*s};
}
constexpr DualVec3 operator/(const DualVec3& a,Dual s) {
    return {a.x/s,a.y/s,a.z/s};
}
constexpr DualVec3& operator+=(DualVec3& a,const DualVec3& b) {
    a=a+b; return a;
}
constexpr Dual dotDual(const DualVec3& a,const DualVec3& b) {
    return a.x*b.x+a.y*b.y+a.z*b.z;
}
constexpr DualVec3 crossDual(const DualVec3& a,const DualVec3& b) {
    return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
constexpr Dual squaredNormDual(const DualVec3& a) { return dotDual(a,a); }
inline Dual normDual(const DualVec3& a) { return sqrtDual(squaredNormDual(a)); }

}  // namespace positronium::dual
