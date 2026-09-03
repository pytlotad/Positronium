#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace positronium::objects {

struct Vec3 {
    double x=0.0,y=0.0,z=0.0;

    constexpr Vec3()=default;
    constexpr Vec3(double xValue,double yValue,double zValue)
        :x(xValue),y(yValue),z(zValue) {}
    constexpr Vec3 operator+(const Vec3& value) const {
        return {x+value.x,y+value.y,z+value.z};
    }
    constexpr Vec3 operator-(const Vec3& value) const {
        return {x-value.x,y-value.y,z-value.z};
    }
    constexpr Vec3 operator*(double scale) const {
        return {x*scale,y*scale,z*scale};
    }
    constexpr Vec3 operator/(double scale) const {
        return {x/scale,y/scale,z/scale};
    }
    constexpr Vec3& operator+=(const Vec3& value) {
        x+=value.x;y+=value.y;z+=value.z;return *this;
    }
    constexpr Vec3& operator-=(const Vec3& value) {
        x-=value.x;y-=value.y;z-=value.z;return *this;
    }
    constexpr Vec3 operator-() const {
        return {-x,-y,-z};
    }
    constexpr double squaredNorm() const { return x*x+y*y+z*z; }

    // Squaring costs half the exponent range, so a plain
    // sqrt(x*x+y*y+z*z) cannot see a vector whose components sit below about
    // 1e-154: the sum of squares underflows, and below ~1e-162 it reaches
    // exactly zero, so a nonzero vector reported a zero length.  The mirror
    // case is a component above ~1e+154, where the sum overflows to infinity.
    //
    // The ordinary path is unchanged, deliberately and bit for bit: when the
    // sum of squares is a normal, finite double the square root of it is
    // correctly rounded, and that is the path every quantity this model
    // computes takes -- positions near 1e-11 m, radiated momenta near
    // 1e-40 kg m/s, magnetic moments near 1e-23 J/T are all dozens of orders
    // of magnitude clear of either threshold.  The cost there is two
    // comparisons on a value already in a register.
    //
    // Only when the sum has lost the answer does this rescale by the largest
    // component: the ratios are then in [0,1], their sum of squares in [0,3],
    // and a square root of that can neither overflow nor underflow, so one
    // multiply restores the true magnitude.
    // The rescuing branch is deliberately a separate, cold, never-inlined
    // function.  Leaving it in the body changed how the optimizer inlined
    // norm() into its callers, and with it which multiply-adds got contracted,
    // which moved the last bits of results everywhere -- a chaotic inspiral
    // amplifies that into the fourth digit.  Kept out of line, the hot path
    // compiles to what it always did.
    [[gnu::noinline, gnu::cold]] double rescaledNorm() const {
        if(x!=x||y!=y||z!=z) {                                // NaN in, NaN out
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double largest=std::max({std::abs(x),std::abs(y),std::abs(z)});
        if(largest==0.0) return 0.0;                          // the zero vector
        if(largest>std::numeric_limits<double>::max())
            return largest;                                   // infinite component
        const Vec3 scaled{x/largest,y/largest,z/largest};
        return largest*std::sqrt(scaled.squaredNorm());
    }

    double norm() const {
        const double sumOfSquares=squaredNorm();
        if(sumOfSquares>=std::numeric_limits<double>::min()
           &&sumOfSquares<=std::numeric_limits<double>::max())
            return std::sqrt(sumOfSquares);
        return rescaledNorm();
    }
};

constexpr Vec3 operator*(double scale, const Vec3& value) {
    return value * scale;
}

constexpr double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

} // namespace positronium::objects
