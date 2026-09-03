#pragma once

#include <cmath>

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
    double norm() const { return std::sqrt(squaredNorm()); }
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
