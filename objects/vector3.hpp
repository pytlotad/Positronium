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
    constexpr double squaredNorm() const { return x*x+y*y+z*z; }
    double norm() const { return std::sqrt(squaredNorm()); }
};

} // namespace positronium::objects
