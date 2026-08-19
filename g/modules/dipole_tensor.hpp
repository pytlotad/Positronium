#pragma once

#include "vector3.hpp"

namespace positronium::objects {

// SI decomposition of an antisymmetric polarization-magnetization tensor.
// electric has units C m and magnetic has units A m^2.  Keeping both parts
// is essential: a pure magnetic dipole in one inertial frame is generally a
// mixed electric-magnetic dipole in another.
struct DipoleTensor {
    Vec3 electric;
    Vec3 magnetic;
};

} // namespace positronium::objects
