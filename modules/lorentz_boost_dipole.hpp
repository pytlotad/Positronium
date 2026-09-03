#pragma once

// The covariant electric/magnetic dipole tensor's Lorentz boost between the
// particle's rest frame and the lab frame -- what
// modules/relativistic_kinematics.hpp's synchronizeCovariantDipoles calls
// to turn a rest-frame moment into the lab-frame firstDipole/secondDipole
// (magnetic) and firstElectricDipole/secondElectricDipole (electric) pair.
//
// Self-contained and order-independent.  It names what it needs through
// using-declarations rather than reopening namespace positronium: the header
// is still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.
//
// The Lorentz factor is taken straight from kinematics::strictLorentzFactor
// rather than through the gamma() alias in relativistic_kinematics.hpp.
// gamma() is a one-line forwarder to exactly this function, so the arithmetic
// is unchanged, and going to the source breaks a cycle: that header's
// synchronizeCovariantDipoles calls lorentzBoostDipole in turn.

#include "dipole_tensor.hpp"
#include "physical_constants.hpp"
#include "two_body_kinematics.hpp"
#include "vector3.hpp"

using positronium::objects::Vec3;
using positronium::objects::DipoleTensor;
using positronium::objects::cross;
using positronium::objects::dot;
using positronium::parameters::c;

inline DipoleTensor lorentzBoostDipole(const DipoleTensor& dipole,
                                       const Vec3& frameVelocity) {
    const double speedSquared=frameVelocity.squaredNorm();
    if(speedSquared==0.0) return dipole;
    const double boostGamma=
        positronium::kinematics::strictLorentzFactor(frameVelocity);
    const double longitudinalFactor=boostGamma*boostGamma
        /(boostGamma+1.0)/(c*c);
    return {
        (dipole.electric+cross(frameVelocity,dipole.magnetic)/(c*c))
            *boostGamma
            -frameVelocity*(longitudinalFactor
                *dot(frameVelocity,dipole.electric)),
        (dipole.magnetic-cross(frameVelocity,dipole.electric))*boostGamma
            -frameVelocity*(longitudinalFactor
                *dot(frameVelocity,dipole.magnetic))};
}
