#pragma once

// Relativistic point-particle kinematics (Lorentz factor, momentum, kinetic
// energy, the momentum-to-velocity inverse) and the lab/proper covariant
// dipole synchronization built on top of them.
//
// Self-contained and order-independent.  It names what it needs through
// using-declarations rather than reopening namespace positronium: the header
// is still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.
//
// interpolateState/interpolateVector/interpolateDipole and isFinite stay in
// modules/state_validity_interpolation.hpp rather than moving here.

#include "dipole_tensor.hpp"
#include "lorentz_boost_dipole.hpp"
#include "physical_constants.hpp"
#include "state.hpp"
#include "two_body_kinematics.hpp"
#include "vector3.hpp"

#include <cmath>

namespace two_body = positronium::kinematics;

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::DipoleTensor;
using positronium::parameters::c;

inline double gamma(const Vec3& velocity) {
    return two_body::strictLorentzFactor(velocity);
}

inline Vec3 momentum(const Vec3& velocity, double mass) { return velocity * (gamma(velocity) * mass); }

inline double kineticEnergy(const Vec3& velocity,double mass) {
    return two_body::kineticEnergyFromVelocity(velocity,mass);
}

inline Vec3 velocityFromMomentum(const Vec3& momentum, double mass) {
    const double gammaFromMomentum = std::sqrt(1.0 + momentum.squaredNorm() / (mass*mass*c*c));
    return momentum / (gammaFromMomentum * mass);
}

inline void synchronizeCovariantDipoles(State& state) {
    if(state.firstProperDipole.squaredNorm()==0.0
        &&state.firstDipole.squaredNorm()>0.0)
        state.firstProperDipole=state.firstDipole;
    if(state.secondProperDipole.squaredNorm()==0.0
        &&state.secondDipole.squaredNorm()>0.0)
        state.secondProperDipole=state.secondDipole;
    const DipoleTensor firstLab=lorentzBoostDipole(
        {{},state.firstProperDipole},state.firstVelocity);
    const DipoleTensor secondLab=lorentzBoostDipole(
        {{},state.secondProperDipole},state.secondVelocity);
    state.firstElectricDipole=firstLab.electric;
    state.firstDipole=firstLab.magnetic;
    state.secondElectricDipole=secondLab.electric;
    state.secondDipole=secondLab.magnetic;
}
