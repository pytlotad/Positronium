#pragma once

#include <deque>

#include "vector3.hpp"
#include "dipole_tensor.hpp"

namespace positronium::objects {

struct State {
    Vec3 firstPosition,secondPosition;
    Vec3 firstVelocity,secondVelocity;
    Vec3 firstAcceleration,secondAcceleration;
    Vec3 firstDipole,secondDipole;
    Vec3 firstElectricDipole,secondElectricDipole;
    // Magnetic moments in the instantaneous particle rest frames.  These are
    // the independent degrees of freedom; the laboratory p and mu above are
    // reconstructed from the covariant tensor after every update.
    Vec3 firstProperDipole,secondProperDipole;
    double time=0.0;
    double radiatedEnergy=0.0;
    double dipoleConstraintEnergy=0.0;
    Vec3 radiatedMomentum,radiatedAngularMomentum;
    // Reconstructed bound/interference field reservoir required to close the
    // particle-plus-field conservation laws on the control surface.
    double boundFieldEnergy=0.0;
    Vec3 boundFieldMomentum,boundFieldAngularMomentum;
    // Independent integral of the mismatch between individual LL reaction
    // and the coherent charge-sector far-field flux.
    double reactionEnergyMismatch=0.0;
    Vec3 reactionMomentumMismatch,reactionAngularMomentumMismatch;
};

using StateHistory=std::deque<State>;

} // namespace positronium::objects
