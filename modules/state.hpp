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
    // Same accumulation as radiatedEnergy, but the E1 (charge-charge,
    // orbital) far-zone flux only, with the M1 (magnetic-dipole/spin-
    // precession) channel excluded rather than merged in.  Exists because
    // radiatedEnergy's merge is exactly what makes it useless as an
    // orbital-decay diagnostic near/below the Compton barrier: M1 dominates
    // it by orders of magnitude there (does not recoil the orbit at all),
    // while orbitalRadiatedEnergy is computed from the exact retarded far-
    // zone Poynting integral and does not share conservativeParticleEnergy's
    // Darwin/near-field approximation, so it stays meaningful exactly where
    // that approximation (period comparable to the light-crossing time)
    // breaks down.  Diagnostic-only: nothing currently reads it in
    // production, and it costs nothing extra to accumulate once the flux
    // quadrature (computeOutwardFlux) is already being computed for
    // radiatedEnergy's own sake.
    double orbitalRadiatedEnergy=0.0;
    double dipoleConstraintEnergy=0.0;
    // Accumulated phase of the orbit-following zero-point field, integral of
    // the osculating orbital angular frequency.  Carried in the state because
    // a band that tracks the orbit has a CHANGING mode frequency, and a phase
    // written as -omega(t)*t would jump whenever omega moved, injecting energy
    // out of nowhere.  Integrating d(phase)/dt = omega keeps it continuous.
    // Zero and unused unless --zpf is on.
    double zeroPointPhase=0.0;
    Vec3 radiatedMomentum,radiatedAngularMomentum;
    // Reconstructed bound/interference field reservoir required to close the
    // particle-plus-field conservation laws on the control surface.
    double boundFieldEnergy=0.0;
    Vec3 boundFieldMomentum,boundFieldAngularMomentum;
    // Independent integral of the mismatch between individual LL reaction
    // and the coherent charge-sector far-field flux.
    double reactionEnergyMismatch=0.0;
    // Previous COMMITTED step length and the rates sampled at its start.
    // They exist only to give the radiation bookkeeping trapezoidal weights,
    // and they live in the State so a rejected trial step discards them with
    // everything else.
    double previousStepDt=0.0;
    bool hasPreviousRates=false;
    double previousFluxEnergy=0.0,previousDipoleFluxEnergy=0.0;
    Vec3 previousFluxMomentum,previousFluxAngularMomentum;
    double previousMismatchEnergy=0.0;
    Vec3 previousMismatchMomentum,previousMismatchAngularMomentum;
    Vec3 reactionMomentumMismatch,reactionAngularMomentumMismatch;
};

using StateHistory=std::deque<State>;

} // namespace positronium::objects
