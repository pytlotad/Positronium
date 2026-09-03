#pragma once

// Finiteness checking and time interpolation for State: isFinite() guards
// every integration step and terminal event against a state that has gone
// non-finite (NaN/Inf, or superluminal for a velocity); interpolateState()
// (built on interpolateVector()/interpolateDipole()) reconstructs the state
// at an interior time between two integrated samples, which is how frames,
// terminal events and historical samples are produced without re-running
// the integrator.
//
// Self-contained and order-independent.  It names what it needs through
// using-declarations rather than reopening namespace positronium: the header
// is still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.

#include "dipole_tensor.hpp"
#include "state.hpp"
#include "two_body_kinematics.hpp"
#include "vector3.hpp"

#include <cmath>

namespace two_body = positronium::kinematics;

using positronium::objects::Vec3;
using positronium::objects::State;

inline bool isFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

inline bool isFinite(const State& state) {
    return isFinite(state.firstPosition) && isFinite(state.secondPosition)
        && isFinite(state.firstVelocity) && isFinite(state.secondVelocity)
        && two_body::subluminal(state.firstVelocity)
        && two_body::subluminal(state.secondVelocity)
        && isFinite(state.firstAcceleration) && isFinite(state.secondAcceleration)
        && isFinite(state.firstDipole) && isFinite(state.secondDipole)
        &&isFinite(state.firstElectricDipole)
        &&isFinite(state.secondElectricDipole)
        &&isFinite(state.firstProperDipole)
        &&isFinite(state.secondProperDipole)
        && isFinite(state.radiatedMomentum)
        && isFinite(state.radiatedAngularMomentum)
        && isFinite(state.boundFieldMomentum)
        && isFinite(state.boundFieldAngularMomentum)
        && isFinite(state.previousFluxMomentum)
        && isFinite(state.previousFluxAngularMomentum)
        && isFinite(state.previousMismatchMomentum)
        && isFinite(state.previousMismatchAngularMomentum)
        && isFinite(state.reactionMomentumMismatch)
        && isFinite(state.reactionAngularMomentumMismatch)
        && std::isfinite(state.time) && std::isfinite(state.radiatedEnergy)
        && std::isfinite(state.orbitalRadiatedEnergy)
        && std::isfinite(state.dipoleConstraintEnergy)
        && std::isfinite(state.zeroPointPhase)
        && std::isfinite(state.boundFieldEnergy)
        && std::isfinite(state.reactionEnergyMismatch)
        && std::isfinite(state.previousStepDt)
        && std::isfinite(state.previousFluxEnergy)
        && std::isfinite(state.previousDipoleFluxEnergy)
        && std::isfinite(state.previousMismatchEnergy);
}

inline Vec3 interpolateVector(const Vec3& before, const Vec3& after, double fraction) {
    return before + (after - before) * fraction;
}

inline Vec3 interpolateDipole(const Vec3& before, const Vec3& after, double fraction) {
    Vec3 result = interpolateVector(before, after, fraction);
    const double targetNorm = before.norm() + (after.norm() - before.norm()) * fraction;
    const double resultNorm = result.norm();
    if (resultNorm > 0.0) result = result * (targetNorm / resultNorm);
    return result;
}

inline State interpolateState(const State& before, const State& after, double fraction) {
    // At an endpoint every field, including the committed-step cache, has an
    // exact meaning.  Returning the operand also avoids avoidable round-off.
    if(fraction==0.0) return before;
    if(fraction==1.0) return after;

    // Preserve a complete fallback if State later gains another member.  All
    // current continuous members are explicitly overwritten below; the
    // discrete previous* cache is explicitly invalidated at the end.
    State result=before;
    result.firstPosition = interpolateVector(
        before.firstPosition, after.firstPosition, fraction);
    result.secondPosition = interpolateVector(
        before.secondPosition, after.secondPosition, fraction);
    result.firstVelocity = interpolateVector(
        before.firstVelocity, after.firstVelocity, fraction);
    result.secondVelocity = interpolateVector(
        before.secondVelocity, after.secondVelocity, fraction);
    result.firstAcceleration = interpolateVector(
        before.firstAcceleration, after.firstAcceleration, fraction);
    result.secondAcceleration = interpolateVector(
        before.secondAcceleration, after.secondAcceleration, fraction);
    result.firstDipole = interpolateDipole(
        before.firstDipole, after.firstDipole, fraction);
    result.secondDipole = interpolateDipole(
        before.secondDipole, after.secondDipole, fraction);
    result.firstElectricDipole=interpolateVector(
        before.firstElectricDipole,after.firstElectricDipole,fraction);
    result.secondElectricDipole=interpolateVector(
        before.secondElectricDipole,after.secondElectricDipole,fraction);
    result.firstProperDipole=interpolateDipole(
        before.firstProperDipole,after.firstProperDipole,fraction);
    result.secondProperDipole=interpolateDipole(
        before.secondProperDipole,after.secondProperDipole,fraction);
    result.time = before.time + (after.time - before.time) * fraction;
    result.radiatedEnergy = before.radiatedEnergy
        + (after.radiatedEnergy - before.radiatedEnergy) * fraction;
    result.orbitalRadiatedEnergy=before.orbitalRadiatedEnergy
        +(after.orbitalRadiatedEnergy-before.orbitalRadiatedEnergy)*fraction;
    result.dipoleConstraintEnergy=before.dipoleConstraintEnergy
        +(after.dipoleConstraintEnergy-before.dipoleConstraintEnergy)*fraction;
    result.zeroPointPhase=before.zeroPointPhase
        +(after.zeroPointPhase-before.zeroPointPhase)*fraction;
    result.radiatedMomentum = interpolateVector(
        before.radiatedMomentum, after.radiatedMomentum, fraction);
    result.radiatedAngularMomentum = interpolateVector(
        before.radiatedAngularMomentum, after.radiatedAngularMomentum, fraction);
    result.boundFieldEnergy=before.boundFieldEnergy
        +(after.boundFieldEnergy-before.boundFieldEnergy)*fraction;
    result.boundFieldMomentum=interpolateVector(
        before.boundFieldMomentum,after.boundFieldMomentum,fraction);
    result.boundFieldAngularMomentum=interpolateVector(
        before.boundFieldAngularMomentum,after.boundFieldAngularMomentum,fraction);
    result.reactionEnergyMismatch=before.reactionEnergyMismatch
        +(after.reactionEnergyMismatch-before.reactionEnergyMismatch)*fraction;
    result.reactionMomentumMismatch=interpolateVector(
        before.reactionMomentumMismatch,after.reactionMomentumMismatch,fraction);
    result.reactionAngularMomentumMismatch=interpolateVector(
        before.reactionAngularMomentumMismatch,
        after.reactionAngularMomentumMismatch,fraction);

    // previous* is an atomic cache of the last committed integration substep,
    // used by the trapezoidal flux correction.  Neither endpoint cache nor a
    // linear blend describes an interior event, especially when an adaptive
    // step consists of several accepted half-steps.  Invalidate the whole
    // cache so continuing from an interpolated state safely restarts the
    // quadrature instead of applying a correction from an unrelated step.
    result.previousStepDt=0.0;
    result.hasPreviousRates=false;
    result.previousFluxEnergy=0.0;
    result.previousDipoleFluxEnergy=0.0;
    result.previousFluxMomentum={};
    result.previousFluxAngularMomentum={};
    result.previousMismatchEnergy=0.0;
    result.previousMismatchMomentum={};
    result.previousMismatchAngularMomentum={};
    return result;
}
