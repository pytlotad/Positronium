#pragma once

// Finiteness checking and time interpolation for State: isFinite() guards
// every integration step and terminal event against a state that has gone
// non-finite (NaN/Inf, or superluminal for a velocity); interpolateState()
// (built on interpolateVector()/interpolateDipole()) reconstructs the state
// at an interior time between two integrated samples, which is how frames,
// terminal events and historical samples are produced without re-running
// the integrator.
//
// Extracted verbatim from positronium.cpp (Stage 0 of splitting engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied, so it depends on that namespace
// already having in scope: Vec3, State, two_body::subluminal().  Not yet a
// standalone, order-independent header.

bool isFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) && std::isfinite(value.z);
}

bool isFinite(const State& state) {
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
        && isFinite(state.reactionMomentumMismatch)
        && isFinite(state.reactionAngularMomentumMismatch)
        && std::isfinite(state.time) && std::isfinite(state.radiatedEnergy)
        && std::isfinite(state.dipoleConstraintEnergy)
        && std::isfinite(state.zeroPointPhase)
        && std::isfinite(state.boundFieldEnergy)
        && std::isfinite(state.reactionEnergyMismatch);
}

Vec3 interpolateVector(const Vec3& before, const Vec3& after, double fraction) {
    return before + (after - before) * fraction;
}

Vec3 interpolateDipole(const Vec3& before, const Vec3& after, double fraction) {
    Vec3 result = interpolateVector(before, after, fraction);
    const double targetNorm = before.norm() + (after.norm() - before.norm()) * fraction;
    const double resultNorm = result.norm();
    if (resultNorm > 0.0) result = result * (targetNorm / resultNorm);
    return result;
}

State interpolateState(const State& before, const State& after, double fraction) {
    State result;
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
    result.dipoleConstraintEnergy=before.dipoleConstraintEnergy
        +(after.dipoleConstraintEnergy-before.dipoleConstraintEnergy)*fraction;
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
    return result;
}
