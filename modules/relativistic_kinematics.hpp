// Relativistic point-particle kinematics (Lorentz factor, momentum, kinetic
// energy, the momentum-to-velocity inverse) and the lab/proper covariant
// dipole synchronization built on top of them.
//
// Extracted verbatim from positronium.cpp (Stage 0 of splitting engine,
// experiments and ROOT presentation apart -- see the session notes on why
// and what stays behind; interpolateState/interpolateVector/interpolateDipole
// and isFinite are deliberately left where they are for now, since an
// unrelated in-progress change on this machine touches them).  Textually
// included at the same point inside positronium.cpp's shared anonymous
// namespace it always occupied, so it depends on that namespace already
// having in scope: Vec3, State, DipoleTensor, the c constant, two_body::
// (strictLorentzFactor, kineticEnergyFromVelocity) and lorentzBoostDipole.
// Not yet a standalone, order-independent header.

double gamma(const Vec3& velocity) {
    return two_body::strictLorentzFactor(velocity);
}

Vec3 momentum(const Vec3& velocity, double mass) { return velocity * (gamma(velocity) * mass); }

double kineticEnergy(const Vec3& velocity,double mass) {
    return two_body::kineticEnergyFromVelocity(velocity,mass);
}

Vec3 velocityFromMomentum(const Vec3& momentum, double mass) {
    const double gammaFromMomentum = std::sqrt(1.0 + momentum.squaredNorm() / (mass*mass*c*c));
    return momentum / (gammaFromMomentum * mass);
}

void synchronizeCovariantDipoles(State& state) {
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
