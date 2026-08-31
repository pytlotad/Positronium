#pragma once

// The covariant electric/magnetic dipole tensor's Lorentz boost between the
// particle's rest frame and the lab frame -- what
// modules/relativistic_kinematics.hpp's synchronizeCovariantDipoles calls
// to turn a rest-frame moment into the lab-frame firstDipole/secondDipole
// (magnetic) and firstElectricDipole/secondElectricDipole (electric) pair.
//
// Extracted verbatim from positronium.cpp (Stage 0 of splitting engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied, so it depends on that namespace
// already having in scope: Vec3, DipoleTensor, the c constant, gamma(),
// cross(), dot().  Not yet a standalone, order-independent header.

DipoleTensor lorentzBoostDipole(const DipoleTensor& dipole,
                                const Vec3& frameVelocity) {
    const double speedSquared=frameVelocity.squaredNorm();
    if(speedSquared==0.0) return dipole;
    const double boostGamma=gamma(frameVelocity);
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
