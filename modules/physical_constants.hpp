#pragma once

namespace positronium::parameters {

inline constexpr double pi = 3.14159265358979323846;
inline constexpr double epsilon0 = 8.8541878128e-12;
inline constexpr double speedOfLight = 299792458.0;
inline constexpr double mu0 = 4.0*pi*1.0e-7;
inline constexpr double elementaryCharge = 1.602176634e-19;
inline constexpr double electronMass = 9.1093837139e-31;
inline constexpr double positronMass = electronMass;
inline constexpr double bohrMagneton = 9.2740100657e-24;
// Measured electron g-factor (CODATA 2018, |g_e-|=2.00231930436256(35)),
// used in place of the classical point-dipole value g=1 for the Thomas-BMT
// precession dynamics.  By CPT the positron g-factor has the same magnitude.
inline constexpr double electronGFactor = 2.00231930436256;
// Magnitude of the electron/positron intrinsic magnetic dipole moment,
// mu = (g/2) * mu_B evaluated at the measured g above (~9.28476469e-24 J/T),
// used instead of quoting the bare Bohr magneton as the dipole strength.
inline constexpr double electronMagneticMoment = 0.5*electronGFactor*bohrMagneton;
inline constexpr double coulombConstant = 1.0/(4.0*pi*epsilon0);
inline constexpr double bohrRadius = 5.29177210903e-11;
// Bohr radius of the POSITRONIUM system.  a0 above is hydrogen's, built from
// the electron mass; positronium's reduced mass is mu = m_e/2, so its scale is
//
//     a_Ps = hbar^2/(mu k e^2) = 2 a0 = 1.058 angstrom.
//
// Equivalently: the orbit carrying L = hbar has semi-major axis a_Ps, and its
// binding energy is k e^2/(2 a_Ps) = 6.8 eV -- the measured Ps binding energy.
// The three statements are one condition seen three ways.
//
// The pair used to be prepared at separation a0 with a sub-circular tangential
// speed, which put the sampled orbits at L = 0.51-0.69 hbar, semi-major axis
// ~0.8 a0 and binding energy ~17 eV: INSIDE the ground state the model is
// meant to represent.  Because the classical inspiral time scales as a^3 that
// alone shortened the collapse by more than an order of magnitude.
//
// L = hbar is the Bohr/SED value, not the quantum 1s value, which is L = 0.
// A classical L = 0 orbit is a radial plunge, so it cannot stand in for a
// ground state at all; L = hbar is the classical orbit that reproduces the
// correct binding energy and mean radius.
inline constexpr double positroniumBohrRadius = 2.0*bohrRadius;
inline constexpr double nuclearCutoff = 1.0e-14;
// The smoothing radius of the magnetic dipole field is NOT here: it depends on
// which two particles are being integrated, so it is derived in
// particle_species.hpp as dipoleRegularizationRadius(pair) and exported there
// under the name magneticRegularizationRadius.  The exponent below is a pure
// shape parameter of the regulator and carries no species scale, so it stays.
//
// It used to be the literal 4.722121244442525e-14 sitting in this file, which
// is the correct value for e+e- and for nothing else -- it was built from the
// electron's moment and the electron's rest energy, so for any other pair the
// ceiling it exists to enforce silently stopped binding.
inline constexpr double magneticRegularizationExponent = 6.0;
inline constexpr double hbar = 1.054571817e-34;
inline constexpr double chargeCloudRestRadius = 0.01*bohrRadius;
// Separation at which a trajectory is declared to have collided, and the depth
// at which the collapse estimator stops integrating.
//
// Numerically equal to chargeCloudRestRadius, and that is the point: until now
// ONE constant played both roles, so the reach of the trajectory sector could
// not be changed without also resizing the charge cloud that the Maxwell grid
// deposits, and vice versa.  They are separate questions -- one is the model's
// spatial resolution, the other is where a trajectory is deemed to have ended
// -- and separating them is what allows the reach to be probed on its own.
inline constexpr double collisionBoundaryRadius = chargeCloudRestRadius;

} // namespace positronium::parameters
