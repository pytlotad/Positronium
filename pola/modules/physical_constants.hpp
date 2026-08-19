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
// Smoothing radius of the magnetic dipole field, fixed by a physical ceiling
// rather than by the point-particle boundary:
//
//     a = cbrt( (mu0/4pi) * mu^2 / (m_e c^2) ) = 47.22 fm
//
// is the separation at which the classical dipole-dipole interaction energy
// reaches the electron rest energy.  Because the regularized profile
// w(r)/r^3 = r^3/(r^6+a^6) peaks exactly at r=a, this choice caps the classical
// dipole energy at m_e c^2 / 2 everywhere.  A classical model has no business
// producing interaction energies above the pair-creation threshold.
//
// mu is electronMagneticMoment = (g/2)*mu_B, the moment the dynamics actually
// carries -- not the bare Bohr magneton.  Deriving the radius from mu_B while
// the dipoles carry (g/2)*mu_B left the true peak at 0.5012*m_e c^2, i.e. just
// above the ceiling this constant exists to enforce.  The (g/2)^(2/3) = 1.00077
// correction is numerically inert for every reported orbit but restores the
// guarantee exactly; positronium_validation measures the peak with the same
// moment, so the check now certifies the production configuration.
//
// The previous value, 0.5*nuclearCutoff = 5 fm, was chosen only so the
// smoothing stayed inside the reported domain.  It let the dipole energy reach
// 2.1e8 eV, i.e. 420 electron rest masses, and it left the magnetic sector 106
// times more point-like than the charge sector, which already assumes a finite
// source of radius chargeCloudRestRadius.  Raising the radius to 47 fm leaves
// every reported orbit bit-identical: at tens of picometres the weight differs
// from one by less than 1e-19.
//
// This does NOT remove the short-range dipole barrier.  That barrier sits where
// the unregularized dipole and Coulomb energies cross,
// r* = sqrt((mu0/4pi)*mu^2/(k e^2)) = 193 fm, which is half the reduced
// Compton wavelength and far outside any smoothing radius that leaves the
// reported orbits untouched.
inline constexpr double magneticRegularizationRadius = 4.722121244442525e-14;
inline constexpr double magneticRegularizationExponent = 6.0;
inline constexpr double hbar = 1.054571817e-34;
inline constexpr double chargeCloudRestRadius = 0.01*bohrRadius;

} // namespace positronium::parameters
