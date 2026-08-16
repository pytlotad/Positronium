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
inline constexpr double coulombConstant = 1.0/(4.0*pi*epsilon0);
inline constexpr double bohrRadius = 5.29177210903e-11;
inline constexpr double nuclearCutoff = 1.0e-14;
// Smoothing radius of the magnetic dipole field, fixed by a physical ceiling
// rather than by the point-particle boundary:
//
//     a = cbrt( (mu0/4pi) * mu_B^2 / (m_e c^2) ) = 47.18 fm
//
// is the separation at which the classical dipole-dipole interaction energy
// reaches the electron rest energy.  Because the regularized profile
// w(r)/r^3 = r^3/(r^6+a^6) peaks exactly at r=a, this choice caps the classical
// dipole energy at m_e c^2 / 2 everywhere.  A classical model has no business
// producing interaction energies above the pair-creation threshold.
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
// r* = sqrt((mu0/4pi)*mu_B^2/(k e^2)) = 193 fm, which is half the reduced
// Compton wavelength and far outside any smoothing radius that leaves the
// reported orbits untouched.
inline constexpr double magneticRegularizationRadius = 4.718474089940226e-14;
inline constexpr double magneticRegularizationExponent = 6.0;
inline constexpr double hbar = 1.054571817e-34;
inline constexpr double chargeCloudRestRadius = 0.01*bohrRadius;

} // namespace positronium::parameters
