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
// Numerical boundary of the point-particle model.  Keeping the magnetic
// smoothing transition inside this boundary prevents it from modifying any
// reported trajectory while still protecting trial evaluations near r=0.
inline constexpr double magneticRegularizationRadius = 0.5*nuclearCutoff;
inline constexpr double magneticRegularizationExponent = 6.0;
inline constexpr double hbar = 1.054571817e-34;
inline constexpr double chargeCloudRestRadius = 0.01*bohrRadius;

} // namespace positronium::parameters
