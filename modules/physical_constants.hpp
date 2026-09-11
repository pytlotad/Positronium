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
// Measured electron g-factor (CODATA 2022, |g_e-|=2.00231930436092(36)),
// used in place of the classical point-dipole value g=1 for the Thomas-BMT
// precession dynamics.  By CPT the positron g-factor has the same magnitude.
// Audit finding: this used to be the CODATA 2018 value (...436256), still
// correctly cited as such right here, while particle_species.hpp claimed
// "g-factors from the same source [CODATA 2022]" for this exact number --
// a stale citation, not a wrong constant (both were genuine measured
// values; the 2022 adjustment folded in the 2023 Fan/Gabrielse
// remeasurement).  Verified against physics.nist.gov directly.  The shift
// is 1.6e-12 absolute, 8.2e-13 relative -- far below every tolerance this
// project checks against (1e-8 and up), so no simulated result moves.
inline constexpr double electronGFactor = 2.00231930436092;
// Magnitude of the electron/positron intrinsic magnetic dipole moment,
// mu = (g/2) * mu_B evaluated at the measured g above (~9.28476469e-24 J/T),
// used instead of quoting the bare Bohr magneton as the dipole strength.
inline constexpr double electronMagneticMoment = 0.5*electronGFactor*bohrMagneton;
inline constexpr double coulombConstant = 1.0/(4.0*pi*epsilon0);
inline constexpr double c = speedOfLight;
inline constexpr double eCharge = elementaryCharge;
inline constexpr double coulomb = coulombConstant;
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
// It used to be the literal 4.722121244442525e-14 sitting in this file.  Apart
// from being tied to e+e-, that value bounded only w/r^3; the production
// curl(A) has a larger transverse maximum.  particle_species.hpp now derives
// the pair scale from that actual field maximum.
// Raised from 6 to 12.  The regulator's VALUE at the radii trajectories
// actually visit was already negligible at 6 (w = 0.998 at the Compton
// barrier), but its GRADIENT was not, and the azimuth-averaged dipole energy
// depends on exactly that gradient:
//
//     <U> = -(mu0/4pi) (mu1.mu2) (2/3) w'(r)/r^2,
//
// which is identically zero for an unregularized 1/r^3 and nonzero only
// through w'.  At exponent 6 that left a -58.6 eV channel-signed systematic at
// the barrier -- 0.79% of Coulomb -- that scaled with the arbitrary regulator
// radius (factor 57 and a sign flip under a factor-2 change), i.e. an artefact
// masquerading as the model's only channel-dependent energy effect.  At 12 the
// same term is -2.57 eV, a 23-fold reduction, and the regulator is LESS
// intrusive in value as well (w = 0.999957 against 0.998).  See the README
// section on the dipole-dipole channel difference.
inline constexpr double magneticRegularizationExponent = 12.0;
inline constexpr double hbar = 1.054571817e-34;
// Compton barrier of the e+e- pair: r* = (g/2)*(reduced Compton wavelength)/2
//   = g*hbar/(4*m_e*c) = 193.30 fm.
//
// NOT AN INDEPENDENT IMPORT.  The same number is the electron's magnetic
// moment divided by charge times c:
//
//     mu/(e c) = [(g/2)(e hbar/2m)]/(e c) = g hbar/(4 m c) = r*,
//
// to within 2.6e-16 of magneticMoment(electron), which particle_species.hpp
// pins with a static_assert (it can include this header; this one cannot
// include it back).  The residue is one or two ULPs from the order of the
// multiplications, not a discrepancy: both sides are g*hbar/(4*m*c)
// algebraically.  The expression is left in the form below so the constant
// stays bit-identical to every result already measured with it.
//
// So r* is the LENGTH THE MOMENT CONFIGURATION ITSELF CARRIES -- the only
// length that can be built from mu, e and c -- and not a separately chosen
// stopping radius.
//
// The localization reading is the second, independent one: the separation at
// which localizing the electron requires momenta of order m_e*c, i.e. where
// classical point-particle electrodynamics stops applying.  Two unrelated
// arguments landing on one number is why this is quoted as derived rather
// than as a model parameter.
// r*/r_e = (g/2)/(2*alpha) = 68.60: this is not a regularization choice or
// anything tunable, it follows from the fine structure constant.  Annihilation
// happens 68.6x deeper, at r_e (was misstated as "137x" here -- that number
// is the reduced Compton wavelength's own ratio to r_e, 1/alpha, not r*'s;
// r* is defined above as half that wavelength times g/2, so its own ratio to
// r_e picks up the extra factor of ~2 that brings 137 down to 68.6).  Used
// in crem_collapse.hpp as the depth
// Experiments 1/2 actually resolve the pair's inspiral down to, replacing the
// unrelated chargeCloudRestRadius-derived collisionBoundaryRadius that used
// to bound that measurement -- see README's floor-lowering measurement for
// why: once the retarded-history/step coupling bug was fixed, the mechanical
// integrator turned out to reach this scale reliably, so there was no longer
// a reason to stop three orders of magnitude short of it.
inline constexpr double comptonBarrierRadius =
    electronGFactor*hbar/(4.0*electronMass*speedOfLight);

// The fine-structure constant, DERIVED from the constants already entered
// rather than tabulated, so it cannot drift away from the coupling the
// trajectories actually integrate: alpha = k e^2/(hbar c) = 1/137.035999000.
//
// It is named here because it is the model's ONE dimensionless scale input,
// and naming it makes that checkable.  Every length in the model is the
// reduced Compton wavelength times a power of it -- a_pair = 2 lambda_C/alpha,
// r_e = alpha lambda_C, r* = (g/2) lambda_C/2 -- and particle_species.hpp
// pins each of those with a static_assert.  See the README section on the
// alpha gap for what that does and does not buy.
inline constexpr double fineStructureConstant =
    coulombConstant*elementaryCharge*elementaryCharge/(hbar*speedOfLight);
// hbar/(m_e c), the reduced Compton wavelength.  The one dimensionful scale
// the moment configuration itself carries (see comptonBarrierRadius).
inline constexpr double reducedComptonWavelength =
    hbar/(electronMass*speedOfLight);
// e^2/(4 pi eps0 m c^2), the length at which the electron's electrostatic
// self-energy equals its rest energy.  Not a model parameter: it is quoted
// here only as the second contact scale the first-passage measurement in
// crem_collapse.hpp compares L against, sitting 1/(2 alpha) = 68.5 below
// comptonBarrierRadius above -- see the README's "annihilation happens 137
// times deeper" note, which is the same statement.
inline constexpr double classicalElectronRadius =
    coulombConstant*elementaryCharge*elementaryCharge
    /(electronMass*speedOfLight*speedOfLight);
inline constexpr double chargeCloudRestRadius = 0.01*bohrRadius;
// The separation at which a trajectory counts as collided is NOT here: it
// scales with the pair, so it is derived in particle_species.hpp as
// collisionBoundaryOf(pair) and exported there as collisionBoundaryRadius.
//
// It used to sit here as 0.01*a0 = 529 fm, hydrogen's Bohr radius, while the
// starting separation has long been the PAIR's Bohr radius.  For mu+mu- that
// start is 512 fm and for p+pbar 57.6 fm, i.e. already inside the boundary:
// the trajectory was declared collided before its first step and came back as
// a numerical failure.

} // namespace positronium::parameters
