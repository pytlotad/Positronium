#pragma once

#include "physical_constants.hpp"

#include <array>
#include <string>
#include <string_view>

namespace positronium::parameters {

// A charged spin-1/2 particle, given by the three properties that fix its
// electrodynamics: inertia, coupling to the field, and the size of its
// intrinsic magnetic moment relative to its own magneton.
//
// The moment is DERIVED rather than tabulated.  Given the g-factor,
//
//     magneton  = |q| hbar / (2 m),      mu = (g/2) * magneton,
//
// which reproduces the measured moments of all six species below to better
// than 2e-9 -- so the three numbers here carry the whole magnetic sector and
// cannot drift out of step with each other.  Tabulating mu separately would
// allow exactly that kind of inconsistency, which this project has already
// been bitten by once (a g-factor of 2.0023 sitting beside a moment of mu_B).
//
// The g-factor is emphatically NOT close to 2 for every particle: the proton
// carries g = 5.5857 because it is composite.  Nothing here may assume the
// Dirac value.
struct ParticleSpecies {
    const char* name;
    double mass;    // kg
    double charge;  // C, signed
    double gFactor; // dimensionless
};

constexpr double magnitude(double value) { return value<0.0?-value:value; }

// |q| hbar / (2 m).  Bohr magneton for the electron, nuclear magneton for the
// proton, and the corresponding scale for every other species.
constexpr double magneton(const ParticleSpecies& species) {
    return magnitude(species.charge)*hbar/(2.0*species.mass);
}

// Magnitude of the intrinsic magnetic moment.  The direction relative to spin
// follows the sign of the charge and is carried by the gyromagnetic ratio in
// the precession equations, not here.
constexpr double magneticMoment(const ParticleSpecies& species) {
    return 0.5*species.gFactor*magneton(species);
}

// CODATA 2022 masses; g-factors from the same source.  Antiparticles take the
// particle's mass and g-factor with the charge reversed, which is CPT rather
// than an independent measurement -- for the electron/positron and the
// proton/antiproton pairs that equality is itself experimentally tested to
// high precision, and the model has no mechanism that could distinguish them.
inline constexpr ParticleSpecies electron{
    "electron",   9.1093837139e-31,  -elementaryCharge, 2.00231930436092};
inline constexpr ParticleSpecies positron{
    "positron",   9.1093837139e-31,  +elementaryCharge, 2.00231930436092};
inline constexpr ParticleSpecies muon{
    "muon",       1.883531627e-28,   -elementaryCharge, 2.00233184123};
inline constexpr ParticleSpecies antimuon{
    "antimuon",   1.883531627e-28,   +elementaryCharge, 2.00233184123};
inline constexpr ParticleSpecies proton{
    "proton",     1.67262192595e-27, +elementaryCharge, 5.5856946893};
inline constexpr ParticleSpecies antiproton{
    "antiproton", 1.67262192595e-27, -elementaryCharge, 5.5856946893};

// The two particles a run actually integrates, plus the scales their pairing
// implies.  Everything the dynamics needs about "which particles are these"
// is meant to come from here.
struct ParticlePair {
    ParticleSpecies first;
    ParticleSpecies second;
};

constexpr double reducedMassOf(const ParticlePair& pair) {
    return pair.first.mass*pair.second.mass
        /(pair.first.mass+pair.second.mass);
}

// Product of the charges.  Negative for an attracting pair, which is the only
// case that has bound states.
constexpr double chargeProduct(const ParticlePair& pair) {
    return pair.first.charge*pair.second.charge;
}

constexpr bool isAttracting(const ParticlePair& pair) {
    return chargeProduct(pair)<0.0;
}

// Bohr radius of the pair, hbar^2/(mu k |q1 q2|).  For e+e- this is 2*a0
// (1.058 angstrom) because the reduced mass is m_e/2; for mu+mu- it is 512 fm
// and for p+pbar 57.6 fm.  a0 itself is hydrogen's and belongs to none of
// them.
// Derived from the pair's MAGNETIC MOMENTS rather than posted independently,
// so hbar enters this model in exactly one place -- the tabulated moments --
// instead of two.
//
// The model's electrodynamics already produces one length by itself: the
// radius where the magnetic and Coulomb energies cross,
//
//     r*^2 = (mu0/4pi) mu1 mu2 / K = mu1 mu2 / (c^2 |q1 q2|),
//
// which is the same number as comptonBarrierRadius = g hbar/(4 m c) -- not a
// coincidence but an identity, since mu = (g/2)(q hbar/2m) gives
// r* = mu/(|q| c) directly.  Verified: 193.3035387648 fm against
// 193.3035388174 fm, a relative 2.7e-10 that is CODATA rounding between the
// tabulated moment and hbar, nothing more.
//
// Eliminating hbar between r* and the Bohr relation leaves
//
//     a_pair = 16 (m1+m2) mu1 mu2 / (g1 g2 |q1 q2| K),
//
// with no explicit hbar and no explicit c.  It agrees with hbar^2/(mu_red K)
// to 5.4e-10, again CODATA rounding.
//
// WHAT THIS IS AND IS NOT.  It is not a classical derivation of the Bohr
// radius: mu carries hbar, so this eliminates hbar between two quantities
// that both contain it.  What it does establish is that a_pair and the
// magnetic moment are NOT independent inputs -- the binding energy this model
// reports follows from the measured moment through classical electrodynamics,
// rather than being a second quantum constant entered alongside it.
constexpr double pairBohrRadius(const ParticlePair& pair) {
    return 16.0*(pair.first.mass+pair.second.mass)
        *magneticMoment(pair.first)*magneticMoment(pair.second)
        /(pair.first.gFactor*pair.second.gFactor
          *magnitude(chargeProduct(pair))
          *coulombConstant*magnitude(chargeProduct(pair)));
}

// Binding energy of the orbit carrying L = hbar, k|q1 q2|/(2 a).  6.80 eV for
// positronium, 1.41 keV for true muonium, 12.5 keV for protonium.
constexpr double pairBindingEnergy(const ParticlePair& pair) {
    return coulombConstant*magnitude(chargeProduct(pair))
        /(2.0*pairBohrRadius(pair));
}

// The pair this build integrates unless something selects another one.
inline constexpr ParticlePair defaultPair{electron,positron};

// Cube root usable in a constant expression, which std::cbrt is not.  The
// argument is first scaled into [1,8) by powers of eight -- exact in binary
// floating point, so the scaling costs no precision -- and Newton's iteration
// x <- (2x + v/x^2)/3 then converges quadratically from that seed.  Sixty
// passes are far more than double precision can use; the count is chosen to
// need no convergence argument at all.
constexpr double cubeRoot(double value) {
    if(!(value>0.0)) return 0.0;
    double scaled=value,factor=1.0;
    while(scaled>8.0) { scaled/=8.0; factor*=2.0; }
    while(scaled<1.0) { scaled*=8.0; factor/=2.0; }
    double x=scaled;
    for(int step=0;step<60;++step) x=(2.0*x+scaled/(x*x))/3.0;
    return x*factor;
}

// Smoothing radius of the magnetic dipole field, set by a physical ceiling
// rather than by the point-particle boundary.  The regulator is applied to
// A=C f(r)(mu x r), not directly to B, so its energy ceiling must be derived
// from curl(A).  With x=r/a and p=6 the two eigenvalues of B/(C mu/a^3) are
//
//     longitudinal:  2 x^3/(1+x^6),
//     transverse:     x^3(5-x^6)/(1+x^6)^2.
//
// The absolute transverse eigenvalue is the larger one.  It peaks at
// x^6=9-2 sqrt(19) with the dimensionless value below.  Therefore
//
//     a^3 = peak * (mu0/4pi) mu1 mu2 / (E/2)
//
// caps the ACTUAL energy -mu1.B_reg at E/2 for every separation and
// orientation when both moments have their tabulated proper magnitudes in a
// common rest frame.  Omitting peak and the factor two capped only f(r), while
// the production curl(A) still reached 2.778 E.
//
// E is the rest energy of the LIGHTER of the two, which is the first scale at
// which the classical description of a constituent stops meaning anything.  A
// classical model has no business producing interaction energies above it.
//
// Both moments enter, mu1 mu2 rather than mu^2.  That is the same number only
// for a particle and its own antiparticle; the proton's moment is 658 times
// smaller than the electron's, so the radii remain strongly pair-dependent.
// Squaring one role's moment would have overestimated the singular energy of
// every mixed pair and smoothed it far more than its own physics asks for.
//
// For e+e- this does NOT remove the short-range dipole barrier: the
// unregularized dipole and Coulomb energies cross at
// r* = sqrt((mu0/4pi) mu1 mu2 / k|q1 q2|) = 193 fm, outside the 83.64 fm
// smoothing radius.  For sufficiently asymmetric pairs the E/2 scale can
// overlap that formal point-dipole crossover, but both then lie at or below
// the 10 fm terminal boundary where the reported trajectory already ends.
inline constexpr double dipoleEnergyCeilingFraction=0.5;
// Peak of |U| in units of (mu0/4pi) mu1 mu2 / r_reg^3, for the CURRENT
// regulator exponent.  Derived rather than fitted: the transverse coefficient
// is w[n(1-w)-1]/r^3, and with u = (r_reg/r)^n its extremum solves
//
//     (n-1)(3-n) u^2 + (n^2+4n-6) u - 3 = 0,
//
// which at n=6 gives 5u^2-18u+1=0, u=(9+2sqrt19)/5, hence the (9-2sqrt19)^(1/6)
// peak radius this file used to carry, and at n=12 gives 33u^2-62u+1=0,
// u=(31+4sqrt58)/33.  Evaluating the coefficient there yields the constant
// below; an independent 4e6-point numerical maximization agrees to 1.2e-13.
inline constexpr double regulatorDipoleCurlPeak=2.7783644145278381;
static_assert(magneticRegularizationExponent==12.0,
    "regulatorDipoleCurlPeak must be re-derived when the profile changes");

constexpr double dipoleRegularizationRadius(const ParticlePair& pair) {
    const double lighterMass=pair.first.mass<pair.second.mass
        ?pair.first.mass:pair.second.mass;
    return cubeRoot((mu0/(4.0*pi))
        *magneticMoment(pair.first)*magneticMoment(pair.second)
        *regulatorDipoleCurlPeak
        /(dipoleEnergyCeilingFraction
          *lighterMass*speedOfLight*speedOfLight));
}

// The name the field code has always used, now answering for whichever pair
// the RUN integrates.  Deliberately not constexpr: --pair chooses the pair at
// startup and applyPair reassigns this before anything integrates.  It is
// constant-initialized to the default pair's derived value.
inline double magneticRegularizationRadius=
    dipoleRegularizationRadius(defaultPair);

// True only for e-e+ in either role order.  Not merely "a particle with its
// own antiparticle": the annihilation data the bound-decay experiments compare
// against is positronium's measured lifetimes, and true muonium and protonium
// annihilate through different physics that this model does not carry.
constexpr bool isPositronium(const ParticlePair& pair) {
    return pair.first.mass==electron.mass
        && pair.second.mass==electron.mass
        && pair.first.charge==-pair.second.charge
        && magnitude(pair.first.charge)==elementaryCharge;
}

// Separation at which a trajectory is declared to have collided, and the depth
// at which the collapse estimator stops integrating.
//
// Expressed as a fraction of the PAIR's own Bohr radius, so the experiments
// stay geometrically similar from one pair to the next.  The coefficient is
// chosen to reproduce the historical 0.01*a0 = 529 fm for positronium, whose
// pair radius is 2*a0: 0.005 * 2 a0 = 0.01 a0.  Every other pair now gets the
// same fraction of its own orbit instead of hydrogen's absolute number --
// 2.56 fm for mu+mu-, 0.29 fm for p+pbar, 265 fm for p+e-.
//
// The default pair shifts by 2.6e-9 relative, because pairBohrRadius derives
// the radius from hbar^2/(mu k q^2) while the old constant came from the
// tabulated a0; that is the same rounding gap the assertions below admit.
constexpr double collisionBoundaryOf(const ParticlePair& pair) {
    return 0.005*pairBohrRadius(pair);
}

// The name the trajectory and collapse code use, following whichever pair the
// run integrates.  Not constexpr for the same reason as the regularization
// radius above: --pair reassigns it before anything integrates.
inline double collisionBoundaryRadius=collisionBoundaryOf(defaultPair);

// Separation at which classical point-particle electrodynamics stops
// applying: comptonBarrierRadius (physical_constants.hpp) for e+e-, since
// that scale is derived from the electron's own Compton wavelength and does
// not generalize to other masses, else the pair's collisionBoundaryOf(pair)
// above. Three independent call sites (BeamConfiguration/exp3-4,
// InteractionConfiguration/exp5, and the bound-state collapse-time
// extrapolation exp5 also does for captured pairs) used to spell this
// ternary out separately; they drifted out of sync once already (exp3/4
// still pointed at nuclearCutoff while exp5 had already moved to
// comptonBarrierRadius, reconciled in 340b67b) before being written here
// once so they cannot drift again.
constexpr double pointParticleBoundaryOf(const ParticlePair& pair) {
    return isPositronium(pair) ? comptonBarrierRadius : collisionBoundaryOf(pair);
}

// Everything --pair is allowed to name.  Kept next to the species themselves
// so a new species cannot be added without becoming selectable.
inline constexpr std::array<const ParticleSpecies*,6> selectableSpecies{
    &electron,&positron,&muon,&antimuon,&proton,&antiproton};

inline const ParticleSpecies* speciesByName(std::string_view name) {
    for(const ParticleSpecies* species:selectableSpecies)
        if(name==species->name) return species;
    return nullptr;
}

inline std::string selectableSpeciesList() {
    std::string result;
    for(const ParticleSpecies* species:selectableSpecies) {
        if(!result.empty()) result+=", ";
        result+=species->name;
    }
    return result;
}

// The species table must agree with the standalone electron constants, or the
// two descriptions of the same particle have drifted apart.
//
// The first three compare identical literals and are exact.  The last two
// cannot be: electronMagneticMoment is built from the TABULATED Bohr magneton
// while magneticMoment(electron) derives it as e*hbar/(2 m_e), and
// positroniumBohrRadius is 2*a0 from the tabulated a0 while pairBohrRadius
// derives it from hbar^2/(mu k e^2).  Independently rounded CODATA values do
// not recombine bit-for-bit; the gaps are 6.1e-10 and 2.6e-9.  The 1e-8 band
// admits exactly that and nothing else -- a wrong mass, charge or g-factor
// would miss by 1e-3 or worse.
constexpr bool agreesTo(double left,double right,double tolerance) {
    return magnitude(left-right)<=tolerance*magnitude(right);
}
static_assert(electron.mass==electronMass);
static_assert(positron.mass==positronMass);
static_assert(electron.gFactor==electronGFactor);
static_assert(magnitude(electron.charge)==elementaryCharge);
static_assert(agreesTo(magneticMoment(electron),electronMagneticMoment,1.0e-8));
// Pinned to e+e- explicitly, not to defaultPair: positroniumBohrRadius is
// positronium's scale and stays the right comparison however the default pair
// is set.  Written against defaultPair this assertion stopped being a
// consistency check and became a ban on ever selecting another pair.
static_assert(agreesTo(pairBohrRadius(ParticlePair{electron,positron}),
                       positroniumBohrRadius,1.0e-8));
// Pin the default scale as a guard against dropping the curl(A) peak factor or
// one of the two role moments.  The former scalar-profile radius was 47.22 fm;
// enforcing the same E/2 ceiling on the production field raises it by
// (2*2.778364...)^(1/3)=1.771560... to 83.64 fm.
static_assert(agreesTo(dipoleRegularizationRadius(ParticlePair{electron,positron}),
                       8.363926379319052e-14,1.0e-8));

} // namespace positronium::parameters
