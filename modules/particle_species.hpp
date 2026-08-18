#pragma once

#include "physical_constants.hpp"

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
    "electron",   9.1093837139e-31,  -elementaryCharge, 2.00231930436256};
inline constexpr ParticleSpecies positron{
    "positron",   9.1093837139e-31,  +elementaryCharge, 2.00231930436256};
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
constexpr double pairBohrRadius(const ParticlePair& pair) {
    return hbar*hbar/(reducedMassOf(pair)*coulombConstant
        *magnitude(chargeProduct(pair)));
}

// Binding energy of the orbit carrying L = hbar, k|q1 q2|/(2 a).  6.80 eV for
// positronium, 1.41 keV for true muonium, 12.5 keV for protonium.
constexpr double pairBindingEnergy(const ParticlePair& pair) {
    return coulombConstant*magnitude(chargeProduct(pair))
        /(2.0*pairBohrRadius(pair));
}

// The pair this build integrates unless something selects another one.
inline constexpr ParticlePair defaultPair{electron,positron};

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
static_assert(agreesTo(pairBohrRadius(defaultPair),positroniumBohrRadius,1.0e-8));

} // namespace positronium::parameters
