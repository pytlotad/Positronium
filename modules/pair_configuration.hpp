#pragma once

// The particle pair the engine is configured to simulate: the role-scoped
// globals every force law reads (firstMass/firstCharge/... for "first",
// mirrored for "second", never "electron"/"positron" -- see the comment
// on activePair below for why), the pair-level Coulomb/dipole/reduced-mass
// quantities derived from them, and applyPair/applyPairFromOption, which is
// how --pair changes any of it at startup.
//
// Self-contained and order-independent.  It names what it needs through a
// using-directive on positronium::parameters -- the same one positronium.cpp
// carries, and the readable choice here because nearly every line below
// reaches into that namespace -- rather than reopening namespace positronium.
// That distinction matters: the header is still textually included inside
// positronium.cpp's anonymous namespace, where reopening a named namespace
// would create {anonymous}::positronium and hide the real one from every
// later lookup.
//
// Every global below is `inline`, so including this header from more than one
// translation unit gives the linker one definition rather than several.

#include "particle_species.hpp"
#include "physical_constants.hpp"
#include "vector3.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using positronium::objects::Vec3;
using namespace positronium::parameters;

// The two ROLES the dynamics integrates, as opposed to the species that fill
// them.  Every force law below refers to "first" and "second"; which particles
// those are is decided here and nowhere else.  Naming them electron/positron
// was accurate only for the default pair and actively misleading for any
// other, since the state fields are role slots, not species.
//
// These are variables, not constants, because --pair chooses the pair at
// startup.  Dropping constexpr costs nothing measurable: on a 30-event beam
// run it moved the wall clock by 0.07%, inside the 0.6% spread between
// repeats, because the inner loop is dominated by the retarded-history search
// and the field solver rather than by folding these scalars.  They are
// constant-initialized to the default pair, so a run that never passes --pair
// is bit-identical to the old build.
inline ParticlePair activePair=defaultPair;
inline ParticleSpecies firstSpecies=activePair.first;
inline ParticleSpecies secondSpecies=activePair.second;
inline double firstMass=firstSpecies.mass;
inline double secondMass=secondSpecies.mass;
inline double firstCharge=firstSpecies.charge;
inline double secondCharge=secondSpecies.charge;
inline double firstGFactor=firstSpecies.gFactor;
inline double secondGFactor=secondSpecies.gFactor;

// Gyromagnetic ratio relating the intrinsic moment to the intrinsic spin,
// mu = gamma S with gamma = g q / (2 m).  The g belongs here and was missing
// from four call sites that wrote q/(2m) instead.
//
// The model carries |mu| = (g/2) * magneton, so the correct ratio returns
// S = hbar/2 exactly, as a spin-1/2 particle must.  Dropping g returned
// 1.00116 hbar, i.e. the reported intrinsic angular momentum was inflated by
// g = 2.0023, and the dipole's response to a radiation-reaction torque,
// d(mu) = gamma tau dt, was weaker by the same factor.
inline double gyromagneticRatio(double charge,double mass,double gFactor) {
    return gFactor*charge/(2.0*mass);
}
inline double firstGyromagneticRatioOf() {
    return gyromagneticRatio(firstCharge,firstMass,firstGFactor);
}
inline double secondGyromagneticRatioOf() {
    return gyromagneticRatio(secondCharge,secondMass,secondGFactor);
}
// Magnitudes differ between roles for every pair except a particle and its
// own antiparticle: the proton carries 1.41e-26 J/T against the electron's
// 9.28e-24, four hundred times smaller.  Code that reached for one role's
// moment while dressing the other was correct only by that accident.
inline double firstMagneticMoment=magneticMoment(firstSpecies);
inline double secondMagneticMoment=magneticMoment(secondSpecies);

// Pair-level quantities.  These three replace the e^2 that used to be written
// out wherever two charges met, and they are NOT interchangeable:
//
//   pairChargeProduct   q1 q2, SIGNED.  Negative for an attracting pair, which
//                       is the only case with bound states.  Belongs in the
//                       Coulomb potential and force, where the sign is the
//                       physics.
//   pairCoulombStrength k|q1 q2|, a magnitude.  Belongs in the Kepler
//                       relations (period, semi-major axis, attraction
//                       parameter), which presuppose attraction.
//   pairDipoleCharge    the effective charge of the pair's electric dipole.
//                       With the centre of mass at rest, m1 r1 + m2 r2 = 0, so
//                       d = q1 r1 + q2 r2 = [(q1 m2 - q2 m1)/(m1+m2)] r.  For
//                       e+e- that collapses to -e and |d| = e|r|, which is why
//                       the old code could write e and be right; for a pair
//                       with unequal masses it does not.
inline double pairChargeProduct=firstCharge*secondCharge;
inline double pairCoulombStrength=coulomb*magnitude(pairChargeProduct);
inline double pairDipoleCharge=
    (firstCharge*secondMass-secondCharge*firstMass)/(firstMass+secondMass);
inline double pairReducedMass=firstMass*secondMass/(firstMass+secondMass);

// Point every role constant at a different pair.  Called once, from argument
// parsing, before anything integrates.
//
// The two conditions below used to be static_asserts on the hard-coded pair.
// As runtime checks they are worth strictly more: they now guard a value the
// USER supplies rather than one a developer edited into the header, and they
// are the reason --pair cannot be used to ask for a system this model has no
// business integrating.  A repelling pair has no bound states at all, so
// every experiment built around capture and inspiral would be meaningless
// rather than merely inaccurate.
inline void applyPair(const ParticlePair& pair) {
    if(!isAttracting(pair))
        throw std::invalid_argument(std::string("pair ")+pair.first.name+"+"
            +pair.second.name+" repels; it has no bound states");
    if(magnitude(magnitude(chargeProduct(pair))-eCharge*eCharge)
        >1.0e-12*eCharge*eCharge)
        throw std::invalid_argument(std::string("pair ")+pair.first.name+"+"
            +pair.second.name+" does not carry opposite unit charges");
    activePair=pair;
    firstSpecies=pair.first;
    secondSpecies=pair.second;
    firstMass=firstSpecies.mass;
    secondMass=secondSpecies.mass;
    firstCharge=firstSpecies.charge;
    secondCharge=secondSpecies.charge;
    firstGFactor=firstSpecies.gFactor;
    secondGFactor=secondSpecies.gFactor;
    firstMagneticMoment=magneticMoment(firstSpecies);
    secondMagneticMoment=magneticMoment(secondSpecies);
    pairChargeProduct=firstCharge*secondCharge;
    pairCoulombStrength=coulomb*magnitude(pairChargeProduct);
    pairDipoleCharge=
        (firstCharge*secondMass-secondCharge*firstMass)/(firstMass+secondMass);
    pairReducedMass=firstMass*secondMass/(firstMass+secondMass);
    // Derived in particle_species.hpp from BOTH moments and the lighter mass,
    // so it has to follow the pair rather than stay at the default's 68.47 fm.
    magneticRegularizationRadius=dipoleRegularizationRadius(pair);
    // Same reason: the boundary is a fraction of the pair's own orbit, and a
    // heavy pair starts inside hydrogen's 529 fm.
    collisionBoundaryRadius=collisionBoundaryOf(pair);
}

// Parse "first,second" as the command line spells it.
inline void applyPairFromOption(const std::string& value) {
    const std::size_t comma=value.find(',');
    if(comma==std::string::npos)
        throw std::invalid_argument(
            "--pair needs two species separated by a comma, e.g. proton,electron");
    const ParticleSpecies* first=speciesByName(value.substr(0,comma));
    const ParticleSpecies* second=speciesByName(value.substr(comma+1));
    if(!first||!second)
        throw std::invalid_argument("unknown species in --pair; known species: "
            +selectableSpeciesList());
    applyPair({*first,*second});
}
// Uniform external magnetic field, zero unless the run asks for one.  Zero is
// the historical behaviour and the default: the model's own list of excluded
// effects named external fields, and every committed result was produced
// without one.
//
// It enters in exactly three places, which are the three the pair can feel it
// through: the instantaneous force sum, the retarded force sum, and the local
// field each particle sees, the last of which carries it into Thomas-BMT
// precession for both roles.
inline Vec3 gExternalMagneticField;

// Earth's field is about this, and it is the value the startup question
// offers.  Worth knowing before reading the output: at 50 uT the cyclotron
// rate eB/m is 8.8e6 rad/s against an orbital rate near 3e15 rad/s, so the
// orbit itself is untouched at the ninth decimal.  What the field does reach
// is the dipoles, which precess at the same 8.8e6 rad/s -- about 3e-4 rad over
// a 35 ps collapse.  The effect is real, small, and mostly magnetic.
inline constexpr double earthScaleMagneticField=50.0e-6;
