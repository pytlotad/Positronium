#pragma once

// The particle pair the engine is configured to simulate: the role-scoped
// globals every force law reads (firstMass/firstCharge/... for "first",
// mirrored for "second", never "electron"/"positron" -- see the comment
// on activePair below for why), the pair-level Coulomb/dipole/reduced-mass
// quantities derived from them, and applyPair/applyPairFromOption, which is
// how --pair changes any of it at startup.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied (right after the namespace opens
// and its own using-declarations/c/eCharge/coulomb constants, which stay in
// positronium.cpp since a namespace-opening brace has no business living
// separately from its close), so it depends on that namespace already
// having in scope: Vec3, ParticlePair, ParticleSpecies, defaultPair,
// magnitude(), isAttracting(), chargeProduct(), speciesByName(),
// selectableSpeciesList(), magneticMoment(), dipoleRegularizationRadius(),
// collisionBoundaryOf(), magneticRegularizationRadius,
// collisionBoundaryRadius, the c/eCharge/coulomb constants.  Not yet a
// standalone, order-independent header.

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
ParticlePair activePair=defaultPair;
ParticleSpecies firstSpecies=activePair.first;
ParticleSpecies secondSpecies=activePair.second;
double firstMass=firstSpecies.mass;
double secondMass=secondSpecies.mass;
double firstCharge=firstSpecies.charge;
double secondCharge=secondSpecies.charge;
double firstGFactor=firstSpecies.gFactor;
double secondGFactor=secondSpecies.gFactor;

// Gyromagnetic ratio relating the intrinsic moment to the intrinsic spin,
// mu = gamma S with gamma = g q / (2 m).  The g belongs here and was missing
// from four call sites that wrote q/(2m) instead.
//
// The model carries |mu| = (g/2) * magneton, so the correct ratio returns
// S = hbar/2 exactly, as a spin-1/2 particle must.  Dropping g returned
// 1.00116 hbar, i.e. the reported intrinsic angular momentum was inflated by
// g = 2.0023, and the dipole's response to a radiation-reaction torque,
// d(mu) = gamma tau dt, was weaker by the same factor.
double gyromagneticRatio(double charge,double mass,double gFactor) {
    return gFactor*charge/(2.0*mass);
}
double firstGyromagneticRatioOf() {
    return gyromagneticRatio(firstCharge,firstMass,firstGFactor);
}
double secondGyromagneticRatioOf() {
    return gyromagneticRatio(secondCharge,secondMass,secondGFactor);
}
// Magnitudes differ between roles for every pair except a particle and its
// own antiparticle: the proton carries 1.41e-26 J/T against the electron's
// 9.28e-24, four hundred times smaller.  Code that reached for one role's
// moment while dressing the other was correct only by that accident.
double firstMagneticMoment=magneticMoment(firstSpecies);
double secondMagneticMoment=magneticMoment(secondSpecies);

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
double pairChargeProduct=firstCharge*secondCharge;
double pairCoulombStrength=coulomb*magnitude(pairChargeProduct);
double pairDipoleCharge=
    (firstCharge*secondMass-secondCharge*firstMass)/(firstMass+secondMass);
double pairReducedMass=firstMass*secondMass/(firstMass+secondMass);

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
void applyPair(const ParticlePair& pair) {
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
void applyPairFromOption(const std::string& value) {
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
Vec3 gExternalMagneticField;

// Earth's field is about this, and it is the value the startup question
// offers.  Worth knowing before reading the output: at 50 uT the cyclotron
// rate eB/m is 8.8e6 rad/s against an orbital rate near 3e15 rad/s, so the
// orbit itself is untouched at the ninth decimal.  What the field does reach
// is the dipoles, which precess at the same 8.8e6 rad/s -- about 3e-4 rad over
// a 35 ps collapse.  The effect is real, small, and mostly magnetic.
inline constexpr double earthScaleMagneticField=50.0e-6;
