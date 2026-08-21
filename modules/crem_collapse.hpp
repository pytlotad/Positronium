#pragma once

// CREM collapse estimator: the orbit-averaged secular integration that turns a
// prepared bound state into a measured classical inspiral time, together with
// the closed-form electrodynamic references it is compared against.
//
// Textual module, like the other headers here: it is included once, from
// inside the anonymous namespace of positronium.cpp and inside the production
// #ifndef, so everything it needs (State, the trajectory engine, simulate(),
// splitMix64) is already in scope.  It deliberately contains no ROOT: nothing
// in this file draws anything, and the panels that display these results live
// in positronium.cpp.

struct CremCollapseEstimate {
    double lifetimeSeconds=std::numeric_limits<double>::quiet_NaN();
    double calibrationSeconds=0.0;
    double meanRadiatedPowerWatts=std::numeric_limits<double>::quiet_NaN();
    SimulationOutcome calibrationOutcome=SimulationOutcome::NumericalFailure;
    // Distinguishes the two ways a trajectory can end as ObservationLimit.
    // Without it, "the reaction force is switched off, so this orbit will
    // never decay" is indistinguishable from "this orbit is decaying but ran
    // out of wall clock", and the report tells the user to raise a budget
    // that cannot possibly help.
    bool secularLossAbsent=false;
    // Kepler period of the quasi-closed orbit.  The inspiral is not a closed
    // orbit: T shrinks as a^(3/2) while the pair sinks, so a single number
    // cannot describe it.  Both ends of the run are reported instead, plus
    // the revolution count that the orbit-averaged integrator accumulates
    // between them (measured orbits and analytically skipped ones alike).
    double initialPeriodSeconds=std::numeric_limits<double>::quiet_NaN();
    double finalPeriodSeconds=std::numeric_limits<double>::quiet_NaN();
    double revolutions=std::numeric_limits<double>::quiet_NaN();
    // Osculating elements of the prepared orbit, and the closed-form classical
    // prediction they imply.  These make the run testable against textbook
    // electrodynamics rather than only against itself.
    double initialSemiMajorAxis=std::numeric_limits<double>::quiet_NaN();
    double initialEccentricity=std::numeric_limits<double>::quiet_NaN();
    double analyticCollapseSeconds=std::numeric_limits<double>::quiet_NaN();
    // Mean over checkpoints of (measured orbital energy loss rate) / (Larmor
    // rate for the same osculating orbit).  1 means the engine reproduces
    // coherent electric-dipole radiation exactly.
    double larmorPowerRatio=std::numeric_limits<double>::quiet_NaN();
    // Classical dipole-dipole interaction energy of the prepared pair,
    // expressed as a frequency so it can sit beside the measured o-Ps/p-Ps
    // hyperfine splitting.
    double dipoleCouplingHz=std::numeric_limits<double>::quiet_NaN();
};

struct OsculatingElements { double specificEnergy=0.0; double specificAngularMomentum=0.0; };

// --- Closed-form classical inspiral, used only as an external reference ---
//
// The pair carries electric dipole moment d = e*r, so it radiates at the
// Larmor rate P = |d''|^2/(6 pi eps0 c^3) with |d''| = 2 k e^3/(m r^2).
// Feeding that into dE/dt with E = -k e^2/(2a) gives da/dt = -C/a^2 with
//
//     C = 8 k e^4 / (6 pi eps0 c^3 m^2),
//
// so the time from a_i to a_f is (a_i^3 - a_f^3)/(3 C).  Nothing in the CREM
// engine is used to obtain this; it is the textbook answer the engine is
// being compared against.
// C = 2 q_eff^2 k|q1 q2| / (6 pi eps0 c^3 mu^2), with q_eff the pair's
// effective dipole charge and mu the reduced mass.  For e+e- (q_eff = e,
// mu = m/2) this collapses to the familiar 8 k e^4/(6 pi eps0 c^3 m^2), which
// is what the previous form hard-coded along with the equal-mass assumption.
// The active pair is selected by --pair, so these inputs are runtime globals.
// Keeping this function constexpr is ill-formed even though GCC historically
// accepted it as an extension with -Winvalid-constexpr.
inline double classicalInspiralCoefficient() {
    return 2.0*pairDipoleCharge*pairDipoleCharge*pairCoulombStrength
        /(6.0*pi*epsilon0*c*c*c*pairReducedMass*pairReducedMass);
}

// Orbit-averaged enhancement of DIPOLE radiation on a Kepler ellipse.
// P ~ r^-4, and <r^-4> = (1 + e^2/2)/(a^4 (1-e^2)^(5/2)) by direct
// integration with dt = (r^2/h) dtheta.  This is NOT the Peters factor
// (1 + 73e^2/24 + 37e^4/96)/(1-e^2)^(7/2), which belongs to QUADRUPOLE
// (gravitational) radiation with P ~ r^-6 and is far larger: at e = 0.5 the
// two differ by more than a factor of two.
double dipoleEccentricityFactor(double eccentricity) {
    const double e2=eccentricity*eccentricity;
    if(!(e2<1.0)) return std::numeric_limits<double>::quiet_NaN();
    return (1.0+0.5*e2)/std::pow(1.0-e2,2.5);
}

// Larmor power of the coherent electric dipole, orbit-averaged over the
// Kepler ellipse with the given elements.
double larmorOrbitAveragedPower(double semiMajorAxis,double eccentricity) {
    if(!(semiMajorAxis>0.0)) return std::numeric_limits<double>::quiet_NaN();
    // |d''| = |q_eff| k|q1 q2| / (mu a^2), so P = |d''|^2/(6 pi eps0 c^3).
    const double secondDerivative=magnitude(pairDipoleCharge)
        *pairCoulombStrength/(pairReducedMass*semiMajorAxis*semiMajorAxis);
    const double circular=secondDerivative*secondDerivative
        /(6.0*pi*epsilon0*c*c*c);
    return circular*dipoleEccentricityFactor(eccentricity);
}

// Time for the closed-form inspiral to carry the orbit from a_i to a_f at
// fixed eccentricity.  Radiation actually circularizes the orbit, so holding
// e fixed is an approximation; it is stated on the panel that uses this.
double classicalInspiralSeconds(double initialSemiMajorAxis,
                                double finalSemiMajorAxis,
                                double eccentricity) {
    const double factor=dipoleEccentricityFactor(eccentricity);
    if(!(initialSemiMajorAxis>finalSemiMajorAxis)||!std::isfinite(factor))
        return std::numeric_limits<double>::quiet_NaN();
    return (std::pow(initialSemiMajorAxis,3.0)-std::pow(finalSemiMajorAxis,3.0))
        /(3.0*classicalInspiralCoefficient()*factor);
}

// Periapsis distance of the two-body Kepler ellipse with the given specific
// (per unit reduced mass) energy and angular momentum -- the same formula
// simulate() already uses for its InitialConditions.predictedClosestApproach.
double osculatingPeriapsis(const OsculatingElements& elements,
                           double attractionParameter) {
    if(!(elements.specificAngularMomentum!=0.0)) return 0.0;
    const double eccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*elements.specificAngularMomentum
            *elements.specificAngularMomentum
            /(attractionParameter*attractionParameter)));
    return elements.specificAngularMomentum*elements.specificAngularMomentum
        /(attractionParameter*(1.0+eccentricity));
}

// Apoapsis of the two-body Kepler ellipse -- osculatingPeriapsis's mirror
// root (a(1+e) instead of a(1-e)).  Infinity for an unbound (e>=1) orbit.
double osculatingApoapsis(const OsculatingElements& elements,
                          double attractionParameter) {
    const double L=elements.specificAngularMomentum;
    if(!(L!=0.0)) return 0.0;
    const double eccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*L*L
            /(attractionParameter*attractionParameter)));
    if(!(eccentricity<1.0)) return std::numeric_limits<double>::infinity();
    return L*L/(attractionParameter*(1.0-eccentricity));
}

// Unperturbed Kepler period at the given specific energy; used only to size
// the one-orbit measurement window, not to compute the reported lifetime.
double osculatingPeriod(double specificEnergy,double attractionParameter) {
    const double semiMajorAxis=-attractionParameter/(2.0*specificEnergy);
    return 2.0*pi*std::sqrt(semiMajorAxis*semiMajorAxis*semiMajorAxis
        /attractionParameter);
}

// Radial potential energy (per unit reduced mass) for the Plummer-softened
// force law: the specific radial force is
// F_r(r) = -attractionParameter/(r^2+floor^2), matching the true Coulomb
// force -attractionParameter/r^2 for r >> floor and saturating smoothly
// (rather than diverging) as r -> 0, with clampedSeparationVector in
// positronium.cpp (r_eff = sqrt(r^2+floor^2)) as the force-law-level source
// of this.  Integrating -dU/dr = F_r gives
//
//     U(r) = -(attractionParameter/floor) * atan(floor/r),
//
// with the integration constant fixed so U -> -attractionParameter/r as
// r -> infinity (atan(floor/r) ~ floor/r there), matching the unclamped
// Coulomb potential exactly in that limit.  Valid and smooth for every
// r > 0: unlike the earlier hard-clamp regularization, there is no separate
// branch for r above/below floor here at all -- one formula, everywhere,
// which is what removes the kink in the force that stalled the adaptive
// integrator at the old floor crossing.
double regularizedPotentialEnergy(double r,double attractionParameter,
                                  double floor) {
    return -(attractionParameter/floor)*std::atan(floor/r);
}

// Radius where the softened force above exactly balances the centrifugal
// term for angular momentum L -- the circular-orbit radius under this force
// law.  Setting -dU/dr = L^2/r^3 gives
// attractionParameter/(r^2+floor^2) = L^2/r^3, i.e.
//
//     attractionParameter*r^3 - L^2*r^2 - L^2*floor^2 = 0,
//
// a cubic with exactly one positive real root (Descartes: coefficient signs
// +,-,0,- change sign once), interpolating smoothly between the deep limit
// (floor dominates, r ~ floor) and the standard Kepler circular-orbit radius
// L^2/attractionParameter (the floor -> 0 limit -- the identity this reduces
// to away from the barrier).  No simpler closed form survives adding the
// floor^2 term, so this is found by bisection.  The bracket is a physically
// motivated guess (floor plus the naive circular radius), doubled until it
// brackets a root as cheap insurance against a bad guess rather than a claim
// of precision: bisection's absolute convergence over even a wildly
// oversized bracket reaches machine precision at the relevant scale well
// within the iteration budget below regardless of how tight the guess was.
double criticalRadius(double L,double attractionParameter,double floor) {
    const auto cubic=[&](double r) {
        return attractionParameter*r*r*r-L*L*(r*r+floor*floor);
    };
    double lo=std::numeric_limits<double>::min();
    double hi=2.0*(floor+L*L/attractionParameter);
    for(int i=0;i<200&&!(cubic(hi)>0.0);++i) hi*=2.0;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(lo+hi);
        if(cubic(mid)>0.0) hi=mid; else lo=mid;
    }
    return 0.5*(lo+hi);
}

struct RegularizedTurningPoints {
    double periapsis,apoapsis;
    // False only in the naive-fallback branch below: there periapsis/
    // apoapsis are NOT roots of h(r) (h(r_min) missed even the marginal
    // circular-orbit tolerance, so no such roots exist), they are the plain
    // osculatingPeriapsis/osculatingApoapsis values returned as the
    // least-bad guess for STATE CONSTRUCTION (regularizedPeriapsis, via
    // osculatingPeriapsisState, where any single r is better than none).
    // regularizedPeriod's quadrature integrates h(r) between periapsis and
    // apoapsis, which is only valid when they truly are h's roots -- with
    // the naive fallback pair instead, h is negative at every quadrature
    // node (that is exactly what "not even marginal" means), so every node
    // gets skipped and the sum silently comes out zero.  Measured directly
    // (seed 5, the checkpoint where h(r_min)/localScale missed the 5%
    // margin by roughly a percentage point): regularizedPeriod returned
    // 0 ps before this flag existed to stop it from even trying.
    bool exact;
};

// Periapsis and apoapsis under the Plummer-softened force law
// (regularizedPotentialEnergy), replacing the naive (pure 1/r) formulas
// osculatingPeriapsis/osculatingApoapsis wherever the orbit's own scale
// approaches separationFloor().  h(r) = specificEnergy -
// regularizedPotentialEnergy(r) - L^2/(2r^2) is smooth for every r > 0 --
// no piecewise construction, unlike the hard-clamp regularization this
// replaces, since the force itself is now one formula everywhere -- concave
// with a single interior maximum at r = criticalRadius(...), and
// h -> -infinity as r -> 0 (centrifugal term dominates) or r -> infinity
// (bound orbit, specificEnergy < 0), so it has at most two roots straddling
// that maximum, found by bisecting on each side of it.  Cost is negligible
// (a few hundred evaluations of an elementary function) against the
// mechanical integration this sizes a measurement for, so -- unlike the
// hard-clamp predecessor -- there is no need to shortcut to the naive
// formula when the orbit is far from the barrier: the exact formula already
// reduces to it there (the correction is a part in 10^7 at Bohr-radius
// scale), so always using it is both simpler and correct everywhere.
//
// If h(r_min) itself is not positive, the (E,L) pair has no strictly
// accessible orbit at all under this force law -- but h is smooth in (E,L),
// so a small negative h(r_min) means the pair is almost exactly circular AT
// r_min (h(r_min)=0 is precisely that circular orbit), not that no orbit
// exists.  Comparing |h(r_min)| against the local force scale there
// (attractionParameter*r_min/(r_min^2+floor^2), i.e. force(r_min)*r_min, the
// specific potential step across a distance ~r_min) keeps the margin
// scale-free.  Within that margin both turning points collapse to r_min (a
// genuinely circular orbit); beyond it the naive values are returned as the
// least-bad fallback, since fabricating a root that is not there would be
// worse than saying so.
RegularizedTurningPoints regularizedTurningPoints(
    const OsculatingElements& elements,double attractionParameter,
    double floor) {
    const double naivePeriapsis=
        osculatingPeriapsis(elements,attractionParameter);
    const double naiveApoapsis=
        osculatingApoapsis(elements,attractionParameter);
    const double L=elements.specificAngularMomentum;
    if(!(floor>0.0)||!(L!=0.0)) return {naivePeriapsis,naiveApoapsis,false};
    const double rMin=criticalRadius(L,attractionParameter,floor);
    if(!(rMin>0.0)||!std::isfinite(rMin))
        return {naivePeriapsis,naiveApoapsis,false};
    const auto h=[&](double r) {
        return elements.specificEnergy
            -regularizedPotentialEnergy(r,attractionParameter,floor)
            -L*L/(2.0*r*r);
    };
    const double hAtMin=h(rMin);
    if(!(hAtMin>0.0)) {
        const double localScale=
            attractionParameter*rMin/(rMin*rMin+floor*floor);
        constexpr double marginalRelativeTolerance=0.05;
        if(localScale>0.0&&-hAtMin<=marginalRelativeTolerance*localScale)
            return {rMin,rMin,true};
        return {naivePeriapsis,naiveApoapsis,false};
    }
    double periapsisLo=std::numeric_limits<double>::min();
    double periapsisHi=rMin;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(periapsisLo+periapsisHi);
        if(h(mid)>0.0) periapsisHi=mid; else periapsisLo=mid;
    }
    double apoapsisLo=rMin;
    double apoapsisHi=(std::isfinite(naiveApoapsis)&&naiveApoapsis>rMin)
        ?naiveApoapsis:2.0*rMin;
    for(int i=0;i<200&&!(h(apoapsisHi)<0.0);++i) apoapsisHi*=2.0;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(apoapsisLo+apoapsisHi);
        if(h(mid)>0.0) apoapsisLo=mid; else apoapsisHi=mid;
    }
    return {0.5*(periapsisLo+periapsisHi),0.5*(apoapsisLo+apoapsisHi),true};
}

double regularizedPeriapsis(const OsculatingElements& elements,
                            double attractionParameter,double floor) {
    return regularizedTurningPoints(elements,attractionParameter,floor)
        .periapsis;
}

// Regularized orbital period, replacing osculatingPeriod() once the orbit's
// own scale approaches separationFloor().  osculatingPeriod() sizes the
// one-orbit measurement window in estimateCremCollapse(); once
// regularizedTurningPoints() starts differing materially from the naive
// values, that window is no longer one true period of the actual (softened)
// orbit, and the mismatch compounds into the jump extrapolation the same
// way periapsis's own mismatch did before it was fixed (c909ee8/2bb2751).
//
// No closed form exists here in general -- the softened force does not
// admit the standard eccentric-anomaly substitution -- so the period is
// found by direct quadrature of T = 2 * integral[r_p,r_a] dr / sqrt(2 h(r)),
// using the substitution r(theta) = (r_p+r_a)/2 + (r_a-r_p)/2 * cos(theta),
// theta in (0,pi): this removes BOTH endpoint (turning-point) singularities
// for ANY smooth h with simple roots there, because near either root
// h(r) ~ h'(r_root)*(r-r_root) while r(theta)-r_root ~ -/+(r_a-r_p)/4 *
// theta^2 near theta=0/pi, so h(r(theta)) ~ theta^2 while dr/dtheta ~ theta,
// and the ratio in the integrand stays finite.  h is now the single smooth
// formula everywhere (no piecewise potential to match across a boundary the
// way the hard-clamp version needed), so this integral is, if anything,
// better conditioned than its predecessor: there is no interior kink in
// curvature left for the quadrature to sit near.
//
// Fixed-node midpoint quadrature is used rather than Gauss-Legendre so no
// tabulated node/weight constants need transcribing; every node falls
// strictly inside (0,pi), so the removable endpoint singularities are never
// evaluated at all.  N=200 carries over unchanged from the hard-clamp
// version (checked there to machine precision against the closed form in
// the far-field limit); cost is negligible either way (a few hundred h(r)
// evaluations against the full mechanical integration this sizes the window
// for), and always running the same quadrature -- rather than shortcutting
// to osculatingPeriod() far from the barrier -- keeps this function as
// simple as regularizedTurningPoints() above for the same reason.
double regularizedPeriod(const OsculatingElements& elements,
                         double attractionParameter,double floor) {
    const double naive=osculatingPeriod(
        elements.specificEnergy,attractionParameter);
    if(!(floor>0.0)) return naive;
    const RegularizedTurningPoints turningPoints=
        regularizedTurningPoints(elements,attractionParameter,floor);
    // Naive fallback (turningPoints.exact == false): periapsis/apoapsis are
    // NOT roots of h(r) below, so integrating between them is meaningless --
    // h is negative at every quadrature node (that is exactly what "missed
    // even the marginal tolerance" means), every node gets skipped, and the
    // sum silently comes out zero.  Measured directly: this is precisely
    // what happened before this check existed (seed 5, one checkpoint whose
    // h(r_min)/localScale missed the 5% margin by about a point -- naive is
    // still the least-bad state-construction guess regularizedPeriapsis
    // hands to osculatingPeriapsisState in that case, but it is not a valid
    // integration bound here, so this falls back to osculatingPeriod()
    // instead of the quadrature).
    if(!turningPoints.exact) return naive;
    const double periapsis=turningPoints.periapsis;
    const double apoapsis=turningPoints.apoapsis;
    const double L=elements.specificAngularMomentum;
    if(!(apoapsis>periapsis)) {
        // Degenerate: both turning points collapsed onto (about) the same
        // radius -- regularizedTurningPoints' own marginal branch, an orbit
        // sitting almost exactly at the bottom of the effective-potential
        // well.  The RADIAL period is not just small there, it is the wrong
        // question: for a force law that does not satisfy Bertrand's
        // theorem, radial and angular period differ, and "one orbit" for a
        // near-circular orbit means one REVOLUTION, not one (vanishing)
        // radial oscillation.  phi-dot = L/r^2 by definition of specific
        // angular momentum, so the angular period is 2*pi*r^2/L at the
        // shared radius.
        const double r=0.5*(periapsis+apoapsis);
        return (L!=0.0)?2.0*pi*r*r/std::abs(L):naive;
    }
    const auto h=[&](double r) {
        return elements.specificEnergy
            -regularizedPotentialEnergy(r,attractionParameter,floor)
            -L*L/(2.0*r*r);
    };
    constexpr int quadratureNodes=200;
    double halfPeriod=0.0;
    for(int i=0;i<quadratureNodes;++i) {
        const double theta=(static_cast<double>(i)+0.5)*pi
            /static_cast<double>(quadratureNodes);
        const double r=0.5*(periapsis+apoapsis)
            +0.5*(apoapsis-periapsis)*std::cos(theta);
        const double radialKineticEnergyTimesTwo=2.0*h(r);
        if(!(radialKineticEnergyTimesTwo>0.0)) continue; // guard only; every
            // node sits strictly between the turning points, so this should
            // not trigger, but a stray non-finite h must not poison the sum.
        halfPeriod+=std::sin(theta)/std::sqrt(radialKineticEnergyTimesTwo);
    }
    halfPeriod*=0.5*(apoapsis-periapsis)*pi/static_cast<double>(quadratureNodes);
    return 2.0*halfPeriod;
}

// Fresh State at periapsis for the given osculating elements, carrying the
// supplied dipole vectors over unchanged.  The orbital plane is reset to the
// canonical x-y plane every time: only the (E,L) magnitudes are propagated
// secularly, and the plane orientation is irrelevant to a collapse-time
// estimate, so nothing is lost by not tracking it.
State osculatingPeriapsisState(const OsculatingElements& elements,
                               double attractionParameter,
                               const Vec3& firstDipole,
                               const Vec3& secondDipole) {
    // regularizedPeriapsis(), not the naive formula: below separationFloor()
    // the naive value systematically overshoots how close the pair can
    // actually get (it assumes the force keeps growing as 1/r^2 rather than
    // being pinned), so teleporting a measurement orbit there starts it from
    // a position/speed pair that is not a real turning point of the clamped
    // force law -- exactly the mismatch that made the one-orbit measurement
    // fail outright once periapsis approached the barrier.  Tangential speed
    // is still L/r at any turning point regardless of the force law (that is
    // the definition of specific angular momentum at zero radial velocity),
    // so only the r fed into it needs correcting.
    const double periapsis=regularizedPeriapsis(
        elements,attractionParameter,separationFloor());
    const double tangentialSpeed=periapsis>0.0
        ?elements.specificAngularMomentum/periapsis:0.0;
    const Vec3 relativePosition{periapsis,0.0,0.0};
    const Vec3 relativeVelocity{0.0,tangentialSpeed,0.0};
    State s;
    s.firstPosition=relativePosition*(secondMass/(firstMass+secondMass));
    s.secondPosition=relativePosition*(-firstMass/(firstMass+secondMass));
    s.firstVelocity=relativeVelocity*(secondMass/(firstMass+secondMass));
    s.secondVelocity=relativeVelocity*(-firstMass/(firstMass+secondMass));
    s.firstDipole=firstDipole;
    s.secondDipole=secondDipole;
    return s;
}

// Orbit-averaged secular integration.  Resolving each of the ~1e5-1e6
// individual orbits a classical inspiral needs (the direct approach this
// replaces) took upward of an hour per trajectory.  Instead, one
// representative orbit is fully resolved with the complete CREM engine --
// retardation, Darwin, dipole coupling, and radiation reaction under
// whichever --radiation-reaction model is active -- starting and ending near
// periapsis.  Its ACTUALLY MEASURED energy and angular-momentum loss (read
// directly from the engine's own exact radiatedEnergy / radiatedAngularMomentum
// flux bookkeeping, not a formula) gives the current dissipation rate.  A
// bounded number of further orbits, capped so the implied energy change per
// jump stays small, are then skipped analytically at that rate before the
// next orbit is resolved and the rate re-measured.  This is the standard
// osculating-element / orbit-averaging technique for multi-timescale
// radiative inspirals (the same approximation EMRI gravitational-wave models
// use), and it stays exactly as sensitive to --radiation-reaction as brute
// force would be: with the reaction force disabled, the very first measured
// orbit shows zero energy loss and the estimate correctly reports no
// collapse instead of iterating uselessly.
CremCollapseEstimate estimateCremCollapse(std::uint64_t seed,
                                           int selectedPhenomenon,
                                           double wallClockBudgetSeconds) {
    CremCollapseEstimate result;
    const double reducedMass=firstMass*secondMass
        /(firstMass+secondMass);
    const double attractionParameter=pairCoulombStrength/reducedMass;

    // Sample the same random bound initial condition simulate() would use.
    // The near-zero observation window means nothing actually integrates;
    // frames.front() is the pristine sampled state (pushed before the loop).
    SimulationOptions seedOptions;
    seedOptions.frameCount=2;
    seedOptions.observationTime=1.0e-24;
    const SimulationResult seedRun=simulate(seed,selectedPhenomenon,seedOptions);
    if(seedRun.frames.empty()||!(seedRun.initial.relativeEnergy<0.0))
        return result;

    OsculatingElements elements{seedRun.initial.relativeEnergy/reducedMass,
        seedRun.initial.orbitalAngularMomentum/reducedMass};
    Vec3 firstDipole=seedRun.frames.front().firstDipole;
    Vec3 secondDipole=seedRun.frames.front().secondDipole;

    // Osculating elements of the prepared orbit and the closed-form classical
    // prediction they imply.  The run is stopped when the PERIAPSIS reaches
    // finalApproachMultiple*comptonBarrierRadius, so the reference must be
    // evaluated between the same two orbits, not down to a = 0.
    result.initialSemiMajorAxis=
        -attractionParameter/(2.0*elements.specificEnergy);
    result.initialEccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*elements.specificAngularMomentum
            *elements.specificAngularMomentum
            /(attractionParameter*attractionParameter)));
    // Classical dipole-dipole interaction energy of the prepared pair, as a
    // frequency, for comparison with the measured hyperfine splitting.
    {
        const Vec3 separation=seedRun.frames.front().first
                             -seedRun.frames.front().second;
        const double distance=separation.norm();
        if(distance>0.0) {
            const Vec3 direction=separation/distance;
            const double coupling=(mu0/(4.0*pi))
                *(dot(firstDipole,secondDipole)
                  -3.0*dot(firstDipole,direction)
                      *dot(secondDipole,direction))
                /(distance*distance*distance);
            result.dipoleCouplingHz=coupling/(2.0*pi*hbar);
        }
    }
    double larmorRatioSum=0.0;
    int larmorRatioCount=0;
    // Cached restart-artefact ratio and the osculating energy it was taken at.
    // Measured trade-off on one full trajectory (seed 42), against the
    // un-thinned collapse time of 131.685 ps:
    //
    //     threshold   wall clock   collapse time   error
    //       0.02        16.27 s      131.685 ps    0.000%   (refreshes every
    //       0.05        13.32 s      130.602 ps    0.822%    checkpoint, so it
    //       0.10        11.26 s      127.694 ps    3.031%    reproduces the
    //                                                        un-thinned answer
    //                                                        exactly)
    //
    // 0.05 is the working point: the error stays well inside the 3% per-jump
    // tolerance the frozen rate already accepts, so the background stays a
    // subordinate approximation.  0.10 was tried first on the assumption that
    // it would too; it does not -- there the background error EQUALS the
    // frozen-rate tolerance and would become co-dominant.  Raising it buys
    // speed at a bias that grows faster than linearly.
    constexpr double backgroundRefreshFraction=0.05;
    bool haveBackgroundRatio=false;
    double backgroundEnergyRatio=0.0;
    double backgroundAngularRatio=0.0;
    double energyAtLastBackground=0.0;
    // Read once, thread-safely: runCremCollapseExperiment runs many of these
    // concurrently, so the global must not be re-read (or, worse, mutated)
    // from here on.
    const ChargeRadiationReactionModel activeReactionModel=gRadiationReactionModel;

    const auto wallClockStart=std::chrono::steady_clock::now();
    const auto wallClockSpent=[&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now()-wallClockStart).count();
    };
    // Size of one analytic jump, as the dimensionless s = (3/2) n du0/u0 of
    // the resummed solution below.  s must stay strictly under 1, where the
    // closed form reaches complete collapse.
    //
    // Convergence on one full trajectory (seed 42), collapse time in ps:
    //
    //     s_max    wall clock   collapse time
    //     0.012      31.0 s       124.169
    //     0.045      12.6 s       124.957     <- same step as the old 3% cap
    //     0.10        8.4 s       125.585
    //     0.20        5.7 s       124.986
    //     0.30        4.6 s       124.328
    //
    // Flat to +/-0.6% over a 25x range of step size, so 0.30 is taken.  The
    // frozen-rate hold this replaces returned 130.6-131.7 ps, i.e. 5-6% HIGH,
    // and it was never convergence-tested: at the SAME step size (s = 0.045)
    // the two schemes differ by 4.5%, which is interpolation shape, not step
    // size.  The sign is the expected one -- a zeroth-order hold freezes a
    // loss rate that actually grows as the orbit tightens, so it under-counts
    // the loss and over-states the collapse time.  This change therefore buys
    // accuracy as well as speed.
    constexpr double maximumJumpParameter=0.30;
    constexpr int maxOrbitsSkippedAtOnce=200000;
    constexpr int maxCheckpoints=4000; // ~150 are needed at 3%/jump to cover
        // a factor 100 in |E|; this is a generous safety margin, not a
        // normal stopping point.
    double simulatedTimeTotal=0.0;
    double radiatedEnergyTotal=0.0;
    // Revolutions completed so far.  Every orbit the loop accounts for is
    // counted here exactly once: the resolved measurement orbit is the first
    // of the orbitsToSkip that each checkpoint advances over.
    double revolutionsTotal=0.0;
    // Period of the orbit the pair actually starts on, recorded at the first
    // checkpoint before any secular decay has been applied.
    result.initialPeriodSeconds=osculatingPeriod(
        elements.specificEnergy,attractionParameter);

    for(int checkpoint=0;checkpoint<maxCheckpoints;++checkpoint) {
        // Refreshed every checkpoint so that whichever branch returns below
        // reports the period of the orbit reached at that point.
        // regularizedPeriod(), not the naive formula: below separationFloor()
        // the naive value systematically under- or overestimates the true
        // period (see regularizedPeriod's own comment), and this is a
        // user-facing report, not just internal bookkeeping.
        result.finalPeriodSeconds=regularizedPeriod(
            elements,attractionParameter,separationFloor());
        result.revolutions=revolutionsTotal;
        if(larmorRatioCount>0) {
            result.larmorPowerRatio=larmorRatioSum
                /static_cast<double>(larmorRatioCount);
        }
        // Closed-form prediction for the stretch actually covered so far, so
        // the reference is always evaluated between the same two orbits as
        // the measurement it will be compared with.
        result.analyticCollapseSeconds=classicalInspiralSeconds(
            result.initialSemiMajorAxis,
            -attractionParameter/(2.0*elements.specificEnergy),
            result.initialEccentricity);
        if(const char* debug=std::getenv("CREM_DEBUG");debug) {
            std::cerr<<"checkpoint "<<checkpoint<<" t="<<simulatedTimeTotal*1e12
                     <<"ps E="<<elements.specificEnergy<<" L="
                     <<elements.specificAngularMomentum<<" wall="
                     <<wallClockSpent()<<"s"<<std::endl;
        }
        if(wallClockSpent()>wallClockBudgetSeconds) {
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal;
            return result;
        }
        const double periapsis=osculatingPeriapsis(elements,attractionParameter);
        if(!(periapsis>0.0)||!std::isfinite(periapsis)) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            return result;
        }
        // Stop the orbit-averaged phase exactly AT comptonBarrierRadius --
        // the same physically meaningful scale ("where classical point-
        // particle electrodynamics stops applying") already used as the
        // Collision boundary in experiment 5 (makeInteractionConfiguration),
        // so "collapsed" now means the same separation there and here rather
        // than two different numbers that happened to both be called a
        // stopping point.  This used to be 10x further out (finalApproachMultiple
        // was 10.0): the one-period measurement window this loop relies on
        // stopped being a good approximation approaching the barrier under
        // the OLD hard-clamp force law (a real kink in the force at the
        // clamp boundary, see f477866), so the loop backed off well before
        // it rather than chasing an unreliable measurement for a truncation
        // that (t~a^3 near collapse) only cost 0.1% of the total collapse
        // time anyway.
        //
        // With that force-law kink removed (Plummer softening, f477866) and
        // the energy signal itself no longer riding on the Darwin
        // approximation (orbitalRadiatedEnergy, 0ebc7d2), the 10x margin is
        // no longer needed: measured directly on the same 6 seeds used
        // throughout this investigation, finalApproachMultiple=1.0 (i.e.
        // stopping exactly at the barrier) completes cleanly on all six,
        // zero numerical failures, zero censored.  0.5 (past the barrier)
        // already is not reliably safe -- 2 of 6 seeds hit NumericalFailure
        // or exhausted their wall-clock budget there, the boundary this
        // whole investigation has independently kept finding from the force-
        // law and Darwin-approximation sides -- so 1.0 is not an arbitrary
        // round number, it is measured to be the deepest value that is
        // consistently reachable, and it happens to coincide exactly with
        // where the model already declares Collision elsewhere.
        constexpr double finalApproachMultiple=1.0;
        if(periapsis<=finalApproachMultiple*comptonBarrierRadius) {
            result.lifetimeSeconds=simulatedTimeTotal;
            result.meanRadiatedPowerWatts=simulatedTimeTotal>0.0
                ?radiatedEnergyTotal/simulatedTimeTotal
                :std::numeric_limits<double>::quiet_NaN();
            result.calibrationOutcome=SimulationOutcome::ReachedCutoff;
            result.calibrationSeconds=simulatedTimeTotal;
            return result;
        }

        const State measurementState=osculatingPeriapsisState(
            elements,attractionParameter,firstDipole,secondDipole);
        // regularizedPeriod(), not the naive formula: this sizes the ONE
        // orbit the measurement below actually integrates, and the naive
        // period stops matching a true orbit of the clamped force law
        // exactly where osculatingPeriapsisState's own teleport target does
        // (see regularizedPeriod's comment for the measured symptom this
        // fixes).
        const double period=regularizedPeriod(
            elements,attractionParameter,separationFloor());
        SimulationOptions measureOptions;
        measureOptions.collectFrames=false;
        measureOptions.observationTime=period;
        measureOptions.terminalSeparation=comptonBarrierRadius;
        // The ENERGY side of the secular update reads
        // finalState.orbitalRadiatedEnergy now, not a
        // conservativeParticleEnergy() difference (see below), so the flux
        // quadrature this enables is no longer pure overhead for run -- it
        // is the measurement.  background still only needs the mechanical
        // (position, velocity) endpoint for its angular-momentum
        // subtraction, so its own copy below turns this back off.
        measureOptions.radiatedEnergyBookkeeping=true;
        // The budget is shared with the whole estimate: a single
        // measurement orbit can itself run into the stiff region near the
        // boundary, so it must be interruptible too.
        measureOptions.stopRequested=[&]() {
            return wallClockSpent()>wallClockBudgetSeconds;
        };
        // nuclearCutoff, not comptonBarrierRadius: measurementState's own
        // position is regularizedPeriapsis() (via osculatingPeriapsisState),
        // which is allowed to sit AT or below comptonBarrierRadius once the
        // orbit has circularized deep enough (see the r_min branch there) --
        // passing comptonBarrierRadius here as well would then make this
        // measurement start already past its own trajectoryCutoff, a
        // zero-step degenerate case (measured directly: startSep=189.5fm <
        // cutoff=193.3fm, elapsed=0ps, NumericalFailure).  nuclearCutoff is
        // the universal, pair-independent "genuine point-particle floor"
        // already used everywhere else in the codebase for exactly this
        // role (e.g. lienardWiechertField's "whichever floor is bigger
        // wins"), and it is always strictly below comptonBarrierRadius, so
        // this measurement now always starts with room to run regardless of
        // how deep regularizedPeriapsis has to go.  In the far-from-barrier
        // regime this changes nothing observable: the measured orbit's
        // deepest point stays far above either cutoff there, so which one
        // is configured is moot.
        const MechanicalTrajectoryResult run=runMechanicalTrajectory(
            measurementState,period,nuclearCutoff,measureOptions,
            activeReactionModel);

        // Angular momentum and energy are both read from the measured
        // orbit's actual start/end state rather than from a formula: L is
        // exactly conserved by a central force regardless of its radial
        // profile, and conservativeParticleEnergy() is the same quantity
        // positronium_validation checks to 1e-12, so together they isolate
        // the genuine secular drift with no reaction-model-dependent
        // bookkeeping in the way.  (radiatedEnergy / radiatedAngularMomentum
        // were tried and rejected: they are system-total EM flux bookkeeping
        // that bundles in the M1 magnetic-dipole/spin-precession channel,
        // which does not recoil the orbit at all, and dominated the total by
        // orders of magnitude here.)
        //
        // Even so, restarting each measurement from a fresh, from-rest
        // "periapsis" State leaves a small settling artefact in the
        // reconstructed retarded-field history (causalInitialHistory has no
        // real orbit behind it to extrapolate from) that reads as a spurious
        // secular loss growing linearly with elapsed time -- present even
        // with --radiation-reaction disabled, where the mechanics guarantee
        // exact conservation.  A second, otherwise identical measurement run
        // with the reaction force explicitly forced off calibrates exactly
        // that artefact away by subtraction, leaving only the effect the
        // active --radiation-reaction model actually adds.  It does not
        // touch gRadiationReactionModel, so it stays safe across the worker
        // threads in runCremCollapseExperiment.
        //
        // The background is the second of the two full-orbit integrations
        // every checkpoint used to pay for, i.e. half the cost of the whole
        // estimate.  Measured across a complete trajectory it is a large but
        // SMOOTH fraction of the raw signal: background/raw runs 0.42..0.74
        // and drifts by a median 0.52% per checkpoint.  What is cached here is
        // therefore that RATIO rather than the absolute artefact, which grows
        // with the orbit and would go stale far faster.  The ratio does jump
        // by up to 12% where the orbit changes character, so the refresh is
        // driven by how far the osculating energy has moved since the cached
        // ratio was taken, not by a fixed checkpoint count: near the boundary,
        // where the jumps live, it re-measures more often on its own.
        const bool refreshBackground = !haveBackgroundRatio
            || !(std::abs(elements.specificEnergy-energyAtLastBackground)
                 <= backgroundRefreshFraction*std::abs(energyAtLastBackground));
        MechanicalTrajectoryResult background;
        if(refreshBackground) {
            // Same cutoff as run above, and for the same reason: it must
            // match run's so the pair of integrations stays comparable (the
            // background subtraction assumes both runs cover the same
            // stretch), and it must clear measurementState's own position.
            // Its own copy of measureOptions with radiatedEnergyBookkeeping
            // turned back off: background is only ever read for its
            // mechanical endpoint (angular momentum), never its
            // orbitalRadiatedEnergy, so the flux quadrature run buys the
            // measurement would be pure overhead here.
            SimulationOptions backgroundOptions=measureOptions;
            backgroundOptions.radiatedEnergyBookkeeping=false;
            background=runMechanicalTrajectory(
                measurementState,period,nuclearCutoff,backgroundOptions,
                ChargeRadiationReactionModel::disabled);
        }
        const auto measuredDelta=[&](const State& end) {
            const Vec3 relPos=end.firstPosition-end.secondPosition;
            const Vec3 relVel=end.firstVelocity-end.secondVelocity;
            return OsculatingElements{
                (conservativeParticleEnergy(end)
                    -conservativeParticleEnergy(measurementState))/reducedMass,
                cross(relPos,relVel).z-elements.specificAngularMomentum};
        };

        // Either the freshly measured artefact, or the cached ratio applied to
        // this checkpoint's own raw measurement.
        const auto backgroundFor=[&](const OsculatingElements& realDelta) {
            if(refreshBackground&&isFinite(background.finalState))
                return measuredDelta(background.finalState);
            if(haveBackgroundRatio)
                return OsculatingElements{
                    backgroundEnergyRatio*realDelta.specificEnergy,
                    backgroundAngularRatio*realDelta.specificAngularMomentum};
            return OsculatingElements{};
        };

        if(run.outcome==SimulationOutcome::ReachedCutoff) {
            result.lifetimeSeconds=simulatedTimeTotal+run.elapsedTime;
            // Already in real (not specific) Joules, and already exactly the
            // E1 energy this partial orbit radiated -- no background
            // subtraction needed (see the full explanation where this same
            // substitution is made below, for the ordinary ObservationLimit
            // path this ReachedCutoff path is the truncated-orbit sibling of).
            radiatedEnergyTotal+=run.finalState.orbitalRadiatedEnergy;
            result.meanRadiatedPowerWatts=result.lifetimeSeconds>0.0
                ?radiatedEnergyTotal/result.lifetimeSeconds
                :std::numeric_limits<double>::quiet_NaN();
            result.calibrationOutcome=SimulationOutcome::ReachedCutoff;
            result.calibrationSeconds=result.lifetimeSeconds;
            // The measurement orbit ran into the boundary partway round, so
            // only the fraction of a revolution it actually covered counts.
            result.revolutions=revolutionsTotal
                +(period>0.0?run.elapsedTime/period:0.0);
            return result;
        }
        if(run.outcome!=SimulationOutcome::ObservationLimit
           ||!isFinite(run.finalState)
           ||(refreshBackground&&!isFinite(background.finalState))) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal+run.elapsedTime;
            return result;
        }
        if(wallClockSpent()>wallClockBudgetSeconds) {
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal+run.elapsedTime;
            return result;
        }

        const double measuredElapsed=std::max(run.elapsedTime,1.0e-30);
        // ANGULAR MOMENTUM only, from here on: realDelta/backgroundDelta's
        // specificEnergy field is computed (conservativeParticleEnergy(), a
        // Darwin/near-field approximation) but deliberately not read for the
        // loss below any more -- see deltaEnergyPerOrbit's own comment.
        // backgroundEnergyRatio is still computed in the cache refresh
        // alongside backgroundAngularRatio (the two are cached together),
        // but nothing downstream reads it either now; it is vestigial,
        // left rather than split the cache apart, and costs nothing to carry.
        const OsculatingElements realDelta=measuredDelta(run.finalState);
        const OsculatingElements backgroundDelta=backgroundFor(realDelta);
        // Refresh the cache only from an actual measurement, and only when the
        // raw values are large enough for the quotient to mean anything.
        if(refreshBackground&&isFinite(background.finalState)) {
            if(realDelta.specificEnergy!=0.0
               &&realDelta.specificAngularMomentum!=0.0) {
                backgroundEnergyRatio=backgroundDelta.specificEnergy
                    /realDelta.specificEnergy;
                backgroundAngularRatio=backgroundDelta.specificAngularMomentum
                    /realDelta.specificAngularMomentum;
                haveBackgroundRatio=true;
                energyAtLastBackground=elements.specificEnergy;
            }
        }
        // orbitalRadiatedEnergy, not a conservativeParticleEnergy()
        // difference: conservativeParticleEnergy() is the Darwin (static
        // near-field) approximation, valid only while the orbital period
        // stays many light-crossing times long, and measured directly to
        // fail exactly that way -- on seed 4, as period/(separation/c) fell
        // from ~950 to ~40, the conservativeParticleEnergy-based dE/E grew
        // from 3.6e-6 (matching the documented ~1e-5 integrator noise floor)
        // to OVER 100% per "orbit" (checkpoint 24: 1.73), while
        // orbitalRadiatedEnergy over the SAME stretch stayed smooth and
        // monotonic throughout (0.042 -> 0.272), only faltering mildly at
        // the very end.  It is computed from the exact retarded far-zone
        // Poynting flux (electromagneticFieldFluxRates via
        // particleMultipoleRadiation), with the M1 (magnetic-dipole/spin)
        // channel already excluded (see orbitalRadiatedEnergy's own comment
        // in state.hpp) -- and unlike conservativeParticleEnergy(measure-
        // mentState), which is a near-zero DIFFERENCE of two large numbers,
        // orbitalRadiatedEnergy starts at exactly 0 for a fresh
        // measurementState and accumulates the genuine physical loss
        // directly, so no background subtraction is needed for it: the
        // "settling artefact" this whole background machinery exists to
        // remove is specific to conservativeParticleEnergy's own near-field
        // bootstrap, not to a directly-measured flux.  Sign matches the
        // existing convention (negative = orbit binding tighter).
        const double deltaEnergyPerOrbit=
            -run.finalState.orbitalRadiatedEnergy/reducedMass;
        const double deltaAngularMomentumPerOrbit=
            realDelta.specificAngularMomentum-backgroundDelta.specificAngularMomentum;
        // A secular inspiral cannot shed a large fraction of its own binding
        // energy in ONE orbit; if it could, orbit averaging would not apply in
        // the first place.  Without this guard a single bad measurement is
        // folded straight into the osculating elements and destroys them.
        // Observed with --radiation-reaction coherent: the reaction force is
        // built from a third derivative of the dipole moment divided by
        // 2*h^3, and when the retarded-history spacing degenerates the stencil
        // returns garbage -- dE/E jumped 22 orders of magnitude in a single
        // checkpoint (2.7e-5 -> 1.4e17) while the reaction-free background
        // stayed at 1.9e-5.  The specific energy then sat at -1.2e30, every
        // later measurement read exactly zero, and the run was reported as
        // "not decaying", which named a downstream symptom rather than the
        // numerical failure that had actually occurred.
        constexpr double maxRelativeLossPerOrbit=0.5;
        if(!std::isfinite(deltaEnergyPerOrbit)
           ||std::abs(deltaEnergyPerOrbit)
               >maxRelativeLossPerOrbit*std::abs(elements.specificEnergy)) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal+measuredElapsed;
            return result;
        }
        if(const char* debug=std::getenv("CREM_DEBUG");debug) {
            std::cerr<<"  measured dE/E="
                     <<deltaEnergyPerOrbit/elements.specificEnergy
                     <<" (orbitalRadiatedEnergy="
                     <<run.finalState.orbitalRadiatedEnergy
                     <<"J) dL/L="<<deltaAngularMomentumPerOrbit
                        /elements.specificAngularMomentum<<std::endl;
        }
        // Measured orbital dissipation of this osculating orbit against the
        // Larmor rate for the same orbit.  Both sides are orbit averages over
        // the same period, so the ratio is a direct, parameter-free check of
        // the engine's radiation sector against textbook electrodynamics.
        {
            const double measuredPower=
                -deltaEnergyPerOrbit*reducedMass/measuredElapsed;
            const double semiMajorAxis=
                -attractionParameter/(2.0*elements.specificEnergy);
            const double eccentricity=std::sqrt(std::max(0.0,1.0
                +2.0*elements.specificEnergy*elements.specificAngularMomentum
                    *elements.specificAngularMomentum
                    /(attractionParameter*attractionParameter)));
            const double larmorPower=
                larmorOrbitAveragedPower(semiMajorAxis,eccentricity);
            if(std::isfinite(measuredPower)&&std::isfinite(larmorPower)
               &&larmorPower>0.0&&measuredPower>0.0) {
                larmorRatioSum+=measuredPower/larmorPower;
                ++larmorRatioCount;
            }
        }
        if(!(deltaEnergyPerOrbit<0.0)||!std::isfinite(deltaEnergyPerOrbit)) {
            // No measurable energy loss this orbit: with the reaction force
            // disabled (or below numerical noise), the orbit is not
            // secularly decaying, so there is nothing further to observe.
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal+measuredElapsed;
            result.secularLossAbsent=true;
            return result;
        }

        // --- Skipped orbits, integrated with the known functional form ---
        //
        // Freezing dE/dorbit across the jump is a zeroth-order hold, and it is
        // wrong in a known direction: the orbit tightens, so the true loss per
        // orbit GROWS and the hold always underestimates it.  That is why the
        // jump had to stay within 3% of |E|.
        //
        // With u = -E and dipole radiation, P ~ a^-4 and T ~ a^(3/2) give
        // du/dn = B u^(5/2), whose exact solution is
        //
        //     u(n) = u0 (1 - s)^(-2/3),      s = (3/2) n du0/u0,
        //
        // and expanding it to first order in s returns exactly the frozen-rate
        // answer, so this is a resummation of the same measurement rather than
        // a different model.  B is still taken from THIS checkpoint's measured
        // loss; only the shape of the interpolation between measurements
        // changes.  The period follows T(n) = T0 (1 - s), so the elapsed time
        // integrates to T0 n (1 - s/2) instead of T0 n.
        //
        // Checked against the engine on a full trajectory: fitting
        // log(dE/E) against log(u) gives slope 1.4617 (1.4465 over the first
        // three quarters) where this form predicts 3/2, so the exponent is
        // right to within 3%.  The angular momentum is carried by the measured
        // ratio k = (dL/L)/(du/u), which came out between -0.454 and -0.490
        // over that same trajectory, i.e. far steadier than dL/dorbit itself;
        // L then follows L = L0 (u/u0)^k.
        const double energyMagnitude=std::abs(elements.specificEnergy);
        const double lossPerOrbit=std::abs(deltaEnergyPerOrbit);
        int orbitsToSkip=1;
        if(lossPerOrbit>0.0&&energyMagnitude>0.0) {
            const double requestedOrbits=maximumJumpParameter*energyMagnitude
                /(1.5*lossPerOrbit);
            // Saturate while the value is still floating-point.  For a tiny
            // measured loss requestedOrbits can exceed INT_MAX or become +inf;
            // converting either value to int before clamp has undefined
            // behaviour.  The bounded value is always representable because
            // maxOrbitsSkippedAtOnce itself is an int.
            const double boundedOrbits=std::isfinite(requestedOrbits)
                ?std::clamp(requestedOrbits,1.0,
                    static_cast<double>(maxOrbitsSkippedAtOnce))
                :static_cast<double>(maxOrbitsSkippedAtOnce);
            orbitsToSkip=static_cast<int>(boundedOrbits);
        }
        const double jumpParameter=std::min(
            1.5*static_cast<double>(orbitsToSkip)*lossPerOrbit/energyMagnitude,
            maximumJumpParameter);
        const double energyGrowth=std::pow(1.0-jumpParameter,-2.0/3.0);
        const double updatedEnergyMagnitude=energyMagnitude*energyGrowth;
        // k, in CLOSED FORM rather than measured.  The previous approach
        // (k = (dL/L)/(dE/E) from this checkpoint's own measurement) tracked
        // smoothly while the orbit was still visibly eccentric but became
        // noise-dominated as e -> 0 (dL on one orbit shrinks toward the scale
        // of the background-subtraction artefact), and a bad k there gets
        // raised to a power and multiplied into L across tens of thousands of
        // skipped orbits, so the error compounds instead of averaging out:
        // on seed 42 the discriminant 1+2EL^2/mu^2 (=e^2), which can only
        // fall as the orbit radiates, instead climbed back from 0.0005 to
        // 0.15 over checkpoints 13-23 and went negative (unphysical) by
        // checkpoint 28.  EMA-smoothing the measurement was tried and made
        // it WORSE (negative by checkpoint 4), which is itself informative:
        // the drift is not noise around a fixed point, so damping it with a
        // stale average fights the real k, it does not filter a fake one.
        //
        // The real k is not empirical at all.  Orbit-averaging the dipole
        // reaction force -- a_react = gamma*r'''(t), Abraham-Lorentz on the
        // relative coordinate, gamma folding in q_eff^2/(6 pi eps0 mu c^3) --
        // over an unperturbed Kepler ellipse (same method as Peters 1964, but
        // for the dipole's P~r^-4 in place of the quadrupole's r^-6) gives,
        // with r'''=-(attractionParameter/r^3)[v-3(r-hat.v)r-hat]:
        //
        //   <dE/dt> = -gamma*attractionParameter^2*(1+e^2/2)
        //                 / (a^4 (1-e^2)^(5/2))
        //   <dL/dt> = -gamma*attractionParameter*L / (a^3 (1-e^2)^(3/2))
        //
        // using the standard time-averages <r^-n> = a^-n I_(n-1), I_n =
        // (1/2pi) integral (1-e cosE)^-n dE, generated by the recursion
        // I_n = I_(n-1) + (e/(n-1)) dI_(n-1)/de (from d/de[(1-e cosE)^-(n-1)]
        // integrated over one period).  <dE/dt> above reproduces
        // dipoleEccentricityFactor() exactly -- an independent check that
        // this is the same physics already used for the Larmor comparison,
        // not a new assumption -- and gamma cancels in the ratio, so k does
        // NOT depend on which reaction model (coherent or individual LL)
        // is active, only on the force law's r-dependence:
        //
        //   k = [<dL/dt>/L] / [d|E|/dt / |E|] = -(1-e^2)/(2+e^2)
        //
        // Checked against the (well-conditioned, e>0.02) part of the
        // measured seed-42 trajectory: matches the empirical k to 1-2e-4
        // absolute at every checkpoint (e.g. e^2=0.063: measured -0.454356,
        // formula -0.454154), i.e. the old measurement was already reading
        // this same formula, just with noise on top that grew as e shrank.
        // e^2 here is exactly the periapsis discriminant, so no new state is
        // needed to evaluate it.
        const double eccentricitySquared=std::max(0.0,1.0
            +2.0*elements.specificEnergy*elements.specificAngularMomentum
                *elements.specificAngularMomentum
                /(attractionParameter*attractionParameter));
        const double angularExponent=
            -(1.0-eccentricitySquared)/(2.0+eccentricitySquared);
        radiatedEnergyTotal+=
            (updatedEnergyMagnitude-energyMagnitude)*reducedMass;
        elements.specificEnergy=-updatedEnergyMagnitude;
        elements.specificAngularMomentum*=
            std::pow(energyGrowth,angularExponent);
        simulatedTimeTotal+=measuredElapsed*static_cast<double>(orbitsToSkip)
            *(1.0-0.5*jumpParameter);
        revolutionsTotal+=static_cast<double>(orbitsToSkip);
        firstDipole=run.finalState.firstDipole;
        secondDipole=run.finalState.secondDipole;

        if(!(elements.specificEnergy<0.0)||!std::isfinite(elements.specificEnergy)
           ||!std::isfinite(elements.specificAngularMomentum)) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal;
            return result;
        }
    }
    result.calibrationOutcome=SimulationOutcome::ObservationLimit;
    result.calibrationSeconds=simulatedTimeTotal;
    return result;
}

std::vector<CremCollapseEstimate> runCremCollapseExperiment(
    std::uint64_t masterSeed,int selectedPhenomenon,int runCount,
    double wallClockBudgetSeconds) {
    // Keep the safety boundary valid for non-CLI callers too.  In particular,
    // comparisons against NaN or +infinity never become true and would disable
    // every wall-clock stop inside estimateCremCollapse().
    if(!(wallClockBudgetSeconds>0.0)
       ||!std::isfinite(wallClockBudgetSeconds)) {
        throw std::invalid_argument(
            "--crem-wallclock-budget-s must be finite and positive");
    }
    std::vector<CremCollapseEstimate> estimates(static_cast<size_t>(runCount));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::mutex outputMutex;
    const int workerCount=std::min(runCount,
        static_cast<int>(std::max(1u,std::thread::hardware_concurrency())));
    std::cout<<"Running "<<runCount<<" CREM collapse calibrations on "
             <<workerCount<<" worker"<<(workerCount==1?"":"s")
             <<" (per-event wall-clock budget "<<wallClockBudgetSeconds<<" s).\n";
    const auto worker=[&]() {
        while(true) {
            const int index=nextIndex.fetch_add(1);
            if(index>=runCount) break;
            estimates[static_cast<size_t>(index)]=estimateCremCollapse(
                splitMix64(masterSeed+static_cast<std::uint64_t>(index)),
                selectedPhenomenon,wallClockBudgetSeconds);
            const int done=completed.fetch_add(1)+1;
            if(done%10==0||done==runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout<<"CREM calibrations: "<<done<<"/"<<runCount<<'\n';
            }
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    for(int index=0;index<workerCount;++index) workers.emplace_back(worker);
    for(std::thread& thread:workers) thread.join();
    return estimates;
}
