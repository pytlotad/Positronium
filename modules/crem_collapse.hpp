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
inline constexpr double classicalInspiralCoefficient() {
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

// Unperturbed Kepler period at the given specific energy; used only to size
// the one-orbit measurement window, not to compute the reported lifetime.
double osculatingPeriod(double specificEnergy,double attractionParameter) {
    const double semiMajorAxis=-attractionParameter/(2.0*specificEnergy);
    return 2.0*pi*std::sqrt(semiMajorAxis*semiMajorAxis*semiMajorAxis
        /attractionParameter);
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
    const double periapsis=osculatingPeriapsis(elements,attractionParameter);
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
    constexpr double reducedMass=firstMass*secondMass
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
    // finalApproachMultiple*chargeCloudRestRadius, so the reference must be
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
        result.finalPeriodSeconds=osculatingPeriod(
            elements.specificEnergy,attractionParameter);
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
        // Stop the orbit-averaged phase a bit before the exact boundary
        // rather than chasing it there.  With t~a^3 near collapse (P~a^-4),
        // covering the remaining stretch from 10*collisionRadius down to
        // collisionRadius accounts for only (1/10)^3 = 0.1% of the total
        // collapse time, well under the 3% per-jump tolerance already
        // accepted elsewhere -- but the one-period measurement window this
        // loop relies on stops being a good approximation exactly in that
        // last stretch (eccentricity and the dipole barrier both grow
        // quickly), so pushing further trades a well-understood, negligible
        // truncation for an unreliable measurement.
        constexpr double finalApproachMultiple=10.0;
        if(periapsis<=finalApproachMultiple*chargeCloudRestRadius) {
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
        const double period=osculatingPeriod(
            elements.specificEnergy,attractionParameter);
        SimulationOptions measureOptions;
        measureOptions.collectFrames=false;
        measureOptions.observationTime=period;
        measureOptions.terminalSeparation=chargeCloudRestRadius;
        // The secular update below is driven entirely by mechanical
        // quantities -- conservativeParticleEnergy() and the orbital angular
        // momentum read off the endpoint states -- and radiatedEnergyTotal is
        // accumulated from those same differences.  Nothing here reads the
        // flux bookkeeping, so its quadrature is pure overhead on the two
        // integrations every checkpoint pays for.
        measureOptions.radiatedEnergyBookkeeping=false;
        // The budget is shared with the whole estimate: a single
        // measurement orbit can itself run into the stiff region near the
        // boundary, so it must be interruptible too.
        measureOptions.stopRequested=[&]() {
            return wallClockSpent()>wallClockBudgetSeconds;
        };
        const MechanicalTrajectoryResult run=runMechanicalTrajectory(
            measurementState,period,chargeCloudRestRadius,measureOptions,
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
            background=runMechanicalTrajectory(
                measurementState,period,chargeCloudRestRadius,measureOptions,
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
            const OsculatingElements realDelta=measuredDelta(run.finalState);
            const OsculatingElements backgroundDelta=backgroundFor(realDelta);
            result.lifetimeSeconds=simulatedTimeTotal+run.elapsedTime;
            radiatedEnergyTotal+=
                (backgroundDelta.specificEnergy-realDelta.specificEnergy)*reducedMass;
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
        const double deltaEnergyPerOrbit=
            realDelta.specificEnergy-backgroundDelta.specificEnergy;
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
                     <<" (raw "<<realDelta.specificEnergy/elements.specificEnergy
                     <<", background "<<backgroundDelta.specificEnergy
                        /elements.specificEnergy
                     <<") dL/L="<<deltaAngularMomentumPerOrbit
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
            orbitsToSkip=std::clamp(static_cast<int>(
                maximumJumpParameter*energyMagnitude/(1.5*lossPerOrbit)),
                1,maxOrbitsSkippedAtOnce);
        }
        const double jumpParameter=std::min(
            1.5*static_cast<double>(orbitsToSkip)*lossPerOrbit/energyMagnitude,
            maximumJumpParameter);
        const double energyGrowth=std::pow(1.0-jumpParameter,-2.0/3.0);
        const double updatedEnergyMagnitude=energyMagnitude*energyGrowth;
        // k from this checkpoint's own measurement, guarded against a
        // vanishing denominator.
        const double relativeEnergyStep=lossPerOrbit/energyMagnitude;
        double angularExponent=0.0;
        if(relativeEnergyStep>0.0&&elements.specificAngularMomentum!=0.0) {
            angularExponent=(deltaAngularMomentumPerOrbit
                /elements.specificAngularMomentum)/relativeEnergyStep;
        }
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
