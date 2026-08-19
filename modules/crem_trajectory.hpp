#pragma once

// CREM trajectory layer: the observable Frame taken from a State, the shared
// mechanical integration loop that every experiment drives, and the sampler
// that prepares a positronium-scale initial condition and classifies what it
// turns into.
//
// Textual module, included from inside the anonymous namespace of
// positronium.cpp and inside the production #ifndef.  Contains no ROOT.

Frame makeFrame(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const double firstKinetic = (gamma(s.firstVelocity) - 1.0) * firstMass * c*c;
    const double secondKinetic = (gamma(s.secondVelocity) - 1.0) * secondMass * c*c;
    const double coulombPotential = -pairCoulombStrength * geometry.inverseDistance;
    const double dipolePotential = -dot(s.firstDipole,
        regularizedDipoleField(geometry.firstMinusSecond, s.secondDipole));
    const double darwinEnergy = darwinInteractionEnergy(s);
    const MutualForces forces = allExternalForces(s);
    const Vec3 firstAcceleration = relativisticAcceleration(s.firstVelocity, forces.first, firstMass);
    const Vec3 secondAcceleration = relativisticAcceleration(s.secondVelocity, forces.second, secondMass);
    // Low-velocity diagnostic proxy for the sum of the two individual Schott
    // near-field energies.  Cross terms belong to the mutual retarded field
    // and must not be inserted again through a collective dipole expression.
    const double schottEnergy =
        -(firstCharge*firstCharge*dot(firstAcceleration, s.firstVelocity)
         + secondCharge*secondCharge*dot(secondAcceleration, s.secondVelocity))
        / (6.0 * pi * epsilon0 * c*c*c);
    // Noether energy of the conservative approximate action.  The q-mu term
    // is homogeneous of degree one in both velocities and therefore cancels
    // in the Legendre transform.
    const double conservativeNoetherEnergy=conservativeParticleEnergy(s);
    // One evaluation feeds all three Noether quantities below.
    const CanonicalMomenta canonical = canonicalMomenta(s);
    return {s.firstPosition,s.secondPosition,
            s.firstProperDipole.squaredNorm()>0.0
                ?s.firstProperDipole:s.firstDipole,
            s.secondProperDipole.squaredNorm()>0.0
                ?s.secondProperDipole:s.secondDipole,
            noetherMomentum(canonical),
            noetherAngularMomentum(s, canonical),
            s.radiatedMomentum, s.radiatedAngularMomentum,
            canonicalMomentumScale(canonical),
            s.time, geometry.distance, s.radiatedEnergy,
            conservativeNoetherEnergy,
            schottEnergy,
            firstKinetic + 0.5 * (coulombPotential + dipolePotential
                + darwinEnergy + s.dipoleConstraintEnergy),
            secondKinetic + 0.5 * (coulombPotential + dipolePotential
                + darwinEnergy + s.dipoleConstraintEnergy),
            s.boundFieldEnergy,s.reactionEnergyMismatch,
            s.boundFieldMomentum,s.boundFieldAngularMomentum,
            s.reactionMomentumMismatch,s.reactionAngularMomentumMismatch};
}

const char* phenomenonName(Phenomenon phenomenon) {
    switch (phenomenon) {
        case Phenomenon::DirectCollision: return "Direct collision";
        case Phenomenon::Scattering: return "Scattering";
        case Phenomenon::ParaPositronium: return "Para-positronium";
        case Phenomenon::OrthoPositronium: return "Ortho-positronium";
    }
    return "Unknown";
}

struct MechanicalTrajectoryResult {
    State finalState;
    SimulationOutcome outcome=SimulationOutcome::NumericalFailure;
    double minimumSeparation=std::numeric_limits<double>::quiet_NaN();
    double elapsedTime=0.0;
    double finalRadiatedEnergy=0.0;
    double maximumBeta=0.0;
    std::vector<Frame> frames;
};

// Integrates `initial` forward with the complete CREM engine (retardation,
// Darwin, dipole coupling, radiation reaction under whichever
// --radiation-reaction model is active) until separation reaches
// trajectoryCutoff (ReachedCutoff), observationTime elapses or
// options.stopRequested() fires (both ObservationLimit), or the state turns
// non-finite (NumericalFailure).  This is the shared mechanical core behind
// simulate(), which samples its own random initial condition, and the
// secular CREM collapse estimator, which reconstructs a fresh osculating
// state after each measured orbit instead of sampling one.
MechanicalTrajectoryResult runMechanicalTrajectory(State s,
                                                    double observationTime,
                                                    double trajectoryCutoff,
                                                    const SimulationOptions& options,
                                                    ChargeRadiationReactionModel
                                                        reactionModel=gRadiationReactionModel) {
    const double reducedMass = firstMass * secondMass / (firstMass + secondMass);
    // Establish the covariant dipoles BEFORE anything samples the state.
    // A zero properDipole is the "not initialized yet" sentinel that
    // synchronizeCovariantDipoles() fills in on its first call, which used to
    // happen inside the first integration step.  Any frame interpolated
    // between the pre-step state (proper dipole still zero) and the post-step
    // state therefore came out scaled by the interpolation fraction:
    // interpolateDipole() carries the norm linearly, and a linear ramp from
    // zero is exactly the fraction.  The endpoint diagnostic then reported a
    // spurious |mu| drift of exactly 1 - frameInterval/dt on that one frame.
    // It stayed hidden only while dt happened to be smaller than the frame
    // interval, so the fraction clamped to 1; widening the orbit made dt
    // larger than the frame interval and the artefact appeared at 47%, 30%
    // and 86% for the three observation windows, matching that formula to
    // better than 0.1%.
    synchronizeCovariantDipoles(s);
    const int frameCount = std::max(2, options.frameCount);
    std::vector<Frame> frames;
    if (options.collectFrames) frames.reserve(frameCount + 1);
    const double frameInterval = observationTime / (frameCount - 1);
    double nextFrame = frameInterval;
    double minimumSeparation = separation(s);
    double maximumBeta = std::max(s.firstVelocity.norm(), s.secondVelocity.norm()) / c;
    double elapsedTime = s.time;
    double finalRadiatedEnergy = s.radiatedEnergy;
    SimulationOutcome outcome = SimulationOutcome::NumericalFailure;
    ClassicalTrajectoryEngine trajectory(s,
        {.relativeTolerance=1.0e-5,.maximumDepth=12,
         .reactionModel=reactionModel,
         .computeOutwardFlux=options.radiatedEnergyBookkeeping});
    if (options.collectFrames) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }

    bool reachedObservationCeiling=false;
    while (s.time < observationTime && separation(s) > trajectoryCutoff
           &&!(options.stopRequested&&options.stopRequested())) {
        // Resolve each instantaneous orbit well enough to keep numerical
        // energy drift below the physical radiation loss.
        const State beforeStep = s;
        const double r = separation(beforeStep);
        const double omega = std::sqrt(pairCoulombStrength / (reducedMass * r*r*r));
        const double remaining = observationTime - s.time;
        // A remaining budget below the floating-point resolution of s.time
        // can never actually change it: state.time += dt becomes a no-op,
        // and the loop would spin at a vanishing dt forever (or until an
        // external wall-clock budget intervenes) instead of ending cleanly.
        // Same failure mode, and same fix, as the clock-resolution guard in
        // the Experiment 5 interaction sampler.
        const double clockResolution = std::max(std::abs(s.time),
            std::abs(observationTime)) * std::numeric_limits<double>::epsilon() * 8.0;
        if (remaining <= clockResolution) { reachedObservationCeiling=true; break; }
        const double dt = std::min({5.0e-18, 2.0 * pi / (128.0 * omega), remaining});
        if (!(dt > 0.0) || !std::isfinite(dt)) break;
        if(!trajectory.advance(s,dt)) break;
        const double currentSeparation = separation(s);
        if (!(currentSeparation > 0.0) || !std::isfinite(currentSeparation)) break;
        if (options.stepReady) options.stepReady(s);

        const bool reachedCutoff = currentSeparation <= trajectoryCutoff;
        const double crossingFraction = reachedCutoff
            ? separationCrossingFraction(beforeStep,s,trajectoryCutoff) : 1.0;
        const double validEndTime = beforeStep.time + crossingFraction * dt;
        if (options.collectFrames) {
            while (nextFrame <= validEndTime
                   && static_cast<int>(frames.size()) < frameCount) {
                const double sampleFraction = std::clamp(
                    (nextFrame - beforeStep.time) / dt, 0.0, crossingFraction);
                const State sampledState = interpolateState(beforeStep, s, sampleFraction);
                if (separation(sampledState) > trajectoryCutoff) {
                    frames.push_back(makeFrame(sampledState));
                    if (options.frameReady) options.frameReady(frames.back());
                }
                nextFrame += frameInterval;
            }
        }

        if (reachedCutoff) {
            // Locate the terminal event linearly inside the last, already
            // very small adaptive step.  This avoids reporting a numerical
            // overshoot below the boundary of the point-particle model.
            elapsedTime = beforeStep.time + crossingFraction * dt;
            finalRadiatedEnergy = beforeStep.radiatedEnergy + crossingFraction
                * (s.radiatedEnergy - beforeStep.radiatedEnergy);
            const Vec3 eventFirstVelocity = beforeStep.firstVelocity
                + (s.firstVelocity - beforeStep.firstVelocity) * crossingFraction;
            const Vec3 eventSecondVelocity = beforeStep.secondVelocity
                + (s.secondVelocity - beforeStep.secondVelocity) * crossingFraction;
            minimumSeparation = std::min(minimumSeparation, trajectoryCutoff);
            maximumBeta = std::max(maximumBeta,
                std::max(eventFirstVelocity.norm(), eventSecondVelocity.norm()) / c);
            State eventState=interpolateState(beforeStep,s,crossingFraction);
            if(options.collectFrames) {
                const double balanceEnergy=conservativeParticleEnergy(beforeStep)
                    +beforeStep.radiatedEnergy+beforeStep.boundFieldEnergy;
                const Vec3 balanceMomentum=noetherMomentum(beforeStep)
                    +beforeStep.radiatedMomentum+beforeStep.boundFieldMomentum;
                const Vec3 balanceAngularMomentum=noetherAngularMomentum(beforeStep)
                    +beforeStep.radiatedAngularMomentum
                    +beforeStep.boundFieldAngularMomentum;
                eventState.boundFieldEnergy=balanceEnergy
                    -conservativeParticleEnergy(eventState)
                    -eventState.radiatedEnergy;
                eventState.boundFieldMomentum=balanceMomentum
                    -noetherMomentum(eventState)-eventState.radiatedMomentum;
                eventState.boundFieldAngularMomentum=balanceAngularMomentum
                    -noetherAngularMomentum(eventState)
                    -eventState.radiatedAngularMomentum;
                if(frames.empty()||frames.back().time<eventState.time) {
                    frames.push_back(makeFrame(eventState));
                    if(options.frameReady) options.frameReady(frames.back());
                }
            }
            outcome = SimulationOutcome::ReachedCutoff;
            s = eventState;
            break;
        }

        minimumSeparation = std::min(minimumSeparation, currentSeparation);
        maximumBeta = std::max(maximumBeta,
            std::max(s.firstVelocity.norm(), s.secondVelocity.norm()) / c);
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;
    }
    // Reaching the simulated-time ceiling and an external stop request (an
    // exhausted wall-clock budget in the CREM collapse estimator, or the
    // visual window closing) are both "we stopped watching before the
    // boundary", not a numerical breakdown -- both are reported the same way.
    if (outcome != SimulationOutcome::ReachedCutoff && isFinite(s)
        && (s.time >= observationTime || reachedObservationCeiling
            || (options.stopRequested && options.stopRequested()))) {
        outcome = SimulationOutcome::ObservationLimit;
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;
    }
    if (options.collectFrames && outcome == SimulationOutcome::ObservationLimit
        && (frames.empty() || frames.back().time < s.time)) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }
    return {s, outcome, minimumSeparation, elapsedTime, finalRadiatedEnergy,
            maximumBeta, std::move(frames)};
}

SimulationResult simulate(std::uint64_t seed, int selectedPhenomenon,
                          SimulationOptions options = {}) {
    const double reducedMass = firstMass * secondMass / (firstMass + secondMass);
    // Every scenario starts at the positronium scale, not hydrogen's.  The
    // unbound scenarios inherit it too: their speeds are quoted as multiples
    // of the circular speed at the starting separation, so the trajectories
    // keep their character and only the overall scale moves.
    // Bohr radius OF THE PAIR, hbar^2/(mu k |q1 q2|) -- not positronium's, and
    // certainly not hydrogen's.  105.8 pm for e+e-, 512 fm for mu+mu-, 57.6 fm
    // for p+pbar.  The orbit carrying L = hbar sits here, which is also where
    // the binding energy matches the measured one.
    const double initialSeparation = pairBohrRadius(activePair);
    const double circularSpeed =
        std::sqrt(pairCoulombStrength / (reducedMass * initialSeparation));
    const double escapeSpeed = std::sqrt(2.0) * circularSpeed;
    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> unitRandom(0.0, 1.0);
    std::uniform_real_distribution<double> signedRandom(-1.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, 2.0*pi);
    const auto randomDirection = [&]() {
        const double z = signedRandom(random);
        const double phi = azimuth(random);
        const double radial = std::sqrt(1.0 - z*z);
        return Vec3{radial*std::cos(phi), radial*std::sin(phi), z};
    };

    State s;
    s.firstPosition =
        {initialSeparation * secondMass / (firstMass + secondMass), 0, 0};
    s.secondPosition =
        {-initialSeparation * firstMass / (firstMass + secondMass), 0, 0};

    // The selected branch chooses a physically useful sampling range while
    // all actual initial values inside that range remain random.
    // Menu order is para, ortho, direct collision, scattering. Internally the
    // samplers retain the order direct, scattering, para, ortho.
    const std::array<int, 5> scenarioForMenuChoice = {0, 2, 3, 0, 1};
    const int sampledScenario = scenarioForMenuChoice[selectedPhenomenon];
    double radialSpeed = 0.0;
    double tangentialSpeed = 0.0;
    if (sampledScenario == 0) {
        radialSpeed = -escapeSpeed * (0.35 + 0.40 * unitRandom(random));
        tangentialSpeed = circularSpeed * (0.001 + 0.004 * unitRandom(random));
    } else if (sampledScenario == 1) {
        const double speed = escapeSpeed * (1.05 + 0.35 * unitRandom(random));
        const double tangentialFraction = 0.35 + 0.35 * unitRandom(random);
        tangentialSpeed = speed * tangentialFraction;
        radialSpeed = -speed * std::sqrt(1.0 - tangentialFraction*tangentialFraction);
    } else {
        // Bound scenarios are now centred on the circular orbit at a_Ps, which
        // is the L = hbar state: at this separation L = f * hbar exactly, with
        // f the tangential speed in units of the circular speed.  The old band
        // f in [0.72, 0.97] sat entirely BELOW the ground state, so every
        // sampled orbit was more tightly bound than positronium actually is.
        // The band now straddles it, leaving a comparable spread.
        radialSpeed = circularSpeed * (-0.10 + 0.20 * unitRandom(random));
        tangentialSpeed = circularSpeed * (0.88 + 0.24 * unitRandom(random));
    }
    const Vec3 relativeVelocity{radialSpeed, tangentialSpeed, 0.0};
    s.firstVelocity = relativeVelocity * (secondMass / (firstMass + secondMass));
    s.secondVelocity = relativeVelocity * (-firstMass / (firstMass + secondMass));

    s.firstDipole = randomDirection() * firstMagneticMoment;
    do {
        s.secondDipole = randomDirection() * secondMagneticMoment;
    } while ((sampledScenario == 2 && dot(s.firstDipole, s.secondDipole)
                                  / (firstMagneticMoment*secondMagneticMoment) < 0.5)
          || (sampledScenario == 3 && dot(s.firstDipole, s.secondDipole)
                                  / (firstMagneticMoment*secondMagneticMoment) >= 0.5));

    const double relativeEnergy = 0.5 * reducedMass * relativeVelocity.squaredNorm()
                                - pairCoulombStrength / initialSeparation;
    const double orbitalAngularMomentum =
        reducedMass * cross(Vec3{initialSeparation, 0, 0}, relativeVelocity).norm();
    const double specificEnergy = relativeEnergy / reducedMass;
    const double specificAngularMomentum = orbitalAngularMomentum / reducedMass;
    const double attractionParameter = pairCoulombStrength / reducedMass;
    const double eccentricity = std::sqrt(std::max(0.0, 1.0 + 2.0 * specificEnergy
        * specificAngularMomentum*specificAngularMomentum
        / (attractionParameter*attractionParameter)));
    const double predictedClosestApproach = specificAngularMomentum == 0.0 ? 0.0
        : specificAngularMomentum*specificAngularMomentum
          / (attractionParameter * (1.0 + eccentricity));
    const double dipoleAlignment = dot(s.firstDipole, s.secondDipole)
                                 / (firstMagneticMoment*secondMagneticMoment);
    Phenomenon phenomenon;
    if (radialSpeed < 0.0 && predictedClosestApproach < nuclearCutoff) {
        phenomenon = Phenomenon::DirectCollision;
    } else if (relativeEnergy >= 0.0) {
        phenomenon = Phenomenon::Scattering;
    } else if (dipoleAlignment >= 0.5) {
        phenomenon = Phenomenon::ParaPositronium;
    } else {
        phenomenon = Phenomenon::OrthoPositronium;
    }

    const double defaultObservationTime = phenomenon == Phenomenon::DirectCollision ? 4.0e-16
                                        : phenomenon == Phenomenon::Scattering ? 2.0e-15
                                        : 1.50e-15;
    const double observationTime = options.observationTime > 0.0
                                 ? options.observationTime : defaultObservationTime;
    // Visual mode is a finite-resolution preview. Its looser local tolerance
    // remains substantially smaller than a screen pixel at the displayed
    // scale, while endpoint diagnostics continue to expose accumulated drift.
    const bool directCollision=phenomenon==Phenomenon::DirectCollision;
    // The visual point-particle picture cannot resolve the interior of the
    // configured charge cloud. Stop a direct-collision animation at that
    // model boundary; beam statistics retain the smaller nuclear cutoff.
    const double trajectoryCutoff=options.terminalSeparation>0.0
        ?options.terminalSeparation
        :(directCollision?chargeCloudRestRadius:nuclearCutoff);
    MechanicalTrajectoryResult run = runMechanicalTrajectory(
        s, observationTime, trajectoryCutoff, options);
    const double timeToCutoff = run.outcome == SimulationOutcome::ReachedCutoff
                              ? run.elapsedTime : std::numeric_limits<double>::infinity();
    return {std::move(run.frames), {relativeEnergy, orbitalAngularMomentum,
            predictedClosestApproach, dipoleAlignment, timeToCutoff, phenomenon, seed},
            run.outcome, run.minimumSeparation, run.elapsedTime,
            run.finalRadiatedEnergy, run.maximumBeta};
}
