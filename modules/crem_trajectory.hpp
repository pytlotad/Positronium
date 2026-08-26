#pragma once

// CREM trajectory layer: the observable Frame taken from a State, the shared
// mechanical integration loop that every experiment drives, and the sampler
// that prepares a positronium-scale initial condition and classifies what it
// turns into.
//
// Textual module, included from inside the anonymous namespace of
// positronium.cpp and inside the production #ifndef.  Contains no ROOT.

Frame makeFrame(const State& s) {
    // TRUE geometry: only for the reported radius below, which must say
    // where the pair actually is, not where the force laws pretend it is.
    const PairGeometry geometry = pairGeometry(s);
    const double firstKinetic=kineticEnergy(s.firstVelocity,firstMass);
    const double secondKinetic=kineticEnergy(s.secondVelocity,secondMass);
    // Clamped geometry for the energy terms: firstMechanicalEnergy +
    // secondMechanicalEnergy must sum to conservativeNoetherEnergy below,
    // which already reads the clamped potential through
    // conservativeParticleEnergy() -- using the true geometry here instead
    // would silently break that identity below the barrier.
    const PairGeometry clampedGeometry = clampedPairGeometry(s);
    const double coulombPotential =
        -pairCoulombStrength * clampedGeometry.inverseDistance;
    const double dipolePotential = regularizedDipoleInteractionEnergy(
        clampedGeometry.firstMinusSecond,s.firstDipole,s.secondDipole);
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

// Local copy of positronium.cpp's splitMix64: that one is defined AFTER this
// header is included (line ~1213 vs. this header's ~1174), so it is not yet
// visible here.  Same well-known bit-mixer, just under its own name to avoid
// masking the later declaration.
std::uint64_t stochasticPhotonHash64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// Draws one Exp(1) variate (mean 1) from the running 64-bit stream state,
// advancing it in place.  This is the "next photon's hazard threshold" in
// the standard thinning/inverse-CDF construction of an inhomogeneous
// Poisson process: accumulate hazard = integral(rate dt) and fire whenever
// it crosses a freshly-drawn Exp(1) threshold.
double drawExponentialUnit(std::uint64_t& streamState) {
    streamState = stochasticPhotonHash64(streamState);
    // Upper 53 bits -> a double uniform in (0,1]; excluding 0 keeps log finite.
    const double uniform =
        1.0 - static_cast<double>(streamState >> 11) * (1.0 / 9007199254740992.0);
    return -std::log(uniform);
}

// Draws one uniform variate in [0,1) from the same stream, advancing it in
// place.  Same bit construction as drawExponentialUnit's own intermediate
// uniform, just returned directly instead of being fed through -log().
double drawUniformUnit(std::uint64_t& streamState) {
    streamState = stochasticPhotonHash64(streamState);
    return static_cast<double>(streamState >> 11) * (1.0 / 9007199254740992.0);
}

// Direction of one photon from a circularly rotating E1 dipole.  Relative to
// the orbital normal its density is proportional to 1+cos^2(theta).  The
// cubic inverse below is the same closed-form sampler used by the secular
// collapse path; keeping it here makes the mechanical and skipped-orbit paths
// consume the same physical angular law.
Vec3 sampleRotatingDipolePhotonDirection(const Vec3& orbitalNormal,
                                         std::uint64_t& streamState) {
    Vec3 axis=orbitalNormal;
    const double axisNorm=axis.norm();
    axis=axisNorm>1.0e-300&&std::isfinite(axisNorm)
        ?axis*(1.0/axisNorm):Vec3{0.0,0.0,1.0};
    const double q=4.0-8.0*drawUniformUnit(streamState);
    const double root=std::sqrt(0.25*q*q+1.0);
    const double cosine=std::cbrt(-0.5*q+root)+std::cbrt(-0.5*q-root);
    const double sine=std::sqrt(std::max(0.0,1.0-cosine*cosine));
    const Vec3 seed=std::abs(axis.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
    Vec3 first=cross(axis,seed);
    first=first*(1.0/first.norm());
    const Vec3 second=cross(axis,first);
    const double azimuth=2.0*pi*drawUniformUnit(streamState);
    return axis*cosine
        +(first*std::cos(azimuth)+second*std::sin(azimuth))*sine;
}

struct StochasticPhotonRecoil {
    bool emitted=false;
    double energy=0.0;       // photon energy in the pre-emission pair COM
    Vec3 direction;          // photon direction in that frame
};

// Emit a photon in the instantaneous centre-of-momentum frame and reconstruct
// the remaining two-particle state exactly on both mass shells.  If P is the
// incoming pair four-momentum and k=(E_gamma/c,E_gamma*n/c), the result obeys
// p1'+p2'+k=P (up to floating-point roundoff).  The residual pair has
// invariant energy W'=sqrt(W^2-2 W E_gamma); its internal momentum keeps the
// pre-emission direction, while its COM carries the full recoil -k.
//
// A draw that would put W' below (m1+m2)c^2 is kinematically forbidden and is
// rejected without modifying the state.  This replaces the old capped kick,
// which could remove all relative kinetic energy while still crediting an
// unrelated photon and did not conserve the pair+photon momentum.
StochasticPhotonRecoil applyStochasticDipolePhoton(
    State& s,double photonEnergy,const Vec3& photonDirection) {
    StochasticPhotonRecoil result;
    if(!(photonEnergy>0.0)||!std::isfinite(photonEnergy)) return result;
    const double directionNorm=photonDirection.norm();
    if(!(directionNorm>0.0)||!std::isfinite(directionNorm)) return result;
    const Vec3 direction=photonDirection*(1.0/directionNorm);
    const auto firstLab=two_body::fourMomentumFromVelocity(
        s.firstVelocity,firstMass);
    const auto secondLab=two_body::fourMomentumFromVelocity(
        s.secondVelocity,secondMass);
    if(!firstLab.valid()||!secondLab.valid()) return result;
    const double incomingEnergy=firstLab.energy+secondLab.energy;
    const Vec3 incomingMomentum=firstLab.momentum+secondLab.momentum;
    const Vec3 incomingComVelocity=incomingMomentum*(c*c/incomingEnergy);
    if(!(incomingComVelocity.squaredNorm()<c*c)) return result;
    const auto firstCom=two_body::boostFourMomentum(
        firstLab,incomingComVelocity*(-1.0));
    const auto secondCom=two_body::boostFourMomentum(
        secondLab,incomingComVelocity*(-1.0));
    if(!firstCom.valid()||!secondCom.valid()) return result;

    const long double W=static_cast<long double>(firstCom.energy)
        +static_cast<long double>(secondCom.energy);
    const long double Eg=photonEnergy;
    const long double cL=c;
    const long double rest1=static_cast<long double>(firstMass)*cL*cL;
    const long double rest2=static_cast<long double>(secondMass)*cL*cL;
    const long double W2=W*W-2.0L*W*Eg;
    const long double threshold=(rest1+rest2)*(rest1+rest2);
    if(!(W2>=threshold)||!std::isfinite(static_cast<double>(W2))) return result;
    const long double residualW=std::sqrt(W2);
    const long double firstEnergy=(W2+rest1*rest1-rest2*rest2)
        /(2.0L*residualW);
    const long double momentumEnergySquared=std::max(0.0L,
        (firstEnergy-rest1)*(firstEnergy+rest1));
    const double relativeMomentum=static_cast<double>(
        std::sqrt(momentumEnergySquared)/cL);
    Vec3 relativeDirection=firstCom.momentum-secondCom.momentum;
    const double relativeNorm=relativeDirection.norm();
    if(!(relativeNorm>1.0e-300)||!std::isfinite(relativeNorm)) {
        relativeDirection=cross(direction,Vec3{0,0,1});
        if(relativeDirection.squaredNorm()<1.0e-24)
            relativeDirection=cross(direction,Vec3{0,1,0});
    }
    relativeDirection=relativeDirection*(1.0/relativeDirection.norm());
    const auto firstResidualRest=two_body::fourMomentumFromMomentum(
        relativeDirection*relativeMomentum,firstMass);
    const auto secondResidualRest=two_body::fourMomentumFromMomentum(
        relativeDirection*(-relativeMomentum),secondMass);
    const double residualEnergyInOldCom=static_cast<double>(W-Eg);
    const Vec3 residualMomentumInOldCom=direction*(-photonEnergy/c);
    const Vec3 residualComVelocity=residualMomentumInOldCom
        *(c*c/residualEnergyInOldCom);
    const auto firstOldCom=two_body::boostFourMomentum(
        firstResidualRest,residualComVelocity);
    const auto secondOldCom=two_body::boostFourMomentum(
        secondResidualRest,residualComVelocity);
    const auto firstFinal=two_body::boostFourMomentum(
        firstOldCom,incomingComVelocity);
    const auto secondFinal=two_body::boostFourMomentum(
        secondOldCom,incomingComVelocity);
    if(!firstFinal.valid()||!secondFinal.valid()) return result;
    const Vec3 firstVelocity=two_body::velocityFromFourMomentum(firstFinal);
    const Vec3 secondVelocity=two_body::velocityFromFourMomentum(secondFinal);
    if(!isFinite(firstVelocity)||!isFinite(secondVelocity)) return result;
    s.firstVelocity=firstVelocity;
    s.secondVelocity=secondVelocity;
    result.emitted=true;
    result.energy=photonEnergy;
    result.direction=direction;
    return result;
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
         .compositionOrder=gIntegratorOrder,
         .reactionModel=reactionModel,
         .computeOutwardFlux=options.radiatedEnergyBookkeeping});
    if (options.collectFrames) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }

    // Poisson-process bookkeeping for
    // ChargeRadiationReactionModel::stochasticElectricDipole, unused (and
    // costing nothing) for every other model.  Seeded from bits of the
    // INITIAL state, which is already seed-derived and unique per
    // trajectory/measurement-orbit call upstream -- no new seed parameter
    // needs threading through this function's signature.  Local to this one
    // call: correct, because a stochastic reaction-force run always drives
    // the trajectory through exactly one runMechanicalTrajectory call, never
    // split across several the way the CREM collapse estimator's checkpoint
    // loop splits a full collapse into many short measurement orbits.
    std::uint64_t stochasticPhotonStream=[&]() {
        std::uint64_t bits=0;
        auto mix=[&](double value) {
            std::uint64_t word;
            std::memcpy(&word,&value,sizeof(word));
            bits=stochasticPhotonHash64(bits^word);
        };
        mix(s.firstPosition.x); mix(s.firstPosition.y); mix(s.firstPosition.z);
        mix(s.firstVelocity.x); mix(s.secondVelocity.y); mix(s.time);
        return bits;
    }();
    double stochasticHazard=0.0;
    double stochasticThreshold=drawExponentialUnit(stochasticPhotonStream);

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

        // Bank this step's TOTAL classical radiated power as hazard instead
        // of removing any of it as a continuous force.  One photon stream
        // carries every channel: the E1 charge power, the M1 magnetic dipole
        // power and the E2 charge quadrupole power are summed before the
        // division, so nothing in this mode radiates continuously.  The charge sector's gate sits in
        // particleMultipoleRadiation (chargeReaction is zero there); the
        // magnetic sector's two continuous sinks -- the reaction torque and
        // the dipoleConstraintEnergy drain -- are gated in
        // integrateElectrodynamicStep, or the M1 energy would leave both
        // continuously and as photons.
        //
        // photonEnergy = hbar*omega uses the SAME instantaneous orbital
        // frequency already computed above for step sizing: the natural,
        // already-established characteristic frequency of this orbit, not a
        // new definition invented for this mode.  Applying it to the M1
        // power too is the deliberate cost of a single unified channel: M1
        // is actually emitted near the spin-precession rate, which is
        // roughly four orders of magnitude BELOW omega_orb here, so its
        // quanta come out fewer and individually larger than that rate
        // would give.  What this does NOT change is the mean radiated
        // energy: the hazard rate is P/(hbar*omega) and each photon removes
        // hbar*omega, so the expected energy per unit time is P for ANY
        // choice of omega.  Only the granularity and the shot statistics
        // depend on it.
        //
        // The secular estimator in crem_collapse.hpp stays E1-only, and
        // cannot be otherwise: it integrates the hazard analytically across
        // skipped orbits from the osculating elements alone, which carry no
        // spin state for an M1 power to be reconstructed from.  That is not
        // a gap in practice -- the bound phenomena that path runs measure
        // M1 identically zero, because the pair's two moments enter the
        // coherent M1 amplitude as m1+m2 and cancel.
        if(reactionModel==ChargeRadiationReactionModel::stochasticElectricDipole) {
            const MutualForces stepForces=
                retardedExternalForces(s,trajectory.history());
            const ParticleMultipoleRadiation stepRadiation=
                particleMultipoleRadiation(s,stepForces,trajectory.history(),
                    false,reactionModel);
            const double photonEnergy=hbar*omega;
            if(photonEnergy>0.0) {
                // Every channel the model computes, summed before the
                // division: E1 charge, M1 magnetic dipole, E2 charge
                // quadrupole.  Multipole powers are orthogonal at this
                // order, so this is a sum and not a double count.  E2 is
                // identically zero for every mass-symmetric pair (see
                // needsQuadrupolePower in electrodynamics.hpp) and carries
                // the channel only for asymmetric ones such as p+e-.
                //
                // WHERE THE LADDER STOPS, and why here rather than one rung
                // further.  The next order is M2/E3, not M1/E2: the ORBITAL
                // magnetic dipole carries the same kappa as the quadrupole
                // above, so it vanishes with it for a mass-symmetric pair,
                // and for any pair it is m = (kappa/2)(d x d-dot), which is
                // proportional to the conserved orbital angular momentum --
                // constant, hence m-ddot = 0 and no radiation, to leading
                // order, whatever kappa is.  (The M1 term summed above is
                // the INTRINSIC moment's, a different object.)
                //
                // Which electric multipoles survive at all is fixed by a
                // parity rule, worth stating once because it answers the
                // whole infinite family rather than one rung.  Writing
                // r1 = (m2/M) d and r2 = -(m1/M) d, the l-th moment carries
                //
                //     kappa_l = [ q1 m2^l + (-1)^l q2 m1^l ] / M^l,
                //
                // which for a MASS-SYMMETRIC pair collapses to
                // [q1 + (-1)^l q2]/2^l: proportional to q1+q2 = 0 for even l,
                // and to q1-q2 = -2e for odd l.  So a neutral symmetric pair
                // has NO even electric multipole and every odd one.  Measured
                // for e+e-: kappa_1 = -1.000 e, kappa_2 = 0, kappa_3 =
                // -0.250 e, kappa_4 = 0, kappa_5 = -0.0625 e.
                //
                // E3 therefore does NOT vanish for positronium the way E2
                // does -- it is the leading correction to E1 there, together
                // with M2.  Both sit at beta^4 relative to E1 (E_l scales as
                // beta^(2(l-1)), M_l as beta^(2l)), one order beyond the
                // usual beta^2, M2 paying both the magnetic and the
                // one-l-higher penalty.  At the pair Bohr radius beta = alpha
                // exactly, so that is 2.84e-9 -- against which the model's
                // own radiation measurement has a floor of 6.9e-5 (the
                // far-field quadrature's agreement with analytic Larmor) and,
                // at best, 1.9e-8 (near-field contamination of the production
                // control sphere).  Adding M2 would inject a term four orders
                // below the noise of everything it would be compared against,
                // and cost another set of history stencils to do it.  The
                // same verdict covers E3, and E5 at beta^8 = 8e-18 is not a
                // question anyone needs to ask twice.
                //
                // Nor is it missing from the model's ENERGY bookkeeping: the
                // far-zone Poynting quadrature is not a truncated expansion
                // at all -- directional retardation across the pair keeps E3,
                // M2, toroidal terms and their interference exactly (see
                // electromagneticFieldFluxRates in electrodynamics.hpp).  The
                // truncation is confined to THIS analytic sum, whose job is
                // an instantaneous hazard rate; the flux cannot serve that
                // job, being retarded by the control radius, 1.8e-13 s, which
                // is longer than some runs last.
                const double quantizedPower=
                    stepRadiation.leadingElectricDipolePower
                    +stepRadiation.magneticDipoleFlux.energy
                    +stepRadiation.electricQuadrupolePower;
                stochasticHazard+=quantizedPower/photonEnergy*dt;
                // A while, not an if: a fast step near periapsis can bank
                // more than one photon's worth of hazard at once.
                while(stochasticHazard>=stochasticThreshold) {
                    stochasticHazard-=stochasticThreshold;
                    // Deliberately NOT also credited to s.radiatedEnergy/
                    // orbitalRadiatedEnergy here: those are already fed,
                    // unconditionally and independently of the reaction
                    // model, by integrateElectrodynamicStep's own flux
                    // quadrature (the true Poynting flux of whatever
                    // trajectory results -- see radiatedEnergyIncrement in
                    // electrodynamics.hpp).  That flux does not stop just
                    // because chargeReaction is zero between photons: the
                    // pair is still accelerating under the bare mutual
                    // Lorentz force, and Maxwell's equations do not know
                    // about this mode's bookkeeping.  Adding the kick's
                    // removed energy to the SAME ledger the flux quadrature
                    // already fills would double-count -- the two are
                    // meant to be compared (as with every other model's
                    // "reaction/flux" diagnostics), not merged.  Only the
                    // MECHANICAL removal (the velocity kick itself) is this
                    // function's job.
                    const Vec3 separationVector=
                        s.firstPosition-s.secondPosition;
                    const Vec3 relativeMomentum=
                        momentum(s.firstVelocity,firstMass)
                        * (secondMass/(firstMass+secondMass))
                        -momentum(s.secondVelocity,secondMass)
                        * (firstMass/(firstMass+secondMass));
                    const Vec3 orbitalNormal=cross(
                        separationVector,relativeMomentum);
                    const Vec3 photonDirection=
                        sampleRotatingDipolePhotonDirection(
                            orbitalNormal,stochasticPhotonStream);
                    const StochasticPhotonRecoil recoil=
                        applyStochasticDipolePhoton(
                            s,photonEnergy,photonDirection);
                    if(std::getenv("CREM_DEBUG"))
                        std::cerr<<"  PHOTON t="<<s.time*1e12<<"ps r="
                                 <<separation(s)*1e15<<"fm hbar*omega="
                                 <<photonEnergy<<"J emitted="<<recoil.emitted
                                 <<" n=("<<photonDirection.x<<','
                                 <<photonDirection.y<<','
                                 <<photonDirection.z<<')'<<std::endl;
                    stochasticThreshold=
                        drawExponentialUnit(stochasticPhotonStream);
                }
            }
        }
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
        :(directCollision?collisionBoundaryRadius:nuclearCutoff);
    MechanicalTrajectoryResult run = runMechanicalTrajectory(
        s, observationTime, trajectoryCutoff, options);
    const double timeToCutoff = run.outcome == SimulationOutcome::ReachedCutoff
                              ? run.elapsedTime : std::numeric_limits<double>::infinity();
    return {std::move(run.frames), {relativeEnergy, orbitalAngularMomentum,
            predictedClosestApproach, dipoleAlignment, timeToCutoff, phenomenon, seed},
            run.outcome, run.minimumSeparation, run.elapsedTime,
            run.finalRadiatedEnergy, run.maximumBeta};
}
