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

// The threshold the quantized emission actually fires on.  Exp(1) gives the
// Poisson process; the constant 1 gives the deterministic variant, which
// fires exactly when the continuously accumulated loss reaches one quantum
// (see gDeterministicEmission for why that is the same thing as emitting on
// a level crossing).  The stream is advanced either way, so switching modes
// does not shift every later draw in the run.
double drawEmissionThreshold(std::uint64_t& streamState) {
    const double exponential=drawExponentialUnit(streamState);
    return gDeterministicEmission ? 1.0 : exponential;
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
    SimulationStopReason stopReason=SimulationStopReason::NumericalFailure;
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
    SimulationStopReason stopReason=SimulationStopReason::NumericalFailure;
    ClassicalTrajectoryEngine trajectory(s,
        {.relativeTolerance=1.0e-5,.maximumDepth=12,
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
    double stochasticThreshold=drawEmissionThreshold(stochasticPhotonStream);
    // Open emission window.  A photon is not handed over in one instant: the
    // event opens a window one orbital period long over which the quantum is
    // paid out continuously.  See the firing site below for why an
    // instantaneous kick cannot work at all here.
    double emissionRemaining=0.0;      // J still owed by the active photon
    double emissionRate=0.0;           // W, the quantum spread over the window
    Vec3 emissionDirection;

    bool reachedObservationCeiling=false;
    bool externalStopRequested=false;
    while (s.time < observationTime && separation(s) > trajectoryCutoff) {
        if(options.stopRequested&&options.stopRequested()) {
            externalStopRequested=true;
            break;
        }
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
            stopReason=SimulationStopReason::ReachedCutoff;
            s = eventState;
            break;
        }

        minimumSeparation = std::min(minimumSeparation, currentSeparation);
        maximumBeta = std::max(maximumBeta,
            std::max(s.firstVelocity.norm(), s.secondVelocity.norm()) / c);
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;

        // Pay out an open emission window, before this step's hazard is
        // banked, so a photon in flight is settled on the state it was
        // actually emitted from.  The increment is the quantum's share of
        // this step, and each one is small enough that the kinematic ceiling
        // -- the reason the window exists at all -- never binds: it is
        // bounded by photonEnergy*dt/T against a kinetic energy of order
        // photonEnergy/2, and dt is a small fraction of the orbital period T.
        //
        // A partial payout that runs out of trajectory (the run ends, or the
        // pair reaches the boundary, mid-window) is left unfinished rather
        // than dumped in one lump: that would be the instantaneous kick this
        // whole construction replaces.
        if(emissionRemaining>0.0&&emissionRate>0.0) {
            // Evenly over the window, clamped so the last step cannot overpay
            // -- and clamped again against what the pair can actually afford
            // this step.  The second clamp is what makes the ceiling closed
            // BY CONSTRUCTION rather than merely usually: the nominal rate is
            // hbar*omega/(2 pi/omega), which grows as omega^2 as the orbit
            // tightens, while the kinetic energy it must come out of grows
            // only as the binding energy.  Deep in the collapse the nominal
            // rate therefore outruns the reservoir -- measured before this
            // clamp existed, refused increments ran 2.5x the available
            // kinetic energy on average and up to 13.7x.  Capping the share
            // at a tenth of that energy simply LENGTHENS the window there;
            // the total paid is still exactly one quantum, and no increment
            // can ever be refused.
            // The reservoir is the CM kinetic energy, which is what the
            // guard inside applyStochasticDipolePhoton actually tests -- not
            // the laboratory one.  The two differ once the pair's COM has
            // picked up recoil, and using the laboratory value here left the
            // clamp too generous: measured, it still allowed 1.3% of
            // increments to be refused.  W is the invariant
            // sqrt((sum E)^2 - |sum p c|^2) and W - (m1+m2)c^2 is the energy
            // a photon can actually be paid out of.
            const auto emissionFirst=two_body::fourMomentumFromVelocity(
                s.firstVelocity,firstMass);
            const auto emissionSecond=two_body::fourMomentumFromVelocity(
                s.secondVelocity,secondMass);
            const double emissionEnergy=
                emissionFirst.energy+emissionSecond.energy;
            const Vec3 emissionMomentum=
                emissionFirst.momentum+emissionSecond.momentum;
            const double emissionInvariantSquared=
                emissionEnergy*emissionEnergy
                -emissionMomentum.squaredNorm()*c*c;
            const double emissionInvariant=emissionInvariantSquared>0.0
                ?std::sqrt(emissionInvariantSquared):0.0;
            const double affordable=0.1*(emissionInvariant
                -(firstMass+secondMass)*c*c);
            const double share=std::min({emissionRemaining,emissionRate*dt,
                                         std::max(affordable,0.0)});
            if(share>0.0) {
                const StochasticPhotonRecoil increment=
                    applyStochasticDipolePhoton(s,share,emissionDirection);
                // A refused increment means even this differential does not
                // fit, which the ceiling analysis says should not happen; if
                // it ever does, close the window rather than retry forever.
                emissionRemaining=increment.emitted
                    ?emissionRemaining-share:0.0;
            }
        }

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
        // The secular estimator in crem_collapse.hpp carries M1 too, and by
        // the same channel split as here: it draws each skipped photon's
        // channel from the two hazards' ratio, gives an M1 photon the
        // coherent moment's precession axis instead of the orbital normal,
        // and withholds the Kepler harmonic series from it.  This comment
        // used to say that path was E1-only and "cannot be otherwise",
        // because the osculating elements carry no spin state -- true when
        // written, false since the coupled secular spin-orbit solver began
        // carrying firstDipole/secondDipole through that same checkpoint
        // loop.  It also claimed the bound phenomena measure M1 identically
        // zero through m1+m2 cancelling; that holds only for ortho (S=1,
        // parallel spins, anti-aligned moments).  Para (S=0) has ALIGNED
        // moments -- opposite charges invert the spin-moment relation -- so
        // its M1 amplitude adds rather than cancels.  See
        // coherentMagneticDipoleOrbitAveragedEmission's own comment.
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
                // Ground-state emission floor, the in-orbit half of the
                // same experiment gated in crem_collapse.hpp.  Practically
                // every photon fires on the secular path (the mechanical
                // hazard is ~5.5e-5 per revolution and a --mode 2 run
                // instrumented here recorded zero calls), so this branch is
                // for completeness and for --level runs rather than for
                // production weight.
                const double floorSpecificEnergy=
                    -pairCoulombStrength
                    /(2.0*pairBohrRadius(activePair)*pairReducedMass);
                const bool belowFloor=gGroundStateEmissionFloor
                    &&conservativeParticleEnergy(s)/pairReducedMass
                        <=floorSpecificEnergy;
                stochasticHazard+=belowFloor
                    ?0.0:quantizedPower/photonEnergy*dt;
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
                    // Which channel THIS photon belongs to, drawn with
                    // probability equal to that channel's share of the
                    // combined power quantizedPower sums (E1 charge, M1
                    // magnetic dipole, E2 charge quadrupole -- see the
                    // comment above quantizedPower).  The unified hazard
                    // still fires at the orbital rate for every channel
                    // (that choice is deliberately unchanged: mean radiated
                    // energy is P for any quantum size, see above), but a
                    // photon attributed to M1 should not be drawn from the
                    // orbiting-charge dipole's own angular law -- it comes
                    // from the coherent moment m=m1+m2 precessing, whose
                    // rotation axis is cross(m,mdot), not the orbital
                    // normal.  Falls back to the orbital normal if the
                    // moment is not actually precessing (m and mdot
                    // parallel, or either vanishing): direction is then
                    // irrelevant to first order anyway, since a
                    // non-precessing coherent moment radiates negligibly
                    // (this photon's very existence already means
                    // quantizedPower includes a real M1 share, but a
                    // vanishing cross product only zeroes the AXIS
                    // estimate, not the power upstream).
                    Vec3 photonAxis=orbitalNormal;
                    if(quantizedPower>0.0
                       &&drawUniformUnit(stochasticPhotonStream)
                           <stepRadiation.magneticDipoleFlux.energy
                               /quantizedPower) {
                        const RetardedDipoleKinematics firstMoment=
                            historicalDipoleKinematics(
                                trajectory.history(),s,true,s.time);
                        const RetardedDipoleKinematics secondMoment=
                            historicalDipoleKinematics(
                                trajectory.history(),s,false,s.time);
                        const Vec3 coherentMoment=
                            firstMoment.moment+secondMoment.moment;
                        const Vec3 coherentMomentRate=
                            firstMoment.firstDerivative
                            +secondMoment.firstDerivative;
                        const Vec3 precessionAxis=
                            cross(coherentMoment,coherentMomentRate);
                        if(precessionAxis.norm()>1.0e-300)
                            photonAxis=precessionAxis;
                    }
                    const Vec3 photonDirection=
                        sampleRotatingDipolePhotonDirection(
                            photonAxis,stochasticPhotonStream);
                    // KINEMATIC CEILING, and why nothing here tries to
                    // raise it.  applyStochasticDipolePhoton moves only the
                    // velocities: the positions, and therefore the potential
                    // energy, are untouched at the emission instant.  So the
                    // photon can be paid for out of the CM KINETIC energy and
                    // nothing else, and its own guard refuses anything larger.
                    //
                    // That ceiling bites hard, because hbar*omega_orb is not
                    // small against it.  The virial theorem puts the kinetic
                    // energy of a circular orbit at exactly the binding
                    // energy, while hbar*omega/E_bind = 2 sqrt(a_Ps/a) -- two
                    // at the pair Bohr radius and worse as the orbit tightens.
                    // Measured as a fraction of the orbital PERIOD on which a
                    // photon of n hbar*omega fits at all:
                    //
                    //     e      n=1      n=5      n=16
                    //     0      0        0        0
                    //     0.3    0        0        0
                    //     0.5    0.149    0        0
                    //     0.9    0.113    0.017    0
                    //     0.97   0.098    0.015    0.003
                    //
                    // A circular or mildly eccentric orbit cannot emit the
                    // fundamental ANYWHERE on its period.
                    //
                    // Two "obvious" repairs are both wrong, which is the
                    // reason for this comment rather than a patch.  Keeping
                    // the banked hazard and retrying on a later step
                    // deadlocks outright for e <~ 0.3: the hazard grows
                    // without bound and no photon ever fires.  Relaxing the
                    // guard to the secular path's own condition (W^2-2 W E >
                    // 0, i.e. E < W/2 ~ 511 keV) is not a unification either
                    // -- crem_collapse.hpp can afford that because it removes
                    // the energy from the osculating ELEMENTS, an
                    // orbit-averaged quantity whose reservoir is the binding
                    // energy, whereas an instantaneous kick at fixed position
                    // has only the kinetic energy to draw on.  The two guards
                    // differ because the two emissions differ in kind, not
                    // because one of them is built on the wrong invariant.
                    //
                    // What is left is a real domain limit: the correspondence
                    // prescription hbar*omega_orb and the instantaneous-kick
                    // emission model are mutually inconsistent on a
                    // near-circular orbit, and reconciling them needs a
                    // different emission model -- not a different threshold,
                    // and not a rescaled quantum either.
                    //
                    // The second half of that is worth the arithmetic, since
                    // rescaling is the obvious thing to try.  At n=1 the
                    // virial puts the kinetic energy at the binding energy R
                    // while hbar*omega is 2R, so the ceiling admits a quantum
                    // only below hbar*omega/2.  But the classical hazard is
                    // (1/3) alpha^5 m_e c^2/hbar against Gamma_para's (1/2),
                    // so the quantum that reproduces the measured
                    // positronium lifetime EXACTLY is (2/3) hbar*omega.
                    // Those two demands are incompatible by a factor of 4/3:
                    //
                    //   quantum        lifetime      vs 124.5 ps   ceiling
                    //   hbar*omega     186.74 ps     x1.500        refused
                    //   (2/3)hbar*w    124.49 ps     x1.000        refused
                    //   (1/2)hbar*w     93.37 ps     x0.750        marginal
                    //   (1/3)hbar*w     62.25 ps     x0.500        fits
                    //   m_e c^2          7.01 us     x5.6e4        fits
                    //
                    // The first quantum that actually clears the ceiling
                    // costs a factor of two on the lifetime, worse than the
                    // 1.5 it starts from.  So no rescaling in EITHER
                    // direction buys the ceiling without losing the
                    // comparison, which is a stronger statement than the
                    // earlier "do not raise the quantum" and closes the
                    // direction the other way too.
                    // Production does not meet this: essentially every photon
                    // fires through the secular path (this one's hazard is
                    // ~5.5e-5 per orbit, and a --mode 2 run instrumented at
                    // this call site recorded zero calls).  Refusals are
                    // reported below through recoil.emitted.
                    // OPEN THE WINDOW instead of kicking.  The kinematic
                    // ceiling documented above makes the instantaneous kick
                    // unusable here -- at n <= 1, which is the model's whole
                    // domain, hbar*omega/E_kinetic = 2/n >= 2 and the guard
                    // refuses every photon, silently.  Two repairs were
                    // measured and both fail: keeping the banked hazard and
                    // retrying deadlocks for e <~ 0.3 (the photon fits
                    // nowhere on such an orbit), and no rescaling of the
                    // quantum works in either direction, the ceiling wanting
                    // E < hbar*omega/2 while the measured positronium
                    // lifetime wants (2/3) hbar*omega.
                    //
                    // What does work is changing the emission's KIND rather
                    // than its size: pay the quantum out over a window one
                    // orbital period long.  At no instant is more than a
                    // differential removed, and across the window the
                    // POSITIONS evolve, so the energy comes from the orbit --
                    // from the potential well -- instead of from the frozen
                    // instant's kinetic energy.  That is the same change of
                    // kind the secular path already relies on, carried out
                    // on a resolved trajectory.  Measured: at n=1 the orbit
                    // goes 105.8 -> 35.3 pm, still 184x above the Compton
                    // barrier, so there is always somewhere to fall to.
                    //
                    // The increments are handed to the same
                    // applyStochasticDipolePhoton as before, which is what
                    // keeps the relativistic bookkeeping exact and, as a
                    // bonus, carries the right angular momentum: scaling the
                    // relative momentum's magnitude gives dE/dL = v^2/(v_t r),
                    // equal to omega for a circular orbit -- the ratio a
                    // rotating E1 dipole actually radiates -- and within
                    // 0.25% of it at the measured median emission
                    // eccentricity of 0.05.
                    //
                    // The cost, stated plainly: the trajectory is no longer
                    // exactly conserved between photons, since a continuous
                    // force acts inside the window.  The duty cycle bounds
                    // that at 1.5e-6 -- one photon per 669000 orbits, one
                    // orbit long -- so 99.99985% of the trajectory keeps the
                    // exact conservation this mode is built on.
                    // ACCUMULATE, never overwrite.  At the real hazard --
                    // one photon per 669000 orbits against a window one
                    // orbit long -- two windows overlapping is a 1.5e-6
                    // event and the policy hardly matters; overwriting
                    // would nevertheless discard an unpaid remainder, which
                    // is an energy leak rather than a modelling choice.
                    // BLEND, don't overwrite, the direction the still-unpaid
                    // remainder is paid out along -- overwriting silently
                    // handed the old window's whole unpaid remainder to the
                    // new photon's direction, which conserves the energy (it
                    // is accumulated, not replaced) but not the momentum of
                    // that already-committed-but-unpaid share.  Weighted by
                    // energy, since that is what each direction is worth in
                    // the payout below.
                    const Vec3 blendedDirection=
                        emissionDirection*emissionRemaining
                            +photonDirection*photonEnergy;
                    const double blendedNorm=blendedDirection.norm();
                    emissionDirection=blendedNorm>0.0
                        ?blendedDirection*(1.0/blendedNorm):photonDirection;
                    emissionRemaining+=photonEnergy;
                    // One orbital period: an E1 photon of frequency omega
                    // cannot be assembled from a shorter wave train.
                    emissionRate=emissionRemaining
                        /(2.0*pi/std::max(omega,1.0e-300));
                    const StochasticPhotonRecoil recoil{true,photonEnergy,
                                                        photonDirection};
                    if(std::getenv("CREM_DEBUG"))
                        std::cerr<<"  PHOTON t="<<s.time*1e12<<"ps r="
                                 <<separation(s)*1e15<<"fm hbar*omega="
                                 <<photonEnergy<<"J emitted="<<recoil.emitted
                                 <<" n=("<<photonDirection.x<<','
                                 <<photonDirection.y<<','
                                 <<photonDirection.z<<')'<<std::endl;
                    stochasticThreshold=
                        drawEmissionThreshold(stochasticPhotonStream);
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
            || externalStopRequested)) {
        outcome = SimulationOutcome::ObservationLimit;
        stopReason=externalStopRequested
            ?SimulationStopReason::StopRequested
            :SimulationStopReason::ObservationTimeLimit;
        elapsedTime = s.time;
        finalRadiatedEnergy = s.radiatedEnergy;
    }
    if (options.collectFrames && outcome == SimulationOutcome::ObservationLimit
        && (frames.empty() || frames.back().time < s.time)) {
        frames.push_back(makeFrame(s));
        if (options.frameReady) options.frameReady(frames.back());
    }
    return {s, outcome, stopReason, minimumSeparation, elapsedTime,
            finalRadiatedEnergy, maximumBeta, std::move(frames)};
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
    // a_n = n^2 a_pair.  gInitialPrincipalLevel is 1 unless --level says
    // otherwise, and at 1 this is bit-identical to the historical expression.
    const double initialSeparation = pairBohrRadius(activePair)
        * static_cast<double>(gInitialPrincipalLevel)
        * static_cast<double>(gInitialPrincipalLevel);
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
        //
        // WHAT THIS BAND ACTUALLY PREPARES, worked out rather than assumed.
        // With v_circ^2 = K/r at r = a_Ps, and f_r the radial speed in the
        // same units,
        //
        //     a_0/a_Ps = 1/(2 - f_r^2 - f^2),   n_E = sqrt(a_0/a_Ps),
        //     e_0^2    = 1 + f^2 (f_r^2 + f^2 - 2),
        //
        // so the sampled state carries TWO different Bohr levels: an
        // angular-momentum one, L/hbar = f, and an energy one, n_E.  Setting
        // them equal gives (f^2-1)^2 = 0, i.e. they agree at f = 1 and NOWHERE
        // ELSE in the band.  Everywhere else the state has e_0 up to 0.28 and
        // is not a Bohr state at all -- which matters only because other parts
        // of the model then apply ladder concepts to it: quantumFor reads the
        // ENERGY level, clampAboveGroundStateAngularMomentum reads the
        // ANGULAR-MOMENTUM one.
        //
        // Consequences, measured:
        //   * n_E runs 0.903-1.166 and P(n_E < 1) = 0.493, so about HALF of
        //     all trajectories are prepared below the ground state.  That is
        //     what makes --ground-state-floor unusable at level 1: the floor
        //     declares them settled before they move (see
        //     CremCollapseEstimate::preparedBelowGroundState).
        //   * Half likewise start with L = f*hbar < hbar, below the same
        //     flag's angular-momentum floor, so the first photon would clamp
        //     L UPWARD -- an emission that increases angular momentum.
        //   * The band is centred on f = 1 but NOT on n = 1: 1/(2 - f^2) is
        //     convex, so E[a_0/a_Ps] = 1.0288 and E[n_E] = 1.0117.  The
        //     ensemble sits 2.9% outside the ground-state radius on average.
        //   * It is not a small effect on the headline observable.  Holding f
        //     fixed, the mean collapse time is 108.6 ps at f=0.88, 188.0 at
        //     f=1.00 and 367.2 at f=1.12 -- a factor 3.4 across the band
        //     (classical t ~ a^3 predicts 4.6).  Decomposing the variance,
        //     sigma/mean is 1.127 with the band against ~0.93 at fixed f, so
        //     the band supplies about 32% of the collapse time's variance and
        //     the photon process the other 68%.
        //
        // The WIDTH is inherited, not derived: it was kept "comparable" to the
        // old band's when the centring was corrected.  Deriving it instead --
        // as the width of an SED fluctuation-dissipation equilibrium -- is the
        // --zpf thread, not a constant that can be written down here.  The
        // sharp alternative (f = 1, f_r = 0 exactly, so n_E = L/hbar = 1 and
        // e_0 = 0) is a one-line change that would make the two level
        // definitions agree by construction; it is NOT taken here because it
        // moves the production mean by -12% and that is a decision about what
        // the model claims, not a bug fix.
        //
        // SHARP PREPARATION, taken.  f = 1 and f_r = 0 exactly, so
        // n_E = L/hbar = 1 and e_0 = 0: the two Bohr-level definitions agree
        // BY CONSTRUCTION rather than at one point of a band, and --level n
        // becomes an exact Bohr state (at r = n^2 a_Ps a circular orbit has
        // L = n hbar and E = -R/n^2 identically).  What is given up is a
        // spread whose WIDTH was never derived -- it was inherited from the
        // old [0.72, 0.97] band when the centring was corrected -- and which
        // supplied about 32% of the collapse time's variance while moving its
        // mean by a factor of 3.4 depending on where in the band a trajectory
        // sat.  The remaining spread is the emission process, which is what
        // the collapse time's spread ought to mean.
        //
        // The two random draws are still CONSUMED, not removed: dropping them
        // would reseed every subsequent draw in the stream and make this a
        // wholesale change of trajectory rather than a change of initial
        // condition, destroying comparability with everything measured
        // before.
        const double bandRadial = -0.10 + 0.20 * unitRandom(random);
        const double bandTangential = 0.88 + 0.24 * unitRandom(random);
        // CREM_INITIAL_BAND=1 restores the old sampled band, kept one flag
        // away for regression comparison -- the same opt-out shape the
        // harmonic correction uses, and for the same reason.
        const char* bandEnv = std::getenv("CREM_INITIAL_BAND");
        const bool sampledBand = bandEnv && std::strcmp(bandEnv, "0") != 0;
        radialSpeed = sampledBand ? circularSpeed * bandRadial : 0.0;
        tangentialSpeed = circularSpeed
            * (sampledBand ? bandTangential : 1.0);
    }
    const Vec3 relativeVelocity{radialSpeed, tangentialSpeed, 0.0};
    s.firstVelocity = relativeVelocity * (secondMass / (firstMass + secondMass));
    s.secondVelocity = relativeVelocity * (-firstMass / (firstMass + secondMass));
    // CREM_COM_DRIFT: give the WHOLE pair a common velocity, in units of c,
    // normal to the orbital plane.  Off by default; any nonzero value is a
    // numerical experiment, the same shape as --zpf-scale's "anything other
    // than 1 is not the physical field".
    //
    // The bound initial conditions are otherwise prepared at exactly zero
    // total momentum, so the two velocities are equal and opposite.  A common
    // drift is what makes both share a velocity SENSE, and it is not a
    // cosmetic relabelling of the frame for the dipole sector: the motional
    // electric dipole is p_i = gamma (v_i x mu_i)/c^2, so with v_1 = -v_2 the
    // drift-free pair has para's two p cancel and ortho's add, while a common
    // drift contributes the same v_drift x mu_i to both and therefore reverses
    // which channel adds.  That is the asymmetry this switch exists to expose;
    // it is deliberately applied normal to the plane so the internal orbit is
    // untouched and the drift enters as a boost rather than as a perturbed
    // orbit.
    //
    // SCOPE, measured, and the reason to read this before trusting a result
    // from it: the drift reaches the mechanical seed run and the carried
    // dipoles, but NOT the secular collapse estimator, which is frame-locked
    // to the pair's zero-momentum frame twice over -- estimateCremCollapse
    // holds centreOfMassVelocity at zero by construction ("starts at the true
    // zero and only ever moves because a photon kicked it"), and
    // orbitAveragedBmtAngularVelocities rebuilds every quadrature node's
    // velocities from the orbital elements as +-relativeVelocity*(m/M), with
    // no drift term.  So a common drift is projected straight back out of the
    // collapse integral.  Measured at 0.30 c: the Kaplan-Meier collapse time
    // does not move at all (147.818 ps para and 147.819 -> 147.820 ps ortho,
    // against a 12-run sigma/mean of 1e-5), while the pair's total motional
    // dipole changes by 86x.  A null collapse-time result from this switch is
    // therefore a statement about the estimator's scope, NOT about physics,
    // and must not be reported as one.  Making it a real measurement means
    // seeding centreOfMassVelocity from the initial state and carrying the
    // drift into the quadrature -- both of which touch that documented
    // zero-momentum assumption deliberately.
    if(const char* driftText=std::getenv("CREM_COM_DRIFT")) {
        const double drift=std::atof(driftText)*c;
        if(std::isfinite(drift)&&drift!=0.0) {
            const Vec3 driftVelocity{0.0,0.0,drift};
            s.firstVelocity=s.firstVelocity+driftVelocity;
            s.secondVelocity=s.secondVelocity+driftVelocity;
        }
    }

    s.firstDipole = randomDirection() * firstMagneticMoment;
    // SPIN QUANTIZATION, third floor of --ground-state-floor and the one that
    // gives the channels a dynamical difference instead of a label.
    //
    // Sampled normally, the mutual angle is drawn from a RANGE -- para gets
    // cos >= 0.5, ortho everything below -- and that range is why the one
    // exact mechanism the model owns does not survive.  The M1 power is
    // computed coherently from the total moment m = mu1 + mu2, so at cos = -1
    // it cancels EXACTLY (measured: 0.0000e+00, not merely small) while at
    // cos = +1 it interferes constructively.  That is the same shape as the
    // QED selection rule, where C|n gamma> = (-1)^n against
    // C|Ps> = (-1)^(L+S) closes the leading channel for one spin state and
    // pushes it an order of alpha higher.  But averaged over the sampled
    // range the M1 share came out 1.875e-3 for para against 1.903e-3 for
    // ortho -- a difference consistent with zero.  The mechanism exists at
    // the endpoint and is washed out by the sampling.
    //
    // Quantizing S removes the range: S=0 and S=1 are exact states, so the
    // moments are exactly aligned or exactly anti-aligned, not somewhere in
    // a band.  Opposite charges invert the spin-moment relation, which is why
    // ANTI-parallel spins (para, S=0) give ALIGNED moments -- see the README's
    // Sonda 4, where that reversal is derived and measured.
    //
    // Like the other two floors this is IMPORTED.  It does not derive the
    // photon count: 2 against 3 is a statement about the final state's
    // C-parity, and this model has no annihilation dynamics at all -- no
    // contact channel and no rate, established separately.  What it buys is
    // that the channel difference becomes a property of the configuration
    // rather than of the --phenomenon switch.
    if(gGroundStateEmissionFloor && (sampledScenario==2||sampledScenario==3)) {
        const double sign = sampledScenario==2 ? 1.0 : -1.0;
        s.secondDipole = s.firstDipole
            * (sign*secondMagneticMoment/firstMagneticMoment);
    } else {
        do {
            s.secondDipole = randomDirection() * secondMagneticMoment;
        } while ((sampledScenario == 2 && dot(s.firstDipole, s.secondDipole)
                                      / (firstMagneticMoment*secondMagneticMoment) < 0.5)
              || (sampledScenario == 3 && dot(s.firstDipole, s.secondDipole)
                                      / (firstMagneticMoment*secondMagneticMoment) >= 0.5));
    }

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
            run.outcome, run.stopReason, run.minimumSeparation, run.elapsedTime,
            run.finalRadiatedEnergy, run.maximumBeta};
}
