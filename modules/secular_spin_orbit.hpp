#pragma once

// Coupled orbit-averaged classical spin-orbit transport.
//
// This module is included after electrodynamics.hpp, inside the program's
// anonymous namespace.  It deliberately advances the three angular-momentum
// reservoirs together:
//
//     J = L_orbit + mu_1/gamma_1 + mu_2/gamma_2.
//
// The magnetic moments follow the same orbit-averaged Thomas-BMT angular
// velocities as the resolved mechanical integrator.  The orbital vector is
// then the conservative reaction required by the same reduced model, rather
// than an independently averaged scalar torque.  Torque from a configured
// external/ZPF field is tracked separately and changes the pair's J instead
// of being spuriously fed into L.  Radiation and photon kicks remain separate
// dissipative operators in crem_collapse.hpp.

struct OrbitAveragedBmtAngularVelocities {
    Vec3 first;
    Vec3 second;
    // Part driven by configured external/ZPF fields.  It changes the
    // particles' J because the source of that field is outside the reduced
    // two-body system; only the remainder receives an orbital reaction.
    Vec3 firstExternal;
    Vec3 secondExternal;
    // < |d2(mu1+mu2)/dt2|^2 >, the orbit average of the SQUARED coherent
    // second derivative of the pair's total moment -- the quantity the M1
    // Larmor power is actually proportional to.  It is accumulated here, at
    // each phase node, rather than reconstructed downstream from the averaged
    // rates above, because those are two different numbers: < |mu''|^2 > is
    // not |mu''(<omega>)|^2.  See the accumulation below for why the gap is
    // not a refinement but the leading term.
    double coherentSecondDerivativeSquared=0.0;
    // < osculatingOrbitalFrequency >, the orbit average of the SAME
    // instantaneous frequency the ZPF band is tuned to at every node, and the
    // rate at which the pair's integrated orbital phase actually grows.  It is
    // NOT the mean motion: osculatingOrbitalFrequency goes as r^-3/2, so it
    // swings by ((1+e)/(1-e))^3/2 around an eccentric orbit -- 595x at
    // e^2=0.945 -- and its time average comes out 1.5978 n there, not n.  The
    // two coincide only for a circle.  advanceCoupledSecularSpinOrbit advances
    // zeroPointPhase by this, so the secular path accumulates the same
    // integrated phase the mechanical path does one step at a time.
    double averagedOrbitalFrequency=0.0;
    int phaseNodes=0;
    bool valid=false;
};

struct SecularSpinOrbitState {
    Vec3 orbitalAngularMomentum;
    Vec3 firstDipole;
    Vec3 secondDipole;
    // Mirrors State::zeroPointPhase for the orbit-averaged treatment: the
    // ZPF sampled inside orbitAveragedBmtAngularVelocities needs the
    // pair's REAL accumulated phase (advanced by the mean motion below,
    // the same quantity osculatingOrbitalFrequency*dt advances it by on
    // the mechanical path), not a phase frozen at 0 for every checkpoint.
    double zeroPointPhase=0.0;
    // Unit vector towards periapsis.  L fixes only the orbital plane; for an
    // eccentric orbit the apsidal line is additional physical data because
    // the strongest r^-3 fields are sampled preferentially near periapsis.
    // A zero vector is accepted as a backwards-compatible "unspecified"
    // value and replaced by a deterministic in-plane direction.
    Vec3 periapsisDirection;
};

struct SecularSpinOrbitAdvance {
    SecularSpinOrbitState state;
    bool completed=false;
    int substeps=0;
    double maximumSubstepAngle=0.0;
    double relativeAngularMomentumResidual=
        std::numeric_limits<double>::infinity();
    Vec3 externalAngularMomentumTransfer;
};

Vec3 rotateDipoleByAngularVelocity(const Vec3& dipole,
                                   const Vec3& angularVelocity,
                                   double elapsedTime) {
    if(elapsedTime==0.0) return dipole;
    const double angularSpeed=angularVelocity.norm();
    if(!(angularSpeed>0.0)) return dipole;
    const Vec3 axis=angularVelocity/angularSpeed;
    const double angle=angularSpeed*elapsedTime;
    const double cosine=std::cos(angle);
    const double sine=std::sin(angle);
    return dipole*cosine+cross(axis,dipole)*sine
        +axis*(dot(axis,dipole)*(1.0-cosine));
}

// Project a preferred apsidal direction into the plane normal to L.  Keeping
// this operation in one place prevents the secular quadrature and the
// mechanically resolved checkpoint orbit from silently choosing different
// planes or different periapsides.
Vec3 orbitPlaneDirection(const Vec3& orbitalAngularMomentum,
                         const Vec3& preferredDirection) {
    const double orbitalNorm=orbitalAngularMomentum.norm();
    if(!(orbitalNorm>0.0)||!std::isfinite(orbitalNorm)) return {};
    const Vec3 normal=orbitalAngularMomentum/orbitalNorm;
    if(isFinite(preferredDirection)) {
        const Vec3 projected=preferredDirection
            -normal*dot(preferredDirection,normal);
        const double projectedNorm=projected.norm();
        if(projectedNorm>1.0e-14*std::max(1.0,preferredDirection.norm()))
            return projected/projectedNorm;
    }
    const Vec3 seed=std::abs(normal.x)<0.9
        ?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    const Vec3 fallback=cross(seed,normal);
    const double fallbackNorm=fallback.norm();
    return fallbackNorm>0.0?fallback/fallbackNorm:Vec3{};
}

// Minimal-rotation transport of the apsidal line when spin exchange or a
// photon tilts L.  The reduced model has no independent Runge-Lenz equation;
// parallel transport is the neutral closure that introduces no arbitrary
// rotation about the new normal.  Projection at the end removes round-off.
Vec3 transportOrbitPlaneDirection(const Vec3& direction,
                                   const Vec3& oldAngularMomentum,
                                   const Vec3& newAngularMomentum) {
    const double oldNorm=oldAngularMomentum.norm();
    const double newNorm=newAngularMomentum.norm();
    if(!(oldNorm>0.0)||!(newNorm>0.0))
        return orbitPlaneDirection(newAngularMomentum,direction);
    const Vec3 oldNormal=oldAngularMomentum/oldNorm;
    const Vec3 newNormal=newAngularMomentum/newNorm;
    const Vec3 radial=orbitPlaneDirection(oldAngularMomentum,direction);
    const Vec3 rotationAxisVector=cross(oldNormal,newNormal);
    const double sine=rotationAxisVector.norm();
    const double cosine=std::clamp(dot(oldNormal,newNormal),-1.0,1.0);
    Vec3 transported=radial;
    if(sine>1.0e-15) {
        const Vec3 axis=rotationAxisVector/sine;
        transported=radial*cosine+cross(axis,radial)*sine
            +axis*(dot(axis,radial)*(1.0-cosine));
    }
    return orbitPlaneDirection(newAngularMomentum,transported);
}

OrbitAveragedBmtAngularVelocities orbitAveragedBmtAngularVelocities(
        double semiMajorAxis,const Vec3& orbitalAngularMomentum,
        const Vec3& firstDipole,const Vec3& secondDipole,
        double reducedMass,double accumulatedZeroPointPhase=0.0,
        Vec3 periapsisDirection={},int phaseNodeOverride=0) {
    OrbitAveragedBmtAngularVelocities result;
    if(!(semiMajorAxis>0.0)||!(reducedMass>0.0)
       ||!std::isfinite(semiMajorAxis)||!std::isfinite(reducedMass)
       ||!isFinite(orbitalAngularMomentum)||!isFinite(firstDipole)
       ||!isFinite(secondDipole)||!isFinite(periapsisDirection)) return result;
    const double orbitalNorm=orbitalAngularMomentum.norm();
    if(!(orbitalNorm>0.0)) return result;
    // The secular state carries a and the full orbital L, hence it also
    // carries the eccentricity; treating |L| as orientation-only silently
    // replaced every osculating ellipse by a circle.  For the relative
    // Coulomb problem
    //
    //     L^2 = mu K a (1-e^2),       K = k |q1 q2|.
    //
    // A value above the circular limit is not an ellipse with this energy.
    // The wider collapse estimator already defines that over-circular seam as
    // e=0: spin exchange can change L while this conservative sub-operator
    // holds a fixed, and adding a new spin-energy equation here would be a
    // different physical model.  Preserve that explicit projection only for
    // L>=L_circular; throughout the physical elliptic domain |L| now determines
    // the actual eccentricity instead of being discarded.
    const double circularAngularMomentumSquared=
        reducedMass*pairCoulombStrength*semiMajorAxis;
    if(!(circularAngularMomentumSquared>0.0)
       ||!std::isfinite(circularAngularMomentumSquared)) return result;
    const double oneMinusEccentricitySquared=
        orbitalNorm*orbitalNorm/circularAngularMomentumSquared;
    if(!(oneMinusEccentricitySquared>0.0)
       ||!std::isfinite(oneMinusEccentricitySquared)) return result;
    const double boundedOneMinusEccentricitySquared=
        std::min(1.0,oneMinusEccentricitySquared);
    const double eccentricity=std::sqrt(std::max(
        0.0,1.0-boundedOneMinusEccentricitySquared));
    const double eccentricityComplement=
        std::sqrt(boundedOneMinusEccentricitySquared);
    const Vec3 normal=orbitalAngularMomentum/orbitalNorm;
    const Vec3 radialHat=orbitPlaneDirection(
        orbitalAngularMomentum,periapsisDirection);
    const double radialNorm=radialHat.norm();
    if(!(radialNorm>0.0)) return result;
    const Vec3 tangentialHat=cross(normal,radialHat);
    const double meanMotion=std::sqrt(pairCoulombStrength
        /(reducedMass*semiMajorAxis*semiMajorAxis*semiMajorAxis));
    if(!std::isfinite(meanMotion)||!(meanMotion>0.0)) return result;

    const double totalMassHere=firstMass+secondMass;
    if(!(totalMassHere>0.0)) return result;
    // Uniform eccentric-anomaly nodes form a periodic trapezoidal quadrature.
    // The factor dM/dE=1-e*cos(E) below converts it to a uniform-in-time
    // average.  Sampling E instead of mean anomaly also resolves periapsis:
    // at e^2=0.945 the old 64 midpoint samples in M skipped the narrow peak,
    // whereas this grid contains E=0 exactly.  The analyticity strip narrows
    // as e approaches one, so refine only eccentric orbits.  Nested powers of
    // two avoid the aliasing seen with unrelated node counts, and a circle
    // keeps the historical 64 evaluations.
    //
    // The schedule's CONSTANT is set by the M1 average below, not by the rate
    // average this function is named for, because the two integrands are not
    // equally peaked.  omega itself rises as (1-e cos E)^-3 towards periapsis,
    // but |mu''|^2 carries omega' as well and rises as (1-e cos E)^-8: at
    // e^2=0.945 that is 2.7e12 between apoapsis and periapsis against 4.6e4
    // for the rate.  Measured refinement ladder for the M1 average at that
    // eccentricity (validation's m1-secular-orbit-average prints it):
    //
    //     64      128     256     512     1024    2048   nodes
    //     -70.8%  -29.5%  -8.6%   -2.2%   -0.57%  -0.14% against the limit
    //     0.71    1.97    3.28    3.80    3.95           convergence order
    //
    // so the old 32/sqrt(1-e) constant, which selects 256 here, sat in the
    // pre-asymptotic regime where the order ratio has not yet reached 4 and
    // the value is still climbing by tens of percent per refinement.  Doubling
    // the constant to 64 selects 512 and lands inside the asymptotic regime;
    // the cap is raised to 2048 for the near-parabolic tail.  Circular orbits
    // are unaffected (target 64, unchanged), so the extra cost falls only on
    // the eccentric orbits that need it.
    const double targetPhaseNodes=
        64.0/std::sqrt(std::max(1.0-eccentricity,1.0e-12));
    int phaseNodes=64;
    while(phaseNodes<2048&&phaseNodes<targetPhaseNodes) phaseNodes*=2;
    // Validation refines this deliberately to measure the quadrature's own
    // convergence -- see the m1-secular-orbit-average check.
    if(phaseNodeOverride>=8) phaseNodes=phaseNodeOverride;
    result.phaseNodes=phaseNodes;
    Vec3 firstOmegaSum,secondOmegaSum;
    Vec3 firstExternalOmegaSum,secondExternalOmegaSum;
    double orbitalFrequencySum=0.0;
    // Instantaneous rates kept per node, not just their running sum: the M1
    // power below needs omega(E) itself and its derivative, neither of which
    // survives averaging.
    std::vector<Vec3> firstOmegaNodes(static_cast<std::size_t>(phaseNodes));
    std::vector<Vec3> secondOmegaNodes(static_cast<std::size_t>(phaseNodes));
    std::vector<double> timeWeightNodes(
        static_cast<std::size_t>(phaseNodes));
    for(int node=0;node<phaseNodes;++node) {
        const double eccentricAnomaly=2.0*pi*node
            /static_cast<double>(phaseNodes);
        const double cosine=std::cos(eccentricAnomaly);
        const double sine=std::sin(eccentricAnomaly);
        const double timeWeight=1.0-eccentricity*cosine;
        if(!(timeWeight>0.0)||!std::isfinite(timeWeight)) return result;
        const Vec3 relativePosition=
            radialHat*(semiMajorAxis*(cosine-eccentricity))
            +tangentialHat*(semiMajorAxis*eccentricityComplement*sine);
        const Vec3 relativeVelocity=
            (radialHat*(-semiMajorAxis*meanMotion*sine)
             +tangentialHat*(semiMajorAxis*meanMotion
                 *eccentricityComplement*cosine))/timeWeight;
        State sample{};
        sample.firstPosition=
            relativePosition*(secondMass/totalMassHere);
        sample.secondPosition=
            relativePosition*(-firstMass/totalMassHere);
        sample.firstVelocity=
            relativeVelocity*(secondMass/totalMassHere);
        sample.secondVelocity=
            relativeVelocity*(-firstMass/totalMassHere);
        if(!(sample.firstVelocity.norm()<c)
           ||!(sample.secondVelocity.norm()<c)) return result;
        // firstDipole/secondDipole (the arguments, and what
        // SecularSpinOrbitState/rotateDipoleByAngularVelocity actually
        // transport) are the PROPER moments: a fixed-magnitude vector
        // precessing under pure rotation is only self-consistent with BMT
        // physics if it is the rest-frame moment BMT actually describes.
        // Copying that value into the lab slot directly (as this used to)
        // skipped the boost the node's own nonzero firstVelocity/
        // secondVelocity calls for -- synchronizeCovariantDipoles derives
        // the lab dipole AND the induced electric dipole
        // (dipole~(v x mu)/c^2) that localRelativisticFields below needs to
        // see the same tensor allExternalForces/retardedExternalForces
        // would reconstruct from this same proper moment and velocity.
        sample.firstProperDipole=firstDipole;
        sample.secondProperDipole=secondDipole;
        sample.zeroPointPhase=accumulatedZeroPointPhase;
        synchronizeCovariantDipoles(sample);
        const StateHistory sampleHistory{State{sample}};
        const LocalElectromagneticFields fields=
            localRelativisticFields(sample,sampleHistory);
        ElectromagneticField externalAtFirst{{},gExternalMagneticField};
        ElectromagneticField externalAtSecond{{},gExternalMagneticField};
        if(gZeroPointField.active()) {
            Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
            const double orbitalFrequency=osculatingOrbitalFrequency(sample);
            gZeroPointField.sample(sample.firstPosition,orbitalFrequency,
                accumulatedZeroPointPhase,firstElectric,firstMagnetic);
            gZeroPointField.sample(sample.secondPosition,orbitalFrequency,
                accumulatedZeroPointPhase,secondElectric,secondMagnetic);
            externalAtFirst.electric+=firstElectric;
            externalAtFirst.magnetic+=firstMagnetic;
            externalAtSecond.electric+=secondElectric;
            externalAtSecond.magnetic+=secondMagnetic;
        }
        // The band frequency this node's ZPF was actually sampled at, both
        // here and inside localRelativisticFields (which reads it off the
        // same sample), time-weighted into the orbit average the phase
        // advance needs.
        orbitalFrequencySum+=
            osculatingOrbitalFrequency(sample)*timeWeight;
        // CREM_DEBUG_FIELDSYM: para is mirror-symmetric in position, velocity
        // and moment, so the two local fields should have equal magnitude and
        // the two precession rates with them.  Printed at one node so the
        // symmetry can be checked rather than argued.
        if(std::getenv("CREM_DEBUG_FIELDSYM")&&node==0) {
            static int fieldSymSamples=0;
            if(fieldSymSamples++<3)
                std::cerr<<"FIELDSYM |B1|="<<fields.atFirst.magnetic.norm()
                    <<" |B2|="<<fields.atSecond.magnetic.norm()
                    <<" ratio="<<fields.atFirst.magnetic.norm()
                        /std::max(fields.atSecond.magnetic.norm(),1.0e-300)
                    <<" |E1|="<<fields.atFirst.electric.norm()
                    <<" |E2|="<<fields.atSecond.electric.norm()
                    <<" |p1|="<<sample.firstElectricDipole.norm()
                    <<" |p2|="<<sample.secondElectricDipole.norm()
                    <<" p1.p2/|p||p|="<<dot(sample.firstElectricDipole,
                        sample.secondElectricDipole)
                        /std::max(sample.firstElectricDipole.norm()
                                 *sample.secondElectricDipole.norm(),1.0e-300)
                    <<'\n';
        }
        const Vec3 firstOmega=thomasBmtEffectiveField(
            sample.firstVelocity,fields.atFirst,firstGFactor)
                *(-firstCharge/firstMass);
        const Vec3 secondOmega=thomasBmtEffectiveField(
            sample.secondVelocity,fields.atSecond,secondGFactor)
                *(-secondCharge/secondMass);
        firstOmegaNodes[static_cast<std::size_t>(node)]=firstOmega;
        secondOmegaNodes[static_cast<std::size_t>(node)]=secondOmega;
        timeWeightNodes[static_cast<std::size_t>(node)]=timeWeight;
        firstOmegaSum+=firstOmega*timeWeight;
        secondOmegaSum+=secondOmega*timeWeight;
        firstExternalOmegaSum+=thomasBmtEffectiveField(
            sample.firstVelocity,externalAtFirst,firstGFactor)
                *(-firstCharge/firstMass*timeWeight);
        secondExternalOmegaSum+=thomasBmtEffectiveField(
            sample.secondVelocity,externalAtSecond,secondGFactor)
                *(-secondCharge/secondMass*timeWeight);
    }
    // --- Orbit-averaged coherent M1 second derivative, < |mu''|^2 > ---
    //
    // For a moment of fixed magnitude carried by a precession omega(t),
    // mu' = omega x mu exactly, hence
    //
    //     mu'' = omega' x mu + omega x (omega x mu).
    //
    // Both terms matter, and the FIRST is normally the larger one.  omega
    // varies over the orbit on the orbital timescale, so |omega'| ~ n |omega|
    // with n the mean motion, making the two terms' ratio |omega'x mu| /
    // |omega x (omega x mu)| ~ n/|omega|.  For a spin-orbit precession -- a
    // fine-structure-scale rate -- n/|omega| is large, so dropping omega' does
    // not lose a correction, it loses the leading term.
    //
    // The average that the M1 Larmor power needs is < |mu''|^2 >, taken over
    // the orbit, and NOT |mu''|^2 rebuilt from the orbit-averaged rates above.
    // The two differ for the usual Jensen reason and the gap is not small
    // here: the power is quartic in omega, while omega itself swings by
    // (1+e)^3/(1-e)^3 between apoapsis and periapsis for an r^-3 field, so an
    // average taken before squaring discards precisely the periapsis spike
    // that dominates the emission.  Both particles are summed at each phase
    // node BEFORE squaring, so the coherence between mu1'' and mu2'' -- and
    // its correlation with orbital phase -- is kept; summing two separately
    // averaged powers would discard exactly that cross term.
    //
    // omega' is taken on the same periodic node ring, by central difference in
    // eccentric anomaly and the chain rule dE/dt = n/(1-e cos E): the nodes
    // are uniform in E and the ring is closed, so the difference is spectrally
    // clean and needs no extra field evaluations.
    const double nodeCount=static_cast<double>(phaseNodes);
    const double eccentricAnomalyStep=2.0*pi/nodeCount;
    double secondDerivativeSquaredSum=0.0;
    for(int node=0;node<phaseNodes;++node) {
        const std::size_t here=static_cast<std::size_t>(node);
        const std::size_t next=static_cast<std::size_t>(
            (node+1)%phaseNodes);
        const std::size_t previous=static_cast<std::size_t>(
            (node+phaseNodes-1)%phaseNodes);
        const double timeWeight=timeWeightNodes[here];
        if(!(timeWeight>0.0)) return result;
        const double eccentricAnomalyRate=meanMotion/timeWeight;
        const Vec3 firstOmegaRate=
            (firstOmegaNodes[next]-firstOmegaNodes[previous])
                *(eccentricAnomalyRate/(2.0*eccentricAnomalyStep));
        const Vec3 secondOmegaRate=
            (secondOmegaNodes[next]-secondOmegaNodes[previous])
                *(eccentricAnomalyRate/(2.0*eccentricAnomalyStep));
        const Vec3 totalSecondDerivative=
            cross(firstOmegaRate,firstDipole)
            +cross(firstOmegaNodes[here],
                cross(firstOmegaNodes[here],firstDipole))
            +cross(secondOmegaRate,secondDipole)
            +cross(secondOmegaNodes[here],
                cross(secondOmegaNodes[here],secondDipole));
        if(!isFinite(totalSecondDerivative)) return result;
        secondDerivativeSquaredSum+=
            totalSecondDerivative.squaredNorm()*timeWeight;
    }
    result.coherentSecondDerivativeSquared=
        secondDerivativeSquaredSum/nodeCount;

    result.averagedOrbitalFrequency=
        orbitalFrequencySum/static_cast<double>(phaseNodes);
    result.first=firstOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.second=secondOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.firstExternal=
        firstExternalOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.secondExternal=
        secondExternalOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.valid=isFinite(result.first)&&isFinite(result.second)
        &&isFinite(result.firstExternal)&&isFinite(result.secondExternal)
        &&std::isfinite(result.coherentSecondDerivativeSquared)
        &&result.coherentSecondDerivativeSquared>=0.0
        &&std::isfinite(result.averagedOrbitalFrequency)
        &&result.averagedOrbitalFrequency>0.0;
    return result;
}

SecularSpinOrbitAdvance advanceCoupledSecularSpinOrbit(
        const SecularSpinOrbitState& initial,double semiMajorAxis,
        double reducedMass,double elapsedTime,
        double maximumRotationPerSubstep=0.05,
        int maximumSubsteps=65536) {
    SecularSpinOrbitAdvance result;
    result.state=initial;
    if(isFinite(initial.orbitalAngularMomentum)
       &&isFinite(initial.periapsisDirection)
       &&initial.orbitalAngularMomentum.norm()>0.0) {
        result.state.periapsisDirection=orbitPlaneDirection(
            initial.orbitalAngularMomentum,initial.periapsisDirection);
    }
    if(elapsedTime==0.0) {
        result.completed=isFinite(initial.orbitalAngularMomentum)
            &&isFinite(initial.firstDipole)&&isFinite(initial.secondDipole)
            &&std::isfinite(initial.zeroPointPhase)
            &&isFinite(initial.periapsisDirection)
            &&result.state.periapsisDirection.norm()>0.0;
        result.relativeAngularMomentumResidual=result.completed?0.0:
            std::numeric_limits<double>::infinity();
        return result;
    }
    if(!(elapsedTime>0.0)||!(semiMajorAxis>0.0)||!(reducedMass>0.0)
       ||!(maximumRotationPerSubstep>0.0)||maximumSubsteps<1
       ||!std::isfinite(elapsedTime)||!std::isfinite(semiMajorAxis)
       ||!std::isfinite(reducedMass)
       ||!std::isfinite(maximumRotationPerSubstep)
       ||!isFinite(initial.orbitalAngularMomentum)
       ||!isFinite(initial.firstDipole)||!isFinite(initial.secondDipole)
       ||!std::isfinite(initial.zeroPointPhase)
       ||!isFinite(initial.periapsisDirection)
       ||!(result.state.periapsisDirection.norm()>0.0)) {
        return result;
    }
    const double firstGyromagneticRatio=firstGyromagneticRatioOf();
    const double secondGyromagneticRatio=secondGyromagneticRatioOf();
    if(!(std::abs(firstGyromagneticRatio)>0.0)
       ||!(std::abs(secondGyromagneticRatio)>0.0)) return result;

    const auto spinTotal=[&](const Vec3& first,const Vec3& second) {
        return first/firstGyromagneticRatio
            +second/secondGyromagneticRatio;
    };
    // The ZPF phase advances at the orbit-averaged OSCULATING orbital
    // frequency, which is what the mechanical path integrates one step at a
    // time (electrodynamics.hpp advances State::zeroPointPhase by
    // osculatingOrbitalFrequency(s)*dt), and what every node of the
    // quadrature tunes its own ZPF band to.
    //
    // This used to advance by the MEAN MOTION instead, with a comment
    // asserting that was "exactly the quantity osculatingOrbitalFrequency*dt
    // advances State::zeroPointPhase by on the mechanical path".  It is not,
    // and only a circle hides the difference: osculatingOrbitalFrequency is
    // built from the instantaneous separation, sqrt(K/(mu r^3)), so it goes
    // as r^-3/2 while the mean motion is fixed at a^-3/2.  Around the
    // validation's own eccentric reference (e^2=0.945) the instantaneous
    // value runs from 0.361 n at apoapsis to 214.7 n at periapsis, a swing of
    // 595x, and its TIME AVERAGE is 1.5978 n -- so the secular path was
    // accumulating ZPF phase ~60% slower than the mechanical path does over
    // the same elapsed time, putting the two on different realizations of the
    // same field.  Taking the average from the same orbit walk that samples
    // the field costs nothing and makes them agree by construction; at e=0 it
    // reduces to the mean motion identically, so circular orbits are
    // unchanged.
    // Taken per substep from the rates below rather than computed here: it
    // depends on the eccentricity, which the transported orbital angular
    // momentum changes as the substeps proceed.
    const Vec3 initialAngularMomentum=
        initial.orbitalAngularMomentum
        +spinTotal(initial.firstDipole,initial.secondDipole);
    Vec3 transportedAngularMomentum=initialAngularMomentum;
    const double angularMomentumScale=std::max(
        initial.orbitalAngularMomentum.norm()
            +initial.firstDipole.norm()/std::abs(firstGyromagneticRatio)
            +initial.secondDipole.norm()/std::abs(secondGyromagneticRatio),
        1.0e-300);

    double advanced=0.0;
    while(advanced<elapsedTime&&result.substeps<maximumSubsteps) {
        const double remaining=elapsedTime-advanced;
        const OrbitAveragedBmtAngularVelocities startRates=
            orbitAveragedBmtAngularVelocities(
                semiMajorAxis,result.state.orbitalAngularMomentum,
                result.state.firstDipole,result.state.secondDipole,
                reducedMass,result.state.zeroPointPhase,
                result.state.periapsisDirection);
        if(!startRates.valid) return result;
        const double startSpeed=std::max(
            startRates.first.norm(),startRates.second.norm());
        double dt=remaining;
        if(startSpeed>0.0)
            dt=std::min(dt,maximumRotationPerSubstep/startSpeed);
        if(!(dt>0.0)||!std::isfinite(dt)) return result;

        // Explicit midpoint in the slowly changing precession field, with
        // exact rotations for both the predictor and the accepted update.
        // If the midpoint field rises sharply, shrink and retry before
        // accepting so the angle bound remains a real bound, not an estimate.
        const auto externalTorque=[&](
                const OrbitAveragedBmtAngularVelocities& rates,
                const Vec3& first,const Vec3& second) {
            return cross(rates.firstExternal,first)/firstGyromagneticRatio
                +cross(rates.secondExternal,second)/secondGyromagneticRatio;
        };
        const Vec3 startExternalTorque=
            externalTorque(startRates,result.state.firstDipole,
                            result.state.secondDipole);
        Vec3 firstMid,secondMid,orbitalMid,midpointAngularMomentum;
        OrbitAveragedBmtAngularVelocities midpointRates;
        for(int retry=0;retry<24;++retry) {
            firstMid=rotateDipoleByAngularVelocity(
                result.state.firstDipole,startRates.first,0.5*dt);
            secondMid=rotateDipoleByAngularVelocity(
                result.state.secondDipole,startRates.second,0.5*dt);
            midpointAngularMomentum=transportedAngularMomentum
                +startExternalTorque*(0.5*dt);
            orbitalMid=midpointAngularMomentum-spinTotal(firstMid,secondMid);
            const Vec3 periapsisMid=transportOrbitPlaneDirection(
                result.state.periapsisDirection,
                result.state.orbitalAngularMomentum,orbitalMid);
            midpointRates=orbitAveragedBmtAngularVelocities(
                semiMajorAxis,orbitalMid,firstMid,secondMid,reducedMass,
                result.state.zeroPointPhase
                    +startRates.averagedOrbitalFrequency*(0.5*dt),
                periapsisMid);
            if(!midpointRates.valid) return result;
            const double midpointAngle=dt*std::max(
                midpointRates.first.norm(),midpointRates.second.norm());
            if(midpointAngle<=1.05*maximumRotationPerSubstep) break;
            dt*=0.5;
            if(retry==23||!(dt>0.0)) return result;
        }

        const Vec3 firstBefore=result.state.firstDipole;
        const Vec3 secondBefore=result.state.secondDipole;
        const Vec3 firstAfter=rotateDipoleByAngularVelocity(
            firstBefore,midpointRates.first,dt);
        const Vec3 secondAfter=rotateDipoleByAngularVelocity(
            secondBefore,midpointRates.second,dt);
        if(!isFinite(firstAfter)||!isFinite(secondAfter)) return result;
        const Vec3 midpointExternalTorque=
            externalTorque(midpointRates,firstMid,secondMid);
        const Vec3 angularMomentumAfter=
            transportedAngularMomentum+midpointExternalTorque*dt;
        const Vec3 orbitalAfter=
            angularMomentumAfter-spinTotal(firstAfter,secondAfter);
        if(!isFinite(orbitalAfter)||!(orbitalAfter.norm()>1.0e-300))
            return result;
        const Vec3 periapsisAfter=transportOrbitPlaneDirection(
            result.state.periapsisDirection,
            result.state.orbitalAngularMomentum,orbitalAfter);
        if(!(periapsisAfter.norm()>0.0)) return result;

        // CREM_DEBUG_ALIGN: the mutual angle obeys
        // d(mu1.mu2)/dt = (omega1-omega2).(mu1 x mu2), so an exactly collinear
        // pair is a FIXED POINT of the transport whichever way it is collinear.
        // What separates para from ortho is therefore not which torque acts but
        // whether that fixed point is stable, and the quantity that decides it
        // is the rate DIFFERENCE.  Printed against the rates themselves so a
        // difference that is merely small can be told from one that vanishes.
        if(std::getenv("CREM_DEBUG_ALIGN")) {
            static int alignSamples=0;
            if(alignSamples++<6) {
                const Vec3 rateDifference=
                    midpointRates.first-midpointRates.second;
                const Vec3 momentCross=cross(firstBefore,secondBefore);
                const double firstNorm=firstBefore.norm();
                const double secondNorm=secondBefore.norm();
                std::cerr<<"ALIGN cos="
                    <<dot(firstBefore,secondBefore)
                        /std::max(firstNorm*secondNorm,1.0e-300)
                    <<" |w1|="<<midpointRates.first.norm()
                    <<" |w2|="<<midpointRates.second.norm()
                    <<" |w1-w2|="<<rateDifference.norm()
                    <<" rel="<<rateDifference.norm()
                        /std::max(midpointRates.first.norm(),1.0e-300)
                    <<" |mu1xmu2|/mumu="<<momentCross.norm()
                        /std::max(firstNorm*secondNorm,1.0e-300)
                    <<" dcos/dt="<<dot(rateDifference,momentCross)
                        /std::max(firstNorm*secondNorm,1.0e-300)<<'\n';
            }
        }
        result.state.firstDipole=firstAfter;
        result.state.secondDipole=secondAfter;
        result.state.orbitalAngularMomentum=orbitalAfter;
        result.state.periapsisDirection=periapsisAfter;
        // Midpoint rate, matching how the rotations themselves are accepted.
        result.state.zeroPointPhase+=
            midpointRates.averagedOrbitalFrequency*dt;
        transportedAngularMomentum=angularMomentumAfter;
        result.maximumSubstepAngle=std::max(result.maximumSubstepAngle,
            dt*std::max(midpointRates.first.norm(),
                        midpointRates.second.norm()));
        advanced+=dt;
        ++result.substeps;
        if(elapsedTime-advanced
           <=16.0*std::numeric_limits<double>::epsilon()*elapsedTime) {
            advanced=elapsedTime;
        }
    }

    const Vec3 finalAngularMomentum=result.state.orbitalAngularMomentum
        +spinTotal(result.state.firstDipole,result.state.secondDipole);
    result.relativeAngularMomentumResidual=
        (finalAngularMomentum-transportedAngularMomentum).norm()
        /angularMomentumScale;
    result.externalAngularMomentumTransfer=
        transportedAngularMomentum-initialAngularMomentum;
    result.completed=advanced>=elapsedTime
        &&std::isfinite(result.relativeAngularMomentumResidual)
        &&result.maximumSubstepAngle
            <=1.051*maximumRotationPerSubstep;
    return result;
}
