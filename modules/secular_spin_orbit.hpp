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
        Vec3 periapsisDirection={}) {
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
    // two avoid the aliasing seen with unrelated node counts: a circle keeps
    // the historical 64 evaluations, e^2=0.945 receives 256, and the cost is
    // bounded at 512 near the parabolic limit.
    const double targetPhaseNodes=
        32.0/std::sqrt(std::max(1.0-eccentricity,1.0e-12));
    int phaseNodes=64;
    while(phaseNodes<512&&phaseNodes<targetPhaseNodes) phaseNodes*=2;
    result.phaseNodes=phaseNodes;
    Vec3 firstOmegaSum,secondOmegaSum;
    Vec3 firstExternalOmegaSum,secondExternalOmegaSum;
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
        firstOmegaSum+=thomasBmtEffectiveField(
            sample.firstVelocity,fields.atFirst,firstGFactor)
                *(-firstCharge/firstMass*timeWeight);
        secondOmegaSum+=thomasBmtEffectiveField(
            sample.secondVelocity,fields.atSecond,secondGFactor)
                *(-secondCharge/secondMass*timeWeight);
        firstExternalOmegaSum+=thomasBmtEffectiveField(
            sample.firstVelocity,externalAtFirst,firstGFactor)
                *(-firstCharge/firstMass*timeWeight);
        secondExternalOmegaSum+=thomasBmtEffectiveField(
            sample.secondVelocity,externalAtSecond,secondGFactor)
                *(-secondCharge/secondMass*timeWeight);
    }
    result.first=firstOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.second=secondOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.firstExternal=
        firstExternalOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.secondExternal=
        secondExternalOmegaSum*(1.0/static_cast<double>(phaseNodes));
    result.valid=isFinite(result.first)&&isFinite(result.second)
        &&isFinite(result.firstExternal)&&isFinite(result.secondExternal);
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
    // Same mean motion orbitAveragedBmtAngularVelocities computes
    // internally (semiMajorAxis is fixed for this whole call, so it is
    // fixed here too): the ZPF phase advances at the pair's real orbital
    // frequency, exactly the quantity osculatingOrbitalFrequency*dt
    // advances State::zeroPointPhase by on the mechanical path -- not left
    // at the initial value for the whole checkpoint, let alone frozen at 0.
    const double meanMotion=std::sqrt(pairCoulombStrength
        /(reducedMass*semiMajorAxis*semiMajorAxis*semiMajorAxis));
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
                result.state.zeroPointPhase+meanMotion*(0.5*dt),
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

        result.state.firstDipole=firstAfter;
        result.state.secondDipole=secondAfter;
        result.state.orbitalAngularMomentum=orbitalAfter;
        result.state.periapsisDirection=periapsisAfter;
        result.state.zeroPointPhase+=meanMotion*dt;
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
