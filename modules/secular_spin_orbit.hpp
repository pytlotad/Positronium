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
    bool valid=false;
};

struct SecularSpinOrbitState {
    Vec3 orbitalAngularMomentum;
    Vec3 firstDipole;
    Vec3 secondDipole;
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

OrbitAveragedBmtAngularVelocities orbitAveragedBmtAngularVelocities(
        double semiMajorAxis,const Vec3& orbitalAngularMomentum,
        const Vec3& firstDipole,const Vec3& secondDipole,
        double reducedMass) {
    OrbitAveragedBmtAngularVelocities result;
    if(!(semiMajorAxis>0.0)||!(reducedMass>0.0)
       ||!std::isfinite(semiMajorAxis)||!std::isfinite(reducedMass)
       ||!isFinite(orbitalAngularMomentum)||!isFinite(firstDipole)
       ||!isFinite(secondDipole)) return result;
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
    const Vec3 seed=std::abs(normal.x)<0.9
        ?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    Vec3 radialHat=cross(seed,normal);
    const double radialNorm=radialHat.norm();
    if(!(radialNorm>0.0)) return result;
    radialHat=radialHat/radialNorm;
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
    // as e approaches one, so refine only eccentric orbits: a circle keeps the
    // historical 64 evaluations, e^2=0.945 receives 192, and the cost is
    // bounded near the parabolic limit.
    const int phaseNodes=std::clamp(static_cast<int>(std::ceil(
        128.0/std::sqrt(std::max(1.0-eccentricity,1.0e-12)))),64,512);
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
        sample.firstDipole=firstDipole;
        sample.secondDipole=secondDipole;
        sample.firstProperDipole=firstDipole;
        sample.secondProperDipole=secondDipole;
        const StateHistory sampleHistory{State{sample}};
        const LocalElectromagneticFields fields=
            localRelativisticFields(sample,sampleHistory);
        ElectromagneticField externalAtFirst{{},gExternalMagneticField};
        ElectromagneticField externalAtSecond{{},gExternalMagneticField};
        if(gZeroPointField.active()) {
            Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
            const double orbitalFrequency=osculatingOrbitalFrequency(sample);
            gZeroPointField.sample(sample.firstPosition,orbitalFrequency,
                sample.zeroPointPhase,firstElectric,firstMagnetic);
            gZeroPointField.sample(sample.secondPosition,orbitalFrequency,
                sample.zeroPointPhase,secondElectric,secondMagnetic);
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
    if(elapsedTime==0.0) {
        result.completed=isFinite(initial.orbitalAngularMomentum)
            &&isFinite(initial.firstDipole)&&isFinite(initial.secondDipole);
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
       ||!isFinite(initial.firstDipole)||!isFinite(initial.secondDipole)) {
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
                result.state.firstDipole,result.state.secondDipole,reducedMass);
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
            midpointRates=orbitAveragedBmtAngularVelocities(
                semiMajorAxis,orbitalMid,firstMid,secondMid,reducedMass);
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

        result.state.firstDipole=firstAfter;
        result.state.secondDipole=secondAfter;
        result.state.orbitalAngularMomentum=orbitalAfter;
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
