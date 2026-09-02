#pragma once

// End-to-end numerical validation of the optional Maxwell field backend.

enum class StatisticalValidationProfile { Small, Publication };

int runMaxwellSelfTest(
        StatisticalValidationProfile statisticalProfile=
            StatisticalValidationProfile::Small) {
    const auto benchmarkStart=std::chrono::steady_clock::now();

    // Fixed-seed distribution checks for the stochastic photon sector.  The
    // small profile is cheap enough for every regression run; publication is
    // deliberately much larger and tightens the tolerances approximately as
    // sqrt(N).  Neither profile accepts a user seed: changing the random
    // sample must be an explicit reviewable source change, not ambient state.
    const bool publicationStatistics=
        statisticalProfile==StatisticalValidationProfile::Publication;
    const std::uint64_t statisticalSampleCount=
        publicationStatistics?1000000ULL:4096ULL;
    std::uint64_t directionStream=0x4352454d5f45315fULL;
    std::uint64_t exponentialStream=0x4352454d5f504f49ULL;
    double cosineSum=0.0;
    double cosineSquaredSum=0.0;
    double azimuthCosineSum=0.0;
    double azimuthSineSum=0.0;
    double exponentialSum=0.0;
    double exponentialSquaredSum=0.0;
    for(std::uint64_t sample=0;sample<statisticalSampleCount;++sample) {
        const Vec3 direction=sampleRotatingDipolePhotonDirection(
            Vec3{0.0,0.0,1.0},directionStream);
        cosineSum+=direction.z;
        cosineSquaredSum+=direction.z*direction.z;
        azimuthCosineSum+=direction.x;
        azimuthSineSum+=direction.y;
        const double exponential=drawExponentialUnit(exponentialStream);
        exponentialSum+=exponential;
        exponentialSquaredSum+=exponential*exponential;
    }
    const double inverseStatisticalSampleCount=
        1.0/static_cast<double>(statisticalSampleCount);
    const double photonCosineMean=cosineSum*inverseStatisticalSampleCount;
    const double photonCosineSecondMoment=
        cosineSquaredSum*inverseStatisticalSampleCount;
    const double photonAzimuthVectorMean=std::hypot(
        azimuthCosineSum,azimuthSineSum)*inverseStatisticalSampleCount;
    const double exponentialMean=
        exponentialSum*inverseStatisticalSampleCount;
    const double exponentialVariance=
        exponentialSquaredSum*inverseStatisticalSampleCount
        -exponentialMean*exponentialMean;
    std::uint64_t replayA=0x4352454d5f524550ULL;
    std::uint64_t replayB=replayA;
    bool statisticalReplayExact=true;
    for(int sample=0;sample<1024;++sample) {
        const Vec3 a=sampleRotatingDipolePhotonDirection(
            Vec3{0.0,0.0,1.0},replayA);
        const Vec3 b=sampleRotatingDipolePhotonDirection(
            Vec3{0.0,0.0,1.0},replayB);
        statisticalReplayExact=statisticalReplayExact
            &&a.x==b.x&&a.y==b.y&&a.z==b.z;
    }
    const double cosineMeanTolerance=publicationStatistics?0.004:0.065;
    const double cosineSecondMomentTolerance=
        publicationStatistics?0.002:0.032;
    const double azimuthMeanTolerance=publicationStatistics?0.002:0.035;
    const double exponentialMeanTolerance=publicationStatistics?0.006:0.10;
    const double exponentialVarianceTolerance=
        publicationStatistics?0.012:0.20;
    const bool stochasticStatisticsOk=statisticalReplayExact
        &&std::abs(photonCosineMean)<cosineMeanTolerance
        &&std::abs(photonCosineSecondMoment-0.4)
            <cosineSecondMomentTolerance
        &&photonAzimuthVectorMean<azimuthMeanTolerance
        &&std::abs(exponentialMean-1.0)<exponentialMeanTolerance
        &&std::abs(exponentialVariance-1.0)<exponentialVarianceTolerance;

    // Experiment 5 models a circular 2D Gaussian transverse beam.  Exercise
    // its production sampler directly: b/sigma must have Rayleigh moments
    // E[b/sigma]=sqrt(pi/2), E[(b/sigma)^2]=2 and no preferred azimuth.  The
    // old folded one-dimensional Gaussian fails both radial moments by a wide
    // margin (sqrt(2/pi) and 1 respectively).
    constexpr double impactTestSigma=1.0;
    constexpr double impactTestMaximum=5.0*impactTestSigma;
    std::mt19937_64 impactRandom(0x4352454d5f494d50ULL);
    double normalizedImpactSum=0.0;
    double normalizedImpactSquaredSum=0.0;
    double impactAzimuthCosineSum=0.0;
    double impactAzimuthSineSum=0.0;
    bool impactSamplesValid=true;
    for(std::uint64_t sample=0;sample<statisticalSampleCount;++sample) {
        const IsotropicGaussianImpactSample impact=
            sampleIsotropicGaussianImpact(
                impactRandom,impactTestSigma,impactTestMaximum);
        impactSamplesValid=impactSamplesValid
            &&impact.valid(impactTestMaximum);
        if (!impact.valid(impactTestMaximum)) continue;
        const double normalizedImpact=impact.impactParameter/impactTestSigma;
        normalizedImpactSum+=normalizedImpact;
        normalizedImpactSquaredSum+=normalizedImpact*normalizedImpact;
        if (impact.impactParameter>0.0) {
            impactAzimuthCosineSum+=
                impact.transverseY/impact.impactParameter;
            impactAzimuthSineSum+=
                impact.transverseZ/impact.impactParameter;
        }
    }
    const double normalizedImpactMean=
        normalizedImpactSum*inverseStatisticalSampleCount;
    const double normalizedImpactSecondMoment=
        normalizedImpactSquaredSum*inverseStatisticalSampleCount;
    const double impactAzimuthVectorMean=std::hypot(
        impactAzimuthCosineSum,impactAzimuthSineSum)
        *inverseStatisticalSampleCount;
    std::mt19937_64 impactReplayA(0x4352454d5f495250ULL);
    std::mt19937_64 impactReplayB(0x4352454d5f495250ULL);
    bool impactReplayExact=true;
    for(int sample=0;sample<1024;++sample) {
        const IsotropicGaussianImpactSample a=sampleIsotropicGaussianImpact(
            impactReplayA,impactTestSigma,impactTestMaximum);
        const IsotropicGaussianImpactSample b=sampleIsotropicGaussianImpact(
            impactReplayB,impactTestSigma,impactTestMaximum);
        impactReplayExact=impactReplayExact
            &&a.transverseY==b.transverseY
            &&a.transverseZ==b.transverseZ
            &&a.impactParameter==b.impactParameter;
    }
    const double impactMeanTolerance=publicationStatistics?0.003:0.040;
    const double impactSecondMomentTolerance=
        publicationStatistics?0.006:0.12;
    const double impactAzimuthTolerance=
        publicationStatistics?0.002:0.040;
    const bool impactParameterProfileOk=impactSamplesValid
        &&impactReplayExact
        &&std::abs(normalizedImpactMean-std::sqrt(pi/2.0))
            <impactMeanTolerance
        &&std::abs(normalizedImpactSecondMoment-2.0)
            <impactSecondMomentTolerance
        &&impactAzimuthVectorMean<impactAzimuthTolerance;

    // Cheap production-kinematics regressions.  These call the same builders
    // as experiments 3--5, so they detect role-dependent energy splits and a
    // return to Galilean velocity addition without running a trajectory.
    const double twoBodyTestEnergy=20.0*eCharge;
    const double twoBodyMatchingRadius=1.0e-9;
    const double twoBodyImpact=0.2*twoBodyMatchingRadius;
    const double twoBodyPotentialGain=pairCoulombStrength/twoBodyMatchingRadius;
    const Vec3 twoBodyRadial{-std::sqrt(1.0-0.2*0.2),0.2,0.0};
    const Vec3 twoBodyTangent{0.2,std::sqrt(1.0-0.2*0.2),0.0};
    const two_body::IncomingTwoBodyKinematics protonElectronIncoming=
        two_body::incomingTwoBodyKinematics(
            twoBodyTestEnergy,twoBodyPotentialGain,
            twoBodyImpact,twoBodyMatchingRadius,
            twoBodyRadial,twoBodyTangent,proton.mass,electron.mass);
    const two_body::IncomingTwoBodyKinematics electronProtonIncoming=
        two_body::incomingTwoBodyKinematics(
            twoBodyTestEnergy,twoBodyPotentialGain,
            twoBodyImpact,twoBodyMatchingRadius,
            twoBodyRadial,twoBodyTangent,electron.mass,proton.mass);
    const two_body::FreeTwoBodyKinematics protonElectronFree=
        two_body::freeTwoBodyKinematics(
            protonElectronIncoming.firstVelocity,proton.mass,
            protonElectronIncoming.secondVelocity,electron.mass);
    const two_body::FreeTwoBodyKinematics electronProtonFree=
        two_body::freeTwoBodyKinematics(
            electronProtonIncoming.firstVelocity,electron.mass,
            electronProtonIncoming.secondVelocity,proton.mass);
    const double twoBodyRoleEnergyResidual=std::max(
        std::abs((protonElectronFree.centreOfMomentumKineticEnergy
                  -twoBodyPotentialGain-twoBodyTestEnergy)/twoBodyTestEnergy),
        std::abs((electronProtonFree.centreOfMomentumKineticEnergy
                  -twoBodyPotentialGain-twoBodyTestEnergy)/twoBodyTestEnergy));
    const double twoBodyRoleMomentumResidual=std::max(
        (protonElectronIncoming.firstFrame.momentum
         +protonElectronIncoming.secondFrame.momentum).norm()
            /std::max(2.0*protonElectronIncoming.finiteRadius.momentumMagnitude,
                      1.0e-300),
        (electronProtonIncoming.firstFrame.momentum
         +electronProtonIncoming.secondFrame.momentum).norm()
            /std::max(2.0*electronProtonIncoming.finiteRadius.momentumMagnitude,
                      1.0e-300));
    const double twoBodyRoleSwapResidual=std::max({
        std::abs(protonElectronIncoming.asymptotic.momentumMagnitude
                -electronProtonIncoming.asymptotic.momentumMagnitude)
            /std::max(protonElectronIncoming.asymptotic.momentumMagnitude,
                      1.0e-300),
        std::abs(protonElectronIncoming.finiteRadius.firstSpeed
                -electronProtonIncoming.finiteRadius.secondSpeed)
            /std::max(protonElectronIncoming.finiteRadius.firstSpeed,1.0),
        std::abs(protonElectronIncoming.finiteRadius.secondSpeed
                -electronProtonIncoming.finiteRadius.firstSpeed)
            /std::max(protonElectronIncoming.finiteRadius.secondSpeed,1.0)});

    const double electronRestEnergy=electron.mass*c*c;
    const two_body::HeadOnLabKinematics asymmetricLab=
        two_body::headOnLabKinematics(
            2.0*electronRestEnergy,electron.mass,
            0.01*electronRestEnergy,positron.mass,Vec3{1.0,0.0,0.0});
    const two_body::IncomingTwoBodyKinematics boostedIncoming=
        two_body::incomingTwoBodyKinematics(
            asymmetricLab.pair.centreOfMomentumKineticEnergy,
            twoBodyPotentialGain,twoBodyImpact,twoBodyMatchingRadius,
            twoBodyRadial,twoBodyTangent,electron.mass,positron.mass,
            asymmetricLab.pair.centreOfMomentumVelocity);
    const two_body::FreeTwoBodyKinematics boostedFree=
        two_body::freeTwoBodyKinematics(
            boostedIncoming.firstVelocity,electron.mass,
            boostedIncoming.secondVelocity,positron.mass);
    const two_body::ParticleFourMomentum recoveredFirst=
        two_body::boostFourMomentum(boostedIncoming.firstFrame,
            asymmetricLab.pair.centreOfMomentumVelocity*(-1.0));
    const two_body::ParticleFourMomentum recoveredSecond=
        two_body::boostFourMomentum(boostedIncoming.secondFrame,
            asymmetricLab.pair.centreOfMomentumVelocity*(-1.0));
    const Vec3 unequalFrameVelocity{0.55*c,0.10*c,0.0};
    const two_body::IncomingTwoBodyKinematics unequalBoostedIncoming=
        two_body::incomingTwoBodyKinematics(
            twoBodyTestEnergy,twoBodyPotentialGain,
            twoBodyImpact,twoBodyMatchingRadius,
            twoBodyRadial,twoBodyTangent,proton.mass,electron.mass,
            unequalFrameVelocity);
    const two_body::FreeTwoBodyKinematics unequalBoostedFree=
        two_body::freeTwoBodyKinematics(
            unequalBoostedIncoming.firstVelocity,proton.mass,
            unequalBoostedIncoming.secondVelocity,electron.mass);
    const two_body::ParticleFourMomentum unequalRecoveredFirst=
        two_body::boostFourMomentum(unequalBoostedIncoming.firstFrame,
            unequalFrameVelocity*(-1.0));
    const two_body::ParticleFourMomentum unequalRecoveredSecond=
        two_body::boostFourMomentum(unequalBoostedIncoming.secondFrame,
            unequalFrameVelocity*(-1.0));
    const double twoBodyBoostEnergyResidual=std::max(std::abs(
        boostedFree.centreOfMomentumKineticEnergy
        -boostedIncoming.finiteRadius.kineticEnergy)
        /std::max(boostedIncoming.finiteRadius.kineticEnergy,1.0e-300),
        std::abs(unequalBoostedFree.centreOfMomentumKineticEnergy
                 -unequalBoostedIncoming.finiteRadius.kineticEnergy)
        /std::max(unequalBoostedIncoming.finiteRadius.kineticEnergy,1.0e-300));
    const double twoBodyBoostVelocityResidual=std::max(
        (boostedFree.centreOfMomentumVelocity
         -asymmetricLab.pair.centreOfMomentumVelocity).norm()/c,
        (unequalBoostedFree.centreOfMomentumVelocity
         -unequalFrameVelocity).norm()/c);
    const double twoBodyInverseBoostResidual=std::max({
        (recoveredFirst.momentum
         -boostedIncoming.firstCentreOfMomentum.momentum).norm()
            /std::max(boostedIncoming.finiteRadius.momentumMagnitude,1.0e-300),
        (recoveredSecond.momentum
         -boostedIncoming.secondCentreOfMomentum.momentum).norm()
            /std::max(boostedIncoming.finiteRadius.momentumMagnitude,1.0e-300),
        (unequalRecoveredFirst.momentum
         -unequalBoostedIncoming.firstCentreOfMomentum.momentum).norm()
            /std::max(unequalBoostedIncoming.finiteRadius.momentumMagnitude,
                      1.0e-300),
        (unequalRecoveredSecond.momentum
         -unequalBoostedIncoming.secondCentreOfMomentum.momentum).norm()
            /std::max(unequalBoostedIncoming.finiteRadius.momentumMagnitude,
                      1.0e-300)});
    const double commonNearLightBeta=1.0-std::ldexp(1.0,-40);
    const Vec3 commonNearLightVelocity{commonNearLightBeta*c,0.0,0.0};
    const two_body::FreeTwoBodyKinematics comovingNearLightPair=
        two_body::freeTwoBodyKinematics(
            commonNearLightVelocity,proton.mass,
            commonNearLightVelocity,electron.mass);
    const double comovingInternalEnergyEv=
        comovingNearLightPair.centreOfMomentumKineticEnergy/eCharge;
    const bool allInitialVelocitiesCausal=
        two_body::subluminal(protonElectronIncoming.firstVelocity)
        &&two_body::subluminal(protonElectronIncoming.secondVelocity)
        &&two_body::subluminal(electronProtonIncoming.firstVelocity)
        &&two_body::subluminal(electronProtonIncoming.secondVelocity)
        &&two_body::subluminal(boostedIncoming.firstVelocity)
        &&two_body::subluminal(boostedIncoming.secondVelocity)
        &&two_body::subluminal(unequalBoostedIncoming.firstVelocity)
        &&two_body::subluminal(unequalBoostedIncoming.secondVelocity)
        &&two_body::subluminal(commonNearLightVelocity);
    const double maximumInitialBeta=std::sqrt(std::max({
        protonElectronIncoming.firstVelocity.squaredNorm(),
        protonElectronIncoming.secondVelocity.squaredNorm(),
        electronProtonIncoming.firstVelocity.squaredNorm(),
        electronProtonIncoming.secondVelocity.squaredNorm(),
        boostedIncoming.firstVelocity.squaredNorm(),
        boostedIncoming.secondVelocity.squaredNorm(),
        unequalBoostedIncoming.firstVelocity.squaredNorm(),
        unequalBoostedIncoming.secondVelocity.squaredNorm(),
        commonNearLightVelocity.squaredNorm()}))/c;

    const CovariantExtendedBody covarianceBody{firstCharge,firstMass,
        chargeCloudRestRadius,{},Vec3{0.31*c,-0.07*c,0.04*c},
        Vec3{0,0,firstMagneticMoment}};
    const double electromagneticMassFraction=
        covarianceBody.electromagneticMass()/covarianceBody.physicalMass;
    const double bareMassFraction=
        covarianceBody.bareMatterMass()/covarianceBody.physicalMass;
    MaxwellBlock selfForceBlock(24,0.5*chargeCloudRestRadius,{});
    double maximumSelfForceFraction=0.0;
    const double referenceForce=pairCoulombStrength
        /std::pow(2.5*chargeCloudRestRadius,2);
    for(int ix=0;ix<2;++ix) for(int iy=0;iy<2;++iy) for(int iz=0;iz<2;++iz) {
        const Vec3 offset{(ix?0.31:-0.19)*selfForceBlock.cellSize(),
                          (iy?0.27:-0.23)*selfForceBlock.cellSize(),
                          (iz?0.17:-0.33)*selfForceBlock.cellSize()};
        const auto [selfElectric,selfMagnetic]=selfForceBlock.numericalSelfField(
            {firstCharge,chargeCloudRestRadius},offset,Vec3{1.0e5,2.0e5,-0.5e5},240);
        const Vec3 selfForce=(selfElectric+cross(Vec3{1.0e5,2.0e5,-0.5e5},
                                                  selfMagnetic))*firstCharge;
        maximumSelfForceFraction=std::max(maximumSelfForceFraction,
                                          selfForce.norm()/referenceForce);
    }
    const FourVector restCurrent{c*firstCharge*covarianceBody.properShape({}),{}};
    const FourVector boostedCurrent=covarianceBody.freeFourCurrent({});
    const double currentInvariantResidual=std::abs(
        minkowskiDot(boostedCurrent,boostedCurrent)
        -minkowskiDot(restCurrent,restCurrent))
        /std::max(std::abs(minkowskiDot(restCurrent,restCurrent)),1.0e-300);
    const AntisymmetricTensor polarizationTensor=
        covarianceBody.polarizationMagnetization({});
    double antisymmetryResidual=0.0;
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
        antisymmetryResidual=std::max(antisymmetryResidual,std::abs(
            polarizationTensor.component[i][j]+polarizationTensor.component[j][i]));
    const RankTwoTensor matterTensor=covarianceBody.matterStressEnergy({});
    double stressSymmetryResidual=0.0;
    for(int i=0;i<4;++i) for(int j=0;j<4;++j)
        stressSymmetryResidual=std::max(stressSymmetryResidual,std::abs(
            matterTensor.component[i][j]-matterTensor.component[j][i]));
    MaxwellAmrHierarchy hierarchy({}, 24, 3);
    MaxwellBlock& finest = hierarchy.finest();
    const Vec3 separationVector{2.5*chargeCloudRestRadius, 0.0, 0.0};
    hierarchy.depositPair(separationVector*(-0.5), {},
                          separationVector*0.5, {});
    finest.projectElectricGaussConstraint(400);
    const double depositedCharge = finest.totalCharge();
    const double gaussBefore = finest.maximumElectricGaussResidual();
    const double relativeGaussResidual = gaussBefore
        / std::max(finest.maximumChargeGaussScale(), 1.0);

    finest.clearSources();
    finest.setVacuumPlaneWave(1.0e5);
    const double initialEnergy = finest.fieldEnergy();
    const double dt = 0.45*finest.courantTimeStep();
    for (int step=0; step<40; ++step) finest.advance(dt);
    const double finalEnergy = finest.fieldEnergy();
    const double relativeEnergyDrift = std::abs(finalEnergy-initialEnergy)
                                     / initialEnergy;
    const double magneticDivergence = finest.maximumMagneticDivergence();
    MaxwellAmrHierarchy coupledHierarchy({},24,3);
    MaxwellBlock& coupledField=coupledHierarchy.finest();
    State particles;
    particles.firstPosition=separationVector*(-0.5);
    particles.secondPosition=separationVector*0.5;
    particles.firstVelocity={0.0,2.0e5,0.0};
    particles.secondVelocity={0.0,-2.0e5,0.0};
    particles.firstDipole={0.0,0.0,firstMagneticMoment};
    particles.secondDipole={0.0,0.0,-secondMagneticMoment};
    coupledField.clearSources();
    coupledField.depositCloud({firstCharge,chargeCloudRestRadius},
        particles.firstPosition,particles.firstVelocity);
    coupledField.depositCloud({secondCharge,chargeCloudRestRadius},
        particles.secondPosition,particles.secondVelocity);
    coupledField.depositCovariantDipole(particles.firstPosition,
        particles.firstVelocity,particles.firstDipole);
    coupledField.depositCovariantDipole(particles.secondPosition,
        particles.secondVelocity,particles.secondDipole);
    coupledField.finalizeBoundInstantaneous();
    coupledField.projectElectricGaussConstraint(400);
    const auto [firstSelfElectric,firstSelfMagnetic]=
        coupledField.numericalSelfField({firstCharge,chargeCloudRestRadius},
            particles.firstPosition,particles.firstVelocity);
    const auto [secondSelfElectric,secondSelfMagnetic]=
        coupledField.numericalSelfField({secondCharge,chargeCloudRestRadius},
            particles.secondPosition,particles.secondVelocity);
    const ElectromagneticField firstSelfField{
        firstSelfElectric,firstSelfMagnetic};
    const ElectromagneticField secondSelfField{
        secondSelfElectric,secondSelfMagnetic};
    const DynamicSelfFieldCalibration dynamicSelfCalibration(
        coupledField,eCharge,120);
    const double phaseEpsilon=1.0e-8*coupledField.cellSize();
    const Vec3 phaseLeft=coupledField.centre()
        +Vec3{-0.5*coupledField.cellSize()+phaseEpsilon,
               0.17*coupledField.cellSize(),-0.23*coupledField.cellSize()};
    const Vec3 phaseRight=phaseLeft+Vec3{
        coupledField.cellSize()-2.0*phaseEpsilon,0,0};
    const Vec3 calibrationVelocity{0.12*c,-0.04*c,0.03*c};
    const ElectromagneticField phaseLeftField=dynamicSelfCalibration.field(
        coupledField,phaseLeft,calibrationVelocity,eCharge);
    const ElectromagneticField phaseRightField=dynamicSelfCalibration.field(
        coupledField,phaseRight,calibrationVelocity,eCharge);
    const ElectromagneticField oppositeChargeField=dynamicSelfCalibration.field(
        coupledField,phaseLeft,calibrationVelocity,-eCharge);
    const double selfCalibrationPhaseJump=std::max(
        (phaseLeftField.electric-phaseRightField.electric).norm()
            /std::max(phaseLeftField.electric.norm(),1.0),
        (phaseLeftField.magnetic-phaseRightField.magnetic).norm()
            /std::max(phaseLeftField.magnetic.norm(),1.0e-300));
    const double selfCalibrationChargeSymmetry=std::max(
        (phaseLeftField.electric+oppositeChargeField.electric).norm()
            /std::max(phaseLeftField.electric.norm(),1.0),
        (phaseLeftField.magnetic+oppositeChargeField.magnetic).norm()
            /std::max(phaseLeftField.magnetic.norm(),1.0e-300));
    const ParticleFieldTotals initialCoupledTotals=
        particleFieldTotals(particles,coupledField);
    const double coupledDt=0.175*coupledHierarchy.timeStep();
    double maximumRelativeContinuityResidual=0.0;
    double maximumRelativeLongitudinalCurl=0.0;
    double coupledEscapedEnergy=0.0;
    Vec3 coupledEscapedMomentum,coupledEscapedAngularMomentum;
    const RelativisticChargeCloud gridFirst{firstCharge,chargeCloudRestRadius};
    const RelativisticChargeCloud gridSecond{secondCharge,chargeCloudRestRadius};
    for(int step=0;step<24;++step) {
        const MaxwellBoundaryFlux beforeFlux=coupledField.boundaryFlux();
        const std::vector<double> previousCharge=
            coupledField.chargeDensitySnapshot();
        const std::vector<Vec3> previousPolarization=
            coupledField.polarizationSnapshot();
        State midpointParticles=particles;
        pushStateWithGridField(midpointParticles,coupledField,0.5*coupledDt,true,
            nullptr,firstSelfField,secondSelfField,
            &dynamicSelfCalibration,&dynamicSelfCalibration);
        MaxwellBlock midpointField=coupledField;
        midpointField.clearSources();
        midpointField.depositCloud(gridFirst,midpointParticles.firstPosition,
                                   midpointParticles.firstVelocity);
        midpointField.depositCloud(gridSecond,midpointParticles.secondPosition,
                                   midpointParticles.secondVelocity);
        midpointField.depositCovariantDipole(midpointParticles.firstPosition,
            midpointParticles.firstVelocity,midpointParticles.firstDipole);
        midpointField.depositCovariantDipole(midpointParticles.secondPosition,
            midpointParticles.secondVelocity,midpointParticles.secondDipole);
        midpointField.finalizeBoundCurrent(previousPolarization,0.5*coupledDt);
        midpointField.finalizeChargeConservingCurrent(
            previousCharge,0.5*coupledDt,240);
        midpointField.advance(0.5*coupledDt);

        State finalParticles=particles;
        pushStateWithGridField(finalParticles,midpointField,coupledDt,true,
            &midpointParticles,firstSelfField,secondSelfField,
            &dynamicSelfCalibration,&dynamicSelfCalibration);
        coupledField.clearSources();
        coupledField.depositChargeDensity(gridFirst,finalParticles.firstPosition);
        coupledField.depositChargeDensity(gridSecond,finalParticles.secondPosition);
        coupledField.depositConvectionCurrent(gridFirst,
            midpointParticles.firstPosition,midpointParticles.firstVelocity);
        coupledField.depositConvectionCurrent(gridSecond,
            midpointParticles.secondPosition,midpointParticles.secondVelocity);
        coupledField.depositCovariantDipole(midpointParticles.firstPosition,
            midpointParticles.firstVelocity,midpointParticles.firstDipole);
        coupledField.depositCovariantDipole(midpointParticles.secondPosition,
            midpointParticles.secondVelocity,midpointParticles.secondDipole);
        coupledField.finalizeBoundCurrent(previousPolarization,coupledDt);
        const ChargeConservingDepositResult deposition=
            coupledField.finalizeChargeConservingCurrent(
            previousCharge,coupledDt,240);
        const double continuityScale=std::max(
            coupledField.maximumChargeGaussScale()*epsilon0/coupledDt,1.0);
        maximumRelativeContinuityResidual=std::max(
            maximumRelativeContinuityResidual,
            deposition.continuityResidual/continuityScale);
        maximumRelativeLongitudinalCurl=std::max(
            maximumRelativeLongitudinalCurl,deposition.longitudinalCurl
                /std::max(coupledField.maximumCurrentGradientScale(),1.0));
        coupledField.advance(coupledDt);
        particles=finalParticles;
        const MaxwellBoundaryFlux afterFlux=coupledField.boundaryFlux();
        coupledEscapedEnergy+=0.5*(beforeFlux.energyRate+afterFlux.energyRate)
                            *coupledDt;
        coupledEscapedMomentum+=(beforeFlux.momentumRate+afterFlux.momentumRate)
                              *(0.5*coupledDt);
        coupledEscapedAngularMomentum+=(beforeFlux.angularMomentumRate
            +afterFlux.angularMomentumRate)*(0.5*coupledDt);
    }
    const ParticleFieldTotals finalCoupledTotals=
        particleFieldTotals(particles,coupledField);
    const double coupledEnergyClosure=(finalCoupledTotals.energy
        +coupledEscapedEnergy-initialCoupledTotals.energy)
        /std::max(std::abs(initialCoupledTotals.energy),1.0e-30);
    const double coupledMomentumScale=std::max(
        momentum(particles.firstVelocity,firstMass).norm()
       +momentum(particles.secondVelocity,secondMass).norm(),1.0e-40);
    const double coupledMomentumClosure=((finalCoupledTotals.momentum
        +coupledEscapedMomentum)-initialCoupledTotals.momentum).norm()
        /coupledMomentumScale;
    const double coupledAngularScale=std::max({
        initialCoupledTotals.angularMomentum.norm(),hbar,1.0e-40});
    const double coupledAngularClosure=((finalCoupledTotals.angularMomentum
        +coupledEscapedAngularMomentum)-initialCoupledTotals.angularMomentum).norm()
        /coupledAngularScale;
    const double coupledBeta=std::max(particles.firstVelocity.norm(),
                                      particles.secondVelocity.norm())/c;
    const double dipoleNormResidual=std::max(
        std::abs(particles.firstDipole.norm()/firstMagneticMoment-1.0),
        std::abs(particles.secondDipole.norm()/secondMagneticMoment-1.0));
    MaxwellBlock absorbingTest(32,0.5*chargeCloudRestRadius,{});
    absorbingTest.setPlaneWavePacketX(1.0e5,chargeCloudRestRadius);
    absorbingTest.enableConvolutionalPml(8,1.0e-10);
    const double absorbingInitialEnergy=absorbingTest.fieldEnergy();
    const double absorbingDt=0.4*absorbingTest.courantTimeStep();
    double integratedBoundaryEnergy=0.0;
    double outgoingCharacteristic=0.0,reflectedCharacteristic=0.0;
    Vec3 integratedBoundaryMomentum,integratedBoundaryAngularMomentum;
    for(int step=0;step<160;++step) {
        const auto [outgoingAtMonitor,incomingAtMonitor]=
            absorbingTest.xCharacteristicEnergy(16);
        if(step<45) outgoingCharacteristic+=outgoingAtMonitor*absorbingDt;
        if(step>=80) reflectedCharacteristic+=incomingAtMonitor*absorbingDt;
        const MaxwellBoundaryFlux beforeFlux=absorbingTest.boundaryFlux();
        absorbingTest.advance(absorbingDt);
        const MaxwellBoundaryFlux afterFlux=absorbingTest.boundaryFlux();
        integratedBoundaryEnergy+=0.5*(beforeFlux.energyRate+afterFlux.energyRate)
                                *absorbingDt;
        integratedBoundaryMomentum+=(beforeFlux.momentumRate+afterFlux.momentumRate)
                                  *(0.5*absorbingDt);
        integratedBoundaryAngularMomentum+=(beforeFlux.angularMomentumRate
            +afterFlux.angularMomentumRate)*(0.5*absorbingDt);
    }
    const double absorbedFraction=absorbingTest.absorbedBoundaryEnergy()
                                 /absorbingInitialEnergy;
    const double lateInteriorEnergyFraction=absorbingTest.interiorFieldEnergy(7)
                                           /absorbingInitialEnergy;
    const double cpmlReflection=reflectedCharacteristic
        /std::max(outgoingCharacteristic,1.0e-300);
    MaxwellBlock magnetizationTest(24,0.5*chargeCloudRestRadius,{});
    magnetizationTest.clearSources();
    magnetizationTest.depositMagneticDipole({},Vec3{0,0,bohrMagneton});
    magnetizationTest.finalizeMagnetizationCurrent();
    const double relativeMagnetizationDivergence=
        magnetizationTest.maximumCurrentDivergence()
        /std::max(magnetizationTest.maximumCurrentGradientScale(),1.0);
    const auto gaussResidualAt=[&](double spacing) {
        MaxwellBlock block(24,spacing,{});
        block.clearSources();
        // The pair's own charges, in the pair's own role order.  Writing -e
        // and +e mirrored the DEFAULT pair and inverted both roles for any
        // other; nothing is pushed here so the mirrored field gives the same
        // normalized residual either way, but this is the exact shape of the
        // hard-coding that did bite the Yee coupled continuity test, where a
        // pusher then moved the charge the deposit had contradicted.
        block.depositCloud({firstCharge,chargeCloudRestRadius},
                           separationVector*(-0.5),{});
        block.depositCloud({secondCharge,chargeCloudRestRadius},
                           separationVector*0.5,{});
        block.projectElectricGaussConstraint(400);
        return block.maximumElectricGaussResidual()
            /std::max(block.maximumChargeGaussScale(),1.0);
    };
    const double coarseGaussResidual=gaussResidualAt(chargeCloudRestRadius);
    const double fineGaussResidual=gaussResidualAt(0.5*chargeCloudRestRadius);
    MaxwellBlock boundCurrentTest(24,0.5*chargeCloudRestRadius,{});
    const Vec3 boundVelocity{1.2e5,-0.8e5,0.4e5};
    const double boundDt=0.1*boundCurrentTest.courantTimeStep();
    boundCurrentTest.clearSources();
    boundCurrentTest.depositCovariantDipole({},boundVelocity,
                                             Vec3{0,0,bohrMagneton});
    boundCurrentTest.finalizeBoundInstantaneous();
    const std::vector<double> oldBoundCharge=boundCurrentTest.chargeDensitySnapshot();
    const std::vector<Vec3> oldPolarization=boundCurrentTest.polarizationSnapshot();
    boundCurrentTest.clearSources();
    boundCurrentTest.depositCovariantDipole(boundVelocity*boundDt,boundVelocity,
                                             Vec3{0,0,bohrMagneton});
    boundCurrentTest.finalizeBoundCurrent(oldPolarization,boundDt);
    const double boundContinuityResidual=boundCurrentTest.maximumContinuityResidual(
        oldBoundCharge,boundDt);
    const double boundContinuityScale=std::max(
        boundCurrentTest.maximumCurrentGradientScale(),1.0);
    const double relativeBoundContinuity=boundContinuityResidual/boundContinuityScale;
    State retardedInitialState;
    retardedInitialState.firstPosition={-1.5*chargeCloudRestRadius,0,0};
    retardedInitialState.secondPosition={1.5*chargeCloudRestRadius,0,0};
    retardedInitialState.firstVelocity={0,0.03*c,0};
    retardedInitialState.secondVelocity={0,-0.02*c,0};
    retardedInitialState.firstAcceleration={1.0e20,-0.5e20,0.25e20};
    retardedInitialState.secondAcceleration={-0.7e20,0.4e20,-0.2e20};
    retardedInitialState.firstDipole={0,0,firstMagneticMoment};
    retardedInitialState.secondDipole={0,0,-secondMagneticMoment};
    const StateHistory retardedInitialHistory{retardedInitialState};
    MaxwellBlock retardedInitialField(24,0.5*chargeCloudRestRadius,{});
    initializeRetardedPairFields(retardedInitialField,retardedInitialState,
                                 retardedInitialHistory,1200);
    const double retardedInitialGauss=
        retardedInitialField.maximumElectricGaussResidual()
        /std::max(retardedInitialField.maximumChargeGaussScale(),1.0);
    const double retardedInitialDivB=
        retardedInitialField.maximumMagneticDivergence()
        /std::max(retardedInitialField.maximumMagneticGradientScale(),1.0);
    const Vec3 retardedProbe{0,4.0*chargeCloudRestRadius,0};
    const ElectromagneticField acceleratedProbe=lienardWiechertField(
        retardedProbe,retardedInitialState.time,retardedInitialHistory,
        retardedInitialState,true,firstCharge,chargeCloudRestRadius);
    State uniformInitialState=retardedInitialState;
    uniformInitialState.firstAcceleration={};
    const StateHistory uniformInitialHistory{uniformInitialState};
    const ElectromagneticField uniformProbe=lienardWiechertField(
        retardedProbe,uniformInitialState.time,uniformInitialHistory,
        uniformInitialState,true,firstCharge,chargeCloudRestRadius);
    const double retardedAccelerationSignal=
        (acceleratedProbe.electric-uniformProbe.electric).norm()
        /std::max(acceleratedProbe.electric.norm(),1.0);
    const auto planeWaveError=[&](int cellsPerAxis) {
        const double physicalExtent=12.0*chargeCloudRestRadius;
        MaxwellBlock block(cellsPerAxis,physicalExtent/cellsPerAxis,{});
        const double waveNumber=2.0*pi/physicalExtent;
        block.setVacuumPlaneWaveWavenumber(1.0e5,waveNumber);
        const double targetTime=0.5*physicalExtent/c;
        const int steps=static_cast<int>(std::ceil(
            targetTime/(0.45*block.courantTimeStep())));
        const double step=targetTime/steps;
        for(int index=0;index<steps;++index) block.advance(step);
        return block.planeWaveRelativeError(1.0e5,waveNumber,targetTime);
    };
    const double coarseWaveError=planeWaveError(12);
    const double fineWaveError=planeWaveError(24);
    const double maxwellConvergenceOrder=std::log(coarseWaveError/fineWaveError)
                                        /std::log(2.0);
    State boostedChargeState;
    boostedChargeState.firstVelocity={0.35*c,0,0};
    const StateHistory boostedChargeHistory{boostedChargeState};
    const Vec3 boostProbe{2.0*chargeCloudRestRadius,
                          3.0*chargeCloudRestRadius,
                         -1.0*chargeCloudRestRadius};
    const ElectromagneticField boostedLw=lienardWiechertField(
        boostProbe,0.0,boostedChargeHistory,boostedChargeState,true,eCharge);
    const double boostGamma=gamma(boostedChargeState.firstVelocity);
    const Vec3 boostAxis=boostedChargeState.firstVelocity
                        /boostedChargeState.firstVelocity.norm();
    const Vec3 parallelPrime=boostAxis*(dot(boostProbe,boostAxis)*boostGamma);
    const Vec3 restEventPosition=parallelPrime+(boostProbe
        -boostAxis*dot(boostProbe,boostAxis));
    const Vec3 restElectric=restEventPosition
        *(coulomb*eCharge/std::pow(restEventPosition.norm(),3));
    const Vec3 expectedBoostedElectric=boostAxis*dot(restElectric,boostAxis)
        +(restElectric-boostAxis*dot(restElectric,boostAxis))*boostGamma;
    const Vec3 expectedBoostedMagnetic=cross(
        boostedChargeState.firstVelocity,expectedBoostedElectric)/(c*c);
    const double lorentzFieldResidual=std::max(
        (boostedLw.electric-expectedBoostedElectric).norm()
            /std::max(expectedBoostedElectric.norm(),1.0),
        (boostedLw.magnetic-expectedBoostedMagnetic).norm()
            /std::max(expectedBoostedMagnetic.norm(),1.0e-300));
    const double fieldInvariantResidual=std::max(
        std::abs(dot(boostedLw.electric,boostedLw.magnetic))
            /std::max(boostedLw.electric.norm()*boostedLw.magnetic.norm(),1.0),
        std::abs(boostedLw.electric.squaredNorm()
                -c*c*boostedLw.magnetic.squaredNorm()
                -restElectric.squaredNorm())
            /std::max(restElectric.squaredNorm(),1.0));
    const Vec3 productionFirst{0.5*bohrRadius,0,0};
    const Vec3 productionSecond{-0.5*bohrRadius,0,0};
    MaxwellPairPatchHierarchy productionGeometry(
        productionFirst,productionSecond);
    const bool productionInitialCoverage=productionGeometry.coversPair(
        productionFirst,productionSecond);
    productionGeometry.farField().setUniformField({1,0,0},{});
    productionGeometry.firstPatch().setUniformField({2,0,0},{});
    productionGeometry.secondPatch().setUniformField({3,0,0},{});
    const auto [firstPatchField,unusedFirstMagnetic]=
        productionGeometry.fieldAt(productionFirst);
    const auto [secondPatchField,unusedSecondMagnetic]=
        productionGeometry.fieldAt(productionSecond);
    const Vec3 shiftedProductionFirst=productionFirst
        +Vec3{0,0.25*productionGeometry.firstPatch().extent(),0};
    const std::size_t productionMovedPatches=productionGeometry.follow(
        shiftedProductionFirst,productionSecond);
    const bool productionMovedCoverage=productionGeometry.coversPair(
        shiftedProductionFirst,productionSecond);
    MaxwellPairPatchHierarchy branchCouplingGeometry(
        productionFirst,productionSecond);
    const double branchWaveNumber=2.0*pi
        /branchCouplingGeometry.farField().extent();
    branchCouplingGeometry.farField().setVacuumPlaneWaveWavenumber(
        1.0e5,branchWaveNumber);
    branchCouplingGeometry.firstPatch().setVacuumPlaneWaveWavenumber(
        1.0e5,branchWaveNumber);
    branchCouplingGeometry.secondPatch().setVacuumPlaneWaveWavenumber(
        1.0e5,branchWaveNumber);
    MaxwellPairPatchHierarchy serialBranchGeometry=branchCouplingGeometry;
    const double branchInitialEnergy=
        branchCouplingGeometry.compositeVolumeIntegrals().energy;
    const auto parallelBranchStart=std::chrono::steady_clock::now();
    const MaxwellPairPatchHierarchy::CouplingStep branchCouplingStep=
        branchCouplingGeometry.advanceSubcycled(
            0.2*branchCouplingGeometry.farField().courantTimeStep());
    const double parallelBranchSeconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-parallelBranchStart).count();
    const auto serialBranchStart=std::chrono::steady_clock::now();
    const MaxwellPairPatchHierarchy::CouplingStep serialBranchCouplingStep=
        serialBranchGeometry.advanceSubcycled(
            0.2*serialBranchGeometry.farField().courantTimeStep(),false);
    const double serialBranchSeconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-serialBranchStart).count();
    const double branchRelativeEnergyDefect=std::abs(
        branchCouplingStep.energyDefect)/std::max(branchInitialEnergy,1.0e-300);
    const double branchMomentumScale=branchInitialEnergy/c;
    const double branchRelativeMomentumDefect=
        branchCouplingStep.momentumDefect.norm()
        /std::max(branchMomentumScale,1.0e-300);
    const auto blockDifference=[](const MaxwellBlock& first,
                                  const MaxwellBlock& second) {
        double maximum=0.0;
        for(std::size_t p=0;p<first.cells().size();++p) {
            maximum=std::max({maximum,
                (first.cells()[p].electric-second.cells()[p].electric).norm(),
                c*(first.cells()[p].magnetic-second.cells()[p].magnetic).norm()});
        }
        return maximum;
    };
    const double branchParallelResidual=std::max({
        blockDifference(branchCouplingGeometry.farField(),
                        serialBranchGeometry.farField()),
        blockDifference(branchCouplingGeometry.firstPatch(),
                        serialBranchGeometry.firstPatch()),
        blockDifference(branchCouplingGeometry.secondPatch(),
                        serialBranchGeometry.secondPatch())})/1.0e5;
    const double branchParallelSpeedup=serialBranchSeconds
        /std::max(parallelBranchSeconds,1.0e-12);
    YeeMaxwellBlock yeeTest(24,0.5*chargeCloudRestRadius,{});
    const Vec3 yeeBefore{-0.37*chargeCloudRestRadius,
                          0.21*chargeCloudRestRadius,
                         -0.16*chargeCloudRestRadius};
    const Vec3 yeeAfter=yeeBefore+Vec3{0.43,-0.31,0.27}
        *yeeTest.courantTimeStep()*c;
    const double yeeDt=0.4*yeeTest.courantTimeStep();
    // Keep the displacement subluminal for the selected dt.
    const Vec3 yeeEnd=yeeBefore+(yeeAfter-yeeBefore)
        *(yeeDt/yeeTest.courantTimeStep());
    yeeTest.depositInitialCic(eCharge,yeeBefore);
    yeeTest.clearSources();
    yeeTest.depositEsirkepovCic(eCharge,yeeBefore,yeeEnd,yeeDt);
    const double yeeContinuity=yeeTest.maximumContinuityResidual(yeeDt)
        /yeeTest.continuityScale(yeeDt);
    YeeMaxwellBlock yeeExtendedSources(28,0.5*chargeCloudRestRadius,{});
    const Vec3 cloudVelocityBefore{0.11*c,-0.07*c,0.04*c};
    const Vec3 cloudVelocityAfter{0.16*c,-0.03*c,0.08*c};
    const Vec3 cloudCentreBefore{-0.31*chargeCloudRestRadius,
                                  0.18*chargeCloudRestRadius,
                                 -0.09*chargeCloudRestRadius};
    const double cloudDt=0.32*yeeExtendedSources.courantTimeStep();
    const Vec3 cloudCentreAfter=cloudCentreBefore
        +(cloudVelocityBefore+cloudVelocityAfter)*(0.5*cloudDt);
    const Vec3 yeeRestDipole{0.21*bohrMagneton,-0.17*bohrMagneton,
                             0.93*bohrMagneton};
    yeeExtendedSources.depositInitialGaussian(
        eCharge,chargeCloudRestRadius,cloudCentreBefore,cloudVelocityBefore);
    yeeExtendedSources.depositCovariantDipoleYee(
        cloudCentreBefore,cloudVelocityBefore,yeeRestDipole,
        chargeCloudRestRadius);
    yeeExtendedSources.finalizeInitialBoundSources();
    yeeExtendedSources.clearSources();
    yeeExtendedSources.depositGaussianEsirkepov(
        eCharge,chargeCloudRestRadius,cloudCentreBefore,cloudVelocityBefore,
        cloudCentreAfter,cloudVelocityAfter,cloudDt);
    yeeExtendedSources.depositCovariantDipoleYee(
        cloudCentreAfter,cloudVelocityAfter,yeeRestDipole,
        chargeCloudRestRadius);
    yeeExtendedSources.finalizeBoundSources(cloudDt);
    const double yeeExtendedContinuity=
        yeeExtendedSources.maximumContinuityResidual(cloudDt)
        /yeeExtendedSources.continuityScale(cloudDt);
    const double yeeGaussianCharge=yeeExtendedSources.totalCharge()/eCharge;
    yeeTest.setPlaneWave(1.0e5);
    const double yeeDivBBefore=yeeTest.maximumMagneticDivergence()
        /yeeTest.magneticDivergenceScale();
    for(int step=0;step<8;++step) yeeTest.advance(0.2*yeeTest.courantTimeStep());
    const double yeeDivBAfter=yeeTest.maximumMagneticDivergence()
        /yeeTest.magneticDivergenceScale();
    YeeAmrHierarchy yeeAmrTest(12,chargeCloudRestRadius,{});
    yeeAmrTest.setPlaneWave(1.0e5);
    for(int step=0;step<6;++step)
        yeeAmrTest.advanceSubcycled(0.28*yeeAmrTest.coarse().courantTimeStep());
    const double yeeAmrRestriction=yeeAmrTest.coarse()
        .maximumRestrictionResidual(yeeAmrTest.fine())/1.0e5;
    const double yeeAmrDivB=yeeAmrTest.coarse().maximumMagneticDivergence()
        /yeeAmrTest.coarse().magneticDivergenceScale();
    YeeAmrHierarchy yeeAmrSourceTest(12,chargeCloudRestRadius,{});
    const Vec3 amrSourceBefore{-0.23*chargeCloudRestRadius,
                                0.14*chargeCloudRestRadius,
                                0.19*chargeCloudRestRadius};
    const double amrSourceDt=0.25*yeeAmrSourceTest.coarse().courantTimeStep();
    const Vec3 amrSourceAfter=amrSourceBefore
        +Vec3{0.19*c,-0.08*c,0.06*c}*amrSourceDt;
    yeeAmrSourceTest.fine().depositInitialGaussian(
        eCharge,chargeCloudRestRadius,amrSourceBefore,{0.19*c,-0.08*c,0.06*c});
    yeeAmrSourceTest.reflux();
    yeeAmrSourceTest.fine().clearSources();
    yeeAmrSourceTest.fine().depositGaussianEsirkepov(
        eCharge,chargeCloudRestRadius,amrSourceBefore,{0.19*c,-0.08*c,0.06*c},
        amrSourceAfter,{0.19*c,-0.08*c,0.06*c},amrSourceDt);
    const std::size_t yeeAmrRefluxed=yeeAmrSourceTest.reflux();
    const double yeeAmrContinuity=
        yeeAmrSourceTest.coarse().maximumContinuityResidual(amrSourceDt)
        /yeeAmrSourceTest.coarse().continuityScale(amrSourceDt);
    YeeMaxwellBlock yeeCpmlTest(48,0.5*chargeCloudRestRadius,{});
    YeeMaxwellBlock yeePeriodicReference(48,0.5*chargeCloudRestRadius,{});
    yeeCpmlTest.setPlaneWavePacketX(1.0e5,1.2*chargeCloudRestRadius);
    yeePeriodicReference.setPlaneWavePacketX(1.0e5,1.2*chargeCloudRestRadius);
    yeeCpmlTest.enableConvolutionalPml(8,1.0e-10);
    const double yeeCpmlInitialEnergy=yeeCpmlTest.fieldEnergy();
    const double yeeCpmlDt=0.35*yeeCpmlTest.courantTimeStep();
    for(int step=0;step<180;++step) {
        yeeCpmlTest.advance(yeeCpmlDt);
        yeePeriodicReference.advance(yeeCpmlDt);
    }
    const double yeeCpmlInteriorFraction=yeeCpmlTest.fieldEnergy(8)
        /yeeCpmlInitialEnergy;
    const double yeePeriodicInteriorFraction=yeePeriodicReference.fieldEnergy(8)
        /yeeCpmlInitialEnergy;
    const double yeeCpmlDivB=yeeCpmlTest.maximumMagneticDivergence(8)
        /yeeCpmlTest.magneticDivergenceScale(8);
    const Vec3 reversibleMomentum{0.17*firstMass*c,-0.09*firstMass*c,
                                   0.06*firstMass*c};
    const Vec3 reversibleElectric{2.0e5,-1.0e5,0.7e5};
    const Vec3 reversibleMagnetic{0.13,-0.08,0.05};
    const double reversibleDt=2.0e-20;
    const Vec3 pushedMomentum=relativisticBorisPush(reversibleMomentum,
        firstCharge,firstMass,reversibleElectric,reversibleMagnetic,reversibleDt);
    const Vec3 recoveredMomentum=relativisticBorisPush(pushedMomentum,
        firstCharge,firstMass,reversibleElectric,reversibleMagnetic,-reversibleDt);
    const double yeePusherReversibility=(recoveredMomentum-reversibleMomentum).norm()
        /reversibleMomentum.norm();
    State conservativeReverseStart=particles;
    conservativeReverseStart.firstDipole={};
    conservativeReverseStart.secondDipole={};
    conservativeReverseStart.firstProperDipole={};
    conservativeReverseStart.secondProperDipole={};
    State conservativeReverse=conservativeReverseStart;
    constexpr double conservativeReverseDt=2.0e-22;
    integrateConservativeMidpoint(conservativeReverse,conservativeReverseDt,12);
    integrateConservativeMidpoint(conservativeReverse,-conservativeReverseDt,12);
    const double conservativeReverseResidual=std::max({
        (conservativeReverse.firstPosition
            -conservativeReverseStart.firstPosition).norm()/bohrRadius,
        (conservativeReverse.secondPosition
            -conservativeReverseStart.secondPosition).norm()/bohrRadius,
        (conservativeReverse.firstVelocity
            -conservativeReverseStart.firstVelocity).norm()/c,
        (conservativeReverse.secondVelocity
            -conservativeReverseStart.secondVelocity).norm()/c});
    YeeMaxwellBlock yeeCoupledField(32,0.5*chargeCloudRestRadius,{});
    State yeeCoupledState;
    yeeCoupledState.firstPosition={-2.2*chargeCloudRestRadius,0,0};
    yeeCoupledState.secondPosition={2.2*chargeCloudRestRadius,0,0};
    yeeCoupledState.firstVelocity={0,0.025*c,0};
    yeeCoupledState.secondVelocity={0,-0.025*c,0};
    yeeCoupledState.firstDipole={0,0,firstMagneticMoment};
    yeeCoupledState.secondDipole={0,0,-secondMagneticMoment};
    // The deposited charge has to be the charge the pusher then moves.  These
    // were written as -e and +e, which matches the default pair and inverts
    // both signs for, say, proton+electron: the current deposited by Esirkepov
    // then contradicted the motion and the continuity residual jumped to 1.0
    // against a 1e-11 tolerance.  The tell was that proton+electron failed
    // while antiproton+positron -- the same system under charge conjugation --
    // passed, because there the hard-coded signs happened to line up.
    yeeCoupledField.depositInitialGaussian(firstCharge,chargeCloudRestRadius,
        yeeCoupledState.firstPosition,yeeCoupledState.firstVelocity);
    yeeCoupledField.depositInitialGaussian(secondCharge,chargeCloudRestRadius,
        yeeCoupledState.secondPosition,yeeCoupledState.secondVelocity);
    yeeCoupledField.depositCovariantDipoleYee(yeeCoupledState.firstPosition,
        yeeCoupledState.firstVelocity,yeeCoupledState.firstDipole,
        chargeCloudRestRadius);
    yeeCoupledField.depositCovariantDipoleYee(yeeCoupledState.secondPosition,
        yeeCoupledState.secondVelocity,yeeCoupledState.secondDipole,
        chargeCloudRestRadius);
    yeeCoupledField.finalizeInitialBoundSources();
    const double yeeCoupledDt=0.18*yeeCoupledField.courantTimeStep();
    double yeeCoupledContinuity=0.0;
    for(int step=0;step<4;++step) {
        pushStateWithYeeField(yeeCoupledState,yeeCoupledField,yeeCoupledDt);
        yeeCoupledContinuity=std::max(yeeCoupledContinuity,
            yeeCoupledField.maximumContinuityResidual(yeeCoupledDt)
            /yeeCoupledField.continuityScale(yeeCoupledDt));
    }
    const double yeeCoupledCharge=yeeCoupledField.totalCharge()/eCharge;
    const double yeeCoupledBeta=std::max(yeeCoupledState.firstVelocity.norm(),
        yeeCoupledState.secondVelocity.norm())/c;
    State sharedEngineVisualState=yeeCoupledState;
    State sharedEngineStatisticalState=yeeCoupledState;
    ClassicalTrajectoryEngine visualEngine(sharedEngineVisualState);
    ClassicalTrajectoryEngine statisticalEngine(sharedEngineStatisticalState);
    const double sharedEngineDt=1.0e-20;
    bool sharedEngineAdvanced=true;
    for(int step=0;step<6;++step)
        sharedEngineAdvanced=visualEngine.advance(
            sharedEngineVisualState,sharedEngineDt)
            &&statisticalEngine.advance(
                sharedEngineStatisticalState,sharedEngineDt)
            &&sharedEngineAdvanced;
    const double sharedEngineResidual=std::max({
        (sharedEngineVisualState.firstPosition
            -sharedEngineStatisticalState.firstPosition).norm()
            /chargeCloudRestRadius,
        (sharedEngineVisualState.secondPosition
            -sharedEngineStatisticalState.secondPosition).norm()
            /chargeCloudRestRadius,
        (sharedEngineVisualState.firstVelocity
            -sharedEngineStatisticalState.firstVelocity).norm()/c,
        (sharedEngineVisualState.secondVelocity
            -sharedEngineStatisticalState.secondVelocity).norm()/c});
    const double sharedInitialTotalEnergy=conservativeParticleEnergy(yeeCoupledState)
        +yeeCoupledState.radiatedEnergy+yeeCoupledState.boundFieldEnergy;
    const double sharedFinalTotalEnergy=
        conservativeParticleEnergy(sharedEngineVisualState)
        +sharedEngineVisualState.radiatedEnergy
        +sharedEngineVisualState.boundFieldEnergy;
    const double sharedEnergyBalanceResidual=std::abs(
        sharedFinalTotalEnergy-sharedInitialTotalEnergy)
        /std::max(std::abs(sharedInitialTotalEnergy),1.0e-300);
    const Vec3 sharedInitialTotalMomentum=noetherMomentum(yeeCoupledState)
        +yeeCoupledState.radiatedMomentum+yeeCoupledState.boundFieldMomentum;
    const Vec3 sharedFinalTotalMomentum=noetherMomentum(sharedEngineVisualState)
        +sharedEngineVisualState.radiatedMomentum
        +sharedEngineVisualState.boundFieldMomentum;
    const double sharedMomentumBalanceResidual=
        (sharedFinalTotalMomentum-sharedInitialTotalMomentum).norm()
        /std::max(canonicalMomentumScale(yeeCoupledState),1.0e-300);
    const Vec3 sharedInitialTotalAngular=noetherAngularMomentum(yeeCoupledState)
        +yeeCoupledState.radiatedAngularMomentum
        +yeeCoupledState.boundFieldAngularMomentum;
    const Vec3 sharedFinalTotalAngular=noetherAngularMomentum(
        sharedEngineVisualState)+sharedEngineVisualState.radiatedAngularMomentum
        +sharedEngineVisualState.boundFieldAngularMomentum;
    const double sharedAngularBalanceResidual=
        (sharedFinalTotalAngular-sharedInitialTotalAngular).norm()
        /std::max(sharedInitialTotalAngular.norm(),hbar);
    // Independent closure diagnostic: unlike pair-field dE/dP/dJ it excludes
    // boundField*, which is the explicitly reconstructed missing reservoir.
    const double sharedRawEnergyResidual=std::abs(
        conservativeParticleEnergy(sharedEngineVisualState)
        +sharedEngineVisualState.radiatedEnergy
        -conservativeParticleEnergy(yeeCoupledState)
        -yeeCoupledState.radiatedEnergy)
        /std::max(std::abs(sharedInitialTotalEnergy),1.0e-300);
    const double sharedRawMomentumResidual=(
        noetherMomentum(sharedEngineVisualState)
        +sharedEngineVisualState.radiatedMomentum
        -noetherMomentum(yeeCoupledState)-yeeCoupledState.radiatedMomentum).norm()
        /std::max(canonicalMomentumScale(yeeCoupledState),1.0e-300);
    const double sharedRawAngularResidual=(
        noetherAngularMomentum(sharedEngineVisualState)
        +sharedEngineVisualState.radiatedAngularMomentum
        -noetherAngularMomentum(yeeCoupledState)
        -yeeCoupledState.radiatedAngularMomentum).norm()
        /std::max(sharedInitialTotalAngular.norm(),hbar);

    // A discrete stochastic photon must close the full four-momentum, not
    // merely remove the requested energy from relative motion.  Probe a
    // drifting, unequal-direction state so neither zero-COM symmetry nor an
    // axis-aligned kick can make the check pass accidentally.
    State photonRecoilState=yeeCoupledState;
    photonRecoilState.firstVelocity+=Vec3{0.031*c,-0.007*c,0.004*c};
    photonRecoilState.secondVelocity+=Vec3{0.031*c,-0.007*c,0.004*c};
    const auto photonFirstBefore=two_body::fourMomentumFromVelocity(
        photonRecoilState.firstVelocity,firstMass);
    const auto photonSecondBefore=two_body::fourMomentumFromVelocity(
        photonRecoilState.secondVelocity,secondMass);
    const double photonIncomingEnergy=
        photonFirstBefore.energy+photonSecondBefore.energy;
    const double photonIncomingKinetic=
        kineticEnergy(photonRecoilState.firstVelocity,firstMass)
        +kineticEnergy(photonRecoilState.secondVelocity,secondMass);
    const Vec3 photonIncomingMomentum=
        photonFirstBefore.momentum+photonSecondBefore.momentum;
    const Vec3 photonIncomingComVelocity=photonIncomingMomentum
        *(c*c/photonIncomingEnergy);
    Vec3 photonProbeDirection{0.31,-0.47,0.826};
    photonProbeDirection=photonProbeDirection
        *(1.0/photonProbeDirection.norm());
    // Keep the recoil finite in units of the active pair instead of asking
    // double precision to recover the same fixed 0.25 eV kick from both an
    // electron and a proton rest-energy momentum.  One pair binding energy is
    // still small beside the prepared state's O(0.025c) mechanical energy,
    // remains equally meaningful for every selectable mass ratio, and keeps
    // the mixed-pair recoil resolved below the strict 1e-8 closure bound
    // without relaxing that bound.
    const double photonProbeEnergy=
        pairBindingEnergy(activePair);
    const StochasticPhotonRecoil photonRecoil=applyStochasticDipolePhoton(
        photonRecoilState,photonProbeEnergy,photonProbeDirection);
    const auto photonFirstAfter=two_body::fourMomentumFromVelocity(
        photonRecoilState.firstVelocity,firstMass);
    const auto photonSecondAfter=two_body::fourMomentumFromVelocity(
        photonRecoilState.secondVelocity,secondMass);
    const double photonFinalKinetic=
        kineticEnergy(photonRecoilState.firstVelocity,firstMass)
        +kineticEnergy(photonRecoilState.secondVelocity,secondMass);
    const auto photonLab=two_body::boostFourMomentum(
        {photonProbeEnergy,photonProbeDirection*(photonProbeEnergy/c)},
        photonIncomingComVelocity);
    const double photonFourEnergyResidual=photonRecoil.emitted
        // Rest energies are identical before and after and cancel exactly.
        // Subtracting total energies first discarded the eV-scale recoil in
        // the GeV-scale proton rest energy; the stable kinetic form tests the
        // same four-energy identity without catastrophic cancellation.
        ?std::abs(photonFinalKinetic+photonLab.energy
            -photonIncomingKinetic)/photonProbeEnergy
        :std::numeric_limits<double>::infinity();
    const double photonFourMomentumResidual=photonRecoil.emitted
        ?(photonFirstAfter.momentum+photonSecondAfter.momentum
            +photonLab.momentum-photonIncomingMomentum).norm()
            /(photonProbeEnergy/c)
        :std::numeric_limits<double>::infinity();
    State forbiddenPhotonState=yeeCoupledState;
    const State forbiddenPhotonBefore=forbiddenPhotonState;
    const StochasticPhotonRecoil forbiddenPhoton=applyStochasticDipolePhoton(
        forbiddenPhotonState,1.0e9*eCharge,photonProbeDirection);
    const double forbiddenPhotonMutation=std::max(
        (forbiddenPhotonState.firstVelocity
            -forbiddenPhotonBefore.firstVelocity).norm()/c,
        (forbiddenPhotonState.secondVelocity
            -forbiddenPhotonBefore.secondVelocity).norm()/c);
    const double reactionMismatchFraction=std::abs(
        sharedEngineVisualState.reactionEnergyMismatch)
        /std::max(std::abs(sharedEngineVisualState.radiatedEnergy),1.0e-300);
    const MutualForces sharedFinalExternal=retardedExternalForces(
        sharedEngineVisualState,visualEngine.history());
    // Ask for the automatic model explicitly: the blending gates reported
    // below (smooth/kR/w) are only evaluated by that model, because under the
    // default individual Landau-Lifshitz model they are unused and skipped.
    const ParticleMultipoleRadiation sharedFinalRadiation=
        particleMultipoleRadiation(sharedEngineVisualState,sharedFinalExternal,
            visualEngine.history(),false,
            ChargeRadiationReactionModel::automatic);
    const double llValidity=sharedFinalRadiation.landauLifshitzValidity;
    // The three inputs of the blending gate.  They were printed as
    // "auto smooth/kR/w" and checked nowhere.
    const double blendSmoothness=
        sharedFinalRadiation.coherentDerivativeConsistency;
    const double blendCompactness=sharedFinalRadiation.sourceCompactness;
    const double blendWeight=sharedFinalRadiation.coherentWeight;

    struct ReactionModelBenchmark {
        State finalState;
        bool advanced=false;
        double rawEnergyResidual=0.0;
        double refinedEnergyResidual=0.0;
        double mechanicalEnergyChange=0.0;
        double radiatedEnergyFraction=0.0;
        double reactionFluxResidual=0.0;
        // Same mismatch, normalized by the orbit's own mechanical energy
        // instead of by the radiated energy.  The flux-normalized form is
        // unusable across pairs: over the benchmark's 8e-22 s a proton pair
        // radiates 1.2e-29 J, so any mismatch divided by it explodes -- it
        // reaches 1893 for p+pbar against 0.026 for e+e-, which reads like a
        // catastrophic defect and is nothing but the denominator vanishing.
        // Against the mechanical energy the same quantity is 9.2e-07 for e+e-,
        // 5.6e-11 for mu+mu- and 7.1e-13 for p+pbar, i.e. tiny and DECREASING
        // with mass, which is the opposite conclusion.
        double reactionMechanicalMismatch=0.0;
        double stepConvergence=0.0;
        double wallSeconds=0.0;
    };
    const auto reactionBenchmark=[&](ChargeRadiationReactionModel model) {
        constexpr double proposedStep=2.0e-22;
        constexpr int fullSteps=4;
        const auto evolve=[&](double dt,int steps,double& wallSeconds) {
            State value=yeeCoupledState;
            StateHistory history=causalInitialHistory(value);
            const auto start=std::chrono::steady_clock::now();
            bool advanced=true;
            for(int step=0;step<steps;++step) {
                // CREM_DARWIN_FORCES=1 integrates the matched Coulomb+Darwin
                // force instead of the retarded one, so the Darwin-vs-retarded
                // comparison is reproducible from here.  It does NOT give the
                // clean version of that comparison: switching the force
                // changes the TRAJECTORY, so the two runs are not the same
                // orbit.  Measured on this probe, Darwin sits a constant
                // ~4.05e-06 WORSE, with or without the charge-dipole energy --
                // the opposite of the matched-force expectation, and the
                // signature of two different paths rather than two
                // bookkeepings of one.  The clean comparison is the one
                // recorded at conservativeParticleEnergy: one orbit at a_Ps,
                // tolerance-controlled.
                integrateElectrodynamicStep(value,dt,history,true,model,
                    std::getenv("CREM_DARWIN_FORCES")==nullptr);
                advanced=isFinite(value)&&advanced;
                if(!advanced) break;
                appendStateHistory(history,value);
            }
            wallSeconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-start).count();
            return std::pair<State,bool>{value,advanced};
        };
        ReactionModelBenchmark result;
        auto [full,fullAdvanced]=evolve(
            proposedStep,fullSteps,result.wallSeconds);
        double halfWall=0.0;
        auto [half,halfAdvanced]=evolve(
            0.5*proposedStep,2*fullSteps,halfWall);
        result.finalState=full;
        result.advanced=fullAdvanced&&halfAdvanced;
        const double initialMechanical=conservativeParticleEnergy(yeeCoupledState);
        result.rawEnergyResidual=std::abs(
            conservativeParticleEnergy(full)+full.radiatedEnergy
            -initialMechanical-yeeCoupledState.radiatedEnergy)
            /std::max(std::abs(initialMechanical),1.0e-300);
        // The residual's two parts, kept apart because summing them hides
        // what each says.  rawEnergyResidual adds the RADIATED energy to the
        // change in mechanical energy; for the DISABLED reaction model no
        // force removes that flux from the mechanics, so including it does
        // not test a balance -- it partially cancels the mechanical drift and
        // makes the residual look smaller than the drift actually is.
        // Measured: dE_mech = -2.48e-06 against E_rad = +8.41e-07, summing to
        // the -1.64e-06 that gets reported.  The physically meaningful number
        // for the disabled model is the first one.
        result.mechanicalEnergyChange=
            (conservativeParticleEnergy(full)-initialMechanical)
            /std::max(std::abs(initialMechanical),1.0e-300);
        result.radiatedEnergyFraction=
            (full.radiatedEnergy-yeeCoupledState.radiatedEnergy)
            /std::max(std::abs(initialMechanical),1.0e-300);
        result.reactionFluxResidual=std::abs(full.reactionEnergyMismatch)
            /std::max(std::abs(full.radiatedEnergy),1.0e-300);
        result.reactionMechanicalMismatch=std::abs(full.reactionEnergyMismatch)
            /std::max(std::abs(initialMechanical),1.0e-300);
        // REFINEMENT RATIO of the energy residual, reported because its
        // absence is how a real missing interaction went unseen.
        //
        // The half-step run was already being computed, but only to check
        // that POSITIONS AND VELOCITIES converge (stepConvergence below).
        // The residual itself was never refined, so a standing 4e-06 read as
        // integrator error -- while the project's own long-horizon balance
        // section states the discriminator plainly: a physical imbalance sits
        // still, a discretization error shrinks.  Applied here it is
        // immediate: the ratio measured 0.9997 with the charge-dipole energy
        // MISSING, i.e. the residual did not shrink at all and was therefore
        // physics, not arithmetic.
        //
        // It still is.  After that term was added the residual fell from
        // 3.97e-06 to 1.64e-06 and the ratio is 0.9993 -- flat again, so what
        // remains is also a physical imbalance and not the integrator.
        // Reported rather than gated: a gate here would fail today, and
        // shipping a red build to mark an open question is worse than
        // printing the number that states it.
        result.refinedEnergyResidual=std::abs(
            conservativeParticleEnergy(half)+half.radiatedEnergy
            -initialMechanical-yeeCoupledState.radiatedEnergy)
            /std::max(std::abs(initialMechanical),1.0e-300);
        result.stepConvergence=std::max({
            (full.firstPosition-half.firstPosition).norm()/bohrRadius,
            (full.secondPosition-half.secondPosition).norm()/bohrRadius,
            (full.firstVelocity-half.firstVelocity).norm()/c,
            (full.secondVelocity-half.secondVelocity).norm()/c});
        return result;
    };
    const ReactionModelBenchmark reactionDisabled=reactionBenchmark(
        ChargeRadiationReactionModel::disabled);
    const ReactionModelBenchmark reactionLl=reactionBenchmark(
        ChargeRadiationReactionModel::individualLandauLifshitz);
    const ReactionModelBenchmark reactionCoherent=reactionBenchmark(
        ChargeRadiationReactionModel::coherentElectricDipole);
    const ReactionModelBenchmark reactionAutomatic=reactionBenchmark(
        ChargeRadiationReactionModel::automatic);
    const double coherentCostRatio=reactionCoherent.wallSeconds
        /std::max(reactionLl.wallSeconds,1.0e-300);

    // End-to-end Lorentz covariance test of the charge sector.  The boost is
    // perpendicular to the orbital plane, so both initial particle events
    // lie on the same boosted simultaneity hyperplane.  Endpoint events are
    // compared after an explicit Lorentz transformation, never at equal
    // untransformed coordinate times.
    const Vec3 covarianceBoost{0,0,0.35*c};
    const double covarianceBoostGamma=gamma(covarianceBoost);
    const auto boostEvent=[&](const Vec3& position,double time) {
        const double speedSquared=covarianceBoost.squaredNorm();
        const double projection=dot(position,covarianceBoost);
        const Vec3 transformedPosition=position
            +covarianceBoost*((covarianceBoostGamma-1.0)*projection/speedSquared
                              -covarianceBoostGamma*time);
        const double transformedTime=covarianceBoostGamma
            *(time-projection/(c*c));
        return std::pair<Vec3,double>{transformedPosition,transformedTime};
    };
    const auto boostVelocity=[&](const Vec3& velocity) {
        const double speedSquared=covarianceBoost.squaredNorm();
        const Vec3 parallel=covarianceBoost
            *(dot(velocity,covarianceBoost)/speedSquared);
        const Vec3 perpendicular=velocity-parallel;
        const double denominator=1.0-dot(velocity,covarianceBoost)/(c*c);
        return (parallel-covarianceBoost+perpendicular/covarianceBoostGamma)
            /denominator;
    };
    const auto boostFourVector=[&](const FourVector& value) {
        const double speedSquared=covarianceBoost.squaredNorm();
        const Vec3 parallel=covarianceBoost
            *(dot(value.space,covarianceBoost)/speedSquared);
        const Vec3 perpendicular=value.space-parallel;
        return FourVector{
            covarianceBoostGamma*(value.time
                -dot(covarianceBoost,value.space)/c),
            perpendicular+(parallel-covarianceBoost*(value.time/c))
                *covarianceBoostGamma};
    };
    State covarianceRest;
    covarianceRest.firstPosition={0,0.5*bohrRadius,0};
    covarianceRest.secondPosition={0,-0.5*bohrRadius,0};
    const double covarianceOrbitalSpeed=std::sqrt(
        pairCoulombStrength/(2.0*firstMass*bohrRadius));
    covarianceRest.firstVelocity={covarianceOrbitalSpeed,0,0};
    covarianceRest.secondVelocity={-covarianceOrbitalSpeed,0,0};
    covarianceRest.firstDipole={0,0,1.0e-3*firstMagneticMoment};
    covarianceRest.secondDipole={0,0,-1.0e-3*secondMagneticMoment};
    covarianceRest.firstProperDipole=covarianceRest.firstDipole;
    covarianceRest.secondProperDipole=covarianceRest.secondDipole;
    // Make the rest state exactly self-consistent BEFORE boosting anything.
    // It used to be initialized with firstProperDipole = firstDipole and a
    // zero electric dipole, which is only true for a particle at rest -- these
    // orbit at beta = 5.2e-3, so the motional electric dipole (v x m)/c^2 was
    // missing from the tensor that then got boosted.  The evolution creates it
    // one step later, so the initial state disagreed with the code's own
    // representation by exactly beta_orb.
    synchronizeCovariantDipoles(covarianceRest);
    State covarianceMoving=covarianceRest;
    const auto [boostedFirstPosition,boostedFirstTime]=boostEvent(
        covarianceRest.firstPosition,covarianceRest.time);
    const auto [boostedSecondPosition,boostedSecondTime]=boostEvent(
        covarianceRest.secondPosition,covarianceRest.time);
    covarianceMoving.firstPosition=boostedFirstPosition;
    covarianceMoving.secondPosition=boostedSecondPosition;
    covarianceMoving.firstVelocity=boostVelocity(
        covarianceRest.firstVelocity);
    covarianceMoving.secondVelocity=boostVelocity(
        covarianceRest.secondVelocity);
    // The FULL tensor is boosted, electric part included.  Boosting {0, m}
    // discarded the motional electric dipole the rest frame already carries.
    const DipoleTensor boostedFirstDipole=lorentzBoostDipole(
        {covarianceRest.firstElectricDipole,covarianceRest.firstDipole},
        covarianceBoost);
    const DipoleTensor boostedSecondDipole=lorentzBoostDipole(
        {covarianceRest.secondElectricDipole,covarianceRest.secondDipole},
        covarianceBoost);
    covarianceMoving.firstDipole=boostedFirstDipole.magnetic;
    covarianceMoving.secondDipole=boostedSecondDipole.magnetic;
    covarianceMoving.firstElectricDipole=boostedFirstDipole.electric;
    covarianceMoving.secondElectricDipole=boostedSecondDipole.electric;
    // Derived through the spin FOUR-VECTOR, which transforms simply, instead
    // of boosting the lab tensor back and keeping only its magnetic part.
    // That discard was the whole problem: what it threw away is exactly the
    // Wigner rotation that composing the frame boost with the particle's own
    // orbital velocity produces, so the moving proper dipole came out wrong at
    // first order in beta_boost * beta_orb.  The four-vector route carries the
    // rotation for free and is exact.
    covarianceMoving.firstProperDipole=properDipoleFromFourVector(
        boostFourVector(dipoleFourVector(covarianceRest.firstProperDipole,
            covarianceRest.firstVelocity)),covarianceMoving.firstVelocity);
    covarianceMoving.secondProperDipole=properDipoleFromFourVector(
        boostFourVector(dipoleFourVector(covarianceRest.secondProperDipole,
            covarianceRest.secondVelocity)),covarianceMoving.secondVelocity);
    covarianceMoving.time=0.5*(boostedFirstTime+boostedSecondTime);
    // What the exact boost says the lab tensor should be, kept before the
    // state is forced back onto the manifold the model can actually represent.
    const Vec3 exactBoostedFirstMagnetic=covarianceMoving.firstDipole;
    // The model stores a PROPER dipole and rebuilds the lab tensor as
    // boost({0, m_proper}, v).  That parametrization cannot express every
    // tensor: composing the frame boost with the particle's own orbital
    // velocity is a boost TIMES A WIGNER ROTATION, and the rotated tensor is
    // not of the form boost({0, m}, v) for any m.  Re-synchronizing here puts
    // the moving state exactly on that manifold, so both states are now
    // self-consistent to machine precision and the residual below measures the
    // gap itself instead of hiding it inside an inconsistent initial state.
    synchronizeCovariantDipoles(covarianceMoving);
    const double covarianceRepresentabilityGap=
        (covarianceMoving.firstDipole-exactBoostedFirstMagnetic).norm()
        /std::max(exactBoostedFirstMagnetic.norm(),1.0e-300);
    ClassicalTrajectoryEngine covarianceRestEngine(covarianceRest);
    ClassicalTrajectoryEngine covarianceMovingEngine(covarianceMoving);
    constexpr double covarianceRestStep=1.0e-20;
    bool covarianceTrajectoriesAdvanced=true;
    for(int step=0;step<8;++step) {
        covarianceTrajectoriesAdvanced=covarianceRestEngine.advance(
            covarianceRest,covarianceRestStep)&&covarianceMovingEngine.advance(
            covarianceMoving,covarianceBoostGamma*covarianceRestStep)
            &&covarianceTrajectoriesAdvanced;
    }
    const auto [expectedFirstPosition,expectedFirstTime]=boostEvent(
        covarianceRest.firstPosition,covarianceRest.time);
    const auto [expectedSecondPosition,expectedSecondTime]=boostEvent(
        covarianceRest.secondPosition,covarianceRest.time);
    const Vec3 expectedFirstVelocity=boostVelocity(
        covarianceRest.firstVelocity);
    const Vec3 expectedSecondVelocity=boostVelocity(
        covarianceRest.secondVelocity);
    const double covarianceWorldlineResidual=std::max({
        (covarianceMoving.firstPosition-expectedFirstPosition).norm()/bohrRadius,
        (covarianceMoving.secondPosition-expectedSecondPosition).norm()/bohrRadius,
        std::abs(covarianceMoving.time-expectedFirstTime)
            /(covarianceRestStep*8.0),
        std::abs(covarianceMoving.time-expectedSecondTime)
            /(covarianceRestStep*8.0)});
    const double covarianceVelocityResidual=std::max(
        (covarianceMoving.firstVelocity-expectedFirstVelocity).norm()/c,
        (covarianceMoving.secondVelocity-expectedSecondVelocity).norm()/c);
    const DipoleTensor expectedFirstTensor=lorentzBoostDipole(
        {covarianceRest.firstElectricDipole,covarianceRest.firstDipole},
        covarianceBoost);
    const DipoleTensor expectedSecondTensor=lorentzBoostDipole(
        {covarianceRest.secondElectricDipole,covarianceRest.secondDipole},
        covarianceBoost);
    const double covarianceDipoleEvolutionResidual=std::max({
        (covarianceMoving.firstDipole-expectedFirstTensor.magnetic).norm()
            /std::max(expectedFirstTensor.magnetic.norm(),1.0e-300),
        (covarianceMoving.secondDipole-expectedSecondTensor.magnetic).norm()
            /std::max(expectedSecondTensor.magnetic.norm(),1.0e-300),
        (covarianceMoving.firstElectricDipole
            -expectedFirstTensor.electric).norm()
            /std::max(expectedFirstTensor.magnetic.norm()/c,1.0e-300),
        (covarianceMoving.secondElectricDipole
            -expectedSecondTensor.electric).norm()
            /std::max(expectedSecondTensor.magnetic.norm()/c,1.0e-300)});
    const MutualForces covarianceRestExternal=retardedExternalForces(
        covarianceRest,covarianceRestEngine.history());
    const MutualForces covarianceMovingExternal=retardedExternalForces(
        covarianceMoving,covarianceMovingEngine.history());
    const LocalElectromagneticFields covarianceRestLocal=
        localRelativisticFields(covarianceRest,covarianceRestEngine.history());
    const LocalElectromagneticFields covarianceMovingLocal=
        localRelativisticFields(covarianceMoving,covarianceMovingEngine.history());
    constexpr double covarianceBmtProbeDt=1.0e-22;
    // Evolve one role's dipole in the rest frame and in the boosted frame, and
    // compare the boosted rest-frame result with the moving-frame one.  The
    // covariance being tested is a property of a SINGLE particle, so one probe
    // is enough to exercise advanceCovariantBmt itself.
    //
    // Both roles are probed because the integrator's behaviour depends on the
    // arguments it is handed.  For a symmetric pair the two residuals come out
    // bit-identical -- e+e- has equal masses, equal g and opposite charges, so
    // the second probe is pure redundancy there.  For an asymmetric pair they
    // are not: measured on proton+electron they differ by a factor of 540
    // (8.1e-9 against 4.4e-6), and the SECOND one is the binding constraint.
    // Without it a defect affecting only the heavier-coupling role would go
    // unseen for every pair.
    //
    // What this does NOT catch is a role mix-up at the CALL SITES -- passing
    // one particle's g-factor or charge-to-mass ratio to the other inside
    // applyDipolePrecession().  The check builds its own arguments explicitly,
    // so it never executes those call sites; injecting exactly that swap
    // leaves the residual unchanged to every printed digit, for e+e- and for
    // proton+electron alike.  Guarding call-site role assignment needs a test
    // that drives applyDipolePrecession() itself and compares each role's
    // precession against its own g, which does not exist yet.
    const auto bmtCovarianceResidual=[&](const Vec3& restProperDipole,
                                         const Vec3& restVelocity,
                                         const ElectromagneticField& restField,
                                         const Vec3& movingProperDipole,
                                         const Vec3& movingVelocity,
                                         const ElectromagneticField& movingField,
                                         double chargeToMass,double gFactor) {
        const Vec3 restProbe=advanceCovariantBmt(restProperDipole,restVelocity,
            restField,chargeToMass,covarianceBmtProbeDt,gFactor);
        const Vec3 movingProbe=advanceCovariantBmt(movingProperDipole,
            movingVelocity,movingField,chargeToMass,
            covarianceBoostGamma*covarianceBmtProbeDt,gFactor);
        const FourVector expected=boostFourVector(
            dipoleFourVector(restProbe,restVelocity));
        const FourVector actual=dipoleFourVector(movingProbe,movingVelocity);
        const double scale=std::max(restProbe.norm(),1.0e-300);
        return std::max(std::abs(actual.time-expected.time)/scale,
                        (actual.space-expected.space).norm()/scale);
    };
    const double covarianceBmtResidual=std::max(
        bmtCovarianceResidual(
            covarianceRest.firstProperDipole,covarianceRest.firstVelocity,
            covarianceRestLocal.atFirst,
            covarianceMoving.firstProperDipole,covarianceMoving.firstVelocity,
            covarianceMovingLocal.atFirst,
            firstCharge/firstMass,firstGFactor),
        bmtCovarianceResidual(
            covarianceRest.secondProperDipole,covarianceRest.secondVelocity,
            covarianceRestLocal.atSecond,
            covarianceMoving.secondProperDipole,covarianceMoving.secondVelocity,
            covarianceMovingLocal.atSecond,
            secondCharge/secondMass,secondGFactor));

    // ---------------------------------------------------------------------
    // Two probes on the PRODUCTION dipole-precession path
    // (advanceThomasBmtDipole -> properDipolePrecessionRate ->
    // thomasBmtEffectiveField).  They answer different questions, and the
    // history of this spot is the reason both are needed.
    //
    // History.  Production first used advanceCovariantBmt, a four-vector RK4
    // route.  A probe here measured a real gap against Jackson's 3D lab-time
    // formula: exactly zero at g=2, growing with anomaly and beta, up to 45%
    // for a proton g at beta=0.9.  Root cause: advanceCovariantBmt holds the
    // four-velocity u FIXED for its call (exactly how applyDipolePrecession
    // invokes it), but the covariant equation preserves its own constraint
    // a.u=0 only when u co-evolves; the final renormalization then injected
    // a correction of the same order as the precession, sourced by the
    // (g/2-1) term that vanishes at g=2.
    //
    // The replacement over-corrected.  It expanded properDipole into the LAB
    // dipole, applied Jackson's rate to that, and mapped back through an
    // inverse boost -- but Jackson 11.170 governs the REST-FRAME moment, so
    // the composition was not a rotation of properDipole and quietly failed
    // to conserve its norm at the relative level gamma-1.  The probe that
    // certified it compared d(mu_LAB)/dt against mu_lab x B_eff, which is
    // precisely the identity that implementation was built to satisfy: it
    // read machine zero by construction and could not have reported
    // otherwise.  A check whose subject and reference are the same equation
    // measures nothing.
    //
    // Hence the split below.  The first probe is a wiring check and is
    // frankly near-tautological -- it re-derives the production rate from
    // the same formula production calls, so it can only catch a mis-wired
    // argument (wrong g, wrong q/m, wrong field) or a broken RK4, not a
    // wrong choice of equation.  The second is not: orthogonality of the
    // rate to the moment is an INVARIANT of precession, independent of which
    // formula is wired in, so it stays meaningful however the sector is
    // rewritten.  It is the one gated, because it is the one that has
    // discriminating power -- it reads 4.3e-3 at beta=0.1 and 0.63 at
    // beta=0.9 against the lab-dipole route, and round-off against a genuine
    // precession.
    // Solver against specification.  advanceThomasBmtDipole now returns a
    // closed-form rotation; properDipolePrecessionRate is the equation that
    // rotation is supposed to solve.  Integrating the specification
    // numerically to a FINITE angle and comparing gives a real check with a
    // well-conditioned subtraction -- it catches a wrong rotation axis, a
    // flipped sign or a dropped factor, none of which the old
    // first-order-difference form could see (it took a step so short that
    // the closed form and the reference agreed trivially).
    const auto bmtEffectiveFieldGap=[&](double betaMagnitude,double gFactorProbe) {
        const Vec3 velocityProbe{0.0,0.0,betaMagnitude*c};
        const ElectromagneticField fieldProbe{Vec3{0.0,1.0,0.0},Vec3{1.0,0.0,0.0}};
        const Vec3 properDipoleProbe{0.3,0.5,0.81};
        constexpr double chargeToMassProbe=1.0;
        // Turn through a finite angle, so the comparison is not dominated by
        // round-off in the difference.
        const Vec3 effectiveField=thomasBmtEffectiveField(
            velocityProbe,fieldProbe,gFactorProbe);
        const double angularSpeed=
            (effectiveField*chargeToMassProbe).norm();
        if(!(angularSpeed>0.0)) return 0.0;
        const double totalDt=0.1/angularSpeed;   // 0.1 rad of precession
        const Vec3 closedForm=advanceThomasBmtDipole(properDipoleProbe,
            velocityProbe,fieldProbe,chargeToMassProbe,totalDt,gFactorProbe);
        // Reference: RK4 on the specification in many small substeps, whose
        // O((dt/N)^4) error at N=4096 is far below the quantity compared.
        constexpr int substeps=4096;
        const double substepDt=totalDt/substeps;
        Vec3 reference=properDipoleProbe;
        for(int step=0;step<substeps;++step) {
            const auto derivative=[&](const Vec3& value) {
                return properDipolePrecessionRate(value,velocityProbe,
                    fieldProbe,chargeToMassProbe,gFactorProbe);
            };
            const Vec3 k1=derivative(reference);
            const Vec3 k2=derivative(reference+k1*(0.5*substepDt));
            const Vec3 k3=derivative(reference+k2*(0.5*substepDt));
            const Vec3 k4=derivative(reference+k3*substepDt);
            reference=reference+(k1+k2*2.0+k3*2.0+k4)*(substepDt/6.0);
        }
        return (closedForm-reference).norm()
            /std::max(reference.norm(),1.0e-300);
    };
    // Fraction of the SPECIFICATION's rate that points along the moment
    // instead of across it.  Zero for any precession; this is the quantity
    // the lab-dipole route got wrong, and the one no norm-drift guard could
    // see back when advanceThomasBmtDipole still repaired the norm by
    // rescaling.  It no longer does, so the two probes are independent: this
    // one gates the equation, bmtNormConservation below gates the solver.
    const auto bmtPrecessionOrthogonality=[&](double betaMagnitude,
                                              double gFactorProbe) {
        const Vec3 velocityProbe{0.0,0.0,betaMagnitude*c};
        const ElectromagneticField fieldProbe{Vec3{0.0,1.0,0.0},Vec3{1.0,0.0,0.0}};
        const Vec3 properDipoleProbe{0.3,0.5,0.81};
        const Vec3 rate=properDipolePrecessionRate(properDipoleProbe,
            velocityProbe,fieldProbe,1.0,gFactorProbe);
        const double scale=properDipoleProbe.norm()*rate.norm();
        return scale>0.0 ? std::abs(dot(properDipoleProbe,rate))/scale : 0.0;
    };
    // Norm conservation of the PRODUCTION routine, over a large angle and
    // with no renormalization left anywhere to arrange the answer.  Large on
    // purpose: a rotation is exact at any angle, so anything that is not a
    // rotation has nowhere to hide at 2 rad.
    const auto bmtNormConservation=[&](double betaMagnitude,
                                       double gFactorProbe) {
        const Vec3 velocityProbe{0.0,0.0,betaMagnitude*c};
        const ElectromagneticField fieldProbe{Vec3{0.0,1.0,0.0},Vec3{1.0,0.0,0.0}};
        const Vec3 properDipoleProbe{0.3,0.5,0.81};
        const double angularSpeed=thomasBmtEffectiveField(
            velocityProbe,fieldProbe,gFactorProbe).norm();
        if(!(angularSpeed>0.0)) return 0.0;
        double worst=0.0;
        Vec3 value=properDipoleProbe;
        for(int step=0;step<64;++step) {
            value=advanceThomasBmtDipole(value,velocityProbe,fieldProbe,
                1.0,2.0/angularSpeed,gFactorProbe);
            worst=std::max(worst,std::abs(value.norm()
                /properDipoleProbe.norm()-1.0));
        }
        return worst;
    };
    // Two probes: the active pair's own g (beta~0 at the orbital scale is far
    // too small to be informative, so a representative orbital-ish beta is
    // used instead) and a synthetic high-anomaly/high-beta point, so these
    // stay informative regardless of which --pair the binary is built for.
    const double bmtEffectiveFieldGapActiveG=bmtEffectiveFieldGap(0.1,firstGFactor);
    const double bmtEffectiveFieldGapHighBeta=bmtEffectiveFieldGap(0.9,5.5857);
    const double bmtOrthogonalityActiveG=
        bmtPrecessionOrthogonality(0.1,firstGFactor);
    const double bmtOrthogonalityHighBeta=
        bmtPrecessionOrthogonality(0.9,5.5857);
    const double bmtNormDriftActiveG=bmtNormConservation(0.1,firstGFactor);
    const double bmtNormDriftHighBeta=bmtNormConservation(0.9,5.5857);
    // Round-off bands.  A precession is orthogonal and norm-preserving
    // exactly; the tolerances admit float noise only, and sit six orders
    // below the 4.3e-3 the lab-dipole route produced at the same beta.  The
    // solver-vs-specification gap carries the reference RK4's own truncation
    // as well, hence the looser band there.
    const bool bmtPrecessionOk=bmtOrthogonalityActiveG<1.0e-12
        && bmtOrthogonalityHighBeta<1.0e-12
        && bmtNormDriftActiveG<1.0e-12 && bmtNormDriftHighBeta<1.0e-12
        && bmtEffectiveFieldGapActiveG<1.0e-9
        && bmtEffectiveFieldGapHighBeta<1.0e-9;

    // Coupled secular spin-orbit transport.  This is separate from the BMT
    // check above: exact rotations can preserve each |mu| while still losing
    // the angular momentum transferred between spin and orbit.  Exercise a
    // finite, nonlinear turn at three substep bounds and independently rebuild
    // J=L+mu1/gamma1+mu2/gamma2 from every result.
    const double secularRadius=pairBohrRadius(activePair);
    const Vec3 secularFirstDirection{0.31,-0.47,0.826619};
    const Vec3 secularSecondDirection{-0.62,0.19,0.761249};
    // Keep (a,L) on the bound Kepler sheet.  The former vector had
    // |L|>sqrt(mu*K*a), i.e. an imaginary eccentricity which the circular
    // implementation could not notice because it ignored |L| altogether.
    const SecularSpinOrbitState secularInitial{
        {0.08*hbar,-0.05*hbar,0.78*hbar},
        secularFirstDirection
            *(firstMagneticMoment/secularFirstDirection.norm()),
        secularSecondDirection
            *(secondMagneticMoment/secularSecondDirection.norm())};
    const OrbitAveragedBmtAngularVelocities secularInitialRates=
        orbitAveragedBmtAngularVelocities(
            secularRadius,secularInitial.orbitalAngularMomentum,
            secularInitial.firstDipole,secularInitial.secondDipole,
            pairReducedMass);
    const double secularAngularSpeed=secularInitialRates.valid
        ?std::max(secularInitialRates.first.norm(),
                  secularInitialRates.second.norm()):0.0;
    const double secularElapsed=secularAngularSpeed>0.0
        ?0.2/secularAngularSpeed:0.0;
    const SecularSpinOrbitAdvance secularCoarse=
        advanceCoupledSecularSpinOrbit(
            secularInitial,secularRadius,pairReducedMass,secularElapsed,0.05);
    const SecularSpinOrbitAdvance secularFine=
        advanceCoupledSecularSpinOrbit(
            secularInitial,secularRadius,pairReducedMass,secularElapsed,0.025);
    const SecularSpinOrbitAdvance secularReference=
        advanceCoupledSecularSpinOrbit(
            secularInitial,secularRadius,pairReducedMass,secularElapsed,0.0125);
    // External torque belongs to the source of the configured field, not to
    // the orbital reaction.  Exercise that branch explicitly and restore the
    // process-wide option before any later validation probe observes it.
    const Vec3 savedExternalMagneticField=gExternalMagneticField;
    const Vec3 secularExternalDirection{0.23,-0.17,0.31};
    // A fixed field made this branch effectively zero for protonium: its
    // internal precession set the elapsed time while the external turn fell
    // below 1e-8 rad.  Match the external angular rate to the pair's own
    // internal rate, keeping the test dimensionless across masses and g.
    const double secularExternalFieldScale=secularAngularSpeed
        /std::max(std::abs(firstGyromagneticRatioOf()),
                  std::abs(secondGyromagneticRatioOf()));
    gExternalMagneticField=secularExternalDirection
        *(secularExternalFieldScale/secularExternalDirection.norm());
    const OrbitAveragedBmtAngularVelocities secularExternalRates=
        orbitAveragedBmtAngularVelocities(
            secularRadius,secularInitial.orbitalAngularMomentum,
            secularInitial.firstDipole,secularInitial.secondDipole,
            pairReducedMass);
    const double secularExternalSpeed=secularExternalRates.valid
        ?std::max(secularExternalRates.first.norm(),
                  secularExternalRates.second.norm()):0.0;
    const SecularSpinOrbitAdvance secularExternal=
        advanceCoupledSecularSpinOrbit(
            secularInitial,secularRadius,pairReducedMass,
            secularExternalSpeed>0.0?0.2/secularExternalSpeed:0.0,0.025);
    gExternalMagneticField=savedExternalMagneticField;

    // Independent eccentric-orbit reference.  Production obtains the
    // geometry analytically from eccentric anomaly.  This probe instead
    // resolves Newton's dimensionless Kepler equations with RK4 at uniform
    // laboratory-time intervals and averages the same instantaneous BMT
    // observable along that mechanically evolved orbit.  It therefore catches
    // all three defects of the former circular surrogate: fixed r=a, missing
    // radial velocity and uniform geometric-phase rather than time weighting.
    constexpr double secularReferenceEccentricitySquared=0.945;
    const double secularReferenceEccentricity=
        std::sqrt(secularReferenceEccentricitySquared);
    // Keep the independent Kepler reference outside the production field
    // floor for every species.  At a_pair the pbar-p probe would otherwise
    // put periapsis at 1.6 fm, deep below its 193 fm floor, and compare an
    // unregularized Newton orbit with deliberately frozen electrodynamics.
    const double secularEccentricReferenceAxis=std::max(
        secularRadius,10.0*separationFloor()
            /(1.0-secularReferenceEccentricity));
    const Vec3 secularReferenceNormal{0.0,0.0,1.0};
    const Vec3 secularReferenceAngularMomentum=secularReferenceNormal*std::sqrt(
        pairReducedMass*pairCoulombStrength*secularEccentricReferenceAxis
        *(1.0-secularReferenceEccentricitySquared));
    const OrbitAveragedBmtAngularVelocities secularEccentricAverage=
        orbitAveragedBmtAngularVelocities(
            secularEccentricReferenceAxis,secularReferenceAngularMomentum,
            secularInitial.firstDipole,secularInitial.secondDipole,
            pairReducedMass);
    const Vec3 secularReferenceSeed=std::abs(secularReferenceNormal.x)<0.9
        ?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    Vec3 secularReferenceRadial=cross(
        secularReferenceSeed,secularReferenceNormal);
    secularReferenceRadial=secularReferenceRadial
        /secularReferenceRadial.norm();
    const Vec3 secularReferenceTangential=cross(
        secularReferenceNormal,secularReferenceRadial);
    Vec3 dimensionlessPosition=secularReferenceRadial
        *(1.0-secularReferenceEccentricity);
    Vec3 dimensionlessVelocity=secularReferenceTangential*std::sqrt(
        (1.0+secularReferenceEccentricity)
        /(1.0-secularReferenceEccentricity));
    const double secularReferenceMeanMotion=std::sqrt(pairCoulombStrength
        /(pairReducedMass*secularEccentricReferenceAxis
            *secularEccentricReferenceAxis*secularEccentricReferenceAxis));
    const double totalPairMass=firstMass+secondMass;
    const auto eccentricReferenceRates=[&](const Vec3& position,
                                            const Vec3& velocity) {
        State sample{};
        const Vec3 relativePosition=position*secularEccentricReferenceAxis;
        const Vec3 relativeVelocity=
            velocity*(secularEccentricReferenceAxis
                *secularReferenceMeanMotion);
        sample.firstPosition=
            relativePosition*(secondMass/totalPairMass);
        sample.secondPosition=
            relativePosition*(-firstMass/totalPairMass);
        sample.firstVelocity=
            relativeVelocity*(secondMass/totalPairMass);
        sample.secondVelocity=
            relativeVelocity*(-firstMass/totalPairMass);
        sample.firstDipole=secularInitial.firstDipole;
        sample.secondDipole=secularInitial.secondDipole;
        sample.firstProperDipole=secularInitial.firstDipole;
        sample.secondProperDipole=secularInitial.secondDipole;
        const StateHistory history{State{sample}};
        const LocalElectromagneticFields fields=
            localRelativisticFields(sample,history);
        return std::pair<Vec3,Vec3>{
            thomasBmtEffectiveField(
                sample.firstVelocity,fields.atFirst,firstGFactor)
                *(-firstCharge/firstMass),
            thomasBmtEffectiveField(
                sample.secondVelocity,fields.atSecond,secondGFactor)
                *(-secondCharge/secondMass)};
    };
    const auto dimensionlessKeplerDerivative=[](
            const Vec3& position,const Vec3& velocity) {
        const double radius=position.norm();
        return std::pair<Vec3,Vec3>{
            velocity,position*(-1.0/(radius*radius*radius))};
    };
    constexpr int resolvedEccentricOrbitSteps=32768;
    constexpr int mechanicalSubstepsPerRateSample=32;
    const double resolvedEccentricOrbitStep=
        2.0*pi/static_cast<double>(
            resolvedEccentricOrbitSteps*mechanicalSubstepsPerRateSample);
    Vec3 resolvedFirstRateSum,resolvedSecondRateSum;
    for(int step=0;step<resolvedEccentricOrbitSteps;++step) {
        const auto rates=eccentricReferenceRates(
            dimensionlessPosition,dimensionlessVelocity);
        resolvedFirstRateSum+=rates.first;
        resolvedSecondRateSum+=rates.second;
        for(int mechanicalStep=0;
            mechanicalStep<mechanicalSubstepsPerRateSample;++mechanicalStep) {
            const auto k1=dimensionlessKeplerDerivative(
                dimensionlessPosition,dimensionlessVelocity);
            const auto k2=dimensionlessKeplerDerivative(
                dimensionlessPosition
                    +k1.first*(0.5*resolvedEccentricOrbitStep),
                dimensionlessVelocity
                    +k1.second*(0.5*resolvedEccentricOrbitStep));
            const auto k3=dimensionlessKeplerDerivative(
                dimensionlessPosition
                    +k2.first*(0.5*resolvedEccentricOrbitStep),
                dimensionlessVelocity
                    +k2.second*(0.5*resolvedEccentricOrbitStep));
            const auto k4=dimensionlessKeplerDerivative(
                dimensionlessPosition+k3.first*resolvedEccentricOrbitStep,
                dimensionlessVelocity+k3.second*resolvedEccentricOrbitStep);
            dimensionlessPosition+=(k1.first+k2.first*2.0
                +k3.first*2.0+k4.first)*(resolvedEccentricOrbitStep/6.0);
            dimensionlessVelocity+=(k1.second+k2.second*2.0
                +k3.second*2.0+k4.second)*(resolvedEccentricOrbitStep/6.0);
        }
    }
    const Vec3 resolvedFirstRate=resolvedFirstRateSum
        *(1.0/static_cast<double>(resolvedEccentricOrbitSteps));
    const Vec3 resolvedSecondRate=resolvedSecondRateSum
        *(1.0/static_cast<double>(resolvedEccentricOrbitSteps));
    const double secularEccentricFirstResidual=
        (secularEccentricAverage.first-resolvedFirstRate).norm()
        /std::max(resolvedFirstRate.norm(),1.0e-300);
    const double secularEccentricSecondResidual=
        (secularEccentricAverage.second-resolvedSecondRate).norm()
        /std::max(resolvedSecondRate.norm(),1.0e-300);
    const double secularEccentricResolvedResidual=std::max(
        secularEccentricFirstResidual,secularEccentricSecondResidual);
    const double secularEccentricOrbitClosure=std::max(
        (dimensionlessPosition
            -secularReferenceRadial*(1.0-secularReferenceEccentricity)).norm(),
        (dimensionlessVelocity-secularReferenceTangential*std::sqrt(
            (1.0+secularReferenceEccentricity)
            /(1.0-secularReferenceEccentricity))).norm());
    const double secularFirstNormDrift=secularReference.completed
        ?std::abs(secularReference.state.firstDipole.norm()
            /secularInitial.firstDipole.norm()-1.0)
        :std::numeric_limits<double>::infinity();
    const double secularSecondNormDrift=secularReference.completed
        ?std::abs(secularReference.state.secondDipole.norm()
            /secularInitial.secondDipole.norm()-1.0)
        :std::numeric_limits<double>::infinity();
    const auto secularDifference=[&](const SecularSpinOrbitState& first,
                                     const SecularSpinOrbitState& second) {
        const double firstGyromagneticRatio=firstGyromagneticRatioOf();
        const double secondGyromagneticRatio=secondGyromagneticRatioOf();
        const double scale=secularInitial.orbitalAngularMomentum.norm()
            +secularInitial.firstDipole.norm()
                /std::abs(firstGyromagneticRatio)
            +secularInitial.secondDipole.norm()
                /std::abs(secondGyromagneticRatio);
        return ((first.orbitalAngularMomentum-second.orbitalAngularMomentum).norm()
            +(first.firstDipole-second.firstDipole).norm()
                /std::abs(firstGyromagneticRatio)
            +(first.secondDipole-second.secondDipole).norm()
                /std::abs(secondGyromagneticRatio))/scale;
    };
    const double secularCoarseFineDifference=
        secularDifference(secularCoarse.state,secularFine.state);
    const double secularFineReferenceDifference=
        secularDifference(secularFine.state,secularReference.state);
    const bool secularSpinOrbitIdentityOk=secularInitialRates.valid
        &&secularCoarse.completed&&secularFine.completed
        &&secularReference.completed
        &&secularCoarse.relativeAngularMomentumResidual<1.0e-12
        &&secularFine.relativeAngularMomentumResidual<1.0e-12
        &&secularReference.relativeAngularMomentumResidual<1.0e-12
        &&secularFirstNormDrift<1.0e-12&&secularSecondNormDrift<1.0e-12;
    const bool secularExternalTorqueOk=secularExternalRates.valid
        &&secularExternal.completed
        &&secularExternal.relativeAngularMomentumResidual<1.0e-12
        &&secularExternal.externalAngularMomentumTransfer.norm()
            >1.0e-6*hbar;
    const bool secularSpinOrbitConvergenceOk=secularSpinOrbitIdentityOk
        &&secularExternalTorqueOk
        &&secularFineReferenceDifference
            <0.6*secularCoarseFineDifference
        &&secularFineReferenceDifference<1.0e-3;
    const bool secularEccentricOrbitOk=secularEccentricAverage.valid
        &&std::isfinite(secularEccentricResolvedResidual)
        &&std::isfinite(secularEccentricOrbitClosure)
        &&secularEccentricResolvedResidual<1.0e-5
        &&secularEccentricOrbitClosure<1.0e-7;

    // ---------------------------------------------------------------------
    // Role routing in applyDipolePrecession().
    //
    // The covariance probe above builds its arguments explicitly, so it never
    // executes the call sites and cannot see them hand one particle's g-factor,
    // charge-to-mass ratio or local field to the other.  This check does: it
    // runs the production routine, then recomputes the same step by calling
    // advanceThomasBmtDipole() directly with each role's OWN parameters, and
    // requires the two to agree exactly.  (advanceCovariantBmt is no longer
    // the production formula -- see properDipolePrecessionRate's comment in
    // electrodynamics.hpp -- but the swapped-role failure mode this test
    // guards against is identical regardless of which formula is wired in.)
    //
    // The state is purpose-built rather than borrowed.  The first attempt
    // reused the covariance state, whose dipoles sit parallel to the local
    // field: the precession is then identically zero and the check compared
    // two unchanged vectors, passing against every injected fault.  Here the
    // moments are tilted off every field axis and the step is long enough to
    // turn them by ~1e-5 rad, which is far above round-off and far below any
    // nonlinearity.
    //
    // Scope limit worth stating: where the two roles share a property the swap
    // is unobservable in principle, not merely undetected.  With the default
    // e+e- pair both g-factors are 2.00231930436256, so a g swap yields
    // identical numbers and NO test can flag it; a charge-to-mass swap is
    // caught, the two roles carrying opposite signs.  For proton+electron both
    // become visible.
    State roleRoutingState;
    {
        const double radius=pairBohrRadius(activePair);
        const double firstShare=secondMass/(firstMass+secondMass);
        const double secondShare=firstMass/(firstMass+secondMass);
        const double circular=std::sqrt(pairCoulombStrength
            /(pairReducedMass*radius));
        roleRoutingState.firstPosition={firstShare*radius,0,0};
        roleRoutingState.secondPosition={-secondShare*radius,0,0};
        roleRoutingState.firstVelocity={0,firstShare*circular,0};
        roleRoutingState.secondVelocity={0,-secondShare*circular,0};
    }
    // Tilted off x, y and z alike, so no field component can leave the
    // precession accidentally zero.
    const Vec3 roleRoutingTilt=Vec3{1.0,1.0,1.0}/std::sqrt(3.0);
    roleRoutingState.firstDipole=roleRoutingTilt*firstMagneticMoment;
    roleRoutingState.secondDipole=roleRoutingTilt*secondMagneticMoment;
    ClassicalTrajectoryEngine roleRoutingEngine(roleRoutingState);
    constexpr double roleRoutingDt=1.0e-17;
    const StateHistory& roleRoutingHistory=roleRoutingEngine.history();
    State roleRoutingReference=roleRoutingState;
    synchronizeCovariantDipoles(roleRoutingReference);
    const LocalElectromagneticFields roleRoutingFields=
        localRelativisticFields(roleRoutingReference,roleRoutingHistory);
    const Vec3 expectedFirstProperDipole=advanceThomasBmtDipole(
        roleRoutingReference.firstProperDipole,
        roleRoutingReference.firstVelocity,roleRoutingFields.atFirst,
        firstCharge/firstMass,roleRoutingDt,firstGFactor);
    const Vec3 expectedSecondProperDipole=advanceThomasBmtDipole(
        roleRoutingReference.secondProperDipole,
        roleRoutingReference.secondVelocity,roleRoutingFields.atSecond,
        secondCharge/secondMass,roleRoutingDt,secondGFactor);
    applyDipolePrecession(roleRoutingState,roleRoutingDt,roleRoutingHistory);
    const double roleRoutingScale=std::max({
        expectedFirstProperDipole.norm(),expectedSecondProperDipole.norm(),
        1.0e-300});
    // How far the step actually turned the moments.  A vanishing value would
    // mean the check had gone vacuous again, so it is required to be finite
    // and non-trivial below rather than merely assumed.
    const double roleRoutingTravel=std::max(
        (expectedFirstProperDipole
            -roleRoutingReference.firstProperDipole).norm()/roleRoutingScale,
        (expectedSecondProperDipole
            -roleRoutingReference.secondProperDipole).norm()/roleRoutingScale);
    const double roleRoutingResidual=std::max(
        (roleRoutingState.firstProperDipole
            -expectedFirstProperDipole).norm()/roleRoutingScale,
        (roleRoutingState.secondProperDipole
            -expectedSecondProperDipole).norm()/roleRoutingScale);
    const ParticleMultipoleRadiation covarianceRestRadiation=
        particleMultipoleRadiation(covarianceRest,covarianceRestExternal,
            covarianceRestEngine.history(),false);
    const ParticleMultipoleRadiation covarianceMovingRadiation=
        particleMultipoleRadiation(covarianceMoving,covarianceMovingExternal,
            covarianceMovingEngine.history(),false);
    const FourVector expectedFirstForce=boostFourVector(fourForce(
        covarianceRest.firstVelocity,covarianceRestExternal.first
            +covarianceRestRadiation.chargeReaction.first));
    const FourVector movingFirstForce=fourForce(
        covarianceMoving.firstVelocity,covarianceMovingExternal.first
            +covarianceMovingRadiation.chargeReaction.first);
    const double covarianceForceResidual=
        (movingFirstForce.space-expectedFirstForce.space).norm()
        /std::max(expectedFirstForce.space.norm(),1.0e-300);
    const FourVector expectedRadiation=boostFourVector({
        covarianceRest.radiatedEnergy/c,covarianceRest.radiatedMomentum});
    const FourVector movingRadiation{covarianceMoving.radiatedEnergy/c,
                                    covarianceMoving.radiatedMomentum};
    const double covarianceAccumulatedRadiationResidual=std::max(
        std::abs(movingRadiation.time-expectedRadiation.time)
            /std::max(std::abs(expectedRadiation.time),1.0e-300),
        (movingRadiation.space-expectedRadiation.space).norm()
            /std::max({expectedRadiation.space.norm(),
                       std::abs(expectedRadiation.time),1.0e-300}));
    const FieldFluxRates covarianceRestFlux=electromagneticFieldFluxRates(
        covarianceRest,covarianceRestEngine.history(),
        {194,1.0e6*bohrRadius,true});
    const FieldFluxRates covarianceMovingFlux=electromagneticFieldFluxRates(
        covarianceMoving,covarianceMovingEngine.history(),
        {194,1.0e6*bohrRadius,true});
    FourVector expectedFlux=boostFourVector(
        {covarianceRestFlux.energy/c,covarianceRestFlux.momentum});
    expectedFlux.time/=covarianceBoostGamma;
    expectedFlux.space=expectedFlux.space/covarianceBoostGamma;
    const FourVector movingFlux{covarianceMovingFlux.energy/c,
                                covarianceMovingFlux.momentum};
    const double covarianceRadiationResidual=std::max(
        std::abs(movingFlux.time-expectedFlux.time)
            /std::max(std::abs(expectedFlux.time),1.0e-300),
        (movingFlux.space-expectedFlux.space).norm()
            /std::max({expectedFlux.space.norm(),std::abs(expectedFlux.time),
                       1.0e-300}));
    const DipoleTensor restDipole{{0,0,0},
        {0.31*bohrMagneton,-0.27*bohrMagneton,0.91*bohrMagneton}};
    const DipoleTensor boostedDipole=lorentzBoostDipole(
        restDipole,covarianceBoost);
    const DipoleTensor recoveredDipole=lorentzBoostDipole(
        boostedDipole,covarianceBoost*-1.0);
    const double dipoleTensorRoundtrip=std::max(
        (recoveredDipole.electric-restDipole.electric).norm()
            /std::max(restDipole.magnetic.norm()/c,1.0e-300),
        (recoveredDipole.magnetic-restDipole.magnetic).norm()
            /restDipole.magnetic.norm());
    const double dipoleInvariantScale=restDipole.magnetic.squaredNorm();
    const double dipoleFirstInvariantResidual=std::abs(
        dipoleFirstInvariant(boostedDipole)-dipoleFirstInvariant(restDipole))
        /dipoleInvariantScale;
    const double dipoleSecondInvariantResidual=std::abs(
        dipoleSecondInvariant(boostedDipole)-dipoleSecondInvariant(restDipole))
        /dipoleInvariantScale;
    const Vec3 expectedInducedElectric=cross(
        covarianceBoost,restDipole.magnetic)
        *(covarianceBoostGamma/(c*c));
    const double inducedElectricDipoleResidual=
        (boostedDipole.electric-expectedInducedElectric).norm()
        /std::max(expectedInducedElectric.norm(),1.0e-300);
    State staticElectricDipoleState;
    staticElectricDipoleState.firstPosition={2.0*bohrRadius,0,0};
    staticElectricDipoleState.secondPosition={0,0,0};
    staticElectricDipoleState.secondElectricDipole={
        0,0,eCharge*bohrRadius};
    State staticElectricDipolePast=staticElectricDipoleState;
    staticElectricDipolePast.time=-8.0*bohrRadius/c;
    const ElectromagneticField staticElectricDipoleField=
        retardedElectricDipoleField(
            staticElectricDipoleState.firstPosition,0.0,
            {staticElectricDipolePast,staticElectricDipoleState},
            staticElectricDipoleState,false);
    const Vec3 expectedStaticElectricDipole={0,0,
        -coulomb*eCharge*bohrRadius/std::pow(2.0*bohrRadius,3)};
    const double staticElectricDipoleResidual=
        (staticElectricDipoleField.electric-expectedStaticElectricDipole).norm()
        /expectedStaticElectricDipole.norm();
    State gradientDipoleState;
    gradientDipoleState.firstPosition={2.0*bohrRadius,0,0};
    gradientDipoleState.secondPosition={0,0,0};
    gradientDipoleState.firstProperDipole={0.3*bohrMagneton,
                                              0.2*bohrMagneton,0};
    gradientDipoleState.secondProperDipole={-0.1*bohrMagneton,
                                               0.4*bohrMagneton,0};
    synchronizeCovariantDipoles(gradientDipoleState);
    State gradientDipolePast=gradientDipoleState;
    gradientDipolePast.time=-8.0*bohrRadius/c;
    const StateHistory gradientDipoleHistory{
        gradientDipolePast,gradientDipoleState};
    const Vec3 tensorGradientForce=covariantDipoleGradientForce(
        gradientDipoleState,gradientDipoleHistory,true);
    const Vec3 analyticGradientForce=regularizedDipoleForce(
        gradientDipoleState.firstPosition
            -gradientDipoleState.secondPosition,
        gradientDipoleState.firstProperDipole,
        gradientDipoleState.secondProperDipole);
    const double tensorGradientStaticResidual=
        (tensorGradientForce-analyticGradientForce).norm()
        /std::max(analyticGradientForce.norm(),1.0e-300);
    // The static check above holds the target at rest, where the lab
    // magnetic and electric dipole channels the fix generalized
    // covariantDipoleGradientForce to (mu_lab.B - p_lab.E, see its comment)
    // cannot be told apart from the old rest-frame-only mu_proper.B_rest
    // recipe: p_rest=0 here, so both give the same static number.  Boost the
    // very same configuration with the charge sector's own
    // covarianceBoost/boostEvent/boostVelocity/lorentzBoostDipole and check
    // that the coupling itself -- mu_lab.B(x)-p_lab.E(x), now with a real
    // induced p_lab from the boost -- is the Lorentz invariant it is built
    // to be, i.e. has the SAME value at the same spacetime event in both
    // frames.  This is what actually exercises and pins down the relative
    // sign between the two channels: a plus sign there is off by two orders
    // of magnitude, the minus sign below matches to five significant
    // figures (the small remainder is discretization/retardation-solver
    // noise, not a channel imbalance).  It does not, on its own, validate
    // the fixed-lab-time spatial gradient the force is built from for a
    // moving target -- see covariantDipoleGradientForce's own comment on
    // what remains open there.
    const auto boostGradientDipoleState=[&](const State& source) {
        State boosted=source;
        const auto [firstPos,firstT]=
            boostEvent(source.firstPosition,source.time);
        const auto [secondPos,secondT]=
            boostEvent(source.secondPosition,source.time);
        boosted.firstPosition=firstPos;
        boosted.secondPosition=secondPos;
        boosted.firstVelocity=boostVelocity(source.firstVelocity);
        boosted.secondVelocity=boostVelocity(source.secondVelocity);
        const DipoleTensor firstTensor=lorentzBoostDipole(
            {source.firstElectricDipole,source.firstDipole},covarianceBoost);
        const DipoleTensor secondTensor=lorentzBoostDipole(
            {source.secondElectricDipole,source.secondDipole},covarianceBoost);
        boosted.firstDipole=firstTensor.magnetic;
        boosted.firstElectricDipole=firstTensor.electric;
        boosted.secondDipole=secondTensor.magnetic;
        boosted.secondElectricDipole=secondTensor.electric;
        boosted.firstProperDipole=properDipoleFromFourVector(
            boostFourVector(dipoleFourVector(source.firstProperDipole,
                source.firstVelocity)),boosted.firstVelocity);
        boosted.secondProperDipole=properDipoleFromFourVector(
            boostFourVector(dipoleFourVector(source.secondProperDipole,
                source.secondVelocity)),boosted.secondVelocity);
        boosted.time=0.5*(firstT+secondT);
        return boosted;
    };
    const State boostedGradientState=
        boostGradientDipoleState(gradientDipoleState);
    const State boostedGradientPast=
        boostGradientDipoleState(gradientDipolePast);
    const StateHistory boostedGradientHistory{
        boostedGradientPast,boostedGradientState};
    const ElectromagneticField restCouplingField=fieldFromOtherParticleAt(
        gradientDipoleState.firstPosition,gradientDipoleState.time,
        gradientDipoleState,gradientDipoleHistory,true);
    const double restCoupling=
        dot(gradientDipoleState.firstDipole,restCouplingField.magnetic)
        -dot(gradientDipoleState.firstElectricDipole,
            restCouplingField.electric);
    const auto [boostedFirstEventPosition,boostedFirstEventTime]=
        boostEvent(gradientDipoleState.firstPosition,gradientDipoleState.time);
    const ElectromagneticField boostedCouplingField=fieldFromOtherParticleAt(
        boostedFirstEventPosition,boostedFirstEventTime,
        boostedGradientState,boostedGradientHistory,true);
    const double boostedCoupling=
        dot(boostedGradientState.firstDipole,boostedCouplingField.magnetic)
        -dot(boostedGradientState.firstElectricDipole,
            boostedCouplingField.electric);
    const double dipoleGradientCouplingInvarianceResidual=
        std::abs(boostedCoupling-restCoupling)
        /std::max(std::abs(restCoupling),1.0e-300);
    State quadrupolePair=yeeCoupledState;
    quadrupolePair.firstPosition={-1.7*bohrRadius,0.4*bohrRadius,
                                      -0.2*bohrRadius};
    quadrupolePair.secondPosition={0.9*bohrRadius,-0.8*bohrRadius,
                                      0.6*bohrRadius};
    const ElectricQuadrupole pairQuadrupole=
        electricQuadrupole(quadrupolePair);
    const Vec3 quadrupoleSeparation=quadrupolePair.firstPosition
                                   -quadrupolePair.secondPosition;
    const double quadrupoleScale=eCharge
        *quadrupoleSeparation.squaredNorm();
    // Analytic two-body quadrupole about the centre of mass,
    // Q_ij = kappa (3 d_i d_j - d^2 delta_ij).  Comparing against THIS rather
    // than against zero is what makes the check a measurement again.  The old
    // form demanded that Q vanish, which the equal-mass midpoint origin
    // guaranteed on its own for every neutral pair -- so it could not tell
    // positronium's genuine cancellation apart from a channel that had been
    // deleted outright, and reported the same 0 for p+e- where the quadrupole
    // is at nearly full strength.
    //
    // kappa carries the entire pair dependence and vanishes exactly when the
    // two masses are equal: 0 for e+e- and mu+mu-, -0.9989 e for p+e-.  That
    // is the real statement "positronium has no E2 channel", now derived from
    // the masses instead of assumed by the choice of origin.
    const double quadrupoleCoefficient=
        (firstCharge*secondMass*secondMass
        +secondCharge*firstMass*firstMass)
        /((firstMass+secondMass)*(firstMass+secondMass));
    ElectricQuadrupole expectedQuadrupole;
    {
        const std::array<double,3> offset{quadrupoleSeparation.x,
            quadrupoleSeparation.y,quadrupoleSeparation.z};
        const double offsetSquared=quadrupoleSeparation.squaredNorm();
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
            expectedQuadrupole.component[static_cast<std::size_t>(3*i+j)]=
                quadrupoleCoefficient*(3.0
                    *offset[static_cast<std::size_t>(i)]
                    *offset[static_cast<std::size_t>(j)]
                  -(i==j?offsetSquared:0.0));
    }
    double quadrupoleTrace=pairQuadrupole.component[0]
        +pairQuadrupole.component[4]+pairQuadrupole.component[8];
    double quadrupoleSymmetry=0.0;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        quadrupoleSymmetry=std::max(quadrupoleSymmetry,std::abs(
            pairQuadrupole.component[static_cast<std::size_t>(3*i+j)]
           -pairQuadrupole.component[static_cast<std::size_t>(3*j+i)]));
    const double quadrupoleResidual=std::sqrt(
        (pairQuadrupole-expectedQuadrupole).squaredNorm())/quadrupoleScale;
    const double quadrupoleMagnitude=std::sqrt(
        pairQuadrupole.squaredNorm())/quadrupoleScale;
    quadrupoleTrace=std::abs(quadrupoleTrace)/quadrupoleScale;
    quadrupoleSymmetry/=quadrupoleScale;

    constexpr double dipoleDerivativeStep=1.0e-22;
    const double derivativeCube=dipoleDerivativeStep*dipoleDerivativeStep
        *dipoleDerivativeStep;
    const Vec3 cubicCoefficient{0.15*bohrRadius/derivativeCube,
                                -0.08*bohrRadius/derivativeCube,0};
    StateHistory cubicDipoleHistory;
    for(int index=-8;index<=0;++index) {
        State sample;
        sample.time=0.5*dipoleDerivativeStep*index;
        sample.secondPosition=Vec3{bohrRadius,0,0}
            +cubicCoefficient*(sample.time*sample.time*sample.time);
        cubicDipoleHistory.push_back(sample);
    }
    const State cubicDipolePresent=cubicDipoleHistory.back();
    // Only the SECOND role moves in this probe, so d''' = q2 r2''' and the
    // expected value carries the second charge, not the elementary charge.
    // Writing +e was right only while the second role was the positron; with
    // an electron there the measured third derivative came out exactly
    // negated and the residual sat at 2.
    const Vec3 expectedDipoleThird=cubicCoefficient*(6.0*secondCharge);
    const Vec3 measuredDipoleThird=electricDipoleThirdDerivative(
        cubicDipolePresent,cubicDipoleHistory);
    const double coherentDerivativeResidual=(measuredDipoleThird
        -expectedDipoleThird).norm()/expectedDipoleThird.norm();
    const MutualForces coherentProbe=coherentElectricDipoleReaction(
        cubicDipolePresent,cubicDipoleHistory);
    const double coherentReactionMomentumResidual=(coherentProbe.first
        +coherentProbe.second).norm()
        /std::max(coherentProbe.first.norm(),1.0e-300);

    // Production-path regression for M1 source interference.  Collinear,
    // positive moments keep interpolateDipole() exactly linear.  Thus
    // m_i(t)=100 +/- t^2 gives m1''=2 and m2''=+/-2 through the real history
    // stencil: equal derivatives must radiate 4 P1 and opposite derivatives
    // must radiate nothing.
    const auto magneticInterferenceHistory=[](bool opposite) {
        StateHistory history;
        for(int index=0;index<=4;++index) {
            State sample;
            sample.time=static_cast<double>(index);
            const double square=sample.time*sample.time;
            sample.firstDipole={100.0+square,0,0};
            sample.secondDipole={100.0+(opposite?-square:square),0,0};
            history.push_back(sample);
        }
        return history;
    };
    const StateHistory alignedMagneticHistory=
        magneticInterferenceHistory(false);
    const StateHistory cancellingMagneticHistory=
        magneticInterferenceHistory(true);
    const State& alignedMagneticState=alignedMagneticHistory.back();
    const State& cancellingMagneticState=cancellingMagneticHistory.back();
    const RetardedDipoleKinematics alignedFirstMagnetic=
        historicalDipoleKinematics(alignedMagneticHistory,
            alignedMagneticState,true,alignedMagneticState.time);
    const RetardedDipoleKinematics alignedSecondMagnetic=
        historicalDipoleKinematics(alignedMagneticHistory,
            alignedMagneticState,false,alignedMagneticState.time);
    const RetardedDipoleKinematics cancellingSecondMagnetic=
        historicalDipoleKinematics(cancellingMagneticHistory,
            cancellingMagneticState,false,cancellingMagneticState.time);
    const DipoleRadiationReaction alignedMagneticRadiation=
        dipoleRadiationReaction(alignedMagneticState,alignedMagneticHistory);
    const DipoleRadiationReaction cancellingMagneticRadiation=
        dipoleRadiationReaction(
            cancellingMagneticState,cancellingMagneticHistory);
    constexpr double magneticRadiationCoefficient=mu0/(6.0*pi*c*c*c);
    const double singleMagneticPower=magneticRadiationCoefficient
        *alignedFirstMagnetic.secondDerivative.squaredNorm();
    const double alignedMagneticPowerRatio=
        alignedMagneticRadiation.power/std::max(singleMagneticPower,1.0e-300);
    const double cancellingMagneticPowerRatio=
        cancellingMagneticRadiation.power/std::max(singleMagneticPower,1.0e-300);
    const double magneticDerivativeResidual=std::max({
        (alignedFirstMagnetic.secondDerivative-Vec3{2,0,0}).norm()/2.0,
        (alignedSecondMagnetic.secondDerivative-Vec3{2,0,0}).norm()/2.0,
        (cancellingSecondMagnetic.secondDerivative-Vec3{-2,0,0}).norm()/2.0});

    // Direct algebra probe complements the production-history check above: a
    // derivative carried by the second source must torque the first moment,
    // angular flux uses M' x M'', and a common boost transports the coherent
    // power with the velocity of the pair as a whole.
    RetardedDipoleKinematics magneticAlgebraFirst;
    RetardedDipoleKinematics magneticAlgebraSecond;
    const Vec3 magneticCommonVelocity{0.125*c,-0.04*c,0.02*c};
    const Vec3 magneticSourcePosition{10.0,-4.0,2.0};
    magneticAlgebraFirst.position=magneticSourcePosition;
    magneticAlgebraSecond.position=magneticSourcePosition;
    magneticAlgebraFirst.velocity=magneticCommonVelocity;
    magneticAlgebraSecond.velocity=magneticCommonVelocity;
    magneticAlgebraFirst.moment={1,0,0};
    magneticAlgebraSecond.moment={0,1,0};
    magneticAlgebraFirst.firstDerivative={1,0,0};
    magneticAlgebraFirst.secondDerivative={0,1,0};
    magneticAlgebraSecond.thirdDerivative={0,1,0};
    const DipoleRadiationReaction magneticAlgebra=
        coherentMagneticDipoleRadiationReaction(
            magneticAlgebraFirst,magneticAlgebraSecond);
    const Vec3 expectedMagneticTorque{0,0,magneticRadiationCoefficient};
    const Vec3 expectedMagneticMomentum=magneticCommonVelocity
        *(magneticAlgebra.power/(c*c));
    const Vec3 expectedMagneticAngularFlux=expectedMagneticTorque
        +cross(magneticSourcePosition,expectedMagneticMomentum);
    const double magneticTorqueResidual=std::max(
        (magneticAlgebra.firstTorque-expectedMagneticTorque).norm()
            /magneticRadiationCoefficient,
        magneticAlgebra.secondTorque.norm()/magneticRadiationCoefficient);
    const double magneticAngularFluxResidual=
        (magneticAlgebra.angularMomentumRate-expectedMagneticAngularFlux).norm()
            /std::max(expectedMagneticAngularFlux.norm(),1.0e-300);
    const double magneticMomentumResidual=
        (magneticAlgebra.momentumRate-expectedMagneticMomentum).norm()
            /std::max(expectedMagneticMomentum.norm(),1.0e-300);
    RetardedDipoleKinematics magneticComFirst=magneticAlgebraFirst;
    RetardedDipoleKinematics magneticComSecond=magneticAlgebraSecond;
    const Vec3 magneticComMomentum{
        0.05*std::min(firstMass,secondMass)*c,0,0};
    magneticComFirst.velocity=velocityFromMomentum(
        magneticComMomentum,firstMass);
    magneticComSecond.velocity=velocityFromMomentum(
        magneticComMomentum*-1.0,secondMass);
    const DipoleRadiationReaction magneticComRadiation=
        coherentMagneticDipoleRadiationReaction(
            magneticComFirst,magneticComSecond);
    const double magneticComMomentumResidual=
        magneticComRadiation.momentumRate.norm()
            /std::max(magneticComRadiation.power/c,1.0e-300);

    // Convergence matrix for the production far-zone surface integral.  The
    // high-resolution, acceleration-field-only result is the reference: it
    // excludes finite-radius 1/R^2 fields without changing the source model.
    const auto fluxResidual=[](const FieldFluxRates& value,
                               const FieldFluxRates& reference) {
        const double energyScale=std::max(std::abs(reference.energy),1.0e-300);
        const double momentumScale=std::max({reference.momentum.norm(),
            energyScale/c,1.0e-300});
        const double angularScale=std::max({reference.angularMomentum.norm(),
            energyScale*bohrRadius/c,1.0e-300});
        return std::max({std::abs(value.energy-reference.energy)/energyScale,
            (value.momentum-reference.momentum).norm()/momentumScale,
            (value.angularMomentum-reference.angularMomentum).norm()/angularScale});
    };
    const auto radiativeFluxResidual=[](const FieldFluxRates& value,
                                        const FieldFluxRates& reference) {
        const double energyScale=std::max(std::abs(reference.energy),1.0e-300);
        const double momentumScale=std::max({reference.momentum.norm(),
            energyScale/c,1.0e-300});
        return std::max(std::abs(value.energy-reference.energy)/energyScale,
            (value.momentum-reference.momentum).norm()/momentumScale);
    };
    const FarFieldSampling farReferenceSampling{194,1.0e6*bohrRadius,true};
    const std::vector<SphereQuadraturePoint> lebedevRule=sphereQuadrature(50);
    double lebedevWeight=0.0;
    Vec3 lebedevFirstMoment;
    std::array<double,9> lebedevSecondMoment{};
    for(const SphereQuadraturePoint& point:lebedevRule) {
        lebedevWeight+=point.solidAngleWeight;
        lebedevFirstMoment+=point.direction*point.solidAngleWeight;
        const std::array<double,3> n{point.direction.x,point.direction.y,
                                     point.direction.z};
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
            lebedevSecondMoment[static_cast<std::size_t>(3*i+j)]+=
                n[static_cast<std::size_t>(i)]*n[static_cast<std::size_t>(j)]
                *point.solidAngleWeight;
    }
    const double lebedevWeightResidual=std::abs(lebedevWeight-4.0*pi)/(4.0*pi);
    double lebedevMomentResidual=lebedevFirstMoment.norm()/(4.0*pi);
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        lebedevMomentResidual=std::max(lebedevMomentResidual,std::abs(
            lebedevSecondMoment[static_cast<std::size_t>(3*i+j)]
            -(i==j?4.0*pi/3.0:0.0))/(4.0*pi));
    const FieldFluxRates farReference=electromagneticFieldFluxRates(
        sharedEngineVisualState,visualEngine.history(),farReferenceSampling);
    const FieldFluxRates farFullReference=electromagneticFieldFluxRates(
        sharedEngineVisualState,visualEngine.history(),
        {194,farReferenceSampling.controlRadius,false});
    const std::array<int,3> farDirectionCounts{26,50,98};
    std::array<double,3> farDirectionResiduals{};
    for(std::size_t index=0;index<farDirectionCounts.size();++index) {
        farDirectionResiduals[index]=fluxResidual(
            electromagneticFieldFluxRates(sharedEngineVisualState,
                visualEngine.history(),{farDirectionCounts[index],
                    farReferenceSampling.controlRadius,false}),farFullReference);
    }
    // --- Far-zone quadrature against the closed-form Larmor rate ---
    //
    // Every other far-field check here compares the quadrature with ITSELF at
    // a different radius or direction count, so all of them pass while the
    // absolute normalization is wrong.  This one supplies the missing external
    // reference: a circular e+e- orbit radiates as a pure electric dipole
    // d = e r with the orbit-averaged power
    //
    //     P = |d''|^2/(6 pi eps0 c^3),   |d''| = 2 k e^3/(m r^2),
    //
    // and the charge quadrupole vanishes identically for a symmetric pair, so
    // this is the whole charge-sector radiation with no free parameters.  The
    // dipoles are set to zero so the M1 channel contributes nothing and the
    // flux needs no decomposition.
    double larmorNormalizationRatio=std::numeric_limits<double>::quiet_NaN();
    double larmorProductionRatio=std::numeric_limits<double>::quiet_NaN();
    double larmorAccumulationRatio=std::numeric_limits<double>::quiet_NaN();
    double longHorizonBalanceResidual=
        std::numeric_limits<double>::quiet_NaN();
    double longHorizonRadiatedEnergy=
        std::numeric_limits<double>::quiet_NaN();
    double longHorizonSchottChange=
        std::numeric_limits<double>::quiet_NaN();
    struct BalanceDiagnostic {
        const char* name="";
        double mechanicalOverFlux=std::numeric_limits<double>::quiet_NaN();
        double schottOverFlux=std::numeric_limits<double>::quiet_NaN();
        double signedResidual=std::numeric_limits<double>::quiet_NaN();
        double fluxEnergy=std::numeric_limits<double>::quiet_NaN();
        double schottChange=std::numeric_limits<double>::quiet_NaN();
        bool advanced=false;
    };
    std::array<BalanceDiagnostic,4> balanceDiagnostics{};
    // Same balance, same trajectory, one refinement step finer.  The residual
    // this probe reports is not a physical imbalance but a discretization
    // error amplified by the measurement's own conditioning (see
    // longHorizonBalanceAmplification below), so the number is meaningless
    // without knowing whether it shrinks when the discretization does.
    double longHorizonBalanceRefined=std::numeric_limits<double>::quiet_NaN();
    // |E_mech| / E_far for this probe: the factor by which ANY relative error
    // in the mechanical energy is magnified in the residual.
    double longHorizonBalanceAmplification=
        std::numeric_limits<double>::quiet_NaN();
    // The same balance swept across a factor of 16 in orbit radius, each
    // point measured coarse and refined.  Two things at once.
    //
    // The conditioning claim is not specific to one radius: the amplification
    // is |E_mech|/E_far with E_mech proportional to 1/R and, by Larmor,
    // E_far per orbit proportional to R^(-5/2), so it must grow as R^(3/2).
    // Checking that the MEASURED amplification follows that law is therefore
    // a check on the radiated energy's own scaling, across a decade and a
    // half of radius, and not merely bookkeeping about the probe.
    //
    // And the residual must shrink under refinement at EVERY radius, not
    // just at the one the headline number is quoted from.
    struct BalanceRadiusPoint {
        double factor=0.0,amplification=0.0,coarse=0.0,refined=0.0;
        bool advanced=false;
    };
    std::array<BalanceRadiusPoint,3> balanceRadiusSweep{};
    // Continuous mechanical radiative ENERGY drain that survives in the fully
    // quantized mode.  Must be exactly zero: under stochasticElectricDipole
    // every energy channel leaves as discrete quanta, so no energy may be
    // removed between photons -- neither the charge sector's reaction force
    // nor the magnetic sector's dipoleConstraintEnergy drain.
    //
    // The reaction TORQUE is deliberately NOT on that list any more, and this
    // clause used to say it was.  The stale wording cost a wrong turn during
    // the symmetry audit: reading it, the M1 torque-on/drain-off pairing looks
    // like an unpaired conjugate slot, and the obvious repair is to turn the
    // drain on.  Two versions of that were tried and both take validation to
    // 38/39 -- "always drain" and "drain iff the torque is applied" alike --
    // because this very check forbids any continuous energy removal in the
    // quantized mode.  The torque is a different conjugate slot and was
    // separated from the drain on purpose; see the note at
    // quantizedDipoleTorqueTravel below for the explicit gate probe.
    double quantizedChargeReactionDrain=
        std::numeric_limits<double>::quiet_NaN();
    double quantizedDipoleConstraintDrain=
        std::numeric_limits<double>::quiet_NaN();
    double quantizedDipoleTorqueTravel=
        std::numeric_limits<double>::quiet_NaN();
    double disabledDipoleTorqueTravel=
        std::numeric_limits<double>::quiet_NaN();
    double quantizedDipoleTorqueNormDrift=
        std::numeric_limits<double>::quiet_NaN();
    {
        const double larmorProbeRadius=pairBohrRadius(activePair);
        const double larmorReducedMass=pairReducedMass;
        const double circularRelativeSpeed=std::sqrt(
            pairCoulombStrength/(larmorReducedMass*larmorProbeRadius));
        State larmorState;
        // Mass-weighted, not halved: the analytic reference below is the
        // orbit-averaged power of a pair whose centre of mass is at rest, so
        // the probe state has to be one.  Halving is the equal-mass special
        // case and silently sets the pair drifting for any other, which drove
        // the measured ratio to 8e5 instead of 1 on proton+electron.
        const double firstShare=secondMass/(firstMass+secondMass);
        const double secondShare=firstMass/(firstMass+secondMass);
        larmorState.firstPosition={firstShare*larmorProbeRadius,0,0};
        larmorState.secondPosition={-secondShare*larmorProbeRadius,0,0};
        larmorState.firstVelocity={0,firstShare*circularRelativeSpeed,0};
        larmorState.secondVelocity={0,-secondShare*circularRelativeSpeed,0};
        // Let the engine build a genuine retarded history rather than reading
        // the reconstructed one, so the comparison tests the quadrature and
        // not causalInitialHistory().
        ClassicalTrajectoryEngine larmorEngine(larmorState);
        const double larmorPeriod=2.0*pi*std::sqrt(
            larmorReducedMass*larmorProbeRadius*larmorProbeRadius
            *larmorProbeRadius/(pairCoulombStrength));
        bool larmorAdvanced=true;
        for(int step=0;step<32&&larmorAdvanced;++step)
            larmorAdvanced=larmorEngine.advance(larmorState,larmorPeriod/512.0);
        const FieldFluxRates larmorFlux=electromagneticFieldFluxRates(
            larmorState,larmorEngine.history(),farReferenceSampling);
        // |d''| = |q_eff| k|q1 q2| / (mu a^2), with q_eff the pair's effective
        // dipole charge and mu the reduced mass.  The earlier form,
        // 2 k e^3/(m a^2), is the equal-mass special case: it silently swaps mu
        // for m/2 and q_eff for e, which for proton+electron overstates the
        // denominator by m_p/(2 m_e) = 918 and drove this ratio to 8e5.
        const double dipoleSecondDerivative=magnitude(pairDipoleCharge)
            *pairCoulombStrength
            /(pairReducedMass*larmorProbeRadius*larmorProbeRadius);
        const double analyticLarmorPower=dipoleSecondDerivative
            *dipoleSecondDerivative/(6.0*pi*epsilon0*c*c*c);
        // Same orbit, same instant, but sampled the way the PRODUCTION path
        // samples it: default FarFieldSampling, i.e. 50 directions and
        // radiationFieldOnly=false so the velocity (near-field) terms are
        // included in the Poynting integral.
        const FieldFluxRates larmorProductionFlux=
            electromagneticFieldFluxRates(larmorState,larmorEngine.history());
        if(larmorAdvanced&&std::isfinite(larmorFlux.energy)
           &&analyticLarmorPower>0.0) {
            larmorNormalizationRatio=larmorFlux.energy/analyticLarmorPower;
            larmorProductionRatio=
                larmorProductionFlux.energy/analyticLarmorPower;
        }
        // Second, separate check: the ACCUMULATED radiatedEnergy over a known
        // stretch of the same circular orbit, against the analytic energy for
        // that stretch.  The ratio above tests the quadrature at one instant;
        // this tests the time integration of it, which is where the bookkeeping
        // used to lose 34% by evaluating the flux at a predictor state.  A
        // check of the instantaneous quadrature alone cannot see that.
        // Run this part at the accuracy and step size the STATISTICAL and
        // VISUAL paths actually use.  The size of the predictor-state error
        // depends on dt relative to the wavefront offset sourceExtent/c, so a
        // finer probe than production simply does not exercise it: at
        // period/512 with tolerance 1e-7 the old trapezoid came out only 1.2%
        // low and slipped through, while at production settings it is 34% low.
        State accumulationState=larmorState;
        ClassicalTrajectoryEngine accumulationEngine(accumulationState,
            {.relativeTolerance=1.0e-5,.maximumDepth=12});
        const double accumulationStart=accumulationState.radiatedEnergy;
        const double accumulationTime=accumulationState.time;
        bool accumulationAdvanced=larmorAdvanced;
        for(int step=0;step<64&&accumulationAdvanced;++step)
            accumulationAdvanced=accumulationEngine.advance(
                accumulationState,larmorPeriod/128.0);
        const double accumulationElapsed=accumulationState.time-accumulationTime;
        if(accumulationAdvanced&&accumulationElapsed>0.0
           &&analyticLarmorPower>0.0) {
            larmorAccumulationRatio=
                (accumulationState.radiatedEnergy-accumulationStart)
                /(analyticLarmorPower*accumulationElapsed);
        }

        // Independent long-horizon balance.  Unlike pair-field IDENTITY this
        // never reads boundFieldEnergy: mechanical energy comes from the
        // conservative Noether expression, radiation from the actual retarded
        // far-sphere Poynting integral accumulated by the engine, and the
        // reversible near-field endpoint term from explicitChargeSchottEnergy.
        // Two complete orbital periods are long compared with the one-step
        // reaction probes above while keeping the validator practical.
        double balanceRadiusFactor=0.2;
        const double balanceRadius=0.2*larmorProbeRadius;
        // The coherent reaction needs enough retained history nodes for its
        // third derivative.  period/512 is below the retarded-history spacing
        // ceiling at this radius; the production-like period/128 step leaves
        // only about three nodes and would intentionally gate that force off.
        constexpr int balancePeriods=2;
        constexpr int balanceStepsPerPeriod=512;
        double balanceTolerance=1.0e-6;
        int balancePeriodCount=balancePeriods;
        const auto runBalance=[&](const char* name,
                                  ChargeRadiationReactionModel reaction,
                                  bool retarded) {
            BalanceDiagnostic diagnostic;
            diagnostic.name=name;
            const double sweptRadius=balanceRadiusFactor*larmorProbeRadius;
            const double sweptSpeed=std::sqrt(pairCoulombStrength
                /(larmorReducedMass*sweptRadius));
            const double sweptPeriod=2.0*pi*std::sqrt(larmorReducedMass
                *sweptRadius*sweptRadius*sweptRadius/pairCoulombStrength);
            State state;
            state.firstPosition={firstShare*sweptRadius,0,0};
            state.secondPosition={-secondShare*sweptRadius,0,0};
            state.firstVelocity={0,firstShare*sweptSpeed,0};
            state.secondVelocity={0,-secondShare*sweptSpeed,0};
            ClassicalTrajectoryEngine engine(state,
                {.relativeTolerance=balanceTolerance,.maximumDepth=14,
                 .reactionModel=reaction,
                 .computeOutwardFlux=true,
                 .useRetardedExternalForces=retarded});
            const double mechanicalStart=conservativeParticleEnergy(state);
            const double fluxStart=state.radiatedEnergy;
            const MutualForces initialExternal=retarded
                ?retardedExternalForces(state,engine.history())
                :allExternalForces(state);
            const double schottStart=explicitChargeSchottEnergy(
                state,initialExternal);
            bool advanced=true;
            for(int step=0;step<balancePeriodCount*balanceStepsPerPeriod
                &&advanced;++step) {
                advanced=engine.advance(
                    state,sweptPeriod/balanceStepsPerPeriod);
            }
            diagnostic.advanced=advanced&&isFinite(state);
            if(!diagnostic.advanced) return diagnostic;
            const MutualForces finalExternal=retarded
                ?retardedExternalForces(state,engine.history())
                :allExternalForces(state);
            const double mechanicalChange=
                conservativeParticleEnergy(state)-mechanicalStart;
            const double flux=state.radiatedEnergy-fluxStart;
            const double schott=explicitChargeSchottEnergy(state,finalExternal)
                -schottStart;
            diagnostic.fluxEnergy=flux;
            diagnostic.schottChange=schott;
            diagnostic.mechanicalOverFlux=mechanicalChange
                /std::max(std::abs(flux),1.0e-300);
            diagnostic.schottOverFlux=schott
                /std::max(std::abs(flux),1.0e-300);
            diagnostic.signedResidual=(mechanicalChange+flux+schott)
                /std::max(std::abs(flux),1.0e-300);
            return diagnostic;
        };
        // --------------------------------------------------------------
        // Full-quantization guard.  The suite's other radiation probes all
        // run continuous reaction models, so none of them executes the
        // stochastic path at all: before this check, gating a continuous
        // sink incorrectly (or forgetting to gate one, which double-counts
        // that channel's energy -- once continuously, once as photons) was
        // invisible to every enforced check here.
        {
            State quantizedState;
            quantizedState.firstPosition={firstShare*larmorProbeRadius,0,0};
            quantizedState.secondPosition={-secondShare*larmorProbeRadius,0,0};
            quantizedState.firstVelocity={0,firstShare*circularRelativeSpeed,0};
            quantizedState.secondVelocity={0,-secondShare*circularRelativeSpeed,0};
            // Tilt the moments off every field axis, so a radiation torque
            // would actually move them if one were still being applied.
            quantizedState.firstProperDipole=
                Vec3{0.3,0.5,0.81}*(firstMagneticMoment/0.99929975);
            quantizedState.secondProperDipole=
                Vec3{-0.4,0.62,0.67}*(secondMagneticMoment/1.00344407);
            synchronizeCovariantDipoles(quantizedState);
            ClassicalTrajectoryEngine quantizedEngine(quantizedState,
                {.relativeTolerance=1.0e-8,.maximumDepth=12,
                 .reactionModel=ChargeRadiationReactionModel
                     ::stochasticElectricDipole,
                 .computeOutwardFlux=true,
                 .useRetardedExternalForces=true});
            const double constraintStart=quantizedState.dipoleConstraintEnergy;
            // Exercise the exact model gate used by the production step with
            // a controlled transverse torque.  The physical M1 torque at
            // this radius only separates the two trajectories at round-off;
            // scaling the input here makes the routing decision observable
            // without changing the production dynamics or requiring a full
            // inspiral.  Renormalization must preserve both dipole norms.
            const State torqueStart=quantizedState;
            ParticleMultipoleRadiation torqueProbe;
            constexpr double torqueProbeStep=1.0e-24;
            const Vec3 firstTransverse=cross(
                torqueStart.firstProperDipole,Vec3{0.0,0.0,1.0});
            const Vec3 secondTransverse=cross(
                torqueStart.secondProperDipole,Vec3{0.0,1.0,0.0});
            torqueProbe.firstDipoleTorque=firstTransverse
                *(0.25*torqueStart.firstProperDipole.norm()
                  /(firstTransverse.norm()
                    *std::abs(firstGyromagneticRatioOf())*torqueProbeStep));
            torqueProbe.secondDipoleTorque=secondTransverse
                *(0.25*torqueStart.secondProperDipole.norm()
                  /(secondTransverse.norm()
                    *std::abs(secondGyromagneticRatioOf())*torqueProbeStep));
            State quantizedTorqueState=torqueStart;
            State disabledTorqueState=torqueStart;
            applyDipoleRadiationTorqueForModel(quantizedTorqueState,torqueProbe,
                torqueProbeStep,
                ChargeRadiationReactionModel::stochasticElectricDipole);
            applyDipoleRadiationTorqueForModel(disabledTorqueState,torqueProbe,
                torqueProbeStep,ChargeRadiationReactionModel::disabled);
            const auto dipoleDirectionTravel=[&](const State& state) {
                return (state.firstProperDipole
                        /state.firstProperDipole.norm()
                       -torqueStart.firstProperDipole
                        /torqueStart.firstProperDipole.norm()).norm()
                    +(state.secondProperDipole
                        /state.secondProperDipole.norm()
                       -torqueStart.secondProperDipole
                        /torqueStart.secondProperDipole.norm()).norm();
            };
            quantizedDipoleTorqueTravel=
                dipoleDirectionTravel(quantizedTorqueState);
            disabledDipoleTorqueTravel=
                dipoleDirectionTravel(disabledTorqueState);
            quantizedDipoleTorqueNormDrift=std::max(
                std::abs(quantizedTorqueState.firstProperDipole.norm()
                         /torqueStart.firstProperDipole.norm()-1.0),
                std::abs(quantizedTorqueState.secondProperDipole.norm()
                         /torqueStart.secondProperDipole.norm()-1.0));
            const double quantizedPeriod=2.0*pi*std::sqrt(
                larmorReducedMass*larmorProbeRadius*larmorProbeRadius
                *larmorProbeRadius/pairCoulombStrength);
            bool quantizedAdvanced=true;
            for(int step=0;step<64&&quantizedAdvanced;++step)
                quantizedAdvanced=quantizedEngine.advance(
                    quantizedState,quantizedPeriod/256.0);
            if(quantizedAdvanced&&isFinite(quantizedState)) {
                const MutualForces quantizedForces=retardedExternalForces(
                    quantizedState,quantizedEngine.history());
                const ParticleMultipoleRadiation quantizedRadiation=
                    particleMultipoleRadiation(quantizedState,quantizedForces,
                        quantizedEngine.history(),true,
                        ChargeRadiationReactionModel::stochasticElectricDipole,
                        true);
                // Charge sector: zero continuous reaction force, normalized
                // against the external force it would have been added to.
                quantizedChargeReactionDrain=
                    (quantizedRadiation.chargeReaction.first.norm()
                    +quantizedRadiation.chargeReaction.second.norm())
                    /std::max(quantizedForces.first.norm()
                             +quantizedForces.second.norm(),1.0e-300);
                // Magnetic sector: the internal reservoir must not have been
                // drained, normalized against the energy the same trajectory
                // radiated through the (still continuous, diagnostic) flux.
                quantizedDipoleConstraintDrain=
                    std::abs(quantizedState.dipoleConstraintEnergy
                             -constraintStart)
                    /std::max(std::abs(quantizedState.radiatedEnergy),1.0e-300);
            }
        }
        balanceDiagnostics={{
            runBalance("retarded + no reaction",
                ChargeRadiationReactionModel::disabled,true),
            runBalance("Darwin + coherent reaction",
                ChargeRadiationReactionModel::coherentElectricDipole,false),
            runBalance("retarded + self-only LL",
                ChargeRadiationReactionModel::individualLandauLifshitzSelfOnly,true),
            runBalance("retarded + coherent reaction",
                ChargeRadiationReactionModel::coherentElectricDipole,true)}};
        const BalanceDiagnostic& productionBalance=balanceDiagnostics.back();
        if(productionBalance.advanced) {
            longHorizonBalanceResidual=
                std::abs(productionBalance.signedResidual);
            longHorizonRadiatedEnergy=productionBalance.fluxEnergy;
            longHorizonSchottChange=productionBalance.schottChange;
            // Conditioning of this probe, stated so the residual above can be
            // read at all.  The pair's binding energy at the probe radius is
            // the scale of conservativeParticleEnergy; the residual divides a
            // difference of two such numbers by the far-field energy, which
            // is four orders smaller.  Any relative error in the mechanical
            // energy -- the integrator's own secular drift included -- is
            // magnified by exactly this factor before it is reported.
            const double bindingScale=pairCoulombStrength/(2.0*balanceRadius);
            longHorizonBalanceAmplification=
                bindingScale/std::max(productionBalance.fluxEnergy,1.0e-300);
            // One refinement step.  A physical imbalance would sit still; a
            // discretization error shrinks.  A single period is enough
            // because the residual is period-independent (measured flat to
            // three digits over 1..32 periods), and halves the cost of the
            // tighter tolerance.
            balanceTolerance=1.0e-8;
            balancePeriodCount=1;
            const BalanceDiagnostic refined=runBalance(
                "retarded + coherent reaction, refined",
                ChargeRadiationReactionModel::coherentElectricDipole,true);
            balanceTolerance=1.0e-6;
            balancePeriodCount=balancePeriods;
            if(refined.advanced)
                longHorizonBalanceRefined=std::abs(refined.signedResidual);
        }
        // Radius sweep.  Three points spanning a factor of 16, each measured
        // at both discretizations.  R=0.2 is deliberately one of them, so the
        // headline numbers above appear inside the sweep too.
        {
            // Preserve the 16x dimensionless radius span, but do not place a
            // probe below its own declared spatial resolution.  For
            // protonium 0.05*a_pair is only 2.9 fm, below nuclearCutoff; the
            // resulting third-derivative subtraction is ill-conditioned and
            // refinement measures more cancellation noise, not convergence.
            // Other selectable pairs remain exactly at 0.05/0.2/0.8.
            const double minimumSweptFactor=std::max(
                0.05,1.25*nuclearCutoff/larmorProbeRadius);
            const std::array<double,3> sweptFactors{
                minimumSweptFactor,
                4.0*minimumSweptFactor,
                16.0*minimumSweptFactor};
            for(std::size_t index=0;index<sweptFactors.size();++index) {
                balanceRadiusFactor=sweptFactors[index];
                balanceTolerance=1.0e-6;
                balancePeriodCount=2;
                const BalanceDiagnostic coarse=runBalance("sweep coarse",
                    ChargeRadiationReactionModel::coherentElectricDipole,true);
                balanceTolerance=1.0e-8;
                balancePeriodCount=1;
                const BalanceDiagnostic fine=runBalance("sweep refined",
                    ChargeRadiationReactionModel::coherentElectricDipole,true);
                BalanceRadiusPoint& point=balanceRadiusSweep[index];
                point.factor=sweptFactors[index];
                point.advanced=coarse.advanced&&fine.advanced
                    &&coarse.fluxEnergy>0.0;
                if(!point.advanced) continue;
                point.amplification=
                    (pairCoulombStrength
                     /(2.0*sweptFactors[index]*larmorProbeRadius))
                    /coarse.fluxEnergy;
                point.coarse=std::abs(coarse.signedResidual);
                point.refined=std::abs(fine.signedResidual);
            }
            balanceRadiusFactor=0.2;
            balanceTolerance=1.0e-6;
            balancePeriodCount=balancePeriods;
        }
    }

    // 1% band around unity.  The measured deviation is 7e-5, so this leaves a
    // margin of more than a hundred and still fails outright on any break of
    // the absolute normalization -- which no other far-field check here can
    // see, since they all compare the quadrature with itself.
    // Separate from particle-covariance so a failure names the actual problem:
    // an argument routed to the wrong role, not a broken boost.
    const bool roleRoutingOk=
        std::isfinite(roleRoutingResidual)
        &&roleRoutingResidual<1.0e-12
        // Guards the guard: if the step stopped turning the moments the
        // residual would be trivially zero and the check meaningless.
        &&std::isfinite(roleRoutingTravel)
        &&roleRoutingTravel>1.0e-8;

    const bool larmorNormalizationOk=
        std::isfinite(larmorNormalizationRatio)
        &&std::isfinite(larmorProductionRatio)
        &&std::abs(larmorNormalizationRatio-1.0)<0.01
        &&std::abs(larmorProductionRatio-1.0)<0.01
        &&std::isfinite(larmorAccumulationRatio)
        // Looser than the instantaneous band: the accumulation is first order
        // in dt by construction.  Still far tighter than the 34% it used to be
        // out by.
        &&std::abs(larmorAccumulationRatio-1.0)<0.02;
    // What this gates, and what it deliberately does not.
    //
    // The residual is NOT a measurement of how well the model conserves
    // energy.  It divides a difference of two mechanical energies by a
    // far-field energy four orders smaller, so it magnifies any relative
    // error in the former by longHorizonBalanceAmplification -- about 1.4e4
    // at the probe's own radius.  At that gain the integrator's own secular
    // drift, at the 1e-6 tolerance it runs, accounts for the whole reading:
    // 4.3e-6 of relative mechanical-energy error is exactly the measured
    // 5.9%.
    //
    // Three independent measurements say the same thing.  The
    // Landau-Lifshitz self-force delivers 0.4997-0.4998 of the far-field
    // flux against a theoretical 0.5000 for this pair (a1 = -a2 makes the
    // coherent |q1 a1 + q2 a2|^2 exactly twice the sum of the individual
    // terms), so the reaction sector is right to three digits.  The far-field
    // energy itself is converged to 2e-5 across four decades of tolerance,
    // while the mechanical difference wanders by 12% over the same range.
    // And refining the validated second-order discretization moves the
    // residual from -5.9% at the production setting to about -1.0% at 1e-8,
    // which a real leak could not do.  (Raising the retarded-history
    // interpolation from cubic
    // to quintic Hermite, tried and reverted, moved it the WRONG way, to
    // -17%: the stored accelerations are only as mutually consistent with
    // the position/velocity samples as the second-order integrator makes
    // them, so a C^2 interpolant over-fits them.)
    //
    // So the gate below is a CONVERGENCE statement, not a conservation one:
    // the refined probe must be no worse than the coarse one, and must land
    // inside a band the coarse setting is not expected to reach.  The 10%
    // ceiling on the coarse value is kept as a blunt guard against the
    // probe breaking outright.
    const bool longHorizonBalanceOk=
        std::isfinite(longHorizonBalanceResidual)
        &&std::isfinite(longHorizonRadiatedEnergy)
        &&longHorizonRadiatedEnergy>0.0
        &&std::isfinite(longHorizonSchottChange)
        &&std::isfinite(longHorizonBalanceRefined)
        &&longHorizonBalanceResidual<0.1
        // Refinement must not make things worse -- but only where there is
        // something to improve.  Once both readings are at the probe's own
        // noise (p+e- measures 7.1e-5 coarse against 1.5e-4 refined, both
        // four orders inside the band below) which of them is larger is
        // decided by round-off, and demanding an ordering there fails a
        // perfectly healthy pair.  The floor makes the comparison bite only
        // when the coarse reading is genuinely above the noise.
        &&longHorizonBalanceRefined
            <=std::max(longHorizonBalanceResidual,1.0e-3)
        // The band that actually catches a physical imbalance.  Measured
        // against an injected 10% error in the radiation-reaction
        // coefficient: the residual does NOT converge away, going 10.9% ->
        // 6.0%, and lands here.  Discretization at the production setting
        // goes 5.9% -> 0.98% and clears it fivefold.  So this resolves a
        // reaction-sector error of roughly ten percent and up; it is not a
        // fine measurement of conservation, and nothing here should be read
        // as one.
        &&longHorizonBalanceRefined<0.05
        // Same two requirements across the sweep, so neither conclusion is an
        // accident of the single radius the headline is quoted from.
        &&std::ranges::all_of(balanceRadiusSweep,
            [](const BalanceRadiusPoint& point) {
                return point.advanced
                    &&std::isfinite(point.coarse)
                    &&std::isfinite(point.refined)
                    &&point.refined<=std::max(point.coarse,1.0e-3);
            })
        // And the measured amplification must follow R^(3/2).  E_mech goes as
        // 1/R and Larmor puts E_far per orbit at R^(-5/2), so this is a check
        // on the radiated energy's own scaling across a factor of 16 in
        // radius -- it fails on any break of the inspiral power law, which no
        // single-radius probe here can see.  Two percent is loose against the
        // 0.1% the ratios actually hold to, and tight against the factor of
        // 2.83 per doubling being tested.
        &&[&]{
            for(std::size_t index=1;index<balanceRadiusSweep.size();++index) {
                const BalanceRadiusPoint& previous=balanceRadiusSweep[index-1];
                const BalanceRadiusPoint& current=balanceRadiusSweep[index];
                if(!(previous.amplification>0.0)) return false;
                const double predicted=std::pow(
                    current.factor/previous.factor,1.5);
                const double measured=
                    current.amplification/previous.amplification;
                if(!(std::abs(measured/predicted-1.0)<0.02)) return false;
            }
            return true;
        }();

    const std::array<double,3> farControlRadii{
        1.0e4*bohrRadius,1.0e5*bohrRadius,1.0e6*bohrRadius};
    std::array<double,3> farRadiusResiduals{};
    std::array<double,3> farNearFieldContamination{};
    for(std::size_t index=0;index<farControlRadii.size();++index) {
        const FarFieldSampling radiationSampling{194,farControlRadii[index],true};
        const FieldFluxRates radiationOnly=electromagneticFieldFluxRates(
            sharedEngineVisualState,visualEngine.history(),radiationSampling);
        const FieldFluxRates fullField=electromagneticFieldFluxRates(
            sharedEngineVisualState,visualEngine.history(),
            {194,farControlRadii[index],false});
        const double angularScale=std::max({farFullReference.angularMomentum.norm(),
            std::abs(farFullReference.energy)*bohrRadius/c,1.0e-300});
        farRadiusResiduals[index]=std::max(
            radiativeFluxResidual(radiationOnly,farReference),
            (fullField.angularMomentum-farFullReference.angularMomentum).norm()
                /angularScale);
        farNearFieldContamination[index]=radiativeFluxResidual(
            fullField,radiationOnly);
    }

    const auto trajectoryResidual=[](const State& value,const State& reference) {
        // The old nuclearCutoff floor was an absolute hydrogen-scale length.
        // It dwarfed protonium and made a failed heavy-pair trajectory look
        // artificially accurate.  The active pair's own terminal surface is
        // the smallest meaningful length for this dimensionless comparison.
        const double lengthScale=std::max(
            separation(reference),collisionBoundaryRadius);
        const double speedScale=std::max(
            (reference.firstVelocity-reference.secondVelocity).norm(),
            1.0e-6*c);
        return std::max({
            (value.firstPosition-reference.firstPosition).norm()/lengthScale,
            (value.secondPosition-reference.secondPosition).norm()/lengthScale,
            (value.firstVelocity-reference.firstVelocity).norm()/speedScale,
            (value.secondVelocity-reference.secondVelocity).norm()/speedScale});
    };
    const auto integrateAccuracyCase=[&](double tolerance,double proposedDt,
                                         int steps) {
        State value=yeeCoupledState;
        ClassicalTrajectoryEngine engine(value,{tolerance,14});
        bool advanced=true;
        for(int step=0;step<steps;++step)
            advanced=engine.advance(value,proposedDt)&&advanced;
        if(!advanced) value.time=std::numeric_limits<double>::quiet_NaN();
        return value;
    };
    // Probe step, scaled to the PAIR's own dynamical time instead of nailed to
    // 1e-20 s.  A tolerance only tests the adaptive controller while the step
    // is long enough for it to bind, and a fixed step covers a different
    // fraction of an orbit for every pair.  For mu+mu- the muon is 207 times
    // heavier, 1e-20 s advanced it so little that tol=1e-5 and tol=1e-6 came
    // out bit-equal (residual 1.3374808e-07 twice), and the monotonicity below
    // failed on that tie -- not because the trajectory was inaccurate (it was
    // eighty times better than the 1e-5 bar) but because two of the three
    // tolerances had stopped testing anything.  Loosening the comparison to
    // <= would have hidden that instead of fixing it, and would equally admit
    // a controller that ignores the tolerance outright.
    //
    // The Kepler time sqrt(mu r^3 / k|q1 q2|) is the scale that binds.  Only
    // its RATIO to the e+e- value is needed, so the probe separation cancels
    // and the default pair divides identical operands: the ratio is exactly
    // 1.0 and every e+e- residual below is bit-for-bit what it was.
    constexpr double referenceReducedMass=
        reducedMassOf(ParticlePair{electron,positron});
    constexpr double referenceCoulombStrength=coulombConstant
        *magnitude(chargeProduct(ParticlePair{electron,positron}));
    const double trajectoryProbeDt=1.0e-20*std::sqrt(
        (pairReducedMass/referenceReducedMass)
       *(referenceCoulombStrength/pairCoulombStrength));
    const State trajectoryReference=
        integrateAccuracyCase(1.0e-8,trajectoryProbeDt,4);
    const std::array<double,3> trajectoryTolerances{1.0e-5,1.0e-6,1.0e-7};
    std::array<double,3> trajectoryToleranceResiduals{};
    for(std::size_t index=0;index<trajectoryTolerances.size();++index)
        trajectoryToleranceResiduals[index]=trajectoryResidual(
            integrateAccuracyCase(trajectoryTolerances[index],
                                  trajectoryProbeDt,4),
            trajectoryReference);
    const State trajectoryHalfStep=
        integrateAccuracyCase(1.0e-7,0.5*trajectoryProbeDt,8);
    const double trajectoryStepResidual=trajectoryResidual(
        integrateAccuracyCase(1.0e-7,trajectoryProbeDt,4),trajectoryHalfStep);

    // A deliberately under-resolved circular step exercises the recursion
    // ceiling itself.  maximumDepth=0 must reject it, while the same tolerance
    // with enough depth must resolve it.  Rejection is transactional: neither
    // the caller's state nor the retained causal history may advance by even a
    // half-step.
    State depthLimitInitial;
    const double depthLimitRadius=pairBohrRadius(activePair);
    const double depthLimitFirstShare=secondMass/(firstMass+secondMass);
    const double depthLimitSecondShare=firstMass/(firstMass+secondMass);
    const double depthLimitRelativeSpeed=std::sqrt(
        pairCoulombStrength/(pairReducedMass*depthLimitRadius));
    depthLimitInitial.firstPosition={
        depthLimitFirstShare*depthLimitRadius,0.0,0.0};
    depthLimitInitial.secondPosition={
        -depthLimitSecondShare*depthLimitRadius,0.0,0.0};
    depthLimitInitial.firstVelocity={
        0.0,depthLimitFirstShare*depthLimitRelativeSpeed,0.0};
    depthLimitInitial.secondVelocity={
        0.0,-depthLimitSecondShare*depthLimitRelativeSpeed,0.0};
    synchronizeCovariantDipoles(depthLimitInitial);
    const double depthLimitPeriod=2.0*pi*std::sqrt(
        pairReducedMass*depthLimitRadius*depthLimitRadius*depthLimitRadius
        /pairCoulombStrength);
    const double depthLimitDt=depthLimitPeriod/8.0;
    const StateHistory depthLimitHistory=causalInitialHistory(depthLimitInitial);
    State depthRejectedState=depthLimitInitial;
    ClassicalTrajectoryEngine depthLimitedEngine(depthLimitHistory,
        {.relativeTolerance=1.0e-6,.maximumDepth=0,
         .reactionModel=ChargeRadiationReactionModel::disabled,
         .computeOutwardFlux=false});
    const std::size_t depthHistorySizeBefore=depthLimitedEngine.history().size();
    const double depthHistoryStartBefore=depthLimitedEngine.history().front().time;
    const double depthHistoryEndBefore=depthLimitedEngine.history().back().time;
    const bool depthLimitAccepted=
        depthLimitedEngine.advance(depthRejectedState,depthLimitDt);
    const double depthRollbackResidual=std::max(
        trajectoryResidual(depthRejectedState,depthLimitInitial),
        std::abs(depthRejectedState.time-depthLimitInitial.time)/depthLimitDt);
    const bool depthHistoryUnchanged=
        depthLimitedEngine.history().size()==depthHistorySizeBefore
        &&depthLimitedEngine.history().front().time==depthHistoryStartBefore
        &&depthLimitedEngine.history().back().time==depthHistoryEndBefore;
    State depthResolvedState=depthLimitInitial;
    ClassicalTrajectoryEngine depthResolvingEngine(depthLimitHistory,
        {.relativeTolerance=1.0e-6,.maximumDepth=12,
         .reactionModel=ChargeRadiationReactionModel::disabled,
         .computeOutwardFlux=false});
    const bool depthResolved=
        depthResolvingEngine.advance(depthResolvedState,depthLimitDt);
    const double depthResolvedTimeResidual=std::abs(
        depthResolvedState.time-depthLimitInitial.time-depthLimitDt)
        /depthLimitDt;

    struct HistorySensitivityCase { State early,settled; bool advanced=false; };
    const auto integrateHistoryCase=[&](double spanFactor,int intervals,
                                        int picardIterations) {
        State value=yeeCoupledState;
        ClassicalTrajectoryEngine engine(
            causalInitialHistory(value,spanFactor,intervals,picardIterations),
            {.relativeTolerance=1.0e-7,.maximumDepth=14,
             .reactionModel=ChargeRadiationReactionModel::disabled,
             .computeOutwardFlux=false});
        bool advanced=true;
        advanced=engine.advance(value,trajectoryProbeDt)&&advanced;
        const State early=value;
        // By eight probe steps the retained window has been replaced by
        // accepted positive-time states.  Any remaining difference is a
        // dynamical memory of the preparation, not direct interpolation of
        // different hidden-past nodes.
        for(int step=1;step<8&&advanced;++step)
            advanced=engine.advance(value,trajectoryProbeDt)&&advanced;
        if(!advanced) value.time=std::numeric_limits<double>::quiet_NaN();
        return HistorySensitivityCase{early,value,advanced};
    };
    // Deliberately vary three independent construction choices together:
    // retained span, node density and self-consistency iteration count.  The
    // coarse history uses only the instantaneous-acceleration seed (0 Picard
    // updates); the reference uses four updates, twice the production count.
    const HistorySensitivityCase startupCase4=
        integrateHistoryCase(4.0,32,0);
    const HistorySensitivityCase startupCase8=
        integrateHistoryCase(8.0,64,2);
    const HistorySensitivityCase startupCase16=
        integrateHistoryCase(16.0,128,4);
    const State& startupHistory4=startupCase4.settled;
    const State& startupHistory8=startupCase8.settled;
    const State& startupHistory16=startupCase16.settled;
    const double startupHistoryResidual4=trajectoryResidual(
        startupHistory4,startupHistory16);
    const double startupHistoryResidual8=trajectoryResidual(
        startupHistory8,startupHistory16);
    const double startupEarlySensitivity=std::max(
        trajectoryResidual(startupCase4.early,startupCase16.early),
        trajectoryResidual(startupCase8.early,startupCase16.early));
    const double startupSettledSensitivity=std::max(
        startupHistoryResidual4,startupHistoryResidual8);

    // Parameter-sensitivity scans use actual trajectories rather than merely
    // evaluating the cutoff/profile formulas.  The cutoff probe is a common
    // radial plunge stopped on three terminal surfaces.  The regulator probe
    // advances one common short-range dipolar state with three smoothing
    // radii, restoring the run-global radius immediately afterwards.
    const std::array<double,3> cutoffFractions{0.004,0.005,0.006};
    std::array<double,3> cutoffArrivalTimes{};
    std::array<double,3> cutoffEventResiduals{};
    const double sensitivityOrbitRadius=pairBohrRadius(activePair);
    const double sensitivityStartRadius=0.02*sensitivityOrbitRadius;
    const double sensitivityTimeScale=std::sqrt(
        pairReducedMass*sensitivityOrbitRadius*sensitivityOrbitRadius
            *sensitivityOrbitRadius/pairCoulombStrength);
    const double sensitivityFirstShare=secondMass/(firstMass+secondMass);
    const double sensitivitySecondShare=firstMass/(firstMass+secondMass);
    const double sensitivityRadialSpeed=0.25*std::sqrt(
        pairCoulombStrength/(pairReducedMass*sensitivityStartRadius));
    State cutoffInitial;
    cutoffInitial.firstPosition={sensitivityFirstShare*sensitivityStartRadius,0,0};
    cutoffInitial.secondPosition={-sensitivitySecondShare*sensitivityStartRadius,0,0};
    cutoffInitial.firstVelocity={-sensitivityFirstShare*sensitivityRadialSpeed,0,0};
    cutoffInitial.secondVelocity={sensitivitySecondShare*sensitivityRadialSpeed,0,0};
    // This probe validates event localization at the cutoff.  Carrying the
    // physical dipoles here mixed in a second question: for protonium their
    // short-range barrier correctly prevents the radial plunge from ever
    // reaching any of the three surfaces.  The regulator has its own probe
    // below, so switch dipoles off and test one mechanism at a time.
    cutoffInitial.firstDipole={};
    cutoffInitial.secondDipole={};
    synchronizeCovariantDipoles(cutoffInitial);
    for(std::size_t index=0;index<cutoffFractions.size();++index) {
        const double cutoff=cutoffFractions[index]*sensitivityOrbitRadius;
        State value=cutoffInitial;
        ClassicalTrajectoryEngine engine(value,
            {.relativeTolerance=1.0e-6,.maximumDepth=12,
             .reactionModel=ChargeRadiationReactionModel::disabled,
             .computeOutwardFlux=false});
        bool advanced=true;
        State before=value;
        for(int step=0;step<4096&&separation(value)>cutoff&&advanced;++step) {
            before=value;
            const double radius=separation(value);
            const double omega=std::sqrt(
                pairCoulombStrength/(pairReducedMass*radius*radius*radius));
            advanced=engine.advance(value,2.0*pi/(256.0*omega));
        }
        if(advanced&&separation(value)<=cutoff) {
            const double fraction=separationCrossingFraction(before,value,cutoff);
            cutoffArrivalTimes[index]=(before.time
                +fraction*(value.time-before.time))/sensitivityTimeScale;
            const State event=interpolateState(before,value,fraction);
            cutoffEventResiduals[index]=std::abs(separation(event)-cutoff)/cutoff;
        } else {
            cutoffArrivalTimes[index]=std::numeric_limits<double>::quiet_NaN();
            cutoffEventResiduals[index]=std::numeric_limits<double>::infinity();
        }
    }

    const std::array<double,3> regularizationRadiusFactors{0.5,1.0,2.0};
    std::array<State,3> regularizationScanStates{};
    std::array<bool,3> regularizationScanAdvanced{};
    const double savedRegularizationRadius=magneticRegularizationRadius;
    const double regulatorProbeSeparation=4.0*savedRegularizationRadius;
    const double regulatorProbeSpeed=std::sqrt(
        pairCoulombStrength/(pairReducedMass*regulatorProbeSeparation));
    const double regulatorProbePeriod=2.0*pi*std::sqrt(
        pairReducedMass*regulatorProbeSeparation*regulatorProbeSeparation
            *regulatorProbeSeparation/pairCoulombStrength);
    State regulatorInitial;
    regulatorInitial.firstPosition={
        sensitivityFirstShare*regulatorProbeSeparation,0,0};
    regulatorInitial.secondPosition={
        -sensitivitySecondShare*regulatorProbeSeparation,0,0};
    regulatorInitial.firstVelocity={0,sensitivityFirstShare*regulatorProbeSpeed,0};
    regulatorInitial.secondVelocity={0,-sensitivitySecondShare*regulatorProbeSpeed,0};
    // Proper moments are the independent representation.  Copying them into
    // the lab slots as well is inconsistent once the particles move: the
    // missing boosted electric dipoles produce a finite discontinuity on the
    // first engine step that adaptive subdivision cannot converge away.
    regulatorInitial.firstProperDipole={0,0,firstMagneticMoment};
    regulatorInitial.secondProperDipole={0,0,secondMagneticMoment};
    synchronizeCovariantDipoles(regulatorInitial);
    for(std::size_t index=0;index<regularizationRadiusFactors.size();++index) {
        magneticRegularizationRadius=
            savedRegularizationRadius*regularizationRadiusFactors[index];
        State value=regulatorInitial;
        ClassicalTrajectoryEngine engine(value,
            {.relativeTolerance=1.0e-6,.maximumDepth=14,
             .reactionModel=ChargeRadiationReactionModel::disabled,
             .computeOutwardFlux=false});
        bool advanced=true;
        // This is a local sensitivity probe, not a stability claim for an
        // orbit whose regulator has deliberately been changed by 2x.  A
        // 1/256-orbit window is long enough to separate all three responses
        // from round-off while staying before the intentionally modified
        // short-range potential can turn the comparison into a collision or
        // escape experiment of its own.
        for(int step=0;step<8&&advanced;++step)
            advanced=engine.advance(value,regulatorProbePeriod/2048.0);
        regularizationScanAdvanced[index]=advanced;
        if(!advanced) value.time=std::numeric_limits<double>::quiet_NaN();
        regularizationScanStates[index]=value;
    }
    magneticRegularizationRadius=savedRegularizationRadius;
    const std::array<double,3> regularizationTrajectoryResiduals{
        trajectoryResidual(regularizationScanStates[0],regularizationScanStates[1]),
        0.0,
        trajectoryResidual(regularizationScanStates[2],regularizationScanStates[1])};
    const StateHistory startupHistory=causalInitialHistory(yeeCoupledState);
    const double startupCoverage=(yeeCoupledState.time-startupHistory.front().time)
        /(separation(yeeCoupledState)/c);
    const MutualForces startupForces=retardedExternalForces(
        yeeCoupledState,startupHistory);
    const bool startupFinite=isFinite(startupForces.first)
        &&isFinite(startupForces.second)
        &&isFinite(startupHistory.back());

    const auto interpolationErrors=[](double span) {
        constexpr double orbitRadius=bohrRadius;
        constexpr double angularRate=2.0e17;
        const auto circularState=[](double time) {
            constexpr double radius=bohrRadius;
            constexpr double omega=2.0e17;
            const double phase=omega*time;
            State state;
            state.time=time;
            state.firstPosition={radius*std::cos(phase),
                                    radius*std::sin(phase),0};
            state.firstVelocity={-radius*omega*std::sin(phase),
                                     radius*omega*std::cos(phase),0};
            state.firstAcceleration=state.firstPosition*(-omega*omega);
            return state;
        };
        const State older=circularState(0.0);
        const State newer=circularState(span);
        const double query=0.5*span;
        const State exact=circularState(query);
        const ChargeKinematics linear=linearlyInterpolatedCharge(
            older,newer,true,query);
        const ChargeKinematics hermite=interpolatedCharge(
            older,newer,true,query);
        const auto residual=[&](const ChargeKinematics& value) {
            return std::max({
                (value.position-exact.firstPosition).norm()/orbitRadius,
                (value.velocity-exact.firstVelocity).norm()
                    /(orbitRadius*angularRate),
                (value.acceleration-exact.firstAcceleration).norm()
                    /(orbitRadius*angularRate*angularRate)});
        };
        return std::array<double,2>{residual(linear),residual(hermite)};
    };
    const std::array<double,2> interpolationCoarse=
        interpolationErrors(1.0e-18);
    const std::array<double,2> interpolationFine=
        interpolationErrors(0.5e-18);
    const double hermiteConvergenceOrder=std::log(
        interpolationCoarse[1]/interpolationFine[1])/std::log(2.0);

    // Audit follow-up (README's "audyt fizyki"): the charge-charge field
    // (lienardWiechertField) interpolates its retarded source through
    // interpolatedCharge()'s Hermite stencil, just measured above, but the
    // DIPOLE field path (retardedElectricDipoleField/
    // retardedMagneticDipoleField, via historicalDipoleKinematics ->
    // historicalSource) interpolates position/velocity with a plain LINEAR
    // blend (historicalSource's own interpolateVector calls, no derivative
    // matching).  No existing check measured this path's own geometric
    // accuracy on a curved worldline -- only the derivative STENCILS built
    // on top of it, at exact history nodes.  Reusing the identical circular
    // reference trajectory and span as interpolationErrors above, so the
    // two are directly comparable number for number.
    const auto sourceInterpolationErrors=[](double span) {
        constexpr double orbitRadius=bohrRadius;
        constexpr double angularRate=2.0e17;
        const auto circularState=[](double time) {
            constexpr double radius=bohrRadius;
            constexpr double omega=2.0e17;
            const double phase=omega*time;
            State state;
            state.time=time;
            state.firstPosition={radius*std::cos(phase),
                                    radius*std::sin(phase),0};
            state.firstVelocity={-radius*omega*std::sin(phase),
                                     radius*omega*std::cos(phase),0};
            state.firstAcceleration=state.firstPosition*(-omega*omega);
            return state;
        };
        const State older=circularState(0.0);
        const State newer=circularState(span);
        const double query=0.5*span;
        const State exact=circularState(query);
        const StateHistory sourceHistory{older,newer};
        const RetardedSourceSample linear=historicalSource(
            sourceHistory,newer,true,query);
        return std::max(
            (linear.position-exact.firstPosition).norm()/orbitRadius,
            (linear.velocity-exact.firstVelocity).norm()
                /(orbitRadius*angularRate));
    };
    const double sourceInterpolationCoarse=
        sourceInterpolationErrors(1.0e-18);
    const double sourceInterpolationFine=
        sourceInterpolationErrors(0.5e-18);
    const double sourceInterpolationOrder=std::log(
        sourceInterpolationCoarse/sourceInterpolationFine)/std::log(2.0);

    State staticDipoleState;
    staticDipoleState.firstPosition={2.0*nuclearCutoff,0,0};
    staticDipoleState.secondPosition={0,0,0};
    staticDipoleState.secondDipole={0,0,secondMagneticMoment};
    State staticDipolePast=staticDipoleState;
    staticDipolePast.time=-8.0*nuclearCutoff/c;
    const StateHistory staticDipoleHistory{staticDipolePast,staticDipoleState};
    const ElectromagneticField retardedStaticDipole=
        retardedMagneticDipoleField(staticDipoleState.firstPosition,0.0,
            staticDipoleHistory,staticDipoleState,false);
    // retardedMagneticDipoleField now softens its distance below
    // separationFloor() (the Compton barrier for e+e-, see positronium.cpp:
    // clampedSeparationVector, r_eff = sqrt(r^2+floor^2)), so its static
    // limit no longer matches the RAW regularizedDipoleField at a probe
    // point this deep (2*nuclearCutoff = 20 fm, well under the 193.3 fm
    // floor) -- it matches the equally-softened one, which is what the
    // reference below must compute for this to still test the intended
    // identity.
    const Vec3 directStaticDipole=regularizedDipoleField(
        clampedSeparationVector(staticDipoleState.firstPosition
            -staticDipoleState.secondPosition,separationFloor()),
        staticDipoleState.secondDipole);
    const double staticDipoleScale=std::max(directStaticDipole.norm(),1.0e-300);
    const double retardedStaticLimitResidual=
        (retardedStaticDipole.magnetic-directStaticDipole).norm()
        /staticDipoleScale;
    const Vec3 equatorialPointField=regularizedDipoleField(
        {nuclearCutoff,0,0},staticDipoleState.secondDipole,0.0);
    const std::array<double,3> regulatorRadii{
        0.75*nuclearCutoff,0.50*nuclearCutoff,0.25*nuclearCutoff};
    std::array<double,3> regulatorFieldResiduals{};
    for(std::size_t index=0;index<regulatorRadii.size();++index) {
        const Vec3 field=regularizedDipoleField({nuclearCutoff,0,0},
            staticDipoleState.secondDipole,regulatorRadii[index]);
        regulatorFieldResiduals[index]=(field-equatorialPointField).norm()
            /equatorialPointField.norm();
    }
    // Probe the regulator on its own terms instead of at a fixed absolute
    // radius.  The previous version evaluated it at nuclearCutoff and required
    // the residual to stay small, which only held while the smoothing radius
    // was deliberately kept below the reported boundary; it therefore encoded
    // that choice rather than testing the operator.  What the regulator must
    // actually do is (a) converge to the point dipole well outside a, faster
    // for a steeper exponent, and (b) stay bounded as r goes to zero.
    const std::array<double,3> regulatorExponents{4.0,6.0,8.0};
    const double regulatorFarRadius=10.0*magneticRegularizationRadius;
    const Vec3 farPointField=regularizedDipoleField(
        {regulatorFarRadius,0,0},staticDipoleState.secondDipole,0.0);
    std::array<double,3> regulatorProfileResiduals{};
    for(std::size_t index=0;index<regulatorExponents.size();++index) {
        const Vec3 field=regularizedDipoleField({regulatorFarRadius,0,0},
            staticDipoleState.secondDipole,magneticRegularizationRadius,
            regulatorExponents[index]);
        regulatorProfileResiduals[index]=(field-farPointField).norm()
            /farPointField.norm();
    }
    // Outside the smoothing core the regularized field must agree with the
    // point dipole in direction.  Inside the core it legitimately reverses,
    // exactly as the field inside a uniformly magnetized sphere is parallel to
    // M rather than antiparallel as on the equator of a point dipole, so the
    // direction check has to be made where the two are supposed to agree.
    const Vec3 farRegularizedField=regularizedDipoleField(
        {regulatorFarRadius,0,0},staticDipoleState.secondDipole);
    const double regulatorFarAlignment=dot(farRegularizedField,farPointField);
    // Deep inside the core the regularized field must be strongly suppressed
    // relative to the point dipole, which is what removes the singularity.
    const double regulatorCoreRadius=0.01*magneticRegularizationRadius;
    const double regulatorCoreSuppression=
        regularizedDipoleField({regulatorCoreRadius,0,0},
            staticDipoleState.secondDipole).norm()
        /regularizedDipoleField({regulatorCoreRadius,0,0},
            staticDipoleState.secondDipole,0.0).norm();
    const Vec3 regulatorOriginField=regularizedDipoleField(
        {},staticDipoleState.secondDipole);
    const Vec3 regulatorOriginForce=regularizedDipoleForce(
        {},{firstMagneticMoment,0,0},staticDipoleState.secondDipole);
    // The physical ceiling the radius is chosen for, measured through the same
    // curl(A) field as production.  With the separation along x the field map
    // is diagonal: radial and transverse source/target dipoles span all of its
    // singular values, hence the larger absolute energy is the maximum over
    // every orientation.  The old test sampled only w(r)*mu1*mu2/r^3; it
    // certified 0.5 while this production operator actually reached 2.778.
    const double lighterMass=std::min(firstMass,secondMass);
    double peakDipoleEnergyOverRestEnergy=0.0;
    double peakRadialDipoleEnergyOverRestEnergy=0.0;
    double peakTransverseDipoleEnergyOverRestEnergy=0.0;
    double measuredDipoleCurlPeak=0.0;
    // 16000 samples, not 4000.  The peak of the regularized energy narrows in
    // log-radius as the regulator exponent rises, and at the production
    // exponent 12 a 4000-point grid over six decades (0.35% spacing in radius
    // at the peak) resolved the curl peak only to 2.778064 against the derived
    // 2.7783644, a 3.0e-4 miss against this check's 1e-4 tolerance.  The error
    // is quadratic in the spacing, so four times the samples brings it to
    // about 2e-5.  The constant itself is not in doubt: evaluated at the
    // ANALYTIC peak radius below, the energy matches the 0.5 ceiling to
    // 1.7e-16.
    for(int sample=0;sample<16000;++sample) {
        const double radius=magneticRegularizationRadius
            *std::pow(10.0,-3.0+6.0*sample/15999.0);
        const Vec3 separation{radius,0,0};
        const Vec3 radialTarget{firstMagneticMoment,0,0};
        const Vec3 transverseTarget{0,firstMagneticMoment,0};
        const double radialEnergy=std::abs(
            regularizedDipoleInteractionEnergy(separation,radialTarget,
                {secondMagneticMoment,0,0}));
        const double transverseEnergy=std::abs(
            regularizedDipoleInteractionEnergy(separation,transverseTarget,
                {0,secondMagneticMoment,0}));
        const double energy=std::max(radialEnergy,transverseEnergy);
        peakRadialDipoleEnergyOverRestEnergy=std::max(
            peakRadialDipoleEnergyOverRestEnergy,
            radialEnergy/(lighterMass*c*c));
        peakTransverseDipoleEnergyOverRestEnergy=std::max(
            peakTransverseDipoleEnergyOverRestEnergy,
            transverseEnergy/(lighterMass*c*c));
        peakDipoleEnergyOverRestEnergy=std::max(
            peakDipoleEnergyOverRestEnergy,energy/(lighterMass*c*c));
        measuredDipoleCurlPeak=std::max(measuredDipoleCurlPeak,
            energy*std::pow(magneticRegularizationRadius,3)
                /((mu0/(4.0*pi))*firstMagneticMoment*secondMagneticMoment));
    }
    // Analytic peak of the transverse coefficient w[n(1-w)-1]/r^3.  With
    // u = (r_reg/r)^n the extremum solves (n-1)(3-n)u^2+(n^2+4n-6)u-3 = 0;
    // at the production exponent 12 that is 33u^2-62u+1 = 0, u = (31+4sqrt58)/33,
    // so the peak sits at r_reg*[33/(31+4sqrt58)]^(1/12) = 0.9494927 r_reg.
    // (At the former exponent 6 the same equation gave (9-2sqrt19)^(1/6).)
    const double exactDipolePeakRadius=magneticRegularizationRadius
        *std::pow(33.0/(31.0+4.0*std::sqrt(58.0)),1.0/12.0);
    const Vec3 exactDipolePeakSeparation{exactDipolePeakRadius,0,0};
    const double exactParallelDipoleEnergy=regularizedDipoleInteractionEnergy(
        exactDipolePeakSeparation,{0,firstMagneticMoment,0},
        {0,secondMagneticMoment,0});
    const double exactAntiparallelDipoleEnergy=
        regularizedDipoleInteractionEnergy(exactDipolePeakSeparation,
            {0,-firstMagneticMoment,0},{0,secondMagneticMoment,0});
    const double exactDipolePeakOverRestEnergy=
        exactAntiparallelDipoleEnergy/(lighterMass*c*c);
    // Both moments come from the roles now.  This mixed a bare magneton for
    // the first dipole with the second role's actual moment, which was a
    // half-finished migration rather than a deliberate pairing; the check is
    // only isFinite, so nothing downstream moved either way.
    const Vec3 cutoffDipoleForce=regularizedDipoleForce(
        {nuclearCutoff,0,0},{firstMagneticMoment,0,0},
        staticDipoleState.secondDipole);

    // State integrity is intentionally tested independently of a physical
    // trajectory.  Distinct, finite values in every member make omissions in
    // interpolateState() visible instead of letting default zeroes pass.
    const std::array stateVectorMembers{
        &State::firstPosition,&State::secondPosition,
        &State::firstVelocity,&State::secondVelocity,
        &State::firstAcceleration,&State::secondAcceleration,
        &State::firstDipole,&State::secondDipole,
        &State::firstElectricDipole,&State::secondElectricDipole,
        &State::firstProperDipole,&State::secondProperDipole,
        &State::radiatedMomentum,&State::radiatedAngularMomentum,
        &State::boundFieldMomentum,&State::boundFieldAngularMomentum,
        &State::previousFluxMomentum,&State::previousFluxAngularMomentum,
        &State::previousMismatchMomentum,&State::previousMismatchAngularMomentum,
        &State::reactionMomentumMismatch,&State::reactionAngularMomentumMismatch
    };
    const std::array stateScalarMembers{
        &State::time,&State::radiatedEnergy,&State::orbitalRadiatedEnergy,
        &State::dipoleConstraintEnergy,&State::zeroPointPhase,
        &State::boundFieldEnergy,&State::reactionEnergyMismatch,
        &State::previousStepDt,&State::previousFluxEnergy,
        &State::previousDipoleFluxEnergy,&State::previousMismatchEnergy
    };
    State interpolationBefore,interpolationAfter;
    double stateTestValue=1.0;
    for(const auto member:stateVectorMembers) {
        interpolationBefore.*member={stateTestValue,stateTestValue+0.25,
                                     stateTestValue+0.5};
        interpolationAfter.*member={stateTestValue+20.0,stateTestValue+20.5,
                                    stateTestValue+21.0};
        stateTestValue+=1.0;
    }
    for(const auto member:stateScalarMembers) {
        interpolationBefore.*member=stateTestValue;
        interpolationAfter.*member=stateTestValue+20.0;
        stateTestValue+=1.0;
    }
    interpolationBefore.hasPreviousRates=true;
    interpolationAfter.hasPreviousRates=false;
    const auto sameVector=[](const Vec3& first,const Vec3& second) {
        return first.x==second.x&&first.y==second.y&&first.z==second.z;
    };
    const auto sameState=[&](const State& first,const State& second) {
        return std::ranges::all_of(stateVectorMembers,[&](const auto member) {
                return sameVector(first.*member,second.*member);
            })
            &&std::ranges::all_of(stateScalarMembers,[&](const auto member) {
                return first.*member==second.*member;
            })
            &&first.hasPreviousRates==second.hasPreviousRates;
    };
    const std::array linearStateVectorMembers{
        &State::firstPosition,&State::secondPosition,
        &State::firstVelocity,&State::secondVelocity,
        &State::firstAcceleration,&State::secondAcceleration,
        &State::firstElectricDipole,&State::secondElectricDipole,
        &State::radiatedMomentum,&State::radiatedAngularMomentum,
        &State::boundFieldMomentum,&State::boundFieldAngularMomentum,
        &State::reactionMomentumMismatch,&State::reactionAngularMomentumMismatch
    };
    const std::array dipoleStateMembers{
        &State::firstDipole,&State::secondDipole,
        &State::firstProperDipole,&State::secondProperDipole
    };
    const std::array continuousStateScalarMembers{
        &State::time,&State::radiatedEnergy,&State::orbitalRadiatedEnergy,
        &State::dipoleConstraintEnergy,&State::zeroPointPhase,
        &State::boundFieldEnergy,&State::reactionEnergyMismatch
    };
    constexpr double stateInterpolationFraction=0.375;
    const State interpolatedState=interpolateState(
        interpolationBefore,interpolationAfter,stateInterpolationFraction);
    const bool stateInterpolationOk=
        sameState(interpolateState(interpolationBefore,interpolationAfter,0.0),
                  interpolationBefore)
        &&sameState(interpolateState(interpolationBefore,interpolationAfter,1.0),
                   interpolationAfter)
        &&std::ranges::all_of(linearStateVectorMembers,[&](const auto member) {
            return sameVector(interpolatedState.*member,
                interpolateVector(interpolationBefore.*member,
                                  interpolationAfter.*member,
                                  stateInterpolationFraction));
        })
        &&std::ranges::all_of(dipoleStateMembers,[&](const auto member) {
            return sameVector(interpolatedState.*member,
                interpolateDipole(interpolationBefore.*member,
                                  interpolationAfter.*member,
                                  stateInterpolationFraction));
        })
        &&std::ranges::all_of(continuousStateScalarMembers,
            [&](const auto member) {
                return interpolatedState.*member==interpolationBefore.*member
                    +(interpolationAfter.*member-interpolationBefore.*member)
                        *stateInterpolationFraction;
            })
        &&!interpolatedState.hasPreviousRates
        &&interpolatedState.previousStepDt==0.0
        &&interpolatedState.previousFluxEnergy==0.0
        &&interpolatedState.previousDipoleFluxEnergy==0.0
        &&interpolatedState.previousFluxMomentum.squaredNorm()==0.0
        &&interpolatedState.previousFluxAngularMomentum.squaredNorm()==0.0
        &&interpolatedState.previousMismatchEnergy==0.0
        &&interpolatedState.previousMismatchMomentum.squaredNorm()==0.0
        &&interpolatedState.previousMismatchAngularMomentum.squaredNorm()==0.0
        &&isFinite(interpolatedState);
    const std::array finiteStateScalarMembers{
        &State::orbitalRadiatedEnergy,&State::previousStepDt,
        &State::previousFluxEnergy,&State::previousDipoleFluxEnergy,
        &State::previousMismatchEnergy
    };
    const std::array finiteStateVectorMembers{
        &State::previousFluxMomentum,&State::previousFluxAngularMomentum,
        &State::previousMismatchMomentum,&State::previousMismatchAngularMomentum
    };
    const std::array invalidStateValues{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()
    };
    const bool stateFiniteGuardOk=isFinite(interpolationBefore)
        &&isFinite(interpolationAfter)
        &&std::ranges::all_of(invalidStateValues,[&](double invalidValue) {
            return std::ranges::all_of(finiteStateScalarMembers,
                [&](const auto member) {
                    State invalid=interpolationBefore;
                    invalid.*member=invalidValue;
                    return !isFinite(invalid);
                })
                &&std::ranges::all_of(finiteStateVectorMembers,
                    [&](const auto member) {
                        State invalid=interpolationBefore;
                        invalid.*member={invalidValue,0.0,0.0};
                        return !isFinite(invalid);
                    });
        });
    const bool stateIntegrityOk=stateInterpolationOk&&stateFiniteGuardOk;

    // Synthetic censoring regression with known directions: failures become
    // more likely as energy rises and impact parameter falls, while the
    // administrative gate becomes more likely at high energy and large b.
    // Repeated deterministic strata avoid making the validator itself
    // stochastic and keep the logistic fits away from complete separation.
    std::vector<positronium::statistics::CensoringObservation>
        syntheticCensoringSample;
    for(int energyIndex=0;energyIndex<10;++energyIndex) {
        for(int impactIndex=0;impactIndex<10;++impactIndex) {
            const double failureProbability=0.02+0.02*energyIndex
                +0.015*(9-impactIndex);
            const double administrativeProbability=0.03+0.01*energyIndex
                +0.02*impactIndex;
            for(int replicate=0;replicate<20;++replicate) {
                const double quantile=(replicate+0.5)/20.0;
                auto disposition=
                    positronium::statistics::ObservationDisposition::Observed;
                if(quantile<failureProbability) {
                    disposition=positronium::statistics::
                        ObservationDisposition::NumericalFailure;
                } else if(quantile
                          <failureProbability+administrativeProbability) {
                    disposition=positronium::statistics::
                        ObservationDisposition::AdministrativelyCensored;
                }
                syntheticCensoringSample.push_back({
                    0.5+0.2*energyIndex,2.0+0.4*impactIndex,disposition,
                    (energyIndex+impactIndex+replicate)%3});
            }
        }
    }
    const auto syntheticCensoring=positronium::statistics::analyzeCensoring(
        syntheticCensoringSample,3);
    const auto syntheticFailureByEnergy=positronium::statistics::binnedRate(
        syntheticCensoringSample,true,
        positronium::statistics::BinaryEndpoint::NumericalFailure);
    const double syntheticCategorySum=std::accumulate(
        syntheticCensoring.ipcwCategoryProbability.begin(),
        syntheticCensoring.ipcwCategoryProbability.end(),0.0);
    std::vector<positronium::statistics::CensoringObservation>
        uncensoredSyntheticSample;
    for(int index=0;index<30;++index) {
        uncensoredSyntheticSample.push_back({1.0+0.1*index,0.5+0.2*index,
            positronium::statistics::ObservationDisposition::Observed,index%3});
    }
    const auto uncensoredSynthetic=positronium::statistics::analyzeCensoring(
        uncensoredSyntheticSample,3);
    const bool censoringModelOk=syntheticCensoring.validCount==2000
        &&syntheticCensoring.observedCount>0
        &&syntheticCensoring.administrativelyCensoredCount>0
        &&syntheticCensoring.numericalFailureCount>0
        &&syntheticCensoring.completionModel.fitted
        &&syntheticCensoring.completionModel.converged
        &&syntheticCensoring.failureModel.fitted
        &&syntheticCensoring.failureModel.converged
        &&syntheticCensoring.failureModel.coefficient[1]>0.0
        &&syntheticCensoring.failureModel.coefficient[2]<0.0
        &&syntheticCensoring.completionModel.coefficient[1]<0.0
        &&syntheticCensoring.completionModel.coefficient[2]<0.0
        &&std::abs(syntheticCategorySum-1.0)<1.0e-12
        &&syntheticCensoring.effectiveObservedSampleSize>0.0
        &&!syntheticFailureByEnergy.empty()
        &&uncensoredSynthetic.completionModel.constant
        &&uncensoredSynthetic.completionModel.constantProbability==1.0
        &&uncensoredSynthetic.failureModel.constant
        &&uncensoredSynthetic.failureModel.constantProbability==0.0
        &&uncensoredSynthetic.maximumIpcwWeight==1.0
        &&uncensoredSynthetic.effectiveObservedSampleSize==30.0;
    SimulationOptions visualTimeLimitOptions;
    visualTimeLimitOptions.collectFrames=false;
    visualTimeLimitOptions.observationTime=1.0e-24;
    visualTimeLimitOptions.radiatedEnergyBookkeeping=false;
    const SimulationResult visualTimeLimited=simulate(
        0x19a72ULL,3,visualTimeLimitOptions);
    SimulationOptions visualStopOptions=visualTimeLimitOptions;
    visualStopOptions.observationTime=1.0e-15;
    visualStopOptions.stopRequested=[](){return true;};
    const SimulationResult visualStopped=simulate(
        0x19a72ULL,3,visualStopOptions);
    const bool visualCensoringSemanticsOk=
        visualTimeLimited.outcome==SimulationOutcome::ObservationLimit
        &&visualTimeLimited.stopReason
            ==SimulationStopReason::ObservationTimeLimit
        &&observationDisposition(visualTimeLimited.stopReason)
            ==SimulationObservationDisposition::AdministrativelyCensored
        &&visualStopped.outcome==SimulationOutcome::ObservationLimit
        &&visualStopped.stopReason==SimulationStopReason::StopRequested
        &&observationDisposition(visualStopped.stopReason)
            ==SimulationObservationDisposition::AdministrativelyCensored
        &&observationDisposition(SimulationStopReason::ReachedCutoff)
            ==SimulationObservationDisposition::ObservedEndpoint
        &&observationDisposition(SimulationStopReason::NumericalFailure)
            ==SimulationObservationDisposition::NumericalFailure;
    State crossingBefore,crossingAfter;
    crossingBefore.firstPosition={1.1*nuclearCutoff,0.2*nuclearCutoff,0};
    crossingAfter.firstPosition={0.7*nuclearCutoff,-0.1*nuclearCutoff,0};
    const double exactCrossingFraction=separationCrossingFraction(
        crossingBefore,crossingAfter,nuclearCutoff);
    const State exactCrossingState=interpolateState(
        crossingBefore,crossingAfter,exactCrossingFraction);
    const double cutoffSurfaceResidual=std::abs(
        separation(exactCrossingState)-nuclearCutoff)/nuclearCutoff;
    MaxwellAmrHierarchy synchronizationTest({},24,3);
    synchronizationTest.finest().setVacuumPlaneWave(1.0e5);
    const std::size_t synchronizedCells=
        synchronizationTest.synchronizeFineToCoarse();
    MaxwellAmrHierarchy movingPatchTest({},24,3);
    movingPatchTest.finest().setLocalizedWavePacket(
        1.0e5,chargeCloudRestRadius);
    const double movingEnergyBefore=movingPatchTest.finest().fieldEnergy();
    const Vec3 requestedPatchCentre{0.25*movingPatchTest.finest().extent(),0,0};
    const std::size_t movedLevels=movingPatchTest.follow(requestedPatchCentre);
    const double movingEnergyRetention=movingPatchTest.finest().fieldEnergy()
                                      /movingEnergyBefore;
    MaxwellAmrHierarchy subcycleTest({},16,3);
    const double commonWaveNumber=2.0*pi/subcycleTest.levels().front().extent();
    for(MaxwellBlock& level:subcycleTest.levels())
        level.setVacuumPlaneWaveWavenumber(1.0e5,commonWaveNumber);
    const double subcycleEnergyBefore=
        subcycleTest.compositeVolumeIntegrals().energy;
    const double coarseDt=0.35*subcycleTest.levels().front().courantTimeStep();
    for(int step=0;step<4;++step) subcycleTest.advanceSubcycled(coarseDt);
    const double subcycleEnergyDrift=std::abs(
        subcycleTest.compositeVolumeIntegrals().energy-subcycleEnergyBefore)
        /subcycleEnergyBefore;
    const double coarseFineBoundaryMismatch=
        subcycleTest.levels()[1].maximumBoundaryMismatch(
            subcycleTest.levels().front())/1.0e5;
    const bool twoBodyRoleOk=protonElectronIncoming.valid()
        &&electronProtonIncoming.valid()
        &&protonElectronFree.valid()&&electronProtonFree.valid()
        &&std::isfinite(twoBodyRoleEnergyResidual)
        &&twoBodyRoleEnergyResidual<1.0e-10
        &&std::isfinite(twoBodyRoleMomentumResidual)
        &&twoBodyRoleMomentumResidual<1.0e-13
        &&std::isfinite(twoBodyRoleSwapResidual)
        &&twoBodyRoleSwapResidual<1.0e-13;
    const bool twoBodyBoostOk=asymmetricLab.valid()&&boostedIncoming.valid()
        &&boostedFree.valid()&&recoveredFirst.valid()&&recoveredSecond.valid()
        &&unequalBoostedIncoming.valid()&&unequalBoostedFree.valid()
        &&unequalRecoveredFirst.valid()&&unequalRecoveredSecond.valid()
        &&std::isfinite(twoBodyBoostEnergyResidual)
        &&twoBodyBoostEnergyResidual<1.0e-10
        &&std::isfinite(twoBodyBoostVelocityResidual)
        &&twoBodyBoostVelocityResidual<1.0e-12
        &&std::isfinite(twoBodyInverseBoostResidual)
        // The unequal-mass 0.56c round trip subtracts a proton boost momentum
        // about 1e5 times larger than the internal momentum represented by the
        // residual, so double precision limits this check near 1e-11.
        &&twoBodyInverseBoostResidual<1.0e-10;
    const bool twoBodyCausalOk=comovingNearLightPair.valid()
        &&comovingInternalEnergyEv==0.0
        &&allInitialVelocitiesCausal
        &&std::isfinite(maximumInitialBeta)&&maximumInitialBeta<1.0;
    const bool chargeOk = std::abs(depositedCharge/eCharge) < 1.0e-12;
    const bool gaussOk = std::isfinite(relativeGaussResidual)
                      && relativeGaussResidual < 2.0e-2;
    const bool divergenceOk = magneticDivergence < 1.0e-6;
    const bool energyOk = relativeEnergyDrift < 5.0e-3;
    const bool couplingOk=maximumRelativeContinuityResidual<1.0e-8
                       &&maximumRelativeLongitudinalCurl<1.0e-12
                       && coupledBeta<1.0 && dipoleNormResidual<1.0e-12
                       && isFinite(particles);
    // All three closures are held to the same 1e-3 band, because all three
    // measure the same thing: how well the coupled particle/grid step balances
    // a conserved quantity against the flux that left the box.
    //
    // The momentum bound used to be 1e-10, and it passed -- but it was never
    // measuring momentum conservation.  The probe launches the pair
    // antisymmetrically (+/-r, +/-v, opposite charges), and under that
    // symmetry every momentum contribution cancels against its mirror
    // identically, so the residual sat at 2e-16 no matter how the scheme
    // behaved.  Breaking the symmetry exposes the real number two ways, and
    // they agree: an unequal-mass pair gives 2.9e-5, and the DEFAULT e+e- pair
    // with one velocity merely halved gives 1.4e-4 -- the same order as the
    // energy (3.5e-4) and angular (4.5e-4) closures that were always checked
    // against 1e-3.  Halving the timestep leaves it unmoved (2.8941e-5 ->
    // 2.8935e-5), so this is the scheme's spatial closure floor, not a
    // convergence error a tighter step could reach.
    const bool balanceFinite=std::isfinite(coupledEnergyClosure)
        &&std::isfinite(coupledMomentumClosure)
        &&std::isfinite(coupledAngularClosure)
        &&std::abs(coupledEnergyClosure)<1.0e-3
        &&coupledMomentumClosure<1.0e-3
        // Was 1.0e-3.  Restoring the g in mu = gamma S did not move the
        // numerator by a single bit -- it stays 9.57308412e-38 -- but it
        // halved the denominator, because the |J_initial| scale contains the
        // intrinsic spin, and that was the very term inflated by g.  The
        // scale fell from 2.11399975e-34 to 1.05698206e-34, i.e. from
        // 2.0046 hbar to 1.0023 hbar, so the ratio doubled from 4.53e-04 to
        // 9.06e-04 with no physical change whatsoever.
        //
        // The bound is raised in the same proportion so the check keeps its
        // previous ABSOLUTE strictness.  Leaving it at 1.0e-3 would have made
        // it twice as strict as a side effect of a units correction, which
        // nobody intended, and would have left 9% of margin.
        &&coupledAngularClosure<2.0e-3;
    const bool boundaryOk=std::isfinite(absorbedFraction)
                       && absorbedFraction>0.05
                       &&std::isfinite(cpmlReflection)&&cpmlReflection<1.0e-6
                       && std::isfinite(integratedBoundaryEnergy)
                       && isFinite(integratedBoundaryMomentum)
                       && isFinite(integratedBoundaryAngularMomentum)
                       && synchronizedCells>0;
    const bool validationOk=relativeMagnetizationDivergence<1.0e-12
        &&std::isfinite(lateInteriorEnergyFraction)
        &&lateInteriorEnergyFraction<0.25
        &&fineGaussResidual<coarseGaussResidual;
    const bool movingAmrOk=movedLevels>=1
        &&std::abs(movingEnergyRetention-1.0)<1.0e-4
        &&std::isfinite(subcycleEnergyDrift)&&subcycleEnergyDrift<5.0e-3
        &&coarseFineBoundaryMismatch<1.0e-12;
    const bool covarianceOk=currentInvariantResidual<1.0e-12
        &&antisymmetryResidual<1.0e-12&&stressSymmetryResidual<1.0e-12
        &&std::abs(boostedCurrent.space.x/boostedCurrent.time
                   -covarianceBody.velocity.x/c)<1.0e-12;
    const bool particleCovarianceOk=covarianceTrajectoriesAdvanced
        &&std::isfinite(covarianceWorldlineResidual)
        &&std::isfinite(covarianceVelocityResidual)
        &&std::isfinite(covarianceForceResidual)
        &&std::isfinite(covarianceRadiationResidual)
        &&covarianceWorldlineResidual<5.0e-2
        &&covarianceVelocityResidual<5.0e-2
        // The remaining sub-percent defect measures the still noncovariant
        // gradient-force sector; tensor transformation itself is tested at
        // machine precision below.
        // 3.86e-3, and it is NOT an integration error: it is immune to the
        // step (8x refinement at fixed total time) and to the tolerance (1e-6
        // to 1e-12), constant in elapsed time, and exactly first order in the
        // orbital speed -- halving beta_orb halves it, zeroing it drops the
        // number 288-fold.  With everything else made exact it coincides with
        // the representability gap below to five digits, so the two lines
        // measure the same thing from opposite ends and the bound belongs to
        // the parametrization, not to the integrator.
        &&covarianceDipoleEvolutionResidual<5.0e-3
        &&std::isfinite(covarianceRepresentabilityGap)
        &&covarianceRepresentabilityGap<5.0e-3
        &&covarianceBmtResidual<1.0e-3
        // Exact agreement is expected: both sides run the same integrator on
        // the same state.  The band admits reordering round-off only.
        &&covarianceForceResidual<1.0e-1
        // A coordinate sphere at constant t is not mapped to the coordinate
        // sphere used in the boosted run.  This is therefore a finite-surface
        // diagnostic, not a strict four-vector pass criterion.  It must remain
        // bounded; the local worldline and four-force tests above are strict.
        // Bylo 2.5e-1, czyli prog dopasowany do bledu, nie do fizyki:
        // brakowalo czynnika Dopplera zamieniajacego moc odbierana na
        // emitowana.  Po jego przywroceniu reszta wynosi 7.0e-04, a reszta
        // skumulowanego czteropedu 2.1e-06.
        &&covarianceRadiationResidual<3.0e-3
        &&covarianceAccumulatedRadiationResidual<1.0e-4;
    const bool dipoleTensorCovarianceOk=dipoleTensorRoundtrip<1.0e-12
        &&dipoleFirstInvariantResidual<1.0e-12
        &&dipoleSecondInvariantResidual<1.0e-12
        &&inducedElectricDipoleResidual<1.0e-12
        &&staticElectricDipoleResidual<1.0e-8
        &&tensorGradientStaticResidual<1.0e-5
        &&std::isfinite(dipoleGradientCouplingInvarianceResidual)
        &&dipoleGradientCouplingInvarianceResidual<1.0e-2;
    const bool massAndSelfForceOk=bareMassFraction>0.0
        &&electromagneticMassFraction>0.0
        &&electromagneticMassFraction<0.01
        &&std::isfinite(maximumSelfForceFraction)
        &&std::isfinite(selfCalibrationPhaseJump)
        &&selfCalibrationPhaseJump<1.0e-5
        &&selfCalibrationChargeSymmetry<1.0e-12;
    const bool boundCurrentOk=std::isfinite(relativeBoundContinuity)
        &&relativeBoundContinuity<1.0e-3
        &&std::abs(boundCurrentTest.totalCharge()/eCharge)<1.0e-12;
    const bool retardedInitializationOk=std::isfinite(retardedInitialGauss)
        &&retardedInitialGauss<1.0e-3
        &&std::isfinite(retardedInitialDivB)&&retardedInitialDivB<1.0e-5
        &&std::isfinite(retardedAccelerationSignal)
        &&retardedAccelerationSignal>1.0e-12;
    const bool convergenceAndBoostOk=std::isfinite(maxwellConvergenceOrder)
        &&maxwellConvergenceOrder>1.7
        &&fineWaveError<coarseWaveError
        &&std::isfinite(lorentzFieldResidual)&&lorentzFieldResidual<1.0e-10
        &&std::isfinite(fieldInvariantResidual)&&fieldInvariantResidual<1.0e-10;
    const bool productionGeometryOk=productionInitialCoverage
        &&productionMovedCoverage&&productionMovedPatches>=1
        &&std::abs(firstPatchField.x-2.0)<1.0e-15
        &&std::abs(secondPatchField.x-3.0)<1.0e-15;
    const bool branchCouplingOk=branchCouplingStep.restrictedCells>=2
        &&serialBranchCouplingStep.restrictedCells
            ==branchCouplingStep.restrictedCells
        &&std::isfinite(branchRelativeEnergyDefect)
        &&branchRelativeEnergyDefect<1.0e-2
        &&std::isfinite(branchRelativeMomentumDefect)
        &&branchRelativeMomentumDefect<1.0e-2
        &&branchParallelResidual<1.0e-13;
    const bool yeeAndEsirkepovOk=std::isfinite(yeeContinuity)
        &&yeeContinuity<1.0e-12
        &&std::isfinite(yeeExtendedContinuity)
        &&yeeExtendedContinuity<1.0e-12
        &&std::abs(yeeGaussianCharge-1.0)<1.0e-12
        &&std::isfinite(yeeDivBAfter)
        &&yeeDivBAfter<=yeeDivBBefore+1.0e-12
        &&yeeAmrRefluxed==12U*12U*12U
        &&yeeAmrRestriction<1.0e-15
        &&yeeAmrContinuity<1.0e-12
        &&yeeAmrDivB<1.0e-12
        &&std::isfinite(yeeCpmlInteriorFraction)
        &&yeeCpmlInteriorFraction<0.2*yeePeriodicInteriorFraction
        &&yeeCpmlDivB<1.0e-8
        &&yeePusherReversibility<1.0e-12
        &&conservativeReverseResidual<1.0e-11
        &&yeeCoupledContinuity<1.0e-11
        &&std::abs(yeeCoupledCharge)<1.0e-12
        &&yeeCoupledBeta<1.0&&isFinite(yeeCoupledState);
    const bool sharedClassicalEngineOk=sharedEngineAdvanced
        &&sharedEngineResidual==0.0
        // These three are an IDENTITY, not a conservation test: boundField* is
        // accumulated every step as the residual that closes the balance, so
        // the sum telescopes and has to come out zero.  Measured: breaking the
        // probe's antisymmetry completely, by taking the second velocity to
        // zero, does not move them at all -- they stay at 0 / 2.4e-19 /
        // 5.7e-26.  The 1e-12 bound is therefore a bound on roundoff and
        // certifies nothing beyond it.
        &&sharedEnergyBalanceResidual<1.0e-12
        &&sharedMomentumBalanceResidual<1.0e-12
        &&sharedAngularBalanceResidual<1.0e-12
        // These three are the MEASUREMENT: they exclude boundField*, so they do
        // not close by construction.  Until now they were printed and never
        // checked, so they could have drifted arbitrarily while the suite kept
        // reporting 28/28.
        //
        // 3e-3 for all three, the momentum one included.  The shipped probe is
        // antisymmetric (+/-0.025c), which puts its momentum residual at
        // 1.3e-08 -- five orders better than energy and angular momentum, but
        // that is the symmetry talking, not the accuracy.  Halving the second
        // velocity gives 2.6e-04 and zeroing it gives 8.1e-04, the same order
        // as energy (8.4e-04) and angular momentum (1.8e-04).  A bound fitted
        // to 1.3e-08 would be a bound on the probe's setup, not on the code.
        &&std::isfinite(sharedRawEnergyResidual)
        &&sharedRawEnergyResidual<3.0e-3
        &&std::isfinite(sharedRawMomentumResidual)
        &&sharedRawMomentumResidual<3.0e-3
        &&std::isfinite(sharedRawAngularResidual)
        &&sharedRawAngularResidual<3.0e-3
        &&std::isfinite(photonFourEnergyResidual)
        &&photonFourEnergyResidual<1.0e-8
        &&std::isfinite(photonFourMomentumResidual)
        &&photonFourMomentumResidual<1.0e-8
        &&!forbiddenPhoton.emitted&&forbiddenPhotonMutation==0.0
        &&std::isfinite(reactionMismatchFraction)
        // Was isfinite alone.  llValidity is the Landau-Lifshitz validity
        // parameter and has to be far below one for the reduced-order form to
        // mean anything; it runs 8.0e-05 for e+e- down to 4.8e-08 for p+pbar,
        // so 1e-2 leaves two orders of margin on the worst pair.
        &&std::isfinite(llValidity)&&llValidity<1.0e-2
        // Compactness kR is the validity condition of the electric-dipole
        // approximation the whole radiative sector rests on: 0.049 for e+e-
        // and smaller for every heavier pair.
        &&std::isfinite(blendCompactness)&&blendCompactness<0.1
        &&std::isfinite(blendSmoothness)
        &&std::isfinite(blendWeight)&&blendWeight>=0.0&&blendWeight<=1.0
        // Consistency of the gate with its own threshold, rather than a bound
        // on the smoothness itself, which legitimately blows up: for mu+mu-
        // the two step sizes disagree by 2.611, i.e. the third derivative is
        // numerical noise there.  What must then hold is that the gate SAW it
        // and zeroed the weight.  Checked across the four pairs: 3.97e-04 -> 1,
        // 1.66e-03 -> 1, 2.611 -> 0, 2.79e-02 -> 0.552.
        &&(blendSmoothness<=5.0e-2||blendWeight==0.0)
        &&quadrupoleTrace<1.0e-14
        &&quadrupoleSymmetry<1.0e-14
        &&quadrupoleResidual<1.0e-14
        &&coherentDerivativeResidual<1.0e-12
        &&coherentReactionMomentumResidual<1.0e-14;
    // Until now this checked isfinite and nothing else, so nine numbers in the
    // radiation-reaction sector -- the model's core claim -- could take any
    // value at all and the suite would still report 28/28.  They are bounded
    // now.
    //
    // reactionFluxResidual is deliberately left ungated: its denominator is
    // the energy radiated over 8e-22 s, which is meaningless for heavy pairs,
    // and it ranges from 0.026 to 1893 across the four pairs for that reason
    // alone.  reactionMechanicalMismatch measures the same mismatch against a
    // scale that survives, and is bounded instead.
    const auto finiteReactionBenchmark=[](const ReactionModelBenchmark& value) {
        return value.advanced&&isFinite(value.finalState)
            &&std::isfinite(value.rawEnergyResidual)
            &&value.rawEnergyResidual<1.0e-4
            &&std::isfinite(value.reactionFluxResidual)
            &&std::isfinite(value.reactionMechanicalMismatch)
            &&value.reactionMechanicalMismatch<1.0e-5
            &&std::isfinite(value.stepConvergence)
            &&value.stepConvergence<1.0e-4
            &&std::isfinite(value.wallSeconds)&&value.wallSeconds>0.0;
    };
    const bool reactionModelsOk=finiteReactionBenchmark(reactionDisabled)
        &&finiteReactionBenchmark(reactionLl)
        &&finiteReactionBenchmark(reactionCoherent)
        &&finiteReactionBenchmark(reactionAutomatic)
        &&std::isfinite(coherentCostRatio)&&coherentCostRatio<10.0;
    const bool coherentMagneticDipoleOk=
        std::isfinite(magneticDerivativeResidual)
        &&magneticDerivativeResidual<1.0e-12
        &&std::isfinite(alignedMagneticPowerRatio)
        &&std::abs(alignedMagneticPowerRatio-4.0)<1.0e-12
        &&std::isfinite(cancellingMagneticPowerRatio)
        &&std::abs(cancellingMagneticPowerRatio)<1.0e-12
        &&std::isfinite(magneticTorqueResidual)
        &&magneticTorqueResidual<1.0e-12
        &&std::isfinite(magneticAngularFluxResidual)
        &&magneticAngularFluxResidual<1.0e-12
        &&std::isfinite(magneticMomentumResidual)
        &&magneticMomentumResidual<1.0e-12
        &&std::isfinite(magneticComMomentumResidual)
        &&magneticComMomentumResidual<1.0e-12;
    const bool farFieldConvergenceOk=std::isfinite(farReference.energy)
        &&farReference.energy>=0.0
        &&std::ranges::all_of(farDirectionResiduals,
            [](double value){return std::isfinite(value);})
        &&std::ranges::all_of(farRadiusResiduals,
            [](double value){return std::isfinite(value);})
        &&std::ranges::all_of(farNearFieldContamination,
            [](double value){return std::isfinite(value);})
        // The symmetric 50-node Lebedev rule can outperform the denser
        // non-tabulated Fibonacci diagnostic, so strict monotonicity in N is
        // neither expected nor desirable here.
        &&farDirectionResiduals[1]<farDirectionResiduals[0]
        &&farDirectionResiduals[2]<farDirectionResiduals[0]
        &&farRadiusResiduals[2]<farRadiusResiduals[1]
        &&farRadiusResiduals[1]<farRadiusResiduals[0]
        &&farNearFieldContamination[2]<farNearFieldContamination[1]
        &&farNearFieldContamination[1]<farNearFieldContamination[0]
        &&farDirectionResiduals.back()<5.0e-2
        &&farRadiusResiduals.back()<1.0e-5
        &&farNearFieldContamination.back()<1.0e-2
        &&lebedevRule.size()==50U&&lebedevWeightResidual<1.0e-14
        &&lebedevMomentResidual<1.0e-14;
    const bool trajectoryConvergenceOk=isFinite(trajectoryReference)
        &&std::ranges::all_of(trajectoryToleranceResiduals,
            [](double value){return std::isfinite(value);})
        &&trajectoryToleranceResiduals[2]<trajectoryToleranceResiduals[1]
        &&trajectoryToleranceResiduals[1]<trajectoryToleranceResiduals[0]
        &&trajectoryToleranceResiduals.back()<1.0e-5
        &&std::isfinite(trajectoryStepResidual)
        &&trajectoryStepResidual<1.0e-5;
    const bool adaptiveDepthRejectionOk=!depthLimitAccepted
        &&depthRollbackResidual==0.0&&depthHistoryUnchanged&&depthResolved
        &&std::isfinite(depthResolvedTimeResidual)
        &&depthResolvedTimeResidual<1.0e-13;
    const bool causalStartupOk=startupFinite
        &&startupCoverage>=7.99
        &&std::abs(startupHistory.back().time-yeeCoupledState.time)<1.0e-30
        &&startupHistoryResidual4<1.0e-5
        &&startupHistoryResidual8<1.0e-6;
    const bool historyConstructionSensitivityOk=
        startupCase4.advanced&&startupCase8.advanced&&startupCase16.advanced
        &&std::isfinite(startupEarlySensitivity)
        &&std::isfinite(startupSettledSensitivity)
        // The hidden nodes disappear, but the small phase perturbation they
        // launched need not contract in a Hamiltonian orbit.  Require that it
        // remain bounded rather than falsely demanding dissipation: less than
        // twofold amplification over the replacement interval, plus the
        // independent absolute accuracy ceiling below.
        &&startupSettledSensitivity
            <2.0*std::max(startupEarlySensitivity,1.0e-15)
        &&startupSettledSensitivity<1.0e-5;
    const bool parameterSensitivityScanOk=
        std::ranges::all_of(cutoffArrivalTimes,
            [](double value){return std::isfinite(value)&&value>0.0;})
        &&cutoffArrivalTimes[2]<cutoffArrivalTimes[1]
        &&cutoffArrivalTimes[1]<cutoffArrivalTimes[0]
        &&std::ranges::all_of(cutoffEventResiduals,
            [](double value){return std::isfinite(value)&&value<1.0e-10;})
        &&std::ranges::all_of(regularizationScanStates,
            [](const State& value){return isFinite(value);})
        &&std::ranges::all_of(regularizationScanAdvanced,
            [](bool value){return value;})
        &&std::ranges::all_of(regularizationTrajectoryResiduals,
            [](double value){return std::isfinite(value)&&value<1.0;})
        // Guards against a scan that accidentally stops routing the mutable
        // radius into the force law and returns three bit-identical outcomes.
        &&regularizationTrajectoryResiduals[0]>1.0e-12
        &&regularizationTrajectoryResiduals[2]>1.0e-12;
    const bool retardedInterpolationOk=
        std::isfinite(hermiteConvergenceOrder)
        &&interpolationFine[1]<interpolationFine[0]
        &&interpolationFine[1]<interpolationCoarse[1]
        &&hermiteConvergenceOrder>1.8;
    const bool shortRangeRegularizationOk=
        retardedStaticLimitResidual<1.0e-12
        &&std::ranges::all_of(regulatorFieldResiduals,
            [](double value){return std::isfinite(value);})
        &&regulatorFieldResiduals[2]<regulatorFieldResiduals[1]
        &&regulatorFieldResiduals[1]<regulatorFieldResiduals[0]
        &&regulatorFieldResiduals[1]<0.2
        &&std::ranges::all_of(regulatorProfileResiduals,
            [](double value){return std::isfinite(value)&&value<1.0e-3;})
        // A steeper exponent must converge faster at the same radius.
        &&regulatorProfileResiduals[1]<regulatorProfileResiduals[0]
        &&regulatorProfileResiduals[2]<regulatorProfileResiduals[1]
        &&std::isfinite(regulatorCoreSuppression)
        &&regulatorCoreSuppression<1.0e-9
        &&isFinite(regulatorOriginField)&&regulatorOriginField.norm()==0.0
        &&isFinite(regulatorOriginForce)&&regulatorOriginForce.norm()==0.0
        // Check both the requested cap and the analytic peak of curl(A).  The
        // lower cap bound prevents an accidentally oversized radius from
        // passing, while the measured 2.778 factor prevents the validator from
        // silently returning to the scalar w/r^3 surrogate.
        &&peakDipoleEnergyOverRestEnergy<0.5005
        &&peakDipoleEnergyOverRestEnergy>0.4995
        &&peakTransverseDipoleEnergyOverRestEnergy<0.5005
        &&peakTransverseDipoleEnergyOverRestEnergy>0.4995
        &&peakRadialDipoleEnergyOverRestEnergy<0.2052
        &&peakRadialDipoleEnergyOverRestEnergy>0.2050
        &&std::abs(measuredDipoleCurlPeak-regulatorDipoleCurlPeak)<1.0e-4
        &&exactParallelDipoleEnergy<0.0
        &&exactAntiparallelDipoleEnergy>0.0
        &&std::abs(exactDipolePeakOverRestEnergy
                    -dipoleEnergyCeilingFraction)<1.0e-12
        &&regulatorFarAlignment>0.0
        &&isFinite(cutoffDipoleForce)
        &&cutoffSurfaceResidual<1.0e-14;

    const double benchmarkSeconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-benchmarkStart).count();
    std::size_t hierarchyBytes=0;
    for(const MaxwellBlock& level:hierarchy.levels())
        hierarchyBytes+=level.cells().size()*sizeof(MaxwellCell);
    const double cflStep=hierarchy.timeStep();
    const double estimatedStepsPerPicosecond=1.0e-12/cflStep;
    const double measuredFieldSteps=40.0+24.0*1.5+160.0;
    const double fieldStepsPerSecond=measuredFieldSteps
        /std::max(benchmarkSeconds,1.0e-12);
    const double estimatedSecondsPerPicosecond=estimatedStepsPerPicosecond
        /fieldStepsPerSecond;
    std::cout << std::setprecision(8)
              << "Maxwell-AMR self-test\n"
              << "2-body role E/P:    " << twoBodyRoleEnergyResidual << " / "
              << twoBodyRoleMomentumResidual << '\n'
              << "2-body role swap:   " << twoBodyRoleSwapResidual << '\n'
              << "2-body boost K/V:   " << twoBodyBoostEnergyResidual << " / "
              << twoBodyBoostVelocityResidual << '\n'
              << "2-body inverse:     " << twoBodyInverseBoostResidual << '\n'
              << "2-body comoving K:  " << comovingInternalEnergyEv << " eV\n"
              << "2-body max beta:    " << std::setprecision(16)
              << maximumInitialBeta << std::setprecision(8) << '\n'
              << "levels:             " << hierarchy.levels().size() << '\n'
              << "finest dx/r_c:      " << finest.cellSize()/chargeCloudRestRadius << '\n'
              << "neutral charge/e:   " << depositedCharge/eCharge << '\n'
              << "Gauss residual:     " << gaussBefore << " V/m2\n"
              << "relative Gauss:     " << relativeGaussResidual << '\n'
              << "max |div B|:        " << magneticDivergence << " T/m\n"
              << "vacuum energy drift:" << relativeEnergyDrift*100.0 << "%\n"
              << "continuity residual:" << maximumRelativeContinuityResidual << '\n'
              << "deposit curl:       " << maximumRelativeLongitudinalCurl << '\n'
              << "coupled max beta:   " << coupledBeta << '\n'
              << "dipole norm residual:" << dipoleNormResidual << '\n'
              << "particle-field dE:  " << coupledEnergyClosure << '\n'
              << "particle-field dP:  " << coupledMomentumClosure << '\n'
              << "particle-field dJ:  " << coupledAngularClosure << '\n'
              << "absorbed fraction:  " << absorbedFraction << '\n'
              << "boundary flux/E0:   "
              << integratedBoundaryEnergy/absorbingInitialEnergy << '\n'
              << "late interior E/E0: " << lateInteriorEnergyFraction << '\n'
              << "CPML reflection:    " << cpmlReflection << '\n'
              << "div J_M relative:   " << relativeMagnetizationDivergence << '\n'
              << "Gauss coarse/fine:  " << coarseGaussResidual << " / "
              << fineGaussResidual << '\n'
              << "AMR synchronized:   " << synchronizedCells << " cells\n"
              << "moving levels:      " << movedLevels << '\n'
              << "regrid E retention: " << movingEnergyRetention << '\n'
              << "subcycle dE:        " << subcycleEnergyDrift << '\n'
              << "coarse/fine shell:  " << coarseFineBoundaryMismatch << '\n'
              << "J four-norm drift:  " << currentInvariantResidual << '\n'
              << "M antisymmetry:     " << antisymmetryResidual << '\n'
              << "T symmetry:         " << stressSymmetryResidual << '\n'
              << "m_em/m_physical:    " << electromagneticMassFraction << '\n'
              << "m_bare/m_physical:  " << bareMassFraction << '\n'
              << "grid self-force/F0: " << maximumSelfForceFraction << '\n'
              << "self phase seam:    " << selfCalibrationPhaseJump << '\n'
              << "self charge sign:   " << selfCalibrationChargeSymmetry << '\n'
              << "bound continuity:   " << relativeBoundContinuity << '\n'
              << "bound net charge/e: " << boundCurrentTest.totalCharge()/eCharge << '\n'
              << "retarded Gauss:     " << retardedInitialGauss << '\n'
              << "retarded div B:     " << retardedInitialDivB << '\n'
              << "retarded accel:     " << retardedAccelerationSignal << '\n'
              << "wave coarse/fine:  " << coarseWaveError << " / "
              << fineWaveError << '\n'
              << "Maxwell order:      " << maxwellConvergenceOrder << '\n'
              << "Lorentz field:      " << lorentzFieldResidual << '\n'
              << "field invariants:   " << fieldInvariantResidual << '\n'
              << "particle boost run: " << covarianceTrajectoriesAdvanced << '\n'
              << "particle boost x/v:" << covarianceWorldlineResidual << " / "
              << covarianceVelocityResidual << '\n'
              << "particle boost D:  " << covarianceDipoleEvolutionResidual << '\n'              << "dipole repr gap:   " << covarianceRepresentabilityGap << '\n'
              << "covariant BMT:    " << covarianceBmtResidual << '\n'
              << "BMT vs eff field: " << bmtEffectiveFieldGapActiveG << " / "
                 << bmtEffectiveFieldGapHighBeta
                 << "  (closed form vs integrated specification)\n"
              << "BMT precession:   " << bmtOrthogonalityActiveG << " / "
                 << bmtOrthogonalityHighBeta
                 << "  (mu.dmu/dt, must be 0 for a precession)\n"
              << "quantized drain:  " << quantizedChargeReactionDrain << " / "
                 << quantizedDipoleConstraintDrain
                 << "  (charge force / dipole reservoir; both must be 0)\n"
              << "quantized torque: " << quantizedDipoleTorqueTravel << " / "
                 << disabledDipoleTorqueTravel << " / "
                 << quantizedDipoleTorqueNormDrift
                 << "  (stochastic travel / disabled travel / norm drift)\n"
              << "BMT norm (2 rad): " << bmtNormDriftActiveG << " / "
                 << bmtNormDriftHighBeta
                 << "  (production routine, no renormalization)\n"
              << "secular spin J/norm: "
                 << secularReference.relativeAngularMomentumResidual << " / "
                 << std::max(secularFirstNormDrift,secularSecondNormDrift)
                 << "\n"
              << "secular spin coarse/fine/ref: "
                 << secularCoarseFineDifference << " / "
                 << secularFineReferenceDifference << "  ("
                 << secularCoarse.substeps << "/" << secularFine.substeps
                 << "/" << secularReference.substeps << " substeps)\n"
              << "secular external dJ/hbar: "
                 << secularExternal.externalAngularMomentumTransfer.norm()/hbar
                 << '\n'
              << "secular eccentric BMT/orbit: "
                 << secularEccentricResolvedResidual << " / "
                 << secularEccentricOrbitClosure
                 << "  (e^2=" << secularReferenceEccentricitySquared << ")\n"
              << "role routing:     " << roleRoutingResidual
              << "  (obrot " << roleRoutingTravel << ")\n"
              << "particle boost F/R:" << covarianceForceResidual << " / "
              << covarianceRadiationResidual << '\n'
              << "boost accumulated R:" << covarianceAccumulatedRadiationResidual
              << '\n'
              << "dipole boost back:  " << dipoleTensorRoundtrip << '\n'
              << "dipole invariants:  " << dipoleFirstInvariantResidual << " / "
              << dipoleSecondInvariantResidual << '\n'
              << "induced electric p:" << inducedElectricDipoleResidual << '\n'
              << "electric dipole E: " << staticElectricDipoleResidual << '\n'
              << "tensor grad static:" << tensorGradientStaticResidual << '\n'
              << "dipole coupling inv:" << dipoleGradientCouplingInvarianceResidual << '\n'
              << "pair patch cover:  " << productionInitialCoverage << " / "
              << productionMovedCoverage << '\n'
              << "pair patches moved:" << productionMovedPatches << '\n'
              << "branch restricted:  " << branchCouplingStep.restrictedCells << '\n'
              << "branch dE:          " << branchRelativeEnergyDefect << '\n'
              << "branch dP:          " << branchRelativeMomentumDefect << '\n'
              << "branch parallel:    " << branchParallelResidual << '\n'
              << "branch speedup:     " << branchParallelSpeedup << '\n'
              << "Yee continuity:     " << yeeContinuity << '\n'
              << "Yee Gauss+P/M cont.:" << yeeExtendedContinuity << '\n'
              << "Yee Gaussian q/e:   " << yeeGaussianCharge << '\n'
              << "Yee div B:          " << yeeDivBBefore << " / "
              << yeeDivBAfter << '\n'
              << "Yee AMR restriction:" << yeeAmrRestriction << '\n'
              << "Yee AMR continuity: " << yeeAmrContinuity << '\n'
              << "Yee AMR div B:      " << yeeAmrDivB << '\n'
              << "Yee CPML interior:  " << yeeCpmlInteriorFraction << " / "
              << yeePeriodicInteriorFraction << '\n'
              << "Yee CPML div B:     " << yeeCpmlDivB << '\n'
              << "Yee pusher reverse: " << yeePusherReversibility << '\n'
              << "conservative reverse:" << conservativeReverseResidual << '\n'
              << "Yee coupled cont.:  " << yeeCoupledContinuity << '\n'
              << "Yee coupled q/e:    " << yeeCoupledCharge << '\n'
              << "Yee coupled beta:   " << yeeCoupledBeta << '\n'
              << "Visual/stat engine: " << sharedEngineResidual << '\n'
              << "pair-field IDENTITY:" << sharedEnergyBalanceResidual << " / "
              << sharedMomentumBalanceResidual << " / "
              << sharedAngularBalanceResidual << '\n'
              << "raw dE/dP/dJ:       " << sharedRawEnergyResidual << " / "
              << sharedRawMomentumResidual << " / " << sharedRawAngularResidual
              << "  (IDENTITY closes by construction; raw is the measurement)"
              << '\n'
              << "photon recoil dE/dP: " << photonFourEnergyResidual << " / "
              << photonFourMomentumResidual << '\n'
              << "reaction/flux dE:   " << reactionMismatchFraction << '\n'
              << "LL validity |F|:    " << llValidity << '\n'
              << "reaction raw off/LL/C:" << reactionDisabled.rawEnergyResidual
              << " / " << reactionLl.rawEnergyResidual << " / "
              << reactionCoherent.rawEnergyResidual << '\n'
              << "  of which dE_mech/|E| = "
              << reactionDisabled.mechanicalEnergyChange
              << " and E_rad/|E| = " << reactionDisabled.radiatedEnergyFraction
              << " (reaction off, so no force removes that flux from the "
                 "mechanics)\n"
              << "  refinement ratio (full/half; ~1 means the residual is "
                 "PHYSICAL, not discretization): "
              << (reactionDisabled.refinedEnergyResidual>0.0
                  ?reactionDisabled.rawEnergyResidual
                      /reactionDisabled.refinedEnergyResidual:0.0)
              << " / "
              << (reactionLl.refinedEnergyResidual>0.0
                  ?reactionLl.rawEnergyResidual
                      /reactionLl.refinedEnergyResidual:0.0)
              << " / "
              << (reactionCoherent.refinedEnergyResidual>0.0
                  ?reactionCoherent.rawEnergyResidual
                      /reactionCoherent.refinedEnergyResidual:0.0)
              << '\n'
              << "reaction flux off/LL/C:"
              << reactionDisabled.reactionFluxResidual << " / "
              << reactionLl.reactionFluxResidual << " / "
              << reactionCoherent.reactionFluxResidual << '\n'
              << "reaction mism/E_mech:" << reactionDisabled.reactionMechanicalMismatch
              << " / " << reactionLl.reactionMechanicalMismatch << " / "
              << reactionCoherent.reactionMechanicalMismatch << '\n'
              << "reaction step off/LL/C:"
              << reactionDisabled.stepConvergence << " / "
              << reactionLl.stepConvergence << " / "
              << reactionCoherent.stepConvergence << '\n'
              << "coherent/LL cost:   " << coherentCostRatio << '\n'
              << "automatic raw/flux/step:"
              << reactionAutomatic.rawEnergyResidual << " / "
              << reactionAutomatic.reactionFluxResidual << " / "
              << reactionAutomatic.stepConvergence << '\n'
              << "auto smooth/kR/w:   "
              << sharedFinalRadiation.coherentDerivativeConsistency << " / "
              << sharedFinalRadiation.sourceCompactness << " / "
              << sharedFinalRadiation.coherentWeight << '\n'
              << "E2 vs analytic:      " << quadrupoleResidual << '\n'
              << "E2 |Q|/scale:        " << quadrupoleMagnitude << '\n'
              << "E2 trace/symmetry:   " << quadrupoleTrace << " / "
              << quadrupoleSymmetry << '\n'
              << "coherent E1 d3/Fsum:" << coherentDerivativeResidual << " / "
              << coherentReactionMomentumResidual << '\n'
              << "coherent M1 +/-:     " << alignedMagneticPowerRatio << " / "
              << cancellingMagneticPowerRatio << '\n'
              << "coherent M1 T/J/P/PC:" << magneticTorqueResidual << " / "
              << magneticAngularFluxResidual << " / "
              << magneticMomentumResidual << " / "
              << magneticComMomentumResidual << '\n'
              << "far N=26/50/98:     " << farDirectionResiduals[0] << " / "
              << farDirectionResiduals[1] << " / "
              << farDirectionResiduals[2] << '\n'
              << "Lebedev w/moment:   " << lebedevWeightResidual << " / "
              << lebedevMomentResidual << '\n'
              << "far R=1e4/5/6 a0:   " << farRadiusResiduals[0] << " / "
              << farRadiusResiduals[1] << " / "
              << farRadiusResiduals[2] << '\n'
              << "Larmor norm ratio:  " << larmorNormalizationRatio
              << " (radiation-only) / " << larmorProductionRatio
              << " (production sampling)\n"
              << "Larmor accumulation:" << larmorAccumulationRatio << '\n'
              << "long balance dE/Efar:" << longHorizonBalanceResidual
              << "  (Efar=" << longHorizonRadiatedEnergy/eCharge
              << " eV, dE_S=" << longHorizonSchottChange/eCharge << " eV)\n"
              << "long balance refined:" << longHorizonBalanceRefined
              << "  (tol 1e-8; discretization shrinks, a real leak does not)\n"
              << "balance amplification:" << longHorizonBalanceAmplification
              << "  (|E_mech|/E_far: gain on any dE_mech error)\n"
              << "balance radius sweep (R/a_pair: amplif / coarse / refined):\n"
              << "  " << balanceRadiusSweep[0].factor << ": "
              << balanceRadiusSweep[0].amplification << " / "
              << balanceRadiusSweep[0].coarse << " / "
              << balanceRadiusSweep[0].refined << '\n'
              << "  " << balanceRadiusSweep[1].factor << ": "
              << balanceRadiusSweep[1].amplification << " / "
              << balanceRadiusSweep[1].coarse << " / "
              << balanceRadiusSweep[1].refined << '\n'
              << "  " << balanceRadiusSweep[2].factor << ": "
              << balanceRadiusSweep[2].amplification << " / "
              << balanceRadiusSweep[2].coarse << " / "
              << balanceRadiusSweep[2].refined
              << "   (amplification must follow R^1.5)\n"
              << "balance matrix M/F,S/F,R:\n"
              << "  retarded + none: "
              << balanceDiagnostics[0].mechanicalOverFlux << " / "
              << balanceDiagnostics[0].schottOverFlux << " / "
              << balanceDiagnostics[0].signedResidual << '\n'
              << "  Darwin + coherent:"
              << balanceDiagnostics[1].mechanicalOverFlux << " / "
              << balanceDiagnostics[1].schottOverFlux << " / "
              << balanceDiagnostics[1].signedResidual << '\n'
              << "  retarded + self: "
              << balanceDiagnostics[2].mechanicalOverFlux << " / "
              << balanceDiagnostics[2].schottOverFlux << " / "
              << balanceDiagnostics[2].signedResidual << '\n'
              << "  retarded + coher:"
              << balanceDiagnostics[3].mechanicalOverFlux << " / "
              << balanceDiagnostics[3].schottOverFlux << " / "
              << balanceDiagnostics[3].signedResidual << '\n'
              << "near R=1e4/5/6 a0:  " << farNearFieldContamination[0] << " / "
              << farNearFieldContamination[1] << " / "
              << farNearFieldContamination[2] << '\n'
              << "step tol=1e-5/6/7:  " << trajectoryToleranceResiduals[0]
              << " / " << trajectoryToleranceResiduals[1] << " / "
              << trajectoryToleranceResiduals[2] << '\n'
              << "step dt/(dt/2):     " << trajectoryStepResidual << '\n'
              << "depth reject/rollback/resolve: " << depthLimitAccepted
              << " / " << depthRollbackResidual << " / " << depthResolved
              << '\n'
              << "history span r/c:   " << startupCoverage << '\n'
              << "history 4/8 vs 16: " << startupHistoryResidual4 << " / "
              << startupHistoryResidual8 << '\n'
              << "history construct early/settled: "
              << startupEarlySensitivity << " / "
              << startupSettledSensitivity << '\n'
              << "cutoff scan f=0.004/5/6 t/t_pair: "
              << cutoffArrivalTimes[0] << " / " << cutoffArrivalTimes[1]
              << " / " << cutoffArrivalTimes[2] << '\n'
              << "cutoff event residuals: " << cutoffEventResiduals[0]
              << " / " << cutoffEventResiduals[1] << " / "
              << cutoffEventResiduals[2] << '\n'
              << "reg-radius scan .5/1/2 residual: "
              << regularizationTrajectoryResiduals[0] << " / "
              << regularizationTrajectoryResiduals[1] << " / "
              << regularizationTrajectoryResiduals[2] << '\n'
              << "reg-radius scan advanced .5/1/2: "
              << regularizationScanAdvanced[0] << " / "
              << regularizationScanAdvanced[1] << " / "
              << regularizationScanAdvanced[2] << '\n'
              << "history linear/H:  " << interpolationFine[0] << " / "
              << interpolationFine[1] << '\n'
              << "Hermite order:     " << hermiteConvergenceOrder << '\n'
              << "dipole src linear: " << sourceInterpolationFine
                 << "  (order " << sourceInterpolationOrder << ")\n"
              << "dipole static limit:" << retardedStaticLimitResidual << '\n'
              << "reg a=.75/.5/.25 rc:" << regulatorFieldResiduals[0] << " / "
              << regulatorFieldResiduals[1] << " / "
              << regulatorFieldResiduals[2] << '\n'
              << "reg p=4/6/8 @10a:  " << regulatorProfileResiduals[0] << " / "
              << regulatorProfileResiduals[1] << " / "
              << regulatorProfileResiduals[2] << '\n'
              << "reg core suppress: " << regulatorCoreSuppression << '\n'
              << "reg peak U T/R:    "
              << peakTransverseDipoleEnergyOverRestEnergy << " / "
              << peakRadialDipoleEnergyOverRestEnergy << '\n'
              << "reg curl(A) peak:  " << measuredDipoleCurlPeak << '\n'
              << "reg exact U/m c2:  " << exactDipolePeakOverRestEnergy << '\n'
              << "cutoff surface:     " << cutoffSurfaceResidual << '\n'
              << "CFL dt:             " << cflStep << " s\n"
              << "field storage:      "
              << static_cast<double>(hierarchyBytes)/(1024.0*1024.0) << " MiB\n"
              << "statistics profile: "
              << (publicationStatistics?"publication":"small")
              << " (N=" << statisticalSampleCount << ")\n"
              << "photon <cos>/<cos2>:" << photonCosineMean << " / "
              << photonCosineSecondMoment << " (expected 0 / 0.4)\n"
              << "photon azimuth mean:" << photonAzimuthVectorMean << '\n'
              << "Poisson Exp mean/var:" << exponentialMean << " / "
              << exponentialVariance << " (expected 1 / 1)\n"
              << "statistics replay:  "
              << (statisticalReplayExact?"exact":"MISMATCH") << '\n'
              << "impact <b/s>/<b2/s2>:" << normalizedImpactMean << " / "
              << normalizedImpactSecondMoment
              << " (Rayleigh: " << std::sqrt(pi/2.0) << " / 2)\n"
              << "impact azimuth mean:" << impactAzimuthVectorMean << '\n'
              << "impact replay:      "
              << (impactReplayExact?"exact":"MISMATCH") << '\n'
              << "validation wall:    " << benchmarkSeconds << " s\n"
              << "field-step rate:    " << fieldStepsPerSecond << " steps/s\n"
              << "steps per ps:       " << estimatedStepsPerPicosecond << '\n'
              << "estimated s/ps:     " << estimatedSecondsPerPicosecond << '\n';

    // One source of truth for both the categorized report and exit status.
    // A section is semantic, not cosmetic: identities establish algebraic
    // wiring, regressions protect numerical behaviour, convergence checks
    // compare resolutions, and physical-domain guards only assert that a
    // documented approximation remains inside its declared validity band.
    // In particular a PASS in the last section is not experimental validation
    // of CREM or QED.
    enum class ValidationSection {
        AlgebraicIdentity,
        NumericalRegression,
        Convergence,
        IndependentBalance,
        PhysicalDomain
    };
    struct ValidationCheck {
        ValidationSection section;
        const char* name;
        bool passed;
    };
    // The two energy gates are exact.  The controlled torque must rotate the
    // stochastic state by a clearly resolved amount, leave the disabled
    // reference unchanged, and preserve both proper-dipole norms.
    const bool quantizedRadiationOk=
        std::isfinite(quantizedChargeReactionDrain)
        && std::isfinite(quantizedDipoleConstraintDrain)
        && std::isfinite(quantizedDipoleTorqueTravel)
        && std::isfinite(disabledDipoleTorqueTravel)
        && std::isfinite(quantizedDipoleTorqueNormDrift)
        && quantizedChargeReactionDrain==0.0
        && quantizedDipoleConstraintDrain==0.0
        && quantizedDipoleTorqueTravel>0.1
        && disabledDipoleTorqueTravel==0.0
        && quantizedDipoleTorqueNormDrift<1.0e-12;
    const std::array<ValidationCheck,46> regressionChecks{{
        {ValidationSection::AlgebraicIdentity,"two-body-role-invariance",twoBodyRoleOk},
        {ValidationSection::NumericalRegression,"two-body-lorentz-boost",twoBodyBoostOk},
        {ValidationSection::PhysicalDomain,"two-body-causality",twoBodyCausalOk},
        {ValidationSection::AlgebraicIdentity,"charge",chargeOk},
        {ValidationSection::NumericalRegression,"gauss",gaussOk},
        {ValidationSection::AlgebraicIdentity,"divergence",divergenceOk},
        {ValidationSection::NumericalRegression,"energy",energyOk},
        {ValidationSection::NumericalRegression,"coupling",couplingOk},
        {ValidationSection::NumericalRegression,"balance-finite",balanceFinite},
        {ValidationSection::NumericalRegression,"boundary",boundaryOk},
        {ValidationSection::NumericalRegression,"validation",validationOk},
        {ValidationSection::NumericalRegression,"moving-amr",movingAmrOk},
        {ValidationSection::AlgebraicIdentity,"covariance",covarianceOk},
        {ValidationSection::PhysicalDomain,"particle-covariance",particleCovarianceOk},
        {ValidationSection::AlgebraicIdentity,"dipole-tensor-covariance",dipoleTensorCovarianceOk},
        {ValidationSection::NumericalRegression,"mass-and-self-force",massAndSelfForceOk},
        {ValidationSection::AlgebraicIdentity,"bound-current",boundCurrentOk},
        {ValidationSection::NumericalRegression,"retarded-initialization",retardedInitializationOk},
        {ValidationSection::Convergence,"convergence-and-boost",convergenceAndBoostOk},
        {ValidationSection::AlgebraicIdentity,"production-geometry",productionGeometryOk},
        {ValidationSection::NumericalRegression,"branch-coupling",branchCouplingOk},
        {ValidationSection::AlgebraicIdentity,"yee-and-esirkepov",yeeAndEsirkepovOk},
        {ValidationSection::NumericalRegression,"shared-classical-engine",sharedClassicalEngineOk},
        {ValidationSection::PhysicalDomain,"reaction-models",reactionModelsOk},
        {ValidationSection::AlgebraicIdentity,"coherent-magnetic-dipole",coherentMagneticDipoleOk},
        {ValidationSection::Convergence,"far-field-convergence",farFieldConvergenceOk},
        {ValidationSection::Convergence,"larmor-normalization",larmorNormalizationOk},
        {ValidationSection::IndependentBalance,"long-horizon-radiative-balance",longHorizonBalanceOk},
        {ValidationSection::AlgebraicIdentity,"role-routing",roleRoutingOk},
        {ValidationSection::AlgebraicIdentity,"bmt-precession-invariant",bmtPrecessionOk},
        {ValidationSection::AlgebraicIdentity,"coupled-secular-spin-orbit",secularSpinOrbitIdentityOk},
        {ValidationSection::Convergence,"coupled-secular-spin-orbit-convergence",secularSpinOrbitConvergenceOk},
        {ValidationSection::Convergence,"secular-eccentric-orbit",secularEccentricOrbitOk},
        {ValidationSection::AlgebraicIdentity,"quantized-radiation-gating",quantizedRadiationOk},
        {ValidationSection::Convergence,"trajectory-convergence",trajectoryConvergenceOk},
        {ValidationSection::NumericalRegression,"state-integrity",stateIntegrityOk},
        {ValidationSection::NumericalRegression,"censoring-model",censoringModelOk},
        {ValidationSection::NumericalRegression,"visual-censoring-semantics",visualCensoringSemanticsOk},
        {ValidationSection::NumericalRegression,"adaptive-depth-rejection",adaptiveDepthRejectionOk},
        {ValidationSection::Convergence,"causal-startup",causalStartupOk},
        {ValidationSection::Convergence,"history-construction-sensitivity",historyConstructionSensitivityOk},
        {ValidationSection::PhysicalDomain,"cutoff-and-regularization-scan",parameterSensitivityScanOk},
        {ValidationSection::Convergence,"retarded-interpolation",retardedInterpolationOk},
        {ValidationSection::NumericalRegression,"short-range-regularization",shortRangeRegularizationOk},
        {ValidationSection::NumericalRegression,"interaction-impact-profile",impactParameterProfileOk},
        {ValidationSection::NumericalRegression,"stochastic-distributions",stochasticStatisticsOk}
    }};
    const auto sectionName=[](ValidationSection section) {
        switch(section) {
            case ValidationSection::AlgebraicIdentity:
                return "algebraic identities";
            case ValidationSection::NumericalRegression:
                return "numerical regressions";
            case ValidationSection::Convergence:
                return "convergence checks";
            case ValidationSection::IndependentBalance:
                return "independent balances";
            case ValidationSection::PhysicalDomain:
                return "physical-domain guards";
        }
        return "unknown";
    };
    constexpr std::array sections{
        ValidationSection::AlgebraicIdentity,
        ValidationSection::NumericalRegression,
        ValidationSection::Convergence,
        ValidationSection::IndependentBalance,
        ValidationSection::PhysicalDomain};
    int failedChecks=0;
    std::cout<<"\nCategorized validation verdict\n";
    for(const ValidationSection section:sections) {
        int passed=0,total=0;
        for(const ValidationCheck& check:regressionChecks) {
            if(check.section!=section) continue;
            ++total;
            if(check.passed) ++passed;
        }
        std::cout<<"  "<<sectionName(section)<<": "<<passed<<'/'<<total
                 <<' '<<(passed==total?"PASS":"FAIL")<<'\n';
        for(const ValidationCheck& check:regressionChecks) {
            if(check.section==section&&!check.passed) {
                ++failedChecks;
                std::cout<<"    FAILED: "<<check.name<<'\n';
            }
        }
    }
    std::cout<<"\nInformational diagnostics (no pass/fail criterion)\n"
             <<"  BMT vs effective-field comparison\n"
             <<"  flux-normalized reaction residual across heavy pairs\n"
             <<"  pair-field balance identities closed by boundField*\n"
             <<"  performance, memory and wall-clock estimates\n";
    std::cout<<"\nEnforced checks:     "
             <<(regressionChecks.size()-static_cast<std::size_t>(failedChecks))
             <<'/'<<regressionChecks.size()<<'\n'
             <<"validation verdict: "<<(failedChecks==0?"PASS":"FAIL")<<'\n'
             <<"scope: regression and declared-domain validation; "
                "not a QED validation\n";
    return failedChecks==0?0:1;
}
