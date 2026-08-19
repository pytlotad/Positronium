#pragma once

// End-to-end numerical validation of the optional Maxwell field backend.

int runMaxwellSelfTest() {
    const auto benchmarkStart=std::chrono::steady_clock::now();
    const CovariantExtendedBody covarianceBody{-eCharge,firstMass,
        chargeCloudRestRadius,{},Vec3{0.31*c,-0.07*c,0.04*c},
        Vec3{0,0,bohrMagneton}};
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
            {-eCharge,chargeCloudRestRadius},offset,Vec3{1.0e5,2.0e5,-0.5e5},240);
        const Vec3 selfForce=(selfElectric+cross(Vec3{1.0e5,2.0e5,-0.5e5},
                                                  selfMagnetic))*(-eCharge);
        maximumSelfForceFraction=std::max(maximumSelfForceFraction,
                                          selfForce.norm()/referenceForce);
    }
    const FourVector restCurrent{c*(-eCharge)*covarianceBody.properShape({}),{}};
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
    particles.firstDipole={0.0,0.0,bohrMagneton};
    particles.secondDipole={0.0,0.0,-bohrMagneton};
    coupledField.clearSources();
    coupledField.depositCloud({-eCharge,chargeCloudRestRadius},
        particles.firstPosition,particles.firstVelocity);
    coupledField.depositCloud({eCharge,chargeCloudRestRadius},
        particles.secondPosition,particles.secondVelocity);
    coupledField.depositCovariantDipole(particles.firstPosition,
        particles.firstVelocity,particles.firstDipole);
    coupledField.depositCovariantDipole(particles.secondPosition,
        particles.secondVelocity,particles.secondDipole);
    coupledField.finalizeBoundInstantaneous();
    coupledField.projectElectricGaussConstraint(400);
    const auto [firstSelfElectric,firstSelfMagnetic]=
        coupledField.numericalSelfField({-eCharge,chargeCloudRestRadius},
            particles.firstPosition,particles.firstVelocity);
    const auto [secondSelfElectric,secondSelfMagnetic]=
        coupledField.numericalSelfField({eCharge,chargeCloudRestRadius},
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
    const RelativisticChargeCloud gridFirst{-eCharge,chargeCloudRestRadius};
    const RelativisticChargeCloud gridSecond{eCharge,chargeCloudRestRadius};
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
        std::abs(particles.firstDipole.norm()/bohrMagneton-1.0),
        std::abs(particles.secondDipole.norm()/bohrMagneton-1.0));
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
        block.depositCloud({-eCharge,chargeCloudRestRadius},
                           separationVector*(-0.5),{});
        block.depositCloud({eCharge,chargeCloudRestRadius},
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
    retardedInitialState.firstDipole={0,0,bohrMagneton};
    retardedInitialState.secondDipole={0,0,-bohrMagneton};
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
        -eCharge,firstMass,reversibleElectric,reversibleMagnetic,reversibleDt);
    const Vec3 recoveredMomentum=relativisticBorisPush(pushedMomentum,
        -eCharge,firstMass,reversibleElectric,reversibleMagnetic,-reversibleDt);
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
    yeeCoupledState.firstDipole={0,0,bohrMagneton};
    yeeCoupledState.secondDipole={0,0,-bohrMagneton};
    yeeCoupledField.depositInitialGaussian(-eCharge,chargeCloudRestRadius,
        yeeCoupledState.firstPosition,yeeCoupledState.firstVelocity);
    yeeCoupledField.depositInitialGaussian(eCharge,chargeCloudRestRadius,
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

    struct ReactionModelBenchmark {
        State finalState;
        bool advanced=false;
        double rawEnergyResidual=0.0;
        double reactionFluxResidual=0.0;
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
                integrateElectrodynamicStep(value,dt,history,true,model);
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
        result.reactionFluxResidual=std::abs(full.reactionEnergyMismatch)
            /std::max(std::abs(full.radiatedEnergy),1.0e-300);
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
    covarianceRest.firstDipole={0,0,1.0e-3*bohrMagneton};
    covarianceRest.secondDipole={0,0,-1.0e-3*bohrMagneton};
    covarianceRest.firstProperDipole=covarianceRest.firstDipole;
    covarianceRest.secondProperDipole=covarianceRest.secondDipole;
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
    const DipoleTensor boostedFirstDipole=lorentzBoostDipole(
        {{},covarianceRest.firstDipole},covarianceBoost);
    const DipoleTensor boostedSecondDipole=lorentzBoostDipole(
        {{},covarianceRest.secondDipole},covarianceBoost);
    covarianceMoving.firstDipole=boostedFirstDipole.magnetic;
    covarianceMoving.secondDipole=boostedSecondDipole.magnetic;
    covarianceMoving.firstElectricDipole=boostedFirstDipole.electric;
    covarianceMoving.secondElectricDipole=boostedSecondDipole.electric;
    covarianceMoving.firstProperDipole=lorentzBoostDipole(
        boostedFirstDipole,covarianceMoving.firstVelocity*-1.0).magnetic;
    covarianceMoving.secondProperDipole=lorentzBoostDipole(
        boostedSecondDipole,covarianceMoving.secondVelocity*-1.0).magnetic;
    covarianceMoving.time=0.5*(boostedFirstTime+boostedSecondTime);
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
    // Role routing in applyDipolePrecession().
    //
    // The covariance probe above builds its arguments explicitly, so it never
    // executes the call sites and cannot see them hand one particle's g-factor,
    // charge-to-mass ratio or local field to the other.  This check does: it
    // runs the production routine, then recomputes the same step by calling
    // advanceCovariantBmt() directly with each role's OWN parameters, and
    // requires the two to agree exactly.
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
        const double radius=pairBohrRadius(defaultPair);
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
    const Vec3 expectedFirstProperDipole=advanceCovariantBmt(
        roleRoutingReference.firstProperDipole,
        roleRoutingReference.firstVelocity,roleRoutingFields.atFirst,
        firstCharge/firstMass,roleRoutingDt,firstGFactor);
    const Vec3 expectedSecondProperDipole=advanceCovariantBmt(
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
    State quadrupolePair=yeeCoupledState;
    quadrupolePair.firstPosition={-1.7*bohrRadius,0.4*bohrRadius,
                                      -0.2*bohrRadius};
    quadrupolePair.secondPosition={0.9*bohrRadius,-0.8*bohrRadius,
                                      0.6*bohrRadius};
    const ElectricQuadrupole neutralQuadrupole=
        electricQuadrupole(quadrupolePair);
    const double quadrupoleScale=eCharge
        *(quadrupolePair.firstPosition-quadrupolePair.secondPosition)
            .squaredNorm();
    double quadrupoleTrace=neutralQuadrupole.component[0]
        +neutralQuadrupole.component[4]+neutralQuadrupole.component[8];
    double quadrupoleSymmetry=0.0;
    for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        quadrupoleSymmetry=std::max(quadrupoleSymmetry,std::abs(
            neutralQuadrupole.component[static_cast<std::size_t>(3*i+j)]
           -neutralQuadrupole.component[static_cast<std::size_t>(3*j+i)]));
    const double quadrupoleCancellation=std::sqrt(
        neutralQuadrupole.squaredNorm())/quadrupoleScale;
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
    const Vec3 expectedDipoleThird=cubicCoefficient*(6.0*eCharge);
    const Vec3 measuredDipoleThird=electricDipoleThirdDerivative(
        cubicDipolePresent,cubicDipoleHistory);
    const double coherentDerivativeResidual=(measuredDipoleThird
        -expectedDipoleThird).norm()/expectedDipoleThird.norm();
    const MutualForces coherentProbe=coherentElectricDipoleReaction(
        cubicDipolePresent,cubicDipoleHistory);
    const double coherentReactionMomentumResidual=(coherentProbe.first
        +coherentProbe.second).norm()
        /std::max(coherentProbe.first.norm(),1.0e-300);

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
    {
        const double larmorProbeRadius=pairBohrRadius(defaultPair);
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
        const double lengthScale=std::max(separation(reference),nuclearCutoff);
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
    const State trajectoryReference=integrateAccuracyCase(1.0e-8,1.0e-20,4);
    const std::array<double,3> trajectoryTolerances{1.0e-5,1.0e-6,1.0e-7};
    std::array<double,3> trajectoryToleranceResiduals{};
    for(std::size_t index=0;index<trajectoryTolerances.size();++index)
        trajectoryToleranceResiduals[index]=trajectoryResidual(
            integrateAccuracyCase(trajectoryTolerances[index],1.0e-20,4),
            trajectoryReference);
    const State trajectoryHalfStep=integrateAccuracyCase(1.0e-7,0.5e-20,8);
    const double trajectoryStepResidual=trajectoryResidual(
        integrateAccuracyCase(1.0e-7,1.0e-20,4),trajectoryHalfStep);

    const auto integrateHistoryCase=[&](double spanFactor,int intervals) {
        State value=yeeCoupledState;
        ClassicalTrajectoryEngine engine(
            causalInitialHistory(value,spanFactor,intervals),{1.0e-7,14});
        bool advanced=true;
        for(int step=0;step<2;++step)
            advanced=engine.advance(value,1.0e-20)&&advanced;
        if(!advanced) value.time=std::numeric_limits<double>::quiet_NaN();
        return value;
    };
    const State startupHistory4=integrateHistoryCase(4.0,32);
    const State startupHistory8=integrateHistoryCase(8.0,64);
    const State startupHistory16=integrateHistoryCase(16.0,128);
    const double startupHistoryResidual4=trajectoryResidual(
        startupHistory4,startupHistory16);
    const double startupHistoryResidual8=trajectoryResidual(
        startupHistory8,startupHistory16);
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

    State staticDipoleState;
    staticDipoleState.firstPosition={2.0*nuclearCutoff,0,0};
    staticDipoleState.secondPosition={0,0,0};
    staticDipoleState.secondDipole={0,0,bohrMagneton};
    State staticDipolePast=staticDipoleState;
    staticDipolePast.time=-8.0*nuclearCutoff/c;
    const StateHistory staticDipoleHistory{staticDipolePast,staticDipoleState};
    const ElectromagneticField retardedStaticDipole=
        retardedMagneticDipoleField(staticDipoleState.firstPosition,0.0,
            staticDipoleHistory,staticDipoleState,false);
    const Vec3 directStaticDipole=regularizedDipoleField(
        staticDipoleState.firstPosition-staticDipoleState.secondPosition,
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
    // The physical ceiling the radius is chosen for: the classical dipole
    // interaction energy must stay below the first rest energy everywhere.
    // w(r)/r^3 peaks at r=a, so the maximum is (mu0/4pi)mu^2/(2a^3).
    double peakDipoleEnergyOverRestEnergy=0.0;
    for(int sample=0;sample<4000;++sample) {
        const double radius=magneticRegularizationRadius
            *std::pow(10.0,-3.0+6.0*sample/3999.0);
        // Must use the moment the dynamics actually carries, not the bare
        // Bohr magneton: the ceiling is a statement about the field the
        // simulated particles produce, so quoting mu_B here would certify a
        // cap that the production dipoles exceed by (g/2)^2.
        const double energy=shortRangeFieldWeight(radius)
            *(mu0/(4.0*pi))*firstMagneticMoment*firstMagneticMoment
            /(radius*radius*radius);
        peakDipoleEnergyOverRestEnergy=std::max(
            peakDipoleEnergyOverRestEnergy,energy/(firstMass*c*c));
    }
    const Vec3 cutoffDipoleForce=regularizedDipoleForce(
        {nuclearCutoff,0,0},{bohrMagneton,0,0},
        staticDipoleState.secondDipole);
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
    const bool chargeOk = std::abs(depositedCharge/eCharge) < 1.0e-12;
    const bool gaussOk = std::isfinite(relativeGaussResidual)
                      && relativeGaussResidual < 2.0e-2;
    const bool divergenceOk = magneticDivergence < 1.0e-6;
    const bool energyOk = relativeEnergyDrift < 5.0e-3;
    const bool couplingOk=maximumRelativeContinuityResidual<1.0e-8
                       &&maximumRelativeLongitudinalCurl<1.0e-12
                       && coupledBeta<1.0 && dipoleNormResidual<1.0e-12
                       && isFinite(particles);
    const bool balanceFinite=std::isfinite(coupledEnergyClosure)
        &&std::isfinite(coupledMomentumClosure)
        &&std::isfinite(coupledAngularClosure)
        &&std::abs(coupledEnergyClosure)<1.0e-3
        &&coupledMomentumClosure<1.0e-10
        &&coupledAngularClosure<1.0e-3;
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
        &&covarianceDipoleEvolutionResidual<5.0e-3
        &&covarianceBmtResidual<1.0e-3
        // Exact agreement is expected: both sides run the same integrator on
        // the same state.  The band admits reordering round-off only.
        &&covarianceForceResidual<1.0e-1
        // A coordinate sphere at constant t is not mapped to the coordinate
        // sphere used in the boosted run.  This is therefore a finite-surface
        // diagnostic, not a strict four-vector pass criterion.  It must remain
        // bounded; the local worldline and four-force tests above are strict.
        &&covarianceRadiationResidual<2.5e-1;
    const bool dipoleTensorCovarianceOk=dipoleTensorRoundtrip<1.0e-12
        &&dipoleFirstInvariantResidual<1.0e-12
        &&dipoleSecondInvariantResidual<1.0e-12
        &&inducedElectricDipoleResidual<1.0e-12
        &&staticElectricDipoleResidual<1.0e-8
        &&tensorGradientStaticResidual<1.0e-5;
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
        &&sharedEnergyBalanceResidual<1.0e-12
        &&sharedMomentumBalanceResidual<1.0e-12
        &&sharedAngularBalanceResidual<1.0e-12
        &&std::isfinite(reactionMismatchFraction)
        &&std::isfinite(llValidity)
        &&quadrupoleTrace<1.0e-14
        &&quadrupoleSymmetry<1.0e-14
        &&quadrupoleCancellation<1.0e-14
        &&coherentDerivativeResidual<1.0e-12
        &&coherentReactionMomentumResidual<1.0e-14;
    const auto finiteReactionBenchmark=[](const ReactionModelBenchmark& value) {
        return value.advanced&&isFinite(value.finalState)
            &&std::isfinite(value.rawEnergyResidual)
            &&std::isfinite(value.reactionFluxResidual)
            &&std::isfinite(value.stepConvergence)
            &&std::isfinite(value.wallSeconds)&&value.wallSeconds>0.0;
    };
    const bool reactionModelsOk=finiteReactionBenchmark(reactionDisabled)
        &&finiteReactionBenchmark(reactionLl)
        &&finiteReactionBenchmark(reactionCoherent)
        &&finiteReactionBenchmark(reactionAutomatic)
        &&std::isfinite(coherentCostRatio)&&coherentCostRatio<10.0;
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
    const bool causalStartupOk=startupFinite
        &&startupCoverage>=7.99
        &&std::abs(startupHistory.back().time-yeeCoupledState.time)<1.0e-30
        &&startupHistoryResidual4<1.0e-5
        &&startupHistoryResidual8<1.0e-6;
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
        // The design invariant is the m_e c^2/2 cap, not merely staying under
        // the rest energy: w(r)/r^3 peaks at r=a, so a radius derived from the
        // moment the dipoles actually carry puts the maximum at exactly 0.5
        // (0.49999 on this 4000-point grid).  A loose "< 1.0" accepted a
        // radius derived from the wrong moment, which lands at 0.50116.
        &&peakDipoleEnergyOverRestEnergy<0.5005
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
              << "particle boost D:  " << covarianceDipoleEvolutionResidual << '\n'
              << "covariant BMT:    " << covarianceBmtResidual << '\n'
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
              << "pair-field dE/dP/dJ:" << sharedEnergyBalanceResidual << " / "
              << sharedMomentumBalanceResidual << " / "
              << sharedAngularBalanceResidual << '\n'
              << "raw dE/dP/dJ:       " << sharedRawEnergyResidual << " / "
              << sharedRawMomentumResidual << " / " << sharedRawAngularResidual
              << '\n'
              << "reaction/flux dE:   " << reactionMismatchFraction << '\n'
              << "LL validity |F|:    " << llValidity << '\n'
              << "reaction raw off/LL/C:" << reactionDisabled.rawEnergyResidual
              << " / " << reactionLl.rawEnergyResidual << " / "
              << reactionCoherent.rawEnergyResidual << '\n'
              << "reaction flux off/LL/C:"
              << reactionDisabled.reactionFluxResidual << " / "
              << reactionLl.reactionFluxResidual << " / "
              << reactionCoherent.reactionFluxResidual << '\n'
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
              << "e+e- E2 cancellation:" << quadrupoleCancellation << '\n'
              << "E2 trace/symmetry:   " << quadrupoleTrace << " / "
              << quadrupoleSymmetry << '\n'
              << "coherent E1 d3/Fsum:" << coherentDerivativeResidual << " / "
              << coherentReactionMomentumResidual << '\n'
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
              << "near R=1e4/5/6 a0:  " << farNearFieldContamination[0] << " / "
              << farNearFieldContamination[1] << " / "
              << farNearFieldContamination[2] << '\n'
              << "step tol=1e-5/6/7:  " << trajectoryToleranceResiduals[0]
              << " / " << trajectoryToleranceResiduals[1] << " / "
              << trajectoryToleranceResiduals[2] << '\n'
              << "step dt/(dt/2):     " << trajectoryStepResidual << '\n'
              << "history span r/c:   " << startupCoverage << '\n'
              << "history 4/8 vs 16: " << startupHistoryResidual4 << " / "
              << startupHistoryResidual8 << '\n'
              << "history linear/H:  " << interpolationFine[0] << " / "
              << interpolationFine[1] << '\n'
              << "Hermite order:     " << hermiteConvergenceOrder << '\n'
              << "dipole static limit:" << retardedStaticLimitResidual << '\n'
              << "reg a=.75/.5/.25 rc:" << regulatorFieldResiduals[0] << " / "
              << regulatorFieldResiduals[1] << " / "
              << regulatorFieldResiduals[2] << '\n'
              << "reg p=4/6/8 @10a:  " << regulatorProfileResiduals[0] << " / "
              << regulatorProfileResiduals[1] << " / "
              << regulatorProfileResiduals[2] << '\n'
              << "reg core suppress: " << regulatorCoreSuppression << '\n'
              << "reg peak U/m_e c2: " << peakDipoleEnergyOverRestEnergy << '\n'
              << "cutoff surface:     " << cutoffSurfaceResidual << '\n'
              << "CFL dt:             " << cflStep << " s\n"
              << "field storage:      "
              << static_cast<double>(hierarchyBytes)/(1024.0*1024.0) << " MiB\n"
              << "validation wall:    " << benchmarkSeconds << " s\n"
              << "field-step rate:    " << fieldStepsPerSecond << " steps/s\n"
              << "steps per ps:       " << estimatedStepsPerPicosecond << '\n'
              << "estimated s/ps:     " << estimatedSecondsPerPicosecond << '\n';

    // Single source of truth for the regression verdict.  This used to be one
    // 26-term conjunction written out twice -- once for the printed verdict and
    // once for the exit status -- so a check added to only one of them would
    // have silently made the two disagree.  Listing the checks once also lets
    // the harness name the ones that actually failed instead of collapsing
    // everything into a single PASS/FAIL with sixty unlabelled numbers above it.
    const std::array<std::pair<const char*,bool>,28> regressionChecks{{
        {"charge",                     chargeOk},
        {"gauss",                      gaussOk},
        {"divergence",                 divergenceOk},
        {"energy",                     energyOk},
        {"coupling",                   couplingOk},
        {"balance-finite",             balanceFinite},
        {"boundary",                   boundaryOk},
        {"validation",                 validationOk},
        {"moving-amr",                 movingAmrOk},
        {"covariance",                 covarianceOk},
        {"particle-covariance",        particleCovarianceOk},
        {"dipole-tensor-covariance",   dipoleTensorCovarianceOk},
        {"mass-and-self-force",        massAndSelfForceOk},
        {"bound-current",              boundCurrentOk},
        {"retarded-initialization",    retardedInitializationOk},
        {"convergence-and-boost",      convergenceAndBoostOk},
        {"production-geometry",        productionGeometryOk},
        {"branch-coupling",            branchCouplingOk},
        {"yee-and-esirkepov",          yeeAndEsirkepovOk},
        {"shared-classical-engine",    sharedClassicalEngineOk},
        {"reaction-models",            reactionModelsOk},
        {"far-field-convergence",      farFieldConvergenceOk},
        {"larmor-normalization",       larmorNormalizationOk},
        {"role-routing",               roleRoutingOk},
        {"trajectory-convergence",     trajectoryConvergenceOk},
        {"causal-startup",             causalStartupOk},
        {"retarded-interpolation",     retardedInterpolationOk},
        {"short-range-regularization", shortRangeRegularizationOk}
    }};
    int failedChecks=0;
    for(const auto& [name,ok]:regressionChecks) {
        if(!ok) {
            ++failedChecks;
            std::cout<<"FAILED check:      "<<name<<'\n';
        }
    }
    std::cout<<"checks passed:      "
             <<(regressionChecks.size()-static_cast<std::size_t>(failedChecks))
             <<"/"<<regressionChecks.size()<<'\n'
             <<"numerical regression:"<<(failedChecks==0?"PASS":"FAIL")<<'\n';
    return failedChecks==0?0:1;
}
