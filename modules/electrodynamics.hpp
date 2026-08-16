#pragma once

// Particle-particle interactions, radiation bookkeeping and the shared
// relativistic integration step. This header is included inside the
// implementation namespace after the field-solver types are defined.

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
struct ParticleFieldTotals {
    double energy=0.0;
    Vec3 momentum, angularMomentum;
};

ParticleFieldTotals particleFieldTotals(const State& state,
                                         const MaxwellBlock& field) {
    const Vec3 electronMomentum=momentum(state.electronVelocity,electronMass);
    const Vec3 positronMomentum=momentum(state.positronVelocity,positronMass);
    const MaxwellVolumeIntegrals fieldTotals=field.volumeIntegrals();
    constexpr double electronGyromagneticRatio=-eCharge/(2.0*electronMass);
    constexpr double positronGyromagneticRatio=eCharge/(2.0*positronMass);
    // Rest energies are constant and omitted from the diagnostic difference.
    return {(gamma(state.electronVelocity)-1.0)*electronMass*c*c
          +(gamma(state.positronVelocity)-1.0)*positronMass*c*c
          +fieldTotals.energy+state.dipoleConstraintEnergy,
            electronMomentum+positronMomentum+fieldTotals.momentum,
            cross(state.electronPosition,electronMomentum)
          +cross(state.positronPosition,positronMomentum)
          +state.electronDipole/electronGyromagneticRatio
          +state.positronDipole/positronGyromagneticRatio
          +fieldTotals.angularMomentum};
}
#endif

// Effective field in the Thomas-BMT equation, expressed in laboratory time:
// d(mu)/dt = (q/m) mu x B_BMT.  The classical orbital value g=1 is retained;
// choosing the measured electron g factor would introduce a QED input.
Vec3 thomasBmtEffectiveField(const Vec3& velocity,
                             const ElectromagneticField& field) {
    constexpr double classicalGFactor = 1.0;
    constexpr double anomaly = 0.5 * (classicalGFactor - 2.0);
    const double relativisticGamma = gamma(velocity);
    const Vec3 beta = velocity / c;
    return field.magnetic * (anomaly + 1.0 / relativisticGamma)
         - beta * (anomaly * relativisticGamma
                   / (relativisticGamma + 1.0)
                   * dot(beta, field.magnetic))
         - cross(beta, field.electric)
             * ((anomaly + 1.0 / (relativisticGamma + 1.0)) / c);
}

Vec3 relativisticAcceleration(const Vec3& velocity, const Vec3& force, double mass) {
    const double velocityForce = dot(velocity, force);
    return (force - velocity * (velocityForce / (c*c))) / (gamma(velocity) * mass);
}

Vec3 lorentzForce(double charge, const Vec3& velocity, const ElectromagneticField& field) {
    return (field.electric + cross(velocity, field.magnetic)) * charge;
}

struct MutualForces { Vec3 electron, positron; };

struct FieldFluxRates {
    double energy = 0.0;
    Vec3 momentum, angularMomentum;
};

struct FarFieldSampling {
    int directionCount = 50;
    double controlRadius = 1.0e6 * bohrRadius;
    bool radiationFieldOnly = false;
};

struct SphereQuadraturePoint {
    Vec3 direction;
    double solidAngleWeight=0.0;
};

// Degree-11, 50-node Lebedev rule.  Its octahedral symmetry integrates the
// low electromagnetic multipoles without the small preferred-axis bias of a
// finite Fibonacci lattice.  We retain the latter as a diagnostic fallback
// for non-tabulated direction counts used by convergence tests.
std::vector<SphereQuadraturePoint> sphereQuadrature(int directionCount) {
    std::vector<SphereQuadraturePoint> result;
    if(directionCount==50) {
        const auto add=[&](double x,double y,double z,double normalizedWeight) {
            result.push_back({{x,y,z},4.0*pi*normalizedWeight});
        };
        constexpr double axisWeight=0.012698412698412698;
        for(int axis=0;axis<3;++axis) for(double sign:{-1.0,1.0}) {
            Vec3 n;
            if(axis==0) n.x=sign; else if(axis==1) n.y=sign; else n.z=sign;
            add(n.x,n.y,n.z,axisWeight);
        }
        constexpr double edge=0.70710678118654752440;
        constexpr double edgeWeight=0.022574955908289241;
        for(int zeroAxis=0;zeroAxis<3;++zeroAxis)
            for(double first:{-edge,edge}) for(double second:{-edge,edge}) {
                Vec3 n;
                if(zeroAxis==0) { n.y=first; n.z=second; }
                else if(zeroAxis==1) { n.x=first; n.z=second; }
                else { n.x=first; n.y=second; }
                add(n.x,n.y,n.z,edgeWeight);
            }
        constexpr double corner=0.57735026918962576451;
        constexpr double cornerWeight=0.02109375;
        for(double x:{-corner,corner}) for(double y:{-corner,corner})
            for(double z:{-corner,corner}) add(x,y,z,cornerWeight);
        constexpr double small=0.30151134457776362265;
        constexpr double large=0.90453403373329086794;
        constexpr double mixedWeight=0.020173335537918871;
        for(int largeAxis=0;largeAxis<3;++largeAxis)
            for(double sx:{-1.0,1.0}) for(double sy:{-1.0,1.0})
                for(double sz:{-1.0,1.0}) {
                    Vec3 n{sx*small,sy*small,sz*small};
                    if(largeAxis==0) n.x=sx*large;
                    else if(largeAxis==1) n.y=sy*large;
                    else n.z=sz*large;
                    add(n.x,n.y,n.z,mixedWeight);
                }
        return result;
    }
    if(directionCount<1) return result;
    constexpr double goldenAngle=pi*(3.0-2.2360679774997896964);
    const double weight=4.0*pi/directionCount;
    result.reserve(static_cast<std::size_t>(directionCount));
    for(int index=0;index<directionCount;++index) {
        const double z=1.0-2.0*(index+0.5)/directionCount;
        const double transverse=std::sqrt(std::max(0.0,1.0-z*z));
        const double azimuth=goldenAngle*index;
        result.push_back({{transverse*std::cos(azimuth),
                           transverse*std::sin(azimuth),z},weight});
    }
    return result;
}

ElectromagneticField farZoneChargeField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool electron,double charge,bool radiationFieldOnly=false) {
    double emissionTime=wavefrontTime;
    ChargeKinematics source=historicalCharge(
        history,present,electron,emissionTime);
    // Far-zone light cone: t_emit=t_wave+n.(x_source-x_centre)/c.
    // Iterating this small correction avoids cancellation between R/c terms.
    for(int iteration=0;iteration<5;++iteration) {
        const double refined=wavefrontTime
            +dot(normal,source.position-centre)/c;
        if(std::abs(refined-emissionTime)<=1.0e-30
            +1.0e-14*std::abs(emissionTime)) {
            emissionTime=refined;
            break;
        }
        emissionTime=refined;
        source=historicalCharge(history,present,electron,emissionTime);
    }
    source=historicalCharge(history,present,electron,emissionTime);
    const Vec3 displacement=observationPosition-source.position;
    const double distance=displacement.norm();
    const Vec3 direction=displacement/distance;
    const Vec3 beta=source.velocity/c;
    const double kappa=std::max(1.0e-12,1.0-dot(direction,beta));
    const Vec3 velocityField=(direction-beta)*((1.0-beta.squaredNorm())
        /(kappa*kappa*kappa*distance*distance));
    const Vec3 accelerationField=cross(direction,
        cross(direction-beta,source.acceleration))
        /(c*c*kappa*kappa*kappa*distance);
    const Vec3 electric=(radiationFieldOnly?accelerationField
        :velocityField+accelerationField)*(coulomb*charge);
    return {electric,cross(direction,electric)/c};
}

FieldFluxRates electromagneticFieldFluxRates(
    const State& state, const StateHistory& history,
    FarFieldSampling sampling={}) {
    if(sampling.directionCount<1||!(sampling.controlRadius>0.0)
        ||!std::isfinite(sampling.controlRadius)) return {};
    const std::vector<SphereQuadraturePoint> quadrature=
        sphereQuadrature(sampling.directionCount);
    const Vec3 centre=(state.electronPosition+state.positronPosition)*0.5;
    const double sourceExtent=std::max(
        (state.electronPosition-centre).norm(),
        (state.positronPosition-centre).norm());
    // The observation event is shifted back by the source radius so every
    // direction samples only the available causal history. Directional
    // retardation across the pair then retains E3, M2, toroidal terms and all
    // their interference without assigning convention-dependent pieces.
    const double wavefrontTime=state.time-sourceExtent/c;
    FieldFluxRates rates;
    for(const SphereQuadraturePoint& point:quadrature) {
        const Vec3 normal=point.direction;
        const Vec3 observationPosition = centre+normal*sampling.controlRadius;
        const ElectromagneticField electronField=farZoneChargeField(
            observationPosition,normal,wavefrontTime,centre,history,state,
            true,-eCharge,sampling.radiationFieldOnly);
        const ElectromagneticField positronField=farZoneChargeField(
            observationPosition,normal,wavefrontTime,centre,history,state,
            false,eCharge,sampling.radiationFieldOnly);
        const Vec3 electric = electronField.electric + positronField.electric;
        const Vec3 magnetic = electronField.magnetic + positronField.magnetic;
        const Vec3 poynting = cross(electric, magnetic) / mu0;
        const double areaWeight=sampling.controlRadius
            *sampling.controlRadius*point.solidAngleWeight;
        rates.energy += dot(poynting, normal) * areaWeight;

        // Outward momentum flux is -T.n for the Maxwell stress convention
        // T_ij=eps0(E_iE_j-E^2 delta_ij/2)+...
        const Vec3 stressOnNormal =
            (electric * dot(electric, normal)
             - normal * (0.5 * electric.squaredNorm())) * epsilon0
          + (magnetic * dot(magnetic, normal)
             - normal * (0.5 * magnetic.squaredNorm())) / mu0;
        const Vec3 momentumFlux = stressOnNormal * (-areaWeight);
        rates.momentum += momentumFlux;
        rates.angularMomentum += cross(observationPosition, momentumFlux);
    }
    return rates;
}

struct LocalElectromagneticFields {
    ElectromagneticField atElectron, atPositron;
};

double shortRangeFieldWeight(double distance,
    double regularizationRadius=magneticRegularizationRadius,
    double exponent=magneticRegularizationExponent) {
    if(!(distance>0.0)) return 0.0;
    const double ratio = regularizationRadius / distance;
    return 1.0/(1.0+std::pow(ratio,exponent));
}

struct MagneticRadialProfile {
    double vectorPotentialFactor;
    double firstDerivative;
    double secondDerivative;
};

struct RetardedDipoleKinematics {
    Vec3 position, velocity, moment, firstDerivative, secondDerivative,
         thirdDerivative;
};

struct RetardedElectricDipoleKinematics {
    Vec3 position,velocity,moment,firstDerivative,secondDerivative;
};

struct DipoleDerivatives { Vec3 electron, positron; };
DipoleDerivatives thomasBmtDipoleDerivatives(
    const State& s,const StateHistory& history);

State historicalState(const StateHistory& history, const State& present,
                      double time) {
    if (history.empty()) return present;
    const State& earliest = history.front();
    if (history.size() == 1) {
        State extrapolated = earliest;
        const double offset = time - earliest.time;
        extrapolated.electronPosition += earliest.electronVelocity*offset;
        extrapolated.positronPosition += earliest.positronVelocity*offset;
        extrapolated.time = time;
        return extrapolated;
    }
    const auto newer=std::lower_bound(history.begin(),history.end(),time,
        [](const State& state,double requestedTime) {
            return state.time<requestedTime;
        });
    if(newer==history.begin()) return *newer;
    if(newer!=history.end()) {
        const State& older=*std::prev(newer);
        return interpolateState(older,*newer,
            (time-older.time)/std::max(newer->time-older.time,1.0e-300));
    }
    const State& older=history.back();
    return interpolateState(older,present,
        (time-older.time)/std::max(present.time-older.time,1.0e-300));
}

// Largest step a backward stencil may use if it must stay inside the retained
// history.  historicalState() clamps to history.front() outside the retained
// window, so a stencil reaching past it silently differentiates constant data:
// the coefficients still sum to zero, but the result is no longer a derivative
// and is inflated by the missing powers of h.  appendStateHistory() keeps only
// 4*r/c, which at production step sizes is 3-6 nodes, while an unbounded
// 5-point third-derivative stencil asks for 8x the last node spacing.
// Measured on a bound para orbit, 99.9% of calls used to overrun the window.
//
// Returning zero means the requested derivative order is not resolvable from
// the retained history at all; callers must then report zero rather than
// noise, because the alternative is a plausible-looking wrong number.
// evaluationTime is the centre of the stencil, which for retarded sources is
// earlier than present.time; measuring the span from present.time instead
// would still let a retarded stencil run off the front of the deque.
double boundedDerivativeStep(const StateHistory& history,
                             double evaluationTime,
                             double requestedStep, int stencilReach) {
    if(history.empty()||stencilReach<1) return 0.0;
    const double span=evaluationTime-history.front().time;
    if(!(span>0.0)||!std::isfinite(span)) return 0.0;
    const double usable=std::min(requestedStep,
        span/static_cast<double>(stencilReach));
    return usable>0.0?usable:0.0;
}

RetardedDipoleKinematics historicalDipoleKinematics(
    const StateHistory& history, const State& present, bool sourceIsElectron,
    double time) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    // The widest branch below samples middle .. middle-4h.
    derivativeStep=boundedDerivativeStep(history,time,derivativeStep,4);
    if(!(derivativeStep>0.0)) {
        const State only=historicalState(history,present,time);
        return {sourceIsElectron?only.electronPosition:only.positronPosition,
                sourceIsElectron?only.electronVelocity:only.positronVelocity,
                sourceIsElectron?only.electronDipole:only.positronDipole,
                {},{},{}};
    }
    const State middle=historicalState(history,present,time);
    const auto moment=[&](const State& state) -> const Vec3& {
        return sourceIsElectron?state.electronDipole:state.positronDipole;
    };
    Vec3 first,second,third;
    if(time+derivativeStep>present.time) {
        const State before=historicalState(
            history,present,time-derivativeStep);
        const State twiceBefore=historicalState(
            history,present,time-2.0*derivativeStep);
        const State threeBefore=historicalState(
            history,present,time-3.0*derivativeStep);
        const State fourBefore=historicalState(
            history,present,time-4.0*derivativeStep);
        first=(moment(middle)*3.0-moment(before)*4.0+moment(twiceBefore))
             /(2.0*derivativeStep);
        second=(moment(middle)-moment(before)*2.0+moment(twiceBefore))
              /(derivativeStep*derivativeStep);
        third=(moment(middle)*5.0-moment(before)*18.0
              +moment(twiceBefore)*24.0-moment(threeBefore)*14.0
              +moment(fourBefore)*3.0)
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    } else if(time-derivativeStep<history.front().time) {
        const State after=historicalState(history,present,time+derivativeStep);
        const State twiceAfter=historicalState(
            history,present,time+2.0*derivativeStep);
        const State threeAfter=historicalState(
            history,present,time+3.0*derivativeStep);
        const State fourAfter=historicalState(
            history,present,time+4.0*derivativeStep);
        first=(moment(after)*4.0-moment(middle)*3.0-moment(twiceAfter))
             /(2.0*derivativeStep);
        second=(moment(twiceAfter)-moment(after)*2.0+moment(middle))
              /(derivativeStep*derivativeStep);
        third=(moment(middle)*-5.0+moment(after)*18.0
              -moment(twiceAfter)*24.0+moment(threeAfter)*14.0
              -moment(fourAfter)*3.0)
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    } else {
        const State before=historicalState(
            history,present,time-derivativeStep);
        const State after=historicalState(
            history,present,time+derivativeStep);
        const State twiceBefore=historicalState(
            history,present,time-2.0*derivativeStep);
        const State twiceAfter=historicalState(
            history,present,time+2.0*derivativeStep);
        first=(moment(after)-moment(before))/(2.0*derivativeStep);
        second=(moment(after)-moment(middle)*2.0+moment(before))
              /(derivativeStep*derivativeStep);
        third=(moment(twiceAfter)-moment(after)*2.0
              +moment(before)*2.0-moment(twiceBefore))
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    }
    return {sourceIsElectron?middle.electronPosition:middle.positronPosition,
            sourceIsElectron?middle.electronVelocity:middle.positronVelocity,
            moment(middle),first,second,third};
}

RetardedElectricDipoleKinematics historicalElectricDipoleKinematics(
    const StateHistory& history,const State& present,bool sourceIsElectron,
    double time) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    // The widest branch below samples middle .. middle-2h.
    derivativeStep=boundedDerivativeStep(history,time,derivativeStep,2);
    const auto sample=[&](double sampleTime) {
        const State state=historicalState(history,present,sampleTime);
        return sourceIsElectron?state.electronElectricDipole
                               :state.positronElectricDipole;
    };
    const State middle=historicalState(history,present,time);
    if(!(derivativeStep>0.0)) {
        return {sourceIsElectron?middle.electronPosition:middle.positronPosition,
                sourceIsElectron?middle.electronVelocity:middle.positronVelocity,
                sample(time),{},{}};
    }
    const Vec3 moment=sample(time);
    Vec3 first,second;
    if(time+derivativeStep>present.time) {
        const Vec3 before=sample(time-derivativeStep);
        const Vec3 twiceBefore=sample(time-2.0*derivativeStep);
        first=(moment*3.0-before*4.0+twiceBefore)/(2.0*derivativeStep);
        second=(moment-before*2.0+twiceBefore)
            /(derivativeStep*derivativeStep);
    } else if(time-derivativeStep<history.front().time) {
        const Vec3 after=sample(time+derivativeStep);
        const Vec3 twiceAfter=sample(time+2.0*derivativeStep);
        first=(after*4.0-moment*3.0-twiceAfter)/(2.0*derivativeStep);
        second=(twiceAfter-after*2.0+moment)
            /(derivativeStep*derivativeStep);
    } else {
        const Vec3 before=sample(time-derivativeStep);
        const Vec3 after=sample(time+derivativeStep);
        first=(after-before)/(2.0*derivativeStep);
        second=(after-moment*2.0+before)/(derivativeStep*derivativeStep);
    }
    return {sourceIsElectron?middle.electronPosition:middle.positronPosition,
            sourceIsElectron?middle.electronVelocity:middle.positronVelocity,
            moment,first,second};
}

struct DipoleRadiationReaction {
    double power=0.0;
    Vec3 momentumRate,angularMomentumRate;
    Vec3 electronTorque,positronTorque;
};

struct ElectricQuadrupole {
    // Symmetric tensor Q_ij = sum_a q_a(3 x_i x_j-r^2 delta_ij).
    std::array<double,9> component{};
    ElectricQuadrupole operator+(const ElectricQuadrupole& other) const {
        ElectricQuadrupole result;
        for(std::size_t i=0;i<component.size();++i)
            result.component[i]=component[i]+other.component[i];
        return result;
    }
    ElectricQuadrupole operator-(const ElectricQuadrupole& other) const {
        ElectricQuadrupole result;
        for(std::size_t i=0;i<component.size();++i)
            result.component[i]=component[i]-other.component[i];
        return result;
    }
    ElectricQuadrupole operator*(double factor) const {
        ElectricQuadrupole result;
        for(std::size_t i=0;i<component.size();++i)
            result.component[i]=component[i]*factor;
        return result;
    }
    double squaredNorm() const {
        double result=0.0;
        for(double value:component) result+=value*value;
        return result;
    }
};

ElectricQuadrupole electricQuadrupole(const State& state) {
    // Use the equal-mass centre as origin.  This keeps the truncated
    // multipole decomposition reproducible while the complete neutral-source
    // radiation remains origin independent.
    const Vec3 origin=(state.electronPosition+state.positronPosition)*0.5;
    ElectricQuadrupole result;
    const auto accumulate=[&](double charge,const Vec3& position) {
        const Vec3 x=position-origin;
        const std::array<double,3> coordinate{x.x,x.y,x.z};
        const double radiusSquared=x.squaredNorm();
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
            result.component[static_cast<std::size_t>(3*i+j)] += charge
                *(3.0*coordinate[static_cast<std::size_t>(i)]
                    *coordinate[static_cast<std::size_t>(j)]
                  -(i==j?radiusSquared:0.0));
    };
    accumulate(-eCharge,state.electronPosition);
    accumulate(eCharge,state.positronPosition);
    return result;
}

ElectricQuadrupole electricQuadrupoleThirdDerivative(
    const State& state,const StateHistory& history) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        8.0*(history.back().time-history[history.size()-2].time));
    // Five-point backward stencil: samples state.time .. state.time-4h.
    derivativeStep=boundedDerivativeStep(history,state.time,derivativeStep,4);
    if(!(derivativeStep>0.0)) return {};
    const auto at=[&](double offset) {
        return electricQuadrupole(
            historicalState(history,state,state.time+offset));
    };
    // Five-point backward derivative is causal and matches the stencil used
    // for radiation from the time-dependent magnetic moments.
    return (at(0.0)*5.0-at(-derivativeStep)*18.0
           +at(-2.0*derivativeStep)*24.0
           -at(-3.0*derivativeStep)*14.0
           +at(-4.0*derivativeStep)*3.0)
        *(1.0/(2.0*derivativeStep*derivativeStep*derivativeStep));
}

Vec3 electricDipoleMoment(const State& state) {
    return (state.positronPosition-state.electronPosition)*eCharge;
}

Vec3 electricDipoleThirdDerivativeAtStep(
    const State& state,const StateHistory& history,double derivativeStep) {
    // Five-point backward stencil: samples state.time .. state.time-4h.
    derivativeStep=boundedDerivativeStep(history,state.time,derivativeStep,4);
    if(!(derivativeStep>0.0)) return {};
    const auto at=[&](double offset) {
        return electricDipoleMoment(
            historicalState(history,state,state.time+offset));
    };
    return (at(0.0)*5.0-at(-derivativeStep)*18.0
           +at(-2.0*derivativeStep)*24.0
           -at(-3.0*derivativeStep)*14.0
           +at(-4.0*derivativeStep)*3.0)
        /(2.0*derivativeStep*derivativeStep*derivativeStep);
}

// Base step for the fine/coarse convergence probe in particleMultipoleRadiation.
// The coarse probe doubles it and its 5-point stencil then reaches 8h, so the
// base must fit eight times into the retained history for the two probes to be
// genuinely different rather than both saturating at the same clamped value.
double electricDipoleDerivativeStep(const StateHistory& history,
                                    const State& present) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    return boundedDerivativeStep(history,present.time,derivativeStep,8);
}

Vec3 electricDipoleThirdDerivative(const State& state,
                                   const StateHistory& history) {
    return electricDipoleThirdDerivativeAtStep(
        state,history,electricDipoleDerivativeStep(history,state));
}

MutualForces coherentElectricDipoleReaction(
    const State& state,const StateHistory& history) {
    // The common Abraham-Lorentz radiation field of a compact neutral
    // source.  Its work is -|p_ddot|^2/(6 pi eps0 c^3) up to the reversible
    // Schott derivative, so dynamics and the leading E1 flux use one model.
    const Vec3 reactionField=electricDipoleThirdDerivative(state,history)
        /(6.0*pi*epsilon0*c*c*c);
    return {reactionField*(-eCharge),reactionField*eCharge};
}

enum class ChargeRadiationReactionModel {
    automatic,
    coherentElectricDipole,
    individualLandauLifshitz,
    disabled
};

double electricQuadrupoleRadiatedPower(
    const State& state,const StateHistory& history) {
    const ElectricQuadrupole third=
        electricQuadrupoleThirdDerivative(state,history);
    // SI coefficient for Q_ij=sum q(3 x_i x_j-r^2 delta_ij).
    return third.squaredNorm()/(180.0*pi*epsilon0*c*c*c*c*c);
}

MutualForces individualLandauLifshitzSelfForces(
    const State& state, const MutualForces& external,
    const StateHistory& history);

DipoleRadiationReaction dipoleRadiationReaction(
    const State& state,const StateHistory& history) {
    constexpr double coefficient=mu0/(6.0*pi*c*c*c);
    const RetardedDipoleKinematics electron=historicalDipoleKinematics(
        history,state,true,state.time);
    const RetardedDipoleKinematics positron=historicalDipoleKinematics(
        history,state,false,state.time);
    DipoleRadiationReaction result;
    const auto accumulate=[&](const RetardedDipoleKinematics& dipole,
                              Vec3& reactionTorque) {
        const double power=coefficient*dipole.secondDerivative.squaredNorm();
        // Abraham-Lorentz reaction torque of a localized magnetic dipole.
        reactionTorque=cross(dipole.moment,dipole.thirdDerivative)*coefficient;
        result.power+=power;
        result.momentumRate+=dipole.velocity*(power/(c*c));
        result.angularMomentumRate+=reactionTorque*-1.0;
    };
    accumulate(electron,result.electronTorque);
    accumulate(positron,result.positronTorque);
    return result;
}

// One response object is the sole interface between particle dynamics and
// radiation bookkeeping.  The charge sector retains the relativistic
// Lienard-Wiechert angular flux and order-reduced Landau-Lifshitz force; the
// first explicit multipole is the time-dependent magnetic dipole.  Keeping
// forces, torques and outward flux together prevents the integrator from
// silently enabling only one side of a radiation channel.
struct ParticleMultipoleRadiation {
    MutualForces chargeReaction;
    Vec3 electronDipoleTorque, positronDipoleTorque;
    FieldFluxRates outwardFlux;
    double leadingElectricDipolePower = 0.0;
    double magneticDipolePower = 0.0;
    double electricQuadrupolePower = 0.0;
    double landauLifshitzValidity = 0.0;
    double coherentDerivativeConsistency = 0.0;
    double sourceCompactness = 0.0;
    double coherentWeight = 0.0;
    bool coherentSelected = false;
};

ParticleMultipoleRadiation particleMultipoleRadiation(
    const State& state, const MutualForces& externalForces,
    const StateHistory& history,bool computeOutwardFlux=true,
    ChargeRadiationReactionModel reactionModel=
        ChargeRadiationReactionModel::individualLandauLifshitz) {
    ParticleMultipoleRadiation result;
    // Only the automatic model consults the blending gates, and only the two
    // coherent models need the coherent reaction force itself.  Under the
    // default individual Landau-Lifshitz model both used to be evaluated and
    // thrown away, at a cost of four five-point history stencils per force
    // evaluation -- and there are six force evaluations per integration step.
    const bool needsCoherentReaction=
        reactionModel==ChargeRadiationReactionModel::coherentElectricDipole
        ||reactionModel==ChargeRadiationReactionModel::automatic;
    const bool needsBlendingGates=
        reactionModel==ChargeRadiationReactionModel::automatic;
    const MutualForces ll=individualLandauLifshitzSelfForces(
        state,externalForces,history);
    result.landauLifshitzValidity=std::max(
        ll.electron.norm()/std::max(externalForces.electron.norm(),1.0e-300),
        ll.positron.norm()/std::max(externalForces.positron.norm(),1.0e-300));
    if(computeOutwardFlux)
        result.outwardFlux=electromagneticFieldFluxRates(state,history);

    const Vec3 electronAcceleration = relativisticAcceleration(
        state.electronVelocity, externalForces.electron, electronMass);
    const Vec3 positronAcceleration = relativisticAcceleration(
        state.positronVelocity, externalForces.positron, positronMass);
    const Vec3 electricDipoleSecondDerivative =
        (positronAcceleration - electronAcceleration) * eCharge;
    result.leadingElectricDipolePower =
        electricDipoleSecondDerivative.squaredNorm()
        / (6.0*pi*epsilon0*c*c*c);
    if(needsBlendingGates) {
        result.electricQuadrupolePower =
            electricQuadrupoleRadiatedPower(state,history);
    }

    const DipoleRadiationReaction magnetic =
        dipoleRadiationReaction(state, history);
    result.electronDipoleTorque = magnetic.electronTorque;
    result.positronDipoleTorque = magnetic.positronTorque;
    result.magneticDipolePower = magnetic.power;
    if(needsBlendingGates) {
        const double derivativeStep=electricDipoleDerivativeStep(history,state);
        const Vec3 dipoleThirdFine=electricDipoleThirdDerivativeAtStep(
            state,history,derivativeStep);
        const Vec3 dipoleThirdCoarse=electricDipoleThirdDerivativeAtStep(
            state,history,2.0*derivativeStep);
        result.coherentDerivativeConsistency=
            (dipoleThirdFine-dipoleThirdCoarse).norm()
            /std::max(dipoleThirdFine.norm(),1.0e-300);
        const double dipoleNorm=std::max(electricDipoleMoment(state).norm(),
                                         eCharge*nuclearCutoff);
        const double angularRate=std::sqrt(
            electricDipoleSecondDerivative.norm()/dipoleNorm);
        result.sourceCompactness=angularRate*separation(state)/c;
        const double nonElectricPower=result.magneticDipolePower
                                      +result.electricQuadrupolePower;
        const auto decreasingGate=[](double value,double full,double zero) {
            return std::clamp((zero-value)/(zero-full),0.0,1.0);
        };
        const double dominance=result.leadingElectricDipolePower
            /std::max(nonElectricPower,1.0e-300);
        const double smoothGate=decreasingGate(
            result.coherentDerivativeConsistency,1.0e-2,5.0e-2);
        const double compactGate=decreasingGate(
            result.sourceCompactness,5.0e-2,1.0e-1);
        const double llGate=decreasingGate(
            result.landauLifshitzValidity,5.0e-3,1.0e-2);
        const double dominanceGate=std::clamp((dominance-10.0)/10.0,0.0,1.0);
        result.coherentWeight=smoothGate*compactGate*llGate*dominanceGate;
    }
    result.coherentSelected=
        reactionModel==ChargeRadiationReactionModel::coherentElectricDipole
        ||(reactionModel==ChargeRadiationReactionModel::automatic
           &&result.coherentWeight>0.0);
    if(reactionModel!=ChargeRadiationReactionModel::disabled) {
        if(!needsCoherentReaction) result.chargeReaction=ll;
        else {
            const MutualForces coherent=
                coherentElectricDipoleReaction(state,history);
            if(reactionModel
                ==ChargeRadiationReactionModel::coherentElectricDipole) {
                result.chargeReaction=coherent;
            } else {
                result.chargeReaction={
                    ll.electron*(1.0-result.coherentWeight)
                        +coherent.electron*result.coherentWeight,
                    ll.positron*(1.0-result.coherentWeight)
                        +coherent.positron*result.coherentWeight};
            }
        }
    }
    if(computeOutwardFlux) {
        result.outwardFlux.energy += magnetic.power;
        result.outwardFlux.momentum += magnetic.momentumRate;
        result.outwardFlux.angularMomentum += magnetic.angularMomentumRate;
    }
    return result;
}

bool finiteRadiationResponse(const ParticleMultipoleRadiation& response) {
    return std::isfinite(response.outwardFlux.energy)
        && std::isfinite(response.leadingElectricDipolePower)
        && std::isfinite(response.magneticDipolePower)
        && std::isfinite(response.electricQuadrupolePower)
        && std::isfinite(response.landauLifshitzValidity)
        && std::isfinite(response.coherentDerivativeConsistency)
        && std::isfinite(response.sourceCompactness)
        && std::isfinite(response.coherentWeight)
        && isFinite(response.chargeReaction.electron)
        && isFinite(response.chargeReaction.positron)
        && isFinite(response.electronDipoleTorque)
        && isFinite(response.positronDipoleTorque)
        && isFinite(response.outwardFlux.momentum)
        && isFinite(response.outwardFlux.angularMomentum);
}

void applyDipoleRadiationTorque(State& state,
                                const ParticleMultipoleRadiation& reaction,
                                double dt) {
    constexpr double electronGyromagneticRatio=-eCharge/(2.0*electronMass);
    constexpr double positronGyromagneticRatio=eCharge/(2.0*positronMass);
    synchronizeCovariantDipoles(state);
    const double electronNorm=state.electronProperDipole.norm();
    const double positronNorm=state.positronProperDipole.norm();
    state.electronProperDipole+=reaction.electronDipoleTorque
                         *(electronGyromagneticRatio*dt);
    state.positronProperDipole+=reaction.positronDipoleTorque
                         *(positronGyromagneticRatio*dt);
    if(state.electronProperDipole.norm()>0.0)
        state.electronProperDipole=state.electronProperDipole
            *(electronNorm/state.electronProperDipole.norm());
    if(state.positronProperDipole.norm()>0.0)
        state.positronProperDipole=state.positronProperDipole
            *(positronNorm/state.positronProperDipole.norm());
    synchronizeCovariantDipoles(state);
}

// f(r)=w(r)/r^3 and its derivatives.  Defining the regulator at the
// vector-potential level ensures that every use of the dipole field employs
// the same B=curl(A), including all derivatives of w.
MagneticRadialProfile magneticRadialProfile(double distance,
    double regularizationRadius=magneticRegularizationRadius,
    double exponent=magneticRegularizationExponent) {
    if(!(distance>0.0)) return {};
    const double weight = shortRangeFieldWeight(
        distance,regularizationRadius,exponent);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;
    const double inverseDistanceCubed = inverseDistanceSquared * inverseDistance;
    const double inverseDistanceFourth = inverseDistanceCubed * inverseDistance;
    const double inverseDistanceFifth = inverseDistanceFourth * inverseDistance;
    const double logarithmicWeightDerivative=exponent*weight*(1.0-weight);
    const double firstNumerator=logarithmicWeightDerivative-3.0*weight;
    const double logarithmicFirstNumeratorDerivative=
        exponent*weight*(1.0-weight)
        *(exponent*(1.0-2.0*weight)-3.0);
    return {weight*inverseDistanceCubed,
        firstNumerator*inverseDistanceFourth,
        (logarithmicFirstNumeratorDerivative-4.0*firstNumerator)
            *inverseDistanceFifth};
}

Vec3 regularizedDipoleField(const Vec3& sourceToTarget,
                            const Vec3& sourceDipole,
                            double regularizationRadius=magneticRegularizationRadius,
                            double exponent=magneticRegularizationExponent);

ElectromagneticField retardedElectricDipoleField(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsElectron) {
    ChargeKinematics source=historicalCharge(
        history,present,sourceIsElectron,observationTime);
    double retardedTime=observationTime
        -(observationPosition-source.position).norm()/c;
    for(int iteration=0;iteration<16;++iteration) {
        source=historicalCharge(history,present,sourceIsElectron,retardedTime);
        const Vec3 displacement=observationPosition-source.position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) return {};
        const Vec3 n=displacement/distance;
        const double residual=retardedTime+distance/c-observationTime;
        const double derivative=std::max(1.0e-8,
            1.0-dot(n,source.velocity/c));
        const double refined=retardedTime-residual/derivative;
        if(std::abs(refined-retardedTime)<=1.0e-30
            +1.0e-14*std::abs(retardedTime)) {
            retardedTime=refined; break;
        }
        retardedTime=refined;
    }
    const RetardedElectricDipoleKinematics dipole=
        historicalElectricDipoleKinematics(
            history,present,sourceIsElectron,retardedTime);
    const Vec3 displacement=observationPosition-dipole.position;
    const double distance=displacement.norm();
    if(!(distance>std::numeric_limits<double>::min())) return {};
    const Vec3 n=displacement/distance;
    const double inverseDistance=1.0/distance;
    const double weight=shortRangeFieldWeight(distance);
    const auto longitudinalPattern=[&](const Vec3& moment) {
        return n*(3.0*dot(n,moment))-moment;
    };
    constexpr double electricCoefficient=1.0/(4.0*pi*epsilon0);
    const Vec3 electric=(longitudinalPattern(dipole.moment)
            *(inverseDistance*inverseDistance*inverseDistance)
        +longitudinalPattern(dipole.firstDerivative)
            *(inverseDistance*inverseDistance/c)
        +cross(n,cross(n,dipole.secondDerivative))
            *(inverseDistance/(c*c)))*(electricCoefficient*weight);
    const Vec3 magnetic=(cross(dipole.firstDerivative,n)
            *(inverseDistance*inverseDistance)
        +cross(dipole.secondDerivative,n)*(inverseDistance/c))
        *(mu0/(4.0*pi)*weight);
    return {electric,magnetic};
}

// Retarded field generated by a time-dependent magnetic dipole. It follows
// directly from A=mu0/(4pi)[m(t_r)x n/r^2 + mdot(t_r)x n/(c r)].
// Consequently it contains the 1/r^3 near field, 1/r^2 induction field and
// 1/r radiation field. The same smooth short-range form factor as the static
// dipole sector is applied to prevent a point-dipole singularity.
ElectromagneticField retardedMagneticDipoleField(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsElectron) {
    ChargeKinematics source=historicalCharge(
        history,present,sourceIsElectron,observationTime);
    double retardedTime=observationTime
        -(observationPosition-source.position).norm()/c;
    for(int iteration=0;iteration<16;++iteration) {
        source=historicalCharge(history,present,sourceIsElectron,retardedTime);
        const Vec3 displacement=observationPosition-source.position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) return {};
        const Vec3 direction=displacement/distance;
        const double residual=retardedTime+distance/c-observationTime;
        const double derivative=std::max(1.0e-8,
            1.0-dot(direction,source.velocity/c));
        const double refined=retardedTime-residual/derivative;
        if(std::abs(refined-retardedTime)<=1.0e-30
            +1.0e-14*std::abs(retardedTime)) {
            retardedTime=refined; break;
        }
        retardedTime=refined;
    }
    const RetardedDipoleKinematics dipole=historicalDipoleKinematics(
        history,present,sourceIsElectron,retardedTime);
    const Vec3 displacement=observationPosition-dipole.position;
    const double distance=displacement.norm();
    if(!(distance>std::numeric_limits<double>::min())) return {};
    const Vec3 n=displacement/distance;
    const double inverseDistance=1.0/distance;
    const double weight=shortRangeFieldWeight(distance);
    constexpr double coefficient=mu0/(4.0*pi);
    const auto transversePattern=[&](const Vec3& moment) {
        return n*(3.0*dot(n,moment))-moment;
    };
    // The static term is exactly curl(A_reg), as in the conservative sector.
    // Retardation adds the induction and radiation pieces of the same
    // potential.  Thus mdot=mddot=0 reproduces regularizedDipoleField bit for
    // bit instead of the former, inconsistent w*B_point approximation.
    Vec3 magnetic=regularizedDipoleField(displacement,dipole.moment)
                 +(transversePattern(dipole.firstDerivative)
                      *(inverseDistance*inverseDistance/c)
                  +cross(n,cross(n,dipole.secondDerivative))
                      *(inverseDistance/(c*c)))*(coefficient*weight);
    Vec3 electric=(cross(n,dipole.firstDerivative)
                      *(inverseDistance*inverseDistance)
                  +cross(n,dipole.secondDerivative)
                      *(inverseDistance/c))*(coefficient*weight);
    return {electric,magnetic};
}

Vec3 regularizedDipoleVectorPotential(const Vec3& sourceToTarget,
                                      const Vec3& sourceDipole) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const MagneticRadialProfile profile = magneticRadialProfile(sourceToTarget.norm());
    return cross(sourceDipole, sourceToTarget)
         * (magneticConstant * profile.vectorPotentialFactor);
}

// B=curl(A) for A=mu0/(4 pi) f(r) mu x r.  In particular, this is not
// merely w times the unregularized point-dipole field.
Vec3 regularizedDipoleField(const Vec3& sourceToTarget,const Vec3& sourceDipole,
                            double regularizationRadius,double exponent) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const double distance = sourceToTarget.norm();
    const Vec3 n = sourceToTarget / distance;
    const MagneticRadialProfile profile = magneticRadialProfile(
        distance,regularizationRadius,exponent);
    const double transverseCoefficient = 2.0 * profile.vectorPotentialFactor
                                       + distance * profile.firstDerivative;
    const double radialCoefficient = -distance * profile.firstDerivative;
    return (sourceDipole * transverseCoefficient
          + n * (dot(sourceDipole, n) * radialCoefficient)) * magneticConstant;
}

// F=-grad(U), U=-mu_target.B_source, using the same regularized field above.
Vec3 regularizedDipoleForce(const Vec3& sourceToTarget,
                            const Vec3& targetDipole,const Vec3& sourceDipole,
                            double regularizationRadius=magneticRegularizationRadius,
                            double exponent=magneticRegularizationExponent) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const double distance = sourceToTarget.norm();
    const Vec3 n = sourceToTarget / distance;
    const MagneticRadialProfile profile = magneticRadialProfile(
        distance,regularizationRadius,exponent);
    const double transverseCoefficientDerivative =
        3.0 * profile.firstDerivative + distance * profile.secondDerivative;
    const double radialCoefficient = -distance * profile.firstDerivative;
    const double radialCoefficientDerivative =
        -profile.firstDerivative - distance * profile.secondDerivative;
    const double targetRadial = dot(targetDipole, n);
    const double sourceRadial = dot(sourceDipole, n);
    const double dipoleDot = dot(targetDipole, sourceDipole);
    return (n * (transverseCoefficientDerivative * dipoleDot
               + radialCoefficientDerivative * targetRadial * sourceRadial)
          + (targetDipole * sourceRadial + sourceDipole * targetRadial
             - n * (2.0 * targetRadial * sourceRadial))
            * (radialCoefficient / distance)) * magneticConstant;
}

MutualForces coulombForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 electron = geometry.electronMinusPositron
                        * (-coulomb * eCharge * eCharge * geometry.inverseDistanceCubed);
    return {electron, electron * -1.0};
}

MutualForces mutualForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const MutualForces electrostatic = coulombForces(s);
    const Vec3 dipoleOnElectron = regularizedDipoleForce(
        geometry.electronMinusPositron, s.electronDipole, s.positronDipole);
    return {electrostatic.electron + dipoleOnElectron,
            electrostatic.positron - dipoleOnElectron};
}

Vec3 darwinForceOnFirst(const Vec3& firstVelocity, const Vec3& secondVelocity,
                        const Vec3& secondLeadingAcceleration,
                        const Vec3& firstMinusSecond, double chargeProduct) {
    const double distance = firstMinusSecond.norm();
    const Vec3 n = firstMinusSecond / distance;
    const double firstRadial = dot(firstVelocity, n);
    const double secondRadial = dot(secondVelocity, n);
    const double radialRate = firstRadial - secondRadial;
    const Vec3 relativeVelocity = firstVelocity - secondVelocity;
    const Vec3 nRate = (relativeVelocity - n * radialRate) / distance;
    const double secondRadialRate = dot(secondLeadingAcceleration, n)
                                  + dot(secondVelocity, nRate);
    const double coefficient = coulomb * chargeProduct / (2.0 * c*c * distance);
    const double velocityProduct = dot(firstVelocity, secondVelocity);

    // Euler-Lagrange force from
    // L_D = C [v1.v2 + (v1.n)(v2.n)]. Accelerations inside d/dt(dL/dv)
    // are reduced to their leading Coulomb values, consistently through
    // order v^2/c^2.
    const Vec3 spatialDerivative =
        (secondVelocity * firstRadial + firstVelocity * secondRadial
       - n * (velocityProduct + 3.0 * firstRadial * secondRadial))
        * (coefficient / distance);
    const Vec3 momentumDerivative =
        (secondVelocity + n * secondRadial)
            * (-coefficient * radialRate / distance)
      + (secondLeadingAcceleration + nRate * secondRadial
       + n * secondRadialRate) * coefficient;
    return spatialDerivative - momentumDerivative;
}

MutualForces darwinForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const MutualForces leading = coulombForces(s);
    const Vec3 electronLeadingAcceleration = leading.electron / electronMass;
    const Vec3 positronLeadingAcceleration = leading.positron / positronMass;
    constexpr double chargeProduct = -eCharge * eCharge;
    return {
        darwinForceOnFirst(s.electronVelocity, s.positronVelocity,
                           positronLeadingAcceleration,
                           geometry.electronMinusPositron, chargeProduct),
        darwinForceOnFirst(s.positronVelocity, s.electronVelocity,
                           electronLeadingAcceleration,
                           geometry.electronMinusPositron * -1.0, chargeProduct)
    };
}

double darwinInteractionEnergy(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 n = geometry.electronMinusPositron * geometry.inverseDistance;
    constexpr double chargeProduct = -eCharge * eCharge;
    return coulomb * chargeProduct * geometry.inverseDistance / (2.0 * c*c)
         * (dot(s.electronVelocity, s.positronVelocity)
          + dot(s.electronVelocity, n) * dot(s.positronVelocity, n));
}

LocalElectromagneticFields localRelativisticFields(
    const State& s, const StateHistory& history) {
    ElectromagneticField atElectron = lienardWiechertField(
        s.electronPosition, s.time, history, s, false, eCharge);
    ElectromagneticField atPositron = lienardWiechertField(
        s.positronPosition, s.time, history, s, true, -eCharge);
    const ElectromagneticField positronDipole=retardedMagneticDipoleField(
        s.electronPosition,s.time,history,s,false);
    const ElectromagneticField electronDipole=retardedMagneticDipoleField(
        s.positronPosition,s.time,history,s,true);
    const ElectromagneticField positronElectricDipole=
        retardedElectricDipoleField(
            s.electronPosition,s.time,history,s,false);
    const ElectromagneticField electronElectricDipole=
        retardedElectricDipoleField(
            s.positronPosition,s.time,history,s,true);
    atElectron.electric+=positronDipole.electric
        +positronElectricDipole.electric;
    atElectron.magnetic+=positronDipole.magnetic
        +positronElectricDipole.magnetic;
    atPositron.electric+=electronDipole.electric
        +electronElectricDipole.electric;
    atPositron.magnetic+=electronDipole.magnetic
        +electronElectricDipole.magnetic;
    return {atElectron, atPositron};
}

DipoleDerivatives thomasBmtDipoleDerivatives(
    const State& s, const StateHistory& history) {
    const LocalElectromagneticFields fields = localRelativisticFields(s, history);
    const Vec3 electronEffectiveField = thomasBmtEffectiveField(
        s.electronVelocity, fields.atElectron);
    const Vec3 positronEffectiveField = thomasBmtEffectiveField(
        s.positronVelocity, fields.atPositron);
    return {
        cross(s.electronDipole, electronEffectiveField)
            * (-eCharge / electronMass),
        cross(s.positronDipole, positronEffectiveField)
            * (eCharge / positronMass)
    };
}

FourVector dipoleFourVector(const Vec3& properDipole,const Vec3& velocity) {
    const double relativisticGamma=gamma(velocity);
    const double projection=dot(velocity,properDipole);
    return {relativisticGamma*projection/c,
        properDipole+velocity*(relativisticGamma*relativisticGamma
            /(relativisticGamma+1.0)*projection/(c*c))};
}

Vec3 properDipoleFromFourVector(const FourVector& dipole,
                                const Vec3& velocity) {
    const double relativisticGamma=gamma(velocity);
    return dipole.space-velocity*(relativisticGamma
        /(relativisticGamma+1.0)*dipole.time/c);
}

FourVector electromagneticTensorAction(const ElectromagneticField& field,
                                        const FourVector& vector) {
    return {dot(field.electric,vector.space)/c,
        field.electric*(vector.time/c)+cross(vector.space,field.magnetic)};
}

Vec3 advanceCovariantBmt(const Vec3& properDipole,const Vec3& velocity,
                         const ElectromagneticField& field,
                         double chargeToMass,double laboratoryDt) {
    const double targetNorm=properDipole.norm();
    if(targetNorm==0.0||laboratoryDt==0.0) return properDipole;
    constexpr double classicalG=1.0;
    constexpr double halfG=0.5*classicalG;
    const FourVector fourVelocityValue{
        gamma(velocity)*c,velocity*gamma(velocity)};
    const FourVector fieldOnVelocity=electromagneticTensorAction(
        field,fourVelocityValue);
    const auto derivative=[&](const FourVector& dipole) {
        const FourVector fieldOnDipole=electromagneticTensorAction(field,dipole);
        const double contraction=minkowskiDot(dipole,fieldOnVelocity)/(c*c);
        return FourVector{
            chargeToMass*(halfG*fieldOnDipole.time
                +(halfG-1.0)*fourVelocityValue.time*contraction),
            (fieldOnDipole.space*halfG
                +fourVelocityValue.space*((halfG-1.0)*contraction))
                *chargeToMass};
    };
    const double properDt=laboratoryDt/gamma(velocity);
    const auto add=[](const FourVector& value,const FourVector& increment,
                      double scale) {
        return FourVector{value.time+increment.time*scale,
                          value.space+increment.space*scale};
    };
    FourVector spin=dipoleFourVector(properDipole,velocity);
    const FourVector k1=derivative(spin);
    const FourVector k2=derivative(add(spin,k1,0.5*properDt));
    const FourVector k3=derivative(add(spin,k2,0.5*properDt));
    const FourVector k4=derivative(add(spin,k3,properDt));
    spin.time+=properDt*(k1.time+2.0*k2.time+2.0*k3.time+k4.time)/6.0;
    spin.space+=(k1.space+k2.space*2.0+k3.space*2.0+k4.space)
        *(properDt/6.0);
    // Remove round-off along u and restore the invariant spacelike norm.
    const double parallel=minkowskiDot(spin,fourVelocityValue)/(c*c);
    spin.time-=fourVelocityValue.time*parallel;
    spin.space=spin.space-fourVelocityValue.space*parallel;
    Vec3 result=properDipoleFromFourVector(spin,velocity);
    if(result.norm()>0.0) result=result*(targetNorm/result.norm());
    return result;
}

void applyDipolePrecession(State& s, double dt,
                           const StateHistory& history) {
    synchronizeCovariantDipoles(s);
    const LocalElectromagneticFields fields = localRelativisticFields(s, history);
    const Vec3 electronDipole=advanceCovariantBmt(s.electronProperDipole,
        s.electronVelocity,fields.atElectron,-eCharge/electronMass,dt);
    const Vec3 positronDipole=advanceCovariantBmt(s.positronProperDipole,
        s.positronVelocity,fields.atPositron,eCharge/positronMass,dt);
    // Update simultaneously so neither particle sees an already-updated peer.
    s.electronProperDipole=electronDipole;
    s.positronProperDipole=positronDipole;
    synchronizeCovariantDipoles(s);
}

struct ChargeDipolePairForces { Vec3 onCharge, onDipole; };

// Euler-Lagrange forces generated by
// L_qmu = q (v_charge-v_dipole).A_mu,
// A_mu = mu0/(4 pi) w(r) mu x R/r^3,
// where R points from the dipole to the charge.  Returning both variations of
// the same translationally and boost-invariant low-velocity interaction gives
// equal and opposite mechanical forces.
ChargeDipolePairForces chargeDipolePairForces(
    const Vec3& relativeVelocity, double charge, const Vec3& sourceDipole,
    const Vec3& sourceDipoleDerivative, const Vec3& sourceToCharge) {
    const double distance = sourceToCharge.norm();
    const MagneticRadialProfile profile = magneticRadialProfile(distance);
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const Vec3 magneticField = regularizedDipoleField(sourceToCharge, sourceDipole);
    const Vec3 inducedElectricTerm = cross(sourceDipoleDerivative, sourceToCharge)
        * (magneticConstant * profile.vectorPotentialFactor);
    const Vec3 onCharge =
        (cross(relativeVelocity, magneticField) - inducedElectricTerm) * charge;
    return {onCharge, onCharge * -1.0};
}

MutualForces chargeDipoleForces(const State& s, const StateHistory& history) {
    const DipoleDerivatives derivatives = thomasBmtDipoleDerivatives(s, history);
    const Vec3 electronMinusPositron = s.electronPosition - s.positronPosition;
    const Vec3 electronMinusPositronVelocity =
        s.electronVelocity - s.positronVelocity;

    const ChargeDipolePairForces electronChargePositronDipole = chargeDipolePairForces(
        electronMinusPositronVelocity, -eCharge, s.positronDipole,
        derivatives.positron, electronMinusPositron);
    const ChargeDipolePairForces positronChargeElectronDipole = chargeDipolePairForces(
        electronMinusPositronVelocity * -1.0, eCharge, s.electronDipole,
        derivatives.electron, electronMinusPositron * -1.0);

    return {electronChargePositronDipole.onCharge + positronChargeElectronDipole.onDipole,
            positronChargeElectronDipole.onCharge + electronChargePositronDipole.onDipole};
}

MutualForces allExternalForces(const State& s) {
    const MutualForces positionForces = mutualForces(s);
    const MutualForces velocityForces = darwinForces(s);
    const StateHistory localHistory{State{s}};
    const MutualForces mixedMagneticForces = chargeDipoleForces(s, localHistory);
    return {positionForces.electron + velocityForces.electron + mixedMagneticForces.electron,
            positionForces.positron + velocityForces.positron + mixedMagneticForces.positron};
}

ElectromagneticField fieldFromOtherParticleAt(
    const Vec3& observationPosition,const State& state,
    const StateHistory& history,bool targetIsElectron) {
    const bool sourceIsElectron=!targetIsElectron;
    const double sourceCharge=sourceIsElectron?-eCharge:eCharge;
    ElectromagneticField field=lienardWiechertField(
        observationPosition,state.time,history,state,sourceIsElectron,
        sourceCharge);
    const ElectromagneticField magneticDipole=retardedMagneticDipoleField(
        observationPosition,state.time,history,state,sourceIsElectron);
    const ElectromagneticField electricDipole=retardedElectricDipoleField(
        observationPosition,state.time,history,state,sourceIsElectron);
    field.electric+=magneticDipole.electric+electricDipole.electric;
    field.magnetic+=magneticDipole.magnetic+electricDipole.magnetic;
    return field;
}

Vec3 magneticFieldInRestFrame(const ElectromagneticField& field,
                              const Vec3& velocity) {
    const double relativisticGamma=gamma(velocity);
    return (field.magnetic-cross(velocity,field.electric)/(c*c))
            *relativisticGamma
        -velocity*(relativisticGamma*relativisticGamma
            /(relativisticGamma+1.0)*dot(velocity,field.magnetic)/(c*c));
}

Vec3 covariantDipoleGradientForce(const State& state,
                                  const StateHistory& history,
                                  bool targetIsElectron) {
    const Vec3 targetPosition=targetIsElectron?state.electronPosition
                                               :state.positronPosition;
    const Vec3 targetVelocity=targetIsElectron?state.electronVelocity
                                               :state.positronVelocity;
    const Vec3 properDipole=targetIsElectron?state.electronProperDipole
                                             :state.positronProperDipole;
    if(properDipole.squaredNorm()==0.0) return {};
    const double gradientStep=std::max(
        1.0e-6*separation(state),1.0e-3*nuclearCutoff);
    const auto coupling=[&](const Vec3& point) {
        return dot(properDipole,magneticFieldInRestFrame(
            fieldFromOtherParticleAt(point,state,history,targetIsElectron),
            targetVelocity));
    };
    Vec3 gradient;
    for(int axis=0;axis<3;++axis) {
        Vec3 offset;
        if(axis==0) offset.x=gradientStep;
        else if(axis==1) offset.y=gradientStep;
        else offset.z=gradientStep;
        const double derivative=(coupling(targetPosition+offset)
                                -coupling(targetPosition-offset))
                               /(2.0*gradientStep);
        if(axis==0) gradient.x=derivative;
        else if(axis==1) gradient.y=derivative;
        else gradient.z=derivative;
    }
    // Spatial component of the covariant gradient divided by gamma gives
    // the laboratory three-force.  At rest this reduces to grad(mu.B).
    return gradient/gamma(targetVelocity);
}

MutualForces retardedExternalForces(const State& s,
                                    const StateHistory& history) {
    const ElectromagneticField positronField = lienardWiechertField(
        s.electronPosition, s.time, history, s, false, eCharge);
    const ElectromagneticField electronField = lienardWiechertField(
        s.positronPosition, s.time, history, s, true, -eCharge);
    const ElectromagneticField positronDipoleField=retardedMagneticDipoleField(
        s.electronPosition,s.time,history,s,false);
    const ElectromagneticField electronDipoleField=retardedMagneticDipoleField(
        s.positronPosition,s.time,history,s,true);
    const ElectromagneticField positronElectricDipoleField=
        retardedElectricDipoleField(
            s.electronPosition,s.time,history,s,false);
    const ElectromagneticField electronElectricDipoleField=
        retardedElectricDipoleField(
            s.positronPosition,s.time,history,s,true);
    const MutualForces chargeCharge{
        lorentzForce(-eCharge,s.electronVelocity,
            {positronField.electric+positronDipoleField.electric
                +positronElectricDipoleField.electric,
             positronField.magnetic+positronDipoleField.magnetic
                +positronElectricDipoleField.magnetic}),
        lorentzForce(eCharge,s.positronVelocity,
            {electronField.electric+electronDipoleField.electric
                +electronElectricDipoleField.electric,
             electronField.magnetic+electronDipoleField.magnetic
                +electronElectricDipoleField.magnetic})};

    const MutualForces tensorGradient{
        covariantDipoleGradientForce(s,history,true),
        covariantDipoleGradientForce(s,history,false)};
    return {chargeCharge.electron+tensorGradient.electron,
            chargeCharge.positron+tensorGradient.positron};
}

struct CanonicalMomenta { Vec3 electron, positron; };

// Canonical momenta obtained from the conservative low-velocity action.  The
// q-mu contributions occur with opposite signs on the charge and on the
// dipole carrier; the Darwin terms do not cancel particle by particle.
CanonicalMomenta canonicalMomenta(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 n = geometry.electronMinusPositron * geometry.inverseDistance;
    constexpr double chargeProduct = -eCharge * eCharge;
    const double darwinCoefficient = coulomb * chargeProduct
        * geometry.inverseDistance / (2.0 * c*c);

    Vec3 electron = momentum(s.electronVelocity, electronMass)
        + (s.positronVelocity + n * dot(s.positronVelocity, n)) * darwinCoefficient;
    Vec3 positron = momentum(s.positronVelocity, positronMass)
        + (s.electronVelocity + n * dot(s.electronVelocity, n)) * darwinCoefficient;

    const Vec3 positronPotentialAtElectron = regularizedDipoleVectorPotential(
        geometry.electronMinusPositron, s.positronDipole);
    const Vec3 electronPotentialAtPositron = regularizedDipoleVectorPotential(
        geometry.electronMinusPositron * -1.0, s.electronDipole);
    const Vec3 mixedCanonicalContribution =
        positronPotentialAtElectron * (-eCharge)
      - electronPotentialAtPositron * eCharge;
    electron += mixedCanonicalContribution;
    positron += mixedCanonicalContribution * -1.0;
    return {electron, positron};
}

// Overloads taking an already-computed CanonicalMomenta.  Callers that need
// more than one of these quantities for the same state used to recompute the
// canonical momenta once per call: three times per frame and four times per
// accepted integration step, each redoing the pair geometry, the relativistic
// momenta and two regularized dipole vector potentials.
Vec3 noetherMomentum(const CanonicalMomenta& canonical) {
    return canonical.electron + canonical.positron;
}

Vec3 noetherMomentum(const State& s) {
    return noetherMomentum(canonicalMomenta(s));
}

double canonicalMomentumScale(const CanonicalMomenta& canonical) {
    return canonical.electron.norm() + canonical.positron.norm();
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Production callers already hold the canonical momenta and use the overload
// above; only the Yee-coupling tests start from a bare state.
double canonicalMomentumScale(const State& s) {
    return canonicalMomentumScale(canonicalMomenta(s));
}
#endif

Vec3 noetherAngularMomentum(const State& s,
                            const CanonicalMomenta& canonical) {
    constexpr double electronGyromagneticRatio = -eCharge / (2.0 * electronMass);
    constexpr double positronGyromagneticRatio = eCharge / (2.0 * positronMass);
    const Vec3 orbital = cross(s.electronPosition, canonical.electron)
                       + cross(s.positronPosition, canonical.positron);
    const Vec3 intrinsic = s.electronDipole / electronGyromagneticRatio
                         + s.positronDipole / positronGyromagneticRatio;
    return orbital + intrinsic;
}

Vec3 noetherAngularMomentum(const State& s) {
    return noetherAngularMomentum(s, canonicalMomenta(s));
}

double conservativeParticleEnergy(const State& state) {
    const PairGeometry geometry=pairGeometry(state);
    const double kinetic=(gamma(state.electronVelocity)-1.0)*electronMass*c*c
        +(gamma(state.positronVelocity)-1.0)*positronMass*c*c;
    const double coulombPotential=-coulomb*eCharge*eCharge
        *geometry.inverseDistance;
    const double dipolePotential=-dot(state.electronDipole,
        regularizedDipoleField(
            geometry.electronMinusPositron,state.positronDipole));
    return kinetic+coulombPotential+dipolePotential
        +darwinInteractionEnergy(state)+state.dipoleConstraintEnergy;
}

FourVector fourVelocity(const Vec3& velocity) {
    const double relativisticGamma = gamma(velocity);
    return {relativisticGamma * c, velocity * relativisticGamma};
}

FourVector fourForce(const Vec3& velocity, const Vec3& force) {
    const double relativisticGamma = gamma(velocity);
    return {relativisticGamma * dot(force, velocity) / c,
            force * relativisticGamma};
}

Vec3 reducedOrderSelfForce(const Vec3& velocity,
                           const FourVector& beforeForce,
                           const FourVector& afterForce,
                           double laboratoryDerivativeStep,
                           double mass) {
    const double relativisticGamma = gamma(velocity);
    const double properDerivativeScale = relativisticGamma
                                       / (2.0 * laboratoryDerivativeStep);
    const FourVector forceDerivative{
        (afterForce.time - beforeForce.time) * properDerivativeScale,
        (afterForce.space - beforeForce.space) * properDerivativeScale};
    const FourVector velocityFour = fourVelocity(velocity);
    const double parallelCoefficient =
        minkowskiDot(velocityFour, forceDerivative) / (c*c);
    const FourVector orthogonalDerivative{
        forceDerivative.time - velocityFour.time * parallelCoefficient,
        forceDerivative.space - velocityFour.space * parallelCoefficient};
    const double characteristicTime = eCharge * eCharge
        / (6.0 * pi * epsilon0 * mass * c*c*c);
    // The spatial component of a four-force is gamma times the ordinary
    // laboratory three-force.
    return orthogonalDerivative.space * (characteristicTime / relativisticGamma);
}

MutualForces individualLandauLifshitzSelfForces(
    const State& s, const MutualForces& external,
    const StateHistory& history) {
    const Vec3 electronAcceleration = relativisticAcceleration(s.electronVelocity, external.electron,
                                                                electronMass);
    const Vec3 positronAcceleration = relativisticAcceleration(s.positronVelocity, external.positron,
                                                                positronMass);
    const double relativeSpeed = (s.electronVelocity - s.positronVelocity).norm();
    const double derivativeStep = std::max(1.0e-24,
        1.0e-4 * std::min(separation(s) / std::max(relativeSpeed, 1.0),
                          separation(s) / c));
    State before = s;
    State after = s;
    before.electronPosition = s.electronPosition - s.electronVelocity * derivativeStep
        + electronAcceleration * (0.5 * derivativeStep * derivativeStep);
    before.positronPosition = s.positronPosition - s.positronVelocity * derivativeStep
        + positronAcceleration * (0.5 * derivativeStep * derivativeStep);
    after.electronPosition = s.electronPosition + s.electronVelocity * derivativeStep
        + electronAcceleration * (0.5 * derivativeStep * derivativeStep);
    after.positronPosition = s.positronPosition + s.positronVelocity * derivativeStep
        + positronAcceleration * (0.5 * derivativeStep * derivativeStep);
    before.electronVelocity = velocityFromMomentum(
        momentum(s.electronVelocity, electronMass)
            - external.electron * derivativeStep, electronMass);
    before.positronVelocity = velocityFromMomentum(
        momentum(s.positronVelocity, positronMass)
            - external.positron * derivativeStep, positronMass);
    after.electronVelocity = velocityFromMomentum(
        momentum(s.electronVelocity, electronMass)
            + external.electron * derivativeStep, electronMass);
    after.positronVelocity = velocityFromMomentum(
        momentum(s.positronVelocity, positronMass)
            + external.positron * derivativeStep, positronMass);
    before.time -= derivativeStep;
    after.time += derivativeStep;
    const MutualForces beforeForces = retardedExternalForces(before, history);
    const MutualForces afterForces = retardedExternalForces(after, history);
    return {
        reducedOrderSelfForce(s.electronVelocity,
            fourForce(before.electronVelocity, beforeForces.electron),
            fourForce(after.electronVelocity, afterForces.electron),
            derivativeStep, electronMass),
        reducedOrderSelfForce(s.positronVelocity,
            fourForce(before.positronVelocity, beforeForces.positron),
            fourForce(after.positronVelocity, afterForces.positron),
            derivativeStep, positronMass)
    };
}

// Relativistic predictor-corrector update with mutually retarded fields and
// an individual, reduced-order Landau-Lifshitz self-force for each charge.
void integrateElectrodynamicStep(State& s, double dt,
                                 const StateHistory& history,
                                 bool computeOutwardFlux=true,
                                 ChargeRadiationReactionModel reactionModel=
                                    ChargeRadiationReactionModel::individualLandauLifshitz) {
    const State balanceStart=s;
    const double initialMechanicalEnergy=conservativeParticleEnergy(balanceStart);
    const CanonicalMomenta initialCanonical=canonicalMomenta(balanceStart);
    const Vec3 initialMechanicalMomentum=noetherMomentum(initialCanonical);
    const Vec3 initialMechanicalAngularMomentum=
        noetherAngularMomentum(balanceStart,initialCanonical);
    applyDipolePrecession(s, 0.5 * dt, history);
    MutualForces forces = retardedExternalForces(s, history);
    const ParticleMultipoleRadiation radiation =
        particleMultipoleRadiation(
            s,forces,history,computeOutwardFlux,reactionModel);
    if (!finiteRadiationResponse(radiation)) {
        s.time = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    applyDipoleRadiationTorque(s,radiation,0.5*dt);
    Vec3 electronMomentum = momentum(s.electronVelocity, electronMass)
        + (forces.electron + radiation.chargeReaction.electron) * (0.5 * dt);
    Vec3 positronMomentum = momentum(s.positronVelocity, positronMass)
        + (forces.positron + radiation.chargeReaction.positron) * (0.5 * dt);
    State trial = s;
    trial.time += dt;
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    electronMomentum = momentum(trial.electronVelocity, electronMass);
    positronMomentum = momentum(trial.positronVelocity, positronMass);
    trial.electronPosition += trial.electronVelocity * dt;
    trial.positronPosition += trial.positronVelocity * dt;
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, forces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, forces.positron, positronMass);

    const MutualForces trialForces = retardedExternalForces(trial, history);
    const ParticleMultipoleRadiation trialRadiation =
        particleMultipoleRadiation(
            trial,trialForces,history,computeOutwardFlux,reactionModel);
    if (!finiteRadiationResponse(trialRadiation)) {
        s.time = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    electronMomentum += (trialForces.electron
        + trialRadiation.chargeReaction.electron) * (0.5 * dt);
    positronMomentum += (trialForces.positron
        + trialRadiation.chargeReaction.positron) * (0.5 * dt);
    trial.electronVelocity = velocityFromMomentum(electronMomentum, electronMass);
    trial.positronVelocity = velocityFromMomentum(positronMomentum, positronMass);
    applyDipolePrecession(trial, 0.5 * dt, history);
    applyDipoleRadiationTorque(trial,trialRadiation,0.5*dt);
    trial.electronAcceleration = relativisticAcceleration(trial.electronVelocity, trialForces.electron, electronMass);
    trial.positronAcceleration = relativisticAcceleration(trial.positronVelocity, trialForces.positron, positronMass);
    const double radiatedEnergyIncrement = 0.5
        * (radiation.outwardFlux.energy + trialRadiation.outwardFlux.energy)*dt;
    trial.radiatedEnergy += radiatedEnergyIncrement;
    // Magnetic-dipole damping changes the constrained internal dipole sector.
    // The charge part is already represented by the particle self-force and
    // its near-field (Schott) term.
    const double dipoleRadiatedEnergy=0.5*(radiation.magneticDipolePower
        +trialRadiation.magneticDipolePower)*dt;
    trial.dipoleConstraintEnergy-=dipoleRadiatedEnergy;
    trial.radiatedMomentum += (radiation.outwardFlux.momentum
        + trialRadiation.outwardFlux.momentum) * (0.5*dt);
    trial.radiatedAngularMomentum +=
        (radiation.outwardFlux.angularMomentum
        + trialRadiation.outwardFlux.angularMomentum) * (0.5*dt);

    if(computeOutwardFlux) {
        // Exact discrete world-tube balance.  This reservoir contains bound
        // near-field, retardation/interference energy and the discretization
        // remainder; it is kept separate from the independently measured
        // mismatch of the LL force and coherent far radiation below.
        const double mechanicalEnergyChange=
            conservativeParticleEnergy(trial)-initialMechanicalEnergy;
        const CanonicalMomenta trialCanonical=canonicalMomenta(trial);
        const Vec3 mechanicalMomentumChange=
            noetherMomentum(trialCanonical)-initialMechanicalMomentum;
        const Vec3 mechanicalAngularChange=
            noetherAngularMomentum(trial,trialCanonical)
            -initialMechanicalAngularMomentum;
        trial.boundFieldEnergy=balanceStart.boundFieldEnergy
            -mechanicalEnergyChange-radiatedEnergyIncrement;
        trial.boundFieldMomentum=balanceStart.boundFieldMomentum
            -mechanicalMomentumChange
            -(trial.radiatedMomentum-balanceStart.radiatedMomentum);
        trial.boundFieldAngularMomentum=balanceStart.boundFieldAngularMomentum
            -mechanicalAngularChange
            -(trial.radiatedAngularMomentum
              -balanceStart.radiatedAngularMomentum);

        const auto chargeMismatchRates=[](
            const State& state,const ParticleMultipoleRadiation& response) {
            FieldFluxRates mismatch=response.outwardFlux;
            mismatch.energy-=response.magneticDipolePower;
            mismatch.energy+=dot(response.chargeReaction.electron,
                                  state.electronVelocity)
                            +dot(response.chargeReaction.positron,
                                  state.positronVelocity);
            const Vec3 reactionForce=response.chargeReaction.electron
                                    +response.chargeReaction.positron;
            mismatch.momentum+=reactionForce;
            mismatch.angularMomentum+=cross(state.electronPosition,
                                             response.chargeReaction.electron)
                +cross(state.positronPosition,response.chargeReaction.positron);
            return mismatch;
        };
        const FieldFluxRates initialMismatch=
            chargeMismatchRates(balanceStart,radiation);
        const FieldFluxRates finalMismatch=
            chargeMismatchRates(trial,trialRadiation);
        trial.reactionEnergyMismatch=balanceStart.reactionEnergyMismatch
            +0.5*(initialMismatch.energy+finalMismatch.energy)*dt;
        trial.reactionMomentumMismatch=balanceStart.reactionMomentumMismatch
            +(initialMismatch.momentum+finalMismatch.momentum)*(0.5*dt);
        trial.reactionAngularMomentumMismatch=
            balanceStart.reactionAngularMomentumMismatch
            +(initialMismatch.angularMomentum+finalMismatch.angularMomentum)
                *(0.5*dt);
    }
    s = trial;
}

void appendStateHistory(StateHistory& history, const State& state) {
    history.push_back(state);
    const double retentionTime = std::max(1.0e-20, 4.0 * separation(state) / c);
    const double earliestNeeded = state.time - retentionTime;
    while (history.size() > 2 && history[1].time < earliestNeeded) {
        history.pop_front();
    }
    // Cap the sample count inside the retarded window.  The window length is
    // set by the light-crossing time of the pair, but the number of samples
    // needed inside it is set by how fast the source moves, not by the
    // integration step.  Without a cap the count is window/step, which grows
    // like r^{-1/2} and reached ten thousand entries at short range.  Ordinary
    // orbits hold three to six nodes, so this never triggers there and their
    // results are bit-identical.
    constexpr std::size_t maximumHistoryNodes = 128;
    if (history.size() > maximumHistoryNodes) {
        StateHistory thinned;
        const std::size_t keepRecent = maximumHistoryNodes/2;
        const std::size_t oldCount = history.size() - keepRecent;
        for (std::size_t index = 0; index < oldCount; index += 2) {
            thinned.push_back(history[index]);
        }
        for (std::size_t index = oldCount; index < history.size(); ++index) {
            thinned.push_back(history[index]);
        }
        history = std::move(thinned);
    }
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Self-adjoint reference map for the instantaneous conservative model.  It is
// deliberately separate from the causal retarded/radiating engine: applying
// the same converged midpoint solve with -dt must recover the initial state.
// Only the time-reversibility test uses it, so it stays out of production.
void integrateConservativeMidpoint(
    State& state,double dt,int iterations=8) {
    if(dt==0.0||!std::isfinite(dt)) return;
    const State start=state;
    const Vec3 electronMomentum0=momentum(start.electronVelocity,electronMass);
    const Vec3 positronMomentum0=momentum(start.positronVelocity,positronMass);
    State endpoint=start;
    for(int iteration=0;iteration<std::max(iterations,1);++iteration) {
        State midpoint=start;
        midpoint.time=start.time+0.5*dt;
        midpoint.electronPosition=(start.electronPosition
                                   +endpoint.electronPosition)*0.5;
        midpoint.positronPosition=(start.positronPosition
                                   +endpoint.positronPosition)*0.5;
        const Vec3 endpointElectronMomentum=momentum(
            endpoint.electronVelocity,electronMass);
        const Vec3 endpointPositronMomentum=momentum(
            endpoint.positronVelocity,positronMass);
        midpoint.electronVelocity=velocityFromMomentum(
            (electronMomentum0+endpointElectronMomentum)*0.5,electronMass);
        midpoint.positronVelocity=velocityFromMomentum(
            (positronMomentum0+endpointPositronMomentum)*0.5,positronMass);
        midpoint.electronDipole=(start.electronDipole
                                 +endpoint.electronDipole)*0.5;
        midpoint.positronDipole=(start.positronDipole
                                 +endpoint.positronDipole)*0.5;
        const MutualForces force=allExternalForces(midpoint);
        endpoint.electronVelocity=velocityFromMomentum(
            electronMomentum0+force.electron*dt,electronMass);
        endpoint.positronVelocity=velocityFromMomentum(
            positronMomentum0+force.positron*dt,positronMass);
        const Vec3 midpointElectronMomentum=(electronMomentum0
            +momentum(endpoint.electronVelocity,electronMass))*0.5;
        const Vec3 midpointPositronMomentum=(positronMomentum0
            +momentum(endpoint.positronVelocity,positronMass))*0.5;
        endpoint.electronPosition=start.electronPosition
            +velocityFromMomentum(midpointElectronMomentum,electronMass)*dt;
        endpoint.positronPosition=start.positronPosition
            +velocityFromMomentum(midpointPositronMomentum,positronMass)*dt;
    }
    endpoint.time=start.time+dt;
    endpoint.electronAcceleration=(endpoint.electronVelocity
                                  -start.electronVelocity)/dt;
    endpoint.positronAcceleration=(endpoint.positronVelocity
                                  -start.positronVelocity)/dt;
    state=endpoint;
}
#endif

// The single classical trajectory engine shared by Visual and Statistical.
// Sampling policies and stopping conditions belong to their callers; the
// equations of motion, retarded history and radiation bookkeeping do not.
