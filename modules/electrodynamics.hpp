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
    const Vec3 firstMomentum=momentum(state.firstVelocity,firstMass);
    const Vec3 secondMomentum=momentum(state.secondVelocity,secondMass);
    const MaxwellVolumeIntegrals fieldTotals=field.volumeIntegrals();
    const double firstGyromagneticRatio=firstCharge/(2.0*firstMass);
    const double secondGyromagneticRatio=secondCharge/(2.0*secondMass);
    // Rest energies are constant and omitted from the diagnostic difference.
    return {(gamma(state.firstVelocity)-1.0)*firstMass*c*c
          +(gamma(state.secondVelocity)-1.0)*secondMass*c*c
          +fieldTotals.energy+state.dipoleConstraintEnergy,
            firstMomentum+secondMomentum+fieldTotals.momentum,
            cross(state.firstPosition,firstMomentum)
          +cross(state.secondPosition,secondMomentum)
          +state.firstDipole/firstGyromagneticRatio
          +state.secondDipole/secondGyromagneticRatio
          +fieldTotals.angularMomentum};
}
#endif

// Effective field in the Thomas-BMT equation, expressed in laboratory time:
// d(mu)/dt = (q/m) mu x B_BMT.  The g-factor is a PER-PARTICLE argument, not a
// constant: it is 2.0023 for a lepton but 5.5857 for a proton, so a pair need
// not share one.  Passing the measured value rather than the classical
// point-dipole g=1 makes the anomalous precession a QED input.
Vec3 thomasBmtEffectiveField(const Vec3& velocity,
                             const ElectromagneticField& field,
                             double gFactor) {
    const double anomaly = 0.5 * (gFactor - 2.0);
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

struct MutualForces { Vec3 first, second; };

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
    bool first,double charge,bool radiationFieldOnly=false) {
    double emissionTime=wavefrontTime;
    ChargeKinematics source=historicalCharge(
        history,present,first,emissionTime);
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
        source=historicalCharge(history,present,first,emissionTime);
    }
    source=historicalCharge(history,present,first,emissionTime);
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
    const Vec3 centre=(state.firstPosition+state.secondPosition)*0.5;
    const double sourceExtent=std::max(
        (state.firstPosition-centre).norm(),
        (state.secondPosition-centre).norm());
    // The observation event is shifted back by the source radius so every
    // direction samples only the available causal history. Directional
    // retardation across the pair then retains E3, M2, toroidal terms and all
    // their interference without assigning convention-dependent pieces.
    const double wavefrontTime=state.time-sourceExtent/c;
    FieldFluxRates rates;
    for(const SphereQuadraturePoint& point:quadrature) {
        const Vec3 normal=point.direction;
        const Vec3 observationPosition = centre+normal*sampling.controlRadius;
        const ElectromagneticField firstField=farZoneChargeField(
            observationPosition,normal,wavefrontTime,centre,history,state,
            true,firstCharge,sampling.radiationFieldOnly);
        const ElectromagneticField secondField=farZoneChargeField(
            observationPosition,normal,wavefrontTime,centre,history,state,
            false,secondCharge,sampling.radiationFieldOnly);
        const Vec3 electric = firstField.electric + secondField.electric;
        const Vec3 magnetic = firstField.magnetic + secondField.magnetic;
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
    ElectromagneticField atFirst, atSecond;
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

struct DipoleDerivatives { Vec3 first, second; };
DipoleDerivatives thomasBmtDipoleDerivatives(
    const State& s,const StateHistory& history);

State historicalState(const StateHistory& history, const State& present,
                      double time) {
    if (history.empty()) return present;
    const State& earliest = history.front();
    if (history.size() == 1) {
        State extrapolated = earliest;
        const double offset = time - earliest.time;
        extrapolated.firstPosition += earliest.firstVelocity*offset;
        extrapolated.secondPosition += earliest.secondVelocity*offset;
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

// Position, velocity and the ONE lab dipole a retarded source sample needs.
// historicalState() rebuilds an entire State for this, and interpolateState()
// runs interpolateDipole() on all four dipole vectors -- three of which are
// discarded here.  With five stencil samples per call that was the single
// largest cost in the profile (interpolateDipole alone: 25% of runtime, three
// square roots per invocation).  The interpolation performed on the fields
// that ARE used is identical to interpolateState's, so results are unchanged.
struct RetardedSourceSample { Vec3 position, velocity, moment; };

RetardedSourceSample historicalSource(const StateHistory& history,
                                      const State& present,
                                      bool sourceIsFirst, double time) {
    const auto pick=[sourceIsFirst](const State& state) {
        return RetardedSourceSample{
            sourceIsFirst?state.firstPosition:state.secondPosition,
            sourceIsFirst?state.firstVelocity:state.secondVelocity,
            sourceIsFirst?state.firstDipole:state.secondDipole};
    };
    if(history.empty()) return pick(present);
    const State& earliest=history.front();
    if(history.size()==1) {
        RetardedSourceSample sample=pick(earliest);
        sample.position+=sample.velocity*(time-earliest.time);
        return sample;
    }
    const auto newer=std::lower_bound(history.begin(),history.end(),time,
        [](const State& state,double requestedTime) {
            return state.time<requestedTime;
        });
    if(newer==history.begin()) return pick(*newer);
    const State& older=(newer!=history.end())?*std::prev(newer):history.back();
    const State& target=(newer!=history.end())?*newer:present;
    const double fraction=(time-older.time)
        /std::max(target.time-older.time,1.0e-300);
    const RetardedSourceSample a=pick(older);
    const RetardedSourceSample b=pick(target);
    return {interpolateVector(a.position,b.position,fraction),
            interpolateVector(a.velocity,b.velocity,fraction),
            interpolateDipole(a.moment,b.moment,fraction)};
}

RetardedDipoleKinematics historicalDipoleKinematics(
    const StateHistory& history, const State& present, bool sourceIsFirst,
    double time) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    // The widest branch below samples middle .. middle-4h.
    derivativeStep=boundedDerivativeStep(history,time,derivativeStep,4);
    if(!(derivativeStep>0.0)) {
        const RetardedSourceSample only=
            historicalSource(history,present,sourceIsFirst,time);
        return {only.position,only.velocity,only.moment,{},{},{}};
    }
    const RetardedSourceSample middle=
        historicalSource(history,present,sourceIsFirst,time);
    const auto moment=[](const RetardedSourceSample& sample) -> const Vec3& {
        return sample.moment;
    };
    Vec3 first,second,third;
    if(time+derivativeStep>present.time) {
        const RetardedSourceSample before=historicalSource(
            history,present,sourceIsFirst,time-derivativeStep);
        const RetardedSourceSample twiceBefore=historicalSource(
            history,present,sourceIsFirst,time-2.0*derivativeStep);
        const RetardedSourceSample threeBefore=historicalSource(
            history,present,sourceIsFirst,time-3.0*derivativeStep);
        const RetardedSourceSample fourBefore=historicalSource(
            history,present,sourceIsFirst,time-4.0*derivativeStep);
        first=(moment(middle)*3.0-moment(before)*4.0+moment(twiceBefore))
             /(2.0*derivativeStep);
        second=(moment(middle)-moment(before)*2.0+moment(twiceBefore))
              /(derivativeStep*derivativeStep);
        third=(moment(middle)*5.0-moment(before)*18.0
              +moment(twiceBefore)*24.0-moment(threeBefore)*14.0
              +moment(fourBefore)*3.0)
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    } else if(time-derivativeStep<history.front().time) {
        const RetardedSourceSample after=historicalSource(
            history,present,sourceIsFirst,time+derivativeStep);
        const RetardedSourceSample twiceAfter=historicalSource(
            history,present,sourceIsFirst,time+2.0*derivativeStep);
        const RetardedSourceSample threeAfter=historicalSource(
            history,present,sourceIsFirst,time+3.0*derivativeStep);
        const RetardedSourceSample fourAfter=historicalSource(
            history,present,sourceIsFirst,time+4.0*derivativeStep);
        first=(moment(after)*4.0-moment(middle)*3.0-moment(twiceAfter))
             /(2.0*derivativeStep);
        second=(moment(twiceAfter)-moment(after)*2.0+moment(middle))
              /(derivativeStep*derivativeStep);
        third=(moment(middle)*-5.0+moment(after)*18.0
              -moment(twiceAfter)*24.0+moment(threeAfter)*14.0
              -moment(fourAfter)*3.0)
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    } else {
        const RetardedSourceSample before=historicalSource(
            history,present,sourceIsFirst,time-derivativeStep);
        const RetardedSourceSample after=historicalSource(
            history,present,sourceIsFirst,time+derivativeStep);
        const RetardedSourceSample twiceBefore=historicalSource(
            history,present,sourceIsFirst,time-2.0*derivativeStep);
        const RetardedSourceSample twiceAfter=historicalSource(
            history,present,sourceIsFirst,time+2.0*derivativeStep);
        first=(moment(after)-moment(before))/(2.0*derivativeStep);
        second=(moment(after)-moment(middle)*2.0+moment(before))
              /(derivativeStep*derivativeStep);
        third=(moment(twiceAfter)-moment(after)*2.0
              +moment(before)*2.0-moment(twiceBefore))
             /(2.0*derivativeStep*derivativeStep*derivativeStep);
    }
    return {middle.position,middle.velocity,
            moment(middle),first,second,third};
}

RetardedElectricDipoleKinematics historicalElectricDipoleKinematics(
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double time) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    // The widest branch below samples middle .. middle-2h.
    derivativeStep=boundedDerivativeStep(history,time,derivativeStep,2);
    const auto sample=[&](double sampleTime) {
        const State state=historicalState(history,present,sampleTime);
        return sourceIsFirst?state.firstElectricDipole
                               :state.secondElectricDipole;
    };
    const State middle=historicalState(history,present,time);
    if(!(derivativeStep>0.0)) {
        return {sourceIsFirst?middle.firstPosition:middle.secondPosition,
                sourceIsFirst?middle.firstVelocity:middle.secondVelocity,
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
    return {sourceIsFirst?middle.firstPosition:middle.secondPosition,
            sourceIsFirst?middle.firstVelocity:middle.secondVelocity,
            moment,first,second};
}

struct DipoleRadiationReaction {
    double power=0.0;
    Vec3 momentumRate,angularMomentumRate;
    Vec3 firstTorque,secondTorque;
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
    // Centre of mass as origin.  For a NEUTRAL pair the quadrupole is origin
    // independent only when the dipole moment also vanishes; this pair carries
    // d = q1 r1 + q2 r2 != 0, so where the expansion is centred is a physical
    // choice and not a convention.  The centre of mass is the origin in which
    // the two-body multipole series is normally written, and it gives
    //
    //     Q_ij = kappa (3 d_i d_j - d^2 delta_ij),
    //     kappa = (q1 m2^2 + q2 m1^2)/(m1+m2)^2,   d = r1 - r2.
    //
    // The equal-mass midpoint (r1+r2)/2 used to stand here, and it silently
    // deleted the whole channel.  About that point x2 = -x1, and 3 x_i x_j -
    // r^2 delta_ij is QUADRATIC in x, so the two contributions add instead of
    // cancelling and the sum collapses to (q1 + q2)(3 x1_i x1_j - x1^2 delta),
    // which is identically zero for any neutral pair -- exactly zero in
    // floating point too, since +e and -e are one literal with a flipped sign.
    //
    // That reproduced the right answer for e+e- for the wrong reason.  The
    // real reason positronium has no E2 channel is that kappa vanishes when
    // the pair is MASS symmetric, which the midpoint origin never consults:
    // it returned zero for p+e- as well, where kappa = -0.9989 e and the
    // channel is at nearly full strength.  Mass-symmetric pairs still give
    // exactly zero here, so e+e- and mu+mu- are unaffected.
    const Vec3 origin=(state.firstPosition*firstMass
                      +state.secondPosition*secondMass)
                     /(firstMass+secondMass);
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
    accumulate(firstCharge,state.firstPosition);
    accumulate(secondCharge,state.secondPosition);
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

// d = q1 r1 + q2 r2.  The old form e*(r2 - r1) is the same thing only when the
// charges are equal and opposite AND the centre of mass is at the origin; it
// is written out here so neither condition is silently assumed.
Vec3 electricDipoleMoment(const State& state) {
    return state.firstPosition*firstCharge+state.secondPosition*secondCharge;
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
    return {reactionField*firstCharge,reactionField*secondCharge};
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
    const RetardedDipoleKinematics first=historicalDipoleKinematics(
        history,state,true,state.time);
    const RetardedDipoleKinematics second=historicalDipoleKinematics(
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
    accumulate(first,result.firstTorque);
    accumulate(second,result.secondTorque);
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
    Vec3 firstDipoleTorque, secondDipoleTorque;
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
        ll.first.norm()/std::max(externalForces.first.norm(),1.0e-300),
        ll.second.norm()/std::max(externalForces.second.norm(),1.0e-300));
    if(computeOutwardFlux)
        result.outwardFlux=electromagneticFieldFluxRates(state,history);

    const Vec3 firstAcceleration = relativisticAcceleration(
        state.firstVelocity, externalForces.first, firstMass);
    const Vec3 secondAcceleration = relativisticAcceleration(
        state.secondVelocity, externalForces.second, secondMass);
    const Vec3 electricDipoleSecondDerivative =
        firstAcceleration*firstCharge+secondAcceleration*secondCharge;
    result.leadingElectricDipolePower =
        electricDipoleSecondDerivative.squaredNorm()
        / (6.0*pi*epsilon0*c*c*c);
    if(needsBlendingGates) {
        result.electricQuadrupolePower =
            electricQuadrupoleRadiatedPower(state,history);
    }

    const DipoleRadiationReaction magnetic =
        dipoleRadiationReaction(state, history);
    result.firstDipoleTorque = magnetic.firstTorque;
    result.secondDipoleTorque = magnetic.secondTorque;
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
                                         magnitude(pairDipoleCharge)*nuclearCutoff);
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
        // Measured, not assumed: correcting the quadrupole origin to the
        // centre of mass moves this gate's INPUT by ten orders of magnitude
        // for an unequal-mass pair, and changes its VERDICT nowhere.
        //
        // For p+e- the electric quadrupole went from 1.6e-30 W (the old
        // midpoint origin, which forced Q == 0 identically and left only
        // rounding) to 1.4e-11 W, dropping dominance from 3.5e13 to 3.3e3.
        // The gate saturates at 1 for anything above 20, so both sit deep in
        // the saturated region and E1 stays overwhelmingly dominant either
        // way.  A 30-event experiment 5 run under --radiation-reaction
        // automatic is byte-identical with and without the correction.
        //
        // For E2 to matter here it would have to come within a factor of 20
        // of E1, which needs a far more compact and relativistic source than
        // this model integrates.  The correction is therefore right and
        // currently inert -- worth stating so nobody re-derives it.
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
                    ll.first*(1.0-result.coherentWeight)
                        +coherent.first*result.coherentWeight,
                    ll.second*(1.0-result.coherentWeight)
                        +coherent.second*result.coherentWeight};
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
        && isFinite(response.chargeReaction.first)
        && isFinite(response.chargeReaction.second)
        && isFinite(response.firstDipoleTorque)
        && isFinite(response.secondDipoleTorque)
        && isFinite(response.outwardFlux.momentum)
        && isFinite(response.outwardFlux.angularMomentum);
}

void applyDipoleRadiationTorque(State& state,
                                const ParticleMultipoleRadiation& reaction,
                                double dt) {
    const double firstGyromagneticRatio=firstCharge/(2.0*firstMass);
    const double secondGyromagneticRatio=secondCharge/(2.0*secondMass);
    synchronizeCovariantDipoles(state);
    const double firstNorm=state.firstProperDipole.norm();
    const double secondNorm=state.secondProperDipole.norm();
    state.firstProperDipole+=reaction.firstDipoleTorque
                         *(firstGyromagneticRatio*dt);
    state.secondProperDipole+=reaction.secondDipoleTorque
                         *(secondGyromagneticRatio*dt);
    if(state.firstProperDipole.norm()>0.0)
        state.firstProperDipole=state.firstProperDipole
            *(firstNorm/state.firstProperDipole.norm());
    if(state.secondProperDipole.norm()>0.0)
        state.secondProperDipole=state.secondProperDipole
            *(secondNorm/state.secondProperDipole.norm());
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
    const StateHistory& history,const State& present,bool sourceIsFirst) {
    ChargeKinematics source=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    double retardedTime=observationTime
        -(observationPosition-source.position).norm()/c;
    for(int iteration=0;iteration<16;++iteration) {
        source=historicalCharge(history,present,sourceIsFirst,retardedTime);
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
            history,present,sourceIsFirst,retardedTime);
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
    const StateHistory& history,const State& present,bool sourceIsFirst) {
    ChargeKinematics source=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    double retardedTime=observationTime
        -(observationPosition-source.position).norm()/c;
    for(int iteration=0;iteration<16;++iteration) {
        source=historicalCharge(history,present,sourceIsFirst,retardedTime);
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
        history,present,sourceIsFirst,retardedTime);
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
    const Vec3 first = geometry.firstMinusSecond
                        * (-pairCoulombStrength * geometry.inverseDistanceCubed);
    return {first, first * -1.0};
}

MutualForces mutualForces(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const MutualForces electrostatic = coulombForces(s);
    const Vec3 dipoleOnFirst = regularizedDipoleForce(
        geometry.firstMinusSecond, s.firstDipole, s.secondDipole);
    return {electrostatic.first + dipoleOnFirst,
            electrostatic.second - dipoleOnFirst};
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
    const Vec3 firstLeadingAcceleration = leading.first / firstMass;
    const Vec3 secondLeadingAcceleration = leading.second / secondMass;
    const double chargeProduct = pairChargeProduct;
    return {
        darwinForceOnFirst(s.firstVelocity, s.secondVelocity,
                           secondLeadingAcceleration,
                           geometry.firstMinusSecond, chargeProduct),
        darwinForceOnFirst(s.secondVelocity, s.firstVelocity,
                           firstLeadingAcceleration,
                           geometry.firstMinusSecond * -1.0, chargeProduct)
    };
}

double darwinInteractionEnergy(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 n = geometry.firstMinusSecond * geometry.inverseDistance;
    const double chargeProduct = pairChargeProduct;
    return coulomb * chargeProduct * geometry.inverseDistance / (2.0 * c*c)
         * (dot(s.firstVelocity, s.secondVelocity)
          + dot(s.firstVelocity, n) * dot(s.secondVelocity, n));
}

LocalElectromagneticFields localRelativisticFields(
    const State& s, const StateHistory& history) {
    ElectromagneticField atFirst = lienardWiechertField(
        s.firstPosition, s.time, history, s, false, secondCharge);
    ElectromagneticField atSecond = lienardWiechertField(
        s.secondPosition, s.time, history, s, true, firstCharge);
    const ElectromagneticField secondDipole=retardedMagneticDipoleField(
        s.firstPosition,s.time,history,s,false);
    const ElectromagneticField firstDipole=retardedMagneticDipoleField(
        s.secondPosition,s.time,history,s,true);
    const ElectromagneticField secondElectricDipole=
        retardedElectricDipoleField(
            s.firstPosition,s.time,history,s,false);
    const ElectromagneticField firstElectricDipole=
        retardedElectricDipoleField(
            s.secondPosition,s.time,history,s,true);
    atFirst.electric+=secondDipole.electric
        +secondElectricDipole.electric;
    atFirst.magnetic+=secondDipole.magnetic
        +secondElectricDipole.magnetic;
    atSecond.electric+=firstDipole.electric
        +firstElectricDipole.electric;
    atSecond.magnetic+=firstDipole.magnetic
        +firstElectricDipole.magnetic;
    // The external field is uniform, so both roles see the same addition and
    // no gradient force follows from it.  Adding it here rather than only in
    // the force sums is what carries it into Thomas-BMT precession, which for
    // a field this weak is the channel that actually does something.
    atFirst.magnetic+=gExternalMagneticField;
    atSecond.magnetic+=gExternalMagneticField;
    return {atFirst, atSecond};
}

DipoleDerivatives thomasBmtDipoleDerivatives(
    const State& s, const StateHistory& history) {
    const LocalElectromagneticFields fields = localRelativisticFields(s, history);
    const Vec3 firstEffectiveField = thomasBmtEffectiveField(
        s.firstVelocity, fields.atFirst, firstGFactor);
    const Vec3 secondEffectiveField = thomasBmtEffectiveField(
        s.secondVelocity, fields.atSecond, secondGFactor);
    return {
        cross(s.firstDipole, firstEffectiveField)
            * (firstCharge/firstMass),
        cross(s.secondDipole, secondEffectiveField)
            * (secondCharge/secondMass)
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
                         double chargeToMass,double laboratoryDt,
                         double gFactor) {
    const double targetNorm=properDipole.norm();
    if(targetNorm==0.0||laboratoryDt==0.0) return properDipole;
    const double halfG=0.5*gFactor;
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
    const Vec3 firstDipole=advanceCovariantBmt(s.firstProperDipole,
        s.firstVelocity,fields.atFirst,firstCharge/firstMass,dt,firstGFactor);
    const Vec3 secondDipole=advanceCovariantBmt(s.secondProperDipole,
        s.secondVelocity,fields.atSecond,secondCharge/secondMass,dt,
        secondGFactor);
    // Update simultaneously so neither particle sees an already-updated peer.
    s.firstProperDipole=firstDipole;
    s.secondProperDipole=secondDipole;
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
    const Vec3 firstMinusSecond = s.firstPosition - s.secondPosition;
    const Vec3 firstMinusSecondVelocity =
        s.firstVelocity - s.secondVelocity;

    const ChargeDipolePairForces firstChargeSecondDipole = chargeDipolePairForces(
        firstMinusSecondVelocity, firstCharge, s.secondDipole,
        derivatives.second, firstMinusSecond);
    const ChargeDipolePairForces secondChargeFirstDipole = chargeDipolePairForces(
        firstMinusSecondVelocity * -1.0, secondCharge, s.firstDipole,
        derivatives.first, firstMinusSecond * -1.0);

    return {firstChargeSecondDipole.onCharge + secondChargeFirstDipole.onDipole,
            secondChargeFirstDipole.onCharge + firstChargeSecondDipole.onDipole};
}

MutualForces allExternalForces(const State& s) {
    const MutualForces positionForces = mutualForces(s);
    const MutualForces velocityForces = darwinForces(s);
    const StateHistory localHistory{State{s}};
    const MutualForces mixedMagneticForces = chargeDipoleForces(s, localHistory);
    const MutualForces externalField{
        lorentzForce(firstCharge, s.firstVelocity, {{}, gExternalMagneticField}),
        lorentzForce(secondCharge, s.secondVelocity, {{}, gExternalMagneticField})};
    return {positionForces.first + velocityForces.first
                + mixedMagneticForces.first + externalField.first,
            positionForces.second + velocityForces.second
                + mixedMagneticForces.second + externalField.second};
}

ElectromagneticField fieldFromOtherParticleAt(
    const Vec3& observationPosition,const State& state,
    const StateHistory& history,bool targetIsFirst) {
    const bool sourceIsFirst=!targetIsFirst;
    const double sourceCharge=sourceIsFirst?firstCharge:secondCharge;
    ElectromagneticField field=lienardWiechertField(
        observationPosition,state.time,history,state,sourceIsFirst,
        sourceCharge);
    const ElectromagneticField magneticDipole=retardedMagneticDipoleField(
        observationPosition,state.time,history,state,sourceIsFirst);
    const ElectromagneticField electricDipole=retardedElectricDipoleField(
        observationPosition,state.time,history,state,sourceIsFirst);
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
                                  bool targetIsFirst) {
    const Vec3 targetPosition=targetIsFirst?state.firstPosition
                                               :state.secondPosition;
    const Vec3 targetVelocity=targetIsFirst?state.firstVelocity
                                               :state.secondVelocity;
    const Vec3 properDipole=targetIsFirst?state.firstProperDipole
                                             :state.secondProperDipole;
    if(properDipole.squaredNorm()==0.0) return {};
    const double gradientStep=std::max(
        1.0e-6*separation(state),1.0e-3*nuclearCutoff);
    const auto coupling=[&](const Vec3& point) {
        return dot(properDipole,magneticFieldInRestFrame(
            fieldFromOtherParticleAt(point,state,history,targetIsFirst),
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
    const ElectromagneticField secondField = lienardWiechertField(
        s.firstPosition, s.time, history, s, false, secondCharge);
    const ElectromagneticField firstField = lienardWiechertField(
        s.secondPosition, s.time, history, s, true, firstCharge);
    const ElectromagneticField secondDipoleField=retardedMagneticDipoleField(
        s.firstPosition,s.time,history,s,false);
    const ElectromagneticField firstDipoleField=retardedMagneticDipoleField(
        s.secondPosition,s.time,history,s,true);
    const ElectromagneticField secondElectricDipoleField=
        retardedElectricDipoleField(
            s.firstPosition,s.time,history,s,false);
    const ElectromagneticField firstElectricDipoleField=
        retardedElectricDipoleField(
            s.secondPosition,s.time,history,s,true);
    const MutualForces chargeCharge{
        lorentzForce(firstCharge,s.firstVelocity,
            {secondField.electric+secondDipoleField.electric
                +secondElectricDipoleField.electric,
             secondField.magnetic+secondDipoleField.magnetic
                +secondElectricDipoleField.magnetic}),
        lorentzForce(secondCharge,s.secondVelocity,
            {firstField.electric+firstDipoleField.electric
                +firstElectricDipoleField.electric,
             firstField.magnetic+firstDipoleField.magnetic
                +firstElectricDipoleField.magnetic})};

    const MutualForces tensorGradient{
        covariantDipoleGradientForce(s,history,true),
        covariantDipoleGradientForce(s,history,false)};
    // Same uniform external field as in the instantaneous sum.  It is not
    // retarded because it is not sourced by either particle.
    const MutualForces externalField{
        lorentzForce(firstCharge,s.firstVelocity,{{},gExternalMagneticField}),
        lorentzForce(secondCharge,s.secondVelocity,{{},gExternalMagneticField})};
    return {chargeCharge.first+tensorGradient.first+externalField.first,
            chargeCharge.second+tensorGradient.second+externalField.second};
}

struct CanonicalMomenta { Vec3 first, second; };

// Canonical momenta obtained from the conservative low-velocity action.  The
// q-mu contributions occur with opposite signs on the charge and on the
// dipole carrier; the Darwin terms do not cancel particle by particle.
CanonicalMomenta canonicalMomenta(const State& s) {
    const PairGeometry geometry = pairGeometry(s);
    const Vec3 n = geometry.firstMinusSecond * geometry.inverseDistance;
    const double chargeProduct = pairChargeProduct;
    const double darwinCoefficient = coulomb * chargeProduct
        * geometry.inverseDistance / (2.0 * c*c);

    Vec3 first = momentum(s.firstVelocity, firstMass)
        + (s.secondVelocity + n * dot(s.secondVelocity, n)) * darwinCoefficient;
    Vec3 second = momentum(s.secondVelocity, secondMass)
        + (s.firstVelocity + n * dot(s.firstVelocity, n)) * darwinCoefficient;

    const Vec3 secondPotentialAtFirst = regularizedDipoleVectorPotential(
        geometry.firstMinusSecond, s.secondDipole);
    const Vec3 firstPotentialAtSecond = regularizedDipoleVectorPotential(
        geometry.firstMinusSecond * -1.0, s.firstDipole);
    const Vec3 mixedCanonicalContribution =
        secondPotentialAtFirst * firstCharge
      + firstPotentialAtSecond * secondCharge;
    first += mixedCanonicalContribution;
    second += mixedCanonicalContribution * -1.0;
    return {first, second};
}

// Overloads taking an already-computed CanonicalMomenta.  Callers that need
// more than one of these quantities for the same state used to recompute the
// canonical momenta once per call: three times per frame and four times per
// accepted integration step, each redoing the pair geometry, the relativistic
// momenta and two regularized dipole vector potentials.
Vec3 noetherMomentum(const CanonicalMomenta& canonical) {
    return canonical.first + canonical.second;
}

Vec3 noetherMomentum(const State& s) {
    return noetherMomentum(canonicalMomenta(s));
}

double canonicalMomentumScale(const CanonicalMomenta& canonical) {
    return canonical.first.norm() + canonical.second.norm();
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
    const double firstGyromagneticRatio = firstCharge / (2.0 * firstMass);
    const double secondGyromagneticRatio = secondCharge / (2.0 * secondMass);
    const Vec3 orbital = cross(s.firstPosition, canonical.first)
                       + cross(s.secondPosition, canonical.second);
    const Vec3 intrinsic = s.firstDipole / firstGyromagneticRatio
                         + s.secondDipole / secondGyromagneticRatio;
    return orbital + intrinsic;
}

Vec3 noetherAngularMomentum(const State& s) {
    return noetherAngularMomentum(s, canonicalMomenta(s));
}

double conservativeParticleEnergy(const State& state) {
    const PairGeometry geometry=pairGeometry(state);
    const double kinetic=(gamma(state.firstVelocity)-1.0)*firstMass*c*c
        +(gamma(state.secondVelocity)-1.0)*secondMass*c*c;
    const double coulombPotential=-pairCoulombStrength
        *geometry.inverseDistance;
    const double dipolePotential=-dot(state.firstDipole,
        regularizedDipoleField(
            geometry.firstMinusSecond,state.secondDipole));
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

// chargeCoupling is q_i^2 for a self term and q_i q_j for the mutual one; it
// used to be a hard-wired e^2, which made both correct only for unit charges.
Vec3 reducedOrderSelfForce(const Vec3& velocity,
                           const FourVector& beforeForce,
                           const FourVector& afterForce,
                           double laboratoryDerivativeStep,
                           double mass,
                           double chargeCoupling) {
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
    const double characteristicTime = chargeCoupling
        / (6.0 * pi * epsilon0 * mass * c*c*c);
    // The spatial component of a four-force is gamma times the ordinary
    // laboratory three-force.
    return orthogonalDerivative.space * (characteristicTime / relativisticGamma);
}

MutualForces individualLandauLifshitzSelfForces(
    const State& s, const MutualForces& external,
    const StateHistory& history) {
    const Vec3 firstAcceleration = relativisticAcceleration(s.firstVelocity, external.first,
                                                                firstMass);
    const Vec3 secondAcceleration = relativisticAcceleration(s.secondVelocity, external.second,
                                                                secondMass);
    const double relativeSpeed = (s.firstVelocity - s.secondVelocity).norm();
    const double derivativeStep = std::max(1.0e-24,
        1.0e-4 * std::min(separation(s) / std::max(relativeSpeed, 1.0),
                          separation(s) / c));
    State before = s;
    State after = s;
    before.firstPosition = s.firstPosition - s.firstVelocity * derivativeStep
        + firstAcceleration * (0.5 * derivativeStep * derivativeStep);
    before.secondPosition = s.secondPosition - s.secondVelocity * derivativeStep
        + secondAcceleration * (0.5 * derivativeStep * derivativeStep);
    after.firstPosition = s.firstPosition + s.firstVelocity * derivativeStep
        + firstAcceleration * (0.5 * derivativeStep * derivativeStep);
    after.secondPosition = s.secondPosition + s.secondVelocity * derivativeStep
        + secondAcceleration * (0.5 * derivativeStep * derivativeStep);
    before.firstVelocity = velocityFromMomentum(
        momentum(s.firstVelocity, firstMass)
            - external.first * derivativeStep, firstMass);
    before.secondVelocity = velocityFromMomentum(
        momentum(s.secondVelocity, secondMass)
            - external.second * derivativeStep, secondMass);
    after.firstVelocity = velocityFromMomentum(
        momentum(s.firstVelocity, firstMass)
            + external.first * derivativeStep, firstMass);
    after.secondVelocity = velocityFromMomentum(
        momentum(s.secondVelocity, secondMass)
            + external.second * derivativeStep, secondMass);
    before.time -= derivativeStep;
    after.time += derivativeStep;
    const MutualForces beforeForces = retardedExternalForces(before, history);
    const MutualForces afterForces = retardedExternalForces(after, history);
    // Self terms: F_i = (q_i^2/(6 pi eps0 c^3)) a_i-dot, in reduced order.
    const Vec3 firstSelf=reducedOrderSelfForce(s.firstVelocity,
        fourForce(before.firstVelocity, beforeForces.first),
        fourForce(after.firstVelocity, afterForces.first),
        derivativeStep, firstMass, firstCharge*firstCharge);
    const Vec3 secondSelf=reducedOrderSelfForce(s.secondVelocity,
        fourForce(before.secondVelocity, beforeForces.second),
        fourForce(after.secondVelocity, afterForces.second),
        derivativeStep, secondMass, secondCharge*secondCharge);

    // MUTUAL terms: each charge also sits in the radiation field of the other,
    // F_i = (q_i q_j/(6 pi eps0 c^3)) a_j-dot.  Decomposing the coherent power
    //
    //     |q_1 a_1 + q_2 a_2|^2 = |q_1 a_1|^2 + |q_2 a_2|^2 + 2 q_1 q_2 a_1.a_2
    //
    // the self forces above reproduce only the first two terms.  For e+e- the
    // interference term equals their sum, so leaving it out cost exactly a
    // factor of two in the orbital energy loss -- measured independently three
    // ways before this was added: the orbit-averaged power against Larmor
    // (0.4984), the reaction work against Larmor (0.4944), and the
    // energy-balance residual |dE_LL-vs-flux|/E_rad (0.501).
    //
    // The partner's contribution is the same object evaluated with the
    // PARTNER's four-force and mass, and with q_i q_j as its charge coupling
    // instead of q_i^2.  The orthogonal projection is taken against the velocity of the
    // particle being acted on, which is what keeps the four-force orthogonal
    // to its own four-velocity; at these speeds the distinction is O(beta^2).
    const Vec3 firstMutual=reducedOrderSelfForce(s.firstVelocity,
        fourForce(before.secondVelocity, beforeForces.second),
        fourForce(after.secondVelocity, afterForces.second),
        derivativeStep, secondMass, pairChargeProduct);
    const Vec3 secondMutual=reducedOrderSelfForce(s.secondVelocity,
        fourForce(before.firstVelocity, beforeForces.first),
        fourForce(after.firstVelocity, afterForces.first),
        derivativeStep, firstMass, pairChargeProduct);

    return {firstSelf+firstMutual, secondSelf+secondMutual};
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
    Vec3 firstMomentum = momentum(s.firstVelocity, firstMass)
        + (forces.first + radiation.chargeReaction.first) * (0.5 * dt);
    Vec3 secondMomentum = momentum(s.secondVelocity, secondMass)
        + (forces.second + radiation.chargeReaction.second) * (0.5 * dt);
    State trial = s;
    trial.time += dt;
    trial.firstVelocity = velocityFromMomentum(firstMomentum, firstMass);
    trial.secondVelocity = velocityFromMomentum(secondMomentum, secondMass);
    firstMomentum = momentum(trial.firstVelocity, firstMass);
    secondMomentum = momentum(trial.secondVelocity, secondMass);
    trial.firstPosition += trial.firstVelocity * dt;
    trial.secondPosition += trial.secondVelocity * dt;
    trial.firstAcceleration = relativisticAcceleration(trial.firstVelocity, forces.first, firstMass);
    trial.secondAcceleration = relativisticAcceleration(trial.secondVelocity, forces.second, secondMass);

    const MutualForces trialForces = retardedExternalForces(trial, history);
    const ParticleMultipoleRadiation trialRadiation =
        particleMultipoleRadiation(
            trial,trialForces,history,computeOutwardFlux,reactionModel);
    if (!finiteRadiationResponse(trialRadiation)) {
        s.time = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    firstMomentum += (trialForces.first
        + trialRadiation.chargeReaction.first) * (0.5 * dt);
    secondMomentum += (trialForces.second
        + trialRadiation.chargeReaction.second) * (0.5 * dt);
    trial.firstVelocity = velocityFromMomentum(firstMomentum, firstMass);
    trial.secondVelocity = velocityFromMomentum(secondMomentum, secondMass);
    applyDipolePrecession(trial, 0.5 * dt, history);
    applyDipoleRadiationTorque(trial,trialRadiation,0.5*dt);
    trial.firstAcceleration = relativisticAcceleration(trial.firstVelocity, trialForces.first, firstMass);
    trial.secondAcceleration = relativisticAcceleration(trial.secondVelocity, trialForces.second, secondMass);
    // Flux bookkeeping is accumulated from the COMMITTED state only, never
    // from `trial`.  At this point `trial` is a predictor: its velocity has
    // taken a half kick and its position was advanced with that new velocity,
    // so the (x,v) pair is not a consistent sample of the trajectory.
    // interpolatedCharge() reconstructs the source acceleration as the second
    // derivative of a cubic Hermite through the two bracketing samples, and
    // feeding it that inconsistent pair makes the reconstruction wrong by a
    // factor (4-6s)/2 -- at the wavefront offset these evaluations actually
    // use, -0.554.  The far-zone field is linear in that acceleration and the
    // Poynting flux quadratic, so the trial-state flux came out at 0.31 of the
    // true value and the trapezoid under-reported the radiated energy by 34%,
    // uniformly, from the first step.  Measured against the integral of the
    // flux over the same trajectory: the trapezoid gave 5.470e-05 eV where
    // the correct answer is 8.268e-05 eV; the committed-state form below gives
    // 8.274e-05 eV, i.e. 0.08%.
    //
    // This is first order in dt rather than second, but advanceAdaptive()
    // evaluates it at the start of each half-step, so a full step samples the
    // start and the midpoint.  None of it feeds back into the trajectory.
    const double radiatedEnergyIncrement = radiation.outwardFlux.energy*dt;
    trial.radiatedEnergy += radiatedEnergyIncrement;
    // Magnetic-dipole damping changes the constrained internal dipole sector.
    // The charge part is already represented by the particle self-force and
    // its near-field (Schott) term.
    const double dipoleRadiatedEnergy=radiation.magneticDipolePower*dt;
    trial.dipoleConstraintEnergy-=dipoleRadiatedEnergy;
    trial.radiatedMomentum += radiation.outwardFlux.momentum*dt;
    trial.radiatedAngularMomentum += radiation.outwardFlux.angularMomentum*dt;

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
            mismatch.energy+=dot(response.chargeReaction.first,
                                  state.firstVelocity)
                            +dot(response.chargeReaction.second,
                                  state.secondVelocity);
            const Vec3 reactionForce=response.chargeReaction.first
                                    +response.chargeReaction.second;
            mismatch.momentum+=reactionForce;
            mismatch.angularMomentum+=cross(state.firstPosition,
                                             response.chargeReaction.first)
                +cross(state.secondPosition,response.chargeReaction.second);
            return mismatch;
        };
        // Committed state only, for the same reason as the flux above: the
        // trial-state rate is built from the same corrupted reconstruction.
        const FieldFluxRates initialMismatch=
            chargeMismatchRates(balanceStart,radiation);
        trial.reactionEnergyMismatch=balanceStart.reactionEnergyMismatch
            +initialMismatch.energy*dt;
        trial.reactionMomentumMismatch=balanceStart.reactionMomentumMismatch
            +initialMismatch.momentum*dt;
        trial.reactionAngularMomentumMismatch=
            balanceStart.reactionAngularMomentumMismatch
            +initialMismatch.angularMomentum*dt;
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
    const Vec3 firstMomentum0=momentum(start.firstVelocity,firstMass);
    const Vec3 secondMomentum0=momentum(start.secondVelocity,secondMass);
    State endpoint=start;
    for(int iteration=0;iteration<std::max(iterations,1);++iteration) {
        State midpoint=start;
        midpoint.time=start.time+0.5*dt;
        midpoint.firstPosition=(start.firstPosition
                                   +endpoint.firstPosition)*0.5;
        midpoint.secondPosition=(start.secondPosition
                                   +endpoint.secondPosition)*0.5;
        const Vec3 endpointFirstMomentum=momentum(
            endpoint.firstVelocity,firstMass);
        const Vec3 endpointSecondMomentum=momentum(
            endpoint.secondVelocity,secondMass);
        midpoint.firstVelocity=velocityFromMomentum(
            (firstMomentum0+endpointFirstMomentum)*0.5,firstMass);
        midpoint.secondVelocity=velocityFromMomentum(
            (secondMomentum0+endpointSecondMomentum)*0.5,secondMass);
        midpoint.firstDipole=(start.firstDipole
                                 +endpoint.firstDipole)*0.5;
        midpoint.secondDipole=(start.secondDipole
                                 +endpoint.secondDipole)*0.5;
        const MutualForces force=allExternalForces(midpoint);
        endpoint.firstVelocity=velocityFromMomentum(
            firstMomentum0+force.first*dt,firstMass);
        endpoint.secondVelocity=velocityFromMomentum(
            secondMomentum0+force.second*dt,secondMass);
        const Vec3 midpointFirstMomentum=(firstMomentum0
            +momentum(endpoint.firstVelocity,firstMass))*0.5;
        const Vec3 midpointSecondMomentum=(secondMomentum0
            +momentum(endpoint.secondVelocity,secondMass))*0.5;
        endpoint.firstPosition=start.firstPosition
            +velocityFromMomentum(midpointFirstMomentum,firstMass)*dt;
        endpoint.secondPosition=start.secondPosition
            +velocityFromMomentum(midpointSecondMomentum,secondMass)*dt;
    }
    endpoint.time=start.time+dt;
    endpoint.firstAcceleration=(endpoint.firstVelocity
                                  -start.firstVelocity)/dt;
    endpoint.secondAcceleration=(endpoint.secondVelocity
                                  -start.secondVelocity)/dt;
    state=endpoint;
}
#endif

// The single classical trajectory engine shared by Visual and Statistical.
// Sampling policies and stopping conditions belong to their callers; the
// equations of motion, retarded history and radiation bookkeeping do not.
