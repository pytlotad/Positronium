#pragma once

// Particle-particle interactions, radiation bookkeeping and the shared
// relativistic integration step.
//
// Self-contained and order-independent.  It names what it needs through a
// using-directive on positronium::parameters and using-declarations for the
// object types, rather than reopening namespace positronium: the header is
// still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.
//
// The MaxwellBlock in particleFieldTotals below is why this header includes
// the validation backend, and is also why the grid-coupled pushers do NOT
// live there: they need thomasBmtEffectiveField from here, which would close
// a cycle.  They sit in maxwell_validation.hpp, above both.

#include "dipole_tensor.hpp"
#include "dual_number.hpp"
#include "pair_configuration.hpp"
#include "pair_geometry.hpp"
#include "physical_constants.hpp"
#include "relativistic_field_types.hpp"
#include "relativistic_kinematics.hpp"
#include "retarded_charge_kinematics.hpp"
#include "state.hpp"
#include "state_validity_interpolation.hpp"
#include "two_body_kinematics.hpp"
#include "vector3.hpp"
#include "zero_point_field.hpp"

// Only particleFieldTotals below needs the Maxwell backend, and it exists only
// in the validation build.  Pulling the whole Yee/AMR/CPML backend into the
// production translation unit, which has never seen it, would be a change of
// substance dressed as an include cleanup.
#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
#include "maxwell_validation_backend.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <span>
#include <vector>

namespace two_body = positronium::kinematics;

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::StateHistory;
using positronium::objects::DipoleTensor;
using positronium::objects::cross;
using positronium::objects::dot;
using namespace positronium::parameters;

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
struct ParticleFieldTotals {
    double energy=0.0;
    Vec3 momentum, angularMomentum;
};

inline ParticleFieldTotals particleFieldTotals(const State& state,
                                         const MaxwellBlock& field) {
    const Vec3 firstMomentum=momentum(state.firstVelocity,firstMass);
    const Vec3 secondMomentum=momentum(state.secondVelocity,secondMass);
    const MaxwellVolumeIntegrals fieldTotals=field.volumeIntegrals();
    const double firstGyromagneticRatio=firstGyromagneticRatioOf();
    const double secondGyromagneticRatio=secondGyromagneticRatioOf();
    // Rest energies are constant and omitted from the diagnostic difference.
    return {kineticEnergy(state.firstVelocity,firstMass)
          +kineticEnergy(state.secondVelocity,secondMass)
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
inline Vec3 thomasBmtEffectiveField(const Vec3& velocity,
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

inline Vec3 relativisticAcceleration(const Vec3& velocity, const Vec3& force, double mass) {
    const double velocityForce = dot(velocity, force);
    return (force - velocity * (velocityForce / (c*c))) / (gamma(velocity) * mass);
}

inline Vec3 lorentzForce(double charge, const Vec3& velocity, const ElectromagneticField& field) {
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

// Lebedev rules: spherical quadratures invariant under the octahedral group.
// Each is a sum of complete orbits of that group, which is what integrates
// the low electromagnetic multipoles without the small preferred-axis bias a
// finite Fibonacci lattice carries.  The lattice is retained as a diagnostic
// fallback for the non-tabulated direction counts the convergence tests use.
//
//   nodes  exact through  orbits
//      26      degree  7  a1 + a2 + a3
//      50      degree 11  a1 + a2 + a3 + b
//     110      degree 17  a1 + a3 + 3b + c
//     194      degree 23  a1 + a2 + a3 + 4b + c + d
//     302      degree 29  a1 + a3 + 6b + 2c + 2d
//
// Before 194 was tabulated here the flux probes' own REFERENCE was the
// Fibonacci lattice, so the 50-node rule was being validated against a grid
// with no exactness degree at all; that is what the higher orders fix.
//
// The node parameters below are not transcribed constants taken on trust.
// They are the solution of each rule's defining moment equations, and
// sphereQuadratureExactDegree re-derives the degree at run time from the
// nodes themselves -- a mistyped digit cannot survive a degree-23 exactness
// test, and the validation suite asserts every rule's degree.
namespace sphere_quadrature {

// Lebedev's orbit names.  axis is the 6 points (1,0,0), edge the 12 points
// (1,1,0)/sqrt2, corner the 8 points (1,1,1)/sqrt3, doubled the 24 points
// (l,l,m) with m=sqrt(1-2l^2), planar the 24 points (p,q,0) with
// q=sqrt(1-p^2), and general the 48 points (r,s,t) with t=sqrt(1-r^2-s^2).
enum class Orbit { axis, edge, corner, doubled, planar, general };

struct Generator {
    Orbit orbit;
    // Unused for the three parameter-free orbits; `second` only for general.
    double first=0.0, second=0.0;
    // Per node, normalized so the whole rule sums to 1 over the sphere.
    double normalizedWeight=0.0;
};

// Expand one orbit: all distinct coordinate permutations of the seed, each
// with all eight sign choices.  The duplicate test is what makes the count
// come out right for the seeds that have a repeated or zero coordinate --
// (l,l,m) gives 24 rather than 48, (p,q,0) likewise -- so the orbit sizes
// are produced by the symmetry rather than asserted alongside it.
//
// It compares only against THIS orbit's own nodes, from `firstOfOrbit` on.
// An orbit is a set, so removing its internal repeats is the right
// operation; two different generators landing on a common node is not, and
// deduplicating across orbits would silently drop the second one's weight
// instead of adding it.  Left this way such a collision shows up as a node
// count that does not match the rule's name, which the degree test then
// fails on.
inline void appendOrbit(std::vector<SphereQuadraturePoint>& into,
                        const Generator& generator) {
    std::array<double,3> seed{};
    switch(generator.orbit) {
    case Orbit::axis: seed={1.0,0.0,0.0}; break;
    case Orbit::edge: {
        constexpr double s=0.70710678118654752440;
        seed={s,s,0.0};
        break;
    }
    case Orbit::corner: {
        constexpr double s=0.57735026918962576451;
        seed={s,s,s};
        break;
    }
    case Orbit::doubled: {
        const double l=generator.first;
        seed={l,l,std::sqrt(std::max(0.0,1.0-2.0*l*l))};
        break;
    }
    case Orbit::planar: {
        const double p=generator.first;
        seed={p,std::sqrt(std::max(0.0,1.0-p*p)),0.0};
        break;
    }
    case Orbit::general: {
        const double r=generator.first,s=generator.second;
        seed={r,s,std::sqrt(std::max(0.0,1.0-r*r-s*s))};
        break;
    }
    }
    const std::size_t firstOfOrbit=into.size();
    std::sort(seed.begin(),seed.end());
    do {
        for(const double signX:{-1.0,1.0})
            for(const double signY:{-1.0,1.0})
                for(const double signZ:{-1.0,1.0}) {
                    const Vec3 node{signX*seed[0],signY*seed[1],
                                    signZ*seed[2]};
                    const bool duplicate=std::any_of(
                        into.begin()+static_cast<std::ptrdiff_t>(firstOfOrbit),
                        into.end(),
                        [&](const SphereQuadraturePoint& existing) {
                            return (existing.direction-node).norm()<1.0e-13;
                        });
                    if(!duplicate)
                        into.push_back(
                            {node,4.0*pi*generator.normalizedWeight});
                }
    } while(std::next_permutation(seed.begin(),seed.end()));
}

inline std::vector<SphereQuadraturePoint> expand(
        std::initializer_list<Generator> generators) {
    std::vector<SphereQuadraturePoint> rule;
    for(const Generator& generator:generators) appendOrbit(rule,generator);
    return rule;
}

// The tabulated rules, built once on first use.  Returning a pointer to the
// function-local static is what keeps the flux loop off the allocator: the
// old sphereQuadrature rebuilt a std::vector on the heap at every call, and
// electromagneticFieldFluxRates calls it once per integration step.
inline const std::vector<SphereQuadraturePoint>* tabulated(int nodeCount) {
    switch(nodeCount) {
    case 26: {
        // Degree 7.  The weights are exact rationals: 1/21, 4/105, 9/280.
        static const std::vector<SphereQuadraturePoint> rule=expand({
            {Orbit::axis,0.0,0.0,0.047619047619047619048},
            {Orbit::edge,0.0,0.0,0.038095238095238095238},
            {Orbit::corner,0.0,0.0,0.032142857142857142857}});
        return &rule;
    }
    case 50: {
        // Degree 11.  Also exact rationals: 4/315, 64/2835, 27/1280,
        // 14641/725760, with the doubled orbit at l=sqrt(1/11).
        static const std::vector<SphereQuadraturePoint> rule=expand({
            {Orbit::axis,0.0,0.0,0.012698412698412698413},
            {Orbit::edge,0.0,0.0,0.022574955908289241623},
            {Orbit::corner,0.0,0.0,0.02109375},
            {Orbit::doubled,0.30151134457776362265,0.0,
             0.020173335537918871252}});
        return &rule;
    }
    case 110: {
        static const std::vector<SphereQuadraturePoint> rule=expand({
            {Orbit::axis,0.0,0.0,0.0038282704950077357},
            {Orbit::corner,0.0,0.0,0.0097937375124621910},
            {Orbit::doubled,0.18511563534526035,0.0,0.0082117372831938609},
            {Orbit::doubled,0.39568947305626379,0.0,0.0095954713360612946},
            {Orbit::doubled,0.69042104838229246,0.0,0.0099428148911751418},
            {Orbit::planar,0.47836902881218757,0.0,0.0096949963616637085}});
        return &rule;
    }
    case 194: {
        static const std::vector<SphereQuadraturePoint> rule=expand({
            {Orbit::axis,0.0,0.0,0.0017823407029467122},
            {Orbit::edge,0.0,0.0,0.0057169059469172594},
            {Orbit::corner,0.0,0.0,0.0055733831531803515},
            {Orbit::doubled,0.44469331839394804,0.0,0.0055187714487162347},
            {Orbit::doubled,0.28924656407574989,0.0,0.0051582376862526888},
            {Orbit::doubled,0.67129734427151899,0.0,0.0056087040757898882},
            {Orbit::doubled,0.12993354723470232,0.0,0.0041067770239256495},
            {Orbit::planar,0.34577021999433044,0.0,0.0050518460683153729},
            {Orbit::general,0.15904171061980538,0.52511857256697980,
             0.0055302489150390339}});
        return &rule;
    }
    case 302: {
        static const std::vector<SphereQuadraturePoint> rule=expand({
            {Orbit::axis,0.0,0.0,0.00085441695053510157},
            {Orbit::corner,0.0,0.0,0.0035992209029283508},
            {Orbit::doubled,0.70117674239278927,0.0,0.0036500744781326022},
            {Orbit::doubled,0.65663372553267974,0.0,0.0036049168619254235},
            {Orbit::doubled,0.47290369278945227,0.0,0.0035768246127088682},
            {Orbit::doubled,0.35155963103670884,0.0,0.0034498282076958005},
            {Orbit::doubled,0.22195835915596318,0.0,0.0031089106132587258},
            {Orbit::doubled,0.09617723338258305,0.0,0.0023518863705014117},
            {Orbit::planar,0.57189344218093918,0.0,0.0036008355448544334},
            {Orbit::planar,0.26440477849926736,0.0,0.0029822269762529285},
            {Orbit::general,0.25100032220726348,0.80007584910134055,
             0.0035715979404095407},
            {Orbit::general,0.12335110110856259,0.41276604486477253,
             0.0033923112909537076}});
        return &rule;
    }
    default: return nullptr;
    }
}

}  // namespace sphere_quadrature

// A view of the rule, not a copy of it.  The tabulated orders cost nothing
// per call; the Fibonacci fallback is cached per thread and rebuilt only
// when the requested count changes, so the hot loop never allocates either.
inline std::span<const SphereQuadraturePoint> sphereQuadratureView(
        int directionCount) {
    if(const std::vector<SphereQuadraturePoint>* rule=
           sphere_quadrature::tabulated(directionCount))
        return *rule;
    thread_local std::vector<SphereQuadraturePoint> fallback;
    thread_local int cachedCount=0;
    if(directionCount<1) return {};
    if(directionCount!=cachedCount) {
        constexpr double goldenAngle=pi*(3.0-2.2360679774997896964);
        const double weight=4.0*pi/directionCount;
        fallback.clear();
        fallback.reserve(static_cast<std::size_t>(directionCount));
        for(int index=0;index<directionCount;++index) {
            const double z=1.0-2.0*(index+0.5)/directionCount;
            const double transverse=std::sqrt(std::max(0.0,1.0-z*z));
            const double azimuth=goldenAngle*index;
            fallback.push_back({{transverse*std::cos(azimuth),
                                 transverse*std::sin(azimuth),z},weight});
        }
        cachedCount=directionCount;
    }
    return fallback;
}

// Owning form, for the callers that keep a rule alive across other work.
inline std::vector<SphereQuadraturePoint> sphereQuadrature(int directionCount) {
    const std::span<const SphereQuadraturePoint> view=
        sphereQuadratureView(directionCount);
    return {view.begin(),view.end()};
}

// The degree through which a rule is exact, re-derived from its own nodes.
// The integral of x^i y^j z^k over the unit sphere, divided by 4 pi, is
// zero when any exponent is odd and (i-1)!!(j-1)!!(k-1)!!/(i+j+k+1)!!
// otherwise, so the rule's claimed degree is a measurement and not a label.
// This is the guard on the tables above: the parameters were solved
// numerically, and nothing here depends on them having been typed correctly.
inline int sphereQuadratureExactDegree(
        std::span<const SphereQuadraturePoint> rule,
        double tolerance=2.0e-14,int limit=32) {
    if(rule.empty()) return -1;
    const auto doubleFactorial=[](int n) {
        double value=1.0;
        for(int factor=n;factor>1;factor-=2) value*=factor;
        return value;
    };
    const auto exactMoment=[&](int i,int j,int k) {
        if(i%2||j%2||k%2) return 0.0;
        return doubleFactorial(i-1)*doubleFactorial(j-1)*doubleFactorial(k-1)
              /doubleFactorial(i+j+k+1);
    };
    double totalWeight=0.0;
    for(const SphereQuadraturePoint& point:rule)
        totalWeight+=point.solidAngleWeight;
    if(!(std::abs(totalWeight-4.0*pi)<=tolerance*4.0*pi)) return -1;
    int degree=0;
    for(int order=1;order<=limit;++order) {
        for(int i=0;i<=order;++i) for(int j=0;i+j<=order;++j) {
            const int k=order-i-j;
            double sum=0.0;
            for(const SphereQuadraturePoint& point:rule)
                sum+=point.solidAngleWeight
                    *std::pow(point.direction.x,i)
                    *std::pow(point.direction.y,j)
                    *std::pow(point.direction.z,k);
            if(std::abs(sum/(4.0*pi)-exactMoment(i,j,k))>tolerance)
                return degree;
        }
        degree=order;
    }
    return degree;
}

// Defined below historicalDipoleKinematics/historicalElectricDipoleKinematics
// (which they share with retardedMagneticDipoleField/
// retardedElectricDipoleField), forward-declared here so
// electromagneticFieldFluxRates can fold their contribution into the SAME
// Poynting/stress integral as the charge field -- see that function's own
// comment on why this is not just an addition.
ElectromagneticField farZoneMagneticDipoleField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool sourceIsFirst,bool radiationFieldOnly=false);
ElectromagneticField farZoneElectricDipoleField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool sourceIsFirst,bool radiationFieldOnly=false);

inline ElectromagneticField farZoneChargeField(
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
    // Same guard every near-field solver in this file already carries: a
    // control sphere drawn through the source leaves the direction undefined
    // and the division returns NaN, which the flux integral then sums into
    // the whole Poynting/stress result rather than losing one direction.  It
    // is not reachable at the control radii anything here uses (1e4-1e6 a0
    // against a source at r*), but the radius is a caller's parameter and
    // this was the one member of the family without the check.
    if(!(distance>std::numeric_limits<double>::min())) return {};
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

inline FieldFluxRates electromagneticFieldFluxRates(
    const State& state, const StateHistory& history,
    FarFieldSampling sampling={}) {
    if(sampling.directionCount<1||!(sampling.controlRadius>0.0)
        ||!std::isfinite(sampling.controlRadius)) return {};
    const std::span<const SphereQuadraturePoint> quadrature=
        sphereQuadratureView(sampling.directionCount);
    const Vec3 centre=(state.firstPosition+state.secondPosition)*0.5;
    const double sourceExtent=std::max(
        (state.firstPosition-centre).norm(),
        (state.secondPosition-centre).norm());
    // The observation event is shifted back by the source radius so every
    // direction samples only the available causal history. Directional
    // retardation across the pair then retains E3, M2, toroidal terms and all
    // their interference without assigning convention-dependent pieces.
    const double wavefrontTime=state.time-sourceExtent/c;
    // Velocity of the pair as a whole.  It is exactly zero whenever the centre
    // of mass is at rest, which is every configuration this model integrates
    // in production, so the Doppler factor below is then identically one and
    // nothing downstream moves.
    const Vec3 centreVelocity=(state.firstVelocity*firstMass
        +state.secondVelocity*secondMass)/(firstMass+secondMass);
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
        // Both particles' magnetic dipole fields, folded into the SAME
        // point-by-point field sum the charge fields go into, before the
        // Poynting vector and Maxwell stress are built from it -- not added
        // as separate power/momentum numbers afterward.  E1-M1 (and every
        // other pair) interference is exactly what the cross terms of
        // (E1+Edipole)x(B1+Bdipole) and the corresponding stress tensor
        // carry; summing pre-computed E1 and M1 flux totals instead, as
        // particleMultipoleRadiation used to, discards it identically,
        // since it can shift directionality (and hence momentum/recoil)
        // while contributing nothing at leading order to the total power.
        //
        // No separate farZoneElectricDipoleField call: firstElectricDipole/
        // secondElectricDipole are not an independent moment, they are
        // exactly the motional electric dipole lorentzBoostDipole produces
        // from boosting this SAME rest-frame magnetic moment (see
        // synchronizeCovariantDipoles).  farZoneMagneticDipoleField's exact
        // relativistic construction already carries that motional channel
        // -- it is a single boosted current loop's field, not two
        // independent point-multipole sources -- so folding in the electric
        // solver here would double-count it.  It stayed load-bearing back
        // when the magnetic channel was only the low-velocity
        // approximation and genuinely lacked that piece; now that it is the
        // exact field, adding it back is exactly the bug this replaces.
        // The standalone electric-dipole solver itself is untouched, for
        // its own tests and for a possible future genuine (non-motional)
        // electric dipole moment.
        const ElectromagneticField firstMagneticDipoleField=
            farZoneMagneticDipoleField(observationPosition,normal,
                wavefrontTime,centre,history,state,true,
                sampling.radiationFieldOnly);
        const ElectromagneticField secondMagneticDipoleField=
            farZoneMagneticDipoleField(observationPosition,normal,
                wavefrontTime,centre,history,state,false,
                sampling.radiationFieldOnly);
        const Vec3 electric = firstField.electric + secondField.electric
            + firstMagneticDipoleField.electric
            + secondMagneticDipoleField.electric;
        const Vec3 magnetic = firstField.magnetic + secondField.magnetic
            + firstMagneticDipoleField.magnetic
            + secondMagneticDipoleField.magnetic;
        if(std::getenv("POSITRONIUM_DEBUG_FIELDS")) {
            static int farDebugSamples=0;
            if(farDebugSamples++<20)
                std::cerr<<"FAR_DEBUG t="<<state.time
                    <<" chargeE="<<(firstField.electric+secondField.electric).norm()
                    <<" mE="<<(firstMagneticDipoleField.electric
                        +secondMagneticDipoleField.electric).norm()
                    <<" totalE="<<electric.norm()<<'\n';
        }
        const Vec3 poynting = cross(electric, magnetic) / mu0;
        const double areaWeight=sampling.controlRadius
            *sampling.controlRadius*point.solidAngleWeight;
        // Doppler factor kappa = 1 - n.beta, converting the RECEIVED rate per
        // unit observer time into the EMITTED rate per unit emitter time.
        // |E_rad|^2 carries kappa^-6, while dP_emit/dOmega carries kappa^-5,
        // and the missing kappa is exactly the difference.
        //
        // Without it the flux from a moving source came out wrong at first
        // order in beta: at beta = 0.35 the radiated power was 1.1955 times
        // the rest-frame value, which the invariance of Larmor's formula
        // forbids, and |p|/(E/c) came out 0.467 where the four-vector
        // structure demands beta = 0.350.  Restoring kappa drops the boost
        // covariance residual from 0.2084 to 6.99e-04, a factor of 298, and
        // brings |p|/(E/c) to 0.35002 against a predicted 0.35003.
        //
        // One kappa is used for both roles, built from the centre-of-mass
        // velocity.  The two particles differ from it by their orbital
        // velocity, so the approximation costs O(beta_orb^2) -- 2.7e-05 for
        // the default pair -- and it is the only choice available for the
        // interference term, which belongs to neither particle alone.
        const double dopplerFactor=1.0-dot(normal,centreVelocity)/c;
        rates.energy += dot(poynting, normal) * areaWeight * dopplerFactor;

        // Outward momentum flux is -T.n for the Maxwell stress convention
        // T_ij=eps0(E_iE_j-E^2 delta_ij/2)+...
        const Vec3 stressOnNormal =
            (electric * dot(electric, normal)
             - normal * (0.5 * electric.squaredNorm())) * epsilon0
          + (magnetic * dot(magnetic, normal)
             - normal * (0.5 * magnetic.squaredNorm())) / mu0;
        const Vec3 momentumFlux = stressOnNormal * (-areaWeight*dopplerFactor);
        rates.momentum += momentumFlux;
        rates.angularMomentum += cross(observationPosition, momentumFlux);
    }
    return rates;
}

struct LocalElectromagneticFields {
    ElectromagneticField atFirst, atSecond;
};

inline double shortRangeFieldWeight(double distance,
    double regularizationRadius=magneticRegularizationRadius,
    double exponent=magneticRegularizationExponent) {
    if(!(distance>0.0)) return 0.0;
    const double ratio = regularizationRadius / distance;
    return 1.0/(1.0+std::pow(ratio,exponent));
}

// dw/dr for shortRangeFieldWeight's w(r)=1/(1+(regularizationRadius/r)^n):
// w'(r)=n w(1-w)/r (the same combination magneticRadialProfile already
// derives inline for w/r^3's own derivative, exposed on its own for
// retardedMagneticDipoleField's induction-term correction below, which needs
// dw/dr directly rather than the profile's d(w/r^3)/dr).
inline double shortRangeFieldWeightDerivative(double distance,
    double regularizationRadius=magneticRegularizationRadius,
    double exponent=magneticRegularizationExponent) {
    if(!(distance>0.0)) return 0.0;
    const double weight=shortRangeFieldWeight(
        distance,regularizationRadius,exponent);
    return exponent*weight*(1.0-weight)/distance;
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

inline State historicalState(const StateHistory& history, const State& present,
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
    // A non-positive span means the two samples carry the same instant, so
    // there is nothing between them to interpolate and the endpoint IS the
    // answer.  The floor this used to divide by, 1.0e-300, did not say that:
    // it turned a zero span into an extrapolation factor of up to infinity,
    // and interpolateVector then scaled the sampled dipoles by ~1e244.  Their
    // sum of squares overflowed, interpolateDipole's rescale divided
    // infinity by infinity, and the far-zone dipole field received a NaN
    // moment.  Measured on a bound para orbit: 120 of 6575645 calls.
    if(newer!=history.end()) {
        const State& older=*std::prev(newer);
        const double span=newer->time-older.time;
        if(!(span>0.0)) return *newer;
        return interpolateState(older,*newer,(time-older.time)/span);
    }
    const State& older=history.back();
    const double span=present.time-older.time;
    if(!(span>0.0)) return present;
    return interpolateState(older,present,(time-older.time)/span);
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
// Retarded derivative stencils scale with the history's GENUINE node
// spacing: the mean gap between its nodes, leaving out the leading one.
//
// The former rule, twice the gap between the last two stored nodes, is not
// a grid property.  appendStateHistory() overwrites the leading node until a
// full spacing has passed, so that gap is the time since the last genuine
// node: it grows from 0 to the spacing inside every node interval, and it
// differs between the adaptive engine's coarse step and its two half steps
// by exactly the step.  Every retarded dipole field therefore depended on the
// integration path, linearly in dt.  Measured on para released radially at
// 3 r* (floor 0.25 r*), same end state evaluated with the coarse and the
// half-step history: the gradient force differed by 1.75e-3 and the dipole
// field by 38% at the failing step, both halving with each halving of dt,
// while the Lienard-Wiechert force (no stencil) differed by exactly 0.  The
// step-doubling velocity error fell only 4x per halving instead of 8x, and at
// zero relative velocity, where it is normalized by a speed of 3.8e3 m/s,
// the run failed at t = 3.8e-24 s at maximum depth.  With the mean genuine
// spacing both paths see the same stencil (only a node push or pop inside a
// step changes it, by about one part in the node count), the difference is
// exactly 0, and that run completes its observation window.
//
// A stencil from the grid the history is laid down on (4 r/c / 128) was tried
// first and rejected: it assumes physical time scales, and on a history whose
// nodes are a second apart it rounds away to nothing.
inline double historyDerivativeStep(const StateHistory& history,
                                    double spacings) {
    double spacing=0.0;
    if(history.size()>=3)
        spacing=(history[history.size()-2].time-history.front().time)
            /static_cast<double>(history.size()-2);
    else if(history.size()==2)
        spacing=history.back().time-history.front().time;
    return std::max(1.0e-24,spacings*spacing);
}

inline double boundedDerivativeStep(const StateHistory& history,
                             double evaluationTime,
                             double requestedStep, int stencilReach) {
    if(history.empty()||stencilReach<1) return 0.0;
    const double span=evaluationTime-history.front().time;
    if(!(span>0.0)||!std::isfinite(span)) return 0.0;
    const double usable=std::min(requestedStep,
        span/static_cast<double>(stencilReach));
    // Arbitrarily shrinking h next to the oldest retained node turns a
    // derivative into subtraction of nearly identical interpolants divided
    // by an arbitrarily small number.  Once less than one millionth of the
    // requested stencil fits, the history contains no numerically resolvable
    // derivative at that event; report the documented zero fallback.
    if(!(usable>0.0)||usable<1.0e-6*requestedStep) return 0.0;
    return usable;
}

// Position, velocity and the ONE lab dipole a retarded source sample needs.
// historicalState() rebuilds an entire State for this, and interpolateState()
// runs interpolateDipole() on all four dipole vectors -- three of which are
// discarded here.  With five stencil samples per call that was the single
// largest cost in the profile (interpolateDipole alone: 25% of runtime, three
// square roots per invocation).  The interpolation performed on the fields
// that ARE used is identical to interpolateState's, so results are unchanged.
struct RetardedSourceSample { Vec3 position, velocity, moment; };

inline RetardedSourceSample historicalSource(const StateHistory& history,
                                      const State& present,
                                      bool sourceIsFirst, double time) {
    // THE MOMENT HERE IS THE LABORATORY ONE, state.firstDipole, and NOT
    // state.firstProperDipole.  That is deliberate: everything downstream of
    // this sampler -- the retarded dipole fields and the coherent M1 power --
    // is written in the laboratory frame, so the moment entering them must be
    // too.
    //
    // It is worth a warning because the two are nearly impossible to tell
    // apart by inspection and very easy to tell apart by consequence.  At an
    // orbital beta of 2.3e-02 they agree in MAGNITUDE to 1.000039, so a check
    // that compares |mu| will pass whichever one you picked; but the boost
    // relating them turns with the velocity, so its own derivative enters at
    // the same order as the precession's, and their first and second time
    // derivatives stand at 0.8666 and 0.7036 of each other.  Audit 192 built
    // a whole comparison against the proper moment, found the channel
    // "wrong" by a factor of two and three orders at one end, and published
    // it; audit 193 withdrew the lot.  Differentiate the same member the
    // channel does, or the disagreement you measure is your own.
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
    // Same degeneracy as in historicalState above: a non-positive span means
    // the two samples share an instant, and dividing by a 1.0e-300 floor
    // turned that into an extrapolation factor of up to infinity instead of
    // saying so.  The endpoint is the answer.
    const double span=target.time-older.time;
    if(!(span>0.0)) return pick(target);
    const double fraction=(time-older.time)/span;
    const RetardedSourceSample a=pick(older);
    const RetardedSourceSample b=pick(target);
    return {interpolateVector(a.position,b.position,fraction),
            interpolateVector(a.velocity,b.velocity,fraction),
            interpolateDipole(a.moment,b.moment,fraction)};
}

inline RetardedDipoleKinematics historicalDipoleKinematics(
    const StateHistory& history, const State& present, bool sourceIsFirst,
    double time) {
    double derivativeStep=historyDerivativeStep(history,2.0);
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
        // Second order, matching first and third beside it.  This used to
        // read (m - 2m(-h) + m(-2h))/h^2, whose Taylor expansion is
        // f'' - h f''' + O(h^2) -- FIRST order, alone in this block, while
        // the first and third derivatives here are both second order.  On a
        // moment whose dominant motion is the local field swinging round at
        // the ORBITAL frequency, f''' ~ omega f'', so that -h f''' term is a
        // relative error of order omega*h, and it does not cancel between
        // the two particles: audit 192 measured the resulting M1 power at
        // half to two thirds of the converged value through the middle of
        // the opening-angle scan and three to four ORDERS too large at the
        // parallel end, where the true power comes from a near-exact
        // cancellation the first-order error swamps.
        second=(moment(middle)*2.0-moment(before)*5.0
               +moment(twiceBefore)*4.0-moment(threeBefore))
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
        // The mirror of the backward branch above, and second order for the
        // same reason.
        second=(moment(middle)*2.0-moment(after)*5.0
               +moment(twiceAfter)*4.0-moment(threeAfter))
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

inline RetardedElectricDipoleKinematics historicalElectricDipoleKinematics(
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double time) {
    double derivativeStep=historyDerivativeStep(history,2.0);
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

// Moment entering a localized source integral, as distinct from the spatial
// component of the covariant dipole tensor stored in State.  A constant-time
// volume element contracts by gamma, hence P_source=p_lab/gamma and
// M_source=mu_lab/gamma (Sautbekov's moving-moment transformation).  Taking
// the derivatives after this conversion is essential for accelerated motion:
// simply dividing p_dot or mu_dot by gamma would omit derivatives of gamma.
inline RetardedElectricDipoleKinematics historicalIntegratedDipoleKinematics(
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double time,bool electricMoment) {
    double derivativeStep=historyDerivativeStep(history,2.0);
    derivativeStep=boundedDerivativeStep(history,time,derivativeStep,2);
    const auto sample=[&](double sampleTime) {
        const State state=historicalState(history,present,sampleTime);
        const Vec3 velocity=sourceIsFirst
            ?state.firstVelocity:state.secondVelocity;
        const Vec3 laboratoryMoment=electricMoment
            ?(sourceIsFirst?state.firstElectricDipole
                           :state.secondElectricDipole)
            :(sourceIsFirst?state.firstDipole:state.secondDipole);
        return laboratoryMoment/gamma(velocity);
    };
    const State middle=historicalState(history,present,time);
    const Vec3 moment=sample(time);
    if(!(derivativeStep>0.0)) {
        return {sourceIsFirst?middle.firstPosition:middle.secondPosition,
                sourceIsFirst?middle.firstVelocity:middle.secondVelocity,
                moment,{},{}};
    }
    Vec3 first,second;
    if(derivativeStep>0.0) {
        // The central and backward stencils meet at time + h = present, where
        // their second derivatives differ at O(h).  A hard switch there made
        // every retarded dipole field jump in time whenever its source sat
        // h c from the observer: once h stops shrinking with r (the history
        // floor of appendStateHistory), that is ~0.25 r* for e+e-, and one
        // pole of one gradient probe crossing it moved m'' by 3x and the
        // moment force by 50% within 1e-40 s (audit section 107).  Between
        // reach 1 and 2 the two are blended with a C1 weight; outside that
        // band, and next to the history front, nothing changes.  The sibling
        // stencils (historicalDipoleKinematics, the electric one) are not on
        // the force path measured there and keep the hard switch.
        const double reach=(present.time-time)/derivativeStep;
        const bool nearFront=time-derivativeStep<history.front().time;
        if(reach<1.0||(reach<2.0&&!nearFront)) {
            const Vec3 before=sample(time-derivativeStep);
            const Vec3 twiceBefore=sample(time-2.0*derivativeStep);
            first=(moment*3.0-before*4.0+twiceBefore)
                /(2.0*derivativeStep);
            second=(moment-before*2.0+twiceBefore)
                /(derivativeStep*derivativeStep);
            if(reach>1.0) {
                const Vec3 after=sample(time+derivativeStep);
                const double x=reach-1.0;
                const double weight=x*x*(3.0-2.0*x);
                first=first*(1.0-weight)
                    +((after-before)/(2.0*derivativeStep))*weight;
                second=second*(1.0-weight)
                    +((after-moment*2.0+before)
                      /(derivativeStep*derivativeStep))*weight;
            }
        } else if(nearFront) {
            const Vec3 after=sample(time+derivativeStep);
            const Vec3 twiceAfter=sample(time+2.0*derivativeStep);
            first=(after*4.0-moment*3.0-twiceAfter)
                /(2.0*derivativeStep);
            second=(twiceAfter-after*2.0+moment)
                /(derivativeStep*derivativeStep);
        } else {
            const Vec3 before=sample(time-derivativeStep);
            const Vec3 after=sample(time+derivativeStep);
            first=(after-before)/(2.0*derivativeStep);
            second=(after-moment*2.0+before)
                /(derivativeStep*derivativeStep);
        }
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

inline ElectricQuadrupole electricQuadrupole(const State& state) {
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

inline ElectricQuadrupole electricQuadrupoleThirdDerivative(
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
inline Vec3 electricDipoleMoment(const State& state) {
    return state.firstPosition*firstCharge+state.secondPosition*secondCharge;
}

inline Vec3 electricDipoleThirdDerivativeAtStep(
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
inline double electricDipoleDerivativeStep(const StateHistory& history,
                                    const State& present) {
    double derivativeStep=1.0e-24;
    if(history.size()>=2) derivativeStep=std::max(derivativeStep,
        2.0*(history.back().time-history[history.size()-2].time));
    return boundedDerivativeStep(history,present.time,derivativeStep,8);
}

inline Vec3 electricDipoleThirdDerivative(const State& state,
                                   const StateHistory& history) {
    return electricDipoleThirdDerivativeAtStep(
        state,history,electricDipoleDerivativeStep(history,state));
}

inline MutualForces coherentElectricDipoleReaction(
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
    individualLandauLifshitzSelfOnly,
    disabled,
    // Same E1 dipole power already computed below as
    // leadingElectricDipolePower, but instead of removing it from the orbit
    // as a continuous drag force, it is banked as Poissonian HAZARD
    // (integral of power/hbar*omega_orb) and paid out in discrete,
    // momentum-conserving kicks of size hbar*omega_orb -- see
    // applyStochasticDipolePhoton() in crem_trajectory.hpp, the only place
    // that reads this tag.  Deliberately carries ZERO continuous reaction
    // force here (same as disabled): between photons the pair moves on the
    // bare mutual Lorentz-force trajectory, exactly conserving energy, the
    // same way a real emitter does not lose energy continuously between
    // quantum jumps.
    //
    // EXPOSURE GAP, found and closed.  The hazard accumulator above only
    // sees a single measured orbit per checkpoint; crem_collapse.hpp's
    // secular estimator analytically SKIPS up to 200000 more between
    // checkpoints, where nearly all of the real hazard lives (measured at
    // a=3.1pm: ~5.5e-5 per single orbit vs ~11 photons of expected hazard
    // across a typical 200000-orbit skip).  A second, independent hazard
    // accumulator in crem_collapse.hpp now integrates the SAME
    // power/(hbar*omega) rate analytically across the skipped span itself
    // (closed form against the classical u(n)=u0(1-Jx)^(-2/3) envelope the
    // skip's own jump already assumes; verified against brute-force
    // quadrature to 1e-12 relative -- see README), firing photons directly
    // into the osculating elements (energy kick plus the same closed-form
    // k=-(1-e^2)/(2+e^2) ratio the continuum models' own angular-momentum
    // update already uses) since no mechanical State exists during a
    // skipped span.
    //
    // This did more than close the gap -- it surfaced a genuine physical
    // effect.  hbar*omega_orb at positronium's OWN scale is comparable to
    // or larger than its total binding energy throughout the collapse
    // (~9-13 eV per photon near the start, against ~6.8 eV of binding), so
    // individual photons are not a small perturbation on the classical
    // curve the way the (u0(1-Jx)^-2/3) derivation's own small-hazard
    // regime assumes -- they are large, rare, discrete jumps, and the
    // PATH between them is exactly conserved (see this reaction model's
    // top comment).  Measured directly (seed 42): the classical/continuum
    // models give ~36-40 ps; this one gives 665 ps for the same seed and a
    // median 149 ps / mean 276+-91 ps (sigma/mean~1.0, i.e. genuinely
    // wide, not a tight distribution around the continuum answer) over 10
    // seeds -- a real, order-of-magnitude difference from discretizing
    // emission this way, not a rounding artefact.
    //
    // This is the production model, and the difference above is the reason
    // rather than a caveat against it.  Unlike --zpf's SED probe it shows a
    // qualitative, physically interpretable departure from the continuum:
    // photon-sized discretization of a classical orbit does not reduce
    // smoothly to the classical answer once the photon scale stops being
    // small against the system's own energy scale.  At the pair Bohr radius
    // hbar*omega is 13.6 eV against 6.8 eV of binding, so that scale
    // separation is absent from the outset and the continuum answer is the
    // approximation, not this one.
    //
    // The direction this serves is a DETERMINISTIC determination of the
    // quantum parameters: the photon energy is fixed by the orbit
    // (correspondence, or the level spacing where a ladder exists -- see
    // crem_collapse.hpp), and --emission deterministic fixes its timing too.
    // What remains stochastic here is the emission threshold, which is the
    // conventional description of spontaneous emission and stays the
    // default; it is not the part the model is committed to.
    //
    // DIRECTION, added.  Every kick above was originally isotropic in
    // effect (only |E| and |L| moved; OsculatingElements itself has no
    // orientation).  Real E1 radiation from an orbiting pair is a
    // ROTATING, not linearly oscillating, dipole, whose angular pattern is
    // dP/dOmega proportional to (1+cos^2(theta)) from the orbital angular
    // momentum axis (maximal along the axis, not the sin^2(theta)
    // in-plane-maximal pattern a single linear dipole would give) -- see
    // crem_collapse.hpp's photon-firing loop for the closed-form (Cardano)
    // sampling of theta and the orbit-averaged treatment of the resulting
    // plane tilt (the true anomaly at firing time is unknowable in this
    // elements-only representation, so the tilt's azimuth is drawn
    // uniformly at random instead of guessed -- an isotropic random walk
    // of the orbital plane, not a systematic precession, which is the
    // correct phase-averaged picture).  Not tied to spin: positronium's
    // quantum spin (S=0 para, S=1 ortho) already governs dipole
    // alignment/annihilation selection rules elsewhere in this model, and
    // is a different axis from the classical orbital angular momentum
    // this emission pattern actually depends on.
    //
    // LINEAR MOMENTUM, fixed.  Every kick above balanced energy exactly but
    // let the photon's own momentum hbar*omega/c vanish -- unlike the
    // continuous models, whose self-force is Newton's-third-law consistent
    // with the field it radiates by construction.  Measured directly (not
    // asserted): this omission is NOT the "many orders of magnitude
    // negligible" simplification the mechanical-loop version of this kick
    // (crem_trajectory.hpp) still describes -- the ratio p_photon/p_orbital
    // equals exactly (reduced Compton wavelength of the pair)/a, ~0.007 at
    // the starting radius but already ~0.25 by a=3.1pm and >1 below
    // ~0.77pm, i.e. genuinely significant for most of the depth this model
    // reaches, well inside comptonBarrierRadius.  Fixed in
    // crem_collapse.hpp: the photon's now-sampled full 3D direction gives a
    // real recoil, applied as a uniform velocity shift to BOTH particles
    // (centreOfMassVelocity) -- correct because CREM's bound initial
    // conditions are prepared at EXACTLY zero total momentum and every
    // continuous model keeps that true, so momentum conservation reduces
    // to "give the whole system a common kick", not a differential one.
    // The resulting CM kinetic-energy change (v_cm.p_photon +
    // p_photon^2/(2M), computed exactly rather than dropped as
    // second-order) is charged to the orbital budget on top of the
    // photon's own energy, so total energy still drops by exactly the
    // photon energy credited to radiatedEnergyTotal.  The mechanical-loop
    // kick in crem_trajectory.hpp was NOT given the same fix (see its own
    // updated comment) -- it lacks the position/orientation information
    // this fix needs and is rarely exercised in production regardless.
    //
    // ANGULAR MOMENTUM, audited -- and found not fixable the same way.  An
    // earlier version of the block above also tilted the tracked orbital
    // plane per photon, treating the recoil as an r x delta-v kick to the
    // RELATIVE motion (r ~ semi-major axis, delta-v ~ photonEnergy/(c*
    // reducedMass)).  That double-spent the same photon momentum: the
    // linear-momentum fix above already balances it as a uniform kick to
    // BOTH particles, and Sum m_i*r_i=0 about the centre of mass means that
    // kick alone contributes EXACTLY zero torque (checked numerically for
    // e+e- and for 1836:1 and 3:1 mass ratios, residual 1e-17) -- a uniform
    // push through a system's own centre of mass cannot torque it.  The
    // tilt was removed rather than reweighted, because the genuine
    // remaining source of torque is not classical at all: a real photon's
    // angular momentum is dominated by SPIN, +-hbar along its own
    // propagation direction, an exact universal fact for any massless
    // spin-1 boson, not an orbit-averaged estimate.  It is not modelled
    // because CREM's own starting point is the Bohr/SED value L=hbar (see
    // physical_constants.hpp), so |L_photon|/|L_orbital| = hbar/hbar = 1
    // exactly at the very first photon, by construction of the initial
    // condition -- not a measured smallness the way p_photon/p_orbital
    // above was, but the opposite: an O(1) disruption every time, which no
    // classical adiabatic bookkeeping (the k ratio included) is built to
    // absorb.  Left as an open, documented limitation: the orbital-plane
    // direction is carried forward unchanged by photon recoil, which is
    // the self-consistent choice given the linear-momentum fix, not a
    // simplification of a mechanism that was ever actually working.
    //
    // ANGULAR MOMENTUM MAGNITUDE, k determined -- and it was silently
    // wrong for nearly every production photon.  k=-(1-e^2)/(2+e^2) itself
    // is unchanged and still exactly correct here: it is a property of the
    // classical dipole force law's r-dependence, not of whether the
    // resulting loss is booked continuously or in photon-sized lumps (same
    // argument already established above for reaction-model-independence).
    // What was wrong is HOW crem_collapse.hpp applied it to a single
    // photon: L*=(E_after/E_before)^k, the frozen-k power law, is only the
    // first-order/small-jump approximation to the differential relation
    // d(lnL)/d(ln|E|)=k(e) it comes from.  The bulk/deterministic branch
    // uses the same approximation but its jump per checkpoint is capped at
    // jumpParameter<=0.30 (energy ratio <~1.33), where the error is mild;
    // a single stochastic photon can carry a much larger multiple of the
    // current orbital energy in one shot (measured up to ~18.5x elsewhere
    // in this file). Measured directly: at that scale the frozen-k formula
    // is not a little off, it is qualitatively wrong -- checked across a
    // real production run (seed 42, o-Ps, 5 photons), EVERY SINGLE ONE
    // gave a negative post-kick e^2, silently clamped to zero by this
    // file's own std::max(0.0, ...) guard, erasing real eccentricity
    // information at essentially every photon event and then propagating
    // that error into every subsequent kHere evaluation.  Fixed by solving
    // the differential relation exactly instead of extrapolating it:
    // separating variables in s=1-e^2 and x=ln|E| gives a closed-form
    // first integral, (1-e^2)^3/(e^4 |E|^3) = const, exactly conserved
    // along any trajectory obeying dL/L=k(e) dE/E for this k(e). Verified
    // two ways: algebraic self-consistency with the eccentricity formula
    // it was derived from, and independently against direct RK4
    // integration of the ODE (1e-12 agreement, RK4's own discretization
    // error, not the closed form's). Solved per photon by bisection on the
    // resulting monotonic cubic in crem_collapse.hpp (Cardano's formula
    // would also close it, but bisection sidesteps casus irreducibilis --
    // multiple real roots needing branch selection -- and photon events
    // are rare enough per trajectory that 80 bisection steps cost
    // nothing). Statistically close to unchanged in aggregate (10 seeds,
    // o-Ps: median 155.7 ps vs the 147.8-151.6 ps this file already
    // recorded pre-fix, 0 numerical failures either way) -- STAMP: level 2
    // (the default at the time), floor on, s_max 0.30, deterministic
    // emission, 10 seeds, budget not recorded; the two figures are a
    // BEFORE/AFTER pair from one scan and are comparable with each other
    // and with nothing else -- even though the
    // per-photon eccentricity value itself moved by O(1) at nearly every
    // event: collapse TIME is set mainly by the energy/hazard integral,
    // which this fix does not touch, while eccentricity feeds the
    // periapsis distance and the period/light-crossing-ratio exit
    // condition -- and both the exact and the old (clamped-to-zero)
    // values were already small (near-circular) at every measured event in
    // this run, so the two differ by less than 0.01 in absolute e^2 even
    // where the old, pre-clamp value was formally unphysical (negative).
    // elements.specificAngularMomentum is carried forward across
    // checkpoints, so this is not a one-off correction -- every subsequent
    // eccentricity, k evaluation and periapsis check downstream of a
    // photon now inherits the corrected value instead of the wrong one,
    // for the rest of that trajectory.
    //
    // ANGULAR MOMENTUM, superseded by photon spin -- the O(1) disruption
    // flagged above as unabsorbable is now implemented anyway, on request,
    // and measured to be survivable.  Both sections above are now
    // superseded: the "direction carried forward unchanged" claim and the
    // k-ratio magnitude update are BOTH replaced in crem_collapse.hpp by a
    // single real vector kick, L_pair -= h*hbar*photonDirection (h=+-1,
    // the photon's helicity), applied to the actual angular-momentum
    // vector (not the specific/magnitude-only representation used
    // elsewhere) and then re-decomposed into magnitude and direction.  Why
    // this replaces rather than adds to the k-ratio result: k comes from
    // the classical reaction TORQUE, and by angular-momentum conservation
    // applied to a classical field, "torque integrated on the orbit" and
    // "what the continuous field carries away" are the same quantity, not
    // two contributions -- so once emission is quantized into a real
    // photon with a real, known spin, applying k on top of it would
    // double-count exactly the way the removed tilt double-counted linear
    // momentum.  h is drawn from the standard Delta-m=+-1 dipole
    // transition's conditional helicity distribution given the emission
    // angle theta already sampled above, P(h=+1|theta)=(1+cos theta)^2 /
    // [2(1+cos^2 theta)], P(h=-1|theta)=(1-cos theta)^2 / [2(1+cos^2
    // theta)] -- not an independent assumption, these two
    // probabilities sum to exactly the (1+cos^2 theta) pattern theta was
    // already drawn from.  What this still cannot capture is the photon's
    // ORBITAL angular momentum relative to the pair (needs the true
    // anomaly this representation does not carry, the same gap that ruled
    // out the tilt): checked directly, the axial expectation <h cos
    // theta> over the full angular distribution is 1/2, not 1, so spin
    // alone recovers only half of the Delta-m=1 selection rule on average.
    // Measured, not just argued, to be survivable: a production batch (o-Ps,
    // 80 trajectories across two seeds) gave 1 numerical failure against 0
    // failures but 3 wall-clock censorings for the same batch size under
    // the pre-spin (k-ratio-applied) code.  Root-caused rather than
    // guessed at (an earlier version of this comment guessed "leaves the
    // bound regime outright", i.e. specificEnergy>=0 -- checked directly
    // and that guess was wrong): reproduced the exact failing trajectory
    // (crem_collapse.hpp's own per-index seed derivation, splitMix64
    // (masterSeed+index), makes any trajectory in a batch replayable
    // singly) and traced it to crem_collapse.hpp's maxRelativeLossPerOrbit
    // guard, not to specificEnergy going non-negative: a spin kick had
    // just driven the orbit to e^2~0.945, and the SINGLE mechanically-
    // measured orbit that follows (the same one every reaction model
    // relies on for its dE/dt, drawn from the real retarded-field flux,
    // independent of which reaction model is active) then genuinely
    // radiated 56% of the binding energy in that one orbit -- over the
    // 50% guard this file documents elsewhere as "a secular inspiral
    // cannot shed a large fraction of its own binding energy in ONE
    // orbit; if it could, orbit averaging would not apply in the first
    // place".  That guard is doing its job in general, but was checked
    // again rather than left alone on that verdict, and the second look
    // found it was over-firing here specifically: for isStochastic,
    // deltaEnergyPerOrbit (the quantity this guard tests) is used for
    // NOTHING else in this function but the guard itself and the
    // Larmor-ratio diagnostic just below it (separately finite-guarded,
    // never fed back into the state) -- expectedLossPerOrbit, the
    // Larmor-orbit-averaged rate at the CURRENT osculating elements, sizes
    // orbitsToSkip instead (see its own comment).  So the failure mode
    // this guard exists to prevent -- a corrupted single-orbit measurement
    // extrapolated across many skipped orbits -- cannot happen through
    // this path for this model regardless of how large deltaEnergyPerOrbit
    // reads.  Nor is the measurement itself suspect the way the guard's
    // own motivating case was: that one came from the coherent model's
    // third-derivative-of-dipole-moment reaction-force stencil degenerating,
    // which stochasticElectricDipole never evaluates (no continuous
    // reaction force is applied for it at all) -- orbitalRadiatedEnergy
    // instead comes from the flux integral
    // (electromagneticFieldFluxRates via particleMultipoleRadiation),
    // structurally unrelated to that failure mode.  Fixed in
    // crem_collapse.hpp: the guard's magnitude half (not its isfinite half,
    // which still applies unconditionally to every model) is gated off for
    // isStochastic.  Checked, not assumed: the exact trajectory that
    // failed above (seed 107) now completes (172.78 ps, matching the scale
    // already on record) with zero guard trips, and a re-run of the same
    // 80-trajectory two-seed batch plus two more batches (100 and 50
    // further trajectories) gives 0 numerical failures across all 230,
    // where the unfixed guard gave 1/80.  Collapse time itself barely
    // moved either way (ortho median 100-155 ps across several seeds,
    // versus the 147.8-151.6 ps already on record from before either
    // (same stamp as above: level 2, floor on, s_max 0.30, deterministic,
    // 10 seeds -- a before/after pair, not a citable lifetime)
    // change), because collapse time is set by the energy/hazard integral,
    // which none of this touches; eccentricity -- now genuinely reaching
    // values like e^2=0.9-0.95 that the old k-ratio path could never
    // produce -- affects which exit condition a trajectory hits and how
    // soon, and, it turned out, whether an unrelated guard mistakes a
    // real result for a corrupted one.
    //
    // ENERGY, a genuine leak found while auditing multipole completeness
    // (asked: is E1 alone, or E1+M2, the right description of this
    // system's radiation?  Checked directly: for equal-mass, opposite-
    // charge pairs the orbital magnetic dipole and electric quadrupole
    // moments vanish EXACTLY -- Sum q_i(m_j/M)^2(r x v) and Sum q_i
    // (3rr-r^2) both collapse to a factor of q1+q2 or q1-q2 that is zero
    // for this mass ratio, checked numerically -- so M2/E3, not the
    // usually-assumed M1/E2, are the genuine next multipole order after
    // E1 here.  Estimated at beta^4 relative to E1 (one order beyond the
    // usual beta^2 for M1/E2, since M2 pays both the "magnetic" and the
    // "one l higher" penalty): ~3e-9 at a_Ps -- genuinely negligible, not
    // fixed, and the search that ruled it out surfaced something that was
    // not negligible.)
    //
    // RE-AUDITED.  The structure holds exactly and is worth stating in
    // closed form: for equal masses the multipole of order n goes as
    // (r/2)^n [q1 + q2(-1)^n], which for q2 = -q1 vanishes for EVEN n and
    // survives for odd; the magnetic series picks up one further sign
    // reversal from v2 = -v1 and so has the OPPOSITE parity.  Hence E1
    // survives, E2 vanishes, E3 survives; M1 (orbital) vanishes, M2
    // survives.  The next non-vanishing order after E1 is {E3, M2}, both
    // beta^4 in power, exactly as claimed.
    //
    // The MAGNITUDE claim was wrong.  This comment used to add that beta^4
    // grows "to order unity right at the Compton barrier".  It does not.
    // A circular orbit at r* has |E| = K/2r* = (2/g) alpha m c^2 = 1/2 mu
    // v^2, giving beta^2 = 4 alpha and
    //
    //     beta^4 = 8.52e-04 at the barrier, not 1 -- wrong by 1174x,
    //     beta^4 = 3.08e-06 at the period/light-crossing threshold of 150,
    //     beta^4 = 2.84e-09 at a_Ps (the ~3e-9 above, which was right).
    //
    // The error is in the SAFE direction: the neglected order is three
    // orders smaller than the comment claimed, so the decision not to
    // implement it is better founded than its own justification said.
    //
    // The ENDPOINT claim that followed here is out of date.  It said
    // production stops at the ground state, at beta^4 = 2.8e-9.  That holds
    // only under --ground-state-floor, which is off by default.  The default
    // run stops on the period/light-crossing limit of 150 (CollapseStopCause::
    // RetardationLimit); trace 42/2 ends on a circular orbit at a = 4.9 r*,
    // i.e. beta^2 = 4 alpha/4.9 = 5.9e-3 and beta^4 = 3.5e-5.  Still
    // negligible for the lifetime, which is set at large radii.
    //
    // MEASURED (audit section 78): full far-zone flux against the secular
    // formula's Coulomb-only Larmor power on the same engine trajectory,
    // circular orbits, no moments.  flux/Larmor - 1 = -3.2e-5, -3.2e-4,
    // -1.09e-3, -3.45e-3, -7.44e-3 at 548, 55, 16, 5.5 and 2.7 r*, i.e.
    // -0.60 beta^2 with a beta^4 coefficient near -9 at the deepest point.
    // So what the analytic ladder misses is dominated by the relativistic
    // correction to E1 itself, not by E3/M2.  Larger still, and not a
    // multipole at all: the flux follows the ACTUAL relative acceleration,
    // which includes the dipole-dipole force.  With moments along L that
    // shifts the E1 power by -1.2e-3 (para) at 55 r* and by -16% (para) /
    // +46% (ortho) at 5.5 r*, which Coulomb-only Larmor cannot see.  The
    // one mechanically-measured orbit every checkpoint performs
    // (runMechanicalTrajectory, feeding deltaEnergyPerOrbit/
    // orbitalRadiatedEnergy above) genuinely radiates via the real
    // retarded-field flux -- true regardless of reaction model -- but for
    // isStochastic that measured loss was going nowhere: elements.
    // specificEnergy there is touched only inside the photon while-loop,
    // whose hazard integral is scoped explicitly to the orbitsToSkip
    // orbits AFTER this one, never to this one's own already-measured
    // loss.  Measured directly, not assumed small the way the M2/E3 gap
    // above was: on a shallow trajectory (seed 42) the total discarded
    // loss was ~2e-5 of the energy actually credited via photons over the
    // whole run, but on one that spent time at high eccentricity (seed
    // 107, the same trajectory the guard fix above investigated) it
    // reached 38.5% -- a real violation of energy conservation, not a
    // rounding error, and one the guard fix directly above makes easier
    // to reach in practice (a trajectory that previously failed outright
    // at high eccentricity can now run on and accumulate this leak).
    // Fixed in crem_collapse.hpp: run.finalState.orbitalRadiatedEnergy is
    // now credited to elements.specificEnergy and radiatedEnergyTotal
    // once per checkpoint for isStochastic, additively with (not instead
    // of) the photon credits below it -- they cover disjoint spans (this
    // one measured orbit vs. the orbitsToSkip that follow), so this is
    // closing a bookkeeping gap, not double-counting the way earlier
    // fixes in this file had to guard against.  For the deterministic
    // branch this was never a gap: its own energy update is built
    // directly from this same measurement (lossPerOrbit=
    // |deltaEnergyPerOrbit| seeds energyGrowth), so crediting it again
    // here is correctly restricted to isStochastic only.  Verified:
    // clean compile, positronium_validation 33/33, the same seed-107
    // trajectory investigated above still completes at the same 172.78 ps
    // once given enough wall-clock budget (crediting real energy costs
    // real compute: more checkpoints are needed once eccentricity is
    // high, since more of the energy budget is now correctly spent per
    // checkpoint), and a re-run of the seed-99/30-trajectory batch is
    // unchanged (29/30 reach the boundary, 0 numerical failures, median
    // 100.653 ps, identical to before this fix).
    stochasticElectricDipole
};

inline double electricQuadrupoleRadiatedPower(
    const State& state,const StateHistory& history) {
    const ElectricQuadrupole third=
        electricQuadrupoleThirdDerivative(state,history);
    // SI coefficient for Q_ij=sum q(3 x_i x_j-r^2 delta_ij).  Gaussian-unit
    // reference (Landau & Lifshitz, Classical Theory of Fields) is
    // P=D_ij_dddot^2/(180 c^5); converting to SI charges the same way the
    // neighbouring E1/M1 terms in this file do (q_Gauss^2 -> q_SI^2/4*pi*
    // epsilon0, e.g. electricCoefficient=1/(4*pi*epsilon0) in
    // retardedElectricDipoleField) gives 1/(720*pi*epsilon0*c^5), not
    // 1/(180*pi*epsilon0*c^5) -- audit found this coefficient short a
    // factor of 4 (missing the 4*pi every other Coulomb-type constant here
    // carries).  Confirmed inert in production: this power only feeds the
    // saturated dominanceGate in particleMultipoleRadiation() (automatic
    // reaction-model blending), never the energy/momentum bookkeeping, and
    // the measured dominance ratios (~3.5e13 for e+e-, ~3.3e3 for p+e-) are
    // orders of magnitude clear of the gate's [10,20] threshold either way.
    return third.squaredNorm()/(720.0*pi*epsilon0*c*c*c*c*c);
}

MutualForces individualLandauLifshitzSelfForces(
    const State& state, const MutualForces& external,
    const StateHistory& history,bool includeMutual=true);

inline DipoleRadiationReaction coherentMagneticDipoleRadiationReaction(
    const RetardedDipoleKinematics& first,
    const RetardedDipoleKinematics& second) {
    constexpr double coefficient=mu0/(6.0*pi*c*c*c);
    // In the long-wavelength M1 approximation the pair is one localized
    // source with m=m1+m2.  Radiation amplitudes add before they are squared;
    // summing the two individual powers would discard the interference term
    // 2*m1''*m2'' and would radiate even when the two derivatives cancel.
    const Vec3 totalFirstDerivative=
        first.firstDerivative+second.firstDerivative;
    const Vec3 totalSecondDerivative=
        first.secondDerivative+second.secondDerivative;
    const Vec3 totalThirdDerivative=
        first.thirdDerivative+second.thirdDerivative;
    DipoleRadiationReaction result;
    result.power=coefficient*totalSecondDerivative.squaredNorm();

    // Both moments see the same Abraham-Lorentz reaction field generated by
    // the total source.  This is the torque counterpart of the power cross
    // term and makes exact destructive interference torque-free as well.
    result.firstTorque=
        cross(first.moment,totalThirdDerivative)*coefficient;
    result.secondTorque=
        cross(second.moment,totalThirdDerivative)*coefficient;
    const Vec3 intrinsicAngularMomentumRate=cross(
        totalFirstDerivative,totalSecondDerivative)*coefficient;

    // A localized M1 pattern has no net recoil in the pair COM frame.  Boost
    // its coherent power with the velocity c^2 P_pair/E_pair instead of
    // assigning separate momenta to powers which no longer exist separately.
    const double firstEnergyOverC2=gamma(first.velocity)*firstMass;
    const double secondEnergyOverC2=gamma(second.velocity)*secondMass;
    const Vec3 totalMomentum=first.velocity*firstEnergyOverC2
                           +second.velocity*secondEnergyOverC2;
    const double totalEnergyOverC2=firstEnergyOverC2+secondEnergyOverC2;
    const double totalEnergy=totalEnergyOverC2*c*c;
    result.momentumRate=totalMomentum*(result.power/totalEnergy);
    const Vec3 centreOfEnergyPosition=
        (first.position*firstEnergyOverC2
         +second.position*secondEnergyOverC2)/totalEnergyOverC2;
    // The intrinsic far-zone flux differs instantaneously from minus the
    // reaction torque by the magnetic-dipole Schott term.  The second term is
    // the orbital angular momentum of the boosted radiation about the global
    // origin; it vanishes in the pair COM frame used by production runs.
    result.angularMomentumRate=intrinsicAngularMomentumRate
        +cross(centreOfEnergyPosition,result.momentumRate);
    return result;
}

inline DipoleRadiationReaction dipoleRadiationReaction(
    const State& state,const StateHistory& history) {
    const RetardedDipoleKinematics first=historicalDipoleKinematics(
        history,state,true,state.time);
    const RetardedDipoleKinematics second=historicalDipoleKinematics(
        history,state,false,state.time);
    return coherentMagneticDipoleRadiationReaction(first,second);
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
    FieldFluxRates outwardFlux, magneticDipoleFlux;
    double leadingElectricDipolePower = 0.0;
    double electricQuadrupolePower = 0.0;
    double landauLifshitzValidity = 0.0;
    double coherentDerivativeConsistency = 0.0;
    double sourceCompactness = 0.0;
    double coherentWeight = 0.0;
    bool coherentSelected = false;
};

inline ParticleMultipoleRadiation particleMultipoleRadiation(
    const State& state, const MutualForces& externalForces,
    const StateHistory& history,bool computeOutwardFlux=true,
    ChargeRadiationReactionModel reactionModel=
        ChargeRadiationReactionModel::individualLandauLifshitz,
    bool mutualRadiationAlreadyRetarded=false) {
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
    // The E2 power is needed by two unrelated consumers: the automatic
    // model's dominance gate, and the quantized mode's hazard, which sums
    // every channel before dividing by hbar*omega.  It stays off for the
    // continuous models, where nothing reads it and the history stencils it
    // costs would be thrown away.
    //
    // For the default pair this is free in the only sense that matters: the
    // pair quadrupole is Q = kappa (3 d d - d^2 delta) with
    // kappa = (q1 m2^2 + q2 m1^2)/(m1+m2)^2, which is EXACTLY zero for any
    // mass-symmetric neutral pair -- e+e-, mu+mu- and p+pbar alike, and
    // exactly in floating point, the two masses being the same literal.  The
    // channel is real only for an asymmetric pair, where the centre of mass
    // and the centre of charge differ: kappa = -0.9989 e for p+e-, i.e.
    // essentially full strength, and until now that pair radiated no E2 at
    // all in any mode.
    const bool needsQuadrupolePower=needsBlendingGates
        ||reactionModel
            ==ChargeRadiationReactionModel::stochasticElectricDipole;
    // disabled and stochasticElectricDipole never read the Landau-Lifshitz
    // force or its validity ratio (the ratio feeds only the automatic
    // model's gate), yet it costs two more retarded force sums per call,
    // three times the whole step, and a NaN in it rejected the step as
    // non-finite (audit section 107).  landauLifshitzValidity stays 0 there.
    const bool needsLandauLifshitz=
        reactionModel!=ChargeRadiationReactionModel::disabled
        &&reactionModel
            !=ChargeRadiationReactionModel::stochasticElectricDipole;
    const MutualForces ll=needsLandauLifshitz
        ?individualLandauLifshitzSelfForces(
            state,externalForces,history,
            !mutualRadiationAlreadyRetarded
                &&reactionModel
                    !=ChargeRadiationReactionModel::individualLandauLifshitzSelfOnly)
        :MutualForces{};
    if(needsLandauLifshitz)
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
    if(needsQuadrupolePower) {
        result.electricQuadrupolePower =
            electricQuadrupoleRadiatedPower(state,history);
    }

    const DipoleRadiationReaction magnetic =
        dipoleRadiationReaction(state, history);
    result.firstDipoleTorque = magnetic.firstTorque;
    result.secondDipoleTorque = magnetic.secondTorque;
    result.magneticDipoleFlux={magnetic.power,magnetic.momentumRate,
                               magnetic.angularMomentumRate};
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
        const double nonElectricPower=result.magneticDipoleFlux.energy
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
    // stochasticElectricDipole carries zero continuous force here, same as
    // disabled: its whole point is that nothing drags the orbit between
    // photons.  leadingElectricDipolePower above is still computed
    // unconditionally, which is all applyStochasticDipolePhotons() needs.
    if(reactionModel!=ChargeRadiationReactionModel::disabled
       &&reactionModel!=ChargeRadiationReactionModel::stochasticElectricDipole) {
        // retardedExternalForces already carries the partner's acceleration
        // field.  Adding the mutual part of LL or q_i*d''' on top counts that
        // radiative exchange twice.  In that configuration every continuous
        // model therefore contributes only the particle self terms.  The full
        // coherent reaction remains correct and available with the
        // non-retarded Coulomb-Darwin external force.
        if(mutualRadiationAlreadyRetarded||!needsCoherentReaction)
            result.chargeReaction=ll;
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
    // No longer added here: electromagneticFieldFluxRates now folds both
    // particles' dipole fields into its own Poynting/stress integral (see
    // its comment), so outwardFlux already carries magneticDipoleFlux's
    // contribution -- with the E1-M1 (and dipole-dipole) interference terms
    // this addition-of-separate-totals could never produce -- and adding it
    // again here would double-count.  magneticDipoleFlux itself is
    // unchanged and still computed unconditionally: the continuous
    // reaction torque and the quantized channel's hazard both read it
    // directly, independent of computeOutwardFlux.
    return result;
}

inline bool finiteRadiationResponse(const ParticleMultipoleRadiation& response) {
    return std::isfinite(response.outwardFlux.energy)
        && std::isfinite(response.leadingElectricDipolePower)
        && std::isfinite(response.magneticDipoleFlux.energy)
        && std::isfinite(response.electricQuadrupolePower)
        && std::isfinite(response.landauLifshitzValidity)
        && std::isfinite(response.coherentDerivativeConsistency)
        && std::isfinite(response.sourceCompactness)
        && std::isfinite(response.coherentWeight)
        && isFinite(response.chargeReaction.first)
        && isFinite(response.chargeReaction.second)
        && isFinite(response.firstDipoleTorque)
        && isFinite(response.secondDipoleTorque)
        && isFinite(response.magneticDipoleFlux.momentum)
        && isFinite(response.magneticDipoleFlux.angularMomentum)
        && isFinite(response.outwardFlux.momentum)
        && isFinite(response.outwardFlux.angularMomentum);
}

inline void applyDipoleRadiationTorque(State& state,
                                const ParticleMultipoleRadiation& reaction,
                                double dt) {
    const double firstGyromagneticRatio=firstGyromagneticRatioOf();
    const double secondGyromagneticRatio=secondGyromagneticRatioOf();
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

inline bool dipoleRadiationTorqueEnabled(ChargeRadiationReactionModel model) {
    // A stochastic photon supplies the linear recoil, but it does not rotate
    // the magnetic moment which emitted it.  Keep that orientation reaction
    // active in the quantized model; only the fully disabled reference turns
    // every radiation-reaction channel off.
    //
    // CREM_NO_M1_TORQUE suppresses just this channel, leaving the charge
    // sector untouched, so that "is the M1 reaction live in production" is a
    // measurement rather than a reading of this function (audit 136).  It is
    // a diagnostic ablation and is off by default.
    static const bool suppressed=std::getenv("CREM_NO_M1_TORQUE")!=nullptr;
    if(suppressed) return false;
    return model!=ChargeRadiationReactionModel::disabled;
}

inline void applyDipoleRadiationTorqueForModel(
    State& state,const ParticleMultipoleRadiation& reaction,double dt,
    ChargeRadiationReactionModel model) {
    if(dipoleRadiationTorqueEnabled(model))
        applyDipoleRadiationTorque(state,reaction,dt);
}

// f(r)=w(r)/r^3 and its derivatives.  Defining the regulator at the
// vector-potential level ensures that every use of the dipole field employs
// the same B=curl(A), including all derivatives of w.
inline MagneticRadialProfile magneticRadialProfile(double distance,
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

// Low-velocity baseline: finds the retarded position of the source but
// then applies the potential of a dipole essentially at rest in the lab
// frame (missing the kappa=1-n.beta convection factor, direction aberration,
// and the velocity-coupled pieces of the acceleration-order terms that an
// arbitrarily moving/accelerated point dipole's exact field has -- see
// Sautbekov for the closed form).  Kept as the regularized baseline used by
// retardedElectricDipoleField below: production adds the exact relativistic
// point-dipole correction and retains this model only inside the declared
// short-range smoothing core.
inline ElectromagneticField retardedElectricDipoleFieldLowVelocity(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    bool regularized=true) {
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
    // Magnitude only: n stays the true direction, but every 1/distance power
    // below reads the floored value, so the field this returns stops
    // strengthening once the true separation dips under the barrier -- same
    // "freeze at the barrier" rule as the conservative sector.
    const double fieldDistance=std::max(distance,separationFloor());
    const double inverseDistance=1.0/fieldDistance;
    const double weight=regularized?shortRangeFieldWeight(fieldDistance):1.0;
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

// Low-velocity baseline, with the same limitation as
// retardedElectricDipoleFieldLowVelocity above.  It follows directly from
// A=mu0/(4pi)[m(t_r)x n/r^2 + mdot(t_r)x n/(c r)], the potential of a
// dipole essentially at rest in the lab frame, missing kappa, aberration,
// and the velocity-coupled acceleration-order terms.  See
// retardedMagneticDipoleFieldExact below for the relativistic correction.
inline ElectromagneticField retardedMagneticDipoleFieldLowVelocity(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    bool regularized=true) {
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
    // Same "freeze at the barrier" floor as retardedElectricDipoleField;
    // clampedDisplacement additionally feeds the static term below, which
    // (unlike the induction/radiation terms here) reads a vector, not a
    // bare magnitude.
    const Vec3 clampedDisplacement=
        clampedSeparationVector(displacement,separationFloor());
    const double fieldDistance=clampedDisplacement.norm();
    const double inverseDistance=1.0/fieldDistance;
    const double weight=regularized?shortRangeFieldWeight(fieldDistance):1.0;
    constexpr double coefficient=mu0/(4.0*pi);
    const auto transversePattern=[&](const Vec3& moment) {
        return n*(3.0*dot(n,moment))-moment;
    };
    // The static term is exactly curl(A_reg), as in the conservative sector.
    // Retardation adds the induction and radiation pieces of the same
    // potential.  Thus mdot=mddot=0 reproduces regularizedDipoleField bit for
    // bit instead of the former, inconsistent w*B_point approximation.
    //
    // The declared potential is A_reg=w(r)[m(t_r)xn/r^2+mdot(t_r)xn/(cr)], so
    // B=curl(A_reg) must expand as w curl(A)+grad(w) x A (product rule), not
    // just w times the unregularized curl.  Working that expansion through
    // (m(t_r) contributes via -mdot(t_r)/c per unit r from the chain rule
    // through the retarded time itself) shows the radiation term needs no
    // correction -- its w' pieces cancel exactly -- but the induction term
    // is missing exactly this transverse piece.  Confirmed numerically
    // against a finite-difference curl of A_reg directly at
    // r=magneticRegularizationRadius (where dw/dr itself peaks): the two
    // agree to ~1e-10 with the correction and were off by ~37% without it
    // (see regularizedInductionCurlResidual in maxwell_validation.hpp).
    const Vec3 missingInductionCurlTerm=
        (dipole.firstDerivative-n*dot(n,dipole.firstDerivative))
        *((regularized?shortRangeFieldWeightDerivative(fieldDistance):0.0)
            /(c*fieldDistance))
        *coefficient;
    Vec3 magnetic=regularizedDipoleField(clampedDisplacement,dipole.moment,
                    regularized?magneticRegularizationRadius:0.0)
                 +(transversePattern(dipole.firstDerivative)
                      *(inverseDistance*inverseDistance/c)
                  +cross(n,cross(n,dipole.secondDerivative))
                      *(inverseDistance/(c*c)))*(coefficient*weight)
                 +missingInductionCurlTerm;
    Vec3 electric=(cross(n,dipole.firstDerivative)
                      *(inverseDistance*inverseDistance)
                  +cross(n,dipole.secondDerivative)
                      *(inverseDistance/c))*(coefficient*weight);
    return {electric,magnetic};
}

// Far-zone form of the same two-pole limit used by the local field below.
// Its retarded-time equation is written in wavefront form so it never
// subtracts two O(controlRadius/c) numbers.  Each pole retains its own
// kappa, aberrated direction and acceleration field; radiationFieldOnly
// removes only the 1/R^2 velocity field, exactly as farZoneChargeField does.
// The same retention measure for the FAR-ZONE two-pole construction, which
// had none until audit 164.  Kept separate from the near-field global: the
// two constructions run at different pole separations and on different
// geometries, and one overwriting the other's reading would hide exactly the
// case the detector exists for.
inline thread_local double gFarPoleCancellationRatio=0.0;

template<class MomentSampler>
inline ElectromagneticField farZoneTwoChargeLimitDipoleField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool sourceIsFirst,bool radiationFieldOnly,
    double poleSeparationFraction,
    MomentSampler&& momentAt) {
    const double controlRadius=(observationPosition-centre).norm();
    if(!(controlRadius>std::numeric_limits<double>::min())) return {};
    Vec3 referenceMoment,referenceFirst,referenceSecond;
    momentAt(wavefrontTime,referenceMoment,referenceFirst,referenceSecond);
    if(!isFinite(referenceMoment)||!isFinite(referenceFirst)
       ||!isFinite(referenceSecond)) return {};
    // The pole separation must stay a fixed small fraction of the SOURCE's
    // own scale, not of the control radius.  The near-field
    // twoChargeLimitDipoleField above can use the observer distance for
    // this because observer and source sit at comparable range there; a
    // far-zone control sphere is deliberately parked many wavelengths out
    // to isolate the radiation field, and tying eps to that arbitrary
    // radius would grow the two-charge construction's own O(eps^2)
    // multipole truncation error right along with it -- worse at larger,
    // supposedly "farther into the far zone", radii, which is exactly
    // backwards.  The source's own distance from the pair's centre is the
    // right scale: it is what actually limits how large a pole separation
    // can stay while still standing in for a point dipole.
    const ChargeKinematics wavefrontCharge=historicalCharge(
        history,present,sourceIsFirst,wavefrontTime);
    const double sourceScale=std::max(
        (wavefrontCharge.position-centre).norm(),separationFloor());
    const double referenceTime=sourceScale/c;
    const double momentScale=std::max({referenceMoment.norm(),
        referenceFirst.norm()*referenceTime,
        referenceSecond.norm()*referenceTime*referenceTime});
    if(!(poleSeparationFraction>0.0)
       ||!std::isfinite(poleSeparationFraction)) return {};
    const double poleCharge=momentScale
        /(poleSeparationFraction*sourceScale);
    if(!(poleCharge>0.0)||!std::isfinite(poleCharge)) return {};
    const auto poleKinematics=[&](double sign,double time,
            Vec3& position,Vec3& velocity,Vec3& acceleration) {
        const ChargeKinematics charge=historicalCharge(
            history,present,sourceIsFirst,time);
        Vec3 moment,first,second;
        momentAt(time,moment,first,second);
        const double inversePole=sign/(2.0*poleCharge);
        position=charge.position+moment*inversePole;
        velocity=charge.velocity+first*inversePole;
        acceleration=charge.acceleration+second*inversePole;
    };
    bool polesValid=true;
    const auto poleField=[&](double sign) {
        double emissionTime=wavefrontTime;
        Vec3 position,velocity,acceleration;
        for(int iteration=0;iteration<8;++iteration) {
            poleKinematics(sign,emissionTime,position,velocity,acceleration);
            const double refined=wavefrontTime
                +dot(normal,position-centre)/c;
            if(std::abs(refined-emissionTime)<=1.0e-30
               +1.0e-14*std::abs(emissionTime)) {
                emissionTime=refined;
                break;
            }
            emissionTime=refined;
        }
        poleKinematics(sign,emissionTime,position,velocity,acceleration);
        if(!isFinite(position)||!isFinite(velocity)||!isFinite(acceleration)) {
            polesValid=false;
            return ElectromagneticField{};
        }
        const Vec3 displacement=observationPosition-position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) {
            polesValid=false;
            return ElectromagneticField{};
        }
        const Vec3 direction=displacement/distance;
        const Vec3 beta=velocity/c;
        const double betaSquared=beta.squaredNorm();
        const double kappa=1.0-dot(direction,beta);
        // The two poles are an auxiliary representation of a physical point
        // dipole.  A reconstructed pole that is actually superluminal is not
        // a valid Lienard-Wiechert source, but a valid source at 0.999c or
        // faster must retain its velocity: clipping it to 0.999c changes the
        // physical field by orders of magnitude at high gamma.
        if(!(betaSquared<1.0) || !(kappa>0.0) || !std::isfinite(kappa)) {
            polesValid=false;
            return ElectromagneticField{};
        }
        const Vec3 velocityField=(direction-beta)
            *((1.0-betaSquared)
                /(kappa*kappa*kappa*distance*distance));
        const Vec3 accelerationField=cross(direction,
            cross(direction-beta,acceleration))
            /(c*c*kappa*kappa*kappa*distance);
        const Vec3 electric=(radiationFieldOnly?accelerationField
            :velocityField+accelerationField)*(coulomb*sign*poleCharge);
        if(!isFinite(electric)) {
            polesValid=false;
            return ElectromagneticField{};
        }
        return ElectromagneticField{electric,cross(direction,electric)/c};
    };
    const ElectromagneticField positive=poleField(1.0);
    const ElectromagneticField negative=poleField(-1.0);
    if(!polesValid) { gFarPoleCancellationRatio=0.0; return {}; }
    const ElectromagneticField total{positive.electric+negative.electric,
                                     positive.magnetic+negative.magnetic};
    // HOW MUCH OF A SINGLE POLE SURVIVED.  By construction the subtraction
    // should retain a fraction poleSeparationFraction of each pole, so this
    // ratio is ~1 when the construction behaves and large when it does not.
    // Both components are measured because the two wrappers read different
    // ones: the magnetic dipole field is built from `magnetic`, the electric
    // one from `electric`, and a failure in either is a failure.
    const double electricPole=std::max(positive.electric.norm(),
                                       negative.electric.norm());
    const double magneticPole=std::max(positive.magnetic.norm(),
                                       negative.magnetic.norm());
    const double electricRatio=electricPole>0.0
        ?(total.electric.norm()/electricPole)/poleSeparationFraction:0.0;
    const double magneticRatio=magneticPole>0.0
        ?(total.magnetic.norm()/magneticPole)/poleSeparationFraction:0.0;
    gFarPoleCancellationRatio=std::max(electricRatio,magneticRatio);
    return total;
}

// Detect and retreat, the far-zone twin of the one
// covariantDipoleGradientForce and dipoleCouplingMaterialRate already apply
// to the near-field construction
// (audit 164).  Without it the far-zone dipole field returned a value up to
// forty times oversized on isolated directions -- about five in 200000 at one
// measured geometry -- which the energy integral barely feels, being a sum of
// positive terms, and which scatters the MOMENTUM integral, whose net is only
// about 2% of the same scale and so is dominated by any outlier.
//
// The retreat is two decades in one step and then verifies itself, for the
// reasons written out at covariantDipoleGradientForce: the error is not
// monotone in the pole separation, so refining gradually can land somewhere
// worse, and if the fallback is also anomalous the original value is kept
// rather than a silently different field being substituted.
// THE LIMIT IS NOT THE NEAR FIELD'S 30.  That constant is set against a
// healthy near-field ratio of ~1; the far-zone construction's healthy ratio
// is an order of magnitude smaller -- over 40000 directions at one measured
// geometry, median 0.068, p99 0.38, p99.9 1.13 -- so 30 never fires and the
// detector sat dead when it was first written this way.
//
// Set from the two distributions instead.  Scanning 200000 directions found
// five where the production separation returns a field 43 to 225 times the
// converged one, and their ratios are 4.74, 6.06, 19.8, 24.9 and 24.9.  The
// healthy and broken ranges do OVERLAP -- a correct field was seen at 5.4 --
// so this cannot be a clean classifier, and it does not have to be: where
// the retreat is unnecessary it is also harmless, because the two
// separations then agree to about 1e-06.  A limit of 1.0 therefore catches
// every broken direction observed while retreating on 0.12% of the healthy
// ones at a cost of one extra evaluation each.
inline constexpr double farPoleCancellationLimit=1.0;
inline constexpr double farPoleRetreatFraction=1.0e-7;

inline ElectromagneticField farZoneMagneticDipoleField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool sourceIsFirst,bool radiationFieldOnly) {
    const auto evaluate=[&](double fraction) {
        return farZoneTwoChargeLimitDipoleField(
            observationPosition,normal,wavefrontTime,centre,history,present,
            sourceIsFirst,radiationFieldOnly,fraction,
            [&](double time,Vec3& moment,Vec3& first,Vec3& second) {
                const RetardedElectricDipoleKinematics dipole=
                    historicalIntegratedDipoleKinematics(
                        history,present,sourceIsFirst,time,false);
                moment=dipole.moment/(c*c);
                first=dipole.firstDerivative/(c*c);
                second=dipole.secondDerivative/(c*c);
            });
    };
    ElectromagneticField dual=evaluate(1.0e-5);
    if(gFarPoleCancellationRatio>farPoleCancellationLimit) {
        const ElectromagneticField retreated=
            evaluate(farPoleRetreatFraction);
        if(gFarPoleCancellationRatio<=farPoleCancellationLimit)
            dual=retreated;
    }
    return {dual.magnetic*(-c*c),dual.electric};
}

inline ElectromagneticField farZoneElectricDipoleField(
    const Vec3& observationPosition,const Vec3& normal,double wavefrontTime,
    const Vec3& centre,const StateHistory& history,const State& present,
    bool sourceIsFirst,bool radiationFieldOnly) {
    const auto evaluate=[&](double fraction) {
        return farZoneTwoChargeLimitDipoleField(
            observationPosition,normal,wavefrontTime,centre,history,present,
            sourceIsFirst,radiationFieldOnly,fraction,
            [&](double time,Vec3& moment,Vec3& first,Vec3& second) {
                const RetardedElectricDipoleKinematics dipole=
                    historicalIntegratedDipoleKinematics(
                        history,present,sourceIsFirst,time,true);
                moment=dipole.moment;
                first=dipole.firstDerivative;
                second=dipole.secondDerivative;
            });
    };
    const ElectromagneticField direct=evaluate(1.0e-5);
    if(gFarPoleCancellationRatio>farPoleCancellationLimit) {
        const ElectromagneticField retreated=
            evaluate(farPoleRetreatFraction);
        if(gFarPoleCancellationRatio<=farPoleCancellationLimit)
            return retreated;
    }
    return direct;
}

// Exact retarded field of a point dipole moving and accelerating
// arbitrarily, built as the eps->0 limit of two point charges +-q separated
// by moment(t)/q and each moving exactly as that separation dictates
// (translation, rotation and internal stretching of the moment all fall out
// automatically). Differentiating the ALREADY EXACT Lienard-Wiechert
// point-charge field this way -- rather than hand-deriving the moving
// dipole's kappa/aberration/acceleration-order terms symbolically, the
// Sautbekov route -- reuses machinery this codebase already trusts
// (identical velocityField/accelerationField/kappa structure to
// lienardWiechertField) instead of introducing a second, independent
// relativistic formula that would need its own from-scratch verification.
// momentAt(time, moment, firstDerivative, secondDerivative) supplies the
// dipole moment and its first two time derivatives; retardedElectricDipoleField
// below feeds it the electric moment directly, retardedMagneticDipoleField
// feeds it moment/c^2 and applies the vacuum duality substitution
// (E,B)->(cB,-E/c) afterwards -- verified against this file's existing
// low-velocity formulas, whose retardation terms match under exactly that
// substitution (electricCoefficient=mu0/4pi*c^2, so p->m/c^2 carries E_p's
// coefficient onto B_m's exactly, and B_p(m/c^2)=-E_m/c^2 the same way).
//
// eps is set from poleSeparationFraction*referenceDistance rather than a
// fixed length.  The two poles are symmetric, so the first omitted multipole
// is O(eps^2), while cancellation contributes O(machine-epsilon/eps).
// Balancing those errors calls for the cube root of machine epsilon, not its
// square root: 1e-5 retains about eleven useful relative digits in the
// uniform-motion probes and is substantially less noisy inside the outer
// spatial gradient used by covariantDipoleGradientForce.
// CREM_DEBUG_GRAD reads this: the last few twoChargeLimitDipoleField calls,
// so an outlier probe of covariantDipoleGradientForce's stencil can be traced
// to the field evaluation that produced it.  Diagnostic only.
struct TwoChargeLimitTrace {
    int centralIterations=0;
    double centralResidual=0.0, poleCharge=0.0, momentScale=0.0;
    // The construction is a DIFFERENCE of two large point-charge fields, so
    // how much of each survives the subtraction is the number that says
    // whether its output can be trusted at a given geometry.
    double poleMagnitude=0.0, sumMagnitude=0.0;
    bool earlyReturn=false;
};
inline thread_local TwoChargeLimitTrace gTwoChargeTrace[8];
inline thread_local int gTwoChargeTraceCount=0;
// Cancellation quality of the most recent twoChargeLimitDipoleField call,
// as |sum|/|pole| divided by the pole separation fraction it was built with.
// The construction is designed so the subtraction retains a fraction eps/d of
// each pole, so this ratio is ~1 whenever the geometry is benign and large
// when it is not.  Measured: 1.09 on healthy probes, 197 on the degenerate
// one.  Callers that can afford to re-evaluate use it to notice that the
// value they just received is not trustworthy.
inline thread_local double gPoleCancellationRatio=0.0;
template<class MomentSampler>
inline ElectromagneticField twoChargeLimitDipoleField(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double poleSeparationFraction,double softening,MomentSampler&& momentAt) {
    const ChargeKinematics nowCharge=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    double centralRetardedTime=observationTime
        -(observationPosition-nowCharge.position).norm()/c;
    // Tracing is diagnostic and this is a hot path -- several thread_local
    // writes per call, one of them inside the Newton loop -- so it is gated.
    // Leaving it ungated cost a measurable slowdown: a para batch that
    // completed 12 of 12 at a 45 s budget dropped to 6 of 12.
    static const bool traceTwoCharge=
        std::getenv("CREM_DEBUG_GRAD")!=nullptr;
    static thread_local TwoChargeLimitTrace discardedTrace;
    TwoChargeLimitTrace& trace=traceTwoCharge
        ?gTwoChargeTrace[gTwoChargeTraceCount++ & 7]:discardedTrace;
    if(traceTwoCharge) trace=TwoChargeLimitTrace{};
    for(int iteration=0;iteration<16;++iteration) {
        if(traceTwoCharge) trace.centralIterations=iteration+1;
        const ChargeKinematics charge=historicalCharge(
            history,present,sourceIsFirst,centralRetardedTime);
        const Vec3 displacement=observationPosition-charge.position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) { trace.earlyReturn=true; return {}; }
        const Vec3 direction=displacement/distance;
        const double residual=centralRetardedTime+distance/c-observationTime;
        if(traceTwoCharge) trace.centralResidual=residual;
        const double derivative=std::max(
            1.0e-8,1.0-dot(direction,charge.velocity/c));
        const double refined=centralRetardedTime-residual/derivative;
        if(std::abs(refined-centralRetardedTime)<=1.0e-30
           +1.0e-14*std::abs(centralRetardedTime)) {
            centralRetardedTime=refined;
            break;
        }
        centralRetardedTime=refined;
    }
    const ChargeKinematics referenceCharge=historicalCharge(
        history,present,sourceIsFirst,centralRetardedTime);
    Vec3 referenceMoment,referenceFirst,referenceSecond;
    momentAt(centralRetardedTime,referenceMoment,
        referenceFirst,referenceSecond);
    if(!isFinite(referenceMoment)||!isFinite(referenceFirst)
       ||!isFinite(referenceSecond)) { trace.earlyReturn=true; return {}; }
    if(referenceMoment.squaredNorm()==0.0
       &&referenceFirst.squaredNorm()==0.0
       &&referenceSecond.squaredNorm()==0.0) { trace.earlyReturn=true; return {}; }
    const double rawReferenceDistance=
        (observationPosition-referenceCharge.position).norm();
    const double referenceDistance=std::max(rawReferenceDistance,
        softening);
    // One short-range regularization for BOTH poles, set by the dipole's
    // own retarded distance.  Clamping each pole's distance separately
    // (max(distance, cutoff, floor) per pole) broke the pair whenever the
    // retarded distance sat within one pole separation of the floor: one
    // pole was clamped and the other was not, and their difference -- the
    // dipole field itself -- picked up the unbalanced 1/d^2 of a single
    // pole, amplified by 1/poleSeparationFraction.  Measured at the tilted
    // para failure at 0.33 r* (floor 0.25 r*): retarded distance 1.46e-20 m
    // above the floor, pole distances spread by 4.1e-20 m; the field ramped
    // by 2e-4 of itself within 4e-31 s (slope 5e26 1/s against ~5e21 1/s
    // physical) and sat 3% off its converged value, while pole separations
    // 1e-6 and 1e-7 agreed with each other to 1e-6.  The adaptive integrator
    // failed on it at maximum depth.  Scaling both pole distances by the
    // same factor keeps the pair intact on both sides and is continuous at
    // the floor, where the factor is exactly 1.
    const double poleDistanceScale=
        rawReferenceDistance>std::numeric_limits<double>::min()
            ?std::max(1.0,std::max(nuclearCutoff,softening)
                           /rawReferenceDistance)
            :1.0;
    if(!(poleSeparationFraction>0.0)
       ||!std::isfinite(poleSeparationFraction)) { trace.earlyReturn=true; return {}; }
    const double poleSeparation=poleSeparationFraction*referenceDistance;
    const double referenceTime=referenceDistance/c;
    const double momentScale=std::max({referenceMoment.norm(),
        referenceFirst.norm()*referenceTime,
        referenceSecond.norm()*referenceTime*referenceTime});
    const double poleCharge=momentScale/poleSeparation;
    if(traceTwoCharge) { trace.poleCharge=poleCharge;
                         trace.momentScale=momentScale; }
    if(!(poleCharge>0.0)||!std::isfinite(poleCharge)) {
        trace.earlyReturn=true; return {}; }
    // Both poles read the CHARGE from one cubic segment -- the one holding the
    // central retarded time -- and expand it by their own small time offset
    // (at most separation/c).  Reading it at each pole's own time instead
    // let the two poles straddle a history node, where the C1 Hermite history
    // makes the acceleration jump, and that jump is multiplied by a pole
    // charge growing as 1/separation: see historicalChargeWithJerk.  Inside
    // the segment this is identical to the old per-pole read, because the
    // cubic is its own third-order Taylor series.  The moment derivatives are
    // still sampled at each pole's own time; they were measured continuous
    // through the same node.
    const ChargeKinematicsWithJerk centralCharge=historicalChargeWithJerk(
        history,present,sourceIsFirst,centralRetardedTime);
    const auto poleKinematics=[&](double sign,double time,
            Vec3& position,Vec3& velocity,Vec3& acceleration) {
        const double offset=time-centralRetardedTime;
        const Vec3 jerkStep=centralCharge.jerk*offset;
        position=centralCharge.position+centralCharge.velocity*offset
            +centralCharge.acceleration*(0.5*offset*offset)
            +jerkStep*(offset*offset/6.0);
        velocity=centralCharge.velocity+centralCharge.acceleration*offset
            +jerkStep*(0.5*offset);
        acceleration=centralCharge.acceleration+jerkStep;
        Vec3 moment,first,second;
        momentAt(time,moment,first,second);
        const double inversePole=sign/(2.0*poleCharge);
        position+=moment*inversePole;
        velocity+=first*inversePole;
        acceleration+=second*inversePole;
    };
    bool polesValid=true;
    const auto poleField=[&](double sign) -> ElectromagneticField {
        Vec3 position,velocity,acceleration;
        double retardedTime=centralRetardedTime;
        for(int iteration=0;iteration<16;++iteration) {
            poleKinematics(sign,retardedTime,position,velocity,acceleration);
            const Vec3 displacement=observationPosition-position;
            const double distance=displacement.norm();
            if(!(distance>std::numeric_limits<double>::min())) {
                polesValid=false;
                trace.earlyReturn=true;
                return {};
            }
            const Vec3 direction=displacement/distance;
            const double residual=retardedTime+distance/c-observationTime;
            const double derivative=std::max(1.0e-8,
                1.0-dot(direction,velocity/c));
            const double refined=retardedTime-residual/derivative;
            if(std::abs(refined-retardedTime)<=1.0e-30
                +1.0e-14*std::abs(retardedTime)) {
                retardedTime=refined; break;
            }
            retardedTime=refined;
        }
        poleKinematics(sign,retardedTime,position,velocity,acceleration);
        const Vec3 displacement=observationPosition-position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) {
            polesValid=false;
            trace.earlyReturn=true;
            return {};
        }
        const Vec3 direction=displacement/distance;
        const Vec3 beta=velocity/c;
        // The pole velocity contains a small derivative correction.  If an
        // ill-conditioned reconstruction makes that auxiliary pole
        // superluminal, reject the construction instead of altering a valid
        // physical source velocity.  In particular, 0.9995c and 0.9999c are
        // valid inputs and must be evaluated at their actual Lorentz factor.
        const double betaSquared=beta.squaredNorm();
        const double kappa=1.0-dot(direction,beta);
        if(!(betaSquared<1.0) || !(kappa>0.0) || !std::isfinite(kappa)) {
            polesValid=false;
            trace.earlyReturn=true;
            return {};
        }
        // With a floor each pole carries lienardWiechertField's covariant
        // Plummer softening, so the pair of poles is the Plummer dipole of
        // pairDipoleField in the static limit (section 66).  The per-pole
        // scaling is smooth, so the straddle of section 56 cannot recur.
        double fieldDistance=distance*poleDistanceScale;
        double plummerScale=1.0;
        if(const double floor=softening; floor>0.0) {
            fieldDistance=distance;
            const double restDistance=distance*kappa
                /std::sqrt(1.0-betaSquared);
            const double ratio=restDistance
                /std::sqrt(restDistance*restDistance+floor*floor);
            plummerScale=ratio*ratio*ratio;
        }
        const Vec3 velocityField=(direction-beta)*((1.0-betaSquared)
            /(kappa*kappa*kappa*fieldDistance*fieldDistance));
        const Vec3 accelerationField=cross(direction,
            cross(direction-beta,acceleration))
            /(c*c*kappa*kappa*kappa*fieldDistance);
        const Vec3 electric=(velocityField+accelerationField)
            *(coulomb*sign*poleCharge*plummerScale);
        if(!isFinite(electric)) {
            polesValid=false;
            trace.earlyReturn=true;
            return {};
        }
        return {electric,cross(direction,electric)/c};
    };
    const ElectromagneticField positivePole=poleField(1.0);
    const ElectromagneticField negativePole=poleField(-1.0);
    if(!polesValid) return {};
    const ElectromagneticField total{
        positivePole.electric+negativePole.electric,
        positivePole.magnetic+negativePole.magnetic};
    const double poleMagnitude=std::max(positivePole.magnetic.norm(),
                                        negativePole.magnetic.norm());
    const double sumMagnitude=total.magnetic.norm();
    gPoleCancellationRatio=
        (poleMagnitude>0.0&&poleSeparationFraction>0.0)
            ?(sumMagnitude/poleMagnitude)/poleSeparationFraction
            :0.0;
    if(traceTwoCharge) { trace.poleMagnitude=poleMagnitude;
                         trace.sumMagnitude=sumMagnitude; }
    if(!isFinite(total.electric)||!isFinite(total.magnetic)) { trace.earlyReturn=true; return {}; }
    return total;
}

// Electric field of the two-pole dipole construction AND its material
// derivative D/Ds along a target worldline x(s)=x0+v s, t(s)=t0+s.
//
// WHY THIS EXISTS.  hiddenMomentumRateForce needs D(mu x E)/Dt.  Taking it
// by differencing the field divides the two-pole construction's cancellation
// floor, ~1e-10 of the dipole field, by a step of about 1/2000 of the
// field's own time scale, and audit 161 measured the result at 1.4e-04 for
// this sector against 1.1e-08 for the charge field.  Carrying the derivative
// through the same arithmetic instead leaves the floor at the cancellation
// itself.
//
// HOW THE RETARDED TIMES ARE DIFFERENTIATED.  Not by differentiating the
// Newton iteration.  Each retarded time is a root of
// F(tau,s) = tau + |x(s)-X(tau)|/c - t(s), so the implicit function theorem
// gives it in closed form from the CONVERGED root:
//     dtau/ds = (1 - nhat.beta_target) / (1 - nhat.beta_source),
// the numerator because the observation event moves at v, the denominator
// the usual Doppler kappa.  The iteration is run in plain doubles exactly as
// twoChargeLimitDipoleField runs it, so both land on the same root.
//
// THE VALUE IS THE TEST.  Everything below mirrors twoChargeLimitDipoleField
// statement for statement, so `electric` must reproduce that function's
// output.  A divergence means the mirror has drifted and the derivative is
// describing a different expression than production evaluates.
//
// STATUS: NOT WIRED IN.  hiddenMomentumRateForce still uses its stencil.
// The mirror is verified -- it reproduces twoChargeLimitDipoleField to
// 5.4e-11 (electric) and 1.0e-10 (magnetic), which is that construction's
// own pole-cancellation floor -- but the DERIVATIVE disagrees with the
// converged stencil it would replace by 1.5e-03, and the stencil is itself
// good to 1e-04 here.  Shipping it in that state would be a regression on a
// term carrying 3.5% of the total force, so it stays unused until the
// discrepancy is closed (audit section 162).
//
// WHAT HAS ALREADY BEEN EXCLUDED, so it is not re-checked:
//   dtau_central/ds   the implicit-function formula, verified against a
//                     numerically differentiated converged root to 4e-12
//   snap              CREM_DUAL_SNAP gives the jerk a differenced rate
//                     instead of the exact zero one cubic segment implies;
//                     the answer does not move, confirming the segment is a
//                     cubic and its jerk constant
//   the third moment derivative's step   CREM_DUAL_MOMENT_STEP sweeps it;
//                     the answer is flat from 1e-02 to 1e-04 of the
//                     light-crossing time.  Zeroing the term outright
//                     (CREM_DUAL_NO_THIRD) moves the error from 1.5e-03 to
//                     6.5e-03,
//                     so the term is needed and is being computed sanely
//   the explicit worldline term in dtau_pole/ds   derived and included
//                     below; it is real but worth only 3e-06 here
//
// WHAT THE REMAINING 1.5e-03 IS.  The moment sampler is not self-consistent
// as a function of the time it is sampled at: differencing `moment`
// numerically disagrees with the `first` it returns by 2.9%, and
// differencing `first` disagrees with `second` by 2.5%, both converged over
// three decades of step.  This chain rule needs a consistent triple, and so
// does the pole reconstruction below, which adds moment, first and second to
// a pole's position, velocity and acceleration as though they were each
// other's time derivatives.  Whether that is a defect of the sampler or a
// different parameterization it is documenting is the open question, and it
// has to be answered before an analytic derivative of this construction can
// mean anything.
struct DipoleElectricWithRate {
    Vec3 electric, electricRate;
    Vec3 magnetic, magneticRate;
    bool valid=false;
};

template<typename MomentSampler>
inline DipoleElectricWithRate twoChargeLimitDipoleElectricWithRate(
    const Vec3& observationPosition,const Vec3& observationVelocity,
    double observationTime,const StateHistory& history,const State& present,
    bool sourceIsFirst,double poleSeparationFraction,double softening,
    MomentSampler&& momentAt) {
    using positronium::dual::Dual;
    using positronium::dual::DualVec3;
    using positronium::dual::dotDual;
    using positronium::dual::crossDual;
    using positronium::dual::normDual;
    using positronium::dual::maxDual;
    DipoleElectricWithRate result;
    const ChargeKinematics nowCharge=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    double centralRetardedTime=observationTime
        -(observationPosition-nowCharge.position).norm()/c;
    for(int iteration=0;iteration<16;++iteration) {
        const ChargeKinematics charge=historicalCharge(
            history,present,sourceIsFirst,centralRetardedTime);
        const Vec3 displacement=observationPosition-charge.position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) return result;
        const Vec3 direction=displacement/distance;
        const double residual=centralRetardedTime+distance/c-observationTime;
        const double derivative=std::max(
            1.0e-8,1.0-dot(direction,charge.velocity/c));
        const double refined=centralRetardedTime-residual/derivative;
        if(std::abs(refined-centralRetardedTime)<=1.0e-30
           +1.0e-14*std::abs(centralRetardedTime)) {
            centralRetardedTime=refined;
            break;
        }
        centralRetardedTime=refined;
    }
    const ChargeKinematics referenceCharge=historicalCharge(
        history,present,sourceIsFirst,centralRetardedTime);
    Vec3 referenceMoment,referenceFirst,referenceSecond;
    momentAt(centralRetardedTime,referenceMoment,
        referenceFirst,referenceSecond);
    if(!isFinite(referenceMoment)||!isFinite(referenceFirst)
       ||!isFinite(referenceSecond)) return result;
    if(referenceMoment.squaredNorm()==0.0
       &&referenceFirst.squaredNorm()==0.0
       &&referenceSecond.squaredNorm()==0.0) return result;
    const double rawReferenceDistance=
        (observationPosition-referenceCharge.position).norm();
    if(!(rawReferenceDistance>std::numeric_limits<double>::min()))
        return result;
    if(!(poleSeparationFraction>0.0)
       ||!std::isfinite(poleSeparationFraction)) return result;
    const Vec3 targetBeta=observationVelocity/c;

    // dtau_central/ds from the converged central root.
    const Vec3 centralDirection=
        (observationPosition-referenceCharge.position)/rawReferenceDistance;
    const double centralKappa=std::max(1.0e-8,
        1.0-dot(centralDirection,referenceCharge.velocity/c));
    const double centralTimeRate=
        (1.0-dot(centralDirection,targetBeta))/centralKappa;

    const DualVec3 observation{
        Dual{observationPosition.x,observationVelocity.x},
        Dual{observationPosition.y,observationVelocity.y},
        Dual{observationPosition.z,observationVelocity.z}};
    const Dual centralTime{centralRetardedTime,centralTimeRate};

    // The moment sampler gives the moment and two derivatives.  A dual time
    // needs one more: d(second)/dt.  That quantity comes from the moment
    // interpolant, which carries no pole cancellation, so differencing it is
    // safe in a way differencing the FIELD is not -- the whole point of this
    // function.  The step is a fixed small fraction of the light-crossing
    // time of the reference distance, well inside the pinned segment.
    double momentStep=1.0e-3*rawReferenceDistance/c;
    if(const char* o=std::getenv("CREM_DUAL_MOMENT_STEP"))
        momentStep=std::atof(o)*rawReferenceDistance/c;
    // EACH OF THE THREE IS DIFFERENTIATED AS ITSELF, not as the next one's
    // integral.  The sampler's `first` and `second` are not d(moment)/dt and
    // d(first)/dt: they are finite-difference stencils of the moment taken
    // at historyDerivativeStep, a HISTORY-NODE scale, so they carry that
    // stencil's truncation.  Measured against a fine difference of the same
    // sampler, `first` sits 2.9% off d(moment)/dt and `second` 2.5% off
    // d(first)/dt, converged over three decades of step (audit 162f).
    //
    // Chaining them as though they were exact derivatives is what left this
    // function's first version 1.5e-03 away from the stencil it replaces.
    // The pole reconstruction below consumes all three as the field defines
    // them, so the derivative has to differentiate the functions the field
    // actually calls -- which means differencing all three, finely.  That is
    // safe for the same reason it was unsafe for the FIELD: the moment
    // interpolant carries no pole cancellation, so a fine step costs nothing
    // here while it cost everything there.
    const auto momentDualAt=[&](Dual time,DualVec3& moment,
                                DualVec3& first,DualVec3& second) {
        Vec3 m,f,s;
        momentAt(time.value,m,f,s);
        Vec3 mAhead,fAhead,sAhead,mBehind,fBehind,sBehind;
        momentAt(time.value+momentStep,mAhead,fAhead,sAhead);
        momentAt(time.value-momentStep,mBehind,fBehind,sBehind);
        const double inverse=1.0/(2.0*momentStep);
        Vec3 momentRate=(mAhead-mBehind)*inverse;
        Vec3 firstRate=(fAhead-fBehind)*inverse;
        Vec3 secondRate=(sAhead-sBehind)*inverse;
        if(std::getenv("CREM_DUAL_CHAINED_MOMENT")) {
            // The refuted version of audit 162, kept so its 1.5e-03 can be
            // reproduced: moment and first chained onto the next quantity,
            // only `second` differenced as itself.
            momentRate=f; firstRate=s;
        }
        if(std::getenv("CREM_DUAL_NO_THIRD")) secondRate=Vec3{};
        const double rate=time.derivative;
        moment={Dual{m.x,momentRate.x*rate},Dual{m.y,momentRate.y*rate},
                Dual{m.z,momentRate.z*rate}};
        first={Dual{f.x,firstRate.x*rate},Dual{f.y,firstRate.y*rate},
               Dual{f.z,firstRate.z*rate}};
        second={Dual{s.x,secondRate.x*rate},Dual{s.y,secondRate.y*rate},
                Dual{s.z,secondRate.z*rate}};
    };

    DualVec3 refMoment,refFirst,refSecond;
    momentDualAt(centralTime,refMoment,refFirst,refSecond);
    const DualVec3 referencePosition{
        Dual{referenceCharge.position.x,referenceCharge.velocity.x
             *centralTimeRate},
        Dual{referenceCharge.position.y,referenceCharge.velocity.y
             *centralTimeRate},
        Dual{referenceCharge.position.z,referenceCharge.velocity.z
             *centralTimeRate}};
    const Dual rawReference=normDual(observation-referencePosition);
    const Dual referenceDistanceDual=maxDual(rawReference,softening);
    const Dual poleDistanceScaleDual=maxDual(
        Dual{1.0,0.0},
        Dual{std::max(nuclearCutoff,softening),0.0}/rawReference);
    const Dual poleSeparationDual=
        Dual{poleSeparationFraction,0.0}*referenceDistanceDual;
    const Dual referenceTimeDual=referenceDistanceDual/Dual{c,0.0};
    const Dual momentScaleDual=maxDual(maxDual(
        normDual(refMoment),
        normDual(refFirst)*referenceTimeDual),
        normDual(refSecond)*referenceTimeDual*referenceTimeDual);
    if(!(momentScaleDual.value>0.0)) return result;
    const Dual poleChargeDual=momentScaleDual/poleSeparationDual;
    if(!(poleChargeDual.value>0.0)||!std::isfinite(poleChargeDual.value))
        return result;

    const ChargeKinematicsWithJerk centralCharge=historicalChargeWithJerk(
        history,present,sourceIsFirst,centralRetardedTime);
    // Inside one pinned cubic segment the jerk is constant, so its own rate
    // is exactly zero.  That is a property of the interpolant, not an
    // approximation: historicalChargeWithJerk reads one Hermite cubic and a
    // cubic's third derivative does not vary along it.
    const auto asDual=[&](const Vec3& value,const Vec3& rate,double scale) {
        return DualVec3{Dual{value.x,rate.x*scale},Dual{value.y,rate.y*scale},
                        Dual{value.z,rate.z*scale}};
    };
    const DualVec3 centralPosition=
        asDual(centralCharge.position,centralCharge.velocity,centralTimeRate);
    const DualVec3 centralVelocity=asDual(centralCharge.velocity,
        centralCharge.acceleration,centralTimeRate);
    const DualVec3 centralAcceleration=asDual(centralCharge.acceleration,
        centralCharge.jerk,centralTimeRate);
    DualVec3 centralJerk=asDual(centralCharge.jerk,Vec3{},0.0);
    if(const char* snap=std::getenv("CREM_DUAL_SNAP")) {
        // Diagnostic: give the jerk a numerically differenced rate instead of
        // the exact zero a single cubic segment implies.
        const double snapStep=1.0e-3*rawReferenceDistance/c;
        const ChargeKinematicsWithJerk ahead=historicalChargeWithJerk(
            history,present,sourceIsFirst,centralRetardedTime+snapStep);
        const ChargeKinematicsWithJerk behind=historicalChargeWithJerk(
            history,present,sourceIsFirst,centralRetardedTime-snapStep);
        const Vec3 snapRate=(ahead.jerk-behind.jerk)*(1.0/(2.0*snapStep));
        centralJerk=asDual(centralCharge.jerk,snapRate,centralTimeRate);
        (void)snap;
    }

    // Plain-double pole kinematics, byte-identical to the production
    // lambda, used only to converge each pole's retarded time.
    const auto poleKinematics=[&](double sign,double time,
            Vec3& position,Vec3& velocity,Vec3& acceleration) {
        const double offset=time-centralRetardedTime;
        const Vec3 jerkStep=centralCharge.jerk*offset;
        position=centralCharge.position+centralCharge.velocity*offset
            +centralCharge.acceleration*(0.5*offset*offset)
            +jerkStep*(offset*offset/6.0);
        velocity=centralCharge.velocity+centralCharge.acceleration*offset
            +jerkStep*(0.5*offset);
        acceleration=centralCharge.acceleration+jerkStep;
        Vec3 moment,first,second;
        momentAt(time,moment,first,second);
        const double inversePole=sign/(2.0*poleChargeDual.value);
        position+=moment*inversePole;
        velocity+=first*inversePole;
        acceleration+=second*inversePole;
    };

    bool polesValid=true;
    struct PoleField { DualVec3 electric, magnetic; };
    const auto poleRate=[&](double sign) -> PoleField {
        Vec3 position,velocity,acceleration;
        double retardedTime=centralRetardedTime;
        for(int iteration=0;iteration<16;++iteration) {
            poleKinematics(sign,retardedTime,position,velocity,acceleration);
            const Vec3 displacement=observationPosition-position;
            const double distance=displacement.norm();
            if(!(distance>std::numeric_limits<double>::min())) {
                polesValid=false; return {};
            }
            const Vec3 direction=displacement/distance;
            const double residual=retardedTime+distance/c-observationTime;
            const double derivative=std::max(1.0e-8,
                1.0-dot(direction,velocity/c));
            const double refined=retardedTime-residual/derivative;
            if(std::abs(refined-retardedTime)<=1.0e-30
                +1.0e-14*std::abs(retardedTime)) {
                retardedTime=refined; break;
            }
            retardedTime=refined;
        }
        poleKinematics(sign,retardedTime,position,velocity,acceleration);
        const Vec3 displacement=observationPosition-position;
        const double distance=displacement.norm();
        if(!(distance>std::numeric_limits<double>::min())) {
            polesValid=false; return {};
        }
        // dtau_pole/ds is NOT just the usual Doppler ratio here.  The pole's
        // worldline depends on s in its own right, not only through its
        // retarded time: it is reconstructed off the CENTRAL charge, whose
        // own retarded time moves with s, and it is displaced by the moment
        // over a pole charge that moves with s too.  Writing the light cone
        // as F(tau_p,s) = tau_p + |x(s) - X_p(tau_p,s)|/c - t(s) = 0,
        //     dtau_p/ds = [1 - nhat.beta_target + nhat.(dX_p/ds)|_tau /c]
        //                 / (1 - nhat.beta_pole).
        // Dropping the middle term is what a textbook derivation gives for a
        // source that does not move with the observation event, and it is
        // wrong by 1.5e-03 of the rate here -- measured against the stencil
        // this function replaces, which is itself good to 1e-04.
        //
        // The explicit part comes from the reconstruction itself: build it
        // once with the retarded time held FIXED, and the derivative that
        // falls out is dX_p/ds at constant tau_p.
        const Vec3 poleDirection=displacement/distance;
        const double poleKappa=std::max(1.0e-8,
            1.0-dot(poleDirection,velocity/c));
        const DualVec3 poleDirectionDual{Dual{poleDirection.x,0.0},
            Dual{poleDirection.y,0.0},Dual{poleDirection.z,0.0}};
        DualVec3 positionDual,velocityDual,accelerationDual;
        const auto reconstruct=[&](Dual poleTime) {
            const Dual offset=poleTime-centralTime;
            const DualVec3 jerkStep=centralJerk*offset;
            positionDual=centralPosition+centralVelocity*offset
                +centralAcceleration*(offset*offset*Dual{0.5,0.0})
                +jerkStep*(offset*offset/Dual{6.0,0.0});
            velocityDual=centralVelocity+centralAcceleration*offset
                +jerkStep*(offset*Dual{0.5,0.0});
            accelerationDual=centralAcceleration+jerkStep;
            DualVec3 moment,first,second;
            momentDualAt(poleTime,moment,first,second);
            const Dual inversePole=
                Dual{sign,0.0}/(Dual{2.0,0.0}*poleChargeDual);
            positionDual+=moment*inversePole;
            velocityDual+=first*inversePole;
            accelerationDual+=second*inversePole;
        };
        reconstruct(Dual{retardedTime,0.0});
        const Dual explicitShift=dotDual(poleDirectionDual,positionDual);
        const double poleTimeRate=
            (1.0-dot(poleDirection,targetBeta)
             +explicitShift.derivative/c)/poleKappa;
        const Dual poleTime{retardedTime,poleTimeRate};
        reconstruct(poleTime);
        const DualVec3 displacementDual=observation-positionDual;
        const Dual distanceDual=normDual(displacementDual);
        const DualVec3 directionDual=displacementDual/distanceDual;
        const DualVec3 betaDual=velocityDual/Dual{c,0.0};
        const Dual betaSquaredDual=dotDual(betaDual,betaDual);
        const Dual kappaDual=Dual{1.0,0.0}-dotDual(directionDual,betaDual);
        if(!(betaSquaredDual.value<1.0)||!(kappaDual.value>0.0)
           ||!std::isfinite(kappaDual.value)) {
            polesValid=false; return {};
        }
        Dual fieldDistanceDual=distanceDual*poleDistanceScaleDual;
        Dual plummerScaleDual{1.0,0.0};
        if(softening>0.0) {
            fieldDistanceDual=distanceDual;
            const Dual restDistance=distanceDual*kappaDual
                /positronium::dual::sqrtDual(Dual{1.0,0.0}-betaSquaredDual);
            const Dual ratio=restDistance/positronium::dual::sqrtDual(
                restDistance*restDistance
                +Dual{softening*softening,0.0});
            plummerScaleDual=ratio*ratio*ratio;
        }
        const Dual kappaCubed=kappaDual*kappaDual*kappaDual;
        const DualVec3 velocityFieldDual=(directionDual-betaDual)
            *((Dual{1.0,0.0}-betaSquaredDual)
              /(kappaCubed*fieldDistanceDual*fieldDistanceDual));
        const DualVec3 accelerationFieldDual=crossDual(directionDual,
            crossDual(directionDual-betaDual,accelerationDual))
            /(Dual{c*c,0.0}*kappaCubed*fieldDistanceDual);
        const DualVec3 electricDual=(velocityFieldDual+accelerationFieldDual)
            *(Dual{coulomb*sign,0.0}*poleChargeDual*plummerScaleDual);
        return PoleField{electricDual,
            crossDual(directionDual,electricDual)/Dual{c,0.0}};
    };
    const PoleField positivePole=poleRate(1.0);
    const PoleField negativePole=poleRate(-1.0);
    if(!polesValid) return result;
    const DualVec3 totalElectric=positivePole.electric+negativePole.electric;
    const DualVec3 totalMagnetic=positivePole.magnetic+negativePole.magnetic;
    result.electric={totalElectric.x.value,totalElectric.y.value,
                     totalElectric.z.value};
    result.electricRate={totalElectric.x.derivative,
                         totalElectric.y.derivative,
                         totalElectric.z.derivative};
    result.magnetic={totalMagnetic.x.value,totalMagnetic.y.value,
                     totalMagnetic.z.value};
    result.magneticRate={totalMagnetic.x.derivative,
                         totalMagnetic.y.derivative,
                         totalMagnetic.z.derivative};
    if(!isFinite(result.electric)||!isFinite(result.electricRate)
       ||!isFinite(result.magnetic)||!isFinite(result.magneticRate))
        return result;
    result.valid=true;
    return result;
}

// Exact point-dipole field before the model's short-range smoothing.  This
// remains separately named so validation can inspect the regularized
// production wrapper against its unsmoothed limit.
inline bool historicalDipoleSourceIsStatic(const StateHistory& history,
    const State& present,bool sourceIsFirst,bool electricMoment) {
    const Vec3 position=sourceIsFirst
        ?present.firstPosition:present.secondPosition;
    const Vec3 moment=electricMoment
        ?(sourceIsFirst?present.firstElectricDipole
                       :present.secondElectricDipole)
        :(sourceIsFirst?present.firstDipole:present.secondDipole);
    const Vec3 velocity=sourceIsFirst
        ?present.firstVelocity:present.secondVelocity;
    const Vec3 acceleration=sourceIsFirst
        ?present.firstAcceleration:present.secondAcceleration;
    if(velocity.squaredNorm()!=0.0||acceleration.squaredNorm()!=0.0)
        return false;
    return std::ranges::all_of(history,[&](const State& state) {
        const Vec3 historicalPosition=sourceIsFirst
            ?state.firstPosition:state.secondPosition;
        const Vec3 historicalVelocity=sourceIsFirst
            ?state.firstVelocity:state.secondVelocity;
        const Vec3 historicalAcceleration=sourceIsFirst
            ?state.firstAcceleration:state.secondAcceleration;
        const Vec3 historicalMoment=electricMoment
            ?(sourceIsFirst?state.firstElectricDipole
                           :state.secondElectricDipole)
            :(sourceIsFirst?state.firstDipole:state.secondDipole);
        return (historicalPosition-position).squaredNorm()==0.0
            &&historicalVelocity.squaredNorm()==0.0
            &&historicalAcceleration.squaredNorm()==0.0
            &&(historicalMoment-moment).squaredNorm()==0.0;
    });
}

// Plummer-softened dipole field; see pairDipoleField (section 66).
// Softening length of a MAGNETIC moment's field.  A magnetic moment is a
// current loop, and the smallest loop that carries mu with the particle's
// charge at speeds not above c has radius R = 2 mu/(|q| c) = 2 r* for e+e-
// (N1: mu/(e c) = r* exactly).  Its magnetization is a uniform disk of that
// radius.  The model smears moments with a Plummer profile, so the loop is
// represented by the Plummer length with the same contact weight: for the
// r^-1/2 density of an orbit ensemble at small r (section 85) the disk gives
// <r^-1/2> = (4/3) R^-1/2 and Plummer 1.5 B(5/4,5/4) eps^-1/2, hence
//
//     eps = [ (3/4) 1.5 B(5/4,5/4) ]^2 R = 0.48341 R = 0.96682 r*.
//
// It replaces the numerical separation floor (0.05 r*) for the magnetic
// dipole sector only; charge fields keep the floor.  One length for the
// whole moment field, not only its contact part: a magnetization term on a
// different length than the rest of the field would give div B != 0.
// CREM_MAGNETIC_RADIUS_SCALE (units of r*) overrides it; 0.05 reproduces
// the floor-softened moment field of sections 66-84.  Audit section 86.
inline double magneticDipoleRadius() {
    const double floor=separationFloor();
    if(!(floor>0.0)) return 0.0;
    static const double overrideScale=[] {
        const char* text=std::getenv("CREM_MAGNETIC_RADIUS_SCALE");
        const double value=text?std::atof(text):0.0;
        return (std::isfinite(value)&&value>0.0)?value:0.0;
    }();
    if(overrideScale>0.0) return overrideScale*comptonBarrierRadius;
    static const double diskPlummerFraction=[] {
        const double beta=std::tgamma(1.25)*std::tgamma(1.25)
            /std::tgamma(2.5);
        const double factor=0.75*1.5*beta;
        return factor*factor;
    }();
    const double loopRadius=2.0*std::abs(firstMagneticMoment)
        /(std::abs(firstCharge)*c);
    return diskPlummerFraction*loopRadius;
}

inline Vec3 plummerDipoleField(const Vec3& sourceToTarget,
                               const Vec3& sourceDipole,double softening) {
    constexpr double magneticConstant=mu0/(4.0*pi);
    const double rhoSquared=sourceToTarget.squaredNorm()+softening*softening;
    if(!(rhoSquared>std::numeric_limits<double>::min())) return {};
    const double inverseRho=1.0/std::sqrt(rhoSquared);
    const double inverseRhoCubed=inverseRho*inverseRho*inverseRho;
    const double inverseRhoFifth=inverseRhoCubed*inverseRho*inverseRho;
    return (sourceToTarget*(3.0*dot(sourceToTarget,sourceDipole)
                *inverseRhoFifth)
           -sourceDipole*inverseRhoCubed)*magneticConstant;
}

// The field above is the Plummer-softened dipole of two separated POLES --
// right for an electric dipole, whose field integrates to -p/(3 eps0) over
// space.  A magnetic moment is a current loop, and a current loop's field
// integrates to +(2 mu0/3) m: that difference is the Fermi contact term,
// which attracts parallel moments (para) and repels antiparallel ones
// (ortho).  With the poles' field used for the magnetic moment (00f0071 to
// section 84) the contact term had the opposite sign and half the strength
// (audit section 83).
//
// The loop field is the curl of the model's own vector potential
// A = (mu0/4 pi) m x r / rho^3, and it differs from the poles' field by a
// purely LOCAL term:
//
//     B_loop = B_poles + mu0 m n(r),     n(r) = 3 eps^2 / (4 pi rho^5),
//
// n being the normalized Plummer density.  mu0 m n is the field of the
// moment's magnetization itself (B = mu0 (H + M)); it does not propagate.
inline Vec3 plummerMagnetizationField(const Vec3& sourceToTarget,
                                      const Vec3& sourceDipole,
                                      double softening) {
    const double rhoSquared=sourceToTarget.squaredNorm()+softening*softening;
    if(!(rhoSquared>std::numeric_limits<double>::min())) return {};
    const double density=3.0*softening*softening
        /(4.0*pi*std::pow(rhoSquared,2.5));
    return sourceDipole*(mu0*density);
}

inline Vec3 plummerMagneticDipoleField(const Vec3& sourceToTarget,
                                       const Vec3& sourceDipole,
                                       double softening) {
    return plummerDipoleField(sourceToTarget,sourceDipole,softening)
        +plummerMagnetizationField(sourceToTarget,sourceDipole,softening);
}

// The magnetization term of a MOVING moment at an observation event, in
// the lab.  In the moment's rest frame it is the static field mu0 m n(r')
// with r' the rest-frame displacement; being local, it carries no
// retardation.  The lab fields follow from the Lorentz transformation of a
// pure rest-frame B:
//
//     B = gamma B' - gamma^2/(gamma+1) v (v.B')/c^2,     E = -gamma v x B'.
//
// r' is taken from the displacement at equal lab time, r'_par = gamma r_par,
// as for a uniformly moving source (the term reaches only a few floors).
inline ElectromagneticField movingMagnetizationField(
    const Vec3& sourceToTarget,const Vec3& sourceVelocity,
    const Vec3& properDipole,double softening) {
    const double speedSquared=sourceVelocity.squaredNorm();
    if(!(speedSquared>0.0)) {
        return {{},plummerMagnetizationField(sourceToTarget,properDipole,
                                             softening)};
    }
    const double lorentzFactor=1.0/std::sqrt(1.0-speedSquared/(c*c));
    const Vec3 axis=sourceVelocity/std::sqrt(speedSquared);
    const double parallel=dot(sourceToTarget,axis);
    const Vec3 restDisplacement=sourceToTarget
        +axis*((lorentzFactor-1.0)*parallel);
    const Vec3 restField=plummerMagnetizationField(restDisplacement,
        properDipole,softening);
    const Vec3 magnetic=restField*lorentzFactor
        -sourceVelocity*(lorentzFactor*lorentzFactor/(lorentzFactor+1.0)
            *dot(sourceVelocity,restField)/(c*c));
    const Vec3 electric=cross(sourceVelocity,restField)*(-lorentzFactor);
    return {electric,magnetic};
}

inline ElectromagneticField retardedElectricDipoleFieldExact(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double poleSeparationFraction=1.0e-5) {
    if(historicalDipoleSourceIsStatic(
            history,present,sourceIsFirst,true)) {
        if(const double floor=separationFloor(); floor>0.0) {
            const Vec3 position=sourceIsFirst
                ?present.firstPosition:present.secondPosition;
            const Vec3 moment=sourceIsFirst
                ?present.firstElectricDipole:present.secondElectricDipole;
            // 1/(4 pi epsilon0) = c^2 mu0/(4 pi)
            return {plummerDipoleField(observationPosition-position,
                                       moment,floor)*(c*c),{}};
        }
        return retardedElectricDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,false);
    }
    return twoChargeLimitDipoleField(observationPosition,observationTime,
        history,present,sourceIsFirst,poleSeparationFraction,separationFloor(),
        [&](double time,Vec3& moment,Vec3& first,Vec3& second) {
            const RetardedElectricDipoleKinematics dipole=
                historicalIntegratedDipoleKinematics(
                    history,present,sourceIsFirst,time,true);
            moment=dipole.moment;
            first=dipole.firstDerivative;
            second=dipole.secondDerivative;
        });
}

// Duality: feeding the magnetic moment scaled by 1/c^2 into the same
// two-charge construction produces the field of an electric dipole
// p=m/c^2; the actual magnetic-dipole field is that field's vacuum dual,
// (E,B)_m = (-c^2 B_p, E_p) -- see twoChargeLimitDipoleField's comment for
// the cross-check against this file's low-velocity formulas that pins down
// this exact sign and scale.
inline ElectromagneticField retardedMagneticDipoleFieldExact(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double poleSeparationFraction=1.0e-5) {
    if(historicalDipoleSourceIsStatic(
            history,present,sourceIsFirst,false)) {
        if(const double floor=separationFloor(); floor>0.0) {
            const Vec3 position=sourceIsFirst
                ?present.firstPosition:present.secondPosition;
            const Vec3 moment=sourceIsFirst
                ?present.firstDipole:present.secondDipole;
            return {{},plummerMagneticDipoleField(
                observationPosition-position,moment,magneticDipoleRadius())};
        }
        return retardedMagneticDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,false);
    }
    const ElectromagneticField dual=twoChargeLimitDipoleField(
        observationPosition,observationTime,history,present,sourceIsFirst,
        poleSeparationFraction,magneticDipoleRadius(),
        [&](double time,Vec3& moment,Vec3& first,Vec3& second) {
            const RetardedElectricDipoleKinematics dipole=
                historicalIntegratedDipoleKinematics(
                    history,present,sourceIsFirst,time,false);
            moment=dipole.moment/(c*c);
            first=dipole.firstDerivative/(c*c);
            second=dipole.secondDerivative/(c*c);
        });
    ElectromagneticField field{dual.magnetic*(-c*c),dual.electric};
    // The poles give the moving dipole's field of separated charges; a
    // magnetic moment is a loop, so its local magnetization term is added
    // (see plummerMagnetizationField).
    //
    // READ AT THE RETARDED TIME, like the rest of the field.  The first
    // version of this term took the source kinematics at the OBSERVATION
    // time and the proper moment from the present state, while the two-pole
    // construction beside it reads both at the retarded time.  The two
    // recipes respond differently when the observation event is displaced,
    // which left the coupling U with different left and right time
    // derivatives -- a one-sided kink of 1.65 in the term alone and 2.12 in
    // the total field (audit section 88).  That kink is what floored the
    // step-doubling error of the whole retarded sector at 1e-8 and made the
    // engine reject deep steps: with the term removed the same state
    // converges cleanly at 4.0x per halving down to 2e-10, and the poles
    // alone show left/right slopes equal to 1.0000.
    if(const double floor=separationFloor(); floor>0.0) {
        double magnetizationTime=observationTime;
        ChargeKinematics source=historicalCharge(
            history,present,sourceIsFirst,magnetizationTime);
        for(int iteration=0;iteration<16;++iteration) {
            const Vec3 displacement=observationPosition-source.position;
            const double distance=displacement.norm();
            const double residual=magnetizationTime+distance/c
                -observationTime;
            const double derivative=std::max(1.0e-8,
                1.0-dot(displacement/std::max(distance,
                    std::numeric_limits<double>::min()),
                    source.velocity/c));
            const double refined=magnetizationTime-residual/derivative;
            const bool converged=std::abs(refined-magnetizationTime)
                <=1.0e-30+1.0e-14*std::abs(magnetizationTime);
            magnetizationTime=refined;
            source=historicalCharge(
                history,present,sourceIsFirst,magnetizationTime);
            if(converged) break;
        }
        // The PROPER moment at the retarded time.  movingMagnetizationField
        // boosts a rest-frame magnetization, so it needs the proper moment;
        // historicalIntegratedDipoleKinematics returns mu_lab/gamma, the
        // moment of a source integral, which differs at second order in beta
        // and left the boosted loop off by 1.2e-2 at beta 0.3 (section 103).
        const State retardedState=
            historicalState(history,present,magnetizationTime);
        Vec3 properMoment=sourceIsFirst?retardedState.firstProperDipole
                                       :retardedState.secondProperDipole;
        if(properMoment.squaredNorm()==0.0)
            properMoment=sourceIsFirst?present.firstProperDipole
                                      :present.secondProperDipole;
        // CENTRED ON THE EXTRAPOLATED PRESENT POSITION, not on the retarded
        // one (audit section 103).  Reading the kinematics at the retarded
        // time is what removed the kink of section 88, but a boosted static
        // field is centred where the source IS at the observation time, as the
        // Lienard-Wiechert velocity field is.  Centring it on the retarded
        // position shifted it back by beta R along the motion: against the
        // boosted loop in uniform motion it erred by 8% at beta 0.05 and 58%
        // at beta 0.3 near 1 r*, and on engine trajectories the tangential
        // force it left did work, pumping energy into para and draining ortho
        // (sections 101e-102).
        //
        // The retarded event is extrapolated to the observation time to
        // SECOND order, with the retarded acceleration the pole fields beside
        // it already use.  Three recipes were measured (section 103):
        //   kinematics at the observation time: the ledger is exact again, but
        //     the one-sided slope ratio returns to 5.84 -- section 88's kink;
        //   first-order extrapolation: smooth (1.0000) and exact in uniform
        //     motion, but on a circular orbit it misses by a^2 dt^2/2 and
        //     turns the velocity by omega dt, leaving the ledger 1e3 off at
        //     2 r* (3.0e-3 against 2.3e-6);
        //   second order: smooth (1.0000), exact in uniform motion, and the
        //     ledger back to the pre-regression digits (2.31e-6 at 2 r*).
        //
        // The velocity is continued through the FOUR-velocity u = gamma v,
        // whose lab-time derivative is gamma a + gamma^3 v (v.a)/c^2: the same
        // slope at lag 0 as v + a lag, but it cannot reach c.  v + a lag did,
        // in deep plunges where a r/c is comparable to c (beta 1.06 measured
        // at 0.13 r*), and the Lorentz factor below turned the force into NaN
        // -- the "non-finite" engine failures of section 106 (section 107).
        const double lag=observationTime-magnetizationTime;
        const Vec3 extrapolatedSourcePosition=source.position
            +source.velocity*lag+source.acceleration*(0.5*lag*lag);
        const double sourceGamma=1.0/std::sqrt(std::max(
            1.0-source.velocity.squaredNorm()/(c*c),1.0e-300));
        const Vec3 extrapolatedFourVelocity=source.velocity*sourceGamma
            +(source.acceleration*sourceGamma
              +source.velocity*(sourceGamma*sourceGamma*sourceGamma
                  *dot(source.velocity,source.acceleration)/(c*c)))*lag;
        const Vec3 extrapolatedSourceVelocity=extrapolatedFourVelocity
            /std::sqrt(1.0+extrapolatedFourVelocity.squaredNorm()/(c*c));
        const ElectromagneticField magnetization=movingMagnetizationField(
            observationPosition-extrapolatedSourcePosition,
            extrapolatedSourceVelocity,
            properMoment,magneticDipoleRadius());
        field.electric+=magnetization.electric;
        field.magnetic+=magnetization.magnetic;
    }
    return field;
}

// Production fields: the full moving-point-dipole result everywhere the
// short-range regulator is inactive, and a smooth correction to the existing
// regularized potential inside its core.  Written as
//
//   F_reg = F_low,reg + w(r) [F_exact,point - F_low,point]
//
// so the static curl(A_reg) and its induction product-rule term are preserved
// exactly at short range, while kappa, aberration, convective and motional
// channels approach the exact result with the same declared weight.  In the
// physical exterior w is one to round-off and the wrappers return the exact
// field directly, avoiding two unnecessary low-velocity evaluations.
inline ElectromagneticField retardedElectricDipoleField(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst) {
    const ChargeKinematics current=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    const double rawDistance=(observationPosition-current.position).norm();
    const double modelFloor=separationFloor();
    if(modelFloor>0.0)
        return retardedElectricDipoleFieldExact(observationPosition,
            observationTime,history,present,sourceIsFirst);
    const double transition=modelFloor>0.0?std::clamp(
        (rawDistance-modelFloor)/modelFloor,0.0,1.0):1.0;
    const double domainWeight=transition*transition*(3.0-2.0*transition);
    const double weight=shortRangeFieldWeight(rawDistance)*domainWeight;
    const ElectromagneticField exact=retardedElectricDipoleFieldExact(
        observationPosition,observationTime,history,present,sourceIsFirst);
    if(weight>=1.0-1.0e-14) return exact;
    const ElectromagneticField regularizedLow=
        retardedElectricDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,true);
    const ElectromagneticField pointLow=
        retardedElectricDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,false);
    return {regularizedLow.electric
                +(exact.electric-pointLow.electric)*weight,
            regularizedLow.magnetic
                +(exact.magnetic-pointLow.magnetic)*weight};
}

inline ElectromagneticField retardedMagneticDipoleField(
    const Vec3& observationPosition,double observationTime,
    const StateHistory& history,const State& present,bool sourceIsFirst,
    double poleSeparationFraction=1.0e-5) {
    const ChargeKinematics current=historicalCharge(
        history,present,sourceIsFirst,observationTime);
    const double rawDistance=(observationPosition-current.position).norm();
    const double modelFloor=separationFloor();
    // With a floor the exact construction is itself Plummer-softened
    // (section 66), so no low-velocity core and no blend.
    if(modelFloor>0.0)
        return retardedMagneticDipoleFieldExact(observationPosition,
            observationTime,history,present,sourceIsFirst,
            poleSeparationFraction);
    const double transition=modelFloor>0.0?std::clamp(
        (rawDistance-modelFloor)/modelFloor,0.0,1.0):1.0;
    const double domainWeight=transition*transition*(3.0-2.0*transition);
    const double weight=shortRangeFieldWeight(rawDistance)*domainWeight;
    const ElectromagneticField exact=retardedMagneticDipoleFieldExact(
        observationPosition,observationTime,history,present,sourceIsFirst,
        poleSeparationFraction);
    if(weight>=1.0-1.0e-14) return exact;
    const ElectromagneticField regularizedLow=
        retardedMagneticDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,true);
    const ElectromagneticField pointLow=
        retardedMagneticDipoleFieldLowVelocity(observationPosition,
            observationTime,history,present,sourceIsFirst,false);
    return {regularizedLow.electric
                +(exact.electric-pointLow.electric)*weight,
            regularizedLow.magnetic
                +(exact.magnetic-pointLow.magnetic)*weight};
}

inline Vec3 regularizedDipoleVectorPotential(const Vec3& sourceToTarget,
                                      const Vec3& sourceDipole) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const MagneticRadialProfile profile = magneticRadialProfile(sourceToTarget.norm());
    return cross(sourceDipole, sourceToTarget)
         * (magneticConstant * profile.vectorPotentialFactor);
}

// B=curl(A) for A=mu0/(4 pi) f(r) mu x r.  In particular, this is not
// merely w times the unregularized point-dipole field.
inline Vec3 regularizedDipoleField(const Vec3& sourceToTarget,const Vec3& sourceDipole,
                            double regularizationRadius,double exponent) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const double distance = sourceToTarget.norm();
    if(!(distance>std::numeric_limits<double>::min())) return {};
    const Vec3 n = sourceToTarget / distance;
    const MagneticRadialProfile profile = magneticRadialProfile(
        distance,regularizationRadius,exponent);
    const double transverseCoefficient = 2.0 * profile.vectorPotentialFactor
                                       + distance * profile.firstDerivative;
    const double radialCoefficient = -distance * profile.firstDerivative;
    return (sourceDipole * transverseCoefficient
          + n * (dot(sourceDipole, n) * radialCoefficient)) * magneticConstant;
}

inline double regularizedDipoleInteractionEnergy(const Vec3& sourceToTarget,
    const Vec3& targetDipole,const Vec3& sourceDipole,
    double regularizationRadius=magneticRegularizationRadius,
    double exponent=magneticRegularizationExponent) {
    return -dot(targetDipole,regularizedDipoleField(
        sourceToTarget,sourceDipole,regularizationRadius,exponent));
}

// F=-grad(U), U=-mu_target.B_source, using the same regularized field above.
inline Vec3 regularizedDipoleForce(const Vec3& sourceToTarget,
                            const Vec3& targetDipole,const Vec3& sourceDipole,
                            double regularizationRadius=magneticRegularizationRadius,
                            double exponent=magneticRegularizationExponent) {
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const double distance = sourceToTarget.norm();
    if(!(distance>std::numeric_limits<double>::min())) return {};
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

// THE PAIR'S DIPOLE SECTOR, ONE POTENTIAL (audit section 66).
//
// With a separation floor the dipole field, the dipole-dipole energy and the
// dipole-dipole force are those of a Plummer-softened dipole, rho^2 = r^2 +
// floor^2 on the TRUE separation:
//
//     B = mu0/(4 pi) [3 r (r.m)/rho^5 - m/rho^3],   U = -m_t . B,
//
// which is exactly what two poles carrying the Plummer-softened charge field
// of lienardWiechertField produce as their separation goes to zero.  The
// retarded two-charge construction therefore reduces to these in the static
// limit, and the energies conservativeParticleEnergy counts are the
// potentials of the forces.  Before, the field was a point dipole at the
// CLAMPED separation blended with one at the true separation through a
// position-dependent weight, and unsoftened beyond two floors while the
// energy was softened: at the default floor its gradient force did about
// +0.33 k/r0 of work per passage that no energy term accounted for (65d).
// Pairs without a floor keep the magnetic regularization profile.

inline Vec3 pairDipoleField(const Vec3& sourceToTarget,
                            const Vec3& sourceDipole) {
    if(const double floor=separationFloor(); floor>0.0)
        return plummerMagneticDipoleField(sourceToTarget,sourceDipole,
                                          magneticDipoleRadius());
    return regularizedDipoleField(sourceToTarget,sourceDipole,
        magneticRegularizationRadius,magneticRegularizationExponent);
}

inline double pairDipoleInteractionEnergy(const Vec3& sourceToTarget,
    const Vec3& targetDipole,const Vec3& sourceDipole) {
    return -dot(targetDipole,pairDipoleField(sourceToTarget,sourceDipole));
}

// -grad of pairDipoleInteractionEnergy with respect to the target position.
inline Vec3 pairDipoleForce(const Vec3& sourceToTarget,
    const Vec3& targetDipole,const Vec3& sourceDipole) {
    const double floor=separationFloor();
    if(!(floor>0.0))
        return regularizedDipoleForce(sourceToTarget,targetDipole,sourceDipole);
    constexpr double magneticConstant=mu0/(4.0*pi);
    const double softening=magneticDipoleRadius();
    const double rhoSquared=sourceToTarget.squaredNorm()+softening*softening;
    const double inverseRho=1.0/std::sqrt(rhoSquared);
    const double inverseRhoFifth=inverseRho*inverseRho*inverseRho
        *inverseRho*inverseRho;
    const double inverseRhoSeventh=inverseRhoFifth*inverseRho*inverseRho;
    const double targetRadial=dot(targetDipole,sourceToTarget);
    const double sourceRadial=dot(sourceDipole,sourceToTarget);
    const double dipoleDot=dot(targetDipole,sourceDipole);
    // Last term: -grad of the magnetization energy
    // -(mu0/4 pi) 3 eps^2 (mu_t.mu_s)/rho^5, the contact attraction of
    // parallel moments (see plummerMagnetizationField).
    return (sourceToTarget*(3.0*dipoleDot*inverseRhoFifth
                -15.0*targetRadial*sourceRadial*inverseRhoSeventh)
           +(targetDipole*sourceRadial+sourceDipole*targetRadial)
                *(3.0*inverseRhoFifth)
           -sourceToTarget*(15.0*softening*softening*dipoleDot
                *inverseRhoSeventh))
        *magneticConstant;
}

inline MutualForces coulombForces(const State& s) {
    const PairGeometry geometry = clampedPairGeometry(s);
    const Vec3 first = clampedSeparationForceToTruePositions(
        s.firstPosition - s.secondPosition, separationFloor(),
        geometry.firstMinusSecond
            * (-pairCoulombStrength * geometry.inverseDistanceCubed));
    return {first, first * -1.0};
}

// CREM_NO_DIPOLE_FORCE: ablate the mutual dipole-dipole FORCE while leaving
// everything else about the moments in place -- they still precess, still
// radiate M1, still enter the annihilation invariant.  This exists to answer
// one question: with the spins quantized, para and ortho at the same seed
// differ ONLY by the sign of mu2 (the quantized branch derives secondDipole
// from firstDipole and consumes no extra draws), so every channel difference
// the model produces has to come from something that reads that sign.  The
// coherent M1 channel is bounded at 1e-12 of the radiated power and cannot
// move a 1e-5 observable; this force is the next candidate, and turning it
// off is how the candidacy gets tested rather than asserted.
inline bool gDipoleForceEnabled=std::getenv("CREM_NO_DIPOLE_FORCE")==nullptr;

inline MutualForces mutualForces(const State& s) {
    const MutualForces electrostatic = coulombForces(s);
    if(!gDipoleForceEnabled) return electrostatic;
    const Vec3 dipoleOnFirst = pairDipoleForce(
        s.firstPosition - s.secondPosition, s.firstDipole, s.secondDipole);
    return {electrostatic.first + dipoleOnFirst,
            electrostatic.second - dipoleOnFirst};
}

inline Vec3 darwinForceOnFirst(const Vec3& firstVelocity, const Vec3& secondVelocity,
                        const Vec3& secondLeadingAcceleration,
                        const Vec3& firstMinusSecond, double chargeProduct) {
    const double distance = firstMinusSecond.norm();
    // Both callers pass clampedPairGeometry's separation, which is floored at
    // separationFloor() -- but that returns 0 for every pair that is not
    // positronium, and the clamp is then the identity.  Coincidence stays
    // unreachable behind the terminal surface; this keeps the function's own
    // contract from depending on which pair is active.
    if(!(distance>std::numeric_limits<double>::min())) return Vec3{};
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

inline MutualForces darwinForces(const State& s) {
    const PairGeometry geometry = clampedPairGeometry(s);
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

inline double darwinInteractionEnergy(const State& s) {
    const PairGeometry geometry = clampedPairGeometry(s);
    const Vec3 n = geometry.firstMinusSecond * geometry.inverseDistance;
    const double chargeProduct = pairChargeProduct;
    return coulomb * chargeProduct * geometry.inverseDistance / (2.0 * c*c)
         * (dot(s.firstVelocity, s.secondVelocity)
          + dot(s.firstVelocity, n) * dot(s.secondVelocity, n));
}

// Osculating orbital angular frequency of the pair, sqrt(k|q1 q2|/(mu r^3)).
// This is the rail the orbit-following zero-point band rides on: as the orbit
// tightens the frequency rises and the whole band rises with it, instead of
// being left behind at the value the orbit had when the run started.
inline double osculatingOrbitalFrequency(const State& s) {
    const double radius=separation(s);
    if(!(radius>0.0)) return 0.0;
    return std::sqrt(pairCoulombStrength
        /(pairReducedMass*radius*radius*radius));
}

// Time derivative of the same radius-defined frequency.  This is needed by
// the optional chirped ZPF modes, whose fields are derived from a vector
// potential and therefore require both omega and d omega/dt.
inline double osculatingOrbitalFrequencyDerivative(const State& s) {
    const Vec3 separationVector=s.firstPosition-s.secondPosition;
    const double radius=separationVector.norm();
    if(!(radius>0.0)) return 0.0;
    const double radialVelocity=dot(separationVector,
        s.firstVelocity-s.secondVelocity)/radius;
    return -1.5*osculatingOrbitalFrequency(s)*radialVelocity/radius;
}

// THE FREQUENCY THE ZERO-POINT BAND RIDES.
//
// Production pins it to the pair's OSCULATING orbital frequency, so every
// mode stays in resonance the whole way down.  That choice has a history:
// with a fixed band the orbit outran it -- the measured period falls by a
// factor of 105 over a collapse -- and the lifetime then depended strongly
// on where the upper edge was cut.
//
// It also has a consequence audit 170 measured.  A band pinned to omega_orb
// carries an amplitude going as omega_orb^2, so the power it feeds in
// inherits the orbit's own r^-4; the radiated power carries the same r^-4,
// their ratio is flat, and the stochastic-electrodynamics balance
// P_abs = P_rad stops being an equation for a radius at all.  It becomes a
// condition on the field's normalization, satisfied at every radius at once
// or at none, which is the opposite of the scale-selecting mechanism the
// plan's Step 4 was looking for.
//
// THE RIDE IS A ONE-PARAMETER FAMILY, and the parameter is what audit 172
// scans.  Write
//     omega_band = omega_ref * (omega_orb/omega_ref)^p,
// so p = 1 is the production band, riding exactly, and p = 0 is a band
// frozen at omega_ref.  CREM_ZPF_RIDE_EXPONENT sets p and
// CREM_ZPF_REFERENCE_FREQUENCY sets omega_ref in rad/s;
// CREM_ZPF_FIXED_FREQUENCY remains as the p = 0 shorthand it was.
//
// WHY THE FAMILY IS THE RIGHT OBJECT.  Audits 170 and 171 measured the two
// endpoints and found them to fail for opposite reasons: p = 1 absorbs but
// inherits r^-4 so the ratio is flat, p = 0 breaks the inheritance but sits
// off resonance and absorbs nothing.  Between them two things compete.  The
// band has a finite width, [0.3, 3] times its centre, so it overlaps the
// orbit while omega_orb/omega_band stays inside about a factor of three;
// that ratio goes as omega_orb^(1-p), and omega_orb moves by 190 over the
// radii measured, so resonance survives only for |1-p| <~ 0.21.  Inside that
// window the amplitude still follows omega_band^2 but with a different
// exponent, which is what could restore a radial dependence to the balance.
//
// Unset, all three variables leave the production band bit for bit.
struct ZeroPointDrive { double frequency=0.0, derivative=0.0; };
inline ZeroPointDrive zeroPointDriveFrequency(const State& s) {
    struct Settings { double exponent, reference; };
    static const Settings settings=[]{
        const auto read=[](const char* name,double fallback){
            const char* text=std::getenv(name);
            if(!text) return fallback;
            const double value=std::atof(text);
            return std::isfinite(value)?value:fallback;
        };
        const double frozen=read("CREM_ZPF_FIXED_FREQUENCY",0.0);
        if(frozen>0.0) return Settings{0.0,frozen};
        const double reference=read("CREM_ZPF_REFERENCE_FREQUENCY",0.0);
        return Settings{read("CREM_ZPF_RIDE_EXPONENT",1.0),reference};
    }();
    const double orbital=osculatingOrbitalFrequency(s);
    const double rate=osculatingOrbitalFrequencyDerivative(s);
    if(settings.exponent==1.0) return {orbital,rate};
    if(!(settings.reference>0.0)||!(orbital>0.0)) return {orbital,rate};
    // omega_band = ref * (orbital/ref)^p, and its rate follows by the chain
    // rule: d(omega_band)/dt = p * (omega_band/omega_orb) * d(omega_orb)/dt.
    const double band=settings.reference
        *std::pow(orbital/settings.reference,settings.exponent);
    if(!(band>0.0)||!std::isfinite(band)) return {orbital,rate};
    return {band,settings.exponent*(band/orbital)*rate};
}

inline LocalElectromagneticFields localRelativisticFields(
    const State& s, const StateHistory& history) {
    ElectromagneticField atFirst = lienardWiechertField(
        s.firstPosition, s.time, history, s, false, secondCharge);
    ElectromagneticField atSecond = lienardWiechertField(
        s.secondPosition, s.time, history, s, true, firstCharge);
    const ElectromagneticField secondDipole=retardedMagneticDipoleField(
        s.firstPosition,s.time,history,s,false);
    const ElectromagneticField firstDipole=retardedMagneticDipoleField(
        s.secondPosition,s.time,history,s,true);
    // These electric components are generated by the exact moving magnetic
    // dipoles themselves.  State's electric dipoles are the boosted
    // components of those same proper magnetic moments, not extra sources;
    // adding them here separately would count the motional channel twice.
    // CREM_SPIN_RATE_SECTOR: diagnostic ablation of the precession rate
    // (audit section 116).  |S1+S2| is conserved exactly when omega1 =
    // omega2, so "which sector makes the two rates differ" is answered by
    // keeping one of them: "charge" drops the partner's dipole field and
    // leaves its Lienard-Wiechert charge field, "dipole" does the reverse.
    // It reaches the Thomas-BMT precession and the moment derivatives that
    // chargeDipoleForces reads, and nothing else: the Lorentz and gradient
    // forces assemble their own fields.  Unset (the default) leaves every
    // number bit-identical.
    // "isotropic" and "tensor" (audit section 117) split the partner's dipole
    // field the way the energy splits: the part proportional to the partner's
    // moment, which carries mu1.mu2 and commutes with S^2, against the part
    // along the separation, which carries (mu1.n)(mu2.n) and does not.  Both
    // are the STATIC Plummer forms rather than the retarded field, an
    // approximation worth 1e-5 at the orbital speeds here, and they are
    // diagnostics: only the precession sees them.
    static const int rateSector=[]{
        const char* text=std::getenv("CREM_SPIN_RATE_SECTOR");
        if(!text) return 0;
        if(std::strcmp(text,"charge")==0) return 1;
        if(std::strcmp(text,"dipole")==0) return 2;
        if(std::strcmp(text,"isotropic")==0) return 3;
        if(std::strcmp(text,"tensor")==0) return 4;
        return 0;
    }();
    if(rateSector>=2) { atFirst=ElectromagneticField{};
                        atSecond=ElectromagneticField{}; }
    if(rateSector<=0||rateSector==2) {
        atFirst.electric+=secondDipole.electric;
        atFirst.magnetic+=secondDipole.magnetic;
        atSecond.electric+=firstDipole.electric;
        atSecond.magnetic+=firstDipole.magnetic;
    } else if(rateSector>=3) {
        const double softening=magneticDipoleRadius();
        constexpr double magneticConstant=mu0/(4.0*pi);
        const auto part=[&](const Vec3& sourceToTarget,const Vec3& moment) {
            const double rhoSquared=sourceToTarget.squaredNorm()
                +softening*softening;
            if(!(rhoSquared>std::numeric_limits<double>::min())) return Vec3{};
            const double inverseRho=1.0/std::sqrt(rhoSquared);
            const double inverseRhoCubed=inverseRho*inverseRho*inverseRho;
            if(rateSector==4)
                return sourceToTarget*(3.0*dot(sourceToTarget,moment)
                    *inverseRhoCubed*inverseRho*inverseRho*magneticConstant);
            return moment*(-inverseRhoCubed*magneticConstant)
                +plummerMagnetizationField(sourceToTarget,moment,softening);
        };
        atFirst.magnetic+=part(s.firstPosition-s.secondPosition,s.secondDipole);
        atSecond.magnetic+=part(s.secondPosition-s.firstPosition,s.firstDipole);
    }
    // The external field is uniform, so both roles see the same addition and
    // no gradient force follows from it.  Adding it here rather than only in
    // the force sums is what carries it into Thomas-BMT precession, which for
    // a field this weak is the channel that actually does something.
    atFirst.magnetic+=gExternalMagneticField;
    atSecond.magnetic+=gExternalMagneticField;
    // The zero-point field, unlike the uniform one, differs between the two
    // positions and changes with time, so each role is sampled separately.
    if(gZeroPointField.active()) {
        Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
        const ZeroPointDrive drive=zeroPointDriveFrequency(s);
        const double orbitalFrequency=drive.frequency;
        const double orbitalFrequencyDerivative=drive.derivative;
        gZeroPointField.sample(s.firstPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,firstElectric,firstMagnetic);
        gZeroPointField.sample(s.secondPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,secondElectric,secondMagnetic);
        atFirst.electric+=firstElectric;
        atFirst.magnetic+=firstMagnetic;
        atSecond.electric+=secondElectric;
        atSecond.magnetic+=secondMagnetic;
    }
    return {atFirst, atSecond};
}

inline DipoleDerivatives thomasBmtDipoleDerivatives(
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

inline FourVector dipoleFourVector(const Vec3& properDipole,const Vec3& velocity) {
    const double relativisticGamma=gamma(velocity);
    const double projection=dot(velocity,properDipole);
    return {relativisticGamma*projection/c,
        properDipole+velocity*(relativisticGamma*relativisticGamma
            /(relativisticGamma+1.0)*projection/(c*c))};
}

inline Vec3 properDipoleFromFourVector(const FourVector& dipole,
                                const Vec3& velocity) {
    const double relativisticGamma=gamma(velocity);
    return dipole.space-velocity*(relativisticGamma
        /(relativisticGamma+1.0)*dipole.time/c);
}

inline FourVector electromagneticTensorAction(const ElectromagneticField& field,
                                        const FourVector& vector) {
    return {dot(field.electric,vector.space)/c,
        field.electric*(vector.time/c)+cross(vector.space,field.magnetic)};
}

// No longer on the production path (see properDipolePrecessionRate's
// comment below for why) -- kept for the self-consistency checks in
// maxwell_validation.hpp (POSITRONIUM_ENABLE_FIELD_VALIDATION only), hence
// unused in the plain production build.
[[maybe_unused]] inline Vec3 advanceCovariantBmt(const Vec3& properDipole,
                         const Vec3& velocity,
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

// The lab-frame magnetic moment obtained by boosting a purely-magnetic rest
// dipole (lorentzBoostDipole with zero electric input) scales only the
// component PERPENDICULAR to velocity by gamma; the component PARALLEL to
// velocity passes through unchanged.  Verified by direct expansion of
// lorentzBoostDipole's magnetic branch (README's "audyt fizyki" BMT
// section).
//
// Nothing calls this any more.  properDipolePrecessionRate below used to
// route Jackson's rate through this inverse boost, and that routing was
// exactly the bug (see its comment); with the rate applied to the rest-frame
// moment directly, no inverse boost is needed anywhere.  Kept, rather than
// deleted, only because the gamma-perpendicular/unity-parallel property it
// states is what makes synchronizeCovariantDipoles's forward boost readable,
// and it is the cheapest place to record it.  Delete it freely if that stops
// being worth a function.
[[maybe_unused]] inline Vec3 inverseTensorBoostMagnetic(const Vec3& labVector,const Vec3& velocity) {
    const double speedSquared=velocity.squaredNorm();
    if(!(speedSquared>0.0)) return labVector;
    const Vec3 axis=velocity*(1.0/std::sqrt(speedSquared));
    const double parallelComponent=dot(labVector,axis);
    const Vec3 parallelPart=axis*parallelComponent;
    return (labVector-parallelPart)*(1.0/gamma(velocity))+parallelPart;
}

// Rate of change of the PROPER (rest-frame) dipole, which is what Jackson's
// Thomas-BMT equation governs:
//
//     d(mu_proper)/dt = (q/m) mu_proper x B_eff,   B_eff as in 11.170,
//
// with B_eff supplied by thomasBmtEffectiveField and t the LABORATORY time.
// Jackson 11.170 states the equation for the rest-frame spin s; mu_proper is
// (g q / 2m) s with all three factors constant, so it obeys the same
// equation.  state.firstProperDipole/secondProperDipole ARE that rest-frame
// moment (see state.hpp), so the rate applies to them directly and nothing
// needs boosting either way.
//
// Two earlier routes are gone, for two different reasons.
//
// advanceCovariantBmt's four-vector route: its final projection, needed
// because a.u=0 is not preserved when u is held frozen for a sub-step (as it
// always is here -- applyDipolePrecession below is called twice per full
// step, each time with velocity fixed for that half), injected a spurious
// correction of the SAME order as the genuine precession whenever the
// anomalous (g/2-1) term is active: exactly zero at g=2, growing with both
// the anomaly and beta otherwise (up to 45% for a proton at beta=0.9).
//
// Its replacement then over-corrected in the other direction: it expanded
// mu_proper into the LAB dipole, applied Jackson's rate to THAT, and mapped
// the result back through inverseTensorBoostMagnetic above.  Jackson's mu is
// the rest-frame moment, not the lab one, so that composition L^-1[L(mu) x
// B_eff] is not a rotation of mu_proper at all: mu_proper.d(mu_proper)/dt
// came out non-zero at the relative level gamma-1 (2.3e-5 at positronium's
// own beta, 4.3e-3 at beta=0.1, 0.63 at beta=0.9), independent of g -- so
// unlike the four-vector bug it did NOT vanish at g=2.  The norm drift it
// produced was absorbed step by step by advanceThomasBmtDipole's
// renormalization and misread there as RK4 truncation; measured with the
// renormalization removed, the drift per unit time is CONSTANT at 7.16e8/s
// across dt from 1e-18 down to 1e-22, i.e. first order in dt, where an
// O(dt^5) truncation term would fall by 1e4 per decade.
//
// The form below is a literal cross product, so it is a pure precession by
// construction: it conserves |mu_proper| exactly, as a rest-frame moment's
// magnitude must be conserved under precession in ANY frame, and it needs no
// inverse boost, no four-vector and no projection.
// Not called by production any more -- advanceThomasBmtDipole below solves
// this equation in closed form rather than integrating it.  It remains the
// SPECIFICATION of what that closed form must solve, and the validation
// suite checks the solver against it, so it is deliberately kept rather
// than folded into the solver.
[[maybe_unused]] inline Vec3 properDipolePrecessionRate(const Vec3& properDipole,const Vec3& velocity,
                                const ElectromagneticField& field,
                                double chargeToMass,double gFactor) {
    return cross(properDipole,
                 thomasBmtEffectiveField(velocity,field,gFactor))
        *chargeToMass;
}

// Closed-form solution of the precession over the sub-step, not a numerical
// integration of it.
//
// applyDipolePrecession holds velocity and field fixed for the whole call
// (two half-steps per full step, each at the velocity current for that
// half), so B_eff = thomasBmtEffectiveField(v,field,g) is a CONSTANT vector
// here.  The equation is then d(mu)/dt = omega x mu with a fixed
// omega = -(q/m) B_eff, whose exact solution is a rigid rotation of mu about
// omega by the angle |omega|*dt.  Rodrigues' formula gives that rotation
// directly, for any step size, with no truncation of any order.
//
// This replaces an RK4 pass followed by an explicit renormalization.  RK4's
// own truncation was never the problem -- at the step angles production
// actually takes (|omega|*dt ~ 2e-5 rad) its O(theta^5) error sits below
// round-off, measured identical to the closed form to 1e-16 -- but the
// renormalization that followed it was, for two reasons.
//
// Round-off: rescaling by targetNorm/result.norm() costs a sqrt and a
// division every sub-step, and their round-off accumulates over the ~1e6
// sub-steps of a trajectory faster than the rotation's own does.  Measured
// over 2e6 steps at production parameters: max |mu| drift 4.24e-11 for
// RK4-plus-rescale against 2.79e-12 for the rotation below, a factor of 15.
//
// Structure, which matters more.  A rescale REPAIRS the norm, and a repair
// is indistinguishable from there being nothing to repair.  That is exactly
// how this sector's previous bug hid: while properDipolePrecessionRate
// routed the rate through the lab dipole it was not a rotation of
// properDipole at all, and the rescale silently absorbed a first-order,
// physical-size error every step (see that function's comment).  A rotation
// conserves the norm BY CONSTRUCTION, so nothing is left that could mask a
// future structural defect: any such defect now has to show up as drift
// instead of being quietly removed.  The validation suite's
// bmt-precession-invariant check guards the same property from the other
// side.
inline Vec3 advanceThomasBmtDipole(const Vec3& properDipole,const Vec3& velocity,
                            const ElectromagneticField& field,
                            double chargeToMass,double laboratoryDt,
                            double gFactor) {
    if(laboratoryDt==0.0) return properDipole;
    // Angular velocity of the precession.  d(mu)/dt = (q/m) mu x B_eff is
    // omega x mu with omega = -(q/m) B_eff; the sign lives here so that the
    // rotation below is a rotation about omega in the usual sense.
    const Vec3 angularVelocity=
        thomasBmtEffectiveField(velocity,field,gFactor)*(-chargeToMass);
    const double angularSpeed=angularVelocity.norm();
    if(!(angularSpeed>0.0)) return properDipole;
    const Vec3 axis=angularVelocity/angularSpeed;
    const double angle=angularSpeed*laboratoryDt;
    const double cosine=std::cos(angle),sine=std::sin(angle);
    // (1.0-cosine) is safe here for the reason written out above
    // rotateDipoleByAngularVelocity in secular_spin_orbit.hpp: Sterbenz makes
    // the subtraction exact, and the cosine's own error cancels against
    // properDipole*cosine in the axis-parallel component.
    return properDipole*cosine+cross(axis,properDipole)*sine
        +axis*(dot(axis,properDipole)*(1.0-cosine));
}

inline void applyDipolePrecession(State& s, double dt,
                           const StateHistory& history) {
    synchronizeCovariantDipoles(s);
    const LocalElectromagneticFields fields = localRelativisticFields(s, history);
    const Vec3 firstDipole=advanceThomasBmtDipole(s.firstProperDipole,
        s.firstVelocity,fields.atFirst,firstCharge/firstMass,dt,firstGFactor);
    const Vec3 secondDipole=advanceThomasBmtDipole(s.secondProperDipole,
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
// where R points from the dipole to the charge.
//
// The RELATIVE velocity here is exact at this order, not a symmetrization
// imposed to make the two forces balance.  A magnetic dipole moving at v_mu
// carries a motional ELECTRIC dipole p = (v_mu x mu)/c^2, and that dipole's
// scalar potential contributes to the charge's Lagrangian
//
//   -q phi_p = -q/(4 pi eps0) w(r) (p.R)/r^3
//            = -q mu0/(4 pi) w(r) (v_mu x mu).R/r^3
//            = -q mu0/(4 pi) w(r) v_mu.(mu x R)/r^3   [(axb).c = a.(bxc)]
//            = -q v_mu.A_mu,
//
// using 1/(4 pi eps0 c^2) = mu0/4pi.  So q v_q.A_mu - q phi_p IS
// q (v_q - v_mu).A_mu: the single-particle coupling q v_q.A with the moving
// dipole's own electric field already in it, not an approximation of it.
// Newton's third law then FOLLOWS -- L depends only on R and on v_q - v_mu,
// so translation invariance forces the two variations to be equal and
// opposite -- rather than being the reason the form was chosen.
//
// Because the motional electric dipole is already accounted for here, this
// term must not be combined with retardedElectricDipoleField, which
// represents the same physics through state.firstElectricDipole.  It is not:
// allExternalForces uses this one and retardedExternalForces uses that one,
// and the two force sums are alternatives, never summed together.
inline ChargeDipolePairForces chargeDipolePairForces(
    const Vec3& relativeVelocity, double charge, const Vec3& sourceDipole,
    const Vec3& sourceDipoleDerivative, const Vec3& sourceToCharge) {
    const double distance = sourceToCharge.norm();
    const double floor = separationFloor();
    constexpr double magneticConstant = mu0 / (4.0 * pi);
    const Vec3 magneticField = pairDipoleField(sourceToCharge, sourceDipole);
    const double radialFactor = floor > 0.0
        ? 1.0/std::pow(distance*distance
              + magneticDipoleRadius()*magneticDipoleRadius(), 1.5)
        : magneticRadialProfile(distance).vectorPotentialFactor;
    const Vec3 inducedElectricTerm = cross(sourceDipoleDerivative, sourceToCharge)
        * (magneticConstant * radialFactor);
    const Vec3 onCharge =
        (cross(relativeVelocity, magneticField) - inducedElectricTerm) * charge;
    return {onCharge, onCharge * -1.0};
}

// Orbital back-reaction of Thomas precession.
//
// chargeDipolePairForces derives the charge-dipole force from
// L_qmu = q (v_charge - v_dipole).A_mu, which couples the moment only to the
// magnetic field in its rest frame.  The moments themselves follow
// Thomas-BMT, whose precession also contains the kinematic Thomas term
// omega_T = -(v x a)/(2c^2).  A force law without its counterpart does not
// conserve J = L + S: measured on Kepler orbits (e=0.3, 0.5; a_pair and
// 0.1 a_pair) the orbit-averaged torque of the instantaneous sum was 4/3 of
// -dS/dt for e+e- and 2x for the electron of p+e-, while the retarded
// assembly balanced it to 3e-4.  4/3 and 2 are exactly the Barker-O'Connell
// spin-orbit weights with and without the Thomas term: for spin p in the
// Coulomb field of o the orbit couples through
//
//     g_p/(2 m_p m_o) + (g_p - 1)/(2 m_p^2)      (with Thomas)
//     g_p/(2 m_p m_o) + g_p/(2 m_p^2)            (magnetic only),
//
// i.e. 3/2 against 2 for equal masses and 1/2 against 1 for m_o >> m_p.
//
// The missing term is the Lagrangian of Thomas precession,
//
//     L_T = -omega_T.S_p = (q_p/(2 m_p c^2)) v_p.(E_o(r_p) x S_p)
//         = v_p.A_T,   A_T = C f(d) (d x S_p),   C = k q_o q_p/(2 m_p c^2),
//
// d = r_p - r_o and f the same radial profile the charge's field uses
// (Plummer rho^-3 when the floor is on).  Euler-Lagrange with L_T depending
// on r_p - r_o only:
//
//     F_p =  grad_d(v_p.A_T) - dA_T/dt,     F_o = -grad_d(v_p.A_T),
//     dA_T/dt = C [ (grad f.v_rel)(d x S) + f (v_rel x S) + f (d x S') ].
//
// Checked against the retarded assembly (audit section 80): orbit-averaged
// J balance 4e-6..4e-4, apsidal rate equal to 5e-4 (e+e-) and 3e-3 (p+e-).
// Node by node the two differ by total time derivatives that average out.
// L_T is linear in the velocities, so it adds nothing to the energy function
// at fixed S.  CREM_NO_THOMAS_BACKREACTION removes it, for comparison.
inline MutualForces thomasBackReactionForces(const State& s,
                                             const DipoleDerivatives& derivatives) {
    const Vec3 firstMinusSecond=s.firstPosition-s.secondPosition;
    const double distance=firstMinusSecond.norm();
    const double floor=separationFloor();
    if(!(distance>0.0)&&!(floor>0.0)) return {};
    double radialFactor=0.0,gradientFactor=0.0;   // f, and grad f = g d
    if(floor>0.0) {
        const double rhoSquared=distance*distance+floor*floor;
        radialFactor=1.0/std::pow(rhoSquared,1.5);
        gradientFactor=-3.0/std::pow(rhoSquared,2.5);
    } else {
        const MagneticRadialProfile profile=magneticRadialProfile(distance);
        radialFactor=profile.vectorPotentialFactor;
        gradientFactor=profile.firstDerivative/distance;
    }
    const auto contribution=[&](bool spinFirst) {
        const double gyromagnetic=spinFirst?firstGyromagneticRatioOf()
                                           :secondGyromagneticRatioOf();
        if(!(std::abs(gyromagnetic)>0.0)) return MutualForces{};
        const double spinCharge=spinFirst?firstCharge:secondCharge;
        const double spinMass=spinFirst?firstMass:secondMass;
        const double otherCharge=spinFirst?secondCharge:firstCharge;
        const Vec3 spinVelocity=spinFirst?s.firstVelocity:s.secondVelocity;
        const Vec3 otherVelocity=spinFirst?s.secondVelocity:s.firstVelocity;
        const Vec3 spin=(spinFirst?s.firstDipole:s.secondDipole)/gyromagnetic;
        const Vec3 spinRate=
            (spinFirst?derivatives.first:derivatives.second)/gyromagnetic;
        const Vec3 d=spinFirst?firstMinusSecond:firstMinusSecond*-1.0;
        const double coefficient=coulomb*otherCharge*spinCharge
            /(2.0*spinMass*c*c);
        const Vec3 spinCrossVelocity=cross(spin,spinVelocity);
        const Vec3 gradient=(spinCrossVelocity*radialFactor
            +d*(gradientFactor*dot(d,spinCrossVelocity)))*coefficient;
        const Vec3 relativeVelocity=spinVelocity-otherVelocity;
        const Vec3 potentialRate=(cross(d,spin)
                *(gradientFactor*dot(d,relativeVelocity))
            +cross(relativeVelocity,spin)*radialFactor
            +cross(d,spinRate)*radialFactor)*coefficient;
        const Vec3 onSpin=gradient-potentialRate;
        const Vec3 onOther=gradient*-1.0;
        return spinFirst?MutualForces{onSpin,onOther}
                        :MutualForces{onOther,onSpin};
    };
    const MutualForces first=contribution(true);
    const MutualForces second=contribution(false);
    return {first.first+second.first,first.second+second.second};
}

inline MutualForces chargeDipoleForces(const State& s, const StateHistory& history) {
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

    static const bool thomasBackReaction=
        std::getenv("CREM_NO_THOMAS_BACKREACTION")==nullptr;
    const MutualForces thomas=thomasBackReaction
        ?thomasBackReactionForces(s,derivatives):MutualForces{};
    return {firstChargeSecondDipole.onCharge + secondChargeFirstDipole.onDipole
                + thomas.first,
            secondChargeFirstDipole.onCharge + firstChargeSecondDipole.onDipole
                + thomas.second};
}

inline MutualForces allExternalForces(const State& s) {
    const MutualForces positionForces = mutualForces(s);
    const MutualForces velocityForces = darwinForces(s);
    // A one-state history for the Thomas-BMT moment derivatives inside
    // chargeDipoleForces, with the stored accelerations ZEROED.  That history
    // extrapolates the partner back by r/c, and inside a step the stored
    // acceleration is the previous step's value, so with it kept the force
    // depended on the step path and the step-doubling velocity error halved
    // per halving of dt (a local error of order dt, 3.1e-9 at a captured
    // tilted para state at 0.40 r*).  Without chargeDipoleForces, or with the
    // accelerations zeroed, it fell 4x per halving (3.8e-11): the
    // extrapolation is then uniform motion, fixed by position and velocity.
    // Neither the dipole-dipole force nor the precession mattered.
    State localState{s};
    localState.firstAcceleration={};
    localState.secondAcceleration={};
    const StateHistory localHistory{localState};
    const MutualForces mixedMagneticForces = chargeDipoleForces(s, localHistory);
    MutualForces externalField{
        lorentzForce(firstCharge, s.firstVelocity, {{}, gExternalMagneticField}),
        lorentzForce(secondCharge, s.secondVelocity, {{}, gExternalMagneticField})};
    if(gZeroPointField.active()) {
        Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
        const ZeroPointDrive drive=zeroPointDriveFrequency(s);
        const double orbitalFrequency=drive.frequency;
        const double orbitalFrequencyDerivative=drive.derivative;
        gZeroPointField.sample(s.firstPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,firstElectric,firstMagnetic);
        gZeroPointField.sample(s.secondPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,secondElectric,secondMagnetic);
        externalField.first=externalField.first+lorentzForce(firstCharge,
            s.firstVelocity,{firstElectric,firstMagnetic})
            +gZeroPointField.gradientForce(s.firstPosition,orbitalFrequency,
                s.zeroPointPhase,s.firstDipole);
        externalField.second=externalField.second+lorentzForce(secondCharge,
            s.secondVelocity,{secondElectric,secondMagnetic})
            +gZeroPointField.gradientForce(s.secondPosition,orbitalFrequency,
                s.zeroPointPhase,s.secondDipole);
    }
    return {positionForces.first + velocityForces.first
                + mixedMagneticForces.first + externalField.first,
            positionForces.second + velocityForces.second
                + mixedMagneticForces.second + externalField.second};
}

// ONE SOFTENING LENGTH FOR BOTH SIDES OF THE MOMENT-CHARGE INTERACTION
// (audit 126), on by default; CREM_NO_MATCHED_SOFTENING restores the old
// mismatch.  The moment is a loop of radius magneticDipoleRadius(), so it
// samples the partner's charge field averaged over its own extent, exactly as
// the charge samples the moment's field softened at that same length.  Audit
// 125 traced the pair's self-propulsion to the mismatch -- 0.05 r* acting on
// the moment against 0.9668 r* acting back on the charge -- and audit 126
// measured what matching them does: in the synthetic circular test at 1 r*
// the net force falls 868 times, on a real orbit the gained momentum falls
// from 1.1714 to 0.0061 m c with the plunge intact, and in the exact
// uniform-motion case both sides then deviate from the point-dipole answer by
// the SAME 1.05% while their net matches the analytic net to four digits.
// Production is untouched: it runs at 15 r* and beyond, where a 0.9668 r*
// softening is nothing -- the collapse time moves by 0.0003 ps and a cascade
// lifetime by 1e-9 ps.
inline double matchedMomentSoftening() {
    static const bool enabled=
        std::getenv("CREM_NO_MATCHED_SOFTENING")==nullptr;
    return enabled?magneticDipoleRadius():0.0;
}

inline ElectromagneticField fieldFromOtherParticleAt(
    const Vec3& observationPosition,double observationTime,
    const State& state,const StateHistory& history,bool targetIsFirst,
    double poleSeparationFraction=1.0e-5,
    double chargePlummerFloor=0.0) {
    const bool sourceIsFirst=!targetIsFirst;
    const double sourceCharge=sourceIsFirst?firstCharge:secondCharge;
    ElectromagneticField field=lienardWiechertField(
        observationPosition,observationTime,history,state,sourceIsFirst,
        sourceCharge,0.0,chargePlummerFloor);
    const ElectromagneticField magneticDipole=retardedMagneticDipoleField(
        observationPosition,observationTime,history,state,sourceIsFirst,
        poleSeparationFraction);
    field.electric+=magneticDipole.electric;
    field.magnetic+=magneticDipole.magnetic;
    return field;
}

// Field derivative DU/Dt = partial_t U + v.grad U for
// U = mu_lab.B + p_lab.E, with the target tensor held fixed at the event.
// Sample the instantaneous TANGENT worldline, not the target's historical
// positions or moments: this is a directional derivative of the field, not
// the torque/BMT evolution of the particle.  Both samples precede the present,
// and use the same causal source-history extension as the spatial stencil.
inline double dipoleCouplingMaterialRate(const State& state,
                                  const StateHistory& history,
                                  bool targetIsFirst) {
    const Vec3 position=targetIsFirst?state.firstPosition:state.secondPosition;
    const Vec3 velocity=targetIsFirst?state.firstVelocity:state.secondVelocity;
    const Vec3 magnetic=targetIsFirst?state.firstDipole:state.secondDipole;
    const Vec3 electric=targetIsFirst?state.firstElectricDipole
                                      :state.secondElectricDipole;
    const auto couplingAt=[&](double offset) {
        const Vec3 point=position+velocity*offset;
        ElectromagneticField field=fieldFromOtherParticleAt(
            point,state.time+offset,state,history,targetIsFirst,1.0e-5,
            matchedMomentSoftening());
        // Match the spatial stencil's pole-cancellation recovery.
        if(gPoleCancellationRatio>30.0) {
            const ElectromagneticField retreated=fieldFromOtherParticleAt(
                point,state.time+offset,state,history,targetIsFirst,1.0e-7);
            if(gPoleCancellationRatio<=30.0) field=retreated;
        }
        return dot(magnetic,field.magnetic)+dot(electric,field.electric);
    };
    // Resolve the same spacetime scale as the spatial gradient.  Using two
    // history-node spacings here gives a finite travel interval, potentially
    // comparable to r/c, instead of a local derivative for a fast target.
    const double derivativeStep=std::max(
        1.0e-4*separation(state),1.0e-3*nuclearCutoff)/c;
    const double now=couplingAt(0.0);
    const double before=couplingAt(-derivativeStep);
    const double twiceBefore=couplingAt(-2.0*derivativeStep);
    return (3.0*(now-before)-(before-twiceBefore))/(2.0*derivativeStep);
}

// Rate of the hidden momentum a current loop carries in an external electric
// field, as a force: F = -(1/c^2) d(mu x E)/dt along the particle's own
// worldline (audit section 123).
//
// Why it belongs in the force sum.  A magnetic moment realized as a current
// loop carries momentum (mu x E)/c^2 that is not in gamma m v, so Newton's
// law for the MECHANICAL momentum reads dp/dt = F_gradient - d/dt (mu x E)/c^2.
// Audit 122 measured the model's own momentum imbalance in a clean synthetic
// setting -- first order in v/c, 92% of it in the force on the moments in the
// partner's charge field -- and found this rate to be 46% of it, with the same
// scaling.  The model had no such term.
//
// The derivative is taken the way dipoleCouplingMaterialRate takes its own:
// along the instantaneous tangent worldline, with a backward three-point
// stencil at the same spacetime step, so the two see the same field history.
// d(mu x E)/dt = (dmu/dt) x E + mu x (DE/Dt), with dmu/dt the Thomas-BMT
// precession from this history.
inline Vec3 hiddenMomentumRateForce(const State& state,
                              const StateHistory& history,
                              bool targetIsFirst) {
    const Vec3 moment=targetIsFirst?state.firstDipole:state.secondDipole;
    if(moment.squaredNorm()==0.0) return {};
    const Vec3 position=targetIsFirst?state.firstPosition:state.secondPosition;
    const Vec3 velocity=targetIsFirst?state.firstVelocity:state.secondVelocity;
    const double derivativeStep=std::max(
        1.0e-4*separation(state),1.0e-3*nuclearCutoff)/c;
    const auto electricAt=[&](double offset) {
        return fieldFromOtherParticleAt(position+velocity*offset,
            state.time+offset,state,history,targetIsFirst,1.0e-5,
            matchedMomentSoftening()).electric;
    };
    // THE THREE SAMPLES READ ONE HISTORY SEGMENT (audit 151).  This is the
    // analytic rate that audits 149g and 150e asked for, and it needs no new
    // algebra: RetardedSegmentPin exists so that a group of retarded reads
    // sees the same polynomial, and a cubic is its own third-order Taylor
    // series, so one segment continued past its ends reproduces the
    // interpolant exactly and smoothly.  Differencing inside it is therefore
    // differentiating a smooth function, not stepping across the C0 jump that
    // audit 149 measured at 4.01% of the field and 15.40 times the Coulomb
    // force.  covariantDipoleGradientForce pins its six spatial probes for
    // the identical reason and against the identical symptom.
    //
    // The samples sit at observation times 0, -h and -2h and at positions
    // shifted by v times those offsets, so their retarded times span about
    // 2h(1+beta)/(1-beta); eight step lengths cover that for any speed this
    // integrator reaches, which is the same margin the gradient force uses.
    const RetardedSegmentPinGuard rateSegment(retardedSegmentPinAt(
        history,state,!targetIsFirst,position,state.time,
        8.0*derivativeStep));
    // THE POLE SECTOR'S SHARE OF THE RATE IS TAKEN ANALYTICALLY (audit 163).
    // Audit 161 separated this rate by sector: the charge field differences
    // cleanly at 1.1e-08, the two-pole dipole field at 1.4e-04, because the
    // latter's ~1e-10 cancellation floor is divided by a step of about
    // 1/2000 of the field's own time scale.  So the pole part is computed
    // with its derivative carried through the same arithmetic
    // (twoChargeLimitDipoleElectricWithRate), and only the REMAINDER -- the
    // charge field and the magnetization term, neither of which cancels --
    // is still differenced.  Where the pole construction is not the path
    // production takes, the correction is identically zero and this
    // degenerates to the stencil it replaces.
    const bool sourceIsFirstHere=!targetIsFirst;
    const auto momentSampler=[&](double time,Vec3& moment,Vec3& first,
                                 Vec3& second) {
        const RetardedElectricDipoleKinematics dipole=
            historicalIntegratedDipoleKinematics(
                history,state,sourceIsFirstHere,time,false);
        moment=dipole.moment/(c*c);
        first=dipole.firstDerivative/(c*c);
        second=dipole.secondDerivative/(c*c);
    };
    static const bool analyticPoleRate=
        std::getenv("CREM_NO_ANALYTIC_POLE_RATE")==nullptr;
    bool poleSectorActive=analyticPoleRate&&separationFloor()>0.0
        &&!historicalDipoleSourceIsStatic(history,state,sourceIsFirstHere,
                                          false);
    const auto poleElectricAt=[&](double offset)->Vec3 {
        if(!poleSectorActive) return {};
        const DipoleElectricWithRate pole=
            twoChargeLimitDipoleElectricWithRate(position+velocity*offset,
                velocity,state.time+offset,history,state,sourceIsFirstHere,
                1.0e-5,magneticDipoleRadius(),momentSampler);
        return pole.valid?pole.magnetic*(-c*c):Vec3{};
    };
    // The offset-zero construction is evaluated ONCE: it carries both the
    // value the remainder needs and the rate that replaces the stencil's
    // pole share.  Three two-pole evaluations are added in total, not five.
    Vec3 analyticPoleContribution, poleElectricNow;
    if(poleSectorActive) {
        const DipoleElectricWithRate pole=
            twoChargeLimitDipoleElectricWithRate(position,velocity,
                state.time,history,state,sourceIsFirstHere,1.0e-5,
                magneticDipoleRadius(),momentSampler);
        // An invalid construction leaves the remainder whole and the
        // correction zero, which is the stencil's own answer.
        if(pole.valid) {
            analyticPoleContribution=pole.magneticRate*(-c*c);
            poleElectricNow=pole.magnetic*(-c*c);
        } else poleSectorActive=false;
    }
    const Vec3 now=electricAt(0.0)-poleElectricNow;
    const Vec3 before=electricAt(-derivativeStep)
        -poleElectricAt(-derivativeStep);
    const Vec3 twiceBefore=electricAt(-2.0*derivativeStep)
        -poleElectricAt(-2.0*derivativeStep);
    // LIMITED STENCIL (audit 150).  The step above is about 1/1500 of the
    // history's own node spacing, and the retarded field is only C0 across a
    // node: audit 149 measured a +4.01% jump in the dipole field's electric
    // part, which divided by this step becomes 15.4 times the Coulomb force
    // and killed every tilted orbit inside 2.5 r* at step 6 of 256.
    //
    // The two half-slopes carry the information needed to tell the two cases
    // apart.  Where the field is smooth they differ by O(h), and
    // (3 d1 - d2)/2 -- which is exactly the three-point form this replaces --
    // is the second-order derivative.  Where the stencil straddles a jump,
    // one half is the jump divided by h and the two disagree by order one.
    // So each component keeps the second-order value while the halves agree
    // and falls back to the smaller half (minmod, zero on opposite signs)
    // when they do not.  That refuses to differentiate a discontinuity
    // instead of amplifying it.
    //
    // It is a numerical treatment and not a physical one.  The honest repair
    // is an analytic d(mu x E)/dt, expressible from the moment's precession
    // and the field's own retarded derivatives; this is the cheap one that
    // makes the term integrable.  CREM_HIDDEN_RATE_UNLIMITED restores the
    // bare three-point form for reproducing earlier numbers.
    //
    // AUDIT 161 MEASURED BOTH HALVES OF THAT AND NARROWED THE TARGET.  The
    // limiter below is dormant in the production range -- it fired 0 of 8110
    // stencil evaluations over an inspiral from 2 r* to 1.2 r* -- because
    // audit 151's pin removed the node crossings outright.  So nothing here
    // is being clipped, and 149's failure no longer reproduces: that orbit
    // now completes 256/256 at tolerances 1e-08 and 1e-09.
    //
    // The repair is still owed, and its target is not what the paragraph
    // above implies.  (dmu/dt) x E is ALREADY analytic, from
    // thomasBmtDipoleDerivatives below.  Of the remaining mu x DE/Dt, the
    // charge field's share is clean: separated by sector at 0.974 r*, the
    // rate's worst relative scatter near h is 1.1e-08 for the
    // Lienard-Wiechert field against 1.4e-04 for the pole-limit dipole
    // field.  The residue is entirely the two-pole construction's
    // cancellation floor, ~1e-10 of the field, which 1/h multiplies by about
    // 2000 -- an amplification, not a lost order, so a higher-order stencil
    // would not touch it and Richardson extrapolation would make it worse.
    //
    // What would: a per-pole analytic Lienard-Wiechert derivative combined by
    // the SAME subtraction twoChargeLimitDipoleField already performs, which
    // moves the floor from 1e-10/h to 1e-10 and needs no new field algebra.
    // It is not small work -- the path runs through two nested retarded-time
    // Newton solves, a Taylor reconstruction off one pinned cubic, and two
    // std::max clamps on the pole distance that are not differentiable at
    // exactly the near-floor geometry they exist for.
    //
    // It is worth doing: the term carries 1.4% of the total force at 5 r* and
    // 4.4% at 0.5 r*, so the rate's 2e-05 of scatter is about 7e-07 of the
    // total force at 1 r*, already past the engine's production relative
    // tolerance of 1e-07, and both factors grow inward.
    // The escape hatch evaluates the ORIGINAL expression verbatim rather than
    // an algebraically equal rearrangement: (3 d1 - d2)/2 and
    // (3 now - 4 before + twiceBefore)/(2h) agree in exact arithmetic and not
    // in double, and the difference moves a chaotic cascade in its fourth
    // digit.  Reproducing earlier numbers is the hatch's whole purpose, so it
    // has to be bit-faithful.
    static const bool unlimitedRate=
        std::getenv("CREM_HIDDEN_RATE_UNLIMITED")!=nullptr;
    const Vec3 fieldRate=[&]() -> Vec3 {
        if(unlimitedRate)
            return (now*3.0-before*4.0+twiceBefore)
                *(1.0/(2.0*derivativeStep));
        // The straddle test is on the FIELD's own scale, not on the slopes:
        // two slopes can disagree by round-off when both are tiny, and an
        // earlier version of this limiter engaged on exactly that noise and
        // moved a smooth synthetic case in its third digit.  A node crossing
        // instead shows up as a second difference that is a PERCENT of the
        // field, against O((h/tau)^2) ~ 1e-09 where the field is smooth, so
        // the two cases are ten decades apart and the threshold is not
        // delicate.
        const Vec3 secondDifference=now-before*2.0+twiceBefore;
        const double fieldScale=std::max(now.norm(),1.0e-300);
        const Vec3 firstHalf=(now-before)*(1.0/derivativeStep);
        const Vec3 secondHalf=(before-twiceBefore)*(1.0/derivativeStep);
        if(secondDifference.norm()<=1.0e-4*fieldScale)
            return (now*3.0-before*4.0+twiceBefore)
                *(1.0/(2.0*derivativeStep));
        const auto minmod=[](double first,double second) {
            if(first*second<=0.0) return 0.0;
            return std::abs(first)<std::abs(second)?first:second;
        };
        return Vec3{minmod(firstHalf.x,secondHalf.x),
                    minmod(firstHalf.y,secondHalf.y),
                    minmod(firstHalf.z,secondHalf.z)};
    }();
    const DipoleDerivatives derivatives=
        thomasBmtDipoleDerivatives(state,history);
    const Vec3 momentRate=targetIsFirst?derivatives.first:derivatives.second;
    // `now` is the remainder, so the FIELD in the first cross product has to
    // be put back together before it is used.
    const Vec3 fullField=now+poleElectricNow;
    const Vec3 totalFieldRate=fieldRate+analyticPoleContribution;
    return (cross(momentRate,fullField)+cross(moment,totalFieldRate))
        *(-1.0/(c*c));
}

inline Vec3 covariantDipoleGradientForce(const State& state,
                                  const StateHistory& history,
                                  bool targetIsFirst) {
    const Vec3 targetPosition=targetIsFirst?state.firstPosition
                                               :state.secondPosition;
    const Vec3 targetVelocity=targetIsFirst?state.firstVelocity
                                               :state.secondVelocity;
    // U = mu_lab.B + p_lab.E is the scalar coupling in this tensor
    // convention.  Differentiate the field with the target tensor fixed;
    // project its four-gradient onto the target's rest space below.
    const Vec3 labMagneticDipole=targetIsFirst?state.firstDipole
                                                 :state.secondDipole;
    const Vec3 labElectricDipole=targetIsFirst?state.firstElectricDipole
                                                  :state.secondElectricDipole;
    if(labMagneticDipole.squaredNorm()==0.0
       &&labElectricDipole.squaredNorm()==0.0) return {};
    // The moving-dipole field is evaluated as a symmetric pole limit.  Its
    // own cancellation error is ~1e-10; differentiating it with the former
    // 1e-6*r step amplified that floor to ~1e-4 and made the adaptive force
    // estimate history-dependent for heavy pairs.  A 1e-4*r stencil leaves
    // O(h^2)=1e-8 spatial truncation while limiting that amplification to
    // ~1e-6.  The independent boosted coupling and trajectory convergence
    // checks below/within validation guard both sides of this balance.
    const double gradientStep=std::max(
        1.0e-4*separation(state),1.0e-3*nuclearCutoff);
    // DETECT AND RETREAT.  The moving-dipole field is built as the eps->0
    // limit of two point charges, so the subtraction that produces it retains
    // only a fraction eps/d of each pole.  At most geometries that is exactly
    // what happens and the result is independent of eps to nine or ten digits.
    // At a rare few it is not: measured, one probe of this very stencil
    // retained 180x too much and returned a field 7.5x oversized, which is
    // what makes this force discontinuous in the integrator's step size and
    // kills the trajectory with reason=accuracy.
    //
    // gPoleCancellationRatio is that retention divided by its design value, so
    // ~1 when the construction is behaving (measured 1.09) and large when it
    // is not (measured 197).  It costs two norms of vectors already computed,
    // so the healthy path -- which the scan below found to be 466 of 466
    // sampled geometries -- pays nothing at all.
    //
    // On detection the separation retreats by two decades rather than being
    // refined gradually, and that is not a free choice: the error is NOT
    // monotone in eps.  Sweeping the degenerate probe gave 892k, 2.16M, 6.79M,
    // 23.4M for eps/d = 1e-4, 3e-5, 1e-5, 3e-6 before collapsing to the
    // correct 907.8k at 1e-6 and 1e-7.  Halving from production would land on
    // 3e-6, four times WORSE than where it started, so anything that refines
    // in small steps makes this defect worse.
    //
    // The retreat then verifies itself: if the fallback's own cancellation is
    // also anomalous the original value is kept, which leaves today's
    // behaviour (a reported numerical failure) rather than substituting a
    // silently different force.
    constexpr double poleCancellationLimit=30.0;
    constexpr double retreatFraction=1.0e-7;
    // CREM_DEBUG_PROBES: log ALL six stencil probes, not only the ones that
    // trip the detector, so healthy and degenerate can be compared within the
    // SAME stencil -- same instant, same source geometry, same moment.
    // Without that control an E.B sign change at the bad probe means nothing,
    // because E.B = 0 surfaces are generic in a dipole field.
    static const bool traceProbes=std::getenv("CREM_DEBUG_PROBES")!=nullptr;
    double probeRatio[6]={},probeNMu[6]={},probeCosEB[6]={},probeEoverB[6]={};
    // The probes differ ONLY in which axis they step along, so if a dependence
    // survives it is on the step direction, not the line of sight.  Logged
    // signed: the + and - probes of one axis share |cos| with the pole axis,
    // and only one of them degenerates, so an unsigned angle cannot explain it.
    double probeDMu[6]={},probeDN[6]={};
    int probeIndex=0;
    const auto coupling=[&](const Vec3& point) {
        ElectromagneticField field=fieldFromOtherParticleAt(
            point,state.time,state,history,targetIsFirst,1.0e-5,
            matchedMomentSoftening());
        // Captured BEFORE the retreat, which overwrites the global with its
        // own (healthy) value -- reading it afterwards reports 1.09 for every
        // probe and hides exactly what is being looked for.
        const double asEvaluatedRatio=gPoleCancellationRatio;
        if(gPoleCancellationRatio>poleCancellationLimit) {
            const double flagged=gPoleCancellationRatio;
            const ElectromagneticField retreated=fieldFromOtherParticleAt(
                point,state.time,state,history,targetIsFirst,retreatFraction);
            if(gPoleCancellationRatio<=poleCancellationLimit)
                field=retreated;
            static const bool traceRetreat=
                std::getenv("CREM_DEBUG_RETREAT")!=nullptr;
            if(traceRetreat) {
                // Does the degenerate geometry have a signature?  Each pole
                // alone is photon-like by construction (B = n x E / c), so
                // E.B and E^2-c^2B^2 vanish for it identically; the SUM of
                // the two need not, which makes those invariants meaningful
                // here.  Printed with the line-of-sight/moment angle, so a
                // null-field or perpendicular-axis coincidence would show.
                const bool sourceIsFirst=!targetIsFirst;
                const Vec3 sourcePosition=sourceIsFirst?state.firstPosition
                                                       :state.secondPosition;
                const Vec3 sourceMoment=sourceIsFirst?state.firstDipole
                                                     :state.secondDipole;
                const Vec3 sight=point-sourcePosition;
                const double sightNorm=sight.norm();
                const double momentNorm=sourceMoment.norm();
                const double eNorm=field.electric.norm();
                const double bNorm=field.magnetic.norm();
                std::cerr<<"RETREAT ratio="<<flagged
                    <<" after="<<gPoleCancellationRatio
                    <<(gPoleCancellationRatio<=poleCancellationLimit
                        ?" accepted":" REJECTED")
                    <<" n.mu="<<(sightNorm>0.0&&momentNorm>0.0
                        ?dot(sight,sourceMoment)/(sightNorm*momentNorm):0.0)
                    <<" cos(E,B)="<<(eNorm>0.0&&bNorm>0.0
                        ?dot(field.electric,field.magnetic)/(eNorm*bNorm):0.0)
                    <<" |E|/c|B|="<<(bNorm>0.0?eNorm/(c*bNorm):0.0)<<'\n';
            }
        }
        // RELATIVE PLUS between the two channels.  This is the potential the
        // force is the gradient of, so the elementary requirement settles it:
        // a magnetic dipole feels grad(m.B) and an electric one grad(p.E), and
        // one potential generating both is m.B + p.E.
        //
        // It carried a MINUS until the boost measurement below caught it, on
        // the argument that m^{ab}F_{ab}/2 carries the same relative minus
        // that F^{ab}F_{ab}=2(B^2-E^2/c^2) does.  That argument is about a
        // particular index convention for the moment tensor, not about which
        // combination generates the force, and the minus does not survive a
        // direct test: transport the SAME proper moment into another frame the
        // way a moving particle actually carries it (rebuild the lab tensor at
        // that frame's own velocity, which is what synchronizeCovariantDipoles
        // does) and boost the field, then compare.  Measured over two boosts
        // and two particle velocities, m.B+p.E holds to 1.000000000 while
        // m.B-p.E comes out at -55.9, -71.9, -1.17e4 and -1.51e4.
        //
        // The old cross-check certified the minus because it transported the
        // tensor the OTHER way -- carrying an actively boosted tensor rather
        // than rebuilding it -- and the two errors cancelled exactly.  Either
        // pairing is self-consistent, which is why that probe reads 5.0e-6
        // for both and can never tell them apart; it has been corrected to
        // rebuild the tensor, and the discriminating check is
        // dipole-gradient-force-covariance, which boosts the assembled FORCE.
        if(traceProbes&&probeIndex<6) {
            const bool probeSourceIsFirst=!targetIsFirst;
            const Vec3 probeSource=probeSourceIsFirst?state.firstPosition
                                                     :state.secondPosition;
            const Vec3 probeMoment=probeSourceIsFirst?state.firstDipole
                                                     :state.secondDipole;
            const Vec3 sight=point-probeSource;
            const double sightNorm=sight.norm();
            const double momentNorm=probeMoment.norm();
            const double eNorm=field.electric.norm();
            const double bNorm=field.magnetic.norm();
            probeRatio[probeIndex]=asEvaluatedRatio;
            probeNMu[probeIndex]=(sightNorm>0.0&&momentNorm>0.0)
                ?dot(sight,probeMoment)/(sightNorm*momentNorm):0.0;
            probeCosEB[probeIndex]=(eNorm>0.0&&bNorm>0.0)
                ?dot(field.electric,field.magnetic)/(eNorm*bNorm):0.0;
            probeEoverB[probeIndex]=bNorm>0.0?eNorm/(c*bNorm):0.0;
            const Vec3 stepVector=point-targetPosition;
            const double stepNorm=stepVector.norm();
            probeDMu[probeIndex]=(stepNorm>0.0&&momentNorm>0.0)
                ?dot(stepVector,probeMoment)/(stepNorm*momentNorm):0.0;
            probeDN[probeIndex]=(stepNorm>0.0&&sightNorm>0.0)
                ?dot(stepVector,sight)/(stepNorm*sightNorm):0.0;
            ++probeIndex;
        }
        return dot(labMagneticDipole,field.magnetic)
              +dot(labElectricDipole,field.electric);
    };
    // All six probes read the source on ONE history segment -- the one that
    // holds the retarded time of the unshifted target -- continued across its
    // ends; see RetardedSegmentPin.  Probe retarded times lie within
    // gradientStep/(c(1-beta)) of the central one, so eight step lengths over
    // c cover them for any speed this integrator reaches.
    const RetardedSegmentPinGuard probeSegment(retardedSegmentPinAt(
        history,state,!targetIsFirst,targetPosition,state.time,
        8.0*gradientStep/c));
    Vec3 gradient;
    double probePlus[3]={},probeMinus[3]={};
    double channelPlusMagnetic[3]={},channelPlusElectric[3]={};
    double lwPlusMagnetic[3]={},dipPlusMagnetic[3]={};
    double channelMinusMagnetic[3]={},channelMinusElectric[3]={};
    for(int axis=0;axis<3;++axis) {
        Vec3 offset;
        if(axis==0) offset.x=gradientStep;
        else if(axis==1) offset.y=gradientStep;
        else offset.z=gradientStep;
        const double plus=coupling(targetPosition+offset);
        const double minus=coupling(targetPosition-offset);
        static const bool traceChannels=
            std::getenv("CREM_DEBUG_GRAD")!=nullptr;
        if(traceChannels) {
            const bool sourceIsFirst=!targetIsFirst;
            const ElectromagneticField lwPlus=lienardWiechertField(
                targetPosition+offset,state.time,history,state,sourceIsFirst,
                sourceIsFirst?firstCharge:secondCharge);
            const ElectromagneticField dipPlus=retardedMagneticDipoleField(
                targetPosition+offset,state.time,history,state,sourceIsFirst);
            lwPlusMagnetic[axis]=lwPlus.magnetic.norm();
            dipPlusMagnetic[axis]=dipPlus.magnetic.norm();
            const ElectromagneticField fp=fieldFromOtherParticleAt(
                targetPosition+offset,state.time,state,history,targetIsFirst);
            channelPlusMagnetic[axis]=dot(labMagneticDipole,fp.magnetic);
            channelPlusElectric[axis]=dot(labElectricDipole,fp.electric);
            const ElectromagneticField fm=fieldFromOtherParticleAt(
                targetPosition-offset,state.time,state,history,targetIsFirst);
            channelMinusMagnetic[axis]=dot(labMagneticDipole,fm.magnetic);
            channelMinusElectric[axis]=dot(labElectricDipole,fm.electric);
        }
        probePlus[axis]=plus;
        probeMinus[axis]=minus;
        const double derivative=(plus-minus)/(2.0*gradientStep);
        if(axis==0) gradient.x=derivative;
        else if(axis==1) gradient.y=derivative;
        else gradient.z=derivative;
    }
    if(traceProbes&&probeIndex==6) {
        bool anyDegenerate=false;
        for(int slot=0;slot<6;++slot)
            if(probeRatio[slot]>poleCancellationLimit) anyDegenerate=true;
        static thread_local int probeDumps=0;
        if(anyDegenerate&&probeDumps++<3) {
            std::cerr<<"PROBES r="<<separation(state)<<'\n';
            for(int slot=0;slot<6;++slot)
                std::cerr<<"  ax"<<(slot/2)<<(slot%2?'-':'+')
                    <<" ratio="<<probeRatio[slot]
                    <<" n.mu="<<probeNMu[slot]
                    <<" cos(E,B)="<<probeCosEB[slot]
                    <<" d.mu="<<probeDMu[slot]
                    <<" d.n="<<probeDN[slot]
                    <<(probeRatio[slot]>poleCancellationLimit
                        ?"   <-- DEGENERATE":"")<<'\n';
        }
    }
    // CREM_DEBUG_GRAD: dump the six stencil couplings when the assembled
    // gradient comes out anomalously large.  The gradient is a difference of
    // near-equal couplings, so either one probe point is an outlier (a
    // retarded-solve or formula-branch flip at that point) or all six are
    // sane and the blow-up is the cancellation losing its conditioning.
    // These two are distinguishable only by looking at the six.
    // CREM_DEBUG_EPS_SCAN: is a two-decade shrink of the pole separation
    // safe GENERALLY, or only at the one geometry measured so far?  The error
    // structure has a peak, so a fixed retreat target could land on an
    // equivalent peak elsewhere.  Sampled sparsely along a real collapse,
    // which sweeps four decades of separation on its own.
    static const bool scanEps=
        std::getenv("CREM_DEBUG_EPS_SCAN")!=nullptr;
    if(scanEps) {
        // Bucketed by DECADE of separation, not by call count: sampling every
        // Nth call spends the whole budget in the outer orbit, which is
        // exactly where nothing goes wrong.  At most 24 samples per decade.
        static thread_local int scanPerDecade[16]={};
        const double scanSeparation=separation(state);
        const int scanDecade=scanSeparation>0.0
            ?std::clamp(static_cast<int>(-std::log10(scanSeparation)),0,15):0;
        if(scanPerDecade[scanDecade]<24) {
            ++scanPerDecade[scanDecade];
            const bool sourceIsFirst=!targetIsFirst;
            static const double scanFractions[]={1.0e-4,3.0e-5,1.0e-5,3.0e-6,
                                                 1.0e-6,3.0e-7,1.0e-7,1.0e-8};
            // All THREE axes: the degeneracy that started this appeared on
            // y, and sampling only x misses it by construction.
            for(int scanAxis=0;scanAxis<3;++scanAxis) {
                Vec3 scanOffset;
                if(scanAxis==0) scanOffset.x=gradientStep;
                else if(scanAxis==1) scanOffset.y=gradientStep;
                else scanOffset.z=gradientStep;
                std::cerr<<"EPSSCAN r="<<separation(state);
                for(const double fraction:scanFractions) {
                    const ElectromagneticField swept=
                        retardedMagneticDipoleFieldExact(
                            targetPosition+scanOffset,state.time,history,state,
                            sourceIsFirst,fraction);
                    std::cerr<<" "<<swept.magnetic.norm();
                }
                std::cerr<<'\n';
            }
        }
    }
    static const bool traceAnomaly=
        std::getenv("CREM_DEBUG_GRAD")!=nullptr;
    if(traceAnomaly&&gradient.norm()/gamma(targetVelocity)>1.0e-3) {
        std::cerr<<std::setprecision(12)<<"GRAD anomaly |grad|="
            <<gradient.norm()/gamma(targetVelocity)
            <<" step="<<gradientStep<<" r="<<separation(state);
        for(int axis=0;axis<3;++axis)
            std::cerr<<"  ax"<<axis<<"[+"<<probePlus[axis]
                <<" -"<<probeMinus[axis]
                <<" d"<<(probePlus[axis]-probeMinus[axis])<<"]";
        // Does the degenerate probe have a GEOMETRIC signature?  Each pole
        // alone is photon-like by construction (B = n x E / c), so E.B and
        // E^2 - c^2 B^2 both vanish for it identically; the SUM of the two
        // poles need not, which makes those invariants meaningful for the
        // assembled dipole field.  Printed with the line-of-sight/moment
        // angle, so a null-field or perpendicular-axis coincidence at the
        // bad probe would show up directly.
        {
            const bool sourceIsFirst=!targetIsFirst;
            const Vec3 sourcePosition=sourceIsFirst?state.firstPosition
                                                   :state.secondPosition;
            const Vec3 sourceMoment=sourceIsFirst?state.firstDipole
                                                 :state.secondDipole;
            for(int axis=0;axis<3;++axis) {
                Vec3 offset;
                if(axis==0) offset.x=gradientStep;
                else if(axis==1) offset.y=gradientStep;
                else offset.z=gradientStep;
                const Vec3 probe=targetPosition+offset;
                const ElectromagneticField dipoleField=
                    retardedMagneticDipoleFieldExact(probe,state.time,
                        history,state,sourceIsFirst);
                const Vec3 sight=probe-sourcePosition;
                const double sightNorm=sight.norm();
                const double momentNorm=sourceMoment.norm();
                const double electricNorm=dipoleField.electric.norm();
                const double magneticNorm=dipoleField.magnetic.norm();
                const double eDotB=dot(dipoleField.electric,
                                       dipoleField.magnetic);
                const double invariant=electricNorm*electricNorm
                    -c*c*magneticNorm*magneticNorm;
                std::cerr<<"\n  geo"<<axis
                    <<" n.mu="<<(sightNorm>0.0&&momentNorm>0.0
                        ?dot(sight,sourceMoment)/(sightNorm*momentNorm):0.0)
                    <<" cos(E,B)="<<(electricNorm>0.0&&magneticNorm>0.0
                        ?eDotB/(electricNorm*magneticNorm):0.0)
                    <<" (E2-c2B2)/(E2+c2B2)="
                    <<invariant/std::max(electricNorm*electricNorm
                        +c*c*magneticNorm*magneticNorm,1.0e-300)
                    <<" |E|/c|B|="<<(magneticNorm>0.0
                        ?electricNorm/(c*magneticNorm):0.0);
            }
        }
        // eps SWEEP at each probe.  The two-charge construction is a device
        // for the eps->0 limit, so the question that decides what to do about
        // a degenerate probe is whether shrinking eps recovers the answer
        // there (then eps is simply too coarse for that geometry) or whether
        // it stays erratic (then the degeneracy is directional and no scalar
        // eps fixes it).  Measured here rather than assumed.
        {
            const bool sourceIsFirst=!targetIsFirst;
            static const double fractions[]={1.0e-4,3.0e-5,1.0e-5,
                                             3.0e-6,1.0e-6,1.0e-7};
            for(int axis=0;axis<3;++axis) {
                Vec3 offset;
                if(axis==0) offset.x=gradientStep;
                else if(axis==1) offset.y=gradientStep;
                else offset.z=gradientStep;
                std::cerr<<"\n  eps"<<axis;
                for(const double fraction:fractions) {
                    const ElectromagneticField swept=
                        retardedMagneticDipoleFieldExact(
                            targetPosition+offset,state.time,history,state,
                            sourceIsFirst,fraction);
                    std::cerr<<" "<<swept.magnetic.norm();
                }
            }
        }
        // The two channels separately: U is their SUM, so if each is large
        // and they nearly cancel, the outlier is a cancellation failure
        // rather than a bad field evaluation.
        for(int axis=0;axis<3;++axis)
            std::cerr<<"\n  src"<<axis<<" +[|B_LW| "<<lwPlusMagnetic[axis]
                <<" |B_dip| "<<dipPlusMagnetic[axis]<<"]";
        for(int axis=0;axis<3;++axis)
            std::cerr<<"\n  ch"<<axis
                <<" +[muB "<<channelPlusMagnetic[axis]
                <<" pE "<<channelPlusElectric[axis]<<"]"
                <<" -[muB "<<channelMinusMagnetic[axis]
                <<" pE "<<channelMinusElectric[axis]<<"]";
        // The six field evaluations behind those couplings, most recent last.
        std::cerr<<"\n  traces";
        for(int slot=0;slot<8;++slot) {
            const TwoChargeLimitTrace& t=
                gTwoChargeTrace[(gTwoChargeTraceCount+slot)&7];
            std::cerr<<" [it"<<t.centralIterations
                <<" |pole|"<<t.poleMagnitude
                <<" |sum|"<<t.sumMagnitude
                <<" keep"<<(t.poleMagnitude>0.0
                    ?t.sumMagnitude/t.poleMagnitude:0.0)
                <<(t.earlyReturn?" EARLY":"")<<"]";
        }
        std::cerr<<std::setprecision(6)<<'\n';
    }
    // Fixed-mass completion of the model's rest-frame force grad(U).
    // With metric (+---), partial^a U = (partial_t U/c, -grad U):
    //   f^a = -(eta^{ab} - u^a u^b/c^2) partial_b U.
    // Thus u.f = 0, f_rest = (0, grad_rest U), and the lab three-force is
    //   F = grad U/gamma + gamma v (partial_t U + v.grad U)/c^2.
    // The PLUS sign follows from the metric and the positive rest gradient.
    // grad U/gamma alone fails for a target moving along the gradient; the
    // earlier transverse-orbit boost probe had almost zero DU/Dt and could
    // not detect that missing component (see the 2026-09-09 audit).
    const double targetGamma=gamma(targetVelocity);
    if(targetVelocity.squaredNorm()==0.0) return gradient;
    const double materialRate=dipoleCouplingMaterialRate(
        state,history,targetIsFirst);
    return gradient/targetGamma
        +targetVelocity*(targetGamma*materialRate/(c*c));
}

// Weight of the RETARDED dipole sector at this state: 1 by default, and the
// smoothstep that CREM_DIPOLE_INSTANTANEOUS_BELOW blends the retarded sector
// with the instantaneous one below N..2N floors.  Shared by the force blend
// and by the energy ledger, whose charge-dipole term belongs to the retarded
// sector only (see chargeDipoleInteractionEnergy).
inline double retardedDipoleSectorWeight(const State& s) {
    static const double instantaneousDipoleFloors=[] {
        const char* text=std::getenv("CREM_DIPOLE_INSTANTANEOUS_BELOW");
        const double value=text?std::atof(text):0.0;
        return (std::isfinite(value)&&value>0.0)?value:0.0;
    }();
    if(!gDipoleForceEnabled||!(instantaneousDipoleFloors>0.0)
       ||!(separationFloor()>0.0)) return 1.0;
    const double inner=instantaneousDipoleFloors*separationFloor();
    const double transition=std::clamp(
        (separation(s)-inner)/inner,0.0,1.0);
    return transition*transition*(3.0-2.0*transition);
}

inline MutualForces retardedExternalForces(const State& s,
                                    const StateHistory& history) {
    const ElectromagneticField secondField = lienardWiechertField(
        s.firstPosition, s.time, history, s, false, secondCharge);
    const ElectromagneticField firstField = lienardWiechertField(
        s.secondPosition, s.time, history, s, true, firstCharge);
    // CREM_NO_DIPOLE_FORCE also gates the RETARDED dipole field, not just
    // the instantaneous term in mutualForces.  Ablating only the latter was
    // an incomplete experiment: this is the path the integrator actually
    // takes (crem_engine.hpp calls retardedExternalForces every step), so
    // leaving it in meant the moments still steered the orbit while the
    // ablation claimed they did not.
    const ElectromagneticField secondDipoleField=gDipoleForceEnabled
        ?retardedMagneticDipoleField(s.firstPosition,s.time,history,s,false)
        :ElectromagneticField{};
    const ElectromagneticField firstDipoleField=gDipoleForceEnabled
        ?retardedMagneticDipoleField(s.secondPosition,s.time,history,s,true)
        :ElectromagneticField{};
    if(std::getenv("POSITRONIUM_DEBUG_FIELDS")) {
        static int debugSamples=0;
        if(debugSamples++<24) {
            const ElectromagneticField secondMagneticLow=
                retardedMagneticDipoleFieldLowVelocity(
                    s.firstPosition,s.time,history,s,false,true);
            std::cerr<<"DIPOLE_DEBUG t="<<s.time
                <<" r="<<separation(s)
                <<" beta="<<s.firstVelocity.norm()/c<<"/"
                <<s.secondVelocity.norm()/c
                <<" Bexact/low="<<secondDipoleField.magnetic.norm()<<"/"
                <<secondMagneticLow.magnetic.norm()
                <<" Eexact/low="<<secondDipoleField.electric.norm()<<"/"
                <<secondMagneticLow.electric.norm()<<'\n';
        }
    }
    MutualForces chargeCharge{
        lorentzForce(firstCharge,s.firstVelocity,
            {secondField.electric+secondDipoleField.electric,
             secondField.magnetic+secondDipoleField.magnetic}),
        lorentzForce(secondCharge,s.secondVelocity,
            {firstField.electric+firstDipoleField.electric,
             firstField.magnetic+firstDipoleField.magnetic})};

    // Third and last path by which the moments steer the orbit, gated with
    // the other two so that CREM_NO_DIPOLE_FORCE means what it says.
    MutualForces tensorGradient=gDipoleForceEnabled
        ?MutualForces{covariantDipoleGradientForce(s,history,true),
                      covariantDipoleGradientForce(s,history,false)}
        :MutualForces{};
    // CREM_DIPOLE_INSTANTANEOUS_BELOW=N (in separation floors; unset or 0 =
    // off, the default, which leaves every force above bit-identical).
    // Below N floors the WHOLE dipole sector is the instantaneous one of
    // allExternalForces -- regularizedDipoleForce plus chargeDipoleForces,
    // the forces whose energies conservativeParticleEnergy counts, both
    // carried back through the separation clamp -- and
    // between N and 2N floors it blends into the retarded sector with a
    // smoothstep in the present separation.  The charge-charge
    // Lienard-Wiechert force stays retarded throughout.
    //
    // Why an option exists (audit section 58): in a deep tilted para passage
    // the retarded gradient force follows the partner's RETARDED position and
    // moment, whose direction turns by up to 0.56 rad during r/c at
    // beta 0.5-0.8, while the mechanical energy is instantaneous.  Its work
    // swung by -263 and +302 k/r0 across 1-4 floors against a dipole-energy
    // change of +-30 k/r0, and that is where the 40-400 k/r0 energy errors of
    // every deep run came from; the integrator's own share was below 7 k/r0.
    // Whether that exchange is physics or an artefact of point dipoles at
    // these radii is a modelling choice, so it is a switch, not a fix.
    {
        const double retardedWeight=retardedDipoleSectorWeight(s);
        if(retardedWeight<1.0) {
            const MutualForces chargeOnly{
                lorentzForce(firstCharge,s.firstVelocity,secondField),
                lorentzForce(secondCharge,s.secondVelocity,firstField)};
            const Vec3 dipoleOnFirst=pairDipoleForce(
                s.firstPosition-s.secondPosition,s.firstDipole,s.secondDipole);
            // The Thomas-BMT moment derivatives inside chargeDipoleForces
            // read the fields from the REAL retarded history.
            // allExternalForces passes a one-state history instead, which
            // extrapolates the partner back by r/c with the acceleration
            // stored in the state -- and inside a step that is the previous
            // step's value.  Measured at a switched tilted para failure at
            // 0.40 r*: with the one-state history the step-doubling velocity
            // error halved per halving of dt (a local error of order dt,
            // 3.1e-9 at the failing step); removing chargeDipoleForces, or
            // giving it this history, restored 4x per halving (3.7e-11).
            // Neither the dipole-dipole force nor the precession mattered.
            const MutualForces mixed=chargeDipoleForces(s,history);
            const MutualForces instantaneous{dipoleOnFirst+mixed.first,
                                             mixed.second-dipoleOnFirst};
            const double w=retardedWeight;
            chargeCharge={
                chargeOnly.first+(chargeCharge.first-chargeOnly.first)*w,
                chargeOnly.second+(chargeCharge.second-chargeOnly.second)*w};
            tensorGradient={
                tensorGradient.first*w+instantaneous.first*(1.0-w),
                tensorGradient.second*w+instantaneous.second*(1.0-w)};
        }
    }
    // CREM_RETARDED_THOMAS: the orbital back-reaction of Thomas precession,
    // which allExternalForces carries inside chargeDipoleForces and this sum
    // has never had (audit sections 105g, 120).  Off by default, so
    // production is bit-identical without it; the moment derivatives come
    // from the same retarded history the rest of this sum uses.
    // The rate of the loops' own hidden momentum (audits 122-124), ON by
    // decision after audit 123 measured it: a current loop carries
    // (mu x E)/c^2 outside gamma m v, and the retarded sum has to carry its
    // rate or it is not Newton's law for the mechanical momentum.  It removes
    // 46% of the sum's measured imbalance at every speed and leaves the
    // rigidly drifting pair exactly balanced.  What it does NOT do is improve
    // a deep passage, because the other 54% is still unsourced (123c): the
    // orbit at 1 r* ends with twice the momentum.  CREM_NO_HIDDEN_MOMENTUM_
    // FORCE removes it for comparison.
    static const bool hiddenMomentumForce=
        std::getenv("CREM_NO_HIDDEN_MOMENTUM_FORCE")==nullptr;
    MutualForces hiddenRate;
    if(hiddenMomentumForce&&gDipoleForceEnabled)
        hiddenRate={hiddenMomentumRateForce(s,history,true),
                    hiddenMomentumRateForce(s,history,false)};
    static const bool retardedThomas=
        std::getenv("CREM_RETARDED_THOMAS")!=nullptr;
    MutualForces thomas;
    if(retardedThomas)
        thomas=thomasBackReactionForces(s,
            thomasBmtDipoleDerivatives(s,history));
    // Same uniform external field as in the instantaneous sum.  It is not
    // retarded because it is not sourced by either particle.
    MutualForces externalField{
        lorentzForce(firstCharge,s.firstVelocity,{{},gExternalMagneticField}),
        lorentzForce(secondCharge,s.secondVelocity,{{},gExternalMagneticField})};
    if(gZeroPointField.active()) {
        Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
        const ZeroPointDrive drive=zeroPointDriveFrequency(s);
        const double orbitalFrequency=drive.frequency;
        const double orbitalFrequencyDerivative=drive.derivative;
        gZeroPointField.sample(s.firstPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,firstElectric,firstMagnetic);
        gZeroPointField.sample(s.secondPosition,orbitalFrequency,
                               orbitalFrequencyDerivative,
                               s.zeroPointPhase,secondElectric,secondMagnetic);
        externalField.first=externalField.first+lorentzForce(firstCharge,
            s.firstVelocity,{firstElectric,firstMagnetic})
            +gZeroPointField.gradientForce(s.firstPosition,orbitalFrequency,
                s.zeroPointPhase,s.firstDipole);
        externalField.second=externalField.second+lorentzForce(secondCharge,
            s.secondVelocity,{secondElectric,secondMagnetic})
            +gZeroPointField.gradientForce(s.secondPosition,orbitalFrequency,
                s.zeroPointPhase,s.secondDipole);
    }
    return {chargeCharge.first+tensorGradient.first+externalField.first
                +thomas.first+hiddenRate.first,
            chargeCharge.second+tensorGradient.second+externalField.second
                +thomas.second+hiddenRate.second};
}

struct CanonicalMomenta { Vec3 first, second; };

// Canonical momenta obtained from the conservative low-velocity action.  The
// q-mu contributions occur with opposite signs on the charge and on the
// dipole carrier; the Darwin terms do not cancel particle by particle.
inline CanonicalMomenta canonicalMomenta(const State& s) {
    const PairGeometry geometry = clampedPairGeometry(s);
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
inline Vec3 noetherMomentum(const CanonicalMomenta& canonical) {
    return canonical.first + canonical.second;
}

inline Vec3 noetherMomentum(const State& s) {
    return noetherMomentum(canonicalMomenta(s));
}

inline double canonicalMomentumScale(const CanonicalMomenta& canonical) {
    return canonical.first.norm() + canonical.second.norm();
}

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
// Production callers already hold the canonical momenta and use the overload
// above; only the Yee-coupling tests start from a bare state.
inline double canonicalMomentumScale(const State& s) {
    return canonicalMomentumScale(canonicalMomenta(s));
}
#endif

inline Vec3 noetherAngularMomentum(const State& s,
                            const CanonicalMomenta& canonical) {
    const double firstGyromagneticRatio = firstGyromagneticRatioOf();
    const double secondGyromagneticRatio = secondGyromagneticRatioOf();
    const Vec3 orbital = cross(s.firstPosition, canonical.first)
                       + cross(s.secondPosition, canonical.second);
    const Vec3 intrinsic = s.firstDipole / firstGyromagneticRatio
                         + s.secondDipole / secondGyromagneticRatio;
    return orbital + intrinsic;
}

inline Vec3 noetherAngularMomentum(const State& s) {
    return noetherAngularMomentum(s, canonicalMomenta(s));
}

// The Darwin energy: exact through O(v^2/c^2) and no further.  It is the
// potential of the Coulomb+Darwin force (allExternalForces), NOT of the
// retarded Lienard-Wiechert force production actually integrates, and the
// gap between the two is the reason this quantity has a floor no tolerance
// can lower.  Measured over one orbit at a_Ps with radiation reaction off:
// integrating the matched Coulomb+Darwin force, the drift falls 5.13e-7 ->
// 1.28e-7 -> 6.42e-8 -> 3.21e-8 as the tolerance goes 1e-6 -> 1e-9, cleanly
// converging; integrating the retarded force it sits at 1.3-2.3e-6 and does
// not converge at all.  Swept over beta (a_Ps down to a_Ps/16) that floor
// scales as beta^2.94 +- 0.15 with floor/beta^3 constant to 7% at 3.8-4.1 --
// i.e. it is the O(beta^3) retardation term, the first order Darwin omits
// and the order at which the two-body interaction stops being derivable
// from a potential at all.
//
// Consequences for callers: differences of this function are meaningful only
// against a reference carrying the same truncation at the same phase, which
// is exactly what crem_collapse.hpp's background run is for.  See the README
// section on the deterministic threshold in the energy balance.
// CHARGE-DIPOLE (SPIN-ORBIT) INTERACTION ENERGY -- the motional electric
// dipole of each moving moment in the Coulomb field of the other charge,
//
//     U = -( p_1 . E_2(r_1) + p_2 . E_1(r_2) ),   p_i = lab electric dipole,
//
// p_i being the tensor partner lorentzBoostDipole gives the proper magnetic
// moment at the particle's OWN velocity (p = gamma v x mu/c^2 at leading
// order), with E the Plummer (floor on) or regularized field of the charge.
//
// WHICH TERM, AND FOR WHICH FORCES.  Measured along the model's own flow
// (audit section 82: x' = v, p' = F, mu' = BMT; dE/dt of the full ledger by
// central difference, spin sector = with moments minus without):
//
//                               no term     old -mu.B(v_other)   -p.E
//   e+e- retarded, ortho 0.1    3.4e-3 W    3.9e-6 W             3.9e-6 W
//   p+e- retarded, para 0.1     2.7e-2 W    2.7e-2 W             9.8e-6 W
//   e+e- instantaneous, ortho   9.3e-7 W    3.4e-3 W             3.4e-3 W
//   p+e- instantaneous, para    1.4e-8 W    2.7e-5 W             2.7e-2 W
//
// The retarded assembly couples each moment through U = mu.B + p.E, so its
// ledger needs -p.E.  The old form, each moment in the Biot-Savart field of
// the OTHER charge's motion, is identical for equal masses (v_2 = -v_1) and
// wrong for any other pair: for p+e- it missed the electron's own motional
// dipole entirely.  The instantaneous sum derives from Lagrangians linear in
// the velocities (q (v_q - v_mu).A and the Thomas term), which contribute
// nothing to the energy function, so it needs NO term.  conservative-
// ParticleEnergy therefore takes the retarded-sector weight: 1 for the
// retarded assembly, 0 for allExternalForces, the smoothstep in between
// under CREM_DIPOLE_INSTANTANEOUS_BELOW.
//
// History: introduced by 8c85397 as -mu.B(v_other).  Its supporting
// measurement, "the balance residual halves", was an artefact of an
// unsynchronized validation start state (section 82): from a synchronized
// start the residual is 1e-8 with or without the term over that window, and
// the flow check above is what actually tests it.  For e+e- in the centre-of-
// mass frame the two forms agree, so production is unchanged.
inline double chargeDipoleInteractionEnergy(const State& state) {
    const Vec3 firstMinusSecond = state.firstPosition - state.secondPosition;
    const double distance = firstMinusSecond.norm();
    // This term is the energy of the force audit 126 rebuilt -- the motional
    // dipole in the partner's charge field -- so in principle it should carry
    // the same softening length, or the force is not minus the gradient of
    // the ledger.  Measured (audit 127): giving it that length moves the
    // one-orbit ledger drift from -1.084 to -1.129 |U_dd| at 1 r* and from
    // -0.230 to -0.246 at 0.5 r*, i.e. 4-7% the WRONG way, and leaves the
    // momentum untouched because no force reads this term.  So the residual
    // drift is not this mismatch, and the ledger keeps its own floor;
    // CREM_MATCHED_ENERGY switches it over for comparison.
    static const bool matchedEnergy=
        std::getenv("CREM_MATCHED_ENERGY")!=nullptr;
    const double floor = (matchedEnergy&&matchedMomentSoftening()>0.0)
        ? matchedMomentSoftening() : separationFloor();
    if(!(distance > 0.0) && !(floor > 0.0)) return 0.0;
    // The Plummer charge's field (section 66): k q d / rho^3.
    const double radialFactor = floor > 0.0
        ? 1.0/std::pow(distance*distance + floor*floor, 1.5)
        : magneticRadialProfile(distance).vectorPotentialFactor;
    const Vec3 fieldAtFirst = firstMinusSecond
        * (coulomb * secondCharge * radialFactor);
    const Vec3 fieldAtSecond = firstMinusSecond
        * (-coulomb * firstCharge * radialFactor);
    // p is rebuilt from the moments and velocities rather than read from
    // state.firstElectricDipole: callers (validation probes, the secular
    // turning point) build states that set a magnetic moment without ever
    // synchronizing its electric partner, and a stale zero there reads as a
    // jump of the whole term on the first engine step (it moved the shared-
    // engine raw energy residual from 1.8e-6 to 8.4e-3).  synchronized is the
    // same boost synchronizeCovariantDipoles applies, so on a synchronized
    // state the result is identical.
    State synchronized = state;
    synchronizeCovariantDipoles(synchronized);
    const double energy = -dot(synchronized.firstElectricDipole, fieldAtFirst)
                          -dot(synchronized.secondElectricDipole, fieldAtSecond);
    return std::isfinite(energy) ? energy : 0.0;
}

// chargeDipoleWeight: the retarded-sector weight of the charge-dipole term,
// see chargeDipoleInteractionEnergy.  The default 1 is the production
// (retarded) ledger; allExternalForces-driven callers pass 0.
inline double conservativeParticleEnergy(const State& state,
                                         double chargeDipoleWeight=1.0) {
    const PairGeometry geometry=clampedPairGeometry(state);
    const double kinetic=kineticEnergy(state.firstVelocity,firstMass)
        +kineticEnergy(state.secondVelocity,secondMass);
    const double coulombPotential=-pairCoulombStrength
        *geometry.inverseDistance;
    const double dipolePotential=pairDipoleInteractionEnergy(
        state.firstPosition-state.secondPosition,
        state.firstDipole,state.secondDipole);
    return kinetic+coulombPotential+dipolePotential
        +chargeDipoleWeight*chargeDipoleInteractionEnergy(state)
        +darwinInteractionEnergy(state)+state.dipoleConstraintEnergy;
}

// THE CONSERVED ENERGY, as opposed to the radial one (audit section 145).
// conservativeParticleEnergy above sums every interaction the model carries,
// including the charge-dipole term.  That term is -p.E with p = (v x mu)/c^2,
// homogeneous of degree one in both velocities, so it cancels in the Legendre
// transform H = sum v.(dL/dv) - L and does NOT belong to the conserved
// energy: it generates a force -- the spin-orbit force, 88% of the moment
// sector by audit 105 -- and that force does no work.
//
// Measured before this was split out: with instantaneous forces, where an
// instantaneous ledger could be exact, including the term inflates the
// ledger's range over one orbit by a factor of 24 to 104 across four
// configurations, while the other two candidate corrections (the Darwin sign,
// the proper moments) make it worse.  With the model's own retarded forces
// the term is masked by the retardation, which contributes about 1.13 |U_dd|
// per orbit of its own and with the opposite sign.
//
// The same energy IS the right radial potential at FIXED angular momentum,
// where it is a function of r and L alone; that is how dipoleAwarePeriapsis
// uses chargeDipoleInteractionEnergy, and nothing here changes it.  The two
// quantities were one function and are now two.
inline double conservedParticleEnergy(const State& state) {
    return conservativeParticleEnergy(state,0.0);
}

inline FourVector fourVelocity(const Vec3& velocity) {
    const double relativisticGamma = gamma(velocity);
    return {relativisticGamma * c, velocity * relativisticGamma};
}

inline FourVector fourForce(const Vec3& velocity, const Vec3& force) {
    const double relativisticGamma = gamma(velocity);
    return {relativisticGamma * dot(force, velocity) / c,
            force * relativisticGamma};
}

// chargeCoupling is q_i^2 for a self term and q_i q_j for the mutual one; it
// used to be a hard-wired e^2, which made both correct only for unit charges.
inline Vec3 reducedOrderSelfForce(const Vec3& velocity,
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

inline MutualForces individualLandauLifshitzSelfForces(
    const State& s, const MutualForces& external,
    const StateHistory& history,bool includeMutual) {
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
    // the self forces above reproduce only the first two terms.  These explicit
    // mutual terms are required only when the external force is the
    // non-retarded Coulomb-Darwin approximation.  With full Lienard-Wiechert
    // mutual fields the partner's radiation field is already present and
    // includeMutual must be false; the independent two-period balance measures
    // 0.88% residual for Darwin+coherent, 5.94% for retarded+self-only, and
    // 51.0% for the incorrect retarded+full-coherent combination.
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

    return includeMutual
        ?MutualForces{firstSelf+firstMutual,secondSelf+secondMutual}
        :MutualForces{firstSelf,secondSelf};
}

// Explicit reversible near-field energy paired with the coherent charge
// radiation reaction.  For d=sum(q_i r_i), the reduced-order reaction work is
//
//   Wdot = C d''' . d' = C d/dt(d'' . d') - C |d''|^2,
//
// so E_S=-C d'' . d' is the Schott energy that makes
// dE_mech/dt + P_far + dE_S/dt vanish in the electric-dipole limit.  This is
// evaluated directly from the external-force accelerations; it does not read
// boundFieldEnergy or any residual accumulator and can therefore serve as an
// independent balance diagnostic.
[[maybe_unused]] inline double explicitChargeSchottEnergy(
    const State& state,const MutualForces& externalForces) {
    const Vec3 firstAcceleration=relativisticAcceleration(
        state.firstVelocity,externalForces.first,firstMass);
    const Vec3 secondAcceleration=relativisticAcceleration(
        state.secondVelocity,externalForces.second,secondMass);
    const Vec3 dipoleVelocity=state.firstVelocity*firstCharge
        +state.secondVelocity*secondCharge;
    const Vec3 dipoleAcceleration=firstAcceleration*firstCharge
        +secondAcceleration*secondCharge;
    return -dot(dipoleAcceleration,dipoleVelocity)
        /(6.0*pi*epsilon0*c*c*c);
}

// Relativistic predictor-corrector update with mutually retarded fields and
// an individual, reduced-order Landau-Lifshitz self-force for each charge.
inline void integrateElectrodynamicStep(State& s, double dt,
                                 const StateHistory& history,
                                 bool computeOutwardFlux=true,
                                 ChargeRadiationReactionModel reactionModel=
                                    ChargeRadiationReactionModel::individualLandauLifshitz,
                                 bool useRetardedExternalForces=true) {
    const State balanceStart=s;
    // The M1 sector has two back-reactions, and they are NOT the same
    // decision, which is why they now sit behind two flags rather than one.
    //
    // The ENERGY drain (dipoleConstraintEnergy) is switched off for
    // stochasticElectricDipole, because every radiation channel leaves as
    // discrete quanta there (crem_trajectory.hpp banks the TOTAL E1+M1+E2
    // power as one photon stream).  Leaving it on would remove the magnetic
    // dipole's radiated energy twice: once continuously here and once as
    // photon energy there.
    //
    // The ORIENTATION back-reaction (the reaction torque) has no such
    // duplicate.  applyStochasticDipolePhoton rebuilds the four-velocities
    // and nothing else -- the photon carries linear recoil only, and leaves
    // the moment that radiated it pointing exactly where it did before.
    // Gating the torque alongside the drain therefore did not prevent a
    // double count; it silently deleted the one back-reaction the quantized
    // path had no other route to.  Measured on the continuous side, over 96
    // random orientation pairs integrated along the whole inspiral, the
    // omission was worth 2.33e-2 +- 0.24e-2 rad, i.e. 1.33 +- 0.14 degrees
    // of moment rotation per collapse.  The torque now runs in the quantized
    // mode too.
    //
    // No energy is double counted by that: the torque's own energetic cost
    // is the change it makes to the dipole-dipole interaction energy, which
    // conservativeParticleEnergy already carries through
    // regularizedDipoleInteractionEnergy.  The photon accounts for the
    // radiated far-zone energy; the torque accounts for the orientation.
    //
    // BOTH are off for disabled, whose whole point is that NOTHING drags the
    // orbit -- it is the mechanical reference crem_collapse.hpp integrates
    // as the background whose delta is subtracted from the real run's.  That
    // half used to be missing, and the omission was not academic: the
    // background then drained the same M1 energy as the run it was
    // subtracted from, so in the continuous models the M1 sink cancelled
    // itself out (measured: 4.0e-5 of it survived the subtraction), while in
    // the production stochastic mode the background carried a sink the real
    // run did not have at all.  The charge sector has always excluded both
    // disabled and stochastic, at the corresponding gate inside
    // particleMultipoleRadiation; the drain flag now matches it.
    const bool quantizedRadiation=
        reactionModel==ChargeRadiationReactionModel::stochasticElectricDipole
        ||reactionModel==ChargeRadiationReactionModel::disabled;
    // The ledger must match the force model the step integrates.
    const auto ledgerWeight=[&](const State& state) {
        return useRetardedExternalForces
            ?retardedDipoleSectorWeight(state):0.0;
    };
    const double initialMechanicalEnergy=conservativeParticleEnergy(
        balanceStart,ledgerWeight(balanceStart));
    const CanonicalMomenta initialCanonical=canonicalMomenta(balanceStart);
    const Vec3 initialMechanicalMomentum=noetherMomentum(initialCanonical);
    const Vec3 initialMechanicalAngularMomentum=
        noetherAngularMomentum(balanceStart,initialCanonical);
    // THIS IS THE OPERATOR SPLITTING, and it is worth naming because it does
    // not look like one: the half-precession here, the mechanical step below,
    // and the second half-precession before the accelerations are recomputed
    // are a Strang splitting of the spin sector against the orbital one.  The
    // spin half-steps are not integrated, they are EXACT -- see
    // advanceThomasBmtDipole, a closed-form Rodrigues rotation about a B_eff
    // held fixed across the half-step -- so the spin sub-step carries no
    // step-size restriction at all.  Measured: at omega*dt from 1.3e-11 up to
    // 1.3e+19 radians it stays norm-conserving to 1.1e-16 and two half-steps
    // reproduce one full step to 2.4e-16.
    //
    // So the spin sector is NOT what drives dt down near the barrier, and it
    // is not the stiff one either.  omega_spin/omega_orbit RISES as the orbit
    // tightens but SATURATES well below one, and then turns over:
    //   r/r*     10     5      3      2     1.5     1.2      1     0.5    0.3
    //   ratio  .0035  .0097  .019   .028   .036   .042   .047   .054   .054
    // At r* the precession is nineteen times SLOWER than the orbit.  It never
    // becomes the fast mode, so there is no stiffness here for an implicit or
    // exponential integrator to buy back.
    //
    // What does drive dt down is the DIPOLE FORCE sector, and it is a C0
    // discontinuity rather than stiffness -- no integrator order helps, this
    // splitting included.  Ablating that sector at r=1.2 r* restores textbook
    // second order (3.98-4.03x per halving, out to 8e-15); leaving it in
    // gives 6.55x, 1.08x, 1.94x on the last three halvings.  Audit section
    // 149 pins the mechanism: a 4% jump in the retarded magnetic-dipole
    // field's ELECTRIC component, turned by the hidden-momentum term's
    // finite-difference d(mu x E)/dt into a force 15.4 times Coulomb.
    applyDipolePrecession(s, 0.5 * dt, history);
    if(std::getenv("POSITRONIUM_DEBUG_DIPOLE")&&!isFinite(s))
        std::cerr<<"STEP_DEBUG nonfinite-after-precession t="<<s.time
            <<" dt="<<dt<<'\n';
    MutualForces forces = useRetardedExternalForces
        ?retardedExternalForces(s,history):allExternalForces(s);
    if(std::getenv("POSITRONIUM_DEBUG_DIPOLE")
       &&(!isFinite(forces.first)||!isFinite(forces.second)))
        std::cerr<<"STEP_DEBUG nonfinite-forces t="<<s.time
            <<" dt="<<dt<<'\n';
    const ParticleMultipoleRadiation radiation =
        particleMultipoleRadiation(
            s,forces,history,computeOutwardFlux,reactionModel,
            useRetardedExternalForces);
    if (!finiteRadiationResponse(radiation)) {
        if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
            std::cerr<<"STEP_DEBUG nonfinite-radiation t="<<s.time
                <<" dt="<<dt<<" flux="<<radiation.outwardFlux.energy
                <<" lead="<<radiation.leadingElectricDipolePower
                <<" m="<<radiation.magneticDipoleFlux.energy
                <<" ll="<<radiation.landauLifshitzValidity<<'\n';
        s.time = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    applyDipoleRadiationTorqueForModel(
        s,radiation,0.5*dt,reactionModel);
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

    const MutualForces trialForces = useRetardedExternalForces
        ?retardedExternalForces(trial,history):allExternalForces(trial);
    if(std::getenv("POSITRONIUM_DEBUG_DIPOLE")
       &&(!isFinite(trialForces.first)||!isFinite(trialForces.second)))
        std::cerr<<"STEP_DEBUG nonfinite-trial-forces t="<<trial.time
            <<" dt="<<dt<<'\n';
    const ParticleMultipoleRadiation trialRadiation =
        particleMultipoleRadiation(
            trial,trialForces,history,computeOutwardFlux,reactionModel,
            useRetardedExternalForces);
    if (!finiteRadiationResponse(trialRadiation)) {
        if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
            std::cerr<<"STEP_DEBUG nonfinite-trial-radiation t="<<trial.time
                <<" dt="<<dt<<" flux="<<trialRadiation.outwardFlux.energy
                <<" lead="<<trialRadiation.leadingElectricDipolePower
                <<" m="<<trialRadiation.magneticDipoleFlux.energy
                <<" ll="<<trialRadiation.landauLifshitzValidity<<'\n';
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
    applyDipoleRadiationTorqueForModel(
        trial,trialRadiation,0.5*dt,reactionModel);
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
    // The committed-state RATE above is kept for exactly that reason; what
    // changed since is only the QUADRATURE applied to it, which used to be a
    // left rule and is now trapezoidal -- see below.  None of it feeds back
    // into the trajectory.
    // Trapezoidal weight for the radiation bookkeeping.  The rate is sampled
    // at the step START, so weighting every sample by its own dt is a LEFT
    // Riemann sum: first order, and near collapse the power roughly doubles
    // per step, so that error is about half the power itself.  Measured on
    // the collapse trajectory against the identity
    // W = -dE_S^coh - P_E1 dt, sampled every fifteenth step: the left rule
    // leaves -0.298, -0.435, -0.492, -0.503 of P dt as r/r_coll runs
    // 104 -> 4.9, while matched trapezoidal quadrature leaves +0.026,
    // +0.041, +0.052, +0.064.  e+e- and mu+mu- agree to four digits at equal
    // r/r_coll, so this is quadrature, not a pair-dependent effect.
    //
    // Summing f_k*(dt_{k-1}+dt_k)/2 IS the trapezoid rule regrouped by
    // sample, and it needs only the previous step length -- no second force
    // evaluation, so the cost is one double.  The earlier attempt used the
    // trial state's rate as the right endpoint and under-reported by 34%
    // because that reconstruction is inconsistent; this form never evaluates
    // anything off the committed trajectory.
    //
    // Concretely, after step k the accumulator holds every COMPLETED interval
    // as an exact trapezoid plus a provisional left-rule value for the
    // current one, which the next step corrects:
    //
    //     increment = f_k*dt_k + (f_k - f_{k-1})*dt_{k-1}/2.
    //
    // Weighting each sample by (dt_{k-1}+dt_k)/2 instead is the same rule
    // regrouped, but it drops the final half interval outright; that cost
    // 0.34% on the suite's Larmor accumulation check (0.99998 -> 0.99656),
    // because the deficit there is half a step's VALUE.  The form above
    // leaves only half a step's CHANGE, which vanishes whenever the power is
    // slowly varying.
    const auto trapezoid=[&](double rate,double previousRate) {
        return rate*dt + (s.hasPreviousRates
            ? 0.5*(rate-previousRate)*s.previousStepDt : 0.0);
    };
    const auto trapezoidVector=[&](const Vec3& rate,const Vec3& previousRate) {
        return rate*dt + (s.hasPreviousRates
            ? (rate-previousRate)*(0.5*s.previousStepDt) : Vec3{});
    };
    trial.previousStepDt = dt;
    trial.hasPreviousRates = true;
    trial.previousFluxEnergy = radiation.outwardFlux.energy;
    trial.previousFluxMomentum = radiation.outwardFlux.momentum;
    trial.previousFluxAngularMomentum = radiation.outwardFlux.angularMomentum;
    trial.previousDipoleFluxEnergy = radiation.magneticDipoleFlux.energy;
    const double radiatedEnergyIncrement =
        trapezoid(radiation.outwardFlux.energy,s.previousFluxEnergy);
    trial.radiatedEnergy += radiatedEnergyIncrement;
    // Unlike the bookkeeping around it, this one DOES feed back: it is the
    // clock the orbit-following zero-point band reads.  Accumulating it here
    // means a rejected trial step discards its phase along with everything
    // else, and the adaptive error estimate sees the phase difference between
    // a full step and two half steps.  The current controller deliberately
    // bases its tolerance on the mechanical (x,v) trajectory; the phase is
    // nevertheless transactional and can never leak out of a rejected step.
    if(gZeroPointField.active())
        trial.zeroPointPhase += zeroPointDriveFrequency(s).frequency*dt;
    // Magnetic-dipole damping changes the constrained internal dipole sector.
    // The charge part is already represented by the particle self-force and
    // its near-field (Schott) term.
    const double dipoleRadiatedEnergy=trapezoid(
        radiation.magneticDipoleFlux.energy,s.previousDipoleFluxEnergy);
    // Diagnostic split below still subtracts the M1 rate unconditionally --
    // orbitalRadiatedEnergy is a decomposition of the measured Poynting flux
    // and is not a sink.  Only the mechanical drain is gated.
    // SYMMETRY AUDIT.  Under a QUANTIZED reaction model the M1 channel exerts
    // its reaction TORQUE (dipoleRadiationTorqueEnabled is true) while its energy
    // drain is gated OFF, because the quantized channel is supposed to remove
    // the energy in photons instead.  The photon hazard is built from
    // larmorOrbitAveragedPower, which is the E1 power and contains no M1, so
    // the M1 energy is radiated on paper and removed from nothing.
    //
    // THIS IS THE DESIGN, not a gap, and the audit's first reading of it as an
    // unpaired conjugate slot was wrong.  The gate is enforced: turning the
    // drain on fails quantized-radiation-gating, whose
    // whole content is that no energy may be removed continuously in this
    // mode.  Both repairs were tried -- unconditional draining, and draining
    // iff the torque is applied -- and both fail it.  So closing the gap
    // really would mean giving the photon hazard an M1 share, i.e. changing
    // the emission prescription, which is what that entry claims.
    //
    // Measured, as a fraction of the total radiated power:
    //
    //   production default (a_Ps, spin quantized):  para 9.24e-18, ortho 0
    //   barrier mode (cos sampled, deeper radii):   para 1.96e-15, ortho 1.11e-15
    //   --radiation-reaction coherent:              drain ON, no gap
    //
    // Ortho is EXACTLY zero under the default because spin quantization pins
    // cos = -1, where the coherent moment mu1 + mu2 cancels exactly -- so the
    // gap is a pure para effect, which is what makes it worth naming.  At
    // 1e-18 to 1e-15 of the total it is far below anything the model resolves,
    // so it is documented rather than fixed: draining it would require the
    // photon hazard to carry an M1 share, which is a change to the emission
    // prescription for an effect seventeen orders down.
    if(!quantizedRadiation) trial.dipoleConstraintEnergy-=dipoleRadiatedEnergy;
    // E1-only, M1 excluded: valid by linearity of the trapezoid rule in its
    // rate argument, so this is exactly trapezoid(outwardFlux-
    // magneticDipoleFlux, previousFluxEnergy-previousDipoleFluxEnergy)
    // without a third trapezoid() call.  See orbitalRadiatedEnergy's own
    // comment in state.hpp for why this is kept separate from radiatedEnergy.
    trial.orbitalRadiatedEnergy += radiatedEnergyIncrement-dipoleRadiatedEnergy;
    trial.radiatedMomentum += trapezoidVector(
        radiation.outwardFlux.momentum,s.previousFluxMomentum);
    trial.radiatedAngularMomentum += trapezoidVector(
        radiation.outwardFlux.angularMomentum,s.previousFluxAngularMomentum);

    if(computeOutwardFlux) {
        // Exact discrete world-tube balance.  This reservoir contains bound
        // near-field, retardation/interference energy and the discretization
        // remainder; it is kept separate from the independently measured
        // mismatch of the LL force and coherent far radiation below.
        const double mechanicalEnergyChange=
            conservativeParticleEnergy(trial,ledgerWeight(trial))
            -initialMechanicalEnergy;
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
            mismatch.energy-=response.magneticDipoleFlux.energy;
            mismatch.momentum=mismatch.momentum
                -response.magneticDipoleFlux.momentum;
            mismatch.angularMomentum=mismatch.angularMomentum
                -response.magneticDipoleFlux.angularMomentum;
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
            +trapezoid(initialMismatch.energy,s.previousMismatchEnergy);
        trial.reactionMomentumMismatch=balanceStart.reactionMomentumMismatch
            +trapezoidVector(initialMismatch.momentum,s.previousMismatchMomentum);
        trial.previousMismatchEnergy=initialMismatch.energy;
        trial.previousMismatchMomentum=initialMismatch.momentum;
        trial.previousMismatchAngularMomentum=initialMismatch.angularMomentum;
        trial.reactionAngularMomentumMismatch=
            balanceStart.reactionAngularMomentumMismatch
            +trapezoidVector(initialMismatch.angularMomentum,
                             s.previousMismatchAngularMomentum);
    }
    s = trial;
}

inline void appendStateHistory(StateHistory& history, const State& state) {
    const double retentionTime = std::max(1.0e-20, 4.0 * separation(state) / c);
    // Node SPACING follows the retarded window, not the integration step.
    //
    // Storing one node per accepted step ties the spacing to dt, and the
    // retarded derivatives are differenced over h = 2*spacing.  The Schott
    // sector needs a third derivative, whose stencil divides by h^3, so
    // halving dt multiplies its round-off amplification by eight.  Refining
    // the step therefore made the velocity channel WORSE: measured against
    // the collapse trajectory, the position residual fell as dt^2 from
    // 5.0e-09 to 4.6e-13 while the velocity residual rose from 1.25e-07 to
    // 4.31e-07 over the same three halvings, and the adaptive engine then
    // exhausted its subdivision depth and reported a numerical failure at
    // 1.1-2.2 collision radii.  That is why the tolerance could not be
    // tightened past 1e-6 at all.
    //
    // Pinning the spacing to retentionTime/targetHistoryNodes breaks that
    // feedback: h stops shrinking with dt, and the node count is bounded by
    // construction for every pair rather than only for e+e-, which leaves
    // the decimation below as a safety net that no longer fires in ordinary
    // running.  Raising maximumDepth instead was measured at 445x the run
    // time and still failed, and a fixed absolute h floor produced a
    // nonsense 2.3e+08 eV radiated energy at fine steps.
    //
    // This buys robustness, NOT a closed radiative balance: with the wall
    // removed, tightening the tolerance to 1e-7 still moves the radiated
    // energy from 0.412 to 0.345 eV, and the reaction mismatch stays near
    // 2.4 eV against 0.41 eV of radiation.
    //
    // That leftover is NOT the reduced-order expansion breaking down, and it
    // is worth recording why, because the opposite is the natural guess.  At
    // the e+e- collision boundary the LL parameter is tau*omega = 3.7e-4
    // (tau = 2 r_e/3c = 6.26e-24 s, omega = 5.9e19 rad/s) and beta = 0.10, so
    // the self-force expansion is nowhere near its limit.
    //
    // It is NOT a missing interference term either.  The mutual retardation
    // already carries that exchange.  Adding firstMutual on top was later
    // isolated by the independent long-horizon balance as double counting;
    // production now retains only firstSelf/secondSelf with retarded fields.
    //
    // Decomposed along the collapse (e+e-, seed 12345, eV): true retarded
    // flux 0.422, coherent E1 Larmor 0.527, sum of individual 0.268,
    // interference 0.259, reaction work +2.071.  The interference is nine
    // times too small to account for a 2.43 eV gap, and the work comes out
    // POSITIVE, i.e. dominated by the reversible Schott term rather than by
    // damping.  E_S at t=0 was measured, not assumed, at -1.2e-06 eV.
    //
    // Structurally the gap tracks -dE_S^coh: integral(flux) - dE_S^coh -
    // integral(P_E1) = 0.422 + 2.667 - 0.527 = 2.562 against a measured
    // 2.432.  The Schott boundary term, which is reversible near-field
    // energy rather than a conservation violation, is thus the bulk of it,
    // and the genuine radiation part (flux - P_E1 = -0.105 eV) is about the
    // retardation correction expected at beta=0.1.
    //
    // individualLandauLifshitzSelfForces itself is SOUND.  Tested locally,
    // step by step, its identity work = -dE_S^coh - P_E1 dt closes to
    // 2.5-10%, and e+e- and mu+mu- agree to four digits at equal r/r_coll
    // (0.0255 vs 0.0254), so the reaction sector is scale invariant as the
    // shared tau*omega = 3.7e-4 requires.
    //
    // Two earlier readings of this were measurement artefacts, recorded so
    // they are not rediscovered: comparing the work INTEGRAL against an
    // endpoint value of E_S makes the identity look 75-97% wrong for
    // mu+mu-, because E_S varies fastest exactly at the endpoint; and
    // evaluating the work at the step start while taking dE_S exactly across
    // the step leaves a spurious -0.5 P dt.  Neither survives a matched
    // quadrature.  Lowering the 1e-24 derivativeStep floor a thousandfold
    // moves the muon result by 0.6%, so that floor is not implicated either.
    //
    // One concrete trap found while checking this: on an ORDINARY BOUND
    // ORBIT the retarded window max(1e-20, 4r/c) is SHORTER than the
    // integration step, so the history holds about three nodes -- too few for
    // the third-derivative stencil, which then correctly returns zero.  The
    // coherent route is thus unavailable there by construction, and comparing
    // it against the LL route at that point compares against a near-zero
    // vector.  That is a measurement trap, not a second bug.
    constexpr std::size_t targetHistoryNodes = 128;
    // Anchor the freshness check on the last COMMITTED node
    // (history[size-2]), not on history.back() itself: back() is exactly the
    // field this branch overwrites, so comparing against it made the gap
    // reset to the current step's dt on every call and never accumulate.
    // Under sustained sub-threshold steps -- the deep, stiff part of a
    // collapsing trajectory, where recursive step halving pushes dt below
    // retentionTime/128 while the retentionTime floor stops shrinking -- that
    // meant no new node was ever committed and stale-sample eviction below
    // (which only runs past this early return) never ran either, while the
    // spacing between the last two genuine nodes grew unboundedly instead of
    // staying pinned near the target.
    if (history.size() >= 2
        && state.time - history[history.size()-2].time
               < retentionTime / static_cast<double>(targetHistoryNodes)) {
        // Keep the leading edge current without densifying the grid.
        history.back() = state;
        return;
    }
    history.push_back(state);
    const double earliestNeeded = state.time - retentionTime;
    while (history.size() > 2 && history[1].time < earliestNeeded) {
        history.pop_front();
    }
    // Backstop for the sample count inside the retarded window.  Before the
    // spacing rule above, the count was window/step, which grows like
    // r^{-1/2} and reached ten thousand entries at short range; this cap was
    // what kept that bounded.  The spacing rule now holds the count near
    // targetHistoryNodes by construction, so the decimation below is reached
    // only if a caller appends out of order.  It stays as a safety net.
    constexpr std::size_t maximumHistoryNodes = 2*targetHistoryNodes;
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
inline void integrateConservativeMidpoint(
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
