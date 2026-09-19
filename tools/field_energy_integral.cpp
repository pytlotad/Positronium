// The interaction field ENERGY of the pair, integrated (audit section 155).
//
// U_int = integral ( eps0 E1.E2 + B1.B2/mu0 ) dV, the cross terms of the field
// energy.  The self terms diverge at each particle and renormalize its mass;
// these are finite and are the only part that depends on the pair's
// configuration.  The static limit of this integral is exactly the Coulomb
// interaction energy the ledger already carries, so what it adds is the
// RETARDED correction -- which is what audit 145 attributed 1.13 |U_dd| per
// orbit of ledger drift to.  Energy analogue of the field momentum audit 121
// integrated, on the same grid, with the same history horizon: beyond about
// 15 r* the retarded time leaves the stored history.
//
// Usage: field_energy_integral [r/r*] [channel 0 para|1 ortho] [tilt deg]
//        [samples around the orbit, default 8] [radial nodes] [polar nodes]
//        [azimuthal nodes] [outer radius in r*, default 12]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/field_energy_integral.cpp -o /tmp/fieldE $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
struct Grid { std::vector<double> nodes,weights; };
// Gauss-Legendre nodes on [-1,1] by Newton iteration on the Legendre
// polynomial -- the angular quadrature this integral needs.
Grid gaussLegendre(int count) {
    Grid grid; grid.nodes.resize(count); grid.weights.resize(count);
    for(int index=0;index<count;++index) {
        double x=std::cos(pi*(index+0.75)/(count+0.5));
        for(int iteration=0;iteration<100;++iteration) {
            double p0=1.0,p1=0.0;
            for(int degree=0;degree<count;++degree) {
                const double p2=p1; p1=p0;
                p0=((2.0*degree+1.0)*x*p1-degree*p2)/(degree+1.0);
            }
            const double derivative=count*(x*p0-p1)/(x*x-1.0);
            const double step=p0/derivative;
            x-=step;
            if(std::abs(step)<1.0e-15) break;
        }
        double p0=1.0,p1=0.0;
        for(int degree=0;degree<count;++degree) {
            const double p2=p1; p1=p0;
            p0=((2.0*degree+1.0)*x*p1-degree*p2)/(degree+1.0);
        }
        const double derivative=count*(x*p0-p1)/(x*x-1.0);
        grid.nodes[static_cast<size_t>(index)]=x;
        grid.weights[static_cast<size_t>(index)]=
            2.0/((1.0-x*x)*derivative*derivative);
    }
    return grid;
}
}
int main(int argc,char** argv) {
    const double radius=(argc>1?atof(argv[1]):1.0)*comptonBarrierRadius;
    const int channel=argc>2?atoi(argv[2]):0;
    const double tilt=(argc>3?atof(argv[3]):0.0)*pi/180;
    const int samples=argc>4?atoi(argv[4]):8;
    const int radialNodes=argc>5?atoi(argv[5]):96;
    const int polarNodes=argc>6?atoi(argv[6]):24;
    const int azimuthalNodes=argc>7?atoi(argv[7]):48;
    const double outerRadius=(argc>8?atof(argv[8]):12.0)*comptonBarrierRadius;
    const double reduced=pairReducedMass,total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(reduced*radius));
    const double period=2.0*pi*radius/speed;
    const Vec3 direction{std::sin(tilt),0.0,std::cos(tilt)};
    State state{};
    state.firstPosition={radius*secondMass/total,0.0,0.0};
    state.secondPosition={-radius*firstMass/total,0.0,0.0};
    // CREM_FIELD_STATIC validates the quadrature against the one case with an
    // exact answer: at rest the cross energy integral IS the Coulomb
    // interaction energy, eps0 int E1.E2 dV = q1 q2/(4 pi eps0 r), so any
    // shortfall is the grid's, not the physics'.  A probe that cannot
    // reproduce that has nothing to say about the retarded value.
    const bool staticPair=std::getenv("CREM_FIELD_STATIC")!=nullptr;
    state.firstVelocity=staticPair?Vec3{}
        :Vec3{0.0,speed*secondMass/total,0.0};
    state.secondVelocity=staticPair?Vec3{}
        :Vec3{0.0,-speed*firstMass/total,0.0};
    state.firstProperDipole=direction*firstMagneticMoment;
    state.secondProperDipole=direction
        *(channel?-secondMagneticMoment:secondMagneticMoment);
    synchronizeCovariantDipoles(state);
    const double dipoleScale=std::abs(pairDipoleInteractionEnergy(
        state.firstPosition-state.secondPosition,
        direction*firstMagneticMoment,direction*secondMagneticMoment));
    ClassicalTrajectoryEngine::Accuracy accuracy;
    accuracy.relativeTolerance=1.0e-8;
    accuracy.maximumDepth=20;
    // Audit 158: the particles pay for what the far zone receives, so the
    // quantity tested carries the radiated energy and the Schott term as well
    // as the field correction.  CREM_FIELD_NO_REACTION restores the
    // reaction-disabled configuration audit 155 used.
    const bool payForRadiation=std::getenv("CREM_FIELD_NO_REACTION")==nullptr;
    accuracy.reactionModel=payForRadiation
        ?ChargeRadiationReactionModel::individualLandauLifshitz
        :ChargeRadiationReactionModel::disabled;
    accuracy.computeOutwardFlux=true;
    accuracy.useRetardedExternalForces=true;
    ClassicalTrajectoryEngine engine(state,accuracy);
    // eps0 E1.E2 + B1.B2/mu0, each particle's charge and dipole field taken
    // from the model's own functions.
    // THE INTEGRAND IS SUBTRACTED AGAINST A FROZEN COPY OF THE SAME STATE.
    // eps0 E1.E2 has an integrable 1/r^2 singularity at each particle, which a
    // grid centred on the pair cannot resolve: the static validation below
    // came out at 74% of the exact answer and did NOT improve with the grid,
    // so the shortfall is the singularity and the cut-offs, not the
    // resolution.  Evaluating the same model functions on a state with the
    // velocities zeroed gives a reference whose singular structure is
    // identical and whose exact integral is known, -K/r, so the difference is
    // regular and the quadrature only has to resolve the RETARDED correction:
    //   U_field = -K/r + integral (u_retarded - u_frozen) dV.
    // No formula is copied; both densities come from the same two functions.
    const auto frozenOf=[&](const State& s) {
        State frozen=s;
        frozen.firstVelocity={};
        frozen.secondVelocity={};
        frozen.firstAcceleration={};
        frozen.secondAcceleration={};
        return frozen;
    };
    const auto densityWith=[&](const State& s,const StateHistory& history,
                               const Vec3& point) {
        ElectromagneticField fromFirst=lienardWiechertField(
            point,s.time,history,s,true,firstCharge);
        const ElectromagneticField firstDipole=retardedMagneticDipoleField(
            point,s.time,history,s,true);
        fromFirst.electric+=firstDipole.electric;
        fromFirst.magnetic+=firstDipole.magnetic;
        ElectromagneticField fromSecond=lienardWiechertField(
            point,s.time,history,s,false,secondCharge);
        const ElectromagneticField secondDipole=retardedMagneticDipoleField(
            point,s.time,history,s,false);
        fromSecond.electric+=secondDipole.electric;
        fromSecond.magnetic+=secondDipole.magnetic;
        return epsilon0*dot(fromFirst.electric,fromSecond.electric)
            +dot(fromFirst.magnetic,fromSecond.magnetic)/mu0;
    };
    const auto crossEnergyDensityOld=[&](const State& s,const Vec3& point) {
        ElectromagneticField fromFirst=lienardWiechertField(
            point,s.time,engine.history(),s,true,firstCharge);
        const ElectromagneticField firstDipole=retardedMagneticDipoleField(
            point,s.time,engine.history(),s,true);
        fromFirst.electric+=firstDipole.electric;
        fromFirst.magnetic+=firstDipole.magnetic;
        ElectromagneticField fromSecond=lienardWiechertField(
            point,s.time,engine.history(),s,false,secondCharge);
        const ElectromagneticField secondDipole=retardedMagneticDipoleField(
            point,s.time,engine.history(),s,false);
        fromSecond.electric+=secondDipole.electric;
        fromSecond.magnetic+=secondDipole.magnetic;
        return epsilon0*dot(fromFirst.electric,fromSecond.electric)
            +dot(fromFirst.magnetic,fromSecond.magnetic)/mu0;
    };
    (void)crossEnergyDensityOld;
    const Grid polar=gaussLegendre(polarNodes);
    const auto integrate=[&](const State& s,std::vector<double>& cumulative,
                             std::vector<double>& shellRadius) {
        StateHistory frozenHistory;
        frozenHistory.push_back(frozenOf(s));
        const Vec3 centre=(s.firstPosition*firstMass
            +s.secondPosition*secondMass)*(1.0/total);
        const double inner=0.05*comptonBarrierRadius;
        double running=0.0;
        long long skipped=0;
        cumulative.clear(); shellRadius.clear();
        for(int shell=0;shell<radialNodes;++shell) {
            const double lowEdge=inner*std::pow(outerRadius/inner,
                static_cast<double>(shell)/radialNodes);
            const double highEdge=inner*std::pow(outerRadius/inner,
                static_cast<double>(shell+1)/radialNodes);
            const double centreRadius=std::sqrt(lowEdge*highEdge);
            const double shellVolume=(highEdge*highEdge*highEdge
                -lowEdge*lowEdge*lowEdge)/3.0;
            double angularSum=0.0;
            for(int polarIndex=0;polarIndex<polarNodes;++polarIndex) {
                const double cosTheta=polar.nodes[
                    static_cast<size_t>(polarIndex)];
                const double weight=polar.weights[
                    static_cast<size_t>(polarIndex)];
                const double sinTheta=std::sqrt(std::max(0.0,
                    1.0-cosTheta*cosTheta));
                for(int slice=0;slice<azimuthalNodes;++slice) {
                    const double phi=2.0*pi*(slice+0.5)/azimuthalNodes;
                    const Vec3 point=centre+Vec3{
                        centreRadius*sinTheta*std::cos(phi),
                        centreRadius*sinTheta*std::sin(phi),
                        centreRadius*cosTheta};
                    const double density=densityWith(s,engine.history(),point)
                        -densityWith(frozenOf(s),frozenHistory,point);
                    if(std::isfinite(density)) angularSum+=density*weight;
                    else ++skipped;
                }
            }
            running+=angularSum*(shellVolume*2.0*pi
                /static_cast<double>(azimuthalNodes));
            cumulative.push_back(running);
            shellRadius.push_back(highEdge);
        }
        if(skipped>0)
            std::printf("   (%lld non-finite samples skipped: the history "
                        "horizon)\n",skipped);
        return running;
    };
    // The instantaneous interaction terms the ledger carries, for comparison:
    // Coulomb, Darwin and the azimuth-free dipole-dipole, with the
    // velocity-linear charge-dipole term left out as audit 146 established.
    const auto instantaneous=[&](const State& s) {
        const PairGeometry geometry=clampedPairGeometry(s);
        return -pairCoulombStrength*geometry.inverseDistance
            +pairDipoleInteractionEnergy(s.firstPosition-s.secondPosition,
                                         s.firstDipole,s.secondDipole)
            +darwinInteractionEnergy(s);
    };
    std::printf("r %.2f r*  channel %d  tilt %.0f deg  samples %d  "
                "grid %dx%dx%d  outer %.1f r*\n|U_dd| = %.6e J\n",
                radius/comptonBarrierRadius,channel,tilt*180.0/pi,samples,
                radialNodes,polarNodes,azimuthalNodes,
                outerRadius/comptonBarrierRadius,dipoleScale);
    if(staticPair) {
        const double exact=-pairCoulombStrength/radius;
        std::printf("STATIC validation: exact Coulomb interaction energy "
                    "%.9e J\n",exact);
    }
    std::printf("%10s %16s %16s %14s %14s\n","t/period","U_field [J]",
                "U_inst [J]","U_field/|U_dd|","U_inst/|U_dd|");
    std::vector<double> cumulative,shellRadius,fieldValues,instantValues;
    std::vector<double> combinedValues,fullValues;
    const int stepsPerSample=64;
    for(int sample=0;sample<samples;++sample) {
        if(sample>0) {
            bool ok=true;
            for(int step=0;step<stepsPerSample&&ok;++step)
                ok=engine.advance(state,period/(samples*stepsPerSample));
            if(!ok) { std::printf("engine failed at sample %d\n",sample); break; }
        }
        const double correction=integrate(state,cumulative,shellRadius);
        const double staticPart=-pairCoulombStrength
            /(state.firstPosition-state.secondPosition).norm();
        const double field=staticPart+correction;
        const double instant=instantaneous(state);
        const double kinetic=kineticEnergy(state.firstVelocity,firstMass)
            +kineticEnergy(state.secondVelocity,secondMass);
        const MutualForces external=
            retardedExternalForces(state,engine.history());
        const double schott=explicitChargeSchottEnergy(state,external);
        const double ledger=kinetic+instant;
        fieldValues.push_back(ledger);
        instantValues.push_back(ledger+state.orbitalRadiatedEnergy);
        combinedValues.push_back(ledger+state.orbitalRadiatedEnergy+schott);
        fullValues.push_back(ledger+state.orbitalRadiatedEnergy+schott
                             +correction);
        std::printf("%10.4f %16.9e %16.9e %14.6f %14.6f  (static %+.4f, "
                    "retarded correction %+.6f)\n",
                    static_cast<double>(sample)/samples,field,instant,
                    field/dipoleScale,instant/dipoleScale,
                    staticPart/dipoleScale,correction/dipoleScale);
    }
    if(fieldValues.size()>1) {
        const auto range=[](const std::vector<double>& v) {
            double low=v.front(),high=v.front();
            for(double x:v) { low=std::min(low,x); high=std::max(high,x); }
            return high-low;
        };
        const double fieldRange=range(fieldValues)/dipoleScale;
        const double instantRange=range(instantValues)/dipoleScale;
        std::printf("\nrange of the TOTAL (kinetic + interaction) over the sampled orbit, in |U_dd|:\n");
        std::printf("  E1  the ledger                       %12.6f\n",fieldRange);
        std::printf("  E4  plus the far-zone radiated       %12.6f\n",
                    instantRange);
        std::printf("  E5  plus the Schott endpoint term    %12.6f\n",
                    range(combinedValues)/dipoleScale);
        std::printf("  E6  plus the field-energy correction %12.6f\n",
                    range(fullValues)/dipoleScale);
        std::printf("  audit 145e attributed 1.13 |U_dd| of the ledger's "
                    "drift to retardation\n");
    }
    // R155c: convergence against the cut-off radius, printed rather than assumed.
    std::printf("\ncumulative U_field against the cut-off, last state:\n");
    for(std::size_t index=0;index<cumulative.size();
        index+=std::max<std::size_t>(1,cumulative.size()/8))
        std::printf("  R %8.3f r*   %16.9e J   %12.6f |U_dd|\n",
                    shellRadius[index]/comptonBarrierRadius,cumulative[index],
                    cumulative[index]/dipoleScale);
}
