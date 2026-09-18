// The interaction field momentum of the pair, integrated (audit section 121).
//
// P_int = eps0 integral (E1 x B2 + E2 x B1) dV, the cross terms of the
// Poynting momentum: the self terms belong to each particle separately (they
// renormalize its mass and diverge at its own position), while these are
// finite and are the only part that depends on the pair's configuration.
// Sections 104-106 measured a pair that accelerates itself; 119 and 120 ruled
// out hidden momentum and the Thomas back-reaction as the missing term.  This
// probe computes the remaining candidate directly, from the model's own field
// functions, on a spherical grid around the centre of mass, and prints the
// cumulative integral against the cut-off radius so its convergence is
// visible rather than assumed.
//
// Usage: field_momentum_integral [r/r*] [channel 0 para|1 ortho] [tilt deg]
//        [steps of the orbit to integrate at: 0 = start only, N = start and
//         after N steps] [radial nodes] [polar nodes] [azimuthal nodes]
//        [outer radius in r*, default 30]
//
// The outer radius is bounded by the model itself: appendStateHistory keeps
// max(1e-20 s, 4r/c) of history, so a field point further than about 15 r*
// asks for a retarded time the history does not cover and the evaluation
// returns garbage or NaN.  The cross-term density falls as 1/r^4 in the near
// zone, so the integral converges well inside that horizon; the printed
// cumulative shows it.
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/field_momentum_integral.cpp -o /tmp/field $(root-config --libs)
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
    const int advanceSteps=argc>4?atoi(argv[4]):0;
    const int radialNodes=argc>5?atoi(argv[5]):96;
    const int polarNodes=argc>6?atoi(argv[6]):24;
    const int azimuthalNodes=argc>7?atoi(argv[7]):48;
    const double outerRadius=(argc>8?atof(argv[8]):30.0)*comptonBarrierRadius;
    const double reduced=pairReducedMass,total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(reduced*radius));
    const double period=2.0*pi*radius/speed;
    const Vec3 direction{std::sin(tilt),0.0,std::cos(tilt)};
    State state{};
    state.firstPosition={radius*secondMass/total,0.0,0.0};
    state.secondPosition={-radius*firstMass/total,0.0,0.0};
    state.firstVelocity={0.0,speed*secondMass/total,0.0};
    state.secondVelocity={0.0,-speed*firstMass/total,0.0};
    state.firstProperDipole=direction*firstMagneticMoment;
    state.secondProperDipole=direction
        *(channel?-secondMagneticMoment:secondMagneticMoment);
    synchronizeCovariantDipoles(state);
    ClassicalTrajectoryEngine::Accuracy accuracy;
    accuracy.relativeTolerance=1.0e-8;
    accuracy.maximumDepth=20;
    accuracy.reactionModel=ChargeRadiationReactionModel::disabled;
    accuracy.computeOutwardFlux=false;
    accuracy.useRetardedExternalForces=true;
    ClassicalTrajectoryEngine engine(state,accuracy);
    const double scale=firstMass*c;
    const auto mechanical=[&](const State& s) {
        return momentum(s.firstVelocity,firstMass)
            +momentum(s.secondVelocity,secondMass);
    };
    // eps0 (E1 x B2 + E2 x B1), with each particle's own charge and dipole
    // field taken from the model.
    const auto crossDensity=[&](const State& s,const Vec3& point) {
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
        return (cross(fromFirst.electric,fromSecond.magnetic)
               +cross(fromSecond.electric,fromFirst.magnetic))*epsilon0;
    };
    const Grid polar=gaussLegendre(polarNodes);
    const auto integrate=[&](const State& s,std::vector<Vec3>& cumulative,
                             std::vector<double>& shellRadius) {
        const Vec3 centre=(s.firstPosition*firstMass
            +s.secondPosition*secondMass)*(1.0/total);
        const double inner=0.05*comptonBarrierRadius;
        const double outer=outerRadius;
        Vec3 running;
        long long skipped=0;
        cumulative.clear(); shellRadius.clear();
        for(int shell=0;shell<radialNodes;++shell) {
            // Logarithmic shells, midpoint rule in log r.
            const double lowEdge=inner*std::pow(outer/inner,
                static_cast<double>(shell)/radialNodes);
            const double highEdge=inner*std::pow(outer/inner,
                static_cast<double>(shell+1)/radialNodes);
            const double centreRadius=std::sqrt(lowEdge*highEdge);
            const double shellVolume=(highEdge*highEdge*highEdge
                -lowEdge*lowEdge*lowEdge)/3.0;
            Vec3 angularSum;
            double angularWeight=0.0;
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
                    const Vec3 density=crossDensity(s,point);
                    if(isFinite(density)) angularSum+=density*weight;
                    else ++skipped;
                    angularWeight+=weight;
                }
            }
            // integral f dV over the shell = ((R^3-r^3)/3) * integral f dOmega,
            // and integral f dOmega = sum(w_polar * f) * 2 pi / azimuthalNodes.
            (void)angularWeight;
            running+=angularSum*(shellVolume*2.0*pi
                /static_cast<double>(azimuthalNodes));
            cumulative.push_back(running);
            shellRadius.push_back(highEdge);
        }
        if(skipped>0)
            std::printf("   (%lld non-finite samples skipped: the history "
                        "horizon)\n",skipped);
    };
    std::vector<Vec3> cumulative; std::vector<double> shellRadius;
    integrate(state,cumulative,shellRadius);
    const Vec3 mechanicalStart=mechanical(state);
    const Vec3 fieldStart=cumulative.back();
    std::printf("start: |P_mech| %.6e mc, |P_field(cross)| %.6e mc\n",
                mechanicalStart.norm()/scale,fieldStart.norm()/scale);
    std::printf("cumulative field momentum against the cut-off:\n");
    for(size_t index=0;index<cumulative.size();index+=radialNodes/12)
        std::printf("   R %10.3e r*   |P_field| %.6e mc\n",
                    shellRadius[index]/comptonBarrierRadius,
                    cumulative[index].norm()/scale);
    if(advanceSteps>0) {
        const double step=period/256.0;
        for(int index=0;index<advanceSteps;++index)
            if(!engine.advance(state,step)) {
                std::printf("FAILED after %d steps\n",index); return 1;
            }
        integrate(state,cumulative,shellRadius);
        const Vec3 mechanicalEnd=mechanical(state);
        const Vec3 fieldEnd=cumulative.back();
        std::printf("\nafter %d steps (r %.4f r*): |P_mech| %.6f mc, "
                    "|P_field| %.6e mc\n",advanceSteps,
                    separation(state)/comptonBarrierRadius,
                    mechanicalEnd.norm()/scale,fieldEnd.norm()/scale);
        std::printf("changes: |dP_mech| %.6f mc, |dP_field| %.6e mc, "
                    "|d(P_mech+P_field)| %.6f mc, ratio %.4f\n",
                    (mechanicalEnd-mechanicalStart).norm()/scale,
                    (fieldEnd-fieldStart).norm()/scale,
                    ((mechanicalEnd+fieldEnd)
                        -(mechanicalStart+fieldStart)).norm()/scale,
                    (mechanicalEnd-mechanicalStart).norm()
                        /std::max(((mechanicalEnd+fieldEnd)
                            -(mechanicalStart+fieldStart)).norm(),1.0e-300));
    }
}
