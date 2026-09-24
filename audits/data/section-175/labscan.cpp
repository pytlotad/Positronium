#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
int main(int argc,char** argv){
    const bool retarded=(argc>1&&std::atoi(argv[1])!=0);
    const double total=firstMass+secondMass;
    const double tilt=40.0*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    std::printf("%s forces, reaction off, one orbit, tilt 40\n",
                retarded?"RETARDED":"instantaneous");
    std::printf("%10s %8s %13s %13s %13s %13s\n","r","beta",
                "U_dd lab rng","U_dd prop rng","lab-prop diff","E_cons rng");
    std::printf("%10s %8s %13s %13s %13s %13s\n","","","[k/r0]","[k/r0]",
                "[k/r0]","[k/r0]");
    for(double f:{1.0,0.1,0.01,0.003,0.001}){
        const double r0=f*pairBohrRadius(activePair);
        const double speed=std::sqrt(pairCoulombStrength
            /(pairReducedMass*r0));
        const double period=2.0*pi*r0/speed;
        const double unit=pairCoulombStrength/r0;
        State s{};
        s.firstPosition={r0*secondMass/total,0,0};
        s.secondPosition={-r0*firstMass/total,0,0};
        s.firstVelocity={0,speed*secondMass/total,0};
        s.secondVelocity={0,-speed*firstMass/total,0};
        s.firstProperDipole=dir*firstMagneticMoment;
        s.secondProperDipole=dir*secondMagneticMoment;
        synchronizeCovariantDipoles(s);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1e-10; acc.maximumDepth=22;
        acc.computeOutwardFlux=false;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.useRetardedExternalForces=retarded;
        ClassicalTrajectoryEngine engine(s,acc);
        double lo=1e300,hi=-1e300,loP=1e300,hiP=-1e300;
        double loE=1e300,hiE=-1e300,worstDiff=0.0;
        bool ok=true;
        for(int i=0;i<=2000;++i){
            if(i&&!engine.advance(s,period/2000)){ ok=false; break; }
            const Vec3 sep=s.firstPosition-s.secondPosition;
            const double a=pairDipoleInteractionEnergy(sep,
                s.firstDipole,s.secondDipole);
            const double b=pairDipoleInteractionEnergy(sep,
                s.firstProperDipole,s.secondProperDipole);
            const double e=conservedParticleEnergy(s);
            lo=std::min(lo,a); hi=std::max(hi,a);
            loP=std::min(loP,b); hiP=std::max(hiP,b);
            loE=std::min(loE,e); hiE=std::max(hiE,e);
            worstDiff=std::max(worstDiff,std::abs(a-b));
        }
        if(!ok){ std::printf("%10.3g   (engine stopped)\n",f); continue; }
        std::printf("%10.3g %8.4f %13.4e %13.4e %13.4e %13.4e\n",
            f,speed/c,(hi-lo)/unit,(hiP-loP)/unit,worstDiff/unit,
            (hiE-loE)/unit);
    }
}
