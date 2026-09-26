// Audit 225: is the projection of the magnetic moment on the orbital
// radius quantized, for para and for ortho?
//
// Three projections are followed over one orbit, for both channels
// under --spin-quantization, where 213b measured mu1.mu2 = +-1.00000
// exactly:  mu1.rhat,  mu2.rhat,  and (mu1+mu2).rhat.
// The seed's only freedom is the moment direction, drawn uniformly on
// the sphere (213b), so delta = angle(mu1, orbit normal) is the one
// parameter and rhat sweeps the equator beneath it.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    std::printf("%8s %6s %11s %11s %11s %11s %13s\n",
                "delta","chan","min m1.rh","max m1.rh","mean m1.rh",
                "sin(delta)","max|(m1+m2).rh|");
    for(double deltaDeg:{0.0,30.0,60.0,90.0,135.0}){
        const double th=deltaDeg*pi/180.0;
        const Vec3 unit{std::sin(th),0.0,std::cos(th)};
        for(int ch=0;ch<2;++ch){
            const Vec3 m1=unit*firstMagneticMoment;
            const Vec3 m2=unit*((ch==0?+1.0:-1.0)*secondMagneticMoment);
            const State start=osculatingPeriapsisState(
                el,k,m1,m2,Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            ClassicalTrajectoryEngine engine(start,acc);
            State s=start;
            double lo=1e300,hi=-1e300,sum=0.0,worstSum=0.0;
            const int n=512; bool ok=true;
            for(int i=0;i<n&&ok;++i){
                ok=engine.advance(s,period/n);
                const Vec3 sep=s.firstPosition-s.secondPosition;
                const Vec3 rh=sep*(1.0/sep.norm());
                const double p1=dot(s.firstDipole,rh)/firstMagneticMoment;
                const double p2=dot(s.secondDipole,rh)/secondMagneticMoment;
                lo=std::min(lo,p1); hi=std::max(hi,p1); sum+=p1;
                worstSum=std::max(worstSum,std::abs(p1+p2));
            }
            if(!ok){ std::printf("%8.1f %6s failed\n",deltaDeg,
                     ch==0?"para":"ortho"); continue; }
            std::printf("%8.1f %6s %11.6f %11.6f %11.3e %11.6f %13.3e\n",
                deltaDeg,ch==0?"para":"ortho",lo,hi,sum/n,
                std::sin(th),worstSum);
        }
    }
    std::printf("\n# m1.rhat is normalised to |mu1|, so a quantized\n");
    std::printf("# projection would show discrete values; a swept one\n");
    std::printf("# fills [-sin(delta), +sin(delta)] continuously.\n");
}
