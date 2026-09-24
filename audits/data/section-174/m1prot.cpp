// Is |mu1+mu2| = 0 dynamically protected in exactly anti-aligned ortho,
// or only true at preparation?
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
int main(){
    const double total=firstMass+secondMass;
    const double r=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r));
    const double period=2.0*pi*r/speed;
    const double tilt=40.0*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    for(int channel=0;channel<2;++channel){
        State s{};
        s.firstPosition={r*secondMass/total,0,0};
        s.secondPosition={-r*firstMass/total,0,0};
        s.firstVelocity={0,speed*secondMass/total,0};
        s.secondVelocity={0,-speed*firstMass/total,0};
        s.firstProperDipole=dir*firstMagneticMoment;
        // channel 0 = para (aligned moments), 1 = ortho (exactly opposite)
        s.secondProperDipole=dir*(channel? -secondMagneticMoment
                                         :  secondMagneticMoment);
        synchronizeCovariantDipoles(s);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1e-10; acc.maximumDepth=20;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(s,acc);
        const double scale=firstMagneticMoment+secondMagneticMoment;
        double worst=0.0, worstRate=0.0;
        std::printf("%s, exact preparation:\n",channel?"ORTHO":"PARA ");
        for(int i=0;i<=4000;++i){
            if(i&&!engine.advance(s,period/4000)) break;
            const double m=(s.firstDipole+s.secondDipole).norm()/scale;
            const DipoleDerivatives d=
                thomasBmtDipoleDerivatives(s,engine.history());
            const double rate=(d.first+d.second).norm()
                /std::max((d.first.norm()+d.second.norm()),1e-300);
            worst=std::max(worst,m);
            worstRate=std::max(worstRate,rate);
            if(i%1000==0)
                std::printf("   orbit %.2f  |m|/mu_sum=%.6e"
                            "  |d(m)/dt|/sum=%.6e\n",
                            i/4000.0,m,rate);
        }
        std::printf("   worst over the orbit: |m|=%.6e   |dm/dt|=%.6e\n\n",
                    worst,worstRate);
    }
}
