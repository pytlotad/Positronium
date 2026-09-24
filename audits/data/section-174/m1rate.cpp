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
    std::printf("%-6s %14s %14s %16s %16s\n","chan","|w1|","|w2|",
                "|w1-w2|/|w1|","M1 flux [W]");
    for(int channel=0;channel<2;++channel){
        State s{};
        s.firstPosition={r*secondMass/total,0,0};
        s.secondPosition={-r*firstMass/total,0,0};
        s.firstVelocity={0,speed*secondMass/total,0};
        s.secondVelocity={0,-speed*firstMass/total,0};
        s.firstProperDipole=dir*firstMagneticMoment;
        s.secondProperDipole=dir*(channel? -secondMagneticMoment
                                         :  secondMagneticMoment);
        synchronizeCovariantDipoles(s);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1e-10; acc.maximumDepth=20;
        acc.computeOutwardFlux=true;
        ClassicalTrajectoryEngine engine(s,acc);
        for(int i=0;i<200;++i) if(!engine.advance(s,period/2000)) break;
        // precession rates of the two moments
        const DipoleDerivatives d=thomasBmtDipoleDerivatives(s,engine.history());
        const double w1=d.first.norm()/std::max(s.firstDipole.norm(),1e-300);
        const double w2=d.second.norm()/std::max(s.secondDipole.norm(),1e-300);
        // the model's own M1 flux
        const DipoleRadiationReaction m1r=dipoleRadiationReaction(s,engine.history());
        std::printf("%-6s %14.6e %14.6e %16.3e %16.6e\n",
            channel?"ortho":"para",w1,w2,
            std::abs(w1-w2)/std::max(w1,1e-300),m1r.power);
    }
}
