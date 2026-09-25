// 181j/181k: the residual over six orbits, with the drift off so the
// moment-azimuth condition is the only one in play.

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
#include <initializer_list>
namespace {
void run(double rf,double bd,int orbits){
    const double total=firstMass+secondMass;
    const double r0=rf*pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double period=2.0*pi*r0/speed;
    const double th=40.0*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity=Vec3{bd*c,speed*secondMass/total,0};
    s.secondVelocity=Vec3{bd*c,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-10; acc.maximumDepth=20;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    std::printf("  r = %-5g a_pair   beta_orb = %.4f   beta_drift = %-6g\n",
                rf,speed/c,bd);
    std::printf("  %8s %14s %14s\n","t/T","worst this orbit","at whole turn");
    const int perOrbit=400;
    double worst=0.0;
    for(int i=0;i<=perOrbit*orbits;++i){
        if(i&&!engine.advance(s,period/perOrbit)) break;
        const Vec3 sep=s.firstPosition-s.secondPosition;
        const LocalElectromagneticFields f=
            localRelativisticFields(s,engine.history());
        const Vec3 b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                              firstGFactor);
        const Vec3 b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                              secondGFactor);
        const Vec3 nh=sep/sep.norm(), zh{0,0,1};
        const Vec3 perp=cross(zh,nh);
        const Vec3 target=(b1-perp*(2.0*dot(b1,perp)))*-1.0;
        const double r=(b2-target).norm()/std::max(b1.norm(),1e-300);
        worst=std::max(worst,r);
        if(i%perOrbit==0){
            std::printf("  %8.1f %14.4e %14.4e\n",double(i)/perOrbit,worst,r);
            worst=0.0;
        }
    }
    std::printf("\n");
}
}
int main(){
    std::printf("Drift ZERO: the condition then holds at EVERY\n"
                "orbit phase; any residual is history alone.\n\n");
    run(1.0,0.0,6);
    run(0.1,0.0,6);
    std::printf("With the drift (the angle confound is back):"
                "\n\n");
    run(1.0,0.01,6);
}
