// 181j: is the growth the identity degrading, or the orbit
// not being back?  The phase offset is measured alongside
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        phasecontrol.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
int main(){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double period=2.0*pi*r0/speed;
    const double th=40.0*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity={0,speed*secondMass/total,0};
    s.secondVelocity={0,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-10; acc.maximumDepth=20;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    std::printf("Drift zero.  At each whole turn the geometry is meant to\n"
                "return to the symmetric one.  Is the residual there the\n"
                "identity degrading, or just the orbit not being back?\n"
                "The drift-free slope at this configuration is 6.116e-03\n"
                "per degree (audit 180h).\n\n");
    std::printf("%6s %14s %14s %16s %12s\n","turn","angle off [deg]",
                "residual","residual/angle","vs slope");
    const int perOrbit=400;
    for(int i=0;i<=perOrbit*6;++i){
        if(i&&!engine.advance(s,period/perOrbit)) break;
        if(i%perOrbit) continue;
        const Vec3 sep=s.firstPosition-s.secondPosition;
        double a=std::atan2(sep.y,sep.x)*180.0/pi;
        while(a>90.0) a-=180.0; while(a<-90.0) a+=180.0;
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
        const double q=(std::abs(a)>0.0)?r/std::abs(a):0.0;
        std::printf("%6d %14.6f %14.4e %16.4e %12.3f\n",i/perOrbit,a,r,q,
                    q/6.116e-03);
    }
}
