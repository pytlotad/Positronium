// 181g: is the preparation residual physical or round-off?  Scaling
// the orbital speed to zero at fixed geometry settles it.

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
#include <initializer_list>
namespace {
Vec3 polar(const Vec3& v){ return { v.x,-v.y, v.z}; }
Vec3 axial(const Vec3& v){ return {-v.x, v.y,-v.z}; }
struct R { double e,b,accel; };
R at(double rf,double orbitScale){
    const double total=firstMass+secondMass;
    const double r0=rf*pairBohrRadius(activePair);
    // orbitScale multiplies the orbital speed: 0 kills the acceleration
    // while leaving the geometry and the drift untouched.
    const double speed=orbitScale
        *std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=40.0*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity={0.01*c,speed*secondMass/total,0};
    s.secondVelocity={0.01*c,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 e1=f.atFirst.electric,e2=f.atSecond.electric;
    const Vec3 b1=f.atFirst.magnetic,b2=f.atSecond.magnetic;
    const double acc=pairCoulombStrength/(pairReducedMass*r0*r0);
    return {(e2-polar(e1)).norm()/e1.norm(),
            (b2-axial(b1)).norm()/b1.norm(),acc};
}
}
int main(){
    std::printf("Hypothesis: exact for the VELOCITY part of the\n"
                "field, broken only by the ACCELERATION part, so\n"
                "the residual should track ~ r^-2.\n\n");
    std::printf("%8s %12s %12s %12s %12s %12s\n","r/a_pair","accel",
                "E residual","B residual","E/accel","B/accel");
    double e0=0,b0=0,a0=0;
    for(double rf:{1.0,0.3,0.1,0.03,0.01}){
        const R r=at(rf,1.0);
        if(a0==0){ e0=r.e; b0=r.b; a0=r.accel; }
        std::printf("%8.3g %12.4e %12.4e %12.4e %12.4f %12.4f\n",rf,r.accel,
            r.e,r.b,(r.e/e0)/(r.accel/a0),(r.b/b0)/(r.accel/a0));
    }
    std::printf("\nDirect test: scale the orbital speed down at fixed\n"
                "geometry, which scales the acceleration as its square.\n");
    std::printf("%12s %12s %12s %12s\n","orbit scale","accel share",
                "E residual","B residual");
    for(double q:{1.0,0.3,0.1,0.03,0.0}){
        const R r=at(1.0,q);
        std::printf("%12.3g %12.4e %12.4e %12.4e\n",q,q*q,r.e,r.b);
    }
}
