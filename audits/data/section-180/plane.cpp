// 180d: which plane is the mirror -- 179b named it wrongly
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        plane.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
Vec3 mirrorPseudo(const Vec3& b,const Vec3& n){ return n*(2.0*dot(b,n))-b; }
// vr = radial velocity (along the separation), vz = out of the orbit plane
double residual(double vr,double vz,const Vec3& n,double r0=0.0){
    const double a=pairBohrRadius(activePair);
    const double r=(r0>0.0)?r0:a, total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r));
    const double th=40.0*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    const Vec3 v1{vr*c,speed*secondMass/total,vz*c};
    const Vec3 v2{vr*c,-speed*firstMass/total,vz*c};
    State s{};
    s.firstPosition={r*secondMass/total,0,0};
    s.secondPosition={-r*firstMass/total,0,0};
    s.firstVelocity=v1; s.secondVelocity=v2;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 b1=thomasBmtEffectiveField(v1,f.atFirst,firstGFactor);
    const Vec3 b2=thomasBmtEffectiveField(v2,f.atSecond,secondGFactor);
    return (b2-mirrorPseudo(b1,n)).norm()/std::max(b1.norm(),1e-300);
}
}
int main(){
    const Vec3 ny{0,1,0}, nx{1,0,0};
    std::printf("Which plane is the mirror?  Add velocity components that\n"
                "separate the two candidate names.  "
                "(theta 40, r = a_pair)\n\n");
    std::printf("%-34s %14s %14s\n","velocity added","normal = y",
                "normal = sep");
    struct Row { const char* name; double vr,vz; };
    for(const Row& t:{Row{"none (circular)",0.0,0.0},
                      Row{"radial 0.01 c (along the sep)",0.01,0.0},
                      Row{"radial 0.05 c",0.05,0.0},
                      Row{"out of plane 0.01 c (along z)",0.0,0.01},
                      Row{"radial 0.03 + out of plane 0.01",0.03,0.01}})
        std::printf("%-34s %14.4e %14.4e\n",t.name,
                    residual(t.vr,t.vz,ny),residual(t.vr,t.vz,nx));
}
