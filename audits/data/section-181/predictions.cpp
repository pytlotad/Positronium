// 181c/181d/181e: the three preconditions the operation predicts --
// the drift axis, the moment azimuth, and equal moment magnitudes.

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
// E and v are polar: pattern (+x,-y,+z).  B is axial: (-x,+y,-z).
Vec3 polar(const Vec3& v){ return { v.x,-v.y, v.z}; }
Vec3 axial(const Vec3& v){ return {-v.x, v.y,-v.z}; }
struct Res { double eRes,bRes,effRes; };
Res at(const Vec3& drift,double thDeg,double phiDeg,double r0){
    const double total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thDeg*pi/180, ph=phiDeg*pi/180;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),
                   std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity=Vec3{0,speed*secondMass/total,0}+drift;
    s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+drift;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 e1=f.atFirst.electric,e2=f.atSecond.electric;
    const Vec3 b1=f.atFirst.magnetic,b2=f.atSecond.magnetic;
    const Vec3 k1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                          firstGFactor);
    const Vec3 k2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                          secondGFactor);
    return {(e2-polar(e1)).norm()/std::max(e1.norm(),1e-300),
            (b2-axial(b1)).norm()/std::max(b1.norm(),1e-300),
            (k2-axial(k1)).norm()/std::max(k1.norm(),1e-300)};
}
}
int main(){
    const double a=pairBohrRadius(activePair);
    std::printf("P1.  The drift must have NO component along y\n"
                "(the orbital velocity); 180d tried only x and z.\n");
    std::printf("%-30s %11s %11s %11s\n","drift direction","E","B","B_eff");
    struct D { const char* n; Vec3 v; };
    for(const D& d:{D{"none",{0,0,0}},
                    D{"0.01 c along x (separation)",{0.01*c,0,0}},
                    D{"0.01 c along z (orbit normal)",{0,0,0.01*c}},
                    D{"0.01 c along y (orbital vel.)",{0,0.01*c,0}},
                    D{"0.001 c along y",{0,0.001*c,0}},
                    D{"1e-5 c along y",{0,1e-5*c,0}}}){
        const Res r=at(d.v,40.0,0.0,a);
        std::printf("%-30s %11.3e %11.3e %11.3e\n",d.n,
                    r.eRes,r.bRes,r.effRes);
    }
    std::printf("\nP2.  It needs the moments out of y, which is the azimuth\n"
                "condition -- now a prediction, not a finding.\n");
    std::printf("%-30s %11s %11s %11s\n","moment azimuth","E","B","B_eff");
    for(double p:{0.0,0.001,0.01,0.1,1.0,45.0}){
        const Res r=at(Vec3{0.01*c,0,0},40.0,p,a);
        char buf[32]; std::snprintf(buf,sizeof buf,"phi = %g deg",p);
        std::printf("%-30s %11.3e %11.3e %11.3e\n",buf,
                    r.eRes,r.bRes,r.effRes);
    }
    std::printf("\nP3.  It needs charge conjugation, so it needs the two\n"
                "moments EQUAL in magnitude.  Try an unequal pair.\n");
    std::printf("%-30s %11s %11s %11s\n","pair","E","B","B_eff");
    {   const Res r=at(Vec3{0.01*c,0,0},40.0,0.0,a);
        std::printf("%-30s %11.3e %11.3e %11.3e\n","electron + positron",
                    r.eRes,r.bRes,r.effRes); }
    applyPair(ParticlePair{electron,antimuon});
    {   const double a2=pairBohrRadius(activePair);
        const Res r=at(Vec3{0.01*c,0,0},40.0,0.0,a2);
        std::printf("%-30s %11.3e %11.3e %11.3e\n","electron + antimuon",
                    r.eRes,r.bRes,r.effRes);
        std::printf("   moment ratio |mu1|/|mu2| = %.4f\n",
                    firstMagneticMoment/secondMagneticMoment); }
}
