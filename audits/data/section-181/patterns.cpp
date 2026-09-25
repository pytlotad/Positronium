// 181b: all eight sign patterns with the components, to show
// the +-y discrimination is real once the drift is on
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        patterns.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
struct Out { Vec3 e1,e2,b1,b2,v1,v2; };
Out probe(double bd,double thDeg,double phiDeg,double r0){
    const double total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thDeg*pi/180, ph=phiDeg*pi/180;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),
                   std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity={bd*c,speed*secondMass/total,0};
    s.secondVelocity={bd*c,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    return {f.atFirst.electric,f.atSecond.electric,
            f.atFirst.magnetic,f.atSecond.magnetic,
            s.firstVelocity,s.secondVelocity};
}
Vec3 apply(const Vec3& v,int p){
    return {(p&1)?-v.x:v.x,(p&2)?-v.y:v.y,(p&4)?-v.z:v.z};
}
const char* tag(int p){
    static char s[8];
    s[0]=(p&1)?'-':'+'; s[1]='x'; s[2]=(p&2)?'-':'+'; s[3]='y';
    s[4]=(p&4)?'-':'+'; s[5]='z'; s[6]=0; return s;
}
void table(const char* nm,const Vec3& a,const Vec3& b){
    std::printf("  %s\n    1 = (%+.6e, %+.6e, %+.6e)\n"
                "    2 = (%+.6e, %+.6e, %+.6e)\n",nm,a.x,a.y,a.z,b.x,b.y,b.z);
    std::printf("   ");
    for(int p=0;p<8;++p)
        std::printf(" %s %.1e",tag(p),
            (b-apply(a,p)).norm()/std::max(a.norm(),1e-300));
    std::printf("\n");
}
}
int main(){
    const double a=pairBohrRadius(activePair);
    for(double bd:{0.0,0.01}){
        std::printf("=== drift %.3g c, moments in the plane (phi 0, "
                    "theta 40) ===\n",bd);
        const Out o=probe(bd,40.0,0.0,a);
        table("E",o.e1,o.e2); table("B",o.b1,o.b2); table("v",o.v1,o.v2);
        std::printf("\n");
    }
}
