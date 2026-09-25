// 178f item 3: do not guess the symmetry operation -- measure the relation
// between the two B_eff vectors and see which one it actually is.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
void report(const char* name,const Vec3& drift,double phiDeg){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=40.0*pi/180, ph=phiDeg*pi/180;
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
    const Vec3 b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                          firstGFactor);
    const Vec3 b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                          secondGFactor);
    // candidate relations between b2 and b1
    const Vec3 mirrorY{b1.x,-b1.y,b1.z};      // reflect in y=0
    const Vec3 mirrorX{-b1.x,b1.y,b1.z};      // reflect in x=0
    const Vec3 negMirrorY{-b1.x,b1.y,-b1.z};  // -(reflect in y=0)
    const auto rel=[&](const Vec3& cand){
        return (b2-cand).norm()/std::max(b1.norm(),1e-300); };
    std::printf("%-30s\n",name);
    std::printf("   B1 = (%+.6e, %+.6e, %+.6e)\n",b1.x,b1.y,b1.z);
    std::printf("   B2 = (%+.6e, %+.6e, %+.6e)\n",b2.x,b2.y,b2.z);
    std::printf("   |B2|/|B1|-1            %+.3e\n",b2.norm()/b1.norm()-1.0);
    std::printf("   B2 vs  mirror_y(B1)     %.3e\n",rel(mirrorY));
    std::printf("   B2 vs  mirror_x(B1)     %.3e\n",rel(mirrorX));
    std::printf("   B2 vs -mirror_y(B1)     %.3e\n",rel(negMirrorY));
    std::printf("   B2 vs  B1               %.3e\n",rel(b1));
    std::printf("   B2 vs -B1               %.3e\n\n",rel(b1*-1.0));
}
}
int main(){
    report("no drift, phi=0",           Vec3{},0.0);
    report("drift 0.01c along x, phi=0",Vec3{0.01*c,0,0},0.0);
    report("drift 0.01c along y, phi=0",Vec3{0,0.01*c,0},0.0);
    report("drift 0.01c along x, phi=45",Vec3{0.01*c,0,0},45.0);
}
