// Is the exact w1 = w2 at the preparation instant protected by a residual
// reflection in the y = 0 plane?  At t = 0 the positions lie along x and the
// orbital velocities along +-y, so y -> -y exchanges the two orbital
// velocities.  If that is the protection, taking the MOMENTS out of the
// x-z plane, or turning the drift out of x, must break it.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
struct Out { double w1,w2,rel,m; };
Out probe(double thetaDeg,double phiDeg,const Vec3& drift,double phaseDeg){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thetaDeg*pi/180, ph=phiDeg*pi/180, ps=phaseDeg*pi/180;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),
                   std::cos(th)};
    // orbital phase: rotate the separation and the orbital velocity together
    const Vec3 rhat{std::cos(ps),std::sin(ps),0.0};
    const Vec3 vhat{-std::sin(ps),std::cos(ps),0.0};
    State s{};
    s.firstPosition=rhat*(r0*secondMass/total);
    s.secondPosition=rhat*(-r0*firstMass/total);
    s.firstVelocity=vhat*(speed*secondMass/total)+drift;
    s.secondVelocity=vhat*(-speed*firstMass/total)+drift;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                          firstGFactor);
    const Vec3 b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                          secondGFactor);
    const double w1=b1.norm()*std::abs(firstCharge/firstMass);
    const double w2=b2.norm()*std::abs(secondCharge/secondMass);
    return {w1,w2,std::abs(w1-w2)/std::max(w1,1e-300),
            (s.firstDipole+s.secondDipole).norm()
            /(firstMagneticMoment+secondMagneticMoment)};
}
}
int main(){
    const double b=0.01*c;
    std::printf("ORTHO, a_pair, preparation instant, drift 0.01c unless "
                "stated\n");
    std::printf("%-46s %12s %12s\n","configuration","|w1-w2|/w1","|m|/mu_sum");
    struct Case { const char* name; double th,ph,phase; Vec3 drift; };
    const Case cases[]={
      {"177's case: moments in x-z, drift along x",  40,  0,  0,{0.01*c,0,0}},
      {"moments out of plane, phi = 30",             40, 30,  0,{0.01*c,0,0}},
      {"moments out of plane, phi = 90",             40, 90,  0,{0.01*c,0,0}},
      {"moments along y exactly (theta 90, phi 90)", 90, 90,  0,{0.01*c,0,0}},
      {"drift along y (parallel to v_orb)",          40,  0,  0,{0,0.01*c,0}},
      {"drift along z (normal to the orbit)",        40,  0,  0,{0,0,0.01*c}},
      {"drift skew (x+y+z)/sqrt3",                   40,  0,  0,
                    {0.01*c/1.7320508,0.01*c/1.7320508,0.01*c/1.7320508}},
      {"orbital phase 30 deg, drift along x",        40,  0, 30,{0.01*c,0,0}},
      {"no drift, phi = 30 (control)",               40, 30,  0,{}},
    };
    for(const Case& k:cases){
        const Out o=probe(k.th,k.ph,k.drift,k.phase);
        std::printf("%-46s %12.4e %12.4e\n",k.name,o.rel,o.m);
    }
    (void)b;
}
