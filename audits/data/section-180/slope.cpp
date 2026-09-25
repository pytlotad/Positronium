// 180h/180i: the per-degree slope across radius, polar angle and drift
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        slope.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
Vec3 rotZ(const Vec3& v,double a){
    const double c_=std::cos(a),s_=std::sin(a);
    return {c_*v.x-s_*v.y,s_*v.x+c_*v.y,v.z};
}
Vec3 mirrorPseudo(const Vec3& b,const Vec3& n){ return n*(2.0*dot(b,n))-b; }
double residual(double psiDeg,double bd,double r0,double thDeg){
    const double total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thDeg*pi/180, psi=psiDeg*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    const Vec3 sep=rotZ(Vec3{1,0,0},psi), orb=rotZ(Vec3{0,1,0},psi);
    State s{};
    s.firstPosition=sep*(r0*secondMass/total);
    s.secondPosition=sep*(-r0*firstMass/total);
    s.firstVelocity=orb*(speed*secondMass/total)+Vec3{bd*c,0,0};
    s.secondVelocity=orb*(-speed*firstMass/total)+Vec3{bd*c,0,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                          firstGFactor);
    const Vec3 b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                          secondGFactor);
    return (b2-mirrorPseudo(b1,rotZ(Vec3{0,1,0},psi))).norm()
           /std::max(b1.norm(),1e-300);
}
double slope(double bd,double r0,double th){
    return residual(2.0,bd,r0,th)/2.0;
}
}
int main(){
    const double a=pairBohrRadius(activePair);
    std::printf("drift-free floor of the slope (beta_drift = 0)\n");
    std::printf("%10s %9s %14s %14s\n","r/a_pair","beta_orb","theta 20",
                "theta 70");
    for(double rf:{1.0,0.1,0.01,0.003}){
        const double bo=std::sqrt(pairCoulombStrength
            /(pairReducedMass*rf*a))/c;
        std::printf("%10.3g %9.4f %14.4e %14.4e\n",rf,bo,
                    slope(0.0,rf*a,20.0),slope(0.0,rf*a,70.0));
    }
    std::printf("\nslope minus floor against the speed ratio "
                "(theta = 40)\n");
    std::printf("%10s %10s %12s %12s %12s\n","r/a_pair","beta_dr",
                "ratio","slope-floor","/ratio^2");
    for(double rf:{1.0,0.1,0.01}){
        const double bo=std::sqrt(pairCoulombStrength
            /(pairReducedMass*rf*a))/c;
        const double f0=slope(0.0,rf*a,40.0);
        for(double bd:{0.001,0.003,0.01,0.03}){
            const double q=bd/bo, d=slope(bd,rf*a,40.0)-f0;
            std::printf("%10.3g %10.3g %12.4f %12.4e %12.4e\n",rf,bd,q,d,
                        d/(q*q));
        }
    }
}
