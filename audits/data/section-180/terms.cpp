// 180c: does B2 = -M(B1) hold separately in each Thomas-BMT term?
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        terms.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
Vec3 bmtTerm(const Vec3& v,const ElectromagneticField& f,double g,int term){
    const double a=0.5*(g-2.0), gam=gamma(v);
    const Vec3 beta=v/c;
    if(term==0) return f.magnetic*(a+1.0/gam);
    if(term==1) return beta*(-a*gam/(gam+1.0)*dot(beta,f.magnetic));
    return cross(beta,f.electric)*(-(a+1.0/(gam+1.0))/c);
}
struct Pair { Vec3 b1,b2; };
Pair fields(double phiDeg,const Vec3& drift,double r0,int term=-1){
    const double total=firstMass+secondMass;
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
    if(term<0)
        return {thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                        firstGFactor),
                thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                        secondGFactor)};
    return {bmtTerm(s.firstVelocity,f.atFirst,firstGFactor,term),
            bmtTerm(s.secondVelocity,f.atSecond,secondGFactor,term)};
}
double residual(const Pair& p){
    const Vec3 target{-p.b1.x,p.b1.y,-p.b1.z};   // -M_y(B1)
    return (p.b2-target).norm()/std::max(p.b1.norm(),1e-300);
}
}
int main(){
    const double a=pairBohrRadius(activePair);
    const Vec3 dx{0.01*c,0,0};
    std::printf("ITEM 3: do the NORMS even move with the moment azimuth?\n");
    std::printf("%8s %18s %18s %13s %13s\n","phi","|B1|","|B2|",
                "|B2|/|B1|-1","identity res");
    for(double p:{0.0,15.0,30.0,45.0,60.0,90.0}){
        const Pair q=fields(p,dx,a);
        std::printf("%8.0f %18.12e %18.12e %13.3e %13.3e\n",p,
            q.b1.norm(),q.b2.norm(),q.b2.norm()/q.b1.norm()-1.0,residual(q));
    }
    std::printf("\nITEM 1: does B2 = -M(B1) hold TERM BY TERM in the BMT\n"
                "decomposition?  (drift in the mirror, phi = 0)\n");
    const char* names[3]={"magnetic  (a+1/gamma)B",
                          "longitudinal beta(beta.B)",
                          "motional  beta x E"};
    for(int t=0;t<3;++t){
        const Pair q=fields(0.0,dx,a);
        const Pair qt=fields(0.0,dx,a,t);
        std::printf("   %-26s |term|/|B_eff| %8.3e   residual %10.3e\n",
            names[t],qt.b1.norm()/std::max(q.b1.norm(),1e-300),residual(qt));
    }
    std::printf("\nITEM 2: is the identity residual's scale general?\n");
    std::printf("%10s %10s %15s\n","r/a_pair","beta_drift",
                "residual at phi=45");
    for(double rf:{1.0,0.1}) for(double bd:{0.003,0.01,0.03}){
        const Pair q=fields(45.0,Vec3{bd*c,0,0},rf*a);
        std::printf("%10.3g %10.3g %15.4e\n",rf,bd,residual(q));
    }
}
