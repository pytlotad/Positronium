// 176g's three open items, measured.
//  (1) is |w1-w2|/|w1| large because one precession stalls?  Print w1 and w2.
//  (2) what is the actual beta dependence?  Scan it.
//  (3) is the motional v x E term the mechanism?  Rebuild B_eff without it.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
// thomasBmtEffectiveField, with each piece switchable.
Vec3 bmt(const Vec3& velocity,const ElectromagneticField& f,double g,
         bool keepMagnetic,bool keepLongitudinal,bool keepMotional){
    const double anomaly=0.5*(g-2.0);
    const double gam=gamma(velocity);
    const Vec3 beta=velocity/c;
    Vec3 out{};
    if(keepMagnetic)     out+=f.magnetic*(anomaly+1.0/gam);
    if(keepLongitudinal) out-=beta*(anomaly*gam/(gam+1.0)
                                    *dot(beta,f.magnetic));
    if(keepMotional)     out-=cross(beta,f.electric)
                              *((anomaly+1.0/(gam+1.0))/c);
    return out;
}
struct Rates { double w1,w2; };
Rates rates(const State& s,const StateHistory& h,
            bool mag,bool lon,bool mot){
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    const Vec3 b1=bmt(s.firstVelocity,f.atFirst,firstGFactor,mag,lon,mot);
    const Vec3 b2=bmt(s.secondVelocity,f.atSecond,secondGFactor,mag,lon,mot);
    return {b1.norm()*std::abs(firstCharge/firstMass),
            b2.norm()*std::abs(secondCharge/secondMass)};
}
}
int main(){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double tilt=40.0*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    std::printf("ORTHO, a_pair, tilt 40, at the preparation instant.\n");
    std::printf("Rates in 1/s; 'full' is the production B_eff.\n\n");
    std::printf("%8s %13s %13s %11s | %13s %13s %11s\n","beta_CM",
        "w1 full","w2 full","|dw|/w1","w1 no-motional","w2 no-motional",
        "|dw|/w1");
    for(double b:{0.0,0.001,0.003,0.01,0.03,0.1,0.3}){
        State s{};
        s.firstPosition={r0*secondMass/total,0,0};
        s.secondPosition={-r0*firstMass/total,0,0};
        const Vec3 com{b*c,0,0};
        s.firstVelocity=Vec3{0,speed*secondMass/total,0}+com;
        s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+com;
        s.firstProperDipole=dir*firstMagneticMoment;
        s.secondProperDipole=dir*(-secondMagneticMoment);
        synchronizeCovariantDipoles(s);
        const StateHistory h=causalInitialHistory(s);
        const Rates full=rates(s,h,true,true,true);
        const Rates noMot=rates(s,h,true,true,false);
        std::printf("%8.3g %13.5e %13.5e %11.4f | %13.5e %13.5e %11.4e\n",
            b,full.w1,full.w2,
            std::abs(full.w1-full.w2)/std::max(full.w1,1e-300),
            noMot.w1,noMot.w2,
            std::abs(noMot.w1-noMot.w2)/std::max(noMot.w1,1e-300));
    }
}
