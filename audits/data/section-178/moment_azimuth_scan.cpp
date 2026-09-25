#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
int main(){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=40.0*pi/180;
    
    std::printf("Preparation instant, ortho, drift along x (in the mirror\n"
                "plane).  The MOMENT azimuth phi is scanned.  mu is a\n"
                "PSEUDOvector: under reflection in y=0 it is odd when it lies\n"
                "in the x-z plane (phi=0,180) and even when it lies along y\n"
                "(phi=90).  Mixed phi has no definite parity; |sin 2phi| is\n"
                "the mixing.\n\n");
    std::printf("%8s %14s %14s %14s\n","phi deg","|sin 2phi|",
                "|w1-w2|/w1","|m|/mu_sum");
    for(double phi:{0.0,10.0,30.0,45.0,60.0,80.0,90.0,135.0,180.0}){
        const double p=phi*pi/180;
        const Vec3 drift{0.01*c,0.0,0.0};
        const Vec3 dir{std::sin(th)*std::cos(p),std::sin(th)*std::sin(p),
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
        const double w1=b1.norm()*std::abs(firstCharge/firstMass);
        const double w2=b2.norm()*std::abs(secondCharge/secondMass);
        std::printf("%8.0f %14.4f %14.4e %14.4e\n",phi,
            std::abs(std::sin(2.0*p)),
            std::abs(w1-w2)/std::max(w1,1e-300),
            (s.firstDipole+s.secondDipole).norm()
            /(firstMagneticMoment+secondMagneticMoment));
    }
}
