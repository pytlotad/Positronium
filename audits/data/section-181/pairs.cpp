// 181e: the moment-magnitude condition across three pairs
// -- the residual reproduces |mu1|/|mu2| - 1
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        pairs.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
namespace { Vec3 axial(const Vec3& v){ return {-v.x,v.y,-v.z}; } }
int main(){
    auto show=[&](const char* nm){
        const double total=firstMass+secondMass;
        const double r0=pairBohrRadius(activePair);
        const double speed=std::sqrt(pairCoulombStrength
            /(pairReducedMass*r0));
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
        const Vec3 b1=f.atFirst.magnetic,b2=f.atSecond.magnetic;
        std::printf("%-22s |B2|/|B1| = %12.5e   mu1/mu2 = %12.5e\n"
                    "%-22s residual  = %12.5e   ratio-1 = %12.5e\n",
            nm,b2.norm()/b1.norm(),firstMagneticMoment/secondMagneticMoment,
            "",(b2-axial(b1)).norm()/b1.norm(),
            firstMagneticMoment/secondMagneticMoment-1.0);
    };
    show("electron + positron");
    applyPair(ParticlePair{electron,antimuon});
    show("electron + antimuon");
    applyPair(ParticlePair{proton,antiproton});
    show("proton + antiproton");
}
