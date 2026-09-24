// Is the exchange parity that forces w1 = w2 in ortho EXACT, or an artifact
// of the symmetric preparation?  Break the preparation's symmetry in ways
// that leave mu2 = -mu1 intact and see whether w1 = w2 survives.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
struct Result { double worstM, worstRate, worstOmega; bool ok; };
Result run(double tiltDeg,double eccentricity,double comBeta,
           const Vec3& externalField,int orbits=1){
    const double total=firstMass+secondMass;
    const double a=pairBohrRadius(activePair);
    // periapsis start of an ellipse with the given eccentricity
    const double r0=a*(1.0-eccentricity);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*a)
        *(1.0+eccentricity)/(1.0-eccentricity));
    const double period=2.0*pi*std::sqrt(pairReducedMass*a*a*a
        /pairCoulombStrength);
    const double tilt=tiltDeg*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    const Vec3 com{comBeta*c,0.0,0.0};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity=Vec3{0,speed*secondMass/total,0}+com;
    s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+com;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);   // ORTHO, exactly
    synchronizeCovariantDipoles(s);
    gExternalMagneticField=externalField;
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-9; acc.maximumDepth=20;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    const double scale=firstMagneticMoment+secondMagneticMoment;
    Result out{0,0,0,true};
    const int steps=500*orbits;
    for(int i=0;i<=steps;++i){
        if(i&&!engine.advance(s,period*orbits/steps)){ out.ok=false; break; }
        out.worstM=std::max(out.worstM,
            (s.firstDipole+s.secondDipole).norm()/scale);
        const DipoleDerivatives d=
            thomasBmtDipoleDerivatives(s,engine.history());
        out.worstRate=std::max(out.worstRate,(d.first+d.second).norm()
            /std::max(d.first.norm()+d.second.norm(),1e-300));
        const double w1=d.first.norm()
            /std::max(s.firstDipole.norm(),1e-300);
        const double w2=d.second.norm()
            /std::max(s.secondDipole.norm(),1e-300);
        out.worstOmega=std::max(out.worstOmega,
            std::abs(w1-w2)/std::max(w1,1e-300));
    }
    gExternalMagneticField=Vec3{};
    return out;
}
}
int main(){
    std::printf("ORTHO prepared exactly antiparallel; worst value over the "
                "orbit\n");
    std::printf("%-34s %13s %13s %13s\n","preparation","|w1-w2|/|w1|",
                "|m|/mu_sum","|dm/dt| norm");
    struct Case { const char* name; double tilt,ecc,com; Vec3 field; };
    const Case cases[]={
        {"baseline: tilt 40, circular",        40.0,0.0,0.0,{}},
        {"tilt 0 (moments in the plane)",       0.0,0.0,0.0,{}},
        {"tilt 90 (moments along L)",          90.0,0.0,0.0,{}},
        {"tilt 17 (generic)",                  17.0,0.0,0.0,{}},
        {"eccentric e=0.3",                    40.0,0.3,0.0,{}},
        {"eccentric e=0.7",                    40.0,0.7,0.0,{}},
        {"centre of mass beta=0.01",           40.0,0.0,0.01,{}},
        {"centre of mass beta=0.10",           40.0,0.0,0.10,{}},
        {"POSITIVE CONTROL: external B 50 uT", 40.0,0.0,0.0,{0,0,50e-6}},
    };
    for(const Case& c:cases){
        const Result r=run(c.tilt,c.ecc,c.com,c.field);
        if(!r.ok){ std::printf("%-34s   (engine stopped)\n",c.name);
                      std::fflush(stdout); continue; }
        std::printf("%-34s %13.4e %13.4e %13.4e\n",
            c.name,r.worstOmega,r.worstM,r.worstRate);
        std::fflush(stdout);
    }
}
