// Does the LAB-frame dipole energy drift because the force does not
// differentiate its velocity dependence?  Measured in units of k/r0.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
int main(){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double period=2.0*pi*r0/speed;
    const double unit=pairCoulombStrength/r0;          // k/r0
    const double tilt=40.0*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    std::printf("k/r0 = %.6e J = %.6f eV,  beta = %.6e\n\n",
                unit,unit/eCharge,speed/c);
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity={0,speed*secondMass/total,0};
    s.secondVelocity={0,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*secondMagneticMoment;
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-11; acc.maximumDepth=22;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.useRetardedExternalForces=false;   // instantaneous: an exact ledger
    ClassicalTrajectoryEngine engine(s,acc);
    const auto uLab=[&](const State& x){
        return pairDipoleInteractionEnergy(
            x.firstPosition-x.secondPosition,x.firstDipole,x.secondDipole); };
    const auto uProper=[&](const State& x){
        return pairDipoleInteractionEnergy(
            x.firstPosition-x.secondPosition,
            x.firstProperDipole,x.secondProperDipole); };
    double lo=1e300,hi=-1e300,loP=1e300,hiP=-1e300;
    double loE=1e300,hiE=-1e300;
    double properDrift=0.0;
    const Vec3 p0=s.firstProperDipole;
    for(int i=0;i<=4000;++i){
        if(i&&!engine.advance(s,period/4000)) break;
        const double a=uLab(s), b=uProper(s);
        const double e=conservedParticleEnergy(s);
        lo=std::min(lo,a); hi=std::max(hi,a);
        loP=std::min(loP,b); hiP=std::max(hiP,b);
        loE=std::min(loE,e); hiE=std::max(hiE,e);
        properDrift=std::max(properDrift,
            (s.firstProperDipole-p0).norm()/p0.norm());
    }
    std::printf("over one orbit, instantaneous forces, reaction off:\n");
    std::printf("  U_dd with LAB moments    range %.6e J = %.6e k/r0\n",
                hi-lo,(hi-lo)/unit);
    std::printf("  U_dd with PROPER moments range %.6e J = %.6e k/r0\n",
                hiP-loP,(hiP-loP)/unit);
    std::printf("  conservedParticleEnergy  range %.6e J = %.6e k/r0\n",
                hiE-loE,(hiE-loE)/unit);
    std::printf("  proper moment rotated by %.6e of itself "
                "(precession is ON here)\n",properDrift);
    std::printf("  |U_dd| itself            %.6e J = %.6e k/r0\n",
                std::abs(uLab(s)),std::abs(uLab(s))/unit);
}
