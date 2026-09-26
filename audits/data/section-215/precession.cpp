// Audit 215: the probe holds delta fixed along the whole inspiral.  The
// estimator does not -- its own debug prints "coupled spin-orbit half
// substeps=69 maxAngle=0.05" per checkpoint.  How far does the moment
// actually turn in ONE orbit, and how far does delta itself move?
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=pairBohrRadius({electron,positron});
    const double th=60.0*3.14159265358979324/180.0;
    const Vec3 unit{std::sin(th),0.0,std::cos(th)};
    const auto ang=[](const Vec3& u,const Vec3& v){
        return std::acos(std::max(-1.0,std::min(1.0,
            dot(u,v)/(u.norm()*v.norm()))))*180.0/3.14159265358979324; };
    std::printf("%10s %14s %14s %14s\n","r/a_pair","turn/orbit deg",
                "delta_end deg","orbits to 90deg");
    for(double f:{1.0,0.5,0.25,0.15,0.10,0.05}){
        const double r=f*a;
        const OsculatingElements el{-k/(2.0*r),std::sqrt(k*r)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        const State start=osculatingPeriapsisState(el,k,
            unit*firstMagneticMoment,unit*secondMagneticMoment,
            Vec3{0,0,1},Vec3{1,0,0},0.0);
        SimulationOptions o; o.frameCount=8;
        o.radiatedEnergyBookkeeping=false;
        const MechanicalTrajectoryResult run=runMechanicalTrajectory(
            start,period,nuclearCutoff,o,gRadiationReactionModel);
        if(!isFinite(run.finalState)){ std::printf("%10.4f failed\n",f);
            continue; }
        const double turn=ang(start.firstDipole,run.finalState.firstDipole);
        std::printf("%10.4f %14.6f %14.4f %14.3g\n",f,turn,
            ang(run.finalState.firstDipole,Vec3{0,0,1}),
            turn>0.0?90.0/turn:0.0);
    }
}
