// Audit 215: is dP/P converged, or is it a start transient?  The same
// radius over one, two and four periods.  orbitalRadiatedEnergy over
// elapsed time must be independent of the window if the orbit is
// stationary; any drift is the settling artefact crem_collapse warns
// about for the conservativeParticleEnergy channel.
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
    std::printf("%9s %8s %16s %16s %14s\n",
                "r/a_pair","periods","P_para [W]","P_ortho [W]","dP/P");
    for(double f:{0.50,0.35,0.25,0.15}){
        const double r=f*a;
        const OsculatingElements el{-k/(2.0*r),std::sqrt(k*r)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        for(double n:{1.0,2.0,4.0}){
            double p[2]={0,0}; bool ok=true;
            for(int ch=0;ch<2;++ch){
                const State s=osculatingPeriapsisState(el,k,
                    unit*firstMagneticMoment,
                    unit*((ch==0?1.0:-1.0)*secondMagneticMoment),
                    Vec3{0,0,1},Vec3{1,0,0},0.0);
                SimulationOptions o; o.frameCount=64;
                o.radiatedEnergyBookkeeping=true;
                const MechanicalTrajectoryResult run=runMechanicalTrajectory(
                    s,n*period,nuclearCutoff,o,gRadiationReactionModel);
                if(!isFinite(run.finalState)||!(run.elapsedTime>0.0)){
                    ok=false; break; }
                p[ch]=run.finalState.orbitalRadiatedEnergy/run.elapsedTime;
            }
            if(!ok){ std::printf("%9.4f %8.0f failed\n",f,n); continue; }
            std::printf("%9.4f %8.0f %16.8e %16.8e %14.6e\n",
                        f,n,p[0],p[1],(p[1]-p[0])/p[0]);
        }
    }
}
