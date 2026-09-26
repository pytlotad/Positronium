// Audit 215 / point 2: WHERE along the inspiral the para-ortho
// difference accumulates, using the orbit-averaged quantities 212
// introduced instead of 201d's instantaneous ones.
//
// The integrand of the collapse time is the energy lost per orbit over
// the orbital period.  crem_collapse.hpp:3213 sets
//   deltaEnergyPerOrbit = -run.finalState.orbitalRadiatedEnergy/mu
// -- the E1 far-zone flux with M1 already excluded, and NOT the
// background-subtracted conservativeParticleEnergy difference that the
// measuredDelta lambda just above it computes.  measuredDelta's result
// is used for the ANGULAR MOMENTUM component only.  A first version of
// this probe used measuredDelta for both and came out a factor 2.06
// low in energy while matching dL/L to six digits; that is what pinned
// the substitution down.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <cstdlib>

int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double deltaDeg=argc>1?std::atof(argv[1]):60.0;
    const double reducedMass=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/reducedMass;
    const double a=pairBohrRadius({electron,positron});
    const double th=deltaDeg*3.14159265358979324/180.0;
    const Vec3 unit{std::sin(th),0.0,std::cos(th)};   // delta from +z
    const Vec3 mu1=unit*firstMagneticMoment;
    const auto lz=[](const State& s){
        return cross(s.firstPosition-s.secondPosition,
                     s.firstVelocity-s.secondVelocity).z; };
    std::printf("# delta = %.2f deg from the orbit normal\n",deltaDeg);
    std::printf("# r/a_pair  P_para  P_ortho  dP/P  null  E_para  E_ortho"
                "  dEoverE_para  dLoverL_para\n");
    for(double f:{1.0,0.7,0.5,0.35,0.25,0.15,0.10,0.07,0.05,0.035,
                  0.025,0.015,0.010,0.007,0.005,0.0035,0.0025,0.0018}){
        const double r=f*a;
        const OsculatingElements el{-k/(2.0*r),std::sqrt(k*r)};
        const double period=osculatingPeriod(el.specificEnergy,k);
        double power[3]={0,0,0},ener[3]={0,0,0},dEoE=0.0,dLoL=0.0;
        bool ok=true;
        for(int ch=0;ch<3;++ch){    // ch 2 duplicates ch 0: a null
            const Vec3 mu2=unit*((ch==1?-1.0:+1.0)*secondMagneticMoment);
            const State start=osculatingPeriapsisState(
                el,k,mu1,mu2,Vec3{0.0,0.0,1.0},Vec3{1.0,0.0,0.0},0.0);
            SimulationOptions on; on.frameCount=64;
            on.radiatedEnergyBookkeeping=true;
            const MechanicalTrajectoryResult run=
                runMechanicalTrajectory(start,period,nuclearCutoff,on,
                                        gRadiationReactionModel);
            if(!isFinite(run.finalState)||!(run.elapsedTime>0.0)){
                ok=false; break; }
            power[ch]=run.finalState.orbitalRadiatedEnergy/run.elapsedTime;
            ener[ch]=conservativeParticleEnergy(start);
            if(ch!=0) continue;
            // Validation against the estimator's own checkpoint 0, which
            // prints dE/E and dL/L: the angular-momentum channel IS
            // background-subtracted, so its background run is kept.
            dEoE=-(run.finalState.orbitalRadiatedEnergy/reducedMass)
                 /(el.specificEnergy);
            SimulationOptions off=on;
            off.radiatedEnergyBookkeeping=false;
            const MechanicalTrajectoryResult bg=
                runMechanicalTrajectory(start,period,nuclearCutoff,off,
                    ChargeRadiationReactionModel::disabled);
            dLoL=isFinite(bg.finalState)
                ?((lz(run.finalState)-lz(start))
                  -(lz(bg.finalState)-lz(start)))
                 /el.specificAngularMomentum:0.0;
        }
        if(!ok){ std::printf("%10.5f failed\n",f); continue; }
        std::printf("%10.5f %16.8e %16.8e %14.6e %11.3e %20.12e %20.12e"
                    " %13.5e %13.5e\n",f,power[0],power[1],
            power[0]!=0.0?(power[1]-power[0])/power[0]:0.0,
            power[0]!=0.0?(power[2]-power[0])/power[0]:0.0,
            ener[0],ener[1],dEoE,dLoL);
    }
}
