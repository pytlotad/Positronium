// Audit 226 control: with useRetardedExternalForces = false, is the
// DIPOLE force still acting?  226's drift is the same at tilt 0 and
// tilt 90 although U_dd differs by 2300x, which is what one would see
// either if the missing velocity term contributes nothing or if the
// dipole force is simply absent in this mode.  The two are told apart
// by scaling the moments down and looking at the TRAJECTORY.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 unit{1.0,0.0,0.0};            // tilt 90, in plane
    std::printf("%10s %8s %20s %20s\n",
                "moments","retarded","final |r| [m]","d|r| vs full");
    for(int ret=0;ret<2;++ret){
        double reference=0.0;
        for(double s0:{1.0,1.0e-6}){
            const State start=osculatingPeriapsisState(
                el,k,unit*(firstMagneticMoment*s0),
                unit*(secondMagneticMoment*s0),
                Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=1.0e-12; acc.maximumDepth=26;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            acc.useRetardedExternalForces=(ret!=0);
            ClassicalTrajectoryEngine engine(start,acc);
            State s=start; bool ok=true;
            for(int i=0;i<512&&ok;++i) ok=engine.advance(s,period/512);
            if(!ok){ std::printf("%10.0e %8d  failed\n",s0,ret);
                     continue; }
            const double r=(s.firstPosition-s.secondPosition).norm();
            if(s0==1.0) reference=r;
            std::printf("%10.0e %8d %20.12e %20.4e\n",
                        s0,ret,r,r-reference);
        }
    }
    std::printf("\n# A nonzero difference means the dipole force acts.\n");
}
