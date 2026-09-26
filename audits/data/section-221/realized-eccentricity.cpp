// Audit 221: is the imposed eccentricity realized?  Below e = beta^2
// the k=2 harmonic does not move off the circular floor AT ALL -- not
// even by the sqrt(2) that adding an equal contribution in quadrature
// would give.  The obvious explanation is that the orbit does not
// actually carry the eccentricity asked of it, because the model's own
// non-Keplerian distortion is itself of order beta^2.  Measured here
// directly from the separation over one radial period.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double A=pairBohrRadius({electron,positron});
    std::printf("%9s %11s %12s %12s %10s\n",
                "a/a_pair","beta_rel^2","e imposed","e realized","ratio");
    for(double f:{1.0,0.01}){
        const double a=f*A;
        const double b2=pairCoulombStrength/(mu*c*c*a);
        for(double scale:{0.0,0.25,0.5,1.0,2.0,4.0,10.0}){
            const double e=scale*b2/2.0;
            const OsculatingElements el{-k/(2.0*a),
                std::sqrt(k*a*(1.0-e*e))};
            const double period=osculatingPeriod(el.specificEnergy,k);
            const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*1.0e-6);
            const State start=osculatingPeriapsisState(el,k,m1,
                m1*(secondMagneticMoment/firstMagneticMoment),
                Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            ClassicalTrajectoryEngine engine(start,acc);
            State s=start;
            double lo=1.0e300,hi=0.0; bool ok=true;
            for(int i=0;i<512&&ok;++i){
                ok=engine.advance(s,period/512);
                const double r=(s.firstPosition-s.secondPosition).norm();
                lo=std::min(lo,r); hi=std::max(hi,r);
            }
            if(!ok){ std::printf("%9.4f  failed\n",f); continue; }
            const double realized=(hi-lo)/(hi+lo);
            std::printf("%9.4f %11.4e %12.5e %12.5e %10.3f\n",
                f,b2,e,realized,e>0.0?realized/e:realized/b2);
        }
    }
    std::printf("\n# ratio column is realized/imposed, except on the\n");
    std::printf("# e = 0 row where it is realized/beta^2.\n");
}
