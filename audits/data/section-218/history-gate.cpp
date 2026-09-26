// Audit 218: the branch label, quantitatively.  The history commits a
// node only when at least retentionTime/targetHistoryNodes has elapsed
// since the last committed one (electrodynamics.hpp:5743), with
// retentionTime = max(1e-20, 4 separation/c) and targetHistoryNodes
// = 128; above maximumHistoryNodes = 256 the history is DECIMATED.
// So the retarded grid is pinned to the separation, not to the
// integration step, and refining the step past the gate cannot refine
// it -- while crossing the cap actively coarsens it.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <vector>

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double r0=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*r0),std::sqrt(k*r0)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const State rest=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    const double span=0.25*period;
    const double retention=std::max(1.0e-20,4.0*r0/c);
    const double gate=retention/128.0;
    std::printf("# retentionTime = 4 r/c = %.4e s, gate = /128 ="
                " %.4e s\n",retention,gate);
    std::printf("# span = quarter period = %.4e s\n",span);
    std::printf("# window/gate = %.1f nodes fit in the retarded window\n",
                retention/gate);
    std::printf("\n%7s %12s %10s %9s %16s\n",
                "steps","dt [s]","dt/gate","history","|dP/dt| [W]");
    for(int steps:{200,400,800,1600,3200,6400,12800}){
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-8; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(rest,acc);
        State s=rest;
        bool ok=true;
        for(int i=0;i<steps&&ok;++i) ok=engine.advance(s,span/steps);
        if(!ok){ std::printf("%7d  failed\n",steps); continue; }
        FarFieldSampling sampling; sampling.directionCount=50;
        const FieldFluxRates f=
            electromagneticFieldFluxRates(s,engine.history(),sampling);
        std::printf("%7d %12.4e %10.3f %9zu %16.8e\n",
                    steps,span/steps,(span/steps)/gate,
                    engine.history().size(),f.energy);
    }
    std::printf("\n# A history that stops growing once dt falls under\n");
    std::printf("# the gate, and shrinks once the 256 cap decimates,\n");
    std::printf("# is the branch structure of audit 217f.\n");
}
