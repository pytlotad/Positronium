// Audit 214: the 24-seed sweep repeated at s_max = 0.075, the step at
// which 213f showed the discreteness of the para column has shrunk by a
// factor 3.5.  Seed range on argv so four chunks can run in parallel.
// The wall-clock budget is generous: a censored row would be a hole in
// the sweep, and at this step each estimate takes minutes.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    const std::uint64_t lo=argc>1?std::strtoull(argv[1],nullptr,10):1;
    const std::uint64_t hi=argc>2?std::strtoull(argv[2],nullptr,10):24;
    for(std::uint64_t seed=lo;seed<=hi;++seed){
        const CremCollapseEstimate p=estimateCremCollapse(seed,1,3000.0);
        const CremCollapseEstimate o=estimateCremCollapse(seed,2,3000.0);
        const bool ok=
            p.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&o.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&std::isfinite(p.lifetimeSeconds)
            &&std::isfinite(o.lifetimeSeconds);
        if(!ok){ std::printf("%6llu censored %d %d\n",
            (unsigned long long)seed,
            static_cast<int>(p.calibrationOutcome),
            static_cast<int>(o.calibrationOutcome)); continue; }
        std::printf("%6llu %20.12e %20.12e %14.6e\n",
            (unsigned long long)seed,p.lifetimeSeconds,o.lifetimeSeconds,
            o.lifetimeSeconds/p.lifetimeSeconds-1.0);
    }
}
