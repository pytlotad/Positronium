// Point 1 / audit 213: is 211's para-ortho split step-independent?
// Same three seeds at three values of maximumJumpParameter (s_max),
// chosen so that seed 1 (largest ratio in 211) and seed 5 (smallest)
// are both present.  What matters is the RATIO column: the absolute
// times are known to move with s_max by ~0.4%.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    const char* tag=argc>1?argv[1]:"?";
    // Run twice: first with {1,2,5}, then with {15,8,13} (audit 213f).
    for(std::uint64_t seed:{15ull,8ull,13ull}){
        const CremCollapseEstimate p=estimateCremCollapse(seed,1,600.0);
        const CremCollapseEstimate o=estimateCremCollapse(seed,2,600.0);
        const bool ok=
            p.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&o.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&std::isfinite(p.lifetimeSeconds)
            &&std::isfinite(o.lifetimeSeconds);
        if(!ok){ std::printf("%8s %6llu censored\n",tag,
            (unsigned long long)seed); continue; }
        std::printf("%8s %6llu %20.12e %20.12e %14.6e\n",tag,
            (unsigned long long)seed,p.lifetimeSeconds,o.lifetimeSeconds,
            o.lifetimeSeconds/p.lifetimeSeconds-1.0);
    }
}
