#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("%6s %20s %20s %14s\n","seed","para [s]","ortho [s]",
                "ortho/para-1");
    for(std::uint64_t seed=1;seed<=24;++seed){
        const CremCollapseEstimate p=estimateCremCollapse(seed,1,300.0);
        const CremCollapseEstimate o=estimateCremCollapse(seed,2,300.0);
        const bool ok=
            p.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&o.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&std::isfinite(p.lifetimeSeconds)
            &&std::isfinite(o.lifetimeSeconds);
        if(!ok){ std::printf("%6llu %20s %20s %14s\n",
            (unsigned long long)seed,"censored","censored","-"); continue; }
        std::printf("%6llu %20.12e %20.12e %14.6e\n",
            (unsigned long long)seed,p.lifetimeSeconds,o.lifetimeSeconds,
            o.lifetimeSeconds/p.lifetimeSeconds-1.0);
    }
}
