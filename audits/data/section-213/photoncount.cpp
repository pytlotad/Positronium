// Audit 213: the para column takes discrete values.  Is the quantum one
// stochastic photon?  Seeds chosen to span the observed levels
// (3.0659, 3.0708, 3.0745, 3.0825 e-11 s), para channel only.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("%6s %20s %8s %16s\n","seed","para [s]","nphoton","<P> [W]");
    for(std::uint64_t seed:{4ull,1ull,15ull,22ull,16ull,10ull,13ull}){
        const CremCollapseEstimate p=estimateCremCollapse(seed,1,600.0);
        std::printf("%6llu %20.12e %8zu %16.8e\n",(unsigned long long)seed,
            p.lifetimeSeconds,p.labFramePhotons.size(),
            p.meanRadiatedPowerWatts);
    }
}
