// Audit 214: four seeds (1, 8, 17, 20) keep an ortho lifetime ~9e-04
// above the cos-delta curve at both steps, and seed 13 sits below it.
// If they terminated at a different radius their total radiated energy
// <P> * t would differ from the rest.  Both channels, s_max = 0.30.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    std::printf("%6s %4s %20s %16s %18s\n",
                "seed","ch","lifetime [s]","<P> [W]","E=<P>t [J]");
    for(std::uint64_t seed:{1ull,8ull,17ull,20ull,13ull,
                            2ull,3ull,9ull,21ull})
        for(int ch=1;ch<=2;++ch){
            const CremCollapseEstimate e=
                estimateCremCollapse(seed,ch,900.0);
            std::printf("%6llu %4s %20.12e %16.8e %18.10e\n",
                (unsigned long long)seed,ch==1?"para":"orth",
                e.lifetimeSeconds,e.meanRadiatedPowerWatts,
                e.meanRadiatedPowerWatts*e.lifetimeSeconds);
        }
}
