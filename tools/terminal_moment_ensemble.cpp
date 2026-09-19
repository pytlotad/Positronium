// What the dynamics does to the ensemble averages of audit 153 (section 154).
//
// Section 153 measured prepared states; this runs the production path to
// annihilation and marks each trajectory on stderr so the terminal moments can
// be read from the last CREM_DEBUG_SPINPAIR line of its block.  The reaction
// model is left at the production default (continuous Landau-Lifshitz, photons
// off since audit 109), which is both the configuration the question is about
// and the fast one.
//
// Usage: terminal_moment_ensemble <phenomenon> [count] [seed base] [budget s]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/terminal_moment_ensemble.cpp -o /tmp/terminal $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const int phenomenon=argc>1?atoi(argv[1]):1;
    const int count=argc>2?atoi(argv[2]):12;
    const unsigned long long base=argc>3?strtoull(argv[3],nullptr,10):42ULL;
    const double budget=argc>4?atof(argv[4]):300.0;
    gMeasureCollapseTransit=false;
    for(int index=0;index<count;++index) {
        const unsigned long long seed=base+static_cast<unsigned long long>(index);
        std::fprintf(stderr,"TRAJ %d seed %llu phenomenon %d\n",
                     index,seed,phenomenon);
        std::fflush(stderr);
        const auto result=runCremCollapseExperiment(seed,phenomenon,1,budget);
        if(result.empty()) { std::printf("RESULT seed %llu EMPTY\n",seed); continue; }
        std::printf("RESULT seed %llu phenomenon %d t_ps %.9g outcome %d "
                    "stop %d\n",seed,phenomenon,
                    result[0].lifetimeSeconds*1e12,
                    static_cast<int>(result[0].calibrationOutcome),
                    static_cast<int>(result[0].stopCause));
        std::fflush(stdout);
    }
    std::printf("DONE\n");
}
