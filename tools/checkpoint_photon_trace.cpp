// One trajectory with both clocks, photons and wall time (audit sections 96-97).
// With CREM_SKIP_CENSUS=1 CREM_PHOTON_LADDER=1 the stderr carries the checkpoint
// census and the photons, which tools/checkpoint_photon_trace.py compares.
//
// Section 96 run (pairs.txt lists "short_seed long_seed" per line):
//   CREM_SKIP_CENSUS=1 CREM_PHOTON_LADDER=1 /tmp/trace <seed> 300 \
//     > trace_<seed>.out 2> trace_<seed>.err
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/checkpoint_photon_trace.cpp -o /tmp/trace $(root-config --libs)

// Both clocks, photons and wall time for one trajectory.  argv: seed budget
#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  const unsigned long long seed=strtoull(argv[1],nullptr,10);
  const double budget=argc>2?atof(argv[2]):900.0;
  const auto t0=std::chrono::steady_clock::now();
  const auto r=runCremCollapseExperiment(seed,1,1,budget);
  const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  std::printf("RESULT seed %llu lab_ps %.9f proper_ps %.9f lab_minus_proper_fs %.6f outcome %d photons %zu wall_s %.1f\n",
    seed,r[0].lifetimeSecondsLab*1e12,r[0].lifetimeSeconds*1e12,
    (r[0].lifetimeSecondsLab-r[0].lifetimeSeconds)*1e15,
    (int)r[0].calibrationOutcome,r[0].labFramePhotons.size(),wall);
}
