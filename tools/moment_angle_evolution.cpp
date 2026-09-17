// Does the inspiral select the mutual moment angle when none is imposed?
// (audit section 94).  One trajectory per call, TRAJ markers on stderr, the
// angle at every spin-transport substep through CREM_DEBUG_ALIGN.
//
// Section 94 run, four processes, then tools/moment_angle_evolution.py <dir>:
//   for w in 0 1 2 3; do CREM_DEBUG_ALIGN=100000000 /tmp/cosevo $((50000+w*100)) \
//     10 300 > cosevo_$w.out 2> cosevo_$w.err & done; wait
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/moment_angle_evolution.cpp -o /tmp/cosevo $(root-config --libs)

// Does the dynamics select mutual moment angles when none is imposed?
// One trajectory per runCremCollapseExperiment call, so the CREM_DEBUG_ALIGN
// lines of different trajectories cannot interleave; a TRAJ marker on stderr
// separates them.  Result lines on stdout carry lifetime, outcome, photons.
// argv: baseSeed count budget
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  // Measures the PHOTON CASCADE, so it selects the quantized channel
  // explicitly: production stopped defaulting to it in audit section 109,
  // and without this the tool would silently measure the continuous
  // inspiral instead and stop being comparable with sections 94-108.
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  // The collapse-transit run of audit section 108 would only double the cost.
  gMeasureCollapseTransit=false;
  const unsigned long long base=argc>1?strtoull(argv[1],nullptr,10):1000ULL;
  const int count=argc>2?atoi(argv[2]):10;
  const double budget=argc>3?atof(argv[3]):300.0;
  std::printf("gSpinQuantization=%d (0 = kat losowany swobodnie)\n",
              gSpinQuantization?1:0);
  for(int i=0;i<count;++i){
    const unsigned long long seed=base+(unsigned long long)i;
    std::fprintf(stderr,"TRAJ %d seed %llu\n",i,seed);
    std::fflush(stderr);
    const auto r=runCremCollapseExperiment(seed,1,1,budget);
    if(r.empty()){ std::printf("RESULT seed %llu EMPTY\n",seed); continue; }
    std::printf("RESULT seed %llu t_ps %.9g outcome %d photons %zu\n",
      seed,r[0].lifetimeSeconds*1e12,(int)r[0].calibrationOutcome,
      r[0].labFramePhotons.size());
    std::fflush(stdout);
  }
  std::printf("DONE\n");
}
