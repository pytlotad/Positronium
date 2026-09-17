// Lab lifetime against the start angle over one free-draw ensemble (audit
// section 95).  The start angle is read from the prepared state exactly as
// estimateCremCollapse prepares it; "prep" mode prints it without integrating.
//
// Section 95 run, then tools/lifetime_vs_start_angle.py <dir>:
//   for w in 0 1 2 3; do /tmp/step2 run $((60000+w*20)) 20 300 > step2_$w.out & done; wait
// NOTE: section 95 printed lifetimeSeconds (proper clock); lab and proper
// differ by 1e-11 relative (section 97c).
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/lifetime_vs_start_angle.cpp -o /tmp/step2 $(root-config --libs)

// Step 2: lifetime and photon count against the start angle, one ensemble.
// The start angle is read from the PREPARED state exactly as
// estimateCremCollapse prepares it: simulate(splitMix64(masterSeed+index),
// ...) with a 1e-24 s window, first frame.  No integration, no diagnostics.
// argv: mode(prep|run) baseSeed count budget
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
static double startCos(unsigned long long masterSeed){
  SimulationOptions o; o.frameCount=2; o.observationTime=1.0e-24;
  const SimulationResult r=simulate(splitMix64(masterSeed),1,o);
  if(r.frames.empty()) return std::numeric_limits<double>::quiet_NaN();
  const Vec3 a=r.frames.front().firstDipole, b=r.frames.front().secondDipole;
  return dot(a,b)/(a.norm()*b.norm());
}
int main(int argc,char** argv){
  // Measures the PHOTON CASCADE, so it selects the quantized channel
  // explicitly: production stopped defaulting to it in audit section 109,
  // and without this the tool would silently measure the continuous
  // inspiral instead and stop being comparable with sections 94-108.
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  // The collapse-transit run of audit section 108 would only double the cost.
  gMeasureCollapseTransit=false;
  const bool run=argc>1&&std::strcmp(argv[1],"run")==0;
  const unsigned long long base=argc>2?strtoull(argv[2],nullptr,10):60000ULL;
  const int count=argc>3?atoi(argv[3]):20;
  const double budget=argc>4?atof(argv[4]):300.0;
  for(int i=0;i<count;++i){
    const unsigned long long seed=base+(unsigned long long)i;
    const double c0=startCos(seed);
    if(!run){ std::printf("PREP seed %llu cos0 %.6f\n",seed,c0); continue; }
    const auto r=runCremCollapseExperiment(seed,1,1,budget);
    if(r.empty()){ std::printf("RESULT seed %llu cos0 %.6f EMPTY\n",seed,c0); continue; }
    std::printf("RESULT seed %llu cos0 %.6f t_ps %.9f outcome %d photons %zu\n",
      seed,c0,r[0].lifetimeSeconds*1e12,(int)r[0].calibrationOutcome,
      r[0].labFramePhotons.size());
    std::fflush(stdout);
  }
  if(run) std::printf("DONE\n");
}
