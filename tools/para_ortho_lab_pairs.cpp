// Para against ortho on the LAB clock (audit sections 98-100).
//   quant  --spin-quantization pairs at a common seed: identical orbit and
//          first moment, cos exactly +1 (para) and -1 (ortho)
//   free   default free draw: start cos and lab lifetime
//   prep   check the pairing premise without integrating
//   tilt   first moment's cos to the orbital normal (section 100b)
//
// Section 98 runs, then tools/para_ortho_lab_pairs.py <dir>:
//   A: for w in 0 1 2 3; do /tmp/pairs quant $((70000+w*5)) 5 900 > poA_$w.out & done; wait
//   B: as A with a binary built by tools/build_variant.sh --max-jump 0.15
//      (poB_*.out, seeds 70000-70009)
//   C: for w in 0 1 2 3; do /tmp/pairs free $((71000+w*10)) 10 900 > poC_$w.out & done; wait
//   tilt: /tmp/pairs tilt 70000 20 > tilt.txt
// Section 99 ablation: same binary with CREM_MAGNETIC_RADIUS_SCALE=0.05,
// CREM_NO_DIPOLE_FORCE=1 or CREM_NO_THOMAS_BACKREACTION=1.  Section 100: a
// binary from tools/build_variant.sh --freeze-spin-orbit, run with
// CREM_FREEZE_SPIN_ORBIT=1.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/para_ortho_lab_pairs.cpp -o /tmp/pairs $(root-config --libs)

// Para against ortho on the LAB clock.
// argv: mode(quant|free|prep) baseSeed count budget
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
static void prepared(unsigned long long seed,int ph,double& c,Vec3& m1){
  SimulationOptions o; o.frameCount=2; o.observationTime=1.0e-24;
  const SimulationResult r=simulate(splitMix64(seed),ph,o);
  m1=r.frames.front().firstDipole;
  const Vec3 m2=r.frames.front().secondDipole;
  c=dot(m1,m2)/(m1.norm()*m2.norm());
}
int main(int argc,char** argv){
  // Measures the PHOTON CASCADE, so it selects the quantized channel
  // explicitly: production stopped defaulting to it in audit section 109,
  // and without this the tool would silently measure the continuous
  // inspiral instead and stop being comparable with sections 94-108.
  gRadiationReactionModel=ChargeRadiationReactionModel::stochasticElectricDipole;
  // The collapse-transit run of audit section 108 would only double the cost.
  gMeasureCollapseTransit=false;
  const char* mode=argv[1];
  const unsigned long long base=strtoull(argv[2],nullptr,10);
  const int count=atoi(argv[3]);
  const double budget=argc>4?atof(argv[4]):900.0;
  gSpinQuantization=std::strcmp(mode,"free")!=0;
  for(int i=0;i<count;++i){
    const unsigned long long seed=base+i;
    if(!std::strcmp(mode,"tilt")){
      // First moment's projection on the orbital normal z (section 100b).
      double c; Vec3 m; prepared(seed,2,c,m);
      std::printf("TILT %llu cosz %+.6f\n",seed,m.z/m.norm());
      continue;
    }
    if(!std::strcmp(mode,"prep")){
      double cp,co; Vec3 mp,mo; prepared(seed,1,cp,mp); prepared(seed,2,co,mo);
      std::printf("PREP seed %llu cos_para %+.9f cos_orto %+.9f |mu1_para-mu1_orto| %.3e\n",
        seed,cp,co,(mp-mo).norm()/mp.norm());
      continue;
    }
    if(!std::strcmp(mode,"quant")){
      const auto p=runCremCollapseExperiment(seed,1,1,budget);
      const auto o=runCremCollapseExperiment(seed,2,1,budget);
      std::printf("PAIR seed %llu para_lab_ps %.9f orto_lab_ps %.9f para_out %d orto_out %d para_ph %zu orto_ph %zu\n",
        seed,p[0].lifetimeSecondsLab*1e12,o[0].lifetimeSecondsLab*1e12,
        (int)p[0].calibrationOutcome,(int)o[0].calibrationOutcome,
        p[0].labFramePhotons.size(),o[0].labFramePhotons.size());
    } else {
      double c; Vec3 m; prepared(seed,1,c,m);
      const auto r=runCremCollapseExperiment(seed,1,1,budget);
      std::printf("FREE seed %llu cos0 %+.6f lab_ps %.9f out %d ph %zu\n",
        seed,c,r[0].lifetimeSecondsLab*1e12,(int)r[0].calibrationOutcome,
        r[0].labFramePhotons.size());
    }
    std::fflush(stdout);
  }
  std::printf("DONE\n");
}
