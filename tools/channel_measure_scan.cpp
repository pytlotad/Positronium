// Would a different MEASURE over the mutual angle give the 1:3 branching
// ratio (audit section 112)?
//
// The weight w = (1+cos)/2 is fixed by the two forced ends (audit 91), so
// E[w] = (1+<cos>)/2 and 1:3 needs <cos> = -1/2 exactly.  This probe measures
// <cos> and E[w] under the measures the model itself can produce: the
// preparation measure it draws from, the terminal measure its own dynamics
// produces (CREM_CHANNEL_AT_ANNIHILATION), and the same terminal sample
// reweighted by the collapse time, which is what an ensemble observed AT
// decay would be.  It also prints the Boltzmann tilt exp(-U_dd/E) for the
// energy scales the model has, and the scale that would be needed.
//
// Usage: channel_measure_scan <trajectories> [seed base, default 42]
//        [phenomenon, default 1] [budget s, default 600]
//        [quant to impose the mutual angle] [stochastic for the photon cascade]
// The last two reproduce the configuration the terminal-moment numbers of
// README ("Kanal wybiera teraz dynamika") were taken in: audit 91 withdrew
// the imposition and audit 109 switched the default channel to continuous.
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/channel_measure_scan.cpp -o /tmp/measure $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc,char** argv) {
  const int count=argc>1?atoi(argv[1]):20;
  const unsigned long long base=argc>2?strtoull(argv[2],nullptr,10):42ULL;
  const int phenomenon=argc>3?atoi(argv[3]):1;
  const double budget=argc>4?atof(argv[4]):600.0;
  bool imposedAngle=false,quantizedChannel=false;
  for(int index=5;index<argc;++index) {
    if(std::string(argv[index])=="quant") imposedAngle=true;
    if(std::string(argv[index])=="stochastic") quantizedChannel=true;
  }
  gSpinQuantization=imposedAngle;
  if(quantizedChannel)
    gRadiationReactionModel=
        ChargeRadiationReactionModel::stochasticElectricDipole;
  // The weight must come from the TERMINAL moments here: that is the measure
  // the dynamics produces, as opposed to the one the sampler draws.
  setenv("CREM_CHANNEL_AT_ANNIHILATION","1",1);
  gMeasureCollapseTransit=false;
  const std::vector<CremCollapseEstimate> estimates=
      runCremCollapseExperiment(base,phenomenon,count,budget);
  std::vector<double> preparedCos,terminalWeight,lifetimes;
  int drawnTwoPhoton=0,emittedTwoPhoton=0,labelPara=0,labelAgrees=0;
  std::uint64_t drawStream=splitMix64(base^0x5851f42d4c957f2dULL);
  for(int index=0;index<count;++index) {
    const CremCollapseEstimate& estimate=estimates[static_cast<size_t>(index)];
    SimulationOptions options;
    options.frameCount=2;
    options.observationTime=1.0e-24;
    const SimulationResult prepared=simulate(
        splitMix64(base+static_cast<unsigned long long>(index)),
        phenomenon,options);
    if(prepared.frames.empty()) continue;
    const Vec3 first=prepared.frames.front().firstDipole;
    const Vec3 second=prepared.frames.front().secondDipole;
    const double cosine=dot(first,second)/(first.norm()*second.norm());
    if(!std::isfinite(estimate.annihilationTwoPhotonWeight)
       ||!std::isfinite(estimate.lifetimeSecondsLab)) {
      std::printf("ROW %d cos %+.6f  (no terminal weight: outcome %d)\n",
                  index,cosine,static_cast<int>(estimate.calibrationOutcome));
      continue;
    }
    preparedCos.push_back(cosine);
    terminalWeight.push_back(estimate.annihilationTwoPhotonWeight);
    lifetimes.push_back(estimate.lifetimeSecondsLab);
    // What the run itself emitted, and what the terminal weight would draw.
    const std::size_t emitted=estimate.annihilationPhotonEnergies.size();
    if(emitted==2) ++emittedTwoPhoton;
    const bool draws=drawUniformUnit(drawStream)
        <estimate.annihilationTwoPhotonWeight;
    if(draws) ++drawnTwoPhoton;
    const bool para=cosine>=0.5;
    if(para) ++labelPara;
    if(para==(emitted==2)) ++labelAgrees;
    std::printf("ROW %d cos_prepared %+.6f w_prepared %.6f w_terminal %.6f "
                "cos_terminal %+.6f photons %zu t_ps %.6f\n",
                index,cosine,0.5*(1.0+cosine),
                estimate.annihilationTwoPhotonWeight,
                2.0*estimate.annihilationTwoPhotonWeight-1.0,
                emitted,estimate.lifetimeSecondsLab*1e12);
  }
  const int total=static_cast<int>(terminalWeight.size());
  if(total==0) { std::printf("no completed trajectories\n"); return 1; }
  const auto average=[](const std::vector<double>& values,
                        const std::vector<double>& weights) {
    double sum=0.0,norm=0.0;
    for(size_t index=0;index<values.size();++index) {
      const double weight=weights.empty()?1.0:weights[index];
      sum+=weight*values[index]; norm+=weight;
    }
    return norm>0.0?sum/norm:NAN;
  };
  const std::vector<double> none;
  std::vector<double> inverseLifetime;
  for(double value:lifetimes)
    inverseLifetime.push_back(value>0.0?1.0/value:0.0);
  std::vector<double> preparedWeight;
  for(double value:preparedCos) preparedWeight.push_back(0.5*(1.0+value));
  std::printf("\ncompleted %d trajectories (seed base %llu, phenomenon %d, "
              "imposed angle %d, quantized channel %d)\n",
              total,base,phenomenon,imposedAngle?1:0,quantizedChannel?1:0);
  {
    std::vector<double> sortedWeight=terminalWeight;
    std::sort(sortedWeight.begin(),sortedWeight.end());
    std::printf("terminal weight: min %.4f max %.4f mean %.4f\n",
                sortedWeight.front(),sortedWeight.back(),
                average(terminalWeight,none));
    std::printf("emitted 2 gamma %d/%d = %.1f%%; a draw on the terminal "
                "weight would give %d/%d = %.1f%%\n",
                emittedTwoPhoton,total,100.0*emittedTwoPhoton/total,
                drawnTwoPhoton,total,100.0*drawnTwoPhoton/total);
    std::printf("label cos>=0.5 %d/%d; label agrees with the multiplicity "
                "%d/%d = %.1f%%\n",labelPara,total,labelAgrees,total,
                100.0*labelAgrees/total);
  }
  std::printf("M1 preparation measure : <cos> %+.4f  E[w] %.4f\n",
              average(preparedCos,none),average(preparedWeight,none));
  std::printf("M2 terminal measure    : <cos> %+.4f  E[w] %.4f\n",
              2.0*average(terminalWeight,none)-1.0,
              average(terminalWeight,none));
  std::printf("M3 weighted by t       : E[w] %.4f;  by 1/t: E[w] %.4f\n",
              average(terminalWeight,lifetimes),
              average(terminalWeight,inverseLifetime));
  // M4: Boltzmann tilt exp(-U_dd(cos)/E) with U_dd = -(mu0/4pi) mu^2 cos/r^3
  // at r = a_pair, azimuth-averaged so only the cos term survives.
  const double momentScale=firstMagneticMoment*secondMagneticMoment;
  const double separation=pairBohrRadius(activePair);
  const double coupling=1.0e-7*momentScale
      /(separation*separation*separation);   // mu0/(4 pi) mu1 mu2 / r^3
  const auto tiltedMean=[&](double energyScale) {
    // <w> under p(cos) ~ exp(+coupling*cos/energyScale) (aligned moments are
    // the LOW-energy end for the tensor term at this geometry), cos uniform.
    const double lambda=energyScale>0.0?coupling/energyScale:0.0;
    const int samples=200001; double sum=0.0,norm=0.0;
    for(int index=0;index<samples;++index) {
      const double cosine=-1.0+2.0*index/(samples-1.0);
      const double density=std::exp(lambda*cosine);
      sum+=density*0.5*(1.0+cosine); norm+=density;
    }
    return sum/norm;
  };
  const double binding=pairBindingEnergy(activePair);
  const double orbitalQuantum=2.0*binding;   // hbar omega_orb at a_pair
  std::printf("M4 Boltzmann tilt, coupling U_dd(a_pair) = %.3e J = %.3e eV\n",
              coupling,coupling/1.602176634e-19);
  std::printf("   E = U_dd itself      : E[w] %.4f\n",tiltedMean(coupling));
  std::printf("   E = binding %.3e eV : E[w] %.4f\n",
              binding/1.602176634e-19,tiltedMean(binding));
  std::printf("   E = hbar omega %.3e eV: E[w] %.4f\n",
              orbitalQuantum/1.602176634e-19,tiltedMean(orbitalQuantum));
  // The scale that WOULD give <cos> = -1/2, i.e. E[w] = 1/4, by bisection on
  // lambda with the sign that favours anti-alignment.
  double low=0.0,high=1.0e4;
  for(int iteration=0;iteration<200;++iteration) {
    const double middle=0.5*(low+high);
    const int samples=20001; double sum=0.0,norm=0.0;
    for(int index=0;index<samples;++index) {
      const double cosine=-1.0+2.0*index/(samples-1.0);
      const double density=std::exp(-middle*cosine);
      sum+=density*0.5*(1.0+cosine); norm+=density;
    }
    if(sum/norm>0.25) low=middle; else high=middle;
  }
  const double lambdaNeeded=0.5*(low+high);
  std::printf("   needed: p(cos) ~ exp(-lambda cos) with lambda %.4f, i.e. an "
              "energy scale %.3e eV\n     = %.2e of U_dd(a_pair) and %.2e of "
              "the binding\n",lambdaNeeded,
              coupling/lambdaNeeded/1.602176634e-19,
              1.0/lambdaNeeded,coupling/lambdaNeeded/binding);
  std::printf("rule: a measure gives 1:3 when E[w] is in [0.24, 0.26]\n");
}
