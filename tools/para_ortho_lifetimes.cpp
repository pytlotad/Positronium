// Paired para/ortho collapse-time comparison.
//
// runCremCollapseExperiment seeds trajectory i from splitMix64(masterSeed+i)
// regardless of the channel, and the quantized-spin branch of the sampler
// derives the second moment from the first WITHOUT consuming an extra draw.
// Running phenomenon 1 and phenomenon 2 at the same master seed therefore
// gives trajectory pairs with identical initial energy, angular momentum and
// first moment, differing ONLY in the sign of the second moment -- para has
// the moments aligned, ortho anti-aligned.  The pairing matters because the
// initial conditions dominate the spread in collapse time and cancel in the
// per-pair difference.
//
// Reports the paired difference, the geometric ratio with a confidence
// interval, an exact sign test, and how often the two channels disagree on
// the PHOTON COUNT -- the last being the observable that separated them most
// strongly in earlier measurements.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . tools/para_ortho_lifetimes.cpp -o /tmp/probe
#include "modules/crem_collapse.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
const char* stopCauseName(CollapseStopCause cause) {
    switch(cause) {
        case CollapseStopCause::ComptonBarrier:   return "barrier";
        case CollapseStopCause::RetardationLimit: return "retard";
        case CollapseStopCause::GroundStateFloor: return "floor";
        default:                                  return "none";
    }
}
// Two-sided exact binomial probability at p=1/2, summing every outcome no
// more likely than the observed one.  With p=1/2 the distribution is
// symmetric, so that is twice the lower tail.
double signTestProbability(int successes,int trials) {
    if(trials<=0) return 1.0;
    const int extreme=std::min(successes,trials-successes);
    long double tail=0.0L, term=1.0L;   // term = C(trials,k)/2^trials at k=0
    for(int k=0;k<trials;++k) term/=2.0L;
    for(int k=0;k<=extreme;++k) {
        tail+=term;
        term=term*static_cast<long double>(trials-k)
                 /static_cast<long double>(k+1);
    }
    const double probability=static_cast<double>(2.0L*tail);
    return probability>1.0?1.0:probability;
}
}

int main(int argc,char** argv) {
    const int runCount=argc>1?std::atoi(argv[1]):20;
    const std::uint64_t masterSeed=
        argc>2?std::strtoull(argv[2],nullptr,10):7ULL;
    const double budget=argc>3?std::atof(argv[3]):90.0;
    // Emission model, because it decides whether this comparison has any
    // resolving power at all.  Deterministic emission (the default) fires at
    // a fixed accumulated-hazard threshold, so two paired trajectories carry
    // no shot noise and a sub-percent channel difference is measurable.
    // Poisson emission draws the threshold from Exp(1) and its spread --
    // sigma/mean near 0.9 -- swamps any such difference.
    if(argc>4&&std::string(argv[4])=="poisson") gDeterministicEmission=false;
    // Starting separation and the ground-state floor, because the answer
    // depends on both: level 2 with the floor measures a cascade to n=1,
    // level 1 without it measures an inspiral to the Compton barrier, and
    // they are different physical questions with different answers.
    if(argc>5) gInitialPrincipalLevel=std::atoi(argv[5]);
    if(argc>6&&std::string(argv[6])=="floor") gGroundStateEmissionFloor=true;
    // Branch diagnosis.  At level 1 the ortho channel takes one of two
    // discrete collapse times rather than scattering about one, so the
    // question "what picks the branch" needs the terminal state of each
    // trajectory, not just its duration.
    const bool diagnoseBranch=argc>7&&std::string(argv[7])=="branch";
    std::cout<<"paired para/ortho: "<<runCount<<" trajectories per channel, "
             <<"master seed "<<masterSeed<<", budget "<<budget<<" s, emission "
             <<(gDeterministicEmission?"deterministic":"poisson")
             <<", level "<<gInitialPrincipalLevel<<", floor "
             <<(gGroundStateEmissionFloor?"on":"off")<<"\n";
    const auto para=runCremCollapseExperiment(masterSeed,1,runCount,budget);
    const auto ortho=runCremCollapseExperiment(masterSeed,2,runCount,budget);

    std::vector<double> differencePs, logRatio;
    int orthoLonger=0, photonCountDiffers=0, bothComplete=0;
    double paraSum=0.0, orthoSum=0.0;
    std::cout<<std::setprecision(9);
    std::cout<<"\n idx    t_para(ps)     t_ortho(ps)    ratio    "
             <<"photons para/ortho\n";
    for(int index=0;index<runCount;++index) {
        const auto& p=para[static_cast<size_t>(index)];
        const auto& o=ortho[static_cast<size_t>(index)];
        const bool complete=
            p.calibrationOutcome==SimulationOutcome::ReachedCutoff
            &&o.calibrationOutcome==SimulationOutcome::ReachedCutoff;
        // The photon count is collected from every trajectory, complete or
        // not: a censored one still fired whatever it fired.
        if(p.emittedPhotonCount!=o.emittedPhotonCount) ++photonCountDiffers;
        if(!complete) continue;
        const double paraPs=p.lifetimeSecondsLab*1.0e12;
        const double orthoPs=o.lifetimeSecondsLab*1.0e12;
        if(!(paraPs>0.0)||!(orthoPs>0.0)) continue;
        ++bothComplete;
        paraSum+=paraPs; orthoSum+=orthoPs;
        differencePs.push_back(orthoPs-paraPs);
        logRatio.push_back(std::log(orthoPs/paraPs));
        if(orthoPs>paraPs) ++orthoLonger;
        std::cout<<std::setw(4)<<index<<std::setw(15)<<paraPs
                 <<std::setw(16)<<orthoPs
                 <<std::setw(11)<<orthoPs/paraPs
                 <<std::setw(9)<<p.emittedPhotonCount<<" /"
                 <<std::setw(6)<<o.emittedPhotonCount<<"\n";
    }
    if(diagnoseBranch) {
        std::cout<<"\n branch diagnosis (both channels, every trajectory)\n"
                 <<" idx ch  t(ps)         stop      peri/barrier"
                 <<"  period/lightcross  a_term(m)      bind(eV)"
                 <<"   revolutions   photons  E_emitted(J)\n";
        for(int index=0;index<runCount;++index) {
            for(int channel=0;channel<2;++channel) {
                const auto& e=(channel==0?para:ortho)[
                    static_cast<size_t>(index)];
                std::cout<<std::setw(4)<<index
                         <<(channel==0?"  p":"  o")
                         <<std::setw(14)<<e.calibrationSecondsLab*1.0e12
                         <<std::setw(10)<<stopCauseName(e.stopCause)
                         <<std::setw(14)<<e.terminalPeriapsisOverBarrier
                         <<std::setw(19)<<e.terminalPeriodToLightCrossing
                         <<std::setw(15)<<e.terminalSemiMajorAxis
                         <<std::setw(12)<<e.terminalBindingEnergy/eCharge
                         <<std::setw(14)<<e.revolutions
                         <<std::setw(9)<<e.emittedPhotonCount
                         <<std::setw(15)<<e.quantizedEmittedEnergyJoules
                         <<"\n";
            }
        }
    }
    if(bothComplete<2) {
        std::cout<<"\nonly "<<bothComplete
                 <<" complete pairs; raise the budget or the run count\n";
        return 1;
    }
    const auto summarize=[](const std::vector<double>& sample,
                            double& mean,double& standardError) {
        mean=0.0;
        for(double value:sample) mean+=value;
        mean/=static_cast<double>(sample.size());
        double variance=0.0;
        for(double value:sample) variance+=(value-mean)*(value-mean);
        variance/=static_cast<double>(sample.size()-1);
        standardError=std::sqrt(variance/static_cast<double>(sample.size()));
    };
    double meanDifference=0.0, differenceError=0.0;
    double meanLog=0.0, logError=0.0;
    summarize(differencePs,meanDifference,differenceError);
    summarize(logRatio,meanLog,logError);
    std::cout<<"\ncomplete pairs:        "<<bothComplete<<" / "<<runCount<<"\n"
             <<"mean t_para:           "<<paraSum/bothComplete<<" ps\n"
             <<"mean t_ortho:          "<<orthoSum/bothComplete<<" ps\n"
             <<"DIFFERENCE ortho-para: "<<meanDifference<<" +/- "
             <<differenceError<<" ps   ("
             <<(differenceError>0?std::abs(meanDifference)/differenceError:0.0)
             <<" sigma)\n"
             <<"GEOMETRIC RATIO:       "<<std::exp(meanLog)
             <<",  95% CI ["<<std::exp(meanLog-1.96*logError)
             <<", "<<std::exp(meanLog+1.96*logError)<<"]\n"
             <<"SIGN TEST:             ortho longer in "<<orthoLonger
             <<" of "<<bothComplete<<",  exact p = "
             <<signTestProbability(orthoLonger,bothComplete)<<"\n"
             <<"PHOTON COUNT DIFFERS:  "<<photonCountDiffers<<" of "
             <<runCount<<"  ("
             <<100.0*photonCountDiffers/static_cast<double>(runCount)<<"%)\n";
    return 0;
}
