// First-passage time of the orbital angular momentum through contact.
//
// A classical Kepler orbit reaches separation r_c only if
// L <= sqrt(mu k r_c (1+e)); on the circular orbits this model produces that
// is sqrt(mu k r_c), which for e+e- is 0.0427 hbar at the Compton barrier and
// 0.00516 hbar at the classical electron radius.  Each photon removes
// angular momentum, so L walks down and "when can the pair first touch" is a
// first-passage time -- a rate scale the model produces itself, with no
// annihilation width imported from QED.
//
// WHAT THIS IS NOT.  Reaching contact is a necessary condition, not a rate.
// The probability of annihilating once there is a cross-section this model
// does not carry, so the number below is "how long until the classical orbit
// stops being the obstacle", not a lifetime.  It is also measured with
// CREM_RETARDATION_LIMIT lowered, because the crossing happens one photon
// after the default stop; below 150 the pair ends up inside the Compton
// barrier, where classical point-particle electrodynamics does not apply.
// Every number here is a probe of the model's own bookkeeping.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/contact_first_passage.cpp -o /tmp/probe $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {
void summarize(const char* label,std::vector<double> sample,int attempted) {
    std::cout<<"  "<<std::left<<std::setw(26)<<label<<std::right;
    if(sample.empty()) {
        std::cout<<"reached by 0 of "<<attempted<<"\n";
        return;
    }
    std::sort(sample.begin(),sample.end());
    const std::size_t count=sample.size();
    const double median=count%2
        ?sample[count/2]
        :0.5*(sample[count/2-1]+sample[count/2]);
    double mean=0.0;
    for(double value:sample) mean+=value;
    mean/=static_cast<double>(count);
    double variance=0.0;
    for(double value:sample) variance+=(value-mean)*(value-mean);
    variance=count>1?variance/static_cast<double>(count-1):0.0;
    const double deviation=std::sqrt(variance);
    std::cout<<"reached by "<<count<<" of "<<attempted
             <<"   median "<<median<<" ps   mean "<<mean
             <<" +/- "<<deviation/std::sqrt(static_cast<double>(count))
             <<" ps   sigma/mean "<<(mean!=0.0?deviation/mean:0.0)<<"\n";
}
}

int main(int argc,char** argv) {
    const int runCount=argc>1?std::atoi(argv[1]):24;
    const std::uint64_t masterSeed=
        argc>2?std::strtoull(argv[2],nullptr,10):42ULL;
    const double budget=argc>3?std::atof(argv[3]):200.0;
    const std::string retardationLimit=argc>4?argv[4]:"5";
    // Emission model, because it decides whether this measurement has a
    // SHAPE at all.  Deterministic emission fires at a fixed accumulated
    // hazard, so the first passage is a near-delta and the model produces a
    // characteristic time rather than a rate.  Poisson emission draws each
    // threshold from Exp(1), so the passage is the SUM of the waiting times
    // for however many photons it takes -- a Gamma distribution of that
    // shape, with sigma/mean = 1/sqrt(count).  Annihilation is exponential,
    // i.e. Gamma of shape one, so this is where the mechanism's distribution
    // can be compared against the law it would have to reproduce.
    if(argc>5&&std::string(argv[5])=="poisson") gDeterministicEmission=false;
    setenv("CREM_RETARDATION_LIMIT",retardationLimit.c_str(),1);
    const double contactAtBarrier=
        std::sqrt(reducedMassOf(activePair)*pairCoulombStrength
                  *comptonBarrierRadius)/hbar;
    const double contactAtElectronRadius=
        std::sqrt(reducedMassOf(activePair)*pairCoulombStrength
                  *classicalElectronRadius)/hbar;
    std::cout<<std::setprecision(6)
             <<"first passage of L through contact, "<<runCount
             <<" trajectories per channel, master seed "<<masterSeed
             <<", budget "<<budget<<" s\n"
             <<"retardation limit lowered to "<<retardationLimit
             <<" (default 150) -- outside the declared domain, see the tool's"
             <<" own comment\n"
             <<"contact L: "<<contactAtBarrier<<" hbar at the Compton barrier,"
             <<" "<<contactAtElectronRadius
             <<" hbar at the classical electron radius\n"
             <<"emission: "
             <<(gDeterministicEmission?"deterministic":"poisson")<<"\n\n";
    for(int phenomenon=1;phenomenon<=2;++phenomenon) {
        const auto estimates=runCremCollapseExperiment(
            masterSeed,phenomenon,runCount,budget);
        std::vector<double> barrier,electronRadius,terminalL;
        int completed=0;
        for(const auto& estimate:estimates) {
            if(estimate.calibrationOutcome==SimulationOutcome::ReachedCutoff)
                ++completed;
            if(std::isfinite(estimate.contactPassageAtBarrierSeconds))
                barrier.push_back(
                    estimate.contactPassageAtBarrierSeconds*1.0e12);
            if(std::isfinite(
                   estimate.contactPassageAtElectronRadiusSeconds))
                electronRadius.push_back(
                    estimate.contactPassageAtElectronRadiusSeconds*1.0e12);
            if(std::isfinite(estimate.terminalAngularMomentum))
                terminalL.push_back(estimate.terminalAngularMomentum);
        }
        std::cout<<(phenomenon==1?"para":"ortho")<<": "<<completed
                 <<" of "<<runCount<<" reached a stopping condition\n";
        summarize("through Compton barrier",barrier,runCount);
        summarize("through electron radius",electronRadius,runCount);
        if(!terminalL.empty()) {
            double mean=0.0;
            for(double value:terminalL) mean+=value;
            mean/=static_cast<double>(terminalL.size());
            std::cout<<"  terminal L mean            "<<mean<<" hbar = "
                     <<mean/contactAtBarrier<<" x contact\n";
        }
        std::cout<<'\n';
    }
    return 0;
}
