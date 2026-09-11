// Checkpoint-by-checkpoint trace of ONE collapse trajectory, to answer what
// selects between the two discrete collapse times the ortho channel takes at
// --level 1.
//
// runCremCollapseExperiment seeds trajectory i from splitMix64(masterSeed+i),
// so naming the master seed and the index here reproduces exactly the
// trajectory that experiment ran, and running a single one keeps the census
// output free of the interleaving four workers would produce.
//
// The census itself is CREM_SKIP_CENSUS, which already exists for this
// question: it prints the pre-truncation float beside the integer number of
// orbits the checkpoint skips, which is the only discrete quantity left in
// the secular path, together with the eccentricity, the loss per orbit, the
// measured period, the net spin along the orbital axis and the orbital
// angular momentum.  Two trajectories that land in different branches must
// differ in one of those before they differ in duration.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/para_ortho_branch_trace.cpp -o /tmp/trace $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    if(argc<4) {
        std::cerr<<"usage: "<<argv[0]<<" <master seed> <index> <para|ortho>"
                 <<" [budget seconds]\n";
        return 2;
    }
    const std::uint64_t masterSeed=std::strtoull(argv[1],nullptr,10);
    const int index=std::atoi(argv[2]);
    const std::string channel=argv[3];
    const double budget=argc>4?std::atof(argv[4]):90.0;
    // Level and floor, because the two configurations behave differently:
    // with the floor the cascade stops at n=1, far above the depth where one
    // photon outweighs the remaining binding; without it the inspiral runs
    // into exactly that region.
    if(argc>5) gInitialPrincipalLevel=std::atoi(argv[5]);
    if(argc>6&&std::string(argv[6])=="floor") gGroundStateEmissionFloor=true;
    const int phenomenon=channel=="para"?1:2;
    // Set before the call: the census is read through getenv at the
    // checkpoint site, so this reaches it without a wrapper script.
    setenv("CREM_SKIP_CENSUS","1",1);
    const std::uint64_t trajectorySeed=
        splitMix64(masterSeed+static_cast<std::uint64_t>(index));
    std::cerr<<"# trace master="<<masterSeed<<" index="<<index
             <<" channel="<<channel<<" seed="<<trajectorySeed
             <<" level="<<gInitialPrincipalLevel<<" floor="
             <<(gGroundStateEmissionFloor?"on":"off")<<"\n";
    const CremCollapseEstimate estimate=
        estimateCremCollapse(trajectorySeed,phenomenon,budget);
    std::cout<<std::setprecision(12)
             <<"channel "<<channel<<"  index "<<index
             <<"  t="<<estimate.calibrationSecondsLab*1.0e12<<" ps"
             <<"  outcome="
             <<static_cast<int>(estimate.calibrationOutcome)
             <<"  stop="<<static_cast<int>(estimate.stopCause)
             <<"  peri/barrier="<<estimate.terminalPeriapsisOverBarrier
             <<"  P/LC="<<estimate.terminalPeriodToLightCrossing
             <<"  revolutions="<<estimate.revolutions
             <<"  photons="<<estimate.emittedPhotonCount<<"\n";
    // Terminal orbital elements.  L is recorded by the estimator itself;
    // L_contact is what it would have to reach for a Kepler orbit to touch
    // the Compton barrier, sqrt(mu k r_c)/hbar on a circular orbit.
    const double terminalA=estimate.terminalSemiMajorAxis;
    const double contactAngularMomentum=
        std::sqrt(pairCoulombStrength*reducedMassOf(activePair)
                  *comptonBarrierRadius)/hbar;
    std::cout<<"  terminal: a="<<terminalA<<" m"
             <<"  a/barrier="<<terminalA/comptonBarrierRadius
             <<"  binding="<<estimate.terminalBindingEnergy/eCharge<<" eV"
             <<"  L="<<estimate.terminalAngularMomentum<<" hbar"
             <<"  L/L_contact="
             <<estimate.terminalAngularMomentum/contactAngularMomentum
             <<"  h^2/(Aa)="<<estimate.terminalKeplerConsistency<<"\n";
    // What the dynamics chose, and the |J| it chose it from.
    std::cout<<"  channel:  |J|="
             <<estimate.annihilationTotalAngularMomentum<<" hbar"
             <<"  ->  J="<<estimate.annihilationTotalAngularMomentumQuantum
             <<"  ->  "<<estimate.annihilationPhotonEnergies.size()
             <<" photons  (prepared as "<<channel<<")\n";
    return 0;
}
