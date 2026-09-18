// Where does the collapse stop, in units of the orbital action (audit 131)?
//
// For the Kepler problem the Bohr-Sommerfeld actions are J_phi = 2 pi L and
// J_r = 2 pi (K sqrt(mu/2|E|) - L), so J_r + J_phi = n h with
// n = (K/hbar) sqrt(mu/(2|E|)).  The model's emission rule is
// E_photon = hbar omega_orb = 2R/n^3, which gives the exact landing law
//   n -> n/sqrt(1+2/n),   Delta n = n [1 - (1+2/n)^(-1/2)],
// tending to one quantum of action only for n >> 1.  This probe prints that
// ladder and the action at each of the model's own stopping conditions, so
// "does the inspiral stop at an integer action" is a table rather than an
// argument.
//
// Usage: orbital_action_ladder [retardation margin, default 150]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/orbital_action_ladder.cpp -o /tmp/ladder $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const double margin=argc>1?atof(argv[1]):150.0;
    const double axis=pairBohrRadius(activePair);
    const double reducedMass=reducedMassOf(activePair);
    const double attraction=pairCoulombStrength;
    // n from the two osculating elements, exactly as the CREM_ACTION_TRACE
    // added to the estimator computes it.
    const auto levelFromSemiMajorAxis=[&](double a) {
        return std::sqrt(reducedMass*attraction*a)/hbar;
    };
    const auto landing=[](double n) { return n/std::sqrt(1.0+2.0/n); };
    std::printf("a_pair = %.6e m = %.3f r*; n = sqrt(mu K a)/hbar, so "
                "n = sqrt(a/a_pair)\n",axis,axis/comptonBarrierRadius);
    std::printf("identity check, n from a against J_r/h + J_phi/h on a "
                "circular orbit:\n");
    for(double fraction:{1.0,0.5,0.1,0.01}) {
        const double a=axis*fraction;
        const double angular=std::sqrt(reducedMass*attraction*a);  // circular
        const double level=levelFromSemiMajorAxis(a);
        std::printf("  a/a_pair %6.3f  n %.9f  J_phi/h %.9f  J_r/h %.9f\n",
                    fraction,level,angular/hbar,level-angular/hbar);
    }
    std::printf("\nthe ladder the emission rule produces, from the prepared "
                "state:\n");
    double level=1.0;
    for(int step=1;step<=6;++step) {
        const double next=landing(level);
        std::printf("  photon %d: n %.6f -> %.6f   Delta n %.6f   "
                    "a/a_pair %.6f -> %.6f\n",step,level,next,level-next,
                    level*level,next*next);
        level=next;
    }
    std::printf("  and for reference Delta n at n = 10, 100, 1000: "
                "%.6f %.6f %.6f (the correspondence limit is 1)\n",
                10.0-landing(10.0),100.0-landing(100.0),
                1000.0-landing(1000.0));
    // The model's three stopping conditions, in action units.
    // 1. The ground-state floor, both of its clamps.
    const double floorLevel=(attraction/hbar)
        *std::sqrt(reducedMass/(2.0*(-groundStateSpecificEnergy()
                                     *reducedMass)));
    const double floorAngular=
        groundStateSpecificAngularMomentum()*reducedMass/hbar;
    // 2. The retardation safety margin: period/light-crossing = margin on a
    //    circular orbit gives 2 pi c sqrt(mu a / K) = margin, hence
    //    n = margin * K / (2 pi hbar c) = margin * alpha / (2 pi).
    const double fineStructure=attraction/(hbar*c);
    const double retardationLevel=margin*fineStructure/(2.0*pi);
    // 3. The Compton barrier, and the transit endpoint of audit 108.
    const double barrierLevel=levelFromSemiMajorAxis(comptonBarrierRadius);
    const double transitLevel=levelFromSemiMajorAxis(0.005*axis);
    std::printf("\nstopping conditions in action units:\n");
    std::printf("  --ground-state-floor   energy clamp n = %.9f, angular "
                "clamp J_phi/h = %.9f\n",floorLevel,floorAngular);
    std::printf("  retardation margin %-5.0f n = %.6f  (= margin alpha / 2 pi,"
                " alpha = %.8f)\n",margin,retardationLevel,fineStructure);
    std::printf("  Compton barrier        n = %.6f  (a = r* = %.4e m)\n",
                barrierLevel,comptonBarrierRadius);
    std::printf("  transit end 0.005 a    n = %.6f  (audit 108's endpoint)\n",
                transitLevel);
    // Which rung the retardation margin actually catches, starting from the
    // prepared state: the first landing below it.
    level=1.0;
    int steps=0;
    while(level>retardationLevel&&steps<40) { level=landing(level); ++steps; }
    std::printf("\nstarting from n = 1, the first rung below the retardation "
                "margin is n = %.6f after %d photons\n",level,steps);
    std::printf("so the stop is set by the margin, not by the action: it is "
                "the ladder's first rung below %.6f, and nothing about it is "
                "an integer\n",retardationLevel);
}
