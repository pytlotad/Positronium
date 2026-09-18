// Does the M1 reaction torque remove the energy the M1 flux carries away
// (audit section 136)?
//
// On a real trajectory M1/E1 is of order 1e-15, so the balance cannot be read
// off the energy budget in double precision.  It can be read off a controlled
// configuration: a moment precessing uniformly about an axis at rate Omega
// with cone angle theta.  There
//   mu''  = -mu Omega^2 sin(theta) (cos, sin, 0),   P = k |mu''|^2,
//   mu''' =  mu Omega^3 sin(theta) (sin, -cos, 0),
//   N = k (mu x mu''') = k mu^2 Omega^3 sin(theta)
//            (cos(theta) cos, cos(theta) sin, -sin(theta)),
// so the work the torque does on the precession, N . Omega_vector, is
// -k mu^2 Omega^4 sin^2(theta) = -P exactly, at every instant and not only on
// average.  The transverse part is the cone damping.  This probe builds that
// history, calls the model's own dipoleRadiationReaction on it, and checks
// both statements against the model rather than against the algebra.
//
// Usage: m1_channel_closure [cone angle deg, default 30] [samples, default 9]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/m1_channel_closure.cpp -o /tmp/closure $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const double coneDegrees=argc>1?atof(argv[1]):30.0;
    const int samples=argc>2?atoi(argv[2]):9;
    const double cone=coneDegrees*pi/180.0;
    const double moment=firstMagneticMoment;
    // The production precession rate at a_pair, audit 134's |w1|.
    const double rate=4.13417e11;
    const double period=2.0*pi/rate;
    const double coefficient=mu0/(6.0*pi*c*c*c);
    const Vec3 axis{0.0,0.0,1.0};
    // A history of the FIRST moment precessing uniformly; the second is left
    // at zero so the coherent sum is the first moment alone and the balance
    // is read without the interference term of the validation suite.
    const auto historyAt=[&](double centre,double step) {
        StateHistory history;
        for(int index=-(samples-1);index<=0;++index) {
            State sample;
            sample.time=centre+index*step;
            const double phase=rate*sample.time;
            sample.firstDipole={moment*std::sin(cone)*std::cos(phase),
                                moment*std::sin(cone)*std::sin(phase),
                                moment*std::cos(cone)};
            sample.secondDipole={0.0,0.0,0.0};
            sample.firstProperDipole=sample.firstDipole;
            sample.secondProperDipole=sample.secondDipole;
            history.push_back(sample);
        }
        return history;
    };
    std::printf("cone %.2f deg, precession %.6e rad/s (period %.6e s), "
                "moment %.6e J/T\n",coneDegrees,rate,period,moment);
    std::printf("%12s %16s %16s %14s %16s %12s\n","step/period","P flux [W]",
                "-N.Omega [W]","ratio","P analytic [W]","flux/analytic");
    for(double fraction:{1.0e-2,3.0e-3,1.0e-3,3.0e-4,1.0e-4}) {
        const double step=fraction*period;
        const StateHistory history=historyAt(0.0,step);
        const State& state=history.back();
        const DipoleRadiationReaction reaction=
            dipoleRadiationReaction(state,history);
        const double workRate=dot(reaction.firstTorque,axis)*rate;
        const double analytic=coefficient*moment*moment*std::pow(rate,4)
            *std::sin(cone)*std::sin(cone);
        std::printf("%12.1e %16.9e %16.9e %14.9f %16.9e %12.9f\n",
                    fraction,reaction.power,-workRate,
                    -workRate/std::max(reaction.power,1.0e-300),
                    analytic,reaction.power/analytic);
    }
    // The cone damping the transverse torque implies, and the time it takes.
    {
        const double step=1.0e-3*period;
        const StateHistory history=historyAt(0.0,step);
        const State& state=history.back();
        const DipoleRadiationReaction reaction=
            dipoleRadiationReaction(state,history);
        const double gyromagnetic=firstGyromagneticRatioOf();
        // d(mu)/dt from the reaction alone, as applyDipoleRadiationTorque
        // applies it, projected onto the cone angle.
        const Vec3 momentNow=state.firstDipole;
        const Vec3 rateVector=reaction.firstTorque*gyromagnetic;
        const double cosine=dot(momentNow,axis)/momentNow.norm();
        const double coneRate=-(dot(rateVector,axis)
            -cosine*dot(rateVector,momentNow)/momentNow.norm())
            /(momentNow.norm()*std::sin(cone));
        std::printf("\ncone damping: d(theta)/dt = %.6e rad/s, time constant "
                    "theta/|d(theta)/dt| = %.6e s\n",
                    coneRate,std::abs(cone/coneRate));
        std::printf("against the collapse of 30.672 ps that is %.3e collapse "
                    "times\n",std::abs(cone/coneRate)/30.672e-12);
        std::printf("energy in the cone, mu B (1-cos) with B from the same "
                    "precession: rate/gyromagnetic = %.6e T\n",
                    rate/std::abs(gyromagnetic));
    }
}
