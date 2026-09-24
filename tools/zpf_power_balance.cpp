// PLAN STEP 4: the stochastic-electrodynamics power balance.
//
// The plan (audits/2026-09-11-plan-uwolnienia-modelu.md) calls this the
// single most important measurement it contains, because it is the only
// remaining candidate for a scale the model does not put in by hand.  The
// question is where the power the zero-point field pumps IN equals the power
// the orbit radiates OUT; that radius would be chosen by the dynamics rather
// than by the L = hbar preparation.
//
// WHAT THE PLAN ALREADY WARNS.  Its own note says the FIELD-equality estimate
// -- ZPF E_rms going as omega^2 ~ r^-3 against a binding field going as r^-2,
// crossing at 1.05e-13 m -- is NOT this measurement: the balance is about
// POWER and carries different exponents.  This tool measures the power.
//
// HOW.  P_abs is the work the zero-point field does on the two charges,
// <q1 E.v1 + q2 E.v2>, accumulated along a trajectory the field is actually
// acting on (gZeroPointField enters the Lorentz force, so the orbit responds
// to it and the correlation that carries the absorption is present).  P_rad
// is the coherent electric-dipole power |d''|^2/(6 pi eps0 c^3) of the same
// state, which is the channel the model's inspiral runs on.
//
// WHAT IT IS NOT, and this is the whole difficulty.  <q E.v> has expectation
// zero for a field whose phases are uncorrelated with the motion, so a short
// average measures its VARIANCE, not its mean.  Audit 169 measured a factor
// of 21 between seeds at one radius: the per-seed number is not a
// measurement of the balance, and only the sign was systematic (six of six
// positive).  Averaging that down is the work this tool does not do.
//
// IS THE FIELD A PERTURBATION AT ALL.  The band rides the orbital frequency,
// so its modes are RESONANT with the orbit by construction, and over tens of
// orbits a resonant mode can pump the pair apart instead of holding it.  When
// that happens the radiated power is computed from an orbit that no longer
// exists and the ratio means nothing: audit 170 measured <P_rad> jumping
// from 0.78 W to 5.5e+07 W between scale 0.3 and 1 at 0.01 a_pair, purely
// from disrupted trajectories.  So every row reports the ratio of the final
// separation to the initial one, and a row whose orbit did not survive is
// not a measurement of a balance.
//
// SAMPLING.  The step is tied to the FASTEST zero-point mode, not to the
// orbit, and the steps-per-cycle column reports it.  The step census in
// positronium.cpp warns why: a band whose upper edge sits above the orbital
// frequency can be aliased by a step that still satisfies the trajectory
// tolerance, and aliasing would fake exactly the net work being looked for.
//
// A FIXED BAND.  CREM_ZPF_FIXED_FREQUENCY (rad/s) freezes the band instead of
// letting it ride the orbital frequency; this tool honours it everywhere it
// matters -- the phase accumulation, the sampling step and the reported
// E_rms all follow the frozen frequency, so the measurement is of a band
// that genuinely does not ride.  Audit 171 is that experiment.
//
// Usage: zpf_power_balance [scale] [modes] [seed] [steps per ZPF cycle]
//                          [orbits] [r/a_pair, or 0 for the default scan]
// With a single radius it prints one machine-readable row, which is how the
// seed ensemble of audit 170 is built: run the seeds in parallel and
// aggregate, rather than averaging inside one process.
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/zpf_power_balance.cpp -o /tmp/zpf $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
#include <vector>

int main(int argc,char** argv){
    const double scale=(argc>1?atof(argv[1]):1.0);
    const int modes=(argc>2?atoi(argv[2]):16);
    const std::uint64_t seed=(argc>3?strtoull(argv[3],nullptr,10):42);
    const double perCycle=(argc>4?atof(argv[4]):32.0);
    const double orbits=(argc>5?atof(argv[5]):1.0);
    const double singleRadius=(argc>6?atof(argv[6]):0.0);
    constexpr double lo=0.3, hi=3.0;
    const double total=firstMass+secondMass;
    gZeroPointField=makeZeroPointField(lo,hi,modes,scale,seed);
    std::printf("ZPF scale %.3g, %d modes over [%.1f, %.1f] x omega_orb, "
                "seed %llu, %.0f steps per fastest cycle\n",
                scale,modes,lo,hi,(unsigned long long)seed,perCycle);
    if(!(singleRadius>0.0))
        std::printf("%9s %13s %13s %13s %12s %9s\n","r/a_pair","P_abs [W]",
                    "P_rad [W]","P_abs/P_rad","E_rms [V/m]","steps/cyc");
    const double aPair=pairBohrRadius(activePair);
    const std::vector<double> scan=singleRadius>0.0
        ?std::vector<double>{singleRadius}
        :std::vector<double>{1.0,0.3,0.1,0.03,0.01};
    for(double fraction:scan){
        const double r=fraction*aPair;
        const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r));
        const double omega=speed/r;
        const double period=2.0*pi/omega;
        State s{};
        s.firstPosition={r*secondMass/total,0,0};
        s.secondPosition={-r*firstMass/total,0,0};
        s.firstVelocity={0,speed*secondMass/total,0};
        s.secondVelocity={0,-speed*firstMass/total,0};
        s.firstProperDipole={0,0,firstMagneticMoment};
        s.secondProperDipole={0,0,secondMagneticMoment};
        synchronizeCovariantDipoles(s);
        ClassicalTrajectoryEngine::Accuracy accuracy;
        accuracy.relativeTolerance=1.0e-7;
        accuracy.maximumDepth=18;
        accuracy.computeOutwardFlux=false;
        accuracy.reactionModel=ChargeRadiationReactionModel::disabled;
        accuracy.useRetardedExternalForces=false;
        ClassicalTrajectoryEngine engine(s,accuracy);
        // Every rate below follows the band's OWN frequency, which is the
        // orbital one unless CREM_ZPF_FIXED_FREQUENCY froze it.
        const double drive=zeroPointDriveFrequency(s).frequency;
        // Resolve the FASTER of the two clocks.  With the riding band the
        // fastest mode always beats the orbit, so this is the mode; with a
        // frozen band far below the orbit it is the orbit instead, and
        // basing the step on the band alone would leave the work integral
        // sampling a few points per revolution -- which is the aliasing the
        // step census warns about, arriving from the other side.
        const double fastestPeriod=
            std::min(2.0*pi/(hi*drive),period);
        const double step=fastestPeriod/perCycle;
        const int stepCount=static_cast<int>(orbits*period/step);
        double work=0.0,elapsed=0.0,phase=0.0;
        const double initialSeparation=separation(s);
        int taken=0;
        for(int index=0;index<stepCount;++index){
            const Vec3 firstPosition=s.firstPosition;
            const Vec3 secondPosition=s.secondPosition;
            const Vec3 firstVelocity=s.firstVelocity;
            const Vec3 secondVelocity=s.secondVelocity;
            const ZeroPointDrive band=zeroPointDriveFrequency(s);
            const double orbital=band.frequency;
            const double orbitalRate=band.derivative;
            Vec3 firstElectric,firstMagnetic,secondElectric,secondMagnetic;
            gZeroPointField.sample(firstPosition,orbital,orbitalRate,
                s.zeroPointPhase,firstElectric,firstMagnetic);
            gZeroPointField.sample(secondPosition,orbital,orbitalRate,
                s.zeroPointPhase,secondElectric,secondMagnetic);
            const double instantaneous=
                firstCharge*dot(firstElectric,firstVelocity)
               +secondCharge*dot(secondElectric,secondVelocity);
            if(!engine.advance(s,step)) break;
            work+=instantaneous*step;
            elapsed+=step;
            phase+=orbital*step;
            s.zeroPointPhase=phase;
            ++taken;
        }
        if(taken<10||!(elapsed>0.0)){
            std::printf("%9.3g   (engine stopped after %d steps)\n",
                        fraction,taken);
            continue;
        }
        const Vec3 dipoleSecond=s.firstAcceleration*firstCharge
                               +s.secondAcceleration*secondCharge;
        const double radiated=
            dipoleSecond.squaredNorm()/(6.0*pi*epsilon0*c*c*c);
        const double absorbed=work/elapsed;
        const double rootMeanSquare=gZeroPointField.amplitudeCoefficient
            *drive*drive*std::sqrt(modes/2.0);
        if(singleRadius>0.0)
            // ROW <seed> <scale> <r/a_pair> <P_abs> <P_rad> <orbits>
            //     <final/initial separation> <band frequency>
            std::printf("ROW %llu %.6g %.6g %.9e %.9e %.4f %.6g %.6e\n",
                (unsigned long long)seed,scale,fraction,absorbed,radiated,
                elapsed/period,separation(s)/initialSeparation,drive);
        else
            std::printf("%9.3g %13.4e %13.4e %13.4e %12.3e %9.1f\n",
                fraction,absorbed,radiated,
                absorbed/std::max(std::abs(radiated),1.0e-300),
                rootMeanSquare,perCycle);
    }
    gZeroPointField=ZeroPointField{};
}
