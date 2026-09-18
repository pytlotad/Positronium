// Does ONE spin content give both 1:3 and 203.3941 GHz (audit section 130)?
//
// The two observables constrain different things.  The annihilation ratio is a
// statement about the MEASURE over the mutual angle: E[w] = (1+<cos>)/2, so
// 1:3 <=> <cos> = -1/2 and nothing else matters (audit 112b).  The splitting
// is a statement about an absolute ENERGY scale.  In quantum mechanics one
// object supplies both with no free parameter: 1:3 is the degeneracy of the
// four levels and 203.3941 GHz is the gap between the same four.  A classical
// continuum model has to get the ratio from a polarization, which needs a
// scale of its own, so it has one parameter more than the data it explains --
// unless that scale comes from inside the model.
//
// Section 129 supplies it: the mutual angle's own oscillation has an action
// J = 0.0985 hbar, adiabatically invariant to a factor 1.40 along the whole
// collapse, so E_lib = J |omega1 - omega2| is an internal energy for exactly
// the degree of freedom the measure lives on.  The measure tested here is
// therefore p(cos) ~ exp(-U(cos)/E_lib) with U from the model's own
// azimuthAveragedDipoleEnergy and omega from its own orbit-averaged BMT rates.
//
// Both observables then become functions of ONE model quantity, the orbital
// angular momentum at fixed semi-major axis, and the test is whether they
// pick the same L.
//
// Usage: two_observable_consistency [action J/hbar, default 0.0985]
//        [semi-major axis in a_pair, default 1]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/two_observable_consistency.cpp -o /tmp/consistency $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>

// The model has no Boltzmann constant of its own, because it has no
// reservoir; this one exists only to quote the measure's fitted scale in
// kelvin, as section 112b did.
constexpr double boltzmannForQuoting=1.380649e-23;

int main(int argc,char** argv) {
    const double action=(argc>1?atof(argv[1]):0.0985)*hbar;
    const double axis=(argc>2?atof(argv[2]):1.0)*pairBohrRadius(activePair);
    const double planck=2.0*pi*hbar;
    const double measuredSplittingGHz=203.3941;
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 first=normal*firstMagneticMoment;
    const Vec3 second=normal*secondMagneticMoment;
    // A(r): half the splitting analogue, i.e. the amplitude of the model's own
    // angle-dependent energy, read off the model's function on the two
    // collinear configurations rather than from a formula.  Its SIGN is kept:
    // positive means the aligned configuration is the higher one, which is the
    // sign that polarizes the angle anti-aligned and is what 1:3 needs.
    const auto amplitude=[&](double separation) {
        const double up=azimuthAveragedDipoleEnergy(
            separation,first,second,normal);
        const double down=azimuthAveragedDipoleEnergy(
            separation,first,second*(-1.0),normal);
        return 0.5*(up-down);
    };
    // Orbit average at fixed semi-major axis, over the eccentric anomaly with
    // the time weight (1 - e cos E), as sections 116d and 117c did.
    const auto orbitAverage=[&](double angular,auto&& quantity) {
        const double eccentricity=
            std::sqrt(std::max(0.0,1.0-angular*angular));
        const int samples=20001;
        double sum=0.0,weight=0.0;
        for(int index=0;index<samples;++index) {
            const double anomaly=2.0*pi*index/(samples-1.0);
            const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
            const double timeWeight=1.0-eccentricity*std::cos(anomaly);
            sum+=quantity(radius)*timeWeight;
            weight+=timeWeight;
        }
        return sum/weight;
    };
    const auto averagedAmplitude=[&](double angular) {
        return orbitAverage(angular,amplitude);
    };
    // The libration's angular velocity: the model's own orbit-averaged BMT
    // rates for the two moments, differenced as vectors.  The mutual angle
    // responds to the RELATIVE precession, so |omega1 - omega2| is its rate;
    // section 129 measured the resulting frequency on a full trace and rule
    // R130d checks this against it.
    const auto librationRate=[&](double angular) {
        const Vec3 orbital=normal*(angular*hbar);
        const OrbitAveragedBmtAngularVelocities rates=
            orbitAveragedBmtAngularVelocities(axis,orbital,first,second,
                                              pairReducedMass);
        if(!rates.valid) return 0.0;
        return (rates.first-rates.second).norm();
    };
    // Langevin: <cos> = coth(x) - 1/x for p ~ exp(x cos), and here
    // x = -A/E_lib, so a POSITIVE A (aligned configuration higher) gives a
    // negative <cos>, which is the direction 1:3 needs.
    const auto langevin=[](double x) {
        if(std::abs(x)<1.0e-8) return x/3.0;
        return 1.0/std::tanh(x)-1.0/x;
    };
    struct Row { double angular,amplitudeEv,splitting,librationGHz,
                 energyEv,x,cosine,weight; };
    const auto evaluate=[&](double angular) {
        Row row{};
        row.angular=angular;
        const double amplitudeJoule=averagedAmplitude(angular);
        row.amplitudeEv=amplitudeJoule/elementaryCharge;
        row.splitting=std::abs(2.0*amplitudeJoule)/planck*1.0e-9;
        const double rate=librationRate(angular);
        row.librationGHz=rate/(2.0*pi)*1.0e-9;
        row.energyEv=action*rate/elementaryCharge;
        row.x=row.energyEv>0.0?-row.amplitudeEv/row.energyEv:0.0;
        row.cosine=langevin(row.x);
        row.weight=0.5*(1.0+row.cosine);
        return row;
    };
    std::printf("semi-major axis %.4f a_pair; libration action J = %.4f hbar "
                "(audit 129)\n",axis/pairBohrRadius(activePair),action/hbar);
    std::printf("O1 asks for E[w] = 0.2500, O2 for a splitting of %.4f GHz\n\n",
                measuredSplittingGHz);
    std::printf("%8s %14s %14s %14s %14s %9s %9s %8s\n","L/hbar",
                "A [eV]","2A [GHz]","f_lib [GHz]","E_lib [eV]","x","<cos>",
                "E[w]");
    for(double angular:{1.0,0.7,0.5,0.3,0.2,0.1,0.05,0.03,0.01}) {
        const Row row=evaluate(angular);
        std::printf("%8.2f %+14.6e %14.4f %14.4f %+14.6e %9.3f %+9.4f %8.4f\n",
                    row.angular,row.amplitudeEv,row.splitting,
                    row.librationGHz,row.energyEv,row.x,row.cosine,row.weight);
    }
    // R130d: the calibration against the measured libration of section 129.
    const Row production=evaluate(1.0);
    std::printf("\nR130d calibration: secular f_lib at L = hbar is %.4f GHz "
                "against the traced %.4f GHz, ratio %.3f\n",
                production.librationGHz,90.551,
                production.librationGHz/90.551);
    // R130c: which way the model's own energy polarizes the angle.
    std::printf("R130c sign at the production geometry: A = %+.6e eV, so the "
                "measure polarizes the angle %s\n",production.amplitudeEv,
                production.amplitudeEv>0.0?"ANTI-ALIGNED (what 1:3 needs)"
                                          :"ALIGNED (the wrong way)");
    // Locate each observable's L by bisection.  Both A and f_lib rise as L
    // falls, so both E[w] (for the anti-aligning sign) and the splitting are
    // monotone there; the search brackets are checked before use.
    const auto bisect=[&](auto&& value,double target,double low,double high) {
        const double atLow=value(low),atHigh=value(high);
        if((atLow-target)*(atHigh-target)>0.0) return -1.0;
        for(int iteration=0;iteration<80;++iteration) {
            const double middle=0.5*(low+high);
            const bool sameSideAsLow=(value(middle)-target)*(atLow-target)>0.0;
            if(sameSideAsLow) low=middle; else high=middle;
        }
        return 0.5*(low+high);
    };
    const double weightL=bisect([&](double a){ return evaluate(a).weight; },
                                0.25,0.005,1.0);
    const double splittingL=bisect(
        [&](double a){ return evaluate(a).splitting; },
        measuredSplittingGHz,0.005,1.0);
    std::printf("\nL(O1), where E[w] = 0.2500:      ");
    if(weightL>0.0) {
        const Row row=evaluate(weightL);
        std::printf("L = %.4f hbar; there 2A = %.4f GHz (factor %.3f from "
                    "%.4f)\n",weightL,row.splitting,
                    row.splitting>0.0?std::max(row.splitting,
                        measuredSplittingGHz)/std::min(row.splitting,
                        measuredSplittingGHz):0.0,measuredSplittingGHz);
    } else std::printf("not crossed in 0.005..1 hbar\n");
    std::printf("L(O2), where 2A = %.4f GHz: ",measuredSplittingGHz);
    if(splittingL>0.0) {
        const Row row=evaluate(splittingL);
        std::printf("L = %.4f hbar; there E[w] = %.4f (factor %.3f from "
                    "0.2500)\n",splittingL,row.weight,
                    row.weight>0.0?std::max(row.weight,0.25)
                        /std::min(row.weight,0.25):0.0);
    } else std::printf("not crossed in 0.005..1 hbar\n");
    if(weightL>0.0&&splittingL>0.0)
        std::printf("R130b L(O1)/L(O2) = %.3f\n",weightL/splittingL);
    // WHY x IS L-INDEPENDENT, checked rather than argued: the precession rate
    // is omega = gamma B and the energy amplitude is A = mu B in the SAME
    // field, so hbar omega / A = hbar gamma / mu = 2 exactly up to the
    // geometry of the two terms, and the orbit cancels out.  Printed here as
    // the dimensionless ratio, and as the spread of E[w] over the outer
    // decade of L, so that the cancellation is visible in the numbers.
    std::printf("\nhbar omega_lib / A at L = 1.00, 0.50, 0.20, 0.10: "
                "%.4f %.4f %.4f %.4f\n",
                hbar*2.0*pi*evaluate(1.0).librationGHz*1.0e9
                    /elementaryCharge/evaluate(1.0).amplitudeEv,
                hbar*2.0*pi*evaluate(0.5).librationGHz*1.0e9
                    /elementaryCharge/evaluate(0.5).amplitudeEv,
                hbar*2.0*pi*evaluate(0.2).librationGHz*1.0e9
                    /elementaryCharge/evaluate(0.2).amplitudeEv,
                hbar*2.0*pi*evaluate(0.1).librationGHz*1.0e9
                    /elementaryCharge/evaluate(0.1).amplitudeEv);
    {
        double lowest=1.0,highest=0.0;
        for(int index=0;index<=90;++index) {
            const double angular=0.1+0.9*index/90.0;
            const double weight=evaluate(angular).weight;
            lowest=std::min(lowest,weight); highest=std::max(highest,weight);
        }
        std::printf("E[w] over L = 0.1..1.0 (91 nodes): %.4f..%.4f, spread "
                    "%.1e -- O1 does not select an orbit at all\n",
                    lowest,highest,highest-lowest);
    }
    // ROBUSTNESS.  Two measured actions (the first and the last libration
    // cycle of 129) against two definitions of the libration frequency (the
    // secular rate used above, and the one traced in 129, which R130d found
    // to be 2.06 times larger).  All four are parameter-free once J is taken
    // from 129; the box shows how much of the verdict rides on the choice.
    std::printf("\nrobustness box, E[w] (all L-independent):\n");
    std::printf("%22s %16s %16s\n","","secular omega","traced omega");
    for(double actionRatio:{0.0985,0.1383}) {
        const double reference=evaluate(1.0).x*0.0985/actionRatio;
        const double traced=reference/(90.551/evaluate(1.0).librationGHz);
        std::printf("  J = %.4f hbar %*s %16.4f %16.4f\n",actionRatio,4,"",
                    0.5*(1.0+langevin(reference)),
                    0.5*(1.0+langevin(traced)));
    }
    // The band of actions that satisfies O1 within the factor 1.5 of R130a,
    // so that the measured J can be placed inside or outside it.
    {
        const double baseX=evaluate(1.0).x;
        const auto weightForAction=[&](double actionRatio) {
            return 0.5*(1.0+langevin(baseX*0.0985/actionRatio));
        };
        const auto solve=[&](double target) {
            double low=0.01,high=3.0;
            for(int iteration=0;iteration<80;++iteration) {
                const double middle=0.5*(low+high);
                if((weightForAction(middle)-target)
                   *(weightForAction(low)-target)>0.0) low=middle;
                else high=middle;
            }
            return 0.5*(low+high);
        };
        std::printf("O1 within a factor 1.5 (E[w] in 0.1667..0.3750) needs "
                    "J in %.4f..%.4f hbar; 129 measured %.4f..%.4f\n",
                    solve(0.375),solve(1.0/6.0),0.0985,0.1383);
        std::printf("O1 exactly (E[w] = 0.2500) needs J = %.4f hbar\n",
                    solve(0.25));
    }

    // The three candidates of 112b, for the free-parameter count of R130e.
    std::printf("\n112b candidates, with the scale each one needs:\n");
    std::printf("  two points cos = +-1 weighted 1:3   <cos> = -1/2 by "
                "construction; the weights are degeneracies, so the measure "
                "implies NO splitting at all\n");
    std::printf("  p ~ (1-cos)^2 = |S_total|^4         scale-free, so it also "
                "implies no splitting\n");
    const double langevinScale=1.7965;
    const double requiredAmplitude=
        0.5*planck*measuredSplittingGHz*1.0e9/elementaryCharge;
    std::printf("  p ~ exp(-1.7965 cos)                needs kT = A/1.7965; "
                "with the model's own A(L = hbar) = %.6e eV that is %.6e eV "
                "= %.4f K\n",production.amplitudeEv,
                std::abs(production.amplitudeEv)/langevinScale,
                std::abs(production.amplitudeEv)/langevinScale
                    *elementaryCharge/boltzmannForQuoting);
    std::printf("  the same measure fitted to O2 instead needs A = %.6e eV, "
                "i.e. %.2f times the model's own at L = hbar\n",
                requiredAmplitude,
                requiredAmplitude/std::abs(production.amplitudeEv));
}
