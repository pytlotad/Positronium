// What depends on the moment radius (audit section 139)?
//
// The model uses magneticDipoleRadius() = 0.9668 r*, derived from a current
// loop rather than fitted, and since audit 126 the SAME length also softens
// the charge side of every moment-charge interaction (matchedMomentSoftening).
// Section 118 measured that the contact term would need 452 r* instead.  This
// probe does not test either number: it maps which of the model's results move
// when the constant does, so that any future change to it has a known blast
// radius.
//
// CREM_MAGNETIC_RADIUS_SCALE sets the radius in units of r*, and it is read
// once into a static, so one process per value.
//
// Usage: moment_radius_sensitivity [budget s, default 300] [seed, default 42]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/moment_radius_sensitivity.cpp -o /tmp/radius $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const double budget=argc>1?atof(argv[1]):300.0;
    const unsigned long long seed=argc>2?strtoull(argv[2],nullptr,10):42ULL;
    const double axis=pairBohrRadius(activePair);
    const double planck=2.0*pi*hbar;
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 first=normal*firstMagneticMoment;
    const Vec3 second=normal*secondMagneticMoment;
    // The model's own splitting analogue at a_pair and at 2 r*, where the
    // softening actually bites.
    const auto splitting=[&](double separation) {
        const double up=azimuthAveragedDipoleEnergy(
            separation,first,second,normal);
        const double down=azimuthAveragedDipoleEnergy(
            separation,first,second*(-1.0),normal);
        return std::abs(up-down)/planck*1.0e-9;
    };
    const auto started=std::chrono::steady_clock::now();
    const auto result=runCremCollapseExperiment(seed,1,1,budget);
    const double wall=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();
    std::printf("RADIUS scale=%.6f r*  matched_softening=%.6f r*  "
                "split_a_pair=%.6f GHz  split_2rstar=%.6e GHz  "
                "split_1rstar=%.6e GHz  transit_ps=%.9f  cascade_ps=%.9f  "
                "outcome=%d  stop=%d  photons=%zu  revolutions=%.6e  wall=%.1f\n",
                magneticDipoleRadius()/comptonBarrierRadius,
                matchedMomentSoftening()/comptonBarrierRadius,
                splitting(axis),splitting(2.0*comptonBarrierRadius),
                splitting(comptonBarrierRadius),
                result[0].collapseTransitSeconds*1.0e12,
                result[0].lifetimeSecondsLab*1.0e12,
                static_cast<int>(result[0].calibrationOutcome),
                static_cast<int>(result[0].stopCause),
                result[0].labFramePhotons.size(),
                result[0].revolutions,wall);
}
