// Is omega1 = omega2 in the charge sector a property of the SECTOR or of the
// e+e- pair (audit section 134)?
//
// Section 116b measured that keeping only the partner's charge field in the
// precession conserves |S1+S2| exactly, and built an argument on it -- on one
// pair, with equal masses and opposite charges.  In the centre-of-mass frame
// that pairing makes v1 = -v2 and r_hat12 = -r_hat21, so the two fields are
// opposite, the two q_i/m_i are opposite, and the rates coincide.  The
// equality should therefore survive for mu+mu- and fail for p+e-, whose mass
// ratio is 1836.  This tool runs one trajectory per pair so the distinction
// is measured rather than argued.
//
// Read it with CREM_DEBUG_PROJECTION on and, for the ablation,
// CREM_SPIN_RATE_SECTOR=charge; the PROJ lines carry Spair_hbar, |w1|, |w2|
// and cos12, which is what the tables of 116b and 117d were built from.
//
// Usage: spin_rate_symmetry <pair, e.g. muon,antimuon> [budget s] [seed]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/spin_rate_symmetry.cpp -o /tmp/symmetry $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc,char** argv) {
    const std::string pair=argc>1?argv[1]:"positron,electron";
    const double budget=argc>2?atof(argv[2]):300.0;
    const unsigned long long seed=argc>3?strtoull(argv[3],nullptr,10):42ULL;
    applyPairFromOption(pair);
    // One trajectory in the trace, not two: the collapse-transit measurement
    // of audit 108 would run the whole thing a second time and interleave its
    // PROJ lines with the first run's.
    // CREM_COLLAPSE_TRANSIT=1 restores it, which is what the production path
    // does and therefore what the traces of audit 116b were taken under:
    // there the estimator runs the trajectory a SECOND time under continuous
    // reaction, and the PROJ lines of both runs land in one stream.  Rule
    // R134e checks the two settings against each other.
    if(!std::getenv("CREM_COLLAPSE_TRANSIT")) gMeasureCollapseTransit=false;
    const auto started=std::chrono::steady_clock::now();
    const auto result=runCremCollapseExperiment(seed,1,1,budget);
    const double wall=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-started).count();
    const char* sector=std::getenv("CREM_SPIN_RATE_SECTOR");
    std::printf("RESULT pair %s sector %s seed %llu a_pair %.6e m "
                "mass_ratio %.4f moment_ratio %.6f outcome %d stop_cause %d "
                "photons %zu wall_s %.1f\n",
                pair.c_str(),sector?sector:"both",seed,
                pairBohrRadius(activePair),
                firstMass/secondMass,
                firstMagneticMoment/secondMagneticMoment,
                static_cast<int>(result[0].calibrationOutcome),
                static_cast<int>(result[0].stopCause),
                result[0].labFramePhotons.size(),wall);
}
