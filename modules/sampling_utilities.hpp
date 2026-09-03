#pragma once

// Small, unrelated-to-each-other utilities shared by the statistical
// experiments: a SplitMix64 seed-mixer (for deriving independent RNG
// streams from one seed) and a histogram bin-count heuristic.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>

constexpr std::uint64_t splitMix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// The default ceiling of 30 suits the trajectory experiments, whose sample
// sizes are limited by integration cost.  The annihilation generator is cheap
// enough to run millions of events, so its kinematics histograms pass a larger
// ceiling instead of throwing that resolution away.
inline int histogramBins(size_t count, int maximumBins = 30) {
    return std::clamp(static_cast<int>(std::lround(2.0 * std::sqrt(count))),
                      6, maximumBins);
}

// A circular, isotropic 2D Gaussian beam has independent Gaussian transverse
// coordinates.  Its radial impact parameter is therefore Rayleigh-distributed
// (the area Jacobian contributes the factor b), not half-normal.
struct IsotropicGaussianImpactSample {
    double transverseY = std::numeric_limits<double>::quiet_NaN();
    double transverseZ = std::numeric_limits<double>::quiet_NaN();
    double impactParameter = std::numeric_limits<double>::quiet_NaN();

    bool valid(double maximumImpactParameter) const {
        return std::isfinite(transverseY)
            && std::isfinite(transverseZ)
            && std::isfinite(impactParameter)
            && impactParameter >= 0.0
            && impactParameter <= maximumImpactParameter;
    }
};

inline IsotropicGaussianImpactSample sampleIsotropicGaussianImpact(
        std::mt19937_64& random, double transverseSigma,
        double maximumImpactParameter) {
    if (!(transverseSigma > 0.0)
        || !std::isfinite(transverseSigma)
        || !(maximumImpactParameter >= 0.0)
        || !std::isfinite(maximumImpactParameter)) {
        return {};
    }
    std::normal_distribution<double> transverseGaussian(0.0,transverseSigma);
    for(int attempt=0;attempt<1000;++attempt) {
        const double transverseY=transverseGaussian(random);
        const double transverseZ=transverseGaussian(random);
        const double impactParameter=std::hypot(transverseY,transverseZ);
        if (impactParameter<=maximumImpactParameter) {
            return {transverseY,transverseZ,impactParameter};
        }
    }
    return {};
}
