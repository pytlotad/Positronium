#pragma once

// Streaming (Welford) maximum-likelihood Gaussian fit: mean and population
// sigma from a sample, skipping non-finite values, in one pass with no
// separate accumulation of the raw values.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.

struct GaussianFitSummary {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double sigma = std::numeric_limits<double>::quiet_NaN();
    std::size_t count = 0;
};

GaussianFitSummary gaussianMaximumLikelihood(const std::vector<double>& values) {
    GaussianFitSummary result;
    result.mean = 0.0;
    double secondMoment = 0.0;
    for (double value : values) {
        if (!std::isfinite(value)) continue;
        ++result.count;
        const double delta = value - result.mean;
        result.mean += delta/static_cast<double>(result.count);
        const double updatedDelta = value - result.mean;
        secondMoment += delta*updatedDelta;
    }
    if (result.count == 0) {
        result.mean = std::numeric_limits<double>::quiet_NaN();
        return result;
    }
    result.sigma = std::sqrt(std::max(0.0, secondMoment)
        / static_cast<double>(result.count));
    return result;
}
