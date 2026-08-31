#pragma once

// Constrained unbinned MLE for the second-Legendre anisotropy coefficient
// a2 of a sample of cos(theta) values, under the likelihood
// product_i [1 + a2 P2(cos(theta_i))], -1 <= a2 <= 2.  Validation-only: used
// by the annihilation-generator unit test.  The production panels plot the
// exact isotropic reference instead of fitting a sampled one.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes; this
// slice is validation-only rather than engine/experiment/presentation, so
// it is entirely POSITRONIUM_ENABLE_FIELD_VALIDATION-gated, same as in the
// original).  Textually included at the same point inside positronium.cpp's
// shared anonymous namespace it always occupied.  Not yet a standalone,
// order-independent header.

#ifdef POSITRONIUM_ENABLE_FIELD_VALIDATION
struct LegendreFitSummary {
    double anisotropy = std::numeric_limits<double>::quiet_NaN();
    double standardError = std::numeric_limits<double>::quiet_NaN();
    std::size_t count = 0;
    bool atBoundary = false;
};

LegendreFitSummary fitSecondLegendreAnisotropy(
    const std::vector<double>& cosines) {
    LegendreFitSummary result;
    std::vector<double> secondLegendre;
    secondLegendre.reserve(cosines.size());
    for (double cosine : cosines) {
        if (!std::isfinite(cosine)) continue;
        secondLegendre.push_back(0.5*(3.0*cosine*cosine - 1.0));
    }
    result.count = secondLegendre.size();
    if (result.count == 0) return result;

    // The event-level likelihood is proportional to
    // product_i [1 + a2 P2(cos(theta_i))], with the physical range
    // -1 <= a2 <= 2.  Its derivative is monotone, so bisection gives the
    // constrained unbinned MLE without treating the two back-to-back photons
    // as independent observations.
    constexpr double lower = -1.0;
    constexpr double upper = 2.0;
    constexpr double guard = 64.0*std::numeric_limits<double>::epsilon();
    const auto derivative = [&](double anisotropy) {
        double value = 0.0;
        for (double legendre : secondLegendre) {
            value += legendre/(1.0 + anisotropy*legendre);
        }
        return value;
    };
    const double safeLower = lower + guard;
    const double safeUpper = upper - guard;
    if (derivative(safeLower) <= 0.0) {
        result.anisotropy = lower;
        result.atBoundary = true;
    } else if (derivative(safeUpper) >= 0.0) {
        result.anisotropy = upper;
        result.atBoundary = true;
    } else {
        double left = safeLower;
        double right = safeUpper;
        for (int iteration = 0; iteration < 100; ++iteration) {
            const double middle = 0.5*(left + right);
            if (derivative(middle) > 0.0) left = middle;
            else right = middle;
        }
        result.anisotropy = 0.5*(left + right);
        double information = 0.0;
        for (double legendre : secondLegendre) {
            const double scaled = legendre
                /(1.0 + result.anisotropy*legendre);
            information += scaled*scaled;
        }
        if (result.count >= 2 && information > 0.0) {
            result.standardError = 1.0/std::sqrt(information);
        }
    }
    return result;
}
#endif
