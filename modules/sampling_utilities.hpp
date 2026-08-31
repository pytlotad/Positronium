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

std::uint64_t splitMix64(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

// The default ceiling of 30 suits the trajectory experiments, whose sample
// sizes are limited by integration cost.  The annihilation generator is cheap
// enough to run millions of events, so its kinematics histograms pass a larger
// ceiling instead of throwing that resolution away.
int histogramBins(size_t count, int maximumBins = 30) {
    return std::clamp(static_cast<int>(std::lround(2.0 * std::sqrt(count))),
                      6, maximumBins);
}
