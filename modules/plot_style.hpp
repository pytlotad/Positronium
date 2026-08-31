#pragma once

// The colour convention every plotted panel follows (what kind of statement
// a line/marker/legend entry makes, encoded by colour so it survives a
// glance and, via the Okabe-Ito palette, colour-vision deficiency and
// greyscale printing too), plus the shared ROOT histogram styling built on
// top of it.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.

// ---------------------------------------------------------------------------
// Plot colour convention.
//
// Every line, marker set and legend entry is coloured by WHERE ITS NUMBERS COME
// FROM, so one glance at any panel says what kind of statement it makes:
//
//   CREM simulation     blue         trajectories integrated by this program
//   sampled input       orange       quantities drawn from a random distribution
//   experimental        bluish green published measurements
//   theory / analytic   vermillion   closed-form reference curves
//
// The hues are the Okabe-Ito colour-vision-deficiency-safe palette rather than
// plain red/green/blue.  Plain red and green are precisely the pair that the
// roughly 8% of men with deuteranopia cannot separate, and here they would have
// carried the two most important comparisons.  Line and marker STYLE repeats
// the same information, so the panels survive greyscale printing too.
//
// Fill colour is reserved for classification (experiment 5) and uses pale tints
// that cannot be mistaken for the saturated provenance outlines.
namespace plot_style {

int crem()         { return TColor::GetColor("#0072B2"); }
int sampled()      { return TColor::GetColor("#E69F00"); }
int experimental() { return TColor::GetColor("#009E73"); }
int theory()       { return TColor::GetColor("#D55E00"); }
int qed()          { return TColor::GetColor("#CC79A7"); }

// Pale classification fills, in the order of InteractionOutcome.
int classificationFill(std::size_t slot) {
    static const std::array<const char*,6> tints{
        "#E8A0A0",  // Collision
        "#8FC3E8",  // Scattering
        "#8FD1B0",  // Para-Positronium
        "#C9A0DC",  // Ortho-Positronium
        "#BEBEBE",  // Unresolved
        "#EBC17A"}; // NumericalFailure
    return TColor::GetColor(tints[std::min(slot, tints.size() - 1)]);
}

// Legend text naming only the provenances a given panel actually contains, so
// the key stays short and never claims something the panel does not show.
std::string key(bool hasCrem, bool hasSampled, bool hasExperimental,
                bool hasTheory) {
    std::string result = "colours: ";
    bool first = true;
    const auto add = [&](bool present, const char* text) {
        if (!present) return;
        if (!first) result += ", ";
        result += text;
        first = false;
    };
    add(hasCrem, "blue CREM");
    add(hasSampled, "orange sampled");
    add(hasExperimental, "green measured");
    add(hasTheory, "vermillion theory");
    return result;
}

} // namespace plot_style

void styleHistogram(TH1D& histogram, int color) {
    histogram.SetDirectory(nullptr);
    histogram.SetLineColor(color);
    histogram.SetFillColorAlpha(color, 0.45);
    histogram.SetLineWidth(2);
}

// Prepares a histogram to have its own bin content drawn as an integer label
// above each occupied bar, via ROOT's "TEXT0" draw option (appended by the
// caller alongside "HIST"; TEXT0 skips empty bins so a wide, mostly-empty
// axis does not fill up with zeroes). Used on every histogram in this file
// whose bins are literal trajectory/event counts -- not on the handful that
// hold a derived physical quantity (differential cross sections, weighted
// spectra), where a bin's "count" is not the number being reported.
void styleBinCounts(TH1D& histogram) {
    histogram.SetMarkerSize(1.3);
    // TH1 has no per-histogram text format; ROOT reads it from TStyle at
    // paint time, so this is a (harmless, idempotent) global style edit.
    // No leading '%': ROOT prepends its own, and "%.0f" here would paint the
    // literal string "%.0f" instead of formatting anything (found by looking
    // at the actual rendered PDF, not just trusting the call compiled).
    gStyle->SetPaintTextFormat(".0f");
}
