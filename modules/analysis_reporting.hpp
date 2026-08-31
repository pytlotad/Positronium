#pragma once

// The analysis-panel reporting helpers used across the statistical
// experiments: a compact number formatter, the colour-tinted analysis-box
// text lines and the ROOT TPaveText they draw into, a Gaussian MLE overlay
// curve for a histogram, and the two "did the save succeed" reporters for
// plot and archive exports.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.
// Depends on that namespace already having in scope: GaussianFitSummary
// (modules/gaussian_fit.hpp, included earlier), root_export::ExportResult,
// statistics_archive::OperationResult.

std::string compactNumber(double value, int precision = 5) {
    std::ostringstream output;
    output << std::setprecision(precision) << value;
    return output.str();
}

// One line of an analysis box, tinted to match the provenance colour
// convention (plot_style::crem/sampled/experimental/theory) of whichever
// series it describes -- so a reader can tell what a number is FOR without
// tracing it back to a legend.  Converts implicitly from a bare string
// literal, defaulting to kBlack, so a line with no single-series meaning
// (counts, captions, warnings) needs no change at the call site.
struct AnalysisLine {
    std::string text;
    int color;
    AnalysisLine(std::string t, int c = kBlack) : text(std::move(t)), color(c) {}
    // A bare string literal is a const char*, and reaching AnalysisLine from
    // one needs a SINGLE user-defined conversion -- going through the
    // std::string overload above would be two (const char*->string->
    // AnalysisLine) chained implicit conversions, which C++ does not allow,
    // and every brace-enclosed lines list below mixes literals with coloured
    // AnalysisLine(...) entries.
    AnalysisLine(const char* t, int c = kBlack) : text(t), color(c) {}
};

TPaveText* drawAnalysisBox(std::vector<std::unique_ptr<TPaveText>>& storage,
                           double x1, double y1, double x2, double y2,
                           const std::vector<AnalysisLine>& lines,
                           double textSize = 0.027) {
    auto box = std::make_unique<TPaveText>(x1, y1, x2, y2, "NDC");
    box->SetFillColorAlpha(kWhite, 0.84);
    box->SetTextAlign(12);
    box->SetTextFont(42);
    box->SetTextSize(textSize);
    for (const AnalysisLine& line : lines) {
        box->AddText(line.text.c_str())->SetTextColor(line.color);
    }
    box->Draw();
    TPaveText* result = box.get();
    storage.push_back(std::move(box));
    return result;
}

std::unique_ptr<TF1> gaussianMleOverlay(const std::string& name,
                                        const GaussianFitSummary& fit,
                                        const TH1D& histogram, int color) {
    if (fit.count < 2 || !(fit.sigma > 0.0) || !std::isfinite(fit.sigma)) {
        return nullptr;
    }
    const double lower = histogram.GetXaxis()->GetXmin();
    const double upper = histogram.GetXaxis()->GetXmax();
    auto function = std::make_unique<TF1>(name.c_str(),
        "[0]*exp(-0.5*((x-[1])/[2])^2)", lower, upper);
    const double amplitude = static_cast<double>(fit.count)
        * histogram.GetBinWidth(1)/(std::sqrt(2.0*pi)*fit.sigma);
    function->SetParameters(amplitude, fit.mean, fit.sigma);
    function->SetLineColor(color);
    function->SetLineWidth(2);
    function->SetLineStyle(2);
    function->Draw("SAME");
    return function;
}

void reportExports(const std::vector<root_export::ExportResult>& results) {
    for (const root_export::ExportResult& result : results) {
        if (result) {
            std::cout << "Saved plot: " << result.path.string() << '\n';
        } else {
            std::cerr << "Warning: could not save " << result.path.string()
                      << ": " << result.error << '\n';
        }
    }
}

bool reportArchiveOperation(const statistics_archive::OperationResult& result,
                            const char* description) {
    if (result) {
        std::cout << "Saved " << description << ": "
                  << result.path.string() << '\n';
        return true;
    }
    std::cerr << "Warning: could not save " << description << ' '
              << result.path.string() << ": " << result.error << '\n';
    return false;
}
