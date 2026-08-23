#pragma once

#include <TCanvas.h>
#include <TError.h>
#include <TFrame.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

// Small, header-only helpers for exporting ROOT pads.  They deliberately write
// a complete PDF beside the destination first and only then replace the old
// file.  A failed render therefore does not truncate the previous good image.
namespace root_export {

struct ExportResult {
    std::filesystem::path path;
    std::string error;

    [[nodiscard]] bool succeeded() const noexcept { return error.empty(); }
    explicit operator bool() const noexcept { return succeeded(); }
};

struct NamedPad {
    TVirtualPad* pad = nullptr; // non-owning; must outlive saveStatisticalPlots()
    int screenNumber = 0;
    // 'a' for a plot of INPUT/sampled data (properties of the prepared or
    // drawn state, before any dynamics act on it -- e.g. a dipole coupling
    // computed from the prepared pair's moments, or a beam parameter drawn
    // for the trial), 'b' for a plot of OUTPUT/result data (anything the
    // simulation itself produced or measured -- collapse times, cross
    // sections, conservation-law residuals, photon spectra).  Required, not
    // defaulted: every call site names it explicitly so a plot's
    // classification is a decision made at the point it is added, not an
    // silently-inherited default.
    char inputOrOutput = '\0';
    int padNumber = 0;
    std::string_view plotName;
};

inline std::string_view phenomenonNameAscii(int selectedPhenomenon) noexcept {
    switch (selectedPhenomenon) {
        case 1: return "para_positronium";
        case 2: return "ortho_positronium";
        case 3: return "direct_collision";
        case 4: return "scattering";
        case 5: return "interactions";
        default: return {};
    }
}

inline bool isAsciiFileStem(std::string_view value) noexcept {
    if (value.empty()) return false;
    for (const unsigned char character : value) {
        const bool lower = character >= 'a' && character <= 'z';
        const bool digit = character >= '0' && character <= '9';
        if (!lower && !digit && character != '_') return false;
    }
    return value.front() != '_' && value.back() != '_';
}

namespace detail {

inline std::filesystem::path temporaryPdfPath(
    const std::filesystem::path& destination) {
    static std::atomic<std::uint64_t> sequence{0};
    const auto tick = std::chrono::steady_clock::now().time_since_epoch().count();
    const auto serial = sequence.fetch_add(1, std::memory_order_relaxed);
    const std::string name = "." + destination.stem().string() + ".tmp_"
        + std::to_string(tick) + "_" + std::to_string(serial) + ".pdf";
    return destination.parent_path() / name;
}

inline bool hasPdfSignature(const std::filesystem::path& path) {
    constexpr std::array<unsigned char, 5> expected{
        0x25, 0x50, 0x44, 0x46, 0x2d};
    std::array<unsigned char, expected.size()> bytes{};
    std::ifstream input(path, std::ios::binary);
    input.read(reinterpret_cast<char*>(bytes.data()),
               static_cast<std::streamsize>(bytes.size()));
    return input.gcount() == static_cast<std::streamsize>(bytes.size())
        && bytes == expected;
}

inline std::string temporaryRootObjectName(std::string_view kind) {
    static std::atomic<std::uint64_t> sequence{0};
    const auto serial = sequence.fetch_add(1, std::memory_order_relaxed);
    return "root_export_" + std::string(kind) + "_" + std::to_string(serial);
}

inline void preparePdfExporter() {}

inline bool renderPadToPdf(TVirtualPad& pad,
                           const std::filesystem::path& temporary,
                           std::string& error) {
    pad.Modified();
    pad.Update();

    // Print the original pad directly: cloning a TPad can lose its painted
    // primitives.  The explicit PDF backend avoids SaveAs's format dialog.
    struct RootMessageLevelGuard {
        Int_t previousLevel = gErrorIgnoreLevel;
        RootMessageLevelGuard() { gErrorIgnoreLevel = kWarning; }
        ~RootMessageLevelGuard() { gErrorIgnoreLevel = previousLevel; }
    } messageGuard;
    struct PaperSizeGuard {
        Float_t previousWidth = 0.0F;
        Float_t previousHeight = 0.0F;
        PaperSizeGuard() {
            gStyle->GetPaperSize(previousWidth, previousHeight);
            gStyle->SetPaperSize(27.0F, 19.0F);
        }
        ~PaperSizeGuard() {
            gStyle->SetPaperSize(previousWidth, previousHeight);
        }
    } paperSizeGuard;
    if (dynamic_cast<TCanvas*>(&pad)) {
        pad.Print(temporary.string().c_str(), "pdf");
        return true;
    }
    if (!pad.GetCanvas()) {
        error = "pad has no parent canvas";
        return false;
    }
    // Statistical export runs in ROOT's batch backend. Clone only the selected
    // pad into a temporary full-size canvas: printing its parent canvas would
    // also paint the other overlapping page and all four diagnostic panels.
    const std::string canvasName=temporaryRootObjectName("canvas");
    TCanvas exportCanvas(canvasName.c_str(),pad.GetTitle(),900,700);
    exportCanvas.cd();
    TObject* cloneObject=pad.DrawClone();
    auto* clone=dynamic_cast<TPad*>(cloneObject);
    if(!clone) {
        error="cannot clone ROOT pad for PDF export";
        return false;
    }
    clone->SetPad(0.0,0.0,1.0,1.0);
    clone->SetBit(kCanDelete);
    exportCanvas.Modified();
    exportCanvas.Update();
    exportCanvas.Print(temporary.string().c_str(), "pdf");
    return true;
}

inline ExportResult savePadPdf(TVirtualPad& pad,
                               const std::filesystem::path& destination) {
    ExportResult result{destination, {}};
    try {
        const std::filesystem::path directory = destination.parent_path().empty()
            ? std::filesystem::path(".") : destination.parent_path();
        std::error_code error;
        std::filesystem::create_directories(directory, error);
        if (error || !std::filesystem::is_directory(directory)) {
            result.error = "cannot create output directory: "
                + (error ? error.message() : std::string("path is not a directory"));
            return result;
        }

        const std::filesystem::path temporary = temporaryPdfPath(destination);
        struct TemporaryCleanup {
            std::filesystem::path path;
            ~TemporaryCleanup() {
                std::error_code ignored;
                std::filesystem::remove(path, ignored);
            }
        } cleanup{temporary};

        // TContext restores the caller's gPad after rendering the requested pad.
        TVirtualPad::TContext padContext(&pad);
        if (!renderPadToPdf(pad, temporary, result.error)) return result;

        if (!hasPdfSignature(temporary)) {
            result.error = "ROOT did not produce a valid PDF";
            return result;
        }

        // On POSIX this atomically replaces an existing destination.  If a
        // platform cannot replace it, keep the previous image and report the
        // error instead of deleting the known-good file first.
        std::filesystem::rename(temporary, destination, error);
        if (error) {
            result.error = "cannot replace destination: " + error.message();
            return result;
        }
    } catch (const std::exception& exception) {
        result.error = exception.what();
    }
    return result;
}

} // namespace detail

// Must be called after TApplication construction and before the first GUI
// TCanvas.  This avoids ROOT 6.40's unsafe batch-to-X11 canvas transition.
inline void preparePdfExporter() { detail::preparePdfExporter(); }

inline ExportResult saveVisualScreenshot(
    TVirtualPad& canvas, int selectedPhenomenon,
    const std::filesystem::path& outputDirectory = "distributions") {
    const std::string_view phenomenon = phenomenonNameAscii(selectedPhenomenon);
    if (phenomenon.empty()) {
        return {outputDirectory, "phenomenon number must be from 1 to 4"};
    }
    return detail::savePadPdf(canvas,
        outputDirectory / (std::to_string(selectedPhenomenon)
            + "_1_1_visual_simulation.pdf"));
}

inline std::vector<ExportResult> saveStatisticalPlots(
    int selectedPhenomenon, std::span<const NamedPad> plots,
    const std::filesystem::path& outputDirectory = "distributions") {
    std::vector<ExportResult> results;
    results.reserve(plots.size());
    const std::string_view phenomenon = phenomenonNameAscii(selectedPhenomenon);
    if (phenomenon.empty()) {
        results.push_back(
            {outputDirectory, "phenomenon number must be from 1 to 4"});
        return results;
    }

    for (const NamedPad& plot : plots) {
        if (!plot.pad) {
            results.push_back({outputDirectory, "cannot export a null ROOT pad"});
            continue;
        }
        if (plot.screenNumber < 1 || plot.padNumber < 1) {
            results.push_back({outputDirectory,
                "screen and pad numbers must be positive"});
            continue;
        }
        if (plot.inputOrOutput != 'a' && plot.inputOrOutput != 'b') {
            results.push_back({outputDirectory,
                "inputOrOutput must be 'a' (input/sampled data) or 'b' "
                "(output/result data), set explicitly at the call site"});
            continue;
        }
        if (!isAsciiFileStem(plot.plotName)) {
            results.push_back({outputDirectory,
                "plot name must contain only lowercase ASCII letters, digits, and underscores"});
            continue;
        }
        const std::string filename = std::to_string(selectedPhenomenon) + "_"
            + std::to_string(plot.screenNumber) + "_"
            + std::string(1, plot.inputOrOutput) + "_"
            + std::to_string(plot.padNumber) + "_"
            + std::string(plot.plotName) + ".pdf";
        results.push_back(
            detail::savePadPdf(*plot.pad, outputDirectory / filename));
    }
    return results;
}

inline std::vector<ExportResult> saveStatisticalPlots(
    int selectedPhenomenon, std::initializer_list<NamedPad> plots,
    const std::filesystem::path& outputDirectory = "distributions") {
    return saveStatisticalPlots(selectedPhenomenon,
        std::span<const NamedPad>(plots.begin(), plots.size()), outputDirectory);
}

} // namespace root_export
