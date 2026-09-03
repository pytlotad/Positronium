#pragma once

// Visual-mode GUI callbacks and the flags they and the animation loop share.
// Deliberately at GLOBAL scope, not inside the engine's anonymous
// namespace: ROOT's TButton invokes ToggleSimulation through the Cling
// interpreter by name, which needs ordinary external linkage.
//
// Self-contained and order-independent: it pulls in root_export.hpp and the
// ROOT headers for TButton/TCanvas/TPaveText/TApplication itself rather than
// relying on positronium.cpp's top-of-file block.
//
// The flags and callbacks stay at namespace scope with EXTERNAL linkage, which
// `inline` preserves -- Cling resolves ToggleSimulation and ExitSimulation by
// name when TButton fires them, and cannot see an internal-linkage symbol.

#include "root_export.hpp"

#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <Rtypes.h>

#include <iostream>

// These functions are intentionally global: ROOT's TButton invokes its action
// through the interpreter while the animation loop observes these flags.
inline bool gSimulationPaused = false;
inline bool gExitRequested = false;
inline bool gVisualSimulationComplete = false;
inline TButton* gStopButton = nullptr;
inline TCanvas* gVisualCanvas = nullptr;
inline TPaveText* gVisualObservationBox = nullptr;
inline int gVisualPhenomenon = 0;
inline bool gVisualExitSaveAttempted = false;

inline void SetVisualObservationStatus(const char* headline,
                                       const char* detail, int color) {
    if(!gVisualObservationBox) return;
    gVisualObservationBox->Clear();
    gVisualObservationBox->SetTextColor(color);
    gVisualObservationBox->AddText(headline);
    gVisualObservationBox->AddText(detail);
}

inline void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
inline void ExitSimulation() {
    gExitRequested = true;
    if(gVisualCanvas&&gVisualPhenomenon>=1&&gVisualPhenomenon<=4) {
        // EXIT censors an observation only while integration is still in
        // progress.  Once the engine has returned, this button merely closes
        // an already classified result and must not overwrite it.
        if(!gVisualSimulationComplete) {
            SetVisualObservationStatus(
                "Observation: ADMINISTRATIVELY CENSORED",
                "User EXIT before the terminal endpoint was observed",kOrange+7);
        }
        gVisualExitSaveAttempted=true;
        gVisualCanvas->Modified();
        gVisualCanvas->Update();
        const root_export::ExportResult screenshot=
            root_export::saveVisualScreenshot(
                *gVisualCanvas,gVisualPhenomenon);
        if(screenshot) {
            std::cout<<"Saved visual screenshot: "
                     <<screenshot.path.string()<<'\n';
        } else {
            std::cerr<<"Warning: could not save visual screenshot "
                     <<screenshot.path.string()<<": "
                     <<screenshot.error<<'\n';
        }
    }
    if (gApplication) gApplication->Terminate(0);
}
