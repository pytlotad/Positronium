#pragma once

// Visual-mode GUI callbacks and the flags they and the animation loop share.
// Deliberately at GLOBAL scope, not inside the engine's anonymous
// namespace: ROOT's TButton invokes ToggleSimulation through the Cling
// interpreter by name, which needs ordinary external linkage.
//
// Extracted verbatim from positronium.cpp (the presentation side of
// splitting engine/experiments/ROOT presentation apart -- see the session
// notes; this is the first slice that is not "engine").  Included at the
// same point positronium.cpp always had it, before the anonymous namespace
// opens, so it depends on root_export.hpp (already included earlier) and
// ROOT's TButton/TCanvas/TApplication (gApplication) being visible, which
// they are from the top-of-file #include block.

// These functions are intentionally global: ROOT's TButton invokes its action
// through the interpreter while the animation loop observes these flags.
bool gSimulationPaused = false;
bool gExitRequested = false;
TButton* gStopButton = nullptr;
TCanvas* gVisualCanvas = nullptr;
int gVisualPhenomenon = 0;
bool gVisualExitSaveAttempted = false;

void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
void ExitSimulation() {
    gExitRequested = true;
    if(gVisualCanvas&&gVisualPhenomenon>=1&&gVisualPhenomenon<=4) {
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
