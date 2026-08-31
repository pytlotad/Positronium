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
bool gVisualSimulationComplete = false;
TButton* gStopButton = nullptr;
TCanvas* gVisualCanvas = nullptr;
TPaveText* gVisualObservationBox = nullptr;
int gVisualPhenomenon = 0;
bool gVisualExitSaveAttempted = false;

void SetVisualObservationStatus(const char* headline,const char* detail,
                                int color) {
    if(!gVisualObservationBox) return;
    gVisualObservationBox->Clear();
    gVisualObservationBox->SetTextColor(color);
    gVisualObservationBox->AddText(headline);
    gVisualObservationBox->AddText(detail);
}

void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
void ExitSimulation() {
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
