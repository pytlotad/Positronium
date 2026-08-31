#pragma once

// The shared interface types between the mechanical engine and everything
// that drives or reports on it (visual rendering, statistical experiment
// orchestration): what one accepted step or terminal state looks like
// (Frame), what a run was asked to do and how (SimulationOptions), what
// kind of phenomenon/outcome/visual style is in play, and what a completed
// run returned (SimulationResult, InitialConditions).
//
// Extracted verbatim from positronium.cpp (Stage 0 of splitting engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied, so it depends on that namespace
// already having in scope: Vec3, State, std::function, std::vector,
// std::uint64_t.  Not yet a standalone, order-independent header.

struct Frame {
    Vec3 first, second, firstDipole, secondDipole;
    Vec3 noetherMomentum, noetherAngularMomentum;
    Vec3 radiatedMomentum, radiatedAngularMomentum;
    double canonicalMomentumScale;
    double time, radius, radiatedEnergy, mechanicalEnergy, schottEnergy;
    double firstMechanicalEnergy, secondMechanicalEnergy;
    double boundFieldEnergy=0.0,reactionEnergyMismatch=0.0;
    Vec3 boundFieldMomentum,boundFieldAngularMomentum;
    Vec3 reactionMomentumMismatch,reactionAngularMomentumMismatch;
};

enum class Phenomenon { DirectCollision, Scattering, ParaPositronium, OrthoPositronium };
enum class SimulationOutcome { ReachedCutoff, ObservationLimit, NumericalFailure };
// SimulationOutcome intentionally remains the coarse contract used by the
// statistical collapse machinery.  This companion records WHY the shared
// ObservationLimit outcome occurred, which matters to a visual observation:
// a fixed window and an explicit user stop are both administrative censoring,
// but they should not be presented as the same action.
enum class SimulationStopReason {
    ReachedCutoff,
    ObservationTimeLimit,
    StopRequested,
    NumericalFailure
};
enum class SimulationObservationDisposition {
    ObservedEndpoint,
    AdministrativelyCensored,
    NumericalFailure
};

constexpr SimulationObservationDisposition observationDisposition(
        SimulationStopReason reason) noexcept {
    switch(reason) {
        case SimulationStopReason::ReachedCutoff:
            return SimulationObservationDisposition::ObservedEndpoint;
        case SimulationStopReason::ObservationTimeLimit:
        case SimulationStopReason::StopRequested:
            return SimulationObservationDisposition::AdministrativelyCensored;
        case SimulationStopReason::NumericalFailure:
            return SimulationObservationDisposition::NumericalFailure;
    }
    return SimulationObservationDisposition::NumericalFailure;
}
enum class VisualStyle { Unselected, Line, Dot };

struct InitialConditions {
    double relativeEnergy;
    double orbitalAngularMomentum;
    double predictedClosestApproach;
    double dipoleAlignment;
    double timeToCutoff;
    Phenomenon phenomenon;
    std::uint64_t seed;
};

struct SimulationResult {
    std::vector<Frame> frames;
    InitialConditions initial;
    SimulationOutcome outcome;
    SimulationStopReason stopReason=SimulationStopReason::NumericalFailure;
    double minimumSeparation;
    double elapsedTime;
    double finalRadiatedEnergy;
    double maximumBeta;
};

struct SimulationOptions {
    bool collectFrames = true;
    int frameCount = 1200;
    double observationTime = 0.0; // zero selects the phenomenon's visual window
    // Positive values override the visual point-particle cutoff. Statistical
    // CREM lifetime studies (experiments 1/2) instead pass comptonBarrierRadius
    // directly as runMechanicalTrajectory's trajectoryCutoff argument, so this
    // field is not what bounds them -- see crem_collapse.hpp.
    double terminalSeparation = 0.0;
    std::function<void(const Frame&)> frameReady;
    // Called after every accepted integration step.  Unlike frameReady this
    // is independent of physical-time sampling and is therefore suitable for
    // keeping an interactive visualization responsive.
    std::function<void(const State&)> stepReady;
    std::function<bool()> stopRequested;
    // Set false by callers that never read radiatedEnergy / radiatedMomentum /
    // radiatedAngularMomentum / boundField* from the result.  It switches off
    // the far-zone Poynting quadrature, which is the dominant per-step cost
    // and does not influence the trajectory (see
    // ClassicalTrajectoryEngine::Accuracy::computeOutwardFlux).
    bool radiatedEnergyBookkeeping = true;
};
