#pragma once

// Small text-formatting helpers used across the reporting/plotting side of
// the codebase: a fixed-precision table value, a picosecond time label that
// switches to scientific notation for very small values, and a dipole-spin
// label (up/down arrows) for a rendered Frame.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes; this is
// the first of several themed slices from a block that is, as a whole,
// presentation/statistics rather than engine).  Textually included at the
// same point inside positronium.cpp's shared anonymous namespace, itself
// already inside an #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region that
// positronium.cpp opens well before this #include and does not close until
// much later -- so this header needs no guard of its own; it inherits that
// one.  Depends on that namespace already having in scope: Frame.

std::string formatTableValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

std::string cutoffTimeLabel(double timeToCutoff) {
    if (std::isinf(timeToCutoff)) return "not reached";
    const double picoseconds = timeToCutoff * 1.0e12;
    std::ostringstream out;
    if (picoseconds < 1.0e-3) out << std::scientific << std::setprecision(2);
    else out << std::fixed << std::setprecision(3);
    out << picoseconds << " ps";
    return out.str();
}

std::string spinLabel(const Frame& frame) {
    const char* firstArrow = frame.firstDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    const char* secondArrow = frame.secondDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    return std::string(firstArrow) + " " + secondArrow;
}
