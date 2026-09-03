#pragma once

// Small text-formatting helpers used across the reporting/plotting side of
// the codebase: a fixed-precision table value, a picosecond time label that
// switches to scientific notation for very small values, and a dipole-spin
// label (up/down arrows) for a rendered Frame.
//
// Self-contained and order-independent.  Frame arrives with
// simulation_interface.hpp rather than from the surrounding namespace.
// Note this header is included inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region of positronium.cpp: the
// #pragma once above is its own guard, unrelated to that one.

#include "simulation_interface.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

inline std::string formatTableValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

inline std::string cutoffTimeLabel(double timeToCutoff) {
    if (std::isinf(timeToCutoff)) return "not reached";
    const double picoseconds = timeToCutoff * 1.0e12;
    std::ostringstream out;
    if (picoseconds < 1.0e-3) out << std::scientific << std::setprecision(2);
    else out << std::fixed << std::setprecision(3);
    out << picoseconds << " ps";
    return out.str();
}

inline std::string spinLabel(const Frame& frame) {
    const char* firstArrow = frame.firstDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    const char* secondArrow = frame.secondDipole.z >= 0.0 ? "#uparrow" : "#downarrow";
    return std::string(firstArrow) + " " + secondArrow;
}
