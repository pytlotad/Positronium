#include <TApplication.h>
#include <TButton.h>
#include <TCanvas.h>
#include <TInterpreter.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TView.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

// These functions are intentionally global: ROOT's TButton invokes its action
// through the interpreter while the animation loop observes these flags.
bool gSimulationPaused = false;
bool gExitRequested = false;
TButton* gStopButton = nullptr;

void ToggleSimulation() {
    gSimulationPaused = !gSimulationPaused;
    // Update immediately on the click, not only on the next animation frame.
    if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
}
void ExitSimulation() {
    gExitRequested = true;
    if (gApplication) gApplication->Terminate(0);
}

// Classical, non-relativistic two-body electrodynamics in SI units.
// The orbital energy lost through Larmor radiation is removed from relative
// motion via a stable, order-reduced radiation-reaction approximation.
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double epsilon0 = 8.8541878128e-12;
constexpr double c = 299792458.0;
constexpr double mu0 = 4.0 * pi * 1.0e-7;
constexpr double eCharge = 1.602176634e-19;
constexpr double electronMass = 9.1093837139e-31;
constexpr double positronMass = electronMass;
constexpr double bohrMagneton = 9.2740100657e-24;
constexpr double nuclearMagneton = 5.0507837461e-27;
constexpr double coulomb = 1.0 / (4.0 * pi * epsilon0);
constexpr double bohrRadius = 5.29177210903e-11;
constexpr double nuclearCutoff = 1.0e-14; // point-particle theory has failed below this scale

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return {x + v.x, y + v.y, z + v.z}; }
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    double squaredNorm() const { return x*x + y*y + z*z; }
    double norm() const { return std::sqrt(squaredNorm()); }
};

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
Vec3 unit(const Vec3& v) { return v / v.norm(); }

struct ElectromagneticField { Vec3 electric, magnetic; };

struct State {
    Vec3 electronPosition, positronPosition;
    Vec3 electronVelocity, positronVelocity;
    Vec3 electronAcceleration, positronAcceleration;
    Vec3 electronDipole, positronDipole; // classical magnetic moments, in J/T
    double time = 0;
    double radiatedEnergy = 0;
};

struct Frame {
    Vec3 electron, positron, electronDipole, positronDipole;
    double time, radius, radiatedEnergy, mechanicalEnergy;
    double electronMechanicalEnergy, positronMechanicalEnergy;
};

double separation(const State& s) { return (s.electronPosition - s.positronPosition).norm(); }

double dot(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }

// Mutual, retarded Lienard-Wiechert field of a moving point charge.  The
// source at t-R/c is reconstructed locally from its velocity and acceleration.
ElectromagneticField lienardWiechertField(const Vec3& observationPosition,
                                          const Vec3& sourcePosition,
                                          const Vec3& sourceVelocity,
                                          const Vec3& sourceAcceleration,
                                          double sourceCharge) {
    double delay = (observationPosition - sourcePosition).norm() / c;
    Vec3 retardedPosition, retardedVelocity;
    for (int iteration = 0; iteration < 3; ++iteration) {
        retardedPosition = sourcePosition - sourceVelocity * delay + sourceAcceleration * (0.5 * delay * delay);
        retardedVelocity = sourceVelocity - sourceAcceleration * delay;
        delay = (observationPosition - retardedPosition).norm() / c;
    }
    const Vec3 displacement = observationPosition - retardedPosition;
    const double distance = displacement.norm();
    const Vec3 direction = displacement / distance;
    const Vec3 beta = retardedVelocity / c;
    const double betaSquared = beta.squaredNorm();
    const double kappa = std::max(1.0e-8, 1.0 - dot(direction, beta));
    const Vec3 velocityField = (direction - beta) * ((1.0 - betaSquared) /
                              (kappa*kappa*kappa * distance*distance));
    const Vec3 accelerationField = cross(direction, cross(direction - beta, sourceAcceleration)) /
                                   (c * kappa*kappa*kappa * distance);
    const Vec3 electric = (velocityField + accelerationField) * (coulomb * sourceCharge);
    return {electric, cross(direction, electric) / c};
}

// Static magnetic field of a point dipole.  'at' is the field point.
Vec3 dipoleField(const Vec3& at, const Vec3& sourcePosition, const Vec3& sourceDipole) {
    const Vec3 r = at - sourcePosition;
    const double rLength = r.norm();
    const Vec3 n = r / rLength;
    return (n * (3.0 * sourceDipole.x*n.x + 3.0 * sourceDipole.y*n.y + 3.0 * sourceDipole.z*n.z) - sourceDipole)
           * (mu0 / (4.0 * pi * rLength*rLength*rLength));
}

// F = grad(m_target . B_source), with r directed from source to target.
Vec3 dipoleForce(const Vec3& targetPosition, const Vec3& sourcePosition,
                 const Vec3& targetDipole, const Vec3& sourceDipole) {
    const Vec3 r = targetPosition - sourcePosition;
    const double rLength = r.norm();
    const Vec3 n = r / rLength;
    const double mnTarget = targetDipole.x*n.x + targetDipole.y*n.y + targetDipole.z*n.z;
    const double mnSource = sourceDipole.x*n.x + sourceDipole.y*n.y + sourceDipole.z*n.z;
    const double dots = targetDipole.x*sourceDipole.x + targetDipole.y*sourceDipole.y + targetDipole.z*sourceDipole.z;
    return (sourceDipole * mnTarget + targetDipole * mnSource + n * (dots - 5.0*mnTarget*mnSource))
           * (3.0 * mu0 / (4.0 * pi * rLength*rLength*rLength*rLength));
}

void precessDipole(Vec3& dipole, const Vec3& field, double gyromagneticRatio, double dt) {
    const double fieldMagnitude = field.norm();
    if (fieldMagnitude == 0.0) return;
    const Vec3 axis = field / fieldMagnitude;
    const double angle = gyromagneticRatio * fieldMagnitude * dt;
    // Rodrigues rotation preserves the magnitude of the magnetic moment.
    const Vec3 rotated = dipole * std::cos(angle) + cross(axis, dipole) * std::sin(angle)
                       + axis * ((axis.x*dipole.x + axis.y*dipole.y + axis.z*dipole.z) * (1.0 - std::cos(angle)));
    dipole = rotated * (dipole.norm() / rotated.norm());
}

Vec3 coulombAcceleration(const Vec3& targetPosition, const Vec3& sourcePosition, double coefficient) {
    const Vec3 displacement = targetPosition - sourcePosition;
    const double distanceSquared = displacement.squaredNorm();
    if (distanceSquared <= 0.0) return {};
    const double distanceCubed = distanceSquared * std::sqrt(distanceSquared);
    return displacement * (coefficient * coulomb / distanceCubed);
}

// Mechanical energy of both particles: translational kinetic energy, Coulomb
// potential and magnetic dipole-dipole potential.  It deliberately excludes
// energy that has already escaped as electromagnetic radiation.
double mechanicalEnergy(const State& s) {
    const double r = separation(s);
    const double kinetic = 0.5 * electronMass * s.electronVelocity.squaredNorm()
                         + 0.5 * positronMass * s.positronVelocity.squaredNorm();
    const double coulombPotential = -coulomb * eCharge*eCharge / r;
    const double dipolePotential = -dot(s.electronDipole,
                                        dipoleField(s.electronPosition, s.positronPosition, s.positronDipole));
    return kinetic + coulombPotential + dipolePotential;
}

double electronMechanicalEnergy(const State& s) {
    const double r = separation(s);
    const double kinetic = 0.5 * electronMass * s.electronVelocity.squaredNorm();
    const double coulombShare = -0.5 * coulomb * eCharge*eCharge / r;
    const double dipoleShare = -0.5 * dot(s.electronDipole,
                                          dipoleField(s.electronPosition, s.positronPosition, s.positronDipole));
    return kinetic + coulombShare + dipoleShare;
}

double positronMechanicalEnergy(const State& s) {
    const double r = separation(s);
    const double kinetic = 0.5 * positronMass * s.positronVelocity.squaredNorm();
    const double coulombShare = -0.5 * coulomb * eCharge*eCharge / r;
    const double dipoleShare = -0.5 * dot(s.positronDipole,
                                          dipoleField(s.positronPosition, s.electronPosition, s.electronDipole));
    return kinetic + coulombShare + dipoleShare;
}

// Coulomb motion plus the force and torque of two classical magnetic dipoles.
void advance(State& s, double dt) {
    Vec3 ae = coulombAcceleration(s.electronPosition, s.positronPosition, -eCharge*eCharge/electronMass);
    Vec3 ap = coulombAcceleration(s.positronPosition, s.electronPosition, -eCharge*eCharge/positronMass);
    ae += dipoleForce(s.electronPosition, s.positronPosition, s.electronDipole, s.positronDipole) / electronMass;
    ap += dipoleForce(s.positronPosition, s.electronPosition, s.positronDipole, s.electronDipole) / positronMass;
    s.electronVelocity += ae * (0.5 * dt);
    s.positronVelocity += ap * (0.5 * dt);
    s.electronPosition += s.electronVelocity * dt;
    s.positronPosition += s.positronVelocity * dt;
    ae = coulombAcceleration(s.electronPosition, s.positronPosition, -eCharge*eCharge/electronMass);
    ap = coulombAcceleration(s.positronPosition, s.electronPosition, -eCharge*eCharge/positronMass);
    ae += dipoleForce(s.electronPosition, s.positronPosition, s.electronDipole, s.positronDipole) / electronMass;
    ap += dipoleForce(s.positronPosition, s.electronPosition, s.positronDipole, s.electronDipole) / positronMass;
    s.electronVelocity += ae * (0.5 * dt);
    s.positronVelocity += ap * (0.5 * dt);

    const Vec3 fieldAtElectron = dipoleField(s.electronPosition, s.positronPosition, s.positronDipole);
    const Vec3 fieldAtPositron = dipoleField(s.positronPosition, s.electronPosition, s.electronDipole);
    // The real dipole precession is extremely slow on an orbital timescale.
    // This factor advances only the displayed spin dynamics so its changing
    // 3D direction is observable; the mechanical dipole force above remains
    // at its physical, unamplified value.
    constexpr double visibleSpinTimeScale = 3.0e4;
    precessDipole(s.electronDipole, fieldAtElectron, -1.76085963023e11, dt * visibleSpinTimeScale);
    precessDipole(s.positronDipole, fieldAtPositron, 1.76085963023e11, dt * visibleSpinTimeScale);

    // Larmor power of both accelerated charges: P = q^2 a^2/(6 pi eps0 c^3).
    // The corresponding energy is removed from relative motion only, preserving
    // the centre-of-mass momentum of the isolated atom.
    const double radiatedPower = eCharge*eCharge * (ae.squaredNorm() + ap.squaredNorm()) /
                                 (6.0 * pi * epsilon0 * c*c*c);
    const Vec3 relativeVelocity = s.electronVelocity - s.positronVelocity;
    const double reducedMass = electronMass * positronMass / (electronMass + positronMass);
    const double relativeKineticEnergy = 0.5 * reducedMass * relativeVelocity.squaredNorm();
    const double removedEnergy = std::min(radiatedPower * dt, 0.02 * relativeKineticEnergy);
    if (relativeKineticEnergy > 0.0 && removedEnergy > 0.0) {
        const double scale = std::sqrt(1.0 - removedEnergy / relativeKineticEnergy);
        const Vec3 centreVelocity = (s.electronVelocity * electronMass + s.positronVelocity * positronMass) /
                                    (electronMass + positronMass);
        const Vec3 newRelativeVelocity = relativeVelocity * scale;
        s.electronVelocity = centreVelocity + newRelativeVelocity * (positronMass / (electronMass + positronMass));
        s.positronVelocity = centreVelocity - newRelativeVelocity * (electronMass / (electronMass + positronMass));
        s.radiatedEnergy += removedEnergy;
    }

    s.time += dt;
}

std::vector<Frame> simulate() {
    const double reducedMass = electronMass * positronMass / (electronMass + positronMass);
    const double circularSpeed = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * bohrRadius));
    State s;
    // Centre of mass is at rest; both bodies therefore move, with the positron's
    // small reflex orbit included.
    s.electronPosition = {bohrRadius * positronMass / (electronMass + positronMass), 0, 0};
    s.positronPosition = {-bohrRadius * electronMass / (electronMass + positronMass), 0, 0};
    s.electronVelocity = {0, circularSpeed * positronMass / (electronMass + positronMass), 0};
    s.positronVelocity = {0, -circularSpeed * electronMass / (electronMass + positronMass), 0};
    // Fixed seed means the initial, randomly distributed directions are
    // repeatable between runs, which makes the simulation comparable.
    std::mt19937_64 random(0x484944524f47454eULL);
    std::uniform_real_distribution<double> distributed(-1.0, 1.0);
    std::uniform_real_distribution<double> azimuth(0.0, 2.0*pi);
    const auto randomDirection = [&]() {
        const double z = distributed(random);
        const double phi = azimuth(random);
        const double radial = std::sqrt(1.0 - z*z);
        return Vec3{radial*std::cos(phi), radial*std::sin(phi), z};
    };
    s.electronDipole = randomDirection() * bohrMagneton;
    s.positronDipole = randomDirection() * bohrMagneton;

    constexpr int frameCount = 15000;
    constexpr double displayedLifetime = 1.50e-10;
    std::vector<Frame> frames;
    frames.reserve(frameCount);
    double nextFrame = 0.0;
    const double frameInterval = displayedLifetime / (frameCount - 1);

    while (frames.size() < frameCount && separation(s) > nuclearCutoff) {
        if (s.time >= nextFrame) {
            frames.push_back({s.electronPosition, s.positronPosition, s.electronDipole, s.positronDipole,
                              s.time, separation(s), s.radiatedEnergy, mechanicalEnergy(s),
                              electronMechanicalEnergy(s), positronMechanicalEnergy(s)});
            nextFrame += frameInterval;
        }
        // At least 80 steps per instantaneous orbit.  This shortens with the
        // radius and keeps the Coulomb integration stable through the plunge.
        const double r = separation(s);
        const double omega = std::sqrt(coulomb * eCharge*eCharge / (reducedMass * r*r*r));
        advance(s, std::min(2.0e-18, 2.0 * pi / (80.0 * omega)));
    }
    if (frames.empty()) throw std::runtime_error("No simulation frames were produced");
    return frames;
}

std::string labelFor(const Frame& f) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << "t = " << f.time * 1.0e12 << " ps"
        << "     E_{e+p} = " << f.mechanicalEnergy / eCharge << " eV"
        << "     E_{rad} = " << f.radiatedEnergy / eCharge << " eV";
    return out.str();
}

std::string deltaLabelFor(const Frame& current, const Frame& initial) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << "#Deltat = " << (current.time - initial.time) * 1.0e12 << " ps"
        << "     #DeltaE_{e+p} = " << (current.mechanicalEnergy - initial.mechanicalEnergy) / eCharge << " eV"
        << "     #DeltaE_{rad} = " << (current.radiatedEnergy - initial.radiatedEnergy) / eCharge << " eV";
    return out.str();
}

std::string formatTableValue(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2) << value;
    return out.str();
}

void setDipoleArrow(TPolyLine3D& shaft, TPolyLine3D& leftHead, TPolyLine3D& rightHead,
                    const Vec3& centre, const Vec3& dipole) {
    const Vec3 direction = unit(dipole);
    const Vec3 tip = centre + direction * 0.31;
    const Vec3 helper = std::abs(direction.z) < 0.8 ? Vec3{0, 0, 1} : Vec3{0, 1, 0};
    const Vec3 side = unit(cross(direction, helper));
    const Vec3 base = tip - direction * 0.085;
    shaft.SetPoint(0, centre.x, centre.y, centre.z);
    shaft.SetPoint(1, tip.x, tip.y, tip.z);
    leftHead.SetPoint(0, tip.x, tip.y, tip.z);
    leftHead.SetPoint(1, base.x + side.x*0.055, base.y + side.y*0.055, base.z + side.z*0.055);
    rightHead.SetPoint(0, tip.x, tip.y, tip.z);
    rightHead.SetPoint(1, base.x - side.x*0.055, base.y - side.y*0.055, base.z - side.z*0.055);
}
} // namespace

int main(int argc, char** argv) {
    std::vector<Frame> frames = simulate();
    std::cout << "Classical radiative collapse simulated for " << frames.back().time << " s; "
              << frames.size() << " animation frames.\n";

    TApplication app("classical_hydrogen", &argc, argv);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(kBlack);
    gStyle->SetPadColor(kBlack);

    TCanvas canvas("atom", "Classical model of positronium", 1100, 820);
    canvas.SetFillColor(kBlack);
    canvas.SetSupportGL(kTRUE);
    TPad scene("scene", "Simulation", 0.0, 0.0, 1.0, 0.86);
    TPad controls("controls", "Controls", 0.0, 0.86, 1.0, 1.0);
    scene.SetFillColor(kBlack);
    controls.SetFillColor(kBlack);
    scene.Draw();
    controls.Draw();

    scene.cd();
    // ROOT rotates this 3D view when the user drags the mouse in the scene.
    // A non-zero Z span preserves perspective for the initially planar orbit.
    TView* view = TView::CreateView(1);
    view->SetRange(-1.10, -1.10, -0.50, 1.10, 1.10, 0.50);
    gPad->SetTheta(70);
    gPad->SetPhi(25);

    const double scale = 1.0 / bohrRadius;
    std::vector<double> x(frames.size()), y(frames.size()), z(frames.size());
    for (size_t i = 0; i < frames.size(); ++i) {
        x[i] = frames[i].electron.x * scale;
        y[i] = frames[i].electron.y * scale;
        z[i] = frames[i].electron.z * scale;
    }
    TPolyLine3D path(static_cast<int>(frames.size()), x.data(), y.data(), z.data());
    path.SetLineColor(kOrange + 7);
    path.SetLineWidth(2);
    path.Draw();

    TPolyMarker3D electron(1), positron(1);
    electron.SetMarkerStyle(20); electron.SetMarkerSize(1.8); electron.SetMarkerColor(kAzure + 1);
    positron.SetMarkerStyle(20); positron.SetMarkerSize(2.5); positron.SetMarkerColor(kRed + 1);
    electron.Draw("same");
    positron.Draw("same");

    TPolyLine3D electronDipoleShaft(2), electronDipoleLeft(2), electronDipoleRight(2);
    TPolyLine3D positronDipoleShaft(2), positronDipoleLeft(2), positronDipoleRight(2);
    for (TPolyLine3D* arrow : {&electronDipoleShaft, &electronDipoleLeft, &electronDipoleRight}) {
        arrow->SetLineColor(kAzure + 1); arrow->SetLineWidth(3); arrow->Draw("same");
    }
    for (TPolyLine3D* arrow : {&positronDipoleShaft, &positronDipoleLeft, &positronDipoleRight}) {
        arrow->SetLineColor(kRed + 1); arrow->SetLineWidth(3); arrow->Draw("same");
    }

    const Frame& initialFrame = frames.front();
    constexpr double bottomHeaderY = 0.125;
    constexpr double bottomInitialY = 0.090;
    constexpr double bottomCurrentY = 0.055;
    constexpr double bottomDeltaY = 0.020;
    const std::array<double, 5> bottomXs = {0.05, 0.21, 0.34, 0.47, 0.61};
    const std::array<const char*, 6> bottomHeaders = {
        "stan", "t [ps]", "E_{e} [eV]", "E_{p} [eV]", "E_{sum} [eV]", "E_{rad} [eV]"
    };
    std::array<TLatex, 5> bottomHeaderLabels;
    for (size_t i = 0; i < bottomHeaderLabels.size(); ++i) {
        bottomHeaderLabels[i].SetNDC(); bottomHeaderLabels[i].SetTextColor(kWhite);
        bottomHeaderLabels[i].SetTextSize(0.024); bottomHeaderLabels[i].SetTextFont(62);
        bottomHeaderLabels[i].SetText(bottomXs[i], bottomHeaderY, bottomHeaders[i]);
        bottomHeaderLabels[i].Draw();
    }
    TLatex bottomRadHeader;
    bottomRadHeader.SetNDC(); bottomRadHeader.SetTextColor(kWhite);
    bottomRadHeader.SetTextSize(0.024); bottomRadHeader.SetTextFont(62);
    bottomRadHeader.SetText(0.75, bottomHeaderY, bottomHeaders.back());
    bottomRadHeader.Draw();
    std::array<TLatex, 5> initialBottomRow;
    std::array<TLatex, 5> currentBottomRow;
    std::array<TLatex, 5> deltaBottomRow;
    TLatex initialBottomRad;
    TLatex currentBottomRad;
    TLatex deltaBottomRad;
    const auto initializeBottomRow = [&](std::array<TLatex, 5>& row, int color, double textSize) {
        for (TLatex& cell : row) {
            cell.SetNDC(); cell.SetTextColor(color);
            cell.SetTextSize(textSize); cell.SetTextFont(62);
            cell.SetText(0, 0, "");
            cell.Draw();
        }
    };
    initializeBottomRow(initialBottomRow, kCyan + 1, 0.024);
    initializeBottomRow(currentBottomRow, kYellow + 1, 0.024);
    initializeBottomRow(deltaBottomRow, kGreen + 2, 0.024);
    initialBottomRad.SetNDC(); initialBottomRad.SetTextColor(kCyan + 1);
    initialBottomRad.SetTextSize(0.024); initialBottomRad.SetTextFont(62);
    currentBottomRad.SetNDC(); currentBottomRad.SetTextColor(kYellow + 1);
    currentBottomRad.SetTextSize(0.024); currentBottomRad.SetTextFont(62);
    deltaBottomRad.SetNDC(); deltaBottomRad.SetTextColor(kGreen + 2);
    deltaBottomRad.SetTextSize(0.024); deltaBottomRad.SetTextFont(62);

    const auto drawBottomRow = [&](std::array<TLatex, 5>& row, double y, const std::array<std::string, 5>& values) {
        for (size_t i = 0; i < row.size(); ++i) {
            row[i].SetText(bottomXs[i], y, values[i].c_str());
            row[i].Draw();
        }
    };

    drawBottomRow(initialBottomRow, bottomInitialY, std::array<std::string, 5>{
        "initial",
        formatTableValue(initialFrame.time * 1.0e12),
        formatTableValue(initialFrame.electronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.positronMechanicalEnergy / eCharge),
        formatTableValue(initialFrame.mechanicalEnergy / eCharge)
    });
    initialBottomRad.SetText(0.75, bottomInitialY, formatTableValue(initialFrame.radiatedEnergy / eCharge).c_str());
    initialBottomRad.Draw();
    drawBottomRow(currentBottomRow, bottomCurrentY, std::array<std::string, 5>{"current", "0.00", "0.00", "0.00", "0.00"});
    currentBottomRad.SetText(0.75, bottomCurrentY, "0.00");
    currentBottomRad.Draw();
    drawBottomRow(deltaBottomRow, bottomDeltaY, std::array<std::string, 5>{"delta", "0.00", "0.00", "0.00", "0.00"});
    deltaBottomRad.SetText(0.75, bottomDeltaY, "0.00");
    deltaBottomRad.Draw();

    controls.cd();
    gInterpreter->Declare("void ToggleSimulation(); void ExitSimulation();");
    TButton stopButton("STOP", "ToggleSimulation();", 0.73, 0.25, 0.85, 0.78);
    stopButton.SetFillColor(kOrange + 7);
    stopButton.SetTextFont(62);
    stopButton.Draw();
    gStopButton = &stopButton;
    TButton exitButton("EXIT", "ExitSimulation();", 0.87, 0.25, 0.97, 0.78);
    exitButton.SetFillColor(kRed + 1);
    exitButton.SetTextFont(62);
    exitButton.Draw();
    const auto syncStopButton = [&]() {
        if (gStopButton) gStopButton->SetTitle(gSimulationPaused ? "START" : "STOP");
        controls.Modified();
    };
    syncStopButton();
    controls.Modified();
    controls.Update();

    const auto updateBottomRow = [&](const Frame& f) {
        const std::array<std::string, 5> currentValues = {
            "current",
            formatTableValue(f.time * 1.0e12),
            formatTableValue(f.electronMechanicalEnergy / eCharge),
            formatTableValue(f.positronMechanicalEnergy / eCharge),
            formatTableValue(f.mechanicalEnergy / eCharge)
        };
        drawBottomRow(currentBottomRow, bottomCurrentY, currentValues);
        currentBottomRad.SetText(0.75, bottomCurrentY, formatTableValue(f.radiatedEnergy / eCharge).c_str());
        currentBottomRad.Draw();

        const std::array<std::string, 5> deltaValues = {
            "delta",
            formatTableValue((f.time - initialFrame.time) * 1.0e12),
            formatTableValue((f.electronMechanicalEnergy - initialFrame.electronMechanicalEnergy) / eCharge),
            formatTableValue((f.positronMechanicalEnergy - initialFrame.positronMechanicalEnergy) / eCharge),
            formatTableValue((f.mechanicalEnergy - initialFrame.mechanicalEnergy) / eCharge)
        };
        drawBottomRow(deltaBottomRow, bottomDeltaY, deltaValues);
        deltaBottomRad.SetText(0.75, bottomDeltaY,
                               formatTableValue((f.radiatedEnergy - initialFrame.radiatedEnergy) / eCharge).c_str());
        deltaBottomRad.Draw();
    };

    updateBottomRow(initialFrame);

    for (size_t i = 0; i < frames.size(); ++i) {
        while (gSimulationPaused && !gExitRequested) {
            syncStopButton();
            controls.Update();
            gSystem->ProcessEvents();
            syncStopButton();
            gSystem->Sleep(20);
        }
        if (gExitRequested) return 0;
        syncStopButton();
        const Frame& f = frames[i];
        electron.SetPoint(0, f.electron.x * scale, f.electron.y * scale, f.electron.z * scale);
        positron.SetPoint(0, f.positron.x * scale, f.positron.y * scale, f.positron.z * scale);
        const Vec3 electronPosition = f.electron * scale;
        const Vec3 positronPosition = f.positron * scale;
        setDipoleArrow(electronDipoleShaft, electronDipoleLeft, electronDipoleRight,
                       electronPosition, f.electronDipole);
        setDipoleArrow(positronDipoleShaft, positronDipoleLeft, positronDipoleRight,
                       positronPosition, f.positronDipole);
        updateBottomRow(f);
        if (i % 4 == 0 || i + 1 == frames.size()) {
            scene.Modified();
            canvas.Modified();
            canvas.Update();
            gSystem->ProcessEvents();
            syncStopButton();
            gSystem->Sleep(8);
        }
    }
    controls.cd();
    canvas.Modified();
    canvas.Update();
    app.Run();
    return 0;
}
