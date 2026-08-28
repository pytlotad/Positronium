#pragma once

// CREM collapse estimator: the orbit-averaged secular integration that turns a
// prepared bound state into a measured classical inspiral time, together with
// the closed-form electrodynamic references it is compared against.
//
// Textual module, like the other headers here: it is included once, from
// inside the anonymous namespace of positronium.cpp and inside the production
// #ifndef, so everything it needs (State, the trajectory engine, simulate(),
// splitMix64) is already in scope.  It deliberately contains no ROOT: nothing
// in this file draws anything, and the panels that display these results live
// in positronium.cpp.

// One stochasticElectricDipole photon (modules/electrodynamics.hpp), recorded
// as it would be measured by a fixed, distant lab observer rather than in the
// pair's own instantaneous rest frame S' the emission physics is computed in.
// Diagnostic/plotting only -- see estimateCremCollapse's own comment at the
// point these are filled for why the internal energy/momentum/angular-
// momentum bookkeeping needs no boost (it is already self-consistent in S')
// and this record exists purely to answer the separate, external question
// "what spectrum/angles would a lab detector see", per README point L.
struct LabFramePhoton {
    double energyJoules=0.0;
    // Angle between the photon's LAB-FRAME direction and the recoil axis
    // (direction of centreOfMassVelocity at the instant of emission), radians
    // in [0,pi] -- the only frame-connecting axis available (there is no
    // fixed external reference direction in this model otherwise).  NaN when
    // the recoil speed at emission was too small for that axis to be
    // meaningful (beta below 1e-9): direction aberration is undefined at
    // exactly zero recoil, not zero itself.
    double angleFromRecoilAxisRadians=std::numeric_limits<double>::quiet_NaN();
    // |centreOfMassVelocity|/c at the instant of emission (pre-kick), kept
    // for weighting/diagnostics -- e.g. restricting a histogram to the
    // deep-cascade photons where the lab/source-frame difference is largest.
    double sourceBeta=0.0;
};

struct CremCollapseEstimate {
    double lifetimeSeconds=std::numeric_limits<double>::quiet_NaN();
    double calibrationSeconds=0.0;
    double meanRadiatedPowerWatts=std::numeric_limits<double>::quiet_NaN();
    // Lab-frame counterparts of the two fields above.  lifetimeSeconds/
    // meanRadiatedPowerWatts are accumulated in S', the pair's own
    // instantaneous rest frame -- correct and self-consistent for internal
    // bookkeeping (see LabFramePhoton's own comment), but centreOfMassVelocity
    // genuinely grows over a stochastic trajectory's photon cascade, and a
    // fixed lab observer's clock runs slow relative to S' by gamma(beta)
    // while that recoil speed is nonzero.  These two fields integrate that
    // dilation through however many discrete recoil kicks occurred, so that
    // "what a lab clock/power-meter would read" (README point N) does not
    // have to be re-derived by consumers of this struct from beta alone --
    // NaN/0 for the continuous (non-stochastic) radiation-reaction models
    // and the mechanical trajectory path, where centreOfMassVelocity never
    // moves and these are identical to the S'-frame fields anyway.
    double lifetimeSecondsLab=std::numeric_limits<double>::quiet_NaN();
    double meanRadiatedPowerWattsLab=std::numeric_limits<double>::quiet_NaN();
    // Lab-frame counterpart of calibrationSeconds, tracked at every exit
    // point that one is (including the censored/failed ones): the
    // right-censored survival sample needs its "still bound as of" bound in
    // the SAME frame convention as the completed-collapse observations
    // (lifetimeSecondsLab) it is combined with, or the Kaplan-Meier estimate
    // would silently mix frames.
    double calibrationSecondsLab=0.0;
    SimulationOutcome calibrationOutcome=SimulationOutcome::NumericalFailure;
    // Every stochasticElectricDipole photon this trajectory fired, converted
    // to lab-frame observables.  Empty for continuous (non-stochastic)
    // radiation-reaction models, and for the mechanical trajectory path.
    std::vector<LabFramePhoton> labFramePhotons;
    // Distinguishes the two ways a trajectory can end as ObservationLimit.
    // Without it, "the reaction force is switched off, so this orbit will
    // never decay" is indistinguishable from "this orbit is decaying but ran
    // out of wall clock", and the report tells the user to raise a budget
    // that cannot possibly help.
    bool secularLossAbsent=false;
    // Kepler period of the quasi-closed orbit.  The inspiral is not a closed
    // orbit: T shrinks as a^(3/2) while the pair sinks, so a single number
    // cannot describe it.  Both ends of the run are reported instead, plus
    // the revolution count that the orbit-averaged integrator accumulates
    // between them (measured orbits and analytically skipped ones alike).
    double initialPeriodSeconds=std::numeric_limits<double>::quiet_NaN();
    double finalPeriodSeconds=std::numeric_limits<double>::quiet_NaN();
    double revolutions=std::numeric_limits<double>::quiet_NaN();
    // Osculating elements of the prepared orbit, and the closed-form classical
    // prediction they imply.  These make the run testable against textbook
    // electrodynamics rather than only against itself.
    double initialSemiMajorAxis=std::numeric_limits<double>::quiet_NaN();
    double initialEccentricity=std::numeric_limits<double>::quiet_NaN();
    double analyticCollapseSeconds=std::numeric_limits<double>::quiet_NaN();
    // Mean over checkpoints of (measured orbital energy loss rate) / (Larmor
    // rate for the same osculating orbit).  1 means the engine reproduces
    // coherent electric-dipole radiation exactly.
    double larmorPowerRatio=std::numeric_limits<double>::quiet_NaN();
    // Classical dipole-dipole interaction energy of the prepared pair,
    // expressed as a frequency so it can sit beside the measured o-Ps/p-Ps
    // hyperfine splitting.
    double dipoleCouplingHz=std::numeric_limits<double>::quiet_NaN();
    // ANNIHILATION TIED TO THE TERMINAL RADIUS.  The project's own
    // annihilation generator is a quantum prescription deliberately
    // independent of the trajectory: it assumes a pair AT REST and puts the
    // full 2 m_e c^2 into the photons, so its 2-gamma line is a compile-time
    // constant.  These fields are the other thing one can ask -- what the
    // pair this model actually integrated has left when it stops.
    //
    // The classical orbit is bound, so its invariant energy is BELOW the
    // rest-mass sum by the binding it has accumulated:
    //
    //     W = (m1+m2) c^2 + mu * epsilon,   epsilon < 0,
    //
    // and W, not 2 m_e c^2, is what the final-state photons have to share.
    // Filled only where a trajectory genuinely reaches the terminal radius;
    // NaN/empty for censored and failed runs, which have no terminal state
    // to annihilate from.
    double terminalSemiMajorAxis=std::numeric_limits<double>::quiet_NaN();
    double terminalBindingEnergy=std::numeric_limits<double>::quiet_NaN();
    // Dipole-dipole interaction energy at the terminal configuration, the
    // one term that distinguishes para from ortho in the final-state
    // invariant.  Signed: negative for the orientation that binds.
    double terminalDipoleEnergy=std::numeric_limits<double>::quiet_NaN();
    double annihilationInvariantEnergy=std::numeric_limits<double>::quiet_NaN();
    // Two entries for para (2 gamma), three for ortho (3 gamma), in the
    // pair's own centre-of-momentum frame and summing to
    // annihilationInvariantEnergy exactly.
    std::vector<double> annihilationPhotonEnergies;
};

// Final-state photon energies for a pair annihilating with invariant energy
// W, in its own centre-of-momentum frame.
//
// Para (2 gamma): back to back, W/2 each, forced by momentum conservation
// alone.  At W = 2 m_e c^2 this is the textbook 511 keV line; tied to a
// bound terminal state it sits BELOW it, by exactly half the binding energy
// the orbit accumulated.
//
// Ortho (3 gamma): the Ore-Powell spectrum, whose only scale is the maximum
// single-photon energy W/2, so it rescales with W rather than being pinned
// to m_e c^2.  Sampled by the same rejection rule the reference curve uses
// (acceptance proportional to the normalized density below), then the three
// energies are closed onto W exactly -- momentum conservation for three
// massless quanta requires them to sum to W and to be constructible as a
// closed triangle, which x1+x2+x3 = 2 with each x <= 1 guarantees.
std::vector<double> annihilationPhotonEnergiesFor(
        double invariantEnergy,bool para,std::uint64_t& stream) {
    if(!(invariantEnergy>0.0)||!std::isfinite(invariantEnergy)) return {};
    if(para) return {0.5*invariantEnergy,0.5*invariantEnergy};
    // Ore-Powell density in x = E/(W/2), normalized to its own maximum.
    const auto density=[](double x) {
        if(!(x>0.0)||x>=1.0) return 0.0;
        const double u=1.0-x, d=2.0-x;
        return x*u/(d*d)-2.0*u*u*std::log(u)/(d*d*d)
              +d/x+2.0*u*std::log(u)/(x*x);
    };
    // The density is maximal at the endpoint; 1.0 there bounds it.
    const double bound=std::max(density(0.999),1.0);
    double x1=0.0,x2=0.0;
    for(int attempt=0;attempt<10000;++attempt) {
        x1=drawUniformUnit(stream);
        x2=drawUniformUnit(stream);
        const double x3=2.0-x1-x2;
        // Physical three-photon region: every energy positive and at most W/2.
        if(!(x1>0.0&&x2>0.0&&x3>0.0&&x1<=1.0&&x2<=1.0&&x3<=1.0)) continue;
        if(drawUniformUnit(stream)*bound<=density(x1)) break;
    }
    const double x3=2.0-x1-x2;
    const double half=0.5*invariantEnergy;
    if(!(x1>0.0&&x2>0.0&&x3>0.0)) {
        // Sampling never landed in the physical region: fall back on the
        // symmetric configuration rather than returning something that does
        // not sum to W.
        return {2.0*half/3.0,2.0*half/3.0,2.0*half/3.0};
    }
    return {x1*half,x2*half,x3*half};
}

struct OsculatingElements { double specificEnergy=0.0; double specificAngularMomentum=0.0; };

// --- Closed-form classical inspiral, used only as an external reference ---
//
// The pair carries electric dipole moment d = e*r, so it radiates at the
// Larmor rate P = |d''|^2/(6 pi eps0 c^3) with |d''| = 2 k e^3/(m r^2).
// Feeding that into dE/dt with E = -k e^2/(2a) gives da/dt = -C/a^2 with
//
//     C = 8 k e^4 / (6 pi eps0 c^3 m^2),
//
// so the time from a_i to a_f is (a_i^3 - a_f^3)/(3 C).  Nothing in the CREM
// engine is used to obtain this; it is the textbook answer the engine is
// being compared against.
// C = 2 q_eff^2 k|q1 q2| / (6 pi eps0 c^3 mu^2), with q_eff the pair's
// effective dipole charge and mu the reduced mass.  For e+e- (q_eff = e,
// mu = m/2) this collapses to the familiar 8 k e^4/(6 pi eps0 c^3 m^2), which
// is what the previous form hard-coded along with the equal-mass assumption.
// The active pair is selected by --pair, so these inputs are runtime globals.
// Keeping this function constexpr is ill-formed even though GCC historically
// accepted it as an extension with -Winvalid-constexpr.
inline double classicalInspiralCoefficient() {
    return 2.0*pairDipoleCharge*pairDipoleCharge*pairCoulombStrength
        /(6.0*pi*epsilon0*c*c*c*pairReducedMass*pairReducedMass);
}

// Orbit-averaged enhancement of DIPOLE radiation on a Kepler ellipse.
// P ~ r^-4, and <r^-4> = (1 + e^2/2)/(a^4 (1-e^2)^(5/2)) by direct
// integration with dt = (r^2/h) dtheta.  This is NOT the Peters factor
// (1 + 73e^2/24 + 37e^4/96)/(1-e^2)^(7/2), which belongs to QUADRUPOLE
// (gravitational) radiation with P ~ r^-6 and is far larger: at e = 0.5 the
// two differ by more than a factor of two.
double dipoleEccentricityFactor(double eccentricity) {
    const double e2=eccentricity*eccentricity;
    if(!(e2<1.0)) return std::numeric_limits<double>::quiet_NaN();
    return (1.0+0.5*e2)/std::pow(1.0-e2,2.5);
}

// Larmor power of the coherent electric dipole, orbit-averaged over the
// Kepler ellipse with the given elements.
double larmorOrbitAveragedPower(double semiMajorAxis,double eccentricity) {
    if(!(semiMajorAxis>0.0)) return std::numeric_limits<double>::quiet_NaN();
    // |d''| = |q_eff| k|q1 q2| / (mu a^2), so P = |d''|^2/(6 pi eps0 c^3).
    const double secondDerivative=magnitude(pairDipoleCharge)
        *pairCoulombStrength/(pairReducedMass*semiMajorAxis*semiMajorAxis);
    const double circular=secondDerivative*secondDerivative
        /(6.0*pi*epsilon0*c*c*c);
    return circular*dipoleEccentricityFactor(eccentricity);
}

// --- Eccentric-orbit harmonic content, CREM_HARMONIC ---
//
// electricDipoleRadiatedPower/larmorOrbitAveragedPower above give the exact
// TOTAL power for an eccentric Kepler dipole (dipoleEccentricityFactor is
// the correct orbit average, not a circular-orbit approximation).  What they
// do not address is how that power is DISTRIBUTED in frequency: an eccentric
// dipole moment is not a pure sinusoid at the orbital frequency omega_orb,
// it has a full harmonic series n*omega_orb (n=1,2,3,...), because the
// motion is fast and radiates strongly near periapsis and slow (radiating
// weakly) near apoapsis -- exactly the classical "impulsive burst" mechanism
// that broadens a time-domain pulse into a wide frequency-domain spectrum.
// stochasticElectricDipole's hazard/photon-energy machinery, by contrast,
// has always used ONLY the fundamental (photonEnergyReference=hbar*omega_orb)
// -- correct in the n=1-dominated near-circular limit this code started in,
// but the spin-based angular-momentum fix (see stochasticElectricDipole's
// own comment) now routinely drives orbits to e^2~0.9-0.95, where n=1 is far
// from dominant.
//
// Measured directly (Python, Kepler's equation solved numerically, DFT of
// the resulting r(t)e^{i nu(t)} against mean anomaly, N=8192 samples): at
// e=0.9, only 0.37% of radiated POWER is in the fundamental; the harmonic
// carrying the MOST power is n=16, the median is n=27, 90% of power needs
// n<=62.  At e=0.945 (the eccentricity actually reached investigating the
// guard fix above): peak n=39, median n=66.  The harmonic where power PEAKS
// fits, to 3 significant figures across e=0.5..0.97 (checked against
// e=0.9/0.945/0.97 independently, residual <0.6%):
//
//     n_peak(e) ~= 0.504 * (1-e)^-1.5
//
// This changes two separate things, not one, and they are NOT the same
// correction:
//
// (1) THE EVENT RATE.  A photon's own quantum ties it to ONE specific
//     harmonic, so treating each harmonic n as its own independent Poisson
//     channel (rate = power_n/(hbar*n*omega_orb), the same Bohr-
//     correspondence logic already applied to n=1 alone, just generalized)
//     gives a total hazard rate hazard = (P_total/(hbar*omega_orb)) * S(e),
//     where S(e) = sum_n (power fraction at harmonic n)/n.  Measured
//     (same numerical decomposition): S(e) drops from 1 at e=0 to 0.061 at
//     e=0.9, 0.022 at e=0.945, 0.0010 at e=0.99 -- far fewer total photon
//     events than the n=1-only estimate, because dividing by n suppresses
//     the very high harmonics that carry most of the POWER but would, if
//     counted at n=1's rate, wildly over-count events.
//
// (2) WHICH HARMONIC A GIVEN EVENT IS.  Given that an event fires, the
//     probability it belongs to harmonic n is proportional to (power
//     fraction at n)/n -- NOT to the power fraction alone: a rare, huge
//     event still only happens once, however much power it carries.
//     Checked: even count-weighted, the modal single harmonic is still
//     n=1 at every tested e (dividing by n keeps low n individually most
//     likely) -- but the DISTRIBUTION is heavy-tailed enough that the
//     MEDIAN event is still far from n=1 (median n=11 at e=0.9, n=27 at
//     e=0.945, n=64 at e=0.97): only ~6% of events are literally n=1 at
//     e=0.9, dropping to ~1.4% at e=0.97.  So "many small photons instead
//     of few big ones" was the right qualitative worry, just imprecise on
//     the exact split between count and size.
//
// What is NOT re-derived: the (1+cos^2 theta) angular pattern and the
// P(h=+-1|theta) helicity split (stochasticElectricDipole's own comment)
// are kept UNCHANGED for every harmonic.  Justified, not assumed: for
// PLANAR Kepler motion, x(t)+iy(t) decomposed into z=x+iy has a
// "prograde" component (co-rotating with the orbit) and a "retrograde"
// one; a pure rotating dipole (n=1, e=0) is 100% prograde, which is
// exactly what makes (1+cos^2 theta)/Delta-m=+1 exact for it.  Checked
// directly at the harmonics that actually matter (peak/median/90%-power n,
// e=0.9/0.945/0.97): retrograde power is 0.10-0.95%, i.e. these harmonics
// are themselves >99% circularly polarized in the SAME sense as the
// orbit -- so reusing the n=1 angular/helicity machinery at the SAMPLED
// harmonic's frequency is a <1%-level approximation, not a new one.
//
// S(e) below is tabulated directly per eccentricity from the same
// numerical decomposition -- no closed form was found (or, unlike the
// Cardano cubic elsewhere in this file, safely re-derivable under this
// session's time budget) for the harmonic content of an eccentric
// Kepler dipole, so this is a verified numerical table, not an
// unverified guess dressed as one.
//
// eccentricOrbitHarmonicNumber (further below) used to be built the
// same way but on a SEPARATE, scale-invariance assumption: fit a single
// "peak harmonic" n_peak(e)~=0.504(1-e)^-1.5 and one universal quantile
// shape in x=n/n_peak, checked only across e=0.75-0.98.  Checked again,
// on request, against a finer, direct per-quantile comparison at
// e=0.2-0.75 rather than left on that first check: the picture was more
// nuanced than "breaks down below e~0.6-0.7" -- the MEDIAN and 90%
// quantile were already good (within 0/+1 harmonic) all the way down to
// e=0.2, it was specifically the 99%+ tail that the universal-shape
// extrapolation got wrong outside its calibrated range (e.g. e=0.2:
// real 99% is n=2, the extrapolated table implied n=5).  A hard gate at
// e>=0.6 (the first fix) was therefore both too conservative (discarding
// real median/90% structure already present at e=0.3-0.5) and not fully
// sufficient (residual +1/+2-harmonic tail error persisted even above
// the gate, at e=0.6-0.75).  Replaced below with a directly tabulated
// 2D table (eccentricity x quantile -> harmonic, no scale-invariance
// assumption anywhere), removing the gate entirely: every entry is a
// real, independently computed number, not an extrapolation, so there
// is no boundary left to be conservative about.
double interpolateMonotonicTable(
        const double* keys,const double* values,int count,double key) {
    if(key<=keys[0]) return values[0];
    if(key>=keys[count-1]) return values[count-1];
    for(int i=1;i<count;++i) {
        if(key<=keys[i]) {
            const double span=keys[i]-keys[i-1];
            const double t=span>0.0?(key-keys[i-1])/span:0.0;
            return values[i-1]+t*(values[i]-values[i-1]);
        }
    }
    return values[count-1];
}

// S(e): hazard-rate suppression from spreading power over many harmonics
// instead of crediting it all to n=1 (see derivation above).
double eccentricOrbitHazardSuppression(double eccentricity) {
    static constexpr double keys[]={0.0000,0.0500,0.1000,0.1500,0.2000,
        0.2500,0.3000,0.3500,0.4000,0.4500,0.5000,0.5500,0.6000,0.6500,
        0.7000,0.7500,0.8000,0.8500,0.9000,0.9300,0.9500,0.9700,0.9900};
    static constexpr double values[]={1.000000,0.995009,0.980140,0.955704,
        0.922208,0.880335,0.830924,0.774946,0.713472,0.647646,0.578663,
        0.507742,0.436110,0.364990,0.295604,0.229184,0.167007,0.110485,
        0.061352,0.036433,0.022231,0.010466,0.000853};
    return interpolateMonotonicTable(keys,values,
        static_cast<int>(std::size(keys)),
        std::clamp(eccentricity,0.0,0.999));
}

// Count-weighted (power_n/n) quantile of harmonic number, tabulated
// DIRECTLY over (eccentricity, cumulative probability) -- see the
// comment above interpolateMonotonicTable for why this replaced an
// earlier, scale-invariance-based version.  18 eccentricities x 16
// quantile levels, each entry independently computed (same DFT-based
// decomposition as eccentricOrbitHazardSuppression, N=8192 samples,
// nmax up to 6000, i.e. this is not an extrapolation of a smaller
// calibration set the way the replaced version was); bilinearly
// interpolated in both dimensions.  The final quantile column (u=1)
// caps each row at 1.5x its own 0.999 entry -- a deliberately
// conservative ceiling for the near-never-hit u->1 edge of a continuous
// draw, not a truncation artifact the way an earlier draft's naive
// summation cutoff was.
double eccentricOrbitHarmonicNumber(double eccentricity,double uniformDraw) {
    static constexpr double eccentricityGrid[]={0.0000,0.1000,0.2000,
        0.3000,0.4000,0.5000,0.6000,0.6500,0.7000,0.7500,0.8000,0.8500,
        0.9000,0.9300,0.9450,0.9600,0.9700,0.9800};
    static constexpr double quantileGrid[]={0.0000,0.0500,0.1000,0.2000,
        0.3000,0.4000,0.5000,0.6000,0.7000,0.8000,0.9000,0.9500,0.9900,
        0.9950,0.9990,1.0000};
    static constexpr double harmonicTable[][16]={
        {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,3},
        {1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,4},
        {1,1,1,1,1,1,1,1,1,1,2,2,3,3,4,6},
        {1,1,1,1,1,1,1,1,1,2,2,3,4,5,6,9},
        {1,1,1,1,1,1,1,2,2,2,3,4,6,6,8,12},
        {1,1,1,1,1,1,2,2,3,3,4,6,8,10,12,18},
        {1,1,1,1,1,2,2,3,3,4,6,7,10,12,15,22},
        {1,1,1,1,2,2,3,3,4,5,7,9,13,15,20,30},
        {1,1,1,1,2,3,3,4,5,7,9,12,18,20,26,39},
        {1,1,1,2,3,3,4,6,7,9,13,17,25,29,38,57},
        {1,1,1,2,4,5,6,8,11,14,20,26,40,46,59,88},
        {1,1,2,4,6,8,11,15,19,26,37,48,74,85,111,166},
        {1,2,3,6,10,14,19,25,33,44,63,82,127,146,191,286},
        {1,2,4,9,14,20,27,35,47,63,91,118,183,210,275,412},
        {1,3,6,13,21,31,42,56,75,101,146,191,295,340,445,668},
        {1,4,9,19,32,46,64,86,115,155,224,294,455,525,687,1030},
        {1,7,15,34,56,84,116,157,209,284,411,540,840,970,1286,1929},
    };
    constexpr int eCount=static_cast<int>(std::size(eccentricityGrid));
    constexpr int qCount=static_cast<int>(std::size(quantileGrid));
    const double e=std::clamp(eccentricity,0.0,eccentricityGrid[eCount-1]);
    const double u=std::clamp(uniformDraw,0.0,1.0);
    int ei=1;
    while(ei<eCount-1&&e>eccentricityGrid[ei]) ++ei;
    int qi=1;
    while(qi<qCount-1&&u>quantileGrid[qi]) ++qi;
    const double eSpan=eccentricityGrid[ei]-eccentricityGrid[ei-1];
    const double eT=eSpan>0.0?(e-eccentricityGrid[ei-1])/eSpan:0.0;
    const double qSpan=quantileGrid[qi]-quantileGrid[qi-1];
    const double qT=qSpan>0.0?(u-quantileGrid[qi-1])/qSpan:0.0;
    const double lowELowQ=harmonicTable[ei-1][qi-1];
    const double lowEHighQ=harmonicTable[ei-1][qi];
    const double highELowQ=harmonicTable[ei][qi-1];
    const double highEHighQ=harmonicTable[ei][qi];
    const double lowE=lowELowQ+qT*(lowEHighQ-lowELowQ);
    const double highE=highELowQ+qT*(highEHighQ-highELowQ);
    return lowE+eT*(highE-lowE);
}

// Time for the closed-form inspiral to carry the orbit from a_i to a_f at
// fixed eccentricity.  Radiation actually circularizes the orbit, so holding
// e fixed is an approximation; it is stated on the panel that uses this.
double classicalInspiralSeconds(double initialSemiMajorAxis,
                                double finalSemiMajorAxis,
                                double eccentricity) {
    const double factor=dipoleEccentricityFactor(eccentricity);
    if(!(initialSemiMajorAxis>finalSemiMajorAxis)||!std::isfinite(factor))
        return std::numeric_limits<double>::quiet_NaN();
    return (std::pow(initialSemiMajorAxis,3.0)-std::pow(finalSemiMajorAxis,3.0))
        /(3.0*classicalInspiralCoefficient()*factor);
}

// Periapsis distance of the two-body Kepler ellipse with the given specific
// (per unit reduced mass) energy and angular momentum -- the same formula
// simulate() already uses for its InitialConditions.predictedClosestApproach.
double osculatingPeriapsis(const OsculatingElements& elements,
                           double attractionParameter) {
    if(!(elements.specificAngularMomentum!=0.0)) return 0.0;
    const double eccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*elements.specificAngularMomentum
            *elements.specificAngularMomentum
            /(attractionParameter*attractionParameter)));
    return elements.specificAngularMomentum*elements.specificAngularMomentum
        /(attractionParameter*(1.0+eccentricity));
}

// Apoapsis of the two-body Kepler ellipse -- osculatingPeriapsis's mirror
// root (a(1+e) instead of a(1-e)).  Infinity for an unbound (e>=1) orbit.
double osculatingApoapsis(const OsculatingElements& elements,
                          double attractionParameter) {
    const double L=elements.specificAngularMomentum;
    if(!(L!=0.0)) return 0.0;
    const double eccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*L*L
            /(attractionParameter*attractionParameter)));
    if(!(eccentricity<1.0)) return std::numeric_limits<double>::infinity();
    return L*L/(attractionParameter*(1.0-eccentricity));
}

// Unperturbed Kepler period at the given specific energy; used only to size
// the one-orbit measurement window, not to compute the reported lifetime.
double osculatingPeriod(double specificEnergy,double attractionParameter) {
    const double semiMajorAxis=-attractionParameter/(2.0*specificEnergy);
    return 2.0*pi*std::sqrt(semiMajorAxis*semiMajorAxis*semiMajorAxis
        /attractionParameter);
}

// Radial potential energy (per unit reduced mass) for the Plummer-softened
// force law: the specific radial force is
// F_r(r) = -attractionParameter/(r^2+floor^2), matching the true Coulomb
// force -attractionParameter/r^2 for r >> floor and saturating smoothly
// (rather than diverging) as r -> 0, with clampedSeparationVector in
// positronium.cpp (r_eff = sqrt(r^2+floor^2)) as the force-law-level source
// of this.  Integrating -dU/dr = F_r gives
//
//     U(r) = -(attractionParameter/floor) * atan(floor/r),
//
// with the integration constant fixed so U -> -attractionParameter/r as
// r -> infinity (atan(floor/r) ~ floor/r there), matching the unclamped
// Coulomb potential exactly in that limit.  Valid and smooth for every
// r > 0: unlike the earlier hard-clamp regularization, there is no separate
// branch for r above/below floor here at all -- one formula, everywhere,
// which is what removes the kink in the force that stalled the adaptive
// integrator at the old floor crossing.
double regularizedPotentialEnergy(double r,double attractionParameter,
                                  double floor) {
    return -(attractionParameter/floor)*std::atan(floor/r);
}

// Radius where the softened force above exactly balances the centrifugal
// term for angular momentum L -- the circular-orbit radius under this force
// law.  Setting -dU/dr = L^2/r^3 gives
// attractionParameter/(r^2+floor^2) = L^2/r^3, i.e.
//
//     attractionParameter*r^3 - L^2*r^2 - L^2*floor^2 = 0,
//
// a cubic with exactly one positive real root (Descartes: coefficient signs
// +,-,0,- change sign once), interpolating smoothly between the deep limit
// (floor dominates, r ~ floor) and the standard Kepler circular-orbit radius
// L^2/attractionParameter (the floor -> 0 limit -- the identity this reduces
// to away from the barrier).  No simpler closed form survives adding the
// floor^2 term, so this is found by bisection.  The bracket is a physically
// motivated guess (floor plus the naive circular radius), doubled until it
// brackets a root as cheap insurance against a bad guess rather than a claim
// of precision: bisection's absolute convergence over even a wildly
// oversized bracket reaches machine precision at the relevant scale well
// within the iteration budget below regardless of how tight the guess was.
double criticalRadius(double L,double attractionParameter,double floor) {
    const auto cubic=[&](double r) {
        return attractionParameter*r*r*r-L*L*(r*r+floor*floor);
    };
    double lo=std::numeric_limits<double>::min();
    double hi=2.0*(floor+L*L/attractionParameter);
    for(int i=0;i<200&&!(cubic(hi)>0.0);++i) hi*=2.0;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(lo+hi);
        if(cubic(mid)>0.0) hi=mid; else lo=mid;
    }
    return 0.5*(lo+hi);
}

struct RegularizedTurningPoints {
    double periapsis,apoapsis;
    // False only in the naive-fallback branch below: there periapsis/
    // apoapsis are NOT roots of h(r) (h(r_min) missed even the marginal
    // circular-orbit tolerance, so no such roots exist), they are the plain
    // osculatingPeriapsis/osculatingApoapsis values returned as the
    // least-bad guess for STATE CONSTRUCTION (regularizedPeriapsis, via
    // osculatingPeriapsisState, where any single r is better than none).
    // regularizedPeriod's quadrature integrates h(r) between periapsis and
    // apoapsis, which is only valid when they truly are h's roots -- with
    // the naive fallback pair instead, h is negative at every quadrature
    // node (that is exactly what "not even marginal" means), so every node
    // gets skipped and the sum silently comes out zero.  Measured directly
    // (seed 5, the checkpoint where h(r_min)/localScale missed the 5%
    // margin by roughly a percentage point): regularizedPeriod returned
    // 0 ps before this flag existed to stop it from even trying.
    bool exact;
};

// Periapsis and apoapsis under the Plummer-softened force law
// (regularizedPotentialEnergy), replacing the naive (pure 1/r) formulas
// osculatingPeriapsis/osculatingApoapsis wherever the orbit's own scale
// approaches separationFloor().  h(r) = specificEnergy -
// regularizedPotentialEnergy(r) - L^2/(2r^2) is smooth for every r > 0 --
// no piecewise construction, unlike the hard-clamp regularization this
// replaces, since the force itself is now one formula everywhere -- concave
// with a single interior maximum at r = criticalRadius(...), and
// h -> -infinity as r -> 0 (centrifugal term dominates) or r -> infinity
// (bound orbit, specificEnergy < 0), so it has at most two roots straddling
// that maximum, found by bisecting on each side of it.  Cost is negligible
// (a few hundred evaluations of an elementary function) against the
// mechanical integration this sizes a measurement for, so -- unlike the
// hard-clamp predecessor -- there is no need to shortcut to the naive
// formula when the orbit is far from the barrier: the exact formula already
// reduces to it there (the correction is a part in 10^7 at Bohr-radius
// scale), so always using it is both simpler and correct everywhere.
//
// If h(r_min) itself is not positive, the (E,L) pair has no strictly
// accessible orbit at all under this force law -- but h is smooth in (E,L),
// so a small negative h(r_min) means the pair is almost exactly circular AT
// r_min (h(r_min)=0 is precisely that circular orbit), not that no orbit
// exists.  Comparing |h(r_min)| against the local force scale there
// (attractionParameter*r_min/(r_min^2+floor^2), i.e. force(r_min)*r_min, the
// specific potential step across a distance ~r_min) keeps the margin
// scale-free.  Within that margin both turning points collapse to r_min (a
// genuinely circular orbit); beyond it the naive values are returned as the
// least-bad fallback, since fabricating a root that is not there would be
// worse than saying so.
RegularizedTurningPoints regularizedTurningPoints(
    const OsculatingElements& elements,double attractionParameter,
    double floor) {
    const double naivePeriapsis=
        osculatingPeriapsis(elements,attractionParameter);
    const double naiveApoapsis=
        osculatingApoapsis(elements,attractionParameter);
    const double L=elements.specificAngularMomentum;
    if(!(floor>0.0)||!(L!=0.0)) return {naivePeriapsis,naiveApoapsis,false};
    const double rMin=criticalRadius(L,attractionParameter,floor);
    if(!(rMin>0.0)||!std::isfinite(rMin))
        return {naivePeriapsis,naiveApoapsis,false};
    const auto h=[&](double r) {
        return elements.specificEnergy
            -regularizedPotentialEnergy(r,attractionParameter,floor)
            -L*L/(2.0*r*r);
    };
    const double hAtMin=h(rMin);
    if(!(hAtMin>0.0)) {
        const double localScale=
            attractionParameter*rMin/(rMin*rMin+floor*floor);
        constexpr double marginalRelativeTolerance=0.05;
        if(localScale>0.0&&-hAtMin<=marginalRelativeTolerance*localScale)
            return {rMin,rMin,true};
        return {naivePeriapsis,naiveApoapsis,false};
    }
    double periapsisLo=std::numeric_limits<double>::min();
    double periapsisHi=rMin;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(periapsisLo+periapsisHi);
        if(h(mid)>0.0) periapsisHi=mid; else periapsisLo=mid;
    }
    double apoapsisLo=rMin;
    double apoapsisHi=(std::isfinite(naiveApoapsis)&&naiveApoapsis>rMin)
        ?naiveApoapsis:2.0*rMin;
    for(int i=0;i<200&&!(h(apoapsisHi)<0.0);++i) apoapsisHi*=2.0;
    for(int i=0;i<200;++i) {
        const double mid=0.5*(apoapsisLo+apoapsisHi);
        if(h(mid)>0.0) apoapsisLo=mid; else apoapsisHi=mid;
    }
    return {0.5*(periapsisLo+periapsisHi),0.5*(apoapsisLo+apoapsisHi),true};
}

double regularizedPeriapsis(const OsculatingElements& elements,
                            double attractionParameter,double floor) {
    return regularizedTurningPoints(elements,attractionParameter,floor)
        .periapsis;
}

// Regularized orbital period, replacing osculatingPeriod() once the orbit's
// own scale approaches separationFloor().  osculatingPeriod() sizes the
// one-orbit measurement window in estimateCremCollapse(); once
// regularizedTurningPoints() starts differing materially from the naive
// values, that window is no longer one true period of the actual (softened)
// orbit, and the mismatch compounds into the jump extrapolation the same
// way periapsis's own mismatch did before it was fixed (c909ee8/2bb2751).
//
// No closed form exists here in general -- the softened force does not
// admit the standard eccentric-anomaly substitution -- so the period is
// found by direct quadrature of T = 2 * integral[r_p,r_a] dr / sqrt(2 h(r)),
// using the substitution r(theta) = (r_p+r_a)/2 + (r_a-r_p)/2 * cos(theta),
// theta in (0,pi): this removes BOTH endpoint (turning-point) singularities
// for ANY smooth h with simple roots there, because near either root
// h(r) ~ h'(r_root)*(r-r_root) while r(theta)-r_root ~ -/+(r_a-r_p)/4 *
// theta^2 near theta=0/pi, so h(r(theta)) ~ theta^2 while dr/dtheta ~ theta,
// and the ratio in the integrand stays finite.  h is now the single smooth
// formula everywhere (no piecewise potential to match across a boundary the
// way the hard-clamp version needed), so this integral is, if anything,
// better conditioned than its predecessor: there is no interior kink in
// curvature left for the quadrature to sit near.
//
// Fixed-node midpoint quadrature is used rather than Gauss-Legendre so no
// tabulated node/weight constants need transcribing; every node falls
// strictly inside (0,pi), so the removable endpoint singularities are never
// evaluated at all.  N=200 carries over unchanged from the hard-clamp
// version (checked there to machine precision against the closed form in
// the far-field limit); cost is negligible either way (a few hundred h(r)
// evaluations against the full mechanical integration this sizes the window
// for), and always running the same quadrature -- rather than shortcutting
// to osculatingPeriod() far from the barrier -- keeps this function as
// simple as regularizedTurningPoints() above for the same reason.
double regularizedPeriod(const OsculatingElements& elements,
                         double attractionParameter,double floor) {
    const double naive=osculatingPeriod(
        elements.specificEnergy,attractionParameter);
    if(!(floor>0.0)) return naive;
    const RegularizedTurningPoints turningPoints=
        regularizedTurningPoints(elements,attractionParameter,floor);
    // Naive fallback (turningPoints.exact == false): periapsis/apoapsis are
    // NOT roots of h(r) below, so integrating between them is meaningless --
    // h is negative at every quadrature node (that is exactly what "missed
    // even the marginal tolerance" means), every node gets skipped, and the
    // sum silently comes out zero.  Measured directly: this is precisely
    // what happened before this check existed (seed 5, one checkpoint whose
    // h(r_min)/localScale missed the 5% margin by about a point -- naive is
    // still the least-bad state-construction guess regularizedPeriapsis
    // hands to osculatingPeriapsisState in that case, but it is not a valid
    // integration bound here, so this falls back to osculatingPeriod()
    // instead of the quadrature).
    if(!turningPoints.exact) return naive;
    const double periapsis=turningPoints.periapsis;
    const double apoapsis=turningPoints.apoapsis;
    const double L=elements.specificAngularMomentum;
    if(!(apoapsis>periapsis)) {
        // Degenerate: both turning points collapsed onto (about) the same
        // radius -- regularizedTurningPoints' own marginal branch, an orbit
        // sitting almost exactly at the bottom of the effective-potential
        // well.  The RADIAL period is not just small there, it is the wrong
        // question: for a force law that does not satisfy Bertrand's
        // theorem, radial and angular period differ, and "one orbit" for a
        // near-circular orbit means one REVOLUTION, not one (vanishing)
        // radial oscillation.  phi-dot = L/r^2 by definition of specific
        // angular momentum, so the angular period is 2*pi*r^2/L at the
        // shared radius.
        const double r=0.5*(periapsis+apoapsis);
        return (L!=0.0)?2.0*pi*r*r/std::abs(L):naive;
    }
    const auto h=[&](double r) {
        return elements.specificEnergy
            -regularizedPotentialEnergy(r,attractionParameter,floor)
            -L*L/(2.0*r*r);
    };
    constexpr int quadratureNodes=200;
    double halfPeriod=0.0;
    for(int i=0;i<quadratureNodes;++i) {
        const double theta=(static_cast<double>(i)+0.5)*pi
            /static_cast<double>(quadratureNodes);
        const double r=0.5*(periapsis+apoapsis)
            +0.5*(apoapsis-periapsis)*std::cos(theta);
        const double radialKineticEnergyTimesTwo=2.0*h(r);
        if(!(radialKineticEnergyTimesTwo>0.0)) continue; // guard only; every
            // node sits strictly between the turning points, so this should
            // not trigger, but a stray non-finite h must not poison the sum.
        halfPeriod+=std::sin(theta)/std::sqrt(radialKineticEnergyTimesTwo);
    }
    halfPeriod*=0.5*(apoapsis-periapsis)*pi/static_cast<double>(quadratureNodes);
    return 2.0*halfPeriod;
}

// Dipole-dipole interaction energy at separation r, averaged over the ONE
// thing the osculating representation genuinely does not know: where in its
// orbital plane the pair is.
//
// osculatingPeriapsisState resets the plane to canonical x-y and puts the
// separation along +x, and says so -- that is sound for a collapse-time
// estimate, which cares only about |E| and |L|.  It is NOT sound for the
// angular factor mu1.mu2 - 3(mu1.n)(mu2.n), which depends on the direction
// of n relative to the moments.  Using that convention's n = x gives an
// answer that is not merely imprecise but can carry the WRONG SIGN: checked
// against the correct average on random configurations, -0.268 against
// +0.484 in one of four trials.
//
// The true anomaly at termination is unknowable here for the same reason the
// photon-firing code states for its own emission azimuth: only the elements
// survive the secular step, not the phase.  What IS tracked is the orbital
// plane normal.  So the azimuth is averaged out rather than guessed, using
//
//     <n_i n_j> = (delta_ij - L_i L_j)/2   for n uniform in the plane,
//
// which turns the angular factor into
//
//     <mu1.mu2 - 3(mu1.n)(mu2.n)> = -(mu1.mu2)/2 + 3(mu1.L)(mu2.L)/2,
//
// verified against brute-force azimuth integration to 1e-16.  No 1/r^3
// weighting enters: r is fixed at the terminal periapsis, so the average is
// uniform in azimuth at that one radius, not an orbit average.
//
// The radial profile is the regularized one the field code uses, so this
// agrees with regularizedDipoleInteractionEnergy term by term and reduces to
// the point-dipole result where the regulator is inactive.
double azimuthAveragedDipoleEnergy(double separation,
                                   const Vec3& firstDipole,
                                   const Vec3& secondDipole,
                                   const Vec3& orbitNormal) {
    if(!(separation>0.0)) return 0.0;
    const MagneticRadialProfile profile=magneticRadialProfile(separation);
    const double transverse=2.0*profile.vectorPotentialFactor
        +separation*profile.firstDerivative;
    const double radial=-separation*profile.firstDerivative;
    const double moments=dot(firstDipole,secondDipole);
    const double alongNormal=dot(firstDipole,orbitNormal)
        *dot(secondDipole,orbitNormal);
    // U = -mu1 . B_reg(mu2), with <(mu1.n)(mu2.n)> substituted for the
    // unknown azimuth.
    const double averagedRadial=0.5*(moments-alongNormal);
    return -(mu0/(4.0*pi))*(moments*transverse+averagedRadial*radial);
}

// Fresh State at periapsis for the given osculating elements, carrying the
// supplied dipole vectors over unchanged.  The orbital plane is reset to the
// canonical x-y plane every time: only the (E,L) magnitudes are propagated
// secularly, and the plane orientation is irrelevant to a collapse-time
// estimate, so nothing is lost by not tracking it.
State osculatingPeriapsisState(const OsculatingElements& elements,
                               double attractionParameter,
                               const Vec3& firstDipole,
                               const Vec3& secondDipole) {
    // regularizedPeriapsis(), not the naive formula: below separationFloor()
    // the naive value systematically overshoots how close the pair can
    // actually get (it assumes the force keeps growing as 1/r^2 rather than
    // being pinned), so teleporting a measurement orbit there starts it from
    // a position/speed pair that is not a real turning point of the clamped
    // force law -- exactly the mismatch that made the one-orbit measurement
    // fail outright once periapsis approached the barrier.  Tangential speed
    // is still L/r at any turning point regardless of the force law (that is
    // the definition of specific angular momentum at zero radial velocity),
    // so only the r fed into it needs correcting.
    const double periapsis=regularizedPeriapsis(
        elements,attractionParameter,separationFloor());
    const double tangentialSpeed=periapsis>0.0
        ?elements.specificAngularMomentum/periapsis:0.0;
    const Vec3 relativePosition{periapsis,0.0,0.0};
    const Vec3 relativeVelocity{0.0,tangentialSpeed,0.0};
    State s;
    s.firstPosition=relativePosition*(secondMass/(firstMass+secondMass));
    s.secondPosition=relativePosition*(-firstMass/(firstMass+secondMass));
    s.firstVelocity=relativeVelocity*(secondMass/(firstMass+secondMass));
    s.secondVelocity=relativeVelocity*(-firstMass/(firstMass+secondMass));
    s.firstDipole=firstDipole;
    s.secondDipole=secondDipole;
    return s;
}

// Orbit-averaged secular integration.  Resolving each of the ~1e5-1e6
// individual orbits a classical inspiral needs (the direct approach this
// replaces) took upward of an hour per trajectory.  Instead, one
// representative orbit is fully resolved with the complete CREM engine --
// retardation, Darwin, dipole coupling, and radiation reaction under
// whichever --radiation-reaction model is active -- starting and ending near
// periapsis.  Its ACTUALLY MEASURED energy and angular-momentum loss (read
// directly from the engine's own exact radiatedEnergy / radiatedAngularMomentum
// flux bookkeeping, not a formula) gives the current dissipation rate.  A
// bounded number of further orbits, capped so the implied energy change per
// jump stays small, are then skipped analytically at that rate before the
// next orbit is resolved and the rate re-measured.  This is the standard
// osculating-element / orbit-averaging technique for multi-timescale
// radiative inspirals (the same approximation EMRI gravitational-wave models
// use), and it stays exactly as sensitive to --radiation-reaction as brute
// force would be: with the reaction force disabled, the very first measured
// orbit shows zero energy loss and the estimate correctly reports no
// collapse instead of iterating uselessly.
CremCollapseEstimate estimateCremCollapse(std::uint64_t seed,
                                           int selectedPhenomenon,
                                           double wallClockBudgetSeconds) {
    CremCollapseEstimate result;
    const double reducedMass=firstMass*secondMass
        /(firstMass+secondMass);
    const double attractionParameter=pairCoulombStrength/reducedMass;

    // Sample the same random bound initial condition simulate() would use.
    // The near-zero observation window means nothing actually integrates;
    // frames.front() is the pristine sampled state (pushed before the loop).
    SimulationOptions seedOptions;
    seedOptions.frameCount=2;
    seedOptions.observationTime=1.0e-24;
    const SimulationResult seedRun=simulate(seed,selectedPhenomenon,seedOptions);
    if(seedRun.frames.empty()||!(seedRun.initial.relativeEnergy<0.0))
        return result;

    OsculatingElements elements{seedRun.initial.relativeEnergy/reducedMass,
        seedRun.initial.orbitalAngularMomentum/reducedMass};
    Vec3 firstDipole=seedRun.frames.front().firstDipole;
    Vec3 secondDipole=seedRun.frames.front().secondDipole;
    // Orbital plane orientation, kept ONLY for stochasticElectricDipole (see
    // its own use below): OsculatingElements itself carries no direction,
    // only magnitudes, which is exactly right for the continuous models
    // (an in-plane reaction force never tilts the plane, so the direction
    // would just sit there unread) but not for individual photon kicks,
    // whose emission angle relative to this axis is physical, not
    // arbitrary.
    //
    // FIXED: used to be seedRun.frames.front().noetherAngularMomentum,
    // normalized, on the claim (in this very comment) that "the dipole/spin
    // contribution to it is already documented elsewhere as ~1e-5 of the
    // orbital term".  That claim was never actually checked here and was
    // wrong: it confused two different ~1e-5-ish numbers (the DIPOLE-DIPOLE
    // INTERACTION ENERGY's ratio to the Coulomb energy, which genuinely is
    // that small -- see stochasticElectricDipole's own comment -- against
    // the SPIN CONTRIBUTION TO ANGULAR MOMENTUM, which is not).  The spin
    // angular momentum per particle is mu/gyromagneticRatio ~ hbar/2 -- the
    // SAME order as the orbital L this state starts at (L=hbar, the Bohr/SED
    // value -- see physical_constants.hpp), not five orders smaller.  Found
    // investigating why para and ortho positronium (same seed, same orbital
    // sampling, differing only in dipole alignment) gave measurably
    // different photon emission directions from the very first photon:
    // instrumented angularMomentumDirection directly and found it was a
    // COMPLETELY different vector between the two channels (not a small
    // tilt), present before any photon had fired, because
    // noetherAngularMomentum includes the mu_e/gamma_e+mu_p/gamma_p spin
    // term, and para/ortho sample genuinely different (not just
    // sign-flipped) dipole vectors under the rejection condition that
    // defines them.  Every stochastic-model trajectory, not only the
    // para/ortho comparison, was sampling photon direction relative to a
    // spin-contaminated axis instead of the true orbital plane normal.
    //
    // Fixed by computing the orbital angular momentum direction directly,
    // the way the comment this replaces said it was avoiding: cross(r0,r1)
    // from the two frames seedOptions above already requests (frameCount=2,
    // observationTime=1e-24s, i.e. two position samples a near-instant
    // apart on the true trajectory).  For r1=r0+v*dt at small dt,
    // cross(r0,r1)=cross(r0,r0)+dt*cross(r0,v)=dt*cross(r0,v) exactly to
    // leading order -- proportional to the true orbital angular momentum
    // direction, with no spin contamination at all, and no dependence on
    // any particular sampling convention (unlike hard-coding the z-axis,
    // which happens to be right for how this file's own bound-state sampler
    // places the initial separation along x with velocity confined to the
    // xy-plane, but would silently break for any other caller). Checked
    // directly (this exact seed): normalizes to (0, 2e-13, 1) -- the z-axis
    // recovered to floating-point noise, not the O(1) difference the spin
    // contamination was producing.
    const Vec3 firstRelativePosition=
        seedRun.frames.front().first-seedRun.frames.front().second;
    const Vec3 secondRelativePosition=
        seedRun.frames.back().first-seedRun.frames.back().second;
    Vec3 angularMomentumDirection=
        cross(firstRelativePosition,secondRelativePosition);
    angularMomentumDirection=angularMomentumDirection.squaredNorm()>0.0
        ?angularMomentumDirection*(1.0/angularMomentumDirection.norm())
        :Vec3{0.0,0.0,1.0};

    // Osculating elements of the prepared orbit and the closed-form classical
    // prediction they imply.  The run is stopped when the PERIAPSIS reaches
    // finalApproachMultiple*comptonBarrierRadius, so the reference must be
    // evaluated between the same two orbits, not down to a = 0.
    result.initialSemiMajorAxis=
        -attractionParameter/(2.0*elements.specificEnergy);
    result.initialEccentricity=std::sqrt(std::max(0.0,1.0
        +2.0*elements.specificEnergy*elements.specificAngularMomentum
            *elements.specificAngularMomentum
            /(attractionParameter*attractionParameter)));
    // Classical dipole-dipole interaction energy of the prepared pair, as a
    // frequency, for comparison with the measured hyperfine splitting.
    {
        const Vec3 separation=seedRun.frames.front().first
                             -seedRun.frames.front().second;
        const double distance=separation.norm();
        if(distance>0.0) {
            const Vec3 direction=separation/distance;
            const double coupling=(mu0/(4.0*pi))
                *(dot(firstDipole,secondDipole)
                  -3.0*dot(firstDipole,direction)
                      *dot(secondDipole,direction))
                /(distance*distance*distance);
            result.dipoleCouplingHz=coupling/(2.0*pi*hbar);
        }
    }
    double larmorRatioSum=0.0;
    int larmorRatioCount=0;
    // Cached restart-artefact ratio and the osculating energy it was taken at.
    // Measured trade-off on one full trajectory (seed 42), against the
    // un-thinned collapse time of 131.685 ps:
    //
    //     threshold   wall clock   collapse time   error
    //       0.02        16.27 s      131.685 ps    0.000%   (refreshes every
    //       0.05        13.32 s      130.602 ps    0.822%    checkpoint, so it
    //       0.10        11.26 s      127.694 ps    3.031%    reproduces the
    //                                                        un-thinned answer
    //                                                        exactly)
    //
    // 0.05 is the working point: the error stays well inside the 3% per-jump
    // tolerance the frozen rate already accepts, so the background stays a
    // subordinate approximation.  0.10 was tried first on the assumption that
    // it would too; it does not -- there the background error EQUALS the
    // frozen-rate tolerance and would become co-dominant.  Raising it buys
    // speed at a bias that grows faster than linearly.
    constexpr double backgroundRefreshFraction=0.05;
    bool haveBackgroundRatio=false;
    double backgroundEnergyRatio=0.0;
    double backgroundAngularRatio=0.0;
    double energyAtLastBackground=0.0;
    // Read once, thread-safely: runCremCollapseExperiment runs many of these
    // concurrently, so the global must not be re-read (or, worse, mutated)
    // from here on.
    const ChargeRadiationReactionModel activeReactionModel=gRadiationReactionModel;

    // Poisson-process bookkeeping for stochasticElectricDipole, unused (and
    // costing nothing) for every other model.  Unlike the analogous state in
    // runMechanicalTrajectory (reset fresh per call, correct there because a
    // stochastic MECHANICAL trajectory never spans more than one call), this
    // one must persist across the WHOLE checkpoint loop below: each checkpoint
    // measures only one orbit before analytically SKIPPING up to
    // maxOrbitsSkippedAtOnce more, and it is that skipped span the hazard
    // needs to see (see this reaction model's own comment in
    // electrodynamics.hpp for the measurement that this exposure gap is real).
    // Seeded from this trajectory's own seed, distinctly from anything else
    // splitMix64(seed) is used for elsewhere in this file.
    std::uint64_t stochasticSkipStream=
        splitMix64(seed^0x506f6973736f6e5fULL);
    double stochasticSkipHazard=0.0;
    double stochasticSkipThreshold=drawEmissionThreshold(stochasticSkipStream);
    // Recoil bookkeeping, same lifetime/gating as the hazard state above.
    // CREM's bound initial conditions are always prepared at EXACTLY zero
    // total momentum (crem_trajectory.hpp splits the sampled relative
    // velocity by mass ratio: firstMass*v1+secondMass*v2=0 identically),
    // and every continuous reaction-force model keeps it that way -- a
    // classical self-force is, by construction, Newton's-third-law
    // consistent with its own radiated field. A photon carrying real
    // momentum away and NOT being balanced by an equal and opposite recoil
    // on the pair is exactly the gap the note in crem_trajectory.hpp's
    // applyStochasticDipolePhoton documents (and that comment's own
    // "negligible" claim, measured directly, was wrong for most of this
    // model's depth range).  centreOfMassVelocity starts at the true
    // zero and only ever moves because a photon kicked it.
    const double totalMass=firstMass+secondMass;
    Vec3 centreOfMassVelocity{0.0,0.0,0.0};

    const auto wallClockStart=std::chrono::steady_clock::now();
    const auto wallClockSpent=[&]() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now()-wallClockStart).count();
    };
    // gamma(beta), reused everywhere below that converts an S'-frame
    // duration into what a fixed lab clock reads for the same stretch of
    // proper time.  Same 1e-9 floor as the LabFramePhoton block below: pure
    // noise from a subnormal 1-beta^2 is worse than just calling it 1.
    const auto gammaFromBeta=[](double beta) {
        return beta>1.0e-9
            ?1.0/std::sqrt(std::max(1.0e-300,1.0-beta*beta))
            :1.0;
    };
    // Size of one analytic jump, as the dimensionless s = (3/2) n du0/u0 of
    // the resummed solution below.  s must stay strictly under 1, where the
    // closed form reaches complete collapse.
    //
    // Convergence on one full trajectory (seed 42), collapse time in ps:
    //
    //     s_max    wall clock   collapse time
    //     0.012      31.0 s       124.169
    //     0.045      12.6 s       124.957     <- same step as the old 3% cap
    //     0.10        8.4 s       125.585
    //     0.20        5.7 s       124.986
    //     0.30        4.6 s       124.328
    //
    // Flat to +/-0.6% over a 25x range of step size, so 0.30 is taken.  The
    // frozen-rate hold this replaces returned 130.6-131.7 ps, i.e. 5-6% HIGH,
    // and it was never convergence-tested: at the SAME step size (s = 0.045)
    // the two schemes differ by 4.5%, which is interpolation shape, not step
    // size.  The sign is the expected one -- a zeroth-order hold freezes a
    // loss rate that actually grows as the orbit tightens, so it under-counts
    // the loss and over-states the collapse time.  This change therefore buys
    // accuracy as well as speed.
    constexpr double maximumJumpParameter=0.30;
    constexpr int maxOrbitsSkippedAtOnce=200000;
    constexpr int maxCheckpoints=4000; // ~150 are needed at 3%/jump to cover
        // a factor 100 in |E|; this is a generous safety margin, not a
        // normal stopping point.
    double simulatedTimeTotal=0.0;
    // Lab-frame counterpart of simulatedTimeTotal (README point N):
    // accumulated with the SAME per-checkpoint increments, each weighted by
    // gamma(beta) for whichever constant-beta segment it covers -- see the
    // single increment point below (the stochastic photon while-loop can
    // fire several kicks within one checkpoint, splitting its elapsed proper
    // time at each one's own sAtPhoton) for the actual integration. Equal to
    // simulatedTimeTotal whenever centreOfMassVelocity never leaves zero
    // (every continuous/non-stochastic model), diverging only once photon
    // recoil has actually accumulated a nonzero pair velocity.
    double labFrameTimeTotal=0.0;
    double radiatedEnergyTotal=0.0;
    // Revolutions completed so far.  Every orbit the loop accounts for is
    // counted here exactly once: the resolved measurement orbit is the first
    // of the orbitsToSkip that each checkpoint advances over.
    double revolutionsTotal=0.0;
    // Period of the orbit the pair actually starts on, recorded at the first
    // checkpoint before any secular decay has been applied.
    result.initialPeriodSeconds=osculatingPeriod(
        elements.specificEnergy,attractionParameter);

    for(int checkpoint=0;checkpoint<maxCheckpoints;++checkpoint) {
        // Refreshed every checkpoint so that whichever branch returns below
        // reports the period of the orbit reached at that point.
        // regularizedPeriod(), not the naive formula: below separationFloor()
        // the naive value systematically under- or overestimates the true
        // period (see regularizedPeriod's own comment), and this is a
        // user-facing report, not just internal bookkeeping.
        result.finalPeriodSeconds=regularizedPeriod(
            elements,attractionParameter,separationFloor());
        result.revolutions=revolutionsTotal;
        if(larmorRatioCount>0) {
            result.larmorPowerRatio=larmorRatioSum
                /static_cast<double>(larmorRatioCount);
        }
        // Closed-form prediction for the stretch actually covered so far, so
        // the reference is always evaluated between the same two orbits as
        // the measurement it will be compared with.
        result.analyticCollapseSeconds=classicalInspiralSeconds(
            result.initialSemiMajorAxis,
            -attractionParameter/(2.0*elements.specificEnergy),
            result.initialEccentricity);
        if(const char* debug=std::getenv("CREM_DEBUG");debug) {
            std::cerr<<"checkpoint "<<checkpoint<<" t="<<simulatedTimeTotal*1e12
                     <<"ps E="<<elements.specificEnergy<<" L="
                     <<elements.specificAngularMomentum<<" wall="
                     <<wallClockSpent()<<"s"<<std::endl;
        }
        if(std::getenv("CREM_DEBUG_PRECISE")) {
            std::cerr<<std::setprecision(17)
                     <<"  PRECISE checkpoint "<<checkpoint<<": E="
                     <<elements.specificEnergy<<" L="
                     <<elements.specificAngularMomentum
                     <<std::setprecision(6)<<std::endl;
        }
        if(wallClockSpent()>wallClockBudgetSeconds) {
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal;
            result.calibrationSecondsLab=labFrameTimeTotal;
            return result;
        }
        const double periapsis=osculatingPeriapsis(elements,attractionParameter);
        if(!(periapsis>0.0)||!std::isfinite(periapsis)) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            return result;
        }

        const State measurementState=osculatingPeriapsisState(
            elements,attractionParameter,firstDipole,secondDipole);
        // regularizedPeriod(), not the naive formula: this sizes the ONE
        // orbit the measurement below actually integrates, and the naive
        // period stops matching a true orbit of the clamped force law
        // exactly where osculatingPeriapsisState's own teleport target does
        // (see regularizedPeriod's comment for the measured symptom this
        // fixes).
        const double period=regularizedPeriod(
            elements,attractionParameter,separationFloor());
        // Stop the orbit-averaged phase when EITHER of two independent
        // physical validity limits is reached, whichever comes first.
        //
        // (a) periapsis <= comptonBarrierRadius: the point-particle model-
        // validity limit ("where classical point-particle electrodynamics
        // stops applying"), independent of orbit shape -- also the Collision
        // boundary experiments 3/4/5 already declare (898dc95, 340b67b), so
        // "collapsed" means the same separation everywhere in the model.
        //
        // (b) period/lightCrossingTime <= minimumPeriodToLightCrossingRatio:
        // the Darwin/retardation-approximation validity limit UNDERLYING the
        // energy signal itself (conservativeParticleEnergy inside
        // regularizedPeriod's own turning points, and orbitalRadiatedEnergy
        // via the far-zone flux quadrature) -- and, unlike (a), NOT a fixed
        // radius: it is set by how fast the orbit's own geometry changes
        // relative to the light-crossing time at periapsis, which depends on
        // L (specific angular momentum), not on r alone.  A fixed stopping
        // radius cannot capture this: measured directly (N=100 at a fixed
        // periapsis<=comptonBarrierRadius stopping rule, i.e. condition (a)
        // alone) every one of 23 genuine NumericalFailures had this ratio
        // between 37.8 and 71.6 at the point of failure, and roughly half of
        // them (periapsis 204-676 fm) were still well ABOVE
        // comptonBarrierRadius when that happened -- low-L (more radial)
        // orbits reach the ratio limit at a LARGER periapsis than high-L
        // ones, so condition (a) alone let them run straight into it.
        // minimumPeriodToLightCrossingRatio=150 keeps a factor ~2 safety
        // margin over the observed failure range's own upper end.
        const double periapsisSeparation=(measurementState.firstPosition
            -measurementState.secondPosition).norm();
        const double lightCrossingTime=periapsisSeparation/c;
        const double periodToLightCrossingRatio=(lightCrossingTime>0.0)
            ?period/lightCrossingTime:std::numeric_limits<double>::infinity();
        constexpr double minimumPeriodToLightCrossingRatio=150.0;
        if(periapsis<=comptonBarrierRadius
           ||periodToLightCrossingRatio<=minimumPeriodToLightCrossingRatio) {
            result.lifetimeSeconds=simulatedTimeTotal;
            result.meanRadiatedPowerWatts=simulatedTimeTotal>0.0
                ?radiatedEnergyTotal/simulatedTimeTotal
                :std::numeric_limits<double>::quiet_NaN();
            result.lifetimeSecondsLab=labFrameTimeTotal;
            result.meanRadiatedPowerWattsLab=labFrameTimeTotal>0.0
                ?radiatedEnergyTotal/labFrameTimeTotal
                :std::numeric_limits<double>::quiet_NaN();
            result.calibrationOutcome=SimulationOutcome::ReachedCutoff;
            result.calibrationSeconds=simulatedTimeTotal;
            result.calibrationSecondsLab=labFrameTimeTotal;
            // Annihilation from the state the trajectory actually reached,
            // rather than from a pair at rest.  What is left to share is the
            // invariant W = (m1+m2)c^2 + mu*epsilon, below the rest-mass sum
            // by the binding this orbit accumulated on its way down.
            result.terminalSemiMajorAxis=
                -attractionParameter/(2.0*elements.specificEnergy);
            result.terminalBindingEnergy=
                reducedMass*std::abs(elements.specificEnergy);
            // DIPOLE-DIPOLE, which the Kepler element does not carry.
            //
            // elements.specificEnergy is an OSCULATING KEPLER element, not a
            // general energy: eight places read it as a = -K/(2 eps),
            // e^2 = 1 + 2 eps l^2/K^2 and the period.  Folding a 1/r^3 term
            // into it would feed every one of those a quantity they were not
            // derived for, and the damage would peak exactly here, at the
            // terminal radius, where the dipole sector is 56% of Coulomb.
            // So it stays a Kepler element and the dipole energy is added
            // where it physically belongs instead: to the invariant the
            // final-state photons actually share.
            //
            // Nothing about the DYNAMICS is being corrected -- the
            // trajectory already feels this force (mutualForces carries
            // regularizedDipoleForce) and the mechanical ledger already
            // carries the energy (conservativeParticleEnergy adds
            // dipolePotential).  Only the secular bookkeeping's energy LABEL
            // omitted it.
            //
            // WHAT THIS DOES AND DOES NOT DO.  An earlier version of this
            // comment claimed the term makes W channel-dependent, since the
            // angular factor changes sign between aligned and anti-aligned
            // moments and that alignment is what separates para from ortho
            // here.  That is wrong, and the reason is worth keeping.
            //
            // What is averaged is the orientation of the RIGID pair relative
            // to L, at FIXED mutual angle cos = mu1^.mu2^.  That angle is
            // what separates the channels and is not averaged over; the
            // identity below was checked separately at cos = -0.9, -0.4, 0,
            // +0.4 and +0.9, so the vanishing happens INSIDE each class.
            //
            // The isotropy that assumes is not merely an initial condition:
            // it survives the dynamics.  Opposite charges give opposite q/m,
            // which suggests the perpendicular components counter-rotate and
            // wash the mutual angle out -- measured, they do not.  cos runs
            // 0.885 -> 0.893, 0.552 -> 0.552, -0.765 -> -0.765 from seed to
            // termination, the largest excursion seen being 0.15.  Opposite
            // charges reverse the spin-moment relation along WITH q/m, so
            // both moments precess the same way about a common B_BMT.
            //
            // WHY THEY DO NOT ALIGN, corrected.  An earlier version argued
            // that B_BMT is dominated by the motional v x E/c^2, which lies
            // along L, so precession conserves mu.L.  That holds only at
            // large separation: the partner's own dipole field is 0.25 of
            // the motional one at a_Ps, equals it at 6630 fm, and is 5.07x
            // it at the terminal radius -- so exactly where this term
            // matters, the precession is about the partner's field, not
            // about L.  The argument failed where it was needed.
            //
            // Nor is "there is no time" the answer: the dipole precession
            // period at a_Ps is 45.6 ps against a collapse of about 133 ps,
            // some three full turns.
            //
            // Nor is "there is nothing to lose energy to" the answer, which
            // is what an earlier pass at this comment said.  The channel
            // EXISTS: M1 enters the quantum beside E1 and E2 (see
            // quantizedPower in crem_trajectory.hpp), and integrated over
            // the collapse it carries 3.98 eV against 3777 eV of E1 -- a
            // share of 1.05e-3, small but not zero.
            //
            // What settles it is the INTEGRATED rotation.  The M1 reaction
            // torque, the same one the integrator applies as
            // mu += gamma*T*dt, turns the moment at 0.85 rad/s at a_Ps and
            // at 1.0e16 rad/s at the terminal radius -- 2.1e-12 and 1.2e-3
            // of the BMT precession rate respectively, so at the bottom it
            // is emphatically not weak, and the alignment time there,
            // 9.9e-17 s, is fourteen orders below the collapse.  But the
            // pair does not LIVE there.  And the accumulation is not even
            // secular: the torque points along mu x m'''_total, and
            // m'''_total turns with the orbit, so successive parts of a
            // period cancel and what survives is an OSCILLATION of amplitude
            // ~Omega_reaction/omega_orb, not a drift ~Omega_reaction*t.
            //
            // Measured directly, which is the only way to see that: a
            // stochastic run and a disabled run from the same seed differ
            // ONLY by this torque (until a photon fires the trajectories are
            // bit-identical), so |d(mu-hat)| between them isolates it.  At
            // a_Ps/200, over 1/2/4/8/16 orbits: 2.16, 2.14, 2.17, 1.84,
            // 2.72e-6 -- bounded, no linear trend across a sixteenfold span.
            // At a_Ps/50: 3.02e-9 at two orbits, 2.77e-9 at eight.
            //
            // So the omitted reorientation is ~Omega_reaction/omega_orb:
            // 4e-17 rad at a_Ps, 2e-6 at a_Ps/200, and at worst ~6e-5 rad
            // (0.003 degrees) at the terminal radius.  An earlier pass here
            // reported 2.33e-2 rad (1.33 degrees) from integrating
            // Omega_reaction dt over the inspiral; that integral sums the
            // MAGNITUDE of the instantaneous rate and so presumes the drift
            // the measurement above rules out.  It is withdrawn.
            //
            // (This supersedes the 1.14e-44 ratio recorded in the README's
            // Sonda 6, which evaluated the torque on an unfilled history and
            // is wrong by some thirty-two orders of magnitude.  It also
            // supersedes an earlier 0.30 degrees written here, which was one
            // unrepresentative orientation rather than an average.)
            //
            // THE GAP, taken apart, and since closed.  The quantized mode
            // used to gate two things at one flag: the reaction torque and
            // the dipoleConstraintEnergy drain.  Only the first was a real
            // omission, and electrodynamics.hpp now splits the two -- the
            // drain stays off (the photon carries that energy), the torque
            // runs.  What follows is the measurement that motivated it.
            //
            //   Energy.  There was a real defect here, since fixed in
            //   electrodynamics.hpp: both M1 sinks were gated for
            //   stochasticElectricDipole only, while the charge sector has
            //   always excluded disabled as well -- and disabled is the model
            //   this file integrates as the BACKGROUND.  The background
            //   therefore drained the same reservoir as the run it was
            //   subtracted from (4.0e-5 of the M1 energy survived), and in
            //   the quantized mode it carried a sink the real run never had,
            //   so between photons run and background were not bit-identical
            //   when they should have been.  Gating disabled too fixes both.
            //
            //   It moves no numbers, and the reason is worth stating: the
            //   reservoir never fed the collapse anyway.  measuredDelta's
            //   specificEnergy is computed but deliberately unread (see the
            //   "ANGULAR MOMENTUM only" comment below); the secular loss is
            //   deltaEnergyPerOrbit = -orbitalRadiatedEnergy/reducedMass, and
            //   orbitalRadiatedEnergy is the flux with M1 already subtracted
            //   out; and in the production path even that is bypassed, since
            //   lossPerOrbit takes expectedLossPerOrbit, the analytic E1
            //   Larmor power.  M1's one entry into the dynamics is the hazard
            //   inside the single measured orbit per checkpoint, against an
            //   analytic skip of up to 200000 E1-only orbits -- roughly 1e-8
            //   of the collapse once the 1.9e-3 M1 share is weighted by the
            //   ~5e-6 of the hazard that measured orbit carries.
            //
            //   Angular momentum.  M1's share of the radiated angular
            //   momentum flux runs 2.2e-15 at a_Ps to 9.3e-5 at the terminal
            //   radius -- negligible once weighted by dwell time.
            //
            //   Affordability.  Cumulative M1 never exceeds the dipole
            //   orientation energy available at the same radius; the ratio
            //   tops out at 3.9e-3 (3.70 eV against 950 eV of U_dd at the
            //   barrier).  Neither account is overdrawn.
            //
            // So the omission was the orientation back-reaction alone, worth
            // at most 0.003 degrees, and it is the one now restored: the
            // supplies linear recoil only and leaves the radiating moment
            // unturned, so nothing was being double counted by gating it --
            // the gate simply deleted it.  The rotation's own energetic cost
            // is the dipole-dipole term conservativeParticleEnergy already
            // carries, so restoring it does not double count either.
            //
            // It did not separate the channels: the
            // destructive M1 interference is exact only at cos = -1, and
            // averaged over orientations the M1 share of the quantized
            // stream is 1.875e-3 +- 3.7e-4 for para against 1.903e-3 +-
            // 2.3e-4 for ortho -- a difference consistent with zero.
            //
            // Measured directly over 144 components from 72 trajectories,
            // mu.L has SD 0.601 at the seed and 0.596 at termination against
            // 1/sqrt(3) = 0.577 for an isotropic distribution, with a median
            // excursion of 0.019 across the whole collapse.  The ungated
            // torque adds at most ~6e-5 to that -- three hundred times less,
            // so the isotropy verdict has three orders of margin.  (An
            // earlier pass put that addition at 2.33e-2 and withdrew the
            // claim of a comfortable margin on the strength of it; the
            // withdrawal rested on the discarded integral and is itself
            // withdrawn.)
            //
            // And had there been dissipation it would not have produced
            // repulsion in ortho anyway: the system would seek the MINIMUM
            // of the dipole energy, whose configuration is attractive in
            // both channels.
            //
            // (Which channel carries the ALIGNED moments is the opposite of
            // the naive guess, and README's own Sonda 4 says so: para is
            // forced to cos >= 0.5 and ortho to cos < 0.5, because opposite
            // charges invert the spin-moment correlation.  Ortho's condition
            // is one-sided -- everything but strongly aligned -- not
            // "anti-aligned".)
            //
            // After the azimuth average the factor depends only on mu1.mu2
            // and (mu1.L)(mu2.L).  For isotropically drawn moments at fixed
            // cos = mu1^.mu2^, <(mu1^.L)(mu2^.L)> = cos/3 -- verified over
            // 2e6 samples -- so
            //
            //     <factor> = -cos/2 + (3/2)(cos/3) = 0
            //
            // IDENTICALLY, for every cos, hence for either channel.  The
            // term therefore has zero expectation and cannot separate para
            // from ortho at all.  Measured over four seeds at fifteen
            // trajectories: para -0.061 +/- 0.219 keV, ortho +0.157 +/-
            // 0.147 keV, difference -0.217 +/- 0.263 keV, every one
            // consistent with zero and with mixed signs inside each channel.
            //
            // What the term does contribute is SCATTER, and a large one: at
            // the terminal radius the dipole sector is 56% of Coulomb, so
            // individual trajectories carry order +/-1.5 keV of it.  The
            // prediction it supports is a BROADENED line, not a shifted one
            // -- roughly 0.15% of 511 keV in width, against the 0.27% by
            // which the binding shifts the centre.
            // Evaluated through the azimuth average, NOT by building a
            // periapsis State: that helper resets the orbital plane to
            // canonical x-y and puts the separation along +x, which is fine
            // for a collapse-time estimate and wrong here, where the answer
            // depends on the direction of n relative to the moments.  See
            // azimuthAveragedDipoleEnergy for the measurement that shows the
            // convention can flip the sign.
            result.terminalDipoleEnergy=azimuthAveragedDipoleEnergy(
                regularizedPeriapsis(elements,attractionParameter,
                                     separationFloor()),
                firstDipole,secondDipole,angularMomentumDirection);
            result.annihilationInvariantEnergy=
                (firstMass+secondMass)*c*c-result.terminalBindingEnergy
                +result.terminalDipoleEnergy;
            result.annihilationPhotonEnergies=annihilationPhotonEnergiesFor(
                result.annihilationInvariantEnergy,
                selectedPhenomenon==1,stochasticSkipStream);
            if(std::getenv("CREM_DEBUG"))
                std::cerr<<std::setprecision(12)
                         <<"  FINAL (periapsis/light-crossing cutoff): t_S'="
                         <<simulatedTimeTotal*1e12<<"ps t_lab="
                         <<labFrameTimeTotal*1e12<<"ps ratio="
                         <<(simulatedTimeTotal>0.0
                             ?labFrameTimeTotal/simulatedTimeTotal:1.0)
                         <<" P_S'="<<result.meanRadiatedPowerWatts
                         <<"W P_lab="<<result.meanRadiatedPowerWattsLab
                         <<"W"<<std::setprecision(6)<<std::endl;
            return result;
        }
        SimulationOptions measureOptions;
        measureOptions.collectFrames=false;
        measureOptions.observationTime=period;
        measureOptions.terminalSeparation=comptonBarrierRadius;
        // The ENERGY side of the secular update reads
        // finalState.orbitalRadiatedEnergy now, not a
        // conservativeParticleEnergy() difference (see below), so the flux
        // quadrature this enables is no longer pure overhead for run -- it
        // is the measurement.  background still only needs the mechanical
        // (position, velocity) endpoint for its angular-momentum
        // subtraction, so its own copy below turns this back off.
        measureOptions.radiatedEnergyBookkeeping=true;
        // The budget is shared with the whole estimate: a single
        // measurement orbit can itself run into the stiff region near the
        // boundary, so it must be interruptible too.
        measureOptions.stopRequested=[&]() {
            return wallClockSpent()>wallClockBudgetSeconds;
        };
        // nuclearCutoff, not comptonBarrierRadius: measurementState's own
        // position is regularizedPeriapsis() (via osculatingPeriapsisState),
        // which is allowed to sit AT or below comptonBarrierRadius once the
        // orbit has circularized deep enough (see the r_min branch there) --
        // passing comptonBarrierRadius here as well would then make this
        // measurement start already past its own trajectoryCutoff, a
        // zero-step degenerate case (measured directly: startSep=189.5fm <
        // cutoff=193.3fm, elapsed=0ps, NumericalFailure).  nuclearCutoff is
        // the universal, pair-independent "genuine point-particle floor"
        // already used everywhere else in the codebase for exactly this
        // role (e.g. lienardWiechertField's "whichever floor is bigger
        // wins"), and it is always strictly below comptonBarrierRadius, so
        // this measurement now always starts with room to run regardless of
        // how deep regularizedPeriapsis has to go.  In the far-from-barrier
        // regime this changes nothing observable: the measured orbit's
        // deepest point stays far above either cutoff there, so which one
        // is configured is moot.
        const MechanicalTrajectoryResult run=runMechanicalTrajectory(
            measurementState,period,nuclearCutoff,measureOptions,
            activeReactionModel);

        // Angular momentum and energy are both read from the measured
        // orbit's actual start/end state rather than from a formula: L is
        // exactly conserved by a central force regardless of its radial
        // profile, and conservativeParticleEnergy() is the same quantity
        // positronium_validation checks to 1e-12, so together they isolate
        // the genuine secular drift with no reaction-model-dependent
        // bookkeeping in the way.  (radiatedEnergy / radiatedAngularMomentum
        // were tried and rejected: they are system-total EM flux bookkeeping
        // that bundles in the M1 magnetic-dipole/spin-precession channel,
        // which does not recoil the orbit at all, and dominated the total by
        // orders of magnitude here.)
        //
        // Even so, restarting each measurement from a fresh, from-rest
        // "periapsis" State leaves a small settling artefact in the
        // reconstructed retarded-field history (causalInitialHistory has no
        // real orbit behind it to extrapolate from) that reads as a spurious
        // secular loss growing linearly with elapsed time -- present even
        // with --radiation-reaction disabled, where the mechanics guarantee
        // exact conservation.  A second, otherwise identical measurement run
        // with the reaction force explicitly forced off calibrates exactly
        // that artefact away by subtraction, leaving only the effect the
        // active --radiation-reaction model actually adds.  It does not
        // touch gRadiationReactionModel, so it stays safe across the worker
        // threads in runCremCollapseExperiment.
        //
        // The background is the second of the two full-orbit integrations
        // every checkpoint used to pay for, i.e. half the cost of the whole
        // estimate.  Measured across a complete trajectory it is a large but
        // SMOOTH fraction of the raw signal: background/raw runs 0.42..0.74
        // and drifts by a median 0.52% per checkpoint.  What is cached here is
        // therefore that RATIO rather than the absolute artefact, which grows
        // with the orbit and would go stale far faster.  The ratio does jump
        // by up to 12% where the orbit changes character, so the refresh is
        // driven by how far the osculating energy has moved since the cached
        // ratio was taken, not by a fixed checkpoint count: near the boundary,
        // where the jumps live, it re-measures more often on its own.
        const bool refreshBackground = !haveBackgroundRatio
            || !(std::abs(elements.specificEnergy-energyAtLastBackground)
                 <= backgroundRefreshFraction*std::abs(energyAtLastBackground));
        MechanicalTrajectoryResult background;
        if(refreshBackground) {
            // Same cutoff as run above, and for the same reason: it must
            // match run's so the pair of integrations stays comparable (the
            // background subtraction assumes both runs cover the same
            // stretch), and it must clear measurementState's own position.
            // Its own copy of measureOptions with radiatedEnergyBookkeeping
            // turned back off: background is only ever read for its
            // mechanical endpoint (angular momentum), never its
            // orbitalRadiatedEnergy, so the flux quadrature run buys the
            // measurement would be pure overhead here.
            SimulationOptions backgroundOptions=measureOptions;
            backgroundOptions.radiatedEnergyBookkeeping=false;
            background=runMechanicalTrajectory(
                measurementState,period,nuclearCutoff,backgroundOptions,
                ChargeRadiationReactionModel::disabled);
        }
        const auto measuredDelta=[&](const State& end) {
            const Vec3 relPos=end.firstPosition-end.secondPosition;
            const Vec3 relVel=end.firstVelocity-end.secondVelocity;
            return OsculatingElements{
                (conservativeParticleEnergy(end)
                    -conservativeParticleEnergy(measurementState))/reducedMass,
                cross(relPos,relVel).z-elements.specificAngularMomentum};
        };

        // Either the freshly measured artefact, or the cached ratio applied to
        // this checkpoint's own raw measurement.
        const auto backgroundFor=[&](const OsculatingElements& realDelta) {
            if(refreshBackground&&isFinite(background.finalState))
                return measuredDelta(background.finalState);
            if(haveBackgroundRatio)
                return OsculatingElements{
                    backgroundEnergyRatio*realDelta.specificEnergy,
                    backgroundAngularRatio*realDelta.specificAngularMomentum};
            return OsculatingElements{};
        };

        if(run.outcome==SimulationOutcome::ReachedCutoff) {
            result.lifetimeSeconds=simulatedTimeTotal+run.elapsedTime;
            // Already in real (not specific) Joules, and already exactly the
            // E1 energy this partial orbit radiated -- no background
            // subtraction needed (see the full explanation where this same
            // substitution is made below, for the ordinary ObservationLimit
            // path this ReachedCutoff path is the truncated-orbit sibling of).
            radiatedEnergyTotal+=run.finalState.orbitalRadiatedEnergy;
            result.meanRadiatedPowerWatts=result.lifetimeSeconds>0.0
                ?radiatedEnergyTotal/result.lifetimeSeconds
                :std::numeric_limits<double>::quiet_NaN();
            // This partial measurement orbit fires no photon (the recoil
            // while-loop lives in the skip/jump branch below, not here), so
            // centreOfMassVelocity -- and therefore beta/gamma -- is exactly
            // constant across the whole of run.elapsedTime.
            result.lifetimeSecondsLab=labFrameTimeTotal
                +gammaFromBeta(centreOfMassVelocity.norm()/c)*run.elapsedTime;
            result.meanRadiatedPowerWattsLab=result.lifetimeSecondsLab>0.0
                ?radiatedEnergyTotal/result.lifetimeSecondsLab
                :std::numeric_limits<double>::quiet_NaN();
            if(std::getenv("CREM_DEBUG"))
                std::cerr<<"  FINAL (mid-orbit cutoff): t_S'="
                         <<result.lifetimeSeconds*1e12<<"ps t_lab="
                         <<result.lifetimeSecondsLab*1e12<<"ps ratio="
                         <<(result.lifetimeSeconds>0.0
                             ?result.lifetimeSecondsLab
                                 /result.lifetimeSeconds:1.0)
                         <<" P_S'="<<result.meanRadiatedPowerWatts
                         <<"W P_lab="<<result.meanRadiatedPowerWattsLab
                         <<"W"<<std::endl;
            result.calibrationOutcome=SimulationOutcome::ReachedCutoff;
            result.calibrationSeconds=result.lifetimeSeconds;
            result.calibrationSecondsLab=result.lifetimeSecondsLab;
            // The measurement orbit ran into the boundary partway round, so
            // only the fraction of a revolution it actually covered counts.
            result.revolutions=revolutionsTotal
                +(period>0.0?run.elapsedTime/period:0.0);
            return result;
        }
        if(run.outcome!=SimulationOutcome::ObservationLimit
           ||!isFinite(run.finalState)
           ||(refreshBackground&&!isFinite(background.finalState))) {
            if(std::getenv("CREM_DEBUG"))
                std::cerr<<"  DIAG mechanical-measurement trip: run.outcome="
                         <<static_cast<int>(run.outcome)<<" finiteFinal="
                         <<isFinite(run.finalState)<<" specificEnergy="
                         <<elements.specificEnergy<<" specificAngularMomentum="
                         <<elements.specificAngularMomentum<<std::endl;
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal+run.elapsedTime;
            // No photon fires during this single measurement orbit (the
            // recoil while-loop lives in the skip/jump branch below), so
            // beta is exactly centreOfMassVelocity's current, constant value
            // across the whole of run.elapsedTime -- same reasoning as the
            // ReachedCutoff sibling of this exit, above.
            result.calibrationSecondsLab=labFrameTimeTotal
                +gammaFromBeta(centreOfMassVelocity.norm()/c)*run.elapsedTime;
            return result;
        }
        if(wallClockSpent()>wallClockBudgetSeconds) {
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal+run.elapsedTime;
            result.calibrationSecondsLab=labFrameTimeTotal
                +gammaFromBeta(centreOfMassVelocity.norm()/c)*run.elapsedTime;
            return result;
        }

        const double measuredElapsed=std::max(run.elapsedTime,1.0e-30);
        // ANGULAR MOMENTUM only, from here on: realDelta/backgroundDelta's
        // specificEnergy field is computed (conservativeParticleEnergy(), a
        // Darwin/near-field approximation) but deliberately not read for the
        // loss below any more -- see deltaEnergyPerOrbit's own comment.
        // backgroundEnergyRatio is still computed in the cache refresh
        // alongside backgroundAngularRatio (the two are cached together),
        // but nothing downstream reads it either now; it is vestigial,
        // left rather than split the cache apart, and costs nothing to carry.
        const OsculatingElements realDelta=measuredDelta(run.finalState);
        const OsculatingElements backgroundDelta=backgroundFor(realDelta);
        // Refresh the cache only from an actual measurement, and only when the
        // raw values are large enough for the quotient to mean anything.
        if(refreshBackground&&isFinite(background.finalState)) {
            if(realDelta.specificEnergy!=0.0
               &&realDelta.specificAngularMomentum!=0.0) {
                backgroundEnergyRatio=backgroundDelta.specificEnergy
                    /realDelta.specificEnergy;
                backgroundAngularRatio=backgroundDelta.specificAngularMomentum
                    /realDelta.specificAngularMomentum;
                haveBackgroundRatio=true;
                energyAtLastBackground=elements.specificEnergy;
            }
        }
        // orbitalRadiatedEnergy, not a conservativeParticleEnergy()
        // difference: conservativeParticleEnergy() is the Darwin (static
        // near-field) approximation, valid only while the orbital period
        // stays many light-crossing times long, and measured directly to
        // fail exactly that way -- on seed 4, as period/(separation/c) fell
        // from ~950 to ~40, the conservativeParticleEnergy-based dE/E grew
        // from 3.6e-6 (matching the documented ~1e-5 integrator noise floor)
        // to OVER 100% per "orbit" (checkpoint 24: 1.73), while
        // orbitalRadiatedEnergy over the SAME stretch stayed smooth and
        // monotonic throughout (0.042 -> 0.272), only faltering mildly at
        // the very end.  It is computed from the exact retarded far-zone
        // Poynting flux (electromagneticFieldFluxRates via
        // particleMultipoleRadiation), with the M1 (magnetic-dipole/spin)
        // channel already excluded (see orbitalRadiatedEnergy's own comment
        // in state.hpp) -- and unlike conservativeParticleEnergy(measure-
        // mentState), which is a near-zero DIFFERENCE of two large numbers,
        // orbitalRadiatedEnergy starts at exactly 0 for a fresh
        // measurementState and accumulates the genuine physical loss
        // directly, so no background subtraction is needed for it: the
        // "settling artefact" this whole background machinery exists to
        // remove is specific to conservativeParticleEnergy's own near-field
        // bootstrap, not to a directly-measured flux.  Sign matches the
        // existing convention (negative = orbit binding tighter).
        const double deltaEnergyPerOrbit=
            -run.finalState.orbitalRadiatedEnergy/reducedMass;
        const double deltaAngularMomentumPerOrbit=
            realDelta.specificAngularMomentum-backgroundDelta.specificAngularMomentum;
        // Read once, here, rather than at its previous location further
        // down: this guard needs it too, see below.
        const bool isStochastic=activeReactionModel
            ==ChargeRadiationReactionModel::stochasticElectricDipole;
        // A secular inspiral cannot shed a large fraction of its own binding
        // energy in ONE orbit; if it could, orbit averaging would not apply in
        // the first place.  Without this guard a single bad measurement is
        // folded straight into the osculating elements and destroys them.
        // Observed with --radiation-reaction coherent: the reaction force is
        // built from a third derivative of the dipole moment divided by
        // 2*h^3, and when the retarded-history spacing degenerates the stencil
        // returns garbage -- dE/E jumped 22 orders of magnitude in a single
        // checkpoint (2.7e-5 -> 1.4e17) while the reaction-free background
        // stayed at 1.9e-5.  The specific energy then sat at -1.2e30, every
        // later measurement read exactly zero, and the run was reported as
        // "not decaying", which named a downstream symptom rather than the
        // numerical failure that had actually occurred.
        //
        // The magnitude half of this guard is gated OFF for isStochastic,
        // checked rather than assumed: traced one actual trip (seed 107,
        // trajectory index 8 of a seed-99 batch) and found deltaEnergyPerOrbit
        // itself was NOT garbage -- a spin kick (electrodynamics.hpp's own
        // comment on stochasticElectricDipole) had pushed the orbit to
        // e^2~0.945, and the following single measured orbit genuinely
        // radiated 56% of the binding energy in one periapsis passage, a
        // real, finite number from the flux integral
        // (electromagneticFieldFluxRates via particleMultipoleRadiation),
        // not the fragile third-derivative stencil the failure this guard
        // was built for came from -- stochasticElectricDipole never even
        // evaluates that stencil (particleMultipoleRadiation applies no
        // continuous reaction force for it).  More directly: for
        // isStochastic, deltaEnergyPerOrbit is used for NOTHING but this
        // guard and the Larmor-ratio diagnostic just below (itself
        // separately finite-guarded) -- expectedLossPerOrbit, the
        // Larmor-orbit-averaged rate at the CURRENT osculating elements,
        // sizes orbitsToSkip instead (see its own comment: this measured
        // value is "noise-dominated and unusable for sizing the skip" in
        // the USUAL near-circular regime, and simply irrelevant, not just
        // noisy, in this one) -- so a large-but-finite value here cannot
        // propagate into an extrapolated skip the way the guard's own
        // justification describes.  isfinite is still checked for every
        // model: that half catches genuine corruption regardless of cause.
        constexpr double maxRelativeLossPerOrbit=0.5;
        if(!std::isfinite(deltaEnergyPerOrbit)
           ||(!isStochastic&&std::abs(deltaEnergyPerOrbit)
               >maxRelativeLossPerOrbit*std::abs(elements.specificEnergy))) {
            if(std::getenv("CREM_DEBUG"))
                std::cerr<<"  DIAG maxRelativeLossPerOrbit trip: deltaEnergyPerOrbit="
                         <<deltaEnergyPerOrbit<<" specificEnergy="
                         <<elements.specificEnergy<<" ratio="
                         <<deltaEnergyPerOrbit/elements.specificEnergy
                         <<" specificAngularMomentum="
                         <<elements.specificAngularMomentum<<std::endl;
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal+measuredElapsed;
            result.calibrationSecondsLab=labFrameTimeTotal
                +gammaFromBeta(centreOfMassVelocity.norm()/c)*measuredElapsed;
            return result;
        }
        if(const char* debug=std::getenv("CREM_DEBUG");debug) {
            std::cerr<<"  measured dE/E="
                     <<deltaEnergyPerOrbit/elements.specificEnergy
                     <<" (orbitalRadiatedEnergy="
                     <<run.finalState.orbitalRadiatedEnergy
                     <<"J) dL/L="<<deltaAngularMomentumPerOrbit
                        /elements.specificAngularMomentum<<std::endl;
        }
        // CREDIT THE MEASURED ORBIT ITSELF, for isStochastic only.  This one
        // orbit was genuinely, mechanically integrated, and genuinely
        // radiated deltaEnergyPerOrbit via the real retarded-field flux
        // (orbitalRadiatedEnergy) -- that is Maxwell's equations acting on
        // the true trajectory, true regardless of which reaction model is
        // switched on.  For the deterministic branch this measurement seeds
        // the whole envelope extrapolation below (energyGrowth, hence
        // updatedEnergyMagnitude, is built directly from lossPerOrbit=
        // |deltaEnergyPerOrbit|), so it is already accounted for there.  For
        // isStochastic it was not: elements.specificEnergy there is touched
        // ONLY inside the photon while-loop below, and that loop's hazard
        // integral is explicitly scoped to the orbitsToSkip orbits AFTER
        // this one (see "n_skip" in its own derivation) -- it never covered
        // this orbit's own, already-measured loss, which was simply
        // discarded.  Measured directly, not assumed to be negligible: on a
        // shallow trajectory (seed 42, 37 checkpoints) the total discarded
        // loss came to ~2e-5 of the energy actually credited via photons
        // over the whole run -- but on one that spent time at high
        // eccentricity (seed 107, the trajectory investigated for the
        // guard fix above) it reached 38.5%, a real violation of energy
        // conservation, not a rounding error, and one this guard fix makes
        // easier to reach (a trajectory that would previously have failed
        // outright can now run on and accumulate it).  Applied here, once
        // per checkpoint, additively with (not instead of) the photon
        // credits below: deltaEnergyPerOrbit already carries this loss's
        // sign by convention (see its own comment), so += is correct, not
        // -=, and matches the ReachedCutoff branch's use of the same
        // orbitalRadiatedEnergy field elsewhere in this function.
        if(isStochastic) {
            elements.specificEnergy+=deltaEnergyPerOrbit;
            radiatedEnergyTotal+=run.finalState.orbitalRadiatedEnergy;
        }
        // Measured orbital dissipation of this osculating orbit against the
        // Larmor rate for the same orbit.  Both sides are orbit averages over
        // the same period, so the ratio is a direct, parameter-free check of
        // the engine's radiation sector against textbook electrodynamics.
        {
            const double measuredPower=
                -deltaEnergyPerOrbit*reducedMass/measuredElapsed;
            const double semiMajorAxis=
                -attractionParameter/(2.0*elements.specificEnergy);
            const double eccentricity=std::sqrt(std::max(0.0,1.0
                +2.0*elements.specificEnergy*elements.specificAngularMomentum
                    *elements.specificAngularMomentum
                    /(attractionParameter*attractionParameter)));
            const double larmorPower=
                larmorOrbitAveragedPower(semiMajorAxis,eccentricity);
            if(std::isfinite(measuredPower)&&std::isfinite(larmorPower)
               &&larmorPower>0.0&&measuredPower>0.0) {
                larmorRatioSum+=measuredPower/larmorPower;
                ++larmorRatioCount;
            }
        }
        // Closed-form Larmor-orbit-averaged loss rate for the CURRENT
        // osculating orbit -- the same quantity larmorOrbitAveragedPower()
        // just above already computed for the Larmor-ratio diagnostic, not
        // a new formula.  Needed only by stochasticElectricDipole, whose
        // single measured orbit almost never contains a real photon (see
        // this reaction model's own comment in electrodynamics.hpp: hazard
        // per SINGLE orbit near a=3pm measures ~5.5e-5), so
        // deltaEnergyPerOrbit from that one orbit is noise-dominated and
        // unusable for sizing the skip -- exactly the envelope the photon
        // hazard integral below is itself derived against, so this is the
        // right substitute, not an approximation invented for convenience.
        // (isStochastic itself is declared above, by the guard that needs
        // it first.)
        const double expectedLossPerOrbit=isStochastic
            ?larmorOrbitAveragedPower(
                  -attractionParameter/(2.0*elements.specificEnergy),
                  std::sqrt(std::max(0.0,1.0+2.0*elements.specificEnergy
                      *elements.specificAngularMomentum
                      *elements.specificAngularMomentum
                      /(attractionParameter*attractionParameter))))
                 *period/reducedMass
            :0.0;
        if(isStochastic
           ?!(expectedLossPerOrbit>0.0)
           :(!(deltaEnergyPerOrbit<0.0)||!std::isfinite(deltaEnergyPerOrbit))) {
            // No measurable energy loss this orbit: with the reaction force
            // disabled (or below numerical noise), the orbit is not
            // secularly decaying, so there is nothing further to observe.
            result.calibrationOutcome=SimulationOutcome::ObservationLimit;
            result.calibrationSeconds=simulatedTimeTotal+measuredElapsed;
            result.calibrationSecondsLab=labFrameTimeTotal
                +gammaFromBeta(centreOfMassVelocity.norm()/c)*measuredElapsed;
            result.secularLossAbsent=true;
            return result;
        }

        // --- Skipped orbits, integrated with the known functional form ---
        //
        // Freezing dE/dorbit across the jump is a zeroth-order hold, and it is
        // wrong in a known direction: the orbit tightens, so the true loss per
        // orbit GROWS and the hold always underestimates it.  That is why the
        // jump had to stay within 3% of |E|.
        //
        // With u = -E and dipole radiation, P ~ a^-4 and T ~ a^(3/2) give
        // du/dn = B u^(5/2), whose exact solution is
        //
        //     u(n) = u0 (1 - s)^(-2/3),      s = (3/2) n du0/u0,
        //
        // and expanding it to first order in s returns exactly the frozen-rate
        // answer, so this is a resummation of the same measurement rather than
        // a different model.  B is still taken from THIS checkpoint's measured
        // loss; only the shape of the interpolation between measurements
        // changes.  The period follows T(n) = T0 (1 - s), so the elapsed time
        // integrates to T0 n (1 - s/2) instead of T0 n.
        //
        // Checked against the engine on a full trajectory: fitting
        // log(dE/E) against log(u) gives slope 1.4617 (1.4465 over the first
        // three quarters) where this form predicts 3/2, so the exponent is
        // right to within 3%.  The angular momentum is carried by the measured
        // ratio k = (dL/L)/(du/u), which came out between -0.454 and -0.490
        // over that same trajectory, i.e. far steadier than dL/dorbit itself;
        // L then follows L = L0 (u/u0)^k.
        const double energyMagnitude=std::abs(elements.specificEnergy);
        const double lossPerOrbit=isStochastic
            ?expectedLossPerOrbit:std::abs(deltaEnergyPerOrbit);
        int orbitsToSkip=1;
        if(lossPerOrbit>0.0&&energyMagnitude>0.0) {
            const double requestedOrbits=maximumJumpParameter*energyMagnitude
                /(1.5*lossPerOrbit);
            // Saturate while the value is still floating-point.  For a tiny
            // measured loss requestedOrbits can exceed INT_MAX or become +inf;
            // converting either value to int before clamp has undefined
            // behaviour.  The bounded value is always representable because
            // maxOrbitsSkippedAtOnce itself is an int.
            const double boundedOrbits=std::isfinite(requestedOrbits)
                ?std::clamp(requestedOrbits,1.0,
                    static_cast<double>(maxOrbitsSkippedAtOnce))
                :static_cast<double>(maxOrbitsSkippedAtOnce);
            orbitsToSkip=static_cast<int>(boundedOrbits);
        }
        const double jumpParameter=std::min(
            1.5*static_cast<double>(orbitsToSkip)*lossPerOrbit/energyMagnitude,
            maximumJumpParameter);
        const double energyGrowth=std::pow(1.0-jumpParameter,-2.0/3.0);
        const double updatedEnergyMagnitude=energyMagnitude*energyGrowth;
        // k, in CLOSED FORM rather than measured.  The previous approach
        // (k = (dL/L)/(dE/E) from this checkpoint's own measurement) tracked
        // smoothly while the orbit was still visibly eccentric but became
        // noise-dominated as e -> 0 (dL on one orbit shrinks toward the scale
        // of the background-subtraction artefact), and a bad k there gets
        // raised to a power and multiplied into L across tens of thousands of
        // skipped orbits, so the error compounds instead of averaging out:
        // on seed 42 the discriminant 1+2EL^2/mu^2 (=e^2), which can only
        // fall as the orbit radiates, instead climbed back from 0.0005 to
        // 0.15 over checkpoints 13-23 and went negative (unphysical) by
        // checkpoint 28.  EMA-smoothing the measurement was tried and made
        // it WORSE (negative by checkpoint 4), which is itself informative:
        // the drift is not noise around a fixed point, so damping it with a
        // stale average fights the real k, it does not filter a fake one.
        //
        // The real k is not empirical at all.  Orbit-averaging the dipole
        // reaction force -- a_react = gamma*r'''(t), Abraham-Lorentz on the
        // relative coordinate, gamma folding in q_eff^2/(6 pi eps0 mu c^3) --
        // over an unperturbed Kepler ellipse (same method as Peters 1964, but
        // for the dipole's P~r^-4 in place of the quadrupole's r^-6) gives,
        // with r'''=-(attractionParameter/r^3)[v-3(r-hat.v)r-hat]:
        //
        //   <dE/dt> = -gamma*attractionParameter^2*(1+e^2/2)
        //                 / (a^4 (1-e^2)^(5/2))
        //   <dL/dt> = -gamma*attractionParameter*L / (a^3 (1-e^2)^(3/2))
        //
        // using the standard time-averages <r^-n> = a^-n I_(n-1), I_n =
        // (1/2pi) integral (1-e cosE)^-n dE, generated by the recursion
        // I_n = I_(n-1) + (e/(n-1)) dI_(n-1)/de (from d/de[(1-e cosE)^-(n-1)]
        // integrated over one period).  <dE/dt> above reproduces
        // dipoleEccentricityFactor() exactly -- an independent check that
        // this is the same physics already used for the Larmor comparison,
        // not a new assumption -- and gamma cancels in the ratio, so k does
        // NOT depend on which reaction model (coherent or individual LL)
        // is active, only on the force law's r-dependence:
        //
        //   k = [<dL/dt>/L] / [d|E|/dt / |E|] = -(1-e^2)/(2+e^2)
        //
        // Checked against the (well-conditioned, e>0.02) part of the
        // measured seed-42 trajectory: matches the empirical k to 1-2e-4
        // absolute at every checkpoint (e.g. e^2=0.063: measured -0.454356,
        // formula -0.454154), i.e. the old measurement was already reading
        // this same formula, just with noise on top that grew as e shrank.
        // e^2 here is exactly the periapsis discriminant, so no new state is
        // needed to evaluate it.
        const double eccentricitySquared=std::max(0.0,1.0
            +2.0*elements.specificEnergy*elements.specificAngularMomentum
                *elements.specificAngularMomentum
                /(attractionParameter*attractionParameter));
        const double angularExponent=
            -(1.0-eccentricitySquared)/(2.0+eccentricitySquared);
        // (sAtPhoton, beta immediately BEFORE that photon's own kick) for
        // every photon fired within this checkpoint's while loop below,
        // read back after the loop to gamma-weight labFrameTimeTotal's own
        // increment per constant-beta segment (README point N).  Stays
        // empty -- and the increment below collapses to the simple
        // gamma(beta-entering-this-checkpoint) case -- for checkpoints that
        // fire no photon, and for the deterministic (non-stochastic) branch
        // entirely.
        std::vector<std::pair<double,double>> photonTimingsThisCheckpoint;
        if(isStochastic) {
            // Replace the deterministic bulk jump below with the sum of
            // individually fired, Poisson-timed photons over this same
            // span -- see stochasticElectricDipole's own comment
            // (electrodynamics.hpp) for why applying BOTH here would
            // double the very quantity the photons are meant to
            // discretize.  jumpParameter/energyGrowth/angularExponent
            // above still describe the ASSUMED classical envelope this
            // span's hazard integral (and each photon's individual
            // energy) is measured against -- inputs to the sampling
            // below, never applied directly for this model.
            //
            // skipHazard = integral(power/(hbar*omega) dt) over the whole
            // skip, in closed form against that same u(n)=u0(1-Jx)^(-2/3)
            // envelope (J=jumpParameter, x=n/orbitsToSkip in [0,1]),
            // T(n)=T0(1-Jx), power(n)=power0*(u(n)/u0)^4, omega(n)=
            // omega0*(u(n)/u0)^1.5 -- verified against brute-force
            // numerical quadrature to 1e-12 relative across the full
            // jumpParameter range this code ever produces (see README).
            // Photon energy: the LEVEL DIFFERENCE where the Bohr ladder has
            // a rung below, and hbar*omega where it has not.
            //
            // hbar*omega is the correspondence-principle value and it is the
            // right one asymptotically -- measured, dE(n->n-1)/hbar*omega is
            // 1.0152 at n=100, 1.0523 at n=30, 1.1728 at n=10 -- but it is
            // wrong by a factor of three at n=2, where the ladder's spacing is
            // nothing like the local orbital frequency.  Using the spacing
            // itself is what the emission physically is: the pair drops a
            // level and the photon carries what the drop released.
            //
            // Below n=2 there is no lower rung and the expression would
            // diverge as n->1, so the fallback is the historical hbar*omega.
            // That is not a patch over an awkward limit: the model's own
            // domain sits at n<=1.09 (measured across all four phenomena),
            // i.e. AT or BELOW the ground state, where a level difference is
            // not defined at all and the classical continuum is the only
            // description available.  The default --level 1 therefore never
            // reaches this branch and is bit-identical to before.
            const double photonEnergyReference=[&]{
                const double classical=hbar*(2.0*pi/period);
                const double bindingScale=pairBindingEnergy(activePair);
                const double orbitalEnergy=
                    reducedMass*std::abs(elements.specificEnergy);
                if(!(orbitalEnergy>0.0)||!(bindingScale>0.0)) return classical;
                const double level=std::sqrt(bindingScale/orbitalEnergy);
                if(!(level>=2.0)) return classical;
                // R(1/(n-1)^2 - 1/n^2), finite for every n >= 2 (0.75 R at
                // n=2 and falling), evaluated at the continuous level rather
                // than snapped to an integer: the orbital energy here is a
                // continuum, and pretending otherwise would quantize the
                // dynamics, which this mode deliberately does not do.
                const double lower=level-1.0;
                return bindingScale*(1.0/(lower*lower)-1.0/(level*level));
            }();
            // CREM_HARMONIC: see eccentricOrbitHazardSuppression's own
            // comment for the full derivation.  hazardReference (NOT
            // photonEnergyReference itself) drives the skip-hazard
            // integral below, so the EVENT RATE reflects S(e) -- fewer
            // events than the n=1-only estimate once power is genuinely
            // spread across many harmonics.  photonEnergyReference is left
            // untouched for use as the n=1 unit each fired photon's own
            // sampled harmonic multiplies below.
            // Production default since this commit: the harmonic
            // correction is now ON unless explicitly disabled
            // (CREM_HARMONIC=0), the same on-by-default/opt-out shape as
            // the stochasticElectricDipole promotion itself.  The env var
            // is kept, not removed, specifically so the pre-harmonic
            // behaviour stays one flag away for regression comparison --
            // see this file's own README section (point K) for the
            // measurement that justified promoting it.
            //
            // SIZE OF THE EFFECT, remeasured.  That section's number came
            // from a single batch: seed 42, fifteen trajectories, 60 s
            // budget, giving RMST 0.304 ns against 0.284 ns, i.e. +7.0%.
            // The comment here used to quote it as +7.7%, which is not what
            // the README's own two numbers give.
            //
            // The batch reproduces exactly -- rerun today it gives 306.8 ps
            // against 288.9 ps, with the median matching the recorded
            // 0.198748 ns to five significant figures -- but it is far too
            // small to size the effect: the two RMSTs carry +/-77 and +/-70
            // ps against a difference of 18 ps, so a single batch measures
            // the effect four times less precisely than the effect itself.
            //
            // Measured properly instead, paired across sixteen seeds at
            // twenty trajectories each (the pairing matters: the initial
            // conditions dominate the 57% spread in collapse time and are
            // shared between the two configurations at a given seed, so
            // they cancel in the difference):
            //
            //     +4.82%, 95% CI [+3.64%, +6.01%], sign positive in 16/16
            //
            // The per-seed spread is 2.41%, which puts the original +7.0%
            // at 0.90 sigma -- an ordinary draw, not a discrepancy.  The
            // median is unchanged, as that section says: the collapse time's
            // median is set by the energy budget, which this correction does
            // not touch.
            const char* harmonicEnv=std::getenv("CREM_HARMONIC");
            const bool harmonicCorrection=
                !harmonicEnv||std::strcmp(harmonicEnv,"0")!=0;
            const double eccentricityHere=std::sqrt(eccentricitySquared);
            const double hazardSuppression=harmonicCorrection
                ?eccentricOrbitHazardSuppression(eccentricityHere):1.0;
            const double hazardReference=photonEnergyReference
                /std::max(hazardSuppression,1.0e-12);
            if(hazardReference>0.0) {
                const double integralFactor=jumpParameter>1.0e-12
                    ?(3.0/jumpParameter)
                        *(1.0-std::pow(1.0-jumpParameter,1.0/3.0))
                    :1.0;
                const double skipHazard=lossPerOrbit*reducedMass
                    /hazardReference
                    *static_cast<double>(orbitsToSkip)*integralFactor;
                double hazardConsumedThisSkip=0.0;
                stochasticSkipHazard+=skipHazard;
                int photonCountDebug=0;
                // Whether elements.specificEnergy/specificAngularMomentum
                // have already been moved by an EARLIER photon within this
                // same while loop pass -- see the refresh block at the top
                // of the loop body, README point E4, for why this matters:
                // period/eccentricityHere/photonEnergyReference above are
                // all evaluated ONCE, from the state at THIS checkpoint's
                // own start, before any photon of this cascade has fired.
                bool cascadeStateAlreadyMoved=false;
                if(std::getenv("CREM_DEBUG"))
                    std::cerr<<"  SKIP orbitsToSkip="<<orbitsToSkip
                             <<" jumpParameter="<<jumpParameter
                             <<" skipHazard="<<skipHazard
                             <<" cumHazard="<<stochasticSkipHazard
                             <<" threshold="<<stochasticSkipThreshold
                             <<" photonEnergyReference="
                             <<photonEnergyReference<<"J"
                             <<" hazardSuppression="<<hazardSuppression
                             <<std::endl;
                while(stochasticSkipHazard>=stochasticSkipThreshold) {
                    // TEMPORARY, investigation-only: full round-trip
                    // precision (setprecision(17)) state dump, gated on a
                    // separate env var so it never disturbs CREM_DEBUG's
                    // existing output. Exists purely to trace the exact
                    // bit-level origin of the para/ortho divergence
                    // documented in README point L, sondy 9-11.
                    if(std::getenv("CREM_DEBUG_PRECISE")) {
                        std::cerr<<std::setprecision(17)
                                 <<"    PRECISE pre-draw: stream="
                                 <<stochasticSkipStream
                                 <<" hazard="<<stochasticSkipHazard
                                 <<" threshold="<<stochasticSkipThreshold
                                 <<" E="<<elements.specificEnergy
                                 <<" L="<<elements.specificAngularMomentum
                                 <<" Ldir=("<<angularMomentumDirection.x
                                 <<","<<angularMomentumDirection.y
                                 <<","<<angularMomentumDirection.z<<")"
                                 <<" vcm=("<<centreOfMassVelocity.x
                                 <<","<<centreOfMassVelocity.y
                                 <<","<<centreOfMassVelocity.z<<")"
                                 <<std::setprecision(6)<<std::endl;
                    }
                    // README point E4: photonEnergyReference/eccentricityHere
                    // above were evaluated ONCE, from this checkpoint's OWN
                    // starting state, before any photon of this cascade had
                    // fired -- correct for the first photon (nothing has
                    // moved yet) but stale for every later one, because
                    // elements.specificEnergy/specificAngularMomentum below
                    // ARE updated per photon within this same while loop.
                    // Measured effect of not refreshing: photons 2 and 3 of
                    // a 3-photon cascade came out 5.8x and 13.5x too small.
                    // Fix: for every photon after the first, recompute the
                    // period/eccentricity/energy-quantum reference from the
                    // CURRENT state directly (energyRatio=1 -- an exact
                    // value needs no envelope extrapolation, unlike the
                    // first photon, which still uses the checkpoint's own
                    // sAtPhoton/energyRatio machinery below because for it
                    // period/eccentricityHere ARE still the current state).
                    // orbitsToSkip/jumpParameter/skipHazard themselves are
                    // NOT re-derived here: they describe how many orbits
                    // this checkpoint's hazard integral spans and how the
                    // Poisson threshold-crossing bookkeeping already in
                    // flight (stochasticSkipHazard/stochasticSkipThreshold)
                    // is being consumed, independent of which reference a
                    // given photon's OWN energy/harmonic draw uses.
                    double effectivePhotonEnergyReference=photonEnergyReference;
                    double effectiveEccentricityHere=eccentricityHere;
                    bool useExactEnergyRatio=false;
                    if(cascadeStateAlreadyMoved) {
                        const OsculatingElements currentElements{
                            elements.specificEnergy,
                            elements.specificAngularMomentum};
                        const double refreshedPeriod=regularizedPeriod(
                            currentElements,attractionParameter,
                            separationFloor());
                        effectivePhotonEnergyReference=
                            hbar*(2.0*pi/refreshedPeriod);
                        const double refreshedEccentricitySquared=
                            std::max(0.0,1.0
                                +2.0*elements.specificEnergy
                                    *elements.specificAngularMomentum
                                    *elements.specificAngularMomentum
                                    /(attractionParameter
                                        *attractionParameter));
                        effectiveEccentricityHere=
                            std::sqrt(refreshedEccentricitySquared);
                        useExactEnergyRatio=true;
                    }
                    stochasticSkipHazard-=stochasticSkipThreshold;
                    hazardConsumedThisSkip+=stochasticSkipThreshold;
                    // Position within [0,1] of this skip where the running
                    // hazard integral reaches this photon's threshold,
                    // found by inverting the same closed form above.  A
                    // threshold consumed mostly out of hazard CARRIED IN
                    // from the previous skip (hFraction saturating at the
                    // clamp) is attributed to this skip's own start
                    // instead -- a boundary approximation, not exact.
                    // Measured, not just flagged (README point E4): in a
                    // 3-photon cascade within one skip call, all three got
                    // an IDENTICAL photonEnergy despite |E| growing by an
                    // order of magnitude between them, because the
                    // threshold was large enough that hFraction saturated
                    // for every one of them -- sAtPhoton pins to
                    // jumpParameter (the approximation's own ceiling)
                    // regardless of how much further the true state has
                    // already moved.  Photons 2 and 3 in that cascade came
                    // out 5.8x and 13.5x too small relative to
                    // photonEnergyReference recomputed from their own,
                    // already-updated E -- growing without bound deeper
                    // into a cascade.  FIXED below (README point E4): this
                    // paragraph still describes the approximation correctly
                    // for the FIRST photon of a cascade (period/
                    // eccentricityHere/photonEnergyReference above ARE the
                    // current state for it, nothing has moved yet) -- the
                    // cascadeStateAlreadyMoved block a few lines down
                    // refreshes them from the true current state for every
                    // photon after the first, instead -- so sAtPhoton/
                    // energyRatio below now only ever apply to that first
                    // photon, for which they are exact, not an
                    // approximation.
                    const double hFraction=skipHazard>0.0
                        ?std::clamp(
                            hazardConsumedThisSkip/skipHazard,0.0,1.0)
                        :1.0;
                    const double base=jumpParameter>1.0e-12
                        ?1.0-std::pow(1.0-jumpParameter,1.0/3.0):1.0;
                    const double x=jumpParameter>1.0e-12
                        ?(1.0-std::pow(1.0-hFraction*base,3.0))
                            /jumpParameter
                        :hFraction;
                    const double sAtPhoton=
                        jumpParameter*std::clamp(x,0.0,1.0);
                    const double energyRatio=
                        std::pow(1.0-sAtPhoton,-2.0/3.0);
                    // CREM_HARMONIC: which harmonic of omega_orb this
                    // specific event actually belongs to -- see
                    // eccentricOrbitHazardSuppression's own comment.
                    // Sampled independently of the hazard-rate correction
                    // above (that one only rescales HOW OFTEN events fire,
                    // this one decides, given that one just fired, which
                    // harmonic it is), from the directly-tabulated 2D
                    // quantile table (see eccentricOrbitHarmonicNumber's
                    // own comment -- no eccentricity gate needed: unlike
                    // the scale-invariance-based version this replaced,
                    // every table entry is independently computed at its
                    // own eccentricity, not extrapolated from a
                    // higher-eccentricity calibration set).  Uses
                    // effectiveEccentricityHere, NOT the once-per-checkpoint
                    // eccentricityHere the hazard-rate/hazardSuppression
                    // integral above is still (correctly) frozen at -- that
                    // integral genuinely describes the whole skip in
                    // aggregate, but THIS specific photon's own harmonic
                    // should reflect the orbit it actually fired from,
                    // refreshed above for every photon after the first in
                    // the same cascade (README point E4).
                    const int harmonicNumber=harmonicCorrection
                        ?std::max(1,static_cast<int>(std::lround(
                            eccentricOrbitHarmonicNumber(
                                effectiveEccentricityHere,
                                drawUniformUnit(stochasticSkipStream)))))
                        :1;
                    const double photonEnergy=effectivePhotonEnergyReference
                        *std::pow(useExactEnergyRatio?1.0:energyRatio,1.5)
                        *static_cast<double>(harmonicNumber);
                    // Photon emission angle relative to the CURRENT orbital
                    // plane normal, drawn from the actual angular pattern of
                    // a rotating (not linearly oscillating) E1 dipole:
                    // dP/dOmega proportional to (1+cos^2(theta)), theta from
                    // the angular-momentum axis -- maximal along the axis,
                    // half that in the plane, not the sin^2(theta)
                    // (in-plane-maximal) pattern of a single linear dipole.
                    // Inverted in closed form via Cardano's formula for the
                    // depressed cubic mu^3+3mu+(4-8u)=0 that solving
                    // CDF(mu)=u for this distribution reduces to (mu=cos
                    // theta); verified against the CDF numerically to
                    // 1e-15 absolute over the full [0,1) range of u.
                    const double cardanoQ=4.0-8.0*drawUniformUnit(
                        stochasticSkipStream);
                    const double cardanoDiscriminant=
                        (cardanoQ*cardanoQ)/4.0+1.0;
                    const double cardanoSqrt=std::sqrt(cardanoDiscriminant);
                    const auto signedCbrt=[](double value) {
                        return std::copysign(std::cbrt(std::abs(value)),value);
                    };
                    const double cosThetaFromAxis=
                        signedCbrt(-cardanoQ/2.0+cardanoSqrt)
                        +signedCbrt(-cardanoQ/2.0-cardanoSqrt);
                    // Shared orthonormal basis (e1,e2) perpendicular to the
                    // current orbital-plane normal -- reused below for both
                    // the photon's own emission azimuth (real recoil) and
                    // the tilt axis (unknown-phase random walk).  Nothing
                    // physical ties a "reference" azimuth to either use, so
                    // one basis serves both.
                    const Vec3 seedAxis=
                        std::abs(angularMomentumDirection.z)<0.9
                            ?Vec3{0.0,0.0,1.0}:Vec3{1.0,0.0,0.0};
                    Vec3 inPlaneFirst=cross(angularMomentumDirection,seedAxis);
                    inPlaneFirst=inPlaneFirst*(1.0/inPlaneFirst.norm());
                    const Vec3 inPlaneSecond=
                        cross(angularMomentumDirection,inPlaneFirst);
                    const double sinThetaFromAxis=std::sqrt(std::max(0.0,
                        1.0-cosThetaFromAxis*cosThetaFromAxis));
                    // Linear-momentum recoil, the fix this whole block
                    // exists for.  CREM's bound initial conditions are
                    // prepared at EXACTLY zero total momentum (see this
                    // function's own comment on centreOfMassVelocity), and
                    // every continuous model keeps that true by
                    // construction; only this discrete-photon model had
                    // been letting the photon carry momentum away
                    // unbalanced.  The photon's FULL 3D direction (not just
                    // its angle from the axis) is needed here, unlike for
                    // the tilt below, because the recoil is a real vector
                    // kick, not an averaged magnitude.
                    const double photonAzimuth=
                        2.0*pi*drawUniformUnit(stochasticSkipStream);
                    const Vec3 photonDirection=
                        angularMomentumDirection*cosThetaFromAxis
                        +(inPlaneFirst*std::cos(photonAzimuth)
                          +inPlaneSecond*std::sin(photonAzimuth))
                            *sinThetaFromAxis;
                    if(std::getenv("CREM_DEBUG_PRECISE")) {
                        std::cerr<<std::setprecision(17)
                                 <<"    PRECISE draws: harmonicNumber="
                                 <<harmonicNumber
                                 <<" cardanoQ="<<cardanoQ
                                 <<" cosThetaFromAxis="<<cosThetaFromAxis
                                 <<" photonAzimuth="<<photonAzimuth
                                 <<" photonEnergy="<<photonEnergy
                                 <<" photonDirection=("<<photonDirection.x
                                 <<","<<photonDirection.y
                                 <<","<<photonDirection.z<<")"
                                 <<std::setprecision(6)<<std::endl;
                    }
                    // LAB-FRAME PHOTON KINEMATICS, diagnostic only.
                    // Everything above (photonEnergy, photonDirection) and
                    // everything below (the recoil kick, the orbital-energy/
                    // angular-momentum updates) is computed correctly in S',
                    // the pair's instantaneous rest frame at this emission --
                    // self-consistent internal bookkeeping, unaffected by
                    // this block.  What follows answers a SEPARATE question,
                    // "what would a fixed, distant lab observer measure for
                    // this photon", via the standard special-relativistic
                    // Doppler/aberration boost from S' (moving with
                    // centreOfMassVelocity -- the PRE-kick recoil speed
                    // accumulated from any earlier photons this same
                    // trajectory already fired) to the lab.  Investigated at
                    // length before being written (README point L,
                    // "charakterystyki fotonu").  The exact four-momentum
                    // subtraction below independently uses this same S'
                    // photon and boost to update the pair; this block records
                    // the already-determined lab observables for plots.
                    {
                        const double sourceSpeed=centreOfMassVelocity.norm();
                        const double sourceBeta=sourceSpeed/c;
                        // Recorded PRE-kick, same convention as labPhoton
                        // below: this segment of the checkpoint's elapsed
                        // time -- up to and including this photon's own
                        // sAtPhoton position -- happened at the recoil speed
                        // the pair had BEFORE this photon's momentum kick.
                        photonTimingsThisCheckpoint.emplace_back(
                            sAtPhoton,sourceBeta);
                        LabFramePhoton labPhoton;
                        labPhoton.sourceBeta=sourceBeta;
                        labPhoton.energyJoules=photonEnergy;
                        if(sourceBeta>1.0e-9) {
                            const double gammaFactor=1.0/std::sqrt(
                                std::max(1.0e-300,
                                    1.0-sourceBeta*sourceBeta));
                            const Vec3 betaHat=
                                centreOfMassVelocity*(1.0/sourceSpeed);
                            const double cosThetaSource=
                                dot(photonDirection,betaHat);
                            labPhoton.energyJoules=gammaFactor*photonEnergy
                                *(1.0+sourceBeta*cosThetaSource);
                            const double photonMomentumSource=photonEnergy/c;
                            const double pParallelLab=gammaFactor
                                *photonMomentumSource
                                *(cosThetaSource+sourceBeta);
                            const Vec3 pPerpLab=
                                (photonDirection-betaHat*cosThetaSource)
                                    *photonMomentumSource;
                            const Vec3 pLab=betaHat*pParallelLab+pPerpLab;
                            const double pLabNorm=pLab.norm();
                            if(pLabNorm>1.0e-300) {
                                const double cosThetaLabFromRecoil=
                                    dot(pLab,betaHat)/pLabNorm;
                                labPhoton.angleFromRecoilAxisRadians=
                                    std::acos(std::clamp(
                                        cosThetaLabFromRecoil,-1.0,1.0));
                            }
                        }
                        result.labFramePhotons.push_back(labPhoton);
                        if(std::getenv("CREM_DEBUG"))
                            std::cerr<<"    LAB FRAME: beta="<<sourceBeta
                                     <<" E_source="<<photonEnergy
                                     <<"J E_lab="<<labPhoton.energyJoules
                                     <<"J (E_lab/E_source-1)="
                                     <<(labPhoton.energyJoules
                                         /photonEnergy-1.0)
                                     <<" angleFromRecoilAxis="
                                     <<labPhoton.angleFromRecoilAxisRadians
                                     <<std::endl;
                    }
                    // Exact composite four-momentum recoil.  The osculating
                    // path does not retain the two instantaneous particle
                    // momenta, but it does retain exactly what is needed for
                    // the pair as a bound composite: its invariant energy
                    // W=M c^2+E_orbit and its lab velocity.  In the incoming
                    // rest frame P'=(W/c,0) and k'=(E_gamma/c,E_gamma*n/c),
                    // hence W_f^2=W^2-2 W E_gamma.  Boost P' and k' to the
                    // lab, subtract them, and read the new COM velocity from
                    // c^2 p_f/E_f.  This remains exact at every accumulated
                    // recoil speed and replaces the former Newtonian
                    // -p_gamma/M velocity increment plus kinetic correction.
                    const double restEnergy=totalMass*c*c;
                    const double invariantEnergyBefore=
                        restEnergy+reducedMass*elements.specificEnergy;
                    const long double invariantSquaredAfter=
                        static_cast<long double>(invariantEnergyBefore)
                            *invariantEnergyBefore
                        -2.0L*invariantEnergyBefore*photonEnergy;
                    if(!(invariantSquaredAfter>0.0L)) {
                        stochasticSkipThreshold=
                            drawEmissionThreshold(stochasticSkipStream);
                        continue;
                    }
                    const double invariantEnergyAfter=static_cast<double>(
                        std::sqrt(invariantSquaredAfter));
                    const two_body::ParticleFourMomentum pairRestBefore{
                        invariantEnergyBefore,{0,0,0}};
                    const two_body::ParticleFourMomentum photonRest{
                        photonEnergy,photonDirection*(photonEnergy/c)};
                    const auto pairLabBefore=two_body::boostFourMomentum(
                        pairRestBefore,centreOfMassVelocity);
                    const auto photonLabFour=two_body::boostFourMomentum(
                        photonRest,centreOfMassVelocity);
                    const two_body::ParticleFourMomentum pairLabAfter{
                        pairLabBefore.energy-photonLabFour.energy,
                        pairLabBefore.momentum-photonLabFour.momentum};
                    const Vec3 recoiledVelocity=
                        two_body::velocityFromFourMomentum(pairLabAfter);
                    if(!isFinite(recoiledVelocity)) {
                        stochasticSkipThreshold=
                            drawEmissionThreshold(stochasticSkipStream);
                        continue;
                    }
                    centreOfMassVelocity=recoiledVelocity;
                    const double recoiledSpecificEnergy=
                        (invariantEnergyAfter-restEnergy)/reducedMass;
                    // ANGULAR MOMENTUM DIRECTION: NOT tilted here, and this
                    // is a finding, not an omission.  An earlier version of
                    // this block tilted angularMomentumDirection by
                    // treating the recoil as an r x delta-v kick to the
                    // RELATIVE motion (delta-v ~ photonEnergy/(c*reducedMass),
                    // r ~ semi-major axis).  That is inconsistent with the
                    // linear-momentum fix directly above: this system's
                    // bound initial conditions are prepared at r1 = r_rel*
                    // m2/M, r2 = -r_rel*m1/M (CM at the origin), so
                    // Sum m_i r_i = 0 identically, and the uniform
                    // centre-of-mass kick applied above -- m_i*deltaV_cm for
                    // both particles, the SAME vector -- therefore satisfies
                    // Sum r_i x (m_i*deltaV_cm) = deltaV_cm x Sum(m_i r_i)
                    // = 0 EXACTLY, for any mass ratio (checked numerically:
                    // e+e-, and 1836:1 and 3:1 mass ratios, residual
                    // 1e-17).  A uniform push through the system's own
                    // centre of mass cannot torque it.  The removed tilt
                    // was therefore spending the SAME photon momentum
                    // twice: once as the centre-of-mass recoil that
                    // conserves linear momentum (above), and again,
                    // uncoordinated with it, as an independent kick to the
                    // relative motion that conservation does not call for.
                    //
                    // The genuine remaining source of orbital-plane torque
                    // is the photon's own angular momentum, and for a real
                    // photon that is dominated by SPIN (+-hbar along its own
                    // propagation direction, an exact, universal fact for
                    // any spin-1 massless boson -- not an orbit-averaged
                    // estimate the way the orbital piece r x p would be).
                    // It is not modelled here because it is not a small
                    // correction: CREM's own starting point is the Bohr/SED
                    // value L = hbar (see physical_constants.hpp), so
                    // |L_photon|/|L_orbital| = hbar/hbar = 1 exactly at the
                    // very first photon, by construction of the initial
                    // condition, not as a result depending on any computed
                    // quantity.  Removing one photon's worth of hbar from an
                    // angular momentum that starts at one hbar is an O(1)
                    // disruption, every time, not a perturbation the
                    // orbit-averaged k ratio below (or any classical
                    // adiabatic treatment) is built to absorb -- this model
                    // sits exactly at the quantum limit where a real bound
                    // system's angular momentum is a discrete hbar-spaced
                    // ladder, not a continuous classical vector, so no
                    // classical bookkeeping trick closes this gap.  Left
                    // unmodelled and documented rather than forced in: the
                    // orbital-plane direction is carried forward unchanged
                    // by photon recoil (self-consistent with the linear-
                    // momentum fix above), and the model's magnitude-only
                    // representation is the honest limit of what this
                    // architecture can support.
                    // Eccentricity/k evaluated BEFORE this photon's own
                    // kick, mirroring the bulk formula's own convention
                    // (its eccentricitySquared above is likewise the
                    // pre-jump value).
                    const double eccentricitySquaredHere=std::max(0.0,1.0
                        +2.0*elements.specificEnergy
                            *elements.specificAngularMomentum
                            *elements.specificAngularMomentum
                            /(attractionParameter*attractionParameter));
                    const double energyBeforeKick=
                        std::abs(elements.specificEnergy);
                    elements.specificEnergy=recoiledSpecificEnergy;
                    const double energyAfterKick=
                        std::abs(elements.specificEnergy);
                    // From here on, any further photon drawn within this
                    // same while loop pass must refresh its own energy/
                    // eccentricity reference above rather than reuse this
                    // checkpoint's stale, pre-cascade one (README point E4).
                    cascadeStateAlreadyMoved=true;
                    // CLASSICAL magnitude prediction from k(e), exactly
                    // integrated over this finite jump (see this file's own
                    // history/README for the derivation and its own
                    // verification).  History: for a long time this was
                    // kept ONLY as a diagnostic comparison, never applied to
                    // the state, on the reasoning that k comes from the
                    // classical reaction TORQUE, and "torque integrated on
                    // the orbit" and "what the continuous field carries
                    // away" are the SAME quantity by angular-momentum
                    // conservation -- so using this number ON TOP OF an
                    // independent, ALSO-magnitude-setting photon spin
                    // vector would double-count it.  That reasoning still
                    // holds for exactly the failure mode it describes
                    // (adding this to a full vector kick that separately
                    // sets magnitude, tried and rejected below, see (2) and
                    // (3)) -- but it does NOT forbid using this value to
                    // REPLACE the magnitude outright, the same
                    // "replace, don't add" relationship photonEnergy
                    // already has to the classical power rate.  Now used
                    // that way: see (4) below, in the photon-spin block,
                    // for the derivation and the numerical verification
                    // that motivated the switch (README point L).
                    double classicalEccentricitySquared=eccentricitySquaredHere;
                    if(eccentricitySquaredHere>1.0e-12
                       &&energyBeforeKick>0.0) {
                        const double energyRatio=
                            energyAfterKick/energyBeforeKick;
                        const double e0Fourth=eccentricitySquaredHere
                            *eccentricitySquaredHere;
                        const double rhsScale=
                            energyRatio*energyRatio*energyRatio
                            *(1.0-eccentricitySquaredHere)
                            *(1.0-eccentricitySquaredHere)
                            *(1.0-eccentricitySquaredHere);
                        double lo=0.0, hi=1.0-1.0e-15;
                        for(int iter=0; iter<80; ++iter) {
                            const double mid=0.5*(lo+hi);
                            const double lhs=
                                (1.0-mid)*(1.0-mid)*(1.0-mid)*e0Fourth;
                            const double rhs=rhsScale*mid*mid;
                            if(lhs>rhs) lo=mid; else hi=mid;
                        }
                        classicalEccentricitySquared=0.5*(lo+hi);
                    }
                    const double classicalAngularMomentumMagnitude=
                        std::sqrt(std::max(0.0,
                            attractionParameter*attractionParameter
                            *(1.0-classicalEccentricitySquared)
                            /(2.0*energyAfterKick)));
                    // PHOTON SPIN, applied for real: a photon is a massless
                    // spin-1 boson, so it carries EXACTLY +-hbar of angular
                    // momentum along its own propagation direction
                    // (photonDirection, already sampled above) -- not an
                    // orbit-averaged estimate, an exact per-photon fact,
                    // same status as photonEnergy=hbar*omega.  The
                    // conditional helicity distribution given the emission
                    // angle theta (cosThetaFromAxis, already sampled) is
                    // the standard result for a Delta-m=+-1 (circular)
                    // dipole transition: P(+|theta)=(1+cos theta)^2 /
                    // [2(1+cos^2 theta)], P(-|theta)=(1-cos theta)^2 /
                    // [2(1+cos^2 theta)] -- these sum to exactly the
                    // (1+cos^2 theta) pattern the emission angle itself was
                    // already drawn from, so this is not a new assumption
                    // layered on top, it is the polarization-resolved
                    // refinement of the SAME distribution (checked: at
                    // theta=0 and theta=pi the two poles give a helicity
                    // that, combined with photonDirection, ALWAYS points
                    // the actual angular-momentum vector h*hbar*n along the
                    // orbital axis, matching Delta-m=+1 exactly; at
                    // theta=pi/2 it is an even 50/50 mix in the orbital
                    // plane).  What this does NOT capture: the photon's
                    // ORBITAL angular momentum relative to the pair (r x p,
                    // r being the true anomaly this elements-only
                    // representation does not carry) -- checked by taking
                    // the expectation of the axial component over the full
                    // (1+cos^2 theta) distribution, <h cos theta> = 1/2,
                    // not 1, so spin alone recovers exactly half of the
                    // Delta-m=1 selection rule on average, the other half
                    // being the same true-anomaly-dependent piece the
                    // removed tilt could not evaluate either.  Left as an
                    // acknowledged, unavoidable shortfall of this
                    // architecture rather than patched with an arbitrary
                    // factor of 2, which would misrepresent an unverified
                    // guess as a measured correction.
                    //
                    // THREE attempts made at closing this gap, all tried
                    // and REJECTED before touching the state below (see
                    // README point L, "orbitalny moment pędu fotonu" for
                    // the full, numbered writeup):
                    //
                    // (1) Literal r x p from a sampled true anomaly
                    // (uniform in mean anomaly, Kepler-solved r(E)=a(1-e
                    // cosE)): parametrically too small.  For dipole
                    // radiation from a point-like source, |r x p_photon| ~
                    // a*(hbar*omega_photon/c) = hbar*n*(v_orbit/c), the
                    // standard multipole suppression of an orbital
                    // contribution relative to the intrinsic-hbar scale.
                    // v_orbit/c ~ 1e-3 for positronium, three orders of
                    // magnitude too small to close a factor-of-2 gap,
                    // whatever r is sampled to.
                    //
                    // (2) Uniformly doubling the spin vector
                    // (photonDirection*2*helicity*hbar), calibrated to
                    // exactly reproduce the classical
                    // k(e)=-(1-e^2)/(2+e^2) secular rate in ENSEMBLE
                    // AVERAGE (closed form: required average axial removal
                    // per photon is 2*|k(e)|*sqrt(1-e^2)/S(e), which the
                    // doubled term matches to within 7% for e in [0,0.97]
                    // -- not numerology, the textbook 50/50 spin/orbital
                    // split for a uniformly rotating classical dipole).
                    // IMPLEMENTED, MEASURED, REVERTED: on 10 independent
                    // seeds the resulting |L| after the first photon came
                    // out 2.3x-4.5x the classical target, not closer to it
                    // -- WORSE than the unpatched spin-only version.  Root
                    // cause, found and proven exactly (law of cosines):
                    // |L_after| = |L_before|*sqrt(1+x^2-2xy), x = R/|
                    // L_before|, y = cos(theta)*helicity in [-1,1].  R =
                    // 2*hbar/reducedMass is, by the SAME Bohr/SED
                    // convention this model's own initial condition uses
                    // (L_initial = hbar/reducedMass), almost exactly TWICE
                    // the model's own characteristic angular-momentum
                    // scale -- so x sits at 1.8-2.2 on every trajectory
                    // tested, not by chance but by construction.  At
                    // x>=2 the minimum of sqrt(1+x^2-2xy) over ALL possible
                    // y is |x-1|>=1: a magnitude DECREASE is geometrically
                    // impossible for every single photon draw, not merely
                    // unlikely.  Any uniform rescaling of a single ~hbar-
                    // scale vector kick hits this same wall, because
                    // positronium's own L lives at the hbar scale by
                    // construction, for the entire trajectory, not just
                    // near the first photon.
                    //
                    // (3) Solve exactly (not calibrate in average) for the
                    // magnitude R along the ALREADY-SAMPLED photonDirection
                    // that hits classicalAngularMomentumMagnitude exactly
                    // for THIS one photon: |L_before + R*n_hat|^2 = target^2
                    // is a quadratic in R with a real solution only if
                    // target >= |L_before|*sin(theta) (the geometric floor
                    // reachable by moving along a FIXED direction).  Tested
                    // on the same 10 seeds using each one's own sampled
                    // theta: 10/10 UNREACHABLE (target strictly below the
                    // floor in every case, because k(e) typically demands
                    // a ~45-50% magnitude drop per photon -- these are
                    // large discrete jumps, same reason photonEnergy itself
                    // jumps energy by O(1) -- while only ~19% of
                    // (1+cos^2 theta)-weighted directions land close enough
                    // to axis to reach that big a drop from any magnitude
                    // at all).  This is the decisive result: it is not a
                    // question of finding the right R, exact or otherwise
                    // -- a single scalar-times-fixed-direction vector has
                    // one degree of freedom and is being asked to satisfy
                    // two generically incompatible constraints (hit the
                    // classical target magnitude; point along a direction
                    // fixed by unrelated energy/polarization physics).
                    // (4) IMPLEMENTED: decouple magnitude from direction,
                    // instead of asking one vector to set both.  The root
                    // cause common to (2) and (3) was never the SIZE of the
                    // kick, it was that a single scalar-times-fixed-
                    // direction vector has one degree of freedom serving
                    // two generically incompatible jobs (hit the classical
                    // target magnitude; point along a direction fixed by
                    // unrelated energy/polarization physics).  Splitting
                    // the jobs removes the conflict instead of trying to
                    // satisfy it:
                    //   magnitude <- classicalAngularMomentumMagnitude
                    //     (the exact-ODE-integrated k(e) result computed
                    //     above, independently verified to 1-2e-4 against a
                    //     measured trajectory -- REPLACES the magnitude,
                    //     exactly the same "replace, don't add" relationship
                    //     photonEnergy already has to the classical power
                    //     rate, so this is not the double-counting the
                    //     comment above (1)-(3) warns against: that warning
                    //     was about ADDING k(e) on top of an independent
                    //     full vector that ALSO sets magnitude, not about
                    //     assigning k(e) the magnitude role outright.
                    //   direction <- the same spin-vector subtraction as
                    //     before, used ONLY for its direction (normalized),
                    //     never for its norm.  Peters-Mathews-style k(e) is
                    //     a magnitude-only result by construction (derived
                    //     orbit-averaged with the orbital plane held fixed)
                    //     -- it says nothing about direction, so assigning
                    //     the photon's real, Stokes-V-verified spin
                    //     structure to fill exactly that gap does not
                    //     double-count anything either.
                    // This sidesteps all three dead ends above: no r x p
                    // reconstruction (1), no single-vector geometry forcing
                    // a magnitude INCREASE for x>=2 (2), no reachability
                    // constraint on a fixed direction (3) -- because
                    // magnitude is now assigned, not solved for or read off
                    // a subtraction's norm.  Side effect, verified below:
                    // this also removes the spurious e^2->0 clamp
                    // (README point L) that the old magnitude-from-norm
                    // approach caused whenever the spin kick left |L| too
                    // large for the now-more-negative E, since (E,L) after
                    // this photon are consistent with each other by
                    // construction (both trace back to the same exact-ODE
                    // integration).
                    const double helicityPlusProbability=std::clamp(
                        (1.0+cosThetaFromAxis)*(1.0+cosThetaFromAxis)
                        /(2.0*(1.0+cosThetaFromAxis*cosThetaFromAxis)),
                        0.0,1.0);
                    const double helicity=
                        drawUniformUnit(stochasticSkipStream)
                            <helicityPlusProbability?1.0:-1.0;
                    const Vec3 photonSpinAngularMomentum=
                        photonDirection*(helicity*hbar);
                    const Vec3 orbitalAngularMomentumBefore=
                        angularMomentumDirection
                            *(elements.specificAngularMomentum*reducedMass);
                    const Vec3 directionTrial=
                        orbitalAngularMomentumBefore-photonSpinAngularMomentum;
                    const double directionTrialNorm=directionTrial.norm();
                    if(directionTrialNorm>1.0e-300) {
                        angularMomentumDirection=
                            directionTrial*(1.0/directionTrialNorm);
                    }
                    // else: keep the previous direction rather than an
                    // undefined one -- the same fallback the old norm-based
                    // branch used (leaving specificAngularMomentum at its
                    // prior value in that branch); here the magnitude line
                    // below still runs regardless, so only direction is
                    // affected by this edge case.
                    elements.specificAngularMomentum=
                        classicalAngularMomentumMagnitude;
                    radiatedEnergyTotal+=photonEnergy;
                    ++photonCountDebug;
                    if(std::getenv("CREM_DEBUG"))
                        std::cerr<<"    PHOTON #"<<photonCountDebug<<" x="<<x
                                 <<" harmonic="<<harmonicNumber
                                 <<" photonEnergy="<<photonEnergy
                                 <<"J W_after="<<invariantEnergyAfter
                                 <<"J |v_cm|="<<centreOfMassVelocity.norm()
                                 <<" energyBeforeKick="<<energyBeforeKick
                                 <<" newE="<<elements.specificEnergy
                                 <<" e0^2="<<eccentricitySquaredHere
                                 <<" helicity="<<helicity
                                 <<" L_spec(magnitude=k(e))="
                                 <<elements.specificAngularMomentum
                                 <<" L_spec(classical k)="
                                 <<classicalAngularMomentumMagnitude
                                 <<" |L_dir|="<<angularMomentumDirection.norm()
                                 <<std::endl;
                    stochasticSkipThreshold=
                        drawEmissionThreshold(stochasticSkipStream);
                }
            }
        } else {
            radiatedEnergyTotal+=
                (updatedEnergyMagnitude-energyMagnitude)*reducedMass;
            elements.specificEnergy=-updatedEnergyMagnitude;
            elements.specificAngularMomentum*=
                std::pow(energyGrowth,angularExponent);
        }
        simulatedTimeTotal+=measuredElapsed*static_cast<double>(orbitsToSkip)
            *(1.0-0.5*jumpParameter);
        // Lab-frame counterpart (README point N): the checkpoint's elapsed
        // proper time as a function of position s in [0,jumpParameter]
        // through its envelope (period shrinks linearly over the
        // orbitsToSkip orbits, so integrating the per-orbit period from 0 to
        // s/jumpParameter gives this closed form -- it reduces to the
        // formula immediately above at s=jumpParameter exactly).  Guarded at
        // jumpParameter<=1e-12 for the same reason sAtPhoton's own formula
        // is (README point E4): in that regime sAtPhoton is forced to ~0 for
        // every photon regardless of hFraction, so "time elapsed up to
        // there" is correctly ~0 too, and the whole checkpoint collapses to
        // one segment at the LATEST beta below.
        const auto properTimeUpToS=[&](double s) {
            if(jumpParameter<=1.0e-12) return 0.0;
            return measuredElapsed*static_cast<double>(orbitsToSkip)
                *(s/jumpParameter)*(1.0-0.5*s);
        };
        {
            double labIncrement=0.0;
            double previousS=0.0;
            for(const auto& timing:photonTimingsThisCheckpoint) {
                const double sAtThisPhoton=timing.first;
                const double betaBeforeThisPhoton=timing.second;
                labIncrement+=gammaFromBeta(betaBeforeThisPhoton)
                    *(properTimeUpToS(sAtThisPhoton)
                        -properTimeUpToS(previousS));
                previousS=sAtThisPhoton;
            }
            // Tail segment: from the last photon fired this checkpoint (or
            // the very start, if none fired) to the checkpoint's own end, at
            // the recoil speed left behind by whichever photons already
            // fired -- exactly centreOfMassVelocity's current value, since
            // nothing else moves it.
            const double totalProperThisCheckpoint=
                measuredElapsed*static_cast<double>(orbitsToSkip)
                *(1.0-0.5*jumpParameter);
            const double tailBeta=centreOfMassVelocity.norm()/c;
            labIncrement+=gammaFromBeta(tailBeta)
                *(totalProperThisCheckpoint-properTimeUpToS(previousS));
            labFrameTimeTotal+=labIncrement;
        }
        revolutionsTotal+=static_cast<double>(orbitsToSkip);
        firstDipole=run.finalState.firstDipole;
        secondDipole=run.finalState.secondDipole;

        if(!(elements.specificEnergy<0.0)||!std::isfinite(elements.specificEnergy)
           ||!std::isfinite(elements.specificAngularMomentum)) {
            result.calibrationOutcome=SimulationOutcome::NumericalFailure;
            result.calibrationSeconds=simulatedTimeTotal;
            result.calibrationSecondsLab=labFrameTimeTotal;
            return result;
        }
    }
    result.calibrationOutcome=SimulationOutcome::ObservationLimit;
    result.calibrationSeconds=simulatedTimeTotal;
    result.calibrationSecondsLab=labFrameTimeTotal;
    return result;
}

std::vector<CremCollapseEstimate> runCremCollapseExperiment(
    std::uint64_t masterSeed,int selectedPhenomenon,int runCount,
    double wallClockBudgetSeconds) {
    // Keep the safety boundary valid for non-CLI callers too.  In particular,
    // comparisons against NaN or +infinity never become true and would disable
    // every wall-clock stop inside estimateCremCollapse().
    if(!(wallClockBudgetSeconds>0.0)
       ||!std::isfinite(wallClockBudgetSeconds)) {
        throw std::invalid_argument(
            "--crem-wallclock-budget-s must be finite and positive");
    }
    std::vector<CremCollapseEstimate> estimates(static_cast<size_t>(runCount));
    std::atomic<int> nextIndex{0};
    std::atomic<int> completed{0};
    std::mutex outputMutex;
    const int workerCount=std::min(runCount,
        static_cast<int>(std::max(1u,std::thread::hardware_concurrency())));
    std::cout<<"Running "<<runCount<<" CREM collapse calibrations on "
             <<workerCount<<" worker"<<(workerCount==1?"":"s")
             <<" (per-event wall-clock budget "<<wallClockBudgetSeconds<<" s).\n";
    const auto worker=[&]() {
        while(true) {
            const int index=nextIndex.fetch_add(1);
            if(index>=runCount) break;
            estimates[static_cast<size_t>(index)]=estimateCremCollapse(
                splitMix64(masterSeed+static_cast<std::uint64_t>(index)),
                selectedPhenomenon,wallClockBudgetSeconds);
            const int done=completed.fetch_add(1)+1;
            if(done%10==0||done==runCount) {
                std::lock_guard<std::mutex> lock(outputMutex);
                std::cout<<"CREM calibrations: "<<done<<"/"<<runCount<<'\n';
            }
        }
    };
    std::vector<std::thread> workers;
    workers.reserve(static_cast<size_t>(workerCount));
    for(int index=0;index<workerCount;++index) workers.emplace_back(worker);
    for(std::thread& thread:workers) thread.join();
    return estimates;
}
