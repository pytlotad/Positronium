#pragma once

// THE CONFIGURATION PANEL: the one place a run is steered from.
//
// Every switch a user may turn on or off is here, next to the reason its
// default is what it is.  Nothing below is a physical constant and nothing
// below is derived: these are the choices, and the rest of the project reads
// them.
//
// THREE WAYS TO STEER A RUN, and they compose in this order, each overriding
// the one before it:
//
//   1. THIS FILE.  Edit a default and rebuild.  This is where a configuration
//      you intend to KEEP belongs: it is version-controlled, a reader sees it
//      without running anything, and it applies to every run of that build.
//   2. THE COMMAND LINE.  Every knob below names its flag.  A flag overrides
//      this file for one run and changes nothing permanently.
//   3. THE IN-PROGRAM MENU.  Four questions -- external field, mode, visual
//      style, experiment.  Each is asked ONLY when neither this file nor the
//      command line has already answered it, which is why a fully specified
//      batch run is never blocked waiting on stdin.
//
// THE MAIN PATH is the configuration every committed result was produced
// with: it is exactly the set of defaults written below, and it is
//
//     electron + positron, CODATA g-factors, started at a_1 = 1 a_pair,
//     stochastic E1-dipole radiation reaction with deterministic emission,
//     spin quantized, NO ground-state floor, NO Bohr-ladder photon energy,
//     NO external field, NO zero-point field.
//
// Anything else is off the main path.  That is not a warning against using
// it -- several of the alternatives exist precisely to be measured against
// the default -- but a result obtained off it has to be reported as such,
// because the numbers quoted in the README and in the audit reports do not
// cover it.  Each knob below states what its non-default settings mean.
//
// WHAT IS DELIBERATELY NOT HERE, so that its absence is not read as
// completeness:
//
//   - PHYSICAL CONSTANTS AND SPECIES DATA are not configuration.  Masses,
//     charges and g-factors are CODATA measurements and live in
//     physical_constants.hpp and particle_species.hpp.  The one exception is
//     the g-factor override below, which exists so the moment can be varied
//     as an experiment; it is off by default and says so loudly.
//   - ENVIRONMENT-VARIABLE SWITCHES are CATALOGUED at the bottom of this file
//     but are not variables here.  They are read at call sites deep inside
//     the engine, and most are one-run probes rather than configurations.
//     The catalogue separates the ones that change a RESULT from the ones
//     that only print.
//   - THE VALIDATION EXECUTABLE DOES NOT READ THIS FILE'S SELECTIONS.
//     positronium_validation stays pinned to the built-in default pair
//     unless its own --pair says otherwise, because its checks carry
//     numeric thresholds measured for that pair.  Editing the panel changes
//     production runs; it must not silently move a regression threshold.

#include "electrodynamics.hpp"
#include "pair_configuration.hpp"
#include "particle_species.hpp"
#include "physical_constants.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

// ===================================================================
// GROUP A -- ENGINE SWITCHES
//
// These six sit at global scope rather than in the configuration namespace
// below, because crem_trajectory.hpp and crem_collapse.hpp read them by
// these names from inside the integration loop.  They are configuration all
// the same, and this is their home.
// ===================================================================

// Selected once from --radiation-reaction and read by every trajectory
// constructed below (visual, beam and interaction experiments alike).
// Defaults to stochasticElectricDipole (radiation ON): the same E1 dipole
// power individualLandauLifshitz would remove as a continuous drag force is
// instead banked as Poissonian hazard and paid out in discrete,
// momentum-conserving photon kicks -- see the enum's own comment in
// electrodynamics.hpp for why (deterministic drag has no photon to report;
// --emission deterministic bypasses the Poisson draw but still uses this
// channel).  estimateCremCollapse now measures the classical inspiral
// mechanically rather than assuming it, so it needs some such channel
// switched on to observe anything; --radiation-reaction individual restores
// the continuous-drag alternative.
// [[maybe_unused]] because the validation executable's main() never reaches
// the trajectory constructors that read this, so GCC sees no use in that
// build.
[[maybe_unused]] inline ChargeRadiationReactionModel gRadiationReactionModel =
    ChargeRadiationReactionModel::stochasticElectricDipole;

// Starting separation the bound scenarios are PREPARED at, a_n = n^2 a_pair
// (a_pair itself comes from the pair's measured magnetic moment -- see
// pairBohrRadius's own comment -- the one length scale this model actually
// has to start an inspiral from).  Not a claim that the pair occupies "the
// n-th energy level" as a physical eigenstate: with the ground-state floor
// and the Bohr-ladder photon energy both off by default (see
// gGroundStateEmissionFloor, gBohrLevelPhotonEnergy below), nothing in the
// default run treats this as a quantized level any more than any other
// starting radius would be -- it is just where the classical inspiral begins.
// n is still an integer and a_n still scales as n^2 purely so --level stays
// a convenient, backward-comparable way to pick a starting separation.
//
// SET TO 1, having been 2.  The old comment justified 2 by the
// instantaneous-kick emission ceiling on the Bohr ladder,
//
//     hbar*omega / E_kinetic = 2/n
//
// (virial: a circular orbit's kinetic energy equals its binding energy),
// which a photon of hbar*omega cannot satisfy at n=1 and fits from n=3 up.
// That argument does not distinguish 1 from 2: the ratio is exactly 1 at
// n=2, so the ceiling excludes n=2 as well.  It never favoured the value it
// was given as the reason for, and 2 was therefore paying a large price for
// nothing.  (The ceiling itself stands, and is documented where it belongs:
// the model never leaves n <= 1.09, so it binds over the whole operating
// range wherever the ladder concept is used at all.)
//
// The price, measured.  The starting radius is what sets the collapse time --
// the floor barely touches it, because t ~ a^3 makes everything below the
// floor cheap (measured: the ground-state floor moves the median by +0.6%,
// and a floor anywhere below a_pair/4 is worth under 2% of the full
// inspiral).  At n=2 the reported median is 6246 ps; at n=1 it is 148 ps,
// a factor 42 (the hazard's own n^5 predicts ~32; the excess is the part of
// the cascade that runs below n=1, where that scaling no longer holds).
//
// Why 1 is the better default under the "classical reproduction attempt"
// scope.  The measured para-Ps lifetime is 124.49 ps.  At n=1 this model
// reports 148 ps -- 19% high, and the closest it gets to any measured
// lifetime.  At n=2 it reports 6246 ps, which is 50x the measurement and is
// a statement about the chosen initial condition rather than about the pair.
// n=1 is also the more robust configuration, not merely the faster one:
// measured on seed 42 with no floor, 8 of 8 trajectories complete cleanly at
// n=1 against 7 of 8 at n=2.
//
// What n=1 does NOT mean.  It is still a starting SEPARATION, a_pair, taken
// from the pair's measured magnetic moment, not a claimed eigenstate -- and
// --ground-state-floor must stay OFF at this level, since the sharp
// preparation sits exactly ON the floor and the floor then declares every
// trajectory settled before it moves (measured: 8 of 8 at t = 0).
//
// a_n = n^2 a_pair and the tangential band is quoted in units of the circular
// speed AT that separation, so the sampled spread in L/(n hbar) is unchanged
// and only the level moves.
inline int gInitialPrincipalLevel = 1;

// Whether the quantized emission draws its next threshold from Exp(1) (a
// genuine Poisson process, available via --emission poisson) or fires
// deterministically as soon as one quantum's worth of energy has
// accumulated -- the default since the measurement below.
//
// The two differ by ONE constant, and the reason is worth stating because it
// makes the deterministic variant nearly free.  The banked hazard is
//
//     integral(rate dt) = integral(P/E_photon dt) = (radiated energy)/E_photon,
//
// i.e. the accumulated loss measured in quanta.  Firing at a fixed threshold
// of 1 therefore means "emit exactly when the orbit has continuously lost one
// photon's worth of energy" -- the orbit still descends continuously, and the
// level crossing itself is the trigger.  Firing at an Exp(1) threshold
// instead is what turns that into spontaneous emission with its shot noise.
//
// Both are production models, not probes.  They serve the same direction --
// a DETERMINISTIC determination of the quantum parameters -- and differ in
// how far along it they go.  The photon ENERGY is already fixed by the orbit
// in both (correspondence, or the level spacing where a ladder exists);
// deterministic emission fixes its TIMING too, leaving nothing about the
// quantum drawn at random.
//
// The departure this makes from the conventional description is worth
// stating plainly rather than hiding: spontaneous emission is normally
// modelled as Poissonian, and that used to be why Poisson was the default.
// What changed the answer is that the "sharp preparation" -- the pair now
// starts on an EXACT circular Bohr orbit (e=0), not a sampled eccentricity
// band -- removed the last non-shot-noise source of collapse-time spread.
// Re-measured after that change (N=100, seed 7, level 2 -> 1): Poisson
// sigma/mean = 0.883 (matches the old 0.894 within sampling noise), while
// deterministic sigma/mean = 4.0e-11 -- floating-point noise, not physics,
// because every trajectory now shares the identical circular Larmor rate and
// there is nothing left for the threshold draw to average over. The mean is
// preserved as claimed: 6398.9+/-568.1 ps (Poisson) against 6206.4 ps
// (deterministic), 0.34 standard errors apart. Poisson's entire reported
// spread is therefore shot noise around that same number, not a competing
// physical prediction, so deterministic is the better default estimator of
// the trajectory itself; --emission poisson remains available for whoever
// specifically wants the spontaneous-emission statistics.
inline bool gDeterministicEmission = true;

// The same standing as --zpf, and aimed at the same gap: CREM is a classical
// radiative inspiral with nothing to halt it at the pair Bohr radius, so left
// alone the collapse runs past the ground state and ends wherever the
// electrodynamics itself stops it -- the Compton barrier, or (measured, the
// large majority of the time -- see the README's Compton-barrier
// re-measurement) a numerical retardation-time safety margin first.
//
// --zpf tried to supply a stopping mechanism from outside, by coupling a
// classical zero-point field and looking for a fluctuation-dissipation
// balance.  That failed: the resonant band moves the collapse time by 0.3%,
// wider bands only pump the orbit toward escape, and the balance condition
// could not even be measured (see the README's ZPF section).
//
// This flag tries the opposite tack.  Instead of a mechanism, it imports ONE
// quantum fact and nothing else: the Bohr ladder terminates, so there is no
// state below n=1 for a photon to leave the pair in.  Emission is refused
// whenever it would bind the pair tighter than the ground state.
//
// It is a CLOSURE, not a derivation.  It does not explain why no lower state
// exists; it asserts it, and any result obtained with it has to be read that
// way.  What it buys over --zpf is that it carries no free parameter, no band
// edge to choose, no mode count to converge, and no cost.  --ground-state-floor
// restores it.
//
// OFF BY DEFAULT, on request, to leave the classical electrodynamics running
// unmodified: this closure, gBohrLevelPhotonEnergy below, and the choice of
// L=n*hbar starting separation (gInitialPrincipalLevel, now documented as a
// starting radius rather than a claimed energy level) are the three places
// this file imports a discrete quantum fact rather than deriving one, and
// this is the first to go.  Measured (N=100x2, seed 7, --no-ground-state-floor
// equivalent, both channels): 97-98% of trajectories stop on the retardation
// safety margin, 2-3% on the Compton barrier, "landed at periapsis" spans a
// smooth 0.26-15.3 r* with no clustering, and para/ortho are statistically
// indistinguishable -- i.e. what the bare electrodynamics produces here is
// continuous scatter, not a rediscovered ladder.
//
// WHAT THE REPORTED TIME MEANS WITH THE FLOOR OFF: a classical inspiral time
// to wherever the electrodynamics actually stops (Compton barrier or
// retardation margin), NOT a cascade to n=1 and NOT an annihilation
// lifetime -- this model has no annihilation dynamics at all, no contact
// channel and no rate, established separately.  --ground-state-floor restores
// the n=1 cascade-time reading documented above.
inline bool gGroundStateEmissionFloor = false;

// SPIN QUANTIZATION, split out of --ground-state-floor and defaulted ON.
//
// S=0 and S=1 are exact states, so the two moments are exactly aligned or
// exactly anti-aligned -- never somewhere in a band.  Opposite charges invert
// the spin-moment relation, so ANTI-parallel spins (para, S=0) give ALIGNED
// moments and |mu1+mu2| = 2mu, while ortho (S=1) gives |mu1+mu2| = 0 exactly.
// That is what makes the coherent M1 channel a real para/ortho asymmetry
// rather than a label.
//
// Sampled from a BAND instead -- para drawn from cos >= 0.5, ortho from
// everything below -- the asymmetry is destroyed, and measurably so, which is
// why this is now the default rather than an option.  Measured over the
// sampled configurations at --level 1: |mu1+mu2| came out 1.775-1.942 mu for
// para against 1.271-1.716 mu for ortho, i.e. ortho carrying 63-86% of para's
// net moment instead of zero, and the resulting M1 shares were comparable
// (ortho's largest, 9.0e-13 of E1, exceeded para's 2.5e-13).  The band also
// puts para's lower edge exactly ON the classification threshold, so half the
// para trajectories precessed across it during the inspiral (6 of 12).
//
// It was previously reachable only through --ground-state-floor, which
// ALSO zeroes the emission hazard at n=1 and therefore stops every trajectory
// before it reaches the collision boundary (measured: 0 of 4 collapses).  The
// two are independent physical assertions and are now independent switches:
// this one fixes the initial mutual angle, that one closes the ladder from
// below.  --no-spin-quantization restores the band sampling.
inline bool gSpinQuantization = true;

// The second of the three imported quantum facts (see
// gGroundStateEmissionFloor's comment).  quantumFor (crem_collapse.hpp) needs
// a photon energy for the secular estimator's hazard bookkeeping; it has
// always had two candidates available, the Bohr LEVEL DIFFERENCE E(n)-E(n-1)
// where the ladder has a rung below, and hbar*omega_orb -- the value the
// orbit's own frequency actually produces, with no ladder concept at all --
// as the fallback everywhere the ladder does not reach (n<2, or now,
// unconditionally, whenever this flag is off).
//
// hbar*omega_orb is not a worse number by construction: it is the
// correspondence-principle value, and measured, dE(n->n-1)/hbar*omega_orb is
// close to 1 (1.0152 at n=100, 1.0523 at n=30) everywhere except close to the
// ground state, where the ladder's spacing stops resembling the local orbital
// frequency at all (a factor of 3 off at n=2) -- see quantumFor's own comment
// for the measurement.  Choosing it unconditionally means the photon energy
// this file reports is always something the electrodynamics itself produces,
// never a level difference imported from the quantum ladder.
//
// OFF BY DEFAULT, on request, alongside gGroundStateEmissionFloor above.
// --bohr-photon-energy restores the level-difference rule (quantumFor's own
// n>=2 branch), for comparison against the historical behaviour.
// [[maybe_unused]] because quantumFor lives in crem_collapse.hpp, which is
// only included outside the validation executable (see its #ifndef just
// below).
[[maybe_unused]] inline bool gBohrLevelPhotonEnergy = false;

// ===================================================================
// GROUP B -- EVERYTHING ELSE
//
// Namespaced, because these are read once at startup rather than from inside
// the engine.  main() copies each into a local, the command line then
// overrides that local, and the menu asks only for what is still unanswered.
// ===================================================================

namespace configuration {

// -------------------------------------------------------------------
// B1.  WHICH PAIR OF PARTICLES THE RUN INTEGRATES
// -------------------------------------------------------------------

// "first,second", the same spelling --pair takes.  Selectable species are
// electron, positron, muon, antimuon, proton and antiproton, and the pair
// must ATTRACT and carry opposite unit charges -- applyPair rejects anything
// else, because a repelling pair has no bound states and every experiment
// here is built around capture and inspiral.
//
// The pair is not a label.  It sets the reduced mass, the pair Bohr radius
// a_pair (derived from the two measured magnetic moments, not imported), the
// binding energy, the dipole regularization radius and the collision
// boundary, all of which follow the pair rather than staying at
// positronium's numbers: 6.80 eV of binding for e+e-, 1.41 keV for true
// muonium, 12.5 keV for protonium.
//
// OFF THE MAIN PATH: anything but electron,positron.  The annihilation
// comparisons in particular are positronium's measured lifetimes and do not
// transfer -- isPositronium() gates them for exactly this reason.
inline std::string pair = "electron,positron";

// -------------------------------------------------------------------
// B2.  THE g-FACTOR
// -------------------------------------------------------------------

// Zero means "use the species' own CODATA value", which is the main path and
// what these should stay at for any run meant to describe real particles.
//
// g is NOT a free parameter of this model.  It is a measured property of the
// species, and it is emphatically not close to 2 for every one of them --
// the proton carries 5.5857 because it is composite -- so nothing anywhere
// may assume the Dirac value.
//
// WHAT AN OVERRIDE ACTUALLY MOVES, measured rather than assumed.  The
// intrinsic moment is mu = (g/2)(|q|hbar/2m), so it scales with g directly,
// and with it everything that reads a moment: the precession rate, the
// mutual dipole-dipole force, the coherent M1 channel, the motional electric
// dipole, and the dipole regularization radius (which goes as the cube root
// of mu1 mu2 -- g1 = 4.0 on the electron measured 8.364e-14 -> 1.053e-13 m).
//
// WHAT IT DOES NOT MOVE, and this is worth knowing before reading a scan:
// the pair Bohr radius and therefore the binding energy and the collision
// boundary.  a_pair = 16 (m1+m2) mu1 mu2 / (g1 g2 |q1 q2|^2 K) contains the
// moments only as mu/g, which is the magneton |q|hbar/2m, so the g's cancel
// identically and the expression reduces to the Bohr relation
// hbar^2/(mu_red |q1 q2| K).  Measured: g1 = 4.0 leaves a_pair at
// 1.0583544e-10 m and the binding at 6.8028466 eV, unchanged to every digit.
// An override is therefore a probe of the MAGNETIC sector alone, not a way
// to rescale the orbit.
//
// OFF THE MAIN PATH: any positive value.  Set one only to answer "how does
// this result depend on the moment", and report it as that.
inline double firstGFactorOverride = 0.0;
inline double secondGFactorOverride = 0.0;

// Resolve B1 and B2 into the role-scoped globals every force law reads.
// Called once from main() after the command line has been parsed and before
// anything integrates or reads a pair-derived scale.
inline void applySelections() {
    const std::size_t comma = pair.find(',');
    if(comma == std::string::npos)
        throw std::invalid_argument("configuration::pair needs two species "
            "separated by a comma, e.g. proton,electron");
    const ParticleSpecies* first = speciesByName(pair.substr(0, comma));
    const ParticleSpecies* second = speciesByName(pair.substr(comma + 1));
    if(!first || !second)
        throw std::invalid_argument("unknown species in configuration::pair; "
            "known species: " + selectableSpeciesList());
    ParticlePair selected{*first, *second};
    // Applied to the COPY, before applyPair derives anything, so every scale
    // that depends on the moment is built from the overridden value rather
    // than being patched afterwards and left inconsistent.
    if(firstGFactorOverride > 0.0) selected.first.gFactor = firstGFactorOverride;
    if(secondGFactorOverride > 0.0) selected.second.gFactor = secondGFactorOverride;
    applyPair(selected);
}

// -------------------------------------------------------------------
// B3.  WHERE THE PAIR STARTS
// -------------------------------------------------------------------
//
// The starting SEPARATION is gInitialPrincipalLevel in group A above,
// a_n = n^2 a_pair; --level sets it and its own comment explains at length
// why n = 1 is the default and what it does and does not claim.
//
// The starting VELOCITY is exactly circular at that separation (e = 0), a
// sharp preparation with no sampled spread.  Restoring the old sampled band
// -- radial in [-0.10, +0.10] and tangential in [0.88, 1.12] of the circular
// speed -- is CREM_INITIAL_BAND=1 in the catalogue below.  That band was
// what deterministic emission's near-zero collapse-time spread depended on
// removing, so the two interact; see gDeterministicEmission's comment.
//
// The pair starts at exactly zero total momentum.  CREM_COM_DRIFT gives it a
// common velocity normal to the orbital plane; read that entry before
// reporting anything from it, because the secular estimator projects the
// drift straight back out.

// -------------------------------------------------------------------
// B4.  RUN IDENTITY: SEED, MODE, EXPERIMENT
// -------------------------------------------------------------------

// Every random draw a run makes comes from this one seed, including the
// external field's orientation and the zero-point field's mode directions,
// so a stated seed reproduces the whole run and not merely the trajectory.
// Left false, the seed is drawn from std::random_device and printed, which
// is the right default for a batch that wants independent events; set it
// true to repeat one exactly.  --seed sets both at once.
inline bool useFixedSeed = false;
inline std::uint64_t fixedSeed = 42;

// 0 asks at startup, 1 = visual, 2 = statistical.  --mode visual|statistical.
// Visual integrates ONE prepared trajectory and draws it; statistical is
// batch-only, never opens a window, and is where every reported number comes
// from.
inline int mode = 0;

// 0 asks at startup.  Which experiment, and the meaning depends on the mode:
//
//   statistical  1 para-Ps CREM collapse + 2 gamma kinematics
//                2 ortho-Ps CREM collapse + 3 gamma kinematics
//                3 e+e- beam, short-range/cutoff channel
//                4 e+e- beam, elastic scattering
//                5 interactions, classify collision outcomes
//   visual       1 para-positronium   2 ortho-positronium
//                3 direct collision   4 scattering
//
// Experiment 5 exists only in statistical mode: visual integrates a single
// trajectory, whereas 5 classifies an ensemble.  --phenomenon.
inline int phenomenon = 0;

// Visual mode only.  0 asks at startup, 1 = line, 2 = dot.  --visual-style.
inline int visualStyle = 0;

// Visual mode only: integrate the trajectory and print the diagnostic dump
// instead of opening a canvas.  Implies visual mode, so a run that also
// names a phenomenon is fully specified and asks nothing.  --diagnose.
inline bool diagnose = false;

// One sample size for every statistical experiment.  The per-experiment
// overrides that used to sit here made the default depend on which channel
// was selected, which is a poor property for a number that appears on every
// plot as "N =".
//
// Ceilings differ because the cost does: experiments 1 and 2 mechanically
// resolve every trajectory to the collision boundary and are capped at 1000,
// while 3, 4 and 5 sample and are capped at 100000.  --runs.
inline int statisticalRuns = 1000;

// -------------------------------------------------------------------
// B5.  BEAM EXPERIMENTS 3 AND 4
// -------------------------------------------------------------------

// Unrelated to the interaction energy in B6 below: the beam channel is a
// scattering measurement whose reported cross sections scale as 1/K_CM^2, so
// this value fixes the axis range of every committed beam plot and must not
// be retuned along with the experiment-5 energy.  --beam-energy-ev.
inline double beamEnergyEv = 20.0;

// Zero by default: unlike experiment 5, experiments 3/4 exist to compare a
// measured cross section against the theoretical Rutherford formula at a
// FIXED energy, and sigma(theta) is itself defined at one K_CM -- a
// zero-width beam is the right default for that comparison, not a gap.
// Nonzero opts into modelling a real beam's finite energy resolution
// instead, at the cost of smearing the theory-vs-model comparison the same
// way real energy spread would smear a real measurement.
// --beam-energy-sigma-ev.
inline double beamEnergySigmaEv = 0.0;

// Smallest scattering angle the analysis bins.  Rutherford diverges as
// theta -> 0, so some floor is required; this one also sets the automatic
// impact-parameter window below.  --theta-min-deg.
inline double thetaMinimumDegrees = 5.0;

// Histogram bins across the analysed angular range.  --angle-bins.
inline int angleBins = 10;

// Widest impact parameter sampled, in picometres.  ZERO SELECTS THE
// AUTOMATIC RULE, which is what the committed plots use: 12x the cutoff
// impact parameter for the short-range channel, else 1.25 Coulomb lengths
// over tan(theta_min/2).  Stating a number instead widens or narrows the
// sampling window, which is a statistics-tuning choice -- how wide to cast
// the net for good counts across every outcome bucket -- and is independent
// of where a collision is declared.  Widening it dilutes every bucket at
// fixed N, not just the one being chased.  --bmax-pm.
inline double impactParameterMaximumPm = 0.0;

// Radius at which the trajectory is matched onto its asymptotic Coulomb
// form, in picometres.  ZERO SELECTS THE AUTOMATIC RULE, the larger of 100x
// bmax and 1000x the Coulomb strength over the energy.  It must exceed bmax
// or the match has nothing to match onto.  --matching-radius-pm.
inline double matchingRadiusPm = 0.0;

// -------------------------------------------------------------------
// B6.  INTERACTION EXPERIMENT 5
// -------------------------------------------------------------------

// Mean and width of the sampled relative mechanical energy, in eV.  A mean
// of 0.6 eV sits well below the 6.8 eV Ps binding energy, so the pair has to
// shed correspondingly less energy to bind and the captured fraction rises.
// --interaction-energy-ev, --interaction-energy-sigma-ev.
inline double interactionEnergyEv = 0.6;
inline double interactionEnergySigmaEv = 0.4;

// Width of the sampled impact parameter b, in picometres.  Fixed in absolute
// terms, NOT tied to the Coulomb length: l_C scales as 1/K_CM, so an auto
// width would widen faster than the capture threshold and lowering the
// energy would produce FEWER bound states, not more.  Zero still selects the
// l_C rule for anyone who wants that scaling.  --interaction-bsigma-pm.
inline double interactionImpactSigmaPm = 60.0;

// -------------------------------------------------------------------
// B7.  HOW LONG ONE CREM TRAJECTORY MAY RUN
// -------------------------------------------------------------------

// Experiments 1 and 2 integrate the full mechanical trajectory to the
// collision boundary instead of extrapolating; a bound orbit near a0 needs a
// huge number of cheap orbits before the radiative loss becomes visible, so
// each event is CENSORED once it spends this long on the wall clock rather
// than left to run indefinitely.
//
// Raised from 20 s with the n=2 cascade default: the cascade takes about
// 7 ns of simulated time against the barrier-limited inspiral's ~150 ps, and
// 45 s per trajectory measured 14 completions out of 16.  90 s leaves margin
// without making a standard batch open-ended.
//
// This is a compute budget, not physics, but it is not free of consequence:
// a censored trajectory enters the Kaplan-Meier estimate as a right-censored
// observation, and a run where everything censors reports a restricted mean
// that identifies neither the median nor the unrestricted mean.  The output
// says so when it happens.  --crem-wallclock-budget-s.
inline double cremWallClockBudgetSeconds = 90.0;

// -------------------------------------------------------------------
// B8.  UNIFORM EXTERNAL MAGNETIC FIELD
// -------------------------------------------------------------------

// Microtesla.  NEGATIVE MEANS "not stated", which is what lets the startup
// question stay silent for a fully specified batch run instead of blocking
// it on stdin; the question offers 0 or Earth's 50 uT.  Setting a value here
// answers it permanently and the question is never asked.
//
// The orientation is isotropic and drawn from the run seed, so it is part of
// what a stated seed reproduces rather than a hidden extra input.
//
// Worth knowing before reading the output: at 50 uT the cyclotron rate eB/m
// is 8.8e6 rad/s against an orbital rate near 3e15 rad/s, so the orbit
// itself is untouched at the ninth decimal.  What the field does reach is the
// dipoles, which precess at the same 8.8e6 rad/s, about 3e-4 rad over a
// 35 ps collapse.  The effect is real, small, and mostly magnetic.
//
// OFF THE MAIN PATH: any nonzero value.  Zero is the historical behaviour,
// the model's own list of excluded effects named external fields, and every
// committed result was produced without one.  --external-field.
inline double externalFieldMicroTesla = -1.0;

// -------------------------------------------------------------------
// B9.  OPTIONAL CLASSICAL ZERO-POINT FIELD
// -------------------------------------------------------------------

// Amplitude of the stochastic-electrodynamics zero-point field, in units of
// the physical level, so 1 is the real thing and the absorbed power scales
// as the square.  OFF BY DEFAULT and firmly off the main path: it is an
// experiment, not part of the model every committed result was produced
// with.
//
// What it is: the FLUCTUATING half of a Langevin pair whose dissipative half
// (radiation reaction) is already present.  It is not a drag -- the omega^3
// spectrum is the unique boost-invariant one, and a drag would single out a
// rest frame.
//
// What it was tried for and what happened: it was meant to supply the
// stopping mechanism a classical radiative inspiral lacks, through a
// fluctuation-dissipation balance.  That failed.  The resonant band moves
// the collapse time by 0.3%, wider bands only pump the orbit toward escape,
// and the balance condition could not be measured at all.  Read the README's
// ZPF section before reporting anything from this.
//
// Anything other than 1 is a numerical experiment on top of an experiment.
// --zpf.
inline double zeroPointScale = 0.0;

// Band edges in units of the pair's OSCULATING orbital angular frequency, so
// the whole band rides up with the orbit as it tightens and stays in
// resonance instead of being left behind.  With a fixed band the orbit
// outran it: the measured period falls from 0.327 fs to 0.0031 fs over a
// collapse, a factor of 105, which is why a fixed-band lifetime depended so
// strongly on where the upper edge was cut.
//
// The upper edge is set by what the trajectory integrator can RESOLVE, not
// by physics.  The step follows the trajectory error probe, which watches the
// orbit and not the field, so a wider band costs a proportionally shorter
// step and, past some width, is sampled only a few times per cycle.  That is
// aliasing, and it fakes exactly the effect a wide band appears to report:
// net work from modes that should average to nothing.  CREM_STEP_CENSUS
// prints steps-per-cycle for the fastest mode; use it whenever this band is
// widened.  --zpf-band lo,hi.
inline double zeroPointBandLow = 0.3;
inline double zeroPointBandHigh = 3.0;

// A convergence knob, not physics: the real spectrum is a continuum and
// these are a finite sampling of it.  Modes carry equal energy, so the
// per-mode amplitude falls as 1/sqrt(N) and the total does not move.
// --zpf-modes.
inline int zeroPointModes = 64;

}  // namespace configuration

// ===================================================================
// GROUP C -- ENVIRONMENT-VARIABLE SWITCHES (catalogue)
//
// Set in the environment, read at the call site, and NOT variables in this
// file: several are one-run probes wired deep inside the integration loop.
// They are listed here so that this file remains the complete answer to
// "what can be turned on or off".  All are off unless set.
//
// C1.  THESE CHANGE A RESULT.  A number produced with any of them set is off
//      the main path and must be reported that way.
//
//   CREM_NO_DIPOLE_FORCE   Ablate the mutual dipole-dipole FORCE while the
//                          moments still precess, still radiate M1 and still
//                          enter the annihilation invariant.  Exists to test
//                          whether that force carries the para/ortho channel
//                          difference, rather than asserting that it does.
//                          (electrodynamics.hpp, gDipoleForceEnabled.)
//   CREM_INITIAL_BAND=1    Restore the sampled initial-velocity band in place
//                          of the sharp circular preparation.  See B3.
//   CREM_COM_DRIFT=<beta>  Give the whole pair a common velocity, in units of
//                          c, normal to the orbital plane.  The secular
//                          estimator is frame-locked to the zero-momentum
//                          frame and projects this straight back out, so a
//                          null collapse-time result here is a statement
//                          about the estimator's scope and NOT about physics.
//   CREM_LADDER_BELOW_2    Let the Bohr-ladder photon energy fire below n=2.
//                          Only meaningful together with
//                          gBohrLevelPhotonEnergy, whose branch otherwise
//                          cuts at n >= 2.
//   CREM_MAX_DEPTH=<n>     Adaptive-step recursion ceiling (default 12).
//   CREM_RETARDATION_LIMIT=<x>
//                          Lower the period/light-crossing safety margin from
//                          its default of 150 so the cascade can be followed
//                          past where it normally stops.  Added for one
//                          measurement, the first passage of L through the
//                          contact value, which happens one photon after the
//                          default stop.  Below 150 the pair ends up inside
//                          the Compton barrier, where classical
//                          point-particle electrodynamics does not apply, and
//                          the retarded-field reconstruction is not validated
//                          either: a probe of the model's own bookkeeping,
//                          never a physical claim.
//   CREM_HARMONIC          Opt out of the harmonic correction in the secular
//                          estimator.
//   CREM_FORCE_M1=<x>      Force the magnetic-dipole share of radiated power
//                          to a stated value instead of the computed one.
//   CREM_AXIAL_SPIN        Constrain the spin to the orbital axis.
//   CREM_SPIN_MAGNITUDE    Hold the spin magnitude fixed under transport.
//   CREM_SPIN_TRIM         Renormalize the spin after each transport step.
//   CREM_PAIR_L_WITH_E     Recompute angular momentum on the same schedule as
//                          energy across a skipped span; this is the change
//                          measured to remove 90% of the L/E separation.
//   CREM_DARWIN_FORCES     Validation build: use Darwin forces where the
//                          retarded ones are the default.
//
// C2.  THESE ONLY PRINT.  Diagnostics and censuses; they do not alter the
//      trajectory, and are safe to leave on while measuring.
//
//   CREM_PROGRESS, CREM_STEP_CENSUS, CREM_SKIP_CENSUS, CREM_QUANTUM_CENSUS,
//   CREM_CAUSALITY, CREM_PHOTON_BALANCE, CREM_M1_SHARE, CREM_TILT,
//   CREM_APSIDAL, CREM_ENERGY_SPLIT, CREM_EMISSION_REACH, CREM_L_UPDATE,
//   CREM_DEBUG, CREM_DEBUG_PRECISE, CREM_DEBUG_CAPTURE, CREM_DEBUG_ORDER,
//   CREM_DEBUG_ORDER_LIVE, CREM_DEBUG_ORDER_TERMS, CREM_DEBUG_AZIMUTH,
//   CREM_DEBUG_CHANNEL, CREM_DEBUG_M1, CREM_DEBUG_GRAD, CREM_DEBUG_PROBES,
//   CREM_DEBUG_RETREAT, CREM_DEBUG_EPS_SCAN, CREM_DEBUG_FIELDSYM,
//   CREM_DEBUG_ALIGN, POSITRONIUM_DEBUG_DIPOLE, POSITRONIUM_DEBUG_FIELDS.
//
// ===================================================================
// GROUP D -- OPTIONS THAT NO LONGER EXIST
//
// Kept documented rather than silently unknown, because an old command line
// that still carries one should say what happened rather than fail with
// "unknown option".  Each of these is rejected with its own message:
//
//   --maxwell-test       Maxwell validation moved to ./positronium_validation.
//   --decay-events       Removed: the photon panels are exact reference
//                        curves and no longer sample the generator.
//   --stat-window-ps     Removed: the CREM calibration window is fixed by the
//                        orbit-averaged collapse estimator.
//   --integrator-order   Removed: the retarded, dissipative production solver
//                        uses its validated adaptive symmetric second-order
//                        step, and there is no longer a choice.
//   --no-gui             Accepted and inert.  Statistical mode is always
//                        batch-only and never opens a window.
// ===================================================================
