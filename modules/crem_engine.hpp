#pragma once

// CREM numerics: reconstruction of a causal retarded history for a freshly
// prepared state, and the adaptive integrator that advances a trajectory with
// it.  This is the step-size machinery -- error probe, subdivision, history
// retention -- as opposed to the force laws in electrodynamics.hpp.
//
// Unlike the other two CREM headers this one sits OUTSIDE the production
// #ifndef: positronium_validation drives the same engine.  Contains no ROOT.
//
// Self-contained and order-independent.  It names what it needs through a
// using-directive on positronium::parameters and using-declarations for the
// object types, rather than reopening namespace positronium: the header is
// still textually included inside positronium.cpp's anonymous namespace,
// where reopening a named namespace would create {anonymous}::positronium and
// hide the real one from every later lookup.

#include "electrodynamics.hpp"
#include "pair_geometry.hpp"
#include "physical_constants.hpp"
#include "retarded_charge_kinematics.hpp"
#include "state.hpp"
#include "state_validity_interpolation.hpp"
#include "vector3.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <iostream>
#include <limits>

namespace two_body = positronium::kinematics;

using positronium::objects::Vec3;
using positronium::objects::State;
using positronium::objects::StateHistory;
using positronium::objects::DipoleTensor;
using positronium::objects::cross;
using positronium::objects::dot;
using namespace positronium::parameters;

// --- Time regularization -------------------------------------------------
//
// A Kepler orbit spends most of its TIME near apoapsis and most of its
// ACCELERATION near periapsis, so a constant dt can only resolve the bottom
// of a deep orbit by paying for the top at the same rate.  Every stepping
// site here instead sets dt from the current separation, a Sundman-type
// regularization:
//
//   localOrbit     dt = 2 pi / (N omega_local) with omega_local =
//                  sqrt(k/(mu r^3)), so dt ~ r^{3/2}: N steps per LOCAL
//                  orbital period.  This is the production law and the
//                  default, and the expression is the one the Experiment 42
//                  loop used inline before it moved here.
//   constantAngle  dt = 2 pi mu r^2 / (N L) = (2 pi / N) / (dphi/dt), i.e. N
//                  steps per full turn of the true anomaly, so the swept
//                  angle per step stays constant all the way down.
//
// What the measurement actually says (audit section 89, apoapsis 3 r*,
// periapsis 0.3 r*, tolerance 1e-8) is narrower than the usual claim for
// regularization, and worth stating exactly, because it is easy to assume
// this buys accuracy or speed and it buys neither.  All three laws -- constant
// dt included -- reach the SAME trajectory: the same periapsis 0.2631 r*, the
// same apoapsis 45 r* afterwards, the same final L of 0.0421 hbar.  Nor is any
// of them appreciably cheaper: 16915 force evaluations for localOrbit, 19021
// for constantAngle, 17385 for constant dt, a spread of 12%.
//
// What regularization buys is SUBDIVISION DEPTH.  At maximumDepth 20 both
// regularized laws carry the window through and constant dt breaks off during
// the first plunge, at 0.523 of the window; constant dt needs depth 26 to
// finish the same passage.  Six levels, a factor of 64 in the smallest step
// the engine can reach, is the whole of the difference.  That is also the
// warning: the binding resource at the bottom of the orbit is the subdivision
// budget, not the outer step, and no choice of law removes the need for
// maximumDepth 20 on this geometry.
//
// THOSE NUMBERS PREDATE AUDIT 103 and were taken with the magnetization term
// centred on the retarded position.  With the corrected field the same
// geometry dips to 0.2095 r*, not 0.2631, and localOrbit at N = 256 needs
// maximumDepth 26 to finish three windows; the periapsis agrees at
// tolerances 1e-8 and 1e-6 (20992 and 2311 force evaluations).  The
// comparison between the three laws has not been repeated on it.
enum class TimeRegularizationLaw { localOrbit, constantAngle };

struct RegularizedStep {
    double stepsPerLocalOrbit=128.0;   // N above
    double ceiling=5.0e-18;            // s, absolute cap on one outer step
    TimeRegularizationLaw law=TimeRegularizationLaw::localOrbit;
};

// CREM_STEP_LAW=angle switches the production loop to the constant-angle law.
// Read once: this sits in the innermost stepping loop.
inline TimeRegularizationLaw configuredTimeRegularizationLaw() {
    static const TimeRegularizationLaw configured=[]{
        const char* requested=std::getenv("CREM_STEP_LAW");
        return requested&&std::strcmp(requested,"angle")==0
            ? TimeRegularizationLaw::constantAngle
            : TimeRegularizationLaw::localOrbit;
    }();
    return configured;
}

// The outer step for one advance, clipped to the budget still owed.  Callers
// keep their own guard on a non-finite or non-positive result: `remaining`
// can legitimately arrive at zero at the end of a window.
inline double regularizedTimeStep(const State& s,double remaining,
                                  const RegularizedStep& rule={}) {
    const double r=separation(s);
    const double reducedMass=firstMass*secondMass/(firstMass+secondMass);
    const double omega=std::sqrt(pairCoulombStrength/(reducedMass*r*r*r));
    double regularized=2.0*pi/(rule.stepsPerLocalOrbit*omega);
    if(rule.law==TimeRegularizationLaw::constantAngle) {
        const double orbitalAngularMomentum=reducedMass
            *cross(s.firstPosition-s.secondPosition,
                   s.firstVelocity-s.secondVelocity).norm();
        // The guard is for the exactly-radial orbit, where a swept angle
        // measures nothing and this law would hand back infinity.  Do not
        // mistake it for protection: a numerical orbit going radial never
        // lands on L == 0 exactly.  Measured on the tilted deep geometry
        // (audit 89e), L falls to 1.7e-4 hbar, this branch fires zero times
        // in 1335 steps, and dt ~ r^2/L meanwhile grows to 44x the
        // local-orbit step -- the law LENGTHENS the step exactly while the
        // orbit is collapsing.  That is why constantAngle is a diagnostic.
        if(orbitalAngularMomentum>0.0)
            regularized=2.0*pi*reducedMass*r*r
                       /(rule.stepsPerLocalOrbit*orbitalAngularMomentum);
    }
    return std::min({rule.ceiling,regularized,remaining});
}

inline StateHistory causalInitialHistory(const State& initial,double spanFactor=8.0,
                                  int intervalCount=64,
                                  int picardIterations=2) {
    State endpoint=initial;
    synchronizeCovariantDipoles(endpoint);
    const double lightCrossingTime=separation(endpoint)/c;
    const double historySpan=std::max(1.0e-24,spanFactor*lightCrossingTime);
    intervalCount=std::max(intervalCount,2);
    const auto buildHistory=[&](const Vec3& firstAcceleration,
                                const Vec3& secondAcceleration) {
        StateHistory result;
        for(int index=intervalCount;index>=0;--index) {
            const double offset=-historySpan*index/intervalCount;
            State sample=endpoint;
            sample.time=initial.time+offset;
            sample.firstPosition=initial.firstPosition
                +initial.firstVelocity*offset
                +firstAcceleration*(0.5*offset*offset);
            sample.secondPosition=initial.secondPosition
                +initial.secondVelocity*offset
                +secondAcceleration*(0.5*offset*offset);
            sample.firstVelocity=initial.firstVelocity
                +firstAcceleration*offset;
            sample.secondVelocity=initial.secondVelocity
                +secondAcceleration*offset;
            sample.firstAcceleration=firstAcceleration;
            sample.secondAcceleration=secondAcceleration;
            // Bookkeeping starts at the requested t=0 event, not in the
            // hidden causal preparation interval.
            sample.radiatedEnergy=initial.radiatedEnergy;
            sample.orbitalRadiatedEnergy=initial.orbitalRadiatedEnergy;
            sample.radiatedMomentum=initial.radiatedMomentum;
            sample.radiatedAngularMomentum=initial.radiatedAngularMomentum;
            synchronizeCovariantDipoles(sample);
            result.push_back(sample);
        }
        return result;
    };
    const MutualForces seedForces=allExternalForces(endpoint);
    Vec3 firstAcceleration=relativisticAcceleration(
        endpoint.firstVelocity,seedForces.first,firstMass);
    Vec3 secondAcceleration=relativisticAcceleration(
        endpoint.secondVelocity,seedForces.second,secondMass);
    StateHistory history=buildHistory(firstAcceleration,secondAcceleration);
    // Two inexpensive Picard updates make the hidden past consistent with
    // the same retarded interaction used at t=0 instead of freezing the
    // instantaneous Coulomb acceleration into the whole preparation span.
    picardIterations=std::max(picardIterations,0);
    for(int iteration=0;iteration<picardIterations;++iteration) {
        const MutualForces retarded=retardedExternalForces(endpoint,history);
        firstAcceleration=relativisticAcceleration(
            endpoint.firstVelocity,retarded.first,firstMass);
        secondAcceleration=relativisticAcceleration(
            endpoint.secondVelocity,retarded.second,secondMass);
        history=buildHistory(firstAcceleration,secondAcceleration);
    }
    return history;
}

// CREM_STEP_CENSUS accumulators -- see the acceptance site inside
// advanceAdaptive for what these are for and why they are unsynchronised.
inline bool gStepCensusEnabled=std::getenv("CREM_STEP_CENSUS")!=nullptr;
inline unsigned long long gStepCensusCount=0;
inline double gStepCensusTotalTime=0.0;
inline double gStepCensusSmallest=std::numeric_limits<double>::infinity();
inline double gStepCensusLargest=0.0;

class ClassicalTrajectoryEngine {
public:
    struct Accuracy {
        double relativeTolerance = 1.0e-7;
        int maximumDepth = 14;
        ChargeRadiationReactionModel reactionModel =
            ChargeRadiationReactionModel::individualLandauLifshitz;
        // Far-zone Poynting quadrature: 50 directions, each solving a
        // retarded time for both charges, on every fine half-step.  It feeds
        // only the flux BOOKKEEPING (radiatedEnergy/Momentum/AngularMomentum,
        // boundField*, reaction mismatches).  The trajectory itself never
        // reads any of them: forces, chargeReaction, magneticDipoleFlux and
        // hence dipoleConstraintEnergy are all computed outside that block,
        // so switching it off leaves positions, velocities and dipoles
        // bit-identical.  Callers that do not report a radiated energy can
        // therefore turn it off outright.
        bool computeOutwardFlux = true;
        // Validation-only diagnostic switch.  Production keeps the default
        // retarded mutual fields; false isolates the conservative
        // Coulomb-Darwin force while leaving the selected reaction and the
        // independently measured far flux unchanged.
        bool useRetardedExternalForces = true;
        // Angular resolution and control radius of the far-zone Poynting
        // quadrature that computeOutwardFlux above switches on and off.
        // It used to be frozen at the FarFieldSampling
        // defaults -- 50 directions, degree 11 -- with no way to reach it
        // from here, and 50 directions is not enough once the source
        // drifts.  The Lienard distribution carries (1 - n.beta)^-k up to
        // k = 6, whose dynamic range over the sphere is
        // ((1+beta)/(1-beta))^k: 80 at beta = 0.35 but 729 at 0.50 and
        // 3.3e+04 at 0.70.  Measured against the exact integral, the
        // 50-direction rule's relative error on that factor is 1.0e-06,
        // 5.7e-05 and 4.7e-04 at those three speeds, and raising the rule
        // to 194 directions cut a measured boost-covariance residual at
        // beta = 0.50 by a factor 22, from 6.8e-05 onto the 3.0e-06 floor
        // (audit sections 219a-219b).  The default is unchanged, so every
        // existing caller keeps the behaviour it was tuned against; a
        // caller that drifts the pair above about beta = 0.4, or that
        // reads the radiated four-momentum rather than only the
        // trajectory, should raise directionCount to 194 or 302.
        FarFieldSampling farFieldSampling{};
    };
    explicit ClassicalTrajectoryEngine(const State& initial)
        :history_(causalInitialHistory(initial)) {}
    ClassicalTrajectoryEngine(const State& initial,Accuracy accuracy)
        :history_(causalInitialHistory(initial)),accuracy_(accuracy) {}
    ClassicalTrajectoryEngine(StateHistory history,Accuracy accuracy)
        :history_(std::move(history)),accuracy_(accuracy) {}
    bool advance(State& state,double dt) {
        if(!(dt>0.0)||!std::isfinite(dt)||!isFinite(state)
            ||!(accuracy_.relativeTolerance>=0.0)
            ||!std::isfinite(accuracy_.relativeTolerance)
            ||accuracy_.maximumDepth<0) return false;

        // Treat a proposed step as a transaction.  Recursive subdivision may
        // finish its first half before the second half discovers that even the
        // deepest allowed step misses the tolerance.  Build the accepted state
        // and history in local outputs, then commit both together only after
        // the complete interval succeeds.
        State accepted;
        StateHistory acceptedHistory;
        if(!advanceAdaptive(state,history_,dt,0,accepted,acceptedHistory)) {
            return false;
        }
        state=std::move(accepted);
        history_=std::move(acceptedHistory);
        return true;
    }
    const StateHistory& history() const { return history_; }
private:
    static double normalizedStepError(const State& coarse,const State& fine) {
        // Error control must follow the selected pair.  The former universal
        // nuclearCutoff floor exceeded a protonium orbit in the inner radius
        // sweep, so a nominal relative tolerance became much looser exactly
        // where the heavy-pair dynamics was hardest.  The pair's own terminal
        // surface is the smallest resolved length in this engine.
        const double lengthScale=std::max(
            separation(fine),collisionBoundaryRadius);
        const double speedScale=std::max(
            (fine.firstVelocity-fine.secondVelocity).norm(),1.0e-6*c);
        return std::max({
            (coarse.firstPosition-fine.firstPosition).norm()/lengthScale,
            (coarse.secondPosition-fine.secondPosition).norm()/lengthScale,
            (coarse.firstVelocity-fine.firstVelocity).norm()/speedScale,
            (coarse.secondVelocity-fine.secondVelocity).norm()/speedScale});
    }

    bool advanceAdaptive(const State& start,const StateHistory& history,
                         double dt,int depth,State& accepted,
                         StateHistory& acceptedHistory) {
        const double requestedEndTime=start.time+dt;
        if(!(requestedEndTime>start.time)||!std::isfinite(requestedEndTime)) {
            if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
                std::cerr<<"ENGINE_DEBUG reason=request-time dt="<<dt
                    <<" t="<<start.time<<" depth="<<depth<<'\n';
            return false;
        }
        const char* rejectionReason="accuracy";
        double rejectionMetric=std::numeric_limits<double>::quiet_NaN();
        const auto subdivide=[&]() {
            if(depth>=accuracy_.maximumDepth) {
                if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
                    std::cerr<<"ENGINE_DEBUG reason="<<rejectionReason
                        <<" metric="<<rejectionMetric<<" depth="<<depth
                        <<" dt="<<dt<<" t="<<start.time
                        <<" r="<<separation(start)<<'\n';
                // CREM_DEBUG_ORDER_LIVE: the same step-doubling sweep
                // CREM_DEBUG_ORDER runs on a synthetic state, but on the REAL
                // one that is failing, with its REAL retarded history, and
                // with the force sectors switched out the same way.
                //
                // This is what identified the ortho collapse's one numerical
                // failure, and it had to be run here rather than on a
                // synthetic state because the synthetic probe gave a FALSE
                // NEGATIVE: with 94 smooth history nodes it showed a clean
                // 4.00x per halving in every sector and exonerated the
                // retarded fields.  The real failure carries 31 nodes and
                // behaves completely differently.  Measured at r=1.056e-12 m,
                // dtFail=5.714e-25 s, sweeping dt from 16x down to 1/32x:
                //
                //   retarded  1.89e-9  9.31e-10  4.67e-10  3.77e-3  9.45e-4
                //             3.19e-13  2.58e-13  1.38e-13  2.30e-14  1.20e-14
                //   Coulomb   4.14e-11 1.04e-11 2.60e-12 6.52e-13 1.63e-13
                //             4.10e-14 1.04e-14 2.55e-15 6.13e-16 2.91e-16
                //
                // Coulomb-only, on the identical state and history, is clean
                // to four digits (3.99x throughout).  The retarded sector
                // spikes SEVEN orders between adjacent step sizes, passes
                // through exactly the 9.45e-4 the engine rejects on, then
                // drops TEN orders.  So it is not a lost convergence order at
                // all, it is a discontinuity in dt.  Radiation reaction is
                // irrelevant: the two retarded rows are bit-identical with it
                // on and off.
                //
                // One mechanism was proposed and REFUTED here, which is worth
                // recording so it is not proposed again: boundedDerivativeStep
                // returning its hard zero fallback would flip derivatives
                // discontinuously, but counting its firings around each trial
                // step gives 0 for every retarded step including the anomalous
                // one, and 108 for every Coulomb step -- constantly firing in
                // the sector that stays smooth, never firing in the sector
                // that jumps.  (The Coulomb count is allExternalForces' own
                // single-node history, where span=0 by construction.)
                //
                // Left narrowed but not pinned: the retarded-time solve, and
                // historicalState's piecewise interpolation across a 31-node
                // history, which is the difference from the 94-node synthetic
                // case that hid the effect.
                if(std::getenv("CREM_DEBUG_ORDER_LIVE")) {
                    static int liveSweeps=0;
                    if(liveSweeps++<1) {
                        struct LiveCase { const char* name; bool retarded;
                                          ChargeRadiationReactionModel model; };
                        const LiveCase liveCases[]={
                            {"retarded+prod",true,accuracy_.reactionModel},
                            {"retarded+noRR",true,
                             ChargeRadiationReactionModel::disabled},
                            {"Coulomb +prod",false,accuracy_.reactionModel},
                            {"Coulomb +noRR",false,
                             ChargeRadiationReactionModel::disabled}};
                        std::cerr<<"LIVE state r="<<separation(start)
                            <<" historyNodes="<<history.size()
                            <<" dtFail="<<dt<<'\n';
                        for(const LiveCase& liveCase:liveCases) {
                            double previous=0.0;
                            std::cerr<<"LIVE "<<liveCase.name;
                            for(int halving=0;halving<10;++halving) {
                                const double sweepStep=
                                    dt*std::pow(2.0,4.0-halving);
                                State coarse=start;
                                integrateElectrodynamicStep(coarse,sweepStep,
                                    history,false,liveCase.model,
                                    liveCase.retarded);
                                State fine=start;
                                StateHistory fineHistory=history;
                                integrateElectrodynamicStep(fine,
                                    0.5*sweepStep,fineHistory,false,
                                    liveCase.model,liveCase.retarded);
                                integrateElectrodynamicStep(fine,
                                    0.5*sweepStep,fineHistory,false,
                                    liveCase.model,liveCase.retarded);
                                const double error=
                                    normalizedStepError(coarse,fine);
                                std::cerr<<"  "<<error;
                                // Which retarded term moves with dt.  Evaluated
                                // at the coarse endpoint: the charge's
                                // Lienard-Wiechert Lorentz force, the retarded
                                // magnetic dipole field's own contribution, and
                                // covariantDipoleGradientForce, which runs its
                                // own six-point spatial stencil and so has the
                                // most machinery to be discontinuous in.
                                if(liveCase.retarded
                                   &&std::getenv("CREM_DEBUG_ORDER_TERMS")) {
                                    const ElectromagneticField lw=
                                        lienardWiechertField(
                                            coarse.firstPosition,coarse.time,
                                            history,coarse,false,secondCharge);
                                    const ElectromagneticField dip=
                                        retardedMagneticDipoleField(
                                            coarse.firstPosition,coarse.time,
                                            history,coarse,false);
                                    const Vec3 grad=
                                        covariantDipoleGradientForce(
                                            coarse,history,true);
                                    std::cerr<<"{LW="<<lw.electric.norm()
                                        <<" dipE="<<dip.electric.norm()
                                        <<" dipB="<<dip.magnetic.norm()
                                        <<" grad="<<grad.norm()<<"}";
                                }
                                if(halving>0&&error>0.0)
                                    std::cerr<<"("<<previous/error<<"x)";
                                previous=error;
                            }
                            std::cerr<<'\n';
                        }
                    }
                }
                return false;
            }
            State midpoint;
            StateHistory midpointHistory;
            if(!advanceAdaptive(start,history,0.5*dt,depth+1,
                                midpoint,midpointHistory)) {
                return false;
            }
            return advanceAdaptive(midpoint,midpointHistory,0.5*dt,depth+1,
                                   accepted,acceptedHistory);
        };

        // The coarse step is a pure error probe: it is discarded on every
        // path, so it never needs the far-zone flux integration.
        State coarse=start;
        integrateElectrodynamicStep(coarse,dt,history,false,
            accuracy_.reactionModel,accuracy_.useRetardedExternalForces);
        if(!isFinite(coarse)) {
            rejectionReason="coarse-nonfinite";
            return subdivide();
        }
        if(!(coarse.time>start.time)) {
            if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
                std::cerr<<"ENGINE_DEBUG reason=coarse-time dt="<<dt
                    <<" t="<<start.time<<" depth="<<depth<<'\n';
            return false;
        }

        // The two half-steps are the path that *becomes* the trajectory when
        // the step is accepted, so they carry the complete bookkeeping from
        // the start and are committed instead of being recomputed.  The flux
        // flag does not feed back into positions or velocities, so the local
        // error estimate below is unchanged by enabling it here.
        State fine=start;
        StateHistory fineHistory=history;
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,
            accuracy_.computeOutwardFlux,accuracy_.reactionModel,
            accuracy_.useRetardedExternalForces,
            accuracy_.farFieldSampling);
        if(!isFinite(fine)) {
            rejectionReason="fine1-nonfinite";
            return subdivide();
        }
        if(!(fine.time>start.time)) {
            if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
                std::cerr<<"ENGINE_DEBUG reason=fine1-time dt="<<dt
                    <<" t="<<start.time<<" depth="<<depth<<'\n';
            return false;
        }
        const double midpointTime=fine.time;
        appendStateHistory(fineHistory,fine);
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,
            accuracy_.computeOutwardFlux,accuracy_.reactionModel,
            accuracy_.useRetardedExternalForces,
            accuracy_.farFieldSampling);
        if(!isFinite(fine)) {
            rejectionReason="fine2-nonfinite";
            return subdivide();
        }
        if(!(fine.time>midpointTime)) {
            if(std::getenv("POSITRONIUM_DEBUG_DIPOLE"))
                std::cerr<<"ENGINE_DEBUG reason=fine2-time dt="<<dt
                    <<" t="<<start.time<<" depth="<<depth<<'\n';
            return false;
        }

        const double error=normalizedStepError(coarse,fine);
        if(!std::isfinite(error)) {
            rejectionReason="error-nonfinite";
            return subdivide();
        }
        if(error<=accuracy_.relativeTolerance) {
            // CREM_STEP_CENSUS: accepted-step statistics, so a run can be
            // asked whether its step actually RESOLVES the forces acting on
            // it.  Written for the --zpf question, where the field band's
            // upper edge is a frequency the trajectory error probe does not
            // necessarily see: a mode far off resonance perturbs the ORBIT
            // very little (so the error stays inside tolerance and the step
            // is never halved) while its own phase advances a great deal per
            // step.  The work integral of such a mode should average to
            // nothing over a cycle and instead gets sampled a few times per
            // cycle, which is the textbook setup for aliasing.  Counting the
            // steps is how that stops being a hypothesis.
            //
            // Deliberately global and unsynchronised: this is a diagnostic
            // for single-trajectory --diagnose runs, not for the threaded
            // statistical mode, and it is off unless the variable is set.
            if(gStepCensusEnabled) {
                ++gStepCensusCount;
                gStepCensusTotalTime+=dt;
                if(dt<gStepCensusSmallest) gStepCensusSmallest=dt;
                if(dt>gStepCensusLargest) gStepCensusLargest=dt;
            }
            appendStateHistory(fineHistory,fine);
            accepted=fine;
            acceptedHistory=std::move(fineHistory);
            return true;
        }
        // Reaching the recursion limit is a failed accuracy contract, not an
        // alternative acceptance rule.  Returning false activates the caller's
        // recovery ladder instead of silently committing an under-resolved
        // trajectory.
        rejectionMetric=error;
        return subdivide();
    }
    StateHistory history_;
    Accuracy accuracy_;
};
