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
            accuracy_.useRetardedExternalForces);
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
            accuracy_.useRetardedExternalForces);
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
