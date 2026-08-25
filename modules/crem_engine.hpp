#pragma once

// CREM numerics: reconstruction of a causal retarded history for a freshly
// prepared state, and the adaptive integrator that advances a trajectory with
// it.  This is the step-size machinery -- error probe, subdivision, history
// retention -- as opposed to the force laws in electrodynamics.hpp.
//
// Textual module, included from inside the anonymous namespace of
// positronium.cpp.  Unlike the other two CREM headers this one sits OUTSIDE
// the production #ifndef: positronium_validation drives the same engine.
// Contains no ROOT.

StateHistory causalInitialHistory(const State& initial,double spanFactor=8.0,
                                  int intervalCount=64) {
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
    for(int iteration=0;iteration<2;++iteration) {
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
        // 2 = the bare symmetric step; 4 = Yoshida composition of three of
        // them.  See composedStep() below for why raising this is the useful
        // knob and symplecticity is not.
        int compositionOrder = 2;
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
        const double lengthScale=std::max(separation(fine),nuclearCutoff);
        const double speedScale=std::max(
            (fine.firstVelocity-fine.secondVelocity).norm(),1.0e-6*c);
        return std::max({
            (coarse.firstPosition-fine.firstPosition).norm()/lengthScale,
            (coarse.secondPosition-fine.secondPosition).norm()/lengthScale,
            (coarse.firstVelocity-fine.firstVelocity).norm()/speedScale,
            (coarse.secondVelocity-fine.secondVelocity).norm()/speedScale});
    }

    // Yoshida's fourth-order composition of the symmetric second-order step,
    //   S4(dt) = S(w1 dt) . S(w0 dt) . S(w1 dt),
    //   w1 = 1/(2 - 2^(1/3)) = 1.3512,   w0 = -2^(1/3) w1 = -1.7024,
    // whose weights sum to one.  It is valid because the base step is exactly
    // time reversible: the suite's "conservative reverse" residual is 0.
    //
    // Why this and not a symplectic rewrite: the base step is ALREADY
    // kick-drift-kick, and the conservative energy error was measured to
    // oscillate with a constant envelope rather than drift -- 4.84, 4.66,
    // 4.62, 4.83e-05 over the four quarters of a run.  There is no secular
    // drift to remove.  What there is, is an amplitude eighteen times the
    // per-orbit physical signal of 2.6e-06, which is why the collapse
    // estimator has to subtract a reference run at all.  Composition attacks
    // exactly that: the bounded error falls as dt^4 instead of dt^2.
    //
    // The middle substep runs BACKWARDS in time, so the intermediate states
    // are deliberately not appended to the history -- doing so would make its
    // time coordinate non-monotonic.  All three substeps therefore read the
    // history as of the step start, which is the same approximation the bare
    // step already makes internally.
    void composedStep(State& state,double dt,const StateHistory& history,
                      bool computeFlux) {
        if(accuracy_.compositionOrder<4) {
            integrateElectrodynamicStep(state,dt,history,computeFlux,
                accuracy_.reactionModel,accuracy_.useRetardedExternalForces);
            return;
        }
        const double cubeRootTwo=std::cbrt(2.0);
        const double outer=1.0/(2.0-cubeRootTwo);
        const double inner=-cubeRootTwo*outer;
        const double weights[3]={outer,inner,outer};
        for(const double weight:weights) {
            integrateElectrodynamicStep(state,weight*dt,history,computeFlux,
                accuracy_.reactionModel,accuracy_.useRetardedExternalForces);
            if(!isFinite(state)) return;
        }
    }

    bool advanceAdaptive(const State& start,const StateHistory& history,
                         double dt,int depth,State& accepted,
                         StateHistory& acceptedHistory) {
        const double requestedEndTime=start.time+dt;
        if(!(requestedEndTime>start.time)||!std::isfinite(requestedEndTime)) {
            return false;
        }
        const auto subdivide=[&]() {
            if(depth>=accuracy_.maximumDepth) return false;
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
        composedStep(coarse,dt,history,false);
        if(!isFinite(coarse)) return subdivide();
        if(!(coarse.time>start.time)) return false;

        // The two half-steps are the path that *becomes* the trajectory when
        // the step is accepted, so they carry the complete bookkeeping from
        // the start and are committed instead of being recomputed.  The flux
        // flag does not feed back into positions or velocities, so the local
        // error estimate below is unchanged by enabling it here.
        State fine=start;
        StateHistory fineHistory=history;
        composedStep(fine,0.5*dt,fineHistory,accuracy_.computeOutwardFlux);
        if(!isFinite(fine)) return subdivide();
        if(!(fine.time>start.time)) return false;
        const double midpointTime=fine.time;
        appendStateHistory(fineHistory,fine);
        composedStep(fine,0.5*dt,fineHistory,accuracy_.computeOutwardFlux);
        if(!isFinite(fine)) return subdivide();
        if(!(fine.time>midpointTime)) return false;

        const double error=normalizedStepError(coarse,fine);
        if(!std::isfinite(error)) return subdivide();
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
        return subdivide();
    }
    StateHistory history_;
    Accuracy accuracy_;
};
