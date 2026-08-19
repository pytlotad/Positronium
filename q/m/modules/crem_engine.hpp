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
        ChargeRadiationReactionModel reactionModel =
            ChargeRadiationReactionModel::individualLandauLifshitz;
        // Far-zone Poynting quadrature: 50 directions, each solving a
        // retarded time for both charges, on every fine half-step.  It feeds
        // only the flux BOOKKEEPING (radiatedEnergy/Momentum/AngularMomentum,
        // boundField*, reaction mismatches).  The trajectory itself never
        // reads any of them: forces, chargeReaction, magneticDipolePower and
        // hence dipoleConstraintEnergy are all computed outside that block,
        // so switching it off leaves positions, velocities and dipoles
        // bit-identical.  Callers that do not report a radiated energy can
        // therefore turn it off outright.
        bool computeOutwardFlux = true;
    };
    explicit ClassicalTrajectoryEngine(const State& initial)
        :history_(causalInitialHistory(initial)) {}
    ClassicalTrajectoryEngine(const State& initial,Accuracy accuracy)
        :history_(causalInitialHistory(initial)),accuracy_(accuracy) {}
    ClassicalTrajectoryEngine(StateHistory history,Accuracy accuracy)
        :history_(std::move(history)),accuracy_(accuracy) {}
    bool advance(State& state,double dt) {
        if(!(dt>0.0)||!std::isfinite(dt)) return false;
        return advanceAdaptive(state,dt,0);
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

    bool advanceAdaptive(State& state,double dt,int depth) {
        const State start=state;

        // The coarse step is a pure error probe: it is discarded on every
        // path, so it never needs the far-zone flux integration.
        State coarse=start;
        integrateElectrodynamicStep(coarse,dt,history_,false,
            accuracy_.reactionModel);

        // The two half-steps are the path that *becomes* the trajectory when
        // the step is accepted, so they carry the complete bookkeeping from
        // the start and are committed instead of being recomputed.  The flux
        // flag does not feed back into positions or velocities, so the local
        // error estimate below is unchanged by enabling it here.
        State fine=start;
        StateHistory fineHistory=history_;
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,
            accuracy_.computeOutwardFlux,accuracy_.reactionModel);
        if(!isFinite(fine)) return false;
        appendStateHistory(fineHistory,fine);
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,
            accuracy_.computeOutwardFlux,accuracy_.reactionModel);
        if(!isFinite(coarse)||!isFinite(fine)) return false;

        const double error=normalizedStepError(coarse,fine);
        if(!std::isfinite(error)) return false;
        if(error<=accuracy_.relativeTolerance||depth>=accuracy_.maximumDepth) {
            appendStateHistory(fineHistory,fine);
            history_=std::move(fineHistory);
            state=fine;
            return true;
        }
        // Nothing above touched state or history_, so the subdivision below
        // restarts from exactly the caller's state.
        return advanceAdaptive(state,0.5*dt,depth+1)
            &&advanceAdaptive(state,0.5*dt,depth+1);
    }
    StateHistory history_;
    Accuracy accuracy_;
};
