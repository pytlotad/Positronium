#pragma once

// Product-limit (Kaplan-Meier) survival estimator with Greenwood standard
// errors, for right-censored collapse-time samples: a trajectory stopped by
// the wall-clock budget has NOT told us nothing, it has told us the
// collapse time exceeds the simulated time it reached, and averaging only
// the completed runs throws that away -- badly, because the budget
// preferentially stops the widest orbits, which are exactly the slowest to
// collapse.
//
// MEASURED, because "informative censoring breaks the independence this
// estimator assumes" is usually where such a note stops and it matters
// whether the estimate actually moves.  An ortho batch at --level 1 was run
// twice from the same seed, once at a 110 s per-event budget leaving 2 of 12
// censored and once at 400 s leaving none:
//
//     110 s   83.3% complete   KM median 147.558 ps   completed-mean 147.664
//     400 s    100% complete   KM median 147.558 ps   completed-mean 147.657
//
// The product-limit median did not move at all, and the completed-run mean
// moved by 5e-5 relative.  So on this sample the estimator absorbed the
// censoring it is warned about, which is the outcome to expect while the
// censored fraction stays small -- not a licence to ignore the warning at
// heavier censoring, where the two ends of the reported range genuinely
// diverge.  One trajectory censored at 110 s turned into a numerical failure
// rather than a completion when given the longer budget.
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace, itself already inside an
// #ifndef POSITRONIUM_VALIDATION_EXECUTABLE region positronium.cpp opens
// well before this #include -- so this header needs no guard of its own.

// One trajectory's contribution to a right-censored survival sample.
// observed=true means the collapse was actually seen at `time`; observed=false
// means the run stopped at `time` with the pair still bound, so all that is
// known is T > time.
struct SurvivalObservation { double time=0.0; bool observed=false; };

struct KaplanMeierPoint {
    double time=0.0;
    double survival=1.0;
    double standardError=0.0; // Greenwood
    int atRisk=0;
    int events=0;
};

struct KaplanMeierEstimate {
    std::vector<KaplanMeierPoint> curve;
    double medianSurvival=std::numeric_limits<double>::quiet_NaN();
    bool medianReached=false;
    // Restricted mean survival time: the area under the curve out to
    // `horizon`.  With censored observations beyond the last collapse the
    // unrestricted mean is not identifiable, so RMST is the honest summary.
    double restrictedMean=std::numeric_limits<double>::quiet_NaN();
    double restrictedMeanError=std::numeric_limits<double>::quiet_NaN();
    double horizon=std::numeric_limits<double>::quiet_NaN();
    double survivalAtHorizon=std::numeric_limits<double>::quiet_NaN();
    double largestCensoredTime=std::numeric_limits<double>::quiet_NaN();
    int eventCount=0;
    int censoredCount=0;
};

// Product-limit estimator with Greenwood standard errors.  This is what the
// bound-decay experiments need: a trajectory stopped by the wall-clock budget
// has NOT told us nothing, it has told us the collapse time exceeds the
// simulated time it reached, and averaging only the completed runs throws that
// away -- badly, because the budget preferentially stops the widest orbits,
// which are exactly the slowest to collapse.
KaplanMeierEstimate kaplanMeier(std::vector<SurvivalObservation> sample) {
    KaplanMeierEstimate result;
    sample.erase(std::remove_if(sample.begin(),sample.end(),
        [](const SurvivalObservation& o) {
            return !std::isfinite(o.time)||o.time<0.0;
        }),sample.end());
    if(sample.empty()) return result;
    // Ties: events are ordered before censorings at the same time, the
    // standard convention -- a run censored at exactly t was still at risk
    // when the collapse at t happened.
    std::sort(sample.begin(),sample.end(),
        [](const SurvivalObservation& a,const SurvivalObservation& b) {
            if(a.time!=b.time) return a.time<b.time;
            return a.observed&&!b.observed;
        });
    for(const SurvivalObservation& o : sample) {
        if(o.observed) ++result.eventCount;
        else {
            ++result.censoredCount;
            result.largestCensoredTime=std::isfinite(result.largestCensoredTime)
                ?std::max(result.largestCensoredTime,o.time):o.time;
        }
    }
    double survival=1.0;
    double greenwoodSum=0.0;
    std::size_t index=0;
    const int total=static_cast<int>(sample.size());
    result.curve.push_back({0.0,1.0,0.0,total,0});
    while(index<sample.size()) {
        const double time=sample[index].time;
        const int atRisk=total-static_cast<int>(index);
        int events=0;
        std::size_t next=index;
        while(next<sample.size()&&sample[next].time==time) {
            if(sample[next].observed) ++events;
            ++next;
        }
        if(events>0) {
            survival*=1.0-static_cast<double>(events)/static_cast<double>(atRisk);
            if(atRisk>events) {
                greenwoodSum+=static_cast<double>(events)
                    /(static_cast<double>(atRisk)
                      *static_cast<double>(atRisk-events));
            } else {
                greenwoodSum=std::numeric_limits<double>::infinity();
            }
            // If the last risk-set member fails, Greenwood's sum is infinite
            // while S(t) is exactly zero.  The literal product 0*sqrt(inf) is
            // NaN; its limiting pointwise standard error is zero.
            const double standardError=survival==0.0
                ?0.0:survival*std::sqrt(greenwoodSum);
            result.curve.push_back(
                {time,survival,standardError,atRisk,events});
        }
        index=next;
    }
    // An all-censored sample is still a valid (degenerate) KM estimate: no
    // observed event lowers S(t), so S(t)=1 through the largest censoring time.
    // Continue through the common integration below to report a finite horizon,
    // RMST=horizon and Greenwood error 0.  The median correctly remains
    // unreached; nothing here extrapolates a lifetime beyond the observed data.
    // Horizon: the last time the estimator actually observed something, so
    // RMST is not extrapolated past the data.
    result.horizon=std::max(result.curve.back().time,
        std::isfinite(result.largestCensoredTime)?result.largestCensoredTime
                                                 :result.curve.back().time);
    result.survivalAtHorizon=result.curve.back().survival;
    for(const KaplanMeierPoint& point : result.curve) {
        if(!result.medianReached&&point.survival<=0.5) {
            result.medianSurvival=point.time;
            result.medianReached=true;
        }
    }
    // RMST = integral of the step function, plus its standard error from the
    // usual sum of squared tail areas weighted by the Greenwood increments.
    double area=0.0;
    for(std::size_t i=0;i+1<result.curve.size();++i) {
        area+=result.curve[i].survival
             *(result.curve[i+1].time-result.curve[i].time);
    }
    area+=result.curve.back().survival
         *(result.horizon-result.curve.back().time);
    result.restrictedMean=area;
    double varianceSum=0.0;
    for(std::size_t i=1;i<result.curve.size();++i) {
        if(result.curve[i].events==0) continue;
        double tailArea=0.0;
        for(std::size_t j=i;j+1<result.curve.size();++j) {
            tailArea+=result.curve[j].survival
                     *(result.curve[j+1].time-result.curve[j].time);
        }
        tailArea+=result.curve.back().survival
                 *(result.horizon-result.curve.back().time);
        const int atRisk=result.curve[i].atRisk;
        const int events=result.curve[i].events;
        if(atRisk>events) {
            varianceSum+=tailArea*tailArea*static_cast<double>(events)
                /(static_cast<double>(atRisk)
                  *static_cast<double>(atRisk-events));
        }
    }
    result.restrictedMeanError=std::sqrt(std::max(0.0,varianceSum));
    return result;
}
