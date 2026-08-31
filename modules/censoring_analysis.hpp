#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <ranges>
#include <vector>

// Outcome-censoring diagnostics for the beam and interaction experiments.
// This is deliberately separate from Kaplan-Meier: KM estimates an event-time
// distribution from right-censored times, whereas the experiments below ask
// which physical channel a trajectory reaches.  An administrative stop hides
// that categorical outcome, and a numerical failure is reported separately.
namespace positronium::statistics {

enum class ObservationDisposition {
    Observed,
    AdministrativelyCensored,
    NumericalFailure
};

struct CensoringObservation {
    double energy=std::numeric_limits<double>::quiet_NaN();
    double impactParameter=std::numeric_limits<double>::quiet_NaN();
    ObservationDisposition disposition=ObservationDisposition::NumericalFailure;
    int physicalCategory=-1;
};

struct CovariateScaling {
    double mean=0.0;
    double standardDeviation=1.0;
    bool varying=false;

    double standardized(double value) const noexcept {
        return varying?(value-mean)/standardDeviation:0.0;
    }
};

struct LogisticModel {
    // logit(p) = beta[0] + beta[1] z_energy + beta[2] z_impact.
    std::array<double,3> coefficient{};
    std::array<double,3> standardError{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    CovariateScaling energyScaling;
    CovariateScaling impactScaling;
    int sampleCount=0;
    int positiveCount=0;
    bool fitted=false;
    bool converged=false;
    bool constant=false;
    double constantProbability=std::numeric_limits<double>::quiet_NaN();
    double brierScore=std::numeric_limits<double>::quiet_NaN();

    double probability(double energy,double impactParameter) const noexcept;
};

struct CensoringAnalysis {
    // completionModel's positive outcome is a physical classification.
    // failureModel's positive outcome is specifically NumericalFailure.
    LogisticModel completionModel;
    LogisticModel failureModel;
    int validCount=0;
    int invalidCovariateCount=0;
    int observedCount=0;
    int administrativelyCensoredCount=0;
    int numericalFailureCount=0;
    // Hajek-normalized inverse-probability-of-censoring weighted physical
    // class distribution.  It estimates the distribution that would have
    // been seen if completion were missing at random conditional on (E,b).
    std::vector<double> ipcwCategoryProbability;
    double effectiveObservedSampleSize=0.0;
    double minimumCompletionProbability=std::numeric_limits<double>::quiet_NaN();
    double maximumIpcwWeight=std::numeric_limits<double>::quiet_NaN();
    // Explicit positivity guard.  Predictions below this value are clipped,
    // so a handful of nearly-unobservable trajectories cannot dominate.
    double probabilityFloor=0.05;
};

enum class BinaryEndpoint {
    NotPhysicallyObserved,
    NumericalFailure
};

struct BinnedRate {
    double covariate=0.0;
    double covariateHalfRange=0.0;
    double rate=0.0;
    double rateErrorLow=0.0;
    double rateErrorHigh=0.0;
    int total=0;
    int positive=0;
};

namespace detail {

inline double logistic(double value) noexcept {
    if(value>=0.0) {
        const double exponential=std::exp(-value);
        return 1.0/(1.0+exponential);
    }
    const double exponential=std::exp(value);
    return exponential/(1.0+exponential);
}

inline CovariateScaling scalingOf(
        const std::vector<CensoringObservation>& observations,
        bool energy) {
    CovariateScaling result;
    if(observations.empty()) return result;
    for(const CensoringObservation& observation:observations) {
        result.mean+=energy?observation.energy:observation.impactParameter;
    }
    result.mean/=static_cast<double>(observations.size());
    double variance=0.0;
    for(const CensoringObservation& observation:observations) {
        const double value=energy?observation.energy:observation.impactParameter;
        const double difference=value-result.mean;
        variance+=difference*difference;
    }
    variance/=static_cast<double>(observations.size());
    result.standardDeviation=std::sqrt(std::max(0.0,variance));
    const double resolution=64.0*std::numeric_limits<double>::epsilon()
        *std::max(1.0,std::abs(result.mean));
    result.varying=result.standardDeviation>resolution;
    if(!result.varying) result.standardDeviation=1.0;
    return result;
}

inline bool solveThreeByThree(std::array<std::array<double,3>,3> matrix,
                              std::array<double,3> right,
                              std::array<double,3>& solution) {
    for(std::size_t column=0;column<3;++column) {
        std::size_t pivot=column;
        for(std::size_t row=column+1;row<3;++row) {
            if(std::abs(matrix[row][column])
               >std::abs(matrix[pivot][column])) pivot=row;
        }
        if(!(std::abs(matrix[pivot][column])>1.0e-18)
           ||!std::isfinite(matrix[pivot][column])) return false;
        if(pivot!=column) {
            std::swap(matrix[pivot],matrix[column]);
            std::swap(right[pivot],right[column]);
        }
        const double divisor=matrix[column][column];
        for(std::size_t entry=column;entry<3;++entry)
            matrix[column][entry]/=divisor;
        right[column]/=divisor;
        for(std::size_t row=0;row<3;++row) {
            if(row==column) continue;
            const double factor=matrix[row][column];
            for(std::size_t entry=column;entry<3;++entry)
                matrix[row][entry]-=factor*matrix[column][entry];
            right[row]-=factor*right[column];
        }
    }
    solution=right;
    return std::ranges::all_of(solution,
        [](double value){return std::isfinite(value);});
}

inline bool invertThreeByThree(
        const std::array<std::array<double,3>,3>& matrix,
        std::array<std::array<double,3>,3>& inverse) {
    for(std::size_t column=0;column<3;++column) {
        std::array<double,3> right{};
        right[column]=1.0;
        std::array<double,3> solution{};
        if(!solveThreeByThree(matrix,right,solution)) return false;
        for(std::size_t row=0;row<3;++row) inverse[row][column]=solution[row];
    }
    return true;
}

template<class Response>
inline LogisticModel fitLogistic(
        const std::vector<CensoringObservation>& observations,
        const CovariateScaling& energyScaling,
        const CovariateScaling& impactScaling,
        Response response) {
    LogisticModel model;
    model.energyScaling=energyScaling;
    model.impactScaling=impactScaling;
    model.sampleCount=static_cast<int>(observations.size());
    for(const CensoringObservation& observation:observations)
        if(response(observation)) ++model.positiveCount;
    if(observations.empty()) return model;
    model.fitted=true;
    const double empirical=static_cast<double>(model.positiveCount)
        /static_cast<double>(model.sampleCount);
    if(model.positiveCount==0||model.positiveCount==model.sampleCount
       ||(!energyScaling.varying&&!impactScaling.varying)) {
        model.constant=true;
        model.constantProbability=empirical;
        model.converged=true;
        model.brierScore=empirical*(1.0-empirical);
        return model;
    }

    // Unit L2 penalty on standardized slopes.  Small interaction samples can
    // be completely separated (for example, every failure below one sampled
    // b threshold), where unpenalized maximum likelihood has infinite slopes
    // and meaningless odds ratios.  lambda=1 is explicit, deterministic and
    // weak relative to a production-size sample while keeping previews finite.
    constexpr double ridge=1.0;
    const double smoothed=(model.positiveCount+0.5)/(model.sampleCount+1.0);
    model.coefficient[0]=std::log(smoothed/(1.0-smoothed));
    const auto objective=[&](const std::array<double,3>& beta) {
        double value=-0.5*ridge*(beta[1]*beta[1]+beta[2]*beta[2]);
        for(const CensoringObservation& observation:observations) {
            const std::array<double,3> x{
                1.0,energyScaling.standardized(observation.energy),
                impactScaling.standardized(observation.impactParameter)};
            const double eta=beta[0]+beta[1]*x[1]+beta[2]*x[2];
            const double logNormalizer=eta>0.0
                ?eta+std::log1p(std::exp(-eta)):std::log1p(std::exp(eta));
            value+=(response(observation)?eta:0.0)-logNormalizer;
        }
        return value;
    };
    for(int iteration=0;iteration<100;++iteration) {
        std::array<double,3> gradient{};
        std::array<std::array<double,3>,3> information{};
        for(const CensoringObservation& observation:observations) {
            const std::array<double,3> x{
                1.0,energyScaling.standardized(observation.energy),
                impactScaling.standardized(observation.impactParameter)};
            const double eta=model.coefficient[0]
                +model.coefficient[1]*x[1]+model.coefficient[2]*x[2];
            const double probability=logistic(eta);
            const double residual=(response(observation)?1.0:0.0)-probability;
            const double weight=probability*(1.0-probability);
            for(std::size_t row=0;row<3;++row) {
                gradient[row]+=x[row]*residual;
                for(std::size_t column=0;column<3;++column)
                    information[row][column]+=weight*x[row]*x[column];
            }
        }
        for(std::size_t index=1;index<3;++index) {
            gradient[index]-=ridge*model.coefficient[index];
            information[index][index]+=ridge;
        }
        std::array<double,3> update{};
        if(!solveThreeByThree(information,gradient,update)) break;
        const double oldObjective=objective(model.coefficient);
        double step=1.0;
        std::array<double,3> candidate=model.coefficient;
        bool accepted=false;
        for(int backtrack=0;backtrack<24;++backtrack) {
            for(std::size_t index=0;index<3;++index)
                candidate[index]=model.coefficient[index]+step*update[index];
            if(objective(candidate)>=oldObjective-1.0e-12) {
                accepted=true;
                break;
            }
            step*=0.5;
        }
        if(!accepted) break;
        model.coefficient=candidate;
        double largestUpdate=0.0;
        for(double value:update)
            largestUpdate=std::max(largestUpdate,std::abs(step*value));
        if(largestUpdate<1.0e-9) {
            model.converged=true;
            break;
        }
    }

    std::array<std::array<double,3>,3> information{};
    double squaredError=0.0;
    for(const CensoringObservation& observation:observations) {
        const std::array<double,3> x{
            1.0,energyScaling.standardized(observation.energy),
            impactScaling.standardized(observation.impactParameter)};
        const double probability=logistic(model.coefficient[0]
            +model.coefficient[1]*x[1]+model.coefficient[2]*x[2]);
        const double weight=probability*(1.0-probability);
        const double residual=(response(observation)?1.0:0.0)-probability;
        squaredError+=residual*residual;
        for(std::size_t row=0;row<3;++row)
            for(std::size_t column=0;column<3;++column)
                information[row][column]+=weight*x[row]*x[column];
    }
    for(std::size_t index=1;index<3;++index)
        information[index][index]+=ridge;
    std::array<std::array<double,3>,3> covariance{};
    if(invertThreeByThree(information,covariance)) {
        model.standardError[0]=std::sqrt(std::max(0.0,covariance[0][0]));
        if(energyScaling.varying)
            model.standardError[1]=std::sqrt(std::max(0.0,covariance[1][1]));
        if(impactScaling.varying)
            model.standardError[2]=std::sqrt(std::max(0.0,covariance[2][2]));
    }
    model.brierScore=squaredError/static_cast<double>(observations.size());
    return model;
}

} // namespace detail

inline double LogisticModel::probability(
        double energy,double impactParameter) const noexcept {
    if(!fitted||!std::isfinite(energy)||!std::isfinite(impactParameter))
        return std::numeric_limits<double>::quiet_NaN();
    if(constant) return constantProbability;
    return detail::logistic(coefficient[0]
        +coefficient[1]*energyScaling.standardized(energy)
        +coefficient[2]*impactScaling.standardized(impactParameter));
}

inline CensoringAnalysis analyzeCensoring(
        const std::vector<CensoringObservation>& sample,
        std::size_t physicalCategoryCount,
        double probabilityFloor=0.05) {
    CensoringAnalysis result;
    result.probabilityFloor=std::clamp(probabilityFloor,1.0e-6,1.0);
    std::vector<CensoringObservation> valid;
    valid.reserve(sample.size());
    for(const CensoringObservation& observation:sample) {
        if(!std::isfinite(observation.energy)
           ||!std::isfinite(observation.impactParameter)) {
            ++result.invalidCovariateCount;
            continue;
        }
        valid.push_back(observation);
        switch(observation.disposition) {
            case ObservationDisposition::Observed: ++result.observedCount; break;
            case ObservationDisposition::AdministrativelyCensored:
                ++result.administrativelyCensoredCount; break;
            case ObservationDisposition::NumericalFailure:
                ++result.numericalFailureCount; break;
        }
    }
    result.validCount=static_cast<int>(valid.size());
    const CovariateScaling energyScaling=detail::scalingOf(valid,true);
    const CovariateScaling impactScaling=detail::scalingOf(valid,false);
    result.completionModel=detail::fitLogistic(
        valid,energyScaling,impactScaling,[](const CensoringObservation& value) {
            return value.disposition==ObservationDisposition::Observed;
        });
    result.failureModel=detail::fitLogistic(
        valid,energyScaling,impactScaling,[](const CensoringObservation& value) {
            return value.disposition==ObservationDisposition::NumericalFailure;
        });

    result.ipcwCategoryProbability.assign(physicalCategoryCount,0.0);
    double weightSum=0.0;
    double squaredWeightSum=0.0;
    double minimumProbability=std::numeric_limits<double>::infinity();
    double maximumWeight=0.0;
    for(const CensoringObservation& observation:valid) {
        const double probability=result.completionModel.probability(
            observation.energy,observation.impactParameter);
        if(std::isfinite(probability))
            minimumProbability=std::min(minimumProbability,probability);
        if(observation.disposition!=ObservationDisposition::Observed) continue;
        const double weight=1.0/std::max(result.probabilityFloor,probability);
        weightSum+=weight;
        squaredWeightSum+=weight*weight;
        maximumWeight=std::max(maximumWeight,weight);
        if(observation.physicalCategory>=0
           &&static_cast<std::size_t>(observation.physicalCategory)
                <physicalCategoryCount) {
            result.ipcwCategoryProbability[
                static_cast<std::size_t>(observation.physicalCategory)]+=weight;
        }
    }
    if(weightSum>0.0) {
        for(double& probability:result.ipcwCategoryProbability)
            probability/=weightSum;
    }
    if(squaredWeightSum>0.0)
        result.effectiveObservedSampleSize=weightSum*weightSum/squaredWeightSum;
    if(std::isfinite(minimumProbability))
        result.minimumCompletionProbability=minimumProbability;
    if(weightSum>0.0) result.maximumIpcwWeight=maximumWeight;
    return result;
}

inline std::vector<BinnedRate> binnedRate(
        const std::vector<CensoringObservation>& sample,
        bool useEnergy,BinaryEndpoint endpoint,int requestedBins=8) {
    std::vector<CensoringObservation> valid;
    valid.reserve(sample.size());
    for(const CensoringObservation& observation:sample) {
        if(std::isfinite(observation.energy)
           &&std::isfinite(observation.impactParameter)) valid.push_back(observation);
    }
    if(valid.empty()) return {};
    const auto covariate=[&](const CensoringObservation& observation) {
        return useEnergy?observation.energy:observation.impactParameter;
    };
    std::sort(valid.begin(),valid.end(),[&](const CensoringObservation& first,
                                            const CensoringObservation& second) {
        return covariate(first)<covariate(second);
    });
    const bool constant=covariate(valid.front())==covariate(valid.back());
    const int targetBins=constant?1:std::clamp(requestedBins,1,
        static_cast<int>(valid.size()));
    const std::size_t targetSize=(valid.size()+targetBins-1)
        /static_cast<std::size_t>(targetBins);
    std::vector<BinnedRate> result;
    constexpr double normal95=1.95996398454005;
    for(std::size_t begin=0;begin<valid.size();) {
        std::size_t end=std::min(valid.size(),begin+targetSize);
        while(end<valid.size()
              &&covariate(valid[end])==covariate(valid[end-1])) ++end;
        BinnedRate point;
        point.total=static_cast<int>(end-begin);
        double covariateSum=0.0;
        double minimum=std::numeric_limits<double>::infinity();
        double maximum=-std::numeric_limits<double>::infinity();
        for(std::size_t index=begin;index<end;++index) {
            const CensoringObservation& observation=valid[index];
            const double value=covariate(observation);
            covariateSum+=value;
            minimum=std::min(minimum,value);
            maximum=std::max(maximum,value);
            const bool positive=endpoint==BinaryEndpoint::NumericalFailure
                ?observation.disposition==ObservationDisposition::NumericalFailure
                :observation.disposition!=ObservationDisposition::Observed;
            if(positive) ++point.positive;
        }
        point.covariate=covariateSum/static_cast<double>(point.total);
        point.covariateHalfRange=std::max(point.covariate-minimum,
                                          maximum-point.covariate);
        point.rate=static_cast<double>(point.positive)
            /static_cast<double>(point.total);
        const double denominator=1.0+normal95*normal95/point.total;
        const double centre=(point.rate+normal95*normal95/(2.0*point.total))
            /denominator;
        const double halfWidth=normal95/denominator*std::sqrt(
            point.rate*(1.0-point.rate)/point.total
            +normal95*normal95/(4.0*point.total*point.total));
        const double lower=std::max(0.0,centre-halfWidth);
        const double upper=std::min(1.0,centre+halfWidth);
        point.rateErrorLow=point.rate-lower;
        point.rateErrorHigh=upper-point.rate;
        result.push_back(point);
        begin=end;
    }
    return result;
}

} // namespace positronium::statistics
