// Audit 217: is 216f's 3-5e-05 covariance floor real non-covariance,
// or the noise of comparing two separately integrated trajectories?
//
// At beta = 0 covariance is satisfied identically, so ANY residual left
// when the same physical problem is integrated twice with different
// step patterns is pure comparison noise.  216f's two runs differed in
// exactly that way -- same step COUNT, step lengths differing by gamma,
// and a drifting trajectory the adaptive controller subdivides
// differently -- so the null has to vary the step pattern, not the
// frame.
//
// The observable and every other setting are 216f's: the far-zone flux
// four-vector rate from electromagneticFieldFluxRates on an orbit at
// 0.25 a_pair evolved a quarter period, tolerance 1e-8, 50 directions,
// and the same relativeGap.  Each row is compared against the 400-step
// run, which is the finest here.
//
// Reading: if these gaps reach 3-5e-05, 216f's floor is this noise and
// not the scheme.  If they sit orders below it, the floor is real.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <array>
#include <vector>

namespace {
struct Four { double time; Vec3 space; };
double relativeGap(const Four& a,const Four& b){
    const double num=std::sqrt(std::pow(a.time-b.time,2)
                              +(a.space-b.space).squaredNorm());
    const double den=std::sqrt(b.time*b.time+b.space.squaredNorm());
    return den>0.0?num/den:std::numeric_limits<double>::quiet_NaN();
}
}  // namespace

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double r0=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*r0),std::sqrt(k*r0)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const State rest=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    const double span=0.25*period;

    const auto rate=[&](int steps,double tol,Four& out,
                        Four& acc_out)->bool{
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(rest,acc);
        State s=rest;
        for(int i=0;i<steps;++i)
            if(!engine.advance(s,span/steps)) return false;
        FarFieldSampling sampling; sampling.directionCount=50;
        const FieldFluxRates f=
            electromagneticFieldFluxRates(s,engine.history(),sampling);
        out=Four{f.energy/c,f.momentum};
        // The ACCUMULATED four-momentum over the same interval: an
        // integral rather than an instantaneous rate, and therefore
        // the candidate for an observable quiet enough to carry a
        // covariance test.
        acc_out=Four{s.radiatedEnergy/c,s.radiatedMomentum};
        return true;
    };

    Four reference,referenceAcc;
    if(!rate(400,1.0e-8,reference,referenceAcc)){ std::printf("reference failed\n");
        return 1; }
    std::printf("# beta = 0 throughout: covariance is exact, so every\n");
    std::printf("# number below is comparison noise.  Reference is the\n");
    std::printf("# 400-step, tol 1e-8 run.  216f's floor was 3.556e-05.\n");
    std::printf("\n%8s %10s %18s %14s %14s\n",
                "steps","tolerance","|dP/dt| [W/c]","gap: rate",
                "gap: integral");
    for(const std::array<double,2>& cfg:std::vector<std::array<double,2>>{
            {400,1e-8},{200,1e-8},{214,1e-8},{201,1e-8},{204,1e-8},
            {231,1e-8},{187,1e-8},{150,1e-8},{800,1e-8},
            {200,1e-7},{200,1e-9},{200,1e-10},{400,1e-9}}){
        const int steps=static_cast<int>(cfg[0]);
        Four v,va;
        if(!rate(steps,cfg[1],v,va)){
            std::printf("%8d %10.0e  failed\n",steps,cfg[1]); continue; }
        std::printf("%8d %10.0e %18.10e %14.6e %14.6e\n",
            steps,cfg[1],
            std::sqrt(v.time*v.time+v.space.squaredNorm()),
            relativeGap(v,reference),relativeGap(va,referenceAcc));
    }
    std::printf("\n# 216f's two runs both used 200 steps with lengths\n");
    std::printf("# differing by gamma = 1.0675, i.e. the 214-step row\n");
    std::printf("# is the closest single analogue of that difference.\n");
}
