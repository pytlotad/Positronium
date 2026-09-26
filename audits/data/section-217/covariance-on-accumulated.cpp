// Audit 217: the covariance test rebuilt on the ACCUMULATED radiated
// four-momentum instead of the instantaneous flux rate.
//
// 216f measured the rate and found a residual of 3-5e-05 that no knob
// moved.  The beta = 0 null shows why: at zero boost, where covariance
// holds identically, the rate itself is only reproducible to 1.7e-05 to
// 1.6e-04 under a change of step pattern, so 216f's floor was the
// noise of its own observable and its verdict has no resolution.  The
// same null puts the ACCUMULATED four-momentum at 1.5e-07 to 3.4e-06 --
// about fifty times quieter, because it is an integral rather than a
// rate.  This rebuilds the test on it.
//
// P^mu accumulated over the interval is a genuine four-vector, so the
// prediction is P'^mu = Lambda^mu_nu P^nu with no 1/gamma: that factor
// in 216f belonged to the rate, not to the total.  The rest run covers
// span, the boosted run gamma*span, which is the same set of events.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <array>
#include <vector>
#include <cstdlib>

namespace {
struct Four { double time; Vec3 space; };
Four boostAlongZ(const Four& v,double beta){
    const double g=1.0/std::sqrt(1.0-beta*beta);
    return Four{g*(v.time-beta*v.space.z),
                Vec3{v.space.x,v.space.y,g*(v.space.z-beta*v.time)}};
}
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

    const auto run=[&](const State& start,double total,int steps,
                       double tol,Four& out)->bool{
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start;
        for(int i=0;i<steps;++i)
            if(!engine.advance(s,total/steps)) return false;
        out=Four{s.radiatedEnergy/c,s.radiatedMomentum};
        return true;
    };

    std::printf("# accumulated radiated four-momentum over a quarter\n");
    std::printf("# period at 0.25 a_pair; boost along the orbit normal.\n");
    std::printf("# beta = 0 step-pattern noise of this observable:\n");
    std::printf("#   1.5e-07 to 3.4e-06 (audit 217's null).\n");
    std::printf("\n%7s %8s %10s %16s %14s\n",
                "beta","steps","tolerance","|P| rest [kg m/s]","residual");
    for(const std::array<double,3>& cfg:std::vector<std::array<double,3>>{
            {0.05,200,1e-8},{0.10,200,1e-8},{0.20,200,1e-8},
            {0.35,200,1e-8},{0.50,200,1e-8},{0.70,200,1e-8},
            {0.35,400,1e-8},{0.35,800,1e-8},
            {0.35,200,1e-7},{0.35,200,1e-9}}){
        const double beta=cfg[0];
        const int steps=static_cast<int>(cfg[1]);
        const double tol=cfg[2];
        const double g=1.0/std::sqrt(1.0-beta*beta);
        State moving=rest;
        moving.firstVelocity=rest.firstVelocity*(1.0/g)+Vec3{0,0,beta*c};
        moving.secondVelocity=rest.secondVelocity*(1.0/g)+Vec3{0,0,beta*c};
        Four pRest,pMoving;
        if(!run(rest,span,steps,tol,pRest)
           ||!run(moving,g*span,steps,tol,pMoving)){
            std::printf("%7.2f %8d %10.0e  failed\n",beta,steps,tol);
            continue; }
        const Four predicted=boostAlongZ(pRest,-beta);
        std::printf("%7.2f %8d %10.0e %16.8e %14.6e\n",beta,steps,tol,
            std::sqrt(pRest.time*pRest.time+pRest.space.squaredNorm()),
            relativeGap(pMoving,predicted));
    }
    std::printf("\n# A residual at or below 3.4e-06 is indistinguishable\n");
    std::printf("# from the null and bounds non-covariance rather than\n");
    std::printf("# measuring it.  Above it, the excess is real.\n");
}
