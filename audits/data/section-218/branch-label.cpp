// Audit 218: what distinguishes the integrator's two solution branches
// that 217f observed.  At beta = 0, quarter period, 0.25 a_pair,
//   (200,1e-8) == (400,1e-8)  to the last digit,
//   (800,1e-8) == (200,1e-9)  to the last digit,
// and the two pairs differ from each other by 30x in the covariance
// residual.  The question is whether the TRAJECTORY differs or only
// the radiated bookkeeping: history density changes the retarded
// field, and the flux is read off that, so a branch could live
// entirely in the sampling with the mechanics bit-identical.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <array>
#include <vector>

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

    struct Result { State end; std::size_t history; double flux;
                    Vec3 fluxMomentum; };
    const auto go=[&](int steps,double tol,Result& out)->bool{
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(rest,acc);
        State s=rest;
        for(int i=0;i<steps;++i)
            if(!engine.advance(s,span/steps)) return false;
        FarFieldSampling sampling; sampling.directionCount=50;
        const FieldFluxRates f=
            electromagneticFieldFluxRates(s,engine.history(),sampling);
        out=Result{s,engine.history().size(),f.energy,f.momentum};
        return true;
    };
    Result ref;
    if(!go(200,1.0e-8,ref)){ std::printf("reference failed\n"); return 1; }
    const auto rel=[](double a,double b){
        return b!=0.0?std::abs(a-b)/std::abs(b):0.0; };

    std::printf("# beta = 0.  Reference is (200 steps, tol 1e-8).\n");
    std::printf("# 'position' and 'velocity' are relative gaps in the\n");
    std::printf("# FINAL mechanical state; 'E_rad' and 'flux' in the\n");
    std::printf("# radiated bookkeeping read off the history.\n");
    std::printf("\n%7s %9s %8s %12s %12s %12s %12s\n",
                "steps","tolerance","history","position","velocity",
                "E_rad","flux");
    for(const std::array<double,2>& cfg:std::vector<std::array<double,2>>{
            {200,1e-8},{400,1e-8},{800,1e-8},{1600,1e-8},
            {200,1e-7},{200,1e-9},{200,1e-10},{400,1e-9},{800,1e-9}}){
        Result r;
        if(!go(static_cast<int>(cfg[0]),cfg[1],r)){
            std::printf("%7.0f %9.0e failed\n",cfg[0],cfg[1]); continue; }
        const Vec3 dp=r.end.firstPosition-ref.end.firstPosition;
        const Vec3 dv=r.end.firstVelocity-ref.end.firstVelocity;
        std::printf("%7.0f %9.0e %8zu %12.4e %12.4e %12.4e %12.4e\n",
            cfg[0],cfg[1],r.history,
            dp.norm()/ref.end.firstPosition.norm(),
            dv.norm()/ref.end.firstVelocity.norm(),
            rel(r.end.radiatedEnergy,ref.end.radiatedEnergy),
            rel(r.flux,ref.flux));
    }
    std::printf("\n# If position and velocity are at roundoff while\n");
    std::printf("# E_rad and flux are not, the branch is in the\n");
    std::printf("# retarded sampling, not in the mechanics.\n");
}
