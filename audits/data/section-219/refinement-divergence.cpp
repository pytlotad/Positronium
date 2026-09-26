// Audit 218: the covariance test with the two frames placed on the
// SAME rung of the history ladder.
//
// The retarded grid's node spacing is pinned to retentionTime/128, not
// to the integration step (electrodynamics.hpp:5743).  The boosted run
// covers gamma*span, so giving it the same step COUNT gives it a step
// gamma times longer and can land it on a different rung -- which is
// audit 217f's straddle.  Matching the coordinate step instead, by
// scaling the count with gamma, puts both frames on the same rung.
// Both choices are reported side by side.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <array>
#include <vector>

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
    const auto go=[&](const State& start,double total,int steps,
                      Four& out,std::size_t& nodes)->bool{
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-8; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start;
        for(int i=0;i<steps;++i)
            if(!engine.advance(s,total/steps)) return false;
        out=Four{s.radiatedEnergy/c,s.radiatedMomentum};
        nodes=engine.history().size();
        return true;
    };
    std::printf("%6s %7s %11s %7s %7s %13s %13s\n",
                "beta","N_rest","N_moving","h_rest","h_mov",
                "residual","matched?");
    for(double beta:{0.35}){
        const double g=1.0/std::sqrt(1.0-beta*beta);
        State moving=rest;
        moving.firstVelocity=rest.firstVelocity*(1.0/g)+Vec3{0,0,beta*c};
        moving.secondVelocity=rest.secondVelocity*(1.0/g)+Vec3{0,0,beta*c};
        for(int N:{800,1200,1600,2000,2400,3200}){
            const int matched=1;
            // matched: scale the count with gamma so the COORDINATE
            // step, and therefore dt/gate, agrees between the frames.
            const int Nm=static_cast<int>(std::lround(N*g));
            Four pr,pm; std::size_t hr=0,hm=0;
            if(!go(rest,span,N,pr,hr)
               ||!go(moving,g*span,Nm,pm,hm)){
                std::printf("%6.2f  failed\n",beta); continue; }
            std::printf("%6.2f %7d %11d %7zu %7zu %13.6e %13s\n",
                beta,N,Nm,hr,hm,
                relativeGap(pm,boostAlongZ(pr,-beta)),
                hr==hm?"same rung":"STRADDLE");
        }
    }
    std::printf("\n# If matching the rung collapses the residual, the\n");
    std::printf("# high-beta excess of 217e was the straddle and not\n");
    std::printf("# the boost.\n");
}
