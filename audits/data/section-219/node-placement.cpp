// Audit 219, item 1: at beta = 0.35, N = 1600 the covariance residual
// jumps 40x although the two frames report the same history node COUNT
// (120/120).  Matching the count is not matching the PLACEMENT.  This
// compares the node times themselves: the moving frame's times are
// divided by gamma to bring them onto the rest frame's clock, and the
// largest relative discrepancy in node spacing is reported.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
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
    const double beta=0.35, g=1.0/std::sqrt(1.0-beta*beta);
    State moving=rest;
    moving.firstVelocity=rest.firstVelocity*(1.0/g)+Vec3{0,0,beta*c};
    moving.secondVelocity=rest.secondVelocity*(1.0/g)+Vec3{0,0,beta*c};
    const double retention=std::max(1.0e-20,4.0*r0/c);
    const double gate=retention/128.0;

    const auto times=[&](const State& start,double total,int steps,
                         std::vector<double>& out)->bool{
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-8; acc.maximumDepth=24;
        ClassicalTrajectoryEngine engine(start,acc);
        State s=start;
        for(int i=0;i<steps;++i)
            if(!engine.advance(s,total/steps)) return false;
        out.clear();
        for(const State& h:engine.history()) out.push_back(h.time);
        return true;
    };
    std::printf("%7s %8s %7s %9s %9s %14s %14s\n",
                "N_rest","N_moving","dt/gate","n_rest","n_mov",
                "max spacing gap","mean gap");
    for(int N:{200,400,800,1600,3200}){
        const int Nm=static_cast<int>(std::lround(N*g));
        std::vector<double> tr,tm;
        if(!times(rest,span,N,tr)||!times(moving,g*span,Nm,tm)){
            std::printf("%7d  failed\n",N); continue; }
        // Bring the moving clock onto the rest clock and compare the
        // spacings from the end backwards, which is what the retarded
        // evaluation at the final instant actually samples.
        double worst=0.0,sum=0.0; std::size_t count=0;
        const std::size_t n=std::min(tr.size(),tm.size());
        for(std::size_t i=1;i<n;++i){
            const double dr=tr[tr.size()-i]-tr[tr.size()-i-1];
            const double dm=(tm[tm.size()-i]-tm[tm.size()-i-1])/g;
            if(!(dr>0.0)) continue;
            const double rel=std::abs(dm-dr)/dr;
            worst=std::max(worst,rel); sum+=rel; ++count;
        }
        std::printf("%7d %8d %7.3f %9zu %9zu %14.6e %14.6e\n",
            N,Nm,(span/N)/gate,tr.size(),tm.size(),worst,
            count?sum/count:0.0);
    }
    std::printf("\n# Equal node COUNTS with unequal spacings mean the\n");
    std::printf("# rungs match but the grids do not.\n");
}
