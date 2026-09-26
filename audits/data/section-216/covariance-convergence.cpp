// Audit 216 / test A: a covariance residual that must fall with BOTH
// the integration step and the angular resolution.
//
// README states that the quantitatively trustworthy results are the
// single-operator unit tests at the stated resolution, and that long-
// time trajectories, radiation and reaction need a convergence and
// sensitivity analysis of their own.  The covariance figures it quotes
// -- 9.89e-16 for one Lienard-Wiechert field, 1.10e-06 for the boosted
// force, 6.99e-04 (later 4.06e-06) for the radiated four-vector -- are
// single numbers at one configuration.  A number is not a convergence
// claim.  This measures the same kind of residual as a SEQUENCE.
//
// The construction.  The observable is the far-zone flux four-vector
// rate dP^mu/dt from electromagneticFieldFluxRates, evaluated on an
// evolved orbit.  It is computed twice: in the pair's rest frame, and
// in a frame boosted along the ORBIT NORMAL.  That axis is chosen for
// an exact reason -- the orbit lies in z = 0, so both charges sit at
// z = 0 at t = 0 and a z-boost maps both to t' = 0.  Simultaneity is
// preserved exactly and the boosted initial data is the honest Lorentz
// transform of the rest-frame data, not an approximation of it.
//
// Under a boost the rate transforms as dP'^mu/dt' = Lambda^mu_nu
// (dP^nu/dt) / gamma, because dt' = gamma dt at the centre of mass.
// The residual reported is |measured - predicted| / |predicted|, over
//   step:       Accuracy::relativeTolerance 1e-5 ... 1e-10
//   resolution: directionCount 26, 50, 194, 302
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

int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double r0=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*r0),std::sqrt(k*r0)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const double beta=argc>1?std::atof(argv[1]):0.35;
    const double g=1.0/std::sqrt(1.0-beta*beta);
    // Moments along the normal, so the boost leaves them invariant and
    // the test isolates the charge/radiation sector.
    const Vec3 m1=Vec3{0,0,1}*firstMagneticMoment;
    const State rest=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    // Exact z-boost of the t = 0 slice: both charges are at z = 0, so
    // both map to t' = 0.  In-plane velocities pick up 1/gamma; the
    // new z-velocity is beta c for both.
    State moving=rest;
    moving.firstVelocity=rest.firstVelocity*(1.0/g)+Vec3{0,0,beta*c};
    moving.secondVelocity=rest.secondVelocity*(1.0/g)+Vec3{0,0,beta*c};

    std::printf("# boost beta = %.2f along the orbit normal, gamma ="
                " %.6f\n",beta,g);
    std::printf("# orbit at 0.25 a_pair, evolved a quarter period\n");
    std::printf("\n%10s %6s %7s %7s %15s %13s\n",
                "tolerance","dirs","radOnly","R/R0","|dP/dt| rest",
                "residual");
    // cfg = {tolerance, directions, radiationFieldOnly, radius/default}
    for(const std::array<double,4>& cfg:
            std::vector<std::array<double,4>>{
            {1e-5, 50,0,1},{1e-6, 50,0,1},{1e-7, 50,0,1},
            {1e-8, 50,0,1},{1e-9, 50,0,1},{1e-10,50,0,1},
            {1e-8, 26,0,1},{1e-8,194,0,1},{1e-8,302,0,1},
            {1e-8, 50,1,1},{1e-8,194,1,1},
            {1e-8, 50,0,100},{1e-8,194,1,100}}){
        if(argc>1&&!(cfg[0]==1e-8&&cfg[1]==50&&cfg[2]<0.5&&cfg[3]==1))
            continue;
        const double tol=cfg[0];
        const int dirs=static_cast<int>(cfg[1]);
        const bool radOnly=cfg[2]>0.5;
        const double radiusScale=cfg[3];
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=tol;
        acc.maximumDepth=24;
        const auto evolve=[&](const State& start,double span,State& out,
                              FieldFluxRates& flux)->bool{
            ClassicalTrajectoryEngine engine(start,acc);
            State s=start;
            const int steps=200;
            for(int i=0;i<steps;++i)
                if(!engine.advance(s,span/steps)) return false;
            out=s;
            FarFieldSampling sampling;
            sampling.directionCount=dirs;
            sampling.radiationFieldOnly=radOnly;
            sampling.controlRadius*=radiusScale;
            flux=electromagneticFieldFluxRates(s,engine.history(),sampling);
            return true;
        };
        State endRest,endMoving; FieldFluxRates fRest,fMoving;
        if(!evolve(rest,0.25*period,endRest,fRest)
           ||!evolve(moving,0.25*g*period,endMoving,fMoving)){
            std::printf("%10.0e %6d  evolution failed\n",tol,dirs);
            continue; }
        const Four rateRest{fRest.energy/c,fRest.momentum};
        const Four rateMoving{fMoving.energy/c,fMoving.momentum};
        // dP'^mu/dt' = Lambda(dP^mu/dt)/gamma
        // The pair was given velocity +beta c, so the frame the rest
        // quantities are carried INTO moves at -beta c: the boost
        // parameter is -beta.  Getting this backwards leaves a flat
        // 0.66 residual that does not move with step or resolution
        // and looks exactly like a covariance violation of the model.
        const Four boosted=boostAlongZ(rateRest,-beta);
        const Four predicted{boosted.time/g,boosted.space*(1.0/g)};
        std::printf("%10.0e %6d %7d %7.0f %15.8e %13.6e\n",tol,dirs,
            radOnly?1:0,radiusScale,
            std::sqrt(rateRest.time*rateRest.time
                     +rateRest.space.squaredNorm()),
            relativeGap(rateMoving,predicted));
    }
    std::printf("\n# A covariant scheme makes this fall with both\n");
    std::printf("# columns.  A floor that does not move is the\n");
    std::printf("# non-covariant part of the construction.\n");
}
