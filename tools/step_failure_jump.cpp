// Where an accuracy failure of the engine comes from (audit section 107).
// Replays ClassicalTrajectoryEngine's adaptive recursion step for step
// (reaction off, retarded forces, flux off, as ledger_engine_drift.cpp) until
// a step misses the tolerance at the maximum depth, then takes the failing
// start state, history and step and scans the endpoint of a single coarse step
// over tau in (0, 2 dt].  At every tau it prints, for both particles, the
// moment force covariantDipoleGradientForce, the material rate, the index of
// the history segment the gradient stencil is pinned to, and the index of the
// segment holding the partner's retarded time, so a jump in the force can be
// matched to the switch that causes it.  Each row also carries the partner's
// charge-field magnitude at particle 1, the retarded lag and the history span.
// It then bisects the first jump of |G1| and prints the six stencil probes on
// both sides: the two-pole cancellation ratio and the coupling U with the
// production pole separation and with the retreat.  With JUMP_SPLIT=1 and
// CREM_DEBUG_GRAD=1, and a binary built by tools/build_variant.sh
// --probe-switches, U is split into the charge, pole and magnetization fields.
// Section 107 additionally printed each pole's moment derivatives and stencil
// branch from a one-off patch inside twoChargeLimitDipoleField (not kept).
//
// Usage: step_failure_jump <r/r*> <channel 0 para|1 ortho> <tilt deg>
//                          <tolerance> <depth> [tau points, default 400]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/step_failure_jump.cpp -o /tmp/jump $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>

namespace {
double normalizedStepError(const State& coarse,const State& fine) {
    const double lengthScale=std::max(separation(fine),collisionBoundaryRadius);
    const double speedScale=std::max(
        (fine.firstVelocity-fine.secondVelocity).norm(),1.0e-6*c);
    return std::max({
        (coarse.firstPosition-fine.firstPosition).norm()/lengthScale,
        (coarse.secondPosition-fine.secondPosition).norm()/lengthScale,
        (coarse.firstVelocity-fine.firstVelocity).norm()/speedScale,
        (coarse.secondVelocity-fine.secondVelocity).norm()/speedScale});
}
struct Failure { bool found=false; State start; StateHistory history; double dt=0; double error=0; };
Failure gFailure;
double gTolerance=1e-8; int gDepth=20;
const auto kModel=ChargeRadiationReactionModel::disabled;

bool adaptive(const State& start,const StateHistory& history,double dt,int depth,
              State& accepted,StateHistory& acceptedHistory) {
    State coarse=start;
    integrateElectrodynamicStep(coarse,dt,history,false,kModel,true);
    State fine=start; StateHistory fineHistory=history;
    bool finite=isFinite(coarse);
    if(finite) {
        integrateElectrodynamicStep(fine,0.5*dt,fineHistory,false,kModel,true);
        finite=isFinite(fine);
        if(finite) {
            appendStateHistory(fineHistory,fine);
            integrateElectrodynamicStep(fine,0.5*dt,fineHistory,false,kModel,true);
            finite=isFinite(fine);
        }
    }
    const double error=finite?normalizedStepError(coarse,fine):NAN;
    if(finite&&error<=gTolerance) {
        appendStateHistory(fineHistory,fine);
        accepted=fine; acceptedHistory=std::move(fineHistory);
        return true;
    }
    if(depth>=gDepth) {
        gFailure={true,start,history,dt,error};
        return false;
    }
    State mid; StateHistory midHistory;
    if(!adaptive(start,history,0.5*dt,depth+1,mid,midHistory)) return false;
    return adaptive(mid,midHistory,0.5*dt,depth+1,accepted,acceptedHistory);
}

std::size_t segmentOf(const StateHistory& h,double t) {
    return static_cast<std::size_t>(std::lower_bound(h.begin(),h.end(),t,
        [](const State& s,double x){ return s.time<x; })-h.begin());
}
double retardedTimeOf(const StateHistory& h,const State& st,bool sourceIsFirst,const Vec3& at) {
    double t=st.time;
    for(int i=0;i<64;++i) {
        const ChargeKinematics k=historicalCharge(h,st,sourceIsFirst,t);
        const Vec3 d=at-k.position; const double r=d.norm();
        const double res=t+r/c-st.time;
        t-=res/std::max(1e-8,1.0-dot(d/r,k.velocity/c));
        if(std::abs(res)<1e-16*std::abs(st.time)) break;
    }
    return t;
}
}

int main(int argc,char** argv) {
    if(argc<6) { std::fprintf(stderr,"usage: r ch tilt tol depth [points]\n"); return 2; }
    const double r=atof(argv[1])*comptonBarrierRadius; const int ch=atoi(argv[2]);
    const double tilt=atof(argv[3])*pi/180; gTolerance=atof(argv[4]); gDepth=atoi(argv[5]);
    const int points=argc>6?atoi(argv[6]):400;
    const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
    const double v=std::sqrt(k/(mu*r)),period=2*pi*r/v;
    const Vec3 dir{std::sin(tilt),0,std::cos(tilt)};
    State s{};
    s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
    s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(ch?-secondMagneticMoment:secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    StateHistory history=causalInitialHistory(s);
    const double dt=period/256;
    int step=0;
    for(;step<256;++step) {
        State next; StateHistory nextHistory;
        if(!adaptive(s,history,dt,0,next,nextHistory)) break;
        s=next; history=std::move(nextHistory);
    }
    if(!gFailure.found) { std::printf("no failure in 256 steps\n"); return 0; }
    const Failure f=gFailure;
    std::printf("failure in outer step %d: r %.5f r*, beta1 %.4f, dt %.4e s, error %.3e, history %zu nodes\n",
        step+1,separation(f.start)/comptonBarrierRadius,f.start.firstVelocity.norm()/c,f.dt,f.error,f.history.size());
    const double gradientStep=std::max(1.0e-4*separation(f.start),1.0e-3*nuclearCutoff);
    std::printf("tau/dt  |G1| |G2| (N)  rate1 rate2  pin1 pin2  seg_ret1 seg_ret2  frac_ret1\n");
    double prevG1=0,prevG2=0; long prevPin1=-1,prevPin2=-1,prevSeg1=-1,prevSeg2=-1;
    for(int i=1;i<=points;++i) {
        const double tau=2.0*f.dt*i/points;
        State e=f.start;
        integrateElectrodynamicStep(e,tau,f.history,false,kModel,true);
        const Vec3 G1=covariantDipoleGradientForce(e,f.history,true);
        const Vec3 G2=covariantDipoleGradientForce(e,f.history,false);
        const double rate1=dipoleCouplingMaterialRate(e,f.history,true);
        const double rate2=dipoleCouplingMaterialRate(e,f.history,false);
        const RetardedSegmentPin p1=retardedSegmentPinAt(f.history,e,false,e.firstPosition,e.time,8.0*gradientStep/c);
        const RetardedSegmentPin p2=retardedSegmentPinAt(f.history,e,true,e.secondPosition,e.time,8.0*gradientStep/c);
        const long pin1=p1.history?static_cast<long>(p1.newerIndex):-1;
        const long pin2=p2.history?static_cast<long>(p2.newerIndex):-1;
        const double tr1=retardedTimeOf(f.history,e,false,e.firstPosition);
        const double tr2=retardedTimeOf(f.history,e,true,e.secondPosition);
        const long seg1=static_cast<long>(segmentOf(f.history,tr1));
        const long seg2=static_cast<long>(segmentOf(f.history,tr2));
        double frac1=NAN;
        if(seg1>0&&seg1<static_cast<long>(f.history.size()))
            frac1=(tr1-f.history[seg1-1].time)/(f.history[seg1].time-f.history[seg1-1].time);
        const double lw1=lienardWiechertField(e.firstPosition,e.time,f.history,e,false,secondCharge).electric.norm();
        const bool jump=i>1&&(std::abs(G1.norm()-prevG1)>1e-3*prevG1||std::abs(G2.norm()-prevG2)>1e-3*prevG2
            ||pin1!=prevPin1||pin2!=prevPin2||seg1!=prevSeg1||seg2!=prevSeg2);
        if(i==1||i==points||jump||i%(points/8)==0)
            std::printf("%.4f %s %.6e %.6e  %+.6e %+.6e  %ld %ld  %ld %ld  %.6f  |E_LW1| %.6e  t-t_ret1 %.6e  hist %zu front %.6e back %.6e\n",tau/f.dt,jump?"*":" ",
                G1.norm(),G2.norm(),rate1,rate2,pin1,pin2,seg1,seg2,frac1,lw1,e.time-tr1,f.history.size(),f.history.front().time,f.history.back().time);
        prevG1=G1.norm(); prevG2=G2.norm(); prevPin1=pin1; prevPin2=pin2; prevSeg1=seg1; prevSeg2=seg2;
    }
    // Bisect the first jump of |G1| and open the six stencil probes on both
    // sides: cancellation ratio of the two-pole field as evaluated, and the
    // coupling with the production separation and with the retreat.
    auto endpoint=[&](double tau){ State e=f.start; integrateElectrodynamicStep(e,tau,f.history,false,kModel,true); return e; };
    auto g1=[&](double tau){ const State e=endpoint(tau); return covariantDipoleGradientForce(e,f.history,true).norm(); };
    double lo=0.0,hi=0.0; double g0=g1(2.0*f.dt/points);
    for(int i=2;i<=points;++i){ const double t=2.0*f.dt*i/points; if(std::abs(g1(t)-g0)>1e-3*g0){ lo=2.0*f.dt*(i-1)/points; hi=t; break; } }
    if(!(hi>0.0)) { std::printf("no jump in |G1|\n"); return 0; }
    for(int it=0;it<60;++it){ const double mid=0.5*(lo+hi); if(std::abs(g1(mid)-g0)>1e-3*g0) hi=mid; else lo=mid; if(hi-lo<=1e-12*f.dt) break; }
    std::printf("\njump of |G1| between tau/dt %.12f and %.12f: %.6e -> %.6e\n",lo/f.dt,hi/f.dt,g1(lo),g1(hi));
    for(double tau:{lo,hi}) {
        const State e=endpoint(tau);
        const RetardedSegmentPinGuard guard(retardedSegmentPinAt(f.history,e,false,e.firstPosition,e.time,8.0*gradientStep/c));
        std::printf("tau/dt %.12f\n",tau/f.dt);
        for(int axis=0;axis<3;++axis) for(int sign:{+1,-1}) {
            Vec3 offset; (axis==0?offset.x:axis==1?offset.y:offset.z)=sign*gradientStep;
            const Vec3 point=e.firstPosition+offset;
            const ElectromagneticField field=fieldFromOtherParticleAt(point,e.time,e,f.history,true);
            const double ratio=gPoleCancellationRatio;
            const ElectromagneticField retreated=fieldFromOtherParticleAt(point,e.time,e,f.history,true,1.0e-7);
            const double ratioRetreated=gPoleCancellationRatio;
            const double U=dot(e.firstDipole,field.magnetic)+dot(e.firstElectricDipole,field.electric);
            const double Ur=dot(e.firstDipole,retreated.magnetic)+dot(e.firstElectricDipole,retreated.electric);
            std::printf("  axis %d%c ratio %.4e retreated %.4e  U %.12e  U_retreated %.12e  %s\n",axis,sign>0?'+':'-',
                ratio,ratioRetreated,U,Ur,ratio>30.0?(ratioRetreated<=30.0?"RETREAT USED":"retreat rejected"):"");
            // With the probe switches (tools/build_variant.sh --probe-switches)
            // split the coupling into the partner's charge, pole and
            // magnetization fields, and read the two-pole construction's trace
            // (needs CREM_DEBUG_GRAD set in the environment).
            if(std::getenv("JUMP_SPLIT")) {
                const auto part=[&](std::initializer_list<const char*> off){
                    for(const char* n:off) setenv(n,"1",1);
                    const ElectromagneticField fp=fieldFromOtherParticleAt(point,e.time,e,f.history,true);
                    for(const char* n:off) unsetenv(n);
                    return dot(e.firstDipole,fp.magnetic)+dot(e.firstElectricDipole,fp.electric); };
                const double uLW=part({"CREM_PROBE_NO_POLES","CREM_PROBE_NO_MAGNETIZATION"});
                const double uPoles=part({"CREM_PROBE_NO_LW_IN_MOMENT","CREM_PROBE_NO_MAGNETIZATION"});
                const TwoChargeLimitTrace tr=gTwoChargeTrace[(gTwoChargeTraceCount-1)&7];
                const double uMag=part({"CREM_PROBE_NO_LW_IN_MOMENT","CREM_PROBE_NO_POLES"});
                std::printf("      U_charge %.12e  U_poles %.12e  U_magnetization %.12e\n",uLW,uPoles,uMag);
                std::printf("      poles: central iterations %d residual %.3e  pole charge %.6e  moment scale %.6e  |pole| %.6e |sum| %.6e early %d\n",
                    tr.centralIterations,tr.centralResidual,tr.poleCharge,tr.momentScale,tr.poleMagnitude,tr.sumMagnitude,tr.earlyReturn?1:0);
            }
        }
    }
}
