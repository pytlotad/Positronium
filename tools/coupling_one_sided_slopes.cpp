// One-sided slopes of the coupling U = mu.B + p.E along the tangent line at
// the state where the engine first refuses a step (audit sections 88 and 103).
// A left/right ratio away from 1.0000 that does not converge with the probe
// step is a kink in time; section 88 found 2.12 and 5.78 before its fix, and
// section 103 found 5.84 when the magnetization term read the observation-
// time kinematics again.
//
// Usage: /tmp/kinkside <apo r*> <peri r*> <tilt deg>
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/coupling_one_sided_slopes.cpp -o /tmp/kinkside $(root-config --libs)
// One-sided slopes of the coupling U along the tangent line: left, right, and per field term.
// Also: which history segment each probe's retarded time falls in.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  const double apo=atof(argv[1])*comptonBarrierRadius, peri=atof(argv[2])*comptonBarrierRadius;
  const double t=atof(argv[3])*pi/180;
  const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
  const double a=0.5*(apo+peri), e=(apo-peri)/(apo+peri), va=std::sqrt(k/mu*(1-e)/apo);
  const double P=2*pi*std::sqrt(mu*a*a*a/k);
  const Vec3 dir{std::sin(t),0,std::cos(t)};
  State s{}; s.firstPosition={apo*secondMass/M,0,0}; s.secondPosition={-apo*firstMass/M,0,0};
  s.firstVelocity={0,va*secondMass/M,0}; s.secondVelocity={0,-va*firstMass/M,0};
  s.firstProperDipole=dir*firstMagneticMoment; s.secondProperDipole=dir*secondMagneticMoment;
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=1e-8; acc.maximumDepth=14;
  acc.reactionModel=ChargeRadiationReactionModel::disabled; acc.computeOutwardFlux=false;
  acc.useRetardedExternalForces=true;
  ClassicalTrajectoryEngine eng(s,acc);
  State x=s; bool failed=false; const int sub=256;
  for(int i=0;i<sub*3&&!failed;++i){ x=s; if(!eng.advance(s,P/sub)) failed=true; }
  const StateHistory h=eng.history();
  const Vec3 pos=x.firstPosition, vel=x.firstVelocity, mMag=x.firstDipole, mEle=x.firstElectricDipole;
  printf("stan %s r=%.4f r*, wezlow %zu\n",failed?"AWARII":"koncowy",separation(x)/comptonBarrierRadius,h.size());
  auto parts=[&](double offset,double& total,double& charge,double& dip,double& tret,int& seg){
    const Vec3 point=pos+vel*offset; const double time=x.time+offset;
    const ElectromagneticField lw=lienardWiechertField(point,time,h,x,false,secondCharge);
    const ElectromagneticField dp=retardedMagneticDipoleField(point,time,h,x,false);
    charge=dot(mMag,lw.magnetic)+dot(mEle,lw.electric);
    dip=dot(mMag,dp.magnetic)+dot(mEle,dp.electric);
    total=charge+dip;
    const Vec3 src=x.secondPosition; tret=time-(point-src).norm()/c;
    seg=-1; for(size_t i=1;i<h.size();++i) if(h[i-1].time<=tret&&tret<=h[i].time){ seg=(int)i; break; } };
  printf("  krok[s]      nachylenie lewe      nachylenie prawe     stosunek | ladunek L/P        dipol L/P\n");
  for(int j=0;j<6;++j){
    const double base=std::max(1.0e-4*separation(x),1.0e-3*nuclearCutoff)/c;
    const double step=base/std::pow(4.0,j);
    double t0,c0,d0,tr0; int s0; parts(0.0,t0,c0,d0,tr0,s0);
    double tm,cm,dm,trm; int sm; parts(-step,tm,cm,dm,trm,sm);
    double tp,cp,dp,trp; int sp; parts(step,tp,cp,dp,trp,sp);
    const double left=(t0-tm)/step, right=(tp-t0)/step;
    const double cl=(c0-cm)/step, cr=(cp-c0)/step, dl=(d0-dm)/step, dr=(dp-d0)/step;
    printf("  %.3e  %+.8e  %+.8e  %.4f | %+.3e/%+.3e  %+.3e/%+.3e\n",
      step,left,right,left/right,cl,cr,dl,dr);
    if(j==0) printf("      segmenty: lewy %d, srodek %d, prawy %d;  t_ret: %.9e %.9e %.9e\n",sm,s0,sp,trm,tr0,trp);
  }
}
