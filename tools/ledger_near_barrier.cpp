// Spin-sector |dE/dt| of the full ledger near r*, relative to |U_dd| n
// (audit sections 86f and 101).  The scale is defined HERE: |U_dd| averaged
// over 16 azimuths of a CIRCULAR orbit, times its angular frequency n; the
// moments point along (0.8 sin t, 0.6 sin t, cos t).  Section 92 compared a
// different probe (Kepler e = 0.3, other moment convention) and wrongly
// concluded this scale was undefined.  Retarded forces, causal history.
//
// Usage: /tmp/deepledger          (prints r = 2, 1, 0.5 r*, para/ortho,
//                                  tilt 0/90)
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/ledger_near_barrier.cpp -o /tmp/deepledger $(root-config --libs)

// Ledger consistency near r*: spin-sector dE/dt along the model's own flow (retarded, causal history).
#include "modules/crem_trajectory.hpp"
#include <cstdio>
static State flow(const State& s,double h){
  State x=s; const StateHistory hist=causalInitialHistory(s);
  const MutualForces F=retardedExternalForces(s,hist);
  const LocalElectromagneticFields fl=localRelativisticFields(s,StateHistory{State{s}});
  x.time=s.time+h; x.firstPosition=s.firstPosition+s.firstVelocity*h; x.secondPosition=s.secondPosition+s.secondVelocity*h;
  x.firstVelocity=velocityFromMomentum(momentum(s.firstVelocity,firstMass)+F.first*h,firstMass);
  x.secondVelocity=velocityFromMomentum(momentum(s.secondVelocity,secondMass)+F.second*h,secondMass);
  const Vec3 w1=thomasBmtEffectiveField(s.firstVelocity,fl.atFirst,firstGFactor)*(-firstCharge/firstMass);
  const Vec3 w2=thomasBmtEffectiveField(s.secondVelocity,fl.atSecond,secondGFactor)*(-secondCharge/secondMass);
  x.firstProperDipole=s.firstProperDipole+cross(w1,s.firstProperDipole)*h; x.secondProperDipole=s.secondProperDipole+cross(w2,s.secondProperDipole)*h;
  synchronizeCovariantDipoles(x); return x;
}
int main(){
  const double k=pairCoulombStrength,mu=pairReducedMass,rs=comptonBarrierRadius,M=firstMass+secondMass;
  printf("# promien momentu %.4f r*\n",magneticDipoleRadius()/rs);
  for(double rr:{2.0,1.0,0.5}) for(int ch:{0,1}) for(double td:{0.0,90.0}){
    const double r=rr*rs, v=std::sqrt(k/(mu*r)); const double t=td*pi/180; const Vec3 dir{std::sin(t)*0.8,std::sin(t)*0.6,std::cos(t)};
    const int N=16; double res=0,scale=0,dd=0;
    for(int i=0;i<N;++i){ double ph=2*pi*i/N; Vec3 rh{std::cos(ph),std::sin(ph),0}, th{-std::sin(ph),std::cos(ph),0};
      State s{}; s.firstPosition=rh*(r*secondMass/M); s.secondPosition=rh*(-r*firstMass/M); s.firstVelocity=th*(v*secondMass/M); s.secondVelocity=th*(-v*firstMass/M);
      State b=s; synchronizeCovariantDipoles(b); s.firstProperDipole=dir*firstMagneticMoment; s.secondProperDipole=dir*(ch?-secondMagneticMoment:secondMagneticMoment); synchronizeCovariantDipoles(s);
      const double h=1e-7*(2*pi*r/v);
      State sp=flow(s,h), sm=flow(s,-h), bp=flow(b,h), bm=flow(b,-h);
      const double dS=(conservativeParticleEnergy(sp)-conservativeParticleEnergy(sm))/(2*h), dB=(conservativeParticleEnergy(bp)-conservativeParticleEnergy(bm))/(2*h);
      res+=std::abs(dS-dB);
      const Vec3 d=s.firstPosition-s.secondPosition; const Vec3 Fdd=pairDipoleForce(d,s.firstDipole,s.secondDipole);
      scale+=std::abs(dot(Fdd,s.firstVelocity-s.secondVelocity)); dd+=std::abs(pairDipoleInteractionEnergy(d,s.firstDipole,s.secondDipole));
    }
    const double n=2*pi*v/(2*pi*r);
    printf("r=%.1f r* %s tilt %2.0f | |dE/dt|_spin %.3e W | skala |U_dd| n %.3e W | wzgl. %.2e\n",rr,ch?"orto":"para",td,res/N,dd/N*n,res/(dd/N*n*N)*1.0);
  }
}
