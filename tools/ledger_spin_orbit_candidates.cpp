// Candidate spin-orbit terms of the energy ledger, checked along the model's
// own flow (audit section 82c).  dE/dt of the ledger by central difference of
// the flow map x' = v, p' = F, mu_proper' = omega_BMT x mu_proper; spin sector
// = the same state with moments minus without.  Compares no term, the old
// -mu.B(v_other) and -p.E, for instantaneous and retarded forces.
//
// Usage: /tmp/ledger <pair|-> <forces: 0 instantaneous, 2 retarded causal>
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/ledger_spin_orbit_candidates.cpp -o /tmp/ledger $(root-config --libs)

// Candidate spin-orbit ledger terms along the model's own flow.
// E_base = conservativeParticleEnergy - U_so (no charge-dipole term).
// candidates: +U_so (current), +U_p and -U_p with U_p = sum p_i . E_j(r_i) (motional electric dipole), 0.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
static int gMode=0; // 0 instantaneous, 2 retarded causal
static State flow(const State& s,double h){
  State x=s; const StateHistory hist=gMode==2?causalInitialHistory(s):StateHistory{State{s}};
  const MutualForces F=gMode?retardedExternalForces(s,hist):allExternalForces(s);
  const LocalElectromagneticFields fl=localRelativisticFields(s,StateHistory{State{s}});
  x.time=s.time+h;
  x.firstPosition=s.firstPosition+s.firstVelocity*h; x.secondPosition=s.secondPosition+s.secondVelocity*h;
  x.firstVelocity=velocityFromMomentum(momentum(s.firstVelocity,firstMass)+F.first*h,firstMass);
  x.secondVelocity=velocityFromMomentum(momentum(s.secondVelocity,secondMass)+F.second*h,secondMass);
  const Vec3 w1=thomasBmtEffectiveField(s.firstVelocity,fl.atFirst,firstGFactor)*(-firstCharge/firstMass);
  const Vec3 w2=thomasBmtEffectiveField(s.secondVelocity,fl.atSecond,secondGFactor)*(-secondCharge/secondMass);
  x.firstProperDipole=s.firstProperDipole+cross(w1,s.firstProperDipole)*h;
  x.secondProperDipole=s.secondProperDipole+cross(w2,s.secondProperDipole)*h;
  synchronizeCovariantDipoles(x); return x;
}
static double motionalDipoleEnergy(const State& s){ // sum_i p_i . E_j(r_i), Plummer Coulomb field
  const Vec3 d=s.firstPosition-s.secondPosition; const double fl=separationFloor();
  const double rho2=d.squaredNorm()+fl*fl; const double f=1.0/std::pow(rho2,1.5);
  const Vec3 E2at1=d*(coulomb*secondCharge*f), E1at2=d*(-coulomb*firstCharge*f);
  return dot(s.firstElectricDipole,E2at1)+dot(s.secondElectricDipole,E1at2);
}
int main(int argc,char** argv){
  if(argc>1&&std::string(argv[1])!="-") applyPairFromOption(argv[1]);
  gMode=argc>2?atoi(argv[2]):0;
  const double k=pairCoulombStrength,mu=pairReducedMass;
  printf("# pair %s  forces %s\n",argc>1?argv[1]:"e-,e+",gMode?"retarded (causal history)":"instantaneous");
  for(double f:{1.0,0.1}) for(int ch:{0,1,2,3}){
    const double e=0.3,a=f*pairBohrRadius(activePair),ec=std::sqrt(1-e*e),n=std::sqrt(k/(mu*a*a*a)),M=firstMass+secondMass;
    const double tt=40*pi/180; const Vec3 dir{std::sin(tt)*std::cos(0.6),std::sin(tt)*std::sin(0.6),std::cos(tt)};
    Vec3 m1=dir*firstMagneticMoment, m2=dir*(ch==1?-secondMagneticMoment:secondMagneticMoment);
    if(ch==2) m2={}; if(ch==3) m1={};
    const char* nm[4]={"oba rownolegle","oba antyrownolegle","tylko pierwszy","tylko drugi"};
    const Vec3 R{1,0,0},T{0,1,0}; const int N=64;
    double r0=0,rSo=0,rPp=0,rPm=0,sc=0;
    for(int i=0;i<N;++i){ double E=2*pi*i/N,co=std::cos(E),si=std::sin(E),w=1-e*co;
      Vec3 r=R*(a*(co-e))+T*(a*ec*si), v=(R*(-a*n*si)+T*(a*n*ec*co))/w;
      State s{}; s.firstPosition=r*(secondMass/M); s.secondPosition=r*(-firstMass/M); s.firstVelocity=v*(secondMass/M); s.secondVelocity=v*(-firstMass/M);
      State b=s; synchronizeCovariantDipoles(b); s.firstProperDipole=m1; s.secondProperDipole=m2; synchronizeCovariantDipoles(s);
      const double h=1e-6*(2*pi/n);
      State sp=flow(s,h), sm=flow(s,-h), bp=flow(b,h), bm=flow(b,-h);
      auto D=[&](double (*g)(const State&)){ return (g(sp)-g(sm))/(2*h); };
      const double dBase=D([](const State& x){return conservativeParticleEnergy(x)-chargeDipoleInteractionEnergy(x);})
                        -(conservativeParticleEnergy(bp)-conservativeParticleEnergy(bm))/(2*h);
      const double dSo=D(chargeDipoleInteractionEnergy), dP=D(motionalDipoleEnergy);
      r0+=std::abs(dBase)*w; rSo+=std::abs(dBase+dSo)*w; rPp+=std::abs(dBase+dP)*w; rPm+=std::abs(dBase-dP)*w; sc+=std::max(std::abs(dSo),std::abs(dP))*w;
    }
    printf("f %-4g %-18s | skala %.3e W | reszta: bez czlonu %.3e  +U_so %.3e  +p.E %.3e  -p.E %.3e\n",f,nm[ch],sc/N,r0/N,rSo/N,rPp/N,rPm/N);
  }
}
