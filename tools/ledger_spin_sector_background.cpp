// Spin-sector |dE/dt| against the moment-free background of the same pair
// (audit section 92c): if the two were comparable, ratios taken from the
// difference would say nothing about the spin sector.
//
// Usage: /tmp/flowbg <periapsis r*> <tilt deg> <retarded 1/0> <nodes>
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/ledger_spin_sector_background.cpp -o /tmp/flowbg $(root-config --libs)

// Absolute scale control: how big is the ledger rate for a pair WITHOUT
// moments, against the spin-sector difference?  If they are comparable, the
// ratios reported from the difference say nothing about the spin sector.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  const double rStar=argc>1?atof(argv[1]):0.5;
  const double tiltDeg=argc>2?atof(argv[2]):0.0;
  const bool retarded=argc>3?atoi(argv[3])!=0:true;
  const int nodes=argc>4?atoi(argv[4]):8;
  const double mu=pairReducedMass,k=pairCoulombStrength,M=firstMass+secondMass;
  const double ecc=0.3, a=rStar*comptonBarrierRadius/(1.0-ecc);
  const double period=2*pi*std::sqrt(mu*a*a*a/k);
  const double tilt=tiltDeg*pi/180.0;
  const Vec3 dir{std::sin(tilt),0,std::cos(tilt)};
  auto advance=[&](const State& base,double step){
    State out=base;
    const StateHistory h=causalInitialHistory(base);
    const MutualForces F=retarded?retardedExternalForces(base,h)
                                 :allExternalForces(base);
    const DipoleDerivatives D=thomasBmtDipoleDerivatives(base,h);
    out.firstPosition=base.firstPosition+base.firstVelocity*step;
    out.secondPosition=base.secondPosition+base.secondVelocity*step;
    out.firstVelocity=velocityFromMomentum(
      momentum(base.firstVelocity,firstMass)+F.first*step,firstMass);
    out.secondVelocity=velocityFromMomentum(
      momentum(base.secondVelocity,secondMass)+F.second*step,secondMass);
    out.firstDipole=base.firstDipole+D.first*step;
    out.secondDipole=base.secondDipole+D.second*step;
    out.time=base.time+step; return out;
  };
  auto rate=[&](const State& s0,double h){
    State p=s0,m=s0;
    for(int i=0;i<32;++i){ p=advance(p,h/32); m=advance(m,-h/32); }
    const double w=retarded?1.0:0.0;
    return (conservativeParticleEnergy(p,w)-conservativeParticleEnergy(m,w))
           /(2.0*h);
  };
  double withSum=0,withoutSum=0,diffSum=0;
  for(int i=0;i<nodes;++i){
    const double E0=2*pi*(i+0.5)/nodes;
    const double x=a*(std::cos(E0)-ecc), y=a*std::sqrt(1-ecc*ecc)*std::sin(E0);
    const double n=std::sqrt(k/(mu*a*a*a));
    const double vx=-a*n*std::sin(E0)/(1-ecc*std::cos(E0));
    const double vy=a*n*std::sqrt(1-ecc*ecc)*std::cos(E0)/(1-ecc*std::cos(E0));
    State withM{},without{};
    for(State* s:{&withM,&without}){
      s->firstPosition={x*secondMass/M,y*secondMass/M,0};
      s->secondPosition={-x*firstMass/M,-y*firstMass/M,0};
      s->firstVelocity={vx*secondMass/M,vy*secondMass/M,0};
      s->secondVelocity={-vx*firstMass/M,-vy*firstMass/M,0};
    }
    withM.firstProperDipole=dir*firstMagneticMoment;
    withM.secondProperDipole=dir*secondMagneticMoment;
    synchronizeCovariantDipoles(withM); synchronizeCovariantDipoles(without);
    const double h=1.0e-6*period;
    const double rw=rate(withM,h), rb=rate(without,h);
    withSum+=std::abs(rw); withoutSum+=std::abs(rb); diffSum+=std::abs(rw-rb);
  }
  printf("peri %.2f r*  kat %3.0f deg  sily %-10s | |dE/dt| z momentami %.4e W  "
         "BEZ momentow %.4e W  roznica (sektor spinowy) %.4e W | "
         "sektor/bez = %.4f\n",
    rStar,tiltDeg,retarded?"opoznione":"chwilowe",
    withSum/nodes,withoutSum/nodes,diffSum/nodes,
    (diffSum/nodes)/std::max(withoutSum/nodes,1e-300));
}
