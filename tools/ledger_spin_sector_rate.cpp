// Spin-sector |dE/dt| along the flow on a Kepler orbit, moment angle measured
// from the orbital axis (audit section 92c-d).  Its scale, |U_dd| n from
// azimuthAveragedDipoleEnergy, is this probe's own and is NOT the scale of
// ledger_near_barrier.cpp / section 86f.
//
// Usage: /tmp/flowrate <periapsis r*> <tilt deg> <ecc> <retarded 1/0>
//                      <nodes> <h as fraction of the period>
// Section 92d: for t in 0 15 30 45 60 75 90; do /tmp/flowrate 0.5 $t 0.3 1 8 1e-6; done
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/ledger_spin_sector_rate.cpp -o /tmp/flowrate $(root-config --libs)

// Section 82b's flow check, rebuilt: is the 1-14% tilt-90 ledger residual of
// 86f still there after e36bf5c moved the magnetization term to the retarded
// time?  dE/dt of the full ledger by central difference of the flow map
// x' = v, p' = F, mu_proper' = omega_BMT x mu_proper.  Spin sector = the same
// state WITH moments minus the same state WITHOUT.  Normalizer |U_dd| is
// azimuthAveragedDipoleEnergy/reducedMass, as crem_collapse.hpp:2091 uses it.
// argv: radius(r*) tilt(deg) ecc retarded(1/0) nodes
#include "modules/crem_collapse.hpp"   // brings crem_trajectory.hpp
                                       // and azimuthAveragedDipoleEnergy
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  const double rStar=argc>1?atof(argv[1]):0.5;
  const double tilt=(argc>2?atof(argv[2]):0.0)*pi/180.0;
  const double ecc=argc>3?atof(argv[3]):0.3;
  const bool retarded=argc>4?atoi(argv[4])!=0:true;
  const int nodes=argc>5?atoi(argv[5]):16;
  const double hFrac=argc>6?atof(argv[6]):1.0e-6;
  const double mu=pairReducedMass, k=pairCoulombStrength, M=firstMass+secondMass;
  const double a=rStar*comptonBarrierRadius/(1.0-ecc);   // periapsis at rStar
  const double period=2*pi*std::sqrt(mu*a*a*a/k);
  // The 40-degree offset of 82b is NOT applied here, and that is deliberate.
  // With it, tilt 0 puts the moments at 40 degrees from the orbital axis and
  // tilt 90 puts them at 130 -- mirror images with the same |cos| = 0.766 --
  // so the two came out nearly equal (0.2687 against 0.2595) by construction
  // rather than by physics, which is useless for asking whether tilt 90 is
  // special.  tilt is the angle from the orbital axis, swept directly.
  const Vec3 dir{std::sin(tilt),0,std::cos(tilt)};

  auto stateAt=[&](double frac,bool withMoments){
    // Kepler node by eccentric anomaly, orbit in the xy plane.
    const double E=2*pi*frac;
    const double r=a*(1.0-ecc*std::cos(E));
    const double x=a*(std::cos(E)-ecc), y=a*std::sqrt(1-ecc*ecc)*std::sin(E);
    const double n=std::sqrt(k/(mu*a*a*a));
    const double vx=-a*n*std::sin(E)/(1-ecc*std::cos(E));
    const double vy=a*n*std::sqrt(1-ecc*ecc)*std::cos(E)/(1-ecc*std::cos(E));
    State s{};
    s.firstPosition={x*secondMass/M,y*secondMass/M,0};
    s.secondPosition={-x*firstMass/M,-y*firstMass/M,0};
    s.firstVelocity={vx*secondMass/M,vy*secondMass/M,0};
    s.secondVelocity={-vx*firstMass/M,-vy*firstMass/M,0};
    if(withMoments){
      s.firstProperDipole=dir*firstMagneticMoment;
      s.secondProperDipole=dir*secondMagneticMoment;
    }
    synchronizeCovariantDipoles(s);
    (void)r;
    return s;
  };

  // One explicit flow step of the ledger, evaluated as a central difference.
  auto ledgerRate=[&](const State& s0,double h){
    // The forces MUST be re-evaluated at the state being advanced, not at the
    // node.  Evaluating them once at `base` and stepping linearly is what made
    // E(+h)-E(-h) exactly proportional to h (measured: 2.037072e-20,
    // 1.018536e-20, 5.092681e-21, 2.546340e-21 on successive halvings, an
    // h-independent quotient of 7.110997e+04 to six digits), because composing
    // N such steps just walks N times along one straight line.  With the
    // evaluation inside, the composition is an Euler transport of the flow and
    // the quotient has something to converge to.
    auto advance=[&](const State& base,double step){
      State out=base;
      const StateHistory history=causalInitialHistory(base);
      const MutualForces F=retarded?retardedExternalForces(base,history)
                                   :allExternalForces(base);
      const DipoleDerivatives D=thomasBmtDipoleDerivatives(base,history);
      out.firstPosition=base.firstPosition+base.firstVelocity*step;
      out.secondPosition=base.secondPosition+base.secondVelocity*step;
      out.firstVelocity=velocityFromMomentum(
        momentum(base.firstVelocity,firstMass)+F.first*step,firstMass);
      out.secondVelocity=velocityFromMomentum(
        momentum(base.secondVelocity,secondMass)+F.second*step,secondMass);
      out.firstDipole=base.firstDipole+D.first*step;
      out.secondDipole=base.secondDipole+D.second*step;
      out.time=base.time+step;
      return out;
    };
    // The earlier version advanced by ONE explicit Euler step each way, which
    // makes E(+h)-E(-h) exactly linear in h: the difference quotient then
    // came out independent of h to every digit (measured: relative change
    // 0.000e+00 on halving), because it was differentiating along a straight
    // line rather than along the flow.  Composing N substeps instead makes
    // the transported state agree with the flow to O(h^2), so the quotient
    // has something to converge to and h means what it should.
    const double w=retarded?retardedDipoleSectorWeight(s0):0.0;
    const int substeps=32;
    auto transport=[&](double total){
        State out=s0;
        for(int i=0;i<substeps;++i) out=advance(out,total/substeps);
        return out;
    };
    const State plus=transport(h), minus=transport(-h);
    return (conservativeParticleEnergy(plus,w)
           -conservativeParticleEnergy(minus,w))/(2.0*h);
  };

  double sumAbs=0.0, maxAbs=0.0, uddSum=0.0;
  for(int i=0;i<nodes;++i){
    const double frac=(i+0.5)/nodes;
    const State withM=stateAt(frac,true), without=stateAt(frac,false);
    const double h=hFrac*period;
    const double spin=ledgerRate(withM,h)-ledgerRate(without,h);
    const Vec3 rel=withM.firstPosition-withM.secondPosition;
    const Vec3 relV=withM.firstVelocity-withM.secondVelocity;
    const Vec3 normal=cross(rel,relV);
    const double udd=std::abs(azimuthAveragedDipoleEnergy(
      rel.norm(),withM.firstDipole,withM.secondDipole,
      normal*(1.0/std::max(normal.norm(),1e-300))));
    sumAbs+=std::abs(spin); maxAbs=std::max(maxAbs,std::abs(spin));
    uddSum+=udd;
  }
  const double meanSpin=sumAbs/nodes, meanUdd=uddSum/nodes;
  // |U_dd| has units of energy; the rate is compared against |U_dd| per
  // orbital period, which is the scale 86f used ("relative to |U_dd| n").
  const double n=2*pi/period;
  printf("r_peri=%.2f r*  przechyl=%3.0f deg  ecc=%.2f  sily=%-12s | "
         "|dE/dt| spin: srednia %.4e W  max %.4e W | |U_dd| n = %.4e W | "
         "stosunek srednia/skala = %.4f  max/skala = %.4f\n",
    rStar,atof(argv[2]?argv[2]:"0"),ecc,retarded?"opoznione":"chwilowe",
    meanSpin,maxAbs,meanUdd*n,
    meanSpin/std::max(meanUdd*n,1e-300),maxAbs/std::max(meanUdd*n,1e-300));
}
