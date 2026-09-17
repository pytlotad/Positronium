// Outer time-step laws at the bottom of a deep orbit, driven through the
// shared regularizedTimeStep (audit section 89): constant dt, dt ~ r^{3/2},
// dt ~ r^2.  Set CREM_STEP_CENSUS=1 in the environment for evaluation counts.
//
// Usage: /tmp/steplaw <apo r*> <peri r*> <tilt deg> <tolerance>
//                     <law: 0 constant, 1.5, 2.0> <steps per orbit> <depth>
//                     <windows>
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/time_step_law_compare.cpp -o /tmp/steplaw $(root-config --libs)

// SHARED-FUNCTION variant: same comparison driven through regularizedTimeStep.
// Constant-time vs Sundman stepping (dt ~ r^p) on the same deep orbit.
// argv: apo peri tilt tol exponent(0=stale, 1.5, 2.0) stepsPerLocalOrbit depth orbits
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <cstdlib>
int main(int argc,char** argv){
  const double apo=atof(argv[1])*comptonBarrierRadius, peri=atof(argv[2])*comptonBarrierRadius;
  const double t=atof(argv[3])*pi/180, tol=atof(argv[4]), expo=atof(argv[5]);
  const double perOrbit=atof(argv[6]); const int depth=atoi(argv[7]), orbits=atoi(argv[8]);
  const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
  const double a=0.5*(apo+peri), e=(apo-peri)/(apo+peri), va=std::sqrt(k/mu*(1-e)/apo);
  const double P=2*pi*std::sqrt(mu*a*a*a/k);
  const Vec3 dir{std::sin(t),0,std::cos(t)};
  State s{}; s.firstPosition={apo*secondMass/M,0,0}; s.secondPosition={-apo*firstMass/M,0,0};
  s.firstVelocity={0,va*secondMass/M,0}; s.secondVelocity={0,-va*firstMass/M,0};
  s.firstProperDipole=dir*firstMagneticMoment; s.secondProperDipole=dir*secondMagneticMoment;
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=tol; acc.maximumDepth=depth;
  acc.reactionModel=ChargeRadiationReactionModel::disabled; acc.computeOutwardFlux=false;
  acc.useRetardedExternalForces=true;
  ClassicalTrajectoryEngine eng(s,acc);
  const auto t0=std::chrono::steady_clock::now();
  const double target=orbits*P; double rmin=1e300; long steps=0; bool ok=true;
  int periapses=0; double previous=separation(s); bool falling=true;
  double smallestRequested=1e300, smallestAt=0, smallestL=0, rmax=0.0;
  auto keplerEnergy=[&](const State& st){
      const Vec3 v=st.firstVelocity-st.secondVelocity;
      return 0.5*pairReducedMass*dot(v,v)-pairCoulombStrength/separation(st); };
  const double energy0=keplerEnergy(s);
  const double angular0=mu*cross(s.firstPosition-s.secondPosition,
                                 s.firstVelocity-s.secondVelocity).norm();
  while(s.time<target&&ok){
    const double r=separation(s);
    // Through the SHARED module function, not the probe's own expression:
    // this is the code the production loop now runs.
    const RegularizedStep rule{perOrbit, 5.0e-18,
        expo>1.75 ? TimeRegularizationLaw::constantAngle
                  : TimeRegularizationLaw::localOrbit};
    double dt = expo<=0.0 ? P/perOrbit
                          : regularizedTimeStep(s, target-s.time, rule);
    if(!(dt>0.0)||!std::isfinite(dt)) { ok=false; break; }
    const double requested=std::min(dt,target-s.time);
    if(requested<smallestRequested){ smallestRequested=requested;
        smallestAt=r/comptonBarrierRadius;
        smallestL=mu*cross(s.firstPosition-s.secondPosition,
                           s.firstVelocity-s.secondVelocity).norm()/hbar; }
    if(!eng.advance(s,requested)) { ok=false; break; }
    ++steps; const double now=separation(s); rmin=std::min(rmin,now); rmax=std::max(rmax,now);
    if(falling&&now>previous){ ++periapses; falling=false; }
    else if(!falling&&now<previous) falling=true;
    previous=now;
  }
  const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  printf("wykladnik %-5s krokow/obieg %5.0f gleb %2d | %-9s krokow %7ld  ukonczone %.3f orbity  perycentrow %d  r_min %.4f r*  czas %.1f s\n",
    expo<=0.0?"staly":(expo>1.75?"r^2":"r^1.5"),perOrbit,depth,ok?"cale okno":"AWARIA",steps,
    s.time/P,periapses,rmin/comptonBarrierRadius,wall);
  const double angular1=mu*cross(s.firstPosition-s.secondPosition,
                                 s.firstVelocity-s.secondVelocity).norm();
  printf("        r_max %.3f r*  dryf E_kepler %+.3e wzgl.  L koncowe %.4f hbar (start %.4f)\n",
    rmax/comptonBarrierRadius,(keplerEnergy(s)-energy0)/std::fabs(energy0),
    angular1/hbar,angular0/hbar);
  printf("        najmniejszy krok ZEWNETRZNY %.3e s przy r=%.4f r*, L=%.4f hbar\n",
    smallestRequested,smallestAt,smallestL);
  printf("        wywolan kroku: %llu, sredni dt %.3e s, najmniejszy %.3e s\n",
    (unsigned long long)gStepCensusCount,
    gStepCensusCount?gStepCensusTotalTime/(double)gStepCensusCount:0.0,
    gStepCensusSmallest);
}
