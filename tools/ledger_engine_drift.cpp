// Energy ledger along a REAL engine trajectory: retarded forces, reaction off,
// circular start, conservativeParticleEnergy tracked over whole orbits and
// normalized by |U_dd| of the moment configuration (audit section 102).
// Unlike ledger_near_barrier.cpp this uses the engine's own accumulated
// history, so it separates a real ledger inconsistency from an artefact of a
// history synthesized from one state.  Build it on two commits to compare,
// e.g. from a `git archive <commit>` tree.
//
// Usage: /tmp/engdrift <r/r*> <channel 0 para|1 ortho|2 no moments> <tilt deg>
//                      <orbits> <tolerance>
// Section 102 grid: r in 1 2 5 10, channel 0 1 2, tilt 0, 4 orbits, 1e-8.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/ledger_engine_drift.cpp -o /tmp/engdrift $(root-config --libs)
// Energy ledger along a real engine trajectory (retarded forces, no reaction):
// before and after e36bf5c.  argv: r/r* channel(0 para,1 ortho,2 bare) tiltDeg orbits tol
#include "modules/crem_trajectory.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  const double r=atof(argv[1])*comptonBarrierRadius; const int ch=atoi(argv[2]);
  const double tilt=atof(argv[3])*pi/180; const int orbits=atoi(argv[4]); const double tol=atof(argv[5]);
  const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
  const double v=std::sqrt(k/(mu*r)), period=2*pi*r/v;
  const Vec3 dir{std::sin(tilt),0,std::cos(tilt)};
  State s{};
  s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
  s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
  if(ch!=2){ s.firstProperDipole=dir*firstMagneticMoment;
             s.secondProperDipole=dir*(ch?-secondMagneticMoment:secondMagneticMoment); }
  synchronizeCovariantDipoles(s);
  // |U_dd| scale of the moment configuration, also used for the bare control.
  const double udd=std::abs(pairDipoleInteractionEnergy(
      s.firstPosition-s.secondPosition,dir*firstMagneticMoment,dir*secondMagneticMoment));
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=tol; acc.maximumDepth=20;
  acc.reactionModel=ChargeRadiationReactionModel::disabled; acc.computeOutwardFlux=false;
  acc.useRetardedExternalForces=true;
  ClassicalTrajectoryEngine eng(s,acc);
  const auto t0=std::chrono::steady_clock::now();
  const double E0=conservativeParticleEnergy(s,1.0); double emin=E0,emax=E0,rmin=1e300,rmax=0;
  const int steps=orbits*256; const double dt=period/256; int done=0; bool ok=true;
  for(int i=0;i<steps;++i){
    if(!eng.advance(s,dt)){ ok=false; break; }
    const double E=conservativeParticleEnergy(s,1.0); emin=std::min(emin,E); emax=std::max(emax,E);
    const double sep=(s.firstPosition-s.secondPosition).norm(); rmin=std::min(rmin,sep); rmax=std::max(rmax,sep);
    ++done;
  }
  const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  std::printf("r %.1f ch %s tilt %.0f tol %.0e | %s %d/%d | zakres E/|U_dd| %.3e | E_koniec-E0 /|U_dd| %+.3e | r %.4f..%.4f r* | %.0f s\n",
    r/comptonBarrierRadius,ch==0?"para":ch==1?"orto":"bez",atof(argv[3]),tol,ok?"ok":"AWARIA",done,steps,
    (emax-emin)/udd,(conservativeParticleEnergy(s,1.0)-E0)/udd,rmin/comptonBarrierRadius,rmax/comptonBarrierRadius,wall);
}
