// Energy ledger along a REAL engine trajectory: retarded forces, reaction off,
// circular start, conservativeParticleEnergy tracked over whole orbits and
// normalized by |U_dd| of the moment configuration (audit section 102).
// Unlike ledger_near_barrier.cpp this uses the engine's own accumulated
// history, so it separates a real ledger inconsistency from an artefact of a
// history synthesized from one state.  Build it on two commits to compare,
// e.g. from a `git archive <commit>` tree.
//
// Usage: /tmp/engdrift <r/r*> <channel 0 para|1 ortho|2 no moments> <tilt deg>
//                      <orbits> <tolerance> [steps per orbit, default 256]
//                      [maximumDepth, default 20] [trace every N steps, 0 = off]
//                      [retarded forces 1/0, default 1]
//                      [rotation of the second moment about y, degrees]
// Trace lines (audit section 104) split the ledger into its terms, in |U_dd|
// units relative to their start values: kinetic, Coulomb, dipole-dipole,
// charge-dipole (-p.E), Darwin, dipole constraint, then the three components of
// the total momentum; a KONCOWY line gives the momentum where the run stopped
// (section 106).
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
             s.secondProperDipole=dir*(ch?-secondMagneticMoment:secondMagneticMoment);
             // Optional 11th argument (section 106): rotate the second moment by
             // delta degrees about y, so ortho (delta 0) turns into para (180).
             if(argc>10){ const double d=atof(argv[10])*pi/180; const Vec3 m=s.secondProperDipole;
               s.secondProperDipole={m.x*std::cos(d)+m.z*std::sin(d),m.y,-m.x*std::sin(d)+m.z*std::cos(d)}; } }
  synchronizeCovariantDipoles(s);
  // |U_dd| scale of the moment configuration, also used for the bare control.
  const double udd=std::abs(pairDipoleInteractionEnergy(
      s.firstPosition-s.secondPosition,dir*firstMagneticMoment,dir*secondMagneticMoment));
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=tol; acc.maximumDepth=argc>7?atoi(argv[7]):20;
  acc.reactionModel=ChargeRadiationReactionModel::disabled; acc.computeOutwardFlux=false;
  acc.useRetardedExternalForces=!(argc>9&&atoi(argv[9])==0);
  ClassicalTrajectoryEngine eng(s,acc);
  const auto t0=std::chrono::steady_clock::now();
  const double E0=conservativeParticleEnergy(s,1.0); double emin=E0,emax=E0,rmin=1e300,rmax=0;
  const int perOrbit=argc>6?atoi(argv[6]):256; const int traceEvery=argc>8?atoi(argv[8]):0;
  const int steps=orbits*perOrbit; const double dt=period/perOrbit;
  struct Terms{double kin,coul,dd,cd,dar,con;};
  auto terms=[&](const State& st){
    const PairGeometry g=clampedPairGeometry(st);
    return Terms{kineticEnergy(st.firstVelocity,firstMass)+kineticEnergy(st.secondVelocity,secondMass),
      -pairCoulombStrength*g.inverseDistance,
      pairDipoleInteractionEnergy(st.firstPosition-st.secondPosition,st.firstDipole,st.secondDipole),
      chargeDipoleInteractionEnergy(st),darwinInteractionEnergy(st),st.dipoleConstraintEnergy};
  };
  const Terms T0=terms(s); int done=0; bool ok=true;
  for(int i=0;i<steps;++i){
    if(!eng.advance(s,dt)){ ok=false; break; }
    const double E=conservativeParticleEnergy(s,1.0); emin=std::min(emin,E); emax=std::max(emax,E);
    const double sep=(s.firstPosition-s.secondPosition).norm(); rmin=std::min(rmin,sep); rmax=std::max(rmax,sep);
    ++done;
    if(traceEvery>0&&done%traceEvery==0){
      const Terms T=terms(s);
      const double beta1=s.firstVelocity.norm()/c;
      // Centre-of-mass motion versus relative motion: a pair that gains
      // energy while its separation barely grows is being pushed as a whole.
      const Vec3 totalMomentum=momentum(s.firstVelocity,firstMass)+momentum(s.secondVelocity,secondMass);
      const Vec3 vcm=(s.firstVelocity*firstMass+s.secondVelocity*secondMass)*(1.0/(firstMass+secondMass));
      const Vec3 vrel=s.firstVelocity-s.secondVelocity;
      const double kcm=0.5*(firstMass+secondMass)*vcm.squaredNorm();
      const double krel=0.5*pairReducedMass*vrel.squaredNorm();
      std::printf("TRACE step %d t/P %.4f r %.4f r* beta %.4f beta_cm %.4f |P|/(m c) %.4f K_cm %.3e K_rel %.3e | dE %+.3e | kin %+.3e coul %+.3e dd %+.3e cd %+.3e dar %+.3e con %+.3e | P/(m c) %+.3e %+.3e %+.3e\n",
        done,s.time/period,sep/comptonBarrierRadius,beta1,vcm.norm()/c,totalMomentum.norm()/(firstMass*c),kcm/udd,krel/udd,(E-E0)/udd,
        (T.kin-T0.kin)/udd,(T.coul-T0.coul)/udd,(T.dd-T0.dd)/udd,(T.cd-T0.cd)/udd,(T.dar-T0.dar)/udd,(T.con-T0.con)/udd,
        totalMomentum.x/(firstMass*c),totalMomentum.y/(firstMass*c),totalMomentum.z/(firstMass*c));
    }
  }
  {
    const Vec3 Pend=momentum(s.firstVelocity,firstMass)+momentum(s.secondVelocity,secondMass);
    std::printf("KONCOWY P/(m c) %+.3e %+.3e %+.3e  |P| %.3e  po %d krokach  r %.4f r*\n",
      Pend.x/(firstMass*c),Pend.y/(firstMass*c),Pend.z/(firstMass*c),Pend.norm()/(firstMass*c),done,separation(s)/comptonBarrierRadius);
  }
  const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
  std::printf("r %.1f ch %s tilt %.0f tol %.0e N %d depth %d | %s %d/%d | zakres E/|U_dd| %.3e | E_koniec-E0 /|U_dd| %+.3e | r %.4f..%.4f r* | %.0f s\n",
    r/comptonBarrierRadius,ch==0?"para":ch==1?"orto":"bez",atof(argv[3]),tol,perOrbit,acc.maximumDepth,ok?"ok":"AWARIA",done,steps,
    (emax-emin)/udd,(conservativeParticleEnergy(s,1.0)-E0)/udd,rmin/comptonBarrierRadius,rmax/comptonBarrierRadius,wall);
}
