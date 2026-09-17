// Work of the dipole sector over one Kepler orbit against the Larmor energy
// of that orbit (audit section 93c).  Small partly by construction: the
// conservative part is a gradient and does no work around a closed orbit.
//
// Usage: /tmp/dipwork [eccentricity, default 0]
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/dipole_work_per_orbit.cpp -o /tmp/dipwork $(root-config --libs)

// Force is not work.  The ZPF balance is a ratio of WORK per orbit, so the
// question is whether the Thomas back-reaction (added after section 75)
// changes the dipole sector's work per orbit enough to matter against the
// pumping share's uncertainty of +-0.152.
// Integrates F_dip . v over one circular orbit at the radii section 75 used,
// and compares it with the radiated (Larmor) energy over the same orbit.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(int argc,char** argv){
  // A CIRCULAR orbit is the one geometry where the dipole force is very
  // nearly perpendicular to the velocity, so it does almost no work by
  // construction.  The ZPF runs start circular but the field makes
  // excursions of 2-4.6x the start radius (75c), so the eccentric case is
  // the one that decides whether the sector can matter.
  const double ecc=argc>1?atof(argv[1]):0.0;
  const double aPair=hbar*hbar/(pairReducedMass*pairCoulombStrength);
  const double M=firstMass+secondMass;
  const int nodes=720;
  printf("%-10s %-14s %-14s %-14s %-12s\n",
    "r/a_pair","W_dip [eV]","W_Larmor [eV]","|W_dip|/Larmor","uwaga");
  for(double f : {1.0,0.3,0.1,0.03}){
    const double a=f*aPair;
    const double r0=a;   // reference radius for the Larmor comparison
    const double period=2*pi*std::sqrt(pairReducedMass*a*a*a
                                       /pairCoulombStrength);
    const double meanMotion=2*pi/period;
    double work=0.0;
    for(int i=0;i<nodes;++i){
      // Kepler node by eccentric anomaly, so the eccentric case is sampled
      // in time correctly rather than uniformly in angle.
      const double E=2*pi*(i+0.5)/nodes;
      const double r=a*(1.0-ecc*std::cos(E));
      const double x=a*(std::cos(E)-ecc);
      const double y=a*std::sqrt(1.0-ecc*ecc)*std::sin(E);
      const double vx=-a*meanMotion*std::sin(E)/(1.0-ecc*std::cos(E));
      const double vy=a*meanMotion*std::sqrt(1.0-ecc*ecc)*std::cos(E)
                      /(1.0-ecc*std::cos(E));
      State s{};
      const Vec3 rel{x,y,0};
      const Vec3 relV{vx,vy,0};
      (void)r;
      s.firstPosition=rel*(secondMass/M); s.secondPosition=rel*(-firstMass/M);
      s.firstVelocity=relV*(secondMass/M); s.secondVelocity=relV*(-firstMass/M);
      s.firstProperDipole={0,0,firstMagneticMoment};
      s.secondProperDipole={0,0,secondMagneticMoment};
      synchronizeCovariantDipoles(s);
      const MutualForces withM=allExternalForces(s);
      State bare{};
      bare.firstPosition=s.firstPosition; bare.secondPosition=s.secondPosition;
      bare.firstVelocity=s.firstVelocity; bare.secondVelocity=s.secondVelocity;
      synchronizeCovariantDipoles(bare);
      const MutualForces none=allExternalForces(bare);
      work+=(dot(withM.first-none.first,s.firstVelocity)
            +dot(withM.second-none.second,s.secondVelocity))*(period/nodes);
    }
    // Larmor over one orbit for the relative motion.
    const double accel=pairCoulombStrength/(pairReducedMass*r0*r0);
    const double larmor=eCharge*eCharge*accel*accel
      /(6.0*pi*epsilon0*c*c*c)*period;
    printf("e=%.1f  %-8.3g %-14.6e %-14.6e %-14.6e %s\n",ecc,
      f,work/eCharge,larmor/eCharge,std::abs(work)/larmor,
      f==0.1?"<- bilans 75c":"");
  }
  return 0;
}
