// Dipole force over Coulomb force on the instantaneous path, circular orbits
// with moments along the normal, at the radii of section 75 (audit 93b).
// Compare CREM_NO_THOMAS_BACKREACTION=1 and CREM_MAGNETIC_RADIUS_SCALE=0.05.
//
// Usage: /tmp/dipforce
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/dipole_force_ratio.cpp -o /tmp/dipforce $(root-config --libs)

// Does the section 86 change of the moment softening radius (0.05 r* -> 0.96682 r*)
// reach the ZPF balance runs at all?  Those runs use 74b's configuration:
// circular orbit, moments along the normal, instantaneous forces.  At the radii
// section 75 measured (1, 0.3, 0.1, 0.03 a_pair) print the dipole force against
// the Coulomb force, with the old and the new softening.
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cstdlib>
int main(){
  const double aPair=hbar*hbar/(pairReducedMass*pairCoulombStrength);
  const double M=firstMass+secondMass;
  printf("%-12s %-12s %-14s %-14s %-12s\n",
    "r/a_pair","r/r*","|F_dip|/|F_C|","|F_dip| [N]","|F_C| [N]");
  for(double f : {1.0,0.3,0.1,0.03,0.01}){
    const double r=f*aPair;
    const double v=std::sqrt(pairCoulombStrength/(pairReducedMass*r));
    State s{};
    s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
    s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
    // 74b: moments along the orbital normal, which is z for this geometry.
    s.firstProperDipole={0,0,firstMagneticMoment};
    s.secondProperDipole={0,0,secondMagneticMoment};
    synchronizeCovariantDipoles(s);
    const MutualForces all=allExternalForces(s);
    State bare=s; bare.firstProperDipole={}; bare.secondProperDipole={};
    State bareFresh{};
    bareFresh.firstPosition=s.firstPosition; bareFresh.secondPosition=s.secondPosition;
    bareFresh.firstVelocity=s.firstVelocity; bareFresh.secondVelocity=s.secondVelocity;
    synchronizeCovariantDipoles(bareFresh);
    const MutualForces none=allExternalForces(bareFresh);
    const double fdip=(all.first-none.first).norm();
    const double fc=pairCoulombStrength/(r*r);
    printf("%-12.4g %-12.4g %-14.6e %-14.6e %-12.6e\n",
      f,r/comptonBarrierRadius,fdip/fc,fdip,fc);
  }
  return 0;
}
