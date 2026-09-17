// What breaks the singlet: |w1-w2| with the partner's dipole field removed
// from the assembled local field, every moment left in place (audit section
// 90a).  Zeroing proper moments does NOT remove the field, because
// synchronizeCovariantDipoles refills the lab moment (section 90f).
//
// Usage: /tmp/singletsplit <channel 0 para|1 ortho> <steps>
//                          [r_low r_high beta] (linear scan in a_pair)
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/singlet_field_split.cpp -o /tmp/singletsplit $(root-config --libs)

// Separates the partner's DIPOLE field from everything else WITHOUT zeroing a
// moment: the dipole contribution is recomputed on its own and subtracted, so
// the motional electric part of the same moment stays where it belongs.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
int main(int argc,char** argv){
  const int channel=argc>1?atoi(argv[1]):0;
  const int steps=argc>2?atoi(argv[2]):9;
  const double aPair=hbar*hbar/(pairReducedMass*pairCoulombStrength);
  const double M=firstMass+secondMass;
  printf("kanal %s\n",channel?"orto":"para");
  printf("%-9s %-11s %-11s %-11s %-11s %-11s\n","r/a_pair",
    "|w1-w2| pel","bez dip.pola","|w1| pel","dw/w pel","dw/w bez dip");
  for(int i=0;i<steps;++i){
    const double rLow=argc>3?atof(argv[3]):0.0, rHigh=argc>4?atof(argv[4]):0.0;
    const double betaFixed=argc>5?atof(argv[5]):0.0;
    const double r=rLow>0.0?aPair*(rLow+(rHigh-rLow)*i/std::max(steps-1,1))
                           :aPair*std::pow(4.0,i-(steps-1)/2.0);
    const double v=betaFixed>0.0?betaFixed*c
                   :std::sqrt(pairCoulombStrength/(pairReducedMass*r));
    State s{};
    s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
    s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
    s.firstProperDipole={0,0,firstMagneticMoment};
    s.secondProperDipole={0,0,channel?-secondMagneticMoment:secondMagneticMoment};
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    LocalElectromagneticFields full=localRelativisticFields(s,h);
    // The partner-dipole contribution, recomputed exactly as the assembly adds it.
    const ElectromagneticField dipAtFirst=
        retardedMagneticDipoleField(s.firstPosition,s.time,h,s,false);
    const ElectromagneticField dipAtSecond=
        retardedMagneticDipoleField(s.secondPosition,s.time,h,s,true);
    LocalElectromagneticFields noDipole=full;
    noDipole.atFirst.electric=full.atFirst.electric-dipAtFirst.electric;
    noDipole.atFirst.magnetic=full.atFirst.magnetic-dipAtFirst.magnetic;
    noDipole.atSecond.electric=full.atSecond.electric-dipAtSecond.electric;
    noDipole.atSecond.magnetic=full.atSecond.magnetic-dipAtSecond.magnetic;
    auto pair=[&](const LocalElectromagneticFields& f){
      const Vec3 w1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                      firstGFactor)*(-firstCharge/firstMass);
      const Vec3 w2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                      secondGFactor)*(-secondCharge/secondMass);
      return std::pair<Vec3,Vec3>{w1,w2};
    };
    const auto F=pair(full); const auto N=pair(noDipole);
    printf("%-9.4g %-11.4e %-11.4e %-11.4e %-11.6f %-11.6f\n",
      r/aPair,(F.first-F.second).norm(),(N.first-N.second).norm(),
      F.first.norm(),
      (F.first-F.second).norm()/std::max(F.first.norm(),1e-300),
      (N.first-N.second).norm()/std::max(N.first.norm(),1e-300));
  }
}
