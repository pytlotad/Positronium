// Contact hyperfine energy of a microcanonical Kepler ensemble at the 1S
// energy, with the model's own pair field (audit section 85, plan step 1).
//
// Usage: /tmp/contact [ensemble semi-major axis in a_pair, default 1]
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/contact_ensemble.cpp -o /tmp/contact $(root-config --libs)

// Step 1: contact hyperfine energy of a microcanonical Kepler ensemble at the 1S energy,
// measured with the model's own pair field, para (parallel moments) and ortho (antiparallel).
// Only the dipole-dipole energy depends on position alone; velocity-linear terms average out
// in an isotropic ensemble.  Radial density rho(r) = sqrt(1/r - 1/2a)/N on r < 2a, N = (pi^2/4)(2a)^{5/2}.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <functional>
int main(int argc,char** argv){
  const double fEns=argc>1?std::atof(argv[1]):1.0;
  const double aPair=pairBohrRadius(activePair), a=fEns*pairBohrRadius(activePair), eps=separationFloor(), rs=comptonBarrierRadius, mu=firstMagneticMoment;
  const std::vector<SphereQuadraturePoint> sphere=sphereQuadrature(194);
  // shell average of U for given moments at radius r
  auto shell=[&](double r,const Vec3& m1,const Vec3& m2){ double s=0,w=0;
    for(const auto& q:sphere){ s+=pairDipoleInteractionEnergy(q.direction*r,m1,m2)*q.solidAngleWeight; w+=q.solidAngleWeight; } return s/w; };
  const Vec3 up{0,0,mu}, down{0,0,-mu};
  auto ensembleAverage=[&](const std::function<double(double)>& rho,double rmax,const Vec3& m1,const Vec3& m2){
    const int N=24000; const double lo=std::log(1e-7*std::min(eps,magneticDipoleRadius())), hi=std::log(rmax); double s=0;
    for(int i=0;i<N;++i){ const double u=lo+(hi-lo)*(i+0.5)/N, r=std::exp(u), dr=r*(hi-lo)/N; s+=rho(r)*shell(r,m1,m2)*4*pi*r*r*dr; }
    return s; };
  const double Nmc=(pi*pi/4)*std::pow(2*a,2.5);
  auto micro=[&](double r){ return r<2*a? std::sqrt(std::max(0.0,1/r-1/(2*a)))/Nmc : 0.0; };
  auto quantum=[&](double r){ return std::exp(-2*r/a)/(pi*a*a*a); };
  const double psi0=1/(pi*a*a*a);
  const double psiPair=1/(pi*aPair*aPair*aPair);
  const double eV=eCharge;
  // QED contact part of the splitting: 4/7 of (7/12) alpha^4 m c^2 = (1/3) alpha^4 m c^2
  const double qedContact=fineStructureConstant*fineStructureConstant*fineStructureConstant*fineStructureConstant*electronMass*c*c/3.0;
  // classical-spin analogue for a given density at 0: dE = 2 (2 mu0/3) mu^2 n(0)
  const double classicalSpinQuantum=2.0*(2.0*mu0/3.0)*mu*mu*psi0;
  const double Upm=ensembleAverage(micro,2*a,up,up), Uom=ensembleAverage(micro,2*a,up,down);
  const double Upq=ensembleAverage(quantum,40*a,up,up), Uoq=ensembleAverage(quantum,40*a,up,down);
  // effective contact density recovered from the energies: dE = 2 (2 mu0/3) mu^2 <n>
  const double nMicro=(Uom-Upm)/(2*(2*mu0/3)*mu*mu), nQuant=(Uoq-Upq)/(2*(2*mu0/3)*mu*mu);
  printf("podloga = %.3f r*, promien momentu = %.5f r*  (a_pair/eps_m = %.1f)\n",eps/rs,magneticDipoleRadius()/rs,a/magneticDipoleRadius());
  printf("  zespol mikrokanoniczny: <U>_para %+.4e eV  <U>_orto %+.4e eV  dE = %.4e eV  <n>/|psi0|^2 = %.4f\n",Upm/eV,Uom/eV,(Uom-Upm)/eV,nMicro/psi0);
  printf("  gestosc kwantowa 1S   : <U>_para %+.4e eV  <U>_orto %+.4e eV  dE = %.4e eV  <n>/|psi0|^2 = %.4f\n",Upq/eV,Uoq/eV,(Uoq-Upq)/eV,nQuant/psi0);
  printf("  QED czlon kontaktowy (4/7 z 204.4 GHz): %.4e eV;  klasyczny spin przy |psi0|^2: %.4e eV\n",qedContact/eV,classicalSpinQuantum/eV);
  printf("  <n>_zespol / |psi0(a_pair)|^2 = %.4e   (a = %.3f a_pair)\n",nMicro/psiPair,fEns);
  printf("  stosunek dE_zespol / QED_kontakt = %.3f ;  dE_zespol / klasyczny_spin_przy_psi0 = %.3f\n",(Uom-Upm)/qedContact,(Uom-Upm)/classicalSpinQuantum);
}
