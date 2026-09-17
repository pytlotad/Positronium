// Does the zero-point field drive L low enough to sample the contact term?
// (audit section 87, plan step 3).  Accumulates <n>, <L>, L_min and the radius
// range along an engine trajectory with the Landau-Lifshitz reaction.
//
// Usage: /tmp/zpfcontact <channel 0 para|1 ortho> <orbits> <seed> <ZPF scale>
//                        <tolerance> <start radius in a_pair> <substeps>
//                        <retarded 0/1>
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/zpf_contact_sampling.cpp -o /tmp/zpfcontact $(root-config --libs)

// Step 3, cheaper and longer: does the ZPF drive L low enough to sample the contact term?
// argv: channel orbits seed scale tol f substeps retarded(0/1)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <chrono>
int main(int argc,char** argv){
  const int ch=atoi(argv[1]), orbits=atoi(argv[2]), seed=atoi(argv[3]);
  const double scale=atof(argv[4]), tol=atof(argv[5]), f=atof(argv[6]);
  const int sub=atoi(argv[7]), retarded=atoi(argv[8]);
  const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
  const double aPair=pairBohrRadius(activePair), a0=f*aPair, rs=comptonBarrierRadius;
  const double eps=magneticDipoleRadius(), psi0=1.0/(pi*aPair*aPair*aPair);
  if(scale>0.0) gZeroPointField=makeZeroPointField(0.3,3.0,64,scale,90000u+seed);
  const double v0=std::sqrt(k/(mu*a0)), P0=2*pi*std::sqrt(mu*a0*a0*a0/k);
  State s{}; s.firstPosition={a0*secondMass/M,0,0}; s.secondPosition={-a0*firstMass/M,0,0};
  s.firstVelocity={0,v0*secondMass/M,0}; s.secondVelocity={0,-v0*firstMass/M,0};
  s.firstProperDipole={0,0,firstMagneticMoment};
  s.secondProperDipole={0,0,ch?-secondMagneticMoment:secondMagneticMoment};
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=tol; acc.maximumDepth=14;
  acc.reactionModel=ChargeRadiationReactionModel::individualLandauLifshitz;
  acc.computeOutwardFlux=false; acc.useRetardedExternalForces=retarded!=0;
  ClassicalTrajectoryEngine eng(s,acc);
  const auto t0=std::chrono::steady_clock::now();
  const double bins[5]={0.5,0.2,0.1,0.05,0.02}; double below[5]={0,0,0,0,0};
  double timeTotal=0,contactSum=0,Lsum=0,Lmin=1e300,rmin=1e300,rmax=0;
  int done=0; bool ok=true; const int every=100;
  for(int o=0;o<orbits&&ok;++o){
    const double E=conservativeParticleEnergy(s); const double aNow=E<0?-k/(2*E):a0;
    double P=2*pi*std::sqrt(mu*std::abs(aNow)*aNow*aNow/k); if(!(P>0)||!std::isfinite(P)) P=P0;
    for(int i=0;i<sub;++i){ const double dt=P/sub;
      if(!eng.advance(s,dt)){ ok=false; break; }
      const Vec3 r=s.firstPosition-s.secondPosition, v=s.firstVelocity-s.secondVelocity;
      const double rn=r.norm(), rho2=rn*rn+eps*eps, L=mu*cross(r,v).norm()/hbar;
      contactSum+=3.0*eps*eps/(4.0*pi*std::pow(rho2,2.5))*dt; Lsum+=L*dt;
      for(int b=0;b<5;++b) if(L<bins[b]) below[b]+=dt;
      Lmin=std::min(Lmin,L); rmin=std::min(rmin,rn/rs); rmax=std::max(rmax,rn/rs); timeTotal+=dt; }
    ++done;
    if(done%every==0||!ok){
      printf("%s f %.2f seed %d skala %g %s orbit %d/%d | <n>/|psi0|^2 %.4e <L> %.4f Lmin %.4f | czas przy L< 0.5/0.2/0.1/0.05/0.02: %.3f %.3f %.3f %.3f %.3f | r %.1f-%.1f r* | wall %.0f s\n",
        ch?"orto":"para",f,seed,scale,retarded?"opozn":"chwil",done,orbits,contactSum/timeTotal/psi0,Lsum/timeTotal,Lmin,
        below[0]/timeTotal,below[1]/timeTotal,below[2]/timeTotal,below[3]/timeTotal,below[4]/timeTotal,rmin,rmax,
        std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count());
      std::fflush(stdout); }
  }
  printf("KONIEC %s f %.2f seed %d skala %g %s: %s po %d/%d orbitach, <n>/|psi0|^2 %.4e, Lmin %.4f, rmin %.2f r*\n",
    ch?"orto":"para",f,seed,scale,retarded?"opozn":"chwil",ok?"pelne":"AWARIA",done,orbits,contactSum/timeTotal/psi0,Lmin,rmin);
  gZeroPointField=ZeroPointField{};
}
