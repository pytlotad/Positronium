// One-off (audit 120): how big is the Thomas back-reaction inside the
// retarded sum, as a force and as a net force on the pair?
#include "modules/crem_trajectory.hpp"
#include <cstdio>
int main(){
  const double r=comptonBarrierRadius,k=pairCoulombStrength,mu=pairReducedMass;
  const double M=firstMass+secondMass,v=std::sqrt(k/(mu*r)),period=2*pi*r/v;
  State s{};
  s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
  s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
  s.firstProperDipole={0,0,firstMagneticMoment};
  s.secondProperDipole={0,0,secondMagneticMoment};
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine::Accuracy a; a.relativeTolerance=1e-8; a.maximumDepth=20;
  a.reactionModel=ChargeRadiationReactionModel::disabled; a.computeOutwardFlux=false;
  a.useRetardedExternalForces=true;
  ClassicalTrajectoryEngine eng(s,a);
  const int steps=256; const double dt=period/steps;
  double worstNet=0,worstForce=0; Vec3 thomasIntegral,totalIntegral;
  for(int i=0;i<steps;++i){
    if(!eng.advance(s,dt)) { std::printf("AWARIA krok %d\n",i); break; }
    const MutualForces total=retardedExternalForces(s,eng.history());
    const MutualForces thomas=thomasBackReactionForces(s,
        thomasBmtDipoleDerivatives(s,eng.history()));
    const Vec3 netTotal=total.first+total.second;
    const Vec3 netThomas=thomas.first+thomas.second;
    thomasIntegral+=netThomas*dt; totalIntegral+=netTotal*dt;
    worstNet=std::max(worstNet,netThomas.norm()/std::max(netTotal.norm(),1e-300));
    worstForce=std::max(worstForce,thomas.first.norm()/std::max(total.first.norm(),1e-300));
    if((i+1)%64==0)
      std::printf("krok %3d r %.4f r*  |F_thomas|/|F_total| %.4e  "
                  "|net Thomas|/|net total| %.4e\n",i+1,
                  separation(s)/comptonBarrierRadius,
                  thomas.first.norm()/std::max(total.first.norm(),1e-300),
                  netThomas.norm()/std::max(netTotal.norm(),1e-300));
  }
  const double mc=firstMass*c;
  std::printf("\nnajwiekszy |F_thomas|/|F_total| %.4e; najwiekszy udzial w "
              "wypadkowej %.4e\ncalka wypadkowej: Thomas %.6f mc, calosc "
              "%.6f mc, stosunek %.4f\n",worstForce,worstNet,
              thomasIntegral.norm()/mc,totalIntegral.norm()/mc,
              thomasIntegral.norm()/std::max(totalIntegral.norm(),1e-300));
}
