// Which part of the retarded dipole sector breaks action and reaction (audit
// section 105).  Follows the section 104 engine trajectory (para, 1 r*, tilt 0,
// one orbit) and at every step evaluates F1+F2 on the SAME state and engine
// history with probe switches, splitting the net force into six contributions
// (charge or moment, in the partner's charge, pole or magnetization field)
// plus the material-derivative term, integrated over the orbit.
//
// Needs the probe switches, so build it with:
//   tools/build_variant.sh tools/dipole_sector_momentum_split.cpp /tmp/split --probe-switches
// Which part of the retarded dipole sector breaks action and reaction?
// Engine trajectory as in section 104 (para, 1 r*, tilt 0, one orbit); at every
// step the net force F1+F2 is evaluated on the SAME state and history with
// probe switches, and each contribution is integrated over the orbit.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
int main(){
  const double r=comptonBarrierRadius; const double k=pairCoulombStrength,mu=pairReducedMass,M=firstMass+secondMass;
  const double v=std::sqrt(k/(mu*r)), period=2*pi*r/v;
  State s{};
  s.firstPosition={r*secondMass/M,0,0}; s.secondPosition={-r*firstMass/M,0,0};
  s.firstVelocity={0,v*secondMass/M,0}; s.secondVelocity={0,-v*firstMass/M,0};
  s.firstProperDipole={0,0,firstMagneticMoment}; s.secondProperDipole={0,0,secondMagneticMoment};
  synchronizeCovariantDipoles(s);
  ClassicalTrajectoryEngine::Accuracy acc; acc.relativeTolerance=1e-8; acc.maximumDepth=20;
  acc.reactionModel=ChargeRadiationReactionModel::disabled; acc.computeOutwardFlux=false; acc.useRetardedExternalForces=true;
  ClassicalTrajectoryEngine eng(s,acc);
  const std::vector<std::string> all={"CREM_PROBE_NO_POLES","CREM_PROBE_NO_MAGNETIZATION","CREM_PROBE_NO_LW_IN_MOMENT",
    "CREM_PROBE_NO_MATERIAL","CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_MOMENT_FORCE"};
  auto net=[&](const std::vector<std::string>& on){
    for(auto& n:all) unsetenv(n.c_str());
    for(auto& n:on) setenv(n.c_str(),"1",1);
    const MutualForces F=retardedExternalForces(s,eng.history());
    for(auto& n:all) unsetenv(n.c_str());
    return F.first+F.second;
  };
  struct Combo{const char* name; std::vector<std::string> on;};
  const std::vector<Combo> combos={
    {"T",{}},
    {"C",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_MOMENT_FORCE"}},
    {"P1_poles",{"CREM_PROBE_NO_MOMENT_FORCE","CREM_PROBE_NO_MAGNETIZATION"}},
    {"P1_mag",{"CREM_PROBE_NO_MOMENT_FORCE","CREM_PROBE_NO_POLES"}},
    {"P2",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE"}},
    {"P2_nomat",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_MATERIAL"}},
    {"P2_lw",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_POLES","CREM_PROBE_NO_MAGNETIZATION"}},
    {"P2_poles",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_LW_IN_MOMENT","CREM_PROBE_NO_MAGNETIZATION"}},
    {"P2_mag",{"CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_LW_IN_MOMENT","CREM_PROBE_NO_POLES"}}};
  std::vector<Vec3> integ(combos.size());
  const int steps=256; const double dt=period/steps; bool ok=true;
  for(int i=0;i<steps;++i){
    if(!eng.advance(s,dt)){ ok=false; break; }
    for(size_t j=0;j<combos.size();++j) integ[j]+=net(combos[j].on)*dt;
    if((i+1)%64==0){
      const Vec3 P=momentum(s.firstVelocity,firstMass)+momentum(s.secondVelocity,secondMass);
      std::printf("t/P %.2f  r %.3f r*  |P|/mc %.4f  |integral T|/mc %.4f\n",(i+1.0)/steps,separation(s)/r,
        P.norm()/(firstMass*c),integ[0].norm()/(firstMass*c));
    }
  }
  const Vec3 P=momentum(s.firstVelocity,firstMass)+momentum(s.secondVelocity,secondMass);
  const Vec3 dir=P*(1.0/std::max(P.norm(),1e-300));
  auto I=[&](const char* n){ for(size_t j=0;j<combos.size();++j) if(std::string(combos[j].name)==n) return integ[j]; return Vec3{}; };
  const Vec3 C=I("C");
  struct Part{const char* label; Vec3 d;};
  const std::vector<Part> parts={
    {"1 ladunek w polu ladunku partnera (LW)",C},
    {"2 ladunek w polu biegunow partnera",I("P1_poles")-C},
    {"3 ladunek w polu magnetyzacji partnera",I("P1_mag")-C},
    {"4 moment w polu ladunku partnera",I("P2_lw")-C},
    {"5 moment w polu biegunow partnera",I("P2_poles")-C},
    {"6 moment w polu magnetyzacji partnera",I("P2_mag")-C},
    {"   w tym czlon materialny (wszystkie pola)",I("P2")-I("P2_nomat")}};
  const double mc=firstMass*c; const double total=dot(I("T"),dir);
  std::printf("\n%s po obiegu; |P| mechaniczny %.4f mc, calka F1+F2 wzdluz P %.4f mc\n",ok?"ok":"AWARIA",P.norm()/mc,total/mc);
  std::printf("%-44s %12s %10s %12s\n","wklad","wzdluz P/mc","udzial","|wektor|/mc");
  Vec3 sum;
  for(int j=0;j<6;++j) sum+=parts[j].d;
  for(auto& p:parts) std::printf("%-44s %+12.4f %9.1f%% %12.4f\n",p.label,dot(p.d,dir)/mc,100*dot(p.d,dir)/total,p.d.norm()/mc);
  std::printf("suma 1-6 wzdluz P %+.4f mc wobec pelnej %+.4f mc (domkniecie %.2e)\n",dot(sum,dir)/mc,total/mc,
    (sum-I("T")).norm()/std::max(I("T").norm(),1e-300));
}
