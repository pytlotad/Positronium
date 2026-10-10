#include "modules/crem_collapse.hpp"
#include <cstdio>
// Audit 404: offline replay of a failing engine step dumped by CREM_DUMP_FAIL.
// argv: dump file, mode ("sweep" | "scan" [points]).
// sweep: the LIVE step-doubling sweep (16 dtFail ... dtFail/32) with |grad| at the coarse endpoint.
// scan: coarse(tau) = one step of length tau from the start, tau log-spaced in [dt/64, 16 dt]; per point the
//       gradient force on each particle, the pinned history segment and the pole-cancellation ratio of each probe.
static State start; static StateHistory hist; static double dtFail=0.0;
static const auto model=ChargeRadiationReactionModel::stochasticElectricDipole;
static void load(const char* path){
  std::FILE* f=std::fopen(path,"rb"); if(!f){std::perror(path); std::exit(1);}
  std::uint64_t nodes=0;
  if(std::fread(&dtFail,sizeof dtFail,1,f)!=1||std::fread(&start,sizeof(State),1,f)!=1
     ||std::fread(&nodes,sizeof nodes,1,f)!=1){std::fprintf(stderr,"short dump\n");std::exit(1);}
  for(std::uint64_t i=0;i<nodes;++i){State s; if(std::fread(&s,sizeof s,1,f)!=1) std::exit(1); hist.push_back(s);}
  std::fclose(f);
}
static State coarseAfter(double tau){ State s=start; integrateElectrodynamicStep(s,tau,hist,false,model,true); return s; }
static State fineAfter(double tau){ State s=start; StateHistory h=hist; integrateElectrodynamicStep(s,0.5*tau,h,false,model,true);
  integrateElectrodynamicStep(s,0.5*tau,h,false,model,true); return s; }
// Per-probe diagnostics of covariantDipoleGradientForce for target `first` at state s (same pin, same probes).
static void probes(const State& s,bool targetIsFirst,double ratio[6],std::size_t& segment){
  const Vec3 p=targetIsFirst?s.firstPosition:s.secondPosition;
  const double step=std::max(1.0e-4*separation(s),1.0e-3*nuclearCutoff);
  const RetardedSegmentPin pin=retardedSegmentPinAt(hist,s,!targetIsFirst,p,s.time,8.0*step/c);
  segment=pin.newerIndex;
  const RetardedSegmentPinGuard guard(pin);
  for(int k=0;k<6;++k){ Vec3 o; double d=(k%2?-step:step); if(k/2==0) o.x=d; else if(k/2==1) o.y=d; else o.z=d;
    (void)fieldFromOtherParticleAt(p+o,s.time,s,hist,targetIsFirst,1.0e-5,matchedMomentSoftening());
    ratio[k]=gPoleCancellationRatio; }
}
int main(int argc,char**argv){
  if(argc<3){std::fprintf(stderr,"usage: h0404 dump sweep|scan [points]\n");return 1;}
  load(argv[1]);
  std::printf("dump: dt=%.6e t=%.9e r=%.6e history=%zu nodes (%.9e .. %.9e)\n",dtFail,start.time,separation(start),
              hist.size(),hist.empty()?0.0:hist.front().time,hist.empty()?0.0:hist.back().time);
  if(std::string(argv[2])=="sweep"){
    for(int h=0;h<10;++h){ const double tau=dtFail*std::pow(2.0,4.0-h);
      const State a=coarseAfter(tau), b=fineAfter(tau);
      const double err=std::max({(a.firstPosition-b.firstPosition).norm()/std::max(separation(b),collisionBoundaryRadius),
        (a.secondPosition-b.secondPosition).norm()/std::max(separation(b),collisionBoundaryRadius),
        (a.firstVelocity-b.firstVelocity).norm()/std::max((b.firstVelocity-b.secondVelocity).norm(),1e-6*c),
        (a.secondVelocity-b.secondVelocity).norm()/std::max((b.firstVelocity-b.secondVelocity).norm(),1e-6*c)});
      std::printf("tau=dt*2^%+d err=%.6e |grad1|=%.6e |grad2|=%.6e\n",4-h,err,
        covariantDipoleGradientForce(a,hist,true).norm(),covariantDipoleGradientForce(a,hist,false).norm()); }
  } else if(std::string(argv[2])=="split"){
    // Which part of the state drives the change of grad1 across tau: take the
    // coarse endpoint at tau and copy ONE group of fields into the start state.
    const int points=argc>3?std::atoi(argv[3]):24;
    for(int i=0;i<points;++i){ const double tau=dtFail*0.25*std::pow(16.0,double(i)/(points-1));
      const State a=coarseAfter(tau);
      State onlyTime=start; onlyTime.time=a.time;
      State onlyPos=start; onlyPos.firstPosition=a.firstPosition; onlyPos.secondPosition=a.secondPosition;
      State onlyVel=start; onlyVel.firstVelocity=a.firstVelocity; onlyVel.secondVelocity=a.secondVelocity;
      onlyVel.firstAcceleration=a.firstAcceleration; onlyVel.secondAcceleration=a.secondAcceleration;
      State onlyDip=start; onlyDip.firstDipole=a.firstDipole; onlyDip.secondDipole=a.secondDipole;
      onlyDip.firstElectricDipole=a.firstElectricDipole; onlyDip.secondElectricDipole=a.secondElectricDipole;
      onlyDip.firstProperDipole=a.firstProperDipole; onlyDip.secondProperDipole=a.secondProperDipole;
      State posTime=onlyPos; posTime.time=a.time;
      std::printf("tau/dt=%.4f full=%.6e time=%.6e pos=%.6e vel=%.6e dip=%.6e pos+time=%.6e start=%.6e\n",tau/dtFail,
        covariantDipoleGradientForce(a,hist,true).norm(),covariantDipoleGradientForce(onlyTime,hist,true).norm(),
        covariantDipoleGradientForce(onlyPos,hist,true).norm(),covariantDipoleGradientForce(onlyVel,hist,true).norm(),
        covariantDipoleGradientForce(onlyDip,hist,true).norm(),covariantDipoleGradientForce(posTime,hist,true).norm(),
        covariantDipoleGradientForce(start,hist,true).norm()); }
  } else if(std::string(argv[2])=="timesplit"){
    // Start state with only the time advanced: gradient of mu.B + p.E split into the
    // Lienard-Wiechert (charge) field and the retarded magnetic-dipole field, same pin.
    const int points=argc>3?std::atoi(argv[3]):24;
    const Vec3 m=start.firstDipole, pe=start.firstElectricDipole, p0=start.firstPosition;
    const double step=std::max(1.0e-4*separation(start),1.0e-3*nuclearCutoff);
    for(int i=0;i<points;++i){ const double tau=dtFail*(0.2+0.8*double(i)/(points-1));
      State s=start; s.time=start.time+tau;
      const RetardedSegmentPin pin=retardedSegmentPinAt(hist,s,false,p0,s.time,8.0*step/c);
      const RetardedSegmentPinGuard guard(pin);
      Vec3 gLW, gDip;
      for(int ax=0;ax<3;++ax){ Vec3 o; if(ax==0) o.x=step; else if(ax==1) o.y=step; else o.z=step;
        const auto lw=[&](const Vec3& x){ const ElectromagneticField f=lienardWiechertField(x,s.time,hist,s,false,secondCharge,0.0,matchedMomentSoftening());
          return dot(m,f.magnetic)+dot(pe,f.electric); };
        const auto dip=[&](const Vec3& x){ const ElectromagneticField f=retardedMagneticDipoleField(x,s.time,hist,s,false,1.0e-5);
          return dot(m,f.magnetic)+dot(pe,f.electric); };
        const double dl=(lw(p0+o)-lw(p0-o))/(2*step), dd=(dip(p0+o)-dip(p0-o))/(2*step);
        if(ax==0){gLW.x=dl;gDip.x=dd;} else if(ax==1){gLW.y=dl;gDip.y=dd;} else {gLW.z=dl;gDip.z=dd;} }
      const ElectromagneticField fd=retardedMagneticDipoleField(p0,s.time,hist,s,false,1.0e-5);
      std::printf("tau/dt=%.4f seg=%zu |gLW|=%.6e gDip=(%.6e %.6e %.6e) |gDip|=%.6e |Bdip|=%.9e |Edip|=%.9e ratio=%.3g full=%.6e\n",
        tau/dtFail,pin.newerIndex,gLW.norm(),gDip.x,gDip.y,gDip.z,gDip.norm(),fd.magnetic.norm(),fd.electric.norm(),gPoleCancellationRatio,
        covariantDipoleGradientForce(s,hist,true).norm()); }
  } else if(std::string(argv[2])=="nodes"){
    // Retarded time of the target's position (source = second particle) for tau in [0, dt], against node times.
    const Vec3 p0=start.firstPosition;
    const auto srcAt=[&](double t){ for(std::size_t k=1;k<hist.size();++k) if(hist[k].time>=t){
        const double f=(t-hist[k-1].time)/(hist[k].time-hist[k-1].time);
        return hist[k-1].secondPosition*(1-f)+hist[k].secondPosition*f; } return hist.back().secondPosition; };
    for(double frac: {0.0,0.35,0.55,0.75,1.0}){ const double tobs=start.time+frac*dtFail; double tr=tobs-separation(start)/c;
      for(int it=0;it<50;++it) tr=tobs-(p0-srcAt(tr)).norm()/c;
      std::size_t k=0; while(k<hist.size()&&hist[k].time<tr) ++k;
      std::printf("tau/dt=%.2f t_ret=%.12e between nodes %zu (%.12e) and %zu (%.12e): from lower %.4f of spacing; step/c=%.3e s\n",frac,tr,k-1,
        hist[k-1].time,k,hist[k].time,(tr-hist[k-1].time)/(hist[k].time-hist[k-1].time),1.0e-4*separation(start)/c); }
    for(std::size_t k=95;k<101&&k<hist.size();++k){ const State& n=hist[k];
      std::printf("node %zu t=%.12e mu2=(%.9e %.9e %.9e) a2=(%.6e %.6e %.6e)\n",k,n.time,n.secondDipole.x,n.secondDipole.y,n.secondDipole.z,
        n.secondAcceleration.x,n.secondAcceleration.y,n.secondAcceleration.z); }
  } else if(std::string(argv[2])=="probesplit"){
    // y-gradient of the dipole-field coupling, split into the two-pole construction and the magnetization term.
    const int points=argc>3?std::atoi(argv[3]):17;
    const Vec3 m=start.firstDipole, pe=start.firstElectricDipole, p0=start.firstPosition;
    const double step=std::max(1.0e-4*separation(start),1.0e-3*nuclearCutoff);
    for(int i=0;i<points;++i){ const double tau=dtFail*(0.2+0.8*double(i)/(points-1));
      State s=start; s.time=start.time+tau;
      const RetardedSegmentPinGuard guard(retardedSegmentPinAt(hist,s,false,p0,s.time,8.0*step/c));
      const auto dual=[&](const Vec3& x){ const ElectromagneticField d=twoChargeLimitDipoleField(x,s.time,hist,s,false,1.0e-5,magneticDipoleRadius(),
          [&](double time,Vec3& mo,Vec3& f1,Vec3& f2){ const RetardedElectricDipoleKinematics k=historicalIntegratedDipoleKinematics(hist,s,false,time,false);
            mo=k.moment/(c*c); f1=k.firstDerivative/(c*c); f2=k.secondDerivative/(c*c); });
        return dot(m,d.magnetic*(-c*c))+dot(pe,d.electric); };
      const auto full=[&](const Vec3& x){ const ElectromagneticField f=retardedMagneticDipoleFieldExact(x,s.time,hist,s,false,1.0e-5);
        return dot(m,f.magnetic)+dot(pe,f.electric); };
      Vec3 o; o.y=step;
      const double gd=(dual(p0+o)-dual(p0-o))/(2*step), gf=(full(p0+o)-full(p0-o))/(2*step);
      std::printf("tau/dt=%.4f gy_dual=%.6e gy_magnetization=%.6e gy_full=%.6e\n",tau/dtFail,gd,gf-gd,gf); }
  } else if(std::string(argv[2])=="fields"){
    // Field magnitudes of the two-pole part and the magnetization part at the target, plus the moment kinematics.
    State s=start;
    const Vec3 p0=s.firstPosition;
    const ElectromagneticField d=twoChargeLimitDipoleField(p0,s.time,hist,s,false,1.0e-5,magneticDipoleRadius(),
          [&](double time,Vec3& mo,Vec3& f1,Vec3& f2){ const RetardedElectricDipoleKinematics k=historicalIntegratedDipoleKinematics(hist,s,false,time,false);
            mo=k.moment/(c*c); f1=k.firstDerivative/(c*c); f2=k.secondDerivative/(c*c); });
    const ElectromagneticField f=retardedMagneticDipoleFieldExact(p0,s.time,hist,s,false,1.0e-5);
    const Vec3 dualB=d.magnetic*(-c*c), dualE=d.electric;
    std::printf("dual |B|=%.6e |E|=%.6e; magnetization |B|=%.6e |E|=%.6e; full |B|=%.6e |E|=%.6e\n",dualB.norm(),dualE.norm(),
      (f.magnetic-dualB).norm(),(f.electric-dualE).norm(),f.magnetic.norm(),f.electric.norm());
    const double tr=s.time-separation(s)/c;
    const RetardedElectricDipoleKinematics k=historicalIntegratedDipoleKinematics(hist,s,false,tr,false);
    std::printf("source moment at t_ret: |m|=%.6e |m'|=%.6e (omega=%.3e 1/s) |m''|=%.6e; point-dipole |B| at r: %.6e\n",k.moment.norm(),
      k.firstDerivative.norm(),k.firstDerivative.norm()/k.moment.norm(),k.secondDerivative.norm(),
      mu0/(4*pi)*2*k.moment.norm()/std::pow(separation(s),3));
    std::printf("target: |mu1|=%.6e |p1|=%.6e (c|p1|/|mu1|=%.3e) v1/c=%.4f v2/c=%.4f\n",s.firstDipole.norm(),s.firstElectricDipole.norm(),
      c*s.firstElectricDipole.norm()/s.firstDipole.norm(),s.firstVelocity.norm()/c,s.secondVelocity.norm()/c);
  } else if(std::string(argv[2])=="mdd"){
    // |m''| and |m'| of the source moment from the history stencil across nodes 95..100, against omega^2 |m|.
    for(int i=0;i<=50;++i){ const double t=hist[95].time+(hist[100].time-hist[95].time)*i/50.0;
      const RetardedElectricDipoleKinematics k=historicalIntegratedDipoleKinematics(hist,start,false,t,false);
      const double w=k.firstDerivative.norm()/k.moment.norm();
      std::printf("t-node95=%.3e (node spacing %.3e) |m'|/|m|=%.4e |m''|=%.4e omega^2|m|=%.4e ratio=%.1f\n",t-hist[95].time,
        hist[96].time-hist[95].time,w,k.secondDerivative.norm(),w*w*k.moment.norm(),k.secondDerivative.norm()/(w*w*k.moment.norm())); }
  } else if(std::string(argv[2])=="bmt"){
    // Thomas-BMT rate at the start state vs the node-to-node change of the stored moments (last 12 nodes).
    const DipoleDerivatives d=thomasBmtDipoleDerivatives(start,hist);
    std::printf("BMT |dmu1/dt|/|mu1|=%.4e  |dmu2/dt|/|mu2|=%.4e 1/s\n",d.first.norm()/start.firstDipole.norm(),d.second.norm()/start.secondDipole.norm());
    for(std::size_t k=hist.size()-12;k+1<hist.size();++k){ const State& a=hist[k]; const State& b=hist[k+1]; const double span=b.time-a.time;
      std::printf("node %zu->%zu span=%.3e  |dmu1|/|mu1|/span=%.4e |dmu2|/|mu2|/span=%.4e  |dmu1_proper|/span=%.4e r=%.4e\n",k,k+1,span,
        (b.firstDipole-a.firstDipole).norm()/a.firstDipole.norm()/span,(b.secondDipole-a.secondDipole).norm()/a.secondDipole.norm()/span,
        (b.firstProperDipole-a.firstProperDipole).norm()/a.firstProperDipole.norm()/span,separation(b)); }
  } else if(std::string(argv[2])=="torque"){
    // M1 radiation-reaction torque rate vs the Thomas-BMT rate at the start state.
    const MutualForces forces=retardedExternalForces(start,hist);
    const ParticleMultipoleRadiation rad=particleMultipoleRadiation(start,forces,hist,false,model,true);
    const DipoleDerivatives d=thomasBmtDipoleDerivatives(start,hist);
    std::printf("BMT rate: %.4e %.4e 1/s; M1 reaction torque rate gamma_g|T|/|mu|: %.4e %.4e 1/s\n",
      d.first.norm()/start.firstProperDipole.norm(),d.second.norm()/start.secondProperDipole.norm(),
      std::abs(firstGyromagneticRatioOf())*rad.firstDipoleTorque.norm()/start.firstProperDipole.norm(),
      std::abs(secondGyromagneticRatioOf())*rad.secondDipoleTorque.norm()/start.secondProperDipole.norm());
  } else {
    const int points=argc>3?std::atoi(argv[3]):400;
    for(int i=0;i<points;++i){ const double tau=dtFail/64.0*std::pow(1024.0,double(i)/(points-1));
      const State a=coarseAfter(tau);
      for(int t=0;t<2;++t){ double ratio[6]; std::size_t seg=0; probes(a,t==0,ratio,seg);
        const Vec3 g=covariantDipoleGradientForce(a,hist,t==0);
        std::printf("%s tau/dt=%.6e t=%.12e g=(%.6e %.6e %.6e) |g|=%.6e seg=%zu ratio=",t==0?"T1":"T2",tau/dtFail,a.time,
                    g.x,g.y,g.z,g.norm(),seg);
        for(int k=0;k<6;++k) std::printf("%.3g%s",ratio[k],k<5?",":"\n"); } }
  }
}
