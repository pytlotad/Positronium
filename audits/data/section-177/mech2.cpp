#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <initializer_list>
int main(int argc,char** argv){
    const double b=(argc>1?atof(argv[1]):0.01);
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double period=2.0*pi*r0/speed;
    const double tilt=40.0*pi/180;
    const Vec3 dir{std::sin(tilt),0.0,std::cos(tilt)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    const Vec3 com{b*c,0,0};
    s.firstVelocity=Vec3{0,speed*secondMass/total,0}+com;
    s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+com;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-9; acc.maximumDepth=20;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    const double scale=firstMagneticMoment+secondMagneticMoment;
    std::printf("beta_CM = %.4g\n",b);
    std::printf("%8s %13s %13s %11s %12s %12s\n","t/T","w1","w2",
                "|dw|/w1","|m|/mu_sum","|dm/dt|");
    double worstOm=0,worstM=0,worstRate=0;
    const int N=500;
    for(int i=0;i<=N;++i){
        if(i&&!engine.advance(s,period/N)) break;
        const DipoleDerivatives d=thomasBmtDipoleDerivatives(s,engine.history());
        const double w1=d.first.norm()/std::max(s.firstDipole.norm(),1e-300);
        const double w2=d.second.norm()/std::max(s.secondDipole.norm(),1e-300);
        const double om=std::abs(w1-w2)/std::max(w1,1e-300);
        const double m=(s.firstDipole+s.secondDipole).norm()/scale;
        const double rate=(d.first+d.second).norm()
            /std::max(d.first.norm()+d.second.norm(),1e-300);
        worstOm=std::max(worstOm,om); worstM=std::max(worstM,m);
        worstRate=std::max(worstRate,rate);
        if(i%100==0)
            std::printf("%8.2f %13.5e %13.5e %11.4e %12.4e %12.4e\n",
                        double(i)/N,w1,w2,om,m,rate);
    }
    std::printf("  worst: |dw|/w1 %.4e   |m| %.4e   |dm/dt| %.4e\n",
                worstOm,worstM,worstRate);
}
