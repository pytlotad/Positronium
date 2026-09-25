// 178f item 1: the mirror is the INSTANTANEOUS separation-normal plane, which
// turns with the orbit.  A drift fixed along x therefore lies in it only when
// the separation is along +-x, i.e. at t/T = 0 and 0.5.  Predicted: the rate
// equality is restored exactly there and nowhere else.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
int main(){
    const double total=firstMass+secondMass;
    const double r0=pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double period=2.0*pi*r0/speed;
    const double th=40.0*pi/180;
    const Vec3 dir{std::sin(th),0.0,std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    const Vec3 drift{0.01*c,0,0};
    s.firstVelocity=Vec3{0,speed*secondMass/total,0}+drift;
    s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+drift;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-10; acc.maximumDepth=20;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    std::printf("drift 0.01c along x.  sepAngle is the angle of the "
                "separation\nfrom x; the drift lies in the mirror when it is "
                "0 or 180.\n\n");
    std::printf("%8s %11s %13s %14s\n","t/T","sepAngle","|w1-w2|/w1",
                "B2 vs -M_y(B1)");
    const int N=1000;
    for(int i=0;i<=N;++i){
        if(i&&!engine.advance(s,period/N)) break;
        const double frac=double(i)/N;
        const bool interesting=(i%125==0)||std::abs(frac-0.5)<0.006;
        if(!interesting) continue;
        const Vec3 sep=s.firstPosition-s.secondPosition;
        const double ang=std::atan2(sep.y,sep.x)*180.0/pi;
        const LocalElectromagneticFields f=
            localRelativisticFields(s,engine.history());
        const Vec3 b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,
                                              firstGFactor);
        const Vec3 b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                              secondGFactor);
        const double w1=b1.norm()*std::abs(firstCharge/firstMass);
        const double w2=b2.norm()*std::abs(secondCharge/secondMass);
        // mirror in the plane of the CURRENT separation and the orbit normal
        const Vec3 nhat=sep/sep.norm();
        const Vec3 zhat{0,0,1};
        const Vec3 perp=cross(zhat,nhat);           // normal to that plane
        const auto mirror=[&](const Vec3& v){
            return v-perp*(2.0*dot(v,perp)); };
        const Vec3 target=mirror(b1)*-1.0;
        std::printf("%8.3f %11.2f %13.4e %14.4e\n",frac,ang,
            std::abs(w1-w2)/std::max(w1,1e-300),
            (b2-target).norm()/std::max(b1.norm(),1e-300));
    }
}
