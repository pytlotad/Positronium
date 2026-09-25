// 181h/181i: the term-by-term check and the drift split across many
// configurations, plus the signed sin(2 phi) fit.

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
#include <initializer_list>
namespace {
Vec3 axial(const Vec3& v){ return {-v.x,v.y,-v.z}; }
Vec3 bmtTerm(const Vec3& v,const ElectromagneticField& f,double g,int term){
    const double a=0.5*(g-2.0), gam=gamma(v);
    const Vec3 beta=v/c;
    if(term==0) return f.magnetic*(a+1.0/gam);
    if(term==1) return beta*(-a*gam/(gam+1.0)*dot(beta,f.magnetic));
    return cross(beta,f.electric)*(-(a+1.0/(gam+1.0))/c);
}
struct R { double res, normRes; };
R at(double bd,double thDeg,double phiDeg,double rf,int term=-1){
    const double total=firstMass+secondMass;
    const double r0=rf*pairBohrRadius(activePair);
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thDeg*pi/180, ph=phiDeg*pi/180;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),
                   std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity={bd*c,speed*secondMass/total,0};
    s.secondVelocity={bd*c,-speed*firstMass/total,0};
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    Vec3 b1,b2;
    if(term<0){
        b1=thomasBmtEffectiveField(s.firstVelocity,f.atFirst,firstGFactor);
        b2=thomasBmtEffectiveField(s.secondVelocity,f.atSecond,secondGFactor);
    } else {
        b1=bmtTerm(s.firstVelocity,f.atFirst,firstGFactor,term);
        b2=bmtTerm(s.secondVelocity,f.atSecond,secondGFactor,term);
    }
    const double n=std::max(b1.norm(),1e-300);
    return {(b2-axial(b1)).norm()/n, b2.norm()/n-1.0};
}
}
int main(){
    std::printf("180j item 2: the term-by-term check across "
                "configurations.\n");
    std::printf("%7s %6s %8s %12s %12s %12s\n","r/a","theta","beta_dr",
                "magnetic","longitudinal","motional");
    for(double rf:{1.0,0.1,0.01}) for(double t:{20.0,70.0})
        for(double bd:{0.001,0.05})
            std::printf("%7.3g %6.0f %8.3g %12.3e %12.3e %12.3e\n",rf,t,bd,
                at(bd,t,0.0,rf,0).res,at(bd,t,0.0,rf,1).res,
                at(bd,t,0.0,rf,2).res);

    std::printf("\n180j item 3: is the vector residual drift-free, and the\n"
                "norm residual exactly linear in the drift, generally?\n");
    std::printf("%7s %6s %10s %14s %14s %12s\n","r/a","phi","beta_dr",
                "vector res","norm res","norm/beta");
    for(double rf:{1.0,0.01}) for(double p:{20.0,45.0,80.0})
        for(double bd:{1e-5,1e-3,1e-1}){
            const R r=at(bd,40.0,p,rf);
            std::printf("%7.3g %6.0f %10.3g %14.6e %14.4e %12.4e\n",
                        rf,p,bd,r.res,r.normRes,r.normRes/bd);
        }

    std::printf("\n180j item 4: the norm residual, SIGNED fit\n"
                "over a dense one-degree scan.\n");
    double sxx=0,sxy=0,syy=0; int n=0;
    double worst=0.0, worstPhi=0.0;
    for(int i=1;i<180;++i){
        const double p=i*1.0;
        const double x=std::sin(2.0*p*pi/180);
        const double y=at(0.01,40.0,p,1.0).normRes;
        sxx+=x*x; sxy+=x*y; syy+=y*y; ++n;
    }
    const double A=sxy/sxx;
    for(int i=1;i<180;++i){
        const double p=i*1.0;
        const double x=std::sin(2.0*p*pi/180);
        const double y=at(0.01,40.0,p,1.0).normRes;
        const double d=std::abs(y-A*x);
        if(d>worst){ worst=d; worstPhi=p; }
    }
    const double rms=std::sqrt(std::max(0.0,syy-A*sxy)/n);
    std::printf("   points %d   A = %.6e\n",n,A);
    std::printf("   rms residual of the fit      %.4e  (%.3f%% of A)\n",
                rms,100.0*rms/A);
    std::printf("   worst deviation              %.4e at phi = %.0f deg\n",
                worst,worstPhi);
}
