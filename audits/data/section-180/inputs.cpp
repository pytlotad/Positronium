// 180b/180e/180f/180g: which relation each INPUT obeys, whether the
// breaking needs the drift, and the chord split
// Build: g++ -std=c++20 -O2 -I <repo root> $(root-config --cflags) \
//        inputs.cpp $(root-config --libs)

#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <initializer_list>
namespace {
struct Out { Vec3 b1,b2,e1,e2,B1,B2,v1,v2; };
Out probe(double phiDeg,const Vec3& drift,double r0,double thDeg=40.0){
    const double total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r0));
    const double th=thDeg*pi/180, ph=phiDeg*pi/180;
    const Vec3 dir{std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),
                   std::cos(th)};
    State s{};
    s.firstPosition={r0*secondMass/total,0,0};
    s.secondPosition={-r0*firstMass/total,0,0};
    s.firstVelocity=Vec3{0,speed*secondMass/total,0}+drift;
    s.secondVelocity=Vec3{0,-speed*firstMass/total,0}+drift;
    s.firstProperDipole=dir*firstMagneticMoment;
    s.secondProperDipole=dir*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    const StateHistory h=causalInitialHistory(s);
    const LocalElectromagneticFields f=localRelativisticFields(s,h);
    return {thomasBmtEffectiveField(s.firstVelocity,f.atFirst,firstGFactor),
            thomasBmtEffectiveField(s.secondVelocity,f.atSecond,
                                    secondGFactor),
            f.atFirst.electric,f.atSecond.electric,
            f.atFirst.magnetic,f.atSecond.magnetic,
            s.firstVelocity,s.secondVelocity};
}
// the eight sign patterns; index bit k set => component k is negated
Vec3 apply(const Vec3& v,int p){
    return {(p&1)?-v.x:v.x,(p&2)?-v.y:v.y,(p&4)?-v.z:v.z};
}
double res(const Vec3& a,const Vec3& b,int p){
    return (b-apply(a,p)).norm()/std::max(a.norm(),1e-300);
}
const char* tag(int p){
    static char s[8];
    s[0]=(p&1)?'-':'+'; s[1]='x'; s[2]=(p&2)?'-':'+'; s[3]='y';
    s[4]=(p&4)?'-':'+'; s[5]='z'; s[6]=0; return s;
}
}
int main(){
    const double a=pairBohrRadius(activePair);
    const Vec3 dx{0.01*c,0,0}, zero{0,0,0};
    std::printf("A. which sign pattern each field obeys (phi=0, drift x)\n");
    {   const Out o=probe(0.0,dx,a);
        const char* nm[4]={"B_eff","E","B","v"};
        const Vec3* aa[4]={&o.b1,&o.e1,&o.B1,&o.v1};
        const Vec3* bb[4]={&o.b2,&o.e2,&o.B2,&o.v2};
        for(int k=0;k<4;++k){
            int best=0; double bv=1e300;
            for(int p=0;p<8;++p){ const double r=res(*aa[k],*bb[k],p);
                if(r<bv){bv=r;best=p;} }
            std::printf("   %-6s best pattern %s  residual %10.3e\n",
                        nm[k],tag(best),bv);
        }
    }
    std::printf("\nB. does the breaking need the drift at all? (phi=45)\n");
    std::printf("%12s %14s %14s\n","beta_drift","identity res","|B2|/|B1|-1");
    for(double bd:{0.0,1e-6,1e-4,0.01,0.1}){
        const Out o=probe(45.0,Vec3{bd*c,0,0},a);
        std::printf("%12.3g %14.4e %14.3e\n",bd,res(o.b1,o.b2,1|4),
                    o.b2.norm()/o.b1.norm()-1.0);
    }
    std::printf("\nC. the residual as a chord: split delta about "
                "the target\n");
    std::printf("%6s %11s %11s %11s %11s\n","phi","|d|/|B1|","par/|d|",
                "|d|/2|B1|","perp/|d|");
    for(double p:{15.0,30.0,45.0,60.0,90.0}){
        const Out o=probe(p,dx,a);
        const Vec3 t=apply(o.b1,1|4);
        const Vec3 d=o.b2-t;
        const Vec3 u=t/t.norm();
        const double par=dot(d,u), perp=(d-u*par).norm();
        std::printf("%6.0f %11.4e %11.4e %11.4e %11.4e\n",p,
            d.norm()/o.b1.norm(),
            par/d.norm(),d.norm()/(2.0*o.b1.norm()),perp/d.norm());
    }
    std::printf("\nD. norm split vs sin(2 phi), normalised at phi=45\n");
    std::printf("%6s %13s %10s %10s\n","phi","|B2|/|B1|-1","ratio","sin2phi");
    const double ref=[&]{ const Out o=probe(45.0,dx,a);
        return o.b2.norm()/o.b1.norm()-1.0; }();
    for(double p:{15.0,30.0,45.0,60.0,75.0,90.0}){
        const Out o=probe(p,dx,a);
        const double v=o.b2.norm()/o.b1.norm()-1.0;
        std::printf("%6.0f %13.4e %10.4f %10.4f\n",p,v,v/ref,
                    std::sin(2.0*p*pi/180));
    }
    std::printf("\nE. is 2.25e-02/deg general?  slope from phi=0 to 2 deg\n");
    std::printf("%9s %9s %9s %14s\n","r/a_pair","theta","beta_dr",
                "slope /deg");
    for(double rf:{1.0,0.1,0.01}) for(double t:{20.0,40.0,70.0})
        for(double bd:{0.0,0.01}){
            const Out o=probe(2.0,Vec3{bd*c,0,0},rf*a,t);
            std::printf("%9.3g %9.0f %9.3g %14.4e\n",rf,t,bd,
                        res(o.b1,o.b2,1|4)/2.0);
        }
}
