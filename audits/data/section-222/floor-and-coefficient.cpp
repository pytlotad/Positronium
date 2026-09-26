// Audit 222: 221h's two open items.
//
// ITEM 1.  221f bracketed the minimum of the realized eccentricity
// with three points and reported 6.141e-06 for k=2 there.  That is
// where the scan landed, not a floor.  Resolved here with a fine scan
// of the launch speed and a tighter tolerance.
//
// ITEM 2.  221f measured k=2 = 3.55 to 3.86 times the realized
// eccentricity, against the 2e of the small-e Bessel limit.  221d
// already contains the counter-case: at e = 5.3251e-04 = 10 beta^2 it
// read k=2 = 1.10169e-03 against 2e = 1.06502e-03, a ratio of 1.034.
// So the coefficient is 2 for a genuine Kepler eccentricity and 3.5
// for the relativistic residual, which would mean the two are not the
// same kind of oscillation.  Measured here as one curve.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

namespace {

struct Point { double ecc; double k2; bool ok; };

Point run(double a,double speedScale,double imposedEcc,int samples,
          double tol){
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double period=osculatingPeriod(-k/(2.0*a),k);
    // Built directly: for a speed above circular the osculating-element
    // helper has a negative discriminant and silently returns another
    // state (audit 221g).  Periapsis speed of an ellipse with the
    // requested eccentricity, times the scan factor.
    const double v=std::sqrt(k/a*(1.0+imposedEcc)/(1.0-imposedEcc))
                  *speedScale;
    const double r0=a*(1.0-imposedEcc);
    const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*1.0e-6);
    State s;
    const double w1=secondMass/(firstMass+secondMass);
    const double w2=-firstMass/(firstMass+secondMass);
    s.firstPosition=Vec3{r0,0,0}*w1;
    s.secondPosition=Vec3{r0,0,0}*w2;
    s.firstVelocity=Vec3{0,v,0}*w1;
    s.secondVelocity=Vec3{0,v,0}*w2;
    s.firstDipole=m1;
    s.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=tol; acc.maximumDepth=26;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(s,acc);
    double lo=1e300,hi=0.0;
    std::vector<Vec3> src; src.reserve(samples);
    for(int i=0;i<samples;++i){
        if(!engine.advance(s,period/samples)) return {0,0,false};
        src.push_back(s.firstAcceleration*firstCharge
                     +s.secondAcceleration*secondCharge);
        const double r=(s.firstPosition-s.secondPosition).norm();
        lo=std::min(lo,r); hi=std::max(hi,r);
    }
    double amp[2]={0,0};
    for(int kk=1;kk<=2;++kk){
        std::complex<double> cx{},cy{},cz{};
        for(int i=0;i<samples;++i){
            const double ph=-2.0*pi*kk*i/samples;
            const std::complex<double> w{std::cos(ph),std::sin(ph)};
            cx+=w*src[i].x; cy+=w*src[i].y; cz+=w*src[i].z;
        }
        amp[kk-1]=std::sqrt(std::norm(cx)+std::norm(cy)+std::norm(cz));
    }
    return {(hi-lo)/(hi+lo),amp[0]>0.0?amp[1]/amp[0]:0.0,true};
}

}  // namespace

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double a=pairBohrRadius({electron,positron});
    const double b2=pairCoulombStrength/(mu*c*c*a);
    std::printf("# a = a_pair, beta_rel^2 = %.6e, moments off,"
                " reaction off\n",b2);

    std::printf("\n# ITEM 1: fine launch scan, tol 1e-12, 512 samples\n");
    std::printf("%14s %14s %10s %14s %10s\n",
                "(v/vC - 1)/b2","e realized","/beta^2","k=2","k2/(2e)");
    for(double x:{-0.30,-0.20,-0.16,-0.14,-0.13,-0.125,-0.12,-0.11,
                  -0.10,-0.05,0.0}){
        const Point p=run(a,1.0+x*b2,0.0,512,1.0e-12);
        if(!p.ok){ std::printf("%14.4f  failed\n",x); continue; }
        std::printf("%14.4f %14.6e %10.4f %14.6e %10.3f\n",
                    x,p.ecc,p.ecc/b2,p.k2,
                    p.ecc>0.0?p.k2/(2.0*p.ecc):0.0);
    }

    std::printf("\n# ITEM 2: the coefficient against eccentricity\n");
    std::printf("%14s %14s %10s %14s %10s\n",
                "e imposed","e realized","/beta^2","k=2","k2/(2e)");
    for(double e:{0.0,1.0e-5,3.0e-5,1.0e-4,3.0e-4,1.0e-3,3.0e-3,
                  1.0e-2,3.0e-2,1.0e-1}){
        const Point p=run(a,1.0,e,512,1.0e-12);
        if(!p.ok){ std::printf("%14.3e  failed\n",e); continue; }
        std::printf("%14.3e %14.6e %10.4f %14.6e %10.3f\n",
                    e,p.ecc,p.ecc/b2,p.k2,
                    p.ecc>0.0?p.k2/(2.0*p.ecc):0.0);
    }
}
