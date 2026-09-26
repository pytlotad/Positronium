// Audit 221: is the circular k=2 floor a property of the model, or of
// the LAUNCH?  The realized eccentricity of a Coulomb-circular launch
// is 1.37269e-05 = 0.258 beta^2, pinned for every imposed e below it.
// If that is a real radial oscillation caused by launching with the
// Coulomb speed into a potential that carries O(beta^2) corrections,
// then the true circular orbit of the FULL potential should have a
// smaller oscillation -- and a smaller k=2 floor with it.  199d read
// the floor as "the relativistic correction and nothing else"; this
// asks whether it is instead radiation from a launch artefact.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=pairBohrRadius({electron,positron});
    const double b2=pairCoulombStrength/(mu*c*c*a);
    const int samples=256;
    std::printf("# a = a_pair, beta_rel^2 = %.4e, moments off\n",b2);
    std::printf("%12s %14s %10s %14s\n",
                "v/v_Coulomb","e realized","/beta^2","k=2");
    for(double s0:{1.0,1.0-0.5*b2,1.0-0.25*b2,1.0-0.125*b2,
                   1.0+0.125*b2,1.0+0.25*b2,1.0+0.5*b2,
                   1.0+1.0*b2,1.0+2.0*b2}){
        // Built directly, NOT through osculatingPeriapsisState: for
        // s0 > 1 the pair {-k/2a, sqrt(ka)*s0} has discriminant
        // 1 - s0^2 < 0, i.e. more angular momentum than that energy
        // allows, and the helper defends itself by returning a state
        // that is not the one asked for.  An earlier version of this
        // probe used it and reported a realized eccentricity pinned at
        // 1.3726e-05 for every s0 > 1 while the k=2 amplitude grew
        // fourfold -- the contradiction that exposed it.
        const double period=osculatingPeriod(-k/(2.0*a),k);
        const double v=std::sqrt(k/a)*s0;
        const Vec3 m1=Vec3{0,0,1}*(firstMagneticMoment*1.0e-6);
        State start;
        const double w1=secondMass/(firstMass+secondMass);
        const double w2=-firstMass/(firstMass+secondMass);
        start.firstPosition=Vec3{a,0,0}*w1;
        start.secondPosition=Vec3{a,0,0}*w2;
        start.firstVelocity=Vec3{0,v,0}*w1;
        start.secondVelocity=Vec3{0,v,0}*w2;
        start.firstDipole=m1;
        start.secondDipole=m1*(secondMagneticMoment/firstMagneticMoment);
        ClassicalTrajectoryEngine::Accuracy acc;
        acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
        acc.reactionModel=ChargeRadiationReactionModel::disabled;
        acc.computeOutwardFlux=false;
        ClassicalTrajectoryEngine engine(start,acc);
        State st=start;
        double lo=1e300,hi=0.0; bool ok=true;
        std::vector<Vec3> src; src.reserve(samples);
        for(int i=0;i<samples&&ok;++i){
            ok=engine.advance(st,period/samples);
            src.push_back(st.firstAcceleration*firstCharge
                         +st.secondAcceleration*secondCharge);
            const double r=(st.firstPosition-st.secondPosition).norm();
            lo=std::min(lo,r); hi=std::max(hi,r);
        }
        if(!ok){ std::printf("%12.9f  failed\n",s0); continue; }
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
        const double ecc=(hi-lo)/(hi+lo);
        std::printf("%12.9f %14.6e %10.4f %14.6e\n",
                    s0,ecc,ecc/b2,amp[0]>0.0?amp[1]/amp[0]:0.0);
    }
}
