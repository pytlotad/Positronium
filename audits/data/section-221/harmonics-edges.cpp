// Audit 221: the knife edge of photon generation.  A circular orbit
// radiates only the fundamental -- 199c/199d measured its k=2 residual
// at 3.3802e-05 with the moments on and 5.2869e-05 with them off, the
// latter being beta^2 and not Kepler content -- and eccentricity
// switches the harmonics on.  Where exactly?
//
// Transform of the E1 source d'' = q1 a1 + q2 a2 over one radial
// period.  200d showed Parseval holds exactly, so |X_k|^2/|X_1|^2 IS
// the power ratio, and 200c validated the model against the Bessel
// closed form to 2-8 parts in 1e4.  Reaction OFF so decay does not
// broaden the peaks, matching 199c; the moments are switched by
// scaling them down by 1e-6, which is 199d's method.
//
// Self-check: at e = 0 and a = a_pair the k=2 entry must reproduce
// 3.3802e-05 (moments on) and 5.2869e-05 (moments off), and at
// e = 0.10 it must reproduce 0.19957 and 0.19947.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>

namespace {

// |X_k|/|X_1| of the E1 source over one radial period.
bool harmonics(double semiMajor,double ecc,bool moments,int samples,
               std::vector<double>& out){
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const OsculatingElements el{-k/(2.0*semiMajor),
        std::sqrt(k*semiMajor*(1.0-ecc*ecc))};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const double scale=moments?1.0:1.0e-6;
    const Vec3 m1=Vec3{0.0,0.0,1.0}*(firstMagneticMoment*scale);
    const State start=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1.0e-10; acc.maximumDepth=24;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    acc.computeOutwardFlux=false;
    ClassicalTrajectoryEngine engine(start,acc);
    State s=start;
    std::vector<Vec3> source; source.reserve(samples);
    // Advance FIRST, then record.  osculatingPeriapsisState leaves
    // firstAcceleration at its default zero, so sampling before the
    // first step injects one null sample -- which is a delta in time
    // and therefore a FLAT pedestal across every k.  Measured, that
    // pedestal read 5.49e-03 at e = 0 where the true value is
    // 3.38e-05, and it sat at 5.49e-03, 5.54e-03, 5.54e-03 for
    // k = 2, 3, 4: flat, which is what gave it away.
    for(int i=0;i<samples;++i){
        if(!engine.advance(s,period/samples)) return false;
        source.push_back(s.firstAcceleration*firstCharge
                        +s.secondAcceleration*secondCharge);
    }
    out.clear();
    double fundamental=0.0;
    for(int kk=1;kk<=6;++kk){
        std::complex<double> cx{},cy{},cz{};
        for(int i=0;i<samples;++i){
            const double ph=-2.0*pi*kk*i/samples;
            const std::complex<double> w{std::cos(ph),std::sin(ph)};
            cx+=w*source[i].x; cy+=w*source[i].y; cz+=w*source[i].z;
        }
        const double a=std::sqrt(std::norm(cx)+std::norm(cy)+std::norm(cz));
        if(kk==1) fundamental=a;
        out.push_back(fundamental>0.0?a/fundamental:0.0);
    }
    return true;
}

}  // namespace

int main(int argc,char** argv){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double A=pairBohrRadius({electron,positron});
    const int samples=argc>1?std::atoi(argv[1]):256;
    const int mode=argc>2?std::atoi(argv[2]):0;
    std::printf("# %d samples per radial period, reaction off\n",samples);
    std::printf("# SELF-CHECK against 199c/199d: e=0 k2 must read\n");
    std::printf("# 3.3802e-05 (moments on) and 5.2869e-05 (off);\n");
    std::printf("# e=0.10 k2 must read 0.19957 and 0.19947.\n");
    std::printf("\n%8s %5s %8s %12s %12s %12s\n",
                "a/a_pair","mom","ecc","k=2","k=3","k=4");
    if(mode==0||mode==1)
    for(double f:{1.0})
        for(int mom=1;mom>=0;--mom)
            for(double e:{0.0,0.10,0.30,0.50,0.70}){
                std::vector<double> h;
                if(!harmonics(f*A,e,mom!=0,samples,h)){
                    std::printf("%8.4f %5d %8.4f  failed\n",f,mom,e);
                    continue; }
                std::printf("%8.4f %5d %8.4f %12.5e %12.5e %12.5e\n",
                            f,mom,e,h[1],h[2],h[3]);
            }

    // EDGE A: where the Kepler harmonic overtakes the relativistic
    // floor.  Moments off throughout, so the floor is pure beta^2.
    std::printf("\n# EDGE A -- k=2 against the circular floor.\n");
    std::printf("%9s %11s %11s %12s %9s\n",
                "a/a_pair","beta_rel^2","ecc","k=2","k2/floor");
    for(double f:{1.0,0.01}){
        if(mode==2&&f!=1.0) continue;
        if(mode==3&&f!=0.01) continue;
        if(mode==1) continue;
        const double mu=firstMass*secondMass/(firstMass+secondMass);
        const double b2=pairCoulombStrength/(mu*c*c*f*A);
        std::vector<double> h0;
        if(!harmonics(f*A,0.0,false,samples,h0)) continue;
        const double floor=h0[1];
        std::printf("%9.4f %11.4e %11s %12.5e %9s\n",
                    f,b2,"0 (floor)",floor,"1.000");
        for(double e:{0.2,0.5,1.0,2.0,4.0,8.0,20.0}){
            const double ecc=e*b2/2.0;
            std::vector<double> h;
            if(!harmonics(f*A,ecc,false,samples,h)) continue;
            std::printf("%9.4f %11s %11.4e %12.5e %9.3f\n",
                        f,"",ecc,h[1],h[1]/floor);
        }
    }

    // EDGE B: where k=2 overtakes k=1.  Closed form predicts 0.540055.
    std::printf("\n# EDGE B -- k=2 against k=1.  Predicted e = 0.540055\n");
    std::printf("%9s %12s %12s\n","ecc","k=2","k=3");
    for(double e:{0.50,0.52,0.53,0.535,0.540,0.545,0.55,0.58,0.60}){
        if(mode!=0&&mode!=4) break;
        std::vector<double> h;
        if(!harmonics(A,e,false,samples,h)) continue;
        std::printf("%9.4f %12.6f %12.6f\n",e,h[1],h[2]);
    }
}
