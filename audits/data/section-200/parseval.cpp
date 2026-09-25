#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
#include <vector>
#include <complex>
#include <initializer_list>
namespace {
void scan(double ecc,int N){
    const double total=firstMass+secondMass;
    const double aSemi=pairBohrRadius(activePair);
    const double period=2.0*pi*std::sqrt(pairReducedMass*aSemi*aSemi*aSemi
                                         /pairCoulombStrength);
    const double r=aSemi*(1.0+ecc);
    const double vsq=(pairCoulombStrength/pairReducedMass)
                     *(2.0/r-1.0/aSemi);
    const double v=std::sqrt(vsq);
    const double th=40.0*pi/180;
    const Vec3 n1{std::sin(th),0.0,std::cos(th)};
    State s{};
    s.firstPosition={r*secondMass/total,0,0};
    s.secondPosition={-r*firstMass/total,0,0};
    s.firstVelocity={0,v*secondMass/total,0};
    s.secondVelocity={0,-v*firstMass/total,0};
    s.firstProperDipole=n1*firstMagneticMoment;
    s.secondProperDipole=n1*(-secondMagneticMoment);
    synchronizeCovariantDipoles(s);
    ClassicalTrajectoryEngine::Accuracy acc;
    acc.relativeTolerance=1e-11; acc.maximumDepth=22;
    acc.computeOutwardFlux=false;
    acc.reactionModel=ChargeRadiationReactionModel::disabled;
    ClassicalTrajectoryEngine engine(s,acc);
    std::vector<Vec3> d2;
    for(int i=0;i<N;++i){
        if(i&&!engine.advance(s,period/N)) break;
        const MutualForces mf=mutualForces(s);
        d2.push_back(relativisticAcceleration(s.firstVelocity,mf.first,
                         firstMass)*firstCharge
                    +relativisticAcceleration(s.secondVelocity,mf.second,
                         secondMass)*secondCharge);
    }
    const int M=int(d2.size());
    if(M<32) return;
    double mean=0;
    for(int i=0;i<M;++i) mean+=d2[size_t(i)].squaredNorm();
    mean/=double(M);
    std::vector<double> amp(size_t(M/2+1),0.0);
    double sumsq=0;
    for(int k=0;k<=M/2;++k){
        std::complex<double> cx{0,0},cy{0,0},cz{0,0};
        for(int i=0;i<M;++i){
            const double ph=-2.0*pi*double(k)*double(i)/double(M);
            const std::complex<double> w{std::cos(ph),std::sin(ph)};
            cx+=w*d2[size_t(i)].x; cy+=w*d2[size_t(i)].y;
            cz+=w*d2[size_t(i)].z;
        }
        const double a2=(std::norm(cx)+std::norm(cy)+std::norm(cz))
                        /(double(M)*double(M));
        amp[size_t(k)]=std::sqrt(a2);
        sumsq+=(k==0||k==M/2)?a2:2.0*a2;
    }
    std::printf("  %5.2f %6d",ecc,M);
    for(int k=1;k<=4;++k)
        std::printf(" %11.5e",amp[size_t(k)]/amp[1]);
    std::printf("  %12.9f\n",sumsq/mean);
    std::fflush(stdout);
}
}
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    std::printf("199f items 2 and 4.  Last column is Parseval: the summed\n"
                "harmonic power over the mean |d''|^2.  If it is 1,\n"
                "the transform of the SOURCE already IS the power spectrum\n"
                "and there is no extra k-weighting to worry about.\n\n");
    std::printf("  %5s %6s %11s %11s %11s %11s %12s\n","ecc","N","k=1","k=2",
                "k=3","k=4","Parseval");
    for(double e:{0.0,0.3,0.7}){ scan(e,256); scan(e,1024); }
}
