// Audit 219: is the sharp onset of 218g's covariance residual between
// beta = 0.35 and 0.50 the fixed 50-direction flux quadrature failing
// on a beamed pattern?
//
// electrodynamics.hpp:1672 calls electromagneticFieldFluxRates(state,
// history) with a DEFAULT FarFieldSampling, i.e. directionCount = 50,
// and ClassicalTrajectoryEngine::Accuracy exposes no way to change it.
// So every accumulated radiated four-momentum in the model is a
// degree-11 Lebedev integral.  The Lienard angular distribution goes
// as (1 - n.beta)^-k with k up to 6, whose dynamic range over the
// sphere is ((1+beta)/(1-beta))^k -- 85 at beta = 0.35, 729 at 0.50
// and 3.2e4 at 0.70.  The exact integral is known:
//   int dOmega/(1-beta mu)^n
//     = 2 pi /((n-1) beta) [ (1-beta)^(1-n) - (1+beta)^(1-n) ],
// so the rule's own error can be measured directly against it.
#include "modules/electrodynamics.hpp"
#include <cstdio>
#include <cmath>

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const auto exactIntegral=[](double beta,int n){
        if(beta<1.0e-12) return 4.0*pi;
        return 2.0*pi/((n-1)*beta)
              *(std::pow(1.0-beta,1-n)-std::pow(1.0+beta,1-n));
    };
    const auto quadrature=[](double beta,int n,int directions){
        double sum=0.0;
        for(const SphereQuadraturePoint& q:
                sphereQuadratureView(directions))
            sum+=q.solidAngleWeight
                *std::pow(1.0-beta*q.direction.z,-n);
        return sum;
    };
    std::printf("# relative error of the sphere rule on (1-beta mu)^-n,\n");
    std::printf("# the Lienard beaming factor.  n = 6 is the Poynting\n");
    std::printf("# integrand's steepest term.\n");
    for(int n:{4,5,6}){
        std::printf("\n# n = %d\n",n);
        std::printf("%7s %12s %14s %14s %14s %14s\n",
                    "beta","range","50 dirs","194 dirs","302 dirs",
                    "218g residual");
        const double ref[]={3.72e-07,1.33e-06,7.25e-07,1.89e-06,
                            6.84e-05,3.55e-04};
        int index=0;
        for(double beta:{0.05,0.10,0.20,0.35,0.50,0.70}){
            const double exact=exactIntegral(beta,n);
            std::printf("%7.2f %12.4g %14.6e %14.6e %14.6e %14.2e\n",
                beta,std::pow((1.0+beta)/(1.0-beta),n),
                std::abs(quadrature(beta,n,50)-exact)/exact,
                std::abs(quadrature(beta,n,194)-exact)/exact,
                std::abs(quadrature(beta,n,302)-exact)/exact,
                ref[index++]);
        }
    }
}
