// Audit 223: the alpha^-5 identity, checked against the code's own
// constants rather than against tables, and against the model's own
// measured collapse time.
//
// Claim: the quasi-circular Larmor inspiral from L = hbar has the
// closed form t = a^3 mu^2 c^3 / (4 k^2) = lambda_C/(2 alpha^5 c),
// while the 2-gamma width gives tau = 2 lambda_C/(alpha^5 c), so the
// same alpha^-5 governs both and the ratio is exactly 4.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double m=electron.mass;
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength;               // = alpha hbar c
    const double a=pairBohrRadius({electron,positron});
    const double lambdaC=hbar/(m*c);                  // reduced Compton
    const double alpha=k/(hbar*c);                    // from the code's k
    std::printf("From the code's own constants:\n");
    std::printf("  pairCoulombStrength k   = %.10e J m\n",k);
    std::printf("  alpha = k/(hbar c)      = %.12f  (1/alpha = %.6f)\n",
                alpha,1.0/alpha);
    std::printf("  lambda_C = hbar/(m c)   = %.10e m\n",lambdaC);
    std::printf("  a_pair                  = %.10e m\n",a);
    std::printf("  2 lambda_C/alpha        = %.10e m   ratio %.12f\n",
                2.0*lambdaC/alpha,a/(2.0*lambdaC/alpha));
    const double a5=std::pow(alpha,5.0);
    // Quasi-circular Larmor inspiral: r^2 dr = -(4/3) k^2/(c^3 mu^2) dt
    const double C=(4.0/3.0)*k*k/(c*c*c*mu*mu);
    const double tClosed=a*a*a/(3.0*C);
    const double tCompton=lambdaC/(2.0*a5*c);
    const double tau2g=2.0*lambdaC/(a5*c);
    std::printf("\nInspiral, two routes to the same number:\n");
    std::printf("  a^3/(3C), C = 4k^2/(3 mu^2 c^3) = %.10e s = %.6f ps\n",
                tClosed,tClosed*1e12);
    std::printf("  lambda_C/(2 alpha^5 c)          = %.10e s = %.6f ps\n",
                tCompton,tCompton*1e12);
    std::printf("  ratio                            = %.14f\n",
                tClosed/tCompton);
    std::printf("\n2-gamma width, leading order:\n");
    std::printf("  2 lambda_C/(alpha^5 c)          = %.10e s = %.6f ps\n",
                tau2g,tau2g*1e12);
    std::printf("  tau_2gamma / t_inspiral          = %.14f\n",
                tau2g/tClosed);
    std::printf("\nAgainst the measured para-Ps lifetime 125.164 ps:\n");
    std::printf("  leading order                    = %.6f ps\n",
                tau2g*1e12);
    std::printf("  shortfall                        = %.4f %%\n",
                100.0*(125.164e-12-tau2g)/tau2g);
    std::printf("  alpha (5 - pi^2/4)/pi            = %.4f %%\n",
                100.0*alpha*(5.0-pi*pi/4.0)/pi);
    // The model's own collapse time, Richardson-extrapolated with the
    // order 0.746 that audit 213f measured on seed 1's three steps.
    const double v1=3.066404061096e-11;   // s_max = 0.30, audit 213
    const double v2=3.099899406385e-11;   // s_max = 0.075, audit 214
    const double p=0.746, ratio=std::pow(0.25,p);
    const double extrap=(v2-v1*ratio)/(1.0-ratio);
    std::printf("\nThe model's own measured collapse time (seed 1, para):\n");
    std::printf("  s_max = 0.30                     = %.6f ps\n",v1*1e12);
    std::printf("  s_max = 0.075                    = %.6f ps\n",v2*1e12);
    std::printf("  Richardson at order %.3f         = %.6f ps\n",p,
                extrap*1e12);
    std::printf("  against the closed form           %.6f ps\n",
                tClosed*1e12);
    std::printf("  gap                              = %.4f %%\n",
                100.0*(extrap-tClosed)/tClosed);
}
