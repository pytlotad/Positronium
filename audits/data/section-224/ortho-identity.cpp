// Audit 224: the same test for ortho.  The 3-gamma width is
//   Gamma_3g = (2/(9 pi)) (pi^2 - 9) alpha^6 m c^2 / hbar,
// so  tau_3g = [9 pi / (2 (pi^2 - 9))] lambda_C/(alpha^6 c),
// one power of alpha HIGHER than the 2-gamma width and therefore one
// power higher than the classical inspiral, which is alpha^-5 and is
// the same orbit for both channels (audit 213b: the estimator prepares
// one circular orbit at a_pair and the seed sets only the moment).
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double m=electron.mass;
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength;
    const double a=pairBohrRadius({electron,positron});
    const double lambdaC=hbar/(m*c);
    const double alpha=k/(hbar*c);
    const double C=(4.0/3.0)*k*k/(c*c*c*mu*mu);
    const double tIn=a*a*a/(3.0*C);
    const double a5=std::pow(alpha,5.0),a6=std::pow(alpha,6.0);
    const double tau2=2.0*lambdaC/(a5*c);
    const double coeff3=9.0*pi/(2.0*(pi*pi-9.0));
    const double tau3=coeff3*lambdaC/(a6*c);
    std::printf("alpha from the code   1/%.6f\n",1.0/alpha);
    std::printf("t_inspiral            %.6f ps   (alpha^-5)\n",tIn*1e12);
    std::printf("tau_2gamma  LO        %.6f ps   (alpha^-5)\n",tau2*1e12);
    std::printf("tau_3gamma  LO        %.6f ns   (alpha^-6)\n",tau3*1e9);
    std::printf("  coefficient 9pi/(2(pi^2-9)) = %.6f\n",coeff3);
    std::printf("\nRatios against the SAME classical inspiral:\n");
    std::printf("  tau_2gamma / t_inspiral = %.10f   (pure number)\n",
                tau2/tIn);
    std::printf("  tau_3gamma / t_inspiral = %.4f\n",tau3/tIn);
    std::printf("  9pi/(pi^2-9) / alpha    = %.4f   ratio %.12f\n",
                (9.0*pi/(pi*pi-9.0))/alpha,
                (tau3/tIn)/((9.0*pi/(pi*pi-9.0))/alpha));
    std::printf("\nAnd against each other:\n");
    std::printf("  tau_3gamma / tau_2gamma = %.4f\n",tau3/tau2);
    std::printf("  9pi/(4(pi^2-9)) / alpha = %.4f\n",
                (9.0*pi/(4.0*(pi*pi-9.0)))/alpha);
    std::printf("  measured 142.05 ns / 125.164 ps = %.4f\n",
                142.05e-9/125.164e-12);
    std::printf("\nQED correction to the 3-gamma width:\n");
    std::printf("  LO           %.4f ns\n",tau3*1e9);
    std::printf("  measured     142.05 ns\n");
    std::printf("  shortfall    %.4f %%\n",100.0*(142.05e-9-tau3)/tau3);
    std::printf("  10.286 alpha/pi = %.4f %%\n",
                100.0*10.286*alpha/pi);
    std::printf("\nWhat the MODEL gives for the same two channels:\n");
    std::printf("  ortho/para - 1, 24 seeds at s_max 0.075 (214b)"
                "  %.4e\n",6.9918e-04);
    std::printf("  required for the physical ratio                 "
                "  %.4e\n",tau3/tau2-1.0);
    std::printf("  shortfall factor                                "
                "  %.4e\n",(tau3/tau2-1.0)/6.9918e-04);
}
