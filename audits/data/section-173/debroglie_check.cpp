#include "modules/crem_trajectory.hpp"
#include <cstdio>
int main(){
    // de Broglie standing wave on a circular orbit: n * lambda = 2 pi r,
    // lambda = h/p, p = mu v, and v = sqrt(k/(mu r)) on the circle.
    // => 2 pi sqrt(mu k r) = n h  =>  r = n^2 hbar^2/(mu k).
    const double k=pairCoulombStrength, mu=pairReducedMass;
    std::printf("n   de Broglie radius [m]   a_pair*n^2 [m]       relative\n");
    for(int n=1;n<=3;++n){
        const double r=n*n*hbar*hbar/(mu*k);
        const double a=pairBohrRadius(activePair)*n*n;
        std::printf("%d   %.15e   %.15e   %.3e\n",n,r,a,std::abs(r-a)/a);
    }
    // and the angular momentum it selects
    const double r1=hbar*hbar/(mu*k);
    const double v1=std::sqrt(k/(mu*r1));
    std::printf("\nL at the n=1 de Broglie circle = %.15e hbar\n",
                mu*v1*r1/hbar);
    std::printf("the model's own preparation uses L = 1 hbar.\n");
}
