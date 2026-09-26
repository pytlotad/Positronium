// Audit 213: E, L and r are seed-independent (prep.cpp).  What DOES the
// seed set?  Orientation: the moment direction and the orbital plane.
// delta is the angle between the first moment and the orbit normal.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    SimulationOptions o; o.frameCount=2; o.observationTime=1.0e-24;
    std::printf("%5s %3s %8s %8s %8s %8s %8s %8s %9s\n",
        "seed","ph","mu1_x","mu1_y","mu1_z","n_x","n_y","n_z","delta");
    for(std::uint64_t seed=1;seed<=24;++seed)
        for(int ph=1;ph<=1;++ph){
            const SimulationResult r=simulate(seed,ph,o);
            if(r.frames.empty()) continue;
            const Frame& s=r.frames.front();
            const Vec3 m=s.firstDipole*(1.0/s.firstDipole.norm());
            const Vec3 L=s.noetherAngularMomentum;
            const Vec3 n=L*(1.0/L.norm());
            std::printf("%5llu %3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f"
                " %9.4f\n",(unsigned long long)seed,ph,
                m.x,m.y,m.z,n.x,n.y,n.z,
                std::acos(std::max(-1.0,std::min(1.0,dot(m,n))))
                    *180.0/3.14159265358979324);
        }
}
