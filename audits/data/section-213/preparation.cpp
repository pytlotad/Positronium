// Audit 213 control: are the para and ortho runs the SAME preparation
// with one sign flipped, or two different orbits?  Near-zero observation
// time, so nothing integrates; frames.front() is the pristine state.
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    gSpinQuantization=true;
    SimulationOptions o; o.frameCount=2; o.observationTime=1.0e-24;
    std::printf("%5s %3s %16s %16s %10s %10s %9s\n",
        "seed","ph","E_rel","L_orb","r0/a_pair","ecc","mu1.mu2");
    for(std::uint64_t seed=1;seed<=8;++seed)
        for(int ph=1;ph<=2;++ph){
            const SimulationResult r=simulate(seed,ph,o);
            if(r.frames.empty()){ std::printf("%5llu %3d empty\n",
                (unsigned long long)seed,ph); continue; }
            const Frame& s=r.frames.front();
            const Vec3 sep=s.first-s.second;
            const double d1=s.firstDipole.norm(),d2=s.secondDipole.norm();
            const double al=d1>0.0&&d2>0.0
                ?dot(s.firstDipole,s.secondDipole)/(d1*d2):0.0;
            const double mu=firstMass*secondMass/(firstMass+secondMass);
            const double eps=r.initial.relativeEnergy/mu;
            const double ell=r.initial.orbitalAngularMomentum/mu;
            const double k=pairCoulombStrength/mu;
            const double disc=1.0+2.0*eps*ell*ell/(k*k);
            std::printf("%5llu %3d %16.8e %16.8e %10.6f %10.6f %9.5f\n",
                (unsigned long long)seed,ph,r.initial.relativeEnergy,
                r.initial.orbitalAngularMomentum,
                sep.norm()/pairBohrRadius({electron,positron}),
                disc>0.0?std::sqrt(disc):-1.0,al);
        }
}
