// Audit 226: is the dipole sector's energy balance a test of physics
// or of numerics?
//
// Section 69 asked whether one Lagrangian gives force, precession and
// ledger together, measured the lab-moment velocity dependence, found
// it <= 2.2e-04 k/r0 over three periods where the model is valid, and
// decided not to implement.  That answers "how big", not "is the test
// empty".  The sharper question is separable by the method 215-218
// established: if the drift of the conserved energy FALLS with the
// integration step it is numerics, and if it FLOORS it is the missing
// term.
//
// The force is made conservative on purpose here:
// useRetardedExternalForces = false isolates the Coulomb-Darwin-dipole
// force, and the reaction is disabled, so a Lagrangian system would
// conserve conservedParticleEnergy exactly.  Any floor is then the
// part of the velocity dependence no force differentiates.
//
// conservedParticleEnergy already carries HALF the Legendre transform:
// its own comment notes the charge-dipole term is -p.E with
// p = (v x mu)/c^2, homogeneous of degree one in the velocities, so it
// cancels in H = sum v.(dL/dv) - L and is excluded by weight 0.  What
// is NOT handled is the lab moments' own velocity dependence, which
// enters pairDipoleInteractionEnergy through the boost.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double a=0.25*pairBohrRadius({electron,positron});
    const double scale=pairCoulombStrength/a;      // k/r0
    const OsculatingElements el{-k/(2.0*a),std::sqrt(k*a)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    std::printf("# 0.25 a_pair, ONE period, retarded OFF, reaction OFF.\n");
    std::printf("# reaction OFF, so a Lagrangian system conserves\n");
    std::printf("# conservedParticleEnergy exactly.  Drifts in k/r0.\n");
    std::printf("\n%7s %11s %13s %13s %13s\n",
                "tilt","tolerance","range E_cons","range E_full",
                "range U_dd");
    for(double tiltDeg:{0.0,90.0}){
        const double th=tiltDeg*pi/180.0;
        const Vec3 unit{std::sin(th),0.0,std::cos(th)};
        for(double tol:{1.0e-12,1.0e-13,1.0e-14}){
            const State start=osculatingPeriapsisState(
                el,k,unit*firstMagneticMoment,
                unit*secondMagneticMoment,
                Vec3{0,0,1},Vec3{1,0,0},0.0);
            ClassicalTrajectoryEngine::Accuracy acc;
            acc.relativeTolerance=tol; acc.maximumDepth=26;
            acc.reactionModel=ChargeRadiationReactionModel::disabled;
            acc.computeOutwardFlux=false;
            acc.useRetardedExternalForces=false;
            ClassicalTrajectoryEngine engine(start,acc);
            State s=start;
            double loC=1e300,hiC=-1e300,loF=1e300,hiF=-1e300;
            double loU=1e300,hiU=-1e300;
            const int n=512; bool ok=true;
            for(int i=0;i<n&&ok;++i){
                ok=engine.advance(s,period/n);
                const double ec=conservedParticleEnergy(s);
                const double ef=conservativeParticleEnergy(s);
                const double ud=pairDipoleInteractionEnergy(
                    s.firstPosition-s.secondPosition,
                    s.firstDipole,s.secondDipole);
                loC=std::min(loC,ec); hiC=std::max(hiC,ec);
                loF=std::min(loF,ef); hiF=std::max(hiF,ef);
                loU=std::min(loU,ud); hiU=std::max(hiU,ud);
            }
            if(!ok){ std::printf("%7.1f %11.0e  failed\n",tiltDeg,tol);
                     continue; }
            std::printf("%7.1f %11.0e %13.5e %13.5e %13.5e\n",
                tiltDeg,tol,(hiC-loC)/scale,(hiF-loF)/scale,
                (hiU-loU)/scale);
        }
    }
    std::printf("\n# A range that falls with the tolerance is the\n");
    std::printf("# integrator.  One that floors is the term no force\n");
    std::printf("# differentiates.  U_dd is the scale it lives on.\n");
}
