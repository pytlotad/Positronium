// The action labels of section 131 against the model's OWN radial potential
// (audit 132).
//
// Section 131 computed n = (K/hbar) sqrt(mu/2|E|) and J_phi/h = L/hbar from
// the Kepler formulas.  For a CIRCULAR orbit the two must agree, since a
// circle has J_r = 0; the trace instead showed them drifting apart, to
// -0.0158 h at the deep end.  The model's radial function is not Keplerian:
// dipoleAwarePeriapsis solves
//
//   h(r) = eps + K/r - u_dd(r) - u_so(r,L) - L^2/(2 r^2),
//
// with u_dd the azimuth-averaged dipole-dipole energy and u_so the
// charge-dipole (spin-orbit) energy evaluated on the turning-point state.
// This probe puts the model's own circular orbit at a given radius -- L from
// dU/dr = 0 in that potential, E from U at that radius -- and prints the
// Kepler energy label against L/hbar, which is the gap the trace reported.
// Running it with the two dipole terms zeroed is the attribution test.
//
// The two energy functions are CALLED, not copied; only the turning-point
// state construction is mirrored from dipoleAwarePeriapsis, because the model
// exposes the potential nowhere else.
//
// Usage: radial_action_full_potential [0 to zero the dipole terms]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/radial_action_full_potential.cpp -o /tmp/radial $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc,char** argv) {
    const bool withDipoles=!(argc>1&&std::string(argv[1])=="0");
    const double axis=pairBohrRadius(activePair);
    const double reducedMass=reducedMassOf(activePair);
    const double attraction=pairCoulombStrength/reducedMass;  // specific
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 first=withDipoles?normal*firstMagneticMoment:Vec3{};
    const Vec3 second=withDipoles?normal*secondMagneticMoment:Vec3{};
    const double totalMass=firstMass+secondMass;
    // The two non-Keplerian pieces of the model's radial potential, specific,
    // exactly as dipoleAwarePeriapsis assembles them.
    const auto dipoleDipole=[&](double r) {
        return azimuthAveragedDipoleEnergy(r,first,second,normal)/reducedMass;
    };
    const auto spinOrbit=[&](double r,double angular) {
        if(!withDipoles||!(r>0.0)||!(totalMass>0.0)) return 0.0;
        const Vec3 radialHat{0.0,-1.0,0.0};      // cross(x, z) normalized
        const Vec3 tangentialHat=cross(normal,radialHat);
        const double tangentialSpeed=angular/r;  // v_r = 0 at the turning point
        State turningPoint{};
        turningPoint.firstPosition=radialHat*(r*secondMass/totalMass);
        turningPoint.secondPosition=radialHat*(-r*firstMass/totalMass);
        turningPoint.firstVelocity=
            tangentialHat*(tangentialSpeed*secondMass/totalMass);
        turningPoint.secondVelocity=
            tangentialHat*(-tangentialSpeed*firstMass/totalMass);
        turningPoint.firstProperDipole=first;
        turningPoint.secondProperDipole=second;
        synchronizeCovariantDipoles(turningPoint);
        return chargeDipoleInteractionEnergy(turningPoint)/reducedMass;
    };
    // Specific potential energy including the centrifugal term.
    const auto effective=[&](double r,double angular) {
        return -attraction/r+dipoleDipole(r)+spinOrbit(r,angular)
            +angular*angular/(2.0*r*r);
    };
    // The model's own circular orbit at radius r: L solves dU/dr = 0 there.
    // u_so depends on L as well, so iterate; the Kepler value seeds it.
    const auto circularAngularMomentum=[&](double r) {
        double angular=std::sqrt(attraction*r);
        for(int iteration=0;iteration<200;++iteration) {
            const double step=1.0e-6*r;
            // d/dr of everything except the centrifugal term, at fixed L.
            const auto withoutCentrifugal=[&](double radius) {
                return -attraction/radius+dipoleDipole(radius)
                    +spinOrbit(radius,angular);
            };
            const double slope=
                (withoutCentrifugal(r+step)-withoutCentrifugal(r-step))
                /(2.0*step);
            const double next=std::sqrt(std::max(0.0,slope*r*r*r));
            if(std::abs(next-angular)<=1.0e-14*std::max(angular,1.0e-300)) {
                angular=next; break;
            }
            angular=next;
        }
        return angular;
    };
    std::printf("dipole terms %s; a_pair = %.3f r*\n",
                withDipoles?"ON (the model as it runs)":"ZEROED (control)",
                axis/comptonBarrierRadius);
    std::printf("%12s %10s %14s %14s %12s %12s %12s\n","a/a_pair","r*",
                "n from E","J_phi/h = L/hbar","gap","|u_dd/u_C|","|u_so/u_C|");
    // The ladder's rungs (n = 1, 0.5774, 0.2733, 0.0947 -> a/a_pair = n^2)
    // and the measured stop.
    for(double fraction:{1.0,0.333333,0.074670,0.009369,0.008976,0.006505}) {
        const double r=axis*fraction;
        const double angular=circularAngularMomentum(r);
        const double energy=effective(r,angular);
        const double level=(attraction*reducedMass/hbar)
            *std::sqrt(reducedMass/(2.0*(-energy*reducedMass)));
        const double azimuthal=angular*reducedMass/hbar;
        const double coulomb=attraction/r;
        std::printf("%12.6f %10.3f %14.9f %14.9f %12.6f %12.3e %12.3e\n",
                    fraction,r/comptonBarrierRadius,level,azimuthal,
                    level-azimuthal,
                    std::abs(dipoleDipole(r)/coulomb),
                    std::abs(spinOrbit(r,angular)/coulomb));
    }
    std::printf("\nthe gap is what CREM_ACTION_TRACE prints as J_r/h on a "
                "circular orbit; the trace measured -0.0158 at the deep end\n");
}
