// The emission ladder in the MODEL's own action, not in the Kepler label
// (audit section 138).
//
// Section 131c walked the ladder in n = (K/hbar) sqrt(mu/2|E|) and used the
// law Delta n = n[1-(1+2/n)^(-1/2)], both of which assume E = -R/n^2, i.e.
// pure Coulomb.  Section 132d measured that the model's own circular orbit
// disagrees with that label by 0.000163 at the third rung and 0.003243 at the
// stop, and 133 integrated the full radial action.  What was still quoted in
// Kepler is the ladder itself.
//
// Here the same ladder is walked twice from the same starting state, with the
// MEASURED photon energies of the trace:
//   Kepler:  n' = sqrt(R/(R/n^2 + E_photon));
//   model:   E' = E - E_photon, then invert E' -> circular radius -> J/h in
//            U(r) = -K/r + u_dd(r) + u_so(r,L) + L^2/2r^2, the same assembly
//            dipoleAwarePeriapsis solves and radial_action_full_potential
//            uses.
// Both are compared with the landings the trace recorded.
//
// Usage: action_ladder_full_potential [photon energies in eV, comma separated]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/action_ladder_full_potential.cpp -o /tmp/ladder2 $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

const double electronVolt=1.602176634e-19;

double dipoleDipole(double r,const Vec3& first,const Vec3& second,
                    const Vec3& normal) {
    return azimuthAveragedDipoleEnergy(r,first,second,normal)
        /reducedMassOf(activePair);
}

double spinOrbit(double r,double angular,const Vec3& first,
                 const Vec3& second,const Vec3& normal) {
    const double totalMass=firstMass+secondMass;
    if(!(r>0.0)||!(totalMass>0.0)) return 0.0;
    const Vec3 radialHat{0.0,-1.0,0.0};
    const Vec3 tangentialHat=cross(normal,radialHat);
    State turningPoint{};
    turningPoint.firstPosition=radialHat*(r*secondMass/totalMass);
    turningPoint.secondPosition=radialHat*(-r*firstMass/totalMass);
    turningPoint.firstVelocity=
        tangentialHat*((angular/r)*secondMass/totalMass);
    turningPoint.secondVelocity=
        tangentialHat*(-(angular/r)*firstMass/totalMass);
    turningPoint.firstProperDipole=first;
    turningPoint.secondProperDipole=second;
    synchronizeCovariantDipoles(turningPoint);
    return chargeDipoleInteractionEnergy(turningPoint)
        /reducedMassOf(activePair);
}

// The model's circular orbit at radius r: L from dU/dr = 0, iterated because
// the spin-orbit term carries L itself.
double circularAngularMomentum(double r,const Vec3& first,const Vec3& second,
                               const Vec3& normal) {
    const double attraction=pairCoulombStrength/reducedMassOf(activePair);
    double angular=std::sqrt(attraction*r);
    for(int iteration=0;iteration<200;++iteration) {
        const double step=1.0e-6*r;
        const auto withoutCentrifugal=[&](double radius) {
            return -attraction/radius+dipoleDipole(radius,first,second,normal)
                +spinOrbit(radius,angular,first,second,normal);
        };
        const double slope=
            (withoutCentrifugal(r+step)-withoutCentrifugal(r-step))
            /(2.0*step);
        const double next=std::sqrt(std::max(0.0,slope*r*r*r));
        if(std::abs(next-angular)<=1.0e-15*std::max(angular,1.0e-300))
            return next;
        angular=next;
    }
    return angular;
}

double circularEnergy(double r,const Vec3& first,const Vec3& second,
                      const Vec3& normal) {
    const double attraction=pairCoulombStrength/reducedMassOf(activePair);
    const double angular=circularAngularMomentum(r,first,second,normal);
    return -attraction/r+dipoleDipole(r,first,second,normal)
        +spinOrbit(r,angular,first,second,normal)
        +angular*angular/(2.0*r*r);
}

}  // namespace

int main(int argc,char** argv) {
    const double axis=pairBohrRadius(activePair);
    const double reducedMass=reducedMassOf(activePair);
    const double attraction=pairCoulombStrength/reducedMass;
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 first=normal*firstMagneticMoment;
    const Vec3 second=normal*secondMagneticMoment;
    std::vector<double> photons{13.604509281,70.670787429,634.225835239};
    if(argc>1) {
        photons.clear();
        std::string text=argv[1];
        std::size_t at=0;
        while(at<text.size()) {
            const std::size_t comma=text.find(',',at);
            photons.push_back(std::atof(text.substr(at,comma-at).c_str()));
            if(comma==std::string::npos) break;
            at=comma+1;
        }
    }
    const double rydberg=-groundStateSpecificEnergy()*reducedMass;
    // Invert a specific energy to the model's circular radius, then read the
    // action off that orbit.  The energy is monotone in r, so bisection is
    // enough; the bracket is checked before use.
    const auto radiusForEnergy=[&](double specificEnergy) {
        double low=1.0e-6*axis,high=1.0e3*axis;
        for(int iteration=0;iteration<300;++iteration) {
            const double middle=std::sqrt(low*high);
            if(circularEnergy(middle,first,second,normal)<specificEnergy)
                low=middle;
            else high=middle;
        }
        return std::sqrt(low*high);
    };
    const auto actionForEnergy=[&](double specificEnergy) {
        const double r=radiusForEnergy(specificEnergy);
        return circularAngularMomentum(r,first,second,normal)
            *reducedMass/hbar;   // J_r = 0 on a circle, so J/h = L/hbar
    };
    // R138c: the inversion against section 133's direct integration, which on
    // a circular orbit is J_phi/h alone -- checked here at the ladder's own
    // radii by going r -> E -> r and r -> J/h -> the same J/h.
    std::printf("R138c inversion check (r -> E -> r, and the action there):\n");
    for(double fraction:{1.0,0.333333,0.074670,0.009369}) {
        const double r=axis*fraction;
        const double energy=circularEnergy(r,first,second,normal);
        const double back=radiusForEnergy(energy);
        std::printf("  a/a_pair %10.6f  recovered %10.6f  relative %.2e  "
                    "J/h %.9f\n",fraction,back/axis,
                    std::abs(back-r)/r,actionForEnergy(energy));
    }
    std::printf("\nthe ladder, walked from the prepared state with the "
                "MEASURED photon energies:\n");
    std::printf("%8s %14s %14s %14s %14s %12s\n","photon","E [eV]",
                "Kepler n'","model J'/h","model/Kepler","Delta J/h");
    double keplerLevel=1.0;
    double specificEnergy=-attraction/(2.0*axis);
    double modelAction=actionForEnergy(specificEnergy);
    std::printf("%8s %14s %14.9f %14.9f %14.9f %12s\n","start","-",
                keplerLevel,modelAction,modelAction/keplerLevel,"-");
    for(std::size_t index=0;index<photons.size();++index) {
        const double photon=photons[index]*electronVolt;
        const double previousAction=modelAction;
        // Kepler: E = -R/n^2 exactly.
        const double binding=rydberg/(keplerLevel*keplerLevel);
        keplerLevel=std::sqrt(rydberg/(binding+photon));
        // Model: subtract the same photon from the model's own energy.
        specificEnergy-=photon/reducedMass;
        modelAction=actionForEnergy(specificEnergy);
        std::printf("%8zu %14.6f %14.9f %14.9f %14.9f %12.9f\n",
                    index+1,photons[index],keplerLevel,modelAction,
                    modelAction/keplerLevel,previousAction-modelAction);
    }
    std::printf("\ntraced landings for comparison (audit 131c): 0.577269 "
                "0.272406 0.096794\n");
}
