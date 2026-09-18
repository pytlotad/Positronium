// Can the classical model carry an S-state analogue: L = 0, sampling r -> 0
// isotropically, with a finite fall time (audit section 118)?
//
// Three parts.  (1) The tensor structure 3(mu1.n)(mu2.n) - mu1.mu2 averaged
// over an ISOTROPIC ensemble of orbit directions, which is what a spherical
// state does and an azimuth average around one normal does not.  (2) The
// isotropic part of the model's own coupling, time-averaged over a RADIAL
// Kepler orbit of the same semi-major axis, against the 116.8 GHz contact
// term -- plus the moment radius that would make the model's contact energy
// at r -> 0 equal the quantum Fermi term.  (3) The free-fall time of that
// radial orbit, which is what "finite fall time" means here.
//
// Usage: radial_state_probe [semi-major axis in a_pair, default 1]
//        [isotropic samples, default 100000]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/radial_state_probe.cpp -o /tmp/radial $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

int main(int argc,char** argv) {
    const double axis=(argc>1?atof(argv[1]):1.0)*pairBohrRadius(activePair);
    const int samples=argc>2?atoi(argv[2]):100000;
    const double planck=2.0*pi*hbar;
    const double moment=firstMagneticMoment;
    const double softening=magneticDipoleRadius();
    // (1) Isotropic average of the tensor structure, in units of mu1.mu2.
    // Deterministic quadrature, not sampling: the average is exactly zero
    // analytically and a Monte Carlo draw only reproduces it to 1/sqrt(N).
    const Vec3 first{0.0,0.0,1.0},second{0.0,0.0,1.0};
    double tensorSum=0.0,tensorWeight=0.0;
    const int polar=2001,azimuthal=64;
    for(int index=0;index<polar;++index) {
        const double z=-1.0+2.0*index/(polar-1.0);
        const double radial=std::sqrt(std::max(0.0,1.0-z*z));
        const double nodeWeight=(index==0||index==polar-1)?0.5:1.0;
        for(int slice=0;slice<azimuthal;++slice) {
            const double phi=2.0*pi*slice/azimuthal;
            const Vec3 direction{radial*std::cos(phi),radial*std::sin(phi),z};
            tensorSum+=nodeWeight*(3.0*dot(first,direction)
                *dot(second,direction)-dot(first,second));
            tensorWeight+=nodeWeight;
        }
    }
    std::printf("P1 isotropic average of 3(mu1.n)(mu2.n)-mu1.mu2 over a "
                "%dx%d quadrature: %.3e of mu1.mu2 (exactly 0 analytically; "
                "%d random draws give %.1e noise)\n",
                polar,azimuthal,tensorSum/tensorWeight,samples,
                1.0/std::sqrt(static_cast<double>(samples)));
    // (2) The isotropic coefficient, time-averaged over a radial orbit.
    // Radial Kepler motion: r = a(1-cos E)/2, t = (T/2pi)(E - sin E)/... use
    // the standard parametrisation r = a sin^2(eta), t ~ eta - sin(2eta)/2.
    const auto isotropicGHz=[&](double separation) {
        // The mu1.mu2 coefficient of the model's own energy, read off by
        // putting both moments perpendicular to the separation.
        const Vec3 along{1.0,0.0,0.0},perpendicular{0.0,0.0,1.0};
        const Vec3 displacement=along*separation;
        const double energy=pairDipoleInteractionEnergy(displacement,
            perpendicular*moment,perpendicular*moment);
        return std::abs(2.0*energy)/planck*1.0e-9;
    };
    const int nodes=200001;
    double weighted=0.0,weight=0.0,deepest=axis;
    for(int index=1;index<nodes;++index) {
        const double eta=pi*index/(nodes-1.0);
        const double radius=axis*std::sin(eta)*std::sin(eta);
        // dt/deta ~ sin^2(eta) for radial Kepler motion.
        const double timeWeight=std::sin(eta)*std::sin(eta);
        weighted+=isotropicGHz(radius)*timeWeight; weight+=timeWeight;
        deepest=std::min(deepest,radius);
    }
    std::printf("P2 isotropic coupling, time-averaged over the radial orbit: "
                "%.1f GHz (target 116.8); deepest sampled %.3e m = %.4f r*\n",
                weighted/weight,deepest,deepest/comptonBarrierRadius);
    std::printf("   the same coupling at r = 0: %.4e GHz; at a_pair: %.3f GHz\n",
                isotropicGHz(0.0),isotropicGHz(axis));
    // The moment radius that would put the model's contact energy at r -> 0
    // on the quantum Fermi term.
    double low=softening,high=1.0e3*softening;
    const double target=116.8;
    for(int iteration=0;iteration<200;++iteration) {
        const double middle=0.5*(low+high);
        const double contact=(mu0/(4.0*pi))*3.0*moment*moment
            /(middle*middle*middle);
        const double ghz=std::abs(2.0*contact)/planck*1.0e-9;
        if(ghz>target) low=middle; else high=middle;
    }
    const double needed=0.5*(low+high);
    std::printf("   moment radius that matches 116.8 GHz at contact: %.4e m "
                "= %.1f r* = %.4f a_pair (the model uses %.4f r*)\n",
                needed,needed/comptonBarrierRadius,
                needed/pairBohrRadius(activePair),
                softening/comptonBarrierRadius);
    // (3) Free-fall time of the radial orbit, and the circular period.
    const double period=2.0*pi*std::sqrt(pairReducedMass*axis*axis*axis
        /pairCoulombStrength);
    std::printf("P3 radial orbit at a = %.4f a_pair: circular period %.4e s, "
                "free fall to the core %.4e s = %.4f fs\n",
                axis/pairBohrRadius(activePair),period,period/(4.0*std::sqrt(2.0)),
                period/(4.0*std::sqrt(2.0))*1.0e15);
    std::printf("   for comparison the measured collapse time is 30.67 ps, "
                "i.e. %.3e times longer\n",
                30.67e-12/(period/(4.0*std::sqrt(2.0))));
    // Does such a state survive its own core passage?  Integrate the model's
    // own E1 power along one radial passage, with the Coulomb force softened
    // at the same separationFloor the trajectory uses.
    {
        const double floor=separationFloor();
        const double binding=pairCoulombStrength/(2.0*axis);
        const double coefficient=pairDipoleCharge*pairDipoleCharge
            /(6.0*pi*epsilon0*c*c*c);
        double radiated=0.0;
        const int steps=2000001;
        for(int index=1;index<steps;++index) {
            const double eta=pi*index/(steps-1.0);
            const double radius=axis*std::sin(eta)*std::sin(eta);
            const double softened=std::sqrt(radius*radius+floor*floor);
            // Relative acceleration of the softened Coulomb pair, and the E1
            // power it radiates: P = (q^2/6 pi eps0 c^3) |a_rel|^2 for a pair
            // whose dipole is q times the relative coordinate.
            const double acceleration=pairCoulombStrength
                /(pairReducedMass*softened*softened);
            const double power=coefficient*acceleration*acceleration;
            const double timeStep=(pi/(steps-1.0))*2.0*axis
                *std::sin(eta)*std::sin(eta)
                /std::sqrt(pairCoulombStrength/(pairReducedMass*axis));
            radiated+=power*timeStep;
        }
        std::printf("P3 radiated in ONE radial passage: %.4e J = %.4e eV "
                    "against a binding of %.4f eV, i.e. %.3e of it\n",
                    radiated,radiated/1.602176634e-19,
                    binding/1.602176634e-19,radiated/binding);
        std::printf("   (Coulomb softened at the model's own separation floor "
                    "%.4f r*)\n",floor/comptonBarrierRadius);
    }
}
