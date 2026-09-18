// How large is the model's own spin-spin coupling, and how does it depend on
// the orbit's angular momentum (audit section 116)?
//
// The measured o-Ps/p-Ps splitting is 203.3941 GHz; in QED 4/7 of it is the
// Fermi contact term (116.8 GHz) and 3/7 virtual annihilation (87.6 GHz),
// which this model does not contain at all.  What it does contain is the
// classical dipole-dipole energy WITH its own contact (magnetization) term,
// azimuthAveragedDipoleEnergy.  The splitting analogue is the difference
// between the two collinear configurations, U(aligned) - U(anti-aligned),
// which is 2U(aligned) because U is odd in the mutual cosine.
//
// The scan is over the angular momentum at fixed semi-major axis: a circular
// L = hbar orbit never approaches the origin, while the physical ground state
// is an S state whose density at contact is what produces the Fermi term.  If
// the model's coupling climbs toward the contact value as L falls, the
// missing factor is the orbit's angular momentum, not the coupling's form.
//
// Usage: spin_spin_coupling_scan [semi-major axis in a_pair, default 1]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/spin_spin_coupling_scan.cpp -o /tmp/scan $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <utility>

int main(int argc,char** argv) {
    const double axis=(argc>1?atof(argv[1]):1.0)*pairBohrRadius(activePair);
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 aligned=normal*firstMagneticMoment;
    const Vec3 partner=normal*secondMagneticMoment;
    const double planck=2.0*pi*hbar;
    // The azimuth-averaged energy is a(r) (mu1.mu2) + b(r) (mu1.n)(mu2.n).
    // The first commutes with S^2 (it is the model's contact/isotropic part),
    // the second does not (audit 117).  Both are read off the model's own
    // function by evaluating it on two moment configurations rather than by
    // copying its formula: aligned along the normal gives a+b, and one moment
    // along the normal with the other in the plane gives a alone.
    const Vec3 inPlane{1.0,0.0,0.0};
    const auto coefficients=[&](double separation) {
        const Vec3 partnerInPlane=inPlane*secondMagneticMoment;
        const double bothNormal=azimuthAveragedDipoleEnergy(
            separation,aligned,partner,normal);
        const double mixed=azimuthAveragedDipoleEnergy(
            separation,aligned,partnerInPlane,normal);
        // mixed has mu1.mu2 = 0 and (mu1.n)(mu2.n) = 0, so it vanishes; the
        // isotropic coefficient comes from the pair with both moments in the
        // plane and parallel: there mu1.mu2 = mu^2 and (mu1.n)(mu2.n) = 0.
        const double bothInPlane=azimuthAveragedDipoleEnergy(
            separation,inPlane*firstMagneticMoment,partnerInPlane,normal);
        (void)mixed;
        const double isotropic=bothInPlane;          // a mu^2
        const double tensor=bothNormal-bothInPlane;  // b mu^2
        return std::pair<double,double>{isotropic,tensor};
    };
    const auto isotropicGHz=[&](double separation) {
        return std::abs(2.0*coefficients(separation).first)/planck*1.0e-9;
    };
    const auto tensorGHz=[&](double separation) {
        return std::abs(2.0*coefficients(separation).second)/planck*1.0e-9;
    };
    const auto splittingGHz=[&](double separation) {
        const double up=azimuthAveragedDipoleEnergy(
            separation,aligned,partner,normal);
        const double down=azimuthAveragedDipoleEnergy(
            separation,aligned,partner*(-1.0),normal);
        return std::abs(up-down)/planck*1.0e-9;
    };
    // The bare point-dipole tensor term the report quotes today: no contact
    // term, no softening, and the separation direction in the orbital plane,
    // which is where the pair actually is.
    const Vec3 separationDirection{1.0,0.0,0.0};
    const auto barePointGHz=[&](double separation) {
        const double moments=dot(aligned,partner);
        const double along=dot(aligned,separationDirection)
            *dot(partner,separationDirection);
        const double up=(mu0/(4.0*pi))*(moments-3.0*along)
            /(separation*separation*separation);
        return std::abs(2.0*up)/planck*1.0e-9;
    };
    std::printf("semi-major axis %.4f a_pair = %.4e m; moment radius %.4f r*\n",
                axis/pairBohrRadius(activePair),axis,
                magneticDipoleRadius()/comptonBarrierRadius);
    std::printf("measured splitting 203.3941 GHz = contact 116.8 + virtual "
                "annihilation 87.6 (the latter has no classical analogue)\n\n");
    std::printf("%8s %12s %14s %16s %16s %14s %14s %10s\n","L/hbar",
                "periapsis","r_p/r*","at periapsis","orbit-averaged",
                "isotropic","tensor","|b/a|");
    for(double angular:{1.0,0.7,0.5,0.3,0.2,0.1,0.05,0.03,0.01}) {
        const double eccentricity=
            std::sqrt(std::max(0.0,1.0-angular*angular));
        const double periapsis=axis*(1.0-eccentricity);
        // Orbit average over the eccentric anomaly: dt = (T/2pi)(1-e cosE) dE
        // and r = a(1 - e cos E).
        const int samples=20001;
        double sum=0.0,weight=0.0;
        for(int index=0;index<samples;++index) {
            const double anomaly=2.0*pi*index/(samples-1.0);
            const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
            const double timeWeight=1.0-eccentricity*std::cos(anomaly);
            sum+=splittingGHz(radius)*timeWeight;
            weight+=timeWeight;
        }
        double isotropicSum=0.0,tensorSum=0.0;
        for(int index=0;index<samples;++index) {
            const double anomaly=2.0*pi*index/(samples-1.0);
            const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
            const double timeWeight=1.0-eccentricity*std::cos(anomaly);
            isotropicSum+=isotropicGHz(radius)*timeWeight;
            tensorSum+=tensorGHz(radius)*timeWeight;
        }
        std::printf("%8.2f %12.4e %14.4f %16.3f %16.3f %14.3f %14.3f %10.3f\n",
                    angular,periapsis,periapsis/comptonBarrierRadius,
                    splittingGHz(periapsis),sum/weight,
                    isotropicSum/weight,tensorSum/weight,
                    (tensorSum/weight)/std::max(isotropicSum/weight,1.0e-300));
    }
    std::printf("\nat L = hbar (the production orbit): model with its contact "
                "term %.3f GHz, bare point dipole %.3f GHz\n",
                splittingGHz(axis),barePointGHz(axis));
    std::printf("the report's dipoleCouplingHz uses the bare form; at this "
                "radius the two agree to %.4f\n",
                splittingGHz(axis)/std::max(barePointGHz(axis),1.0e-300));
    // Where the model's orbit-averaged coupling equals the quantum numbers.
    const auto averaged=[&](double angular) {
        const double eccentricity=
            std::sqrt(std::max(0.0,1.0-angular*angular));
        const int samples=20001;
        double sum=0.0,weight=0.0;
        for(int index=0;index<samples;++index) {
            const double anomaly=2.0*pi*index/(samples-1.0);
            const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
            const double timeWeight=1.0-eccentricity*std::cos(anomaly);
            sum+=splittingGHz(radius)*timeWeight;
            weight+=timeWeight;
        }
        return sum/weight;
    };
    // Where the ISOTROPIC part alone reaches the contact value, and what the
    // tensor part is doing there: the question of audit 117.
    {
        const auto averagedIsotropic=[&](double angular) {
            const double eccentricity=
                std::sqrt(std::max(0.0,1.0-angular*angular));
            const int samples=20001;
            double sum=0.0,weight=0.0;
            for(int index=0;index<samples;++index) {
                const double anomaly=2.0*pi*index/(samples-1.0);
                const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
                const double timeWeight=1.0-eccentricity*std::cos(anomaly);
                sum+=isotropicGHz(radius)*timeWeight; weight+=timeWeight;
            }
            return sum/weight;
        };
        const auto averagedTensor=[&](double angular) {
            const double eccentricity=
                std::sqrt(std::max(0.0,1.0-angular*angular));
            const int samples=20001;
            double sum=0.0,weight=0.0;
            for(int index=0;index<samples;++index) {
                const double anomaly=2.0*pi*index/(samples-1.0);
                const double radius=axis*(1.0-eccentricity*std::cos(anomaly));
                const double timeWeight=1.0-eccentricity*std::cos(anomaly);
                sum+=tensorGHz(radius)*timeWeight; weight+=timeWeight;
            }
            return sum/weight;
        };
        double low=0.005,high=1.0;
        for(int iteration=0;iteration<60;++iteration) {
            const double middle=0.5*(low+high);
            if(averagedIsotropic(middle)>116.8) low=middle; else high=middle;
        }
        const double match=0.5*(low+high);
        std::printf("\nisotropic part = 116.8 GHz at L = %.4f hbar; there the "
                    "tensor part is %.1f GHz, ratio %.2f\n",match,
                    averagedTensor(match),
                    averagedTensor(match)/std::max(averagedIsotropic(match),
                                                   1.0e-300));
        std::printf("ratio tensor/isotropic at L = 1.00, 0.30, 0.10, 0.03: "
                    "%.2f %.2f %.2f %.2f\n",
                    averagedTensor(1.0)/averagedIsotropic(1.0),
                    averagedTensor(0.3)/averagedIsotropic(0.3),
                    averagedTensor(0.1)/averagedIsotropic(0.1),
                    averagedTensor(0.03)/averagedIsotropic(0.03));
    }
    for(double target:{116.8,203.3941}) {
        double low=0.02,high=1.0;
        for(int iteration=0;iteration<60;++iteration) {
            const double middle=0.5*(low+high);
            if(averaged(middle)>target) low=middle; else high=middle;
        }
        std::printf("the model's orbit-averaged coupling equals %.1f GHz at "
                    "L = %.4f hbar\n",target,0.5*(low+high));
    }
}
