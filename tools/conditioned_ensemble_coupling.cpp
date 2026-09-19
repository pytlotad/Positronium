// Sections 117, 118 and the M1 result, restated as ENSEMBLE averages over
// what the default run actually prepares (audit section 153).
//
// Those three are statements about the preparation audit 91 withdrew, the
// moments exactly collinear.  Audit 143 gave the two channels back by
// CONDITIONING -- menu 1 keeps only draws the model's own classifier calls
// para, menu 2 only ortho -- so the same quantities can be averaged over the
// sub-ensembles the default now produces.
//
// What is averaged, per channel, over prepared states only:
//   <cos>, <cos^2>                         the mutual angle's two moments
//   <|m|/mu>, <|m|^2/mu^2>                 the second is what M1 needs, since
//                                          P_M1 ~ |m|^2
//   <mu1.mu2>, <(mu1.n)(mu2.n)>            the two structures of the
//                                          azimuth-averaged energy, which
//                                          audit 117b split as
//                                          a (mu1.mu2) + b (mu1.n)(mu2.n)
// a and b are read off the model's own azimuthAveragedDipoleEnergy on two
// calibration configurations, exactly as 117b read them, so no formula is
// copied here either.
//
// Usage: conditioned_ensemble_coupling [count, default 4000]
//        [seed base, default 42] [semi-major axis in a_pair, default 1]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/conditioned_ensemble_coupling.cpp -o /tmp/ensemble $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc,char** argv) {
    const int count=argc>1?atoi(argv[1]):4000;
    const unsigned long long base=argc>2?strtoull(argv[2],nullptr,10):42ULL;
    const double axis=(argc>3?atof(argv[3]):1.0)*pairBohrRadius(activePair);
    const double planck=2.0*pi*hbar;
    const Vec3 normal{0.0,0.0,1.0};
    const Vec3 inPlane{1.0,0.0,0.0};
    // a and b from the model's own function, as audit 117b took them: both
    // moments in the orbital plane isolates a, both along the normal gives
    // a + b.  Divided by mu^2 so they multiply the structures below directly.
    const double momentProduct=firstMagneticMoment*secondMagneticMoment;
    const double bothInPlane=azimuthAveragedDipoleEnergy(axis,
        inPlane*firstMagneticMoment,inPlane*secondMagneticMoment,normal);
    const double bothNormal=azimuthAveragedDipoleEnergy(axis,
        normal*firstMagneticMoment,normal*secondMagneticMoment,normal);
    const double isotropicCoefficient=bothInPlane/momentProduct;
    const double tensorCoefficient=(bothNormal-bothInPlane)/momentProduct;
    std::printf("a = %.6e J/mu^2, b = %.6e J/mu^2, |b/a| = %.4f  "
                "(a_pair x %.4f)\n",isotropicCoefficient,tensorCoefficient,
                std::abs(tensorCoefficient/isotropicCoefficient),
                axis/pairBohrRadius(activePair));
    std::printf("%8s %7s %10s %10s %11s %12s %12s %12s %12s\n","channel",
                "states","<cos>","<cos^2>","<|m|/mu>","<|m|^2/mu^2>",
                "<iso> [GHz]","<ten> [GHz]","<U> [GHz]");
    struct Row { double cosine,cosineSquared,coherent,coherentSquared,
                 isotropic,tensor; int states; };
    std::vector<Row> rows;
    for(int phenomenon:{1,2,5}) {
        Row row{};
        for(int index=0;index<count;++index) {
            SimulationOptions options;
            options.frameCount=2;
            options.observationTime=1.0e-24;
            const SimulationResult prepared=simulate(
                splitMix64(base+static_cast<unsigned long long>(index)),
                phenomenon,options);
            if(prepared.frames.empty()) continue;
            const Vec3 first=prepared.frames.front().firstDipole;
            const Vec3 second=prepared.frames.front().secondDipole;
            const double norms=first.norm()*second.norm();
            if(!(norms>0.0)) continue;
            const double cosine=dot(first,second)/norms;
            const double coherent=(first+second).norm()/firstMagneticMoment;
            row.cosine+=cosine;
            row.cosineSquared+=cosine*cosine;
            row.coherent+=coherent;
            row.coherentSquared+=coherent*coherent;
            // The two structures, in units of mu^2, with the orbital normal
            // as the model's own preparation uses it.
            row.isotropic+=dot(first,second)/momentProduct;
            row.tensor+=dot(first,normal)*dot(second,normal)/momentProduct;
            ++row.states;
        }
        if(row.states==0) continue;
        const double n=static_cast<double>(row.states);
        const double isotropicEnergy=
            isotropicCoefficient*momentProduct*(row.isotropic/n);
        const double tensorEnergy=
            tensorCoefficient*momentProduct*(row.tensor/n);
        std::printf("%8d %7d %10.4f %10.4f %11.4f %12.4f %12.4f %12.4f "
                    "%12.4f\n",phenomenon,row.states,row.cosine/n,
                    row.cosineSquared/n,row.coherent/n,row.coherentSquared/n,
                    isotropicEnergy/planck*1.0e-9,tensorEnergy/planck*1.0e-9,
                    (isotropicEnergy+tensorEnergy)/planck*1.0e-9);
        rows.push_back(Row{row.cosine/n,row.cosineSquared/n,row.coherent/n,
                           row.coherentSquared/n,isotropicEnergy,tensorEnergy,
                           row.states});
    }
    if(rows.size()>=2) {
        const Row& para=rows[0];
        const Row& ortho=rows[1];
        std::printf("\nM1 power ratio para/ortho = <|m|^2>_para / "
                    "<|m|^2>_ortho = %.4f\n",
                    para.coherentSquared/std::max(ortho.coherentSquared,
                                                  1.0e-300));
        const double separation=
            (para.isotropic+para.tensor-ortho.isotropic-ortho.tensor)
            /planck*1.0e-9;
        std::printf("channel separation in spin-spin energy = %.4f GHz "
                    "against the measured 203.3941, a factor %.2f\n",
                    std::abs(separation),
                    203.3941/std::max(std::abs(separation),1.0e-300));
    }
}
