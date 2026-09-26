// PROBE TEMPLATE -- TWO COLUMNS, ALWAYS.
// ======================================================================
// Copy into audits/data/section-NNN/, rename, and replace measure()
// with whatever this probe is actually for.  Keep the shape: every
// quantity is reported BOTH at the preparation instant and averaged
// over a resolved orbit, side by side, in the same table.  The extra
// cost is one integration per measurement.
//
// Why the shape is mandatory (audits 190, 201-203, 211-212, 215):
//   190c  M1 at the preparation instant was wrong by 0.70x to 736x.
//   203b  the E1 split REVERSES SIGN between the two columns,
//         -9.41e-05 evolved against +9.59e-06 prepared, which withdrew
//         201b-201d outright.
//   212c  the dipole barrier reads +0.20035 instantaneous and -0.11107
//         averaged: the orbit averages it away entirely.
//   215f  211d's sign test failed only because at the radius it chose
//         0.01% of the effect exists yet.
// Four of the five retractions in the 180-215 stretch were this one
// mistake.  A ratio does NOT protect against it: the bias is
// configuration-dependent, so it does not cancel.
//
// Two further habits this template encodes:
//   - A NULL column.  Run one configuration that must give exactly the
//     reference value (here: para against itself).  It reads 0.000e+00
//     or the probe is wrong, and that is checked before anything else.
//   - Measure the quantity the PRODUCTION path consumes, not the one
//     described by the comment beside it.  Audit 215b lost a factor
//     2.06 by differentiating the wrong energy channel; the estimator
//     substitutes orbitalRadiatedEnergy at crem_collapse.hpp:3213 while
//     the lambda above it computes something else, for the angular
//     momentum.  Follow the value, not the prose.
//
// As shipped this template measures the dimensionless dipole-dipole
// coupling factor  mu1.mu2 - 3 (mu1.rhat)(mu2.rhat)  over unit moments,
// for parallel moments (para).  That case has a closed form in both
// columns, so running the template unmodified checks itself:
//     prepared (moment in the r-mu plane, azimuth 0):  1 - 3 sin^2 d
//     orbit-averaged over a circular orbit:            1 - (3/2) sin^2 d
// They differ by up to a factor of ten and cross zero at different
// angles -- which is the whole argument for the second column.
// ======================================================================
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>

namespace {

// REPLACE THIS.  Anything computable from a Frame works; use stepReady
// and a State if the quantity needs accelerations or proper dipoles.
double measure(const Vec3& separation,
               const Vec3& firstDipole,const Vec3& secondDipole) {
    const double r=separation.norm();
    const double a=firstDipole.norm(),b=secondDipole.norm();
    if(!(r>0.0&&a>0.0&&b>0.0)) return 0.0;
    const Vec3 rhat=separation*(1.0/r);
    return dot(firstDipole,secondDipole)/(a*b)
          -3.0*dot(firstDipole,rhat)*dot(secondDipole,rhat)/(a*b);
}

}  // namespace

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double r0=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements elements{-k/(2.0*r0),std::sqrt(k*r0)};
    const double period=osculatingPeriod(elements.specificEnergy,k);

    std::printf("# r = 0.25 a_pair, circular, one resolved orbit.\n");
    std::printf("%7s %12s %12s %10s %12s %12s %10s\n",
                "delta","prepared","evolved","evol/prep",
                "1-3sin^2","1-1.5sin^2","null");
    for(double deltaDeg:{0.0,15.0,30.0,45.0,54.7356,60.0,75.0,90.0}){
        const double d=deltaDeg*3.14159265358979324/180.0;
        const Vec3 unit{std::sin(d),0.0,std::cos(d)};
        const Vec3 mu1=unit*firstMagneticMoment;
        double column[2][2]={{0,0},{0,0}};   // [channel][prepared,evolved]
        bool ok=true;
        // Channel 1 duplicates channel 0 exactly: it is the null.
        for(int channel=0;channel<2;++channel){
            const Vec3 mu2=unit*secondMagneticMoment;
            const State start=osculatingPeriapsisState(
                elements,k,mu1,mu2,Vec3{0.0,0.0,1.0},Vec3{1.0,0.0,0.0},0.0);

            // COLUMN 1 -- the preparation instant.  On its own this is
            // the number that has been wrong five times; it is reported
            // only because column 2 is beside it.
            column[channel][0]=measure(
                start.firstPosition-start.secondPosition,
                start.firstDipole,start.secondDipole);

            // COLUMN 2 -- the same quantity over a resolved orbit.
            // Frames are sampled at uniform physical time, so the plain
            // mean over them is the time average.
            double sum=0.0; long count=0;
            SimulationOptions options;
            options.frameCount=512;
            options.radiatedEnergyBookkeeping=false;
            options.frameReady=[&](const Frame& f){
                sum+=measure(f.first-f.second,f.firstDipole,f.secondDipole);
                ++count;
            };
            const MechanicalTrajectoryResult run=runMechanicalTrajectory(
                start,period,nuclearCutoff,options,gRadiationReactionModel);
            if(!isFinite(run.finalState)||count==0){ ok=false; break; }
            column[channel][1]=sum/static_cast<double>(count);
        }
        if(!ok){ std::printf("%7.2f failed\n",deltaDeg); continue; }
        const double s=std::sin(d);
        std::printf("%7.2f %12.5f %12.5f %10.3f %12.5f %12.5f %10.3e\n",
            deltaDeg,column[0][0],column[0][1],
            column[0][0]!=0.0?column[0][1]/column[0][0]:0.0,
            1.0-3.0*s*s,1.0-1.5*s*s,
            column[0][1]!=0.0
                ?(column[1][1]-column[0][1])/column[0][1]:0.0);
    }
    std::printf("# SELF-CHECK.  The null column must read 0.000e+00\n");
    std::printf("# throughout.  'prepared' must reproduce 1-3sin^2 to\n");
    std::printf("# five decimals and 'evolved' 1-1.5sin^2 to about\n");
    std::printf("# 1e-03 -- that residual is the orbit not closing\n");
    std::printf("# exactly on itself plus finite frame sampling, not an\n");
    std::printf("# error.  The point of the table is the third column:\n");
    std::printf("# at delta = 45 deg the two disagree in SIGN, -0.500\n");
    std::printf("# against +0.249, and at 54.7356 deg by a factor 500.\n");
    std::printf("# One column would have reported either of those as\n");
    std::printf("# the answer.\n");
}
