// Audit 216 / test B: an INDEPENDENT integral of T^{0nu} over a common
// t = const hypersurface.
//
// Why it is needed.  crem_trajectory.hpp:672 defines
//   boundFieldEnergy = (E_particle + E_rad + E_bound)_before
//                      - E_particle(now) - E_rad(now)
// and boundFieldMomentum/AngularMomentum the same way.  They are the
// RESIDUALS that close the balance, not measurements.  The reported
// conservation therefore cannot fail and carries no information about
// whether the residual is the field's energy.  Only an integral of the
// field's own stress-energy over a slice of constant time can test it.
//
// What is integrated: the CROSS (interaction) part, finite without any
// self-energy regulator,
//   u_x = eps0 E1.E2 + B1.B2/mu0          [T^{00}]
//   g_x = eps0 (E1 x B2 + E2 x B1)        [T^{0i}/c]
// with E_i, B_i the retarded fields the model itself produces at an
// arbitrary point.  In the static limit the energy has a closed form
// over a ball of radius R about the pair,
//   int_{|x|<R} u_x d3x = -K (1/d - 1/R),
// so the probe is checked against an exact value before it is asked to
// judge anything.  The O(beta^2) relativistic correction is 5.3e-05
// here and is below the quadrature error, so it is not modelled.
//
// Two constraints the geometry imposes, both of them findings:
//   1. CAUSAL HORIZON.  Evaluating a retarded field at distance R needs
//      history reaching back R/c.  Past that the solver extrapolates
//      and the integral diverges outright -- at s_out = 1e6 d an early
//      version returned +2.9e+02 J for a -8.7e-18 J quantity.  The slice
//      can only be integrated inside c * (history duration).
//   2. The integrand has a 1/s^2 singularity at each charge.  A single
//      grid about the centre of mass cannot resolve them; a Becke
//      partition of unity, one grid per charge, makes the measure
//      s^2 ds cancel each singularity exactly and needs no mask.
#include "modules/crem_trajectory.hpp"
#include "modules/crem_collapse.hpp"
#include <cstdio>
#include <cmath>
#include <array>
#include <vector>

int main(){
    std::setvbuf(stdout,nullptr,_IOLBF,0);
    const double mu=firstMass*secondMass/(firstMass+secondMass);
    const double k=pairCoulombStrength/mu;
    const double r0=0.25*pairBohrRadius({electron,positron});
    const OsculatingElements el{-k/(2.0*r0),std::sqrt(k*r0)};
    const double period=osculatingPeriod(el.specificEnergy,k);
    const Vec3 m1=Vec3{0.5,0.0,0.8660254037844386}*firstMagneticMoment;
    const State start=osculatingPeriapsisState(
        el,k,m1,m1*(secondMagneticMoment/firstMagneticMoment),
        Vec3{0,0,1},Vec3{1,0,0},0.0);

    // Half an orbit, so the causal horizon c*T clears the outer radius
    // by a wide margin, and so the instant is not the symmetric
    // preparation state.
    ClassicalTrajectoryEngine::Accuracy acc;
    ClassicalTrajectoryEngine engine(start,acc);
    State s=start;
    const int steps=400;
    const double span=0.5*period;
    for(int i=0;i<steps;++i)
        if(!engine.advance(s,span/steps)){
            std::printf("advance failed at step %d\n",i); return 1; }
    const StateHistory& hist=engine.history();
    const Vec3 sep=s.firstPosition-s.secondPosition;
    const double d=sep.norm();
    const double horizon=c*span/d;
    std::printf("# separation %.6e m = %.6f a_pair\n",
                d,d/pairBohrRadius({electron,positron}));
    std::printf("# beta %.3e, beta^2 %.3e (below the quadrature error)\n",
                std::max(s.firstVelocity.norm(),s.secondVelocity.norm())/c,
                std::pow(std::max(s.firstVelocity.norm(),
                                  s.secondVelocity.norm())/c,2));
    std::printf("# causal horizon c*T_history = %.4g d\n",horizon);

    const double eps0=1.0/(mu0*c*c);
    const auto partition=[&](const Vec3& x,bool first){
        const double s1=(x-s.firstPosition).squaredNorm();
        const double s2=(x-s.secondPosition).squaredNorm();
        const double p1=s2*s2,p2=s1*s1,sum=p1+p2;
        return sum>0.0?(first?p1:p2)/sum:0.5;
    };
    // A node whose cross energy density exceeds this is not physics: it
    // is a failed retarded-time solve.  The bound is the density of the
    // Coulomb field at one tenth of the inner radius, which no node of
    // the grid ever legitimately reaches.  Rejections are COUNTED and
    // reported, never silently dropped.
    const auto integrate=[&](int radialNodes,int angularNodes,
                             double innerFraction,double outerFraction,
                             double& energy,Vec3& mom,long& rejected){
        const std::span<const SphereQuadraturePoint> sphere=
            sphereQuadratureView(angularNodes);
        const double sin_=innerFraction*d,sout=outerFraction*d;
        const double h=std::log(sout/sin_)/radialNodes;
        const double densityBound=
            pairCoulombStrength/std::pow(0.1*sin_,4)/(4.0*pi*eps0);
        energy=0.0; mom=Vec3{}; rejected=0;
        for(int centre=0;centre<2;++centre){
            const Vec3 origin=centre==0?s.firstPosition:s.secondPosition;
            for(int i=0;i<=radialNodes;++i){
                const double r=sin_*std::exp(i*h);
                const double wr=(i==0||i==radialNodes)?1.0:((i%2)?4.0:2.0);
                const double radial=wr*(h/3.0)*r*r*r;
                for(const SphereQuadraturePoint& q:sphere){
                    const Vec3 x=origin+q.direction*r;
                    const ElectromagneticField f1=
                        fieldFromOtherParticleAt(x,s.time,s,hist,false);
                    const ElectromagneticField f2=
                        fieldFromOtherParticleAt(x,s.time,s,hist,true);
                    if(!isFinite(f1.electric)||!isFinite(f2.electric)
                       ||!isFinite(f1.magnetic)||!isFinite(f2.magnetic)){
                        ++rejected; continue; }
                    const double density=
                        eps0*dot(f1.electric,f2.electric)
                       +dot(f1.magnetic,f2.magnetic)/mu0;
                    if(!(std::abs(density)<densityBound)){
                        ++rejected; continue; }
                    const double w=radial*q.solidAngleWeight
                                  *partition(x,centre==0);
                    energy+=w*density;
                    mom+=(cross(f1.electric,f2.magnetic)
                         +cross(f2.electric,f1.magnetic))*(w*eps0);
                }
            }
        }
    };

    std::printf("\n%7s %8s %8s %8s %16s %16s %9s %8s\n",
                "radial","angular","s_in/d","s_out/d","int u_x [J]",
                "analytic [J]","measured/","rejected");
    std::printf("%7s %8s %8s %8s %16s %16s %9s %8s\n",
                "","","","","","-K(1/d-1/R)","analytic","nodes");
    for(const std::array<double,4>& cfg:std::vector<std::array<double,4>>{
            { 200, 50,1e-4,1e2},{ 400, 50,1e-4,1e2},{ 800, 50,1e-4,1e2},
            {1600, 50,1e-4,1e2},
            { 800, 26,1e-4,1e2},{ 800,194,1e-4,1e2},{ 800,302,1e-4,1e2},
            { 800,194,1e-6,1e2},{ 800,194,1e-4,1e1},{ 800,194,1e-4,1e3}}){
        double e=0.0; Vec3 g; long rej=0;
        integrate(static_cast<int>(cfg[0]),static_cast<int>(cfg[1]),
                  cfg[2],cfg[3],e,g,rej);
        const double exact=-pairCoulombStrength*(1.0/d-1.0/(cfg[3]*d));
        std::printf("%7.0f %8.0f %8.0e %8.0e %16.8e %16.8e %9.5f %8ld\n",
                    cfg[0],cfg[1],cfg[2],cfg[3],e,exact,e/exact,rej);
    }
    std::printf("\n# the model's own residual bookkeeping, same instant:\n");
    std::printf("#   boundFieldEnergy           %16.8e J\n",
                s.boundFieldEnergy);
    std::printf("#   conservativeParticleEnergy %16.8e J\n",
                conservativeParticleEnergy(s));
    std::printf("#   radiatedEnergy             %16.8e J\n",
                s.radiatedEnergy);
    std::printf("#   boundFieldMomentum         %16.8e kg m/s\n",
                s.boundFieldMomentum.norm());
}
