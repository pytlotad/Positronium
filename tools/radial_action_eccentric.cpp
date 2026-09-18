// The full radial action of the model's own potential (audit 133).
//
// Sections 131 and 132 used the Kepler closed form
// J_r = 2 pi (K sqrt(mu/2|E|) - L), which goes NEGATIVE along the collapse --
// impossible for a real radial action.  132 traced that to the dipole term but
// assumed the orbit is circular in the model's potential, which is the same
// statement as the gap being small.  This tool drops the assumption and
// integrates
//
//   J_r = 2 integral_{r1}^{r2} p_r dr,   p_r = mu sqrt(2 h(r)),
//   h(r) = eps + K/r - u_dd(r) - u_so(r,L) - L^2/(2 r^2),
//
// between the model's own turning points, h being the same radial function
// dipoleAwarePeriapsis solves.  The two energy functions are called, not
// copied; the turning-point state construction and the radial frame are
// mirrored from it, including its spin-orbit term, which is exact at the
// turning points and is used as an r-dependent potential in between exactly
// as the estimator uses it.
//
// The sqrt vanishes at both ends, so the integral is taken in
// r = c - b cos(theta) with c, b the midpoint and half-width of the interval:
// there the integrand is smooth, and Gauss-Legendre in theta converges.
//
// Usage: grep ^ACTION trace.err | radial_action_eccentric [nodes, default 400]
//        radial_action_eccentric --selftest    (rule R133a)
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/radial_action_eccentric.cpp -o /tmp/ecc $(root-config --libs)
#include "modules/crem_collapse.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct RadialProblem {
    double specificEnergy=0.0;
    double specificAngularMomentum=0.0;
    Vec3 firstDipole{};
    Vec3 secondDipole{};
    Vec3 normal{0.0,0.0,1.0};
};

// h(r) exactly as dipoleAwarePeriapsis assembles it.
double radialFunction(const RadialProblem& problem,double r) {
    if(!(r>0.0)) return -std::numeric_limits<double>::infinity();
    const double reducedMass=reducedMassOf(activePair);
    const double attraction=pairCoulombStrength/reducedMass;
    const double L=problem.specificAngularMomentum;
    const Vec3 normalHat=problem.normal*(1.0/std::max(problem.normal.norm(),
                                                      1.0e-300));
    const double dipole=azimuthAveragedDipoleEnergy(
        r,problem.firstDipole,problem.secondDipole,normalHat)/reducedMass;
    const Vec3 seedAxis=std::abs(normalHat.x)<0.9
        ?Vec3{1.0,0.0,0.0}:Vec3{0.0,1.0,0.0};
    Vec3 radialHat=cross(seedAxis,normalHat);
    const double radialHatNorm=radialHat.norm();
    if(radialHatNorm>0.0) radialHat=radialHat*(1.0/radialHatNorm);
    const Vec3 tangentialHat=cross(normalHat,radialHat);
    const double totalMass=firstMass+secondMass;
    double spinOrbit=0.0;
    if(firstMass>0.0&&totalMass>0.0&&radialHatNorm>0.0) {
        const double tangentialSpeed=L/r;
        State turningPoint{};
        turningPoint.firstPosition=radialHat*(r*secondMass/totalMass);
        turningPoint.secondPosition=radialHat*(-r*firstMass/totalMass);
        turningPoint.firstVelocity=
            tangentialHat*(tangentialSpeed*secondMass/totalMass);
        turningPoint.secondVelocity=
            tangentialHat*(-tangentialSpeed*firstMass/totalMass);
        turningPoint.firstProperDipole=problem.firstDipole;
        turningPoint.secondProperDipole=problem.secondDipole;
        synchronizeCovariantDipoles(turningPoint);
        spinOrbit=chargeDipoleInteractionEnergy(turningPoint)/reducedMass;
    }
    return problem.specificEnergy+attraction/r-dipole-spinOrbit
        -L*L/(2.0*r*r);
}

struct TurningPoints { double inner=0.0,outer=0.0; bool accessible=false; };

TurningPoints turningPoints(const RadialProblem& problem) {
    TurningPoints result;
    const double axis=pairBohrRadius(activePair);
    // h -> -infinity as r -> 0 (centrifugal) and -> eps < 0 as r -> infinity,
    // so a bound orbit has a single maximum with two roots around it.
    // h need not be unimodal once the dipole terms are in (u_dd goes as
    // r^-3), so the maximum is found by SCANNING a logarithmic grid rather
    // than by a ternary search, which would silently pick the wrong hump and
    // report a state as inaccessible.  The grid is refined around its own
    // best point before the sign test.
    const int gridPoints=20001;
    const double lowest=1.0e-9*axis,highest=1.0e4*axis;
    double peak=lowest,best=-std::numeric_limits<double>::infinity();
    for(int index=0;index<gridPoints;++index) {
        const double radius=lowest*std::pow(highest/lowest,
            static_cast<double>(index)/(gridPoints-1.0));
        const double value=radialFunction(problem,radius);
        if(value>best) { best=value; peak=radius; }
    }
    {
        const double span=std::pow(highest/lowest,2.0/(gridPoints-1.0));
        double low=peak/span,high=peak*span;
        for(int iteration=0;iteration<300;++iteration) {
            const double a=low+(high-low)/3.0,b=high-(high-low)/3.0;
            if(radialFunction(problem,a)<radialFunction(problem,b)) low=a;
            else high=b;
        }
        const double refined=0.5*(low+high);
        if(radialFunction(problem,refined)>best) {
            best=radialFunction(problem,refined); peak=refined;
        }
    }
    if(!(best>0.0)) return result;
    const auto bisect=[&](double inside,double outside) {
        for(int iteration=0;iteration<200;++iteration) {
            const double middle=0.5*(inside+outside);
            if(radialFunction(problem,middle)>0.0) inside=middle;
            else outside=middle;
        }
        return 0.5*(inside+outside);
    };
    result.inner=bisect(peak,1.0e-9*axis);
    result.outer=bisect(peak,1.0e4*axis);
    result.accessible=true;
    return result;
}

// Gauss-Legendre nodes on [-1,1] by Newton on the Legendre polynomial.
void gaussLegendre(int count,std::vector<double>& nodes,
                   std::vector<double>& weights) {
    nodes.assign(count,0.0); weights.assign(count,0.0);
    for(int index=0;index<count;++index) {
        double x=std::cos(pi*(index+0.75)/(count+0.5));
        for(int iteration=0;iteration<100;++iteration) {
            double previous=1.0,current=x,derivative=0.0;
            for(int degree=2;degree<=count;++degree) {
                const double next=((2.0*degree-1.0)*x*current
                                   -(degree-1.0)*previous)/degree;
                previous=current; current=next;
            }
            derivative=count*(x*current-previous)/(x*x-1.0);
            const double step=current/derivative;
            x-=step;
            if(std::abs(step)<1.0e-15) break;
        }
        double previous=1.0,current=x;
        for(int degree=2;degree<=count;++degree) {
            const double next=((2.0*degree-1.0)*x*current
                               -(degree-1.0)*previous)/degree;
            previous=current; current=next;
        }
        const double derivative=count*(x*current-previous)/(x*x-1.0);
        nodes[index]=x;
        weights[index]=2.0/((1.0-x*x)*derivative*derivative);
    }
}

double radialAction(const RadialProblem& problem,const TurningPoints& points,
                    int count) {
    if(!points.accessible) return 0.0;
    std::vector<double> nodes,weights;
    gaussLegendre(count,nodes,weights);
    const double centre=0.5*(points.inner+points.outer);
    const double half=0.5*(points.outer-points.inner);
    const double reducedMass=reducedMassOf(activePair);
    double sum=0.0;
    for(int index=0;index<count;++index) {
        const double angle=0.5*pi*(nodes[index]+1.0);        // [0, pi]
        const double radius=centre-half*std::cos(angle);
        const double value=radialFunction(problem,radius);
        if(value<=0.0) continue;
        sum+=weights[index]*0.5*pi*std::sqrt(2.0*value)*half*std::sin(angle);
    }
    return 2.0*reducedMass*sum;
}

double keplerRadialAction(const RadialProblem& problem) {
    const double reducedMass=reducedMassOf(activePair);
    const double binding=-problem.specificEnergy*reducedMass;
    if(!(binding>0.0)) return 0.0;
    return 2.0*pi*(pairCoulombStrength*std::sqrt(reducedMass/(2.0*binding))
                   -problem.specificAngularMomentum*reducedMass);
}

}  // namespace

int main(int argc,char** argv) {
    const double planck=2.0*pi*hbar;
    const double reducedMass=reducedMassOf(activePair);
    const double attraction=pairCoulombStrength/reducedMass;
    if(argc>1&&std::string(argv[1])=="--selftest") {
        // R133a: with the dipoles zeroed the integral must reproduce the
        // Kepler closed form.
        std::printf("R133a self-test, dipole terms zeroed:\n");
        std::printf("%12s %18s %18s %14s\n","e","J_r/h numerical",
                    "J_r/h closed form","relative");
        for(double eccentricity:{0.1,0.5,0.9}) {
            const double axis=pairBohrRadius(activePair);
            RadialProblem problem;
            problem.specificEnergy=-attraction/(2.0*axis);
            problem.specificAngularMomentum=std::sqrt(
                attraction*axis*(1.0-eccentricity*eccentricity));
            const TurningPoints points=turningPoints(problem);
            const double numerical=radialAction(problem,points,400)/planck;
            const double closed=keplerRadialAction(problem)/planck;
            std::printf("%12.3f %18.12f %18.12f %14.2e\n",eccentricity,
                        numerical,closed,
                        std::abs(numerical-closed)
                            /std::max(std::abs(closed),1.0e-300));
        }
        // And the turning points themselves.
        std::printf("turning points at e = 0.5: inner/outer %.12f %.12f "
                    "a_pair against 0.5 and 1.5\n",
                    [&]{ RadialProblem p; p.specificEnergy=
                         -attraction/(2.0*pairBohrRadius(activePair));
                         p.specificAngularMomentum=std::sqrt(
                             attraction*pairBohrRadius(activePair)*0.75);
                         return turningPoints(p).inner
                             /pairBohrRadius(activePair); }(),
                    [&]{ RadialProblem p; p.specificEnergy=
                         -attraction/(2.0*pairBohrRadius(activePair));
                         p.specificAngularMomentum=std::sqrt(
                             attraction*pairBohrRadius(activePair)*0.75);
                         return turningPoints(p).outer
                             /pairBohrRadius(activePair); }());
        return 0;
    }
    const int nodes=argc>1?atoi(argv[1]):400;
    std::printf("%10s %12s %14s %14s %14s %14s %12s %10s\n","checkpoint",
                "a/a_pair","J_r/h model","J_r/h Kepler","J_phi/h",
                "total model","total Kepler","r_out/r_in");
    std::string line;
    int inaccessible=0,rows=0;
    double firstTotal=0.0,lastTotal=0.0,maximumRadial=0.0;
    while(std::getline(std::cin,line)) {
        if(line.rfind("ACTION",0)!=0) continue;
        const auto number=[&](const std::string& key)->double {
            const std::size_t at=line.find(key);
            if(at==std::string::npos) return 0.0;
            return std::atof(line.c_str()+at+key.size());
        };
        const auto vector=[&](const std::string& key)->Vec3 {
            const std::size_t at=line.find(key);
            if(at==std::string::npos) return Vec3{};
            const char* text=line.c_str()+at+key.size();
            Vec3 value{};
            std::sscanf(text,"%lf,%lf,%lf",&value.x,&value.y,&value.z);
            return value;
        };
        RadialProblem problem;
        problem.specificEnergy=number("eps=");
        problem.specificAngularMomentum=number("L=");
        problem.firstDipole=vector("mu1=");
        problem.secondDipole=vector("mu2=");
        problem.normal=vector("nhat=");
        if(!(problem.specificEnergy<0.0)) continue;
        const int checkpoint=static_cast<int>(number("checkpoint="));
        const double axisOverPair=number("a_over_apair=");
        const TurningPoints points=turningPoints(problem);
        const double azimuthal=
            problem.specificAngularMomentum*reducedMass/hbar;
        const double keplerRadial=keplerRadialAction(problem)/planck;
        if(!points.accessible) {
            ++inaccessible;
            std::printf("%10d %12.6f %14s %14.6f %14.6f %14s %12.6f %10s\n",
                        checkpoint,axisOverPair,"NO REGION",keplerRadial,
                        azimuthal,"-",keplerRadial+azimuthal,"-");
            continue;
        }
        const double coarse=radialAction(problem,points,nodes/2)/planck;
        const double radial=radialAction(problem,points,nodes)/planck;
        (void)coarse;
        ++rows;
        if(rows==1) firstTotal=radial+azimuthal;
        lastTotal=radial+azimuthal;
        maximumRadial=std::max(maximumRadial,radial);
        std::printf("%10d %12.6f %14.6f %14.6f %14.6f %14.6f %12.6f %10.4f\n",
                    checkpoint,axisOverPair,radial,keplerRadial,azimuthal,
                    radial+azimuthal,keplerRadial+azimuthal,
                    points.outer/std::max(points.inner,1.0e-300));
    }
    std::printf("\nrows %d, inaccessible %d, max J_r/h %.6f\n",
                rows,inaccessible,maximumRadial);
    std::printf("total action (J_r+J_phi)/h: first %.9f, last %.9f\n",
                firstTotal,lastTotal);
}
