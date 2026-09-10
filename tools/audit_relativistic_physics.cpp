// Reproduction probes for audits/2026-09-09-relativistic-physics.md.
// These call production functions on prescribed source histories. They print
// measurements, not the pass/fail verdict of the regular validation suite.
// Build from the repository root:
// g++ -std=c++20 -O2 -I . tools/audit_relativistic_physics.cpp -o /tmp/positronium-physics-audit
#include "modules/electrodynamics.hpp"
#include <iostream>
#include <iomanip>
State sample(double time, double boost, double speed) {
    const double R=2*bohrRadius;
    const double g=1/std::sqrt(1-boost*boost);
    State s;
    s.time=time;
    s.firstPosition={(R/g+(speed-boost)*c*time)/(1-boost*speed),0,0};
    s.secondPosition={-boost*c*time,0,0};
    s.firstVelocity={(speed-boost)*c/(1-boost*speed),0,0};
    s.secondVelocity={-boost*c,0,0};
    s.firstProperDipole={0,0,firstMagneticMoment};
    s.secondProperDipole={0,0,secondMagneticMoment};
    synchronizeCovariantDipoles(s);
    return s;
}
int main() {
    // Isolate the dipole force; prescribed inertial source histories do not
    // have to solve the interacting two-body initial-value problem.
    secondCharge=0;
    std::cout<<std::setprecision(12);
    for(double speed:{0.0,0.2,0.6}) {
        double ref=0;
        for(double boost:{0.0,0.35,0.6}) {
            const double time=-boost*(2*bohrRadius)/c/std::sqrt(1-boost*boost);
            State s=sample(time,boost,speed);
            StateHistory h;
            for(int k=128;k>=0;--k) h.push_back(sample(time-k*1e-19,boost,speed));
            Vec3 f=covariantDipoleGradientForce(s,h,true);
            if(boost==0) ref=f.x;
            std::cout<<"force beta="<<speed<<" boost="<<boost<<" Fx="<<f.x<<" Fprime/F="<<f.x/ref<<" expected=1\n";
        }
    }
    for(double beta:{0.0,0.3,0.8,0.999,0.9995,0.9999}) {
        State s;
        s.firstPosition={2*bohrRadius,0,0};
        s.secondVelocity={0,beta*c,0};
        s.secondProperDipole={0,0,secondMagneticMoment};
        synchronizeCovariantDipoles(s);
        StateHistory h;
        for(int k=128;k>=0;--k) {
            State p=s;
            p.time=-k*1e-17;
            p.secondPosition=s.secondVelocity*p.time;
            h.push_back(p);
        }
        const auto f=retardedMagneticDipoleField(s.firstPosition,0,h,s,false);
        const double R=2*bohrRadius;
        const double expectedB=-(mu0/(4*pi))*secondMagneticMoment
            /(R*R*R)/std::sqrt(1-beta*beta);
        std::cout<<"field beta="<<beta<<" Bz="<<f.magnetic.z
                 <<" Bz/boosted_static="<<f.magnetic.z/expectedB<<" expected=1\n";
    }

    ZeroPointField zpf;
    zpf.direction={{0,0,1}};
    zpf.polarization={{1,0,0}};
    zpf.magneticDirection={{0,1,0}};
    zpf.phase={0};
    zpf.frequencyFactor={1};
    const double omega0=1e16, chirp=1e14;
    zpf.amplitudeCoefficient=1/(omega0*omega0);
    for(double h:{1e-18,1e-19,1e-20}) {
        Vec3 ep,bp,em,bm;
        zpf.sample({},omega0*(1+chirp*h),omega0*(h+chirp*h*h/2),ep,bp);
        zpf.sample({},omega0*(1-chirp*h),omega0*(-h+chirp*h*h/2),em,bm);
        Vec3 rightE,rightB,leftE,leftB;
        const double dx=c*h;
        zpf.sample({0,0,dx},omega0,0,rightE,rightB);
        zpf.sample({0,0,-dx},omega0,0,leftE,leftB);
        const double dBydt=(bp.y-bm.y)/(2*h);
        const double curlEy=(rightE.x-leftE.x)/(2*dx);
        std::cout<<"zpf h="<<h<<" Faraday_y="<<dBydt+curlEy
                 <<" expected=0 normalized_to_2chirp_over_c="
                 <<(dBydt+curlEy)/(2*chirp/c)<<"\n";
    }

}
