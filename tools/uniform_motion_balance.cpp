// Does the retarded dipole sector conserve momentum in UNIFORM MOTION
// (audit section 122)?
//
// A pair drifting rigidly -- same constant velocity, fixed separation, no
// acceleration -- carries a field pattern that simply translates, so the
// field momentum is constant and Newton's third law must hold exactly:
// F1 + F2 = 0.  Radiation, retardation transients and hidden momentum are
// all absent or constant, so any net force here is a defect of the force
// construction itself rather than of the bookkeeping (audits 119-121 having
// measured the three candidate carriers and found them 1-2 orders too small).
//
// The history is synthesized as uniform motion, which is what the retarded
// fields need to see; the accelerations are exactly zero.
//
// With "circular" as the sixth argument the history is exact CIRCULAR motion
// instead: the same separation, a chosen speed, and therefore a chosen
// acceleration v^2/r.  Uniform motion must give zero net force; circular
// motion is where audits 104-106 see the imbalance, and sweeping the speed
// gives its power in beta, which says which term of the retarded field
// carries it.
//
// Usage: uniform_motion_balance [beta] [separation in r*] [tilt deg]
//        [history nodes] [history span in r/c] [circular]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/uniform_motion_balance.cpp -o /tmp/uniform $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc,char** argv) {
    const double beta=argc>1?atof(argv[1]):0.05;
    const double separation=(argc>2?atof(argv[2]):1.0)*comptonBarrierRadius;
    const double tilt=(argc>3?atof(argv[3]):0.0)*pi/180;
    const int nodes=argc>4?atoi(argv[4]):512;
    const double span=argc>5?atof(argv[5]):40.0;
    const bool circular=argc>6&&std::string(argv[6])=="circular";
    const Vec3 drift{0.0,beta*c,0.0};
    // CREM_BALANCE_NO_MOMENTS zeroes both magnetic moments, leaving the two
    // charges and their retarded fields alone.  It answers "is the residual
    // in the dipole sector at all", which no ablation of a dipole TERM can
    // answer, since those leave the moment's own field in place (audit 142).
    const double momentScale=[] {
        if(std::getenv("CREM_BALANCE_NO_MOMENTS")) return 0.0;
        // CREM_BALANCE_MOMENT_SCALE multiplies both moments, so the residual's
        // POWER in mu separates a charge-moment interaction (linear) from a
        // moment-moment one (quadratic) without ablating any term.
        const char* text=std::getenv("CREM_BALANCE_MOMENT_SCALE");
        const double value=text?std::atof(text):1.0;
        return std::isfinite(value)?value:1.0;
    }();
    const Vec3 direction{std::sin(tilt),0.0,std::cos(tilt)};
    const double angular=circular?beta*c/(0.5*separation):0.0;
    const auto build=[&](double time) {
        State s{};
        s.time=time;
        if(circular) {
            // Both particles on a circle of radius separation/2 about the
            // centre, opposite phases: exact circular kinematics, so the
            // history carries a genuine acceleration v^2/(separation/2).
            const double phase=angular*time;
            const double armRadius=0.5*separation;
            const Vec3 offset{armRadius*std::cos(phase),
                              armRadius*std::sin(phase),0.0};
            const Vec3 velocity{-armRadius*angular*std::sin(phase),
                                 armRadius*angular*std::cos(phase),0.0};
            const Vec3 acceleration{-armRadius*angular*angular*std::cos(phase),
                                    -armRadius*angular*angular*std::sin(phase),
                                    0.0};
            s.firstPosition=offset; s.secondPosition=offset*-1.0;
            s.firstVelocity=velocity; s.secondVelocity=velocity*-1.0;
            s.firstAcceleration=acceleration;
            s.secondAcceleration=acceleration*-1.0;
            s.firstProperDipole=direction*(firstMagneticMoment*momentScale);
            s.secondProperDipole=direction*(secondMagneticMoment*momentScale);
            synchronizeCovariantDipoles(s);
            return s;
        }
        s.firstPosition=Vec3{0.5*separation,0.0,0.0}+drift*time;
        s.secondPosition=Vec3{-0.5*separation,0.0,0.0}+drift*time;
        s.firstVelocity=drift;
        s.secondVelocity=drift;
        s.firstProperDipole=direction*(firstMagneticMoment*momentScale);
        s.secondProperDipole=direction*(secondMagneticMoment*momentScale);
        synchronizeCovariantDipoles(s);
        return s;
    };
    StateHistory history;
    const double lightTime=separation/c;
    for(int index=0;index<nodes;++index) {
        const double time=-span*lightTime
            *(1.0-static_cast<double>(index)/nodes);
        history.push_back(build(time));
    }
    const State present=build(0.0);
    const MutualForces forces=retardedExternalForces(present,history);
    const Vec3 net=forces.first+forces.second;
    const double reference=std::max(forces.first.norm(),
                                    forces.second.norm());
    std::printf("beta %.4f  separation %.2f r*  tilt %.0f deg  nodes %d  "
                "span %.0f r/c\n",beta,separation/comptonBarrierRadius,
                tilt*180/pi,nodes,span);
    std::printf("  |F1| %.6e  |F2| %.6e  |F1+F2| %.6e  ratio %.3e\n",
                forces.first.norm(),forces.second.norm(),net.norm(),
                net.norm()/std::max(reference,1.0e-300));
    std::printf("  F1 %+.6e %+.6e %+.6e\n  F2 %+.6e %+.6e %+.6e\n",
                forces.first.x,forces.first.y,forces.first.z,
                forces.second.x,forces.second.y,forces.second.z);
    // Is the net force the rate of change of the pair's hidden momentum?
    // A moment in an electric field carries (mu x E)/c^2, and in circular
    // motion that vector rotates with the field, so its rate is first order
    // in beta -- the same order the net force turns out to have.  Both sides
    // are computed from the synthetic kinematics, the rate by a centred
    // difference over the same history (audit 122).
    {
        const auto hiddenAt=[&](double time) {
            const State sample=build(time);
            StateHistory shifted;
            for(int index=0;index<nodes;++index) {
                const double node=time-span*lightTime
                    *(1.0-static_cast<double>(index)/nodes);
                shifted.push_back(build(node));
            }
            const auto electric=[&](bool targetIsFirst) {
                const Vec3 position=targetIsFirst?sample.firstPosition
                                                 :sample.secondPosition;
                return fieldFromOtherParticleAt(position,sample.time,sample,
                    shifted,targetIsFirst).electric;
            };
            return (cross(sample.firstDipole,electric(true))
                   +cross(sample.secondDipole,electric(false)))*(1.0/(c*c));
        };
        const double step=circular&&angular!=0.0
            ?1.0e-3*2.0*pi/std::abs(angular):1.0e-3*lightTime;
        const Vec3 rate=(hiddenAt(step)-hiddenAt(-step))*(1.0/(2.0*step));
        std::printf("  hidden momentum rate d/dt (mu x E)/c^2: %.6e N "
                    "(net force %.6e N, ratio %.4f)\n",rate.norm(),
                    net.norm(),rate.norm()/std::max(net.norm(),1.0e-300));
        std::printf("  net + rate %.6e N; net - rate %.6e N\n",
                    (net+rate).norm(),(net-rate).norm());
    }
    // And the rate of the interaction FIELD momentum, the other term a
    // consistent balance needs: eps0 integral (E1 x B2 + E2 x B1) dV, on the
    // same synthetic history, differenced in time the same way (audit 122).
    if(std::getenv("UNIFORM_FIELD_RATE")) {
        const int radialNodes=48,polarNodes=12,azimuthalNodes=24;
        const double outer=12.0*comptonBarrierRadius;
        const double inner=0.05*comptonBarrierRadius;
        // Gauss-Legendre in cos(theta).
        std::vector<double> nodesCos(polarNodes),weights(polarNodes);
        for(int index=0;index<polarNodes;++index) {
            double x=std::cos(pi*(index+0.75)/(polarNodes+0.5));
            for(int iteration=0;iteration<100;++iteration) {
                double p0=1.0,p1=0.0;
                for(int degree=0;degree<polarNodes;++degree) {
                    const double p2=p1; p1=p0;
                    p0=((2.0*degree+1.0)*x*p1-degree*p2)/(degree+1.0);
                }
                const double derivative=polarNodes*(x*p0-p1)/(x*x-1.0);
                const double stepSize=p0/derivative; x-=stepSize;
                if(std::abs(stepSize)<1.0e-15) break;
            }
            double p0=1.0,p1=0.0;
            for(int degree=0;degree<polarNodes;++degree) {
                const double p2=p1; p1=p0;
                p0=((2.0*degree+1.0)*x*p1-degree*p2)/(degree+1.0);
            }
            const double derivative=polarNodes*(x*p0-p1)/(x*x-1.0);
            nodesCos[static_cast<size_t>(index)]=x;
            weights[static_cast<size_t>(index)]=
                2.0/((1.0-x*x)*derivative*derivative);
        }
        const auto fieldAt=[&](double time) {
            const State sample=build(time);
            StateHistory shifted;
            for(int index=0;index<nodes;++index)
                shifted.push_back(build(time-span*lightTime
                    *(1.0-static_cast<double>(index)/nodes)));
            const Vec3 centre=(sample.firstPosition*firstMass
                +sample.secondPosition*secondMass)/(firstMass+secondMass);
            Vec3 running;
            for(int shell=0;shell<radialNodes;++shell) {
                const double lowEdge=inner*std::pow(outer/inner,
                    static_cast<double>(shell)/radialNodes);
                const double highEdge=inner*std::pow(outer/inner,
                    static_cast<double>(shell+1)/radialNodes);
                const double centreRadius=std::sqrt(lowEdge*highEdge);
                const double shellVolume=(highEdge*highEdge*highEdge
                    -lowEdge*lowEdge*lowEdge)/3.0;
                Vec3 angularSum;
                for(int polarIndex=0;polarIndex<polarNodes;++polarIndex) {
                    const double cosTheta=nodesCos[
                        static_cast<size_t>(polarIndex)];
                    const double weight=weights[
                        static_cast<size_t>(polarIndex)];
                    const double sinTheta=std::sqrt(std::max(0.0,
                        1.0-cosTheta*cosTheta));
                    for(int slice=0;slice<azimuthalNodes;++slice) {
                        const double phi=2.0*pi*(slice+0.5)/azimuthalNodes;
                        const Vec3 point=centre+Vec3{
                            centreRadius*sinTheta*std::cos(phi),
                            centreRadius*sinTheta*std::sin(phi),
                            centreRadius*cosTheta};
                        ElectromagneticField fromFirst=lienardWiechertField(
                            point,sample.time,shifted,sample,true,firstCharge);
                        const ElectromagneticField firstDipole=
                            retardedMagneticDipoleField(point,sample.time,
                                shifted,sample,true);
                        fromFirst.electric+=firstDipole.electric;
                        fromFirst.magnetic+=firstDipole.magnetic;
                        ElectromagneticField fromSecond=lienardWiechertField(
                            point,sample.time,shifted,sample,false,
                            secondCharge);
                        const ElectromagneticField secondDipole=
                            retardedMagneticDipoleField(point,sample.time,
                                shifted,sample,false);
                        fromSecond.electric+=secondDipole.electric;
                        fromSecond.magnetic+=secondDipole.magnetic;
                        const Vec3 density=
                            (cross(fromFirst.electric,fromSecond.magnetic)
                            +cross(fromSecond.electric,fromFirst.magnetic))
                            *epsilon0;
                        if(isFinite(density)) angularSum+=density*weight;
                    }
                }
                running+=angularSum*(shellVolume*2.0*pi
                    /static_cast<double>(azimuthalNodes));
            }
            return running;
        };
        const double stepSize=circular&&angular!=0.0
            ?1.0e-3*2.0*pi/std::abs(angular):1.0e-3*lightTime;
        const Vec3 fieldRate=
            (fieldAt(stepSize)-fieldAt(-stepSize))*(1.0/(2.0*stepSize));
        std::printf("  field momentum rate: %.6e N (net force %.6e N, "
                    "ratio %.4f)\n",fieldRate.norm(),net.norm(),
                    fieldRate.norm()/std::max(net.norm(),1.0e-300));
    }
    // Is the remainder the difference between the GRADIENT form the model uses
    // and the CONVECTIVE form a dipole force has (audit 124)?
    //   grad(mu.B) - (mu.grad)B = mu x (curl B),
    //   grad(p.E)  - (p.grad)E  = p  x (curl E),
    // and for the partner's field curl E = -dB/dt, curl B = (1/c^2) dE/dt, so
    // both are first order in the frequency -- the order the remainder has.
    {
        const double stencil=1.0e-4*separation;
        const auto curls=[&](bool targetIsFirst,Vec3& curlE,Vec3& curlB) {
            const Vec3 position=targetIsFirst?present.firstPosition
                                             :present.secondPosition;
            Vec3 electricPlus[3],electricMinus[3];
            Vec3 magneticPlus[3],magneticMinus[3];
            for(int axis=0;axis<3;++axis) {
                Vec3 offset;
                (axis==0?offset.x:axis==1?offset.y:offset.z)=stencil;
                const ElectromagneticField plus=fieldFromOtherParticleAt(
                    position+offset,present.time,present,history,targetIsFirst);
                const ElectromagneticField minus=fieldFromOtherParticleAt(
                    position-offset,present.time,present,history,targetIsFirst);
                electricPlus[axis]=plus.electric;
                electricMinus[axis]=minus.electric;
                magneticPlus[axis]=plus.magnetic;
                magneticMinus[axis]=minus.magnetic;
            }
            const auto component=[&](const Vec3* plus,const Vec3* minus) {
                const double inverse=1.0/(2.0*stencil);
                return Vec3{
                    (plus[1].z-minus[1].z)*inverse-(plus[2].y-minus[2].y)*inverse,
                    (plus[2].x-minus[2].x)*inverse-(plus[0].z-minus[0].z)*inverse,
                    (plus[0].y-minus[0].y)*inverse-(plus[1].x-minus[1].x)*inverse};
            };
            curlE=component(electricPlus,electricMinus);
            curlB=component(magneticPlus,magneticMinus);
        };
        Vec3 curlEFirst,curlBFirst,curlESecond,curlBSecond;
        curls(true,curlEFirst,curlBFirst);
        curls(false,curlESecond,curlBSecond);
        const Vec3 gradientExcess=
            cross(present.firstDipole,curlBFirst)
           +cross(present.firstElectricDipole,curlEFirst)
           +cross(present.secondDipole,curlBSecond)
           +cross(present.secondElectricDipole,curlESecond);
        const double projection=net.norm()>0.0
            ?dot(gradientExcess,net)/net.norm():0.0;
        std::printf("  gradient-minus-convective excess: |S| %.6e N, "
                    "projection on the net force %.6e N (%.4f of it), "
                    "|net - S| %.6e N (%.4f)\n",
                    gradientExcess.norm(),projection,
                    projection/std::max(net.norm(),1.0e-300),
                    (net-gradientExcess).norm(),
                    (net-gradientExcess).norm()/std::max(net.norm(),1.0e-300));
    }
    // The same net force with single parts of the dipole sector removed, when
    // the binary carries the probe switches of audit 105.
    const char* switches[]={"CREM_PROBE_NO_POLES","CREM_PROBE_NO_MAGNETIZATION",
        "CREM_PROBE_NO_LW_IN_MOMENT","CREM_PROBE_NO_MATERIAL",
        "CREM_PROBE_NO_DIPOLE_ON_CHARGE","CREM_PROBE_NO_MOMENT_FORCE"};
    for(const char* name:switches) {
        setenv(name,"1",1);
        const MutualForces ablated=retardedExternalForces(present,history);
        unsetenv(name);
        const Vec3 ablatedNet=ablated.first+ablated.second;
        std::printf("  without %-28s |F1+F2| %.6e  ratio %.3e\n",name,
                    ablatedNet.norm(),ablatedNet.norm()
                        /std::max(reference,1.0e-300));
    }
}
