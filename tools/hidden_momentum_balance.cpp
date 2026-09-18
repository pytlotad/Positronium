// Does hidden momentum close the pair's momentum balance (audit section 119)?
//
// A magnetic moment in an electric field carries hidden momentum
// (mu x E)/c^2, and section 105 measured that 88% of the pair's unbalanced
// force is the moment sitting in the partner's retarded charge field -- the
// same configuration.  This probe follows the trajectory of sections 104-106
// (para, 1 r*, tilt 0, retarded forces, one orbit) and accumulates, at every
// step, the mechanical momentum and the hidden momentum, for the partner's
// FULL field and for its Lienard-Wiechert charge part alone.  What is tested
// is whether P_mech + P_hidden is conserved where P_mech alone is not.
//
// Usage: hidden_momentum_balance [r/r*, default 1] [channel 0 para|1 ortho]
//        [tilt deg] [orbits] [steps per orbit] [depth]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/hidden_momentum_balance.cpp -o /tmp/hidden $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const double radius=(argc>1?atof(argv[1]):1.0)*comptonBarrierRadius;
    const int channel=argc>2?atoi(argv[2]):0;
    const double tilt=(argc>3?atof(argv[3]):0.0)*pi/180;
    const int orbits=argc>4?atoi(argv[4]):1;
    const int perOrbit=argc>5?atoi(argv[5]):256;
    const int depth=argc>6?atoi(argv[6]):20;
    const double reduced=pairReducedMass,total=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(reduced*radius));
    const double period=2.0*pi*radius/speed;
    const Vec3 direction{std::sin(tilt),0.0,std::cos(tilt)};
    State state{};
    state.firstPosition={radius*secondMass/total,0.0,0.0};
    state.secondPosition={-radius*firstMass/total,0.0,0.0};
    state.firstVelocity={0.0,speed*secondMass/total,0.0};
    state.secondVelocity={0.0,-speed*firstMass/total,0.0};
    state.firstProperDipole=direction*firstMagneticMoment;
    state.secondProperDipole=direction
        *(channel?-secondMagneticMoment:secondMagneticMoment);
    synchronizeCovariantDipoles(state);
    ClassicalTrajectoryEngine::Accuracy accuracy;
    accuracy.relativeTolerance=1.0e-8;
    accuracy.maximumDepth=depth;
    accuracy.reactionModel=ChargeRadiationReactionModel::disabled;
    accuracy.computeOutwardFlux=false;
    accuracy.useRetardedExternalForces=true;
    ClassicalTrajectoryEngine engine(state,accuracy);
    const double scale=firstMass*c;
    const auto mechanical=[&](const State& s) {
        return momentum(s.firstVelocity,firstMass)
            +momentum(s.secondVelocity,secondMass);
    };
    // (mu x E)/c^2 for both particles, with E taken from the partner.
    const auto hidden=[&](const State& s,bool chargeOnly) {
        const auto electric=[&](bool targetIsFirst) {
            const Vec3 position=targetIsFirst?s.firstPosition:s.secondPosition;
            if(chargeOnly)
                return lienardWiechertField(position,s.time,engine.history(),s,
                    !targetIsFirst,targetIsFirst?secondCharge:firstCharge)
                    .electric;
            return fieldFromOtherParticleAt(position,s.time,s,engine.history(),
                targetIsFirst).electric;
        };
        return (cross(s.firstDipole,electric(true))
               +cross(s.secondDipole,electric(false)))*(1.0/(c*c));
    };
    const Vec3 mechanicalStart=mechanical(state);
    const Vec3 hiddenStart=hidden(state,false);
    const Vec3 hiddenChargeStart=hidden(state,true);
    std::printf("start |P_mech| %.4e mc  |P_hidden(full)| %.4e mc  "
                "|P_hidden(charge)| %.4e mc\n",
                mechanicalStart.norm()/scale,hiddenStart.norm()/scale,
                hiddenChargeStart.norm()/scale);
    const int steps=orbits*perOrbit;
    const double step=period/perOrbit;
    for(int index=0;index<steps;++index) {
        if(!engine.advance(state,step)) {
            std::printf("FAILED after %d steps at r %.4f r*\n",index,
                        separation(state)/comptonBarrierRadius);
            break;
        }
        if((index+1)%(perOrbit/4)==0) {
            const Vec3 mech=mechanical(state);
            const Vec3 hid=hidden(state,false);
            const Vec3 hidCharge=hidden(state,true);
            std::printf("t/P %.2f r %.4f r*  |dP_mech| %.6f  "
                        "|d(P_mech+P_hid_full)| %.6f  "
                        "|d(P_mech+P_hid_charge)| %.6f  |P_hid_full| %.6f\n",
                        (index+1.0)/perOrbit,
                        separation(state)/comptonBarrierRadius,
                        (mech-mechanicalStart).norm()/scale,
                        ((mech+hid)-(mechanicalStart+hiddenStart)).norm()/scale,
                        ((mech+hidCharge)
                            -(mechanicalStart+hiddenChargeStart)).norm()/scale,
                        hid.norm()/scale);
        }
    }
    const Vec3 mech=mechanical(state);
    const Vec3 hid=hidden(state,false);
    const Vec3 hidCharge=hidden(state,true);
    const double mechanicalChange=(mech-mechanicalStart).norm()/scale;
    const double fullChange=
        ((mech+hid)-(mechanicalStart+hiddenStart)).norm()/scale;
    const double chargeChange=
        ((mech+hidCharge)-(mechanicalStart+hiddenChargeStart)).norm()/scale;
    std::printf("\nover the run: |dP_mech| %.6f mc; with hidden momentum "
                "(full field) %.6f, ratio %.3f; (charge field only) %.6f, "
                "ratio %.3f\n",mechanicalChange,fullChange,
                mechanicalChange/std::max(fullChange,1.0e-300),chargeChange,
                mechanicalChange/std::max(chargeChange,1.0e-300));
}
