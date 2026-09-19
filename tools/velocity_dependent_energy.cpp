// What does the dynamics actually conserve (audit section 145)?
//
// The model's ledger is T + U_Coulomb + U_dd + U_chargeDipole + U_Darwin +
// U_constraint, and three of those depend on VELOCITY: the charge-dipole term
// is -p.E with p = (v x mu)/c^2, linear in v; the Darwin term goes as v1.v2,
// quadratic; and U_dd is evaluated on the LAB moments, which
// synchronizeCovariantDipoles boosts.  For a Lagrangian L = T - U(r,v) the
// conserved quantity is H = sum v.(dL/dv) - L, not T + U: a velocity-LINEAR
// potential drops out of H entirely and a velocity-quadratic one enters with
// the opposite sign.  So the ledger may be summing quantities that do not
// belong to the conserved energy at all.
//
// Four ledgers are read off the SAME states of one engine orbit, so no
// difference can come from the trajectories diverging:
//   E0  the model's ledger as it stands
//   E1  E0 without the charge-dipole term
//   E2  E1 with the Darwin term's sign flipped (2 U_Darwin subtracted)
//   E3  E2 with U_dd on the PROPER moments instead of the lab ones
//
// Usage: velocity_dependent_energy <r/r*> <channel 0 para|1 ortho|2 none>
//        <tilt deg> <orbits> <tolerance> [steps per orbit, default 256]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/velocity_dependent_energy.cpp -o /tmp/vde $(root-config --libs)
#include "modules/crem_engine.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc,char** argv) {
    const double r=(argc>1?atof(argv[1]):1.0)*comptonBarrierRadius;
    const int channel=argc>2?atoi(argv[2]):0;
    const double tilt=(argc>3?atof(argv[3]):0.0)*pi/180.0;
    const int orbits=argc>4?atoi(argv[4]):1;
    const double tolerance=argc>5?atof(argv[5]):1.0e-8;
    const int perOrbit=argc>6?atoi(argv[6]):256;
    // Rule R145b's fallback candidate, audit 125e: the ledger is an
    // INSTANTANEOUS functional of the state while the forces are RETARDED.
    // Passing 0 here integrates with instantaneous forces instead, which is
    // the configuration the ledger would be exact for.
    const bool retarded=!(argc>7&&atoi(argv[7])==0);
    const double totalMass=firstMass+secondMass;
    const double speed=std::sqrt(pairCoulombStrength/(pairReducedMass*r));
    const double period=2.0*pi*r/speed;
    const Vec3 direction{std::sin(tilt),0.0,std::cos(tilt)};
    State state{};
    state.firstPosition={r*secondMass/totalMass,0.0,0.0};
    state.secondPosition={-r*firstMass/totalMass,0.0,0.0};
    state.firstVelocity={0.0,speed*secondMass/totalMass,0.0};
    state.secondVelocity={0.0,-speed*firstMass/totalMass,0.0};
    if(channel!=2) {
        state.firstProperDipole=direction*firstMagneticMoment;
        state.secondProperDipole=direction
            *(channel?-secondMagneticMoment:secondMagneticMoment);
    }
    synchronizeCovariantDipoles(state);
    const double dipoleScale=std::abs(pairDipoleInteractionEnergy(
        state.firstPosition-state.secondPosition,
        direction*firstMagneticMoment,direction*secondMagneticMoment));
    // The four ledgers, all read off one state.
    struct Ledgers { double e0,e1,e2,e3; };
    const auto ledgers=[&](const State& s) {
        const PairGeometry geometry=clampedPairGeometry(s);
        const double kinetic=kineticEnergy(s.firstVelocity,firstMass)
            +kineticEnergy(s.secondVelocity,secondMass);
        const double coulombPart=-pairCoulombStrength*geometry.inverseDistance;
        const Vec3 separation=s.firstPosition-s.secondPosition;
        const double dipoleLab=pairDipoleInteractionEnergy(
            separation,s.firstDipole,s.secondDipole);
        const double dipoleProper=pairDipoleInteractionEnergy(
            separation,s.firstProperDipole,s.secondProperDipole);
        const double chargeDipole=chargeDipoleInteractionEnergy(s);
        const double darwin=darwinInteractionEnergy(s);
        const double base=kinetic+coulombPart+s.dipoleConstraintEnergy;
        Ledgers out;
        out.e0=base+dipoleLab+chargeDipole+darwin;
        out.e1=base+dipoleLab+darwin;
        out.e2=base+dipoleLab-darwin;
        out.e3=base+dipoleProper-darwin;
        return out;
    };
    ClassicalTrajectoryEngine::Accuracy accuracy;
    accuracy.relativeTolerance=tolerance;
    accuracy.maximumDepth=20;
    // With the reaction disabled the retarded fields still radiate and nothing
    // is debited for it, so part of the ledger's drift is simply unpaid
    // outgoing flux.  The engine accumulates it when asked, and audit 156
    // compares it with what is left unexplained.
    accuracy.reactionModel=ChargeRadiationReactionModel::disabled;
    accuracy.computeOutwardFlux=true;
    accuracy.useRetardedExternalForces=retarded;
    ClassicalTrajectoryEngine engine(state,accuracy);
    const Ledgers start=ledgers(state);
    Ledgers lowest=start,highest=start;
    const int steps=orbits*perOrbit;
    const double step=period/perOrbit;
    int completed=0;
    bool ok=true;
    Ledgers last=start;
    for(int index=0;index<steps;++index) {
        if(!engine.advance(state,step)) { ok=false; break; }
        last=ledgers(state);
        lowest.e0=std::min(lowest.e0,last.e0);
        lowest.e1=std::min(lowest.e1,last.e1);
        lowest.e2=std::min(lowest.e2,last.e2);
        lowest.e3=std::min(lowest.e3,last.e3);
        highest.e0=std::max(highest.e0,last.e0);
        highest.e1=std::max(highest.e1,last.e1);
        highest.e2=std::max(highest.e2,last.e2);
        highest.e3=std::max(highest.e3,last.e3);
        ++completed;
    }
    // radiatedEnergy MERGES the channels and M1 dominates it near the barrier
    // while not recoiling the orbit at all (see state.hpp); orbitalRadiatedEnergy
    // is the exact retarded far-zone Poynting integral of the charge sector,
    // which is the part that must come out of the orbital energy.
    std::printf("radiated over the run: merged %.6e J = %.6f |U_dd|, "
                "orbital (far-zone E1) %.6e J = %.6f |U_dd|\n",
                state.radiatedEnergy,state.radiatedEnergy/dipoleScale,
                state.orbitalRadiatedEnergy,
                state.orbitalRadiatedEnergy/dipoleScale);
    std::printf("r %.2f r*  channel %d  tilt %.0f deg  orbits %d  tol %.0e  "
                "steps %d/%d %s  forces %s\n|U_dd| = %.6e J\n",
                r/comptonBarrierRadius,channel,tilt*180.0/pi,orbits,tolerance,
                completed,steps,ok?"ok":"FAILED",
                retarded?"retarded":"instantaneous",dipoleScale);
    std::printf("%42s %14s %14s %10s\n","ledger","range/|U_dd|",
                "drift/|U_dd|","vs E0");
    const double range0=(highest.e0-lowest.e0)/dipoleScale;
    struct Row { const char* name; double low,high,first,final; };
    for(const Row& entry:{
            Row{"E0  the model's ledger",lowest.e0,highest.e0,start.e0,last.e0},
            Row{"E1  minus the charge-dipole term",
                lowest.e1,highest.e1,start.e1,last.e1},
            Row{"E2  and the Darwin sign flipped",
                lowest.e2,highest.e2,start.e2,last.e2},
            Row{"E3  and U_dd on the proper moments",
                lowest.e3,highest.e3,start.e3,last.e3}}) {
        const double range=(entry.high-entry.low)/dipoleScale;
        const double drift=(entry.final-entry.first)/dipoleScale;
        std::printf("%42s %14.6e %+14.6e %10.4f\n",entry.name,range,drift,
                    range0>0.0?range/range0:0.0);
    }
}
