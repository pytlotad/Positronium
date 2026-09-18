// The model's force on a moment against the exact Amperian force, in the one
// case where everything is analytic (audit section 125).
//
// A charge in UNIFORM motion past a magnetic moment at rest: the fields are
// the Heaviside ellipsoid in closed form, and the force on an Amperian dipole
// is F = grad(m.B) - (1/c^2) m x dE/dt, both computable without the model.
// Audit 122 showed the model's momentum imbalance is first order in v/c and
// audit 124 left 31% of it unexplained; this is the setting where a
// first-order error in the moment force itself must show, because the rigid
// drift case (zero relative velocity) already balances exactly.
//
// Usage: moving_charge_dipole_force [separation in r*] [beta] [tilt deg]
//        [history nodes] [history span in r/c]
// Build from the repository root:
// g++ -std=c++20 -O2 -I . $(root-config --cflags)
//     tools/moving_charge_dipole_force.cpp -o /tmp/moving $(root-config --libs)
#include "modules/crem_trajectory.hpp"
#include <cstdio>
#include <cstdlib>

namespace {
// Exact field of a charge in uniform motion, at the field point, from the
// charge's PRESENT position (Heaviside ellipsoid).
ElectromagneticField uniformChargeField(const Vec3& fieldPoint,
                                        const Vec3& chargePosition,
                                        const Vec3& velocity,double charge) {
    const Vec3 separationVector=fieldPoint-chargePosition;
    const double distance=separationVector.norm();
    if(!(distance>0.0)) return {};
    const double betaSquared=velocity.squaredNorm()/(c*c);
    const Vec3 direction=separationVector/distance;
    const double cosine=velocity.norm()>0.0
        ?dot(direction,velocity)/velocity.norm():0.0;
    const double sineSquared=std::max(0.0,1.0-cosine*cosine);
    const double denominator=std::pow(1.0-betaSquared*sineSquared,1.5);
    const Vec3 electric=separationVector*(charge*(1.0-betaSquared)
        /(4.0*pi*epsilon0*distance*distance*distance*denominator));
    return {electric,cross(velocity,electric)*(1.0/(c*c))};
}
}

int main(int argc,char** argv) {
    const double separation=(argc>1?atof(argv[1]):20.0)*comptonBarrierRadius;
    const double beta=argc>2?atof(argv[2]):0.05;
    const double tilt=(argc>3?atof(argv[3]):0.0)*pi/180;
    const int nodes=argc>4?atoi(argv[4]):512;
    const double span=argc>5?atof(argv[5]):40.0;
    const Vec3 velocity{0.0,beta*c,0.0};
    const Vec3 moment=Vec3{std::sin(tilt),0.0,std::cos(tilt)}*firstMagneticMoment;
    // Particle 1 carries the moment and is at rest; particle 2 is the charge,
    // moving, and carries no moment at all.
    const auto build=[&](double time) {
        State s{};
        s.time=time;
        s.firstPosition={0.5*separation,0.0,0.0};
        s.secondPosition=Vec3{-0.5*separation,0.0,0.0}+velocity*time;
        s.firstVelocity={};
        s.secondVelocity=velocity;
        s.firstProperDipole=moment;
        s.secondProperDipole={};
        synchronizeCovariantDipoles(s);
        s.secondDipole={};
        s.secondElectricDipole={};
        return s;
    };
    StateHistory history;
    const double lightTime=separation/c;
    for(int index=0;index<nodes;++index)
        history.push_back(build(-span*lightTime
            *(1.0-static_cast<double>(index)/nodes)));
    const State present=build(0.0);
    // The model: the moment force, and the hidden-momentum rate beside it.
    const Vec3 modelGradient=covariantDipoleGradientForce(present,history,true);
    const Vec3 modelHidden=hiddenMomentumRateForce(present,history,true);
    // The exact Amperian force, by symmetric stencils on the closed-form
    // fields: grad(m.B) with m fixed, and m x dE/dt at the moment's position.
    const auto analytic=[&](double step) {
        const Vec3 point=present.firstPosition;
        Vec3 gradient;
        for(int axis=0;axis<3;++axis) {
            Vec3 offset;
            (axis==0?offset.x:axis==1?offset.y:offset.z)=step;
            const double plus=dot(moment,uniformChargeField(point+offset,
                present.secondPosition,velocity,secondCharge).magnetic);
            const double minus=dot(moment,uniformChargeField(point-offset,
                present.secondPosition,velocity,secondCharge).magnetic);
            (axis==0?gradient.x:axis==1?gradient.y:gradient.z)
                =(plus-minus)/(2.0*step);
        }
        const double timeStep=step/c;
        const Vec3 electricPlus=uniformChargeField(point,
            present.secondPosition+velocity*timeStep,velocity,
            secondCharge).electric;
        const Vec3 electricMinus=uniformChargeField(point,
            present.secondPosition-velocity*timeStep,velocity,
            secondCharge).electric;
        const Vec3 electricRate=(electricPlus-electricMinus)*(1.0/(2.0*timeStep));
        return gradient-cross(moment,electricRate)*(1.0/(c*c));
    };
    const Vec3 exactCoarse=analytic(1.0e-4*separation);
    const Vec3 exactFine=analytic(0.5e-4*separation);
    const Vec3 modelTotal=modelGradient+modelHidden;
    std::printf("separation %.2f r*  beta %.4f  tilt %.0f deg\n",
                separation/comptonBarrierRadius,beta,tilt*180/pi);
    std::printf("  stencil check: |coarse-fine|/|fine| %.3e\n",
                (exactCoarse-exactFine).norm()
                    /std::max(exactFine.norm(),1.0e-300));
    std::printf("  exact Amperian    %+.6e %+.6e %+.6e  |F| %.6e\n",
                exactFine.x,exactFine.y,exactFine.z,exactFine.norm());
    std::printf("  model gradient    %+.6e %+.6e %+.6e  |F| %.6e\n",
                modelGradient.x,modelGradient.y,modelGradient.z,
                modelGradient.norm());
    std::printf("  model hidden rate %+.6e %+.6e %+.6e  |F| %.6e\n",
                modelHidden.x,modelHidden.y,modelHidden.z,modelHidden.norm());
    std::printf("  model total       %+.6e %+.6e %+.6e  |F| %.6e\n",
                modelTotal.x,modelTotal.y,modelTotal.z,modelTotal.norm());
    // The reaction: the Lorentz force on the moving charge in the moment's
    // field.  The moment is at rest, so its exact field is the static dipole
    // one (the model softens it at the moment radius; at large separations
    // that is negligible, and the printed ratio shows when it is not).
    {
        const Vec3 chargePosition=present.secondPosition;
        const Vec3 displacement=chargePosition-present.firstPosition;
        const double distance=displacement.norm();
        const Vec3 unit=displacement/distance;
        const Vec3 exactMagnetic=(unit*(3.0*dot(unit,moment))-moment)
            *(mu0/(4.0*pi*distance*distance*distance));
        const Vec3 exactReaction=cross(velocity,exactMagnetic)*secondCharge;
        // Only the moment's own field: fieldFromOtherParticleAt would add the
        // partner's CHARGE field, i.e. the Coulomb force, which is not the
        // reaction being compared here.
        const ElectromagneticField modelField=retardedMagneticDipoleField(
            chargePosition,present.time,history,present,true);
        const Vec3 modelReaction=lorentzForce(secondCharge,velocity,modelField);
        std::printf("  exact reaction    %+.6e %+.6e %+.6e  |F| %.6e\n",
                    exactReaction.x,exactReaction.y,exactReaction.z,
                    exactReaction.norm());
        std::printf("  model reaction    %+.6e %+.6e %+.6e  |F| %.6e\n",
                    modelReaction.x,modelReaction.y,modelReaction.z,
                    modelReaction.norm());
        std::printf("  |model-exact| reaction / |exact| %.4e\n",
                    (modelReaction-exactReaction).norm()
                        /std::max(exactReaction.norm(),1.0e-300));
        std::printf("  net of the two, model %.6e N, exact %.6e N "
                    "(|net|/|action| model %.4e exact %.4e)\n",
                    (modelTotal+modelReaction).norm(),
                    (exactFine+exactReaction).norm(),
                    (modelTotal+modelReaction).norm()
                        /std::max(modelTotal.norm(),1.0e-300),
                    (exactFine+exactReaction).norm()
                        /std::max(exactFine.norm(),1.0e-300));
    }
    std::printf("  |model-exact|/|exact| %.4e   (gradient alone %.4e)\n",
                (modelTotal-exactFine).norm()/std::max(exactFine.norm(),1e-300),
                (modelGradient-exactFine).norm()
                    /std::max(exactFine.norm(),1.0e-300));
}
