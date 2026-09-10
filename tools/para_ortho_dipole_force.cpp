// Para/ortho separation in the MUTUAL DIPOLE FORCE, measured directly on
// prescribed states rather than inferred from a collapse-time ensemble.
//
// Why this probe exists.  Para and ortho differ in exactly one place at
// preparation: the quantized spin makes the two magnetic moments exactly
// ALIGNED for para (S=0) and exactly ANTI-aligned for ortho (S=1) -- see the
// spin-quantization branch in crem_trajectory.hpp.  Every channel difference
// the model produces therefore has to be carried by something that reads that
// relative sign, and the mutual dipole force is the leading candidate.
//
// covariantDipoleGradientForce gained a term on 2026-09-10 (the relativistic
// dipole-force fix): the lab three-force is now
//
//     F = grad(U)/gamma + gamma v (partial_t U + v.grad U)/c^2,
//
// where U = mu_lab.B + p_lab.E.  The second piece is new, is proportional to
// the target's own velocity, and is the part the old transverse-orbit test
// geometry could not see.  This probe reports the two pieces separately for
// both channels, so "what did the fix change, and does it change para and
// ortho differently" is a measurement rather than an argument.
//
// Build from the repository root:
// g++ -std=c++20 -O2 -I . tools/para_ortho_dipole_force.cpp -o /tmp/probe
#include "modules/electrodynamics.hpp"
#include <iostream>
#include <iomanip>

// Circular-orbit pair at separation R, moments along z.  sign=+1 is para
// (aligned moments), sign=-1 is ortho (anti-aligned) -- the same convention
// the sampler uses for sampledScenario 2 and 3.
State pairAt(double time, double radius, double sign) {
    const double relativeSpeed =
        std::sqrt(pairCoulombStrength/(pairReducedMass*radius));
    const double firstShare = secondMass/(firstMass+secondMass);
    const double secondShare = -firstMass/(firstMass+secondMass);
    // A GENUINE circular orbit, not a straight tangent.  This matters: the
    // added term is gamma v (partial_t U + v.grad U)/c^2, and on a straight
    // worldline with a straight source history the field the target moves
    // through barely changes, so DU/Dt is suppressed by the prescription
    // rather than by the physics.  Curving both worldlines lets the field
    // rotate at the orbital rate the way it actually does.
    const double angularRate = relativeSpeed/radius;
    const double angle = angularRate*time;
    const Vec3 separationDirection{std::cos(angle), std::sin(angle), 0.0};
    const Vec3 tangentDirection{-std::sin(angle), std::cos(angle), 0.0};
    State s;
    s.time = time;
    s.firstPosition  = separationDirection*(radius*firstShare);
    s.secondPosition = separationDirection*(radius*secondShare);
    s.firstVelocity  = tangentDirection*(relativeSpeed*firstShare);
    s.secondVelocity = tangentDirection*(relativeSpeed*secondShare);
    s.firstProperDipole  = {0, 0, firstMagneticMoment};
    s.secondProperDipole = {0, 0, sign*secondMagneticMoment};
    synchronizeCovariantDipoles(s);
    return s;
}

int main() {
    std::cout<<std::setprecision(10);
    const double aPair = pairBohrRadius(activePair);
    std::cout<<"pair "<<firstSpecies.name<<"+"<<secondSpecies.name
             <<"  a_pair="<<aPair<<" m\n"
             <<"F_grad = grad(U)/gamma (the pre-fix force)\n"
             <<"F_mat  = gamma v (DU/Dt)/c^2 (the term the fix added)\n\n";
    std::cout<<"  R/a_pair   beta_target    channel      F_grad_x"
             <<"          F_mat_y           |F|\n";
    // Down to the two limits the model actually stops at: the collision
    // boundary at 0.005 a_pair and the Compton barrier at 193.3 fm, which is
    // 1.83e-3 a_pair.  Nothing in a production inspiral goes deeper, so the
    // last row bounds the added term over the WHOLE operating range.
    for(double fraction : {1.0, 0.5, 0.2, 0.05, 0.01, 0.005, 0.00183}) {
        const double radius = fraction*aPair;
        Vec3 paraForce, orthoForce;
        double paraTotal=0, orthoTotal=0, paraMat=0, orthoMat=0;
        for(double sign : {1.0, -1.0}) {
            State s = pairAt(0.0, radius, sign);
            StateHistory h;
            // Span the retarded reach of the pair at this separation.
            const double step = 0.02*radius/c;
            for(int k=256;k>=0;--k) h.push_back(pairAt(-k*step,radius,sign));
            const Vec3 total = covariantDipoleGradientForce(s,h,true);
            const double rate = dipoleCouplingMaterialRate(s,h,true);
            const Vec3 material =
                s.firstVelocity*(gamma(s.firstVelocity)*rate/(c*c));
            const Vec3 gradientPart = total-material;
            // The TARGET's own speed, which is what the added term carries;
            // the relative speed is twice this for equal masses.
            std::cout<<std::setw(10)<<fraction
                     <<std::setw(12)<<s.firstVelocity.norm()/c
                     <<std::setw(12)<<(sign>0?"para":"ortho")
                     <<std::setw(18)<<gradientPart.x
                     <<std::setw(18)<<material.y
                     <<std::setw(18)<<total.norm()<<"\n";
            if(sign>0) { paraForce=total; paraTotal=total.norm();
                         paraMat=material.norm(); }
            else       { orthoForce=total; orthoTotal=total.norm();
                         orthoMat=material.norm(); }
        }
        const double mean=0.5*(paraTotal+orthoTotal);
        // The dipole force in context.  Within the dipole sector the two
        // channels differ completely, but that sector is itself a small
        // fraction of the Coulomb force that actually drives the orbit, so
        // both statements have to be quoted together or either one misleads.
        const double coulombForce=pairCoulombStrength/(radius*radius);
        // VECTOR difference, not the gap between the two magnitudes.  At the
        // barrier the two channels have nearly equal |F| but point in
        // OPPOSITE directions, which a magnitude comparison reports as
        // agreement.
        std::cout<<"           separation |F_para-F_ortho|/mean = "
                 <<(mean>0?(paraForce-orthoForce).norm()/mean:0.0)
                 <<"   material-term share para/ortho = "
                 <<(paraTotal>0?paraMat/paraTotal:0.0)<<" / "
                 <<(orthoTotal>0?orthoMat/orthoTotal:0.0)<<"\n"
                 <<"           |F_dipole|/|F_Coulomb| para/ortho = "
                 <<paraTotal/coulombForce<<" / "<<orthoTotal/coulombForce
                 <<"\n\n";
    }
}
