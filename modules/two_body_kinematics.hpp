#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "physical_constants.hpp"
#include "vector3.hpp"

namespace positronium::kinematics {

using objects::Vec3;

inline constexpr double kSpeedOfLight = parameters::speedOfLight;

inline double invalidScalar() {
    return std::numeric_limits<double>::quiet_NaN();
}

inline Vec3 invalidVector() {
    const double value = invalidScalar();
    return {value, value, value};
}

inline bool finiteVector(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y)
        && std::isfinite(value.z);
}

inline long double strictLorentzFactorLongDouble(const Vec3& velocity) {
    const long double inverseC = 1.0L/static_cast<long double>(kSpeedOfLight);
    const long double betaX = static_cast<long double>(velocity.x)*inverseC;
    const long double betaY = static_cast<long double>(velocity.y)*inverseC;
    const long double betaZ = static_cast<long double>(velocity.z)*inverseC;
    const long double betaSquared = betaX*betaX + betaY*betaY + betaZ*betaZ;
    if (!(betaSquared >= 0.0L && betaSquared < 1.0L)) {
        return std::numeric_limits<long double>::quiet_NaN();
    }
    return 1.0L/std::sqrt(1.0L-betaSquared);
}

inline double strictLorentzFactor(const Vec3& velocity) {
    return static_cast<double>(strictLorentzFactorLongDouble(velocity));
}

inline bool subluminal(const Vec3& velocity) {
    return std::isfinite(strictLorentzFactor(velocity));
}

inline double momentumFromKineticEnergy(double kineticEnergy, double mass) {
    if (!(kineticEnergy >= 0.0) || !std::isfinite(kineticEnergy)
        || !(mass > 0.0) || !std::isfinite(mass)) {
        return invalidScalar();
    }
    const double restEnergy = mass*kSpeedOfLight*kSpeedOfLight;
    return std::sqrt(kineticEnergy)
        *std::sqrt(kineticEnergy + 2.0*restEnergy)/kSpeedOfLight;
}

inline double speedFromMomentum(double momentumMagnitude, double mass) {
    if (!(momentumMagnitude >= 0.0) || !std::isfinite(momentumMagnitude)
        || !(mass > 0.0) || !std::isfinite(mass)) {
        return invalidScalar();
    }
    const double restEnergy = mass*kSpeedOfLight*kSpeedOfLight;
    const double momentumEnergy = momentumMagnitude*kSpeedOfLight;
    const double totalEnergy = std::hypot(restEnergy, momentumEnergy);
    return kSpeedOfLight*momentumEnergy/totalEnergy;
}

inline double speedFromKineticEnergy(double kineticEnergy, double mass) {
    return speedFromMomentum(momentumFromKineticEnergy(kineticEnergy, mass), mass);
}

inline double kineticEnergyFromVelocity(const Vec3& velocity,double mass) {
    if (!(mass > 0.0) || !std::isfinite(mass)) return invalidScalar();
    const long double lorentzFactor=strictLorentzFactorLongDouble(velocity);
    if (!std::isfinite(lorentzFactor)) return invalidScalar();
    const long double inverseC=1.0L/static_cast<long double>(kSpeedOfLight);
    const long double betaX=static_cast<long double>(velocity.x)*inverseC;
    const long double betaY=static_cast<long double>(velocity.y)*inverseC;
    const long double betaZ=static_cast<long double>(velocity.z)*inverseC;
    const long double betaSquared=betaX*betaX+betaY*betaY+betaZ*betaZ;
    // gamma-1 = beta^2 gamma^2/(gamma+1) avoids subtracting one from gamma.
    const long double gammaMinusOne=betaSquared*lorentzFactor*lorentzFactor
        /(lorentzFactor+1.0L);
    const long double kineticEnergy=gammaMinusOne
        *static_cast<long double>(mass)*kSpeedOfLight*kSpeedOfLight;
    if (!(kineticEnergy >= 0.0L) || !std::isfinite(kineticEnergy)) {
        return invalidScalar();
    }
    return static_cast<double>(kineticEnergy);
}

struct ParticleFourMomentum {
    double energy = invalidScalar();
    Vec3 momentum = invalidVector();

    bool valid() const {
        return energy > 0.0 && std::isfinite(energy) && finiteVector(momentum);
    }
};

inline ParticleFourMomentum fourMomentumFromMomentum(
    const Vec3& momentum, double mass) {
    if (!finiteVector(momentum) || !(mass > 0.0) || !std::isfinite(mass)) {
        return {};
    }
    const double restEnergy = mass*kSpeedOfLight*kSpeedOfLight;
    const double momentumEnergy = std::hypot(
        std::hypot(momentum.x*kSpeedOfLight, momentum.y*kSpeedOfLight),
        momentum.z*kSpeedOfLight);
    const double energy = std::hypot(restEnergy, momentumEnergy);
    if (!std::isfinite(energy)) return {};
    return {energy, momentum};
}

inline ParticleFourMomentum fourMomentumFromVelocity(
    const Vec3& velocity, double mass) {
    const double lorentzFactor = strictLorentzFactor(velocity);
    if (!std::isfinite(lorentzFactor) || !(mass > 0.0) || !std::isfinite(mass)) {
        return {};
    }
    const double energy = lorentzFactor*mass*kSpeedOfLight*kSpeedOfLight;
    const Vec3 momentum = velocity*(lorentzFactor*mass);
    if (!std::isfinite(energy) || !finiteVector(momentum)) return {};
    return {energy, momentum};
}

inline Vec3 velocityFromFourMomentum(const ParticleFourMomentum& value) {
    if (!value.valid()) return invalidVector();
    const Vec3 velocity = value.momentum
        *(kSpeedOfLight*kSpeedOfLight/value.energy);
    return subluminal(velocity) ? velocity : invalidVector();
}

// Active Lorentz boost from a centre-of-momentum frame to a frame in which
// that centre moves with frameVelocity.  Boosting four-momentum instead of
// adding three-velocities preserves the mass shell and cannot create a
// superluminal particle from physical inputs.
inline ParticleFourMomentum boostFourMomentum(
    const ParticleFourMomentum& value, const Vec3& frameVelocity) {
    if (!value.valid() || !subluminal(frameVelocity)) return {};
    if (frameVelocity.squaredNorm() == 0.0) return value;
    const Vec3 beta = frameVelocity/kSpeedOfLight;
    const double boostGamma = strictLorentzFactor(frameVelocity);
    const double betaMomentum = beta.x*value.momentum.x
        + beta.y*value.momentum.y + beta.z*value.momentum.z;
    const double boostedEnergy = boostGamma
        *(value.energy + kSpeedOfLight*betaMomentum);
    const double longitudinalCoefficient = boostGamma*value.energy/kSpeedOfLight
        + boostGamma*boostGamma/(boostGamma+1.0)*betaMomentum;
    const Vec3 boostedMomentum = value.momentum + beta*longitudinalCoefficient;
    if (!(boostedEnergy > 0.0) || !std::isfinite(boostedEnergy)
        || !finiteVector(boostedMomentum)) {
        return {};
    }
    return {boostedEnergy, boostedMomentum};
}

inline Vec3 boostVelocity(const Vec3& velocity, double mass,
                          const Vec3& frameVelocity) {
    return velocityFromFourMomentum(boostFourMomentum(
        fourMomentumFromVelocity(velocity, mass), frameVelocity));
}

struct CentreOfMomentumKinematics {
    double kineticEnergy = invalidScalar();
    double totalEnergy = invalidScalar();
    double momentumMagnitude = invalidScalar();
    double firstEnergy = invalidScalar();
    double secondEnergy = invalidScalar();
    double firstSpeed = invalidScalar();
    double secondSpeed = invalidScalar();
    double relativeSpeed = invalidScalar();

    bool valid() const {
        return kineticEnergy >= 0.0 && totalEnergy > 0.0
            && momentumMagnitude >= 0.0 && firstEnergy > 0.0
            && secondEnergy > 0.0 && firstSpeed >= 0.0
            && secondSpeed >= 0.0 && relativeSpeed >= 0.0
            && std::isfinite(kineticEnergy) && std::isfinite(totalEnergy)
            && std::isfinite(momentumMagnitude) && std::isfinite(firstEnergy)
            && std::isfinite(secondEnergy) && std::isfinite(firstSpeed)
            && std::isfinite(secondSpeed) && std::isfinite(relativeSpeed)
            && firstSpeed < kSpeedOfLight && secondSpeed < kSpeedOfLight;
    }
};

// Exact two-body centre-of-momentum kinematics.  The factored form avoids the
// loss of precision in W^2-(m1+m2)^2 c^4 when K is measured in eV and one of
// the rest energies is hundreds of MeV.
inline CentreOfMomentumKinematics centreOfMomentumKinematics(
    double kineticEnergy, double firstMass, double secondMass) {
    CentreOfMomentumKinematics result;
    if (!(kineticEnergy >= 0.0) || !std::isfinite(kineticEnergy)
        || !(firstMass > 0.0) || !(secondMass > 0.0)
        || !std::isfinite(firstMass) || !std::isfinite(secondMass)) {
        return result;
    }
    const double firstRestEnergy = firstMass*kSpeedOfLight*kSpeedOfLight;
    const double secondRestEnergy = secondMass*kSpeedOfLight*kSpeedOfLight;
    const double restEnergy = firstRestEnergy + secondRestEnergy;
    const double totalEnergy = restEnergy + kineticEnergy;
    if (!std::isfinite(totalEnergy) || !(totalEnergy > 0.0)) return result;

    // sqrt(A*B)/(2 W c), evaluated as sqrt(A/W)*sqrt(B/W)/(2 c)
    // so neither the eV limit nor a large finite K squares an extreme number.
    const double firstReducedFactor = kineticEnergy
        *((kineticEnergy + 2.0*restEnergy)/totalEnergy);
    const double secondReducedFactor = (kineticEnergy + 2.0*firstRestEnergy)
        *((kineticEnergy + 2.0*secondRestEnergy)/totalEnergy);
    if (!(firstReducedFactor >= 0.0) || !(secondReducedFactor >= 0.0)
        || !std::isfinite(firstReducedFactor)
        || !std::isfinite(secondReducedFactor)) {
        return result;
    }
    const double momentumMagnitude = 0.5*std::sqrt(firstReducedFactor)
        *std::sqrt(secondReducedFactor)/kSpeedOfLight;
    const double momentumEnergy = momentumMagnitude*kSpeedOfLight;
    const double firstEnergy = std::hypot(firstRestEnergy, momentumEnergy);
    const double secondEnergy = std::hypot(secondRestEnergy, momentumEnergy);
    const double firstSpeed = kSpeedOfLight*momentumEnergy/firstEnergy;
    const double secondSpeed = kSpeedOfLight*momentumEnergy/secondEnergy;
    result = {kineticEnergy, totalEnergy, momentumMagnitude,
              firstEnergy, secondEnergy, firstSpeed, secondSpeed,
              firstSpeed+secondSpeed};
    return result.valid() ? result : CentreOfMomentumKinematics{};
}

struct FreeTwoBodyKinematics {
    double laboratoryEnergy = invalidScalar();
    double laboratoryKineticEnergy = invalidScalar();
    Vec3 laboratoryMomentum = invalidVector();
    Vec3 centreOfMomentumVelocity = invalidVector();
    double invariantEnergy = invalidScalar();
    double centreOfMomentumKineticEnergy = invalidScalar();

    bool valid() const {
        return laboratoryEnergy > 0.0 && laboratoryKineticEnergy >= 0.0
            && invariantEnergy > 0.0 && centreOfMomentumKineticEnergy >= 0.0
            && std::isfinite(laboratoryEnergy)
            && std::isfinite(laboratoryKineticEnergy)
            && finiteVector(laboratoryMomentum)
            && finiteVector(centreOfMomentumVelocity)
            && std::isfinite(invariantEnergy)
            && std::isfinite(centreOfMomentumKineticEnergy)
            && subluminal(centreOfMomentumVelocity);
    }
};

// Reconstruct the free-pair invariant from two simultaneous velocities.  The
// relative Lorentz factor is used instead of subtracting E_lab^2-c^2 P_lab^2;
// this retains the small internal energy of a highly boosted heavy pair.
inline FreeTwoBodyKinematics freeTwoBodyKinematics(
    const Vec3& firstVelocity, double firstMass,
    const Vec3& secondVelocity, double secondMass) {
    FreeTwoBodyKinematics result;
    if (!(firstMass > 0.0) || !(secondMass > 0.0)
        || !std::isfinite(firstMass) || !std::isfinite(secondMass)) {
        return result;
    }
    const long double firstGamma = strictLorentzFactorLongDouble(firstVelocity);
    const long double secondGamma = strictLorentzFactorLongDouble(secondVelocity);
    if (!std::isfinite(firstGamma) || !std::isfinite(secondGamma)) {
        return result;
    }
    const long double c = kSpeedOfLight;
    const long double inverseC=1.0L/c;
    const long double firstBetaX=static_cast<long double>(firstVelocity.x)*inverseC;
    const long double firstBetaY=static_cast<long double>(firstVelocity.y)*inverseC;
    const long double firstBetaZ=static_cast<long double>(firstVelocity.z)*inverseC;
    const long double secondBetaX=static_cast<long double>(secondVelocity.x)*inverseC;
    const long double secondBetaY=static_cast<long double>(secondVelocity.y)*inverseC;
    const long double secondBetaZ=static_cast<long double>(secondVelocity.z)*inverseC;
    const long double deltaX=firstBetaX-secondBetaX;
    const long double deltaY=firstBetaY-secondBetaY;
    const long double deltaZ=firstBetaZ-secondBetaZ;
    const long double deltaSquared=deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ;
    const long double deltaDotSecond=deltaX*secondBetaX
        +deltaY*secondBetaY+deltaZ*secondBetaZ;

    // gamma_rel-1 must vanish exactly for two equal velocities, even when both
    // are extremely close to c.  Forming gamma1*gamma2*(1-beta1.beta2)-1
    // subtracts two nearly equal numbers and can invent eV of internal energy
    // for that state.  The positive identity below first computes
    // q=gamma_rel^2-1 and then rationalizes sqrt(1+q)-1:
    //
    // q = gamma1^2 (|beta1-beta2|^2
    //                  + gamma2^2 ((beta1-beta2).beta2)^2).
    const long double relativeGammaSquaredMinusOne=firstGamma*firstGamma
        *(deltaSquared+secondGamma*secondGamma
          *deltaDotSecond*deltaDotSecond);
    if (!(relativeGammaSquaredMinusOne >= 0.0L)
        || !std::isfinite(relativeGammaSquaredMinusOne)) {
        return result;
    }
    const long double relativeGammaMinusOne=relativeGammaSquaredMinusOne
        /(std::sqrt(1.0L+relativeGammaSquaredMinusOne)+1.0L);

    const long double firstMassLong = firstMass;
    const long double secondMassLong = secondMass;
    const long double massSum = firstMassLong+secondMassLong;
    const long double invariantMassIncrementSquared =
        2.0L*firstMassLong*secondMassLong*relativeGammaMinusOne;
    const long double invariantMass = std::sqrt(
        massSum*massSum+invariantMassIncrementSquared);
    const long double internalMass = invariantMassIncrementSquared
        /(invariantMass+massSum);
    const long double internalKineticEnergy = internalMass*c*c;

    const long double firstBetaSquared=firstBetaX*firstBetaX
        +firstBetaY*firstBetaY+firstBetaZ*firstBetaZ;
    const long double secondBetaSquared=secondBetaX*secondBetaX
        +secondBetaY*secondBetaY+secondBetaZ*secondBetaZ;
    const long double firstGammaMinusOne = firstBetaSquared*firstGamma*firstGamma
        /(firstGamma+1.0L);
    const long double secondGammaMinusOne = secondBetaSquared*secondGamma*secondGamma
        /(secondGamma+1.0L);
    const long double firstRestEnergy = firstMassLong*c*c;
    const long double secondRestEnergy = secondMassLong*c*c;
    const long double laboratoryKineticEnergy =
        firstGammaMinusOne*firstRestEnergy
        +secondGammaMinusOne*secondRestEnergy;
    const long double laboratoryEnergy = firstRestEnergy+secondRestEnergy
        +laboratoryKineticEnergy;
    const Vec3 laboratoryMomentum = firstVelocity
        *static_cast<double>(firstGamma*firstMassLong)
        +secondVelocity*static_cast<double>(secondGamma*secondMassLong);
    const Vec3 centreOfMomentumVelocity = laboratoryMomentum
        *static_cast<double>(c*c/laboratoryEnergy);
    const long double invariantEnergy = massSum*c*c+internalKineticEnergy;

    result = {static_cast<double>(laboratoryEnergy),
              static_cast<double>(laboratoryKineticEnergy),
              laboratoryMomentum, centreOfMomentumVelocity,
              static_cast<double>(invariantEnergy),
              static_cast<double>(internalKineticEnergy)};
    return result.valid() ? result : FreeTwoBodyKinematics{};
}

struct HeadOnLabKinematics {
    ParticleFourMomentum first;
    ParticleFourMomentum second;
    FreeTwoBodyKinematics pair;

    bool valid() const { return first.valid() && second.valid() && pair.valid(); }
};

inline HeadOnLabKinematics headOnLabKinematics(
    double firstKineticEnergy, double firstMass,
    double secondKineticEnergy, double secondMass,
    const Vec3& firstDirection) {
    HeadOnLabKinematics result;
    if (!(firstKineticEnergy >= 0.0) || !(secondKineticEnergy >= 0.0)
        || !std::isfinite(firstKineticEnergy)
        || !std::isfinite(secondKineticEnergy)
        || !(firstMass > 0.0) || !(secondMass > 0.0)
        || !std::isfinite(firstMass) || !std::isfinite(secondMass)) {
        return result;
    }
    const double directionNorm = firstDirection.norm();
    if (!(directionNorm > 0.0) || !std::isfinite(directionNorm)) return result;
    const Vec3 direction = firstDirection/directionNorm;
    const double firstMomentum = momentumFromKineticEnergy(
        firstKineticEnergy,firstMass);
    const double secondMomentum = momentumFromKineticEnergy(
        secondKineticEnergy,secondMass);
    if (!std::isfinite(firstMomentum) || !std::isfinite(secondMomentum)) {
        return result;
    }
    result.first = fourMomentumFromMomentum(direction*firstMomentum,firstMass);
    result.second = fourMomentumFromMomentum(direction*(-secondMomentum),secondMass);
    if (!result.first.valid() || !result.second.valid()) return {};

    // For opposing lab beams the invariant can be expanded without subtracting
    // rest-energy squares.  This retains all input digits when, for example,
    // an eV electron beam is paired with a proton rest energy near one GeV.
    const long double c = kSpeedOfLight;
    const long double firstRest = static_cast<long double>(firstMass)*c*c;
    const long double secondRest = static_cast<long double>(secondMass)*c*c;
    const long double firstKinetic = firstKineticEnergy;
    const long double secondKinetic = secondKineticEnergy;
    const long double firstMomentumEnergy =
        static_cast<long double>(firstMomentum)*c;
    const long double secondMomentumEnergy =
        static_cast<long double>(secondMomentum)*c;
    const long double invariantSquaredIncrement = 2.0L*(
        firstRest*secondKinetic + secondRest*firstKinetic
        +firstKinetic*secondKinetic
        +firstMomentumEnergy*secondMomentumEnergy);
    const long double restEnergy = firstRest+secondRest;
    const long double invariantEnergy = std::sqrt(
        restEnergy*restEnergy+invariantSquaredIncrement);
    const long double internalKinetic = invariantSquaredIncrement
        /(invariantEnergy+restEnergy);
    const long double laboratoryKinetic = firstKinetic+secondKinetic;
    const long double laboratoryEnergy = restEnergy+laboratoryKinetic;
    const Vec3 laboratoryMomentum = result.first.momentum+result.second.momentum;
    const Vec3 centreVelocity = laboratoryMomentum
        *static_cast<double>(c*c/laboratoryEnergy);
    result.pair = {static_cast<double>(laboratoryEnergy),
                   static_cast<double>(laboratoryKinetic),
                   laboratoryMomentum,centreVelocity,
                   static_cast<double>(restEnergy+internalKinetic),
                   static_cast<double>(internalKinetic)};
    return result.valid() ? result : HeadOnLabKinematics{};
}

struct IncomingTwoBodyKinematics {
    CentreOfMomentumKinematics asymptotic;
    CentreOfMomentumKinematics finiteRadius;
    double tangentialFraction = invalidScalar();
    Vec3 momentumDirection = invalidVector();
    ParticleFourMomentum firstCentreOfMomentum;
    ParticleFourMomentum secondCentreOfMomentum;
    ParticleFourMomentum firstFrame;
    ParticleFourMomentum secondFrame;
    Vec3 firstVelocity = invalidVector();
    Vec3 secondVelocity = invalidVector();

    bool valid() const {
        return asymptotic.valid() && finiteRadius.valid()
            && tangentialFraction >= 0.0 && tangentialFraction <= 1.0
            && std::isfinite(tangentialFraction)
            && finiteVector(momentumDirection)
            && firstCentreOfMomentum.valid() && secondCentreOfMomentum.valid()
            && firstFrame.valid() && secondFrame.valid()
            && subluminal(firstVelocity) && subluminal(secondVelocity);
    }
};

inline IncomingTwoBodyKinematics incomingTwoBodyKinematics(
    double asymptoticKineticEnergy, double potentialEnergyGain,
    double impactParameter, double matchingRadius,
    const Vec3& radialDirection, const Vec3& tangentDirection,
    double firstMass, double secondMass,
    const Vec3& frameVelocity = {}) {
    IncomingTwoBodyKinematics result;
    if (!(asymptoticKineticEnergy >= 0.0)
        || !(potentialEnergyGain >= 0.0)
        || !(impactParameter >= 0.0) || !(matchingRadius > 0.0)
        || impactParameter > matchingRadius
        || !std::isfinite(asymptoticKineticEnergy)
        || !std::isfinite(potentialEnergyGain)
        || !std::isfinite(impactParameter) || !std::isfinite(matchingRadius)
        || !subluminal(frameVelocity)) {
        return result;
    }
    result.asymptotic = centreOfMomentumKinematics(
        asymptoticKineticEnergy,firstMass,secondMass);
    result.finiteRadius = centreOfMomentumKinematics(
        asymptoticKineticEnergy+potentialEnergyGain,firstMass,secondMass);
    if (!result.asymptotic.valid() || !result.finiteRadius.valid()
        || !(result.finiteRadius.momentumMagnitude > 0.0)) {
        return IncomingTwoBodyKinematics{};
    }

    double fraction = impactParameter*result.asymptotic.momentumMagnitude
        /(matchingRadius*result.finiteRadius.momentumMagnitude);
    const double fractionTolerance = 64.0*std::numeric_limits<double>::epsilon();
    if (fraction > 1.0 && fraction <= 1.0+fractionTolerance) fraction = 1.0;
    if (!(fraction >= 0.0 && fraction <= 1.0) || !std::isfinite(fraction)) {
        return IncomingTwoBodyKinematics{};
    }

    const double radialNorm = radialDirection.norm();
    if (!(radialNorm > 0.0) || !std::isfinite(radialNorm)) {
        return IncomingTwoBodyKinematics{};
    }
    const Vec3 radial = radialDirection/radialNorm;
    const Vec3 tangentProjection = tangentDirection
        -radial*(radial.x*tangentDirection.x
               +radial.y*tangentDirection.y
               +radial.z*tangentDirection.z);
    const double tangentNorm = tangentProjection.norm();
    if (!(tangentNorm > 0.0) || !std::isfinite(tangentNorm)) {
        return IncomingTwoBodyKinematics{};
    }
    const Vec3 tangent = tangentProjection/tangentNorm;
    const Vec3 direction = radial*(-std::sqrt(std::max(0.0,1.0-fraction*fraction)))
                         +tangent*fraction;
    const Vec3 firstMomentum = direction*result.finiteRadius.momentumMagnitude;
    result.tangentialFraction = fraction;
    result.momentumDirection = direction;
    result.firstCentreOfMomentum = fourMomentumFromMomentum(firstMomentum,firstMass);
    result.secondCentreOfMomentum = fourMomentumFromMomentum(
        firstMomentum*(-1.0),secondMass);
    result.firstFrame = boostFourMomentum(
        result.firstCentreOfMomentum,frameVelocity);
    result.secondFrame = boostFourMomentum(
        result.secondCentreOfMomentum,frameVelocity);
    result.firstVelocity = velocityFromFourMomentum(result.firstFrame);
    result.secondVelocity = velocityFromFourMomentum(result.secondFrame);
    return result.valid() ? result : IncomingTwoBodyKinematics{};
}

} // namespace positronium::kinematics
