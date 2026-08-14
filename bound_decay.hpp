#ifndef POSITRONIUM_BOUND_DECAY_HPP
#define POSITRONIUM_BOUND_DECAY_HPP

// Ideal-vacuum annihilation generator for ground-state positronium.
//
// This module is deliberately independent of the classical orbital model in
// positronium.cpp.  Positronium annihilation is a quantum process: its decay
// time and photon final state cannot be obtained from a classical trajectory
// reaching a short-distance cutoff.  The lifetimes below are phenomenological
// vacuum inputs, while the photon kinematics use the leading 2-gamma channel
// for para-positronium and the unpolarized Ore-Powell distribution for the
// leading 3-gamma channel of ortho-positronium.
//
// Detector resolution, material effects (pick-off/quenching), higher photon
// multiplicities, radiative corrections and hyperfine mass splitting are not
// included.  Energies and times are returned in SI units; convenience
// conversion constants are provided for histogram code.
//
// The 3-gamma sampling prescription follows the unpolarized at-rest model
// documented by Geant4 (Ore and Powell, Phys. Rev. 75, 1696 (1949)):
// https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/decay/OrthoPositronium.html

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>

namespace bound_decay {

inline constexpr double pi = 3.141592653589793238462643383279502884;
inline constexpr double speedOfLight = 299792458.0;
inline constexpr double electronMass = 9.1093837139e-31;
inline constexpr double electronVolt = 1.602176634e-19;
inline constexpr double hbar = 1.054571817e-34;
inline constexpr double epsilon0 = 8.8541878128e-12;
inline constexpr double elementaryCharge = 1.602176634e-19;
inline constexpr double coulombConstant = 1.0 / (4.0 * pi * epsilon0);
inline constexpr double fineStructureConstant =
    coulombConstant * elementaryCharge * elementaryCharge
    / (hbar * speedOfLight);

inline constexpr double paraLifetimeSeconds = 125.0e-12;
inline constexpr double orthoLifetimeSeconds = 142.0e-9;

// Ground-state positronium has reduced mass m_e/2, hence one half of the
// hydrogen binding energy.  The same spin-averaged mass is used for p-Ps and
// o-Ps here; their much smaller hyperfine mass difference is outside the
// accuracy of this ideal event generator.
inline constexpr double groundStateBindingEnergyJoules =
    fineStructureConstant * fineStructureConstant * electronMass
    * speedOfLight * speedOfLight / 4.0;
inline constexpr double positroniumRestEnergyJoules =
    2.0 * electronMass * speedOfLight * speedOfLight
    - groundStateBindingEnergyJoules;

inline constexpr double joulesPerKeV = 1.0e3 * electronVolt;

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3 operator+(const Vec3& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    Vec3 operator-() const { return {-x, -y, -z}; }
    Vec3 operator*(double scale) const {
        return {x * scale, y * scale, z * scale};
    }
};

inline double dot(const Vec3& first, const Vec3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

inline Vec3 cross(const Vec3& first, const Vec3& second) {
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

inline double norm(const Vec3& vector) {
    return std::sqrt(dot(vector, vector));
}

inline Vec3 normalized(const Vec3& vector) {
    const double magnitude = norm(vector);
    if (!(magnitude > 0.0)) throw std::runtime_error("cannot normalize zero vector");
    return vector * (1.0 / magnitude);
}

enum class PositroniumState { Para, Ortho };

struct Photon {
    double energyJoules = 0.0;
    Vec3 direction;

    // Photon momentum in kg m/s.
    Vec3 momentum() const {
        return direction * (energyJoules / speedOfLight);
    }
};

struct DecayEvent {
    PositroniumState state = PositroniumState::Para;
    double decayTimeSeconds = 0.0;
    std::array<Photon, 3> photons{};
    std::size_t photonCount = 0;

    // Directly usable histogram observables.  Observables which do not apply
    // to a decay channel contain quiet NaN rather than a fabricated value.
    double paraPhotonCosPolar = std::numeric_limits<double>::quiet_NaN();
    double orthoMaximumEnergyJoules = std::numeric_limits<double>::quiet_NaN();
    double orthoMiddleEnergyJoules = std::numeric_limits<double>::quiet_NaN();
    double orthoLeadingPairAngleRadians = std::numeric_limits<double>::quiet_NaN();
};

inline double energyKeV(double energyJoules) {
    return energyJoules / joulesPerKeV;
}

namespace detail {

inline double clampCosine(double value) {
    return std::clamp(value, -1.0, 1.0);
}

template <class RandomEngine>
double uniformUnit(RandomEngine& random) {
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(random);
}

template <class RandomEngine>
double uniformOpenUnit(RandomEngine& random) {
    // generate_canonical returns [0,1).  Excluding the exactly representable
    // zero keeps log(u) finite without changing the continuous distribution.
    double value = 0.0;
    do {
        value = uniformUnit(random);
    } while (!(value > 0.0));
    return value;
}

template <class RandomEngine>
double exponentialTime(RandomEngine& random, double meanLifetime) {
    return -meanLifetime * std::log(uniformOpenUnit(random));
}

template <class RandomEngine>
Vec3 isotropicDirection(RandomEngine& random) {
    const double cosine = 2.0 * uniformUnit(random) - 1.0;
    const double azimuth = 2.0 * pi * uniformUnit(random);
    const double transverse = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    return {transverse * std::cos(azimuth),
            transverse * std::sin(azimuth), cosine};
}

struct OrthoEnergyPoint {
    std::array<double, 3> fractions{};
    std::array<double, 3> pairCosines{}; // (12), (13), (23)
};

inline OrthoEnergyPoint makeOrthoEnergyPoint(double first, double second) {
    const double third = 2.0 - first - second;
    const double cosine12 = (third * third - first * first - second * second)
                          / (2.0 * first * second);
    const double cosine13 = (second * second - first * first - third * third)
                          / (2.0 * first * third);
    const double cosine23 = (first * first - second * second - third * third)
                          / (2.0 * second * third);
    return {{first, second, third},
            {clampCosine(cosine12), clampCosine(cosine13),
             clampCosine(cosine23)}};
}

inline double orePowellWeight(const OrthoEnergyPoint& point) {
    double weight = 0.0;
    for (double cosine : point.pairCosines) {
        const double difference = 1.0 - cosine;
        weight += difference * difference;
    }
    return weight;
}

inline double pairCosine(const OrthoEnergyPoint& point,
                         std::size_t first, std::size_t second) {
    if (first > second) std::swap(first, second);
    if (first == 0 && second == 1) return point.pairCosines[0];
    if (first == 0 && second == 2) return point.pairCosines[1];
    return point.pairCosines[2];
}

template <class RandomEngine>
OrthoEnergyPoint sampleOrePowellEnergies(RandomEngine& random) {
    // Geant4's unpolarized at-rest prescription samples x1,x2 uniformly in
    // [0,1], takes x3=2-x1-x2, rejects points outside the physical triangle,
    // and accepts the remainder with weight
    //   sum_ij (1-cos(theta_ij))^2.
    // Each term is at most four, so 12 is a conservative global envelope.
    constexpr double weightEnvelope = 12.0;
    constexpr std::size_t maximumAttempts = 1000000;
    for (std::size_t attempt = 0; attempt < maximumAttempts; ++attempt) {
        const double first = uniformOpenUnit(random);
        const double second = uniformOpenUnit(random);
        const double third = 2.0 - first - second;
        if (!(third > 0.0 && third <= 1.0)) continue;

        const OrthoEnergyPoint point = makeOrthoEnergyPoint(first, second);
        const double acceptance = orePowellWeight(point) / weightEnvelope;
        if (uniformUnit(random) < acceptance) return point;
    }
    throw std::runtime_error("Ore-Powell rejection sampler did not converge");
}

template <class RandomEngine>
std::array<Vec3, 3> makeThreePhotonDirections(
    RandomEngine& random, const OrthoEnergyPoint& point) {
    // Derive the direction of the most energetic photon from momentum
    // conservation.  Dividing by the largest, rather than a potentially
    // vanishing soft-photon energy, avoids catastrophic loss of precision at
    // the edges of the Dalitz triangle.
    const std::size_t derivedIndex = static_cast<std::size_t>(
        std::distance(point.fractions.begin(),
                      std::max_element(point.fractions.begin(),
                                       point.fractions.end())));
    std::array<std::size_t, 2> explicitIndices{};
    std::size_t explicitCount = 0;
    for (std::size_t index = 0; index < 3; ++index) {
        if (index != derivedIndex) explicitIndices[explicitCount++] = index;
    }
    const std::size_t firstIndex = explicitIndices[0];
    const std::size_t secondIndex = explicitIndices[1];
    const Vec3 firstDirection = isotropicDirection(random);

    // A random azimuth around photon 1 makes the decay-plane orientation
    // isotropic.  Photon 3 is obtained from momentum conservation rather than
    // from a second independently rounded angle.
    const Vec3 helper = std::abs(firstDirection.z) < 0.9
                      ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
    const Vec3 firstTransverse = normalized(cross(helper, firstDirection));
    const Vec3 secondTransverse = cross(firstDirection, firstTransverse);
    const double planeAzimuth = 2.0 * pi * uniformUnit(random);
    const Vec3 radial = firstTransverse * std::cos(planeAzimuth)
                      + secondTransverse * std::sin(planeAzimuth);
    const double pairAngleCosine = pairCosine(point, firstIndex, secondIndex);
    const double pairAngleSine =
        std::sqrt(std::max(0.0, 1.0 - pairAngleCosine * pairAngleCosine));
    const Vec3 secondDirection =
        firstDirection * pairAngleCosine + radial * pairAngleSine;

    std::array<Vec3, 3> directions{};
    directions[firstIndex] = firstDirection;
    directions[secondIndex] = secondDirection;
    directions[derivedIndex] =
        -(firstDirection * point.fractions[firstIndex]
        + secondDirection * point.fractions[secondIndex])
        * (1.0 / point.fractions[derivedIndex]);
    return directions;
}

} // namespace detail

template <class RandomEngine>
DecayEvent generateParaDecay(RandomEngine& random) {
    DecayEvent event;
    event.state = PositroniumState::Para;
    event.decayTimeSeconds =
        detail::exponentialTime(random, paraLifetimeSeconds);
    event.photonCount = 2;

    const Vec3 direction = detail::isotropicDirection(random);
    const double photonEnergy = 0.5 * positroniumRestEnergyJoules;
    event.photons[0] = {photonEnergy, direction};
    event.photons[1] = {photonEnergy, -direction};
    event.paraPhotonCosPolar = direction.z;
    return event;
}

template <class RandomEngine>
DecayEvent generateOrthoDecay(RandomEngine& random) {
    DecayEvent event;
    event.state = PositroniumState::Ortho;
    event.decayTimeSeconds =
        detail::exponentialTime(random, orthoLifetimeSeconds);
    event.photonCount = 3;

    const detail::OrthoEnergyPoint point =
        detail::sampleOrePowellEnergies(random);
    const std::array<Vec3, 3> directions =
        detail::makeThreePhotonDirections(random, point);
    for (std::size_t index = 0; index < event.photonCount; ++index) {
        event.photons[index] = {
            0.5 * positroniumRestEnergyJoules * point.fractions[index],
            directions[index]
        };
    }

    std::array<std::size_t, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](std::size_t first,
                                               std::size_t second) {
        return event.photons[first].energyJoules
             > event.photons[second].energyJoules;
    });
    const Photon& maximum = event.photons[order[0]];
    const Photon& middle = event.photons[order[1]];
    event.orthoMaximumEnergyJoules = maximum.energyJoules;
    event.orthoMiddleEnergyJoules = middle.energyJoules;
    event.orthoLeadingPairAngleRadians = std::acos(
        detail::clampCosine(dot(maximum.direction, middle.direction)));
    return event;
}

template <class RandomEngine>
DecayEvent generateDecay(PositroniumState state, RandomEngine& random) {
    return state == PositroniumState::Para
         ? generateParaDecay(random) : generateOrthoDecay(random);
}

} // namespace bound_decay

#endif // POSITRONIUM_BOUND_DECAY_HPP
