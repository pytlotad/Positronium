#pragma once

// Classical zero-point field of stochastic electrodynamics: an optional,
// off-by-default physics model (--zpf), not part of the model every
// committed result was produced with -- see gZeroPointField's own use
// sites for how it enters the force sum only when active().
//
// Extracted verbatim from positronium.cpp (continuing the split of engine,
// experiments and ROOT presentation apart -- see the session notes).
// Textually included at the same point inside positronium.cpp's shared
// anonymous namespace it always occupied, so it depends on that namespace
// already having in scope: Vec3, the c/hbar/pi/epsilon0 constants, dot(),
// cross().  Not yet a standalone, order-independent header.

// Classical zero-point field of stochastic electrodynamics: a random,
// homogeneous, isotropic radiation field with spectral energy density
//
//     rho(omega) = hbar omega^3 / (2 pi^2 c^3).
//
// It is NOT a drag.  Radiation reaction is the dissipative half and is already
// in the model; this is the fluctuating half, and it pumps energy INTO the
// orbit.  Together they are a Langevin pair whose stationary state is the SED
// ground state, which is where the L = hbar initial condition this project
// already uses comes from.  The omega^3 spectrum is the unique one invariant
// under boosts, which is exactly why it cannot act as a drag: a drag would
// single out a rest frame.
//
// Represented as equal-energy plane waves.  Sampling omega from the inverse
// CDF of omega^3 over the band, omega = (w1^4 + x(w2^4 - w1^4))^(1/4), makes
// every mode carry the same energy, so one amplitude serves all of them:
// u = eps0 <E^2> and <cos^2> = 1/2 give E = sqrt(2u/(N eps0)).
//
// The band cannot be the whole spectrum.  It is centred on the PAIR's own
// orbital frequency, because that is where the secular exchange happens, and
// its upper edge is what the integrator has to resolve -- the step is set by
// the trajectory error probe, so a wider band buys accuracy in the field at
// the price of a much shorter step.
struct ZeroPointField {
    // Directions and polarizations are fixed; the FREQUENCIES are not.  Each
    // mode carries a dimensionless factor f_n against the pair's osculating
    // orbital frequency, so the whole band rides up with the orbit as it
    // shrinks and stays in resonance instead of being left behind.
    std::vector<Vec3> direction, polarization, magneticDirection;
    std::vector<double> phase, frequencyFactor;
    // amplitude = coefficient * omega_orb^2, because the band energy density
    // integral of hbar omega^3 over [f_lo, f_hi]*omega_orb scales as
    // omega_orb^4 and the amplitude is its square root.  The field therefore
    // strengthens as the orbit tightens, which is what the omega^3 spectrum
    // says it should do.
    double amplitudeCoefficient=0.0;
    bool active() const {
        return amplitudeCoefficient>0.0&&!direction.empty();
    }
    void sample(const Vec3& position,double orbitalFrequency,
                double accumulatedPhase,Vec3& electric,Vec3& magnetic) const {
        electric=Vec3{};
        magnetic=Vec3{};
        if(!(orbitalFrequency>0.0)) return;
        const double amplitude=amplitudeCoefficient
            *orbitalFrequency*orbitalFrequency;
        for(std::size_t n=0;n<direction.size();++n) {
            const double waveNumber=frequencyFactor[n]*orbitalFrequency/c;
            // The temporal phase is f_n times the INTEGRATED orbital phase,
            // so its time derivative is exactly the mode's instantaneous
            // frequency and nothing jumps when the band moves.  The spatial
            // term is kept for completeness but is tiny: k*r is about 3e-3 rad
            // across the pair, which is the long-wavelength limit.
            const double wave=std::cos(
                dot(direction[n],position)*waveNumber
                -frequencyFactor[n]*accumulatedPhase+phase[n]);
            electric+=polarization[n]*(amplitude*wave);
            magnetic+=magneticDirection[n]*(amplitude*wave/c);
        }
    }
    // Translational force on a magnetic dipole from the SPATIAL variation
    // of this field, F=grad(moment.B(r)).  Unlike the mutual dipole-dipole
    // force (whose source moves, so its gradient needs a numerical
    // difference), every mode here is a genuine plane wave with a known
    // closed-form spatial gradient: for B_n=magneticDirection[n]*(amplitude
    // /c)*cos(k_n.r+psi_n), grad(moment.B_n) = -moment.magneticDirection[n]
    // *(amplitude/c)*waveNumber_n*sin(k_n.r+psi_n)*direction[n], since
    // k_n=waveNumber_n*direction[n] and grad(dot(direction[n],r)) is just
    // direction[n]. Called only where the dipole actually sits (each
    // particle's own position), not summed over both the way sample() is
    // used for the mutual-field bookkeeping in localRelativisticFields.
    Vec3 gradientForce(const Vec3& position,double orbitalFrequency,
                       double accumulatedPhase,const Vec3& moment) const {
        Vec3 force;
        if(!(orbitalFrequency>0.0)) return force;
        const double amplitude=amplitudeCoefficient
            *orbitalFrequency*orbitalFrequency;
        for(std::size_t n=0;n<direction.size();++n) {
            const double waveNumber=frequencyFactor[n]*orbitalFrequency/c;
            const double phaseValue=dot(direction[n],position)*waveNumber
                -frequencyFactor[n]*accumulatedPhase+phase[n];
            const double momentDotMagnetic=dot(moment,magneticDirection[n]);
            force=force-direction[n]*(momentDotMagnetic*(amplitude/c)
                *waveNumber*std::sin(phaseValue));
        }
        return force;
    }
};
ZeroPointField gZeroPointField;

// scale multiplies the field AMPLITUDE, so 1 is the full zero-point level and
// the absorbed power scales as scale^2.  Anything other than 1 is a numerical
// experiment, not the physical field.
ZeroPointField makeZeroPointField(double lowFactor,double highFactor,
                                  int modeCount,double scale,
                                  std::uint64_t seed) {
    ZeroPointField field;
    if(!(scale>0.0)||modeCount<=0||!(highFactor>lowFactor)) return field;
    const double lowFourth=std::pow(lowFactor,4.0);
    const double highFourth=std::pow(highFactor,4.0);
    // u(omega_orb) = hbar (f_hi^4 - f_lo^4) omega_orb^4 / (8 pi^2 c^3), and
    // u = eps0 <E^2> with <cos^2> = 1/2 over N equal-energy modes.
    field.amplitudeCoefficient=scale*std::sqrt(
        2.0*hbar*(highFourth-lowFourth)
        /(8.0*pi*pi*c*c*c*static_cast<double>(modeCount)*epsilon0));
    std::mt19937_64 random(seed^0x5851f42d4c957f2dULL);
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    for(int mode=0;mode<modeCount;++mode) {
        const double factor=std::pow(lowFourth
            +uniform(random)*(highFourth-lowFourth),0.25);
        const double cosine=2.0*uniform(random)-1.0;
        const double azimuth=2.0*pi*uniform(random);
        const double transverse=std::sqrt(std::max(0.0,1.0-cosine*cosine));
        const Vec3 propagation{transverse*std::cos(azimuth),
                               transverse*std::sin(azimuth),cosine};
        const Vec3 seedAxis=std::abs(propagation.z)<0.9?Vec3{0,0,1}:Vec3{1,0,0};
        Vec3 first=cross(propagation,seedAxis);
        first=first/first.norm();
        const Vec3 second=cross(propagation,first);
        const double polarizationAngle=2.0*pi*uniform(random);
        const Vec3 polarization=first*std::cos(polarizationAngle)
                               +second*std::sin(polarizationAngle);
        field.direction.push_back(propagation);
        field.polarization.push_back(polarization);
        field.magneticDirection.push_back(cross(propagation,polarization));
        field.phase.push_back(2.0*pi*uniform(random));
        field.frequencyFactor.push_back(factor);
    }
    return field;
}
