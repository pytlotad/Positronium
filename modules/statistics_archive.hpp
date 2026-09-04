#pragma once

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

// Scientific reference catalogue shared by statistical analyses and their
// human-readable export. It deliberately has no dependency on ROOT.
namespace statistics_archive {

inline constexpr int schemaVersion = 1;
inline constexpr unsigned int allPhenomenaMask = 0x0fU;

inline constexpr double paraDecayRatePerMicrosecond = 7990.9;
inline constexpr double paraDecayRateTotalError = 1.7;
inline constexpr double paraExperimentalLifetimePs =
    1.0e6/paraDecayRatePerMicrosecond;
inline constexpr double paraExperimentalLifetimeErrorPs =
    1.0e6*paraDecayRateTotalError
    /(paraDecayRatePerMicrosecond*paraDecayRatePerMicrosecond);

inline constexpr double orthoDecayRatePerMicrosecond = 7.0404;
inline constexpr double orthoDecayRateStatError = 0.0010;
inline constexpr double orthoDecayRateSysError = 0.0008;
inline constexpr double orthoDecayRateTotalError =
    0.0012806248474865698;
inline constexpr double orthoExperimentalLifetimeNs =
    1000.0/orthoDecayRatePerMicrosecond;
inline constexpr double orthoExperimentalLifetimeStatErrorNs =
    1000.0*orthoDecayRateStatError
    /(orthoDecayRatePerMicrosecond*orthoDecayRatePerMicrosecond);
inline constexpr double orthoExperimentalLifetimeSysErrorNs =
    1000.0*orthoDecayRateSysError
    /(orthoDecayRatePerMicrosecond*orthoDecayRatePerMicrosecond);
inline constexpr double orthoExperimentalLifetimeErrorNs =
    1000.0*orthoDecayRateTotalError
    /(orthoDecayRatePerMicrosecond*orthoDecayRatePerMicrosecond);

inline constexpr double pdgPositroniumMassEnergyMeV = 1.021991097;
inline constexpr double pdgPositroniumMassEnergyErrorMeV = 0.000000003;
inline constexpr double pdgPositroniumPhotonEnergyKeV = 510.9955485;
inline constexpr double pdgPositroniumPhotonEnergyErrorKeV = 0.0000015;

struct ScientificSource {
    std::string id;
    std::string authors;
    std::string title;
    std::string publication;
    int year = 0;
    std::string doi;
    std::string url;
    std::string accessedIsoDate;
};

struct ScientificValue {
    std::string id;
    std::string sourceId;
    unsigned int phenomenonMask = 0;
    std::string observableId;
    std::string role;
    std::string quantity;
    std::string symbol;
    double value = 0.0;
    std::string unit;
    bool hasValue = false;
    double statUncertainty = 0.0;
    double systematicUncertainty = 0.0;
    double totalUncertainty = 0.0;
    bool hasStatUncertainty = false;
    bool hasSystematicUncertainty = false;
    bool hasTotalUncertainty = false;
    bool exact = false;
    bool derived = false;
    std::string derivation;
    std::string note;
};

struct ScientificModel {
    std::string id;
    std::string sourceId;
    unsigned int phenomenonMask = 0;
    std::string observableId;
    std::string role;
    std::string formulaAscii;
    std::string domain;
    std::string assumptions;
    std::string implementationNote;
};

inline const std::vector<ScientificSource>& scientificSources() {
    static const std::vector<ScientificSource> sources{
        {"pdg_positronium_2025", "Particle Data Group",
         "Atomic and nuclear properties of positronium (Ps)",
         "Review of Particle Physics 2025 tables", 2025, "",
         "https://pdg.lbl.gov/2025/AtomicNuclearProperties/positronium.html",
         "2026-08-11"},
        {"al_ramadhan_gidley_1994", "A. H. Al-Ramadhan; D. W. Gidley",
         "New precision measurement of the decay rate of singlet positronium",
         "Physical Review Letters 72, 1632-1635", 1994,
         "10.1103/PhysRevLett.72.1632",
         "https://doi.org/10.1103/PhysRevLett.72.1632", "2026-08-11"},
        {"vallery_zitzewitz_gidley_2003",
         "R. S. Vallery; P. W. Zitzewitz; D. W. Gidley",
         "Resolution of the Orthopositronium-Lifetime Puzzle",
         "Physical Review Letters 90, 203402", 2003,
         "10.1103/PhysRevLett.90.203402",
         "https://doi.org/10.1103/PhysRevLett.90.203402", "2026-08-11"},
        {"ore_powell_1949", "A. Ore; J. L. Powell",
         "Three-Photon Annihilation of an Electron-Positron Pair",
         "Physical Review 75, 1696-1699", 1949,
         "10.1103/PhysRev.75.1696",
         "https://doi.org/10.1103/PhysRev.75.1696", "2026-08-11"},
        {"geant4_ore_powell", "Geant4 Collaboration",
         "Decays of orthopositronium at rest to 3 photons",
         "Geant4 Physics Reference Manual 11.4", 2025, "",
         "https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/decay/OrthoPositronium.html",
         "2026-08-11"},
        {"chang_tang_li_1985", "T. C. Chang; C. M. Tang; W. L. Li",
         "Gamma-ray energy spectrum from orthopositronium three-gamma decay",
         "Physics Letters B 157, 357-360", 1985,
         "10.1016/0370-2693(85)90380-6",
         "https://doi.org/10.1016/0370-2693(85)90380-6", "2026-08-11"},
        {"rutherford_1911", "E. Rutherford",
         "The Scattering of Alpha and Beta Particles by Matter and the Structure of the Atom",
         "Philosophical Magazine 21, 669-688", 1911,
         "10.1080/14786440508637080",
         "https://doi.org/10.1080/14786440508637080", "2026-08-11"},
        {"internal_ideal_model", "Positronium simulation",
         "Ideal statistical model assumptions", "In-program model definition",
         2026, "", "", "2026-08-11"},
        // Stochastic electrodynamics.  These back no plotted VALUE -- nothing
        // in scientificValues() cites them -- and are catalogued because the
        // README's --zpf section rests on the SED literature while the
        // bibliography carried nothing from it, so a reader had no way to
        // check the framing or to find where the same question is argued.
        //
        // Metadata verified against Crossref, not written from memory, and
        // the check caught an error worth recording: the Cole and Zou title
        // ends "from classical electrodynamics", not "from classical
        // physics".  What these papers CONCLUDE is deliberately not
        // paraphrased anywhere -- Crossref returns no abstracts for them and
        // they were not read here, so only their existence, authorship and
        // titles are asserted.  The two 2003 entries are published Comments,
        // which their own titles establish; that the Cole-Zou result was
        // contested in print is therefore citable, while who was right is
        // not claimed.
        {"boyer_1975", "Timothy H. Boyer",
         "Random electrodynamics: The theory of classical electrodynamics "
         "with classical electromagnetic zero-point radiation",
         "Physical Review D 11, 790-808", 1975,
         "10.1103/PhysRevD.11.790",
         "https://doi.org/10.1103/PhysRevD.11.790", "2026-08-28"},
        {"cole_zou_2003", "Daniel C. Cole; Yi Zou",
         "Quantum mechanical ground state of hydrogen obtained from "
         "classical electrodynamics",
         "Physics Letters A 317, 14-20", 2003,
         "10.1016/j.physleta.2003.08.022",
         "https://doi.org/10.1016/j.physleta.2003.08.022", "2026-08-28"},
        {"boyer_2003_comment", "Timothy H. Boyer",
         "Comments on Cole and Zou's Calculation of the Hydrogen Ground "
         "State in Classical Physics",
         "Foundations of Physics Letters 16, 613-617", 2003,
         "10.1023/B:FOPL.0000012787.05764.4d",
         "https://doi.org/10.1023/B:FOPL.0000012787.05764.4d", "2026-08-28"},
        {"milonni_2003_comment", "Peter W. Milonni",
         "Comment on Cole and Zou's Classical Computations of the Hydrogen "
         "Ground State",
         "Foundations of Physics Letters 16, 619-621", 2003,
         "10.1023/B:FOPL.0000012788.94843.5b",
         "https://doi.org/10.1023/B:FOPL.0000012788.94843.5b", "2026-08-28"},
        {"nieuwenhuizen_liska_2015", "Theo M. Nieuwenhuizen; Matthew T. P. Liska",
         "Simulation of the hydrogen ground state in stochastic "
         "electrodynamics",
         "Physica Scripta T165, 014006", 2015,
         "10.1088/0031-8949/2015/T165/014006",
         "https://doi.org/10.1088/0031-8949/2015/T165/014006", "2026-08-28"},
        {"nieuwenhuizen_liska_2015b",
         "Theodorus M. Nieuwenhuizen; Matthew T. P. Liska",
         "Simulation of the Hydrogen Ground State in Stochastic "
         "Electrodynamics-2: Inclusion of Relativistic Corrections",
         "Foundations of Physics 45, 1190-1202", 2015,
         "10.1007/s10701-015-9919-0",
         "https://doi.org/10.1007/s10701-015-9919-0", "2026-08-28"},

        // The engine's own physics and numerics.  Until now the catalogue
        // covered only the reference curves the plots are compared against
        // (annihilation spectra, Rutherford, the SED literature); the models
        // the program actually integrates carried their attribution in source
        // comments alone, where nothing checks them.  These are the sources
        // for what runMechanicalTrajectory and its callees implement.
        {"jackson_classical_electrodynamics",
         "John David Jackson", "Classical Electrodynamics",
         "3rd edition, Wiley; chapters 6, 11, 12, 14, 16", 1998, "",
         "https://www.wiley.com/en-us/Classical+Electrodynamics%2C+3rd+Edition-p-9780471309321",
         "2026-09-04"},
        {"landau_lifshitz_classical_fields",
         "Lev D. Landau; Evgeny M. Lifshitz", "The Classical Theory of Fields",
         "Course of Theoretical Physics volume 2, 4th revised English edition, "
         "Butterworth-Heinemann; section 75 and section 76 problem 1", 1975,
         "", "https://archive.org/details/ClassicalTheoryOfFields",
         "2026-09-04"},
        {"dirac_1938",
         "Paul A. M. Dirac", "Classical theory of radiating electrons",
         "Proceedings of the Royal Society A 167, 148-169", 1938,
         "10.1098/rspa.1938.0124",
         "https://doi.org/10.1098/rspa.1938.0124", "2026-09-04"},
        {"bargmann_michel_telegdi_1959",
         "Valentine Bargmann; Louis Michel; Valentine L. Telegdi",
         "Precession of the Polarization of Particles Moving in a "
         "Homogeneous Electromagnetic Field",
         "Physical Review Letters 2, 435-436", 1959,
         "10.1103/PhysRevLett.2.435",
         "https://doi.org/10.1103/PhysRevLett.2.435", "2026-09-04"},
        {"peters_1964",
         "Peter C. Peters",
         "Gravitational Radiation and the Motion of Two Point Masses",
         "Physical Review 136, B1224-B1232", 1964,
         "10.1103/PhysRev.136.B1224",
         "https://doi.org/10.1103/PhysRev.136.B1224", "2026-09-04"},
        {"plummer_1911",
         "Henry C. Plummer",
         "On the Problem of Distribution in Globular Star Clusters: "
         "(Plate 8.)",
         "Monthly Notices of the Royal Astronomical Society 71, 460-470", 1911,
         "10.1093/mnras/71.5.460",
         "https://doi.org/10.1093/mnras/71.5.460", "2026-09-04"},
        {"yee_1966",
         "Kane S. Yee",
         "Numerical solution of initial boundary value problems involving "
         "Maxwell's equations in isotropic media",
         "IEEE Transactions on Antennas and Propagation 14, 302-307", 1966,
         "10.1109/TAP.1966.1138693",
         "https://doi.org/10.1109/TAP.1966.1138693", "2026-09-04"},
        {"roden_gedney_2000",
         "J. Alan Roden; Stephen D. Gedney",
         "Convolution PML (CPML): An efficient FDTD implementation of the "
         "CFS-PML for arbitrary media",
         "Microwave and Optical Technology Letters 27, 334-339", 2000,
         "10.1002/1098-2760(20001205)27:5<334::AID-MOP14>3.0.CO;2-A",
         "https://doi.org/10.1002/1098-2760(20001205)27:5<334::AID-MOP14>"
         "3.0.CO;2-A",
         "2026-09-04"},
        {"boris_1970",
         "Jay P. Boris",
         "Relativistic plasma simulation - optimization of a hybrid code",
         "Proceedings of the Fourth Conference on Numerical Simulation of "
         "Plasmas, Naval Research Laboratory, 3-67", 1970, "",
         "https://apps.dtic.mil/sti/citations/ADA023511", "2026-09-04"},
        {"qin_2013_boris",
         "Hong Qin; Shuangxi Zhang; Jianyuan Xiao; Jian Liu; Yajuan Sun; "
         "William M. Tang",
         "Why is Boris algorithm so good?",
         "Physics of Plasmas 20, 084503", 2013,
         "10.1063/1.4818428",
         "https://doi.org/10.1063/1.4818428", "2026-09-04"},
        {"kaplan_meier_1958",
         "Edward L. Kaplan; Paul Meier",
         "Nonparametric Estimation from Incomplete Observations",
         "Journal of the American Statistical Association 53, 457-481", 1958,
         "10.1080/01621459.1958.10501452",
         "https://doi.org/10.1080/01621459.1958.10501452", "2026-09-04"},
        {"matsumoto_nishimura_1998",
         "Makoto Matsumoto; Takuji Nishimura",
         "Mersenne twister",
         "ACM Transactions on Modeling and Computer Simulation 8, 3-30", 1998,
         "10.1145/272991.272995",
         "https://doi.org/10.1145/272991.272995", "2026-09-04"},
        {"steele_lea_flood_2014",
         "Guy L. Steele Jr.; Doug Lea; Christine H. Flood",
         "Fast splittable pseudorandom number generators",
         "ACM SIGPLAN Notices 49, 453-472 (OOPSLA 2014)", 2014,
         "10.1145/2714064.2660195",
         "https://doi.org/10.1145/2714064.2660195", "2026-09-04"},
        {"tiesinga_codata_2018",
         "Eite Tiesinga; Peter J. Mohr; David B. Newell; Barry N. Taylor",
         "CODATA recommended values of the fundamental physical constants: "
         "2018",
         "Reviews of Modern Physics 93, 025010", 2021,
         "10.1103/RevModPhys.93.025010",
         "https://doi.org/10.1103/RevModPhys.93.025010", "2026-09-04"}
    };
    return sources;
}

inline const std::vector<ScientificValue>& scientificValues() {
    static const std::vector<ScientificValue> values{
        {"para_decay_rate_measurement", "al_ramadhan_gidley_1994", 0x01U,
         "annihilation_time", "experimental_measurement",
         "singlet positronium decay rate", "lambda_p",
         paraDecayRatePerMicrosecond, "us^-1", true,
         0.0, 0.0, paraDecayRateTotalError,
         false, false, true, false, false, "",
         "Published total uncertainty."},
        {"para_lifetime_from_rate", "al_ramadhan_gidley_1994", 0x01U,
         "annihilation_time", "experimental_derived",
         "singlet positronium lifetime", "tau_p",
         paraExperimentalLifetimePs, "ps", true,
         0.0, 0.0, paraExperimentalLifetimeErrorPs,
         false, false, true, false, true,
         "tau_p[ps]=1e6/lambda_p[us^-1]; linear uncertainty propagation",
         "Value displayed in the Experimental section."},
        {"ortho_decay_rate_measurement", "vallery_zitzewitz_gidley_2003", 0x02U,
         "annihilation_time", "experimental_measurement",
         "orthopositronium decay rate", "lambda_o",
         orthoDecayRatePerMicrosecond, "us^-1", true,
         orthoDecayRateStatError, orthoDecayRateSysError,
         orthoDecayRateTotalError, true, true, true, false, false, "",
         "Statistical and systematic uncertainties retained separately."},
        {"ortho_lifetime_from_rate", "vallery_zitzewitz_gidley_2003", 0x02U,
         "annihilation_time", "experimental_derived",
         "orthopositronium lifetime", "tau_o",
         orthoExperimentalLifetimeNs, "ns", true,
         orthoExperimentalLifetimeStatErrorNs,
         orthoExperimentalLifetimeSysErrorNs,
         orthoExperimentalLifetimeErrorNs, true, true, true, false, true,
         "tau_o[ns]=1000/lambda_o[us^-1]; linear uncertainty propagation",
         "Value displayed in the Experimental section."},
        {"pdg_positronium_mass_energy", "pdg_positronium_2025", 0x03U,
         "photon_energy", "evaluated_accepted_value",
         "ground-state positronium mass energy", "m_Ps_c2",
         pdgPositroniumMassEnergyMeV, "MeV", true,
         0.0, 0.0, pdgPositroniumMassEnergyErrorMeV,
         false, false, true, false, false, "",
         "PDG evaluated value."},
        {"pdg_two_photon_energy", "pdg_positronium_2025", 0x03U,
         "photon_energy", "evaluated_derived",
         "photon energy for at-rest two-body split", "E_gamma",
         pdgPositroniumPhotonEnergyKeV, "keV", true,
         0.0, 0.0, pdgPositroniumPhotonEnergyErrorKeV,
         false, false, true, false, true,
         "E_gamma[keV]=500*m_Ps_c2[MeV]",
         "Line energy for p-Ps and kinematic endpoint used for o-Ps comparison."},
        {"ideal_para_anisotropy", "internal_ideal_model", 0x01U,
         "photon_polar_angle", "theory_benchmark",
         "unpolarized p-Ps Legendre anisotropy", "a2",
         0.0, "1", true, 0.0, 0.0, 0.0,
         false, false, false, true, false, "",
         "Exact only within the stated ideal unpolarized spin-zero model."},
        {"rutherford_shape_normalization", "rutherford_1911", 0x0cU,
         "differential_cross_section", "analytic_benchmark",
         "Rutherford-shape normalization", "C_R",
         1.0, "1", true, 0.0, 0.0, 0.0,
         false, false, false, true, false, "",
         "Exact only within the pure instantaneous Coulomb reference model; not experimental data."}
    };
    return values;
}

inline const std::vector<ScientificModel>& scientificModels() {
    static const std::vector<ScientificModel> models{
        {"crem_secular_collapse", "internal_ideal_model", 0x03U,
         "crem_collapse_time", "trajectory_extrapolation",
         "dE/dt=P_CREM; E=-k*e^2/(2a); <P>(a) proportional to a^-4",
         "bound CREM trajectory; terminal radius 0.01*a0",
         "Short full-CREM calibration followed by orbit-averaged secular evolution.",
         "The descriptive exponential fit is not an annihilation law; published lifetimes are comparison curves only."},
        {"ore_powell_three_gamma", "ore_powell_1949", 0x02U,
         "photon_energy,three_photon_dalitz,leading_photon_angle",
         "theory_model",
         "P proportional sum_ij(1-cos(theta_ij))^2 dk1 dk2 dOmega1; "
         "F(x)=x(1-x)/(2-x)^2-2(1-x)^2*ln(1-x)/(2-x)^3+(2-x)/x+2(1-x)*ln(1-x)/x^2",
         "x=E/Emax; 0<=x<=1 (x=0 understood by the limiting value)",
         "Unpolarized o-Ps at rest; leading order; photon polarizations summed.",
         "Sampling follows the Geant4 G4OrePowellAtRestModel prescription."},
        {"geant4_ore_powell_implementation", "geant4_ore_powell", 0x02U,
         "photon_energy,three_photon_dalitz,leading_photon_angle",
         "implementation_reference", "acceptance proportional to sum_ij(1-cos(theta_ij))^2",
         "physical three-photon phase space", "Random polarization.",
         "Reference for the implemented rejection sampler."},
        {"ortho_continuum_validation", "chang_tang_li_1985", 0x02U,
         "photon_energy", "experimental_validation", "",
         "apparatus-specific measured spectrum", "Detector response required.",
         "Published continuum was reported consistent with QED; no numeric overlay is embedded."},
        {"rutherford_cross_section", "rutherford_1911", 0x0cU,
         "differential_cross_section,cumulative_cross_section",
         "analytic_benchmark",
         "dSigma/dOmega=l_C^2/(4*sin(theta/2)^4); sigma(theta>=theta0)=pi*l_C^2*cot(theta0/2)^2",
         "0<theta<=pi", "Classical Coulomb scattering.",
         "The trajectory result fits a single normalization C_R to this fixed shape."},

        // The engine's own models.  Every entry below names the routine that
        // implements it, so the catalogue can be read as a map from the
        // literature into the source rather than as a bibliography.  The
        // phenomenon mask is 0x0f wherever the mechanical engine runs, which
        // is all four.
        {"lienard_wiechert_retarded_field", "jackson_classical_electrodynamics",
         0x0fU, "crem_collapse_time,differential_cross_section",
         "engine_model",
         "E=q/(4 pi eps0)[(n-beta)(1-beta^2)/(kappa^3 R^2) "
         "+ n x ((n-beta) x a)/(c^2 kappa^3 R)]; B=n x E/c; kappa=1-n.beta",
         "point charges outside the declared regularization core; "
         "retarded time from the light-cone equation",
         "Finite retained history; the light-cone solve is a damped Newton "
         "iteration, not exact.",
         "lienardWiechertField in modules/retarded_charge_kinematics.hpp; the "
         "Newton loop and its coincidence guard are in the same function."},
        {"landau_lifshitz_radiation_reaction",
         "landau_lifshitz_classical_fields", 0x0fU, "crem_collapse_time",
         "engine_model",
         "f = (2 q^3 gamma/(3 c^3))[dE/dt + ...] with the ADD acceleration "
         "replaced by its zeroth-order equation-of-motion value",
         "order-reduced; valid while the reaction is a small correction to "
         "the Lorentz force",
         "Order reduction removes the runaway solutions of Abraham-Lorentz-"
         "Dirac at the cost of being perturbative in the reaction.",
         "reducedOrderSelfForce and individualLandauLifshitzSelfForces in "
         "modules/electrodynamics.hpp; selected by --radiation-reaction."},
        {"abraham_lorentz_dirac_origin", "dirac_1938", 0x0fU,
         "crem_collapse_time", "theory_context",
         "m dv/dt = F_ext + (2 q^2/3 c^3)(d2v/dt2 + ...)",
         "the equation the order reduction above starts from",
         "Runaway and pre-accelerating solutions are why the reduced form is "
         "used instead.",
         "Not integrated directly; recorded because the reduction is only "
         "meaningful against it."},
        {"thomas_bmt_spin_transport", "bargmann_michel_telegdi_1959", 0x0fU,
         "crem_collapse_time", "engine_model",
         "d s/dt = (q/m) s x [(g/2-1+1/gamma) B "
         "- (g/2-1) gamma/(gamma+1) (beta.B) beta "
         "- (g/2 - gamma/(gamma+1)) beta x E/c]",
         "classical spin of a particle with an anomalous moment in an "
         "external field",
         "g is a per-particle argument, not a constant: 2.0023 for a lepton, "
         "5.5857 for a proton.",
         "thomasBmtEffectiveField in modules/electrodynamics.hpp and the "
         "orbit-averaged form in modules/secular_spin_orbit.hpp."},
        {"darwin_two_body_lagrangian", "jackson_classical_electrodynamics",
         0x0fU, "crem_collapse_time", "engine_model",
         "L_D = C[v1.v2 + (v1.n)(v2.n)], C = q1 q2/(2 c^2 r)",
         "two charges to order (v/c)^2",
         "Accelerations inside d/dt(dL/dv) are reduced to their leading "
         "Coulomb values, consistently through order v^2/c^2.",
         "darwinForceOnFirst and darwinForces in modules/electrodynamics.hpp."},
        {"moving_point_dipole_transformation",
         "jackson_classical_electrodynamics", 0x0fU, "crem_collapse_time",
         "engine_model",
         "P_source = p_lab/gamma; M_source = mu_lab/gamma, differentiated "
         "AFTER the conversion",
         "a moving, accelerating point dipole",
         "Dividing p_dot or mu_dot by gamma instead would omit the "
         "derivatives of gamma.",
         "historicalIntegratedDipoleKinematics in modules/electrodynamics.hpp; "
         "the closed form is quoted in the source comments as Sautbekov's."},
        {"orbit_averaged_dipole_inspiral", "peters_1964", 0x03U,
         "crem_collapse_time", "engine_model",
         "<dE/dt> = -gamma k^2 (1+e^2/2)/(a^4 (1-e^2)^(5/2)); "
         "<dL/dt> = -gamma k L/(a^3 (1-e^2)^(3/2))",
         "orbit-averaged over an unperturbed Kepler ellipse",
         "Same averaging method as Peters, applied to the DIPOLE's P ~ r^-4 "
         "rather than the quadrupole's r^-6; the eccentricity factors differ "
         "and the Peters factor must not be substituted.",
         "The k(e) used by estimateCremCollapse in modules/crem_collapse.hpp."},
        {"plummer_short_range_softening", "plummer_1911", 0x0fU,
         "crem_collapse_time", "numerical_method",
         "r_eff = sqrt(r^2 + floor^2), direction preserved",
         "close encounters inside the declared regularization radius",
         "Every derivative stays continuous, unlike a hard clamp to the "
         "floor, which puts a measurable kink in the force at r = floor.",
         "clampedSeparationVector in modules/pair_geometry.hpp; the floor is "
         "the Compton barrier for e+e-, not an arbitrary softening length."},
        {"yee_staggered_fdtd", "yee_1966", 0x0fU, "field_validation",
         "numerical_method",
         "curl updates on a staggered E/B grid, second order in space and "
         "time",
         "the optional Maxwell backend of the validation build only",
         "Not a production path: positronium.cpp includes the backend only "
         "under POSITRONIUM_ENABLE_FIELD_VALIDATION.",
         "modules/maxwell_validation_backend.hpp."},
        {"convolutional_pml_boundary", "roden_gedney_2000", 0x0fU,
         "field_validation", "numerical_method",
         "CFS-PML stretched coordinates with a recursive convolution update",
         "outer boundary of the Yee grid above",
         "Absorbs outgoing radiation so a finite grid does not reflect it "
         "back onto the pair.",
         "The CPML sections of modules/maxwell_validation_backend.hpp."},
        {"boris_particle_push", "boris_1970", 0x0fU, "field_validation",
         "numerical_method",
         "electric half-kick, magnetic rotation at the intermediate gamma, "
         "electric half-kick, written in momentum variables",
         "grid-coupled particle advance in the validation build",
         "Time reversible and cannot produce |v| >= c.",
         "relativisticBorisPush in modules/maxwell_validation.hpp."},
        {"boris_phase_space_volume", "qin_2013_boris", 0x0fU,
         "field_validation", "theory_context",
         "the Boris map is volume preserving in phase space",
         "why the pusher above is used rather than a higher-order "
         "non-symplectic one",
         "Volume preservation, not order of accuracy, is what bounds the "
         "long-run energy error.",
         "Recorded as the justification for relativisticBorisPush."},
        {"kaplan_meier_survival", "kaplan_meier_1958", 0x03U,
         "crem_collapse_time", "statistical_method",
         "S(t) = prod_{t_i <= t} (1 - d_i/n_i)",
         "right-censored collapse times",
         "A trajectory stopped by the wall-clock budget is censored, not "
         "missing: it carries the information that collapse had not happened "
         "by the time it reached.",
         "kaplanMeier in modules/kaplan_meier.hpp; the budget preferentially "
         "stops the widest orbits, which are the slowest to collapse, so "
         "averaging only the completed runs would bias the result."},
        {"mersenne_twister_stream", "matsumoto_nishimura_1998", 0x0fU,
         "crem_collapse_time,differential_cross_section",
         "numerical_method",
         "MT19937-64",
         "every stochastic draw in the program",
         "Reproducibility of a run depends on the draw ORDER as well as the "
         "seed.",
         "std::mt19937_64 throughout; seeded per trajectory."},
        {"splitmix64_seed_mixing", "steele_lea_flood_2014", 0x0fU,
         "crem_collapse_time,differential_cross_section",
         "numerical_method",
         "z += 0x9e3779b97f4a7c15; z ^= z>>30; z *= 0xbf58476d1ce4e5b9; "
         "z ^= z>>27; z *= 0x94d049bb133111eb; z ^= z>>31",
         "deriving independent per-trajectory seeds from one run seed",
         "Avoids the correlated low bits a plain counter would hand the "
         "generator.",
         "splitMix64 in modules/sampling_utilities.hpp and its local twin "
         "stochasticPhotonHash64 in modules/crem_trajectory.hpp."},
        {"codata_2018_constants", "tiesinga_codata_2018", 0x0fU,
         "crem_collapse_time,annihilation_time,differential_cross_section",
         "constants_source",
         "",
         "every physical constant the program uses",
         "The 2019 SI redefinition makes c, e and h exact; the electron mass "
         "and the g-factor remain measured.",
         "modules/physical_constants.hpp, with the exact and measured values "
         "kept distinguishable there."},
        {"larmor_dipole_power", "jackson_classical_electrodynamics", 0x0fU,
         "crem_collapse_time", "engine_model",
         "P = q_eff^2 |d2d/dt2|^2/(6 pi eps0 c^3); orbit-averaged for a "
         "Kepler ellipse this gives <P> proportional to a^-4",
         "electric-dipole (E1) radiation of the pair's relative motion",
         "The a^-4 scaling, not the instantaneous formula, is what the "
         "secular estimator integrates.",
         "larmorOrbitAveragedPower in modules/crem_collapse.hpp; it also sets "
         "the photon hazard rate lambda = P/(hbar omega)."},
        {"schott_bound_field_energy", "jackson_classical_electrodynamics",
         0x0fU, "crem_collapse_time", "engine_model",
         "E_Schott = -q^2 (v.a)/(6 pi eps0 c^3), the total-derivative term of "
         "the radiation-reaction power",
         "an accelerating charge with a bound near field",
         "It is a reservoir, not a loss: it returns to the orbit and must not "
         "be counted as radiated energy.",
         "Tracked as schottEnergy through Frame and reported separately from "
         "radiatedEnergy in positronium.cpp."},
        {"poisson_photon_hazard", "internal_ideal_model", 0x03U,
         "crem_collapse_time", "engine_model",
         "P(no photon in dt) = exp(-lambda dt), lambda = P_Larmor/(hbar omega)",
         "the stochastic branch of the radiation-reaction switch",
         "The same E1 power that the deterministic drag would remove "
         "continuously is instead banked as hazard and paid out in discrete, "
         "momentum-conserving kicks, so a photon has an energy and a "
         "direction to report.",
         "stochasticElectricDipole in modules/electrodynamics.hpp; "
         "--emission deterministic bypasses the Poisson draw."},
        {"adaptive_step_doubling", "jackson_classical_electrodynamics", 0x0fU,
         "crem_collapse_time", "numerical_method",
         "one full step against two half steps; accept on the Richardson "
         "error estimate, otherwise subdivide",
         "the mechanical integrator of every trajectory",
         "Measured convergence is second order, 4.00x per halving held for "
         "about ten halvings, with the retarded sector in or out and the "
         "reaction on or off.",
         "ClassicalTrajectoryEngine in modules/crem_engine.hpp; the "
         "measurement is CREM_DEBUG_ORDER."},
        {"rayleigh_impact_parameter", "internal_ideal_model", 0x08U,
         "differential_cross_section", "statistical_method",
         "a circular isotropic 2D Gaussian beam has Rayleigh-distributed "
         "radial impact parameter, p(b) proportional to b exp(-b^2/2 sigma^2)",
         "beam and interaction experiments",
         "The area Jacobian contributes the factor b; sampling |N(0,sigma)| "
         "instead would be half-normal and wrong.",
         "sampleIsotropicGaussianImpact in modules/sampling_utilities.hpp."},
        {"stochastic_electrodynamics_field", "boyer_1975", 0x0fU,
         "crem_collapse_time", "optional_model",
         "rho(omega) = hbar omega^3/(2 pi^2 c^3)",
         "off by default; enabled by --zpf",
         "The fluctuating half of a Langevin pair whose dissipative half is "
         "the radiation reaction already in the model.  The omega^3 spectrum "
         "is the unique Lorentz-invariant one, which is why it cannot act as "
         "a drag.",
         "modules/zero_point_field.hpp; NOT part of the model any committed "
         "result was produced with."}
    };
    return models;
}

inline bool appliesToPhenomenon(unsigned int mask, int phenomenon) {
    return phenomenon >= 1 && phenomenon <= 4
        && (mask & (1U << static_cast<unsigned int>(phenomenon - 1))) != 0U;
}

inline std::string phenomenaForMask(unsigned int mask) {
    std::string result;
    for (int phenomenon = 1; phenomenon <= 4; ++phenomenon) {
        if (!appliesToPhenomenon(mask, phenomenon)) continue;
        if (!result.empty()) result += ',';
        result += std::to_string(phenomenon);
    }
    return result.empty() ? "none" : result;
}

inline std::string formatScientificNumber(double value) {
    char buffer[128];
    const auto conversion = std::to_chars(
        std::begin(buffer), std::end(buffer), value, std::chars_format::general);
    if (conversion.ec == std::errc{}) {
        return std::string(buffer, conversion.ptr);
    }
    std::ostringstream fallback;
    fallback << std::setprecision(std::numeric_limits<double>::max_digits10)
             << value;
    return fallback.str();
}

inline const ScientificValue& scientificValue(std::string_view id) {
    const auto& values = scientificValues();
    const auto found = std::find_if(values.begin(), values.end(),
        [&](const ScientificValue& value) { return value.id == id; });
    if (found == values.end()) throw std::out_of_range("unknown scientific value id");
    return *found;
}

struct OperationResult {
    std::filesystem::path path;
    std::string error;
    explicit operator bool() const noexcept { return error.empty(); }
};

inline std::string uniqueTemporarySuffix() {
    static std::atomic<std::uint64_t> sequence{0};
    const auto ticks = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    return std::to_string(ticks) + "_"
        + std::to_string(sequence.fetch_add(1, std::memory_order_relaxed));
}

inline OperationResult writeScientificReferencesText(
    const std::filesystem::path& destination = "ScientificalReferences.txt") {
    OperationResult result{destination, {}};
    const std::filesystem::path temporary = destination.parent_path()
        / ("." + destination.filename().string() + ".tmp_"
           + uniqueTemporarySuffix());
    struct Cleanup {
        std::filesystem::path path;
        ~Cleanup() {
            std::error_code ignored;
            std::filesystem::remove(path, ignored);
        }
    } cleanup{temporary};
    try {
        std::ofstream output(temporary, std::ios::trunc);
        if (!output) {
            result.error = "cannot open temporary reference file";
            return result;
        }
        output << "# Scientific references used by positronium statistical plots\n"
               << "schema_version=" << schemaVersion << "\n"
               << "generated_from=in-program scientific reference catalogue\n"
               << "numeric_format=shortest round-trip decimal\n\n";
        for (const ScientificSource& source : scientificSources()) {
            output << "[source:" << source.id << "]\n"
                   << "authors=" << source.authors << "\n"
                   << "title=" << source.title << "\n"
                   << "publication=" << source.publication << "\n"
                   << "year=" << source.year << "\n"
                   << "doi=" << source.doi << "\n"
                   << "url=" << source.url << "\n"
                   << "accessed=" << source.accessedIsoDate << "\n\n";
        }
        for (const ScientificValue& value : scientificValues()) {
            output << "[value:" << value.id << "]\n"
                   << "source=" << value.sourceId << "\n"
                   << "phenomenon_mask=" << value.phenomenonMask << "\n"
                   << "phenomena=" << phenomenaForMask(value.phenomenonMask) << "\n"
                   << "observable=" << value.observableId << "\n"
                   << "role=" << value.role << "\n"
                   << "quantity=" << value.quantity << "\n"
                   << "symbol=" << value.symbol << "\n"
                   << "value=";
            if (value.hasValue) output << formatScientificNumber(value.value);
            else output << "n/a";
            output << "\nunit=" << value.unit << "\nstat_uncertainty=";
            if (value.hasStatUncertainty) {
                output << formatScientificNumber(value.statUncertainty);
            }
            else output << "n/a";
            output << "\nsystematic_uncertainty=";
            if (value.hasSystematicUncertainty) {
                output << formatScientificNumber(value.systematicUncertainty);
            }
            else output << "n/a";
            output << "\ntotal_uncertainty=";
            if (value.hasTotalUncertainty) {
                output << formatScientificNumber(value.totalUncertainty);
            }
            else output << "n/a";
            output << "\n"
                   << "exact=" << (value.exact ? "true" : "false") << "\n"
                   << "derived=" << (value.derived ? "true" : "false") << "\n"
                   << "derivation=" << value.derivation << "\n"
                   << "note=" << value.note << "\n\n";
        }
        for (const ScientificModel& model : scientificModels()) {
            output << "[model:" << model.id << "]\n"
                   << "source=" << model.sourceId << "\n"
                   << "phenomenon_mask=" << model.phenomenonMask << "\n"
                   << "phenomena=" << phenomenaForMask(model.phenomenonMask) << "\n"
                   << "observable=" << model.observableId << "\n"
                   << "role=" << model.role << "\n"
                   << "formula=" << model.formulaAscii << "\n"
                   << "domain=" << model.domain << "\n"
                   << "assumptions=" << model.assumptions << "\n"
                   << "implementation=" << model.implementationNote << "\n\n";
        }
        output.close();
        if (!output) {
            result.error = "failed while writing scientific references";
            return result;
        }
        std::error_code error;
        std::filesystem::rename(temporary, destination, error);
        if (error) result.error = "cannot replace reference file: " + error.message();
    } catch (const std::exception& exception) {
        result.error = exception.what();
    }
    return result;
}

} // namespace statistics_archive
