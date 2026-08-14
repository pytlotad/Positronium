#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
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
         2026, "", "", "2026-08-11"}
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
         "The trajectory result fits a single normalization C_R to this fixed shape."}
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
