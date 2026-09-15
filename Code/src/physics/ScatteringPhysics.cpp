#include "transport/physics/ScatteringPhysics.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace transport {
namespace physics {
namespace {

// Port provenance: Code/provenance/reference_physics.json, unit mean_relative_speed.
// Preserve the reference implementation's constant for T02 numerical parity.
constexpr double legacyBoltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;
constexpr double legacyProtonMassGeV = 0.9382720813;
constexpr double referenceDarkMatterMassGeV = 0.1;
constexpr double smallSpeedRatio = 1.0e-4;
constexpr double largeSpeedRatio = 26.0;

double checked_result(double speed_cm_s) {
    if (!std::isfinite(speed_cm_s)) {
        throw std::overflow_error("mean relative speed is not representable");
    }
    return speed_cm_s;
}

void validate_sd_proton_model(const SdProtonModel& model) {
    if (!std::isfinite(model.dark_matter_mass_GeV) ||
        model.dark_matter_mass_GeV != referenceDarkMatterMassGeV ||
        !std::isfinite(model.proton_cross_section_cm2) ||
        model.proton_cross_section_cm2 < 0.0) {
        throw std::invalid_argument(
            "SD-proton MVP requires mass = 0.1 GeV and finite cross section >= 0");
    }
}

void validate_solar_target(const SolarTarget& target) {
    if (target.atomic_number == 0 ||
        target.mass_number < target.atomic_number ||
        !std::isfinite(target.mass_GeV) || target.mass_GeV <= 0.0 ||
        !std::isfinite(target.spin) || target.spin < 0.0 ||
        !std::isfinite(target.proton_spin) ||
        !std::isfinite(target.neutron_spin)) {
        throw std::invalid_argument("invalid solar target metadata");
    }
}

double reference_reduced_mass_GeV(double first_mass_GeV,
                                  double second_mass_GeV) {
    // Preserve libphysica's evaluation order for the physical MVP range.
    const double reduced_mass = first_mass_GeV * second_mass_GeV /
        (first_mass_GeV + second_mass_GeV);
    if (!std::isfinite(reduced_mass) || reduced_mass <= 0.0) {
        throw std::overflow_error("reduced mass is not representable");
    }
    return reduced_mass;
}

}  // namespace

double mean_relative_speed_cm_s(double temperature_K,
                                double target_mass_GeV,
                                double dm_speed_cm_s) {
    if (!std::isfinite(temperature_K) ||
        !std::isfinite(target_mass_GeV) ||
        !std::isfinite(dm_speed_cm_s) || temperature_K <= 0.0 ||
        target_mass_GeV <= 0.0 || dm_speed_cm_s < 0.0) {
        throw std::invalid_argument(
            "mean relative speed requires finite T > 0, mass > 0, and v >= 0");
    }

    const double temperature_mass_ratio = temperature_K / target_mass_GeV;
    if (!std::isfinite(temperature_mass_ratio) ||
        temperature_mass_ratio <= 0.0) {
        throw std::invalid_argument(
            "temperature-to-mass ratio must be finite and positive");
    }
    // Split the square root so a representable thermal speed is not lost when
    // the product inside sqrt would underflow first.
    const double thermal_speed_cm_s = speedOfLightCmS *
        std::sqrt(2.0 * legacyBoltzmannGeVPerK) *
        std::sqrt(temperature_mass_ratio);
    if (!std::isfinite(thermal_speed_cm_s) || thermal_speed_cm_s <= 0.0) {
        throw std::overflow_error("thermal speed is not representable");
    }

    const double x = dm_speed_cm_s / thermal_speed_cm_s;
    const double inverse_sqrt_pi = 1.0 / std::sqrt(std::acos(-1.0));

    if (x < smallSpeedRatio) {
        const double x2 = x * x;
        // Stable expansion through x^6 of
        // exp(-x^2)/sqrt(pi) + (x + 1/(2x))*erf(x).
        const double series = 2.0 * inverse_sqrt_pi *
            (1.0 + x2 * (1.0 / 3.0 +
                         x2 * (-1.0 / 30.0 + x2 / 210.0)));
        return checked_result(thermal_speed_cm_s * series);
    }

    if (x > largeSpeedRatio) {
        // erf(x) is unity to double precision here. This arrangement also
        // avoids squaring the thermal speed or multiplying by a possibly
        // overflowing x.
        return checked_result(dm_speed_cm_s +
                              0.5 * thermal_speed_cm_s / x);
    }

    const double closed_form = std::exp(-x * x) * inverse_sqrt_pi +
        (x + 0.5 / x) * std::erf(x);
    return checked_result(thermal_speed_cm_s * closed_form);
}

double sd_proton_nucleus_cross_section_cm2(const SdProtonModel& model,
                                           const SolarTarget& target) {
    validate_sd_proton_model(model);
    validate_solar_target(target);
    if (target.spin == 0.0 || model.proton_cross_section_cm2 == 0.0) {
        return 0.0;
    }

    // obscura's SD particle uses prefactor 3 when converting sigma_p to fp,
    // followed by 4 mu_A^2/pi * (J+1)/J * (fp Sp + fn Sn)^2.  The MVP fixes
    // fn=0 and keeps this evaluation structure in explicit cm^2/GeV units.
    const double pi = std::acos(-1.0);
    const double proton_reduced_mass = reference_reduced_mass_GeV(
        model.dark_matter_mass_GeV, legacyProtonMassGeV);
    const double proton_coupling =
        std::sqrt(pi * model.proton_cross_section_cm2 / 3.0) /
        proton_reduced_mass;
    if (!std::isfinite(proton_coupling) ||
        (model.proton_cross_section_cm2 > 0.0 && proton_coupling <= 0.0)) {
        throw std::overflow_error("SD proton coupling is not representable");
    }

    const double target_reduced_mass = reference_reduced_mass_GeV(
        model.dark_matter_mass_GeV, target.mass_GeV);
    const double neutron_coupling = 0.0;
    const double spin_response = proton_coupling * target.proton_spin +
        neutron_coupling * target.neutron_spin;
    const double cross_section =
        4.0 * std::pow(target_reduced_mass, 2.0) / pi *
        (target.spin + 1.0) / target.spin *
        std::pow(spin_response, 2.0);
    if (!std::isfinite(cross_section) || cross_section < 0.0) {
        throw std::overflow_error(
            "SD nucleus cross section is not representable");
    }
    return cross_section;
}

SdScatteringRates direct_sd_proton_scattering_rates(
    const SolarBackground& background,
    const SdProtonModel& model,
    double radius_cm,
    double dm_speed_cm_s) {
    validate_sd_proton_model(model);
    if (!std::isfinite(radius_cm) || radius_cm < 0.0 ||
        !std::isfinite(dm_speed_cm_s) || dm_speed_cm_s < 0.0) {
        throw std::invalid_argument(
            "direct scattering rate requires finite radius and speed >= 0");
    }

    SdScatteringRates result;
    result.target_rates.reserve(background.target_count());
    result.total_rate_s_inv = 0.0;
    const bool exterior = radius_cm > background.solar_radius_cm();
    const double temperature_K = exterior
        ? 0.0
        : background.temperature_K(radius_cm);

    for (std::size_t target_index = 0;
         target_index < background.target_count(); ++target_index) {
        const SolarTarget& target = background.target(target_index);
        const double cross_section_cm2 =
            sd_proton_nucleus_cross_section_cm2(model, target);
        const double mean_speed_cm_s = exterior
            ? 0.0
            : mean_relative_speed_cm_s(
                  temperature_K, target.mass_GeV, dm_speed_cm_s);
        const double number_density_cm3 =
            background.number_density_cm3(target_index, radius_cm);
        if (!std::isfinite(number_density_cm3) || number_density_cm3 < 0.0) {
            throw std::runtime_error(
                "solar target number density is invalid");
        }
        const double rate_s_inv =
            number_density_cm3 * cross_section_cm2 * mean_speed_cm_s;
        if (!std::isfinite(rate_s_inv) || rate_s_inv < 0.0) {
            throw std::overflow_error(
                "SD target scattering rate is not representable");
        }
        if (rate_s_inv >
            std::numeric_limits<double>::max() - result.total_rate_s_inv) {
            throw std::overflow_error(
                "total SD scattering rate is not representable");
        }
        result.target_rates.push_back(SdTargetScatteringRate{
            target_index,
            cross_section_cm2,
            mean_speed_cm_s,
            rate_s_inv
        });
        result.total_rate_s_inv += rate_s_inv;
    }
    return result;
}

}  // namespace physics
}  // namespace transport
