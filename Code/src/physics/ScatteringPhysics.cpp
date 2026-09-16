#include "transport/physics/ScatteringPhysics.hpp"

#include <algorithm>
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
constexpr std::size_t referenceSolarTargetCount = 63;
constexpr unsigned long targetVelocityMaxRejectionAttempts = 10000UL;

using Velocity3 = std::array<double, 3>;

double checked_result(double speed_cm_s) {
    if (!std::isfinite(speed_cm_s)) {
        throw std::overflow_error("mean relative speed is not representable");
    }
    return speed_cm_s;
}

double compute_thermal_speed_cm_s(double temperature_K,
                                  double target_mass_GeV) {
    if (!std::isfinite(temperature_K) ||
        !std::isfinite(target_mass_GeV) ||
        temperature_K <= 0.0 || target_mass_GeV <= 0.0) {
        throw std::invalid_argument(
            "thermal speed requires finite T > 0 and mass > 0");
    }
    const double temperature_mass_ratio = temperature_K / target_mass_GeV;
    if (!std::isfinite(temperature_mass_ratio) ||
        temperature_mass_ratio <= 0.0) {
        throw std::invalid_argument(
            "temperature-to-mass ratio must be finite and positive");
    }
    // Split the square root so a representable thermal speed is not lost when
    // the product inside sqrt would underflow first.
    const double speed = speedOfLightCmS *
        std::sqrt(2.0 * legacyBoltzmannGeVPerK) *
        std::sqrt(temperature_mass_ratio);
    if (!std::isfinite(speed) || speed <= 0.0) {
        throw std::overflow_error("thermal speed is not representable");
    }
    return speed;
}

double sample_uniform(std::mt19937& rng,
                      double minimum = 0.0,
                      double maximum = 1.0) {
    std::uniform_real_distribution<double> distribution(minimum, maximum);
    return distribution(rng);
}

double positive_unit_uniform(std::mt19937& rng) {
    return std::max(sample_uniform(rng),
                    std::numeric_limits<double>::min());
}

double vector_norm(const Velocity3& vector) {
    return std::sqrt(vector[0] * vector[0] +
                     vector[1] * vector[1] +
                     vector[2] * vector[2]);
}

Velocity3 vector_cross(const Velocity3& left, const Velocity3& right) {
    return Velocity3{{
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
    }};
}

Velocity3 normalized(const Velocity3& vector) {
    const double norm = vector_norm(vector);
    if (!std::isfinite(norm) || norm <= 0.0) {
        throw std::runtime_error("cannot normalize a zero or nonfinite vector");
    }
    return Velocity3{{
        vector[0] / norm,
        vector[1] / norm,
        vector[2] / norm
    }};
}

Velocity3 unit_vector_at_angle_from_axis(double cosine,
                                         double phi,
                                         const Velocity3& axis) {
    if (!std::isfinite(cosine) || cosine < -1.0 - 1.0e-12 ||
        cosine > 1.0 + 1.0e-12 || !std::isfinite(phi)) {
        throw std::runtime_error("sampled target direction is invalid");
    }
    const Velocity3 e3 = normalized(axis);
    const Velocity3 reference = std::fabs(e3[2]) < 0.9
        ? Velocity3{{0.0, 0.0, 1.0}}
        : Velocity3{{1.0, 0.0, 0.0}};
    const Velocity3 e1 = normalized(vector_cross(reference, e3));
    const Velocity3 e2 = normalized(vector_cross(e3, e1));
    cosine = std::max(-1.0, std::min(1.0, cosine));
    const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    const double cosine_phi = std::cos(phi);
    const double sine_phi = std::sin(phi);
    return Velocity3{{
        cosine * e3[0] + sine * cosine_phi * e1[0] +
            sine * sine_phi * e2[0],
        cosine * e3[1] + sine * cosine_phi * e1[1] +
            sine * sine_phi * e2[1],
        cosine * e3[2] + sine * cosine_phi * e1[2] +
            sine * sine_phi * e2[2]
    }};
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
    if (!std::isfinite(dm_speed_cm_s) || dm_speed_cm_s < 0.0) {
        throw std::invalid_argument(
            "mean relative speed requires finite T > 0, mass > 0, and v >= 0");
    }

    const double thermal_speed_cm_s =
        compute_thermal_speed_cm_s(temperature_K, target_mass_GeV);

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

std::size_t sample_sd_proton_target_index(
    const SdScatteringRates& rates,
    std::mt19937& rng) {
    if (rates.target_rates.size() != referenceSolarTargetCount ||
        !std::isfinite(rates.total_rate_s_inv)) {
        throw std::invalid_argument(
            "target sampling requires 63 rates and a finite total");
    }
    if (rates.total_rate_s_inv <= 0.0) {
        throw std::domain_error(
            "a nuclear target cannot be sampled from a nonpositive total rate");
    }

    double source_order_sum = 0.0;
    for (std::size_t index = 0; index < rates.target_rates.size(); ++index) {
        const SdTargetScatteringRate& entry = rates.target_rates[index];
        if (entry.target_index != index ||
            !std::isfinite(entry.rate_s_inv) || entry.rate_s_inv < 0.0) {
            throw std::invalid_argument(
                "target rates must be finite, nonnegative, and in source order");
        }
        if (entry.rate_s_inv >
            std::numeric_limits<double>::max() - source_order_sum) {
            throw std::overflow_error("target-rate sum is not representable");
        }
        source_order_sum += entry.rate_s_inv;
    }
    const double total_scale = std::max(
        std::fabs(source_order_sum), std::fabs(rates.total_rate_s_inv));
    const double total_tolerance = 16.0 *
        std::numeric_limits<double>::epsilon() * total_scale;
    if (source_order_sum <= 0.0 ||
        std::fabs(source_order_sum - rates.total_rate_s_inv) >
            total_tolerance) {
        throw std::invalid_argument(
            "target-rate total does not match its source-order breakdown");
    }

    const double xi = sample_uniform(rng);
    double cumulative_probability = 0.0;
    for (std::size_t index = 0; index < rates.target_rates.size(); ++index) {
        cumulative_probability +=
            rates.target_rates[index].rate_s_inv / rates.total_rate_s_inv;
        if (cumulative_probability > xi ||
            index + 1 == rates.target_rates.size()) {
            return index;
        }
    }
    throw std::runtime_error("no nuclear target could be sampled");
}

std::array<double, 3> sample_collision_conditioned_target_velocity_cm_s(
    double temperature_K,
    double target_mass_GeV,
    const std::array<double, 3>& dm_velocity_cm_s,
    std::mt19937& rng) {
    for (double component : dm_velocity_cm_s) {
        if (!std::isfinite(component)) {
            throw std::invalid_argument(
                "incoming dark-matter velocity must be finite");
        }
    }
    if (!std::isfinite(temperature_K) || temperature_K <= 0.0 ||
        !std::isfinite(target_mass_GeV) || target_mass_GeV <= 0.0) {
        throw std::invalid_argument(
            "target-velocity sampling requires finite T > 0 and mass > 0");
    }
    const double temperature_GeV = legacyBoltzmannGeVPerK * temperature_K;
    const double kappa = std::sqrt(
        target_mass_GeV / 2.0 / temperature_GeV);
    const Velocity3 dm_velocity_natural{{
        dm_velocity_cm_s[0] / speedOfLightCmS,
        dm_velocity_cm_s[1] / speedOfLightCmS,
        dm_velocity_cm_s[2] / speedOfLightCmS
    }};
    const double dm_speed_natural = vector_norm(dm_velocity_natural);
    if (!std::isfinite(kappa) || kappa <= 0.0) {
        throw std::overflow_error(
            "target inverse thermal speed is not representable");
    }
    if (!std::isfinite(dm_speed_natural) || dm_speed_natural <= 0.0) {
        throw std::invalid_argument(
            "legacy target-velocity sampling requires nonzero DM speed");
    }
    const double y = kappa * dm_speed_natural;
    const double maximum_safe_ratio =
        0.25 * std::sqrt(std::numeric_limits<double>::max());
    if (!std::isfinite(y) || y <= 0.0 || y > maximum_safe_ratio) {
        throw std::overflow_error(
            "dimensionless incoming speed is not representable");
    }

    // Romano-Walsh rejection sampler as used by the fixed reference.  Keep
    // the initial x=y, mu=1 sentinel because its failed acceptance draw is
    // part of the legacy std::mt19937 call order.
    double x = y;
    double cosine = 1.0;
    unsigned long rejection_attempts = 0;
    const auto relative_speed_ratio = [&]() {
        const double x_squared = x * x;
        const double y_squared = y * y;
        const double relative_speed_squared_raw =
            x_squared + y_squared - 2.0 * x * y * cosine;
        const double denominator = x + y;
        const double roundoff_scale = std::max(
            1.0, x_squared + y_squared + 2.0 * std::fabs(x * y));
        if (!std::isfinite(x) || x < 0.0 ||
            !std::isfinite(cosine) || cosine < -1.0 || cosine > 1.0 ||
            !std::isfinite(relative_speed_squared_raw) ||
            relative_speed_squared_raw < -1.0e-12 * roundoff_scale ||
            !std::isfinite(denominator) || denominator <= 0.0) {
            throw std::runtime_error(
                "target-velocity rejection state is invalid");
        }
        const double relative_speed_squared =
            std::max(0.0, relative_speed_squared_raw);
        const double ratio =
            std::sqrt(relative_speed_squared) / denominator;
        if (!std::isfinite(ratio) || ratio < 0.0 ||
            ratio > 1.0 + 1.0e-12) {
            throw std::runtime_error(
                "target-velocity rejection ratio is invalid");
        }
        return std::min(1.0, ratio);
    };

    const double pi = std::acos(-1.0);
    while (relative_speed_ratio() < sample_uniform(rng)) {
        ++rejection_attempts;
        if (rejection_attempts > targetVelocityMaxRejectionAttempts) {
            throw std::runtime_error(
                "target-velocity rejection sampling exceeded maximum attempts");
        }
        cosine = sample_uniform(rng, -1.0, 1.0);
        const double mixture_draw = sample_uniform(rng);
        if (mixture_draw < 2.0 / (std::sqrt(pi) * y + 2.0)) {
            const double draw_2 = positive_unit_uniform(rng);
            const double draw_3 = positive_unit_uniform(rng);
            x = std::sqrt(-std::log(draw_2) - std::log(draw_3));
        } else {
            const double draw_2 = positive_unit_uniform(rng);
            const double draw_3 = sample_uniform(rng);
            const double draw_4 = positive_unit_uniform(rng);
            const double cosine_draw = std::cos(0.5 * pi * draw_3);
            x = std::sqrt(-std::log(draw_2) -
                          cosine_draw * cosine_draw * std::log(draw_4));
        }
    }

    const double target_speed_cm_s = x / kappa * speedOfLightCmS;
    if (!std::isfinite(target_speed_cm_s)) {
        throw std::overflow_error("sampled target speed is not representable");
    }
    const double phi = sample_uniform(rng, 0.0, 2.0 * pi);
    const Velocity3 direction =
        unit_vector_at_angle_from_axis(cosine, phi, dm_velocity_natural);
    const Velocity3 velocity{{
        target_speed_cm_s * direction[0],
        target_speed_cm_s * direction[1],
        target_speed_cm_s * direction[2]
    }};
    if (!std::isfinite(vector_norm(velocity))) {
        throw std::overflow_error(
            "sampled target velocity is not representable");
    }
    return velocity;
}

}  // namespace physics
}  // namespace transport
