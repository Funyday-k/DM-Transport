#include "transport/physics/ScatteringPhysics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using transport::physics::SdProtonModel;
using transport::physics::SdScatteringRates;
using transport::physics::SolarBackground;
using transport::physics::SolarTarget;
using transport::physics::direct_sd_proton_scattering_rates;
using transport::physics::sd_proton_nucleus_cross_section_cm2;

constexpr double protonMassGeV = 0.9382720813;
constexpr double boltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;
constexpr double solarRadiusCm = 6.957e10;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, const std::string& message,
                   double relative_tolerance,
                   double absolute_tolerance = 0.0) {
    require(std::isfinite(actual) && std::isfinite(expected),
            message + ": nonfinite comparison");
    const double scale = std::max(std::abs(actual), std::abs(expected));
    require(std::abs(actual - expected) <=
                std::max(absolute_tolerance, relative_tolerance * scale),
            message + ": values differ");
}

template <typename Exception, typename Callable>
void require_throws(Callable operation, const std::string& message) {
    try {
        operation();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message + ": expected exception was not thrown");
}

long double reduced_mass(long double first, long double second) {
    return first * second / (first + second);
}

long double source_convention_cross_section(
    const SdProtonModel& model,
    const SolarTarget& target) {
    if (target.spin == 0.0) {
        return 0.0L;
    }
    const long double mu_p = reduced_mass(
        model.dark_matter_mass_GeV, protonMassGeV);
    const long double mu_A = reduced_mass(
        model.dark_matter_mass_GeV, target.mass_GeV);
    return static_cast<long double>(model.proton_cross_section_cm2) *
        (4.0L / 3.0L) * (mu_A * mu_A) / (mu_p * mu_p) *
        (target.spin + 1.0L) / target.spin *
        target.proton_spin * target.proton_spin;
}

void test_sd_cross_section_convention(const SolarBackground& background) {
    const SdProtonModel model{0.1, 1.0e-34};
    for (std::size_t index = 0; index < background.target_count(); ++index) {
        const SolarTarget& target = background.target(index);
        const double actual =
            sd_proton_nucleus_cross_section_cm2(model, target);
        const double expected = static_cast<double>(
            source_convention_cross_section(model, target));
        require_close(actual, expected,
                      "fixed obscura SD nucleus convention", 3.0e-14,
                      expected == 0.0 ? 0.0 : 1.0e-300);
    }
    require_close(sd_proton_nucleus_cross_section_cm2(
                      model, background.target(0)),
                  model.proton_cross_section_cm2,
                  "H-1 cross section must equal sigma_p", 3.0e-14);
    require(sd_proton_nucleus_cross_section_cm2(
                model, background.target(1)) == 0.0,
            "spin-zero He-4 must have zero SD cross section");
    require(sd_proton_nucleus_cross_section_cm2(
                model, background.target(5)) == 0.0,
            "N-14 has zero proton-spin response in the fixed nuclear table");
    require(sd_proton_nucleus_cross_section_cm2(
                model, background.target(2)) > 0.0 &&
                sd_proton_nucleus_cross_section_cm2(
                    model, background.target(13)) > 0.0,
            "He-3 and Na-23 must remain active non-hydrogen targets");
}

void test_rate_breakdown(const SolarBackground& background) {
    const SdProtonModel model{0.1, 1.0e-34};
    const SdScatteringRates rates = direct_sd_proton_scattering_rates(
        background, model, 0.0, 0.0);
    require(rates.target_rates.size() == background.target_count(),
            "rate breakdown must retain all 63 targets");

    double direct_sum = 0.0;
    for (std::size_t index = 0; index < rates.target_rates.size(); ++index) {
        const auto& entry = rates.target_rates[index];
        require(entry.target_index == index,
                "rate breakdown target order changed");
        require(std::isfinite(entry.nucleus_cross_section_cm2) &&
                    entry.nucleus_cross_section_cm2 >= 0.0 &&
                    std::isfinite(entry.mean_relative_speed_cm_s) &&
                    entry.mean_relative_speed_cm_s > 0.0 &&
                    std::isfinite(entry.rate_s_inv) &&
                    entry.rate_s_inv >= 0.0,
                "target rate fields must be finite and nonnegative");
        direct_sum += entry.rate_s_inv;
    }
    require(rates.total_rate_s_inv == direct_sum,
            "total rate must be the source-order sum of target rates");
    require(rates.total_rate_s_inv > rates.target_rates[0].rate_s_inv,
            "proton-only coupling must still include non-hydrogen targets");
    require(rates.target_rates[2].rate_s_inv > 0.0 &&
                rates.target_rates[13].rate_s_inv > 0.0,
            "He-3 and Na-23 must contribute to the direct rate");
    require(rates.target_rates[5].rate_s_inv == 0.0,
            "N-14 must expose its zero proton-spin contribution");

    const long double temperature_K = 1.549e7L;
    const long double hydrogen_density_cm3 = 3.2580314157719068e25L;
    const long double thermal_speed_cm_s = speedOfLightCmS * std::sqrt(
        2.0L * boltzmannGeVPerK * temperature_K / protonMassGeV);
    const long double expected_hydrogen_rate =
        hydrogen_density_cm3 * model.proton_cross_section_cm2 *
        2.0L * thermal_speed_cm_s / std::sqrt(std::acos(-1.0L));
    require_close(rates.target_rates[0].rate_s_inv,
                  static_cast<double>(expected_hydrogen_rate),
                  "zero-speed analytic H-1 rate", 5.0e-13);
}

void test_rate_limits_and_scaling(const SolarBackground& background) {
    const SdProtonModel reference{0.1, 1.0e-34};
    const SdProtonModel scaled{0.1, 1.0e-36};
    const double radius_cm = 0.37 * solarRadiusCm;
    const double speed_cm_s = 8.0e7;
    const SdScatteringRates reference_rates =
        direct_sd_proton_scattering_rates(
            background, reference, radius_cm, speed_cm_s);
    const SdScatteringRates scaled_rates = direct_sd_proton_scattering_rates(
        background, scaled, radius_cm, speed_cm_s);
    require_close(reference_rates.total_rate_s_inv,
                  100.0 * scaled_rates.total_rate_s_inv,
                  "total rate must scale linearly with sigma_p", 5.0e-14);
    for (std::size_t index = 0;
         index < reference_rates.target_rates.size(); ++index) {
        const double reference_rate =
            reference_rates.target_rates[index].rate_s_inv;
        const double scaled_rate = scaled_rates.target_rates[index].rate_s_inv;
        if (reference_rate == 0.0) {
            require(scaled_rate == 0.0,
                    "zero target rate must remain zero under scaling");
        } else {
            require_close(reference_rate, 100.0 * scaled_rate,
                          "target rate must scale linearly with sigma_p",
                          5.0e-14);
        }
    }

    const SdScatteringRates zero = direct_sd_proton_scattering_rates(
        background, SdProtonModel{0.1, 0.0}, radius_cm, 0.0);
    require(zero.total_rate_s_inv == 0.0,
            "zero proton cross section must give zero total rate");
    for (const auto& entry : zero.target_rates) {
        require(entry.nucleus_cross_section_cm2 == 0.0 &&
                    entry.rate_s_inv == 0.0 &&
                    entry.mean_relative_speed_cm_s > 0.0,
                "zero-cross-section target diagnostics are inconsistent");
    }

    const SdScatteringRates exterior = direct_sd_proton_scattering_rates(
        background, reference, 1.01 * solarRadiusCm, speed_cm_s);
    require(exterior.total_rate_s_inv == 0.0,
            "exterior total rate must vanish");
    for (const auto& entry : exterior.target_rates) {
        require(entry.rate_s_inv == 0.0 &&
                    entry.mean_relative_speed_cm_s == 0.0,
                "exterior target rate diagnostics must vanish");
    }
}

void test_invalid_inputs(const SolarBackground& background) {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double mass : {0.0, -1.0, 0.5, infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            direct_sd_proton_scattering_rates(
                background, SdProtonModel{mass, 1.0e-34}, 0.0, 0.0);
        }, "invalid dark-matter mass");
    }
    for (double cross_section : {-1.0, infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            direct_sd_proton_scattering_rates(
                background, SdProtonModel{0.1, cross_section}, 0.0, 0.0);
        }, "invalid proton cross section");
    }
    for (double radius : {-1.0, infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            direct_sd_proton_scattering_rates(
                background, SdProtonModel{0.1, 1.0e-34}, radius, 0.0);
        }, "invalid radius");
    }
    for (double speed : {-1.0, infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            direct_sd_proton_scattering_rates(
                background, SdProtonModel{0.1, 1.0e-34}, 0.0, speed);
        }, "invalid speed");
    }
    require_throws<std::overflow_error>([&] {
        sd_proton_nucleus_cross_section_cm2(
            SdProtonModel{0.1, std::numeric_limits<double>::max()},
            background.target(0));
    }, "unrepresentable proton coupling");

    SolarTarget invalid_target = background.target(0);
    invalid_target.spin = -0.5;
    require_throws<std::invalid_argument>([&] {
        sd_proton_nucleus_cross_section_cm2(
            SdProtonModel{0.1, 1.0e-34}, invalid_target);
    }, "invalid target metadata");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "usage: test_scattering_rate SOLAR_MODEL NUCLEAR_DATA");
        }
        const SolarBackground background(argv[1], argv[2]);
        test_sd_cross_section_convention(background);
        test_rate_breakdown(background);
        test_rate_limits_and_scaling(background);
        test_invalid_inputs(background);
        std::cout << "scattering rate tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "scattering rate test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
