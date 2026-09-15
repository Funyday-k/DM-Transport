#include "transport/physics/ScatteringPhysics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using transport::physics::mean_relative_speed_cm_s;

constexpr double legacyBoltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, const std::string& message,
                   double relative_tolerance, double absolute_tolerance = 0.0) {
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

double thermal_speed(double temperature_K, double target_mass_GeV) {
    return speedOfLightCmS * std::sqrt(
        2.0 * legacyBoltzmannGeVPerK * temperature_K / target_mass_GeV);
}

long double independent_closed_form(long double thermal_speed_cm_s,
                                    long double x) {
    const long double pi = std::acos(-1.0L);
    if (x == 0.0L) {
        return 2.0L * thermal_speed_cm_s / std::sqrt(pi);
    }
    return thermal_speed_cm_s *
        (std::exp(-x * x) / std::sqrt(pi) +
         (x + 0.5L / x) * std::erf(x));
}

void test_zero_speed_limit() {
    const double temperature_K = 1.0e7;
    const double target_mass_GeV = 0.9382720813;
    const double v_thermal = thermal_speed(temperature_K, target_mass_GeV);
    const double expected = 2.0 * v_thermal / std::sqrt(std::acos(-1.0));
    require_close(mean_relative_speed_cm_s(temperature_K,
                                           target_mass_GeV, 0.0),
                  expected, "zero-speed Maxwellian mean", 3.0e-15);
}

void test_legacy_golden_value() {
    const double temperature_K = 1.0e7;
    const double proton_mass_GeV = 0.9382720813;
    const double dark_matter_speed_cm_s = 1.0e-3 * speedOfLightCmS;
    const double dimensionless_mean = mean_relative_speed_cm_s(
        temperature_K, proton_mass_GeV, dark_matter_speed_cm_s) /
        speedOfLightCmS;
    require_close(dimensionless_mean, 0.0017928027334128615,
                  "legacy proton golden value", 1.0e-12);
}

void test_representable_subnormal_temperature_mass_ratios() {
    require_close(mean_relative_speed_cm_s(1.0e-310, 1.0, 0.0),
                  1.40435572778300447e-151,
                  "subnormal temperature-to-mass ratio", 2.0e-12);
    require_close(mean_relative_speed_cm_s(1.0e-311, 1.0, 0.0),
                  4.44096274489759081e-152,
                  "smaller subnormal temperature-to-mass ratio", 2.0e-12);
}

void test_branch_boundaries() {
    const double temperature_K = 1.0e7;
    const double target_mass_GeV = 1.0;
    const double v_thermal = thermal_speed(temperature_K, target_mass_GeV);
    const std::vector<double> ratios{
        5.0e-5,
        std::nextafter(1.0e-4, 0.0), 1.0e-4,
        std::nextafter(1.0e-4, 1.0),
        std::nextafter(26.0, 0.0), 26.0,
        std::nextafter(26.0, std::numeric_limits<double>::infinity()),
        30.0
    };
    for (double ratio : ratios) {
        const double actual = mean_relative_speed_cm_s(
            temperature_K, target_mass_GeV, ratio * v_thermal);
        const double expected = static_cast<double>(independent_closed_form(
            static_cast<long double>(v_thermal),
            static_cast<long double>(ratio)));
        require_close(actual, expected, "branch-boundary value", 8.0e-15);
    }
}

void test_invalid_inputs() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (double temperature_K : {0.0, -1.0, infinity, -infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            mean_relative_speed_cm_s(temperature_K, 1.0, 0.0);
        }, "invalid temperature");
    }
    for (double target_mass_GeV : {0.0, -1.0, infinity, -infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            mean_relative_speed_cm_s(1.0, target_mass_GeV, 0.0);
        }, "invalid target mass");
    }
    for (double speed_cm_s : {-1.0, infinity, -infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            mean_relative_speed_cm_s(1.0, 1.0, speed_cm_s);
        }, "invalid dark-matter speed");
    }
    require_throws<std::invalid_argument>([&] {
        mean_relative_speed_cm_s(std::numeric_limits<double>::denorm_min(),
                                 std::numeric_limits<double>::max(), 0.0);
    }, "temperature-to-mass ratio underflow");
    require_throws<std::invalid_argument>([&] {
        mean_relative_speed_cm_s(std::numeric_limits<double>::max(),
                                 std::numeric_limits<double>::denorm_min(),
                                 0.0);
    }, "temperature-to-mass ratio overflow");
}

}  // namespace

int main() {
    try {
        test_zero_speed_limit();
        test_legacy_golden_value();
        test_representable_subnormal_temperature_mass_ratios();
        test_branch_boundaries();
        test_invalid_inputs();
        std::cout << "scattering physics tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "scattering physics test failed: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
