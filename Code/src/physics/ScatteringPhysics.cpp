#include "transport/physics/ScatteringPhysics.hpp"

#include <cmath>
#include <stdexcept>

namespace transport {
namespace physics {
namespace {

// Port provenance: Code/provenance/reference_physics.json, unit mean_relative_speed.
// Preserve the reference implementation's constant for T02 numerical parity.
constexpr double legacyBoltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;
constexpr double smallSpeedRatio = 1.0e-4;
constexpr double largeSpeedRatio = 26.0;

double checked_result(double speed_cm_s) {
    if (!std::isfinite(speed_cm_s)) {
        throw std::overflow_error("mean relative speed is not representable");
    }
    return speed_cm_s;
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

}  // namespace physics
}  // namespace transport
