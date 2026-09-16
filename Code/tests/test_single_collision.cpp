#include "transport/physics/ScatteringPhysics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace {

using transport::physics::CartesianVelocityCmS;
using transport::physics::CollisionSample;
using transport::physics::SdProtonModel;
using transport::physics::SolarBackground;
using transport::physics::direct_sd_proton_scattering_rates;
using transport::physics::elastic_outgoing_dm_velocity_cm_s;
using transport::physics::sample_collision_conditioned_target_velocity_cm_s;
using transport::physics::sample_sd_proton_collision;
using transport::physics::sample_sd_proton_target_index;

using Velocity3 = CartesianVelocityCmS;

constexpr double darkMatterMassGeV = 0.1;
constexpr double protonMassGeV = 0.9382720813;
constexpr double solarRadiusCm = 6.957e10;
constexpr double legacyBoltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
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

template <typename Exception, typename Callable>
void require_throws_without_consuming_rng(Callable operation,
                                          std::mt19937& rng,
                                          const std::string& message) {
    const std::mt19937 initial_state = rng;
    require_throws<Exception>(operation, message);
    require(rng == initial_state,
            message + ": validation consumed random state");
}

double dot(const Velocity3& left, const Velocity3& right) {
    return left[0] * right[0] + left[1] * right[1] +
        left[2] * right[2];
}

double norm(const Velocity3& vector) {
    return std::hypot(std::hypot(vector[0], vector[1]), vector[2]);
}

long double squared_norm_long_double(const Velocity3& vector) {
    long double result = 0.0L;
    for (double component : vector) {
        const long double value = static_cast<long double>(component);
        result += value * value;
    }
    return result;
}

Velocity3 add(const Velocity3& left, const Velocity3& right) {
    return Velocity3{{
        left[0] + right[0],
        left[1] + right[1],
        left[2] + right[2]
    }};
}

Velocity3 subtract(const Velocity3& left, const Velocity3& right) {
    return Velocity3{{
        left[0] - right[0],
        left[1] - right[1],
        left[2] - right[2]
    }};
}

Velocity3 multiply(double scalar, const Velocity3& vector) {
    return Velocity3{{
        scalar * vector[0],
        scalar * vector[1],
        scalar * vector[2]
    }};
}

Velocity3 cross(const Velocity3& left, const Velocity3& right) {
    return Velocity3{{
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0]
    }};
}

Velocity3 normalized(const Velocity3& vector) {
    const double magnitude = norm(vector);
    require(std::isfinite(magnitude) && magnitude > 0.0,
            "test attempted to normalize an invalid vector");
    return multiply(1.0 / magnitude, vector);
}

Velocity3 rotate_cyclic(const Velocity3& vector) {
    // (x,y,z) -> (z,x,y) is a proper 120-degree rotation about (1,1,1).
    return Velocity3{{vector[2], vector[0], vector[1]}};
}

void require_close(double actual,
                   double expected,
                   double relative_tolerance,
                   double absolute_tolerance,
                   const std::string& message) {
    const double scale = std::max(std::fabs(actual), std::fabs(expected));
    require(std::fabs(actual - expected) <=
                absolute_tolerance + relative_tolerance * scale,
            message);
}

void require_vector_close(const Velocity3& actual,
                          const Velocity3& expected,
                          double relative_tolerance,
                          double absolute_tolerance,
                          const std::string& message) {
    for (std::size_t component = 0; component < 3; ++component) {
        require_close(actual[component], expected[component],
                      relative_tolerance, absolute_tolerance,
                      message + " component " + std::to_string(component));
    }
}

Velocity3 reference_outgoing_dm_velocity(
    long double dark_matter_mass_GeV,
    long double target_mass_GeV,
    const Velocity3& incoming_dm_velocity_cm_s,
    const Velocity3& target_velocity_cm_s,
    const Velocity3& outgoing_direction_unit) {
    long double relative_squared = 0.0L;
    for (std::size_t component = 0; component < 3; ++component) {
        const long double relative =
            static_cast<long double>(incoming_dm_velocity_cm_s[component]) -
            static_cast<long double>(target_velocity_cm_s[component]);
        relative_squared += relative * relative;
    }
    const long double relative_speed = std::sqrt(relative_squared);
    const long double total_mass = dark_matter_mass_GeV + target_mass_GeV;
    Velocity3 result{{}};
    for (std::size_t component = 0; component < 3; ++component) {
        const long double center_of_mass_velocity =
            (dark_matter_mass_GeV *
                 static_cast<long double>(
                     incoming_dm_velocity_cm_s[component]) +
             target_mass_GeV *
                 static_cast<long double>(target_velocity_cm_s[component])) /
            total_mass;
        result[component] = static_cast<double>(
            center_of_mass_velocity +
            target_mass_GeV / total_mass * relative_speed *
                static_cast<long double>(
                    outgoing_direction_unit[component]));
    }
    return result;
}

Velocity3 center_of_mass_velocity(double dark_matter_mass_GeV,
                                  double target_mass_GeV,
                                  const Velocity3& dm_velocity_cm_s,
                                  const Velocity3& target_velocity_cm_s) {
    const double total_mass_GeV =
        dark_matter_mass_GeV + target_mass_GeV;
    return Velocity3{{
        (dark_matter_mass_GeV * dm_velocity_cm_s[0] +
         target_mass_GeV * target_velocity_cm_s[0]) / total_mass_GeV,
        (dark_matter_mass_GeV * dm_velocity_cm_s[1] +
         target_mass_GeV * target_velocity_cm_s[1]) / total_mass_GeV,
        (dark_matter_mass_GeV * dm_velocity_cm_s[2] +
         target_mass_GeV * target_velocity_cm_s[2]) / total_mass_GeV
    }};
}

void require_elastic_invariants(double dark_matter_mass_GeV,
                                double target_mass_GeV,
                                const Velocity3& incoming_dm_velocity_cm_s,
                                const Velocity3& target_velocity_cm_s,
                                const Velocity3& outgoing_direction_unit,
                                const Velocity3& outgoing_dm_velocity_cm_s,
                                const std::string& context) {
    const double total_mass_GeV =
        dark_matter_mass_GeV + target_mass_GeV;
    const double relative_speed_cm_s =
        norm(subtract(incoming_dm_velocity_cm_s, target_velocity_cm_s));
    const Velocity3 center = center_of_mass_velocity(
        dark_matter_mass_GeV, target_mass_GeV,
        incoming_dm_velocity_cm_s, target_velocity_cm_s);
    const double expected_dm_cm_speed =
        target_mass_GeV / total_mass_GeV * relative_speed_cm_s;
    require_close(norm(subtract(outgoing_dm_velocity_cm_s, center)),
                  expected_dm_cm_speed, 3.0e-13, 1.0e-6,
                  context + ": outgoing CM radius changed");

    const Velocity3 expected = reference_outgoing_dm_velocity(
        dark_matter_mass_GeV, target_mass_GeV,
        incoming_dm_velocity_cm_s, target_velocity_cm_s,
        outgoing_direction_unit);
    require_vector_close(outgoing_dm_velocity_cm_s, expected,
                         3.0e-13, 1.0e-6,
                         context + ": independent long-double formula mismatch");

    const Velocity3 outgoing_target_velocity_cm_s = add(
        center,
        multiply(-dark_matter_mass_GeV / total_mass_GeV *
                     relative_speed_cm_s,
                 outgoing_direction_unit));
    for (std::size_t component = 0; component < 3; ++component) {
        const long double initial_momentum =
            static_cast<long double>(dark_matter_mass_GeV) *
                static_cast<long double>(
                    incoming_dm_velocity_cm_s[component]) +
            static_cast<long double>(target_mass_GeV) *
                static_cast<long double>(target_velocity_cm_s[component]);
        const long double final_momentum =
            static_cast<long double>(dark_matter_mass_GeV) *
                static_cast<long double>(
                    outgoing_dm_velocity_cm_s[component]) +
            static_cast<long double>(target_mass_GeV) *
                static_cast<long double>(
                    outgoing_target_velocity_cm_s[component]);
        const long double scale = std::max(
            1.0L,
            static_cast<long double>(dark_matter_mass_GeV) *
                    std::fabs(incoming_dm_velocity_cm_s[component]) +
                static_cast<long double>(target_mass_GeV) *
                    std::fabs(target_velocity_cm_s[component]) +
                static_cast<long double>(dark_matter_mass_GeV) *
                    std::fabs(outgoing_dm_velocity_cm_s[component]) +
                static_cast<long double>(target_mass_GeV) *
                    std::fabs(outgoing_target_velocity_cm_s[component]));
        require(std::fabs(final_momentum - initial_momentum) <=
                    2.0e-12L * scale,
                context + ": momentum is not conserved");
    }
    const long double initial_energy_twice =
        static_cast<long double>(dark_matter_mass_GeV) *
            squared_norm_long_double(incoming_dm_velocity_cm_s) +
        static_cast<long double>(target_mass_GeV) *
            squared_norm_long_double(target_velocity_cm_s);
    const long double final_energy_twice =
        static_cast<long double>(dark_matter_mass_GeV) *
            squared_norm_long_double(outgoing_dm_velocity_cm_s) +
        static_cast<long double>(target_mass_GeV) *
            squared_norm_long_double(outgoing_target_velocity_cm_s);
    require(std::fabs(final_energy_twice - initial_energy_twice) <=
                8.0e-13L * std::max(1.0L, std::fabs(initial_energy_twice)),
            context + ": kinetic energy is not conserved");
}

void test_fixed_angle_kinematics() {
    const Velocity3 incoming_direction =
        normalized(Velocity3{{1.0, 2.0, 3.0}});
    const Velocity3 incoming_dm_velocity_cm_s =
        multiply(8.0e7, incoming_direction);
    const Velocity3 target_at_rest{{0.0, 0.0, 0.0}};

    const Velocity3 forward = elastic_outgoing_dm_velocity_cm_s(
        darkMatterMassGeV, protonMassGeV,
        incoming_dm_velocity_cm_s, target_at_rest, incoming_direction);
    require_vector_close(forward, incoming_dm_velocity_cm_s,
                         3.0e-13, 1.0e-6,
                         "forward scatter must preserve the DM velocity");
    require_elastic_invariants(
        darkMatterMassGeV, protonMassGeV,
        incoming_dm_velocity_cm_s, target_at_rest, incoming_direction,
        forward, "forward scatter");

    const Velocity3 backward_direction = multiply(-1.0, incoming_direction);
    const Velocity3 backward = elastic_outgoing_dm_velocity_cm_s(
        darkMatterMassGeV, protonMassGeV,
        incoming_dm_velocity_cm_s, target_at_rest, backward_direction);
    const double speed_factor =
        (darkMatterMassGeV - protonMassGeV) /
        (darkMatterMassGeV + protonMassGeV);
    require_vector_close(backward,
                         multiply(speed_factor, incoming_dm_velocity_cm_s),
                         3.0e-13, 1.0e-6,
                         "backscatter velocity has the wrong mass factor");
    const double loss_fraction = 1.0 -
        dot(backward, backward) /
            dot(incoming_dm_velocity_cm_s, incoming_dm_velocity_cm_s);
    const double maximum_loss_fraction =
        4.0 * darkMatterMassGeV * protonMassGeV /
        std::pow(darkMatterMassGeV + protonMassGeV, 2.0);
    require_close(loss_fraction, maximum_loss_fraction,
                  2.0e-13, 2.0e-15,
                  "backscatter does not attain the analytic maximum loss");
    require_elastic_invariants(
        darkMatterMassGeV, protonMassGeV,
        incoming_dm_velocity_cm_s, target_at_rest, backward_direction,
        backward, "backscatter");
}

void test_moving_target_covariance() {
    const Velocity3 incoming_dm_velocity_cm_s{{7.0e7, -2.0e7, 3.0e7}};
    const Velocity3 target_velocity_cm_s{{-1.0e7, 4.0e7, 2.0e7}};
    const Velocity3 outgoing_direction =
        normalized(Velocity3{{-0.2, 0.9, 0.3}});
    const double target_mass_GeV = 3.2;

    const Velocity3 outgoing = elastic_outgoing_dm_velocity_cm_s(
        darkMatterMassGeV, target_mass_GeV,
        incoming_dm_velocity_cm_s, target_velocity_cm_s,
        outgoing_direction);
    require_elastic_invariants(
        darkMatterMassGeV, target_mass_GeV,
        incoming_dm_velocity_cm_s, target_velocity_cm_s,
        outgoing_direction, outgoing, "moving-target scatter");

    const Velocity3 rotated = elastic_outgoing_dm_velocity_cm_s(
        darkMatterMassGeV, target_mass_GeV,
        rotate_cyclic(incoming_dm_velocity_cm_s),
        rotate_cyclic(target_velocity_cm_s),
        rotate_cyclic(outgoing_direction));
    require_vector_close(rotated, rotate_cyclic(outgoing),
                         4.0e-13, 1.0e-6,
                         "elastic mapping is not rotationally covariant");

    const Velocity3 boost{{1.0e6, -3.0e6, 2.0e6}};
    const Velocity3 boosted = elastic_outgoing_dm_velocity_cm_s(
        darkMatterMassGeV, target_mass_GeV,
        add(incoming_dm_velocity_cm_s, boost),
        add(target_velocity_cm_s, boost), outgoing_direction);
    require_vector_close(boosted, add(outgoing, boost),
                         4.0e-13, 1.0e-6,
                         "elastic mapping is not Galilean covariant");
}

Velocity3 unit_vector_at_angle(double cosine,
                               double phi,
                               const Velocity3& axis) {
    const Velocity3 e3 = normalized(axis);
    const Velocity3 reference = std::fabs(e3[2]) < 0.9
        ? Velocity3{{0.0, 0.0, 1.0}}
        : Velocity3{{1.0, 0.0, 0.0}};
    const Velocity3 e1 = normalized(cross(reference, e3));
    const Velocity3 e2 = normalized(cross(e3, e1));
    const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
    return add(multiply(cosine, e3),
               add(multiply(sine * std::cos(phi), e1),
                   multiply(sine * std::sin(phi), e2)));
}

void test_full_collision_draw_order(const SolarBackground& background) {
    const SdProtonModel model{darkMatterMassGeV, 1.0e-34};
    const double radius_cm = 0.37 * solarRadiusCm;
    const Velocity3 incoming_dm_velocity_cm_s =
        multiply(8.0e7, normalized(Velocity3{{1.0, -2.0, 3.0}}));
    std::mt19937 full_rng(20260926U);
    std::mt19937 manual_rng = full_rng;

    const CollisionSample full = sample_sd_proton_collision(
        background, model, radius_cm, incoming_dm_velocity_cm_s, full_rng);
    const auto rates = direct_sd_proton_scattering_rates(
        background, model, radius_cm, norm(incoming_dm_velocity_cm_s));
    const std::size_t target_index =
        sample_sd_proton_target_index(rates, manual_rng);
    const auto& target = background.target(target_index);
    const Velocity3 target_velocity_cm_s =
        sample_collision_conditioned_target_velocity_cm_s(
            background.temperature_K(radius_cm), target.mass_GeV,
            incoming_dm_velocity_cm_s, manual_rng);
    std::uniform_real_distribution<double> unit_uniform(0.0, 1.0);
    const double cosine = 2.0 * unit_uniform(manual_rng) - 1.0;
    std::uniform_real_distribution<double> phi_uniform(
        0.0, 2.0 * std::acos(-1.0));
    const double phi = phi_uniform(manual_rng);
    const Velocity3 outgoing_direction = unit_vector_at_angle(
        cosine, phi, incoming_dm_velocity_cm_s);
    const Velocity3 outgoing_dm_velocity_cm_s =
        elastic_outgoing_dm_velocity_cm_s(
            model.dark_matter_mass_GeV, target.mass_GeV,
            incoming_dm_velocity_cm_s, target_velocity_cm_s,
            outgoing_direction);

    require(full.target_index == target_index,
            "full collision changed the target-selection order");
    require(full.target_velocity_cm_s == target_velocity_cm_s,
            "full collision changed the target-velocity draw order");
    require_vector_close(
        full.outgoing_dm_velocity_cm_s, outgoing_dm_velocity_cm_s,
        5.0e-13, 1.0e-6,
        "full collision changed the angle or azimuth draw order");
    require(full_rng == manual_rng,
            "full collision consumed an unexpected random variate");
}

void require_binomial_count(std::size_t count,
                            std::size_t sample_count,
                            double probability,
                            const std::string& message) {
    const double expected = sample_count * probability;
    const double standard_deviation = std::sqrt(
        sample_count * probability * (1.0 - probability));
    require(std::fabs(static_cast<double>(count) - expected) <=
                7.0 * standard_deviation + 12.0,
            message);
}

void test_full_collision_physics(const SolarBackground& background) {
    const SdProtonModel model{darkMatterMassGeV, 1.0e-34};
    const double radius_cm = 0.41 * solarRadiusCm;
    const Velocity3 incoming_dm_velocity_cm_s =
        multiply(8.5e7, normalized(Velocity3{{2.0, 1.0, -3.0}}));

    std::mt19937 first(20260927U);
    std::mt19937 replay(20260927U);
    for (std::size_t sample = 0; sample < 256; ++sample) {
        const CollisionSample first_sample = sample_sd_proton_collision(
            background, model, radius_cm,
            incoming_dm_velocity_cm_s, first);
        const CollisionSample replay_sample = sample_sd_proton_collision(
            background, model, radius_cm,
            incoming_dm_velocity_cm_s, replay);
        require(first_sample.target_index == replay_sample.target_index &&
                    first_sample.target_velocity_cm_s ==
                        replay_sample.target_velocity_cm_s &&
                    first_sample.outgoing_dm_velocity_cm_s ==
                        replay_sample.outgoing_dm_velocity_cm_s,
                "same-seed full collision did not replay");
    }

    constexpr std::size_t sampleCount = 12000;
    constexpr std::size_t angularBins = 12;
    std::array<std::size_t, angularBins> cosine_counts{{}};
    std::array<long double, 3> direction_sums{{}};
    std::array<std::array<long double, 3>, 3> direction_products{{}};
    std::mt19937 rng(20260928U);
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        const CollisionSample collision = sample_sd_proton_collision(
            background, model, radius_cm,
            incoming_dm_velocity_cm_s, rng);
        require(collision.target_index < background.target_count(),
                "full collision returned an invalid target index");
        const double target_mass_GeV =
            background.target(collision.target_index).mass_GeV;
        const Velocity3 relative_incoming = subtract(
            incoming_dm_velocity_cm_s, collision.target_velocity_cm_s);
        const double relative_speed_cm_s = norm(relative_incoming);
        require(std::isfinite(relative_speed_cm_s) &&
                    relative_speed_cm_s > 0.0,
                "full collision sampled an invalid relative speed");
        const Velocity3 center = center_of_mass_velocity(
            model.dark_matter_mass_GeV, target_mass_GeV,
            incoming_dm_velocity_cm_s, collision.target_velocity_cm_s);
        const Velocity3 outgoing_cm = subtract(
            collision.outgoing_dm_velocity_cm_s, center);
        const double expected_outgoing_cm_speed =
            target_mass_GeV /
            (model.dark_matter_mass_GeV + target_mass_GeV) *
            relative_speed_cm_s;
        require_close(norm(outgoing_cm), expected_outgoing_cm_speed,
                      4.0e-12, 2.0e-6,
                      "full collision changed the CM speed");
        const Velocity3 outgoing_direction = normalized(outgoing_cm);
        for (std::size_t first_component = 0;
             first_component < 3; ++first_component) {
            const long double first_value = static_cast<long double>(
                outgoing_direction[first_component]);
            direction_sums[first_component] += first_value;
            for (std::size_t second_component = 0;
                 second_component < 3; ++second_component) {
                direction_products[first_component][second_component] +=
                    first_value * static_cast<long double>(
                        outgoing_direction[second_component]);
            }
        }
        require_elastic_invariants(
            model.dark_matter_mass_GeV, target_mass_GeV,
            incoming_dm_velocity_cm_s, collision.target_velocity_cm_s,
            outgoing_direction, collision.outgoing_dm_velocity_cm_s,
            "full collision event");

        const double cosine = std::max(-1.0, std::min(
            1.0, dot(normalized(relative_incoming), outgoing_direction)));
        const std::size_t bin = std::min<std::size_t>(
            angularBins - 1,
            static_cast<std::size_t>(
                0.5 * (cosine + 1.0) * angularBins));
        ++cosine_counts[bin];
    }
    for (std::size_t count : cosine_counts) {
        require_binomial_count(
            count, sampleCount, 1.0 / angularBins,
            "outgoing CM direction is not isotropic");
    }
    const long double sample_count =
        static_cast<long double>(sampleCount);
    const long double mean_tolerance =
        7.0L * std::sqrt(1.0L / (3.0L * sample_count)) + 0.002L;
    const long double diagonal_tolerance =
        7.0L * std::sqrt(4.0L / (45.0L * sample_count)) + 0.002L;
    const long double cross_tolerance =
        7.0L * std::sqrt(1.0L / (15.0L * sample_count)) + 0.002L;
    for (std::size_t component = 0; component < 3; ++component) {
        const long double mean =
            direction_sums[component] / sample_count;
        const long double second_moment =
            direction_products[component][component] / sample_count;
        require(std::fabs(mean) <= mean_tolerance,
                "outgoing CM direction has a nonzero Cartesian mean");
        require(std::fabs(second_moment - 1.0L / 3.0L) <=
                    diagonal_tolerance,
                "outgoing CM direction has an anisotropic second moment");
        for (std::size_t other = component + 1; other < 3; ++other) {
            const long double cross_moment =
                direction_products[component][other] / sample_count;
            require(std::fabs(cross_moment) <= cross_tolerance,
                    "outgoing CM direction has a nonzero cross moment");
        }
    }
}

void test_errors_before_random_draws(const SolarBackground& background) {
    const SdProtonModel model{darkMatterMassGeV, 1.0e-34};
    const Velocity3 valid_velocity{{8.0e7, 0.0, 0.0}};
    const double nan = std::numeric_limits<double>::quiet_NaN();

    std::mt19937 zero_speed_rng(51U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_collision(
            background, model, 0.3 * solarRadiusCm,
            Velocity3{{0.0, 0.0, 0.0}}, zero_speed_rng);
    }, zero_speed_rng, "zero incoming collision speed");

    std::mt19937 nonfinite_rng(52U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_collision(
            background, model, 0.3 * solarRadiusCm,
            Velocity3{{nan, 0.0, 0.0}}, nonfinite_rng);
    }, nonfinite_rng, "nonfinite incoming collision velocity");

    std::mt19937 exterior_rng(53U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_collision(
            background, model, 1.01 * solarRadiusCm,
            valid_velocity, exterior_rng);
    }, exterior_rng, "exterior collision radius");

    std::mt19937 zero_rate_rng(54U);
    require_throws_without_consuming_rng<std::domain_error>([&] {
        sample_sd_proton_collision(
            background, SdProtonModel{darkMatterMassGeV, 0.0},
            0.3 * solarRadiusCm, valid_velocity, zero_rate_rng);
    }, zero_rate_rng, "zero-rate collision model");

    const double overflow_radius_cm = 0.3 * solarRadiusCm;
    const double overflow_temperature_K =
        background.temperature_K(overflow_radius_cm);
    const auto nominal_rates = direct_sd_proton_scattering_rates(
        background, model, overflow_radius_cm, norm(valid_velocity));
    double maximum_positive_rate_kappa = 0.0;
    for (const auto& rate : nominal_rates.target_rates) {
        if (rate.rate_s_inv > 0.0) {
            const double target_mass_GeV =
                background.target(rate.target_index).mass_GeV;
            maximum_positive_rate_kappa = std::max(
                maximum_positive_rate_kappa,
                std::sqrt(target_mass_GeV /
                          (2.0 * legacyBoltzmannGeVPerK *
                           overflow_temperature_K)));
        }
    }
    const std::size_t fallback_target_index =
        background.target_count() - 1;
    require(nominal_rates.target_rates[fallback_target_index].rate_s_inv == 0.0,
            "source-order fallback target unexpectedly has positive rate");
    const double fallback_kappa = std::sqrt(
        background.target(fallback_target_index).mass_GeV /
        (2.0 * legacyBoltzmannGeVPerK * overflow_temperature_K));
    require(fallback_kappa > maximum_positive_rate_kappa,
            "overflow-path precondition lacks a fallback-only domain window");
    const double maximum_safe_ratio =
        0.25 * std::sqrt(std::numeric_limits<double>::max());
    const double overflow_speed_natural = 0.5 * maximum_safe_ratio *
        (1.0 / fallback_kappa + 1.0 / maximum_positive_rate_kappa);
    const Velocity3 deterministic_overflow_velocity{{
        overflow_speed_natural * speedOfLightCmS, 0.0, 0.0}};
    require(maximum_positive_rate_kappa * overflow_speed_natural <
                maximum_safe_ratio &&
                fallback_kappa * overflow_speed_natural > maximum_safe_ratio,
            "overflow speed does not isolate the zero-rate fallback target");
    const auto overflow_rates = direct_sd_proton_scattering_rates(
        background, model, overflow_radius_cm,
        deterministic_overflow_velocity[0]);
    require(std::isfinite(overflow_rates.total_rate_s_inv) &&
                overflow_rates.total_rate_s_inv > 0.0,
            "overflow-path precondition did not reach target sampling");
    std::mt19937 overflow_rng(55U);
    require_throws_without_consuming_rng<std::overflow_error>([&] {
        sample_sd_proton_collision(
            background, model, overflow_radius_cm,
            deterministic_overflow_velocity, overflow_rng);
    }, overflow_rng, "deterministic target-sampler overflow");

    require_throws<std::invalid_argument>([&] {
        elastic_outgoing_dm_velocity_cm_s(
            0.0, protonMassGeV, valid_velocity,
            Velocity3{{0.0, 0.0, 0.0}},
            Velocity3{{1.0, 0.0, 0.0}});
    }, "nonpositive DM mass");
    require_throws<std::invalid_argument>([&] {
        elastic_outgoing_dm_velocity_cm_s(
            darkMatterMassGeV, protonMassGeV, valid_velocity,
            Velocity3{{0.0, 0.0, 0.0}},
            Velocity3{{2.0, 0.0, 0.0}});
    }, "non-unit outgoing direction");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "usage: test_single_collision SOLAR_MODEL NUCLEAR_DATA");
        }
        const SolarBackground background(argv[1], argv[2]);
        test_fixed_angle_kinematics();
        test_moving_target_covariance();
        test_full_collision_draw_order(background);
        test_full_collision_physics(background);
        test_errors_before_random_draws(background);
        std::cout << "single collision tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "single collision test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
