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
#include <vector>

namespace {

using transport::physics::SdProtonModel;
using transport::physics::SdScatteringRates;
using transport::physics::SdTargetScatteringRate;
using transport::physics::SolarBackground;
using transport::physics::direct_sd_proton_scattering_rates;
using transport::physics::sample_collision_conditioned_target_velocity_cm_s;
using transport::physics::sample_sd_proton_target_index;

using Velocity3 = std::array<double, 3>;

constexpr double boltzmannGeVPerK = 8.6173303e-14;
constexpr double speedOfLightCmS = 2.99792458e10;
constexpr double solarRadiusCm = 6.957e10;

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

double norm(const Velocity3& vector) {
    return std::hypot(std::hypot(vector[0], vector[1]), vector[2]);
}

double dot(const Velocity3& left, const Velocity3& right) {
    return left[0] * right[0] + left[1] * right[1] +
        left[2] * right[2];
}

Velocity3 subtract(const Velocity3& left, const Velocity3& right) {
    return Velocity3{{
        left[0] - right[0],
        left[1] - right[1],
        left[2] - right[2]
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
    return Velocity3{{
        vector[0] / magnitude,
        vector[1] / magnitude,
        vector[2] / magnitude
    }};
}

double thermal_speed_cm_s(double temperature_K, double target_mass_GeV) {
    return speedOfLightCmS * std::sqrt(
        2.0 * boltzmannGeVPerK * temperature_K / target_mass_GeV);
}

SdScatteringRates synthetic_rates(double scale = 1.0) {
    SdScatteringRates rates;
    rates.target_rates.reserve(63);
    for (std::size_t index = 0; index < 63; ++index) {
        double rate = 0.0;
        if (index == 1) {
            rate = 1.0 * scale;
        } else if (index == 2) {
            rate = 2.0 * scale;
        } else if (index == 4) {
            rate = 3.0 * scale;
        }
        rates.target_rates.push_back(
            SdTargetScatteringRate{index, 0.0, 0.0, rate});
    }
    rates.total_rate_s_inv = 6.0 * scale;
    return rates;
}

void require_binomial_count(std::size_t count,
                            std::size_t sample_count,
                            double probability,
                            const std::string& message) {
    const double expected = sample_count * probability;
    const double standard_deviation = std::sqrt(
        sample_count * probability * (1.0 - probability));
    require(std::fabs(static_cast<double>(count) - expected) <=
                7.0 * standard_deviation + 8.0,
            message);
}

void test_target_selection_contract(const SolarBackground& background) {
    const SdScatteringRates rates = synthetic_rates();

    std::mt19937 consumption_rng(20260915U);
    std::mt19937 one_uniform_rng = consumption_rng;
    sample_sd_proton_target_index(rates, consumption_rng);
    std::uniform_real_distribution<double> unit_uniform(0.0, 1.0);
    unit_uniform(one_uniform_rng);
    require(consumption_rng == one_uniform_rng,
            "target selection must consume exactly one uniform draw");

    std::mt19937 first(20260916U);
    std::mt19937 replay(20260916U);
    std::mt19937 scaled_rng(20260916U);
    const SdScatteringRates scaled = synthetic_rates(8.0);
    for (std::size_t draw = 0; draw < 4096; ++draw) {
        const std::size_t first_index =
            sample_sd_proton_target_index(rates, first);
        require(first_index == sample_sd_proton_target_index(rates, replay),
                "same-seed target selection did not replay");
        require(first_index ==
                    sample_sd_proton_target_index(scaled, scaled_rng),
                "common rate scaling changed the target sequence");
    }

    constexpr std::size_t sampleCount = 180000;
    std::array<std::size_t, 63> counts{{}};
    std::mt19937 statistical_rng(20260917U);
    for (std::size_t draw = 0; draw < sampleCount; ++draw) {
        ++counts[sample_sd_proton_target_index(rates, statistical_rng)];
    }
    require_binomial_count(counts[1], sampleCount, 1.0 / 6.0,
                           "weight-1 target frequency is inconsistent");
    require_binomial_count(counts[2], sampleCount, 2.0 / 6.0,
                           "weight-2 target frequency is inconsistent");
    require_binomial_count(counts[4], sampleCount, 3.0 / 6.0,
                           "weight-3 target frequency is inconsistent");
    for (std::size_t index = 0; index < counts.size(); ++index) {
        if (index != 1 && index != 2 && index != 4) {
            require(counts[index] == 0,
                    "a zero-rate synthetic target was selected");
        }
    }

    const SdScatteringRates solar_rates = direct_sd_proton_scattering_rates(
        background, SdProtonModel{0.1, 1.0e-34},
        0.37 * solarRadiusCm, 8.0e7);
    std::mt19937 solar_rng(20260918U);
    for (std::size_t draw = 0; draw < 10000; ++draw) {
        const std::size_t index =
            sample_sd_proton_target_index(solar_rates, solar_rng);
        require(index < solar_rates.target_rates.size() &&
                    solar_rates.target_rates[index].target_index == index,
                "solar target selection returned an invalid index");
        require(solar_rates.target_rates[index].rate_s_inv > 0.0 ||
                    index + 1 == solar_rates.target_rates.size(),
                "solar target selection returned a zero-rate non-fallback target");
    }
}

void test_target_selection_errors() {
    SdScatteringRates zero = synthetic_rates(0.0);
    std::mt19937 zero_rng(11U);
    require_throws_without_consuming_rng<std::domain_error>([&] {
        sample_sd_proton_target_index(zero, zero_rng);
    }, zero_rng, "zero total rate");

    SdScatteringRates truncated = synthetic_rates();
    truncated.target_rates.pop_back();
    std::mt19937 truncated_rng(12U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_target_index(truncated, truncated_rng);
    }, truncated_rng, "truncated target breakdown");

    SdScatteringRates wrong_order = synthetic_rates();
    wrong_order.target_rates[3].target_index = 4;
    std::mt19937 wrong_order_rng(13U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_target_index(wrong_order, wrong_order_rng);
    }, wrong_order_rng, "out-of-order target breakdown");

    SdScatteringRates negative = synthetic_rates();
    negative.target_rates[0].rate_s_inv = -1.0;
    std::mt19937 negative_rng(14U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_target_index(negative, negative_rng);
    }, negative_rng, "negative target rate");

    SdScatteringRates nonfinite = synthetic_rates();
    nonfinite.target_rates[0].rate_s_inv =
        std::numeric_limits<double>::infinity();
    std::mt19937 nonfinite_rng(15U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_target_index(nonfinite, nonfinite_rng);
    }, nonfinite_rng, "nonfinite target rate");

    SdScatteringRates wrong_total = synthetic_rates();
    wrong_total.total_rate_s_inv = 7.0;
    std::mt19937 wrong_total_rng(16U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_sd_proton_target_index(wrong_total, wrong_total_rng);
    }, wrong_total_rng, "inconsistent target-rate total");
}

template <typename Integrand>
long double composite_simpson(Integrand integrand,
                              long double lower,
                              long double upper) {
    constexpr int intervals = 20000;
    if (upper <= lower) {
        return 0.0L;
    }
    const long double step = (upper - lower) / intervals;
    long double sum = integrand(lower) + integrand(upper);
    for (int index = 1; index < intervals; ++index) {
        sum += (index % 2 == 0 ? 2.0L : 4.0L) *
            integrand(lower + index * step);
    }
    return sum * step / 3.0L;
}

long double speed_weight(long double target_speed_ratio,
                         long double dm_speed_ratio) {
    if (target_speed_ratio == 0.0L) {
        return 0.0L;
    }
    const long double sum = target_speed_ratio + dm_speed_ratio;
    const long double difference =
        std::fabs(target_speed_ratio - dm_speed_ratio);
    const long double angular_integral =
        (sum * sum * sum - difference * difference * difference) /
        (3.0L * target_speed_ratio * dm_speed_ratio);
    return target_speed_ratio * target_speed_ratio *
        std::exp(-target_speed_ratio * target_speed_ratio) *
        angular_integral;
}

double expected_target_speed_cdf(double threshold, double dm_speed_ratio) {
    const auto integrand = [=](long double target_speed_ratio) {
        return speed_weight(target_speed_ratio, dm_speed_ratio);
    };
    const long double normalization =
        composite_simpson(integrand, 0.0L, 10.0L);
    const long double partial = composite_simpson(
        integrand, 0.0L, std::min(10.0L,
                                  static_cast<long double>(threshold)));
    return static_cast<double>(partial / normalization);
}

long double relative_speed_weight(long double relative_speed_ratio,
                                  long double dm_speed_ratio) {
    const long double lower_exponent =
        -(dm_speed_ratio - relative_speed_ratio) *
        (dm_speed_ratio - relative_speed_ratio);
    const long double upper_exponent =
        -(dm_speed_ratio + relative_speed_ratio) *
        (dm_speed_ratio + relative_speed_ratio);
    return relative_speed_ratio * relative_speed_ratio /
        (2.0L * dm_speed_ratio) *
        (std::exp(lower_exponent) - std::exp(upper_exponent));
}

double expected_relative_speed_cdf(double threshold,
                                   double dm_speed_ratio) {
    const auto integrand = [=](long double relative_speed_ratio) {
        return relative_speed_weight(relative_speed_ratio, dm_speed_ratio);
    };
    const long double upper = dm_speed_ratio + 10.0L;
    const long double normalization =
        composite_simpson(integrand, 0.0L, upper);
    const long double partial = composite_simpson(
        integrand, 0.0L, std::min(upper,
                                  static_cast<long double>(threshold)));
    return static_cast<double>(partial / normalization);
}

double conditional_angle_cdf(double cosine,
                             double target_speed_ratio,
                             double dm_speed_ratio) {
    const double sum = target_speed_ratio + dm_speed_ratio;
    const double difference =
        std::fabs(target_speed_ratio - dm_speed_ratio);
    const double relative_squared = std::max(
        0.0, target_speed_ratio * target_speed_ratio +
            dm_speed_ratio * dm_speed_ratio -
            2.0 * target_speed_ratio * dm_speed_ratio * cosine);
    const double numerator = sum * sum * sum -
        std::pow(relative_squared, 1.5);
    const double denominator = sum * sum * sum -
        difference * difference * difference;
    return std::max(0.0, std::min(1.0, numerator / denominator));
}

void require_empirical_cdf(std::size_t count,
                           std::size_t sample_count,
                           double expected,
                           const std::string& message) {
    const double empirical =
        static_cast<double>(count) / static_cast<double>(sample_count);
    const double standard_error = std::sqrt(
        std::max(0.0, expected * (1.0 - expected) / sample_count));
    require(std::fabs(empirical - expected) <=
                7.0 * standard_error + 3.0e-3,
            message);
}

void validate_conditioned_state(double temperature_K,
                                double target_mass_GeV,
                                double dm_speed_ratio,
                                unsigned int seed) {
    constexpr std::size_t sampleCount = 80000;
    const double thermal_speed =
        thermal_speed_cm_s(temperature_K, target_mass_GeV);
    const Velocity3 axis = normalized(Velocity3{{1.0, 2.0, 3.0}});
    const Velocity3 dm_velocity{{
        dm_speed_ratio * thermal_speed * axis[0],
        dm_speed_ratio * thermal_speed * axis[1],
        dm_speed_ratio * thermal_speed * axis[2]
    }};
    const Velocity3 reference{{0.0, 0.0, 1.0}};
    const Velocity3 transverse_1 = normalized(cross(reference, axis));
    const Velocity3 transverse_2 = normalized(cross(axis, transverse_1));

    const std::vector<double> speed_thresholds{0.5, 1.0, 1.5, 2.0, 3.0};
    std::vector<double> relative_thresholds{
        std::max(0.1, dm_speed_ratio - 2.0),
        std::max(0.2, dm_speed_ratio - 1.0),
        dm_speed_ratio,
        dm_speed_ratio + 1.0,
        dm_speed_ratio + 2.0
    };
    std::sort(relative_thresholds.begin(), relative_thresholds.end());
    relative_thresholds.erase(
        std::unique(relative_thresholds.begin(), relative_thresholds.end()),
        relative_thresholds.end());

    std::vector<std::size_t> speed_counts(speed_thresholds.size(), 0);
    std::vector<std::size_t> relative_counts(relative_thresholds.size(), 0);
    std::array<std::size_t, 10> angle_pit_counts{{}};
    double relative_sum = 0.0;
    double relative_square_sum = 0.0;
    double transverse_1_sum = 0.0;
    double transverse_2_sum = 0.0;
    double transverse_1_square_sum = 0.0;
    double transverse_2_square_sum = 0.0;

    std::mt19937 rng(seed);
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        const Velocity3 target_velocity =
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, target_mass_GeV, dm_velocity, rng);
        require(std::isfinite(norm(target_velocity)),
                "target velocity sample is nonfinite");

        const double target_ratio = norm(target_velocity) / thermal_speed;
        const double relative_ratio =
            norm(subtract(dm_velocity, target_velocity)) / thermal_speed;
        const double cosine = dot(target_velocity, axis) /
            norm(target_velocity);
        const double angle_pit = conditional_angle_cdf(
            cosine, target_ratio, dm_speed_ratio);
        const std::size_t pit_bin = std::min<std::size_t>(
            9, static_cast<std::size_t>(10.0 * angle_pit));
        ++angle_pit_counts[pit_bin];

        for (std::size_t index = 0; index < speed_thresholds.size(); ++index) {
            if (target_ratio <= speed_thresholds[index]) {
                ++speed_counts[index];
            }
        }
        for (std::size_t index = 0;
             index < relative_thresholds.size(); ++index) {
            if (relative_ratio <= relative_thresholds[index]) {
                ++relative_counts[index];
            }
        }

        relative_sum += relative_ratio;
        relative_square_sum += relative_ratio * relative_ratio;
        const double transverse_component_1 =
            dot(target_velocity, transverse_1) / thermal_speed;
        const double transverse_component_2 =
            dot(target_velocity, transverse_2) / thermal_speed;
        transverse_1_sum += transverse_component_1;
        transverse_2_sum += transverse_component_2;
        transverse_1_square_sum +=
            transverse_component_1 * transverse_component_1;
        transverse_2_square_sum +=
            transverse_component_2 * transverse_component_2;
    }

    for (std::size_t index = 0; index < speed_thresholds.size(); ++index) {
        require_empirical_cdf(
            speed_counts[index], sampleCount,
            expected_target_speed_cdf(
                speed_thresholds[index], dm_speed_ratio),
            "conditioned target-speed CDF is inconsistent");
    }
    for (std::size_t index = 0;
         index < relative_thresholds.size(); ++index) {
        require_empirical_cdf(
            relative_counts[index], sampleCount,
            expected_relative_speed_cdf(
                relative_thresholds[index], dm_speed_ratio),
            "conditioned relative-speed CDF is inconsistent");
    }
    for (std::size_t count : angle_pit_counts) {
        require_binomial_count(count, sampleCount, 0.1,
                               "conditional-angle PIT is not uniform");
    }

    const double inverse_sqrt_pi = 1.0 / std::sqrt(std::acos(-1.0));
    const double mean_relative_bath =
        std::exp(-dm_speed_ratio * dm_speed_ratio) * inverse_sqrt_pi +
        (dm_speed_ratio + 0.5 / dm_speed_ratio) *
            std::erf(dm_speed_ratio);
    const double expected_conditioned_relative_mean =
        (dm_speed_ratio * dm_speed_ratio + 1.5) /
        mean_relative_bath;
    const double relative_mean = relative_sum / sampleCount;
    const double relative_variance = std::max(
        0.0, relative_square_sum / sampleCount -
            relative_mean * relative_mean);
    require(std::fabs(relative_mean - expected_conditioned_relative_mean) <=
                7.0 * std::sqrt(relative_variance / sampleCount) + 5.0e-4,
            "conditioned mean relative speed violates the analytic identity");

    const double transverse_mean_1 = transverse_1_sum / sampleCount;
    const double transverse_mean_2 = transverse_2_sum / sampleCount;
    const double transverse_variance_1 = std::max(
        0.0, transverse_1_square_sum / sampleCount -
            transverse_mean_1 * transverse_mean_1);
    const double transverse_variance_2 = std::max(
        0.0, transverse_2_square_sum / sampleCount -
            transverse_mean_2 * transverse_mean_2);
    require(std::fabs(transverse_mean_1) <=
                7.0 * std::sqrt(transverse_variance_1 / sampleCount) + 1.0e-3 &&
            std::fabs(transverse_mean_2) <=
                7.0 * std::sqrt(transverse_variance_2 / sampleCount) + 1.0e-3,
            "target sampler has a transverse directional bias");
    require(std::fabs(transverse_variance_1 - transverse_variance_2) <=
                0.05 * (transverse_variance_1 + transverse_variance_2),
            "target sampler breaks azimuthal symmetry");
}

void test_conditioned_velocity_distribution(
    const SolarBackground& background) {
    const double temperature_K =
        background.temperature_K(0.3 * solarRadiusCm);
    const std::array<std::size_t, 4> target_indices{{0, 2, 13, 0}};
    const std::array<double, 4> speed_ratios{{0.1, 1.0, 3.0, 10.0}};
    for (std::size_t state = 0; state < speed_ratios.size(); ++state) {
        validate_conditioned_state(
            temperature_K,
            background.target(target_indices[state]).mass_GeV,
            speed_ratios[state],
            static_cast<unsigned int>(20260920U + state));
    }
}

void test_zero_speed_limit_and_reproducibility(
    const SolarBackground& background) {
    const double temperature_K =
        background.temperature_K(0.5 * solarRadiusCm);
    const double target_mass_GeV = background.target(0).mass_GeV;
    const double thermal_speed =
        thermal_speed_cm_s(temperature_K, target_mass_GeV);

    std::mt19937 invalid_rng(31U);
    require_throws_without_consuming_rng<std::invalid_argument>([&] {
        sample_collision_conditioned_target_velocity_cm_s(
            temperature_K, target_mass_GeV,
            Velocity3{{0.0, 0.0, 0.0}}, invalid_rng);
    }, invalid_rng, "zero-speed legacy sampler");

    constexpr std::size_t sampleCount = 80000;
    const double one_sided_ratio = 1.0e-6;
    const Velocity3 one_sided_velocity{{
        one_sided_ratio * thermal_speed, 0.0, 0.0}};
    double speed_sum = 0.0;
    double speed_square_sum = 0.0;
    double cosine_sum = 0.0;
    std::mt19937 limit_rng(20260924U);
    for (std::size_t sample = 0; sample < sampleCount; ++sample) {
        const Velocity3 target_velocity =
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, target_mass_GeV,
                one_sided_velocity, limit_rng);
        const double target_ratio = norm(target_velocity) / thermal_speed;
        speed_sum += target_ratio;
        speed_square_sum += target_ratio * target_ratio;
        cosine_sum += target_velocity[0] / norm(target_velocity);
    }
    const double mean_speed = speed_sum / sampleCount;
    const double speed_variance = std::max(
        0.0, speed_square_sum / sampleCount - mean_speed * mean_speed);
    const double expected_limit =
        3.0 * std::sqrt(std::acos(-1.0)) / 4.0;
    require(std::fabs(mean_speed - expected_limit) <=
                7.0 * std::sqrt(speed_variance / sampleCount) + 5.0e-4,
            "one-sided target-speed limit is inconsistent");
    require(std::fabs(cosine_sum / sampleCount) < 1.0e-2,
            "one-sided target direction is not approaching isotropy");

    const Velocity3 positive_velocity{{thermal_speed, 0.0, 0.0}};
    std::mt19937 first(20260925U);
    std::mt19937 replay(20260925U);
    for (std::size_t sample = 0; sample < 256; ++sample) {
        const Velocity3 first_sample =
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, target_mass_GeV,
                positive_velocity, first);
        const Velocity3 replay_sample =
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, target_mass_GeV,
                positive_velocity, replay);
        require(first_sample == replay_sample,
                "same-seed target-velocity sampling did not replay");
    }
}

void test_velocity_errors(const SolarBackground& background) {
    const double temperature_K =
        background.temperature_K(0.5 * solarRadiusCm);
    const double target_mass_GeV = background.target(0).mass_GeV;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();

    for (double temperature : {0.0, -1.0, nan, infinity}) {
        std::mt19937 rng(41U);
        require_throws_without_consuming_rng<std::invalid_argument>([&] {
            sample_collision_conditioned_target_velocity_cm_s(
                temperature, target_mass_GeV,
                Velocity3{{1.0, 0.0, 0.0}}, rng);
        }, rng, "invalid target temperature");
    }
    for (double mass : {0.0, -1.0, nan, infinity}) {
        std::mt19937 rng(42U);
        require_throws_without_consuming_rng<std::invalid_argument>([&] {
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, mass,
                Velocity3{{1.0, 0.0, 0.0}}, rng);
        }, rng, "invalid target mass");
    }
    for (double invalid_component : {nan, infinity}) {
        std::mt19937 rng(43U);
        require_throws_without_consuming_rng<std::invalid_argument>([&] {
            sample_collision_conditioned_target_velocity_cm_s(
                temperature_K, target_mass_GeV,
                Velocity3{{invalid_component, 0.0, 0.0}}, rng);
        }, rng, "invalid incoming velocity");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            throw std::invalid_argument(
                "usage: test_target_sampling SOLAR_MODEL NUCLEAR_DATA");
        }
        const SolarBackground background(argv[1], argv[2]);
        test_target_selection_contract(background);
        test_target_selection_errors();
        test_conditioned_velocity_distribution(background);
        test_zero_speed_limit_and_reproducibility(background);
        test_velocity_errors(background);
        std::cout << "target sampling tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "target sampling test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
