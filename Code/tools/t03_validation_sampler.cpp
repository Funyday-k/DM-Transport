#include "transport/physics/ScatteringPhysics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Velocity = transport::physics::CartesianVelocityCmS;
using Background = transport::physics::SolarBackground;
using Model = transport::physics::SdProtonModel;

constexpr double speedOfLightCmS = 2.99792458e10;
constexpr double boltzmannGeVPerK = 8.6173303e-14;
constexpr double darkMatterMassGeV = 0.1;
constexpr double protonCrossSectionCm2 = 1.0e-34;
constexpr std::uint32_t collisionSeedMask = 0x9e3779b9U;
constexpr std::size_t weakBinCount = 6;
constexpr std::size_t vectorMomentCount = 9;

const std::array<double, weakBinCount> weakSpeedBinLowerEdges{{
    0.0,
    0.7089337109534111,
    0.966738952196692,
    1.2137063222011224,
    1.5234217531739933,
    1.9767053274642596
}};
const std::array<double, 11> outgoingSpeedHistogramEdges{{
    0.0, 0.5, 0.8, 0.95, 0.99, 1.0, 1.01, 1.05, 1.2, 1.5, 2.0
}};
const std::array<double, 9> outgoingCosineHistogramEdges{{
    -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0
}};
const std::array<double, 11> outgoingCmCosineHistogramEdges{{
    -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
    0.2, 0.4, 0.6, 0.8, 1.0
}};
const std::array<double, 9> canonicalComponentCutpoints{{
    -2.0, -1.0, -0.5, -0.25, 0.0,
    0.25, 0.5, 1.0, 2.0
}};
const std::array<const char*, vectorMomentCount> vectorMomentNames{{
    "x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz"
}};

struct Moment {
    long double sum = 0.0L;
    long double sum_sq = 0.0L;

    void add(long double value) {
        if (!std::isfinite(value)) {
            throw std::overflow_error("nonfinite validation observable");
        }
        sum += value;
        sum_sq += value * value;
        if (!std::isfinite(sum) || !std::isfinite(sum_sq)) {
            throw std::overflow_error("validation moment is not representable");
        }
    }
};

double parse_double(const char* text, const std::string& label) {
    const std::string value(text);
    std::size_t parsed = 0;
    double result = 0.0;
    try {
        result = std::stod(value, &parsed);
    } catch (const std::exception&) {
        throw std::invalid_argument(label + " must be a finite number");
    }
    if (parsed != value.size() || !std::isfinite(result)) {
        throw std::invalid_argument(label + " must be a finite number");
    }
    return result;
}

std::uint64_t parse_uint64(const char* text, const std::string& label) {
    const std::string value(text);
    if (value.empty() || value.front() == '-') {
        throw std::invalid_argument(label + " must be an unsigned integer");
    }
    std::size_t parsed = 0;
    std::uint64_t result = 0;
    try {
        result = std::stoull(value, &parsed);
    } catch (const std::exception&) {
        throw std::invalid_argument(label + " must be an unsigned integer");
    }
    if (parsed != value.size()) {
        throw std::invalid_argument(label + " must be an unsigned integer");
    }
    return result;
}

double dot(const Velocity& first, const Velocity& second) {
    return first[0] * second[0] + first[1] * second[1] +
        first[2] * second[2];
}

double norm(const Velocity& velocity) {
    return std::sqrt(dot(velocity, velocity));
}

Velocity subtract(const Velocity& first, const Velocity& second) {
    return Velocity{{
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2]
    }};
}

Velocity parse_direction(const char* text) {
    const std::string value(text);
    const std::size_t first_comma = value.find(',');
    const std::size_t second_comma = value.find(',', first_comma + 1);
    if (first_comma == std::string::npos ||
        second_comma == std::string::npos ||
        value.find(',', second_comma + 1) != std::string::npos) {
        throw std::invalid_argument("direction must be x,y,z");
    }
    const std::string x = value.substr(0, first_comma);
    const std::string y = value.substr(first_comma + 1,
                                       second_comma - first_comma - 1);
    const std::string z = value.substr(second_comma + 1);
    Velocity direction{{
        parse_double(x.c_str(), "direction x"),
        parse_double(y.c_str(), "direction y"),
        parse_double(z.c_str(), "direction z")
    }};
    const double length = norm(direction);
    if (!std::isfinite(length) || length <= 0.0) {
        throw std::invalid_argument("direction must have a finite positive norm");
    }
    for (double& component : direction) {
        component /= length;
    }
    return direction;
}

std::array<Velocity, 3> shortest_rotation_axes(const Velocity& direction) {
    // The columns of the shortest proper rotation R taking z to direction.
    // Canonical components are dot(R e_i, sampled_lab_vector), i.e. R^-1.
    if (1.0 + direction[2] < 1.0e-12) {
        // The shortest rotation is nonunique at -z; choose a pi turn about x.
        return std::array<Velocity, 3>{{
            Velocity{{1.0, 0.0, 0.0}},
            Velocity{{0.0, -1.0, 0.0}},
            Velocity{{0.0, 0.0, -1.0}}
        }};
    }
    const double x = direction[0];
    const double y = direction[1];
    const double denominator = 1.0 + direction[2];
    return std::array<Velocity, 3>{{
        Velocity{{1.0 - x * x / denominator,
                  -x * y / denominator, -x}},
        Velocity{{-x * y / denominator,
                  1.0 - y * y / denominator, -y}},
        direction
    }};
}

template <typename Container>
void print_array(const Container& values) {
    std::cout << '[';
    bool first = true;
    for (const auto& value : values) {
        if (!first) {
            std::cout << ',';
        }
        std::cout << value;
        first = false;
    }
    std::cout << ']';
}

void print_moment(const Moment& moment) {
    std::cout << "{\"sum\":" << moment.sum
              << ",\"sum_sq\":" << moment.sum_sq << '}';
}

template <std::size_t N>
void print_histogram3(
    const std::array<std::array<std::uint64_t, N>, 3>& counts) {
    std::cout << '[';
    for (std::size_t component = 0; component < counts.size(); ++component) {
        if (component != 0) {
            std::cout << ',';
        }
        print_array(counts[component]);
    }
    std::cout << ']';
}

template <std::size_t N>
std::size_t histogram_bin(double value,
                          const std::array<double, N>& lower_edges) {
    if (!std::isfinite(value) || value < lower_edges.front()) {
        throw std::runtime_error("histogram observable out of range");
    }
    const auto boundary = std::upper_bound(
        lower_edges.begin(), lower_edges.end(), value);
    return static_cast<std::size_t>(boundary - lower_edges.begin() - 1);
}

double energy_change_eV(const Velocity& incoming,
                        const Velocity& outgoing) {
    const long double speed_squared_difference =
        static_cast<long double>(dot(outgoing, outgoing)) -
        static_cast<long double>(dot(incoming, incoming));
    const long double result = 0.5L * darkMatterMassGeV * 1.0e9L *
        speed_squared_difference /
        (static_cast<long double>(speedOfLightCmS) * speedOfLightCmS);
    if (!std::isfinite(result)) {
        throw std::overflow_error("energy change is not representable");
    }
    return static_cast<double>(result);
}

void print_common(const char* mode,
                  std::uint64_t sample_count,
                  double radius_Rsun,
                  double radius_cm,
                  double temperature_K,
                  std::uint32_t incoming_seed,
                  std::uint32_t collision_seed,
                  const std::vector<std::uint64_t>& target_counts) {
    std::cout << "{\"schema_version\":1"
              << ",\"producer\":\"DM-Transport/t03_validation_sampler\""
              << ",\"qualification\":\"project_physics_validation_not_independent_oracle\""
              << ",\"mode\":\"" << mode << '\"'
              << ",\"sample_count\":" << sample_count
              << ",\"radius_Rsun\":" << radius_Rsun
              << ",\"radius_cm\":" << radius_cm
              << ",\"bath_temperature_K\":" << temperature_K
              << ",\"dark_matter_mass_GeV\":" << darkMatterMassGeV
              << ",\"proton_cross_section_cm2\":" << protonCrossSectionCm2
              << ",\"seed_incoming\":" << incoming_seed
              << ",\"seed_collision\":" << collision_seed
              << ",\"target_counts\":";
    print_array(target_counts);
}

void run_thermal(const Background& background,
                 double radius_Rsun,
                 double temperature_ratio,
                 std::uint64_t sample_count,
                 std::uint32_t seed) {
    if (!std::isfinite(temperature_ratio) || temperature_ratio <= 0.0) {
        throw std::invalid_argument("Tchi/T must be finite and positive");
    }
    const double radius_cm = radius_Rsun * background.solar_radius_cm();
    const double bath_temperature_K = background.temperature_K(radius_cm);
    const double bath_dm_thermal_speed_cm_s = speedOfLightCmS *
        std::sqrt(2.0 * boltzmannGeVPerK * bath_temperature_K /
                  darkMatterMassGeV);
    const double component_sigma_cm_s = bath_dm_thermal_speed_cm_s *
        std::sqrt(temperature_ratio / 2.0);
    const double kbt_eV = boltzmannGeVPerK * bath_temperature_K * 1.0e9;
    if (!std::isfinite(bath_dm_thermal_speed_cm_s) ||
        bath_dm_thermal_speed_cm_s <= 0.0 ||
        !std::isfinite(component_sigma_cm_s) ||
        component_sigma_cm_s <= 0.0 ||
        !std::isfinite(kbt_eV) || kbt_eV <= 0.0) {
        throw std::overflow_error("thermal validation scale is invalid");
    }
    const Model model{darkMatterMassGeV, protonCrossSectionCm2};
    long double analytic_mb_mean_total_rate_s_inv = 0.0L;
    for (std::size_t target_index = 0;
         target_index < background.target_count(); ++target_index) {
        const auto& target = background.target(target_index);
        const double target_cross_section_cm2 =
            transport::physics::sd_proton_nucleus_cross_section_cm2(
                model, target);
        const double target_number_density_cm3 =
            background.number_density_cm3(target_index, radius_cm);
        const double relative_mean_speed_cm_s = speedOfLightCmS *
            std::sqrt(8.0 / std::acos(-1.0) * boltzmannGeVPerK *
                      (bath_temperature_K * temperature_ratio /
                           darkMatterMassGeV +
                       bath_temperature_K / target.mass_GeV));
        analytic_mb_mean_total_rate_s_inv +=
            static_cast<long double>(target_number_density_cm3) *
            target_cross_section_cm2 * relative_mean_speed_cm_s;
    }
    if (!std::isfinite(analytic_mb_mean_total_rate_s_inv) ||
        analytic_mb_mean_total_rate_s_inv <= 0.0L) {
        throw std::runtime_error("analytic thermal mean rate is invalid");
    }
    const double gamma_ref_s_inv =
        transport::physics::direct_sd_proton_scattering_rates(
            background, model, radius_cm, bath_dm_thermal_speed_cm_s)
            .total_rate_s_inv;
    if (!std::isfinite(gamma_ref_s_inv) || gamma_ref_s_inv <= 0.0) {
        throw std::runtime_error("thermal reference rate must be positive");
    }

    const std::uint32_t collision_seed = seed ^ collisionSeedMask;
    std::mt19937 incoming_rng(seed);
    std::mt19937 collision_rng(collision_seed);
    std::normal_distribution<double> component_normal(0.0,
                                                       component_sigma_cm_s);
    std::array<Moment, weakBinCount> weak_speed_bin_moments;
    Moment rate_ratio;
    Moment weak_energy;
    Moment rate_weighted_delta_energy;
    std::vector<std::uint64_t> target_counts(background.target_count(), 0);

    for (std::uint64_t sample = 0; sample < sample_count; ++sample) {
        const Velocity incoming{{
            component_normal(incoming_rng),
            component_normal(incoming_rng),
            component_normal(incoming_rng)
        }};
        const double incoming_speed = norm(incoming);
        if (!std::isfinite(incoming_speed) || incoming_speed <= 0.0) {
            throw std::runtime_error("sampled incoming speed is invalid");
        }
        const double gamma_s_inv =
            transport::physics::direct_sd_proton_scattering_rates(
                background, model, radius_cm, incoming_speed)
                .total_rate_s_inv;
        const auto collision = transport::physics::sample_sd_proton_collision(
            background, model, radius_cm, incoming, collision_rng);
        if (collision.target_index >= target_counts.size()) {
            throw std::runtime_error("sampled target index is out of range");
        }
        ++target_counts[collision.target_index];
        const double outgoing_speed = norm(collision.outgoing_dm_velocity_cm_s);
        const double normalized_rate = gamma_s_inv / gamma_ref_s_inv;
        rate_ratio.add(normalized_rate);
        const double delta_energy_eV = energy_change_eV(
            incoming, collision.outgoing_dm_velocity_cm_s);
        const std::size_t incoming_bin = histogram_bin(
            incoming_speed / bath_dm_thermal_speed_cm_s,
            weakSpeedBinLowerEdges);
        const std::size_t outgoing_bin = histogram_bin(
            outgoing_speed / bath_dm_thermal_speed_cm_s,
            weakSpeedBinLowerEdges);
        for (std::size_t bin = 0; bin < weakBinCount; ++bin) {
            const int indicator_difference =
                static_cast<int>(outgoing_bin == bin) -
                static_cast<int>(incoming_bin == bin);
            weak_speed_bin_moments[bin].add(
                normalized_rate * indicator_difference);
        }
        weak_energy.add(normalized_rate * delta_energy_eV / kbt_eV);
        rate_weighted_delta_energy.add(gamma_s_inv * delta_energy_eV);
    }

    print_common("thermal", sample_count, radius_Rsun, radius_cm,
                 bath_temperature_K, seed, collision_seed, target_counts);
    std::cout << ",\"tchi_over_t\":" << temperature_ratio
              << ",\"bath_dm_thermal_speed_cm_s\":"
              << bath_dm_thermal_speed_cm_s
              << ",\"gamma_ref_s_inv\":" << gamma_ref_s_inv
              << ",\"analytic_mb_mean_total_rate_s_inv\":"
              << analytic_mb_mean_total_rate_s_inv
              << ",\"kbt_eV\":" << kbt_eV
              << ",\"weak_observable_definition\":\"(Gamma(v)/Gamma_ref)*(phi(v_out)-phi(v_in)); speed bins use indicators; energy phi=E/(kBT_bath)\""
              << ",\"weak_speed_bin_lower_edges_over_bath_dm_thermal_speed\":";
    print_array(weakSpeedBinLowerEdges);
    std::cout << ",\"weak_speed_bin_last_open_ended\":true"
              << ",\"weak_speed_bin_moments\":[";
    for (std::size_t bin = 0; bin < weakBinCount; ++bin) {
        if (bin != 0) {
            std::cout << ',';
        }
        print_moment(weak_speed_bin_moments[bin]);
    }
    std::cout << "]"
              << ",\"rate_ratio\":";
    print_moment(rate_ratio);
    std::cout
              << ",\"weak_energy\":";
    print_moment(weak_energy);
    std::cout << ",\"rate_weighted_delta_energy_eV_per_s\":";
    print_moment(rate_weighted_delta_energy);
    std::cout << "}\n";
}

void run_fixed(const Background& background,
               double radius_Rsun,
               double fixed_speed_over_escape,
               std::uint64_t sample_count,
               std::uint32_t seed,
               const Velocity& direction) {
    if (!std::isfinite(fixed_speed_over_escape) ||
        fixed_speed_over_escape <= 0.0) {
        throw std::invalid_argument("v/vesc must be finite and positive");
    }
    const double radius_cm = radius_Rsun * background.solar_radius_cm();
    const double bath_temperature_K = background.temperature_K(radius_cm);
    const double escape_speed_cm_s = background.escape_speed_cm_s(radius_cm);
    const double canonical_component_scale_cm_s = speedOfLightCmS *
        std::sqrt(2.0 * boltzmannGeVPerK * bath_temperature_K /
                  darkMatterMassGeV);
    const double fixed_speed_cm_s = fixed_speed_over_escape * escape_speed_cm_s;
    if (!std::isfinite(fixed_speed_cm_s) || fixed_speed_cm_s <= 0.0 ||
        !std::isfinite(canonical_component_scale_cm_s) ||
        canonical_component_scale_cm_s <= 0.0) {
        throw std::overflow_error("fixed speed is not representable");
    }
    const std::array<Velocity, 3> canonical_axes =
        shortest_rotation_axes(direction);
    const Velocity incoming{{
        fixed_speed_cm_s * direction[0],
        fixed_speed_cm_s * direction[1],
        fixed_speed_cm_s * direction[2]
    }};
    const Model model{darkMatterMassGeV, protonCrossSectionCm2};
    const double total_rate_s_inv =
        transport::physics::direct_sd_proton_scattering_rates(
            background, model, radius_cm, fixed_speed_cm_s)
            .total_rate_s_inv;
    if (!std::isfinite(total_rate_s_inv) || total_rate_s_inv <= 0.0) {
        throw std::runtime_error("fixed-state scattering rate must be positive");
    }

    const std::uint32_t collision_seed = seed ^ collisionSeedMask;
    std::mt19937 collision_rng(collision_seed);
    std::vector<std::uint64_t> target_counts(background.target_count(), 0);
    std::vector<std::uint64_t> tail_target_counts(background.target_count(), 0);
    std::array<std::uint64_t, outgoingSpeedHistogramEdges.size()>
        outgoing_speed_histogram_counts{{}};
    std::array<std::uint64_t, outgoingCosineHistogramEdges.size() - 1>
        outgoing_cosine_histogram_counts{{}};
    std::array<std::uint64_t, outgoingCmCosineHistogramEdges.size() - 1>
        outgoing_cm_cosine_histogram_counts{{}};
    std::array<std::array<std::uint64_t,
                          canonicalComponentCutpoints.size() + 1>, 3>
        canonical_target_component_histogram_counts{{}};
    std::array<std::array<std::uint64_t,
                          canonicalComponentCutpoints.size() + 1>, 3>
        canonical_outgoing_component_histogram_counts{{}};
    std::array<Moment, vectorMomentCount> vector_moments;
    std::array<Moment, vectorMomentCount> canonical_vector_moments;
    Moment relative_speed_over_escape;
    Moment outgoing_speed_over_escape;
    Moment outgoing_lab_cosine;
    std::uint64_t tail_total_count = 0;
    std::uint64_t tail_first_half_count = 0;
    std::uint64_t domain_overflow_count_1p5_vesc = 0;
    std::uint64_t domain_overflow_count_2p0_vesc = 0;
    std::uint64_t domain_overflow_count_3p0_vesc = 0;
    std::uint64_t undefined_outgoing_lab_cosine_count = 0;

    for (std::uint64_t sample = 0; sample < sample_count; ++sample) {
        const auto collision = transport::physics::sample_sd_proton_collision(
            background, model, radius_cm, incoming, collision_rng);
        if (collision.target_index >= target_counts.size()) {
            throw std::runtime_error("sampled target index is out of range");
        }
        ++target_counts[collision.target_index];
        const Velocity& outgoing = collision.outgoing_dm_velocity_cm_s;
        const double target_mass_GeV =
            background.target(collision.target_index).mass_GeV;
        const double total_mass_GeV = darkMatterMassGeV + target_mass_GeV;
        Velocity outgoing_cm{{}};
        for (std::size_t component = 0; component < 3; ++component) {
            const double center_of_mass_velocity_cm_s =
                (darkMatterMassGeV * incoming[component] +
                 target_mass_GeV * collision.target_velocity_cm_s[component]) /
                total_mass_GeV;
            outgoing_cm[component] =
                outgoing[component] - center_of_mass_velocity_cm_s;
        }
        const double outgoing_cm_speed = norm(outgoing_cm);
        if (!std::isfinite(outgoing_cm_speed) || outgoing_cm_speed <= 0.0) {
            throw std::runtime_error("outgoing CM direction is undefined");
        }
        const double raw_cm_cosine = dot(outgoing_cm, incoming) /
            (outgoing_cm_speed * fixed_speed_cm_s);
        if (!std::isfinite(raw_cm_cosine) ||
            raw_cm_cosine < -1.0 - 1.0e-10 ||
            raw_cm_cosine > 1.0 + 1.0e-10) {
            throw std::runtime_error("outgoing CM angle cosine is invalid");
        }
        const double cm_cosine = std::max(-1.0,
            std::min(1.0, raw_cm_cosine));
        std::size_t cm_cosine_bin = histogram_bin(
            cm_cosine, outgoingCmCosineHistogramEdges);
        cm_cosine_bin = std::min(cm_cosine_bin,
            outgoing_cm_cosine_histogram_counts.size() - 1);
        ++outgoing_cm_cosine_histogram_counts[cm_cosine_bin];
        for (std::size_t component = 0; component < 3; ++component) {
            const double target_canonical =
                dot(canonical_axes[component], collision.target_velocity_cm_s) /
                canonical_component_scale_cm_s;
            const double outgoing_canonical =
                dot(canonical_axes[component], outgoing) /
                canonical_component_scale_cm_s;
            const auto target_boundary = std::upper_bound(
                canonicalComponentCutpoints.begin(),
                canonicalComponentCutpoints.end(), target_canonical);
            const auto outgoing_boundary = std::upper_bound(
                canonicalComponentCutpoints.begin(),
                canonicalComponentCutpoints.end(), outgoing_canonical);
            ++canonical_target_component_histogram_counts[component][
                static_cast<std::size_t>(
                    target_boundary - canonicalComponentCutpoints.begin())];
            ++canonical_outgoing_component_histogram_counts[component][
                static_cast<std::size_t>(
                    outgoing_boundary - canonicalComponentCutpoints.begin())];
        }
        const double outgoing_speed = norm(outgoing);
        const double speed_ratio = outgoing_speed / escape_speed_cm_s;
        const double relative_speed = norm(subtract(
            incoming, collision.target_velocity_cm_s));
        relative_speed_over_escape.add(relative_speed / escape_speed_cm_s);
        outgoing_speed_over_escape.add(speed_ratio);
        ++outgoing_speed_histogram_counts[histogram_bin(
            speed_ratio, outgoingSpeedHistogramEdges)];

        if (outgoing_speed > escape_speed_cm_s) {
            ++tail_total_count;
            if (sample < sample_count / 2) {
                ++tail_first_half_count;
            }
            ++tail_target_counts[collision.target_index];
        }
        if (speed_ratio > 1.5) {
            ++domain_overflow_count_1p5_vesc;
        }
        if (speed_ratio > 2.0) {
            ++domain_overflow_count_2p0_vesc;
        }
        if (speed_ratio > 3.0) {
            ++domain_overflow_count_3p0_vesc;
        }

        if (outgoing_speed > 0.0) {
            const double cosine = std::max(-1.0, std::min(1.0,
                dot(incoming, outgoing) /
                (fixed_speed_cm_s * outgoing_speed)));
            outgoing_lab_cosine.add(cosine);
            std::size_t cosine_bin = histogram_bin(
                cosine, outgoingCosineHistogramEdges);
            cosine_bin = std::min(cosine_bin,
                outgoing_cosine_histogram_counts.size() - 1);
            ++outgoing_cosine_histogram_counts[cosine_bin];
        } else {
            ++undefined_outgoing_lab_cosine_count;
        }

        const long double x = outgoing[0] / escape_speed_cm_s;
        const long double y = outgoing[1] / escape_speed_cm_s;
        const long double z = outgoing[2] / escape_speed_cm_s;
        const std::array<long double, vectorMomentCount> values{{
            x, y, z, x * x, x * y, x * z, y * y, y * z, z * z
        }};
        const long double cx = dot(canonical_axes[0], outgoing) /
            escape_speed_cm_s;
        const long double cy = dot(canonical_axes[1], outgoing) /
            escape_speed_cm_s;
        const long double cz = dot(canonical_axes[2], outgoing) /
            escape_speed_cm_s;
        const std::array<long double, vectorMomentCount> canonical_values{{
            cx, cy, cz, cx * cx, cx * cy, cx * cz,
            cy * cy, cy * cz, cz * cz
        }};
        for (std::size_t index = 0; index < values.size(); ++index) {
            vector_moments[index].add(values[index]);
            canonical_vector_moments[index].add(canonical_values[index]);
        }
    }

    print_common("fixed", sample_count, radius_Rsun, radius_cm,
                 bath_temperature_K, seed, collision_seed, target_counts);
    std::cout << ",\"fixed_speed_over_escape\":" << fixed_speed_over_escape
              << ",\"escape_speed_cm_s\":" << escape_speed_cm_s
              << ",\"direction_unit\":";
    print_array(direction);
    std::cout << ",\"total_rate_s_inv\":" << total_rate_s_inv
              << ",\"canonical_component_scale_cm_s\":"
              << canonical_component_scale_cm_s
              << ",\"canonical_component_histogram_cutpoints\":";
    print_array(canonicalComponentCutpoints);
    std::cout << ",\"canonical_component_histogram_rule\":\"(-inf,c0),[c0,c1),...,[c_last,inf), x/y/z order after inverse shortest proper rotation mapping z to direction\""
              << ",\"canonical_target_component_histogram_counts\":";
    print_histogram3(canonical_target_component_histogram_counts);
    std::cout << ",\"canonical_outgoing_component_histogram_counts\":";
    print_histogram3(canonical_outgoing_component_histogram_counts);
    std::cout
              << ",\"tail_total_count\":" << tail_total_count
              << ",\"tail_first_half_count\":" << tail_first_half_count
              << ",\"tail_first_half_sample_count\":" << sample_count / 2
              << ",\"tail_target_counts\":";
    print_array(tail_target_counts);
    std::cout << ",\"domain_overflow_count_1p5_vesc\":"
              << domain_overflow_count_1p5_vesc
              << ",\"domain_overflow_count_2p0_vesc\":"
              << domain_overflow_count_2p0_vesc
              << ",\"domain_overflow_count_3p0_vesc\":"
              << domain_overflow_count_3p0_vesc
              << ",\"outgoing_speed_over_escape_histogram_edges\":";
    print_array(outgoingSpeedHistogramEdges);
    std::cout << ",\"outgoing_speed_over_escape_histogram_counts\":";
    print_array(outgoing_speed_histogram_counts);
    std::cout << ",\"outgoing_speed_histogram_last_open_ended\":true"
              << ",\"outgoing_lab_cosine_histogram_edges\":";
    print_array(outgoingCosineHistogramEdges);
    std::cout << ",\"outgoing_lab_cosine_histogram_counts\":";
    print_array(outgoing_cosine_histogram_counts);
    std::cout << ",\"outgoing_lab_cosine_right_endpoint_in_last_bin\":true"
              << ",\"cm_scattering_cosine_histogram_edges\":";
    print_array(outgoingCmCosineHistogramEdges);
    std::cout << ",\"cm_scattering_cosine_histogram_counts\":";
    print_array(outgoing_cm_cosine_histogram_counts);
    std::cout << ",\"cm_scattering_cosine_right_endpoint_in_last_bin\":true"
              << ",\"undefined_outgoing_lab_cosine_count\":"
              << undefined_outgoing_lab_cosine_count
              << ",\"relative_speed_over_escape\":";
    print_moment(relative_speed_over_escape);
    std::cout << ",\"outgoing_speed_over_escape\":";
    print_moment(outgoing_speed_over_escape);
    std::cout << ",\"outgoing_lab_cosine\":";
    print_moment(outgoing_lab_cosine);
    std::cout << ",\"outgoing_vector_over_escape_moments\":{";
    for (std::size_t index = 0; index < vector_moments.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << '\"' << vectorMomentNames[index] << "\":";
        print_moment(vector_moments[index]);
    }
    std::cout << "},\"canonical_outgoing_vector_over_escape_moments\":{";
    for (std::size_t index = 0;
         index < canonical_vector_moments.size(); ++index) {
        if (index != 0) {
            std::cout << ',';
        }
        std::cout << '\"' << vectorMomentNames[index] << "\":";
        print_moment(canonical_vector_moments[index]);
    }
    std::cout << "}}\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 8 && argc != 10) {
            throw std::invalid_argument(
                "usage: t03_validation_sampler SOLAR_MODEL NUCLEAR_DATA "
                "thermal|fixed RADIUS_RSUN TCHI_OVER_T|V_OVER_VESC "
                "SAMPLES SEED [--direction x,y,z]");
        }
        const std::string mode(argv[3]);
        if (mode != "thermal" && mode != "fixed") {
            throw std::invalid_argument("mode must be thermal or fixed");
        }
        const double radius_Rsun = parse_double(argv[4], "radius_Rsun");
        if (radius_Rsun < 0.0 || radius_Rsun > 1.0) {
            throw std::invalid_argument("radius_Rsun must be in [0,1]");
        }
        const double parameter = parse_double(argv[5], "mode parameter");
        const std::uint64_t sample_count = parse_uint64(argv[6], "samples");
        const std::uint64_t parsed_seed = parse_uint64(argv[7], "seed");
        if (sample_count == 0) {
            throw std::invalid_argument("samples must be positive");
        }
        if (parsed_seed > std::numeric_limits<std::uint32_t>::max()) {
            throw std::invalid_argument("seed must fit std::mt19937 seed width");
        }
        if (argc == 10 && (mode != "fixed" ||
                           std::string(argv[8]) != "--direction")) {
            throw std::invalid_argument(
                "--direction x,y,z is supported only in fixed mode");
        }
        const Velocity direction = argc == 10
            ? parse_direction(argv[9])
            : Velocity{{0.0, 0.0, 1.0}};
        const Background background(argv[1], argv[2]);
        std::cout << std::setprecision(
            std::numeric_limits<long double>::max_digits10);
        if (mode == "thermal") {
            run_thermal(background, radius_Rsun, parameter,
                        sample_count, static_cast<std::uint32_t>(parsed_seed));
        } else {
            run_fixed(background, radius_Rsun, parameter,
                      sample_count, static_cast<std::uint32_t>(parsed_seed),
                      direction);
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "T03 validation sampler failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
