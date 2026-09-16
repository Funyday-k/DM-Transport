#include "transport/CollisionKernel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>

namespace transport {
namespace {

constexpr std::uint64_t fnv_offset = UINT64_C(14695981039346656037);
constexpr std::uint64_t fnv_prime = UINT64_C(1099511628211);
constexpr std::uint64_t collision_stage = UINT64_C(0x434f4c4c4953494f);

void fingerprint_word(std::uint64_t& fingerprint, std::uint64_t word) {
    // Fixed little-endian byte order, independent of the host byte order.
    for (unsigned int shift = 0; shift < 64; shift += 8) {
        fingerprint ^= (word >> shift) & UINT64_C(0xff);
        fingerprint *= fnv_prime;
    }
}

void fingerprint_double(std::uint64_t& fingerprint, double value) {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    fingerprint_word(fingerprint, bits);
}

std::uint64_t model_grid_fingerprint(const PhaseSpaceGrid& grid,
                                     const physics::SdProtonModel& model) {
    if (sizeof(double) != sizeof(std::uint64_t) ||
        !std::numeric_limits<double>::is_iec559) {
        throw std::runtime_error("collision RNG requires IEEE-754 binary64");
    }
    std::uint64_t result = fnv_offset;
    fingerprint_word(result, UINT64_C(0x444d54524b45524e));
    fingerprint_double(result, model.dark_matter_mass_GeV);
    fingerprint_double(result, model.proton_cross_section_cm2);
    const std::array<const std::vector<double>*, 3> axes{{
        &grid.r_faces(), &grid.v_faces(), &grid.mu_faces()}};
    for (std::size_t axis = 0; axis < axes.size(); ++axis) {
        fingerprint_word(result, static_cast<std::uint64_t>(axis));
        fingerprint_word(result,
                         static_cast<std::uint64_t>(axes[axis]->size()));
        for (double face : *axes[axis]) {
            fingerprint_double(result, face);
        }
    }
    return result;
}

std::uint64_t splitmix64(std::uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

std::mt19937 logical_collision_rng(std::uint64_t master_seed,
                                   std::uint64_t fingerprint,
                                   std::size_t source,
                                   std::size_t node,
                                   std::uint64_t sample) {
    std::uint64_t state = splitmix64(master_seed);
    for (std::uint64_t field : {
             fingerprint, static_cast<std::uint64_t>(source),
             static_cast<std::uint64_t>(node), sample, collision_stage}) {
        state = splitmix64(state ^ splitmix64(field));
    }
    const std::uint64_t first = splitmix64(state);
    const std::uint64_t second = splitmix64(first);
    std::seed_seq seed{
        static_cast<std::uint32_t>(first),
        static_cast<std::uint32_t>(first >> 32),
        static_cast<std::uint32_t>(second),
        static_cast<std::uint32_t>(second >> 32)};
    return std::mt19937(seed);
}

double checked_rate(long double value) {
    if (!std::isfinite(value) || value < 0.0L ||
        value > static_cast<long double>(
                    std::numeric_limits<double>::max())) {
        throw std::overflow_error("collision rate is not representable");
    }
    const double converted = static_cast<double>(value);
    if (value > 0.0L && converted == 0.0) {
        throw std::underflow_error("collision rate rounds to zero");
    }
    return converted;
}

std::size_t deposit_collision_output(const PhaseSpaceGrid& grid,
                                     std::size_t source_ir,
                                     const physics::CartesianVelocityCmS& v) {
    for (double component : v) {
        if (!std::isfinite(component)) {
            throw std::runtime_error("collision returned nonfinite velocity");
        }
    }
    const double speed = std::hypot(std::hypot(v[0], v[1]), v[2]);
    if (!std::isfinite(speed)) {
        throw std::overflow_error("collision speed is not representable");
    }
    if (speed < grid.v_faces().front() || speed > grid.v_faces().back()) {
        return grid.size();  // Explicit numerical velocity-domain overflow.
    }
    double mu = speed == 0.0 ? 0.0 : v[2] / speed;
    const double tolerance = 32.0 * std::numeric_limits<double>::epsilon();
    if (mu < -1.0 - tolerance || mu > 1.0 + tolerance) {
        throw std::runtime_error("collision output has invalid polar cosine");
    }
    mu = std::max(-1.0, std::min(1.0, mu));

    // The collision has no radial displacement. Use the source cell's lower
    // face solely to ask locate_cell for v/mu ownership; an actual quadrature
    // radius can round onto its upper face in a one-ULP-wide radial cell.
    const PhaseSpaceGrid::Index destination = grid.locate_cell(
        grid.r_faces()[source_ir], speed, mu);
    return grid.flatten(source_ir, destination[1], destination[2]);
}

}  // namespace

CollisionKernelRow assemble_collision_kernel_row(
    const PhaseSpaceGrid& grid,
    std::size_t source,
    const std::vector<CollisionNodeHistogram>& nodes,
    std::uint64_t samples_per_node) {
    const std::size_t source_ir = grid.unflatten(source)[0];
    if (nodes.empty() || samples_per_node == 0) {
        throw std::invalid_argument(
            "collision assembly requires nodes and samples per node");
    }

    std::map<std::size_t, long double> off_diagonal;
    long double weight_sum = 0.0L;
    long double mean_rate = 0.0L;
    long double self_rate = 0.0L;
    long double overflow_rate = 0.0L;
    for (const CollisionNodeHistogram& node : nodes) {
        if (!std::isfinite(node.quadrature_weight) ||
            node.quadrature_weight <= 0.0 ||
            !std::isfinite(node.event_rate_s_inv) ||
            node.event_rate_s_inv < 0.0) {
            throw std::invalid_argument("invalid collision node weight or rate");
        }
        weight_sum += node.quadrature_weight;
        const long double weighted_rate =
            static_cast<long double>(node.quadrature_weight) *
            static_cast<long double>(node.event_rate_s_inv);
        mean_rate += weighted_rate;

        if (node.event_rate_s_inv == 0.0) {
            if (node.sample_count != 0 || !node.cell_counts.empty() ||
                node.numerical_velocity_overflow_count != 0) {
                throw std::invalid_argument("zero-rate node has events");
            }
            continue;
        }
        if (node.sample_count != samples_per_node) {
            throw std::invalid_argument(
                "positive-rate node has wrong sample count");
        }
        std::uint64_t counted = node.numerical_velocity_overflow_count;
        if (counted > node.sample_count) {
            throw std::invalid_argument("collision overflow count exceeds samples");
        }
        std::map<std::size_t, std::uint64_t> sorted_counts;
        for (const TransitionCount& bin : node.cell_counts) {
            if (bin.destination >= grid.size() ||
                grid.unflatten(bin.destination)[0] != source_ir) {
                throw std::invalid_argument(
                    "collision destination is invalid or changes radius");
            }
            if (bin.count > node.sample_count - counted) {
                throw std::invalid_argument("collision counts exceed samples");
            }
            counted += bin.count;
            std::uint64_t& merged = sorted_counts[bin.destination];
            if (bin.count > std::numeric_limits<std::uint64_t>::max() - merged) {
                throw std::overflow_error("collision bin count overflows");
            }
            merged += bin.count;
        }
        if (counted != node.sample_count) {
            throw std::invalid_argument("collision counts do not exhaust samples");
        }
        const long double event_weight = weighted_rate /
            static_cast<long double>(node.sample_count);
        for (const auto& bin : sorted_counts) {
            const long double rate = event_weight *
                static_cast<long double>(bin.second);
            if (bin.first == source) {
                self_rate += rate;
            } else {
                off_diagonal[bin.first] += rate;
            }
        }
        overflow_rate += event_weight *
            static_cast<long double>(node.numerical_velocity_overflow_count);
    }
    if (std::abs(weight_sum - 1.0L) >
        64.0L * static_cast<long double>(
                    std::numeric_limits<double>::epsilon())) {
        throw std::invalid_argument("collision quadrature weights do not sum to one");
    }

    CollisionKernelRow row{};
    row.source = source;
    row.mean_event_rate_s_inv = checked_rate(mean_rate);
    row.self_event_rate_s_inv = checked_rate(self_rate);
    row.numerical_velocity_overflow_rate_s_inv = checked_rate(overflow_rate);
    row.samples_per_node = samples_per_node;
    long double leaving_rate =
        static_cast<long double>(row.numerical_velocity_overflow_rate_s_inv);
    for (const auto& entry : off_diagonal) {
        const double rate = checked_rate(entry.second);
        if (rate <= 0.0) {
            throw std::runtime_error("off-diagonal collision rate is not positive");
        }
        row.off_diagonal.push_back({entry.first, rate});
        leaving_rate += static_cast<long double>(rate);
    }
    row.diagonal_s_inv = -checked_rate(leaving_rate);
    row.numerical_velocity_overflow_probability =
        row.mean_event_rate_s_inv == 0.0 ? 0.0 :
        row.numerical_velocity_overflow_rate_s_inv /
            row.mean_event_rate_s_inv;
    if (!std::isfinite(row.numerical_velocity_overflow_probability) ||
        row.numerical_velocity_overflow_probability < 0.0 ||
        row.numerical_velocity_overflow_probability >
            1.0 + 64.0 * std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("invalid collision overflow probability");
    }
    row.numerical_velocity_overflow_probability = std::min(
        1.0, row.numerical_velocity_overflow_probability);
    return row;
}

CollisionKernelRow build_reference_collision_row(
    const PhaseSpaceGrid& grid,
    const physics::SolarBackground& background,
    const physics::SdProtonModel& model,
    std::size_t source,
    CellQuadratureOrder quadrature_order,
    std::uint64_t samples_per_node,
    std::uint64_t master_seed) {
    if (samples_per_node == 0) {
        throw std::invalid_argument("collision samples per node must be positive");
    }
    const std::size_t source_ir = grid.unflatten(source)[0];
    const CellQuadrature rule =
        positive_cell_quadrature(grid, source, quadrature_order);
    const std::uint64_t fingerprint = model_grid_fingerprint(grid, model);
    std::vector<CollisionNodeHistogram> histograms;
    histograms.reserve(rule.size);

    for (std::size_t q = 0; q < rule.size; ++q) {
        const CellQuadratureNode& node = rule.nodes[q];
        const double rate = physics::direct_sd_proton_scattering_rates(
            background, model, node.r_cm, node.v_cm_s).total_rate_s_inv;
        CollisionNodeHistogram histogram{};
        histogram.quadrature_weight = node.weight;
        histogram.event_rate_s_inv = rate;
        if (rate > 0.0) {
            histogram.sample_count = samples_per_node;
            std::map<std::size_t, std::uint64_t> counts;
            const double transverse = node.v_cm_s *
                std::sqrt(std::max(0.0, 1.0 - node.mu * node.mu));
            const physics::CartesianVelocityCmS incoming{{
                transverse, 0.0, node.v_cm_s * node.mu}};
            for (std::uint64_t sample = 0; sample < samples_per_node;
                 ++sample) {
                std::mt19937 rng = logical_collision_rng(
                    master_seed, fingerprint, source, q, sample);
                const physics::CollisionSample collision =
                    physics::sample_sd_proton_collision(
                        background, model, node.r_cm, incoming, rng);
                const std::size_t destination = deposit_collision_output(
                    grid, source_ir, collision.outgoing_dm_velocity_cm_s);
                if (destination == grid.size()) {
                    ++histogram.numerical_velocity_overflow_count;
                } else {
                    ++counts[destination];
                }
            }
            for (const auto& bin : counts) {
                histogram.cell_counts.push_back({bin.first, bin.second});
            }
        }
        histograms.push_back(std::move(histogram));
    }
    return assemble_collision_kernel_row(
        grid, source, histograms, samples_per_node);
}

}  // namespace transport
