#include "transport/CollisionKernel.hpp"

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

using transport::CellQuadratureOrder;
using transport::CollisionKernelRow;
using transport::CollisionNodeHistogram;
using transport::PhaseSpaceGrid;
using transport::physics::SdProtonModel;
using transport::physics::SolarBackground;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected,
                   const std::string& message,
                   double relative_tolerance = 1.0e-13) {
    require(std::isfinite(actual), message + ": nonfinite result");
    const double scale = std::max(1.0, std::abs(expected));
    require(std::abs(actual - expected) <= relative_tolerance * scale,
            message + ": value differs");
}

template <typename Exception, typename Callable>
void require_throws(Callable operation, const std::string& message) {
    try {
        operation();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message + ": expected exception");
}

double rate_to(const CollisionKernelRow& row, std::size_t destination) {
    for (const auto& transition : row.off_diagonal) {
        if (transition.destination == destination) {
            return transition.rate_s_inv;
        }
    }
    return 0.0;
}

void check_ledger(const CollisionKernelRow& row) {
    double leaving = row.numerical_velocity_overflow_rate_s_inv;
    std::size_t previous = 0;
    bool first = true;
    for (const auto& entry : row.off_diagonal) {
        require(entry.destination != row.source,
                "self events must cancel from the generator");
        require(entry.rate_s_inv > 0.0 && std::isfinite(entry.rate_s_inv),
                "off-diagonal rates must be positive and finite");
        require(first || previous < entry.destination,
                "off-diagonal destinations must be strictly sorted");
        first = false;
        previous = entry.destination;
        leaving += entry.rate_s_inv;
    }
    require_close(row.diagonal_s_inv + leaving, 0.0,
                  "extended generator row ledger");
    require_close(row.mean_event_rate_s_inv,
                  row.self_event_rate_s_inv + leaving,
                  "all event channels exhaust the mean rate");
    if (row.mean_event_rate_s_inv > 0.0) {
        require_close(row.numerical_velocity_overflow_probability,
                      row.numerical_velocity_overflow_rate_s_inv /
                          row.mean_event_rate_s_inv,
                      "numerical overflow probability");
    }
}

void test_rate_weighted_assembly_and_overflow() {
    const PhaseSpaceGrid grid({0.0, 1.0}, {0.0, 1.0, 2.0, 3.0},
                              {-1.0, 1.0});
    const std::vector<CollisionNodeHistogram> nodes{
        {0.5, 1.0, 4, {{2, 1}, {0, 2}, {1, 1}}, 0},
        {0.5, 3.0, 4, {{2, 1}, {1, 2}}, 1}};
    const CollisionKernelRow row = transport::assemble_collision_kernel_row(
        grid, 0, nodes, 4);
    require(row.source == 0 && row.samples_per_node == 4,
            "row provenance");
    require(row.off_diagonal.size() == 2,
            "two non-self destinations remain");
    require(row.off_diagonal[0].destination == 1 &&
                row.off_diagonal[1].destination == 2,
            "unsorted event bins are sorted in output");
    require_close(row.mean_event_rate_s_inv, 2.0, "cell-average event rate");
    require_close(rate_to(row, 1), 0.875, "rate-weighted destination one");
    require_close(rate_to(row, 2), 0.5, "rate-weighted destination two");
    require_close(row.self_event_rate_s_inv, 0.25,
                  "self event diagnostic");
    require_close(row.numerical_velocity_overflow_rate_s_inv, 0.375,
                  "numerical overflow rate");
    require_close(row.numerical_velocity_overflow_probability, 0.1875,
                  "numerical overflow probability");
    require_close(row.diagonal_s_inv, -1.75,
                  "diagonal includes off-diagonal and overflow only");
    require(std::abs(rate_to(row, 1) - 2.0 * (0.25 + 0.5) / 2.0) > 0.1,
            "Gamma times unweighted mean P must not replace mean Gamma P");
    check_ledger(row);

    const std::vector<CollisionNodeHistogram> zero_nodes{
        {0.5, 0.0, 0, {}, 0}, {0.5, 0.0, 0, {}, 0}};
    const CollisionKernelRow zero = transport::assemble_collision_kernel_row(
        grid, 0, zero_nodes, 4);
    require(zero.off_diagonal.empty(), "zero-rate row has no transitions");
    require(zero.diagonal_s_inv == 0.0 &&
                zero.mean_event_rate_s_inv == 0.0 &&
                zero.numerical_velocity_overflow_probability == 0.0,
            "zero-rate row and probability are exact zero");
    check_ledger(zero);

    auto incomplete = nodes;
    incomplete[1].numerical_velocity_overflow_count = 0;
    require_throws<std::invalid_argument>([&] {
        transport::assemble_collision_kernel_row(grid, 0, incomplete, 4);
    }, "unaccounted event is rejected");
    auto bad_weights = nodes;
    bad_weights[0].quadrature_weight = 0.4;
    require_throws<std::invalid_argument>([&] {
        transport::assemble_collision_kernel_row(grid, 0, bad_weights, 4);
    }, "quadrature weights must normalize");
}

void test_thin_radial_face_and_unbound_deposition(
    const SolarBackground& background) {
    const double r_lower = 0.3 * background.solar_radius_cm();
    const double r_upper = std::nextafter(
        r_lower, std::numeric_limits<double>::infinity());
    const double escape = background.escape_speed_cm_s(r_lower);
    const PhaseSpaceGrid grid(
        {r_lower, r_upper, 0.4 * background.solar_radius_cm()},
        {0.0, 0.8 * escape, escape, 2.0 * escape},
        {-1.0, 0.0, 1.0});
    const std::size_t source = grid.flatten(0, 1, 0);
    const auto rule = transport::positive_cell_quadrature(
        grid, source, CellQuadratureOrder::two_point);
    bool rounded_to_upper_face = false;
    for (std::size_t q = 0; q < rule.size; ++q) {
        rounded_to_upper_face = rounded_to_upper_face ||
            rule.nodes[q].r_cm == r_upper;
    }
    require(rounded_to_upper_face,
            "test exercises quadrature radius rounded onto upper face");

    const SdProtonModel model{0.1, 1.0e-34};
    const CollisionKernelRow row = transport::build_reference_collision_row(
        grid, background, model, source, CellQuadratureOrder::two_point,
        64, 20260916);
    require(row.source == source && row.mean_event_rate_s_inv > 0.0,
            "real-physics row assembled from eight nodes");
    bool unbound_destination = false;
    for (const auto& transition : row.off_diagonal) {
        const auto index = grid.unflatten(transition.destination);
        require(index[0] == 0,
                "collision deposition retains the source radial cell");
        if (index[1] == 2) {
            unbound_destination = true;
        }
    }
    require(unbound_destination,
            "above-escape output inside v_max remains a grid destination");
    check_ledger(row);

    const CollisionKernelRow repeated = transport::build_reference_collision_row(
        grid, background, model, source, CellQuadratureOrder::two_point,
        64, 20260916);
    require(row.off_diagonal.size() == repeated.off_diagonal.size() &&
                row.mean_event_rate_s_inv == repeated.mean_event_rate_s_inv &&
                row.diagonal_s_inv == repeated.diagonal_s_inv &&
                row.self_event_rate_s_inv == repeated.self_event_rate_s_inv &&
                row.numerical_velocity_overflow_rate_s_inv ==
                    repeated.numerical_velocity_overflow_rate_s_inv,
            "same logical sample set reproduces identical row metadata");
    for (std::size_t i = 0; i < row.off_diagonal.size(); ++i) {
        require(row.off_diagonal[i].destination ==
                    repeated.off_diagonal[i].destination &&
                    row.off_diagonal[i].rate_s_inv ==
                        repeated.off_diagonal[i].rate_s_inv,
                "same logical sample set reproduces identical transitions");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "expected solar-model and nuclear-data paths");
        test_rate_weighted_assembly_and_overflow();
        const SolarBackground background(argv[1], argv[2]);
        test_thin_radial_face_and_unbound_deposition(background);
        std::cout << "collision kernel tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "collision kernel test failed: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
