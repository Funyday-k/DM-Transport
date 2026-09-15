#include "transport/PhaseSpaceProjector.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using transport::PhaseSpaceGrid;
using transport::PhaseSpaceProjector;

static_assert(
    std::is_same<
        decltype(std::declval<const PhaseSpaceProjector&>()
                     .source_particles_s()),
        const std::vector<double>&>::value,
    "projected cell sources must be accessible without mutation or copying");

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

void require_unchanged(const PhaseSpaceProjector& projector,
                       const std::vector<double>& expected_sources,
                       double expected_total, const std::string& message) {
    require(projector.source_particles_s() == expected_sources,
            message + ": cell sources changed");
    require(projector.total_source_particles_s() == expected_total,
            message + ": total source changed");
}

void test_containing_cell_deposition_and_conservation() {
    const PhaseSpaceGrid grid({0.0, 1.0, 4.0}, {0.0, 2.0, 7.0},
                              {-1.0, -0.4, 0.3, 1.0});
    PhaseSpaceProjector projector(grid);
    require(projector.source_particles_s().size() == grid.size(),
            "source vector follows grid size");
    require(std::all_of(projector.source_particles_s().begin(),
                        projector.source_particles_s().end(),
                        [](double source) { return source == 0.0; }),
            "source vector starts at zero");

    projector.deposit(0.25, 0.5, -0.8, 1.25);
    projector.deposit(1.0, 2.0, -0.4, 2.5);
    projector.deposit(4.0, 7.0, 1.0, 4.0);
    projector.deposit(3.0, 5.0, 0.8, 0.75);

    const std::size_t first = grid.flatten(0, 0, 0);
    const std::size_t internal_faces = grid.flatten(1, 1, 1);
    const std::size_t last = grid.flatten(1, 1, 2);
    const auto& source = projector.source_particles_s();
    require(source[first] == 1.25,
            "interior state deposits its full integrated source");
    require(source[internal_faces] == 2.5,
            "internal faces follow the grid's right-cell ownership");
    require(source[last] == 4.75,
            "closed global upper faces and same-cell points accumulate");

    double cell_sum = 0.0;
    for (double value : source) {
        require(std::isfinite(value) && value >= 0.0,
                "projected cell source remains finite and nonnegative");
        cell_sum += value;
    }
    require(cell_sum == 8.5 &&
                projector.total_source_particles_s() == cell_sum,
            "containing-cell deposition conserves total particles/s");
}

void test_reset() {
    const PhaseSpaceGrid grid({0.0, 1.0}, {0.0, 1.0}, {-1.0, 1.0});
    PhaseSpaceProjector projector(grid);
    projector.deposit(0.5, 0.5, 0.0, 3.0);
    projector.reset();
    require(projector.total_source_particles_s() == 0.0,
            "reset clears total source");
    require(std::all_of(projector.source_particles_s().begin(),
                        projector.source_particles_s().end(),
                        [](double source) { return source == 0.0; }),
            "reset clears every cell source");
}

void test_invalid_deposits_are_transactional() {
    const PhaseSpaceGrid grid({0.0, 1.0, 2.0}, {0.0, 1.0, 2.0},
                              {-1.0, 0.0, 1.0});
    PhaseSpaceProjector projector(grid);
    projector.deposit(0.5, 0.5, -0.5, 4.0);
    const std::vector<double> expected_sources =
        projector.source_particles_s();
    const double expected_total = projector.total_source_particles_s();
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();

    for (double invalid_weight : {-1.0, infinity, -infinity, nan}) {
        require_throws<std::invalid_argument>([&] {
            projector.deposit(0.5, 0.5, -0.5, invalid_weight);
        }, "invalid source weight");
        require_unchanged(projector, expected_sources, expected_total,
                          "invalid source weight is transactional");
    }

    const std::vector<std::vector<double>> invalid_states{
        {-0.1, 0.5, -0.5}, {2.1, 0.5, -0.5},
        {0.5, -0.1, -0.5}, {0.5, 2.1, -0.5},
        {0.5, 0.5, -1.1}, {0.5, 0.5, 1.1},
        {nan, 0.5, -0.5}, {0.5, infinity, -0.5},
        {0.5, 0.5, nan}
    };
    for (const auto& state : invalid_states) {
        require_throws<std::out_of_range>([&] {
            projector.deposit(state[0], state[1], state[2], 1.0);
        }, "invalid physical source state");
        require_unchanged(projector, expected_sources, expected_total,
                          "invalid source state is transactional");
    }
}

void test_overflow_is_transactional() {
    const PhaseSpaceGrid grid({0.0, 1.0, 2.0}, {0.0, 1.0},
                              {-1.0, 1.0});
    PhaseSpaceProjector projector(grid);
    const double maximum = std::numeric_limits<double>::max();
    projector.deposit(0.5, 0.5, 0.0, maximum);
    const std::vector<double> expected_sources =
        projector.source_particles_s();
    const double expected_total = projector.total_source_particles_s();

    require_throws<std::overflow_error>([&] {
        // The destination cell would remain finite, but the total would not.
        projector.deposit(1.5, 0.5, 0.0, maximum);
    }, "total source overflow");
    require_unchanged(projector, expected_sources, expected_total,
                      "overflow is transactional");

    PhaseSpaceProjector same_cell_projector(grid);
    same_cell_projector.deposit(0.5, 0.5, 0.0, maximum);
    const std::vector<double> same_cell_sources =
        same_cell_projector.source_particles_s();
    require_throws<std::overflow_error>([&] {
        same_cell_projector.deposit(0.5, 0.5, 0.0, maximum);
    }, "cell source overflow");
    require_unchanged(same_cell_projector, same_cell_sources, maximum,
                      "cell overflow is transactional");
}

}  // namespace

int main() {
    try {
        test_containing_cell_deposition_and_conservation();
        test_reset();
        test_invalid_deposits_are_transactional();
        test_overflow_is_transactional();
        std::cout << "phase-space projector tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "phase-space projector test failed: " << error.what()
                  << '\n';
        return EXIT_FAILURE;
    }
}
