#include "transport/PhaseSpaceCellGeometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using transport::CellQuadrature;
using transport::CellQuadratureOrder;
using transport::EscapeThresholdClass;
using transport::PhaseSpaceGrid;
using transport::physics::SolarBackground;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, const std::string& message,
                   double relative_tolerance = 1.0e-13) {
    require(std::isfinite(actual) && std::isfinite(expected),
            message + ": nonfinite value");
    const double scale = std::max(1.0, std::max(std::abs(actual),
                                                std::abs(expected)));
    require(std::abs(actual - expected) <= relative_tolerance * scale,
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

void check_rule_in_cell(const PhaseSpaceGrid& grid, std::size_t flat,
                        const CellQuadrature& rule, std::size_t expected_size) {
    require(rule.size == expected_size, "quadrature node count");
    const auto bounds = grid.cell_bounds(flat);
    double weight_sum = 0.0;
    for (std::size_t i = 0; i < rule.size; ++i) {
        const auto& node = rule.nodes[i];
        require(node.r_cm >= bounds.r_cm[0] &&
                node.r_cm <= bounds.r_cm[1], "quadrature radius in cell");
        require(node.v_cm_s >= bounds.v_cm_s[0] &&
                node.v_cm_s <= bounds.v_cm_s[1], "quadrature speed in cell");
        require(node.mu >= bounds.mu[0] && node.mu <= bounds.mu[1],
                "quadrature mu in cell");
        require(std::isfinite(node.weight) && node.weight > 0.0,
                "strictly positive quadrature weight");
        weight_sum += node.weight;
    }
    require_close(weight_sum, 1.0, "quadrature weights sum to one", 0.0);
    require_close(weight_sum * grid.cell_phase_measure(flat),
                  grid.cell_phase_measure(flat),
                  "quadrature recovers the full phase measure", 0.0);
}

void test_positive_quadrature() {
    const PhaseSpaceGrid grid({0.0, 2.0}, {0.0, 3.0},
                              {-1.0, -0.4, 1.0});
    for (std::size_t flat = 0; flat < grid.size(); ++flat) {
        const auto one = transport::positive_cell_quadrature(
            grid, flat, CellQuadratureOrder::one_point);
        const auto two = transport::positive_cell_quadrature(
            grid, flat, CellQuadratureOrder::two_point);
        check_rule_in_cell(grid, flat, one, 1);
        check_rule_in_cell(grid, flat, two, 8);
        require_close(one.nodes[0].weight, 1.0, "one-point weight", 0.0);

        // r^3 and v^3 are uniform under the FV measure. Their means are
        // half the upper-face cubes when the lower faces are zero.
        for (const CellQuadrature* rule : {&one, &two}) {
            double r_cubed_mean = 0.0;
            double v_cubed_mean = 0.0;
            double mu_mean = 0.0;
            for (std::size_t i = 0; i < rule->size; ++i) {
                const auto& node = rule->nodes[i];
                r_cubed_mean += node.weight * std::pow(node.r_cm, 3);
                v_cubed_mean += node.weight * std::pow(node.v_cm_s, 3);
                mu_mean += node.weight * node.mu;
            }
            require_close(r_cubed_mean, 4.0, "measure-mean r cubed");
            require_close(v_cubed_mean, 13.5, "measure-mean v cubed");
            const auto bounds = grid.cell_bounds(flat);
            require_close(mu_mean, 0.5 * (bounds.mu[0] + bounds.mu[1]),
                          "uniform-mu midpoint");
        }

        double second_radial_moment = 0.0;
        double second_angular_moment = 0.0;
        double coupled_moment = 0.0;
        const auto bounds = grid.cell_bounds(flat);
        const double expected_mu_squared =
            (bounds.mu[0] * bounds.mu[0] +
             bounds.mu[0] * bounds.mu[1] +
             bounds.mu[1] * bounds.mu[1]) / 3.0;
        for (std::size_t i = 0; i < two.size; ++i) {
            const auto& node = two.nodes[i];
            second_radial_moment += node.weight * std::pow(node.r_cm, 6);
            second_angular_moment += node.weight * node.mu * node.mu;
            coupled_moment += node.weight * std::pow(node.r_cm, 6) *
                std::pow(node.v_cm_s, 6) * node.mu * node.mu;
        }
        require_close(second_radial_moment, 64.0 / 3.0,
                      "two-point rule integrates quadratic r-cube moment");
        require_close(second_angular_moment, expected_mu_squared,
                      "two-point rule integrates quadratic mu moment");
        require_close(coupled_moment,
                      (64.0 / 3.0) * (729.0 / 3.0) * expected_mu_squared,
                      "tensor rule integrates coupled r/v/mu moment");
        require(std::abs(std::pow(one.nodes[0].r_cm, 6) - 64.0 / 3.0) > 1.0,
                "one-point rule has a visible higher-moment error");
        require(std::abs(one.nodes[0].mu * one.nodes[0].mu -
                         expected_mu_squared) > 1.0e-3,
                "one-point rule has a visible angular-moment error");
    }

    const PhaseSpaceGrid shell({1.0, 2.0}, {2.0, 4.0}, {-1.0, 1.0});
    const auto two = transport::positive_cell_quadrature(
        shell, 0, CellQuadratureOrder::two_point);
    check_rule_in_cell(shell, 0, two, 8);
    double radial_sixth_moment = 0.0;
    for (std::size_t i = 0; i < two.size; ++i) {
        radial_sixth_moment += two.nodes[i].weight *
            std::pow(two.nodes[i].r_cm, 6);
    }
    require_close(radial_sixth_moment, 511.0 / 21.0,
                  "nonzero lower-face quadratic measure moment");

    const double huge = std::ldexp(1.0, 500);
    const double tiny = std::ldexp(1.0, -500);
    const PhaseSpaceGrid scaled({0.0, huge}, {0.0, tiny}, {-1.0, 1.0});
    check_rule_in_cell(scaled, 0,
                       transport::positive_cell_quadrature(
                           scaled, 0, CellQuadratureOrder::two_point), 8);
    const double adjacent = std::nextafter(1.0, 2.0);
    const PhaseSpaceGrid thin({1.0, adjacent}, {0.0, 1.0}, {-1.0, 1.0});
    check_rule_in_cell(thin, 0,
                       transport::positive_cell_quadrature(
                           thin, 0, CellQuadratureOrder::two_point), 8);

    require_throws<std::out_of_range>([&] {
        transport::positive_cell_quadrature(
            grid, grid.size(), CellQuadratureOrder::one_point);
    }, "invalid flat cell index");
    require_throws<std::invalid_argument>([&] {
        transport::positive_cell_quadrature(
            grid, 0, static_cast<CellQuadratureOrder>(99));
    }, "invalid quadrature order");
}

EscapeThresholdClass classify(const SolarBackground& background,
                              double r_lo_cm, double r_hi_cm,
                              double v_lo_cm_s, double v_hi_cm_s) {
    const PhaseSpaceGrid grid({r_lo_cm, r_hi_cm},
                              {v_lo_cm_s, v_hi_cm_s}, {-1.0, 1.0});
    return transport::classify_escape_threshold_cell(
        grid, 0, 0, background).classification;
}

void test_escape_threshold(const SolarBackground& background) {
    const double solar_radius_cm = background.solar_radius_cm();
    double previous_escape = background.escape_speed_cm_s(0.0);
    for (int i = 1; i <= 4000; ++i) {
        const double radius_cm = solar_radius_cm * i / 4000.0;
        const double escape = background.escape_speed_cm_s(radius_cm);
        require(escape <= previous_escape,
                "interpolated solar escape profile decreases on a dense mesh");
        previous_escape = escape;
    }
    const double r_lo_cm = 0.3 * solar_radius_cm;
    const double r_hi_cm = 0.7 * solar_radius_cm;
    const double escape_inner = background.escape_speed_cm_s(r_lo_cm);
    const double escape_outer = background.escape_speed_cm_s(r_hi_cm);
    require(escape_inner > escape_outer, "solar escape speed decreases radially");

    require(classify(background, r_lo_cm, r_hi_cm, 0.0, escape_outer) ==
                EscapeThresholdClass::fully_bound,
            "cell touching E=0 at its outer, upper-speed corner is bound");
    require(classify(background, r_lo_cm, r_hi_cm, escape_inner,
                     1.2 * escape_inner) ==
                EscapeThresholdClass::fully_unbound,
            "cell touching E=0 at its inner, lower-speed corner is unbound");
    const double speed_gap = escape_inner - escape_outer;
    require(classify(background, r_lo_cm, r_hi_cm,
                     escape_outer + speed_gap / 3.0,
                     escape_outer + 2.0 * speed_gap / 3.0) ==
                EscapeThresholdClass::threshold_crossing,
            "radial escape boundary crosses a speed cell with no corner at E=0");
    require(classify(background, r_lo_cm, r_hi_cm,
                     0.5 * escape_outer, 1.1 * escape_inner) ==
                EscapeThresholdClass::threshold_crossing,
            "wide cell straddles both energy signs");

    // Check the stated whole-cell interpretation against an independent
    // grid of physical states, including all faces and both solar regimes.
    for (double radius_scale : {0.1, 1.0}) {
        const double inner = radius_scale * solar_radius_cm;
        const double outer = (radius_scale + 0.2) * solar_radius_cm;
        const double e_inner = background.escape_speed_cm_s(inner);
        const double e_outer = background.escape_speed_cm_s(outer);
        for (const auto& speeds : {
                 std::array<double, 2>{{0.0, 0.9 * e_outer}},
                 std::array<double, 2>{{1.1 * e_inner, 1.3 * e_inner}}}) {
            const auto kind = classify(background, inner, outer,
                                       speeds[0], speeds[1]);
            const bool bound = kind == EscapeThresholdClass::fully_bound;
            require(bound || kind == EscapeThresholdClass::fully_unbound,
                    "test cell has a definite threshold side");
            for (int ir = 0; ir <= 20; ++ir) {
                const double r = inner + (outer - inner) * ir / 20.0;
                const double escape = background.escape_speed_cm_s(r);
                for (int iv = 0; iv <= 20; ++iv) {
                    const double v = speeds[0] +
                        (speeds[1] - speeds[0]) * iv / 20.0;
                    require(bound ? v <= escape : v >= escape,
                            "whole-cell threshold label agrees with sampled E sign");
                }
            }
        }
    }

    const PhaseSpaceGrid grid({r_lo_cm, r_hi_cm},
                              {0.0, escape_outer}, {-1.0, 1.0});
    const auto geometry = transport::classify_escape_threshold_cell(
        grid, 0, 0, background);
    require_close(geometry.inner_escape_speed_cm_s, escape_inner,
                  "inner escape-speed metadata", 0.0);
    require_close(geometry.outer_escape_speed_cm_s, escape_outer,
                  "outer escape-speed metadata", 0.0);
    require_throws<std::out_of_range>([&] {
        transport::classify_escape_threshold_cell(grid, 1, 0, background);
    }, "invalid radial cell index");
    require_throws<std::out_of_range>([&] {
        transport::classify_escape_threshold_cell(grid, 0, 1, background);
    }, "invalid speed cell index");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "expected solar-model and nuclear-data paths");
        test_positive_quadrature();
        const SolarBackground background(argv[1], argv[2]);
        test_escape_threshold(background);
        std::cout << "phase-space cell geometry tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "phase-space cell geometry test failed: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
