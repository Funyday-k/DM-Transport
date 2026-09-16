#include "transport/PhaseSpaceCellGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace transport {
namespace {

double measure_quantile(double lower, double upper, long double fraction) {
    // Scale by upper before cubing, so neither r^3 nor v^3 needs to be
    // representable in double. The integral of x^2 dx is uniform in x^3.
    const long double lower_ratio =
        static_cast<long double>(lower) / static_cast<long double>(upper);
    const long double lower_cube =
        lower_ratio * lower_ratio * lower_ratio;
    const long double scaled_cube =
        lower_cube + fraction * (1.0L - lower_cube);
    const double result = static_cast<double>(
        static_cast<long double>(upper) * std::cbrt(scaled_cube));
    // For an interval only one representable double wide, an interior node
    // can round to a face. It still remains inside the closed FV cell.
    return result < lower ? lower : (result > upper ? upper : result);
}

struct RadialNode {
    long double x_fraction;
    long double weight;
    long double bound_y_fraction;
    long double unbound_y_fraction;
};

struct RadialPanel {
    RadialNode nodes[2];
    long double bound_weight;
    long double unbound_weight;
};

struct RootLocation {
    long double lower_x_fraction;
    long double upper_x_fraction;
    long double bracket_width;
};

struct ThresholdRule {
    EscapeThresholdMeasure measure;
    std::vector<RadialNode> radial_nodes;
};

double evaluated_escape(const EscapeSpeedEvaluator& evaluator, double radius) {
    const double escape = evaluator(radius);
    if (!std::isfinite(escape) || escape <= 0.0) {
        throw std::runtime_error("escape-speed evaluator returned a nonpositive or nonfinite value");
    }
    return escape;
}

double escape_at_x_fraction(const PhaseSpaceGrid::CellBounds& bounds,
                            long double x_fraction,
                            const EscapeSpeedEvaluator& evaluator) {
    return evaluated_escape(
        evaluator, measure_quantile(bounds.r_cm[0], bounds.r_cm[1], x_fraction));
}

// Evaluate (v_escape^3-v_lower^3)/(v_upper^3-v_lower^3) without cubing
// physical speeds or subtracting nearly equal cubes.
std::array<long double, 2> clipped_y_fractions(double lower, double upper,
                                                double escape) {
    if (escape <= lower) {
        return {{0.0L, 1.0L}};
    }
    if (escape >= upper) {
        return {{1.0L, 0.0L}};
    }
    const long double lo = static_cast<long double>(lower) / upper;
    const long double hi = static_cast<long double>(escape) / upper;
    const long double bound_numerator =
        (static_cast<long double>(escape - lower) / upper) *
        (hi * hi + hi * lo + lo * lo);
    const long double unbound_numerator =
        (static_cast<long double>(upper - escape) / upper) *
        (1.0L + hi + hi * hi);
    const long double denominator =
        (static_cast<long double>(upper - lower) / upper) *
        (1.0L + lo + lo * lo);
    const long double bound = std::max(0.0L,
        std::min(1.0L, bound_numerator / denominator));
    const long double unbound = std::max(0.0L,
        std::min(1.0L, unbound_numerator / denominator));
    // Form the larger complement from the smaller side, so a thin sliver at
    // either velocity face is not lost to subtraction from one.
    return bound <= unbound ?
        std::array<long double, 2>{{bound, 1.0L - bound}} :
        std::array<long double, 2>{{1.0L - unbound, unbound}};
}

RootLocation threshold_root(const PhaseSpaceGrid::CellBounds& bounds,
                            double inner_escape, double outer_escape,
                            double threshold,
                            const EscapeSpeedEvaluator& evaluator) {
    if (inner_escape <= threshold) {
        return {0.0L, 0.0L, 0.0L};
    }
    if (outer_escape >= threshold) {
        return {1.0L, 1.0L, 0.0L};
    }
    long double lower = 0.0L;
    long double upper = 1.0L;
    for (int iteration = 0; iteration < 96; ++iteration) {
        const long double middle = lower + 0.5L * (upper - lower);
        if (middle == lower || middle == upper) {
            break;
        }
        if (escape_at_x_fraction(bounds, middle, evaluator) >= threshold) {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    return {lower, upper, upper - lower};
}

RadialPanel gauss_panel(const PhaseSpaceGrid::CellBounds& bounds,
                       long double lower, long double upper,
                       const EscapeSpeedEvaluator& evaluator) {
    const long double displacement = std::sqrt(1.0L / 12.0L);
    const long double center = lower + 0.5L * (upper - lower);
    const long double half_width = 0.5L * (upper - lower);
    RadialPanel result{};
    result.bound_weight = 0.0L;
    result.unbound_weight = 0.0L;
    for (std::size_t index = 0; index < 2; ++index) {
        const long double node_x =
            center + (index == 0 ? -1.0L : 1.0L) *
                         2.0L * displacement * half_width;
        const double escape = escape_at_x_fraction(bounds, node_x, evaluator);
        const auto fractions = clipped_y_fractions(
            bounds.v_cm_s[0], bounds.v_cm_s[1], escape);
        result.nodes[index] = {node_x, half_width,
                               fractions[0], fractions[1]};
        result.bound_weight += half_width * fractions[0];
        result.unbound_weight += half_width * fractions[1];
    }
    return result;
}

void append_constant_panel(long double lower, long double upper,
                           double bound_fraction,
                           std::vector<RadialNode>& nodes) {
    if (!(upper > lower)) {
        return;
    }
    const long double displacement = std::sqrt(1.0L / 12.0L);
    const long double center = lower + 0.5L * (upper - lower);
    const long double half_width = 0.5L * (upper - lower);
    nodes.push_back({center - 2.0L * displacement * half_width,
                     half_width, bound_fraction, 1.0L - bound_fraction});
    nodes.push_back({center + 2.0L * displacement * half_width,
                     half_width, bound_fraction, 1.0L - bound_fraction});
}

void append_active_panels(const PhaseSpaceGrid::CellBounds& bounds,
                          long double lower, long double upper,
                          const EscapeSpeedEvaluator& evaluator,
                          unsigned int depth,
                          std::vector<RadialNode>& nodes,
                          long double& estimated_error) {
    if (!(upper > lower)) {
        return;
    }
    const RadialPanel coarse = gauss_panel(bounds, lower, upper, evaluator);
    const long double middle = lower + 0.5L * (upper - lower);
    if (middle == lower || middle == upper) {
        nodes.insert(nodes.end(), coarse.nodes, coarse.nodes + 2);
        estimated_error += upper - lower;
        return;
    }
    const RadialPanel left = gauss_panel(bounds, lower, middle, evaluator);
    const RadialPanel right = gauss_panel(bounds, middle, upper, evaluator);
    // Four times the raw coarse/fine discrepancy is deliberately more
    // conservative than a Richardson-scaled estimate. The positive fine
    // nodes, not an independently computed integral, define the fractions.
    const long double discrepancy = 4.0L * std::max(
        std::abs(left.bound_weight + right.bound_weight - coarse.bound_weight),
        std::abs(left.unbound_weight + right.unbound_weight -
                 coarse.unbound_weight));
    constexpr long double absolute_fraction_tolerance = 1.0e-10L;
    if (discrepancy <= absolute_fraction_tolerance * (upper - lower)) {
        nodes.insert(nodes.end(), left.nodes, left.nodes + 2);
        nodes.insert(nodes.end(), right.nodes, right.nodes + 2);
        estimated_error += discrepancy;
        return;
    }
    if (depth >= 18) {
        throw std::runtime_error("threshold fraction quadrature failed to converge");
    }
    append_active_panels(bounds, lower, middle, evaluator, depth + 1,
                         nodes, estimated_error);
    append_active_panels(bounds, middle, upper, evaluator, depth + 1,
                         nodes, estimated_error);
}

ThresholdRule build_threshold_rule(const PhaseSpaceGrid& grid,
                                   std::size_t ir, std::size_t iv,
                                   const EscapeSpeedEvaluator& evaluator) {
    if (ir >= grid.shape()[0] || iv >= grid.shape()[1]) {
        throw std::out_of_range("threshold cell index is out of range");
    }
    if (!evaluator) {
        throw std::invalid_argument("escape-speed evaluator is empty");
    }
    const PhaseSpaceGrid::CellBounds bounds = grid.cell_bounds(ir, iv, 0);
    const double inner_escape = evaluated_escape(evaluator, bounds.r_cm[0]);
    const double outer_escape = evaluated_escape(evaluator, bounds.r_cm[1]);
    if (inner_escape < outer_escape) {
        throw std::runtime_error("escape-speed profile is not decreasing over the cell");
    }
    ThresholdRule result{};
    if (bounds.v_cm_s[1] <= outer_escape) {
        result.measure = {EscapeThresholdClass::fully_bound, 1.0, 0.0, 0.0};
        append_constant_panel(0.0L, 1.0L, 1.0,
                              result.radial_nodes);
        return result;
    }
    if (bounds.v_cm_s[0] >= inner_escape) {
        result.measure = {EscapeThresholdClass::fully_unbound, 0.0, 1.0, 0.0};
        append_constant_panel(0.0L, 1.0L, 0.0,
                              result.radial_nodes);
        return result;
    }

    const RootLocation fully_bound_end = threshold_root(
        bounds, inner_escape, outer_escape, bounds.v_cm_s[1], evaluator);
    const RootLocation fully_unbound_start = threshold_root(
        bounds, inner_escape, outer_escape, bounds.v_cm_s[0], evaluator);
    // A full-bound plateau ends at the proven bound side of its root; a
    // full-unbound plateau starts at the proven unbound side. The tiny root
    // brackets remain in the actively evaluated region, so no positive
    // plateau node can be placed on the wrong side of E=0.
    const long double full_end = fully_bound_end.lower_x_fraction;
    const long double unbound_start =
        fully_unbound_start.upper_x_fraction;
    if (full_end > unbound_start) {
        throw std::runtime_error("escape-speed profile is not decreasing over the cell");
    }
    append_constant_panel(0.0L, full_end, 1.0,
                          result.radial_nodes);
    long double estimated_error = fully_bound_end.bracket_width +
                                  fully_unbound_start.bracket_width;
    append_active_panels(bounds, full_end,
                         unbound_start, evaluator, 0,
                         result.radial_nodes, estimated_error);
    append_constant_panel(unbound_start, 1.0L, 0.0,
                          result.radial_nodes);

    long double bound_fraction = 0.0L;
    long double unbound_fraction = 0.0L;
    for (const RadialNode& node : result.radial_nodes) {
        bound_fraction += node.weight * node.bound_y_fraction;
        unbound_fraction += node.weight * node.unbound_y_fraction;
    }
    bound_fraction = std::max(0.0L, std::min(1.0L, bound_fraction));
    unbound_fraction = std::max(0.0L, std::min(1.0L, unbound_fraction));
    const double smaller = static_cast<double>(std::min(bound_fraction,
                                                         unbound_fraction));
    const double bound = bound_fraction <= unbound_fraction ?
        smaller : 1.0 - smaller;
    result.measure = {EscapeThresholdClass::threshold_crossing,
                      bound, bound_fraction <= unbound_fraction ?
                          1.0 - smaller : smaller,
                      static_cast<double>(estimated_error)};
    return result;
}

}  // namespace

CellQuadrature positive_cell_quadrature(const PhaseSpaceGrid& grid,
                                        std::size_t flat_index,
                                        CellQuadratureOrder order) {
    const PhaseSpaceGrid::CellBounds bounds = grid.cell_bounds(flat_index);
    CellQuadrature result{};
    if (order == CellQuadratureOrder::one_point) {
        result.nodes[0] = {
            measure_quantile(bounds.r_cm[0], bounds.r_cm[1], 0.5L),
            measure_quantile(bounds.v_cm_s[0], bounds.v_cm_s[1], 0.5L),
            bounds.mu[0] + 0.5 * (bounds.mu[1] - bounds.mu[0]), 1.0};
        result.size = 1;
    } else if (order == CellQuadratureOrder::two_point) {
        const long double displacement = std::sqrt(1.0L / 12.0L);
        const long double fractions[2] = {
            0.5L - displacement, 0.5L + displacement};
        for (std::size_t ir = 0; ir < 2; ++ir) {
            for (std::size_t iv = 0; iv < 2; ++iv) {
                for (std::size_t imu = 0; imu < 2; ++imu) {
                    result.nodes[4 * ir + 2 * iv + imu] = {
                        measure_quantile(bounds.r_cm[0], bounds.r_cm[1],
                                         fractions[ir]),
                        measure_quantile(bounds.v_cm_s[0], bounds.v_cm_s[1],
                                         fractions[iv]),
                        static_cast<double>(
                            static_cast<long double>(bounds.mu[0]) +
                            fractions[imu] *
                                (static_cast<long double>(bounds.mu[1]) -
                                 static_cast<long double>(bounds.mu[0]))),
                        0.125};
                }
            }
        }
        result.size = 8;
    } else {
        throw std::invalid_argument("unsupported positive cell quadrature order");
    }
    return result;
}

EscapeThresholdGeometry classify_escape_threshold_cell(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const physics::SolarBackground& background) {
    if (ir >= grid.shape()[0] || iv >= grid.shape()[1]) {
        throw std::out_of_range("threshold cell index is out of range");
    }
    const double escape_inner =
        background.escape_speed_cm_s(grid.r_faces()[ir]);
    const double escape_outer =
        background.escape_speed_cm_s(grid.r_faces()[ir + 1]);
    if (!std::isfinite(escape_inner) || !std::isfinite(escape_outer) ||
        escape_inner <= 0.0 || escape_outer <= 0.0 ||
        escape_inner < escape_outer) {
        throw std::runtime_error("escape-speed profile is not decreasing over the cell");
    }

    const double speed_lower = grid.v_faces()[iv];
    const double speed_upper = grid.v_faces()[iv + 1];
    EscapeThresholdClass classification =
        EscapeThresholdClass::threshold_crossing;
    if (speed_upper <= escape_outer) {
        classification = EscapeThresholdClass::fully_bound;
    } else if (speed_lower >= escape_inner) {
        classification = EscapeThresholdClass::fully_unbound;
    }
    return {classification, escape_inner, escape_outer};
}

EscapeThresholdMeasure threshold_cell_measure(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const EscapeSpeedEvaluator& escape_speed_cm_s) {
    return build_threshold_rule(grid, ir, iv, escape_speed_cm_s).measure;
}

EscapeThresholdMeasure threshold_cell_measure(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const physics::SolarBackground& background) {
    return threshold_cell_measure(grid, ir, iv,
                                  [&background](double radius) {
                                      return background.escape_speed_cm_s(radius);
                                  });
}

EscapeThresholdQuadrature conditional_escape_threshold_quadrature(
    const PhaseSpaceGrid& grid, std::size_t flat_index,
    EscapeThresholdSide side, const EscapeSpeedEvaluator& escape_speed_cm_s) {
    if (side != EscapeThresholdSide::bound &&
        side != EscapeThresholdSide::unbound) {
        throw std::invalid_argument("unsupported escape threshold side");
    }
    const PhaseSpaceGrid::Index index = grid.unflatten(flat_index);
    const PhaseSpaceGrid::CellBounds bounds = grid.cell_bounds(flat_index);
    const ThresholdRule rule = build_threshold_rule(grid, index[0], index[1],
                                                    escape_speed_cm_s);
    EscapeThresholdQuadrature result{rule.measure, {}};
    const long double displacement = std::sqrt(1.0L / 12.0L);
    const long double fractions[2] = {0.5L - displacement,
                                      0.5L + displacement};
    for (const RadialNode& radial : rule.radial_nodes) {
        const long double conditional_fraction =
            side == EscapeThresholdSide::bound
                ? radial.bound_y_fraction
                : radial.unbound_y_fraction;
        if (conditional_fraction <= 0.0L) {
            continue;
        }
        const double radius = measure_quantile(
            bounds.r_cm[0], bounds.r_cm[1], radial.x_fraction);
        for (long double y_fraction : fractions) {
            const long double full_y_fraction =
                side == EscapeThresholdSide::bound
                    ? conditional_fraction * y_fraction
                    : 1.0L - conditional_fraction * (1.0L - y_fraction);
            const double speed = measure_quantile(
                bounds.v_cm_s[0], bounds.v_cm_s[1], full_y_fraction);
            for (long double mu_fraction : fractions) {
                const double mu = static_cast<double>(
                    static_cast<long double>(bounds.mu[0]) +
                    mu_fraction *
                        (static_cast<long double>(bounds.mu[1]) -
                         static_cast<long double>(bounds.mu[0])));
                const double weight = static_cast<double>(
                    radial.weight * conditional_fraction * 0.25L);
                if (weight > 0.0) {
                    result.nodes.push_back({radius, speed, mu, weight});
                }
            }
        }
    }
    return result;
}

EscapeThresholdQuadrature conditional_escape_threshold_quadrature(
    const PhaseSpaceGrid& grid, std::size_t flat_index,
    EscapeThresholdSide side, const physics::SolarBackground& background) {
    return conditional_escape_threshold_quadrature(
        grid, flat_index, side,
        [&background](double radius) {
            return background.escape_speed_cm_s(radius);
        });
}

}  // namespace transport
