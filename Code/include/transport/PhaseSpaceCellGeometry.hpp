#ifndef TRANSPORT_PHASE_SPACE_CELL_GEOMETRY_HPP
#define TRANSPORT_PHASE_SPACE_CELL_GEOMETRY_HPP

#include "transport/PhaseSpaceGrid.hpp"
#include "transport/physics/SolarBackground.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace transport {

// Nodes integrate the normalized FV measure r^2 v^2 dr dv dmu within a cell.
// The one-point rule uses the midpoint in r^3, v^3, and mu; the two-point
// rule uses two Gauss-Legendre nodes on each of those measure coordinates.
// This is a quadrature rule,
// not a replacement for the finite-volume faces or output deposition rule.
enum class CellQuadratureOrder { one_point, two_point };

struct CellQuadratureNode {
    double r_cm;
    double v_cm_s;
    double mu;
    double weight;  // Nonnegative fraction of the full cell phase measure.
};

struct CellQuadrature {
    std::array<CellQuadratureNode, 8> nodes;
    std::size_t size;
};

CellQuadrature positive_cell_quadrature(const PhaseSpaceGrid& grid,
                                        std::size_t flat_index,
                                        CellQuadratureOrder order);

enum class EscapeThresholdClass {
    fully_bound,
    fully_unbound,
    threshold_crossing
};

struct EscapeThresholdGeometry {
    EscapeThresholdClass classification;
    double inner_escape_speed_cm_s;
    double outer_escape_speed_cm_s;
};

// E = (v^2 - v_escape(r)^2)/2. The classification applies to the entire
// (r,v) cell, independent of mu. A cell touching E=0 only on its boundary
// is assigned to the side occupied by its interior. This relies on the
// decreasing radial escape-speed profile supplied by SolarBackground.
EscapeThresholdGeometry classify_escape_threshold_cell(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const physics::SolarBackground& background);

// The FV measure is uniform in x=r^3, y=v^3 and mu. For a crossing cell,
// bound_fraction is the area under y_escape(x), clipped to the cell's y
// interval, divided by the full x-y rectangle area. Fractions describe a
// cell-integrated occupation; a point collision output must be deposited by
// its actual (r,v,mu) coordinates, never split by these fractions.
struct EscapeThresholdMeasure {
    EscapeThresholdClass classification;
    double bound_fraction;
    double unbound_fraction;
    // Sum of conservative coarse/refined Gauss-rule differences, including
    // root-bracket uncertainty. This is a convergence indicator, not a
    // rigorous mathematical error bound for an arbitrary escape profile.
    double estimated_absolute_error;
};

// The injected evaluator supports analytic test profiles. Its escape speed
// must be finite, positive, continuous, and nonincreasing across the radial
// cell; the adaptive convergence indicator assumes a smooth profile.
using EscapeSpeedEvaluator = std::function<double(double)>;

EscapeThresholdMeasure threshold_cell_measure(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const EscapeSpeedEvaluator& escape_speed_cm_s);
EscapeThresholdMeasure threshold_cell_measure(
    const PhaseSpaceGrid& grid, std::size_t ir, std::size_t iv,
    const physics::SolarBackground& background);

enum class EscapeThresholdSide { bound, unbound };

struct EscapeThresholdQuadrature {
    EscapeThresholdMeasure measure;
    // Positive weights are fractions of the FULL cell FV measure. Their sum
    // equals the requested side's fraction, within floating-point rounding.
    std::vector<CellQuadratureNode> nodes;
};

// On crossing cells, adaptive positive radial quadrature follows the escape
// curve in x=r^3, with conditional two-point Gauss rules in y=v^3 and mu.
// The original global cell index and its point-deposition rule are unchanged.
EscapeThresholdQuadrature conditional_escape_threshold_quadrature(
    const PhaseSpaceGrid& grid, std::size_t flat_index,
    EscapeThresholdSide side, const EscapeSpeedEvaluator& escape_speed_cm_s);
EscapeThresholdQuadrature conditional_escape_threshold_quadrature(
    const PhaseSpaceGrid& grid, std::size_t flat_index,
    EscapeThresholdSide side, const physics::SolarBackground& background);

}  // namespace transport

#endif
