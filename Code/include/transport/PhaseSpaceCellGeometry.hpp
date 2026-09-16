#ifndef TRANSPORT_PHASE_SPACE_CELL_GEOMETRY_HPP
#define TRANSPORT_PHASE_SPACE_CELL_GEOMETRY_HPP

#include "transport/PhaseSpaceGrid.hpp"
#include "transport/physics/SolarBackground.hpp"

#include <array>
#include <cstddef>

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

}  // namespace transport

#endif
