#ifndef TRANSPORT_GRID_SCHEMA_HPP
#define TRANSPORT_GRID_SCHEMA_HPP

#include "transport/PhaseSpaceGrid.hpp"

#include <iosfwd>

namespace transport {

// Minimal versioned metadata for exchanging one finite-volume grid with
// Python. Faces are inlined as JSON numbers, so this is not a bulk-state
// storage format. The quadrature order counts nodes per measure coordinate.
// Throws invalid_argument unless the order is one or two.
void write_phase_space_grid_schema_json(std::ostream& output,
                                        const PhaseSpaceGrid& grid,
                                        int quadrature_order_per_axis);

}  // namespace transport

#endif
