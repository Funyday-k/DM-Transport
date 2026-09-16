#include "transport/GridSchema.hpp"

#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace transport {
namespace {

void write_faces(std::ostream& output, const std::vector<double>& faces) {
    output << '[';
    for (std::size_t i = 0; i < faces.size(); ++i) {
        if (i != 0) {
            output << ',';
        }
        output << faces[i];
    }
    output << ']';
}

}  // namespace

void write_phase_space_grid_schema_json(std::ostream& output,
                                        const PhaseSpaceGrid& grid,
                                        int quadrature_order_per_axis) {
    if (quadrature_order_per_axis != 1 && quadrature_order_per_axis != 2) {
        throw std::invalid_argument("quadrature order per axis must be 1 or 2");
    }
    // A max_digits10 decimal round-trips every binary64 face in Python.
    const std::streamsize original_precision = output.precision();
    const std::ios::fmtflags original_flags = output.flags();
    // Clear inherited flags such as hexfloat and showpos, either of which
    // would make numeric tokens invalid JSON.
    output.flags(std::ios::dec);
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    output << "{\"schema_version\":1,\"r_faces_cm\":";
    write_faces(output, grid.r_faces());
    output << ",\"v_faces_cm_s\":";
    write_faces(output, grid.v_faces());
    output << ",\"mu_faces\":";
    write_faces(output, grid.mu_faces());
    output << ",\"flatten_order\":[\"r\",\"v\",\"mu\"]"
           << ",\"units\":{\"r\":\"cm\",\"v\":\"cm/s\",\"mu\":\"1\","
           << "\"cell_phase_measure\":\"cm^6/s^3\"}"
           << ",\"quadrature_order_per_axis\":" << quadrature_order_per_axis
           << ",\"escape_threshold\":{"
           << "\"energy_definition\":\"0.5*(v_cm_s^2-v_escape_cm_s(r_cm)^2)\","
           << "\"bound_condition\":\"E<0\","
           << "\"unbound_condition\":\"E>0\","
           << "\"zero_surface\":\"zero_phase_measure\","
           << "\"touching_cell\":\"classify_by_interior\","
           << "\"point_deposition\":\"containing_fv_cell_no_fractional_split\"}}";
    output.flags(original_flags);
    output.precision(original_precision);
}

}  // namespace transport
