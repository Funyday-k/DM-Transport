#include "transport/PhaseSpaceCellGeometry.hpp"

#include <cmath>
#include <stdexcept>

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

}  // namespace transport
