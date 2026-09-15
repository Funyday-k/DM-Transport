#ifndef TRANSPORT_PHASE_SPACE_GRID_HPP
#define TRANSPORT_PHASE_SPACE_GRID_HPP

#include <array>
#include <cstddef>
#include <vector>

namespace transport {

// Spherical (r, v, mu) finite-volume geometry with measure
// dGamma = 8*pi^2*r^2*v^2 dr dv dmu; see Proposal section 12.
// Faces use cm, cm/s, and dimensionless mu, respectively.
class PhaseSpaceGrid {
public:
    using Index = std::array<std::size_t, 3>;

    // Every axis needs at least two finite, strictly increasing faces.
    // r and v must be nonnegative; mu must span exactly [-1, 1].
    // Throws invalid_argument for invalid faces, overflow_error for an
    // unrepresentable cell count/volume, and underflow_error if a positive
    // cell volume rounds to zero. Positive subnormal volumes are accepted.
    PhaseSpaceGrid(std::vector<double> r_faces_cm,
                   std::vector<double> v_faces_cm_s,
                   std::vector<double> mu_faces);

    const Index& shape() const noexcept;
    std::size_t size() const noexcept;

    // mu varies fastest, then v, then r. Indices are zero-based; the upper
    // boundary of each axis is excluded. Invalid indices throw out_of_range.
    std::size_t flatten(std::size_t ir, std::size_t iv,
                        std::size_t imu) const;
    Index unflatten(std::size_t flat_index) const;

    // Integral of dGamma over the full cell, in cm^6/s^3.
    double cell_volume(std::size_t ir, std::size_t iv,
                       std::size_t imu) const;
    double cell_volume(std::size_t flat_index) const;

private:
    double volume_unchecked(std::size_t ir, std::size_t iv,
                            std::size_t imu) const;

    std::vector<double> r_faces_cm_;
    std::vector<double> v_faces_cm_s_;
    std::vector<double> mu_faces_;
    Index shape_;
    std::size_t size_;
};

}  // namespace transport

#endif
