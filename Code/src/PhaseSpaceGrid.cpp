#include "transport/PhaseSpaceGrid.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace transport {
namespace {

void validate_faces(const std::vector<double>& faces, const char* name,
                    bool nonnegative) {
    if (faces.size() < 2) {
        throw std::invalid_argument(std::string(name) +
                                    " needs at least two faces");
    }
    for (std::size_t i = 0; i < faces.size(); ++i) {
        if (!std::isfinite(faces[i])) {
            throw std::invalid_argument(std::string(name) +
                                        " faces must be finite");
        }
        if (nonnegative && faces[i] < 0.0) {
            throw std::invalid_argument(std::string(name) +
                                        " faces must be nonnegative");
        }
        if (i > 0 && !(faces[i] > faces[i - 1])) {
            throw std::invalid_argument(std::string(name) +
                                        " faces must strictly increase");
        }
    }
}

std::size_t checked_count_product(std::size_t a, std::size_t b) {
    if (a > std::numeric_limits<std::size_t>::max() / b) {
        throw std::overflow_error("phase-space cell count overflows size_t");
    }
    return a * b;
}

}  // namespace

PhaseSpaceGrid::PhaseSpaceGrid(std::vector<double> r_faces_cm,
                               std::vector<double> v_faces_cm_s,
                               std::vector<double> mu_faces)
    : r_faces_cm_(std::move(r_faces_cm)),
      v_faces_cm_s_(std::move(v_faces_cm_s)),
      mu_faces_(std::move(mu_faces)),
      shape_{{0, 0, 0}},
      size_(0) {
    validate_faces(r_faces_cm_, "r_cm", true);
    validate_faces(v_faces_cm_s_, "v_cm_s", true);
    validate_faces(mu_faces_, "mu", false);
    // These endpoints are exactly representable and define the domain,
    // rather than being approximate numerical comparison results.
    if (mu_faces_.front() != -1.0 || mu_faces_.back() != 1.0) {
        throw std::invalid_argument("mu faces must span exactly [-1, 1]");
    }
    shape_ = {{r_faces_cm_.size() - 1, v_faces_cm_s_.size() - 1,
               mu_faces_.size() - 1}};
    size_ = checked_count_product(
        checked_count_product(shape_[0], shape_[1]), shape_[2]);

    // Reject unusable geometry at construction; retain only the faces,
    // not a second three-dimensional array of cell volumes.
    for (std::size_t ir = 0; ir < shape_[0]; ++ir) {
        for (std::size_t iv = 0; iv < shape_[1]; ++iv) {
            for (std::size_t imu = 0; imu < shape_[2]; ++imu) {
                volume_unchecked(ir, iv, imu);
            }
        }
    }
}

const PhaseSpaceGrid::Index& PhaseSpaceGrid::shape() const noexcept {
    return shape_;
}

std::size_t PhaseSpaceGrid::size() const noexcept {
    return size_;
}

std::size_t PhaseSpaceGrid::flatten(std::size_t ir, std::size_t iv,
                                   std::size_t imu) const {
    if (ir >= shape_[0] || iv >= shape_[1] || imu >= shape_[2]) {
        throw std::out_of_range("phase-space cell index is out of range");
    }
    return (ir * shape_[1] + iv) * shape_[2] + imu;
}

PhaseSpaceGrid::Index PhaseSpaceGrid::unflatten(std::size_t flat_index) const {
    if (flat_index >= size_) {
        throw std::out_of_range("flat phase-space index is out of range");
    }
    const std::size_t imu = flat_index % shape_[2];
    const std::size_t radial_velocity_index = flat_index / shape_[2];
    return {{radial_velocity_index / shape_[1],
             radial_velocity_index % shape_[1], imu}};
}

double PhaseSpaceGrid::cell_volume(std::size_t ir, std::size_t iv,
                                  std::size_t imu) const {
    flatten(ir, iv, imu);
    return volume_unchecked(ir, iv, imu);
}

double PhaseSpaceGrid::cell_volume(std::size_t flat_index) const {
    const Index index = unflatten(flat_index);
    return volume_unchecked(index[0], index[1], index[2]);
}

double PhaseSpaceGrid::volume_unchecked(std::size_t ir, std::size_t iv,
                                       std::size_t imu) const {
    const double r_lo = r_faces_cm_[ir];
    const double r_hi = r_faces_cm_[ir + 1];
    const double v_lo = v_faces_cm_s_[iv];
    const double v_hi = v_faces_cm_s_[iv + 1];
    const double r_ratio = r_lo / r_hi;
    const double v_ratio = v_lo / v_hi;
    const double pi = std::acos(-1.0);

    // b^3-a^3 = (b-a)*b*b*(1+a/b+(a/b)^2), avoiding cancellation
    // of nearly equal cubes. Track binary exponents separately so a
    // representable final volume does not fail on an intermediate cube.
    const std::array<double, 10> factors{{
        8.0 * pi * pi / 9.0,
        r_hi - r_lo, r_hi, r_hi, 1.0 + r_ratio + r_ratio * r_ratio,
        v_hi - v_lo, v_hi, v_hi, 1.0 + v_ratio + v_ratio * v_ratio,
        mu_faces_[imu + 1] - mu_faces_[imu]
    }};
    double mantissa = 1.0;
    int exponent = 0;
    for (double factor : factors) {
        int factor_exponent = 0;
        mantissa *= std::frexp(factor, &factor_exponent);
        exponent += factor_exponent;
        int adjustment = 0;
        mantissa = std::frexp(mantissa, &adjustment);
        exponent += adjustment;
    }
    if (exponent > std::numeric_limits<double>::max_exponent) {
        throw std::overflow_error("phase-space cell volume overflows double");
    }
    const double volume = std::ldexp(mantissa, exponent);
    if (!std::isfinite(volume)) {
        throw std::overflow_error("phase-space cell volume overflows double");
    }
    if (volume == 0.0) {
        throw std::underflow_error("phase-space cell volume rounds to zero");
    }
    return volume;
}

}  // namespace transport
