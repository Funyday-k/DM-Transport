#include "transport/PhaseSpaceProjector.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace transport {

PhaseSpaceProjector::PhaseSpaceProjector(const PhaseSpaceGrid& grid)
    : grid_(grid),
      source_particles_s_(grid_.size(), 0.0),
      total_source_particles_s_(0.0) {}

void PhaseSpaceProjector::deposit(double r_cm, double v_cm_s, double mu,
                                  double weight_particles_s) {
    if (!std::isfinite(weight_particles_s) || weight_particles_s < 0.0) {
        throw std::invalid_argument(
            "source weight must be finite and nonnegative");
    }

    // Finish every operation that can fail before committing either sum.
    const PhaseSpaceGrid::Index index = grid_.locate_cell(r_cm, v_cm_s, mu);
    const std::size_t flat = grid_.flatten(index[0], index[1], index[2]);
    const double next_cell_source =
        source_particles_s_[flat] + weight_particles_s;
    const double next_total_source =
        total_source_particles_s_ + weight_particles_s;
    if (!std::isfinite(next_cell_source) ||
        !std::isfinite(next_total_source)) {
        throw std::overflow_error("source accumulation overflows double");
    }

    source_particles_s_[flat] = next_cell_source;
    total_source_particles_s_ = next_total_source;
}

void PhaseSpaceProjector::reset() noexcept {
    std::fill(source_particles_s_.begin(), source_particles_s_.end(), 0.0);
    total_source_particles_s_ = 0.0;
}

const std::vector<double>&
PhaseSpaceProjector::source_particles_s() const noexcept {
    return source_particles_s_;
}

double PhaseSpaceProjector::total_source_particles_s() const noexcept {
    return total_source_particles_s_;
}

}  // namespace transport
