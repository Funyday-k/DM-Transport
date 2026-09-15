#ifndef TRANSPORT_PHASE_SPACE_PROJECTOR_HPP
#define TRANSPORT_PHASE_SPACE_PROJECTOR_HPP

#include "transport/PhaseSpaceGrid.hpp"

#include <vector>

namespace transport {

// Conservative piecewise-constant source deposition. Each point contributes
// its full particles/s weight to the containing finite-volume cell. This class
// does not apply potential/separatrix cuts, smoothing, or CIC deposition.
class PhaseSpaceProjector {
public:
    explicit PhaseSpaceProjector(const PhaseSpaceGrid& grid);

    // Add one physical source state in cm, cm/s, and dimensionless mu.
    // The weight is an integrated source rate in particles/s, not a
    // phase-space density. Invalid coordinates follow PhaseSpaceGrid's
    // out_of_range contract. Weights must be finite and nonnegative.
    // A failed deposit leaves the accumulated source unchanged.
    void deposit(double r_cm, double v_cm_s, double mu,
                 double weight_particles_s);

    void reset() noexcept;

    const std::vector<double>& source_particles_s() const noexcept;
    double total_source_particles_s() const noexcept;

private:
    // Own a geometry copy so the source layout cannot be invalidated by a
    // caller reassigning or destroying its grid after construction.
    PhaseSpaceGrid grid_;
    std::vector<double> source_particles_s_;
    double total_source_particles_s_;
};

}  // namespace transport

#endif
