#ifndef TRANSPORT_COLLISION_KERNEL_HPP
#define TRANSPORT_COLLISION_KERNEL_HPP

#include "transport/PhaseSpaceCellGeometry.hpp"
#include "transport/physics/ScatteringPhysics.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace transport {

// Q[source,destination] is a rate in s^-1; occupations are column vectors
// advanced by Q^T. A collision staying in its source cell is an event but
// cancels from the generator. Velocity-domain overflow is an unresolved
// numerical channel, never a physical escape sink.
struct TransitionEntry {
    std::size_t destination;
    double rate_s_inv;
};

struct TransitionCount {
    std::size_t destination;
    std::uint64_t count;
};

// One positive FV quadrature node's complete, unweighted event histogram.
// Counts may arrive in any order; assembly sorts and combines destinations.
// A zero-rate node has sample_count=0 and no events. A positive-rate node
// must have the requested samples_per_node events, including self events.
struct CollisionNodeHistogram {
    double quadrature_weight;
    double event_rate_s_inv;
    std::uint64_t sample_count;
    std::vector<TransitionCount> cell_counts;
    std::uint64_t numerical_velocity_overflow_count;
};

struct CollisionKernelRow {
    std::size_t source;
    double mean_event_rate_s_inv;
    std::vector<TransitionEntry> off_diagonal;  // Sorted by destination.
    double diagonal_s_inv;
    double self_event_rate_s_inv;  // Diagnostic; cancels from Q.
    double numerical_velocity_overflow_rate_s_inv;
    double numerical_velocity_overflow_probability;
    std::uint64_t samples_per_node;
};

// Deterministic count-to-rate assembly shared by the reference sampler and
// future distributed count reductions. Each node contributes w_q Gamma_q
// times its own event frequencies: averaging Gamma and P separately is wrong.
// The extended row ledger includes numerical overflow, so
// sum_beta Q[source,beta] + overflow_rate = 0 to floating-point precision.
CollisionKernelRow assemble_collision_kernel_row(
    const PhaseSpaceGrid& grid,
    std::size_t source,
    const std::vector<CollisionNodeHistogram>& nodes,
    std::uint64_t samples_per_node);

// Plain reference Monte Carlo: one independent std::mt19937 per logical
// (model, grid, source, quadrature node, sample, stage) tuple. The seed uses
// a specified integer/IEEE-754 mixing scheme, not std::hash, MPI rank or
// worker scheduling. Floating and library sampling sequences are only
// promised to repeat with the same platform, compiler and C++ library.
// The builder keeps the source radial index because an interior r quadrature
// node may round onto the upper face of an ultra-thin radial cell.
CollisionKernelRow build_reference_collision_row(
    const PhaseSpaceGrid& grid,
    const physics::SolarBackground& background,
    const physics::SdProtonModel& model,
    std::size_t source,
    CellQuadratureOrder quadrature_order,
    std::uint64_t samples_per_node,
    std::uint64_t master_seed);

}  // namespace transport

#endif
