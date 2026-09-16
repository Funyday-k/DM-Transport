#ifndef TRANSPORT_PHYSICS_SCATTERING_PHYSICS_HPP
#define TRANSPORT_PHYSICS_SCATTERING_PHYSICS_HPP

#include <array>
#include <cstddef>
#include <random>
#include <vector>

#include "transport/physics/SolarBackground.hpp"

namespace transport {
namespace physics {

using CartesianVelocityCmS = std::array<double, 3>;

// Port provenance: Code/provenance/reference_physics.json, unit mean_relative_speed.
// Mean speed of a particle moving at dm_speed_cm_s relative to a
// Maxwell-Boltzmann target population. Temperature is in kelvin, target mass
// is in GeV, and the result is in cm/s.
//
// The Maxwell-Boltzmann model is nonrelativistic and is intended for
// k_B T / (m c^2) << 1. This function does not impose a separate
// application-specific cutoff on that ratio.
//
// All inputs must be finite. Temperature and target mass must be strictly
// positive with a finite, positive floating-point ratio, and the dark-matter
// speed must be nonnegative. Invalid inputs throw std::invalid_argument; an
// unrepresentable result throws std::overflow_error.
double mean_relative_speed_cm_s(double temperature_K,
                                double target_mass_GeV,
                                double dm_speed_cm_s);

// Frozen MVP interaction: m_chi = 0.1 GeV constant-contact spin-dependent
// scattering with proton-only coupling.  The cross section is the reference
// proton cross section in cm^2; electrons and neutron coupling are outside
// this type.  Other dark-matter masses are rejected until separately ported
// and validated.
struct SdProtonModel {
    double dark_matter_mass_GeV;
    double proton_cross_section_cm2;
};

struct SdTargetScatteringRate {
    std::size_t target_index;
    double nucleus_cross_section_cm2;
    double mean_relative_speed_cm_s;
    double rate_s_inv;
};

struct SdScatteringRates {
    std::vector<SdTargetScatteringRate> target_rates;
    double total_rate_s_inv;
};

struct CollisionSample {
    std::size_t target_index;
    CartesianVelocityCmS target_velocity_cm_s;
    CartesianVelocityCmS outgoing_dm_velocity_cm_s;
};

// Port of the fixed obscura SD normalization for a_n = 0.  Spin-zero targets
// have zero cross section.  Invalid model or target metadata throws
// std::invalid_argument; an unrepresentable result throws std::overflow_error.
double sd_proton_nucleus_cross_section_cm2(const SdProtonModel& model,
                                           const SolarTarget& target);

// Direct (non-tabulated) solar scattering rate in s^-1.  Every SolarBackground
// isotope remains present and in source order, including zero-rate targets.
// Radius and speed use cm and cm/s.  Exterior rates are zero because the
// background target densities vanish there.  Invalid inputs throw
// std::invalid_argument; unrepresentable rates throw std::overflow_error.
SdScatteringRates direct_sd_proton_scattering_rates(
    const SolarBackground& background,
    const SdProtonModel& model,
    double radius_cm,
    double dm_speed_cm_s);

// Draw one nuclear target from the complete source-order rate breakdown.
// The input must contain a finite, positive total equal to the source-order
// sum of finite, nonnegative target rates.  The fixed MVP has no electron
// target.  Sampling requests exactly one uniform variate from rng.
std::size_t sample_sd_proton_target_index(
    const SdScatteringRates& rates,
    std::mt19937& rng);

// Collision-conditioned thermal target velocity for a constant cross
// section.  This samples f_MB(u) |v_chi-u| rather than the unconditioned
// Maxwell-Boltzmann bath.  Temperature, mass, and all vector components use
// kelvin, GeV, and cm/s.  The legacy sampler requires a nonzero incoming
// velocity; its zero-speed analytic limit is validated separately.
CartesianVelocityCmS sample_collision_conditioned_target_velocity_cm_s(
    double temperature_K,
    double target_mass_GeV,
    const CartesianVelocityCmS& dm_velocity_cm_s,
    std::mt19937& rng);

// Nonrelativistic two-body elastic kinematics in explicit cgs velocity units.
// outgoing_dm_cm_direction_unit is the unit direction of the outgoing dark
// matter velocity relative to the center of mass.  This deterministic helper
// is exposed so conservation laws and fixed-angle limits can be tested without
// coupling them to random sampling.  Inputs must be finite, masses positive,
// and the direction normalized to within floating-point tolerance.
CartesianVelocityCmS elastic_outgoing_dm_velocity_cm_s(
    double dark_matter_mass_GeV,
    double target_mass_GeV,
    const CartesianVelocityCmS& incoming_dm_velocity_cm_s,
    const CartesianVelocityCmS& incoming_target_velocity_cm_s,
    const CartesianVelocityCmS& outgoing_dm_cm_direction_unit);

// Complete local collision for the fixed SD-proton, low-mass contact MVP:
// direct 63-target rates -> source-order target -> collision-conditioned
// thermal target velocity -> isotropic CM direction -> elastic outgoing DM
// velocity.  The routine contains no trajectory, binning, binding-energy,
// escape, or kernel logic.  Radius is in cm and velocities are in cm/s.
CollisionSample sample_sd_proton_collision(
    const SolarBackground& background,
    const SdProtonModel& model,
    double radius_cm,
    const CartesianVelocityCmS& incoming_dm_velocity_cm_s,
    std::mt19937& rng);

}  // namespace physics
}  // namespace transport

#endif
