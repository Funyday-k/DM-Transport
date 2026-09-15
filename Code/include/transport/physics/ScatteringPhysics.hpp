#ifndef TRANSPORT_PHYSICS_SCATTERING_PHYSICS_HPP
#define TRANSPORT_PHYSICS_SCATTERING_PHYSICS_HPP

namespace transport {
namespace physics {

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

}  // namespace physics
}  // namespace transport

#endif
