#ifndef TRANSPORT_PHYSICS_SOLAR_BACKGROUND_HPP
#define TRANSPORT_PHYSICS_SOLAR_BACKGROUND_HPP

#include <cstddef>
#include <memory>
#include <string>

namespace transport {
namespace physics {

struct SolarTarget {
    std::string name;
    unsigned int atomic_number;
    unsigned int mass_number;
    double spin;
    double proton_spin;
    double neutron_spin;
    double mass_GeV;
};

// Project-local AGSS09 solar background in explicit cgs units. Construction
// reads project-owned copies of the solar model and nuclear target tables; it
// never locates or reads a reference checkout. Port details are recorded in
// Code/provenance/reference_physics.json.
class SolarBackground {
public:
    SolarBackground(const std::string& solar_model_path,
                    const std::string& nuclear_data_path);
    ~SolarBackground();

    SolarBackground(SolarBackground&&) noexcept;
    SolarBackground& operator=(SolarBackground&&) noexcept;
    SolarBackground(const SolarBackground&) = delete;
    SolarBackground& operator=(const SolarBackground&) = delete;

    double solar_radius_cm() const noexcept;

    std::size_t target_count() const noexcept;
    const SolarTarget& target(std::size_t target_index) const;

    // All radii are in cm and must be finite and nonnegative. Temperature is
    // defined only through the photosphere. Enclosed mass becomes the total
    // solar mass outside it; target number density becomes zero.
    // Escape speed uses the exterior point-mass solution.
    double temperature_K(double radius_cm) const;
    double mass_enclosed_g(double radius_cm) const;
    double escape_speed_cm_s(double radius_cm) const;
    double number_density_cm3(std::size_t target_index,
                              double radius_cm) const;

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

}  // namespace physics
}  // namespace transport

#endif
