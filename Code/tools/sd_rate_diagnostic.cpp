#include "transport/physics/ScatteringPhysics.hpp"

#include <cmath>
#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

double parse_number(const char* text, const std::string& label) {
    std::size_t parsed = 0;
    const std::string value(text);
    double result = 0.0;
    try {
        result = std::stod(value, &parsed);
    } catch (const std::exception&) {
        throw std::invalid_argument(label + " must be a number");
    }
    if (parsed != value.size() || !std::isfinite(result)) {
        throw std::invalid_argument(label + " must be a finite number");
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 7) {
            throw std::invalid_argument(
                "usage: sd_rate_diagnostic SOLAR_MODEL NUCLEAR_DATA "
                "RADIUS_CM SPEED_CM_S DM_MASS_GEV SIGMA_P_CM2");
        }
        const transport::physics::SolarBackground background(argv[1], argv[2]);
        const double radius_cm = parse_number(argv[3], "radius");
        const double speed_cm_s = parse_number(argv[4], "speed");
        const transport::physics::SdProtonModel model{
            parse_number(argv[5], "dark-matter mass"),
            parse_number(argv[6], "proton cross section")
        };
        const auto rates =
            transport::physics::direct_sd_proton_scattering_rates(
                background, model, radius_cm, speed_cm_s);

        std::cout << std::setprecision(17)
                  << "# producer=DM-Transport/sd_rate_diagnostic\n"
                  << "# qualification=diagnostic_not_independent_oracle\n"
                  << "# interaction=SD_constant_contact_proton_only\n"
                  << "# radius_cm=" << radius_cm << '\n'
                  << "# dm_speed_cm_s=" << speed_cm_s << '\n'
                  << "# dm_mass_GeV=" << model.dark_matter_mass_GeV << '\n'
                  << "# sigma_p_cm2=" << model.proton_cross_section_cm2 << '\n'
                  << "# total_rate_s_inv=" << rates.total_rate_s_inv << '\n'
                  << "target_index\tisotope\tZ\tA\tspin\tproton_spin\t"
                     "neutron_spin\tnumber_density_cm3\t"
                     "nucleus_cross_section_cm2\tmean_relative_speed_cm_s\t"
                     "rate_s_inv\tfraction\n";
        for (const auto& entry : rates.target_rates) {
            const auto& target = background.target(entry.target_index);
            const double fraction = rates.total_rate_s_inv == 0.0
                ? 0.0
                : entry.rate_s_inv / rates.total_rate_s_inv;
            std::cout << entry.target_index << '\t'
                      << target.name << '\t'
                      << target.atomic_number << '\t'
                      << target.mass_number << '\t'
                      << target.spin << '\t'
                      << target.proton_spin << '\t'
                      << target.neutron_spin << '\t'
                      << background.number_density_cm3(
                             entry.target_index, radius_cm) << '\t'
                      << entry.nucleus_cross_section_cm2 << '\t'
                      << entry.mean_relative_speed_cm_s << '\t'
                      << entry.rate_s_inv << '\t'
                      << fraction << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "SD rate diagnostic failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
