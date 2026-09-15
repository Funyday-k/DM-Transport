#include "transport/physics/SolarBackground.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using transport::physics::SolarBackground;
using transport::physics::SolarTarget;

// Horner endpoint evaluation and libm square roots have small compiler-dependent
// rounding differences; this remains far below the source table precision.
constexpr double crossPlatformRelativeTolerance = 1.0e-12;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, const std::string& message,
                   double relative_tolerance,
                   double absolute_tolerance = 0.0) {
    require(std::isfinite(actual) && std::isfinite(expected),
            message + ": nonfinite comparison");
    const double scale = std::max(std::abs(actual), std::abs(expected));
    require(std::abs(actual - expected) <=
                std::max(absolute_tolerance, relative_tolerance * scale),
            message + ": values differ");
}

template <typename Exception, typename Callable>
void require_throws(Callable operation, const std::string& message) {
    try {
        operation();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message + ": expected exception was not thrown");
}

template <typename Exception, typename Callable>
void require_throws_containing(Callable operation,
                               const std::string& expected_text,
                               const std::string& message) {
    try {
        operation();
    } catch (const Exception& error) {
        require(std::string(error.what()).find(expected_text) !=
                    std::string::npos,
                message + ": exception lacks source context");
        return;
    }
    throw std::runtime_error(message + ": expected exception was not thrown");
}

class TemporaryFile {
public:
    TemporaryFile(std::string path, const std::string& contents)
        : path_(std::move(path)) {
        std::ofstream output(path_.c_str(), std::ios::binary);
        if (!output || !(output << contents)) {
            throw std::runtime_error("unable to create loader test fixture");
        }
    }

    ~TemporaryFile() {
        std::remove(path_.c_str());
    }

    const std::string& path() const noexcept {
        return path_;
    }

private:
    std::string path_;
};

std::string minimal_solar_row() {
    std::string row = "0.0000004 0.00150 1 1";
    for (std::size_t column = 4; column < 35; ++column) {
        row += " 0";
    }
    return row + "\n";
}

void test_target_inventory(const SolarBackground& background) {
    const std::array<std::pair<unsigned int, unsigned int>, 63> expected{{
        {1, 1}, {2, 4}, {2, 3}, {6, 12}, {6, 13}, {7, 14}, {7, 15},
        {8, 16}, {8, 17}, {8, 18}, {10, 20}, {10, 21}, {10, 22},
        {11, 23}, {12, 24}, {12, 25}, {12, 26}, {13, 27}, {14, 28},
        {14, 29}, {14, 30}, {15, 31}, {16, 32}, {16, 33}, {16, 34},
        {16, 36}, {17, 35}, {17, 37}, {18, 36}, {18, 38}, {18, 40},
        {19, 39}, {19, 40}, {19, 41}, {20, 40}, {20, 42}, {20, 43},
        {20, 44}, {20, 46}, {20, 48}, {21, 45}, {22, 46}, {22, 47},
        {22, 48}, {22, 49}, {22, 50}, {23, 50}, {23, 51}, {24, 50},
        {24, 52}, {24, 53}, {24, 54}, {25, 55}, {26, 54}, {26, 56},
        {26, 57}, {26, 58}, {27, 59}, {28, 58}, {28, 60}, {28, 61},
        {28, 62}, {28, 64}
    }};
    require(background.target_count() == 63,
            "AGSS09 target inventory must contain 63 isotopes");
    const SolarTarget& hydrogen = background.target(0);
    require(hydrogen.name == "H-1" && hydrogen.atomic_number == 1 &&
                hydrogen.mass_number == 1 && hydrogen.spin == 0.5 &&
                hydrogen.proton_spin == 0.5 &&
                hydrogen.neutron_spin == 0.0 &&
                hydrogen.mass_GeV == 0.9382720813,
            "first target must be the explicit H-1 abundance channel");
    const SolarTarget& helium4 = background.target(1);
    const SolarTarget& helium3 = background.target(2);
    require(helium4.name == "He-4" && helium4.atomic_number == 2 &&
                helium4.mass_number == 4,
            "second target must be the explicit He-4 abundance channel");
    require(helium3.name == "He-3" && helium3.atomic_number == 2 &&
                helium3.mass_number == 3,
            "third target must be the explicit He-3 abundance channel");
    require(background.target(10).name == "Ne-20" &&
                background.target(34).name == "Ca-40",
            "natural-isotope groups must begin in reference order");
    const SolarTarget& last = background.target(62);
    require(last.name == "Ni-64" && last.atomic_number == 28 &&
                last.mass_number == 64,
            "last target must be Ni-64");

    std::set<std::string> names;
    for (std::size_t index = 0; index < background.target_count(); ++index) {
        const SolarTarget& target = background.target(index);
        require(target.atomic_number > 0 &&
                    target.mass_number >= target.atomic_number &&
                    target.spin >= 0.0 && target.mass_GeV > 0.0 &&
                    std::isfinite(target.mass_GeV),
                "solar target metadata is invalid");
        require(target.atomic_number == expected[index].first &&
                    target.mass_number == expected[index].second,
                "solar isotope inventory or order changed");
        require(names.insert(target.name).second,
                "solar target names must be unique");
    }
}

void test_reference_nodes(const SolarBackground& background) {
    constexpr double solar_radius_cm = 6.957e10;
    constexpr double solar_mass_g = 1.98848e33;
    require(background.solar_radius_cm() == solar_radius_cm,
            "legacy solar radius");

    require_close(background.mass_enclosed_g(0.0), 0.0,
                  "central enclosed mass", 0.0, 0.0);
    require_close(background.mass_enclosed_g(0.00150 * solar_radius_cm),
                  0.0000004 * solar_mass_g,
                  "first AGSS09 enclosed-mass node", 2.0e-15);
    require_close(background.mass_enclosed_g(solar_radius_cm), solar_mass_g,
                  "surface enclosed mass", crossPlatformRelativeTolerance);

    require_close(background.temperature_K(0.0), 1.549e7,
                  "central temperature", 2.0e-15);
    require_close(background.temperature_K(0.96950 * solar_radius_cm),
                  1.544e5, "AGSS09 temperature node", 2.0e-15);
    require_close(background.temperature_K(solar_radius_cm), 5800.0,
                  "photosphere temperature", crossPlatformRelativeTolerance);
    require_close(background.number_density_cm3(0, 0.0),
                  3.2580314157719068e25,
                  "central H-1 number density", 2.0e-14);
    require_close(background.number_density_cm3(
                      0, 0.001 * solar_radius_cm),
                  3.2580314157719068e25,
                  "central plateau H-1 number density", 2.0e-14);
    require_close(background.number_density_cm3(0, solar_radius_cm),
                  4.515067196525714e14,
                  "photosphere H-1 number density", 1.0e-9);
    require_close(background.number_density_cm3(
                      0, 0.50025 * solar_radius_cm),
                  5.8042676968996363e23,
                  "off-node H-1 number-density interpolation",
                  crossPlatformRelativeTolerance);
    // The source nuclear table's Ca fractions sum to 1.00003.  This selected
    // density changes if those fractions are silently normalized.
    require_close(background.number_density_cm3(34, 0.0),
                  1.4598770716214873e20,
                  "central Ca-40 number density", 2.0e-14);
}

void test_escape_speed_and_exterior(const SolarBackground& background) {
    constexpr double solar_radius_cm = 6.957e10;
    constexpr double solar_mass_g = 1.98848e33;
    require_close(background.escape_speed_cm_s(0.0),
                  1.384130645174225e8,
                  "central escape speed", 0.0, 1.0e4);
    require_close(background.escape_speed_cm_s(solar_radius_cm),
                  6.176755830414014e7,
                  "surface escape speed", crossPlatformRelativeTolerance);
    require_close(background.escape_speed_cm_s(2.0 * solar_radius_cm),
                  4.3676259334192932e7,
                  "exterior point-mass escape speed",
                  crossPlatformRelativeTolerance);
    require(background.mass_enclosed_g(2.0 * solar_radius_cm) == solar_mass_g,
            "exterior enclosed mass must equal the solar mass");
    require(background.number_density_cm3(
                0, 2.0 * solar_radius_cm) == 0.0,
            "exterior target density must vanish");
}

void test_load_failures(const std::string& solar_model_path,
                        const std::string& nuclear_data_path,
                        const std::string& temporary_directory) {
    require_throws<std::runtime_error>([&] {
        SolarBackground missing_solar(solar_model_path + ".missing",
                                      nuclear_data_path);
        (void)missing_solar;
    }, "missing solar model");
    require_throws<std::runtime_error>([&] {
        SolarBackground missing_nuclear(solar_model_path,
                                        nuclear_data_path + ".missing");
        (void)missing_nuclear;
    }, "missing nuclear data");

    const std::string malformed_path =
        temporary_directory + "/solar_background_malformed.dat";
    const TemporaryFile malformed(
        malformed_path, minimal_solar_row() + "0.0000009 broken\n");
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(malformed.path(), nuclear_data_path);
        (void)invalid;
    }, malformed.path() + ":2:", "malformed solar row");

    const std::string extra_path =
        temporary_directory + "/solar_background_extra_token.dat";
    std::string extra_contents = minimal_solar_row();
    extra_contents.insert(extra_contents.size() - 1, " trailing");
    const TemporaryFile extra(extra_path, extra_contents);
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(extra.path(), nuclear_data_path);
        (void)invalid;
    }, extra.path() + ":1:", "extra solar token");

    const std::string short_solar_path =
        temporary_directory + "/solar_background_short_solar.dat";
    const TemporaryFile short_solar(short_solar_path, minimal_solar_row());
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(short_solar.path(), nuclear_data_path);
        (void)invalid;
    }, "expected 1968 solar-model rows", "short solar table");

    const std::string duplicate_path =
        temporary_directory + "/solar_background_duplicate_nucleus.dat";
    const TemporaryFile duplicate(
        duplicate_path,
        "H 1 1 1 0.5 0.5 0\nH 1 1 1 0.5 0.5 0\n");
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(solar_model_path, duplicate.path());
        (void)invalid;
    }, duplicate.path(), "duplicate nuclear isotope");

    const std::string truncated_nuclear_path =
        temporary_directory + "/solar_background_truncated_nucleus.dat";
    const TemporaryFile truncated_nuclear(
        truncated_nuclear_path, "H 1 1 1 0.5");
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(solar_model_path, truncated_nuclear.path());
        (void)invalid;
    }, truncated_nuclear.path() + ": malformed nuclear record 1",
       "truncated nuclear record");

    const std::string short_nuclear_path =
        temporary_directory + "/solar_background_short_nucleus.dat";
    const TemporaryFile short_nuclear(
        short_nuclear_path, "H 1 1 1 0.5 0.5 0\n");
    require_throws_containing<std::runtime_error>([&] {
        SolarBackground invalid(solar_model_path, short_nuclear.path());
        (void)invalid;
    }, "expected 295 nuclear records", "short nuclear table");
}

void test_invalid_queries(const SolarBackground& background) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    require_throws<std::invalid_argument>([&] {
        background.temperature_K(-1.0);
    }, "negative radius");
    require_throws<std::invalid_argument>([&] {
        background.escape_speed_cm_s(nan);
    }, "nonfinite radius");
    require_throws<std::out_of_range>([&] {
        background.temperature_K(2.0 * background.solar_radius_cm());
    }, "exterior temperature");
    require_throws<std::out_of_range>([&] {
        background.number_density_cm3(background.target_count(), 0.0);
    }, "invalid target index");
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "usage: test_solar_background MODEL_TABLE NUCLEAR_DATA "
                     "TEMPORARY_DIRECTORY\n";
        return EXIT_FAILURE;
    }
    try {
        const SolarBackground background(argv[1], argv[2]);
        test_target_inventory(background);
        test_reference_nodes(background);
        test_escape_speed_and_exterior(background);
        test_invalid_queries(background);
        test_load_failures(argv[1], argv[2], argv[3]);
        std::cout << "solar background tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "solar background test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
