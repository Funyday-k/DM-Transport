#include "transport/GridSchema.hpp"

#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char* argv[]) {
    try {
        if (argc != 2 || (std::strcmp(argv[1], "--schema") != 0 &&
                          std::strcmp(argv[1], "--schema-fixed") != 0 &&
                          std::strcmp(argv[1], "--probe") != 0)) {
            throw std::invalid_argument(
                "usage: t04_grid_schema_probe --schema|--schema-fixed|--probe");
        }
        // Unequal widths on every axis are deliberate cross-language fixtures.
        const transport::PhaseSpaceGrid grid(
            {0.0, 2.0e10, 3.1e10, 6.95e10},
            {0.0, 1.0e6, 4.0e7, 1.4e8},
            {-1.0, -0.37, 0.2, 1.0});
        if (std::strcmp(argv[1], "--schema") == 0 ||
            std::strcmp(argv[1], "--schema-fixed") == 0) {
            if (std::strcmp(argv[1], "--schema-fixed") == 0) {
                std::cout << std::fixed << std::showpos << std::setprecision(2);
            }
            const auto original_flags = std::cout.flags();
            const auto original_precision = std::cout.precision();
            transport::write_phase_space_grid_schema_json(std::cout, grid, 2);
            if (std::cout.flags() != original_flags ||
                std::cout.precision() != original_precision) {
                throw std::runtime_error("schema writer changed caller stream format");
            }
            std::cout << '\n';
            return 0;
        }
        std::cout << "{\"schema\":";
        transport::write_phase_space_grid_schema_json(std::cout, grid, 2);
        std::cout << ",\"cpp_cells\":[";
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (std::size_t ir = 0; ir < grid.shape()[0]; ++ir) {
            for (std::size_t iv = 0; iv < grid.shape()[1]; ++iv) {
                for (std::size_t imu = 0; imu < grid.shape()[2]; ++imu) {
                    const std::size_t flat = grid.flatten(ir, iv, imu);
                    if (flat != 0) {
                        std::cout << ',';
                    }
                    std::cout << "{\"index\":[" << ir << ',' << iv << ',' << imu
                              << "],\"flat\":" << flat
                              << ",\"phase_measure\":"
                              << grid.cell_phase_measure(ir, iv, imu) << '}';
                }
            }
        }
        std::cout << "]}\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
