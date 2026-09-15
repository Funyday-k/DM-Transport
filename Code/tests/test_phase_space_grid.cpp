#include "transport/PhaseSpaceGrid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using transport::PhaseSpaceGrid;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_close(double actual, double expected, const std::string& message,
                   double relative_tolerance = 3e-14) {
    require(std::isfinite(actual) && std::isfinite(expected),
            message + ": nonfinite comparison");
    const double scale = std::max(std::abs(actual), std::abs(expected));
    require(std::abs(actual - expected) <= relative_tolerance * scale,
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

void test_volume_and_indices() {
    const std::vector<double> r_cm{0.0, 0.125, 1.0, 4.0};
    const std::vector<double> v_cm_s{0.0, 0.25, 2.0};
    const std::vector<double> mu{-1.0, -0.9, -0.25, 0.5, 1.0};
    const PhaseSpaceGrid grid(r_cm, v_cm_s, mu);
    require(grid.shape() == PhaseSpaceGrid::Index{{3, 2, 4}}, "shape");
    require(grid.size() == 24, "size");

    const double pi = std::acos(-1.0);
    double total_volume = 0.0;
    std::size_t expected_flat = 0;
    for (std::size_t ir = 0; ir < 3; ++ir) {
        for (std::size_t iv = 0; iv < 2; ++iv) {
            for (std::size_t imu = 0; imu < 4; ++imu) {
                const std::size_t flat = grid.flatten(ir, iv, imu);
                require(flat == expected_flat++, "mu-fastest flatten order");
                require(grid.unflatten(flat) ==
                            PhaseSpaceGrid::Index{{ir, iv, imu}},
                        "index round trip");
                // Independent Cartesian ball/shell volumes, times the
                // fraction of the relative polar solid angle in this cell.
                const double position_shell_cm3 = 4.0 * pi / 3.0 *
                    (std::pow(r_cm[ir + 1], 3) - std::pow(r_cm[ir], 3));
                const double velocity_shell_cm3_s3 = 4.0 * pi / 3.0 *
                    (std::pow(v_cm_s[iv + 1], 3) - std::pow(v_cm_s[iv], 3));
                const double expected = position_shell_cm3 *
                    velocity_shell_cm3_s3 * (mu[imu + 1] - mu[imu]) / 2.0;
                const double volume = grid.cell_volume(ir, iv, imu);
                require(volume > 0.0, "positive cell volume");
                require_close(volume, expected, "nonuniform cell volume");
                require_close(grid.cell_volume(flat), volume,
                              "flat volume agrees");
                total_volume += volume;
            }
        }
    }
    const double position_ball_cm3 = 4.0 * pi * std::pow(4.0, 3) / 3.0;
    const double velocity_ball_cm3_s3 = 4.0 * pi * std::pow(2.0, 3) / 3.0;
    require_close(total_volume, position_ball_cm3 * velocity_ball_cm3_s3,
                  "total phase volume equals product of ball volumes");

    // The radial/velocity domain may begin above zero.
    const PhaseSpaceGrid shell({1.0, 2.0}, {2.0, 3.0}, {-1.0, 1.0});
    require_close(shell.cell_volume(0),
                  (4.0 * pi / 3.0 * 7.0) * (4.0 * pi / 3.0 * 19.0),
                  "nonzero lower domain faces");
}

void test_invalid_faces() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<std::vector<double>> invalid_nonnegative{
        {}, {0.0}, {0.0, 0.0}, {1.0, 0.0}, {-1.0, 1.0},
        {0.0, infinity}, {0.0, nan}, {nan, 1.0}
    };
    for (const auto& faces : invalid_nonnegative) {
        require_throws<std::invalid_argument>([&] {
            PhaseSpaceGrid grid(faces, {0.0, 1.0}, {-1.0, 1.0});
        }, "invalid radial faces");
        require_throws<std::invalid_argument>([&] {
            PhaseSpaceGrid grid({0.0, 1.0}, faces, {-1.0, 1.0});
        }, "invalid speed faces");
    }
    const std::vector<std::vector<double>> invalid_mu{
        {}, {-1.0}, {-1.0, 0.0}, {0.0, 1.0}, {-1.01, 1.0},
        {-1.0, 1.01}, {-1.0, 0.0, 0.0, 1.0}, {-1.0, 0.5, 0.0, 1.0},
        {-1.0, nan, 1.0}, {-infinity, 1.0}, {-1.0, infinity},
        {std::nextafter(-1.0, 0.0), 1.0}
    };
    for (const auto& faces : invalid_mu) {
        require_throws<std::invalid_argument>([&] {
            PhaseSpaceGrid grid({0.0, 1.0}, {0.0, 1.0}, faces);
        }, "invalid mu faces");
    }
}

void test_invalid_indices() {
    const PhaseSpaceGrid grid({0.0, 1.0, 2.0}, {0.0, 1.0},
                              {-1.0, 0.0, 0.5, 1.0});
    const std::size_t huge = std::numeric_limits<std::size_t>::max();
    const std::vector<PhaseSpaceGrid::Index> invalid{
        {{2, 0, 0}}, {{0, 1, 0}}, {{0, 0, 3}},
        {{huge, 0, 0}}, {{0, huge, 0}}, {{0, 0, huge}}
    };
    for (const auto& index : invalid) {
        require_throws<std::out_of_range>([&] {
            grid.flatten(index[0], index[1], index[2]);
        }, "invalid cell index");
        require_throws<std::out_of_range>([&] {
            grid.cell_volume(index[0], index[1], index[2]);
        }, "invalid cell volume index");
    }
    for (std::size_t flat : {grid.size(), huge}) {
        require_throws<std::out_of_range>([&] { grid.unflatten(flat); },
                                          "invalid flat index");
        require_throws<std::out_of_range>([&] { grid.cell_volume(flat); },
                                          "invalid flat volume index");
    }
}

void test_numerical_range() {
    const double pi = std::acos(-1.0);
    // Power-of-two scales make the exact answer independent of intermediate
    // enormous/tiny cubes; the final volume is comfortably representable.
    const double huge = std::ldexp(1.0, 500);
    const double tiny = std::ldexp(1.0, -500);
    const PhaseSpaceGrid scaled({0.0, huge}, {0.0, tiny}, {-1.0, 1.0});
    require_close(scaled.cell_volume(0), 16.0 * pi * pi / 9.0,
                  "intermediate cube range must not reject valid volume");

    const double thin_hi = std::nextafter(1.0, 2.0);
    const double width = thin_hi - 1.0;
    const PhaseSpaceGrid thin({1.0, thin_hi}, {0.0, 1.0}, {-1.0, 1.0});
    // Integral of (1+x)^2 over x in [0, width], avoiding cube subtraction.
    const double thin_integral = width + width * width + width * width * width / 3.0;
    require_close(thin.cell_volume(0), 16.0 * pi * pi / 3.0 * thin_integral,
                  "adjacent representable radial faces");

    require_throws<std::overflow_error>([&] {
        PhaseSpaceGrid grid({0.0, huge}, {0.0, huge}, {-1.0, 1.0});
    }, "unrepresentable final volume");
    require_throws<std::underflow_error>([&] {
        PhaseSpaceGrid grid({0.0, tiny}, {0.0, tiny}, {-1.0, 1.0});
    }, "volume rounded to zero");

    const double small = std::ldexp(1.0, -350);
    const PhaseSpaceGrid subnormal({0.0, small}, {0.0, 1.0}, {-1.0, 1.0});
    require(subnormal.cell_volume(0) > 0.0 &&
            subnormal.cell_volume(0) < std::numeric_limits<double>::min(),
            "positive subnormal volume is retained");
}

}  // namespace

int main() {
    try {
        test_volume_and_indices();
        test_invalid_faces();
        test_invalid_indices();
        test_numerical_range();
        std::cout << "phase-space grid tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "phase-space grid test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
