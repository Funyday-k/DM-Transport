#include "transport/physics/SolarBackground.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace transport {
namespace physics {
namespace {

// Port provenance: Code/provenance/reference_physics.json, unit solar_background.
constexpr double legacySolarMassG = 1.98848e33;
constexpr double legacySolarRadiusCm = 6.957e10;
constexpr double legacyGravitationalConstantCgs = 6.6740836890128835e-8;
constexpr double legacyGeVPerGram = 5.60958884493318e23;
constexpr double legacyProtonMassGeV = 0.9382720813;
constexpr double legacyNucleonMassGeV = 0.932;
constexpr std::size_t solarColumnCount = 35;
constexpr std::size_t abundanceColumnOffset = 6;
constexpr std::size_t expectedSolarRowCount = 1968;
constexpr std::size_t expectedNuclearRecordCount = 295;
constexpr std::size_t expectedTargetCount = 63;

int signum(double value) {
    return (value > 0.0) - (value < 0.0);
}

class MonotoneCubic {
public:
    MonotoneCubic(std::vector<double> x, std::vector<double> y)
        : x_(std::move(x)) {
        if (x_.size() != y.size() || x_.size() < 3) {
            throw std::runtime_error(
                "monotone interpolation requires at least three x-y pairs");
        }
        const std::size_t count = x_.size();
        std::vector<double> width(count - 1);
        std::vector<double> slope(count - 1);
        for (std::size_t i = 0; i + 1 < count; ++i) {
            if (!std::isfinite(x_[i]) || !std::isfinite(y[i]) ||
                !std::isfinite(x_[i + 1]) || !std::isfinite(y[i + 1]) ||
                x_[i + 1] <= x_[i]) {
                throw std::runtime_error(
                    "interpolation coordinates must be finite and increasing");
            }
            width[i] = x_[i + 1] - x_[i];
            slope[i] = (y[i + 1] - y[i]) / width[i];
        }

        std::vector<double> derivative(count);
        const double left_candidate =
            slope[0] * (1.0 + width[0] / (width[0] + width[1])) -
            slope[1] * width[0] / (width[0] + width[1]);
        derivative[0] = (signum(left_candidate) + signum(slope[0])) *
            std::min(std::abs(slope[0]), 0.5 * std::abs(left_candidate));

        const std::size_t last = count - 1;
        const double right_candidate =
            slope[last - 1] *
                (1.0 + width[last - 1] /
                           (width[last - 1] + width[last - 2])) -
            slope[last - 2] * width[last - 1] /
                (width[last - 1] + width[last - 2]);
        derivative[last] =
            (signum(right_candidate) + signum(slope[last - 1])) *
            std::min(std::abs(slope[last - 1]),
                     0.5 * std::abs(right_candidate));

        for (std::size_t i = 1; i < last; ++i) {
            const double candidate =
                (slope[i - 1] * width[i] +
                 slope[i] * width[i - 1]) /
                (width[i - 1] + width[i]);
            derivative[i] = (signum(slope[i - 1]) + signum(slope[i])) *
                std::min(0.5 * std::abs(candidate),
                         std::min(std::abs(slope[i - 1]),
                                  std::abs(slope[i])));
        }

        a_.reserve(last);
        b_.reserve(last);
        c_.reserve(last);
        d_.reserve(last);
        for (std::size_t i = 0; i < last; ++i) {
            const double width_squared = width[i] * width[i];
            const double a = (derivative[i] + derivative[i + 1] -
                              2.0 * slope[i]) / width_squared;
            const double b = (3.0 * slope[i] - 2.0 * derivative[i] -
                              derivative[i + 1]) / width[i];
            if (!std::isfinite(a) || !std::isfinite(b) ||
                !std::isfinite(derivative[i])) {
                throw std::runtime_error(
                    "monotone interpolation coefficients are not finite");
            }
            a_.push_back(a);
            b_.push_back(b);
            c_.push_back(derivative[i]);
            d_.push_back(y[i]);
        }
    }

    double evaluate(double x) const {
        if (!std::isfinite(x) || x < x_.front() || x > x_.back()) {
            throw std::out_of_range("interpolation coordinate is outside its domain");
        }
        std::vector<double>::const_iterator upper =
            std::upper_bound(x_.begin(), x_.end(), x);
        std::size_t interval = upper == x_.end()
            ? x_.size() - 2
            : static_cast<std::size_t>(upper - x_.begin() - 1);
        const double offset = x - x_[interval];
        return ((a_[interval] * offset + b_[interval]) * offset +
                c_[interval]) * offset + d_[interval];
    }

private:
    std::vector<double> x_;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
};

using SolarRow = std::array<double, solarColumnCount>;

std::vector<SolarRow> load_solar_rows(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) {
        throw std::runtime_error("unable to open solar model table: " + path);
    }
    input.imbue(std::locale::classic());

    std::vector<SolarRow> rows;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream fields(line);
        fields.imbue(std::locale::classic());
        SolarRow row{};
        std::size_t parsed_columns = 0;
        for (std::size_t column = 0; column < row.size(); ++column) {
            if (!(fields >> row[column])) {
                break;
            }
            ++parsed_columns;
        }
        if (parsed_columns != row.size()) {
            const std::size_t first = line.find_first_not_of(" \t\r");
            const bool ignorable = first == std::string::npos ||
                line[first] == '#';
            if (rows.empty() || ignorable) {
                continue;
            }
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": expected exactly 35 numeric columns");
        }
        std::string extra;
        if (fields >> extra) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": expected exactly 35 numeric columns");
        }
        for (double value : row) {
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    path + ":" + std::to_string(line_number) +
                    ": nonfinite solar-model value");
            }
        }
        if (row[0] <= 0.0 || row[0] > 1.0 || row[1] <= 0.0 ||
            row[1] >= 1.0 || row[2] <= 0.0 || row[3] < 0.0) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": solar-model value is outside its physical domain");
        }
        if (!rows.empty() &&
            (row[0] <= rows.back()[0] || row[1] <= rows.back()[1])) {
            throw std::runtime_error(
                path + ":" + std::to_string(line_number) +
                ": mass and radius columns must be strictly increasing");
        }
        for (std::size_t column = abundanceColumnOffset;
             column < solarColumnCount; ++column) {
            if (row[column] < 0.0 || row[column] > 1.0) {
                throw std::runtime_error(
                    path + ":" + std::to_string(line_number) +
                    ": abundance columns must lie in [0, 1]");
            }
        }
        rows.push_back(row);
    }
    if (input.bad()) {
        throw std::runtime_error("failed while reading solar model table: " +
                                 path);
    }
    if (rows.size() != expectedSolarRowCount) {
        throw std::runtime_error(
            path + ": expected 1968 solar-model rows, found " +
            std::to_string(rows.size()));
    }
    if (rows.front()[0] != 0.0000004 || rows.front()[1] != 0.00150 ||
        rows.back()[0] != 0.9999930 || rows.back()[1] != 0.98500) {
        throw std::runtime_error(
            path + ": unexpected AGSS09 radial endpoints");
    }
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const double expected_radius =
            0.00150 + 0.00050 * static_cast<double>(index);
        if (std::abs(rows[index][1] - expected_radius) > 1.0e-12) {
            throw std::runtime_error(
                path + ": AGSS09 radius grid is not the expected sequence");
        }
    }

    SolarRow center = rows.front();
    center[0] = 0.0;
    center[1] = 0.0;
    rows.insert(rows.begin(), center);

    SolarRow surface = rows.back();
    surface[0] = 1.0;
    surface[1] = 1.0;
    surface[2] = 5800.0;
    surface[3] = 1.0e-9;
    rows.push_back(surface);
    return rows;
}

struct NuclearRecord {
    std::string symbol;
    unsigned int atomic_number;
    unsigned int mass_number;
    double abundance;
    double spin;
    double proton_spin;
    double neutron_spin;
};

std::vector<NuclearRecord> load_nuclear_records(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) {
        throw std::runtime_error("unable to open nuclear data table: " + path);
    }
    input.imbue(std::locale::classic());
    std::vector<NuclearRecord> records;
    std::size_t record_number = 0;
    while (true) {
        NuclearRecord record{};
        if (!(input >> record.symbol)) {
            if (input.eof()) {
                break;
            }
            throw std::runtime_error(
                "failed while reading nuclear data table: " + path);
        }
        ++record_number;
        if (!(input >> record.atomic_number >> record.mass_number >>
              record.abundance >> record.spin >> record.proton_spin >>
              record.neutron_spin)) {
            throw std::runtime_error(
                path + ": malformed nuclear record " +
                std::to_string(record_number));
        }
        if (record.atomic_number == 0 || record.atomic_number > 92 ||
            record.mass_number < record.atomic_number ||
            record.mass_number > 238 ||
            !std::isfinite(record.abundance) || record.abundance < 0.0 ||
            record.abundance > 1.0 || !std::isfinite(record.spin) ||
            record.spin < 0.0 ||
            !std::isfinite(record.proton_spin) ||
            !std::isfinite(record.neutron_spin)) {
            throw std::runtime_error(
                "nuclear data table contains invalid values: " + path);
        }
        for (const NuclearRecord& existing : records) {
            if (existing.atomic_number == record.atomic_number &&
                existing.mass_number == record.mass_number) {
                throw std::runtime_error(
                    "nuclear data table contains a duplicate isotope: " + path);
            }
        }
        records.push_back(record);
    }
    if (records.size() != expectedNuclearRecordCount) {
        throw std::runtime_error(
            path + ": expected 295 nuclear records, found " +
            std::to_string(records.size()));
    }
    if (records.front().atomic_number != 1 ||
        records.front().mass_number != 1 ||
        records.back().atomic_number != 92 ||
        records.back().mass_number != 238) {
        throw std::runtime_error(
            path + ": unexpected nuclear-table endpoints");
    }
    return records;
}

const NuclearRecord& find_isotope(const std::vector<NuclearRecord>& records,
                                  unsigned int atomic_number,
                                  unsigned int mass_number) {
    for (const NuclearRecord& record : records) {
        if (record.atomic_number == atomic_number &&
            record.mass_number == mass_number) {
            return record;
        }
    }
    throw std::runtime_error("required solar isotope is missing");
}

struct TargetDefinition {
    SolarTarget target;
    std::size_t composition_column;
    double isotopic_fraction;
};

TargetDefinition make_target(const NuclearRecord& record,
                             std::size_t composition_column,
                             double isotopic_fraction) {
    TargetDefinition definition;
    definition.target.name =
        record.symbol + "-" + std::to_string(record.mass_number);
    definition.target.atomic_number = record.atomic_number;
    definition.target.mass_number = record.mass_number;
    definition.target.spin = record.spin;
    definition.target.proton_spin = record.proton_spin;
    definition.target.neutron_spin = record.neutron_spin;
    definition.target.mass_GeV = record.mass_number == 1
        ? legacyProtonMassGeV
        : record.mass_number * legacyNucleonMassGeV;
    definition.composition_column = composition_column;
    definition.isotopic_fraction = isotopic_fraction;
    return definition;
}

std::vector<TargetDefinition> make_solar_targets(
    const std::vector<NuclearRecord>& records) {
    const std::array<unsigned int, 29> atomic_numbers{{
        1, 2, 2, 6, 6, 7, 7, 8, 8, 8, 10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28
    }};
    const std::array<unsigned int, 10> explicit_mass_numbers{{
        1, 4, 3, 12, 13, 14, 15, 16, 17, 18
    }};

    std::vector<TargetDefinition> targets;
    for (std::size_t column = 0; column < explicit_mass_numbers.size();
         ++column) {
        const NuclearRecord& record = find_isotope(
            records, atomic_numbers[column], explicit_mass_numbers[column]);
        targets.push_back(make_target(record, column, 1.0));
    }
    for (std::size_t column = explicit_mass_numbers.size();
         column < atomic_numbers.size(); ++column) {
        const std::size_t before = targets.size();
        for (const NuclearRecord& record : records) {
            if (record.atomic_number == atomic_numbers[column]) {
                targets.push_back(
                    make_target(record, column, record.abundance));
            }
        }
        if (targets.size() == before) {
            throw std::runtime_error("solar element has no nuclear data");
        }
    }
    if (targets.size() != expectedTargetCount) {
        throw std::runtime_error("unexpected number of solar isotope targets");
    }
    return targets;
}

std::vector<double> column_values(const std::vector<SolarRow>& rows,
                                  std::size_t column, double scale) {
    std::vector<double> values;
    values.reserve(rows.size());
    for (const SolarRow& row : rows) {
        values.push_back(row[column] * scale);
    }
    return values;
}

struct GaussNode {
    double location;
    double weight;
};

void evaluate_legendre_30(double location,
                          double& polynomial,
                          double& derivative) {
    polynomial = 1.0;
    double preceding = 0.0;
    for (unsigned int order = 1; order <= 30; ++order) {
        const double older = preceding;
        preceding = polynomial;
        polynomial =
            ((2.0 * order - 1.0) * location * preceding -
             (order - 1.0) * older) / order;
    }
    derivative = 30.0 * (location * polynomial - preceding) /
        (location * location - 1.0);
}

const std::array<GaussNode, 15>& positive_gauss_legendre_30_nodes() {
    static const std::array<GaussNode, 15> nodes = [] {
        std::array<GaussNode, 15> result{};
        const double pi = std::acos(-1.0);
        const double tolerance = 4.0 * std::numeric_limits<double>::epsilon();
        for (std::size_t i = 0; i < result.size(); ++i) {
            double location = std::cos(
                pi * (static_cast<double>(i) + 0.75) / 30.5);
            bool converged = false;
            for (unsigned int iteration = 0; iteration < 100; ++iteration) {
                double polynomial = 0.0;
                double derivative = 0.0;
                evaluate_legendre_30(location, polynomial, derivative);
                const double updated = location - polynomial / derivative;
                if (std::abs(updated - location) <= tolerance) {
                    location = updated;
                    converged = true;
                    break;
                }
                location = updated;
            }
            if (!converged) {
                throw std::runtime_error(
                    "30-point Gauss-Legendre nodes did not converge");
            }
            double polynomial = 0.0;
            double derivative = 0.0;
            evaluate_legendre_30(location, polynomial, derivative);
            result[i] = GaussNode{
                location,
                2.0 / ((1.0 - location * location) *
                       derivative * derivative)
            };
        }
        return result;
    }();
    return nodes;
}

double integrate_mass_over_radius_squared(const MonotoneCubic& mass,
                                          double lower_radius_cm) {
    if (lower_radius_cm == legacySolarRadiusCm) {
        return 0.0;
    }
    const double midpoint =
        0.5 * (lower_radius_cm + legacySolarRadiusCm);
    const double half_width =
        0.5 * (legacySolarRadiusCm - lower_radius_cm);
    double weighted_sum = 0.0;
    for (const GaussNode& node : positive_gauss_legendre_30_nodes()) {
        const double left = midpoint - half_width * node.location;
        const double right = midpoint + half_width * node.location;
        weighted_sum += node.weight *
            (mass.evaluate(left) / (left * left) +
             mass.evaluate(right) / (right * right));
    }
    return half_width * weighted_sum;
}

void require_radius(double radius_cm) {
    if (!std::isfinite(radius_cm) || radius_cm < 0.0) {
        throw std::invalid_argument("solar radius must be finite and nonnegative");
    }
}

}  // namespace

class SolarBackground::Implementation {
public:
    Implementation(const std::string& solar_model_path,
                   const std::string& nuclear_data_path) {
        const std::vector<SolarRow> rows = load_solar_rows(solar_model_path);
        targets = make_solar_targets(load_nuclear_records(nuclear_data_path));

        const std::vector<double> radii =
            column_values(rows, 1, legacySolarRadiusCm);
        mass.reset(new MonotoneCubic(
            radii, column_values(rows, 0, legacySolarMassG)));
        temperature.reset(new MonotoneCubic(
            radii, column_values(rows, 2, 1.0)));

        number_densities.reserve(targets.size());
        for (const TargetDefinition& definition : targets) {
            std::vector<double> values;
            values.reserve(rows.size());
            for (const SolarRow& row : rows) {
                values.push_back(
                    row[abundanceColumnOffset +
                        definition.composition_column] *
                    row[3] * definition.isotopic_fraction *
                    legacyGeVPerGram / definition.target.mass_GeV);
            }
            // Match the reference order of operations: form each isotope's
            // number density at every radial node, then interpolate n_i(r).
            number_densities.emplace_back(radii, std::move(values));
        }

        std::vector<double> escape_speed_squared;
        escape_speed_squared.reserve(radii.size());
        const double surface_scale =
            2.0 * legacyGravitationalConstantCgs * legacySolarMassG /
            legacySolarRadiusCm;
        for (double radius : radii) {
            const double integral =
                integrate_mass_over_radius_squared(*mass, radius);
            const double value = surface_scale *
                (1.0 + legacySolarRadiusCm / legacySolarMassG * integral);
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error(
                    "solar escape-speed profile is not representable");
            }
            escape_speed_squared.push_back(value);
        }
        escape_speed_squared_profile.reset(new MonotoneCubic(
            radii, std::move(escape_speed_squared)));
    }

    std::vector<TargetDefinition> targets;
    std::unique_ptr<MonotoneCubic> mass;
    std::unique_ptr<MonotoneCubic> temperature;
    std::vector<MonotoneCubic> number_densities;
    std::unique_ptr<MonotoneCubic> escape_speed_squared_profile;
};

SolarBackground::SolarBackground(const std::string& solar_model_path,
                                 const std::string& nuclear_data_path)
    : implementation_(new Implementation(solar_model_path,
                                         nuclear_data_path)) {}

SolarBackground::~SolarBackground() = default;
SolarBackground::SolarBackground(SolarBackground&&) noexcept = default;
SolarBackground& SolarBackground::operator=(SolarBackground&&) noexcept = default;

double SolarBackground::solar_radius_cm() const noexcept {
    return legacySolarRadiusCm;
}

std::size_t SolarBackground::target_count() const noexcept {
    return implementation_->targets.size();
}

const SolarTarget& SolarBackground::target(std::size_t target_index) const {
    return implementation_->targets.at(target_index).target;
}

double SolarBackground::temperature_K(double radius_cm) const {
    require_radius(radius_cm);
    if (radius_cm > legacySolarRadiusCm) {
        throw std::out_of_range(
            "solar temperature is undefined outside the photosphere");
    }
    return implementation_->temperature->evaluate(radius_cm);
}

double SolarBackground::mass_enclosed_g(double radius_cm) const {
    require_radius(radius_cm);
    return radius_cm > legacySolarRadiusCm
        ? legacySolarMassG
        : implementation_->mass->evaluate(radius_cm);
}

double SolarBackground::escape_speed_cm_s(double radius_cm) const {
    require_radius(radius_cm);
    if (radius_cm > legacySolarRadiusCm) {
        return std::sqrt(2.0 * legacyGravitationalConstantCgs *
                         legacySolarMassG / radius_cm);
    }
    return std::sqrt(
        implementation_->escape_speed_squared_profile->evaluate(radius_cm));
}

double SolarBackground::number_density_cm3(std::size_t target_index,
                                           double radius_cm) const {
    require_radius(radius_cm);
    const MonotoneCubic& density =
        implementation_->number_densities.at(target_index);
    return radius_cm > legacySolarRadiusCm ? 0.0 : density.evaluate(radius_cm);
}

}  // namespace physics
}  // namespace transport
