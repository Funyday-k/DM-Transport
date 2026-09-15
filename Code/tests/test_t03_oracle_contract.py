"""Validate the frozen import contract for the future independent T03 oracle."""

from __future__ import annotations

import json
import math
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "Code/configs/validation/t03_oracle_contract.json"
REFERENCE_COMMIT = "b5678f5b193aa567ca10715c2a6c764c9e72eec7"
OBSCURA_COMMIT = "4b1b9d87f8a69da4d8081175d322b66e48314328"
LIBPHYSICA_COMMIT = "fefbe47993bbd343277718708407fa7a40629738"


class OracleContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT.read_text(encoding="utf-8"))

    def test_independent_source_boundary(self) -> None:
        self.assertEqual(self.contract["schema_version"], 1)
        self.assertEqual(
            self.contract["contract_id"], "legacy_golden_physics_v1")
        self.assertEqual(
            self.contract["status"], "awaiting_independent_artifact")
        source = self.contract["source"]
        self.assertEqual(source["commit"], REFERENCE_COMMIT)
        self.assertEqual(
            source["dependencies"],
            {
                "obscura_commit": OBSCURA_COMMIT,
                "libphysica_commit": LIBPHYSICA_COMMIT,
            },
        )
        self.assertIn("outside DM-Transport", source["producer"])
        self.assertIn("externally generated", source["import_rule"])
        self.assertEqual(
            set(source["required_producer_metadata"]),
            {
                "source_and_dependency_commits",
                "compiler_and_standard_library",
                "build_type_and_options",
                "applied_patch_state",
                "solar_and_nuclear_data_paths",
                "rng_engine_seed_and_mpi_ranks",
            },
        )
        self.assertEqual(
            set(source["forbidden_generator_components"]),
            {
                "transport_reference_physics",
                "SolarBackground",
                "direct_sd_proton_scattering_rates",
            },
        )

        self.assertEqual(
            self.contract["model"],
            {
                "dark_matter_mass_GeV": 0.1,
                "dark_matter_spin": 0.5,
                "interaction": "SD",
                "form_factor": "Contact",
                "low_mass_mode": True,
                "proton_relative_coupling": 1.0,
                "neutron_relative_coupling": 0.0,
                "electron_cross_section_cm2": 0.0,
                "proton_cross_section_cm2": 1e-34,
                "solar_target_policy":
                    "all_63_isotopes_in_reference_order",
                "direct_rate_without_interpolation": True,
            },
        )

    def test_rate_grid_and_fields(self) -> None:
        rate_grid = self.contract["rate_grid"]
        for axis in ("radius_Rsun", "speed_over_c"):
            values = rate_grid[axis]
            self.assertEqual(len(values), 6)
            self.assertEqual(values, sorted(set(values)))
            self.assertTrue(all(math.isfinite(value) and value >= 0.0
                                for value in values))
        self.assertEqual(
            rate_grid["target_policy"],
            "all_63_isotopes_in_reference_order")
        self.assertEqual(
            set(rate_grid["required_fields"]),
            {
                "temperature_K",
                "mass_enclosed_g",
                "escape_speed_cm_s",
                "target_number_density_cm3",
                "target_cross_section_cm2",
                "target_mean_relative_speed_cm_s",
                "target_rate_s_inv",
                "total_rate_s_inv",
            },
        )

    def test_collision_and_report_layers(self) -> None:
        selection = self.contract["target_selection_reference"]
        self.assertEqual(selection["radius_Rsun"], [0.1, 0.5, 0.9])
        self.assertEqual(selection["speed_over_c"], [0.0001, 0.001, 0.006])
        self.assertEqual(selection["engine"], "std::mt19937")
        self.assertGreaterEqual(
            selection["samples_per_seed_per_state"] *
            len(selection["seeds"]),
            1_000_000,
        )
        self.assertEqual(
            set(selection["required_fields"]),
            {
                "target_rate_s_inv",
                "target_probability",
                "target_selection_count",
                "sample_count",
                "confidence_interval",
            },
        )

        parity = self.contract["collision_reference_parity"]
        self.assertEqual(
            parity["speed_over_target_thermal_speed"],
            [0.1, 1.0, 3.0, 10.0],
        )
        self.assertNotIn(0.0, parity["speed_over_target_thermal_speed"])
        self.assertIn("rejects non-positive", parity[
            "zero_speed_excluded_reason"])
        self.assertEqual(parity["radius_Rsun"], [0.1, 0.5, 0.9])
        self.assertEqual(parity["target_indices"], [0, 2, 13])
        self.assertEqual(parity["incoming_direction_lab"], [0.0, 0.0, 1.0])
        self.assertIn("Cartesian product", parity["state_construction"])
        rng = parity["rng"]
        self.assertEqual(rng["engine"], "std::mt19937")
        self.assertEqual(rng["event_prefix_seed"], 20260915)
        self.assertGreater(rng["event_prefix_per_state"], 0)
        self.assertEqual(len(rng["aggregate_seeds"]), 4)
        self.assertGreaterEqual(
            rng["aggregate_samples_per_seed_per_state"] *
            len(rng["aggregate_seeds"]),
            1_000_000,
        )
        self.assertIn("outgoing_velocity_cm_s_xyz",
                      parity["event_prefix_fields"])
        self.assertIn("confidence_interval", parity["aggregate_fields"])
        self.assertIn("dark_matter_energy_change_eV",
                      parity["definitions"])
        for histogram in parity["histograms"].values():
            self.assertLess(histogram["minimum"], histogram["maximum"])
            self.assertGreater(histogram["bins"], 0)
            self.assertIsInstance(histogram["include_overflow"], bool)

        physics = self.contract["collision_physics_validation"]
        self.assertEqual(
            physics["speed_over_target_thermal_speed"],
            [0.0, 0.1, 1.0, 3.0, 10.0],
        )
        self.assertIn("do not request a zero-speed legacy sample",
                      physics["zero_speed_rule"])
        self.assertIn("thermal_bath_equilibrium",
                      physics["required_checks"])

        comparison = self.contract["comparison_policy"]
        self.assertEqual(comparison["confidence_level"], 0.95)
        self.assertGreater(comparison[
            "deterministic_rate_relative_tolerance"], 0.0)
        self.assertEqual(
            set(comparison["event_prefix_tolerances"]),
            {
                "velocity_relative",
                "velocity_absolute_cm_s",
                "cosine_absolute",
                "energy_relative",
                "energy_absolute_eV",
            },
        )
        self.assertTrue(all(
            value > 0.0
            for value in comparison["event_prefix_tolerances"].values()
        ))
        self.assertIn("Bonferroni", comparison["histogram_interval"])
        self.assertEqual(
            self.contract["validation_report_sections"],
            ["reference_parity", "physics_validation"],
        )
        self.assertIn("no local background", self.contract["unavailable_reason"])


if __name__ == "__main__":
    unittest.main()
