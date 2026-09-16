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
            self.contract["status"],
            "physics_validation_in_progress_reference_parity_not_evaluated")
        report_policy = self.contract["report_status_policy"]
        parity_status = report_policy["reference_parity"]
        self.assertEqual(parity_status["status"], "not_evaluated")
        self.assertEqual(
            parity_status["reason_code"],
            "independent_artifact_unavailable")
        self.assertEqual(parity_status["gate"], "conditional")
        self.assertIn("externally produced frozen artifact",
                      parity_status["activation_condition"])
        self.assertEqual(
            set(parity_status["prohibited_substitutes"]),
            {
                "source_derived_expected_values",
                "project_implementation_outputs",
                "analytic_physics_checks",
                "t01_small_trajectory_samples",
            },
        )
        physics_status = report_policy["physics_validation"]
        self.assertEqual(physics_status["status"], "in_progress")
        self.assertEqual(physics_status["gate"], "required")
        self.assertTrue(physics_status["independent_of_reference_artifact"])
        self.assertIn("reopens G0", report_policy["g0_completion_rule"])
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
        self.assertEqual(
            selection["event_prefix"],
            {
                "seed": 20260915,
                "samples_per_state": 64,
                "fields": ["selected_target_index"],
            },
        )
        self.assertIn("standard-library",
                      selection["distribution_implementation"])
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
        state_groups = parity["state_groups"]
        fixed_target = state_groups["conditioned_target_velocity"]
        full_chain = state_groups["full_chain"]
        self.assertEqual(
            fixed_target["speed_over_target_thermal_speed"],
            [0.1, 1.0, 3.0, 10.0],
        )
        self.assertNotIn(
            0.0, fixed_target["speed_over_target_thermal_speed"])
        self.assertIn("rejects non-positive", parity[
            "zero_speed_excluded_reason"])
        self.assertEqual(fixed_target["radius_Rsun"], [0.1, 0.5, 0.9])
        self.assertEqual(fixed_target["target_indices"], [0, 2, 13])
        self.assertEqual(
            fixed_target["incoming_direction_lab"], [0.0, 0.0, 1.0])
        self.assertIn("Cartesian product", fixed_target["state_construction"])
        self.assertIn("externally fixed", fixed_target["state_construction"])
        self.assertEqual(
            fixed_target["entry_point"],
            "sample_collision_conditioned_target_velocity_cm_s",
        )
        self.assertNotIn(
            "selected_target_index", fixed_target["event_prefix_fields"])
        self.assertEqual(
            fixed_target["event_prefix_fields"],
            ["target_velocity_cm_s_xyz"],
        )
        self.assertIn("no target-selection draw", fixed_target["target_policy"])
        self.assertEqual(
            full_chain["entry_point"], "sample_sd_proton_collision")
        self.assertEqual(
            full_chain["incoming_direction_lab"], [0.0, 0.0, 1.0])
        self.assertIn(
            "target_selection_reference.radius_Rsun",
            full_chain["state_source"],
        )
        self.assertIn(
            "selected_target_index", full_chain["event_prefix_fields"])
        self.assertIn("never force", full_chain["target_policy"])
        self.assertIn(
            "target_selection_count", full_chain["aggregate_fields"])
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
        self.assertIn("confidence_interval", fixed_target["aggregate_fields"])
        self.assertIn("outgoing_velocity_cm_s_xyz",
                      full_chain["event_prefix_fields"])
        self.assertIn("dark_matter_energy_change_eV",
                      parity["definitions"])
        for histogram in parity["histograms"].values():
            self.assertLess(histogram["minimum"], histogram["maximum"])
            self.assertGreater(histogram["bins"], 0)
            self.assertIsInstance(histogram["include_overflow"], bool)
        self.assertEqual(
            set(fixed_target["aggregate_fields"]),
            {
                "sample_count",
                "target_speed_histogram",
                "relative_speed_histogram",
                "target_angle_histogram",
                "confidence_interval",
            },
        )

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
        prefix_policy = comparison["rng_prefix_environment_policy"]
        self.assertEqual(
            prefix_policy["uniform_mapping"],
            "std::uniform_real_distribution<double>",
        )
        self.assertIn(
            "matching_cpp_standard_library_implementation_and_version",
            prefix_policy["exact_prefix_requires"],
        )
        self.assertIn(
            "matching_compiler_and_math_build_options",
            prefix_policy["exact_prefix_requires"],
        )
        self.assertIn("not_evaluated",
                      prefix_policy["standard_library_mismatch"])
        self.assertIn("Bonferroni", comparison["histogram_interval"])
        tiers = self.contract["execution_tiers"]
        self.assertEqual(tiers["fast_ci"]["runner"], "ctest -L fast")
        self.assertLessEqual(
            tiers["fast_ci"]["target_runtime_seconds"], 60)
        self.assertFalse(
            tiers["fast_ci"]["may_complete_t03_physics_validation"])
        self.assertFalse(
            tiers["scientific_validation"]["included_in_default_ctest"])
        self.assertEqual(
            tiers["scientific_validation"]["runner_status"],
            "implemented_scientific_acceptance_pending",
        )
        self.assertEqual(
            tiers["scientific_validation"]["report_generator_status"],
            "implemented_scientific_acceptance_pending",
        )
        self.assertEqual(tiers["scientific_validation"]["runner"],
                         "Code/python/run_t03_validation.py")
        self.assertEqual(
            tiers["scientific_validation"]["intended_execution"],
            "manual_nightly_or_hpc",
        )
        self.assertEqual(
            tiers["scientific_validation"]["report_path"],
            "Output/Result/validation/validation_report.json")
        self.assertIn(
            "thermal_bath_equilibrium_weak_residual",
            tiers["scientific_validation"]["scope"])
        self.assertEqual(
            self.contract["validation_report_sections"],
            ["reference_parity", "physics_validation"],
        )
        self.assertIn("no local background", self.contract["unavailable_reason"])

    def test_scientific_protocol_is_frozen_and_independent(self) -> None:
        protocol = self.contract["scientific_physics_validation"]
        self.assertEqual(protocol["mode"], "full")
        self.assertEqual(protocol["runner"],
                         "Code/python/run_t03_validation.py")
        self.assertEqual(protocol["sampler_executable"],
                         "t03_validation_sampler")
        self.assertEqual(protocol["report_path"],
                         self.contract["execution_tiers"][
                             "scientific_validation"]["report_path"])
        common = protocol["common"]
        self.assertEqual(len(common["seeds"]), 4)
        self.assertEqual(len(set(common["seeds"])), 4)
        self.assertTrue(set(common["seeds"]).isdisjoint(
            common["pilot_seeds_excluded_from_full"]))
        self.assertGreaterEqual(common["samples_per_seed"], 1000000)
        self.assertEqual(common["confidence_level"], 0.95)
        self.assertEqual(common["bonferroni_scalar_count"], 162)
        self.assertIn("disjoint", common["holdout_policy"])
        canonical = common["canonical_source_policy"]
        self.assertTrue(canonical["required_for_g0"])
        self.assertEqual(canonical["noncanonical_status"], "not_evaluated")
        self.assertEqual(canonical["smoke_status"], "not_evaluated")
        self.assertIn("Code/configs/validation/t03_oracle_contract.json",
                      canonical["tracked_clean_paths"])
        self.assertIn("Code/tools/t03_validation_sampler.cpp",
                      canonical["tracked_clean_paths"])
        self.assertIn("clean tracked", canonical["sampler_build_requirement"])
        self.assertEqual(common["smoke_override_policy"]["mode"],
                         "smoke")
        self.assertEqual(common["smoke_override_policy"][
            "physics_validation_status"], "not_evaluated")
        self.assertFalse(common["smoke_override_policy"]["g0_eligible"])

        thermal = protocol["thermal_equilibrium"]
        self.assertEqual(thermal["radius_Rsun"], [0.3, 0.7])
        self.assertEqual(thermal["tchi_over_t"], 1.0)
        edges = thermal["speed_bin_edges_over_bath_dm_thermal_speed"]
        self.assertEqual(len(edges), 7)
        self.assertIsNone(edges[-1])
        self.assertEqual(edges[0], 0.0)
        self.assertEqual(edges[:-1], sorted(set(edges[:-1])))
        def cdf(x: float) -> float:
            return (math.erf(x) -
                    2.0 * x * math.exp(-x*x) / math.sqrt(math.pi))
        for edge, quantile in zip(edges[:-1],
                                  [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]):
            self.assertAlmostEqual(cdf(edge), quantile, places=12)
        self.assertIn("Gamma_total", thermal["estimator"])
        self.assertGreater(thermal["acceptance"][
            "speed_bin_max_absolute_mean_plus_bonferroni_se"], 0.0)
        self.assertGreater(thermal["acceptance"][
            "energy_max_absolute_mean_plus_bonferroni_se"], 0.0)
        self.assertGreater(thermal["acceptance"][
            "analytic_mb_rate_max_relative_difference_plus_bonferroni_se"],
            0.0)

        heating = protocol["heating_cooling"]
        self.assertEqual(heating["tchi_over_t"], [0.5, 1.0, 2.0])
        self.assertTrue(heating["reuse_equilibrium_at_tchi_over_t_one"])
        self.assertGreater(heating["acceptance"][
            "cold_minimum_lower_confidence_bound"], 0.0)
        self.assertLess(heating["acceptance"][
            "hot_maximum_upper_confidence_bound"], 0.0)

        rotation = protocol["full_chain_rotation"]
        self.assertEqual(rotation["radius_Rsun"], [0.3, 0.7])
        self.assertEqual(len(rotation["incoming_directions_lab"]), 2)
        for direction in rotation["incoming_directions_lab"]:
            self.assertAlmostEqual(sum(x*x for x in direction), 1.0)
        self.assertNotEqual(rotation["incoming_directions_lab"][0],
                            rotation["incoming_directions_lab"][1])
        rotated_seeds = {
            seed + rotation["rotated_seed_offset"]
            for seed in common["seeds"]
        }
        self.assertTrue(rotated_seeds.isdisjoint(common["seeds"]))
        self.assertTrue(rotated_seeds.isdisjoint(
            common["pilot_seeds_excluded_from_full"]))
        self.assertEqual(set(rotation["vectors"]),
                         {"target_velocity", "outgoing_dm_velocity"})
        self.assertEqual(rotation[
            "component_cdf_thresholds_over_bath_dm_thermal_speed"],
            [-2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0])
        self.assertEqual(len(rotation["moments"]), 9)
        self.assertEqual(rotation["moment_vector"], "outgoing_dm_velocity")
        cm = protocol["cm_scattering_angle"]
        self.assertEqual(cm["radius_Rsun"], rotation["radius_Rsun"])
        self.assertEqual(cm["speed_over_escape"],
                         rotation["speed_over_escape"])
        self.assertTrue(cm["reuse_full_chain_rotation_base_batches"])
        self.assertEqual(len(cm["cosine_histogram_edges"]), 11)
        self.assertEqual(cm["cosine_histogram_edges"],
                         [round(-1.0 + 0.2 * index, 1)
                          for index in range(11)])
        self.assertIn("incoming lab DM", cm["observable"])
        self.assertGreater(cm["acceptance"][
            "cdf_max_absolute_difference_plus_bonferroni_se"], 0.0)
        self.assertEqual(
            common["bonferroni_scalar_count"],
            len(thermal["radius_Rsun"]) * (6 + 1 + 1) +
            (len(heating["tchi_over_t"]) - 1) +
            len(rotation["radius_Rsun"]) *
            (len(rotation["vectors"]) * 3 * len(rotation[
                "component_cdf_thresholds_over_bath_dm_thermal_speed"]) +
             len(rotation["moments"])) +
            len(cm["radius_Rsun"]) *
            (len(cm["cosine_histogram_edges"]) - 2),
        )
        self.assertIn("inverse", rotation["comparison"])

        tail = protocol["near_escape_tail"]
        self.assertEqual(tail["speed_over_escape"],
                         [0.8, 0.95, 0.99, 1.01])
        self.assertEqual(tail["passive_velocity_ceiling_over_escape"],
                         [1.5, 2.0, 3.0])
        self.assertIn("never truncate", tail["ceiling_rule"])
        self.assertEqual(tail["precision_status"]["role_in_physics_gate"],
                         "diagnostic_only")
        self.assertGreater(tail["precision_status"][
            "heterogeneity_flag_z_threshold"], 0.0)
        self.assertNotIn("near_escape_tail", protocol["physics_gate_checks"])
        self.assertIn("existing_scattering_angle_distribution",
                      protocol["physics_gate_checks"])
        self.assertIn("cm_scattering_angle_distribution",
                      protocol["physics_gate_checks"])


if __name__ == "__main__":
    unittest.main()
