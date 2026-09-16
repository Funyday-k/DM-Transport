"""Check T03 statistics, smoke labeling and the local sampler's data contract."""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "Code/python/run_t03_validation.py"
CONTRACT = ROOT / "Code/configs/validation/t03_oracle_contract.json"
SAMPLER = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else None
sys.argv = [sys.argv[0]]
SPEC = importlib.util.spec_from_file_location("run_t03_validation", RUNNER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load the T03 validation runner")
VALIDATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATION)


class T03RunnerTests(unittest.TestCase):
    """Protect scientific acceptance against vacuous or mislabeled passes."""

    def test_pooled_second_moments_and_equivalence_boundary(self) -> None:
        """Use four known values to validate pooled SE and confidence-band logic."""
        records = [
            {"sample_count": 2, "value": {"sum": 3.0, "sum_sq": 5.0}},
            {"sample_count": 2, "value": {"sum": 7.0, "sum_sq": 25.0}},
        ]
        pooled = VALIDATION._object_moments(records, "value")
        self.assertEqual(pooled["mean"], 2.5)
        self.assertAlmostEqual(pooled["standard_error"], math.sqrt(5.0 / 12.0))
        self.assertEqual(VALIDATION._equivalence(
            {"mean": 0.002, "standard_error": 0.0002},
            critical=3.0, tolerance=0.003)["status"], "pass")
        self.assertEqual(VALIDATION._equivalence(
            {"mean": 0.0026, "standard_error": 0.0002},
            critical=3.0, tolerance=0.003)["status"], "fail")

    def test_inconsistent_second_moment_is_rejected(self) -> None:
        """Malformed sampler moments cannot become falsely certain results."""
        with self.assertRaisesRegex(ValueError, "negative variance"):
            VALIDATION._object_moments([
                {"sample_count": 2,
                 "value": {"sum": 2.0, "sum_sq": 1.0}},
            ], "value")
        with self.assertRaisesRegex(ValueError, "negative variance"):
            VALIDATION._named_vector_moment([
                {"sample_count": 2,
                 "canonical_outgoing_vector_over_escape_moments": {
                     "x": {"sum": 2.0, "sum_sq": 1.0}}},
            ], "x")

    def test_wilson_zero_tail_and_rotation_cdf_bug(self) -> None:
        """Zero observations retain an upper bound; anisotropic CDFs fail."""
        lower, upper = VALIDATION._wilson_interval(0, 1000)
        self.assertEqual(lower, 0.0)
        self.assertGreater(upper, 0.0)
        self.assertLess(upper, 0.01)
        check = VALIDATION._cdf_equivalence(
            0, 1000, 1000, 1000, critical=3.0, tolerance=0.01)
        self.assertEqual(check["status"], "fail")

    def test_cm_angle_uniform_cdf_and_anisotropy(self) -> None:
        """The scientific CM-angle gate distinguishes uniform and biased chains."""
        edges = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
                 0.2, 0.4, 0.6, 0.8, 1.0]
        suite = {"cosine_histogram_edges": edges,
                 "acceptance": {
                     "cdf_max_absolute_difference_plus_bonferroni_se": 0.003}}
        record = {
            "sample_count": 1_000_000,
            "radius_Rsun": 0.3,
            "seed_collision": 123,
            "target_counts": [1_000_000] + [0] * 62,
            "cm_scattering_cosine_histogram_edges": edges,
            "cm_scattering_cosine_histogram_counts": [100_000] * 10,
        }
        self.assertEqual(VALIDATION._cm_angle_distribution(
            [record], suite, critical=3.52)["status"], "pass")
        record["cm_scattering_cosine_histogram_counts"] = [1_000_000] + [0] * 9
        self.assertEqual(VALIDATION._cm_angle_distribution(
            [record], suite, critical=3.52)["status"], "fail")

    def test_frozen_family_and_smoke_label(self) -> None:
        """Count the declared family and never label reduced runs as accepted."""
        suite = json.loads(CONTRACT.read_text(encoding="utf-8"))[
            "scientific_physics_validation"]
        self.assertEqual(VALIDATION._declared_family_size(suite), 162)
        self.assertEqual(suite["common"]["bonferroni_scalar_count"], 162)
        tree = {"status": "pass", "checks": [{"status": "fail"}]}
        VALIDATION._mark_smoke_diagnostics(tree)
        self.assertEqual(tree["status"], "not_evaluated")
        self.assertEqual(tree["checks"][0]["status"], "not_evaluated")
        self.assertEqual(tree["checks"][0]["diagnostic_full_bound_result"],
                         "fail")
        self.assertEqual(VALIDATION._contract_qualification(
            ROOT / "Code/configs/benchmark/mvp.json")["status"],
            "diagnostic_noncanonical_contract")

    @unittest.skipIf(SAMPLER is None, "no built sampler supplied")
    def test_smoke_uses_pilot_streams_and_preserves_fast_evidence(self) -> None:
        """Reduced runs cannot consume holdout seeds or erase actual fast passes."""
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        with mock.patch.object(VALIDATION, "_fast_prerequisite",
                               return_value={"status": "pass"}):
            report = VALIDATION.run_validation(
                contract, SAMPLER,
                {"status": "version_controlled_clean",
                 "path": "Code/configs/validation/t03_oracle_contract.json"},
                {"status": "project_local_clean_fresh_build",
                 "path": "build/t03_validation_sampler"},
                smoke=True, smoke_samples_per_seed=128)
        physics = report["physics_validation"]
        full = set(contract["scientific_physics_validation"]["common"]["seeds"])
        offset = contract["scientific_physics_validation"][
            "full_chain_rotation"]["rotated_seed_offset"]
        self.assertTrue(set(report["run_configuration"]["seeds"]).isdisjoint(
            full | {seed + offset for seed in full}))
        self.assertEqual(physics["status"], "not_evaluated")
        self.assertFalse(physics["g0_eligible"])
        self.assertEqual(physics["existing_elastic_kinematics"]["status"],
                         "pass")
        self.assertEqual(physics["thermal_equilibrium"]["status"],
                         "not_evaluated")

    @unittest.skipIf(SAMPLER is None, "no built sampler supplied")
    def test_project_sampler_aggregates_every_event(self) -> None:
        """Both modes preserve 63-target counts and declared weak/bin schema."""
        for mode, parameter, direction in (
                ("thermal", "1.0", []),
                ("fixed", "0.95", ["--direction", "1,2,3"])):
            command = [
                str(SAMPLER), str(ROOT / "Code/data/solar/model_agss09.dat"),
                str(ROOT / "Code/data/solar/Nuclear_Data.txt"),
                mode, "0.3", parameter, "128", "20261002", *direction,
            ]
            completed = subprocess.run(
                command, cwd=ROOT, text=True, capture_output=True,
                timeout=20, check=True)
            data = json.loads(completed.stdout)
            self.assertEqual(data["sample_count"], 128)
            self.assertEqual(len(data["target_counts"]), 63)
            self.assertEqual(sum(data["target_counts"]), 128)
            if mode == "thermal":
                self.assertEqual(len(data["weak_speed_bin_moments"]), 6)
                total_weak_count = sum(
                    item["sum"] for item in data["weak_speed_bin_moments"])
                self.assertAlmostEqual(total_weak_count, 0.0, places=11)
                self.assertGreater(data["gamma_ref_s_inv"], 0.0)
                self.assertGreater(data["analytic_mb_mean_total_rate_s_inv"],
                                   0.0)
            else:
                self.assertEqual(sum(data["tail_target_counts"]),
                                 data["tail_total_count"])
                self.assertEqual(len(data[
                    "canonical_target_component_histogram_counts"]), 3)
                self.assertEqual(len(data[
                    "canonical_outgoing_component_histogram_counts"]), 3)
                self.assertEqual(len(data[
                    "cm_scattering_cosine_histogram_counts"]), 10)
                self.assertEqual(sum(data[
                    "cm_scattering_cosine_histogram_counts"]), 128)
                for key in ("canonical_target_component_histogram_counts",
                            "canonical_outgoing_component_histogram_counts"):
                    for row in data[key]:
                        self.assertEqual(sum(row), 128)
                        self.assertEqual(len(row), 10)
                self.assertTrue(
                    0 <= data["domain_overflow_count_3p0_vesc"] <=
                    data["domain_overflow_count_2p0_vesc"] <=
                    data["domain_overflow_count_1p5_vesc"] <=
                    data["tail_total_count"] <= 128)


if __name__ == "__main__":
    unittest.main()
