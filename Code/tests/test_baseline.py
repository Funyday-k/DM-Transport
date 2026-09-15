"""Regression checks for read-only validation of the frozen T01 bundle."""

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Code/python"))
import run_baseline
from run_baseline import _project_path, validate_baseline_bundle


def synthetic_bundle(config: dict) -> tuple[dict, dict, dict, bytes]:
    """Build the smallest valid bundle without paths or ignored artifacts."""
    frozen_config_bytes = (
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    analysis = {"method": "in-memory unit-test fixture"}
    cases = []
    performance_cases = []
    for index, (mode, sigma) in enumerate(
        (mode, sigma)
        for mode in ("Capture", "Parameter point")
        for sigma in config["benchmark"]["cross_sections_cm2"]
    ):
        case_id = f"{'capture' if mode == 'Capture' else 'ordinary'}-sigma{sigma:.0e}"
        wall_time = float(index + 1)
        attempted = config["benchmark"]["attempts_per_case"]
        complete_evaporations = None if mode == "Capture" else 0
        cases.append({
            "id": case_id,
            "mode": mode,
            "sigma_cm2": sigma,
            "status": "passed",
            "wall_time_s": wall_time,
            "smoke_summary_complete": True,
            "physical_baseline_qualified": False,
            "contract_checks": {
                "runtime_model_verified": True,
                "artifact_mass_sigma_verified": None if mode == "Capture" else True,
            },
            "metrics": {
                "attempted": attempted,
                "classified": attempted,
                "unresolved_noncaptures": 0,
                "captured": 0,
                "capture_probability_raw": 0.0,
                "capture_probability_classified": 0.0,
                "complete_evaporations": complete_evaporations,
            },
        })
        performance_cases.append({
            "id": case_id,
            "status": "passed",
            "wall_time_s": wall_time,
            "attempted_trajectories": attempted,
            "attempted_trajectories_per_wall_s": attempted / wall_time,
            "complete_evaporations": complete_evaporations,
        })
    manifest = {
        "schema_version": 1,
        "status": "completed_smoke_after_mpi_retry",
        "config": copy.deepcopy(config),
        "config_sha256": hashlib.sha256(frozen_config_bytes).hexdigest(),
        "build_matches_source_head": True,
        "build_matches_legacy_1au": True,
        "total_wall_time_s": 10.0,
        "analysis": analysis,
        "cases": cases,
    }
    validation = {
        "schema_version": 1,
        "status": manifest["status"],
        "build_matches_source_head": True,
        "source_unchanged": True,
        "ctest_all_passed": True,
        "ctest_tests": [{"name": "fixture", "status": "passed"}],
        "ctest_command": {"wall_time_s": 2.0},
        "mvp_physics_acceptance": False,
        "absolute_capture_rate_particles_per_s": None,
        "occupation_validation": None,
        "lifetime_validation": None,
        "finite_age_validation": None,
        "analysis": analysis,
        "cases": copy.deepcopy(cases),
    }
    performance = {
        "schema_version": 1,
        "total_wall_time_s": manifest["total_wall_time_s"],
        "ctest_wall_time_s": 2.0,
        "physical_absolute_capture_rate_particles_per_s": None,
        "collision_kernel_wall_time_s": None,
        "online_solve_wall_time_s": None,
        "analysis": analysis,
        "cases": performance_cases,
    }
    return manifest, validation, performance, frozen_config_bytes


class BaselineBundleTests(unittest.TestCase):
    """Exercise bundle integrity without requiring ignored local evidence."""

    def setUp(self) -> None:
        self.config = json.loads(
            (ROOT / "Code/configs/benchmark/mvp.json").read_text(encoding="utf-8")
        )
        (self.manifest, self.validation, self.performance,
         self.frozen_config_bytes) = synthetic_bundle(self.config)

    def errors(self, *, config=None, manifest=None, validation=None,
               performance=None, frozen_config_bytes=None) -> list[str]:
        return validate_baseline_bundle(
            self.config if config is None else config,
            self.manifest if manifest is None else manifest,
            self.validation if validation is None else validation,
            self.performance if performance is None else performance,
            self.frozen_config_bytes if frozen_config_bytes is None else frozen_config_bytes,
        )

    def test_minimal_bundle_is_internally_consistent(self) -> None:
        self.assertEqual(self.errors(), [])

    def test_bundle_rejects_config_hash_mismatch(self) -> None:
        errors = self.errors(frozen_config_bytes=self.frozen_config_bytes + b" ")
        self.assertTrue(any("SHA-256 differs" in error for error in errors))

    def test_bundle_rejects_embedded_config_mismatch(self) -> None:
        manifest = copy.deepcopy(self.manifest)
        manifest["config"]["source"]["halo"]["local_density_GeV_cm3"] = 0.35
        errors = self.errors(manifest=manifest)
        self.assertTrue(any("differs from the config embedded" in error for error in errors))

    def test_bundle_rejects_cross_report_case_drift(self) -> None:
        performance = copy.deepcopy(self.performance)
        performance["cases"][1]["attempted_trajectories_per_wall_s"] *= 2
        errors = self.errors(performance=performance)
        self.assertTrue(any("throughput is inconsistent" in error for error in errors))

    def test_bundle_rejects_current_config_drift(self) -> None:
        config = copy.deepcopy(self.config)
        config["source"]["halo"]["local_density_GeV_cm3"] = 0.35
        errors = self.errors(config=config)
        self.assertTrue(any("local_density_GeV_cm3 differs" in error for error in errors))

    def test_bundle_preserves_smoke_only_qualification(self) -> None:
        manifest = copy.deepcopy(self.manifest)
        manifest["cases"][0]["physical_baseline_qualified"] = True
        errors = self.errors(manifest=manifest)
        self.assertTrue(any("must remain unqualified" in error for error in errors))

    def test_local_t01_bundle_when_available(self) -> None:
        baseline = ROOT / "Output/Result/baseline/20260915-initial"
        paths = {
            name: baseline / name
            for name in ("config.json", "run_manifest.json", "validation_report.json",
                         "performance_report.json")
        }
        if not all(path.is_file() for path in paths.values()):
            self.skipTest("local ignored T01 bundle is not present")
        manifest = json.loads(paths["run_manifest.json"].read_text(encoding="utf-8"))
        validation = json.loads(paths["validation_report.json"].read_text(encoding="utf-8"))
        performance = json.loads(paths["performance_report.json"].read_text(encoding="utf-8"))
        self.assertEqual(validate_baseline_bundle(
            self.config, manifest, validation, performance, paths["config.json"].read_bytes()
        ), [])


class ProjectPathTests(unittest.TestCase):
    """Keep every CLI read within the selected project root."""

    def test_rejects_absolute_path_outside_project(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary).resolve()
            project = base / "project"
            outside = base / "outside"
            project.mkdir()
            outside.mkdir()
            with patch.object(run_baseline, "PROJECT_ROOT", project.resolve()):
                with self.assertRaisesRegex(ValueError, "must be inside"):
                    _project_path(outside.resolve(), "test path")

    def test_rejects_parent_escape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary).resolve()
            project = base / "project"
            outside = base / "outside"
            project.mkdir()
            outside.mkdir()
            with patch.object(run_baseline, "PROJECT_ROOT", project.resolve()):
                with self.assertRaises(ValueError):
                    _project_path(Path("../outside"), "test path")

    def test_rejects_internal_symlink_to_outside_project(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary).resolve()
            project = base / "project"
            outside = base / "outside"
            project.mkdir()
            outside.mkdir()
            link = project / "escaped-link"
            link.symlink_to(outside, target_is_directory=True)
            with patch.object(run_baseline, "PROJECT_ROOT", project.resolve()):
                with self.assertRaisesRegex(ValueError, "resolves outside"):
                    _project_path(link, "test path")


if __name__ == "__main__":
    unittest.main()
