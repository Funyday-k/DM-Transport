#!/usr/bin/env python3
"""Validate the saved T01 baseline bundle without modifying or rerunning it.

The CLI only reads the frozen config snapshot, three JSON reports, and current
MVP JSON from this repository. It does not inspect source checkouts, invoke
external tools, or write into the baseline directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASELINE_DIR = PROJECT_ROOT / "Output/Result/baseline/20260915-initial"
DEFAULT_CONFIG = PROJECT_ROOT / "Code/configs/benchmark/mvp.json"
BUNDLE_NAMES = ("config.json", "run_manifest.json", "validation_report.json",
                "performance_report.json")


def validate_config(config: dict[str, Any]) -> None:
    """Reject settings outside the supported MVP; masses are GeV, areas cm²."""
    expected = {
        "mass_GeV": 0.1, "interaction": "SD", "spin": 0.5, "light": True,
        "proton_relative_coupling": 1.0, "neutron_relative_coupling": 0.0,
        "electron_cross_section_cm2": 0.0, "form_factor": "Contact",
        "solar_target_policy": "damascus_default_isotopes", "rate_interpolation_points": 0,
    }
    if not isinstance(config, dict) or config.get("schema_version") != 1 or config.get("model") != expected:
        raise ValueError("This tool supports only the schema 1 SD proton-only MVP model")
    boundary = config.get("boundary")
    if not isinstance(boundary, dict):
        raise ValueError("boundary must be an object")
    for key, value in {"policy": "legacy_1au_removal", "escape_radius_Rsun": 1.0,
                       "legacy_matching_radius_Rsun": 1.1, "outer_removal_radius_au": 1.0}.items():
        if boundary.get(key) != value:
            raise ValueError(f"Unsupported boundary setting: {key}")
    benchmark = config.get("benchmark")
    if not isinstance(benchmark, dict):
        raise ValueError("benchmark must be an object")
    if benchmark.get("cross_sections_cm2") != [1e-36, 1e-34] or benchmark.get("mpi_ranks") != 1:
        raise ValueError("This MVP requires two cross sections and one MPI rank")
    source = config.get("source")
    halo = source.get("halo") if isinstance(source, dict) else None
    if (not isinstance(halo, dict) or halo.get("distribution") != "SHM"
            or not isinstance(halo.get("observer_velocity_km_s"), list)
            or len(halo["observer_velocity_km_s"]) != 3):
        raise ValueError("The MVP requires SHM with a three-component observer velocity")
    for key, upper in {"attempts_per_case": 16, "seed": 2147483647}.items():
        value = benchmark.get(key)
        if type(value) is not int or not 1 <= value <= upper:
            raise ValueError(f"{key} must be an integer in [1, {upper}]")
    for key, upper in {"trajectory_wall_timeout_s": 2, "case_timeout_s": 60,
                       "ctest_timeout_s": 60, "ctest_total_timeout_s": 600}.items():
        value = benchmark.get(key)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not 0 < value <= upper:
            raise ValueError(f"{key} must be positive and at most {upper} wall seconds")


def _snapshot_mismatches(frozen: Any, current: Any, path: str = "config") -> list[str]:
    """Compare a frozen JSON value with a current value that may add fields."""
    if isinstance(frozen, dict) and isinstance(current, dict):
        errors = []
        for key, value in frozen.items():
            child = f"{path}.{key}"
            if key not in current:
                errors.append(f"{child} is missing from the current MVP config")
            else:
                errors.extend(_snapshot_mismatches(value, current[key], child))
        return errors
    if isinstance(frozen, list) and isinstance(current, list):
        if len(frozen) != len(current):
            return [f"{path} length differs from the frozen T01 config"]
        errors = []
        for index, (old, new) in enumerate(zip(frozen, current)):
            errors.extend(_snapshot_mismatches(old, new, f"{path}[{index}]"))
        return errors
    if type(frozen) is not type(current) or frozen != current:
        return [f"{path} differs from the frozen T01 config"]
    return []


def _index_cases(report: dict[str, Any], label: str, errors: list[str]) -> dict[str, dict[str, Any]]:
    """Index case objects while reporting malformed or duplicate identifiers."""
    cases = report.get("cases")
    if not isinstance(cases, list):
        errors.append(f"{label}.cases must be a list")
        return {}
    indexed: dict[str, dict[str, Any]] = {}
    for case in cases:
        if not isinstance(case, dict) or not isinstance(case.get("id"), str):
            errors.append(f"{label}.cases contains an entry without a string id")
            continue
        case_id = case["id"]
        if case_id in indexed:
            errors.append(f"{label}.cases contains duplicate id {case_id}")
        indexed[case_id] = case
    return indexed


def _matches_ratio(value: Any, numerator: int, denominator: float) -> bool:
    """Return whether a JSON number equals a defined finite ratio."""
    return (denominator > 0 and isinstance(value, (int, float))
            and not isinstance(value, bool) and math.isfinite(value)
            and math.isclose(value, numerator / denominator, rel_tol=1e-12, abs_tol=0))


def validate_baseline_bundle(config: dict[str, Any], manifest: dict[str, Any],
                             validation: dict[str, Any],
                             performance: dict[str, Any],
                             frozen_config_bytes: bytes) -> list[str]:
    """Return consistency failures for the frozen T01 reports and current MVP."""
    errors: list[str] = []
    documents = {"mvp config": config, "run manifest": manifest,
                 "validation report": validation, "performance report": performance}
    if any(not isinstance(document, dict) for document in documents.values()):
        return [f"{name} must contain a JSON object" for name, document in documents.items()
                if not isinstance(document, dict)]
    for name, document in documents.items():
        if document.get("schema_version") != 1:
            errors.append(f"{name} must use schema_version 1")
    try:
        validate_config(config)
    except (AttributeError, KeyError, TypeError, ValueError) as error:
        errors.append(f"invalid MVP config: {error}")
        return errors

    benchmark = config["benchmark"]
    regression_sigma = benchmark.get("regression_cross_section_cm2")
    if regression_sigma not in benchmark["cross_sections_cm2"]:
        errors.append("benchmark.regression_cross_section_cm2 must select a baseline cross section")
    try:
        frozen_config = json.loads(frozen_config_bytes)
    except (UnicodeError, json.JSONDecodeError) as error:
        frozen_config = None
        errors.append(f"invalid frozen config.json: {error}")
    if not isinstance(frozen_config, dict):
        errors.append("frozen config.json must contain a JSON object")
    elif frozen_config != manifest.get("config"):
        errors.append("frozen config.json differs from the config embedded in the run manifest")
    if hashlib.sha256(frozen_config_bytes).hexdigest() != manifest.get("config_sha256"):
        errors.append("frozen config.json SHA-256 differs from run manifest config_sha256")
    if not isinstance(manifest.get("config"), dict):
        errors.append("run manifest must embed its frozen config")
    if isinstance(frozen_config, dict):
        try:
            validate_config(frozen_config)
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            errors.append(f"invalid frozen T01 config: {error}")
        errors.extend(_snapshot_mismatches(frozen_config, config))

    expected_status = "completed_smoke_after_mpi_retry"
    if manifest.get("status") != expected_status:
        errors.append(f"run manifest status must be {expected_status}")
    if validation.get("status") != manifest.get("status"):
        errors.append("validation status differs from run manifest status")
    for key in ("build_matches_source_head", "build_matches_legacy_1au"):
        if manifest.get(key) is not True:
            errors.append(f"run manifest {key} is not true")
    for key in ("build_matches_source_head", "source_unchanged", "ctest_all_passed"):
        if validation.get(key) is not True:
            errors.append(f"validation report {key} is not true")
    tests = validation.get("ctest_tests")
    if not isinstance(tests, list) or not tests:
        errors.append("validation report has no final CTest results")
    elif any(not isinstance(test, dict) or test.get("status") != "passed" for test in tests):
        errors.append("validation report contains a non-passing final CTest result")
    if validation.get("mvp_physics_acceptance") is not False:
        errors.append("bounded T01 smoke must not claim MVP physics acceptance")
    for key in ("absolute_capture_rate_particles_per_s", "occupation_validation",
                "lifetime_validation", "finite_age_validation"):
        if validation.get(key) is not None:
            errors.append(f"validation report must leave {key} unqualified")

    expected_cases = {
        f"{'capture' if mode == 'Capture' else 'ordinary'}-sigma{sigma:.0e}": (mode, sigma)
        for mode in ("Capture", "Parameter point")
        for sigma in benchmark["cross_sections_cm2"]
    }
    manifest_cases = _index_cases(manifest, "run manifest", errors)
    performance_cases = _index_cases(performance, "performance report", errors)
    if set(manifest_cases) != set(expected_cases):
        errors.append("run manifest case ids do not match the MVP mode/cross-section matrix")
    if validation.get("cases") != manifest.get("cases"):
        errors.append("validation cases differ from the run manifest cases")
    if set(performance_cases) != set(expected_cases):
        errors.append("performance case ids do not match the MVP mode/cross-section matrix")

    for case_id, (mode, sigma) in expected_cases.items():
        case = manifest_cases.get(case_id)
        if case is None:
            continue
        if case.get("mode") != mode or case.get("sigma_cm2") != sigma:
            errors.append(f"{case_id} mode or cross section differs from the MVP")
        if case.get("status") != "passed" or case.get("smoke_summary_complete") is not True:
            errors.append(f"{case_id} is not a completed smoke case")
        if case.get("physical_baseline_qualified") is not False:
            errors.append(f"{case_id} must remain unqualified as a physical baseline")
        contract = case.get("contract_checks")
        if not isinstance(contract, dict) or contract.get("runtime_model_verified") is not True:
            errors.append(f"{case_id} runtime model contract is not verified")
        elif mode == "Parameter point" and contract.get("artifact_mass_sigma_verified") is not True:
            errors.append(f"{case_id} artifact mass/cross-section contract is not verified")
        metrics = case.get("metrics")
        if not isinstance(metrics, dict):
            errors.append(f"{case_id} metrics must be an object")
            continue
        attempted = metrics.get("attempted")
        classified = metrics.get("classified")
        unresolved = metrics.get("unresolved_noncaptures")
        captured = metrics.get("captured")
        counts = (attempted, classified, unresolved, captured)
        if any(type(value) is not int or value < 0 for value in counts):
            errors.append(f"{case_id} has invalid trajectory counts")
        else:
            if attempted != benchmark["attempts_per_case"]:
                errors.append(f"{case_id} attempted count differs from the MVP")
            if classified + unresolved != attempted or captured > classified:
                errors.append(f"{case_id} trajectory accounting does not close")
            raw = metrics.get("capture_probability_raw")
            classified_probability = metrics.get("capture_probability_classified")
            if not _matches_ratio(raw, captured, attempted):
                errors.append(f"{case_id} raw capture probability is inconsistent")
            if not _matches_ratio(classified_probability, captured, classified):
                errors.append(f"{case_id} classified capture probability is inconsistent")

        timing = performance_cases.get(case_id)
        if timing is None:
            continue
        for report_key, case_key in (("status", "status"), ("wall_time_s", "wall_time_s"),
                                     ("attempted_trajectories", "attempted"),
                                     ("complete_evaporations", "complete_evaporations")):
            source = case if case_key in case else metrics
            if timing.get(report_key) != source.get(case_key):
                errors.append(f"{case_id} {report_key} differs between reports")
        wall_time = timing.get("wall_time_s")
        throughput = timing.get("attempted_trajectories_per_wall_s")
        if not isinstance(wall_time, (int, float)) or isinstance(wall_time, bool) or wall_time <= 0:
            errors.append(f"{case_id} has invalid wall time")
        elif type(attempted) is not int or not _matches_ratio(throughput, attempted, wall_time):
            errors.append(f"{case_id} throughput is inconsistent with count and wall time")

    if performance.get("total_wall_time_s") != manifest.get("total_wall_time_s"):
        errors.append("total wall time differs between performance report and run manifest")
    ctest_command = validation.get("ctest_command")
    if isinstance(ctest_command, dict) and performance.get("ctest_wall_time_s") != ctest_command.get("wall_time_s"):
        errors.append("CTest wall time differs between performance and validation reports")
    if manifest.get("analysis") != validation.get("analysis") or manifest.get("analysis") != performance.get("analysis"):
        errors.append("post-run analysis provenance differs between reports")
    for key in ("physical_absolute_capture_rate_particles_per_s", "collision_kernel_wall_time_s",
                "online_solve_wall_time_s"):
        if performance.get(key) is not None:
            errors.append(f"performance report must leave {key} unmeasured")
    return errors


def _project_path(path: Path, kind: str) -> Path:
    """Resolve a CLI input while refusing to leave this repository."""
    candidate = path if path.is_absolute() else PROJECT_ROOT / path
    candidate = candidate.absolute()
    try:
        candidate.relative_to(PROJECT_ROOT)
    except ValueError as error:
        raise ValueError(f"{kind} must be inside {PROJECT_ROOT}") from error
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(PROJECT_ROOT)
    except (OSError, ValueError) as error:
        raise ValueError(f"{kind} is missing or resolves outside {PROJECT_ROOT}") from error
    return resolved


def _read_bytes(path: Path) -> bytes:
    """Read one file without creating or modifying anything."""
    try:
        return path.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read {path}: {error}") from error


def _read_json_object(path: Path) -> dict[str, Any]:
    """Read one JSON object without creating or modifying any file."""
    try:
        value = json.loads(_read_bytes(path))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot parse {path}: {error}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def main(argv: list[str] | None = None) -> int:
    """Read and cross-check the saved T01 JSON bundle; never execute it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR,
                        help="Directory containing config.json and the three saved T01 reports")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                        help="Current MVP JSON to check against the frozen manifest config")
    args = parser.parse_args(argv)
    try:
        baseline_dir = _project_path(args.baseline_dir, "baseline directory")
        if not baseline_dir.is_dir():
            raise ValueError(f"baseline directory is not a directory: {baseline_dir}")
        config_path = _project_path(args.config, "MVP config")
        if not config_path.is_file():
            raise ValueError(f"MVP config is not a file: {config_path}")
        bundle_paths = {name: _project_path(baseline_dir / name, name) for name in BUNDLE_NAMES}
        config = _read_json_object(config_path)
        frozen_config_bytes = _read_bytes(bundle_paths["config.json"])
        reports = {name: _read_json_object(bundle_paths[name]) for name in BUNDLE_NAMES[1:]}
    except ValueError as error:
        parser.error(str(error))
    errors = validate_baseline_bundle(config, reports["run_manifest.json"],
                                      reports["validation_report.json"],
                                      reports["performance_report.json"],
                                      frozen_config_bytes)
    if errors:
        print("T01 baseline bundle validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    case_count = len(reports["run_manifest.json"]["cases"])
    test_count = len(reports["validation_report.json"]["ctest_tests"])
    print(f"T01 baseline bundle is valid ({case_count} cases, {test_count} final CTest results).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
