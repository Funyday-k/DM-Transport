#!/usr/bin/env python3
"""Run project-local T03 scientific checks and write a bounded JSON report.

The runner executes only the sampler built in this repository.  All speeds
are cm/s, temperatures kelvin, rates s^-1, and energy changes eV.  It never
uses the read-only reference source tree as an executable or oracle.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "Code/configs/validation/t03_oracle_contract.json"
DEFAULT_SAMPLER = ROOT / "build/t03_validation_sampler"
DEFAULT_REPORT = ROOT / "Output/Result/validation/validation_report.json"
VECTOR_MOMENTS = ("x", "y", "z", "xx", "xy", "xz", "yy", "yz", "zz")
COLLISION_SEED_MASK = 0x9E3779B9


def _reject_nonfinite_json_constant(token: str) -> None:
    """Reject NaN and infinities, which are not valid scientific JSON data."""
    raise ValueError(f"Nonfinite JSON value: {token}")


def _inside_project(path: Path) -> Path:
    """Resolve a path and reject source, executable or output paths outside this project."""
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(ROOT)
    except ValueError as error:
        raise ValueError(f"Path must remain inside the project: {path}") from error
    return resolved


def _finite_number(value: Any, label: str) -> float:
    """Read one finite JSON number; booleans are not measurements."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _count(value: Any, label: str) -> int:
    """Read one nonnegative integer count."""
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def _sample_variance_numerator(count: int, total: float, total_square: float,
                               label: str) -> float:
    """Reject inconsistent second moments, allowing only roundoff cancellation."""
    if count < 1:
        raise ValueError(f"{label} requires a positive sample count")
    centered = total_square - total * total / count
    roundoff = 64.0 * sys.float_info.epsilon * max(
        abs(total_square), abs(total * total / count), sys.float_info.min)
    if centered < -roundoff:
        raise ValueError(f"{label} has an inconsistent negative variance")
    return max(0.0, centered)


def _object_moments(records: list[dict[str, Any]], key: str,
                    component: int | None = None) -> dict[str, float]:
    """Pool a sampler `{sum,sum_sq}` observable, optionally one array element."""
    total_count = sum(_count(record["sample_count"], "sample_count")
                      for record in records)
    if total_count < 2:
        raise ValueError("At least two samples are required")
    total = 0.0
    total_square = 0.0
    for record in records:
        record_count = _count(record["sample_count"], "sample_count")
        entry = record[key]
        if component is not None:
            entry = entry[component]
        if not isinstance(entry, dict):
            raise ValueError(f"{key} must contain moment objects")
        record_sum = _finite_number(entry["sum"], f"{key}.sum")
        record_square = _finite_number(entry["sum_sq"], f"{key}.sum_sq")
        _sample_variance_numerator(record_count, record_sum, record_square, key)
        total += record_sum
        total_square += record_square
    mean = total / total_count
    variance = _sample_variance_numerator(
        total_count, total, total_square, key) / (total_count - 1)
    return {"sample_count": total_count, "mean": mean,
            "standard_error": math.sqrt(variance / total_count)}


def _normal_critical(familywise_alpha: float, two_sided_tests: int) -> float:
    """Return a Bonferroni normal critical value for declared scalar checks."""
    if not 0.0 < familywise_alpha < 1.0 or two_sided_tests < 1:
        raise ValueError("Invalid familywise alpha or test count")
    return statistics.NormalDist().inv_cdf(
        1.0 - familywise_alpha / (2.0 * two_sided_tests))


def _equivalence(moment: dict[str, float], critical: float,
                 tolerance: float) -> dict[str, Any]:
    """Require the simultaneous confidence interval to fit inside ±tolerance."""
    halfwidth = critical * moment["standard_error"]
    return {
        **moment,
        "simultaneous_halfwidth": halfwidth,
        "absolute_tolerance": tolerance,
        "status": "pass" if abs(moment["mean"]) + halfwidth <= tolerance
        else "fail",
    }


def _run_sampler(sampler: Path, mode: str, radius_Rsun: float,
                 parameter: float, sample_count: int, seed: int,
                 direction: tuple[float, float, float] | None = None,
                 timeout_seconds: int = 600) -> dict[str, Any]:
    """Run one bounded project-local Monte Carlo batch and parse its JSON."""
    command = [
        str(sampler),
        str(ROOT / "Code/data/solar/model_agss09.dat"),
        str(ROOT / "Code/data/solar/Nuclear_Data.txt"),
        mode, repr(radius_Rsun), repr(parameter), str(sample_count), str(seed),
    ]
    if direction is not None:
        command.extend(("--direction", ",".join(repr(x) for x in direction)))
    completed = subprocess.run(
        command, cwd=ROOT, text=True, capture_output=True,
        timeout=timeout_seconds, check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"T03 sampler failed ({mode}, r/Rsun={radius_Rsun}, seed={seed}): "
            f"{completed.stderr.strip()[:1500]}")
    try:
        result = json.loads(completed.stdout,
                            parse_constant=_reject_nonfinite_json_constant)
    except json.JSONDecodeError as error:
        raise ValueError("T03 sampler did not return one JSON object") from error
    if not isinstance(result, dict) or result.get("mode") != mode:
        raise ValueError("T03 sampler returned a mismatched mode")
    if (result.get("schema_version") != 1 or
            result.get("producer") != "DM-Transport/t03_validation_sampler" or
            result.get("qualification") !=
            "project_physics_validation_not_independent_oracle"):
        raise ValueError("T03 sampler provenance or schema is not recognized")
    if _count(result.get("sample_count"), "sample_count") != sample_count:
        raise ValueError("T03 sampler returned a mismatched sample count")
    if not math.isclose(_finite_number(result.get("radius_Rsun"), "radius_Rsun"),
                        radius_Rsun, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("T03 sampler returned a mismatched radius")
    if (result.get("seed_incoming") != seed or
            result.get("seed_collision") != (seed ^ COLLISION_SEED_MASK)):
        raise ValueError("T03 sampler returned a mismatched RNG seed")
    if (result.get("dark_matter_mass_GeV") != 0.1 or
            result.get("proton_cross_section_cm2") != 1e-34):
        raise ValueError("T03 sampler returned a different physics model")
    parameter_key = "tchi_over_t" if mode == "thermal" else "fixed_speed_over_escape"
    if not math.isclose(_finite_number(result.get(parameter_key), parameter_key),
                        parameter, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("T03 sampler returned a mismatched state parameter")
    if mode == "fixed":
        if direction is None:
            raise ValueError("Fixed state requires an incoming direction")
        observed = result.get("direction_unit")
        expected_length = math.sqrt(sum(value * value for value in direction))
        expected = [value / expected_length for value in direction]
        if (not isinstance(observed, list) or len(observed) != 3 or
                any(not math.isclose(
                    _finite_number(actual, "direction_unit"), target,
                    rel_tol=0.0, abs_tol=2e-12)
                    for actual, target in zip(observed, expected))):
            raise ValueError("T03 sampler returned a mismatched direction")
    return result


def _sum_counts(records: list[dict[str, Any]], key: str) -> list[int]:
    """Pool one fixed-length categorical histogram over independent batches."""
    first = records[0].get(key)
    if not isinstance(first, list):
        raise ValueError(f"{key} must be a list")
    pooled = [0] * len(first)
    for record in records:
        values = record.get(key)
        if not isinstance(values, list) or len(values) != len(pooled):
            raise ValueError(f"{key} histogram shapes differ")
        for index, value in enumerate(values):
            pooled[index] += _count(value, f"{key}[{index}]")
    return pooled


def _wilson_interval(successes: int, trials: int,
                     z_value: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson binomial interval for an event probability, including zero tails."""
    if trials <= 0 or not 0 <= successes <= trials:
        raise ValueError("Wilson interval requires 0 <= successes <= trials")
    proportion = successes / trials
    denominator = 1.0 + z_value * z_value / trials
    center = (proportion + z_value * z_value / (2.0 * trials)) / denominator
    spread = z_value / denominator * math.sqrt(
        proportion * (1.0 - proportion) / trials +
        z_value * z_value / (4.0 * trials * trials))
    lower = 0.0 if successes == 0 else max(0.0, center - spread)
    upper = 1.0 if successes == trials else min(1.0, center + spread)
    return lower, upper


def _rotation_components(records: list[dict[str, Any]],
                         key: str) -> list[list[int]]:
    """Pool three canonical Cartesian histograms and check every total."""
    bin_count = len(records[0]["canonical_component_histogram_cutpoints"]) + 1
    pooled = [[0] * bin_count for _ in range(3)]
    for record in records:
        matrix = record.get(key)
        if (not isinstance(matrix, list) or len(matrix) != 3 or
                any(not isinstance(row, list) or len(row) != bin_count
                    for row in matrix)):
            raise ValueError(f"{key} must be a 3x{bin_count} histogram")
        for component, row in enumerate(matrix):
            if sum(_count(value, key) for value in row) != record["sample_count"]:
                raise ValueError(f"{key} component total differs from sample count")
            for bin_index, value in enumerate(row):
                pooled[component][bin_index] += value
    return pooled


def _difference_equivalence(first: dict[str, float],
                            second: dict[str, float], critical: float,
                            tolerance: float) -> dict[str, Any]:
    """Compare two disjoint-sample estimates with simultaneous equivalence."""
    difference = second["mean"] - first["mean"]
    standard_error = math.hypot(first["standard_error"],
                                second["standard_error"])
    halfwidth = critical * standard_error
    return {
        "base_mean": first["mean"], "rotated_mean": second["mean"],
        "difference": difference, "standard_error": standard_error,
        "simultaneous_halfwidth": halfwidth,
        "absolute_tolerance": tolerance,
        "status": "pass" if abs(difference) + halfwidth <= tolerance
        else "fail",
    }


def _cdf_equivalence(first_count: int, first_total: int,
                     second_count: int, second_total: int,
                     critical: float, tolerance: float) -> dict[str, Any]:
    """Compare two independent empirical CDF ordinates."""
    first = first_count / first_total
    second = second_count / second_total
    first_se_sq = first * (1.0 - first) / first_total
    second_se_sq = second * (1.0 - second) / second_total
    return _difference_equivalence(
        {"mean": first, "standard_error": math.sqrt(first_se_sq)},
        {"mean": second, "standard_error": math.sqrt(second_se_sq)},
        critical, tolerance)


def _check_target_totals(records: list[dict[str, Any]]) -> None:
    """Verify that every generated event selected exactly one of 63 targets."""
    for record in records:
        counts = record.get("target_counts")
        if (not isinstance(counts, list) or len(counts) != 63 or
                sum(_count(value, "target_counts") for value in counts) !=
                record["sample_count"]):
            raise ValueError("Target-selection counts do not sum to sample count")


def _thermal_equilibrium(records: list[dict[str, Any]],
                         suite: dict[str, Any], critical: float) -> dict[str, Any]:
    """Test the stationary continuous-time weak form for one solar radius."""
    _check_target_totals(records)
    expected_edges = suite["speed_bin_edges_over_bath_dm_thermal_speed"][:-1]
    if any(record["weak_speed_bin_lower_edges_over_bath_dm_thermal_speed"]
           != expected_edges for record in records):
        raise ValueError("Sampler speed-bin edges differ from the frozen contract")
    tolerances = suite["acceptance"]
    bins = [
        _equivalence(
            _object_moments(records, "weak_speed_bin_moments", index),
            critical,
            _finite_number(tolerances[
                "speed_bin_max_absolute_mean_plus_bonferroni_se"],
                "speed-bin tolerance"))
        for index in range(len(expected_edges))
    ]
    energy = _equivalence(
        _object_moments(records, "weak_energy"), critical,
        _finite_number(tolerances[
            "energy_max_absolute_mean_plus_bonferroni_se"],
            "energy tolerance"))
    rate_ratio = _object_moments(records, "rate_ratio")
    gamma_reference = _finite_number(records[0]["gamma_ref_s_inv"],
                                     "gamma_ref_s_inv")
    analytic_rate = _finite_number(
        records[0]["analytic_mb_mean_total_rate_s_inv"],
        "analytic_mb_mean_total_rate_s_inv")
    if gamma_reference <= 0.0 or analytic_rate <= 0.0:
        raise ValueError("Thermal reference and analytic rates must be positive")
    for record in records[1:]:
        if (not math.isclose(record["gamma_ref_s_inv"], gamma_reference,
                             rel_tol=1e-12) or
                not math.isclose(record[
                    "analytic_mb_mean_total_rate_s_inv"], analytic_rate,
                    rel_tol=1e-12)):
            raise ValueError("Thermal normalization differs across seeds")
    rate_relative = {
        "sample_count": rate_ratio["sample_count"],
        "mean": rate_ratio["mean"] * gamma_reference / analytic_rate - 1.0,
        "standard_error": rate_ratio["standard_error"] *
        gamma_reference / analytic_rate,
    }
    rate_check = _equivalence(
        rate_relative, critical,
        _finite_number(tolerances[
            "analytic_mb_rate_max_relative_difference_plus_bonferroni_se"],
            "analytic MB rate tolerance"))
    rate_check["sampled_mean_rate_s_inv"] = (
        rate_ratio["mean"] * gamma_reference)
    rate_check["analytic_mean_rate_s_inv"] = analytic_rate
    return {
        "radius_Rsun": records[0]["radius_Rsun"],
        "sample_count": sum(record["sample_count"] for record in records),
        "seeds": [record["seed_incoming"] for record in records],
        "speed_bin_weak_residuals": bins,
        "energy_weak_residual": energy,
        "analytic_mb_mean_rate": rate_check,
        "status": "pass" if all(item["status"] == "pass"
                                for item in bins + [energy, rate_check])
        else "fail",
    }


def _heating_cooling(records: list[dict[str, Any]],
                     temperature_ratio: float, suite: dict[str, Any],
                     critical: float) -> dict[str, Any]:
    """Require a resolved rate-weighted drift toward the bath temperature."""
    _check_target_totals(records)
    moment = _object_moments(records, "weak_energy")
    halfwidth = critical * moment["standard_error"]
    acceptance = suite["acceptance"]
    if temperature_ratio < 1.0:
        margin = _finite_number(
            acceptance["cold_minimum_lower_confidence_bound"],
            "cold drift margin")
        passes = moment["mean"] - halfwidth > margin
    elif temperature_ratio > 1.0:
        margin = _finite_number(
            acceptance["hot_maximum_upper_confidence_bound"],
            "hot drift margin")
        passes = moment["mean"] + halfwidth < margin
    else:
        raise ValueError("Tchi=T is handled by thermal equilibrium")
    return {"tchi_over_t": temperature_ratio, **moment,
            "simultaneous_halfwidth": halfwidth, "required_margin": margin,
            "status": "pass" if passes else "fail"}


def _full_chain_rotation(base: list[dict[str, Any]],
                         rotated: list[dict[str, Any]],
                         suite: dict[str, Any], critical: float) -> dict[str, Any]:
    """Compare inverse-rotated full-chain vector CDFs and moments."""
    _check_target_totals(base + rotated)
    if set(record["seed_collision"] for record in base) & set(
            record["seed_collision"] for record in rotated):
        raise ValueError("Rotation samples must use disjoint RNG seeds")
    cutpoints = suite["component_cdf_thresholds_over_bath_dm_thermal_speed"]
    if any(record["canonical_component_histogram_cutpoints"] != cutpoints
           for record in base + rotated):
        raise ValueError("Rotation component cutpoints differ from contract")
    for record in base + rotated:
        if not math.isclose(
                _finite_number(record["canonical_component_scale_cm_s"],
                               "canonical_component_scale_cm_s"),
                _finite_number(base[0]["canonical_component_scale_cm_s"],
                               "canonical_component_scale_cm_s"),
                rel_tol=1e-12):
            raise ValueError("Rotation component scale differs across runs")
    acceptance = suite["acceptance"]
    cdf_tolerance = _finite_number(
        acceptance["component_cdf_max_absolute_difference_plus_bonferroni_se"],
        "rotation CDF tolerance")
    cdf_checks = []
    for vector_name, key in (
            ("target_velocity", "canonical_target_component_histogram_counts"),
            ("outgoing_dm_velocity",
             "canonical_outgoing_component_histogram_counts")):
        base_counts = _rotation_components(base, key)
        rotated_counts = _rotation_components(rotated, key)
        base_total = sum(record["sample_count"] for record in base)
        rotated_total = sum(record["sample_count"] for record in rotated)
        for component in range(3):
            for threshold in range(len(cutpoints)):
                result = _cdf_equivalence(
                    sum(base_counts[component][:threshold + 1]), base_total,
                    sum(rotated_counts[component][:threshold + 1]),
                    rotated_total, critical, cdf_tolerance)
                result.update({"vector": vector_name,
                               "component": "xyz"[component],
                               "threshold_over_bath_dm_thermal_speed":
                               cutpoints[threshold]})
                cdf_checks.append(result)
    moment_checks = []
    for name in VECTOR_MOMENTS:
        first = _named_vector_moment(base, name)
        second = _named_vector_moment(rotated, name)
        tolerance_key = (
            "component_mean_max_absolute_difference_plus_bonferroni_se"
            if len(name) == 1 else
            "second_moment_max_absolute_difference_plus_bonferroni_se")
        check = _difference_equivalence(
            first, second, critical,
            _finite_number(acceptance[tolerance_key], tolerance_key))
        check["observable"] = name
        moment_checks.append(check)
    return {
        "radius_Rsun": base[0]["radius_Rsun"],
        "base_direction": base[0]["direction_unit"],
        "rotated_direction": rotated[0]["direction_unit"],
        "sample_count_per_direction": sum(record["sample_count"] for record in base),
        "component_cdf_checks": cdf_checks,
        "outgoing_vector_moment_checks": moment_checks,
        "target_counts_base": _sum_counts(base, "target_counts"),
        "target_counts_rotated": _sum_counts(rotated, "target_counts"),
        "status": "pass" if all(item["status"] == "pass"
                                for item in cdf_checks + moment_checks)
        else "fail",
    }


def _named_vector_moment(records: list[dict[str, Any]],
                         name: str) -> dict[str, float]:
    """Pool one canonical outgoing-vector moment and its square."""
    key = "canonical_outgoing_vector_over_escape_moments"
    total_count = sum(_count(record["sample_count"], "sample_count")
                      for record in records)
    if total_count < 2:
        raise ValueError("At least two vector samples are required")
    total = 0.0
    total_sq = 0.0
    for record in records:
        record_sum = _finite_number(record[key][name]["sum"],
                                    f"{key}.{name}.sum")
        record_sq = _finite_number(record[key][name]["sum_sq"],
                                   f"{key}.{name}.sum_sq")
        _sample_variance_numerator(record["sample_count"], record_sum,
                                   record_sq, f"{key}.{name}")
        total += record_sum
        total_sq += record_sq
    mean = total / total_count
    variance = _sample_variance_numerator(
        total_count, total, total_sq, f"{key}.{name}") / (total_count - 1)
    return {"sample_count": total_count, "mean": mean,
            "standard_error": math.sqrt(variance / total_count)}


def _cm_angle_distribution(records: list[dict[str, Any]],
                           suite: dict[str, Any],
                           critical: float) -> dict[str, Any]:
    """Test full-chain CM scattering cosines against the isotropic CDF."""
    _check_target_totals(records)
    edges = suite["cosine_histogram_edges"]
    if len(edges) != 11 or edges[0] != -1.0 or edges[-1] != 1.0:
        raise ValueError("CM-angle contract must declare ten bins on [-1,1]")
    for record in records:
        observed_edges = record["cm_scattering_cosine_histogram_edges"]
        if (len(observed_edges) != len(edges) or
                any(not math.isclose(
                    _finite_number(actual, "CM-angle edge"), expected,
                    rel_tol=0.0, abs_tol=1e-14)
                    for actual, expected in zip(observed_edges, edges))):
            raise ValueError("CM-angle histogram edges differ from contract")
        if sum(_count(value, "CM-angle count") for value in record[
                "cm_scattering_cosine_histogram_counts"]) != record[
                    "sample_count"]:
            raise ValueError("CM-angle histogram count differs from samples")
    counts = _sum_counts(records, "cm_scattering_cosine_histogram_counts")
    if len(counts) != len(edges) - 1:
        raise ValueError("CM-angle histogram shape differs from contract")
    total = sum(record["sample_count"] for record in records)
    tolerance = _finite_number(suite["acceptance"][
        "cdf_max_absolute_difference_plus_bonferroni_se"],
        "CM-angle CDF tolerance")
    checks = []
    for boundary in range(1, len(counts)):
        observed = sum(counts[:boundary]) / total
        expected = 0.5 * (edges[boundary] + 1.0)
        moment = {
            "sample_count": total,
            "mean": observed - expected,
            "standard_error": math.sqrt(observed * (1.0 - observed) / total),
        }
        check = _equivalence(moment, critical, tolerance)
        check.update({"cosine_threshold": edges[boundary],
                      "observed_cdf": observed,
                      "uniform_cdf": expected})
        checks.append(check)
    return {
        "radius_Rsun": records[0]["radius_Rsun"],
        "sample_count": total,
        "seeds": [record["seed_collision"] for record in records],
        "cdf_checks": checks,
        "status": "pass" if all(check["status"] == "pass"
                                for check in checks) else "fail",
    }


def _near_escape_tail(records: list[dict[str, Any]],
                      suite: dict[str, Any]) -> dict[str, Any]:
    """Report local threshold crossings, uncertainty and numerical-domain tails."""
    _check_target_totals(records)
    total = sum(record["sample_count"] for record in records)
    crossing = sum(_count(record["tail_total_count"], "tail_total_count")
                   for record in records)
    first_total = sum(_count(record["tail_first_half_sample_count"],
                             "tail_first_half_sample_count")
                      for record in records)
    first_crossing = sum(_count(record["tail_first_half_count"],
                                "tail_first_half_count")
                         for record in records)
    target_counts = _sum_counts(records, "target_counts")
    tail_target_counts = _sum_counts(records, "tail_target_counts")
    if (sum(tail_target_counts) != crossing or crossing > total or
            first_crossing > crossing or first_crossing > first_total):
        raise ValueError("Near-escape tail counts are internally inconsistent")
    overflow = {
        "1.5": sum(_count(record["domain_overflow_count_1p5_vesc"],
                          "overflow 1.5") for record in records),
        "2.0": sum(_count(record["domain_overflow_count_2p0_vesc"],
                          "overflow 2.0") for record in records),
        "3.0": sum(_count(record["domain_overflow_count_3p0_vesc"],
                          "overflow 3.0") for record in records),
    }
    if not (0 <= overflow["3.0"] <= overflow["2.0"] <=
            overflow["1.5"] <= crossing):
        raise ValueError("Passive velocity-ceiling counts are not monotone")
    lower, upper = _wilson_interval(crossing, total)
    proportion = crossing / total
    halfwidth = 0.5 * (upper - lower)
    minimum = _count(suite["precision_status"]["minimum_crossing_count"],
                     "minimum_crossing_count")
    maximum_relative = _finite_number(
        suite["precision_status"]["maximum_wilson_relative_halfwidth"],
        "maximum_wilson_relative_halfwidth")
    heterogeneity_z_threshold = _finite_number(
        suite["precision_status"]["heterogeneity_flag_z_threshold"],
        "heterogeneity_flag_z_threshold")
    precision_sufficient = (crossing >= minimum and
                            halfwidth / proportion <= maximum_relative)
    per_seed = []
    for record in records:
        n = record["sample_count"]
        count = _count(record["tail_total_count"], "tail_total_count")
        interval = _wilson_interval(count, n)
        per_seed.append({"seed": record["seed_collision"],
                         "sample_count": n, "crossing_count": count,
                         "probability": count / n,
                         "wilson_95_interval": interval})
    second_total = total - first_total
    second_crossing = crossing - first_crossing
    first_probability = first_crossing / first_total
    second_probability = second_crossing / second_total
    stability = {
        "first_half": {"sample_count": first_total,
                       "crossing_count": first_crossing,
                       "probability": first_probability,
                       "wilson_95_interval": _wilson_interval(
                           first_crossing, first_total)},
        "second_half": {"sample_count": second_total,
                        "crossing_count": second_crossing,
                        "probability": second_probability,
                        "wilson_95_interval": _wilson_interval(
                            second_crossing, second_total)},
    }
    if crossing >= minimum:
        pooled_variance = proportion * (1.0 - proportion)
        half_standard_error = math.sqrt(
            pooled_variance * (1.0 / first_total + 1.0 / second_total))
        stability["half_sample_z"] = (
            abs(first_probability - second_probability) / half_standard_error
            if half_standard_error > 0.0 else 0.0)
        seed_z = []
        for record in records:
            n = record["sample_count"]
            c = record["tail_total_count"]
            other_n = total - n
            if other_n == 0:
                continue  # One-seed smoke cannot assess seed heterogeneity.
            other_c = crossing - c
            standard_error = math.sqrt(
                pooled_variance * (1.0 / n + 1.0 / other_n))
            seed_z.append(
                abs(c / n - other_c / other_n) / standard_error
                if standard_error > 0.0 else 0.0)
        stability["maximum_seed_vs_rest_z"] = max(seed_z) if seed_z else None
        stability["heterogeneity_z_threshold"] = heterogeneity_z_threshold
        stability["within_declared_z_threshold"] = (
            stability["half_sample_z"] <= heterogeneity_z_threshold and
            stability["maximum_seed_vs_rest_z"] <= heterogeneity_z_threshold
        ) if seed_z else None
    else:
        stability["heterogeneity_z_threshold"] = heterogeneity_z_threshold
        stability["within_declared_z_threshold"] = None
    major_targets = []
    for index in sorted(range(len(target_counts)),
                        key=lambda value: target_counts[value],
                        reverse=True)[:8]:
        selected = target_counts[index]
        tail_selected = tail_target_counts[index]
        major_targets.append({
            "target_index": index, "selected_count": selected,
            "crossing_count": tail_selected,
            "conditional_crossing_probability":
                tail_selected / selected if selected > 0 else None,
            "wilson_95_interval": _wilson_interval(tail_selected, selected)
            if selected > 0 else None,
        })
    return {
        "fixed_speed_over_escape": records[0]["fixed_speed_over_escape"],
        "sample_count": total, "crossing_count": crossing,
        "crossing_probability": proportion,
        "wilson_95_interval": [lower, upper],
        "wilson_relative_halfwidth": halfwidth / proportion
        if crossing > 0 else None,
        "precision_status": "precision_sufficient" if precision_sufficient
        else suite["precision_status"]["insufficient_count_status"],
        "stability": stability,
        "per_seed": per_seed,
        "target_counts": target_counts,
        "tail_target_counts": tail_target_counts,
        "major_target_conditional_probabilities": major_targets,
        "passive_velocity_ceiling_exceedance_counts": overflow,
        "crossing_probability_if_truncated_at_ceiling": {
            ceiling: (crossing - count) / total
            for ceiling, count in overflow.items()},
        "interpretation": "One local collision; not an evaporation rate",
    }


def _write_report(path: Path, report: dict[str, Any]) -> None:
    """Atomically replace a generated report inside this project's Output tree."""
    output_root = ROOT / "Output"
    path = _inside_project(path)
    try:
        path.relative_to(output_root)
    except ValueError as error:
        raise ValueError("T03 report must be written under Output/") from error
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    temporary.replace(path)


def _mark_smoke_diagnostics(value: Any) -> None:
    """Prevent reduced-sample subchecks from masquerading as full acceptance."""
    if isinstance(value, dict):
        if value.get("status") in {"pass", "fail"}:
            value["diagnostic_full_bound_result"] = value["status"]
            value["status"] = "not_evaluated"
        for child in value.values():
            _mark_smoke_diagnostics(child)
    elif isinstance(value, list):
        for child in value:
            _mark_smoke_diagnostics(child)


def _contract_qualification(path: Path) -> dict[str, Any]:
    """Only a clean, tracked canonical contract can produce G0 evidence."""
    relative = path.relative_to(ROOT).as_posix()
    if path != DEFAULT_CONTRACT:
        return {"path": relative, "status": "diagnostic_noncanonical_contract"}
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", "--", relative], cwd=ROOT,
        text=True, capture_output=True, check=False)
    state = subprocess.run(
        ["git", "status", "--porcelain", "--", relative], cwd=ROOT,
        text=True, capture_output=True, check=False)
    if tracked.returncode != 0 or state.returncode != 0:
        raise RuntimeError("Cannot verify the canonical T03 contract is tracked")
    return {
        "path": relative,
        "status": "version_controlled_clean" if not state.stdout.strip()
        else "diagnostic_uncommitted_contract",
    }


def _sampler_qualification(path: Path) -> dict[str, Any]:
    """Record whether this is a fresh local CMake target from clean T03 code."""
    relative = path.relative_to(ROOT).as_posix()
    result: dict[str, Any] = {"path": relative}
    if path.name != "t03_validation_sampler":
        result["status"] = "diagnostic_noncanonical_sampler"
        return result
    cache = path.parent / "CMakeCache.txt"
    if not cache.is_file():
        result["status"] = "diagnostic_missing_project_build_cache"
        return result
    source_lines = [line.partition("=")[2].strip()
                    for line in cache.read_text(encoding="utf-8").splitlines()
                    if line.startswith("CMAKE_HOME_DIRECTORY:INTERNAL=")]
    if len(source_lines) != 1 or Path(source_lines[0]).resolve() != ROOT:
        result["status"] = "diagnostic_nonproject_build"
        return result
    sources = (
        "CMakeLists.txt",
        "Code/tools/t03_validation_sampler.cpp",
        "Code/python/run_t03_validation.py",
        "Code/src/physics/ScatteringPhysics.cpp",
        "Code/src/physics/SolarBackground.cpp",
        "Code/include/transport/physics/ScatteringPhysics.hpp",
        "Code/include/transport/physics/SolarBackground.hpp",
    )
    state = subprocess.run(
        ["git", "status", "--porcelain", "--", *sources], cwd=ROOT,
        text=True, capture_output=True, check=False)
    if state.returncode != 0:
        raise RuntimeError("Cannot inspect the project T03 source state")
    if state.stdout.strip():
        result["status"] = "diagnostic_uncommitted_sampler_source"
        return result
    if path.stat().st_mtime_ns < max((ROOT / source).stat().st_mtime_ns
                                    for source in sources if source !=
                                    "Code/python/run_t03_validation.py"):
        result["status"] = "diagnostic_stale_sampler_binary"
        return result
    result["status"] = "project_local_clean_fresh_build"
    return result


def _fast_prerequisite(sampler: Path) -> dict[str, Any]:
    """Require the existing fast CTest suite in the sampler's local build tree."""
    completed = subprocess.run(
        ["ctest", "--test-dir", str(sampler.parent), "-L", "fast",
         "--output-on-failure"],
        cwd=ROOT, text=True, capture_output=True, timeout=120, check=False)
    combined = completed.stdout + completed.stderr
    if completed.returncode != 0 or "100% tests passed" not in combined:
        raise RuntimeError("Fast prerequisite CTest failed: " +
                           combined.strip()[-2500:])
    return {"status": "pass", "runner": "ctest -L fast",
            "build_directory": sampler.parent.name,
            "result_excerpt": next(
                (line.strip() for line in combined.splitlines()
                 if "100% tests passed" in line), "")}


def _declared_family_size(suite: dict[str, Any]) -> int:
    """Count every predeclared scalar equivalence or sign check once."""
    equilibrium = suite["thermal_equilibrium"]
    rotation = suite["full_chain_rotation"]
    thermal_count = len(equilibrium["radius_Rsun"]) * (
        len(equilibrium["speed_bin_edges_over_bath_dm_thermal_speed"]) - 1
        + 1 + 1)
    heat_count = sum(value != 1.0 for value in suite[
        "heating_cooling"]["tchi_over_t"])
    rotation_count = len(rotation["radius_Rsun"]) * (
        len(rotation["vectors"]) * 3 *
        len(rotation[
            "component_cdf_thresholds_over_bath_dm_thermal_speed"])
        + len(rotation["moments"]))
    angle = suite["cm_scattering_angle"]
    angle_count = len(angle["radius_Rsun"]) * (
        len(angle["cosine_histogram_edges"]) - 2)
    return thermal_count + heat_count + rotation_count + angle_count


def run_validation(contract: dict[str, Any], sampler: Path,
                   contract_qualification: dict[str, Any],
                   sampler_qualification: dict[str, Any],
                   smoke: bool = False,
                   smoke_samples_per_seed: int = 1000) -> dict[str, Any]:
    """Execute the frozen T03 state matrix and assemble a scientific report."""
    suite = contract["scientific_physics_validation"]
    if suite["mode"] != "full" or contract["schema_version"] != 1:
        raise ValueError("Unsupported T03 scientific contract")
    common = suite["common"]
    declared_tests = _count(common["bonferroni_scalar_count"],
                            "bonferroni_scalar_count")
    actual_tests = _declared_family_size(suite)
    if actual_tests != declared_tests:
        raise ValueError(
            f"Declared Bonferroni family has {declared_tests}, "
            f"but state matrix has {actual_tests} scalar checks")
    full_samples = _count(common["samples_per_seed"], "samples_per_seed")
    full_seeds = common["seeds"]
    pilot_seeds = common["pilot_seeds_excluded_from_full"]
    if (full_samples < 2 or not isinstance(full_seeds, list) or
            len(full_seeds) < 2 or len(set(full_seeds)) != len(full_seeds) or
            any(type(seed) is not int or not 0 <= seed <= 2**32 - 1
                for seed in full_seeds) or
            not isinstance(pilot_seeds, list) or not pilot_seeds or
            any(type(seed) is not int or not 0 <= seed <= 2**32 - 1
                for seed in pilot_seeds)):
        raise ValueError("Invalid frozen scientific sample count or seeds")
    rotation_offset = _count(suite["full_chain_rotation"]["rotated_seed_offset"],
                             "rotated_seed_offset")
    full_stream_seeds = set(full_seeds) | {
        seed + rotation_offset for seed in full_seeds}
    pilot_stream_seeds = set(pilot_seeds) | {
        seed + rotation_offset for seed in pilot_seeds}
    if full_stream_seeds & pilot_stream_seeds:
        raise ValueError("Pilot and full RNG seed streams overlap")
    if smoke:
        if not 2 <= smoke_samples_per_seed < full_samples:
            raise ValueError("Smoke samples must be in [2, full samples)")
        sample_count = smoke_samples_per_seed
        seeds = pilot_seeds[-1:]
    else:
        sample_count = full_samples
        seeds = full_seeds
    confidence_level = _finite_number(common["confidence_level"],
                                      "confidence_level")
    critical = _normal_critical(1.0 - confidence_level, declared_tests)
    fast = _fast_prerequisite(sampler)
    started = time.monotonic()
    cache: dict[tuple[Any, ...], dict[str, Any]] = {}

    def get_batch(mode: str, radius: float, parameter: float, seed: int,
                  direction: tuple[float, float, float] | None = None
                  ) -> dict[str, Any]:
        key = (mode, radius, parameter, seed, direction)
        if key not in cache:
            print(f"T03 {mode}: r/Rsun={radius}, parameter={parameter}, "
                  f"seed={seed}, samples={sample_count}",
                  file=sys.stderr, flush=True)
            cache[key] = _run_sampler(
                sampler, mode, radius, parameter, sample_count, seed,
                direction=direction)
        return cache[key]

    equilibrium_suite = suite["thermal_equilibrium"]
    equilibrium_states = []
    for radius in equilibrium_suite["radius_Rsun"]:
        records = [get_batch("thermal", radius,
                             equilibrium_suite["tchi_over_t"], seed)
                   for seed in seeds]
        equilibrium_states.append(_thermal_equilibrium(
            records, equilibrium_suite, critical))

    heating_suite = suite["heating_cooling"]
    heating_states = []
    for temperature_ratio in heating_suite["tchi_over_t"]:
        if temperature_ratio == 1.0:
            matching = [state for state in equilibrium_states
                        if state["radius_Rsun"] == heating_suite["radius_Rsun"]]
            if len(matching) != 1:
                raise ValueError("Equilibrium heating state cannot be reused")
            heating_states.append({
                "tchi_over_t": 1.0,
                "status": matching[0]["energy_weak_residual"]["status"],
                "reused_thermal_equilibrium_radius_Rsun":
                    heating_suite["radius_Rsun"],
            })
        else:
            records = [get_batch("thermal", heating_suite["radius_Rsun"],
                                 temperature_ratio, seed) for seed in seeds]
            heating_states.append(_heating_cooling(
                records, temperature_ratio, heating_suite, critical))

    rotation_suite = suite["full_chain_rotation"]
    directions = [tuple(_finite_number(x, "direction component")
                        for x in direction)
                  for direction in rotation_suite["incoming_directions_lab"]]
    if len(directions) != 2 or any(len(direction) != 3 for direction in directions):
        raise ValueError("Rotation requires two three-component directions")
    offset = rotation_offset
    if set(seeds) & {seed + offset for seed in seeds}:
        raise ValueError("Rotation seed sets overlap")
    rotation_states = []
    for radius in rotation_suite["radius_Rsun"]:
        base = [get_batch("fixed", radius,
                          rotation_suite["speed_over_escape"], seed,
                          directions[0]) for seed in seeds]
        rotated = [get_batch("fixed", radius,
                             rotation_suite["speed_over_escape"], seed + offset,
                             directions[1]) for seed in seeds]
        rotation_states.append(_full_chain_rotation(
            base, rotated, rotation_suite, critical))

    angle_suite = suite["cm_scattering_angle"]
    angle_states = []
    for radius in angle_suite["radius_Rsun"]:
        records = [get_batch("fixed", radius,
                             angle_suite["speed_over_escape"], seed,
                             directions[0]) for seed in seeds]
        angle_states.append(_cm_angle_distribution(
            records, angle_suite, critical))

    tail_suite = suite["near_escape_tail"]
    tail_states = []
    for speed_ratio in tail_suite["speed_over_escape"]:
        records = [get_batch("fixed", tail_suite["radius_Rsun"],
                             speed_ratio, seed, directions[0])
                   for seed in seeds]
        tail_states.append(_near_escape_tail(records, tail_suite))

    scientific_groups_pass = all(
        state["status"] == "pass"
        for state in equilibrium_states + heating_states + rotation_states +
        angle_states)
    qualifying_contract = (
        contract_qualification["status"] == "version_controlled_clean")
    qualifying_sampler = (
        sampler_qualification["status"] == "project_local_clean_fresh_build")
    qualifying_inputs = qualifying_contract and qualifying_sampler
    physics_status = (
        common["smoke_override_policy"]["physics_validation_status"]
        if smoke else "not_evaluated" if not qualifying_inputs
        else "pass" if scientific_groups_pass else "fail")
    parity_policy = contract["report_status_policy"]["reference_parity"]
    if parity_policy["status"] != "not_evaluated":
        raise ValueError("External reference artifact evaluation is not implemented")
    report = {
        "schema_version": 1,
        "contract_id": contract["contract_id"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_mode": "smoke" if smoke else "full",
        "run_configuration": {
            "contract": contract_qualification,
            "sampler": sampler_qualification,
            "seeds": seeds,
            "samples_per_seed_per_state": sample_count,
            "unique_sampler_batches": len(cache),
            "bonferroni_scalar_count": declared_tests,
            "bonferroni_critical_value": critical,
            "confidence_level": confidence_level,
            "elapsed_seconds": time.monotonic() - started,
        },
        "fast_prerequisite": fast,
        "physics_validation": {
            "status": physics_status,
            "g0_eligible": not smoke and qualifying_inputs and
            scientific_groups_pass,
            "existing_elastic_kinematics": {
                "status": fast["status"], "evidence": "fast CTest single_collision"},
            "existing_target_conditioning": {
                "status": fast["status"], "evidence": "fast CTest target_sampling"},
            "existing_scattering_angle_distribution": {
                "status": fast["status"], "evidence": "fast CTest single_collision"},
            "thermal_equilibrium": {
                "status": "pass" if all(state["status"] == "pass"
                                        for state in equilibrium_states) else "fail",
                "states": equilibrium_states,
            },
            "heating_cooling": {
                "status": "pass" if all(state["status"] == "pass"
                                        for state in heating_states) else "fail",
                "states": heating_states,
            },
            "full_chain_rotation": {
                "status": "pass" if all(state["status"] == "pass"
                                        for state in rotation_states) else "fail",
                "states": rotation_states,
            },
            "cm_scattering_angle_distribution": {
                "status": "pass" if all(state["status"] == "pass"
                                        for state in angle_states) else "fail",
                "states": angle_states,
            },
            "near_escape_tail": {
                "status": "reported_diagnostic_only",
                "states": tail_states,
            },
            "interpretation": "Selected continuous-time weak-form equilibrium "
                              "checks, heating/cooling, CM-angle and rotation "
                              "checks; "
                              "not a pointwise stationary-distribution proof "
                              "or external numerical parity",
        },
        "reference_parity": {
            "status": "not_evaluated",
            "reason_code": parity_policy["reason_code"],
            "condition": parity_policy["activation_condition"],
        },
    }
    if smoke:
        for key in ("thermal_equilibrium", "heating_cooling",
                    "full_chain_rotation", "cm_scattering_angle_distribution",
                    "near_escape_tail"):
            _mark_smoke_diagnostics(report["physics_validation"][key])
        report["physics_validation"]["g0_eligible"] = False
    elif not qualifying_inputs:
        report["physics_validation"]["diagnostic_full_bound_result"] = (
            "pass" if scientific_groups_pass else "fail")
    return report


def _arguments(argv: list[str]) -> argparse.Namespace:
    """Parse the T03 CLI; full runs use contract-defined samples and seeds."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--sampler", type=Path, default=DEFAULT_SAMPLER)
    parser.add_argument("--output", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--smoke", action="store_true",
                        help="small wiring run; cannot pass T03 or G0")
    parser.add_argument("--smoke-samples-per-seed", type=int, default=1000)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Execute T03 and return nonzero if required scientific checks fail."""
    args = _arguments(sys.argv[1:] if argv is None else argv)
    try:
        contract_path = _inside_project(args.contract)
        sampler_path = _inside_project(args.sampler)
        if not sampler_path.is_file():
            raise ValueError(f"Build the project-local sampler first: {sampler_path}")
        contract = json.loads(contract_path.read_text(encoding="utf-8"),
                              parse_constant=_reject_nonfinite_json_constant)
        qualification = _contract_qualification(contract_path)
        sampler_state = _sampler_qualification(sampler_path)
        report = run_validation(contract, sampler_path, qualification,
                                sampler_state,
                                args.smoke, args.smoke_samples_per_seed)
        _write_report(args.output, report)
        print(json.dumps({
            "physics_validation": report["physics_validation"]["status"],
            "reference_parity": report["reference_parity"]["status"],
            "report": str(args.output),
        }))
        return 0 if args.smoke or report["physics_validation"]["status"] == "pass" else 1
    except (KeyError, IndexError, TypeError, ValueError, RuntimeError,
            OSError, subprocess.TimeoutExpired) as error:
        print(f"T03 validation failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
