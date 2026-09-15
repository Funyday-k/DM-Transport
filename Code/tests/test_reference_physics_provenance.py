"""Validate the project-local records for the T02 physics ports."""

from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
import tempfile
from typing import Optional
import unittest


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "Code/provenance/reference_physics.json"
UPSTREAMS = {
    "damascus_sun_evap": (
        "https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP.git",
        "b5678f5b193aa567ca10715c2a6c764c9e72eec7",
        "Code/provenance/licenses/DaMaSCUS-SUN-EVAP-MIT.txt",
    ),
    "obscura": (
        "https://github.com/temken/obscura.git",
        "4b1b9d87f8a69da4d8081175d322b66e48314328",
        "Code/provenance/licenses/obscura-MIT.txt",
    ),
    "libphysica": (
        "https://github.com/temken/libphysica.git",
        "fefbe47993bbd343277718708407fa7a40629738",
        "Code/provenance/licenses/libphysica-MIT.txt",
    ),
}
SOURCES = {
    "solar_model": (
        "damascus_sun_evap", "src/Solar_Model.cpp", {
            "solar_model_import_and_tables": (166, 233),
            "solar_target_construction": (251, 289),
            "solar_background_queries": (291, 370),
            "thermal_average_relative_velocity": (616, 649),
        }),
    "agss09_solar_model_data": (
        "damascus_sun_evap", "data/model_agss09.dat",
        {"agss09_model_table": (1, 1988)}),
    "patched_interpolation_evaluation": (
        "damascus_sun_evap", "cmake/PatchLibphysicaInterpolation.cmake",
        {"horner_interpolation_evaluation": (16, 24)}),
    "nuclear_target_model": (
        "obscura", "src/Target_Nucleus.cpp", {
            "isotope_name_and_mass_model": (25, 40),
            "nuclear_data_import": (130, 179),
        }),
    "nuclear_target_data": (
        "obscura", "data/Nuclear_Data.txt",
        {"nuclear_data_records": (1, 295)}),
    "steffen_interpolation": (
        "libphysica", "src/Numerics.cpp", {
            "steffen_coefficients": (21, 70),
            "steffen_construction_and_evaluation": (182, 238),
        }),
    "gauss_legendre_integration": (
        "libphysica", "src/Integration.cpp",
        {"gauss_legendre_30": (95, 107)}),
    "legacy_natural_units": (
        "libphysica", "src/Natural_Units.cpp", {
            "mass_and_length_units": (38, 63),
            "time_and_temperature_units": (78, 106),
            "particle_and_gravity_constants": (135, 164),
            "solar_mass_and_radius": (170, 178),
        }),
}
AGSS09_ATTRIBUTION = {
    "model": "Standard Solar Model AGSS09",
    "authors": "Serenelli et al.",
    "year": 2009,
    "reference": "arXiv:0909.2668",
    "basis": "Upstream table header",
    "licensing_note":
        "The table file states no separate data license; its inclusion is "
        "tracked under the enclosing DaMaSCUS-SUN-EVAP MIT source record.",
}
ALL_FRAGMENTS = {
    fragment
    for _, _, fragments in SOURCES.values()
    for fragment in fragments
}
PORT_FRAGMENTS = {
    "mean_relative_speed": {"thermal_average_relative_velocity"},
    "solar_background": ALL_FRAGMENTS - {"thermal_average_relative_velocity"},
}
DESTINATIONS = {
    "mean_relative_speed": {
        "Code/include/transport/physics/ScatteringPhysics.hpp",
        "Code/src/physics/ScatteringPhysics.cpp",
    },
    "solar_background": {
        "Code/include/transport/physics/SolarBackground.hpp",
        "Code/src/physics/SolarBackground.cpp",
        "Code/data/solar/model_agss09.dat",
        "Code/data/solar/Nuclear_Data.txt",
    },
}
DIFFERENCE_CATEGORIES = {
    "mean_relative_speed": {"units", "scope", "validation"},
    "solar_background": {
        "units", "data_loading", "data_format", "validation_and_domain",
        "input_contract",
    },
}
PRESERVED_BEHAVIORS = {
    "mean_relative_speed": set(),
    "solar_background": {
        "horner_interpolation",
        "escape_profile_numerics",
        "number_density_order",
        "isotopic_abundances",
    },
}
EXCLUSIONS = {
    "mpi_orchestration", "trajectory_propagation", "snapshot_io",
    "parameter_scan", "reference_build_artifacts",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _records(container: dict, key: str) -> list[dict]:
    value = container.get(key)
    _require(isinstance(value, list)
             and all(isinstance(item, dict) for item in value),
             f"{key} must be a list of objects")
    return value


def _index(records: list[dict], label: str) -> dict[str, dict]:
    ids = [record.get("id") for record in records]
    _require(all(isinstance(item, str) and item for item in ids),
             f"{label} ids must be nonempty strings")
    _require(len(ids) == len(set(ids)), f"{label} ids must be unique")
    return {record["id"]: record for record in records}


def _safe_path(value: object) -> bool:
    if not isinstance(value, str) or not value or "\\" in value:
        return False
    path = PurePosixPath(value)
    return (not path.is_absolute() and path.as_posix() == value
            and all(part not in {"", ".", ".."} for part in path.parts))


def _require_project_file(root: Path, value: object, label: str) -> None:
    _require(_safe_path(value), f"{label} must be a safe project-relative path")
    resolved_root = root.resolve()
    try:
        candidate = (resolved_root / PurePosixPath(value)).resolve(strict=True)
        candidate.relative_to(resolved_root)
    except (OSError, ValueError) as error:
        raise ValueError(f"{label} must exist inside the project") from error
    _require(candidate.is_file(), f"{label} must be a project file")


def _validate_differences(port: dict, expected: set[str]) -> None:
    differences = _records(port, "intentional_differences")
    categories = []
    for difference in differences:
        _require(
            set(difference) == {"category", "description", "behavior_change"}
            and isinstance(difference["category"], str)
            and bool(difference["category"].strip())
            and isinstance(difference["description"], str)
            and bool(difference["description"].strip())
            and type(difference["behavior_change"]) is bool,
            "intentional difference record is invalid")
        categories.append(difference["category"])
    _require(len(categories) == len(set(categories))
             and set(categories) == expected,
             f"{port['id']} intentional difference categories changed")


def validate_manifest(manifest: object, project_root: Path) -> None:
    """Check fixed declarations and local files without reading an upstream."""
    _require(isinstance(manifest, dict), "manifest must be an object")
    _require(set(manifest) == {
        "schema_version", "component", "upstreams", "source_files",
        "port_units", "scope_exclusions",
    }, "manifest fields changed")
    _require(manifest.get("schema_version") == 1, "schema_version must be 1")
    _require(manifest.get("component") == "transport_reference_physics",
             "unexpected component")
    exclusions = manifest.get("scope_exclusions")
    _require(isinstance(exclusions, list)
             and len(exclusions) == len(set(exclusions))
             and set(exclusions) == EXCLUSIONS,
             "scope exclusions changed")

    upstreams = _index(_records(manifest, "upstreams"), "upstream")
    _require(set(upstreams) == set(UPSTREAMS), "fixed upstream ids changed")
    for upstream_id, (repository, commit, notice) in UPSTREAMS.items():
        expected_license = {
            "status": "verified",
            "spdx_expression": "MIT",
            "upstream_path": "LICENSE",
            "local_notice_path": notice,
        }
        upstream = upstreams[upstream_id]
        _require(upstream == {
            "id": upstream_id,
            "repository": repository,
            "commit": commit,
            "license": expected_license,
        }, f"{upstream_id} fixed upstream record changed")
        _require(_safe_path(expected_license["upstream_path"]),
                 f"{upstream_id} upstream license path is unsafe")
        _require(PurePosixPath(notice).parts[:3] ==
                 ("Code", "provenance", "licenses"),
                 f"{upstream_id} local license path is outside provenance")
        _require_project_file(
            project_root, notice, f"{upstream_id} local license notice")

    sources = _index(_records(manifest, "source_files"), "source")
    _require(set(sources) == set(SOURCES), "fixed source ids changed")
    fragment_ids: set[str] = set()
    for source_id, (upstream_id, path, expected_fragments) in SOURCES.items():
        source = sources[source_id]
        fields = {"id", "upstream_id", "path", "fragments"}
        if source_id == "agss09_solar_model_data":
            fields.add("attribution")
        _require(set(source) == fields, f"{source_id} source fields changed")
        _require(source.get("upstream_id") == upstream_id,
                 f"{source_id} fixed upstream mapping changed")
        _require(source.get("path") == path and _safe_path(path),
                 f"{source_id} fixed source location changed")
        if source_id == "agss09_solar_model_data":
            _require(source.get("attribution") == AGSS09_ATTRIBUTION,
                     "AGSS09 attribution changed")

        fragments = _index(_records(source, "fragments"), "fragment")
        expected_records = {}
        for fragment_id, bounds in expected_fragments.items():
            if fragment_id == "nuclear_data_records":
                locator = ("start_record", "end_record")
            else:
                locator = ("start_line", "end_line")
            expected_records[fragment_id] = {
                "id": fragment_id,
                locator[0]: bounds[0],
                locator[1]: bounds[1],
            }
        _require(fragments == expected_records,
                 f"{source_id} source fragments changed")
        _require(not fragment_ids.intersection(fragments),
                 "fragment ids must be globally unique")
        fragment_ids.update(fragments)

    ports = _index(_records(manifest, "port_units"), "port")
    _require(set(ports) == set(PORT_FRAGMENTS), "fixed port ids changed")
    used_fragments: set[str] = set()
    for port_id, expected_fragments in PORT_FRAGMENTS.items():
        port = ports[port_id]
        fields = {
            "id", "derived_from", "destinations", "intentional_differences",
        }
        if PRESERVED_BEHAVIORS[port_id]:
            fields.add("preserved_behaviors")
        _require(set(port) == fields, f"{port_id} port fields changed")

        derived_from = port.get("derived_from")
        _require(isinstance(derived_from, list)
                 and len(derived_from) == len(set(derived_from))
                 and set(derived_from) == expected_fragments,
                 f"{port_id} fixed source mapping changed")
        used_fragments.update(derived_from)

        destinations = _records(port, "destinations")
        paths = [record.get("path") for record in destinations]
        _require(all(set(record) == {"path"} for record in destinations),
                 "destination records must contain only a path")
        for path in paths:
            _require(_safe_path(path)
                     and PurePosixPath(path).parts[0] == "Code",
                     "destination must be a safe path below Code")
            _require_project_file(project_root, path, "destination")
        _require(len(paths) == len(set(paths))
                 and set(paths) == DESTINATIONS[port_id],
                 f"{port_id} fixed destinations changed")

        _validate_differences(port, DIFFERENCE_CATEGORIES[port_id])
        expected_preserved = PRESERVED_BEHAVIORS[port_id]
        if expected_preserved:
            preserved = _index(
                _records(port, "preserved_behaviors"), "preserved behavior")
            _require(set(preserved) == expected_preserved,
                     f"{port_id} preserved behaviors changed")
            _require(all(
                set(record) == {"id", "description"}
                and isinstance(record["description"], str)
                and bool(record["description"].strip())
                for record in preserved.values()),
                "preserved behavior record is invalid")

    _require(used_fragments == fragment_ids,
             "every source fragment must map to a port")


def _load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _write_fixture(root: Path, omitted: Optional[str] = None) -> None:
    paths = {upstream[2] for upstream in UPSTREAMS.values()}
    for destinations in DESTINATIONS.values():
        paths.update(destinations)
    for path in paths - {omitted}:
        output = root / path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("fixture\n", encoding="utf-8")


class ProvenanceTests(unittest.TestCase):
    def test_repository_manifest(self) -> None:
        validate_manifest(_load_manifest(), ROOT)

    def test_rejects_destination_escape(self) -> None:
        manifest = _load_manifest()
        port = next(item for item in manifest["port_units"]
                    if item["id"] == "solar_background")
        port["destinations"][0]["path"] = "../outside"
        with self.assertRaisesRegex(ValueError, "safe path below Code"):
            validate_manifest(manifest, ROOT)

    def test_rejects_missing_local_artifacts(self) -> None:
        cases = (
            (UPSTREAMS["obscura"][2], "local license notice"),
            (next(iter(DESTINATIONS["solar_background"])), "destination"),
        )
        for omitted, message in cases:
            with self.subTest(path=omitted):
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    _write_fixture(root, omitted)
                    with self.assertRaisesRegex(ValueError, message):
                        validate_manifest(_load_manifest(), root)

    def test_accepts_equivalent_local_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _write_fixture(root)
            validate_manifest(_load_manifest(), root)


if __name__ == "__main__":
    unittest.main()
