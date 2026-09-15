"""Validate the local records for the T02 physics port."""

from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "Code/provenance/reference_physics.json"
REPOSITORY = "https://github.com/Funyday-k/DaMaSCUS-SUN-EVAP.git"
COMMIT = "b5678f5b193aa567ca10715c2a6c764c9e72eec7"
LICENSE = {
    "status": "verified",
    "spdx_expression": "MIT",
    "upstream_path": "LICENSE",
    "local_notice_path":
        "Code/provenance/licenses/DaMaSCUS-SUN-EVAP-MIT.txt",
}
DESTINATIONS = {
    "Code/include/transport/physics/ScatteringPhysics.hpp",
    "Code/src/physics/ScatteringPhysics.cpp",
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
    _require(isinstance(value, list) and
             all(isinstance(item, dict) for item in value),
             f"{key} must be a list of objects")
    return value


def _record(records: list[dict], record_id: str, label: str) -> dict:
    matches = [record for record in records if record.get("id") == record_id]
    _require(len(matches) == 1, f"expected one {label} record")
    return matches[0]


def _safe_relative_path(value: object) -> bool:
    if not isinstance(value, str) or not value or "\\" in value:
        return False
    path = PurePosixPath(value)
    return (not path.is_absolute() and path.as_posix() == value
            and all(part not in {"", ".", ".."} for part in path.parts))


def _require_project_file(project_root: Path, value: object,
                          label: str) -> None:
    _require(_safe_relative_path(value),
             f"{label} must be a safe project-relative path")
    root = project_root.resolve()
    try:
        candidate = (root / PurePosixPath(value)).resolve(strict=True)
        candidate.relative_to(root)
    except (OSError, ValueError) as error:
        raise ValueError(f"{label} must exist inside the project") from error
    _require(candidate.is_file(), f"{label} must be a project file")


def validate_manifest(manifest: object, project_root: Path) -> None:
    """Check fixed declarations and local files without reading an upstream."""
    _require(isinstance(manifest, dict), "manifest must be an object")
    _require(manifest.get("schema_version") == 1,
             "schema_version must be 1")
    _require(manifest.get("component") == "transport_reference_physics",
             "unexpected component")
    exclusions = manifest.get("scope_exclusions")
    _require(isinstance(exclusions, list)
             and all(isinstance(item, str) for item in exclusions)
             and len(exclusions) == len(set(exclusions))
             and set(exclusions) == EXCLUSIONS,
             "scope exclusions changed")

    upstream = _record(
        _records(manifest, "upstreams"), "damascus_sun_evap", "upstream")
    _require(upstream.get("repository") == REPOSITORY,
             "fixed repository changed")
    _require(upstream.get("commit") == COMMIT, "fixed commit changed")
    _require(upstream.get("license") == LICENSE,
             "fixed license record changed")
    _require_project_file(
        project_root, LICENSE["local_notice_path"], "local license notice")

    source = _record(
        _records(manifest, "source_files"), "solar_model", "source")
    _require(source.get("upstream_id") == "damascus_sun_evap"
             and source.get("path") == "src/Solar_Model.cpp",
             "fixed source location changed")
    fragment = _record(
        _records(source, "fragments"),
        "thermal_average_relative_velocity", "fragment")
    _require(fragment == {
        "id": "thermal_average_relative_velocity",
        "start_line": 616,
        "end_line": 649,
    }, "fixed source fragment changed")

    port = _record(
        _records(manifest, "port_units"), "mean_relative_speed", "port")
    _require(port.get("derived_from") ==
             ["thermal_average_relative_velocity"],
             "fixed source mapping changed")
    destinations = _records(port, "destinations")
    _require(all(set(record) == {"path"} for record in destinations),
             "destination records must contain only a path")
    for destination in destinations:
        path = destination["path"]
        _require(isinstance(path, str) and path.startswith("Code/"),
                 "destination must be below Code")
        _require_project_file(project_root, path, "destination")
    _require(len(destinations) == len(DESTINATIONS)
             and {item["path"] for item in destinations} == DESTINATIONS,
             "fixed destinations changed")

    differences = _records(port, "intentional_differences")
    _require(bool(differences), "intentional differences are required")
    for difference in differences:
        _require(
            set(difference) == {"category", "description", "behavior_change"}
            and isinstance(difference["category"], str)
            and bool(difference["category"].strip())
            and isinstance(difference["description"], str)
            and bool(difference["description"].strip())
            and type(difference["behavior_change"]) is bool,
            "intentional difference record is invalid")


def _load_manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def _write_project_files(root: Path, include_notice: bool = True) -> None:
    paths = set(DESTINATIONS)
    if include_notice:
        paths.add(LICENSE["local_notice_path"])
    for path in paths:
        output = root / path
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("fixture\n", encoding="utf-8")


class ProvenanceTests(unittest.TestCase):
    def test_repository_manifest(self) -> None:
        validate_manifest(_load_manifest(), ROOT)

    def test_rejects_commit_drift(self) -> None:
        manifest = _load_manifest()
        manifest["upstreams"][0]["commit"] = "0" * 40
        with self.assertRaisesRegex(ValueError, "fixed commit"):
            validate_manifest(manifest, ROOT)

    def test_rejects_empty_source_records(self) -> None:
        manifest = _load_manifest()
        manifest["source_files"] = []
        with self.assertRaisesRegex(ValueError, "one source"):
            validate_manifest(manifest, ROOT)

    def test_rejects_source_location_drift(self) -> None:
        manifest = _load_manifest()
        manifest["source_files"][0]["path"] = "src/Other.cpp"
        with self.assertRaisesRegex(ValueError, "source location"):
            validate_manifest(manifest, ROOT)

    def test_rejects_destination_escape(self) -> None:
        manifest = _load_manifest()
        manifest["port_units"][0]["destinations"][0]["path"] = "../outside"
        with self.assertRaisesRegex(ValueError, "below Code"):
            validate_manifest(manifest, ROOT)

    def test_rejects_missing_license_notice(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _write_project_files(root, include_notice=False)
            with self.assertRaisesRegex(ValueError, "license notice"):
                validate_manifest(_load_manifest(), root)

    def test_accepts_equivalent_local_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _write_project_files(root)
            validate_manifest(_load_manifest(), root)


if __name__ == "__main__":
    unittest.main()
