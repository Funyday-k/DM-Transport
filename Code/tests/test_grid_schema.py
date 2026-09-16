"""Cross-language T04 metadata, indexing, and finite-volume measure checks."""

import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from grid_schema import GridSchema, load_grid_schema  # noqa: E402


class GridSchemaTest(unittest.TestCase):
    """Validate C++ witness values against independent Python calculations."""

    @classmethod
    def setUpClass(cls) -> None:
        """Load the standalone metadata and per-cell C++ witness records."""
        if len(sys.argv) != 2:
            raise RuntimeError("usage: test_grid_schema.py PROBE_EXECUTABLE")
        executable = sys.argv[1]
        schema_text = subprocess.run(
            [executable, "--schema"], check=True, capture_output=True, text=True
        ).stdout
        probe_text = subprocess.run(
            [executable, "--probe"], check=True, capture_output=True, text=True
        ).stdout
        cls.metadata = json.loads(schema_text)
        cls.probe = json.loads(probe_text)
        cls.schema = GridSchema.from_mapping(cls.metadata)
        cls.schema_from_fixed_stream = json.loads(subprocess.run(
            [executable, "--schema-fixed"], check=True,
            capture_output=True, text=True
        ).stdout)

    def test_standalone_schema_round_trip(self) -> None:
        """Verify exactly the same metadata travels in the diagnostic wrapper."""
        self.assertEqual(self.probe["schema"], self.metadata)
        self.assertEqual(self.schema_from_fixed_stream, self.metadata)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "grid.json"
            path.write_text(json.dumps(self.metadata), encoding="utf-8")
            self.assertEqual(load_grid_schema(path), self.schema)

    def test_every_cpp_cell_matches_python(self) -> None:
        """Compare all nonuniform cells including the 8*pi^2 cgs measure."""
        nr, nv, nmu = self.schema.shape
        cells = self.probe["cpp_cells"]
        self.assertEqual(len(cells), nr * nv * nmu)
        self.assertEqual(len({tuple(row["index"]) for row in cells}), len(cells))
        self.assertEqual(sorted(row["flat"] for row in cells), list(range(len(cells))))
        r = self.schema.r_faces_cm
        v = self.schema.v_faces_cm_s
        mu = self.schema.mu_faces
        self.assertGreater(len(set(r[i + 1] - r[i] for i in range(nr))), 1)
        self.assertGreater(len(set(v[i + 1] - v[i] for i in range(nv))), 1)
        self.assertGreater(len(set(mu[i + 1] - mu[i] for i in range(nmu))), 1)
        for row in cells:
            ir, iv, imu = row["index"]
            self.assertEqual(row["flat"], self.schema.flatten(ir, iv, imu))
            python_measure = self.schema.cell_phase_measure(ir, iv, imu)
            # Independent shell expression: (8*pi^2/9) Delta(r^3)
            # Delta(v^3) Delta(mu), in cm^6/s^3.
            shell_measure = (8.0 * math.pi * math.pi / 9.0 *
                             (r[ir + 1] ** 3 - r[ir] ** 3) *
                             (v[iv + 1] ** 3 - v[iv] ** 3) *
                             (mu[imu + 1] - mu[imu]))
            self.assertTrue(math.isclose(row["phase_measure"], python_measure,
                                         rel_tol=5e-15, abs_tol=0.0))
            self.assertTrue(math.isclose(row["phase_measure"], shell_measure,
                                         rel_tol=5e-15, abs_tol=0.0))

    def test_invalid_schema_rejected(self) -> None:
        """Reject metadata that silently changes geometry or physics meaning."""
        mutations = (
            ("schema_version", 2),
            ("schema_version", True),
            ("flatten_order", ["mu", "v", "r"]),
            ("units", {**self.metadata["units"], "v": "km/s"}),
            ("quadrature_order_per_axis", 3),
            ("r_faces_cm", [0.0, 2.0, 1.0]),
            ("v_faces_cm_s", [0.0, float("nan")]),
            ("mu_faces", [-1.0, 0.2, 0.99]),
            ("escape_threshold", {**self.metadata["escape_threshold"],
                                  "bound_condition": "E<=0"}),
        )
        for key, value in mutations:
            with self.subTest(key=key, value=value):
                invalid = copy.deepcopy(self.metadata)
                invalid[key] = value
                with self.assertRaises(ValueError):
                    GridSchema.from_mapping(invalid)
        invalid = copy.deepcopy(self.metadata)
        del invalid["units"]
        with self.assertRaises(ValueError):
            GridSchema.from_mapping(invalid)
        with self.assertRaises(IndexError):
            self.schema.flatten(-1, 0, 0)

    def test_duplicate_and_nonfinite_json_rejected(self) -> None:
        """Reject ambiguous keys and non-standard JSON constants on file input."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "grid.json"
            path.write_text('{"schema_version":1,"schema_version":2}',
                            encoding="utf-8")
            with self.assertRaises(ValueError):
                load_grid_schema(path)
            nonfinite = copy.deepcopy(self.metadata)
            nonfinite["r_faces_cm"][0] = float("nan")
            path.write_text(json.dumps(nonfinite), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_grid_schema(path)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
