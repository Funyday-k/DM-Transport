"""Regression checks for safe, correctly labelled legacy baseline execution."""

import copy
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "Code/python"))
from run_baseline import check_case_contract, parse_case_log, render_case_config, run_process


def setting(text: str, name: str) -> str:
    """Read one unambiguous generated libconfig assignment for assertions."""
    matches = re.findall(r"^" + re.escape(name) + r"\s*=\s*(.+);$", text, re.MULTILINE)
    if len(matches) != 1:
        raise AssertionError(f"Expected exactly one setting named {name}")
    return matches[0]


class BaselineTests(unittest.TestCase):
    """Protect coupling semantics, source inputs and wall-time interpretation."""

    def setUp(self) -> None:
        """Load the one maintained MVP configuration (GeV, cm², wall seconds)."""
        self.config = json.loads(
            (ROOT / "Code/configs/benchmark/mvp.json").read_text(encoding="utf-8")
        )

    def test_effective_proton_only_configuration(self) -> None:
        """Do not let the old isospin switch silently replace fn=0 by fn=fp."""
        for mode in ("Capture", "Parameter point"):
            for sigma in self.config["benchmark"]["cross_sections_cm2"]:
                text = render_case_config(self.config, mode, sigma, Path("/tmp/baseline data"))
                self.assertEqual(setting(text, "DM_isospin_conserved"), "false")
                self.assertEqual(setting(text, "DM_relative_couplings"), "(1.0, 0.0)")
                self.assertEqual(float(setting(text, "DM_cross_section_nucleon")), sigma)
                self.assertEqual(float(setting(text, "DM_cross_section_electron")), 0)
                self.assertEqual(setting(text, "DM_light"), "true")
                self.assertEqual(setting(text, "interpolation_points"), "0")
                self.assertEqual(
                    int(setting(text, "max_trajectories")),
                    self.config["benchmark"]["attempts_per_case"],
                )

    def test_halo_and_quoted_paths_follow_inputs(self) -> None:
        """Read halo values from the canonical input and escape output filenames."""
        config = copy.deepcopy(self.config)
        config["source"]["halo"]["local_density_GeV_cm3"] = 0.35
        output = Path('/tmp/baseline "quoted" path')
        text = render_case_config(config, "Capture", 1e-34, output)
        self.assertEqual(float(setting(text, "DM_local_density")), 0.35)
        self.assertEqual(json.loads(setting(text, "output_dir")), str(output.resolve()) + "/")
        with self.assertRaises(ValueError):
            render_case_config(config, "Capture", 1e-34, Path("/tmp/invalid\npath"))

    def test_counts_are_not_absolute_capture_rates(self) -> None:
        """Missing summaries stay unknown and throughput is not the physical C."""
        metrics = parse_case_log(
            "Simulated trajectories:\t16\n"
            "Capture-classified trajectories:\t12\n"
            "Captured count:\t3\n"
            "Numerical failure count:\t1\n"
            "Captured particle rate [1/s]:\t9999\n"
            "Simulation time [s]:\t0.125\n"
        )
        self.assertEqual(metrics["capture_probability_raw"], 3 / 16)
        self.assertEqual(metrics["capture_probability_classified"], 3 / 12)
        self.assertEqual(metrics["legacy_simulation_wall_time_s"], 0.125)
        self.assertIsNone(metrics["complete_evaporations"])
        self.assertNotIn("absolute_capture_rate_particles_per_s", metrics)
        self.assertIsNone(parse_case_log("")["attempted"])

    def test_subprocess_timeout_is_reported(self) -> None:
        """A runaway baseline is terminated and kept as a timeout, not a pass."""
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            result = run_process(
                [sys.executable, "-c", "import time; time.sleep(30)"],
                directory, directory / "timeout.log", 0.1,
            )
        self.assertEqual(result["status"], "timeout")
        self.assertNotEqual(result["returncode"], 0)

    def test_executed_couplings_must_match_the_requested_model(self) -> None:
        """A valid cfg cannot excuse a different model in the actual binary log."""
        log = """git:b5678f5
Mass: 100 MeV
Spin: 0.5
Interaction: Spin-Dependent (SD)
Low mass: [x]
Isospin conservation: [ ]
Sc. rate interpolation: [ ]
Coupling ratio: fn/fp = 0
Sigma_P[cm^2]: 1e-34
Sigma_N[cm^2]: 0
Sigma_E[cm^2]: 0
Trajectory boundary [Rsun]: 1.1
"""
        head = "b5678f5b193aa567ca10715c2a6c764c9e72eec7"
        contract = check_case_contract(parse_case_log(log), self.config, 1e-34, head, {})
        self.assertTrue(contract["runtime_model_verified"])
        self.assertIsNone(contract["artifact_mass_sigma_verified"])
        wrong_coupling = log.replace("fn/fp = 0", "fn/fp = 1")
        contract = check_case_contract(
            parse_case_log(wrong_coupling), self.config, 1e-34, head, {}
        )
        self.assertFalse(contract["runtime_model_verified"])
        contract = check_case_contract(parse_case_log(""), self.config, 1e-34, head, {})
        self.assertFalse(contract["runtime_model_verified"])


if __name__ == "__main__":
    unittest.main()
