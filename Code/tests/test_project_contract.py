"""Check the document inventory and internal consistency of the bootstrap contract.

These checks protect repository conventions, not numerical physics accuracy.
"""

import ast
import hashlib
import json
import math
import os
from pathlib import Path
import re
import unittest
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[2]
ACTIVE_DOCUMENTS = {
    Path("README.md"),
    Path("Note/Proposal.md"),
    Path("Note/Task_Plan.md"),
    Path("CHANGELOG.md"),
}
ARCHIVE = Path("Note/archive/Proposal_2026-09-14_v1.md")


def source_tree_files() -> list[Path]:
    """Return repository source files while ignoring generated local trees."""
    files = []
    for directory, directories, names in os.walk(ROOT, followlinks=False):
        directory_path = Path(directory)
        at_root = directory_path == ROOT
        directories[:] = [
            name for name in directories
            if name not in {".git", ".venv", "Output", "__pycache__", "build"}
            and not (at_root and name.startswith("build"))
        ]
        files.extend(directory_path / name for name in names)
    return files


def dotted_name(node: ast.AST):
    """Return a static dotted callable name, if the AST node has one."""
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = dotted_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return None


class ProjectContractTests(unittest.TestCase):
    """Prevent document proliferation, broken task dependencies and silent gates."""

    def test_fixed_documents_and_frozen_original(self) -> None:
        """Keep four maintained Markdown files and the one frozen original."""
        documents = set()
        for path in ROOT.rglob("*.md"):
            relative = path.relative_to(ROOT)
            top = relative.parts[0]
            if top.startswith(".") or top.startswith("build") or top == "Output":
                continue
            documents.add(relative)
        self.assertEqual(documents, ACTIVE_DOCUMENTS | {ARCHIVE})
        self.assertEqual(
            hashlib.sha256((ROOT / ARCHIVE).read_bytes()).hexdigest(),
            "5b0c75bfe29b8a2b3aca6c5945dec5917e37d0ff65b5c34af4e5bb2cac72810b",
        )

    def test_reference_repository_is_not_a_build_or_execution_input(self) -> None:
        """Keep the EVAP repository outside this project's build and runtime tools."""
        source_files = source_tree_files()
        build_inputs = {
            path for path in source_files
            if path.name == "CMakeLists.txt"
            or path.suffix.lower() == ".cmake"
            or (
                path.suffix.lower() in {".yml", ".yaml"}
                and path.relative_to(ROOT).parts[:2] == (".github", "workflows")
            )
        }
        runtime_python = {
            path for path in source_files
            if path.suffix.lower() == ".py" and "tests" not in path.relative_to(ROOT).parts
        }
        self.assertIn(ROOT / "CMakeLists.txt", build_inputs)
        self.assertIn(ROOT / ".github/workflows/ci.yml", build_inputs)
        self.assertIn(ROOT / "Code/python/run_baseline.py", runtime_python)
        forbidden_build_tokens = {
            "DMTRANSPORT_ENABLE_DAMASCUS_PHYSICS",
            "DMTRANSPORT_DAMASCUS_SOURCE_DIR",
            "DaMaSCUS::Physics",
            "DaMaSCUS-SUN-EVAP",
            "damascus_physics_consumer",
            "--damascus-source",
            "--damascus-build",
        }
        for path in sorted(build_inputs | runtime_python):
            try:
                path.resolve(strict=True).relative_to(ROOT)
            except (OSError, ValueError) as error:
                self.fail(f"Build/runtime source resolves outside the project: {path}: {error}")
            text = path.read_text(encoding="utf-8")
            for token in forbidden_build_tokens:
                self.assertNotIn(token, text, f"{path.relative_to(ROOT)}: {token}")

        runner_path = ROOT / "Code/python/run_baseline.py"
        runner_source = runner_path.read_text(encoding="utf-8")
        tree = ast.parse(runner_source, filename=str(runner_path))
        imported_modules = {
            alias.name.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_modules.update(
            node.module.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module
        )
        safe_imports = {
            "__future__", "argparse", "hashlib", "json", "math", "pathlib",
            "re", "sys", "typing",
        }
        self.assertEqual(imported_modules - safe_imports, set())

        forbidden_calls = {
            "__import__", "compile", "eval", "exec", "open",
            "FileType", "dump",
            "chmod", "hardlink_to", "lchmod", "link_to", "mkdir",
            "rename", "replace", "rmdir", "symlink_to", "touch", "unlink",
            "write_bytes", "write_text",
            "Popen", "call", "check_call", "check_output", "create_subprocess_exec",
            "create_subprocess_shell", "execv", "execve", "execvp", "execvpe",
            "fork", "forkpty", "popen", "posix_spawn", "posix_spawnp", "run",
            "spawnl", "spawnle", "spawnlp", "spawnlpe", "spawnv", "spawnve",
            "spawnvp", "spawnvpe", "system",
            "create_connection", "socket", "urlopen", "urlretrieve",
        }
        violations = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = dotted_name(node.func)
            leaf = name.rsplit(".", 1)[-1] if name else None
            if leaf in forbidden_calls:
                violations.append((node.lineno, name))
        self.assertEqual(violations, [], f"unsafe calls in {runner_path.relative_to(ROOT)}")

    def test_document_links_and_unique_prose(self) -> None:
        """Check project-local links and exact duplicate long prose paragraphs.

        Local links are relative to their document and use explicit HTML anchors.
        External URLs are skipped; unpublished artifacts should use code text.
        Semantic redundancy still requires a review of each change.
        """
        paragraphs = {}
        for relative in sorted(ACTIVE_DOCUMENTS):
            text = (ROOT / relative).read_text(encoding="utf-8")
            self.assertEqual(text.count("```") % 2, 0, str(relative))
            for target in re.findall(
                r"\[[^\]\n]+\]\((<[^>\n]+>|[^)\s]+)\)", text
            ):
                target = target[1:-1] if target.startswith("<") else target
                decoded_target = unquote(target)
                self.assertFalse(
                    decoded_target.startswith(("/", "\\", "~/"))
                    or re.match(r"^[A-Za-z]:[\\/]", decoded_target),
                    f"Local links must be relative to {relative}",
                )
                link = urlsplit(target)
                self.assertNotEqual(
                    link.scheme.lower(), "file",
                    f"Local links must be relative to {relative}",
                )
                if link.scheme:
                    continue
                filename = unquote(link.path)
                document = ROOT / relative
                path = (document.parent / filename).resolve() if filename else document
                self.assertTrue(path.is_file(), f"{relative}: {target}")
                if link.fragment:
                    anchor = unquote(link.fragment)
                    self.assertRegex(
                        path.read_text(encoding="utf-8"),
                        rf"\bid\s*=\s*([\"']){re.escape(anchor)}\1",
                        f"{relative}: {target}",
                    )
            if relative == Path("CHANGELOG.md"):
                continue  # Historical summaries may intentionally repeat facts.
            for block in text.split("\n\n"):
                block = block.strip()
                if len(block) < 100 or "\n" in block or block.startswith(
                    ("|", "#", "```", "$$", "\\[")
                ):
                    continue
                self.assertNotIn(
                    block, paragraphs,
                    f"Duplicate prose in {relative} and {paragraphs.get(block)}",
                )
                paragraphs[block] = relative

    def test_task_dependencies(self) -> None:
        """Every named prerequisite exists and the task graph is acyclic."""
        text = (ROOT / "Note/Task_Plan.md").read_text(encoding="utf-8")
        rows = re.findall(
            r"^\| ([TP]\d{2}) \|.*?\| ([^|]+) \|[^|]+\|$", text, re.MULTILINE
        )
        graph = {task: re.findall(r"[TP]\d{2}", deps) for task, deps in rows}
        expected = {f"T{i:02}" for i in range(21)} | {f"P{i:02}" for i in range(10)}
        self.assertEqual(set(graph), expected)
        self.assertEqual(len(rows), len(graph))
        for task, dependencies in rows:
            self.assertTrue(
                re.fullmatch(r"[TP]\d{2}(?:,[TP]\d{2})*", dependencies.strip())
                or (task == "T00" and dependencies.strip() == "—")
                or (task == "T20" and dependencies.strip() == "对应前序门槛"),
                f"Malformed dependency column: {task}: {dependencies}",
            )
        visited = set()
        active = set()

        def visit(task: str) -> None:
            self.assertIn(task, graph)
            self.assertNotIn(task, active, f"Dependency cycle at {task}")
            if task in visited:
                return
            active.add(task)
            for dependency in graph[task]:
                visit(dependency)
            active.remove(task)
            visited.add(task)

        for task in graph:
            visit(task)

    def test_configured_and_pending_gates(self) -> None:
        """Null tolerances remain explicitly pending instead of passing by default."""
        config = json.loads(
            (ROOT / "Code/configs/benchmark/mvp.json").read_text(encoding="utf-8")
        )
        expected_tolerances = {
            "probability_row_sum_abs", "generator_scaled_row_sum_abs",
            "steady_relative_residual", "source_sink_relative_mismatch",
            "integrated_observable_relative_budget",
            "streaming_fake_escape_absolute_budget", "tail_observable_absolute_budget",
            "negative_occupation_tolerance",
        }
        self.assertEqual(set(config["tolerances"]), expected_tolerances)
        for value in config["tolerances"].values():
            if value is not None:
                self.assertIn(type(value), (float, int))
                self.assertTrue(math.isfinite(value) and value > 0)
        pending = {key for key, value in config["tolerances"].items() if value is None}
        self.assertEqual(pending, set(config["pending_gates"]))
        for entry in config["pending_gates"].values():
            self.assertTrue(entry["reason"])
            self.assertRegex(entry["owner_task"], r"^T\d{2}$")
            self.assertTrue(entry["must_set_before"].startswith(entry["owner_task"]))
            self.assertTrue(entry["must_set_before"])
        cells = config["state"]["grid_cells"]
        self.assertEqual(len(cells), 3)
        self.assertTrue(all(type(n) is int and n > 0 for n in cells))
        faces = config["state"]["mu_faces"]
        self.assertEqual(len(faces), cells[2] + 1)
        self.assertEqual((faces[0], faces[-1]), (-1.0, 1.0))
        self.assertTrue(all(a < b for a, b in zip(faces, faces[1:])))
        self.assertIn(
            config["benchmark"]["regression_cross_section_cm2"],
            config["benchmark"]["cross_sections_cm2"],
        )


if __name__ == "__main__":
    unittest.main()
