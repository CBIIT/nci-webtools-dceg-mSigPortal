"""Unit checks for safe Docker command construction; no Docker daemon required."""
import argparse
from pathlib import Path
import unittest

import run_r


class Commands(unittest.TestCase):
    def args(self, **overrides):
        values = dict(image="backend:fixture", baseline=None, service="backend",
                      server="/deploy/server", timeout=20, case=None, list=False,
                      suite="backend", container=None, staging="/tmp/r-test-staging")
        values.update(overrides)
        return argparse.Namespace(**values)

    def test_image_runs_offline_and_mounts_only_tests_and_reports(self):
        cmd = run_r.command(self.args(), Path("/tmp/path with spaces"))
        self.assertEqual(cmd[:6], ["docker", "run", "--rm", "--network", "none", "--entrypoint"])
        self.assertIn("/tmp/path with spaces:/regression", cmd)
        self.assertIn(f"{run_r.HERE / 'r'}:/opt/msigportal-r-regression:ro", cmd)
        self.assertFalse(any(str(run_r.REPO / "server") in arg for arg in cmd))

    def test_compose_does_not_start_dependencies(self):
        cmd = run_r.command(self.args(image=None), Path("/tmp/reports"))
        self.assertEqual(cmd[:6], ["docker", "compose", "run", "--rm", "--no-deps", "--entrypoint"])
        self.assertIn("backend", cmd)

    def test_baseline_is_read_only_and_selector_is_one_argument(self):
        selector = "msLandscape|association"
        cmd = run_r.command(self.args(baseline=Path("/tmp/baseline"), case=selector), Path("/tmp/reports"))
        self.assertIn(f"{Path('/tmp/baseline').resolve()}:/baseline:ro", cmd)
        self.assertEqual(cmd[cmd.index("--case") + 1], selector)

    def test_existing_container_uses_exec_and_temporary_paths(self):
        cmd = run_r.command(self.args(image=None, container="backend-live"), Path("/tmp/reports"))
        self.assertEqual(cmd[:4], ["docker", "exec", "backend-live", "Rscript"])
        self.assertIn("/tmp/r-test-staging/suite/run.R", cmd)
        self.assertIn("/tmp/r-test-staging/reports", cmd)
        self.assertNotIn("--volume", cmd)

    def test_refitting_uses_separate_entrypoint(self):
        cmd = run_r.command(self.args(suite="refitting", server="/refitting-service"), Path("/tmp/reports"))
        self.assertIn("/opt/msigportal-r-regression/refitting_run.R", cmd)
        self.assertNotIn("/opt/msigportal-r-regression/run.R", cmd)


if __name__ == "__main__":
    unittest.main()
