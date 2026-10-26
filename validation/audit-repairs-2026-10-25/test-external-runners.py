"""Check runner failure/freshness contracts with fake executables, never external engines."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[2]


class RunnerContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="latents-runner-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.env = os.environ.copy()
        self.env["JOBS"] = "1"

    def fake(self, name, source):
        path = self.root / name
        path.write_text("#!/bin/sh\nset -e\n" + source)
        path.chmod(0o755)
        return str(path)

    def run_script(self, path, *args):
        return subprocess.run(["sh", str(path), *args], env=self.env,
                              cwd=self.root, capture_output=True, text=True, timeout=20)

    def test_feature_requires_normal_mplus_termination(self):
        runner = self.root / "features" / "run.sh"
        runner.parent.mkdir()
        shutil.copy2(REPO / "equivalence/features/run.sh", runner)
        kit = runner.parent / "fake" / "kit"
        kit.mkdir(parents=True)
        (kit / "case.inp").write_text("fake input\n")
        (kit / "case.dat").write_text("1 2\n")
        self.env["MPLUS"] = self.fake("mplus", 'printf "%s\\n" "$FAKE_OUTPUT" > "$2"\n')
        self.env["FAKE_OUTPUT"] = "ERROR: invalid model"
        failed = self.run_script(runner, "fake", "mplus")
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("no normally terminated output", failed.stderr)
        self.env["FAKE_OUTPUT"] = "THE MODEL ESTIMATION TERMINATED NORMALLY"
        self.assertEqual(self.run_script(runner, "fake", "mplus").returncode, 0)

    def test_families_cannot_recopy_stale_outputs(self):
        runner = REPO / "equivalence/latentgold-families/run-kit.sh"
        work = self.root / "work"
        work.mkdir()
        (work / "g01.lgs").write_text("fake syntax\n")
        (work / "g01.lst").write_text("Log-likelihood (LL) old\n")
        (work / "g01_posteriors.txt").write_text("old posteriors\n")
        returned = self.root / "returned"
        self.env.update(LATENTGOLD_WORK=str(work), LATENTGOLD_RETURNED=str(returned))
        self.env["WINE"] = self.fake("wine-no-output", "exit 0\n")
        failed = self.run_script(runner)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("missing fresh fitted listing", failed.stderr)
        self.assertFalse((returned / "g01.lst").exists())
        self.env["WINE"] = self.fake("wine-fresh", 'printf "Log-likelihood (LL) -10\\n" > "$5"\n'
                                     'printf "new posteriors\\n" > "${2%.lgs}_posteriors.txt"\n')
        self.assertEqual(self.run_script(runner).returncode, 0)
        self.assertEqual((returned / "g01_posteriors.txt").read_text(), "new posteriors\n")

    def check_montecarlo(self, relative):
        runner = self.root / "run.sh"
        shutil.copy2(REPO / relative, runner)
        kit = self.root / "kit"
        kit.mkdir()
        (kit / "case.inp").write_text("fake input\n")
        (kit / "case.dat").write_text("1 2\n")
        self.env["MPLUS"] = self.fake("mplus", 'printf "THE MODEL ESTIMATION TERMINATED NORMALLY\\n" > "$2"\n')
        result = self.run_script(runner, "mplus")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / "returned" / "case.out").is_file())

    def test_montecarlo_functions_survive_sh_entry(self):
        self.check_montecarlo("equivalence/montecarlo/run.sh")

    def test_extended_functions_survive_sh_entry(self):
        self.check_montecarlo("equivalence/montecarlo/extended/run.sh")


if __name__ == "__main__":
    unittest.main()
