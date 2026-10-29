"""Historical provenance must never turn changed numeric evidence into a run."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).parents[1] / "attest-historical.py"
SPEC = importlib.util.spec_from_file_location("attestation", SCRIPT)
ATTESTATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ATTESTATION)


class HistoricalEvidenceTests(unittest.TestCase):
    original = b'case,fn,latents,agrees\nexample,enumeration,1.25,TRUE\n'

    def validate(self, contents):
        with tempfile.TemporaryDirectory() as directory:
            result = Path(directory) / "comparison.csv"
            result.write_bytes(contents)
            with patch.object(ATTESTATION, "git_bytes", return_value=self.original):
                return ATTESTATION.validate_result(str(result))

    def test_identical_evidence_and_function_reattribution_are_explicit(self):
        self.assertIn("bytes match", self.validate(self.original))
        metadata = self.original.replace(b"enumeration", b"lca")
        note = self.validate(metadata)
        self.assertIn("Function-label metadata", note)
        self.assertIn("not a numerical rerun", note)

    def test_changed_estimate_cannot_receive_a_historical_attestation(self):
        with self.assertRaisesRegex(ValueError, "cannot attest"):
            self.validate(self.original.replace(b"1.25", b"9999"))

    def test_changed_outcome_cannot_receive_a_historical_attestation(self):
        with self.assertRaisesRegex(ValueError, "cannot attest"):
            self.validate(self.original.replace(b"TRUE", b"FALSE"))

    def test_dropped_comparison_cannot_receive_a_historical_attestation(self):
        with self.assertRaisesRegex(ValueError, "cannot attest"):
            self.validate(b"case,fn,latents,agrees\n")


if __name__ == "__main__":
    unittest.main()
