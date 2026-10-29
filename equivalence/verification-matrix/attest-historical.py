"""Validate the 1.0.0 historical evidence; --write adds retrospective sidecars.

This does not run comparisons or pretend a source hash was captured during an
old run. It validates the retained release-fixes SHA256 manifest, both release
commits, and the original 1.0.0 result blobs. Only documented function-label
reattribution may differ in a result CSV. Run from the repository root.
"""
import argparse
import csv
import datetime
import hashlib
import io
import json
from pathlib import Path
import subprocess


def git_bytes(*arguments):
    return subprocess.check_output(["git", *arguments])


def runtime_files():
    return ["NAMESPACE"] + sorted(str(path) for path in Path("R").rglob("*") if path.is_file())


def validate_sources(paths):
    manifest_path = Path("validation/release-fixes-2026-10-28/release/source-manifest-final.json")
    retained = json.loads(manifest_path.read_text())["sha256"]
    expected = {name: digest for name, digest in retained.items()
                if name == "NAMESPACE" or name.startswith("R/")}
    if set(expected) != set(paths):
        raise ValueError("Historical source manifest does not cover the complete current runtime source.")
    for path in paths:
        data = Path(path).read_bytes()
        if hashlib.sha256(data).hexdigest() != expected[path]:
            raise ValueError(f"Runtime source differs from the retained release manifest: {path}")
        for revision in ("344dc7c", "7d0e0df"):
            if git_bytes("show", f"{revision}:{path}") != data:
                raise ValueError(f"Runtime source differs from {revision}: {path}")
    for revision in ("344dc7c", "7d0e0df"):
        saved = set(git_bytes("ls-tree", "-r", "--name-only", revision, "R", "NAMESPACE")
                    .decode().splitlines())
        if saved != set(paths):
            raise ValueError(f"Runtime file inventory differs from {revision}.")


def validate_result(path):
    original = git_bytes("show", f"7d0e0df:{path}")
    current = Path(path).read_bytes()
    if current == original:
        return "Result bytes match the result committed in 7d0e0df."
    old_rows = list(csv.DictReader(io.StringIO(original.decode())))
    new_rows = list(csv.DictReader(io.StringIO(current.decode())))
    def without_function(rows):
        return [{key: value for key, value in row.items() if key != "fn"} for row in rows]
    if not old_rows or without_function(old_rows) != without_function(new_rows):
        raise ValueError(f"Result changed beyond function-label metadata; cannot attest: {path}")
    return ("Function-label metadata was reattributed after review; all other result fields, "
            "including compared values and outcomes, match 7d0e0df. This is not a numerical rerun.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write validated retrospective sidecars.")
    options = parser.parse_args()
    paths = runtime_files()
    validate_sources(paths)
    evidence = [
        "validation/release-fixes-2026-10-28/equivalence/report.csv",
        "validation/release-fixes-2026-10-28/release/external-tests-consolidated.csv",
        "validation/release-fixes-2026-10-28/release/retained-references-results.csv",
        "validation/matrix-refresh-2026-10-28/equivalence-results.csv",
    ]
    with open("equivalence/verification-matrix/kit-evidence.csv", newline="") as handle:
        evidence.extend(row["path"] for row in csv.DictReader(handle))
    # Validate every result before writing anything.
    notes = {path: validate_result(path) for path in evidence}
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    source_rows = [{"kind": "source", "path": path,
                    "md5": hashlib.md5(Path(path).read_bytes()).hexdigest()} for path in paths]
    fields = ["kind", "path", "md5", "basis", "recorded_at", "captured_at", "version", "note"]
    for path, result_note in notes.items():
        note = ("Retrospective attestation, not a rerun or a contemporaneous capture. "
                "All runtime R files and NAMESPACE match the retained release-fixes "
                "source-manifest-final.json SHA256 inventory and commits 344dc7c (0.98.1) "
                "and 7d0e0df (1.0.0). Existing run records supply the original run dates. "
                + result_note)
        rows = source_rows + [{"kind": "result", "path": Path(path).name,
                               "md5": hashlib.md5(Path(path).read_bytes()).hexdigest()}]
        if options.write:
            with open(path + ".provenance.csv", "w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(dict(row, basis="retrospective source-equivalence attestation",
                                      recorded_at=stamp, captured_at=stamp,
                                      version="1.0.0", note=note) for row in rows)
        print(("Attested: " if options.write else "Validated: ") + path)
    print(f"{len(paths)} runtime files and {len(evidence)} results validated; "
          "0.96.1 evidence remains unattested and stale.")


if __name__ == "__main__":
    main()
