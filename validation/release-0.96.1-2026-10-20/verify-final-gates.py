"""Refuse release sign-off unless every recorded final-revision gate passes."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

base = Path('validation/release-0.96.1-2026-10-20')
commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert commit == (base / 'release-commit.txt').read_text().strip()
rows = []
reference_commit = (base / 'outer-final/commit.txt').read_text().strip()
allowed_delta = {'.github/workflows/R-CMD-check.yaml', 'tests/testthat/test-engine-fixes.R'}
actual_delta = set(subprocess.check_output(
    ['git', 'diff', '--name-only', reference_commit, commit], text=True).splitlines())
assert actual_delta == allowed_delta, actual_delta
snapshot = json.loads((base / 'outer-final/snapshot-manifest.json').read_text())
tracked = set(subprocess.check_output(['git', 'ls-files'], text=True).splitlines())
assert tracked == {item['path'] for item in snapshot}
for item in snapshot:
    if item['path'] not in allowed_delta:
        assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256'], item['path']


def count(value):
    if value in ('TRUE', 'FALSE'):
        return int(value == 'TRUE')
    assert value.isdigit(), value
    return int(value)


def check_results(name, relative):
    cases = list(csv.DictReader((base / relative).open()))
    assert cases, name
    totals = {key: sum(count(case[key]) for case in cases)
              for key in ('passed', 'failed', 'error', 'warning', 'skipped')}
    assert not any(totals[key] for key in ('failed', 'error', 'warning', 'skipped')), (name, totals)
    rows.append({'gate': name, 'result': 'passed', 'cases': len(cases), **totals})


for kind in ('unit', 'equivalence', 'references'):
    log = (base / f'{kind}-final.log').read_text()
    execution_commit = commit if kind == 'unit' else reference_commit
    assert f'Commit: {execution_commit}' in log and 'Production source unchanged: TRUE' in log
    fingerprints = list(csv.DictReader((base / f'{kind}-final-source.csv').open()))
    assert len(fingerprints) == len(list(Path('R').glob('*.R')))
    assert all(hashlib.md5(Path(item['file']).read_bytes()).hexdigest() == item['md5']
               for item in fingerprints)
    check_results(kind, f'{kind}-final.csv')

assert (base / 'outer-final/commit.txt').read_text().strip() == reference_commit
for item in csv.DictReader((base / 'outer-final/production-before.csv').open()):
    assert hashlib.md5(Path(item['path']).read_bytes()).hexdigest() == item['md5']
assert 'Scratch production source unchanged: TRUE' in (base / 'outer-final/full-suite.log').read_text()
check_results('complete external equivalence', 'outer-final/cases.csv')

for run_id, expected_jobs in ((39456155958, 6), (39456160271, 1), (39456165092, 1)):
    run = json.loads((base / f'ci-final-d4c9f69/{run_id}-latest.json').read_text())
    assert run['headSha'] == commit and run['conclusion'] == 'success', run_id
    assert len(run['jobs']) == expected_jobs
    assert all(job['conclusion'] == 'success' for job in run['jobs'])
    rows.append({'gate': f"GitHub {run_id}", 'result': 'passed', 'jobs': expected_jobs})

check_log = (base / 'check-final/latents.Rcheck/00check.log').read_text()
assert '* DONE' in check_log and 'Status: 1 NOTE' in check_log
assert 'unable to verify current time' in check_log
assert ' ERROR' not in check_log and ' WARNING' not in check_log
rows.append({'gate': 'native source check', 'result': '0 errors; 0 warnings; 1 time NOTE'})
smoke = (base / 'installed-public-proofs.log').read_text()
assert 'Installed version: 0.96.1' in smoke
assert 'Installed archive: typed worker failure, zero cleanup warnings, no leaked sockets.' in smoke
assert '* DONE (latents)' in (base / 'install-final-archive.log').read_text()
rows.append({'gate': 'installed archive public and worker-death reproductions', 'result': 'passed'})
subprocess.run(['python3', str(base / 'verify-release.py')], check=True)
report = {'commit': commit, 'reference_execution_commit': reference_commit,
          'reference_runtime_inputs_identical': True, 'test_ci_delta': sorted(actual_delta), 'verified_utc': datetime.now(timezone.utc).isoformat(), 'gates': rows}
(base / 'final-validation.json').write_text(json.dumps(report, indent=2) + '\n')
fields = ['gate', 'result', 'cases', 'passed', 'failed', 'error', 'warning', 'skipped', 'jobs']
with (base / 'final-validation.csv').open('w') as destination:
    writer = csv.DictWriter(destination, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(report, indent=2))
