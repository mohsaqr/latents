import csv
import hashlib
import re
import tarfile
from pathlib import Path

root = Path('.')
report = root / 'validation/fixes-2026-10-20/FIXES.md'
assert 'FINAL_VALIDATION_PENDING' not in report.read_text()
links = re.findall(r'\]\(([^)]+)\)', report.read_text())
missing = [target for target in links
           if '://' not in target and not (report.parent / target.split('#')[0]).exists()]
assert not missing, missing
print('Report local links resolve:', len(links))
manifest = list(csv.DictReader((report.parent / 'final-delivery-manifest.csv').open()))
archive_path = root / 'latents_0.96.0.tar.gz'
if not archive_path.exists():
    archive_path = report.parent / 'latents_0.96.0.tar.gz'
checked, excluded = 0, 0

def dcf_fields(content):
    fields, current = {}, None
    for line in content.splitlines():
        if line[:1].isspace():
            fields[current] += line.strip()
        elif ':' in line:
            current, value = line.split(':', 1)
            fields[current] = value.strip()
    return {key: re.sub(r'\s+', '', value) for key, value in fields.items()}

with tarfile.open(archive_path) as archive:
    for item in manifest:
        assert hashlib.md5((root / item['file']).read_bytes()).hexdigest() == item['md5'], item['file']
        if item['file'].startswith('tests/equivalence/'):
            excluded += 1  # Explicit existing .Rbuildignore exclusion; tested from source.
            continue
        payload = archive.extractfile('latents/' + item['file']).read()
        if item['file'] == 'DESCRIPTION':
            original = dcf_fields((root / item['file']).read_text())
            built = dcf_fields(payload.decode())
            assert all(built[key] == value for key, value in original.items())
            assert set(built) - set(original) == {
                'NeedsCompilation', 'Packaged', 'Author', 'Maintainer'}
            checked += 1
            continue  # R CMD build wraps DCF and adds normal package metadata.
        assert hashlib.md5(payload).hexdigest() == item['md5'], item['file']
        checked += 1
print('Working files match retained manifest:', len(manifest))
print('Delivery archive matches:', checked, 'files;', excluded,
      'equivalence source files excluded by .Rbuildignore')
print('Source archive SHA256:', hashlib.sha256(archive_path.read_bytes()).hexdigest())
summary = list(csv.DictReader((report.parent / 'validation-summary.csv').open()))
assert [(x['suite'], x['passed']) for x in summary] == [('unit', '8104'), ('equivalence', '405')]
assert all(all(int(x[key]) == 0 for key in ('failed', 'error', 'warning', 'skipped'))
           for x in summary)
print('Final unit/equivalence summary: verified')
