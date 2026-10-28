"""Verify the checked source archive against the committed release sources."""
import csv
import hashlib
import re
import subprocess
import tarfile
from pathlib import Path

base = Path('validation/release-0.96.1-2026-10-20')
expected_sha = (base / 'release-commit.txt').read_text().strip()
actual_sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip()
assert actual_sha == expected_sha
subprocess.run(['git', 'diff', '--exit-code', '--quiet'], check=True)
subprocess.run(['git', 'diff', '--cached', '--exit-code', '--quiet'], check=True)
index = {}
for line in subprocess.check_output(['git', 'ls-files', '--stage'], text=True).splitlines():
    metadata, path = line.split('\t', 1)
    index[path] = metadata.split()[1]


def dcf_fields(content):
    fields, current = {}, None
    for line in content.splitlines():
        if line[:1].isspace():
            fields[current] += line.strip()
        elif ':' in line:
            current, value = line.split(':', 1)
            fields[current] = value.strip()
    return {key: re.sub(r'\s+', '', value) for key, value in fields.items()}


manifest = list(csv.DictReader((base / 'release-source-manifest.csv').open()))
archive_path = base / 'latents_0.96.1.tar.gz'
checked, excluded, all_tracked = 0, 0, 0
with tarfile.open(archive_path) as archive:
    for item in manifest:
        path = item['file']
        current = Path(path).read_bytes()
        assert hashlib.md5(current).hexdigest() == item['md5'], path
        blob = b'blob ' + str(len(current)).encode() + b'\0' + current
        assert hashlib.sha1(blob).hexdigest() == index[path], path
        if path.startswith(('tests/equivalence/', 'equivalence/', '.github/')):
            excluded += 1  # Explicit existing .Rbuildignore exclusions.
            continue
        packaged = archive.extractfile('latents/' + path).read()
        if path == 'DESCRIPTION':
            original, built = dcf_fields(current.decode()), dcf_fields(packaged.decode())
            assert all(built[key] == value for key, value in original.items())
            assert built['Version'] == '0.96.1'
            assert set(built) - set(original) == {
                'NeedsCompilation', 'Packaged', 'Author', 'Maintainer'}
        else:
            assert packaged == current, path
        checked += 1
    for member in archive.getmembers():
        path = member.name.removeprefix('latents/')
        if not member.isfile() or path not in index:
            continue
        packaged = archive.extractfile(member).read()
        committed = subprocess.check_output(['git', 'show', 'HEAD:' + path])
        if path == 'DESCRIPTION':
            original, built = dcf_fields(committed.decode()), dcf_fields(packaged.decode())
            assert all(built[key] == value for key, value in original.items())
        else:
            assert packaged == committed, path
        all_tracked += 1
print('Release commit:', actual_sha)
print('Tracked working tree and index: clean')
print('Manifest files match committed content:', len(manifest))
print('Archive matches committed package files:', checked)
print('All tracked archive files independently match Git content:', all_tracked)
print('Development-only manifest files excluded by .Rbuildignore:', excluded)
print('Archive SHA256:', hashlib.sha256(archive_path.read_bytes()).hexdigest())
