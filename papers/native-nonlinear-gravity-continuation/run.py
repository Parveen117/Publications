"""I/O only: extract pinned original engines and launch the native verifier."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile

assert sys.version_info[:2] == (3, 12), 'Use Python 3.12 for this packet'
root = Path(__file__).resolve().parent
pins = json.loads((root / 'SOURCE_PINS.json').read_text())
archive = root.parent / 'Native_Phase_Word_Reconstruction_Source.zip'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == pins['native_archive_sha256']
assert len(sys.argv) == 2 and sys.argv[1] in ('--check', '--write')
with tempfile.TemporaryDirectory(prefix='native-nonlinear-gravity-') as tmp:
    with zipfile.ZipFile(archive) as bundle:
        for name, expected in pins['native_engine'].items():
            matches = [n for n in bundle.namelist() if n.endswith('/' + name)]
            assert len(matches) == 1, name
            data = bundle.read(matches[0])
            assert hashlib.sha256(data).hexdigest() == expected, name
            Path(tmp, name).write_bytes(data)
    subprocess.run(['node', str(root / 'verify.cjs'), '--native-root', tmp, sys.argv[1]], check=True)
