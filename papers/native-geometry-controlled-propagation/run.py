"""Python 3.12 I/O only; native JS owns scientific arithmetic and proof replay."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

if sys.version_info[:2] != (3, 12):
    raise SystemExit('Use Python 3.12')
if len(sys.argv) != 2 or sys.argv[1] not in ('--write', '--check'):
    raise SystemExit('Usage: python3 run.py --check|--write')
root = Path(__file__).resolve().parent
pins = json.loads((root / 'SOURCE_PINS.json').read_text())
for relative, expected in pins['files'].items():
    actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit('Source pin mismatch: ' + relative)
subprocess.run(['node', str(root / 'verify.cjs'), sys.argv[1]], check=True)
