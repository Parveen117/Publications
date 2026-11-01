"""Python 3.12 I/O; scientific arithmetic is in the unchanged native engines."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

if sys.version_info[:2] != (3, 12):
    raise SystemExit('Use Python 3.12')
if len(sys.argv) != 2 or sys.argv[1] not in ('--check', '--write'):
    raise SystemExit('Usage: python3.12 run.py --check|--write')
root = Path(__file__).resolve().parent
pins = json.loads((root / 'SOURCE_PINS.json').read_text())
for relative, expected in pins['files'].items():
    if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
        raise SystemExit('Source pin mismatch: ' + relative)
subprocess.run(['node', str(root / 'verify.cjs'), sys.argv[1]], check=True)
