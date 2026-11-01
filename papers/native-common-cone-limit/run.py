"""Replay exact native checks. Python only handles files and process execution."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

if sys.version_info[:2] != (3, 12):
    raise SystemExit('Use Python 3.12')
if len(sys.argv) != 2 or sys.argv[1] not in ('--check', '--write'):
    raise SystemExit('Usage: python3 run.py --check|--write')
root = Path(__file__).resolve().parent
pins = json.loads((root / 'SOURCE_PINS.json').read_text())
for name, expected in pins['packet_files'].items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit('Hash mismatch: ' + name)
subprocess.run(['node', str(root / 'verify.cjs'), sys.argv[1]], check=True)
