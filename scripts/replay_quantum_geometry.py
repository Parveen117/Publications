"""Read-only replay of the integrated quantum geometry research chain."""
from pathlib import Path
import json
import subprocess
import sys

if sys.version_info[:2] != (3, 12):
    raise SystemExit('Use Python 3.12')
root = Path(__file__).resolve().parents[1]
stages = [
    'native-jacobian-cut-certification',
    'native-angular-calibration',
    'native-exponential-gravity-response',
    'native-nonlinear-gravity-continuation',
    'native-cut-record-completeness',
    'native-quantum-gravity-consistency',
    'native-four-dimensional-hinge',
    'native-matter-geometry-holonomy',
    'native-common-cone-limit',
    'native-geometry-controlled-propagation',
    'native-action-geometry-dynamics',
    'native-zero-reference-recovery',
    'native-invariant-ratio-transport',
    'native-cut-local-dynamics',
    'native-relational-cut-space',
]
for stage in stages:
    proc = subprocess.run([sys.executable, str(root / 'papers' / stage / 'run.py'), '--check'], cwd=root, text=True, capture_output=True)
    if proc.returncode:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f'FAILED: {stage}')
    report = json.loads(proc.stdout)
    print(json.dumps({'stage': stage, 'status': report.get('status', 'PASS' if report.get('passed') else 'UNKNOWN'), 'exact_checks': report.get('exact_checks'), 'symbolic_certificates': report.get('symbolic_certificates', 0), 'negative_controls': report.get('negative_controls')}), flush=True)
print(json.dumps({'passed': True, 'stages': len(stages), 'mode': 'read-only replay'}))
