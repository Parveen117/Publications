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
    'native-local-coframe-covariance',
    'native-local-normal-closure',
    'native-relational-event-clock',
    'native-singular-response-transport',
    'native-moving-zero-covariance',
    'native-heat-observer-memory',
    'native-cut-generated-history',
    'native-gravity-signalling-audit',
    'native-homogeneous-quantum-sector',
    'native-local-quantum-anomaly',
    'native-local-physical-sector',
    'native-local-sector-classification',
    'native-complementary-clock-sector',
    'native-critical-origin-sector',
    'native-global-observable-algebra',
    'native-relational-pointer-measurement',
    'native-regular-apparatus-lift',
    'native-material-reference-field',
    'native-material-profile-quantization',
    'native-inverse-gauge-quantization',
    'native-central-gauge-domain',
    'native-residual-gauge-rigging',
    'native-resolved-material-geometry',
    'native-coarse-material-geometry',
    'native-material-detector-dilation',
    'native-material-reference-momentum',
    'native-logarithmic-material-density',
    'native-tvsp-lambda-matrices',
    'native-tvsp-exponential-transport',
    'native-quantum-response-transport',
    'native-density-cut-dynamics',
    'native-phase-current-response',
    'native-wick-calibrated-geometry',
    'native-prime-observer-sectors',
    'native-prime-observation-closure',
    'native-prime-innovation-window',
    'native-infinite-layer-cut',
    'native-layer-continuation-window',
    'native-layer-anchor-transport',
    'native-layer-anchor-atlas',
    'native-layer-depth-compatibility',
    'native-layer-observer-transport',
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
