"""Python 3.12 I/O only; native JS source/baseline and pinned NRQC execution."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PIN = 'fdd39048a6ddaa2dd78aa3197c274d7766f31a0e'


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, capture_output=True, **kwargs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nrqc', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nrqc, out = args.nrqc.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    assert run(['git', 'rev-parse', 'HEAD'], cwd=nrqc).stdout.strip() == PIN
    assert not run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=nrqc).stdout.strip()
    source = run(['node', str(HERE / 'source_intervals.cjs')])
    (out / 'source_intervals.json').write_text(source.stdout)
    words = [{'0': 'X'}, {'1': 'Z', '0': 'X'}, {'0': 'Y'},
             {'1': 'Z', '0': 'Y'}, {'1': 'X'}, {'1': 'X', '0': 'Z'}]
    # Controlled S and S-dagger are genuinely non-Clifford phase circuits.
    # Their explicit final states are independently specified in the JS baseline.
    states = {
        'vacuum': [],
        'bell': [['h', [0]], ['cx', [0, 1]]],
        'controlled_s': [['h', [0]], ['h', [1]], ['t', [0]], ['t', [1]],
                         ['cx', [0, 1]], ['tdg', [1]], ['cx', [0, 1]]],
        'controlled_s_dagger': [['h', [0]], ['h', [1]], ['tdg', [0]], ['tdg', [1]],
                                ['cx', [0, 1]], ['t', [1]], ['cx', [0, 1]]],
    }
    execution = []
    for state, operations in states.items():
        for backend in ['observer', 'clifford_t']:
            stem = state + '_' + backend
            request = {'schema': 'nrqc-calculator-v1', 'n_qubits': 2,
                       'backend': backend, 'operations': operations,
                       'queries': [{'id': 'six_pauli_basis', 'kind': 'geometry',
                                    'generators': [{'paulis': word} for word in words]}]}
            req, res = out / (stem + '.request.json'), out / (stem + '.result.json')
            req.write_text(json.dumps(request, indent=2) + '\n')
            before = time.perf_counter()
            run([sys.executable, '-m', 'native_relational_qc', 'calculate', str(req),
                 '--timeout', '20', '--out', str(res)], cwd=nrqc)
            execution.append({'state': state, 'backend': backend,
                              'command_wall_seconds': time.perf_counter() - before})
    sources = ['source_intervals.cjs', 'verify.cjs', 'run_validation.py',
               '../native-common-cone-limit/native_engine/native_operator.cjs']
    meta = {'nrqc_commit': PIN, 'python': platform.python_version(),
            'node': run(['node', '--version']).stdout.strip(), 'execution': execution,
            'source_sha256': {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                              for name in sources},
            'hardware_called': False, 'benchmark_claim': False,
            'boundary': 'NRQC computes finite Pauli means/covariances; native JS assembles supplied source controls and pulse shears.'}
    (out / 'execution.json').write_text(json.dumps(meta, indent=2) + '\n')
    result = run(['node', str(HERE / 'verify.cjs'), str(out)])
    (out / 'verification.json').write_text(result.stdout)
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(out.glob('*.json')) if p.name != 'manifest.json'}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(result.stdout)


if __name__ == '__main__':
    main()
