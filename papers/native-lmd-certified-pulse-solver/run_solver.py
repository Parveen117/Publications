"""Launch pinned NRQC phase queries; all scientific processing is native JS."""
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
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(result.stderr or result.stdout)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--nrqc', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--request', type=Path, default=HERE / 'request.json')
    args = parser.parse_args()
    nrqc, out = args.nrqc.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    assert run(['git', 'rev-parse', 'HEAD'], cwd=nrqc).stdout.strip() == PIN
    assert not run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=nrqc).stdout.strip()
    (out / 'solver.request.json').write_bytes(args.request.read_bytes())
    executions = []
    # Independent phase probes; at most two qubits per request.
    for block, start in enumerate(range(0, 17, 2)):
        labels = list(range(start, min(start + 2, 17)))
        operations, queries = [], []
        for q, phase in enumerate(labels):
            operations.append(['h', [q]])
            operations.extend([['p64', [q]]] * phase)
            for pauli in ['X', 'Y']:
                queries.append({'id': f'phase_{phase}_{pauli}', 'kind': 'expectation',
                                'paulis': {str(q): pauli}})
        request = {'schema': 'nrqc-calculator-v1', 'n_qubits': len(labels),
                   'backend': 'cyclotomic64', 'decimal_places': 24,
                   'operations': operations, 'queries': queries}
        req, res = out / f'phase_{block}.request.json', out / f'phase_{block}.result.json'
        req.write_text(json.dumps(request, indent=2) + '\n')
        before = time.perf_counter()
        run([sys.executable, '-m', 'native_relational_qc', 'calculate', str(req),
             '--timeout', '30', '--out', str(res)], cwd=nrqc)
        executions.append({'block': block, 'phase_indices': labels,
                           'command_wall_seconds': time.perf_counter() - before})
    source_paths = ['interval.cjs', 'solve.cjs', 'run_solver.py',
                    '../native-common-cone-limit/native_engine/native_operator.cjs',
                    '../native-lmd-pulse-characteristic/results/source_intervals.json']
    start = time.perf_counter()
    run(['node', str(HERE / 'solve.cjs'), str(out)])
    metadata = {'nrqc_commit': PIN, 'python': platform.python_version(),
                'node': run(['node', '--version']).stdout.strip(), 'executions': executions,
                'native_solver_command_seconds': time.perf_counter() - start,
                'source_sha256': {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
                                  for p in source_paths},
                'boundary': 'Native conditional-expectation solver and error bounds; NRQC computes exact phase-readout nodes. No field state preparation or full time-evolution circuit.',
                'hardware_called': False, 'advantage_claimed': False}
    (out / 'execution.json').write_text(json.dumps(metadata, indent=2) + '\n')
    manifest = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(out.glob('*.json')) if p.name != 'manifest.json'}
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    report = json.loads((out / 'verification.json').read_text())
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
