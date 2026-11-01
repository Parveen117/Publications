"""Orchestrate the frozen nonuniform-source native certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

HERE = Path(__file__).resolve().parent


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--protocol-commit', required=True)
    args = parser.parse_args()
    if len(args.protocol_commit) != 40 or any(c not in '0123456789abcdef' for c in args.protocol_commit):
        parser.error('full protocol commit SHA required')
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'status.json', dict(status='running', protocol_commit=args.protocol_commit))
    try:
        pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
        for name, expected in pins['sha256'].items():
            if hashlib.sha256((HERE / name).read_bytes()).hexdigest() != expected:
                raise ValueError('inherited source hash mismatch: ' + name)
        (out / 'request.json').write_bytes((HERE / 'request.json').read_bytes())
        sources = ('source.cjs', 'solve.cjs', 'verify.cjs', 'run.py', 'request.json',
                   'THEOREM.md', 'SOURCE_MAP.md', 'SOURCE_PINS.json')
        hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in sources}
        before = time.perf_counter()
        subprocess.run(['node', str(HERE / 'solve.cjs'), str(out)], check=True, timeout=900)
        elapsed = time.perf_counter() - before
        if hashes != {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in sources}:
            raise ValueError('execution sources changed during evaluation')
        report = json.loads((out / 'verification.json').read_text())
        if report['status'] != 'PASS':
            raise RuntimeError('verification did not pass')
        metadata = dict(protocol_commit=args.protocol_commit, frozen_parent=pins['frozen_parent'],
                        python=platform.python_version(),
                        node=subprocess.check_output(['node', '--version'], text=True).strip(),
                        native_solver_seconds=elapsed, source_sha256=hashes,
                        fresh_nrqc_requests=0, reused_nrqc_requests=9,
                        nrqc_commit='fdd39048a6ddaa2dd78aa3197c274d7766f31a0e',
                        hardware_called=False, ci_dispatched=False,
                        boundary='One nonuniform source and a total-label marginal, with retained relative angle and a complete first-mode error certificate; no full field-state circuit or autonomous gravity.')
        write(out / 'execution.json', metadata)
        write(out / 'status.json', dict(status='completed', verification_passed=True,
                                        tolerance_unmet_cases=report['unresolved_cases']))
        write(out / 'manifest.json', {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in sorted(out.glob('*.json')) if p.name != 'manifest.json'})
        print(json.dumps(dict(status='completed', native_solver_seconds=elapsed)))
    except Exception as error:
        write(out / 'status.json', dict(status='incomplete', verification_passed=False, error=str(error)))
        raise


if __name__ == '__main__':
    main()
