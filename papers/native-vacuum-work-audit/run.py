"""Execute the frozen finite work-accounting calibration on the pinned NRQC."""
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


def write(p, v):
    p.write_text(json.dumps(v, sort_keys=True, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--nrqc-root', type=Path, required=True)
    parser.add_argument('--protocol-commit', required=True)
    args = parser.parse_args()
    if len(args.protocol_commit) != 40 or any(c not in '0123456789abcdef' for c in args.protocol_commit):
        parser.error('full protocol SHA required')
    root, out = args.nrqc_root.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    write(out / 'status.json', dict(status='running', protocol_commit=args.protocol_commit))
    try:
        assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip() == PIN
        assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=root, text=True).strip()
        pins = json.loads((HERE / 'SOURCE_PINS.json').read_text())
        for p, expected in pins['sha256'].items():
            assert hashlib.sha256((HERE / p).read_bytes()).hexdigest() == expected, p
        source_names = ('PROTOCOL.md', 'build.cjs', 'cases.json', 'request.json', 'verify.cjs', 'run.py', 'SOURCE_PINS.json')
        hashes = {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in source_names}
        calc_files = sorted((root / 'native_relational_qc').rglob('*.py'))
        calc_hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in calc_files}
        (out / 'request.json').write_bytes((HERE / 'request.json').read_bytes())
        before = time.perf_counter()
        subprocess.run([sys.executable, '-m', 'native_relational_qc', 'calculate',
                        str(out / 'request.json'), '--out', str(out / 'calculator.result.json')],
                       cwd=root, check=True, timeout=40, stdout=subprocess.DEVNULL)
        subprocess.run(['node', str(HERE / 'verify.cjs'), str(out)], check=True, timeout=40)
        elapsed = time.perf_counter() - before
        assert hashes == {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest() for p in source_names}
        assert calc_hashes == {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in calc_files}
        write(out / 'execution.json', dict(protocol_commit=args.protocol_commit, nrqc_commit=PIN,
                source_sha256=hashes, calculator_source_sha256=calc_hashes,
                python=platform.python_version(), node=subprocess.check_output(['node', '--version'], text=True).strip(),
                wall_seconds=elapsed, fresh_nrqc_requests=1, hardware_called=False, ci_dispatched=False))
        write(out / 'status.json', dict(status='completed', verification_passed=True))
        write(out / 'manifest.json', {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(out.glob('*.json')) if p.name != 'manifest.json'})
        print(json.dumps(dict(status='completed', wall_seconds=elapsed)))
    except Exception as e:
        write(out / 'status.json', dict(status='incomplete', error=str(e)))
        raise


if __name__ == '__main__':
    main()
