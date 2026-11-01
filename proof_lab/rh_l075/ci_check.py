"""Reproduction check for the BRT-1 RH positivity packets.

Reruns the committed certificates in a temporary copy of each packet and
compares the fresh result with the committed evidence field by field
(everything except the interpreter version string must agree), after
checking that every committed result carries the SHA-256 of the engine
that is in the tree. Nothing in the repository is written.

    python3 ci_check.py            # fast packets, both engines (minutes)
    python3 ci_check.py --full     # also the L=9/10, 19/20 and 1 certificates (hours)

Exit status 0 means every rerun reproduced its committed evidence.
"""
import hashlib, json, shutil, subprocess, sys, tempfile, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
L055 = ROOT / 'rh_l055'
FULL = '--full' in sys.argv[1:]

# packet -> (committed result file, engine, extra flags)
ARB_FAST = [
    ('l055_control', 'result-arb-220.json', []),
    ('l0625_route', 'result-arb-160.json', []),
    ('t53_route', 'result-arb-220.json', []),
    ('boundary_route', 'result-arb-220.json', []),
    ('l08_route', 'result-arb-220.json', []),
    ('l085_route', 'result-arb-230.json', []),
    ('near_null', 'rayleigh-arb-160.json', ['--rayleigh']),
    ('near_null_08', 'rayleigh-arb-160.json', ['--rayleigh']),
    ('near_null_085', 'rayleigh-arb-160.json', ['--rayleigh']),
    ('near_null_09', 'rayleigh-arb-160.json', ['--rayleigh']),
    ('near_null_095', 'rayleigh-arb-200.json', ['--rayleigh']),
    ('near_null_1', 'rayleigh-arb-300.json', ['--rayleigh']),
]
ARB_FULL = [
    ('l09_route', 'result-arb-260.json', []),
    ('l095_route', 'result-arb-320.json', []),
    ('l1_route', 'result-arb-420.json', []),
]
DECIMAL_FAST = [
    ('l055_control', 'result-220.json', []),
    ('near_null', 'rayleigh-160.json', ['--rayleigh']),
]
IGNORED_KEYS = {'python'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fail(msg):
    print('FAIL:', msg)
    sys.exit(1)


def precision_of(result_name):
    return result_name.rsplit('-', 1)[1].split('.')[0]


def compare(fresh, committed, label):
    a = json.loads(Path(fresh).read_text())
    b = json.loads(Path(committed).read_text())
    for k in IGNORED_KEYS:
        a.pop(k, None); b.pop(k, None)
    if a == b:
        print(f'  reproduced: {label} ({a.get("status")})')
        return
    keys = sorted(set(a) | set(b))
    diff = [k for k in keys if a.get(k) != b.get(k)]
    fail(f'{label}: fresh result differs from committed evidence in keys {diff}')


def check_hashes():
    arb_hash = sha256(HERE / 'verify_arb.py')
    dec_hash = sha256(HERE / 'verify.py')
    n = 0
    for f in sorted(HERE.glob('*/result*.json')) + sorted(HERE.glob('*/rayleigh*.json')):
        r = json.loads(f.read_text())
        want = arb_hash if '-arb-' in f.name else dec_hash
        if r.get('code_sha256') != want:
            fail(f'{f.relative_to(ROOT)} carries code_sha256 {r.get("code_sha256")[:12]}..., the engine in the tree is {want[:12]}...')
        n += 1
    l055 = json.loads((L055 / 'result-220.json').read_text())
    if l055['code_sha256'] != sha256(L055 / 'verify.py'):
        fail('rh_l055/result-220.json does not carry the hash of rh_l055/verify.py')
    print(f'engine hashes: {n} rh_l075 evidence files and rh_l055/result-220.json carry the hash of the engine in the tree')


def run_packet(engine, packet, result_name, flags):
    src = HERE / packet
    with tempfile.TemporaryDirectory() as tmp:
        for name in ('config.json', 'trial.json'):
            shutil.copy(src / name, Path(tmp) / name)
        cmd = [sys.executable, '-B', str(HERE / engine), precision_of(result_name), '--packet', tmp] + flags
        t = time.time()
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            print(p.stdout[-2000:]); print(p.stderr[-2000:])
            fail(f'{engine} on {packet} exited with status {p.returncode}')
        print(f'  {engine} {packet} ran in {time.time()-t:.0f} s')
        compare(Path(tmp) / result_name, src / result_name, f'{packet}/{result_name}')


def run_l055():
    with tempfile.TemporaryDirectory() as tmp:
        for name in ('verify.py', 'trial.json'):
            shutil.copy(L055 / name, Path(tmp) / name)
        t = time.time()
        p = subprocess.run([sys.executable, '-B', str(Path(tmp) / 'verify.py'), '220'], capture_output=True, text=True, cwd=tmp)
        if p.returncode != 0:
            print(p.stdout[-2000:]); print(p.stderr[-2000:])
            fail('rh_l055/verify.py 220 failed')
        print(f'  rh_l055 verify.py 220 ran in {time.time()-t:.0f} s')
        compare(Path(tmp) / 'result-220.json', L055 / 'result-220.json', 'rh_l055/result-220.json')


def main():
    check_hashes()
    print('rh_l055 (L=11/20, standard-library engine):')
    run_l055()
    print('rh_l075, standard-library engine:')
    for packet, result_name, flags in DECIMAL_FAST:
        run_packet('verify.py', packet, result_name, flags)
    print('rh_l075, ball-arithmetic engine:')
    for packet, result_name, flags in ARB_FAST + (ARB_FULL if FULL else []):
        run_packet('verify_arb.py', packet, result_name, flags)
    if not FULL:
        print('skipped (committed evidence only, rerun with --full): ' + ', '.join(p for p, _, _ in ARB_FULL))
    print('ALL_PACKETS_REPRODUCED')


if __name__ == '__main__':
    main()
