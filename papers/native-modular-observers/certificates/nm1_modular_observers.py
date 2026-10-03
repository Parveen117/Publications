"""NM1 exact native/modular observer controls. Python 3.12 only.

Default and --check are read-only; --write changes only this stage's evidence.
Written general proofs are in THEOREMS.md; these are finite implementation checks.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
FOLDER = 'papers/native-modular-observers/'
PINS = HERE / 'NM1_SOURCE_PINS.json'
RESULT = HERE / 'NM1_RESULT.json'
EXPECTED = HERE / 'EXPECTED_NM1.sha256'
SOURCE = ROOT / 'papers/emk-ugd-algebra/certificates/emk1_determinant_seam_ladder.py'
spec = importlib.util.spec_from_file_location('nm1_native_emk1', SOURCE)
native = importlib.util.module_from_spec(spec)
spec.loader.exec_module(native)


def exact(x):
    if type(x) not in (int, F):
        raise TypeError('integer or Fraction required; no floats or bools')
    return F(x)


def integer(n):
    if type(n) is not int:
        raise TypeError('integer required')
    return n


def mat(A):
    if len(A) != 2 or any(len(row) != 2 for row in A):
        raise ValueError('2 by 2 matrix required')
    return tuple(tuple(exact(x) for x in row) for row in A)


def add(A, B):
    return tuple(tuple(A[i][j]+B[i][j] for j in range(2)) for i in range(2))


def scale(A, t):
    t = exact(t)
    return tuple(tuple(t*x for x in row) for row in A)


def mm(A, B):
    return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def det(A):
    return A[0][0]*A[1][1]-A[0][1]*A[1][0]


def inv(A):
    if not det(A):
        raise ValueError('singular matrix')
    return scale(((A[1][1], -A[0][1]), (-A[1][0], A[0][0])), 1/det(A))


def power(A, n):
    integer(n)
    if n < 0:
        return power(inv(A), -n)
    out = I
    while n:
        if n % 2:
            out = mm(out, A)
        A = mm(A, A); n //= 2
    return out


I = mat(((1, 0), (0, 1)))
R, K = mat(native.R2), mat(native.K2)
L = mm(R, K)
N = scale(add(L, scale(R, -1)), F(1, 2))
T = add(I, N)
S = add(I, R)
U = mat(((1, 1), (0, 1)))


def adapted(A):
    return mm(mm(inv(S), mat(A)), S)


def native_chart(A):
    return mm(mm(S, mat(A)), inv(S))


def sl2(A):
    A = mat(A)
    if any(x.denominator != 1 for row in A for x in row) or det(A) != 1:
        raise ValueError('SL(2,Z) in the declared adapted lattice required')
    return A


def sign(A):
    for row in A:
        for x in row:
            if x:
                return 1 if x > 0 else -1
    raise ValueError('zero has no projective representative')


def section(A):
    A = sl2(A)
    return scale(A, sign(A))


def star(g, h):
    return section(mm(section(g), section(h)))


def cocycle(g, h):
    return sign(mm(section(g), section(h)))


def lift(A):
    A = sl2(A)
    return sign(A), section(A)


def validate_lift(pair):
    epsilon, g = pair
    if type(epsilon) is not int or epsilon not in (-1, 1):
        raise ValueError('a lift requires a sign in {-1,1}')
    g = sl2(g)
    if section(g) != g:
        raise ValueError('canonical projective representative required')
    return epsilon, g


def unlift(pair):
    epsilon, g = validate_lift(pair)
    return scale(g, epsilon)


def compose(x, y):
    epsilon, g = validate_lift(x); eta, h = validate_lift(y)
    return epsilon*eta*cocycle(g, h), star(g, h)


def factor(A):
    """Exact Euclidean word, multiplied in the listed left-to-right order."""
    B = sl2(A); word = []
    while B[1][0]:
        a, c = int(B[0][0]), int(B[1][0])
        r = a % abs(c); q = (a-r)//c
        word.extend((('U', q), ('R', -1)))
        B = mm(R, mm(power(U, -q), B))
        assert B[1][0] == r and abs(r) < abs(c)
    epsilon = int(B[0][0])
    assert epsilon in (-1, 1) and B[1][1] == epsilon
    if epsilon == -1:
        word.append(('R', 2))
    word.append(('U', int(epsilon*B[0][1])))
    return word


def eval_word(word):
    out = I
    for name, n in word:
        if name not in ('R', 'U'):
            raise ValueError('unknown modular generator')
        out = mm(out, power(R if name == 'R' else U, n))
    return out


# Exact Gaussian rational arithmetic; no floating complex numbers.
def zpair(x, y=0):
    return exact(x), exact(y)


def za(z, w):
    return z[0]+w[0], z[1]+w[1]


def zs(z, w):
    return z[0]-w[0], z[1]-w[1]


def zm(z, w):
    return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]


def zi(z):
    norm = z[0]**2+z[1]**2
    if not norm:
        raise ValueError('zero Gaussian denominator')
    return z[0]/norm, -z[1]/norm


def zp(z, n):
    integer(n)
    if n < 0:
        return zp(zi(z), -n)
    out = zpair(1)
    while n:
        if n % 2:
            out = zm(out, z)
        z = zm(z, z); n //= 2
    return out


def upper(z):
    z = zpair(*z)
    if z[1] <= 0:
        raise ValueError('upper-half-plane point required')
    return z


def j(A, z):
    A, z = sl2(A), upper(z)
    return za(zm(zpair(A[1][0]), z), zpair(A[1][1]))


def act(A, z):
    A, z = sl2(A), upper(z)
    return zm(za(zm(zpair(A[0][0]), z), zpair(A[0][1])), zi(j(A, z)))


def slash(f, A, k, z):
    integer(k)
    return zm(zp(j(A, z), -k), f(act(A, z)))


def defect(f, A, k, z):
    return zs(slash(f, A, k, z), f(z))


def sample_group():
    # A deliberately finite set for controls, not an exhaustive infinite proof.
    return sorted({mm(mm(power(U, a), R), power(U, b))
                   for a in range(-2, 3) for b in range(-2, 3)} | {I, scale(I, -1), R})


def controls():
    zero = scale(I, 0)
    assert mm(K, K) == I and mm(R, R) == scale(I, -1)
    assert mm(K, R) == scale(L, -1) and mm(N, N) == zero and N != zero
    finite = {scale(A, e) for A in (I, R, K, L) for e in (-1, 1)}
    assert len(finite) == 8 and all(mm(A, B) in finite for A in finite for B in finite)
    assert T not in finite and adapted(R) == R and adapted(K) == mat(((1, 0), (0, -1)))
    assert adapted(T) == U and power(mm(R, T), 3) == scale(I, -1)
    assert det(S) == 2 and T[0][0].denominator == 2
    for n in range(-12, 13):
        assert power(T, n) == add(I, scale(N, n))
    # Compare the new rational implementation to unchanged native EMK assembly.
    basis = (I, K, R, L)
    for coeffs in itertools.product((-1, 0, 1), repeat=4):
        B = zero
        for t, A in zip(coeffs, basis):
            B = add(B, scale(A, t))
        assert B == mat(native.emk_block(*coeffs))
        for A in basis:
            assert mm(B, A) == mat(native.mmul(B, A))
    group = sample_group()
    factored = 0
    for entries in itertools.product(range(-3, 4), repeat=4):
        A = mat((entries[:2], entries[2:]))
        if det(A) == 1:
            assert eval_word(factor(A)) == A
            assert adapted(native_chart(A)) == A
            factored += 1
    for A, B in itertools.product(group, repeat=2):
        assert unlift(compose(lift(A), lift(B))) == mm(A, B)
        assert section(mm(A, B)) == section(mm(scale(A, -1), B))
    for A, B, C in itertools.product(group[:12], repeat=3):
        assert cocycle(A, B)*cocycle(star(A, B), C) == cocycle(B, C)*cocycle(A, star(B, C))
        assert compose(compose(lift(A), lift(B)), lift(C)) == compose(lift(A), compose(lift(B), lift(C)))
    g = section(R)
    assert star(g, g) == I and cocycle(g, g) == -1
    for e in (-1, 1):
        assert power(scale(g, e), 2) == scale(I, -1)
    continuation = mm(R, U)
    tr = lambda A: A[0][0]+A[1][1]
    assert tr(U) == tr(inv(U)) == 2
    assert (tr(mm(continuation, U)), tr(mm(continuation, inv(U)))) == (2, 0)
    assert power(R, 0) == power(R, 4) and lift(power(R, 0)) == lift(power(R, 4))
    points = [zpair(0, 1), zpair(F(1, 2), F(3, 2)), zpair(-2, F(1, 3))]
    pairs = list(itertools.product(group[:8], repeat=2))
    automorphy_checks = 0
    for A, B in pairs:
        for z in points:
            assert act(A, z)[1] == z[1]/(j(A, z)[0]**2+j(A, z)[1]**2)
            assert act(A, z) == act(scale(A, -1), z)
            assert act(mm(A, B), z) == act(A, act(B, z))
            assert j(mm(A, B), z) == zm(j(A, act(B, z)), j(B, z))
            sg, sh = section(A), section(B)
            left = zm(j(sg, act(sh, z)), j(sh, z))
            right = zm(zpair(cocycle(A, B)), j(star(A, B), z))
            assert left == right
            for k in (0, 1, 2, 3, -2):
                assert zp(left, k) == zm(zpair(F(cocycle(A, B))**k), zp(j(star(A, B), z), k))
            for k in (0, 2):
                f = lambda w: zp(w, 2)
                lhs = defect(f, mm(A, B), k, z)
                rhs = za(slash(lambda w: defect(f, A, k, w), B, k, z), defect(f, B, k, z))
                assert lhs == rhs
            automorphy_checks += 1
    for z in points:
        assert slash(lambda w: w, scale(I, -1), 1, z) == zm(zpair(-1), z)
        for A in group[:8]:
            f = lambda w: w
            for constant in (0, 1):
                correction = lambda w: zs(zpair(constant), w)
                assert za(defect(f, A, 0, z), defect(correction, A, 0, z)) == zpair(0)
    return {'finite_RK_group_order': 8, 'native_T': T, 'lattice_index': 2,
            'native_basis_comparisons': 81, 'integer_matrices_factored': factored,
            'sign_product_pairs': len(group)**2, 'cocycle_triples': 12**3,
            'automorphy_point_pairs': automorphy_checks,
            'sign_lift_minimum_labels_per_fiber': 2,
            'nonsplitting_witness': 'both lifts of [R] square to -I',
            'trace_return_witness': {'before': [2, 2], 'after': [2, 0]},
            'history_limit': 'R^0=R^4 even after the sign is retained',
            'completion_control': 'f(z)=z has corrections -z and 1-z at weight zero'}


SCOPE = {'written_results': ['NM1-T1', 'NM1-T2', 'NM1-T3', 'NM1-T4', 'NM1-T5'],
         'certification': 'written proofs and exact finite controls, not formal or external expert certification',
         'choices': ['rational lattice S Z^2', 'projective observer', 'classical analytic adapter in T4-T5'],
         'open': ['native selection of lattice and observer', 'history generating function',
                  'half-integral multiplier', 'mock theta shadow and completion', 'RH', 'Yang-Mills/Clay'],
         'curvature': 'no identification of section cocycle with connection curvature',
         'novelty': 'native docking and explicit controls; classical theory is credited'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_checks(root=ROOT, pins=None):
    pins = json.loads(PINS.read_text()) if pins is None else pins
    out = {}
    for path, sha in pins['upstream_sha256'].items():
        actual = digest(root/path)
        if actual != sha:
            raise ValueError('source drift: '+path)
        out[path] = actual
    for path in pins['local_inputs']:
        out[path] = digest(root/path)
    if root == ROOT:
        out[str(PINS.relative_to(ROOT))] = digest(PINS)
    return out


def serial(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [serial(v) for v in x]
    return x


def canonical(cert):
    return hashlib.sha256(json.dumps(cert, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run():
    if sys.version_info[:2] != (3, 12) or not __debug__:
        raise RuntimeError('Python 3.12 without optimization required')
    return serial({'certificate_type': 'NM1_NATIVE_MODULAR_OBSERVERS', 'verdict': 'PASS',
                   'source_sha256': source_checks(), 'controls': controls(), 'scope': SCOPE})


def check(cert, result=RESULT, expected=EXPECTED):
    if cert.get('scope') != SCOPE:
        raise ValueError('scope promotion')
    sha = canonical(cert)
    if json.loads(result.read_text()) != cert or expected.read_text().strip() != sha:
        raise ValueError('fresh NM1 result/digest mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true'); group.add_argument('--write', action='store_true')
    args = parser.parse_args(); cert = run(); sha = canonical(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n'); EXPECTED.write_text(sha+'\n')
    else:
        check(cert)
    print('NM1 PASS', sha)
    print('5 written results; native modular sign lift; analytic completion remains OPEN')


if __name__ == '__main__':
    main()
