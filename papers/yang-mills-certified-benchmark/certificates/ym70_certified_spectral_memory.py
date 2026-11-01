"""YM70: native leading Wilson memory and full spectral enclosures.

Python 3.12 only. Exact finite controls accompany written proofs; the
interacting vacuum and its full moment matrices are not computed here.
Default/--check are read-only. --write changes only this result and hash.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym64_vacuum_adapted_tower as previous

native, p, a = previous.native, previous.p, previous.a
model, exact = previous.model, previous.exact
RESULT = HERE/'YM70_RESULT.json'
PIN = HERE/'EXPECTED_YM70.sha256'
SOURCES = HERE/'YM70_SOURCE_PINS.json'


def nonnegative(value):
    value = exact(value)
    if value < 0:
        raise ValueError('Nonnegative exact scalar required')
    return value


def positive(value):
    value = exact(value)
    if value <= 0:
        raise ValueError('Positive exact scalar required')
    return value


def chebyshev(order, value):
    native.width(order, 0)
    value = exact(value)
    left, right = F(1), value
    for _ in range(order):
        left, right = right, 2*value*right-left
    return left


def error_factor(order, damping, gap):
    native.width(order, 0)
    damping, gap = nonnegative(damping), positive(gap)
    if not damping:
        return F(0)
    return 1/chebyshev(order+1, 1+2*gap/damping)**2


def psd_basis(gram):
    """Independent principal columns by exact positive Schur elimination."""
    gram = previous.matrix(gram)
    if len(gram) != len(gram[0]) or gram != a.transpose(gram):
        raise ValueError('Symmetric square Gram matrix required')
    residual = [row[:] for row in gram]
    indices = []
    while True:
        if any(residual[i][i] < 0 for i in range(len(gram))):
            raise ValueError('Gram matrix is not positive')
        pivot = next((i for i in range(len(gram)) if residual[i][i] > 0), None)
        if pivot is None:
            if any(v for row in residual for v in row):
                raise ValueError('Nonzero residual with zero diagonal')
            return tuple(indices)
        indices.append(pivot)
        d = residual[pivot][pivot]
        residual = [[residual[i][j]-residual[i][pivot]*residual[pivot][j]/d
                     for j in range(len(gram))] for i in range(len(gram))]


def quotient_form(gram, pencil, rhs):
    """rhs* pencil^dagger rhs on the fixed Gram quotient, without floats."""
    gram, pencil, rhs = map(previous.matrix, (gram, pencil, rhs))
    if (len(pencil) != len(gram) or len(pencil[0]) != len(gram)
            or len(rhs) != len(gram) or pencil != a.transpose(pencil)):
        raise ValueError('Compatible symmetric pencil and right-hand side required')
    basis = psd_basis(gram)
    if not basis:
        if any(v for row in pencil+rhs for v in row):
            raise ValueError('Nonzero input on a zero Gram space')
        return a.scale(a.identity(len(rhs[0])), 0), 0
    G = previous.sub(gram, basis)
    K = previous.sub(pencil, basis)
    b = previous.sub(rhs, basis, range(len(rhs[0])))
    embed = a.multiply(previous.sub(gram, range(len(gram)), basis), a.inverse(G))
    if (a.multiply(embed, b) != rhs or
            a.multiply(a.multiply(embed, K), a.transpose(embed)) != pencil):
        raise ValueError('Pencil or rhs is inconsistent with the Gram kernel')
    if len(psd_basis(K)) != len(basis):
        raise ValueError('Pencil must be positive on the Gram quotient')
    return a.multiply(a.multiply(a.transpose(b), a.inverse(K)), b), len(basis)


def enclose(moments, order, damping, gap):
    """Bounds under the exact X-moment contract of YM70-T2.

    The finite positivity checks are necessary checks, not a reconstruction
    theorem for an arbitrary supplied infinite moment sequence.
    """
    native.width(order, 0)
    damping, gap = nonnegative(damping), positive(gap)
    if len(moments) < 2*order+3:
        raise ValueError('Moments through index 2N+2 are required')
    M = [previous.matrix(v) for v in moments[:2*order+3]]
    width = len(M[0])
    if any(len(v) != width or len(v[0]) != width or v != a.transpose(v) for v in M):
        raise ValueError('Same-size symmetric moment matrices required')
    G, G1, G2 = [previous.hankel(M, order+1, shift) for shift in (0, 1, 2)]
    for value in (G, G1, G2, a.add(a.scale(G, 1/gap), G1, -1),
                  a.add(a.scale(G1, 1/gap), G2, -1)):
        psd_basis(value)
    b = [row[:] for value in M[:order+1] for row in value]
    b1 = [row[:] for value in M[1:order+2] for row in value]
    lower, rank = quotient_form(G, a.add(G, G1, damping), b)
    weighted, rank1 = quotient_form(G1, a.add(G1, G2, damping), b1)
    upper = a.add(M[0], weighted, -damping)
    return dict(lower=lower, upper=upper, rank=rank, weighted_rank=rank1,
                error_factor=error_factor(order, damping, gap))


def leading_atoms(factor):
    native.width(factor, 3)
    atoms = []
    for k, count in ((1, 16*factor-8), (2, 4)):
        atoms.extend(((F(6*factor+6-3*k, 2), F(count*k*k, 256)),
                      (F(6*factor+6+k, 2), F(count*k*k, 768))))
    return tuple(sorted(atoms))


def atom_moments(atoms, count):
    native.width(count)
    atoms = tuple((positive(rate), positive(weight)) for rate, weight in atoms)
    if not atoms:
        raise ValueError('At least one positive atom required')
    return [[[sum((weight/rate**(n+1) for rate, weight in atoms), F(0))]]
            for n in range(count)]


def atom_response(atoms, damping):
    damping = nonnegative(damping)
    return [[sum((positive(weight)/(positive(rate)+damping)
                  for rate, weight in atoms), F(0))]]


@lru_cache(maxsize=None)
def native_controls():
    rows = []
    for k in (1, 2):
        width = k+2
        shared = tuple((i, 1) for i in range(k))
        left = model.loop_polynomial(shared+((k, 1),), width)
        right = model.loop_polynomial(shared+((k+1, 1),), width)
        product = p.mul(left, right)
        singlet = p.scale(model.loop_polynomial(((k, 1), (k+1, -1)), width), F(1, 4))
        triplet = p.add(product, p.scale(singlet, -1))
        norm = lambda v: native.phi(p.mul(v, v), width)
        assert (norm(product), norm(singlet), norm(triplet)) == (F(1, 16), F(1, 64), F(3, 64))
        assert native.phi(p.mul(singlet, triplet), width) == 0
        response = p.add(*(p.mul(native.derivative(left, e, axis),
                                native.derivative(right, e, axis))
                           for e in range(k) for axis in (1, 2, 3)))
        expected = p.scale(p.add(p.scale(singlet, 3), p.scale(triplet, -1)), F(k, 4))
        assert norm(p.add(response, p.scale(expected, -1))) == 0
        for e in range(width):
            for poly, rate in ((singlet, F(0) if e < k else F(3, 4)),
                               (triplet, F(2) if e < k else F(3, 4))):
                electric = p.scale(p.add(*(native.derivative(native.derivative(poly, e, axis), e, axis)
                                          for axis in (1, 2, 3))), -1)
                assert norm(p.add(electric, p.scale(poly, -rate))) == 0
        assert norm(response) == F(3*k*k, 256)
        rows.append(dict(shared_links=k, singlet_norm_squared=norm(singlet),
                         triplet_norm_squared=norm(triplet), response_norm_squared=norm(response),
                         native_derivative_and_individual_link_Casimir_checks=True))
    geometry_rows = []
    for factor in (3, 4, 5):
        geometry = previous.previous.blocked_geometry(4, factor)
        word = geometry['coarse']['faces'][0]
        active = {e for link, sign in word for e in geometry['paths'][link]}
        touched = [{e for e, sign in face} for face in geometry['fine']['faces']
                   if active.intersection(e for e, sign in face)]
        counts = Counter(len(active & face) for face in touched)
        assert counts == {1: 16*factor-8, 2: 4}
        assert all(len((face-active)-other) > 0 for i, face in enumerate(touched)
                   for other in touched[i+1:])
        response_norm = sum((F(3*k*k*count, 256) for k, count in counts.items()), F(0))
        atoms = leading_atoms(factor)
        assert sum((weight for rate, weight in atoms), F(0)) == F(4, 9)*response_norm
        geometry_rows.append(dict(block_factor=factor, boundary_links=len(active),
                                  one_link_plaquettes=counts[1], two_link_plaquettes=counts[2],
                                  unique_external_link_cross_term_check=True,
                                  response_norm_squared=response_norm, atoms=atoms,
                                  static_coefficient=atom_moments(atoms, 1)[0][0][0]))
    assert geometry_rows[0]['static_coefficient'] == F(113, 4095)
    return dict(native_shared_path_algebra=rows, actual_cubic_lattice_counts=geometry_rows,
                normalization='Half-trace Wilson loop, free squared L2 norm 1/4',
                coefficient_order=2, uniform_perturbative_remainder_established=False)


def assert_enclosure(row, sigma, M0):
    lower, upper, eps = row['lower'], row['upper'], row['error_factor']
    for value in (lower, a.add(sigma, lower, -1), a.add(upper, sigma, -1),
                  a.add(M0, upper, -1),
                  a.add(a.scale(sigma, eps), a.add(sigma, lower, -1), -1),
                  a.add(a.scale(a.add(M0, sigma, -1), eps), a.add(upper, sigma, -1), -1),
                  a.add(a.scale(M0, eps), a.add(upper, lower, -1), -1)):
        assert a.psd(value)


def leading_controls():
    rows = []
    for factor in (3, 4, 5):
        atoms = leading_atoms(factor)
        M = atom_moments(atoms, 11)
        for s in (F(0), F(1, 7), F(9), F(1000)):
            sigma = atom_response(atoms, s)
            old = None
            for order in range(5):
                row = enclose(M, order, s, min(rate for rate, weight in atoms))
                assert row['rank'] == min(order+1, 4) == row['weighted_rank']
                assert_enclosure(row, sigma, M[0])
                if old:
                    assert a.psd(a.add(row['lower'], old['lower'], -1))
                    assert a.psd(a.add(old['upper'], row['upper'], -1))
                if order >= 3:
                    assert row['lower'] == sigma == row['upper']
                old = row
            rows.append(dict(block_factor=factor, damping=s, exact_response=sigma[0][0],
                             exact_at_order=3, singular_order_four_quotient=True))
    # Beyond the nearest leading pole, the static Taylor series gets worse.
    atoms = leading_atoms(3)
    M = atom_moments(atoms, 7)
    s = F(1000)
    sigma = atom_response(atoms, s)[0][0]
    taylor = [sum(((-s)**n*M[n][0][0] for n in range(order+1)), F(0)) for order in (0, 2, 4)]
    assert abs(taylor[0]-sigma) < abs(taylor[1]-sigma) < abs(taylor[2]-sigma)
    return dict(scope='Exact theta-squared coefficient only; not the full interacting memory',
                four_rate_reconstruction=rows,
                outside_Taylor_disk=dict(damping=s, orders=[0, 2, 4], values=taylor,
                                         errors_strictly_increase=True))


def operator_controls():
    H, C, B, D = previous.fixture()
    M = [previous.moment(B, D, n) for n in range(9)]
    gap = F(2)
    commutator = a.add(a.multiply(M[0], M[1]), a.multiply(M[1], M[0]), -1)
    assert any(v for row in commutator for v in row)
    rows = []
    for s in (F(0), F(1, 7), F(2), F(1000)):
        sigma = previous.self_energy(B, D, s)
        old = None
        for order in range(4):
            row = enclose(M, order, s, gap)
            assert row['rank'] == (2 if order == 0 else 3)
            assert_enclosure(row, sigma, M[0])
            if old:
                assert a.psd(a.add(row['lower'], old['lower'], -1))
                assert a.psd(a.add(old['upper'], row['upper'], -1))
            if order >= 1:
                assert row['lower'] == sigma == row['upper']
            old = row
        rows.append(dict(damping=s, rank_at_orders_zero_and_one=[2, 3],
                         exact_at_order_one=True, includes_exact_vacuum_nullspace=True))
    scaling = []
    for kappa in (F(1, 17), F(7, 3), F(100)):
        scaled = [a.scale(value, kappa**(1-n)) for n, value in enumerate(M)]
        for order in (0, 1, 3):
            for s in (F(0), F(1, 3), F(1000)):
                left = enclose(scaled, order, s, kappa*gap)
                right = enclose(M, order, s/kappa, gap)
                assert left['lower'] == a.scale(right['lower'], kappa)
                assert left['upper'] == a.scale(right['upper'], kappa)
                assert left['error_factor'] == right['error_factor']
        scaling.append(dict(kappa=kappa, comparisons=9))
    nested = []
    for s in (F(1, 7), F(2), F(1000)):
        row = enclose(M, 0, s, gap)
        lo, hi = a.add(C, row['upper'], -1), a.add(C, row['lower'], -1)
        for value in (lo, hi):
            assert a.psd(a.add(value, a.diag((0, 1, 1)), -1))
        direct = previous.effective(H, (0, 1), s)
        lower = previous.effective(lo, (0, 1), s)
        upper = previous.effective(hi, (0, 1), s)
        assert a.psd(a.add(direct, lower, -1)) and a.psd(a.add(upper, direct, -1))
        full_effective = previous.effective(H, (0, 1, 2), s)
        sampled = previous.sub(full_effective, (0, 1))
        assert sampled != direct
        nested.append(dict(damping=s, lower=lower, exact=direct, upper=upper,
                           sampled_matrix_is_not_exact_elimination=True))
    return dict(fixture='YM64 six-dimensional positive matrix, not a Yang--Mills truncation',
                noncommuting_moment_commutator=commutator, full_response_enclosures=rows,
                clock_covariance=scaling, nested_same_parameter_enclosures=nested)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = native.digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = native.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = native.digest(SOURCES)
    return checked


def run():
    native.previous.runtime_check()
    return native.previous.encode(dict(
        certificate_type='YM70_NATIVE_FOUR_RATE_COEFFICIENT_AND_CERTIFIED_SPECTRAL_MEMORY',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        native=native_controls(), leading=leading_controls(), operators=operator_controls(),
        spectral_band_budgets=[dict(ratio=s, order=N, relative_width=error_factor(N, s, 1))
                               for s in (F(1), F(10), F(100), F(1000)) for N in (0, 3, 15, 63)],
        model=dict(actual_YM62_operator_and_YM64_adapted_cut_retained=True,
                   actual_leading_Wilson_memory_computed=True,
                   full_interacting_positive_parameter_enclosure_theorem=True,
                   theorem_requires_exact_moment_inputs=True,
                   unbounded_hidden_generator_allowed=True,
                   same_parameter_nested_order=True,
                   full_interacting_moments_numerically_computed=False,
                   four_rates_for_full_interacting_memory=False,
                   uniform_perturbative_remainder=False,
                   enlarged_bare_coupling_window=False,
                   full_spatial_Schur_memory_locality=False,
                   continuum_R4_QFT=False, physical_mass_prediction=False),
        claim_status='EXACT_LEADING_GAUGE_MEMORY__FULL_RESPONSE_APPROXIMATION_THEOREM__CONTINUUM_OPEN',
        proof_scope='Written proofs with native exact coefficients and separate rational matrix controls',
        open='Interacting moments with spatial/scale-uniform errors and a continuum coupling trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM70 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM70 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], leading_static_coefficient='113/4095',
                          full_positive_parameter_enclosure_theorem=True, continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
