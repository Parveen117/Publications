"""YM56: exact Hessian, symmetry and compact-gap controls; Python 3.12 only.

Uses the existing native engines. Default/--check is read-only.
--write updates only this chapter's result and certificate pin.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym54_response_protocol as native
import ym52_energy_observer as energy
import ym55_anisotropic_joint_limit as joint

a, p, q, y = native.a, native.p, native.q, native.y
exact = native.exact
RESULT = HERE / 'YM56_RESULT.json'
PIN = HERE / 'EXPECTED_YM56.sha256'
SOURCES = HERE / 'YM56_SOURCE_PINS.json'
SAMPLES = tuple(map(F, (-2, -1, F(-1,2), 0, F(1,10), F(1,2), 1, 2, 4)))


def runtime_check():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError('YM56 requires Python 3.12')
    if not __debug__:
        raise RuntimeError('Optimized Python is refused')


def diff(poly, axis):
    out = {}
    for exponent, coefficient in poly.items():
        if exponent[axis]:
            new = list(exponent)
            new[axis] -= 1
            out[tuple(new)] = coefficient * exponent[axis]
    return out


def evaluate(poly, point=(F(0), F(0))):
    return sum((c * point[0]**m[0] * point[1]**m[1]
                for m, c in poly.items()), F(0))


def potential(lam):
    lam = exact(lam)
    return {(1,0): F(10), (0,1): F(-10), (2,0): F(1,2),
            (0,2): F(1,2), (3,0): (lam+2)/6, (1,2): lam/2}


def hessian(lam, point=(F(0), F(0))):
    u = potential(lam)
    return [[evaluate(diff(diff(u, i), j), point) for j in range(2)]
            for i in range(2)]


def family(lam):
    lam = exact(lam)
    u = potential(lam)
    derivatives = [[[evaluate(diff(diff(diff(u, i), j), k))
                     for j in range(2)] for i in range(2)] for k in range(2)]
    return native.response_data(hessian(lam), derivatives, native.I)


def rates(lam):
    lam = exact(lam)
    even, odd = min(F(1), lam*lam)/2, (1+lam*lam)/8
    return min(even, odd), even, odd


def witness(axis):
    if isinstance(axis, bool) or axis not in (1, 2, 3):
        raise ValueError('Internal axis must be 1, 2 or 3')
    return p.add(*(p.mul(p.COORD[k], p.COORD[k]) for k in (1,2,3) if k != axis),
                 p.scale(y.RADIUS, F(-1,2)))


def involution(poly):
    # Pullback of quaternion conjugation by e_1.
    return {m: c*(-1)**(m[2]+m[3]) for m, c in poly.items()}


def profile_bounds(ell, bound):
    ell, bound = exact(ell), exact(bound)
    if not 0 < ell <= bound:
        raise ValueError('Require 0 < ell <= B')
    beta = min(F(1), ell*ell)/2
    trace = (1+bound*bound)/2
    return dict(beta=beta, trace_max=trace, theta_ceiling=beta/2400,
                curvature_only_beta=ell*ell/(2*(1+bound*bound)))


def hessian_controls():
    rows = []
    for lam in SAMPLES:
        d = family(lam)
        assert d['X'] == [a.diag((lam+2, lam)), a.scale(native.L, lam)]
        assert d['Y'] == [native.K, a.scale(native.L, lam)]
        assert d['sigma'] == [1+lam, 0]
        assert d['F'] == a.scale(native.R, lam/2)
        assert d['marker'] == lam/2 and d['tau'] == 2*(1+lam*lam)
        assert native.det2(d['G']) == 16*d['marker']**2
        # At H=I the fourth derivatives vanish. Differentiate the two
        # Christoffel matrices independently, then include their product.
        dx, dy = d['derivatives']
        derivative_difference = a.scale(native.comm(dx, dy), F(-1,2))
        curvature = a.add(derivative_difference,
                          native.comm(a.scale(dx,F(1,2)), a.scale(dy,F(1,2))))
        assert curvature == d['F']
        rows.append(dict(lam=lam, marker=d['marker'], tau=d['tau'],
                         scalar_x=d['sigma'][0], full_rate=rates(lam)[0]))
    corners = 0
    for bound in (F(1,2), F(1), F(2), F(4)):
        r = 1/(4*(bound+1))
        for lam, x, z in itertools.product((-bound, 0, bound), (-r,r), (-r,r)):
            H = hessian(lam, (x,z))
            assert a.psd(a.add(H, a.scale(native.I, F(-1,2))))
            corners += 1
    return dict(reference_jets=len(rows), rows=rows, positivity_corners=corners,
                all_parameter_positivity='T1 affine norm bound',
                naked_K_lambdaL_hessian_compatibility_requires_lambda_minus_one=True)


def symmetry_controls():
    actions = brackets = moments = 0
    for lam in SAMPLES:
        d = family(lam)
        C = d['C']
        assert C == a.diag((F(1,2), lam*lam/2, 0))
        assert family(-lam)['C'] == C
        assert family(-lam)['marker'] == -d['marker']
        records = [(1,0,0), (-1,0,0), (0,lam,0), (0,-lam,0)]
        reflected = [(v[0],-v[1],-v[2]) for v in records]
        assert sorted(records) == sorted(reflected)
        direct = [[sum((F(v[i])*F(v[j])/4 for v in records), F(0))
                   for j in range(3)] for i in range(3)]
        assert direct == C
        moments += 1
        for degree in range(1,4):
            for m in y.basis(degree):
                f = {m:F(1)}
                assert involution(involution(f)) == f
                assert involution(q.lap(f,C)) == q.lap(involution(f),C)
                for axis, sign in ((1,1),(2,-1),(3,-1)):
                    assert involution(y.native_generator(involution(f),axis)) == p.scale(
                        y.native_generator(f,axis),sign)
                left = p.add(y.native_generator(y.native_generator(f,2),1),
                             p.scale(y.native_generator(y.native_generator(f,1),2),-1))
                assert p.scale(left,lam) == p.scale(y.native_generator(f,3),lam)
                actions += 1
                brackets += 1
    return dict(polynomial_symmetry_checks=actions, bracket_checks=brackets,
                literal_four_label_moments=moments, degree_range=[1,3],
                sign_of_curvature_lost_by_heat=True)


def gap_controls():
    blocks = sharp = 0
    for lam in SAMPLES:
        C = family(lam)['C']
        full, even, odd = rates(lam)
        assert (full,even,odd) == energy.rates(tuple(sorted((F(1,2),lam*lam/2,F(0)))))
        for axis, rate in ((1,lam*lam/2),(2,F(1,2))):
            w = witness(axis)
            assert q.lap(w,C) == p.scale(w,rate)
            assert y.phi(w) == 0 and energy.variance(w) == F(1,12)
            sharp += 1
        assert q.lap(p.COORD[0],C) == p.scale(p.COORD[0],odd)
        for n in range(1,13):
            squares, G = energy.spin_squares(n)
            op = a.scale(a.add(squares[0],a.scale(squares[1],lam*lam)),F(1,2))
            floor = odd if n%2 else even
            residual = a.multiply(G,a.add(op,a.scale(a.identity(n+1),-floor)))
            assert a.psd(residual)
            blocks += 1
        tau = 2*(1+lam*lam)
        chi2 = 4*lam*lam/(1+lam*lam)**2
        normalized = 16*full/tau
        if chi2 >= F(3,4):
            assert normalized == 1
        else:
            assert normalized < 1 and (1-normalized/2)**2 == 1-chi2
        if lam:
            assert rates(1/lam)[0] == full/(lam*lam)
    assert not q.lap(witness(1),family(0)['C'])
    return dict(exact_spin_inequalities=blocks, tensor_degrees=[1,12],
                sharp_quadratic_modes=sharp, all_content_proof='YM52-T4 plus YM56-T3',
                normalized_curvature_law_checks=len(SAMPLES),
                zero_parameter_stationary_centered_witness=True)


def general_shape_controls():
    count = 0
    # A common 3/5,4/5 rotation gives non-axis-aligned actual response
    # words with known, independent eigenvalues r^2/2,s^2/2.
    for r,s in ((F(1),F(2)),(F(2),F(2)),(F(1),F(10)),(F(3),F(0))):
        d = native.from_shapes(((3*r/5,4*r/5),(-4*s/5,3*s/5)))
        tau, f = d['tau'], d['marker']
        b,c = sorted((r*r/2,s*s/2))
        assert tau == 4*(b+c) and 16*f*f == 16*b*c
        chi2 = 64*f*f/(tau*tau)
        root = (c-b)/(c+b)
        assert root*root == 1-chi2
        actual = energy.rates((F(0),b,c))[0]
        assert 16*actual/tau == min(F(1),2*(1-root))
        count += 1
    return dict(rotated_shape_families=count)


def chain_controls():
    d = profile_bounds(F(1,2),F(2))
    theta = F(1,40000)
    assert d['beta'] == F(1,8) and d['trace_max'] == F(5,2)
    assert d['curvature_only_beta'] == F(1,40)
    assert d['curvature_only_beta']/2400 < theta < d['theta_ceiling']
    matrices = [family(lam)['C'] for lam in (F(-2),F(-1,2),F(1),F(2))]
    rot = q.rotation((F(1,2),)*4)
    matrices += [q.rotate_tensor(C,rot) for C in matrices]
    joint.trace_profile(matrices,d['beta'],d['trace_max'])
    J,R,rho,u = joint.parameters(d['beta'],theta)
    c = joint.constants(d['beta'],theta,d['trace_max'],J,R,rho,u,F(16))
    assert c['gamma'].lo > 0 and c['margin'] > 0
    slow = []
    for n in (2,10,100):
        lam = F(1,n)
        C = family(lam)['C']
        quotient = energy.energy(witness(1),C)/energy.variance(witness(1))
        assert quotient == rates(lam)[0] == F(1,2*n*n)
        slow.append(dict(n=n, quotient=quotient))
    return dict(**d, theta=theta, profile_tensor_checks=len(matrices),
                J=J, R=R, rho=rho, margin=c['margin'],
                chain_gamma_lower=c['gamma'].lo, chain_gamma_upper=c['gamma'].hi,
                fixed_profile_required=True, nonuniform_free_profile=slow,
                physical_mass_identified=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_checks():
    pins = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in pins['upstream_sha256'].items():
        actual = digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in pins['local_inputs']:
        checked[path] = digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = digest(SOURCES)
    return checked


def encode(value):
    if isinstance(value,F):
        return str(value)
    if isinstance(value,dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [encode(v) for v in value]
    return value


def run():
    runtime_check()
    return encode(dict(certificate_type='YM56_HESSIAN_LAMBDA_SYMMETRY_AND_GAP',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        hessian=hessian_controls(), symmetry=symmetry_controls(), gap=gap_controls(),
        general_shape=general_shape_controls(), chain=chain_controls(),
        claim_status='EXACT_SELECTED_FREE_PROTOCOL_AND_INHERITED_FIXED_PROFILE_CHAIN__4D_CLAY_OPEN',
        evidence_scope=dict(general_proof='written, not mechanically formalized or externally reviewed',
                            finite_controls='exact rational native-engine replay',
                            open='physical selection, clock, actual row closure, spatial continuum, 4D Yang-Mills')))


def check(cert, result=RESULT, pin=PIN):
    sha = native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM56 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM56 PASS',sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],
                          hessian_jets=cert['hessian']['reference_jets'],
                          spin_inequalities=cert['gap']['exact_spin_inequalities'],
                          chain_gamma_lower_display=float(F(cert['chain']['chain_gamma_lower']))),sort_keys=True))


if __name__ == '__main__':
    main()
