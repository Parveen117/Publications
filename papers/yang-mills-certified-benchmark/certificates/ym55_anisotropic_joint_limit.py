"""YM55: anisotropic joint local-history controls. Python 3.12 only.

Written proofs establish the infinite result. Finite fixtures test comparison
operations and failure boundaries; they do not discretize SU(2) or prove Clay.
Default/--check never writes. No duplicate native operator or interval engine.
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
import ym47_joint_history as joint

old = native.chain
a = native.a
Iv = old.Iv
exact = native.exact
integer = old.integer
ex = old.refine.ex
sqrt_iv = joint.iv_sqrt
decimals = old.old.decimals
RESULT = HERE/'YM55_RESULT.json'
PIN = HERE/'EXPECTED_YM55.sha256'
SOURCES = HERE/'YM55_SOURCE_PINS.json'


def gate(beta, theta, J, R, rho):
    beta, theta, R, rho = map(exact, (beta, theta, R, rho))
    integer(J, 1)
    if beta <= 0 or min(R, rho) <= 1:
        raise ValueError('Positive floor and R,rho>1 required')
    lam = abs(theta)/beta
    zeta = R*(F(2)/(old.old.EPS*J)+(160+80*rho)*lam*J)
    if zeta >= 1:
        return None
    inherited = old.block_gate(beta, theta, J, R)
    assert inherited is not None
    g = inherited['gamma']/Iv(beta)
    return dict(beta=beta, theta=theta, lam=lam, J=J, R=R, rho=rho,
                zeta=zeta, margin=1-zeta, g=g, gamma=inherited['gamma'])


def parameters(beta, theta):
    beta, theta = exact(beta), exact(theta)
    if beta <= 0 or abs(theta)/beta >= F(1, 2400):
        raise ValueError('Strict sufficient window |theta|/beta<1/2400 required')
    lam, J = abs(theta)/beta, 5
    alpha = F(1, 2)+1200*lam
    rho = 1+(1-alpha)/(160*lam*J) if lam else F(2)
    spatial = F(2)/(old.old.EPS*J)+(160+80*rho)*lam*J
    R = (1+1/spatial)/2
    u = 2*rho/(rho+1)
    assert gate(beta, theta, J, R, rho) is not None and 1 < u and u*u < rho
    return J, R, rho, u


def constants(beta, theta, trace_max, J, R, rho, u, cstar):
    c = gate(beta, theta, J, R, rho)
    if c is None:
        raise ValueError('Spatial comparison gate failed')
    trace_max, u, cstar = map(exact, (trace_max, u, cstar))
    ratio = trace_max/c['beta']
    if ratio < 2 or not 1 < u or u*u >= c['rho'] or cstar <= 0 or cstar*cstar < 12*ratio:
        raise ValueError('Require M>=2 beta, 1<u, u^2<rho and cstar^2>=12 M/beta')
    g = c['g'].lo
    if g <= 0:
        raise ValueError('Interval failed to separate the gap from zero')
    c.update(trace_max=trace_max, trace_ratio=ratio, u=u, cstar=cstar,
             glo=g, d=g/(1+g),
             side=40*c['lam']*J*c['rho']*c['R']*(1+2/g)/c['margin'])
    return c


def trace_profile(tensors, beta, trace_max):
    beta, trace_max = exact(beta), exact(trace_max)
    if beta <= 0 or trace_max < 2*beta or not tensors:
        raise ValueError('Nonempty profile and compatible positive bounds required')
    matrices = [native.q.tensor(C) for C in tensors]
    for C in matrices:
        # C-beta I has at most one negative eigenvalue iff the middle is >=beta.
        inertia = a.y37.inertia(a.add(C, a.scale(a.identity(3), -beta)))
        if inertia[1] > 1 or native.trace(C) > trace_max:
            raise ValueError('Actual tensor violates the declared profile bounds')
    return matrices


def step_budget(step, c):
    step = exact(step)
    if not 0 < step <= 1:
        raise ValueError('Dimensionless step must lie in (0,1]')
    r = old.skeleton(1, step)
    ell = c['J']*r
    q = ex(-c['glo']*step).hi
    if not 0 < q < 1:
        raise ValueError('Rounded temporal ratio is not separated from one')
    D = F(r)/old.old.EPS+8*step*c['lam']*ell*ell
    Htime = (1+q)/(1-q)
    Hspace = (c['rho']+1)/(c['rho']-1)
    return dict(r=r, ell=ell, q=q, D=D,
                side=4*step*c['lam']*ell*c['rho']*c['R']*Htime/c['margin'],
                time=c['R']*D*Hspace/c['margin'],
                half_time=c['R']*(D+1)*Hspace/c['margin'])


def normalized_budget(width, step, duration, c):
    integer(width, 1)
    step, duration = exact(step), exact(duration)
    if not 0 < step <= 1 or duration < 0:
        raise ValueError('Positive step <=1 and nonnegative duration required')
    w = c['lam']*(width-1)
    return c['cstar']*w*duration*sqrt_iv(Iv(step)).hi*ex(2*w*step).hi


def history_inputs(T, count, delta):
    T, delta = exact(T), exact(delta)
    integer(count)
    if T < 0 or delta <= 0 or (count and T < count*delta):
        raise ValueError('Invalid separated observation schedule')
    return T, delta


def history_budget(width, step, T, count, delta, c):
    T, delta = history_inputs(T, count, delta)
    step = exact(step)
    normalized_budget(width, step, T, c)
    if step > delta/4:
        raise ValueError('Step exceeds the separated-time gate')
    w = c['lam']*(width-1)
    B = c['cstar']*(8/c['d']+T+2)
    return B*w*sqrt_iv(Iv(step)).hi*ex(2*w*step).hi+4*count*step/delta


def tail(k, support, T, count, delta, N, c):
    integer(k, 1)
    integer(support, 1)
    integer(N)
    T, delta = history_inputs(T, count, delta)
    q, v = 1/c['u'], c['u']**2/c['rho']
    if joint.power_iv(q*q, N).hi > delta/4:
        raise ValueError('Tail begins before separated-time mesh gate')
    B = c['cstar']*(8/c['d']+T+2)
    Z = ex(2*c['lam']*(k+2+2/(c['u']**2-1))).hi
    H = joint.power_iv(q, N).hi*((k+2+2*N)/(1-q)+2*q/(1-q)**2)
    G2 = joint.power_iv(q*q, N).hi/(1-q*q)
    Gv = joint.power_iv(v, N).hi/(1-v)
    return 2*B*c['lam']*Z*H+8*count/delta*G2+4*support*c['side']*Gv


def tail_radius(k, support, T, count, delta, tolerance, c):
    tolerance = exact(tolerance)
    if tolerance <= 0:
        raise ValueError('Positive target tolerance required')
    history_inputs(T, count, delta)
    integer(k, 1)
    integer(support, 1)
    N = 1
    while joint.power_iv(1/c['u']**2, N).hi > exact(delta)/4:
        N *= 2
    while 3*tail(k, support, T, count, delta, N, c) > tolerance:
        N *= 2
    return N, 3*tail(k, support, T, count, delta, N, c)


def crop_radius(step, rho):
    step, rho = exact(step), exact(rho)
    if not 0 < step < 1 or rho <= 1:
        raise ValueError('Require 0<a<1 and rho>1')
    return joint.crop_radius(step, rho)


def joint_budget(step, left, right, k, support, T, count, delta, c):
    integer(left)
    integer(right)
    integer(k, 1)
    integer(support, 1)
    D = crop_radius(step, c['rho'])
    N = min(left, right, D)
    crop = 8*support*c['side']*exact(step)
    local = history_budget(k+2*D+1, step, T, count, delta, c)
    # +1 makes w=lambda*(k+2D), the stated conservative w_*.
    volume = 3*tail(k, support, T, count, delta, N, c)
    return dict(crop_padding=D, padding_used=N, cropped_width=k+2*D,
                crop_error=crop, refinement_error=local,
                volume_error=volume, total=crop+local+volume)


def rank_two_cell():
    return constants(F(3, 2), F(1, 4096), 3, 8, F(4, 3), F(3, 2), F(9, 8), 5)


def thermo_cell():
    return constants(F(5, 234), F(1, 262144), F(65, 648),
                     6, F(5, 4), F(3, 2), F(9, 8), 8)


def parameter_controls():
    rows, comparisons = [], 0
    for name, c in [('rank_two', rank_two_cell()), ('thermo_source', thermo_cell())]:
        for step in [F(1), F(1, 8), F(7, 256), F(1, 1024)]:
            b = step_budget(step, c)
            actual = c['R']*(2*F(b['r'])/(old.old.EPS*b['ell'])+
                     (16+8*c['rho'])*step*c['lam']*b['ell'])
            assert actual <= c['zeta'] < 1
            assert b['side'] <= c['side']/step
            assert 9 <= step*b['r'] < 10
            comparisons += 1
        N, error = tail_radius(3, 6, 2, 2, 1, F(1, 10**6), c)
        rows.append(dict(name=name, beta=str(c['beta']), theta=str(c['theta']),
                         trace_ratio=str(c['trace_ratio']), margin=str(c['margin']),
                         gamma_lower=decimals(c['gamma'].lo, False),
                         gamma_upper=decimals(c['gamma'].hi),
                         radius_for_1e_minus_6=N, remaining_error_upper=decimals(error)))
    generic = []
    for lam in [F(0), F(1, 100000), F(1, 6144), F(1, 3000), F(999, 2400000)]:
        J, R, rho, u = parameters(1, lam)
        c = constants(1, lam, 2, J, R, rho, u, 5)
        assert c['margin'] > 0 and c['glo'] > 0
        generic.append(dict(lam=str(lam), J=J, R=str(R), rho=str(rho)))
    assert rank_two_cell()['margin'] == F(7, 72)
    assert thermo_cell()['margin'] == F(10249, 98304)
    return dict(cells=rows, fine_step_comparisons=comparisons, entire_open_window_controls=generic)


def incidence_controls():
    layouts = interior = boundary = 0
    tiny, rho, r = F(1, 10**6), F(3, 2), 3
    for m, N, ell in itertools.product([1, 3], [1, 2, 5], [2, 5, 8]):
        labels = old.old.blocks(m, N, ell)
        qt = 1-F(1, 20*ell)
        R = qt**(-ell)
        D = F(r)/old.old.EPS+8*tiny*ell*ell
        targets = [(m//2, 1), (0, N)]
        def weight(i, u):
            return sum((rho**(-abs(i-j))*qt**abs(u-v) for j, v in targets), F(0))
        outgoing = R*(2*F(r)/old.old.EPS+(16+8*rho)*tiny*ell*ell)
        for i in range(m):
            for u in range(1, N+1):
                removed, cost = 0, F(0)
                for j, block in labels:
                    top = max(weight(j, v) for v in block)
                    if i == j and u in block:
                        removed += 1
                    if i == j and u in (block[0]-1, block[-1]+1):
                        cost += top*D
                    if abs(i-j) == 1 and u in block:
                        cost += top*4*tiny*ell
                assert removed == ell and cost <= weight(i, u)*outgoing
                interior += 1
        for i in [-1, m]:
            for u in range(1, N+1):
                cost = sum((max(weight(j, v) for v in block)*4*tiny*ell
                            for j, block in labels if abs(i-j) == 1 and u in block), F(0))
                assert cost <= weight(i, u)*4*tiny*ell*ell*rho*R
                boundary += 1
        for i in range(m):
            for u in [0, N+1]:
                cost = sum((max(weight(j, v) for v in block)*D for j, block in labels
                            if j == i and u in (block[0]-1, block[-1]+1)), F(0))
                assert cost <= weight(i, u)*ell*R*D
                boundary += 1
        layouts += 1
    return dict(layouts=layouts, interior_checks=interior, boundary_checks=boundary,
                short_strips_and_duplicate_labels=True)


def zero_midpoint_controls():
    H = old.path_kernel()
    K = a.multiply(H, H)
    D = a.diag((1, F(3, 2), F(5, 4), 2))
    W = a.multiply(D, D)
    S = a.multiply(a.multiply(H, W), H)
    T = a.multiply(a.multiply(D, K), D)
    assert S != T and a.psd(S) and a.psd(T)
    valid = null = identities = 0
    for x, z in itertools.product(range(4), repeat=2):
        weights = [H[x][y]*H[y][z] for y in range(4)]
        if K[x][z]:
            assert sum(w/K[x][z] for w in weights) == 1
            valid += 1
        else:
            assert weights == [0]*4
            null += 1
    for obs in [(F(1),)*4, (F(-1), F(0), F(1), F(2))]:
        for filler in [F(-7), F(13)]:
            for x, z in itertools.product(range(4), repeat=2):
                direct = sum((S[x][y]*obs[y]*S[y][z] for y in range(4)), F(0))
                lifted = F(0)
                for u, v in itertools.product(range(4), repeat=2):
                    numerator = sum((H[u][y]*obs[y]*H[y][v] for y in range(4)), F(0))
                    readout = numerator/K[u][v] if K[u][v] else filler
                    lifted += H[x][u]*W[u][u]*K[u][v]*readout*W[v][v]*H[v][z]
                assert direct == lifted
                identities += 1
    grouped = old.zero_bridge_controls()
    assert grouped['zero_normalizer_rejected']
    return dict(admissible_midpoint_pairs=valid, null_pairs=null,
                actual_ordered_path_identities=identities, arbitrary_null_fillers_immaterial=True,
                old_grouped_endpoint_control=grouped,
                fixture_scope='Four-state rational control; not a compact-carrier discretization')


def profile_controls():
    C = a.diag((0, F(3, 2), F(3, 2)))
    rotations = [native.q.rotation(q) for q in
                 [(F(1), 0, 0, 0), (F(3, 5), 0, F(4, 5), 0), (F(1, 2),)*4]]
    profile = [native.q.rotate_tensor(C, Q) for Q in rotations]
    trace_profile(profile, F(3, 2), 3)
    assert a.multiply(profile[0], profile[1]) != a.multiply(profile[1], profile[0])
    checks = 0
    for matrix in profile:
        for degree in range(1, 4):
            G = native.y.degree_data(degree)['G']
            L = native.q.operator(degree, matrix)
            assert a.psd(a.multiply(G, L))
            checks += 1
    fixture = native.fixture()
    trace_profile([fixture['C']], F(5, 234), F(65, 648))
    assert fixture['beta_floor'] == F(5, 234) and fixture['tau']/4 == F(65, 648)
    beta = 4*abs(fixture['marker'])**2/fixture['tau']
    assert beta == thermo_cell()['beta']
    assert fixture['tau']/(4*beta) == F(169, 36)
    # Two uniformly admissible but cutoff-alternating free local profiles.
    q0 = native.p.COORD[0]
    assert native.y.phi(native.p.mul(q0, q0)) == F(1, 4)
    for C, rate in [(a.diag((0, 1, 1)), F(1, 2)), (a.diag((0, 1, 2)), F(3, 4))]:
        trace_profile([C], 1, 3)
        assert native.q.lap(q0, C) == native.p.scale(q0, rate)
    correlations = [ex(-rate)/Iv(4) for rate in [F(1, 2), F(3, 4)]]
    assert correlations[0].separated_from(correlations[1])
    return dict(fixed_periodic_profile_sites=3, common_diagonalization_not_used=True,
                native_energy_controls=checks, original_thermo_source_reused=True,
                alternating_profile_failure=[dict(lo=decimals(v.lo, False), hi=decimals(v.hi))
                                             for v in correlations])


def normalized_controls():
    rows = []
    observations = ((F(1), F(-1)), (F(1, 2), F(1)), (F(-1, 3), F(1)))
    for speed, theta in [(F(2), F(1, 32)), (F(5), F(-1, 32))]:
        # Exact two-state exp(t(-speed*L+theta*V)), independent of native profiles.
        U = joint.y44.reference2(speed, theta/speed)
        direct = joint.y44.taylor_reference2(speed, theta/speed, N=90)
        assert all(not U[i][j].separated_from(direct[i][j])
                   for i, j in itertools.product(range(2), repeat=2))
        nu, projector = joint.spectral_pair(U)
        V = joint.normalize2(U, nu)
        target = joint.history2(projector, [V, V], observations)
        c = dict(lam=abs(theta), cstar=F(10), d=F(1, 4))
        for n in [4, 16, 64]:
            step = F(1, n)
            S = joint.y44.split2(speed*step, theta/speed)
            lam, projection = joint.spectral_pair(S)
            An = joint.y44.ivpow2(joint.normalize2(S, lam), n)
            difference = [[An[i][j]-V[i][j] for j in range(2)] for i in range(2)]
            error = joint.y44.row_bound2(difference)
            bound = normalized_budget(2, step, 1, c)
            assert error < bound
            value = joint.history2(projection, [An, An], observations)
            err = max(abs((value-target).lo), abs((value-target).hi))
            hbound = history_budget(2, step, 2, 2, 1, c)
            assert err < hbound
            rows.append(dict(speed=str(speed), steps=n, operator_error_upper=decimals(error),
                             operator_budget_upper=decimals(bound),
                             history_error_upper=decimals(err), history_budget_upper=decimals(hbound)))
    return dict(independent_scaled_two_state_histories=rows,
                control_contract='Finite normalized-refinement identities only; no joint-chain gate claimed for these fixture couplings',
                controls_are_not_the_native_infinite_proof=True)


def geometric_controls():
    checks = 0
    for q, N, length in itertools.product([F(1, 2), F(8, 9), F(27, 32)], [0, 3, 12], [1, 7, 20]):
        k = 3
        H = lambda n: q**n*((k+2+2*n)/(1-q)+2*q/(1-q)**2)
        G = lambda n: q**n/(1-q)
        assert H(N)-H(N+length) == sum((k+2+2*n)*q**n for n in range(N, N+length))
        assert G(N)-G(N+length) == sum(q**n for n in range(N, N+length))
        checks += 2
    paths = 0
    for left, right in itertools.product(range(1, 13), repeat=2):
        n = min(left, right)
        for distance in range(n, max(left, right)):
            width = 3+n+distance+1
            assert width <= 3+2*distance+2
        paths += 1
    return dict(exact_tail_difference_identities=checks, asymmetric_extension_paths=paths)


def cutoff_controls():
    c = rank_two_cell()
    rows = []
    for power in [32, 64, 128]:
        step = F(1, 2**power)
        D = crop_radius(step, c['rho'])
        b = joint_budget(step, D, D*D+3, 3, 6, 2, 2, 1, c)
        assert c['rho']**(-D) <= step*step
        assert D == 1 or c['rho']**(1-D) > step*step
        rows.append(dict(step=str(step), crop_padding=D, asymmetric_outer_padding=D*D+3,
                         total_error_upper=decimals(b['total'])))
    assert F(rows[-1]['total_error_upper']) < F(1, 10**6)
    # A slow spatial exhaustion has a separate tail even at the same tiny step.
    slow = joint_budget(F(1, 2**128), 16, 10**9, 3, 6, 2, 2, 1, c)
    assert slow['padding_used'] == 16 and slow['volume_error'] > 1
    return dict(joint_cutoff_cells=rows, slow_exhaustion_not_silently_dropped=True,
                slow_padding_used=16)


def scale_controls():
    c = rank_two_cell()
    checks = 0
    for factor in [F(1, 10), F(2), F(7, 3)]:
        d = constants(c['beta']*factor, c['theta']*factor, c['trace_max']*factor,
                      c['J'], c['R'], c['rho'], c['u'], c['cstar'])
        for key in ['lam', 'trace_ratio', 'margin', 'glo', 'd', 'side']:
            assert d[key] == c[key]
        assert d['gamma'].lo == c['gamma'].lo*factor
        h = F(1, 32)/c['beta']
        assert (h/factor)*d['beta'] == h*c['beta']
        checks += 1
    return dict(clock_rescaling_checks=checks, physical_clock_selected=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_checks():
    pins = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in pins['upstream_sha256'].items():
        actual = digest(ROOT/path)
        if actual != expected:
            raise ValueError('Upstream source changed: '+path)
        checked[path] = actual
    for path in pins['local_inputs']:
        checked[path] = digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = digest(SOURCES)
    return checked


def run():
    if not __debug__ or sys.version_info[:2] != (3, 12):
        raise RuntimeError('Unoptimized Python 3.12 only')
    return dict(certificate_type='YM55_NATIVE_ANISOTROPIC_JOINT_HISTORY_LIMIT',
                verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
                parameter_gates=parameter_controls(), weighted_incidence=incidence_controls(),
                singular_midpoints=zero_midpoint_controls(), native_profiles=profile_controls(),
                normalized_refinement=normalized_controls(), geometric_tails=geometric_controls(),
                joint_cutoffs=cutoff_controls(), scaling=scale_controls(),
                claim_status='FIXED_PROFILE_JOINT_LOCAL_VACUUM_HISTORY_LIMIT_WITH_CORRELATION_GAP',
                evidence_scope=dict(written_results=7, general_proof='Written, not proof-assistant or expert certified',
                    carrier='Native compact chain; fixed site profile; lower middle eigenvalue and upper trace',
                    limit='Both spatial ends to infinity and heat step to zero at unrestricted relative rates',
                    retained='Positivity, heat-time reflection positivity, and centered time-correlation decay',
                    not_asserted='Continuous-time row-Markov closure; physical state/action/clock selection; '
                                 'moving interacting vacua; spatial-spacing limit; 4D Yang--Mills or Clay'))


def check(cert, result=RESULT, pin=PIN):
    expected = old.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != expected:
        raise ValueError('YM55 fresh certificate/pin mismatch')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    expected = old.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(expected+'\n')
    else:
        check(cert)
    print('YM55 PASS', expected)
    print(json.dumps(dict(written_results=7, cutoff_cells=cert['joint_cutoffs']['joint_cutoff_cells'],
                          physical_clock='OPEN', four_dimensional_ym='OPEN'), sort_keys=True))


if __name__ == '__main__':
    main()
