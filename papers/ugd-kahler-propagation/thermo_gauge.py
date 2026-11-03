"""CP-1--CP-8 exact publication controls; no primitive arithmetic engine."""
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from math import isqrt

import model as m


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def lin(coeffs, basis):
    answer = m.zero(len(basis[0]))
    for c, b in zip(coeffs, basis):
        answer = m.add(answer, m.scale(c, b))
    return answer


def outer(u):
    return [[x*y for y in u] for x in u]


def native_delta(x):
    return 2*(m.trace(x)/4)**2-m.trace(m.mm(x, x))/4


def pure(psi, j):
    return m.add(outer(psi), outer(m.mv(j, psi)))


def cone_and_current_checks():
    _, _, a, j = m.fixture()
    ident = m.eye(4)
    basis = [ident]+a
    units = []
    for p, q in product(range(4), repeat=2):
        e = m.zero(4)
        e[p][q] = Q(1)
        units.append(e)
    columns = [sum(m.add(e, m.scale(-1, m.transpose(e))), [])+sum(m.comm(e, j), []) for e in units]
    require(m.rank(m.transpose(columns)) == 12, 'J-Hermitian sector dimension four')
    require(m.rank([sum(x, []) for x in basis]) == 4, 'native Hermitian basis independent')
    cone_count = 0
    for t, x1, x2, x3 in product(range(-2, 3), (-1, 0, 1), (-1, 0, 1), (-1, 0, 1)):
        x = lin([Q(t), Q(x1), Q(x2), Q(x3)], basis)
        delta = Q(t*t-x1*x1-x2*x2-x3*x3)
        require(native_delta(x) == delta and m.det(x) == delta*delta, 'native quadratic determinant')
        positive = all(m.det([row[:n] for row in x[:n]]) > 0 for n in range(1, 5))
        require(positive == (t > 0 and delta > 0), 'Sylvester positivity versus future cone')
        cone_count += 1

    # A non-orthogonal coframe tests the full current identity, including shift.
    v = [[Q(2), Q(1), Q(0)], [Q(0), Q(1), Q(1)], [Q(0), Q(0), Q(3)]]
    beta = [Q(1, 3), Q(-1, 5), Q(2, 7)]
    spatial_inverse = m.inverse(m.mm(m.transpose(v), v))
    matrices = [lin([beta[i]]+[v[k][i] for k in range(3)], basis) for i in range(3)]

    def current(psi):
        return [m.dot(psi, psi)/2]+[-m.dot(psi, m.mv(x, psi))/2 for x in matrices]

    def norm(cur):
        shifted = [cur[i+1]+beta[i]*cur[0] for i in range(3)]
        return -cur[0]*cur[0]+m.dot(shifted, m.mv(spatial_inverse, shifted))

    phase = m.add(m.scale(Q(3, 5), ident), m.scale(Q(4, 5), j))
    pure_count = 0
    for values in product((-1, 0, 1), repeat=4):
        psi = list(map(Q, values))
        n = m.dot(psi, psi)
        bloch = [m.dot(psi, m.mv(x, psi)) for x in a]
        p = pure(psi, j)
        require(p == lin([n/2]+[b/2 for b in bloch], basis), 'real outer product versus Pauli bilinear')
        require(m.dot(bloch, bloch) == n*n and native_delta(p) == 0, 'pure null identity')
        require(m.rank(p) == (2 if n else 0), 'complex-rank-one response')
        require(norm(current(psi)) == 0, 'current null in shifted anisotropic metric')
        require(pure(m.mv(phase, psi), j) == p, 'central phase cannot create relative coherence')
        pure_count += 1

    pair_count = 0
    for s, t, u, cs in product((Q(1), Q(2)), (Q(-1), Q(0), Q(1)), (Q(1), Q(3)),
                                ((Q(1), Q(0)), (Q(3, 5), Q(4, 5)))):
        co, si = cs
        aa, bb, cc = s*s, s*t, t*t+u*u
        psi1, psi2 = [s, Q(0), t*co, -t*si], [Q(0), Q(0), u, Q(0)]
        h = lin([(aa+cc)/2, bb*co, bb*si, (aa-cc)/2], basis)
        require(m.add(pure(psi1, j), pure(psi2, j)) == h, 'explicit two-response thermo factorization')
        cur = [x+y for x, y in zip(current(psi1), current(psi2))]
        require(norm(cur) == -(aa*cc-bb*bb), 'thermo determinant is summed-current invariant')
        pair_count += 1
    relative = m.add(m.scale(Q(3, 5), ident), m.scale(Q(-4, 5), m.mm(j, a[2])))
    moved = m.mm(m.mm(relative, a[0]), m.transpose(relative))
    require(moved == lin([Q(-7, 25), Q(24, 25)], a[:2]), 'relative phase rotation has the stated sign')
    require(moved != a[0], 'relative and central phases are distinct')
    example_fields = [[Q(2), Q(0), Q(1), Q(0)], [Q(0), Q(0), Q(3), Q(0)]]
    example_current = [sum(m.dot(p, p)/2 for p in example_fields)]+[
        -sum(m.dot(p, m.mv(x, p))/2 for p in example_fields) for x in a]
    example_norm = -example_current[0]**2+m.dot(example_current[1:], example_current[1:])
    require(example_current == [7, -2, 0, 3] and example_norm == -36, 'reported rational current example')
    return dict(hermitian_dimension=4, cone_fixtures=cone_count, null_current_fixtures=pure_count,
                thermo_pair_decompositions=pair_count, relative_phase_control=True,
                rational_thermo_example=dict(hessian=[[4, 2], [2, 10]], current=[7, -2, 0, 3], norm=-36))


@dataclass(frozen=True)
class Jet:
    """Value plus four coordinate derivatives, for independent EL evaluation."""
    value: Q
    grad: tuple

    @staticmethod
    def lift(x):
        return x if isinstance(x, Jet) else Jet(Q(x), (Q(0),)*4)

    def __add__(self, other):
        other = self.lift(other)
        return Jet(self.value+other.value, tuple(x+y for x, y in zip(self.grad, other.grad)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.value, tuple(-x for x in self.grad))

    def __sub__(self, other):
        return self+(-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other)+(-self)

    def __mul__(self, other):
        other = self.lift(other)
        return Jet(self.value*other.value,
                   tuple(x*other.value+self.value*y for x, y in zip(self.grad, other.grad)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.lift(other)
        reciprocal = Jet(1/other.value, tuple(-x/(other.value*other.value) for x in other.grad))
        return self*reciprocal

    def __rtruediv__(self, other):
        return self.lift(other)/self


def square_root(jet):
    n, d = isqrt(jet.value.numerator), isqrt(jet.value.denominator)
    root = Q(n, d)
    require(root*root == jet.value, 'square-root fixture is exact rational')
    return Jet(root, tuple(x/(2*root) for x in jet.grad))


def thermo_chart_checks():
    count = 0
    curvature = None
    poly = {(1, 0): Q(10), (0, 1): Q(-1), (2, 0): Q(1),
            (1, 1): Q(1), (0, 2): Q(5, 2), (2, 1): Q(1, 2)}
    require(m.derivative(poly, [0, 0], (0,)) == 10
            and -m.derivative(poly, [0, 0], (1,)) == 1, 'positive reference temperature and pressure')
    for s0, root in product((Q(-1, 3), Q(0), Q(1, 4)), (Q(2), Q(3), Q(4))):
        v0 = ((1+s0)**2+root**2)/5-2
        s = Jet(s0, (Q(1), Q(0), Q(0), Q(0)))
        v = Jet(v0, (Q(0), Q(1), Q(0), Q(0)))
        a, b = 2+v, 1+s
        delta = 5*a-b*b
        sd = square_root(delta)
        h = [[m.derivative(poly, [s0, v0], (i, k)) for k in range(2)] for i in range(2)]
        require(h == [[a.value, b.value], [b.value, Q(5)]], 'thermo potential Hessian')
        require(m.det(h) == root*root and a.value > 0, 'strict stability of chart fixture')
        alpha_s, alpha_v = 1/(2*sd), -b/(2*a*sd)
        f = alpha_v.grad[0]-alpha_s.grad[1]
        require(f == -5/(4*root**3), 'independently differentiated thermo phase curvature')
        p = Q(1, 6)-1/(2*sd)
        require(-p.grad[1] == f and p.grad[1] > 0, 'Darboux curvature and regular Jacobian')
        chi_s, chi_v = 1/sd-Q(1, 6), -b/(2*a*sd)
        require((p+chi_s).value == alpha_s.value and chi_v.value == alpha_v.value,
                'connection equals P dQ plus analytic gauge gradient')
        require(chi_v.grad[0]-chi_s.grad[1] == 0, 'analytic gauge gradient is closed')
        inverted = ((1+s0)**2+1/(4*(p.value-Q(1, 6))**2))/5-2
        require(inverted == v0, 'explicit thermo Darboux inverse')
        if s0 == 0 and v0 == 0:
            curvature = str(f)
        count += 1
    require(curvature == '-5/108', 'same curved thermo origin as PF-1')
    return dict(exact_chart_fixtures=count, origin_curvature=curvature,
                stable_inverse_and_gauge_gradient=True)


def gauge_checks():
    _, _, a, j = m.fixture()
    ident = m.eye(4)
    unitary = m.add(m.scale(Q(3, 5), ident), m.scale(Q(-4, 5), j))
    count = variations = 0
    g, z, rho = Q(2, 3), Q(5, 7), Q(3, 2)
    for r, s in product(range(-2, 3), range(-1, 2)):
        psi = [Q(r+i*s+1, i+1) for i in range(4)]
        dpsi = [[Q(r+(mu+1)*(i+1)-s, mu+i+2) for i in range(4)] for mu in range(4)]
        potential = [Q(r-s+mu, mu+3) for mu in range(4)]
        dchi = [Q(s+2*r-mu, mu+2) for mu in range(4)]
        matrices = [m.add(a[i], m.scale(Q(s+i, 5), ident)) for i in range(3)]
        # M is locally constant, d_t log rho=2/5 and spatial rho derivatives vanish.
        drift = m.scale(Q(1, 5), ident)

        def equation(field, derivatives, aa):
            cov = [[x+g*aa[mu]*y for x, y in zip(derivatives[mu], m.mv(j, field))] for mu in range(4)]
            spatial = [m.mv(matrices[i], cov[i+1]) for i in range(3)]
            return [cov[0][k]-sum(v[k] for v in spatial)+m.mv(drift, field)[k] for k in range(4)]

        moved_psi = m.mv(unitary, psi)
        moved_derivatives = [m.mv(unitary, [x-g*dchi[mu]*y for x, y in zip(dpsi[mu], m.mv(j, psi))])
                             for mu in range(4)]
        moved_a = [x+y for x, y in zip(potential, dchi)]
        residual = equation(psi, dpsi, potential)
        div_current = z*g*rho*(m.dot(psi, dpsi[0])
                              -sum(m.dot(psi, m.mv(matrices[i], dpsi[i+1])) for i in range(3))
                              +m.dot(psi, psi)/5)
        require(div_current == z*g*rho*m.dot(psi, residual), 'weighted off-shell Noether identity')
        moved_residual = equation(moved_psi, moved_derivatives, moved_a)
        require(moved_residual == m.mv(unitary, residual), 'local gauge covariance with derivative phase')

        def density(field, derivatives, aa):
            return z*rho*m.dot(field, m.mv(j, equation(field, derivatives, aa)))/2

        require(density(psi, dpsi, potential) == density(moved_psi, moved_derivatives, moved_a),
                'local gauge action invariance')
        current = [z*g*rho*m.dot(psi, psi)/2]+[-z*g*rho*m.dot(psi, m.mv(x, psi))/2 for x in matrices]
        for mu in range(4):
            plus, minus = list(potential), list(potential)
            plus[mu] += Q(1, 3)
            minus[mu] -= Q(1, 3)
            require((density(psi, dpsi, plus)-density(psi, dpsi, minus))/Q(2, 3) == -current[mu],
                    'independent potential variation derives source current')
            variations += 1
        count += 1
    return dict(local_gauge_jets=count, independent_source_variations=variations,
                off_shell_noether_identities=count)


def sigma(n):
    result = m.zero(2*n)
    for k in range(n):
        result[2*k][2*k+1] = Q(-1)
        result[2*k+1][2*k] = Q(1)
    return result


def rank_checks():
    table = []
    for n in range(1, 5):
        d = m.zero(2*n, 4)
        if n <= 2:
            for a in range(2*n):
                d[a][a] = Q(1)
        elif n == 3:
            for a, mu in [(0, 0), (1, 1), (2, 2), (4, 3)]:
                d[a][mu] = Q(1)
        else:
            for mu in range(4):
                d[2*mu][mu] = Q(1)
        c = m.mm(sigma(n), d)
        f = m.mm(m.transpose(d), c)
        expected = [(2, 2), (4, 4), (4, 2), (4, 0)][n-1]
        require((m.rank(d), m.rank(f)) == expected, 'sharp channel table')
        require(m.rank(c) == m.rank(d), 'variation rank equals state-map rank')
        table.append(dict(channels=n, map_rank=m.rank(d), curvature_rank=m.rank(f)))
    sampled = 0
    for n, seed in product(range(1, 5), range(1, 8)):
        d = [[Q(((a+1)*(mu+seed+1)+a*a+mu*mu) % 7-3, seed)
              for mu in range(4)] for a in range(2*n)]
        c = m.mm(sigma(n), d)
        f = m.mm(m.transpose(d), c)
        r = m.rank(d)
        require(m.rank(c) == r and m.rank(f) >= max(0, 2*r-2*n), 'symplectic restriction rank bound')
        sampled += 1
    return dict(sharp_channel_table=table, additional_rank_fixtures=sampled)


def canonical_connection(phi, dphi):
    n = len(phi)//2
    potential = [sum(phi[2*k+1]*dphi[2*k][mu] for k in range(n)) for mu in range(4)]
    f = [[sum(dphi[2*k+1][mu]*dphi[2*k][nu]-dphi[2*k+1][nu]*dphi[2*k][mu]
              for k in range(n)) for nu in range(4)] for mu in range(4)]
    return potential, f


def lagrangian(phi, dphi, current, e2):
    potential, f = canonical_connection(phi, dphi)
    eta = [-1, 1, 1, 1]
    kinetic = -sum(eta[mu]*eta[nu]*f[mu][nu]*f[mu][nu]
                   for mu, nu in product(range(4), repeat=2))/(4*e2)
    return kinetic-sum(x*y for x, y in zip(current, potential))


def direct_map_el(phi, dphi, current, e2):
    """Differentiate L in its independent jets, then take coordinate divergence.

    L is at most quadratic in each single derivative entry; the centered
    difference is its exact partial derivative, not a numerical approximation.
    """
    el = []
    variations = 0
    for a in range(len(phi)):
        pp, pm = list(phi), list(phi)
        pp[a], pm[a] = pp[a]+1, pm[a]-1
        value = ((lagrangian(pp, dphi, current, e2)-lagrangian(pm, dphi, current, e2))/2).value
        variations += 1
        for mu in range(4):
            dp, dm = [list(row) for row in dphi], [list(row) for row in dphi]
            dp[a][mu], dm[a][mu] = dp[a][mu]+1, dm[a][mu]-1
            momentum = (lagrangian(phi, dp, current, e2)-lagrangian(phi, dm, current, e2))/2
            value -= momentum.grad[mu]
            variations += 1
        el.append(value)
    return el, variations


def pullback_el_checks():
    eta = [-1, 1, 1, 1]
    cases = variations = nonconserved_cases = 0
    for n, seed in product(range(1, 5), (1, 2)):
        phi, dphi = [], []
        for a in range(2*n):
            first = [Q(((a+1)*(mu+2)+seed) % 5-2, 3) for mu in range(4)]
            second = [[Q(((a+seed)*(mu+nu+1)+mu*nu) % 5-2, 7) for nu in range(4)] for mu in range(4)]
            phi.append(Jet(Q(a-seed, 5), tuple(first)))
            dphi.append([Jet(first[mu], tuple(second[mu])) for mu in range(4)])
        current = [Jet(Q(mu+seed, 3), tuple(Q((mu+1)*(nu+1), 11) if seed == 2 else Q(0)
                                          for nu in range(4))) for mu in range(4)]
        e2 = Q(3, 2)
        _, f = canonical_connection(phi, dphi)
        residual = [sum(eta[nu]*eta[mu]*f[nu][mu].grad[nu] for nu in range(4))/e2-current[mu].value
                    for mu in range(4)]
        div_residual = -sum(current[mu].grad[mu] for mu in range(4))
        c = m.mm(sigma(n), [[z.value for z in row] for row in dphi])
        alpha = [phi[a+1].value if a % 2 == 0 else Q(0) for a in range(2*n)]
        expected = [m.dot(row, residual)-aa*div_residual for row, aa in zip(c, alpha)]
        actual, tested = direct_map_el(phi, dphi, current, e2)
        require(actual == expected, 'independent state-map Euler derivative versus Cartan projection')
        if div_residual:
            require(expected != m.mv(c, residual), 'current-conservation premise cannot be dropped')
            nonconserved_cases += 1
        cases += 1
        variations += tested

    # An actual constant matter solution with a nonzero, conserved null current.
    _, _, a, _ = m.fixture()
    psi = [Q(1), Q(0), Q(0), Q(1)]
    cur = [m.dot(psi, psi)/2]+[-m.dot(psi, m.mv(x, psi))/2 for x in a]
    require(cur == [1, 0, 1, 0], 'spurious coupled solution current')
    controls = []
    for n in (2, 4):
        phi = [Jet.lift(0) for _ in range(2*n)]
        d = m.zero(2*n, 4)
        if n == 2:
            d[0][1] = d[2][3] = Q(1)
        else:
            for mu in range(4):
                d[2*mu][mu] = Q(1)
        # Values at the coordinate origin, with fixed affine map derivatives.
        phi = [Jet(Q(0), tuple(row)) for row in d]
        dphi = [[Jet.lift(x) for x in row] for row in d]
        cj = [Jet.lift(x) for x in cur]
        potential, f = canonical_connection(phi, dphi)
        require(all(x.value == 0 for x in potential) and all(z.value == 0 for row in f for z in row),
                'zero field control')
        el, tested = direct_map_el(phi, dphi, cj, Q(1))
        variations += tested
        if n == 2:
            require(el == [0]*4, 'two-channel map accepts a non-Maxwell zero-field solution')
        else:
            require([el[2*mu+1] for mu in range(4)] == [-x for x in cur],
                    'four-channel variation detects the nonzero Maxwell residual')
        controls.append(dict(channels=n, rank=m.rank(d), map_equations_satisfied=not any(el),
                             full_maxwell_equations_satisfied=False))
    return dict(independent_jet_cases=cases, exact_lagrangian_variations=variations,
                nonconserved_source_controls=nonconserved_cases, zero_field_counterexamples=controls)


def run():
    return dict(protocol='THERMO_GAUGE_COMPLETION_CP1_CP8',
                cone=cone_and_current_checks(), thermo_chart=thermo_chart_checks(),
                gauge=gauge_checks(), ranks=rank_checks(), euler_lagrange=pullback_el_checks(),
                boundaries=dict(continuum_and_metric='SUPPLIED', gauge_kinetic_action='DECLARED',
                                current_mass_interpretation='NOT_CLAIMED',
                                global_or_coupled_well_posedness='NOT_PROVED',
                                physical_constants='NOT_PREDICTED'))


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
