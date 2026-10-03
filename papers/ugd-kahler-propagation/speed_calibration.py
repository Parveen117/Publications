"""SC exact classical characteristic controls; no empirical parameter fitting."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
import json
from pathlib import Path

import model as m
import spin_matter as sm

IX = range(4)
SYMMETRIC = list(combinations_with_replacement(IX, 2))
ETA = [[Q((-1, 1, 1, 1)[i] if i == j else 0) for j in IX] for i in IX]
_, GAMMA, A, J = m.fixture()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def outer(u, v):
    return [[x*y for y in v] for x in u]


def symmetric_vector(h):
    return [h[i][j] for i, j in SYMMETRIC]


def symmetric_basis():
    result = []
    for i, j in SYMMETRIC:
        h = m.zero(4)
        h[i][j] = h[j][i] = Q(1)
        result.append(h)
    return result


def quadratic(gi, p):
    return m.dot(p, m.mv(gi, p))


def maxwell_symbol(gi, p):
    return m.add(m.scale(quadratic(gi, p), m.eye(4)),
                 m.scale(-1, outer(p, m.mv(gi, p))))


def maxwell_from_curvature(gi, p, a):
    """Contract the independently formed exterior derivative F=p wedge a."""
    f = m.add(outer(p, a), m.scale(-1, outer(a, p)))
    pu = m.mv(gi, p)
    return [sum(pu[nu]*f[nu][mu] for nu in IX) for mu in IX]


def einstein_symbol(g, gi, p, h):
    pu = m.mv(gi, p)
    q = m.dot(p, pu)
    s = m.mv(h, pu)
    tr = sum(gi[i][j]*h[i][j] for i in IX for j in IX)
    pp = m.dot(pu, s)
    return [[(p[i]*s[j]+p[j]*s[i]-q*h[i][j]-p[i]*p[j]*tr
              -g[i][j]*(pp-q*tr))/2 for j in IX] for i in IX]


def einstein_from_connection(g, gi, p, h):
    """Different route: linear Christoffel, its Ricci contraction, then trace."""
    connection = [[[sum(gi[r][s]*(p[i]*h[s][j]+p[j]*h[s][i]-p[s]*h[i][j])
                          for s in IX)/2 for j in IX] for i in IX] for r in IX]
    ricci = [[sum(p[r]*connection[r][i][j]-p[j]*connection[r][i][r]
                  for r in IX) for j in IX] for i in IX]
    scalar = sum(gi[i][j]*ricci[i][j] for i in IX for j in IX)
    return m.add(ricci, m.scale(-scalar/2, g))


def cross(u, v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0]]


def transverse(n):
    require(m.dot(n, n) == 1, 'unit spatial covector')
    u = [-n[1], n[0], Q(0)] if n[0] or n[1] else [Q(1), Q(0), Q(0)]
    v = cross(n, u)
    require(m.dot(n, u) == m.dot(n, v) == m.dot(u, v) == 0
            and m.dot(u, u) == m.dot(v, v) > 0, 'two equal-norm transverse vectors')
    return u, v


def characteristic_checks():
    coframes = [m.eye(4),
                [[Q(1), 0, 0, 0], [0, Q(2), 0, 0],
                 [0, 0, Q(3, 2), 0], [0, 0, 0, Q(5, 4)]],
                [[Q(1), 0, 0, 0], [Q(1, 3), Q(2), Q(1, 4), 0],
                 [Q(-1, 5), 0, Q(3, 2), Q(1, 7)],
                 [Q(2, 7), 0, 0, Q(5, 4)]]]
    directions = [list(map(Q, n)) for n in
                  [(1, 0, 0), (0, 0, 1), (Q(3, 5), Q(4, 5), 0),
                   (Q(2, 3), Q(1, 3), Q(2, 3))]]
    off_cone = [list(map(Q, p)) for p in [(1, 0, 0, 0), (0, 1, 0, 0), (2, 1, 1, 0)]]
    null_count = nonnull_count = ricci_columns = maxwell_columns = 0
    for e in coframes:
        et = m.transpose(e)
        ei = m.inverse(e)
        g = m.mm(m.mm(et, ETA), e)
        gi = m.inverse(g)
        for pf in [[Q(1)]+n for n in directions]+off_cone:
            p = m.mv(et, pf)
            pu = m.mv(gi, p)
            q = quadratic(gi, p)
            require(q == quadratic(ETA, pf), 'coframe preserves the scalar symbol')
            em = maxwell_symbol(gi, p)
            for a in m.eye(4):
                require(m.mv(em, a) == maxwell_from_curvature(gi, p, a), 'Maxwell from F')
                maxwell_columns += 1
            require(m.mv(em, p) == [0]*4, 'Maxwell gauge direction')
            grav_columns = []
            for h in symmetric_basis():
                out = einstein_symbol(g, gi, p, h)
                require(out == einstein_from_connection(g, gi, p, h), 'Einstein from Ricci')
                require(m.mv(out, pu) == [0]*4, 'linearized Bianchi symbol identity')
                grav_columns.append(symmetric_vector(out))
                ricci_columns += 1
            grav = m.transpose(grav_columns)
            gauge = [m.add(outer(p, xi), outer(xi, p)) for xi in m.eye(4)]
            gauge_columns = [symmetric_vector(h) for h in gauge]
            require(m.rank(m.transpose(gauge_columns)) == 4, 'four metric gauge directions')
            require(all(m.mv(grav, h) == [0]*10 for h in gauge_columns), 'metric gauge kernel')
            internal = m.mv(m.transpose(ei), p)
            matter = sm.lin(internal, sm.GAMMA)
            require(m.mm(matter, matter) == m.scale(q, m.eye(8))
                    and m.det(matter) == q**4, 'native matter square and determinant')
            if q:
                require((m.rank(em), m.rank(grav), m.rank(matter)) == (3, 6, 8),
                        'off-cone kernels are only gauge')
                nonnull_count += 1
            else:
                require((m.rank(em), m.rank(grav), m.rank(matter)) == (1, 4, 4),
                        'null characteristic ranks')
                u, v = transverse(pf[1:])
                af = [[Q(0)]+u, [Q(0)]+v]
                photons = [m.mv(et, a) for a in af]
                require(all(m.mv(em, a) == [0]*4 for a in photons)
                        and m.rank(m.transpose([p]+photons)) == 3,
                        'two photon polarizations exhaust kernel modulo gauge')
                plus = m.add(outer(af[0], af[0]), m.scale(-1, outer(af[1], af[1])))
                cross_mode = m.add(outer(af[0], af[1]), outer(af[1], af[0]))
                tensors = [symmetric_vector(m.mm(m.mm(et, h), e)) for h in [plus, cross_mode]]
                require(all(m.mv(grav, h) == [0]*10 for h in tensors)
                        and m.rank(m.transpose(gauge_columns+tensors)) == 6,
                        'two metric polarizations exhaust kernel modulo gauge')
                null_count += 1
    return dict(coframes=len(coframes), null_covectors=null_count, nonnull_covectors=nonnull_count,
                independent_ricci_columns=ricci_columns, independent_maxwell_columns=maxwell_columns,
                null_ranks=dict(maxwell=1, einstein=4, native_matter=4),
                nonnull_ranks=dict(maxwell=3, einstein=6, native_matter=8),
                physical_null_polarizations=dict(photon=2, classical_metric=2),
                quantum_graviton='NOT_DERIVED',
                scope='LINEARIZED_EINSTEIN_VACUUM_WITH_ZERO_BACKGROUND_MATTER_AND_FIELD_STRENGTH')


def speed_family_checks():
    speeds = [Q(1, 2), Q(1), Q(3)]
    n = [Q(3, 5), Q(4, 5), Q(0)]
    for speed in speeds:
        matrices = [m.scale(speed, a) for a in A]
        require(all(x == m.transpose(x) and m.comm(J, x) == m.zero(4)
                    and m.transpose(m.mm(J, x)) == m.scale(-1, m.mm(J, x))
                    for x in matrices), 'every positive speed retains norm and action compatibility')
        e = [[Q(i == j)/(speed if i else 1) for j in IX] for i in IX]
        g = m.mm(m.mm(m.transpose(e), ETA), e)
        gi = m.inverse(g)
        for omega in [Q(0), speed, 2*speed]:
            p = [omega]+n
            symbol = m.add(m.scale(omega, m.eye(4)), m.scale(-1, sm.lin(n, matrices)))
            require(m.det(symbol) == (omega**2-speed**2)**2, 'native speed determinant')
            d = sm.lin([omega]+[speed*x for x in n], GAMMA)
            require(m.mm(d, d) == m.scale(-omega**2+speed**2, m.eye(4)), 'squared wave law')
            if omega == speed:
                require(m.rank(symbol) == 2 and m.rank(maxwell_symbol(gi, p)) == 1,
                        'native and photon null frequency')
                gr = m.transpose([symmetric_vector(einstein_symbol(g, gi, p, h))
                                  for h in symmetric_basis()])
                require(m.rank(gr) == 4, 'metric null frequency with the same coframe')
    # The second theory shares the native algebra, but deliberately not the coframe.
    gi_one, gi_two = ETA, [[Q((-1, 4, 4, 4)[i] if i == j else 0) for j in IX] for i in IX]
    p = [Q(1)]+n
    require(m.rank(maxwell_symbol(gi_one, p)) == 1
            and m.rank(maxwell_symbol(gi_two, p)) == 3,
            'native algebra alone does not force two independently chosen metrics to share a cone')
    calibrated = [v*length/time for v, length, time in
                  [(Q(1), Q(1), Q(1)), (Q(1), Q(3), Q(1)), (Q(2, 3), Q(5), Q(2))]]
    require(len(set(calibrated)) == 3, 'unselected length/time adapter changes reported speed')
    # The inherited lattice model also has a free constitutive ratio before calibration.
    lattice_speed_squared = [nu/eps*length**2/time**2 for nu, eps, length, time in
                             [(Q(1), Q(1), Q(1), Q(1)), (Q(4), Q(1), Q(1), Q(1))]]
    require(lattice_speed_squared == [1, 4], 'PF6 fixed-geometry constitutive nonselection')
    return dict(admitted_native_speeds=[str(x) for x in speeds],
                calibration_fixture_values=[str(x) for x in calibrated],
                calibration_fixture_units='DECLARED_ABSTRACT_LENGTH_PER_TIME_NOT_MEASUREMENTS',
                lattice_speed_squared=[str(x) for x in lattice_speed_squared],
                absolute_SI_speed='NOT_DERIVED',
                countermodel_scope='DIFFERENT_CONSTITUTIVE_METRICS; NOT_A_COUNTEREXAMPLE_TO_SHARED_METRIC_THEOREM')


def massive_checks():
    cases = [(Q(3), Q(4), Q(5)), (Q(2), Q(0), Q(2)),
             (Q(1, 2), Q(2, 3), Q(5, 6)), (Q(3), Q(4), Q(6))]
    k = [Q(3, 5), Q(4, 5), Q(0)]
    shell_count = 0
    for speed, mass, omega in cases:
        d = sm.lin([-omega]+[speed*x for x in k], sm.GAMMA)
        t = m.mm(sm.PHASE, d)
        op = m.add(t, m.scale(-mass, m.eye(8)))
        conjugate = m.add(t, m.scale(mass, m.eye(8)))
        dispersion = omega**2-speed**2*m.dot(k, k)-mass**2
        require(m.mm(op, conjugate) == m.scale(dispersion, m.eye(8))
                and m.det(op) == dispersion**4, 'real native phase yields massive dispersion')
        if not dispersion:
            require(m.rank(op) == 4, 'four real on-shell amplitude directions')
            group_squared = speed**4*m.dot(k, k)/omega**2
            require(group_squared <= speed**2 and (group_squared < speed**2) == (mass > 0),
                    'massive group speed is below the wavefront speed')
            shell_count += 1
    return dict(dispersion_cases=len(cases), on_shell_cases=shell_count,
                group_front_distinction='MASS_IS_LOWER_ORDER; FINITE_MOMENTUM_GROUP_SPEED_IS_SMALLER_IF_MASS_POSITIVE',
                units='NATIVE_DIMENSIONLESS_CHART; NO_HBAR_OR_PARTICLE_MASS_CALIBRATION')


def reference_comparison():
    reference = json.loads((Path(__file__).with_name('SPEED_REFERENCE_DATA.json')).read_text())
    observation = reference['gw170817']
    # Derived in SC6 from the physical characteristic kernels, not fitted to this interval.
    ratio = Q(1)
    delta = ratio-1
    lower, upper = map(Q, observation['fractional_speed_interval'])
    require(lower <= delta <= upper, 'SC6 equal-speed result lies inside the published interval')
    observed_delay = Q(observation['gamma_minus_gw_arrival_seconds'])
    delay_error = Q(observation['arrival_uncertainty_seconds'])
    emission_lower, emission_upper = map(Q, observation['intrinsic_emission_endpoints_seconds'])
    require(emission_lower <= observed_delay-delay_error <= observed_delay+delay_error <= emission_upper,
            'zero propagation delay is compatible with, but does not predict, the source lag')
    require(reference['SI']['c_m_per_s'] == 299792458
            and reference['SI']['status'] == 'EXACT_DEFINING_CONSTANT_NOT_AN_INDEPENDENT_FIT_TARGET',
            'SI definition is labelled explicitly')
    return dict(speed_ratio=str(ratio), fractional_difference=str(delta),
                published_fractional_interval=[str(lower), str(upper)],
                comparison='CONSISTENT_WITH_PUBLISHED_EMISSION_MODEL_DEPENDENT_BOUND',
                origin='CONDITIONAL_SHARED_METRIC_EINSTEIN_MAXWELL_REDUCTION; NOT_PRIMITIVE_ONLY_SELECTION',
                SI_value_status=reference['SI']['status'],
                raw_data_reanalysis=False, source_lag_prediction=False,
                independent_discriminator_against_GR=False, experimental_fit_parameters=0)


def run():
    return dict(characteristics=characteristic_checks(), calibration=speed_family_checks(),
                massive=massive_checks(), comparison=reference_comparison())


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
