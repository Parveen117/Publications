"""SM exact classical matter controls; consumes unchanged native fixtures."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product

import model as m
import metric_dynamics as mg
import loop_curvature as lr
from thermo_gauge import Jet

IX, PAIRS, TOP, ETA = range(4), mg.PAIRS, mg.TOP, mg.ETA
require = mg.require
NATIVE, OLD_GAMMA, OLD_A, OLD_J = m.fixture()
GAMMA = [m.kron(m.eye(2), x) for x in OLD_GAMMA]  # UPPER Lorentz index
A = [m.kron(m.eye(2), x) for x in OLD_A]
PHASE = m.kron(NATIVE[1], m.eye(4))
VOLUME = m.kron(m.eye(2), OLD_J)
H = m.scale(-1, m.mm(PHASE, GAMMA[0]))
HG = [m.mm(H, x) for x in GAMMA]
HG_PHASE = [m.mm(x, PHASE) for x in HG]
SPIN = [m.scale(Q(ETA[a]*ETA[b], 2), m.mm(GAMMA[a], GAMMA[b]))
        for a, b in PAIRS]  # independent upper omega^{ab}, a<b
HGS = [[m.mm(x, y) for y in SPIN] for x in HG]
G3 = {abc: (m.scale(mg.sign(abc),
                   m.mm(m.mm(GAMMA[sorted(abc)[0]], GAMMA[sorted(abc)[1]]),
                        GAMMA[sorted(abc)[2]])) if len(set(abc)) == 3 else m.zero(8))
      for abc in product(IX, repeat=3)}
HG3 = {abc: m.mm(H, x) for abc, x in G3.items()}


def val(x):
    return x.value if isinstance(x, Jet) else x


def bil(x, matrix, y=None):
    y = x if y is None else y
    return sum(x[i]*c*y[j] for i, row in enumerate(matrix)
               for j, c in enumerate(row) if c)


def lin(coefficients, basis):
    return [[sum(c*x[i][j] for c, x in zip(coefficients, basis))
             for j in range(len(basis[0]))] for i in range(len(basis[0]))]


def axial(psi):
    return {abc: bil(psi, x) for abc, x in HG3.items()}


def density(e, omega, gauge, psi, dpsi, mass, norm, charge):
    """w Z/2 [bar(psi) gamma^mu D_mu psi - m bar(psi) psi]."""
    inverse = mg.inverse_generic(e)
    spin_bilinear = [[bil(psi, x) for x in row] for row in HGS]
    gauge_bilinear = [bil(psi, x) for x in HG_PHASE]
    kinetic = sum(inverse[mu][a]*(bil(psi, HG[a], dpsi[mu])
                  +sum(omega[p][mu]*spin_bilinear[a][p] for p in range(6))
                  +charge*gauge[mu]*gauge_bilinear[a]) for a, mu in product(IX, repeat=2))
    return norm*mg.determinant_generic(e)*(kinetic-mass*bil(psi, H))/2


def connections(omega, gauge, charge):
    return [m.add(lin([omega[p][mu] for p in range(6)], SPIN),
                  m.scale(charge*gauge[mu], PHASE)) for mu in IX]


def field_fixture(seed):
    ef, om, tf, _ = lr.jet_fixture(seed)
    e = [[ef[a][mu,] for mu in IX] for a in IX]
    omega = [[om[a][b][mu,] for mu in IX] for a, b in PAIRS]
    gauge = [Jet(Q(seed-mu, 13), tuple(Q(seed+mu-k, 17) for k in IX)) for mu in IX]
    psi = [Jet(Q((j*j+seed*j+1) % 11-5, 7),
               tuple(Q(seed+(j+1)*(mu+1), 19) for mu in IX)) for j in range(8)]
    dpsi = [[Jet(psi[j].grad[mu], tuple(Q(seed+j+(mu+1)*(nu+1), 23) for nu in IX))
             for j in range(8)] for mu in IX]
    return e, omega, gauge, psi, dpsi, tf


def quadratic_polynomial(matrix):
    return {(i, j): matrix[i][j]+(matrix[j][i] if i != j else 0)
            for i in range(8) for j in range(i, 8)
            if matrix[i][j]+(matrix[j][i] if i != j else 0)}


def square_polynomial(poly):
    result = {}
    for a, x in poly.items():
        for b, y in poly.items():
            key = tuple(sorted(a+b))
            result[key] = result.get(key, Q(0))+x*y
    return result


def algebra_checks():
    require(m.transpose(H) == H and m.mm(H, H) == m.eye(8), 'nondegenerate symmetric H')
    require(all(m.transpose(x) == m.scale(-1, x) for x in HG), 'kinetic bilinears skew')
    require(all(m.add(m.mm(m.transpose(s), H), m.mm(H, s)) == m.zero(8)
                for s in SPIN+[PHASE]), 'Lorentz and phase bilinear invariance')
    require(all(m.comm(PHASE, x) == m.zero(8) for x in GAMMA), 'mass-compatible native phase')
    require(m.mm(H, GAMMA[0]) == PHASE
            and all(HG[i+1] == m.scale(-1, m.mm(PHASE, A[i])) for i in range(3)),
            'curved action reduces to Ical B_m with NP signs')
    require(all(m.transpose(x) == x for x in HG3.values()), 'only antisymmetric triple survives')
    # A full 64-unknown multiplier problem, with no assumed tensor ansatz.
    columns, skew_columns = [], []
    for i, j in product(range(8), repeat=2):
        unit = m.zero(8)
        unit[i][j] = Q(1)
        columns.append([v for x in A+[GAMMA[0]] for row in m.comm(unit, x) for v in row])
        skew_columns.append(sum(m.add(unit, m.transpose(unit)), []))
    rank = m.rank(m.transpose(columns))
    skew_rank = m.rank(m.transpose(columns)+m.transpose(skew_columns))
    require((rank, skew_rank) == (60, 63), 'doubled full and skew commutant dimensions')
    similarity = GAMMA[0]
    for (a, b), s in zip(PAIRS, SPIN):
        old_s = m.kron(m.eye(2), m.scale(Q(1, 2), m.mm(OLD_GAMMA[a], OLD_GAMMA[b])))
        require(s == m.mm(m.mm(similarity, old_s), m.scale(-1, similarity)),
                'matter spin connection is a fixed similarity of MG')

    # Exhaust every quartic monomial, not just sampled field values.
    scalar = square_polynomial(quadratic_polynomial(H))
    pseudoscalar = square_polynomial(quadratic_polynomial(m.mm(H, VOLUME)))
    vectors = [square_polynomial(quadratic_polynomial(m.scale(-1, x))) for x in HG_PHASE]
    triples = {abc: square_polynomial(quadratic_polynomial(x)) for abc, x in HG3.items()}
    monomials = list(combinations_with_replacement(range(8), 4))
    for key in monomials:
        vv = sum(ETA[a]*vectors[a].get(key, 0) for a in IX)
        aa = sum(ETA[a]*ETA[b]*ETA[c]*p.get(key, 0) for (a, b, c), p in triples.items())
        require(vv+scalar.get(key, 0)+pseudoscalar.get(key, 0) == 0, 'universal current Fierz polynomial')
        require(aa == 6*vv, 'universal axial Fierz polynomial')

    rest = m.mm(PHASE, GAMMA[0])
    p, q = next((p, q) for p, q in product(range(8), repeat=2) if p != q and rest[p][q])
    plus, minus = [Q(0)]*8, [Q(0)]*8
    plus[p] = minus[p] = plus[q] = Q(1)
    minus[q] = Q(-1)
    energies = [-bil(x, rest)/2 for x in (plus, minus)]
    require(energies[0] == -energies[1] != 0 and m.dot(plus, plus) == m.dot(minus, minus),
            'massive positive norm does not make Hamiltonian positive')
    require(m.transpose(m.mm(H, PHASE)) == m.scale(-1, m.mm(H, PHASE)),
            'positive-Lambda Cartan slashed phase is erased in a commuting quadratic action')
    translations = [m.kron(NATIVE[1], x) for x in OLD_GAMMA]
    moved = [m.mm(m.mm(similarity, x), m.scale(-1, similarity)) for x in translations]
    require(all(x == m.scale(-ETA[a], m.mm(PHASE, GAMMA[a])) for a, x in enumerate(moved)),
            'same similarity maps the positive-Lambda translation sector')
    slash = lin([Q(1)]*4, [m.mm(GAMMA[a], moved[a]) for a in IX])
    require(slash == m.scale(-4, PHASE), 'naive Cartan slash has phase, not scalar mass')
    benchmark = list(map(Q, [2, 0, 1, 0, 0, 0, 3, 0]))
    bench_ax = axial(benchmark)
    require([-bil(benchmark, x) for x in HG_PHASE] == [14, -4, 0, 6]
            and bil(benchmark, H) == -12 and bil(benchmark, m.mm(H, VOLUME)) == 0,
            'reported massive-current rational benchmark')
    require({abc: bench_ax[abc] for abc in mg.TRIPLES} ==
            {(0, 1, 2): 0, (0, 1, 3): 12, (0, 2, 3): 0, (1, 2, 3): 0},
            'reported axial rational benchmark')
    return dict(real_components=8, commutant_unknowns=64, commutant_rank=rank,
                skew_commutant_rank=skew_rank, quartic_monomials_per_identity=len(monomials),
                fierz_identities=2, rest_energy_witness=[str(x) for x in energies],
                naive_positive_Lambda_mass_identification_rejected=True,
                benchmark=dict(vector=[14, -4, 0, 6], scalar=-12, pseudoscalar=0,
                               vector_square=-144, axial_square=-864))


def torsion_trace(e, tf):
    inverse = mg.inverse_generic(e)
    return [sum(inverse[nu][a]*mg.sign((mu, nu))*tf[a].get(tuple(sorted((mu, nu))), 0)
                for nu, a in product(IX, repeat=2) if mu != nu) for mu in IX]


def field_and_source_checks():
    mass, norm, charge = Q(3, 5), Q(7, 3), Q(2, 7)
    euler_count = source_count = trace_controls = 0
    for seed in (1, 2):
        e, omega, gauge, psi, dpsi, tf = field_fixture(seed)
        w, inv = mg.determinant_generic(e), mg.inverse_generic(e)
        require(w.value > 0, 'oriented fixture')
        gam = [lin(inv[mu], GAMMA) for mu in IX]
        conn = connections(omega, gauge, charge)
        cov = [[x+y for x, y in zip(dpsi[mu], m.mv(conn[mu], psi))] for mu in IX]
        trace_t = torsion_trace(e, tf)
        expected = [sum(m.mv(gam[mu], cov[mu])[j]
                        +trace_t[mu]*m.mv(gam[mu], psi)[j]/2 for mu in IX)-mass*psi[j]
                    for j in range(8)]
        expected = [val(norm*w*x) for x in m.mv(H, expected)]
        actual = []
        for j in range(8):
            pp, pm = psi[:], psi[:]
            pp[j], pm[j] = pp[j]+1, pm[j]-1
            partial = (density(e, omega, gauge, pp, dpsi, mass, norm, charge)
                       -density(e, omega, gauge, pm, dpsi, mass, norm, charge))/2
            divergence = Q(0)
            for mu in IX:
                dp, dm = [x[:] for x in dpsi], [x[:] for x in dpsi]
                dp[mu][j], dm[mu][j] = dp[mu][j]+1, dm[mu][j]-1
                momentum = (density(e, omega, gauge, psi, dp, mass, norm, charge)
                            -density(e, omega, gauge, psi, dm, mass, norm, charge))/2
                divergence += momentum.grad[mu]
            actual.append(partial.value-divergence)
            euler_count += 1
        require(actual == expected, 'direct variational EL versus torsion-trace Dirac operator')
        trace_correction = [val(x) for x in m.mv(H, [norm*w*sum(trace_t[mu]*m.mv(gam[mu], psi)[j]/2
                                                 for mu in IX) for j in range(8)])]
        require(any(trace_correction), 'omitting off-shell torsion trace changes EL')
        trace_controls += 1
        current = [-norm*charge*bil(psi, m.mm(H, m.mm(gam[mu], PHASE)))/2 for mu in IX]
        for mu in IX:
            ap, am = gauge[:], gauge[:]
            ap[mu], am[mu] = ap[mu]+1, am[mu]-1
            source = (density(e, omega, ap, psi, dpsi, mass, norm, charge)
                      -density(e, omega, am, psi, dpsi, mass, norm, charge))/2
            require(source == -w*current[mu], 'independent gauge source variation including first jets')
            source_count += 1
        divergence = sum((w*current[mu]).grad[mu] for mu in IX)
        require(divergence == charge*m.dot(actual, m.mv(PHASE, [val(x) for x in psi])),
                'off-shell U(1) Noether identity with independently varied EL')
        phase = m.add(m.scale(Q(3, 5), m.eye(8)), m.scale(Q(4, 5), PHASE))
        chi = [Q(seed-mu, 11) for mu in IX]
        moved_psi = m.mv(phase, [val(x) for x in psi])
        moved_dpsi = [m.mv(phase, [val(x)-charge*chi[mu]*y for x, y in
                                   zip(dpsi[mu], m.mv(PHASE, [val(z) for z in psi]))]) for mu in IX]
        ev, ov = lr.matrix_values(e), lr.matrix_values(omega)
        av, pv, dv = [val(x) for x in gauge], [val(x) for x in psi], lr.matrix_values(dpsi)
        require(density(ev, ov, [av[mu]+chi[mu] for mu in IX], moved_psi, moved_dpsi, mass, norm, charge)
                == density(ev, ov, av, pv, dv, mass, norm, charge), 'local phase transformation of action density')
    return dict(independent_euler_components=euler_count, gauge_source_components=source_count,
                off_shell_noether_jets=2, local_phase_covariance_jets=2,
                missing_torsion_trace_failures=trace_controls)


def coframe_and_spin_checks():
    mass, norm, charge = Q(3, 5), Q(7, 3), Q(2, 7)
    coframe_count = spin_count = 0
    for seed in (1, 2):
        ej, oj, aj, pj, dpj, _ = field_fixture(seed)
        e, omega, gauge = lr.matrix_values(ej), lr.matrix_values(oj), [val(x) for x in aj]
        psi, dpsi = [val(x) for x in pj], lr.matrix_values(dpj)
        w, inv = m.det(e), m.inverse(e)
        gam = [lin(inv[mu], GAMMA) for mu in IX]
        conn = connections(omega, gauge, charge)
        cov = [[x+y for x, y in zip(dpsi[mu], m.mv(conn[mu], psi))] for mu in IX]
        lag = density(e, omega, gauge, psi, dpsi, mass, norm, charge)/w
        for a, mu in product(IX, repeat=2):
            ep, em = [row[:] for row in e], [row[:] for row in e]
            ep[a][mu], em[a][mu] = ep[a][mu]+Q(1, 13), em[a][mu]-Q(1, 13)
            actual = (density(ep, omega, gauge, psi, dpsi, mass, norm, charge)
                      -density(em, omega, gauge, psi, dpsi, mass, norm, charge))*Q(13, 2)
            expected = w*(inv[mu][a]*lag-norm*sum(inv[nu][a]*bil(psi, m.mm(H, gam[mu]), cov[nu])
                                                for nu in IX)/2)
            require(actual == expected, 'independent coframe variation versus canonical source')
            coframe_count += 1
        ax = axial(psi)
        for p, (a, b) in enumerate(PAIRS):
            for mu in IX:
                op, om = [row[:] for row in omega], [row[:] for row in omega]
                op[p][mu], om[p][mu] = op[p][mu]+1, om[p][mu]-1
                actual = (density(e, op, gauge, psi, dpsi, mass, norm, charge)
                          -density(e, om, gauge, psi, dpsi, mass, norm, charge))/2
                expected = norm*w*ETA[a]*ETA[b]*sum(inv[mu][c]*ax[c, a, b] for c in IX)/4
                require(actual == expected, 'independent 24-component spin source versus axial three-form')
                spin_count += 1
    return dict(independent_coframe_variations=coframe_count, independent_spin_variations=spin_count)


def lorentz_covariance_checks():
    """Finite rational spin transformations with arbitrary first derivatives."""
    vector_basis = []
    for a, b in PAIRS:
        unit = m.zero(4)
        unit[a][b], unit[b][a] = Q(ETA[b]), Q(-ETA[a])
        vector_basis.append(unit)
    ej, oj, aj, pj, dpj, _ = field_fixture(2)
    e, omega = lr.matrix_values(ej), lr.matrix_values(oj)
    gauge, psi, dpsi = [val(x) for x in aj], [val(x) for x in pj], lr.matrix_values(dpj)
    mass, norm, charge = Q(3, 5), Q(7, 3), Q(2, 7)
    original = density(e, omega, gauge, psi, dpsi, mass, norm, charge)
    cases = [(0, Q(5, 3), Q(4, 3)), (3, Q(3, 5), Q(4, 5))]
    for p, co, si in cases:
        spin = m.add(m.scale(co, m.eye(8)), m.scale(2*si, SPIN[p]))
        si_inverse = m.mm(m.mm(H, m.transpose(spin)), H)
        require(m.mm(spin, si_inverse) == m.eye(8), 'rational Spin transformation')
        moved_gamma = [m.mm(m.mm(spin, x), si_inverse) for x in GAMMA]
        li = [[ETA[b]*m.trace(m.mm(GAMMA[b], moved_gamma[a]))/8 for b in IX] for a in IX]
        lorentz = m.inverse(li)
        require(mg.metric(lorentz) == mg.metric(m.eye(4)) and m.det(lorentz) == 1,
                'induced proper Lorentz frame transformation')
        moved_e, moved_psi = m.mm(lorentz, e), m.mv(spin, psi)
        moved_omega, moved_dp = [[Q(0)]*4 for _ in PAIRS], []
        for mu in IX:
            coefficients = [Q((mu+1)*(k+1)-3, 31) for k in range(6)]
            derivative_spin = lin(coefficients, SPIN)  # (d_mu S) S^{-1}
            derivative_lorentz = lin(coefficients, vector_basis)
            vector_omega = lin([omega[k][mu] for k in range(6)], vector_basis)
            transformed = m.add(m.mm(m.mm(lorentz, vector_omega), li), m.scale(-1, derivative_lorentz))
            for k, (a, b) in enumerate(PAIRS):
                moved_omega[k][mu] = ETA[b]*transformed[a][b]
            moved_dp.append([x+y for x, y in zip(m.mv(spin, dpsi[mu]), m.mv(derivative_spin, moved_psi))])
        require(density(moved_e, moved_omega, gauge, moved_psi, moved_dp, mass, norm, charge) == original,
                'local Lorentz action covariance including inhomogeneous connection')
        inv = m.inverse(e)
        moved_inv = m.inverse(moved_e)
        require(all(lin(moved_inv[mu], GAMMA) == m.mm(m.mm(spin, lin(inv[mu], GAMMA)), si_inverse)
                    for mu in IX), 'world Clifford vector covariance')
    return dict(local_lorentz_jets=2, nonzero_connection_inhomogeneities=8,
                rotation_and_boost_checked=True)


def connection_curvature(omega):
    """Constant-coframe, constant upper-connection curvature, evaluated directly."""
    om = [[{} for _ in IX] for _ in IX]
    for p, (a, b) in enumerate(PAIRS):
        om[a][b] = {(mu,): omega[p][mu] for mu in IX}
        om[b][a] = mg.fscale(-1, om[a][b])
    rf = [[mg.fadd(*(mg.fscale(ETA[c], mg.wedge(om[a][c], om[c][b])) for c in IX))
           for b in IX] for a in IX]
    spin = mg.matrix_form([rf[a][b] for a, b in PAIRS],
                         [m.scale(Q(1, 2), m.mm(OLD_GAMMA[a], OLD_GAMMA[b])) for a, b in PAIRS])
    return spin, om


def contorsion_checks():
    kappa, norm = Q(5, 3), Q(7, 3)
    e, zero_dp, zero_a = m.eye(4), [[Q(0)]*8 for _ in IX], [Q(0)]*4
    zero_om = [[Q(0)]*4 for _ in PAIRS]

    def gravitational(om):
        curvature, _ = connection_curvature(om)
        return mg.gravity_density(e, curvature, Q(0), Q(0), kappa)

    # Quadratic Hessian of ALL 24 connection coefficients, not only the axial ansatz.
    units, diagonal = [], []
    for p, mu in product(range(6), IX):
        om = [row[:] for row in zero_om]
        om[p][mu] = Q(1)
        units.append(om)
        diagonal.append(gravitational(om))
    hessian = [[gravitational(m.add(a, b))-diagonal[i]-diagonal[j]
                for j, b in enumerate(units)] for i, a in enumerate(units)]
    require(m.rank(hessian) == 24, 'full connection quadratic is nonsingular')
    stationary_count = 0
    values = []
    for seed in (1, 2, 3):
        psi = [Q((j*j+seed*j+1) % 11-5, 7) for j in range(8)]
        ax = axial(psi)
        # K_lower_abc=kappa Z A_lower_abc/8; omega_upper_ab,c raises a,b.
        om = [[kappa*norm*ETA[c]*ax[a, b, c]/8 for c in IX] for a, b in PAIRS]
        source = [density(e, u, zero_a, psi, zero_dp, Q(0), norm, Q(0)) for u in units]
        residual = [x+y for x, y in zip(m.mv(hessian, sum(om, [])), source)]
        require(residual == [0]*24, 'closed torsion solution stationary in every connection direction')
        actual = gravitational(om)+density(e, om, zero_a, psi, zero_dp, Q(0), norm, Q(0))
        aa = sum(ETA[a]*ETA[b]*ETA[c]*x*x for (a, b, c), x in ax.items())
        expected = kappa*norm*norm*aa/128
        require(actual == expected, 'direct first-order action elimination versus quartic interaction')
        _, om_forms = connection_curvature(om)
        ef = mg.coframe_forms(e)
        tf = [mg.fadd(*(mg.fscale(ETA[b], mg.wedge(om_forms[a][b], ef[b])) for b in IX)) for a in IX]
        require(torsion_trace(e, tf) == [0]*4, 'eliminated spin torsion has zero trace')
        for c, d in PAIRS:
            left = mg.fadd(*(mg.fscale(mg.sign((a, b, c, d)), mg.wedge(tf[a], ef[b]))
                            for a, b in product(IX, repeat=2)))
            right = {tuple(i for i in IX if i != a):
                     -kappa*norm*ETA[c]*ETA[d]*ax[a, c, d]*(-1)**a/4 for a in IX}
            require(mg.equal_forms(left, right), 'Cartan three-form source sign and normalization')
        for a, b, c in product(IX, repeat=3):
            actual_t = ETA[a]*mg.sign((b, c))*tf[a].get(tuple(sorted((b, c))), Q(0)) if b != c else Q(0)
            require(actual_t == -kappa*norm*ETA[a]*ETA[b]*ETA[c]*ax[a, b, c]/4,
                    'torsion component sign with T=de+omega e')
        conn = connections(om, zero_a, Q(0))
        cubic = [kappa*norm*sum(ETA[a]*ETA[b]*ETA[c]*ax[a, b, c]*m.mv(G3[a, b, c], psi)[j]
                               for a, b, c in product(IX, repeat=3))/32 for j in range(8)]
        require(cubic == [sum(m.mv(GAMMA[mu], m.mv(conn[mu], psi))[j] for mu in IX) for j in range(8)],
                'cubic effective equation equals eliminated covariant connection')
        for j in range(8):
            pjet = [Jet(x, (Q(i == j), Q(0), Q(0), Q(0))) for i, x in enumerate(psi)]
            ajet = axial(pjet)
            contact = kappa*norm*norm*sum(ETA[a]*ETA[b]*ETA[c]*x*x for (a, b, c), x in ajet.items())/128
            require(contact.grad[0] == norm*m.mv(H, cubic)[j], 'exact quartic variation gives cubic source')
        stationary_count += 24
        values.append(str(expected))
    return dict(full_connection_hessian_rank=24, connection_unknowns=24,
                stationary_components=stationary_count, quartic_action_checks=3,
                effective_cubic_variations=24, torsion_trace_zero=True,
                rational_contact_densities=values)


def run():
    return dict(algebra=algebra_checks(), matter=field_and_source_checks(),
                sources=coframe_and_spin_checks(), covariance=lorentz_covariance_checks(),
                torsion=contorsion_checks())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
