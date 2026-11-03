"""MG exact controls: publication exterior algebra, not a native operator engine."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product

import model as m
from thermo_gauge import Jet, square_root

IX = range(4)
PAIRS = list(combinations(IX, 2))
TRIPLES = list(combinations(IX, 3))
ETA = [-1, 1, 1, 1]
TOP = (0, 1, 2, 3)
_, GAMMA, _, J = m.fixture()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sign(indices):
    if len(set(indices)) != len(indices):
        return 0
    return (-1)**sum(indices[i] > indices[k]
                     for i in range(len(indices)) for k in range(i+1, len(indices)))


def wedge(left, right, matrix=False):
    result = {}
    for a, x in left.items():
        for b, y in right.items():
            s = sign(a+b)
            if not s:
                continue
            key = tuple(sorted(a+b))
            value = m.scale(s, m.mm(x, y)) if matrix else s*x*y
            if key in result:
                value = m.add(result[key], value) if matrix else result[key]+value
            result[key] = value
    return result


def fadd(*forms):
    result = {}
    for form in forms:
        for key, value in form.items():
            result[key] = result.get(key, Q(0))+value
    return result


def fscale(c, form):
    return {key: c*value for key, value in form.items()}


def exterior_d(form):
    result = {}
    for key, value in form.items():
        if isinstance(value, Jet):
            for i in IX:
                result = fadd(result, wedge({(i,): value.grad[i]}, {key: Q(1)}))
    return result


def values(form):
    return {key: (value.value if isinstance(value, Jet) else value)
            for key, value in form.items()}


def equal_forms(a, b):
    return all(x == 0 for x in fadd(values(a), fscale(-1, values(b))).values())


def matrix_form(forms, basis):
    result = {}
    for form, coefficient in zip(forms, basis):
        for key, value in form.items():
            result[key] = m.add(result.get(key, m.zero(4)), m.scale(value, coefficient))
    return result


def native_word(word, forms, insertion=None):
    result = {(): m.eye(4)}
    for letter in word:
        result = wedge(result, forms[letter], matrix=True)
    coefficient = result.get(TOP, m.zero(4))
    return m.trace(m.mm(insertion or m.eye(4), coefficient))/4


def coframe_forms(e):
    return [{(mu,): e[a][mu] for mu in IX} for a in IX]


def metric(e):
    return m.mm(m.mm(m.transpose(e), [[Q(ETA[a]) if a == b else Q(0)
                                       for b in IX] for a in IX]), e)


def coframe_fixtures():
    return [m.eye(4),
            [[Q(1), Q(1, 5), Q(0), Q(0)],
             [Q(0), Q(2), Q(1, 3), Q(0)],
             [Q(1, 7), Q(0), Q(1), Q(1, 4)],
             [Q(0), Q(0), Q(0), Q(3, 2)]],
            [[Q(i == k)+Q((i+2*k) % 3-1, 400) for k in IX] for i in IX]]


def trace_and_action_class_checks():
    trace_count = 0
    eta = lambda a, b: ETA[a] if a == b else 0
    for a, b, c, d in product(IX, repeat=4):
        word = m.mm(m.mm(GAMMA[a], GAMMA[b]), m.mm(GAMMA[c], GAMMA[d]))
        require(m.trace(word)/4 == eta(a, b)*eta(c, d)-eta(a, c)*eta(b, d)+eta(a, d)*eta(b, c),
                'ordinary Clifford fourth trace')
        require(m.trace(m.mm(J, word))/4 == -sign((a, b, c, d)), 'oriented native trace')
        trace_count += 2
    bivectors = [m.mm(GAMMA[a], GAMMA[b]) for a, b in PAIRS]
    require(m.rank([[m.trace(m.mm(x, y))/4 for y in bivectors] for x in bivectors]) == 6,
            'bivector trace pairing nondegenerate')
    require(all(m.comm(J, b) == m.zero(4) for b in bivectors), 'J is spin-invariant')
    for a, b in PAIRS:
        for c in IX:
            expected = m.scale(2*eta(b, c), GAMMA[a])
            expected = m.add(expected, m.scale(-2*eta(a, c), GAMMA[b]))
            require(m.comm(m.mm(GAMMA[a], GAMMA[b]), GAMMA[c]) == expected,
                    'bivector action is Lorentz action')

    words = ['EEEE', 'EET', 'ETE', 'TEE', 'EEF', 'EFE', 'FEE', 'TT', 'TF', 'FT', 'FF']
    counts = 0
    for e in coframe_fixtures():
        ef = coframe_forms(e)
        tf = [{pair: Q((a+1)*(mu+2)-nu, 7) for mu, nu in PAIRS for pair in [(mu, nu)]}
              for a in IX]
        rf = [{pair: Q((a+1)*(nu+2)-(b+1)*(mu+1), 11)
               for mu, nu in PAIRS for pair in [(mu, nu)]} for a, b in PAIRS]
        forms = dict(E=matrix_form(ef, GAMMA), T=matrix_form(tf, GAMMA),
                     F=matrix_form(rf, [m.scale(Q(1, 2), b) for b in bivectors]))
        outputs = {(word, parity): native_word(word, forms, J if parity else None)
                   for word, parity in product(words, (False, True))}
        volume = wedge(wedge(wedge(ef[0], ef[1]), ef[2]), ef[3]).get(TOP, Q(0))
        holst = 2*sum(ETA[a]*ETA[b]*wedge(wedge(ef[a], ef[b]), r).get(TOP, Q(0))
                      for (a, b), r in zip(PAIRS, rf))
        palatini = sum(2*sign((a, b, c, d))*wedge(wedge(ef[a], ef[b]), rf[PAIRS.index((c, d))]).get(TOP, Q(0))
                       for a, b in product(IX, repeat=2) for c, d in PAIRS)
        torsion = sum(ETA[a]*wedge(tf[a], tf[a]).get(TOP, Q(0)) for a in IX)
        pontryagin_trace = -sum(ETA[a]*ETA[b]*wedge(r, r).get(TOP, Q(0))
                                for (a, b), r in zip(PAIRS, rf))/4
        euler_trace = -sum(sign((a, b, c, d))*wedge(r, s).get(TOP, Q(0))
                          for (a, b), r in zip(PAIRS, rf)
                          for (c, d), s in zip(PAIRS, rf))/4
        require(outputs['EEEE', True] == -24*volume and outputs['EEEE', False] == 0, 'volume trace')
        require(outputs['EEF', True] == -palatini/4, 'epsilon curvature trace')
        require(outputs['EEF', False] == -holst/2, 'metric curvature trace')
        require(outputs['TT', False] == torsion and outputs['TT', True] == 0, 'torsion traces')
        require(outputs['FF', False] == pontryagin_trace and outputs['FF', True] == euler_trace,
                'characteristic traces versus independent index contractions')
        for word in ['EET', 'ETE', 'TEE', 'TF', 'FT']:
            require(outputs[word, False] == outputs[word, True] == 0, 'odd Clifford traces vanish')
        for parity in (False, True):
            require(outputs['FEE', parity] == outputs['EEF', parity]
                    and outputs['EFE', parity] == (1 if parity else -1)*outputs['EEF', parity],
                    'graded trace cyclicity including J anticommutation with E')
        counts += len(outputs)
    for gamma in (Q(0), Q(2, 3), Q(-5)):
        factor = m.add(J, m.scale(gamma, m.eye(4)))
        reciprocal = m.scale(1/(1+gamma*gamma), m.add(m.scale(gamma, m.eye(4)), m.scale(-1, J)))
        require(m.mm(factor, reciprocal) == m.eye(4), 'real Holst multiplier invertible')
    return dict(fourth_trace_entries=trace_count, spin_vector_actions=24, bivector_pairing_rank=6,
                degree_four_words=len(words), trace_word_controls=counts, real_multiplier_checks=3)


def nieh_yan_checks():
    """Differentiate e_a T^a directly, without replacing dT by a Bianchi formula."""
    count = 0
    for seed in (1, 2, 3):
        e = [[Jet(Q(a == mu)+Q(seed+a-mu, 17),
                  tuple(Q(seed+(a+1)*(mu+1)-k, 23) for k in IX)) for mu in IX] for a in IX]
        ef = coframe_forms(e)
        omega_up = [[{} for _ in IX] for _ in IX]
        for a, b in PAIRS:
            omega_up[a][b] = {(mu,): Jet(Q(seed+a+b-mu, 19),
                                       tuple(Q((a+1)*(k+1)-b+mu+seed, 29) for k in IX)) for mu in IX}
            omega_up[b][a] = fscale(-1, omega_up[a][b])
        tf = [fadd(exterior_d(ef[a]), *(fscale(ETA[b], wedge(omega_up[a][b], ef[b])) for b in IX))
              for a in IX]
        rf = [[fadd(exterior_d(omega_up[a][b]),
                     *(fscale(ETA[c], wedge(omega_up[a][c], omega_up[c][b])) for c in IX))
               for b in IX] for a in IX]
        current = fadd(*(fscale(ETA[a], wedge(ef[a], tf[a])) for a in IX))
        torsion_sq = fadd(*(fscale(ETA[a], wedge(tf[a], tf[a])) for a in IX))
        curvature = fadd(*(fscale(ETA[a]*ETA[b], wedge(wedge(ef[a], ef[b]), rf[a][b]))
                          for a, b in product(IX, repeat=2)))
        require(equal_forms(exterior_d(current), fadd(torsion_sq, fscale(-1, curvature))),
                'Nieh-Yan from independent first jets')
        require(not equal_forms(torsion_sq, curvature), 'nonzero boundary density control')
        count += 1
    return dict(independent_first_jet_checks=count, boundary_term_cannot_be_omitted_pointwise=True)


def torsion_map(e):
    ef = coframe_forms(e)
    columns = []
    for a, pair in product(IX, PAIRS):
        tf = [{} for _ in IX]
        tf[a] = {pair: Q(1)}
        residuals = [fadd(wedge(ef[c], tf[d]), fscale(-1, wedge(ef[d], tf[c]))) for c, d in PAIRS]
        columns.append([form.get(triple, Q(0)) for form in residuals for triple in TRIPLES])
    return m.transpose(columns)


def dq(e, de):
    eta = [[Q(ETA[a]) if a == b else Q(0) for b in IX] for a in IX]
    left = m.mm(m.mm(m.transpose(e), eta), de)
    return m.add(left, m.transpose(left))


def variation_rank_checks():
    rows = []
    symmetric_pairs = [(a, b) for a in IX for b in range(a, 4)]
    for e in coframe_fixtures():
        q = metric(e)
        qi = m.inverse(q)
        require(m.det(q) == -m.det(e)**2 and m.det(e) > 0, 'metric determinant and orientation')
        columns = []
        for a, mu in product(IX, repeat=2):
            de = m.zero(4)
            de[a][mu] = Q(1)
            delta = dq(e, de)
            columns.append([delta[i][k] for i, k in symmetric_pairs])
        require(m.rank(m.transpose(columns)) == 10, 'complete metric variations')
        for i, k in symmetric_pairs:
            h = m.zero(4)
            h[i][k] = h[k][i] = Q(1)
            de = m.scale(Q(1, 2), m.mm(m.mm(e, qi), h))
            require(dq(e, de) == h, 'explicit right inverse of coframe variation')
        kernel = []
        for a, b in PAIRS:
            lorentz = m.zero(4)
            lorentz[a][b], lorentz[b][a] = Q(ETA[a]), Q(-ETA[b])
            de = m.mm(lorentz, e)
            require(dq(e, de) == m.zero(4), 'Lorentz coframe kernel')
            kernel.append(sum(de, []))
        require(m.rank(kernel) == 6, 'six independent Lorentz gauge directions')
        require(m.rank(torsion_map(e)) == 24, 'connection equation removes all torsion')
        rows.append(dict(coframe_determinant=str(m.det(e)), metric_variation_rank=10,
                         lorentz_kernel_dimension=6, torsion_equation_rank=24))
    degenerate = m.eye(4)
    degenerate[0][0] = Q(0)
    degenerate_torsion_rank = m.rank(torsion_map(degenerate))
    require(m.rank(metric(degenerate)) == 3 and degenerate_torsion_rank < 24,
            'degenerate coframe refuses torsion conclusion')
    return dict(nondegenerate_fixtures=rows, right_inverse_controls=30,
                degenerate_metric_rank=3, degenerate_torsion_rank=degenerate_torsion_rank)


def riemann_fixtures():
    eta = lambda a, b: Q(ETA[a]) if a == b else Q(0)
    s = [[Q(a == b)*(a+1)+Q(a+b+1, 7) for b in IX] for a in IX]
    schouten = {(a, b, c, d): eta(a, c)*s[b][d]-eta(a, d)*s[b][c]
                -eta(b, c)*s[a][d]+eta(b, d)*s[a][c] for a, b, c, d in product(IX, repeat=4)}
    constant = {(a, b, c, d): Q(2, 3)*(eta(a, c)*eta(b, d)-eta(a, d)*eta(b, c))
                for a, b, c, d in product(IX, repeat=4)}
    weyl = {indices: Q(0) for indices in product(IX, repeat=4)}
    for (a, b), value in zip(PAIRS, (-2, 1, 1, -1, -1, 2)):
        for c, d, s1 in [(a, b, 1), (b, a, -1)]:
            for f, g, s2 in [(a, b, 1), (b, a, -1)]:
                weyl[c, d, f, g] = Q(s1*s2*value)
    return [('constant_curvature', constant, Q(2)), ('non_einstein', schouten, Q(-1, 3)),
            ('ricci_flat_weyl', weyl, Q(0))]


def curvature_forms(e, riemann):
    ef = coframe_forms(e)
    return [fadd(*(fscale(ETA[a]*ETA[b]*riemann[a, b, c, d], wedge(ef[c], ef[d]))
                   for c, d in PAIRS)) for a, b in PAIRS]


def gravity_density(e, spin_curvature, cosmological, holst, kappa=Q(5, 3)):
    forms = dict(E=matrix_form(coframe_forms(e), GAMMA), F=spin_curvature)
    return -(native_word('EEF', forms, J)+holst*native_word('EEF', forms))/kappa \
        +cosmological*native_word('EEEE', forms, J)/(24*kappa)


def einstein_variation_checks():
    count = holst_count = 0
    witnesses = []
    bivectors = [m.scale(Q(1, 2), m.mm(GAMMA[a], GAMMA[b])) for a, b in PAIRS]
    for name, r, cosmological in riemann_fixtures():
        for a, b, c, d in product(IX, repeat=4):
            require(r[a, b, c, d] == -r[b, a, c, d] == -r[a, b, d, c] == r[c, d, a, b],
                    'algebraic Riemann symmetries')
            require(r[a, b, c, d]+r[a, c, d, b]+r[a, d, b, c] == 0, 'algebraic first Bianchi')
        ric = [[sum(ETA[a]*r[a, b, a, d] for a in IX) for d in IX] for b in IX]
        scalar = sum(ETA[a]*ric[a][a] for a in IX)
        internal_residual = m.add(ric, [[Q(ETA[a])*(cosmological-scalar/2) if a == b else Q(0)
                                       for b in IX] for a in IX])
        if name == 'ricci_flat_weyl':
            curvature_sq = sum(ETA[a]*ETA[b]*ETA[c]*ETA[d]*v*v for (a, b, c, d), v in r.items())
            require(ric == m.zero(4) and curvature_sq == 48, 'nonflat algebraic vacuum curvature')
        require((internal_residual == m.zero(4)) == (name != 'non_einstein'), 'Einstein fixture labels')
        witnesses.append(dict(name=name, scalar=str(scalar), einstein=internal_residual == m.zero(4)))
        for e in coframe_fixtures():
            volume = m.det(e)
            qi = m.inverse(metric(e))
            rc = m.mm(m.mm(m.transpose(e), internal_residual), e)
            raised = m.mm(m.mm(qi, rc), qi)
            expected = m.scale(-volume/Q(5, 3), m.mm(m.mm([[Q(ETA[a]) if a == b else Q(0)
                                                          for b in IX] for a in IX], e), raised))
            spin = matrix_form(curvature_forms(e, r), bivectors)
            for holst in (Q(0), Q(2, 3), Q(-5)):
                require(gravity_density(e, spin, cosmological, holst) == volume*(scalar-2*cosmological)/(2*Q(5, 3)),
                        'native density equals Einstein-Hilbert on torsion-free curvature')
                for a, mu in product(IX, repeat=2):
                    ep, em = [row[:] for row in e], [row[:] for row in e]
                    ep[a][mu] += Q(1, 7)
                    em[a][mu] -= Q(1, 7)
                    # Each coframe entry appears at most linearly in the four-form.
                    actual = (gravity_density(ep, spin, cosmological, holst)
                              -gravity_density(em, spin, cosmological, holst))/Q(2, 7)
                    require(actual == expected[a][mu], 'independent native-action variation versus Einstein tensor')
                    hp = native_word('EEF', dict(E=matrix_form(coframe_forms(ep), GAMMA), F=spin))
                    hm = native_word('EEF', dict(E=matrix_form(coframe_forms(em), GAMMA), F=spin))
                    require(hp-hm == 0, 'Holst coframe variation vanishes on first Bianchi')
                    count += 1
                    holst_count += 1
    return dict(independent_coframe_variations=count, holst_bianchi_variations=holst_count,
                algebraic_curvature_witnesses=witnesses,
                pure_holst_does_not_impose_einstein=True,
                nonflat_weyl_scope='Algebraic pointwise curvature witness, not a constructed global solution')


def determinant_generic(a):
    n = len(a)
    return sum(sign(p)*product_scalar(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))


def product_scalar(items):
    result = Q(1)
    for value in items:
        result *= value
    return result


def inverse_generic(a):
    n, determinant = len(a), determinant_generic(a)
    return [[(-1)**(i+k)*determinant_generic([[a[r][c] for c in range(n) if c != i]
                                             for r in range(n) if r != k])/determinant
             for k in range(n)] for i in range(n)]


def maxwell_stress_checks():
    count = 0
    coupling_sq = Q(7, 5)
    field = [[Q(i-k, 3)+Q(i*k*(i-k), 11) for k in IX] for i in IX]
    for e in coframe_fixtures():
        q, w = metric(e), m.det(e)
        qi = m.inverse(q)
        up = m.mm(m.mm(qi, field), qi)
        field_sq = sum(field[i][k]*up[i][k] for i, k in product(IX, repeat=2))
        stress = m.scale(1/coupling_sq, m.add(m.mm(m.mm(field, qi), m.transpose(field)), m.scale(-field_sq/4, q)))
        require(sum(qi[i][k]*stress[i][k] for i, k in product(IX, repeat=2)) == 0, 'Maxwell stress trace')
        expected = m.scale(w, m.mm(m.mm([[Q(ETA[a]) if a == b else Q(0) for b in IX] for a in IX], e),
                                    m.mm(m.mm(qi, stress), qi)))
        for a, mu in product(IX, repeat=2):
            ej = [[Jet(e[b][nu], (Q((b, nu) == (a, mu)), Q(0), Q(0), Q(0))) for nu in IX] for b in IX]
            qj = metric(ej)
            inverse = inverse_generic(qj)
            root = square_root(-determinant_generic(qj))
            raised = m.mm(m.mm(inverse, field), inverse)
            density = -root*sum(field[i][k]*raised[i][k] for i, k in product(IX, repeat=2))/(4*coupling_sq)
            require(density.grad[0] == expected[a][mu], 'dual-number Maxwell variation versus stress tensor')
            count += 1
    return dict(independent_metric_action_variations=count, trace_free_stress_fixtures=3,
                coupling_square=str(coupling_sq), coupling_is_calibrated=True)


def thermo_coframe_checks():
    e = coframe_fixtures()[2]
    x = [Q(-1, 50), Q(0), Q(1, 50), Q(1, 25)]
    derivatives, inverse_count = [], 0
    for a, mu in product(IX, repeat=2):
        p, s = e[a][mu]-Q(a == mu), x[mu]
        require(p < Q(1, 6), 'Darboux branch')
        v = ((1+s)**2+1/(4*(p-Q(1, 6))**2))/5-2
        sv = Jet(s, (Q(0),)*4)
        vv = Jet(v, (Q(1), Q(0), Q(0), Q(0)))
        delta = 5*(2+vv)-(1+sv)*(1+sv)
        reconstructed = Q(1, 6)-1/(2*square_root(delta))
        require(reconstructed.value == p and reconstructed.grad[0] > 0, 'thermo coframe inverse and derivative')
        require(2+v > 0 and delta.value > 0, 'strict stable response')
        require(10+2*s+v+s*v > 0 and 1-s-5*v-s*s/2 > 0,
                'positive reference-patch temperature and pressure')
        derivatives.append(reconstructed.grad[0])
        inverse_count += 1
    jac = [[derivatives[i] if i == k else Q(0) for k in range(16)] for i in range(16)]
    require(m.rank(jac) == 16, 'sixteen thermo amplitude variations span all coframe variations')
    symmetric = [(a, b) for a in IX for b in range(a, 4)]
    cols = []
    for index, (a, mu) in enumerate(product(IX, repeat=2)):
        de = m.zero(4)
        de[a][mu] = derivatives[index]
        h = dq(e, de)
        cols.append([h[i][k] for i, k in symmetric])
    require(m.rank(m.transpose(cols)) == 10, 'thermo state variations span all metric variations')
    return dict(explicit_stable_channel_inverses=inverse_count, thermo_to_coframe_rank=16,
                thermo_to_metric_rank=10, number_of_channels_claimed_minimal=False)


def restricted_variation_controls():
    metric_poly = [[{} for _ in IX] for _ in IX]
    for a in IX:
        metric_poly[a][a] = {(2, 0, 0, 0): ETA[a]}
    results = []
    for t in (Q(1), Q(2), Q(3, 2)):
        q, ric, scalar = m.ricci(metric_poly, [t, Q(0), Q(0), Q(0)])
        require(scalar == 0 and ric == [[Q(3 if a == 0 else 1)/(t*t) if a == b else Q(0)
                                        for b in IX] for a in IX], 'direct conformal Ricci computation')
        require(ric != m.zero(4), 'trace equation is strictly weaker than vacuum Einstein')
        results.append(dict(t=str(t), scalar=str(scalar), ricci_diagonal=[str(ric[a][a]) for a in IX]))
    # A nonlinear coordinate pullback of Minkowski remains exactly flat.
    flat = [[{} for _ in IX] for _ in IX]
    flat[0][0] = {(0, 0, 0, 0): Q(-1)}
    flat[2][2] = flat[3][3] = {(0, 0, 0, 0): Q(1)}
    # Use X^1=(x^1)^2, hence the power belongs to coordinate 1.
    flat[1][1] = {(0, 2, 0, 0): Q(4)}
    _, ric, scalar = m.ricci(flat, [Q(0), Q(1), Q(0), Q(0)])
    require(ric == m.zero(4) and scalar == 0, 'coordinate-only coframe control')
    return dict(conformal_trace_only_counterexamples=results, nonlinear_coordinate_pullback_ricci_zero=True)


def run():
    return dict(trace_class=trace_and_action_class_checks(), nieh_yan=nieh_yan_checks(),
                variation_ranks=variation_rank_checks(), einstein=einstein_variation_checks(),
                maxwell_stress=maxwell_stress_checks(), thermo_coframe=thermo_coframe_checks(),
                restricted_variations=restricted_variation_controls())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
