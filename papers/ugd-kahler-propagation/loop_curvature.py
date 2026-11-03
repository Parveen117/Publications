"""LR exact loop/Cartan controls; no copied native arithmetic engine."""
from fractions import Fraction as Q
from itertools import combinations, product

import model as m
import metric_dynamics as mg
from thermo_gauge import Jet

IX, PAIRS, TOP = range(4), mg.PAIRS, mg.TOP
ETA, GAMMA, J = mg.ETA, mg.GAMMA, mg.J
require = mg.require


def mfadd(*forms):
    result = {}
    for form in forms:
        for key, value in form.items():
            result[key] = m.add(result[key], value) if key in result else value
    return result


def mfscale(c, form):
    return {key: m.scale(c, value) for key, value in form.items()}


def mform(forms, basis):
    return mfadd(*({key: m.scale(value, coefficient) for key, value in form.items()}
                   for form, coefficient in zip(forms, basis)))


def matrix_values(matrix):
    return [[x.value if isinstance(x, Jet) else x for x in row] for row in matrix]


def mfvalues(form):
    return {key: matrix_values(value) for key, value in form.items()}


def mfd(form):
    result = {}
    for key, matrix in form.items():
        for i in IX:
            s = mg.sign((i,)+key)
            if s:
                derivative = [[x.grad[i] if isinstance(x, Jet) else Q(0) for x in row] for row in matrix]
                result = mfadd(result, {tuple(sorted((i,)+key)): m.scale(s, derivative)})
    return result


def mf_equal(a, b):
    difference = mfvalues(mfadd(a, mfscale(-1, b)))
    return all(x == 0 for matrix in difference.values() for row in matrix for x in row)


def trace_top(form, insertion):
    n = len(insertion)
    return m.trace(m.mm(insertion, form.get(TOP, m.zero(n))))/n


def pair_trace(a, b, insertion):
    return trace_top(mg.wedge(a, b, matrix=True), insertion)


def representations():
    native, gamma, _, j = m.fixture()
    ident, r, k, _ = native
    bivectors = [m.mm(gamma[a], gamma[b]) for a, b in PAIRS]
    return [dict(name='original_real_four', sigma=Q(1), vectors=gamma,
                 bivectors=bivectors, j=j),
            dict(name='extra_native_cut', sigma=Q(1), vectors=[m.kron(k, g) for g in gamma],
                 bivectors=[m.kron(ident, b) for b in bivectors], j=m.kron(ident, j)),
            dict(name='extra_native_return', sigma=Q(-1), vectors=[m.kron(r, g) for g in gamma],
                 bivectors=[m.kron(ident, b) for b in bivectors], j=m.kron(ident, j))]


def coefficient_rank_checks():
    """Solve the 64-unknown Lorentz-vector intertwiner problem directly."""
    bivectors = [m.mm(GAMMA[a], GAMMA[b]) for a, b in PAIRS]
    units = []
    for i, k in product(IX, repeat=2):
        unit = m.zero(4)
        unit[i][k] = Q(1)
        units.append(unit)
    columns = []
    for component, unit in product(IX, units):
        y = [unit if a == component else m.zero(4) for a in IX]
        residual = []
        for (a, b), generator in zip(PAIRS, bivectors):
            for c in IX:
                expected = m.add(m.scale(2*ETA[b] if b == c else 0, y[a]),
                                 m.scale(-2*ETA[a] if a == c else 0, y[b]))
                residual += sum(m.add(m.comm(generator, y[c]), m.scale(-1, expected)), [])
        columns.append(residual)
    rank = m.rank(m.transpose(columns))
    require(rank == 62, 'Lorentz-vector intertwiners have dimension two')
    candidates = [sum((sum(g, []) for g in GAMMA), []),
                  sum((sum(m.mm(J, g), []) for g in GAMMA), [])]
    require(m.rank(candidates) == 2, 'Gamma and J Gamma are independent')
    constraints = m.transpose(columns)
    require(all(x == 0 for candidate in candidates for x in m.mv(constraints, candidate)),
            'both candidates satisfy every covariance constraint')
    phase_checks = 0
    for u, v in [(Q(1), Q(0)), (Q(0), Q(1)), (Q(3, 5), Q(4, 5)), (Q(2), Q(-3))]:
        y = [m.add(m.scale(u, g), m.scale(v, m.mm(J, g))) for g in GAMMA]
        for a, b in product(IX, repeat=2):
            require(m.mm(y[a], y[b]) == m.scale(u*u+v*v, m.mm(GAMMA[a], GAMMA[b])),
                    'real four-component embedding has nonnegative square factor')
            phase_checks += 1
    dimension_controls = []
    for rep in representations():
        vectors, bv, sigma, jj = rep['vectors'], rep['bivectors'], rep['sigma'], rep['j']
        n = len(jj)
        require(m.rank([sum(x, []) for x in bv+vectors]) == 10, 'Cartan Lie algebra dimension ten')
        for (a, b), generator in zip(PAIRS, bv):
            require(m.comm(vectors[a], vectors[b]) == m.scale(2*sigma, generator), 'native square fixes bracket sign')
            for c in IX:
                expected = m.add(m.scale(2*ETA[b] if b == c else 0, vectors[a]),
                                 m.scale(-2*ETA[a] if a == c else 0, vectors[b]))
                require(m.comm(generator, vectors[c]) == expected, 'extended Lorentz covariance')
        require(m.mm(jj, jj) == m.scale(-1, m.eye(n)), 'retained volume square')
        dimension_controls.append(dict(name=rep['name'], real_dimension=n, square_sign=int(sigma), lie_dimension=10))
    return dict(intertwiner_unknowns=64, intertwiner_constraint_rank=rank, intertwiner_dimension=2,
                real_embedding_product_checks=phase_checks, representations=dimension_controls,
                minimality_scope='Two copies suffice and one copy fails, within copies of the fixed real-four module')


def jet_fixture(seed):
    ef = [{(mu,): Jet(Q(a == mu)+Q(seed+a-mu, 31),
                     tuple(Q(seed+(a+1)*(mu+1)-k, 37) for k in IX)) for mu in IX} for a in IX]
    omega = [[{} for _ in IX] for _ in IX]
    for a, b in PAIRS:
        omega[a][b] = {(mu,): Jet(Q(seed+a+b-mu, 41),
                                 tuple(Q((a+1)*(k+1)-b+mu+seed, 43) for k in IX)) for mu in IX}
        omega[b][a] = mg.fscale(-1, omega[a][b])
    tf = [mg.fadd(mg.exterior_d(ef[a]), *(mg.fscale(ETA[b], mg.wedge(omega[a][b], ef[b])) for b in IX))
          for a in IX]
    rf = [[mg.fadd(mg.exterior_d(omega[a][b]),
                   *(mg.fscale(ETA[c], mg.wedge(omega[a][c], omega[c][b])) for c in IX))
           for b in IX] for a in IX]
    return ef, omega, tf, rf


def lift_geometry(rep, ef, omega, tf, rf, scale):
    bv, y = rep['bivectors'], rep['vectors']
    om = mform([omega[a][b] for a, b in PAIRS], [m.scale(Q(1, 2), b) for b in bv])
    ey = mform(ef, y)
    conn = mfadd(om, mfscale(scale, ey))
    direct = mfvalues(mfadd(mfd(conn), mg.wedge(conn, conn, matrix=True)))
    curvature = mfvalues(mform([rf[a][b] for a, b in PAIRS], [m.scale(Q(1, 2), b) for b in bv]))
    ee = mfvalues(mform([mg.wedge(ef[a], ef[b]) for a, b in PAIRS], [m.scale(2, b) for b in bv]))
    ty = mfvalues(mform(tf, y))
    return conn, direct, curvature, ee, ty


def cartan_and_readout_checks():
    counts = seed_counts = 0
    for rep, seed, scale in product(representations(), (1, 2), (Q(1, 4), Q(2, 3))):
        ef, omega, tf, rf = jet_fixture(seed)
        conn, direct, curvature, ee, ty = lift_geometry(rep, ef, omega, tf, rf, scale)
        jj, sigma = rep['j'], rep['sigma']
        ident = m.eye(len(jj))
        rho = sigma*scale*scale
        require(mf_equal(direct, mfadd(curvature, mfscale(scale, ty), mfscale(rho, ee))),
                'direct dA+A^2 versus independent curvature/torsion split')
        tt = sum(ETA[a]*mg.wedge(mg.values(tf[a]), mg.values(tf[a])).get(TOP, Q(0)) for a in IX)
        ny = mg.exterior_d(mg.fadd(*(mg.fscale(ETA[a], mg.wedge(ef[a], tf[a])) for a in IX))).get(TOP, Q(0))
        p = pair_trace(ee, curvature, jj)
        h = pair_trace(ee, curvature, ident)
        v = pair_trace(ee, ee, jj)
        top_i = pair_trace(curvature, curvature, ident)
        top_j = pair_trace(curvature, curvature, jj)
        require(pair_trace(direct, direct, jj) == top_j+2*rho*p+rho*rho*v, 'oriented return action expansion')
        require(pair_trace(direct, direct, ident) == top_i+rho*ny, 'full unoriented trace is characteristic plus Nieh-Yan')
        require(tt+2*h == ny, 'independent boundary form normalization')
        # J conjugation is the admitted even/odd cut of the curvature record.
        jf_j = {key: m.mm(m.mm(jj, value), jj) for key, value in direct.items()}
        even = mfscale(Q(1, 2), mfadd(direct, mfscale(-1, jf_j)))
        odd = mfscale(Q(1, 2), mfadd(direct, jf_j))
        require(mf_equal(even, mfadd(curvature, mfscale(rho, ee))), 'cut retains even curvature')
        require(mf_equal(odd, mfscale(scale, ty)), 'cut isolates torsion component')
        for aa, bb, cc in [(Q(2), Q(-3), Q(2)), (Q(2), Q(-3), Q(0)), (Q(-1), Q(4), Q(5))]:
            actual = aa*pair_trace(even, even, ident)+bb*pair_trace(even, even, jj)+cc*pair_trace(odd, odd, ident)
            bulk = 2*rho*((aa-cc)*h+bb*p)+bb*rho*rho*v
            require(actual == aa*top_i+bb*top_j+cc*rho*ny+bulk, 'projected return family including boundary terms')
            kappa, cosmological, gamma = -1/(2*rho*bb), -12*rho, (aa-cc)/bb
            require(bulk == -(p+gamma*h)/kappa+cosmological*v/(24*kappa), 'MG coefficients recovered exactly')
            require(kappa*cosmological == 6/bb, 'coupling product')
            seed_counts += 1
        counts += 1
    return dict(independent_connection_jets=counts, three_weight_readout_reductions=seed_counts,
                full_trace_has_zero_bulk_holst=True, even_cut_restores_holst_freedom=True)


def series_multiply(a, b, order):
    n = len(a[0])
    return [sum_matrices([m.mm(a[i], b[k-i]) for i in range(k+1)], n) for k in range(order+1)]


def sum_matrices(items, n):
    answer = m.zero(n)
    for item in items:
        answer = m.add(answer, item)
    return answer


def edge_series(a, b, c, order):
    """Exact epsilon series of U'=U[epsilon*a+epsilon^2*(b+t*c)], U(0)=I."""
    n = len(a)
    polynomials = [[m.eye(n)]]
    for k in range(1, order+1):
        derivative = []
        for power, coefficient in enumerate(polynomials[k-1]):
            while len(derivative) <= power:
                derivative.append(m.zero(n))
            derivative[power] = m.add(derivative[power], m.mm(coefficient, a))
        if k >= 2:
            for power, coefficient in enumerate(polynomials[k-2]):
                while len(derivative) <= power+1:
                    derivative.append(m.zero(n))
                derivative[power] = m.add(derivative[power], m.mm(coefficient, b))
                derivative[power+1] = m.add(derivative[power+1], m.mm(coefficient, c))
        polynomials.append([m.zero(n)]+[m.scale(Q(1, power+1), coefficient)
                                        for power, coefficient in enumerate(derivative)])
    return [sum_matrices(polynomial, n) for polynomial in polynomials]


def rectangle_series(connection, mu, nu, order=4):
    coefficients = [connection[(i,)] for i in IX]
    n = len(coefficients[0])
    coordinates = [0]*4
    answer = [m.eye(n)]+[m.zero(n) for _ in range(order)]
    for direction, step in [(mu, 1), (nu, 1), (mu, -1), (nu, -1)]:
        matrix = coefficients[direction]
        a = m.scale(step, matrix_values(matrix))
        b = [[step*sum(coordinates[i]*x.grad[i] for i in IX) for x in row] for row in matrix]
        c = [[x.grad[direction] for x in row] for row in matrix]  # step^2=1
        answer = series_multiply(answer, edge_series(a, b, c, order), order)
        coordinates[direction] += step
    require(coordinates == [0]*4, 'based rectangle closes')
    return answer


def holonomy_checks():
    loop_count = products_checked = 0
    for rep in [representations()[0], representations()[2]]:
        ef, omega, tf, rf = jet_fixture(1)
        connection, curvature, _, _, _ = lift_geometry(rep, ef, omega, tf, rf, Q(1, 4))
        jj, n = rep['j'], len(rep['j'])
        loops = {}
        for mu, nu in PAIRS:
            coefficients = rectangle_series(connection, mu, nu)
            require(coefficients[0] == m.eye(n) and coefficients[1] == m.zero(n), 'loop starts at second order')
            require(coefficients[2] == curvature[mu, nu], 'affine-connection holonomy reconstructs full curvature')
            reverse = rectangle_series(connection, nu, mu)
            identity = series_multiply(coefficients, reverse, 4)
            require(identity == [m.eye(n)]+[m.zero(n) for _ in range(4)], 'reverse path is formal inverse through order four')
            coefficients[0] = m.zero(n)
            loops[mu, nu] = coefficients
            loop_count += 1
        for insertion in [m.eye(n), jj]:
            readout = [m.zero(n) for _ in range(5)]
            for left, right, sign in [((0, 1), (2, 3), 1), ((0, 2), (1, 3), -1), ((0, 3), (1, 2), 1)]:
                lr = series_multiply(loops[left], loops[right], 4)
                rl = series_multiply(loops[right], loops[left], 4)
                readout = [m.add(x, m.scale(sign, m.add(y, z))) for x, y, z in zip(readout, lr, rl)]
            require(all(x == m.zero(n) for x in readout[:4]), 'quadratic return starts at four-volume order')
            require(m.trace(m.mm(insertion, readout[4]))/n == pair_trace(curvature, curvature, insertion),
                    'loop readout coefficient equals curvature-square density')
            products_checked += 1
    return dict(based_affine_connection_rectangles=loop_count, inverse_series_order=4,
                independent_four_volume_readouts=products_checked,
                approximation_scope='Exact Taylor coefficients; continuum convergence is proved in LR-1, not sampled')


def invariant_form_checks():
    records = []
    symmetric = list(combinations(range(10), 2))+[(i, i) for i in range(10)]
    index = {tuple(sorted(pair)): k for k, pair in enumerate(symmetric)}
    for rep in [representations()[0], representations()[2]]:
        basis = rep['bivectors']+rep['vectors']
        n = len(rep['j'])
        gram = [[m.trace(m.mm(x, y))/n for y in basis] for x in basis]
        require(all(gram[i][k] == 0 for i, k in product(range(10), repeat=2) if i != k), 'orthogonal Lie basis')
        require(all(gram[i][i] != 0 for i in range(10)), 'nondegenerate Lie trace pairing')
        constraints = []
        lorentz_rank = None
        for gindex, generator in enumerate(basis):
            adj = m.transpose([[m.trace(m.mm(x, m.comm(generator, y)))/(n*gram[i][i])
                                for i, x in enumerate(basis)] for y in basis])
            for i, k in symmetric:
                row = [Q(0)]*55
                for j in range(10):
                    row[index[tuple(sorted((j, k)))]] += adj[j][i]
                    row[index[tuple(sorted((i, j)))]] += adj[j][k]
                constraints.append(row)
            if gindex == 5:
                lorentz_rank = m.rank(constraints)
                require(lorentz_rank == 52, 'three Lorentz-invariant symmetric bilinear forms')
        full_rank = m.rank(constraints)
        require(full_rank == 54, 'one fully Cartan-invariant symmetric bilinear form')
        records.append(dict(name=rep['name'], symmetric_unknowns=55, lorentz_constraint_rank=lorentz_rank,
                            lorentz_invariant_dimension=3, full_cartan_constraint_rank=full_rank,
                            full_cartan_invariant_dimension=1))
    return dict(representations=records)


def covariant_d(connection, two_form):
    return mfadd(mfd(two_form), mg.wedge(connection, two_form, matrix=True),
                 mfscale(-1, mg.wedge(two_form, connection, matrix=True)))


def cut(form, jj, even=True):
    moved = {key: m.mm(m.mm(jj, value), jj) for key, value in form.items()}
    return mfscale(Q(1, 2), mfadd(form, mfscale(-1 if even else 1, moved)))


def left_j(form, jj):
    return {key: m.mm(jj, value) for key, value in form.items()}


def compact_equation_checks():
    count = 0
    for rep, scale in product([representations()[0], representations()[2]], (Q(1, 4), Q(2, 3))):
        ef, omega, _, _ = jet_fixture(2)
        bv, jj, rho = rep['bivectors'], rep['j'], rep['sigma']*scale*scale
        om = mform([omega[a][b] for a, b in PAIRS], [m.scale(Q(1, 2), b) for b in bv])
        translation = mfscale(scale, mform(ef, rep['vectors']))
        connection = mfadd(om, translation)
        curvature = mfadd(mfd(connection), mg.wedge(connection, connection, matrix=True))
        require(mf_equal(covariant_d(connection, curvature), {}), 'independent affine-connection Bianchi identity')
        even, odd = cut(curvature, jj), cut(curvature, jj, even=False)
        ee = mform([mg.wedge(ef[a], ef[b]) for a, b in PAIRS], [m.scale(2, b) for b in bv])
        dee = covariant_d(om, ee)
        for aa, bb, cc in [(Q(2), Q(-3), Q(2)), (Q(2), Q(-3), Q(0)), (Q(-1), Q(4), Q(5))]:
            response = mfadd(mfscale(aa, even), mfscale(bb, left_j(even, jj)), mfscale(cc, odd))
            equation = covariant_d(connection, response)
            expected_even = mfscale(rho, mfadd(mfscale(aa-cc, dee), mfscale(bb, left_j(dee, jj))))
            require(mf_equal(cut(equation, jj), expected_even), 'compact equation even part is torsion equation')
            active_even = mfadd(mfscale(aa-cc, even), mfscale(bb, left_j(even, jj)))
            expected_odd = mfadd(mg.wedge(translation, active_even, matrix=True),
                                 mfscale(-1, mg.wedge(active_even, translation, matrix=True)))
            require(mf_equal(cut(equation, jj, even=False), expected_odd), 'compact equation odd part is coframe equation')
            if aa == cc:
                anti_j = mfscale(2, left_j(even, jj))
                require(mf_equal(equation, mfscale(bb/2, covariant_d(connection, anti_j))),
                        'unprojected native stationary law is D{J,F}=0')
            count += 1
    return dict(independent_compact_equation_splits=count, bianchi_connection_jets=4,
                includes_full_trace_and_projected_readouts=True)


def seed_variation_checks():
    count = 0
    e = mg.coframe_fixtures()[2]
    scale, kappa = Q(1, 4), Q(5, 3)
    for rep in [representations()[0], representations()[2]]:
        jj, bv, sigma = rep['j'], rep['bivectors'], rep['sigma']
        rho = sigma*scale*scale
        bb, cosmological = -1/(2*kappa*rho), -12*rho
        for _, riemann, _ in mg.riemann_fixtures():
            curvature = mform(mg.curvature_forms(e, riemann), [m.scale(Q(1, 2), b) for b in bv])
            reference_curvature = mg.matrix_form(mg.curvature_forms(e, riemann),
                                                [m.scale(Q(1, 2), m.mm(GAMMA[a], GAMMA[b])) for a, b in PAIRS])

            def density(frame):
                ef = mg.coframe_forms(frame)
                ee = mform([mg.wedge(ef[a], ef[b]) for a, b in PAIRS], [m.scale(2, b) for b in bv])
                full = mfadd(curvature, mfscale(rho, ee))
                return bb*pair_trace(full, full, jj)

            for a, mu in product(IX, repeat=2):
                ep, em = [row[:] for row in e], [row[:] for row in e]
                ep[a][mu] += Q(1, 11)
                em[a][mu] -= Q(1, 11)
                direct = (density(ep)-density(em))/Q(2, 11)
                expected = (mg.gravity_density(ep, reference_curvature, cosmological, Q(0), kappa)
                            -mg.gravity_density(em, reference_curvature, cosmological, Q(0), kappa))/Q(2, 11)
                require(direct == expected, 'independent Cartan-square variation equals source-pinned MG variation')
                count += 1
    return dict(independent_seed_coframe_variations=count, compared_to='Unchanged MG Einstein-tensor-checked action')


def curved_flat_cartan_checks():
    records = []
    for rep, direction in [(representations()[0], 3), (representations()[2], 0)]:
        sigma = rep['sigma']
        require(sigma == ETA[direction], 'spacelike versus timelike conformal coordinate')
        for ell, coordinate in product((Q(1), Q(2)), (Q(1), Q(2, 3))):
            f = Jet(ell/coordinate, tuple(-ell/(coordinate*coordinate) if i == direction else Q(0) for i in IX))
            ef = [{(a,): f} for a in IX]
            omega = [[mg.fscale(sigma/ell, mg.fadd(mg.fscale(Q(a == direction), ef[b]),
                                                  mg.fscale(-Q(b == direction), ef[a]))) for b in IX] for a in IX]
            tf = [mg.fadd(mg.exterior_d(ef[a]), *(mg.fscale(ETA[b], mg.wedge(omega[a][b], ef[b])) for b in IX))
                  for a in IX]
            rf = [[mg.fadd(mg.exterior_d(omega[a][b]),
                           *(mg.fscale(ETA[c], mg.wedge(omega[a][c], omega[c][b])) for c in IX))
                   for b in IX] for a in IX]
            require(all(mg.equal_forms(t, {}) for t in tf), 'explicit conformal coframe torsion zero')
            for a, b in PAIRS:
                require(mg.equal_forms(rf[a][b], mg.fscale(-sigma/(ell*ell), mg.wedge(ef[a], ef[b]))),
                        'direct spin curvature is nonzero constant sectional curvature')
            _, full, curvature, _, _ = lift_geometry(rep, ef, omega, tf, rf, 1/(2*ell))
            require(mf_equal(full, {}) and not mf_equal(curvature, {}), 'flat Cartan transport with curved Levi-Civita geometry')
            polys = [[{} for _ in IX] for _ in IX]
            powers = tuple(-2 if i == direction else 0 for i in IX)
            for a in IX:
                polys[a][a] = {powers: ETA[a]*ell*ell}
            point = [Q(1)]*4
            point[direction] = coordinate
            q, ric, scalar = m.ricci(polys, point)
            cosmological = -3*sigma/(ell*ell)
            require(ric == m.scale(cosmological, q) and scalar == 4*cosmological, 'independent metric Ricci equals selected Lambda')
            records.append(dict(sector='de_Sitter' if sigma == -1 else 'anti_de_Sitter', ell=str(ell),
                                coordinate=str(coordinate), cosmological=str(cosmological), scalar=str(scalar),
                                cartan_curvature_zero=True, metric_curvature_nonzero=True))
    # The seed admits more than the flat-Cartan maximally symmetric subclass.
    weyl = mg.riemann_fixtures()[2][1]
    for rep in [representations()[0], representations()[2]]:
        rho = rep['sigma']/16
        eta = lambda a, b: Q(ETA[a]) if a == b else Q(0)
        r = {(a, b, c, d): weyl[a, b, c, d]-4*rho*(eta(a, c)*eta(b, d)-eta(a, d)*eta(b, c))
             for a, b, c, d in product(IX, repeat=4)}
        ric = [[sum(ETA[a]*r[a, b, a, d] for a in IX) for d in IX] for b in IX]
        require(ric == [[-12*rho*eta(a, b) for b in IX] for a in IX], 'Einstein algebraic curvature with nonzero Weyl part')
        ef = mg.coframe_forms(m.eye(4))
        curvature = mform(mg.curvature_forms(m.eye(4), r), [m.scale(Q(1, 2), b) for b in rep['bivectors']])
        ee = mform([mg.wedge(ef[a], ef[b]) for a, b in PAIRS], [m.scale(2, b) for b in rep['bivectors']])
        require(not mf_equal(mfadd(curvature, mfscale(rho, ee)), {}), 'Einstein does not require all Cartan returns trivial')
    return dict(explicit_metric_checks=records, nonzero_weyl_algebraic_controls=2,
                weyl_scope='Pointwise algebraic controls; the displayed conformal metrics are actual local solutions')


def selection_and_scale_controls():
    f = {(0, 1): m.mm(GAMMA[0], GAMMA[1]), (2, 3): GAMMA[2]}
    generator = GAMMA[3]
    df = {key: m.comm(generator, value) for key, value in f.items()}
    violation = pair_trace(df, f, J)+pair_trace(f, df, J)
    require(violation == 4, 'fixed native J breaks the enlarged Cartan symmetry off shell')
    require(pair_trace(df, f, m.eye(4))+pair_trace(f, df, m.eye(4)) == 0, 'ordinary full trace retains that symmetry')
    for a, b in PAIRS:
        df = {key: m.comm(m.mm(GAMMA[a], GAMMA[b]), value) for key, value in f.items()}
        require(pair_trace(df, f, J)+pair_trace(f, df, J) == 0, 'Lorentz symmetry survives the native cut readout')
    two_plane = {(0, 1): m.mm(GAMMA[0], GAMMA[1]), (2, 3): m.mm(GAMMA[2], GAMMA[3])}
    flipped = dict(two_plane)
    flipped[2, 3] = m.scale(-1, flipped[2, 3])
    densities = [pair_trace(curvature, curvature, J) for curvature in [two_plane, flipped]]
    norms = [sum(m.trace(m.mm(m.transpose(value), value))/4 for value in curvature.values())
             for curvature in [two_plane, flipped]]
    require(densities == [-2, 2] and norms == [2, 2], 'oriented loop action is not a positive defect norm')
    records = []
    rep = representations()[2]
    ef, omega, tf, rf = jet_fixture(1)
    scale, bb = Q(1, 4), Q(8)
    connection, _, _, _, _ = lift_geometry(rep, ef, omega, tf, rf, scale)
    rho = rep['sigma']*scale*scale
    kappa, cosmological = -1/(2*bb*rho), -12*rho
    for dilation in (Q(2), Q(3, 2), Q(1, 3)):
        new_ef = [mg.fscale(dilation, e) for e in ef]
        new_tf = [mg.fscale(dilation, t) for t in tf]
        rescaled, _, _, _, _ = lift_geometry(rep, new_ef, omega, new_tf, rf, scale/dilation)
        require(connection == rescaled, 'all local connection jets and hence returns unchanged by unit dilation')
        new_rho = rho/(dilation*dilation)
        new_kappa, new_cosmological = -1/(2*bb*new_rho), -12*new_rho
        require(new_kappa == dilation*dilation*kappa and new_cosmological == cosmological/(dilation*dilation),
                'dimensionful coefficient ambiguity under identical return records')
        require(new_kappa*new_cosmological == 6/bb, 'coefficient product survives unit dilation')
        records.append(dict(dilation=str(dilation), kappa=str(new_kappa), cosmological=str(new_cosmological),
                            unchanged_product=str(new_kappa*new_cosmological)))
    _, zero_curvature, curvature, _, _ = lift_geometry(rep, ef, omega, tf, rf, Q(0))
    require(mf_equal(zero_curvature, curvature), 'zero soldering scale removes coframe from the seed')
    return dict(nonlorentz_readout_variation=str(violation), oriented_density_sign_control=list(map(str, densities)),
                positive_norms=list(map(str, norms)), unit_dilations=records,
                finite_seed_zero_scale_is_characteristic_only=True,
                seed_and_readout_selection='DECLARED_NOT_PRIMITIVE_DERIVED')


def run():
    return dict(coefficients=coefficient_rank_checks(), invariant_forms=invariant_form_checks(),
                cartan_action=cartan_and_readout_checks(), holonomy=holonomy_checks(),
                compact_equation=compact_equation_checks(),
                action_variation=seed_variation_checks(), solutions=curved_flat_cartan_checks(),
                selection=selection_and_scale_controls())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
