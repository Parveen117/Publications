"""BI exact controls for general homogeneous shear; consumes unchanged SM/CS/PS."""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial

import model as m
import metric_dynamics as mg
import loop_curvature as lr
import spin_matter as sm
import coupled_cosmology as cs
import perturbation_stability as ps
from thermo_gauge import Jet

require = mg.require
IX, PAIRS, ETA = sm.IX, sm.PAIRS, sm.ETA
SPACE = range(3)
SPIN_FORMS = [[sm.HG3[0, i+1, j+1] for j in SPACE] for i in SPACE]


def spin_matrix(psi):
    return [[sm.bil(psi, SPIN_FORMS[i][j]) for j in SPACE] for i in SPACE]


def shear_square(shear):
    return m.trace(m.mm(shear, shear))/2


def spatial_connection(expansion):
    return [sm.lin([-expansion[i][j]/2 for i in SPACE], sm.A) for j in SPACE]


def universal_matrix_checks():
    require(sm.HG3[1, 2, 3] == m.scale(-1, ps.AXIAL), 'axial-time bilinear equals minus B')
    for i, j in combinations(SPACE, 2):
        matrix = SPIN_FORMS[i][j]
        require(m.transpose(matrix) == matrix, 'symmetric spin polarization form')
        require(m.comm(matrix, sm.PHASE) == m.zero(8), 'native phase preserves polarization')
        for generator in [sm.GAMMA[0], m.mm(sm.GAMMA[0], sm.VOLUME)]:
            require(m.comm(matrix, generator) == m.zero(8), 'all homogeneous spin components conserved after volume rescaling')
    count = 0
    for i, j in combinations(range(4), 2):
        if j == 3:
            expansion = m.zero(3)
            expansion[i][i] = Q(1)
        else:
            expansion = m.zero(3)
            expansion[i][j] = expansion[j][i] = Q(1)
        connection = spatial_connection(expansion)
        slash = sm.lin([Q(1)]*3, [m.mm(sm.GAMMA[k+1], connection[k]) for k in SPACE])
        require(slash == m.scale(m.trace(expansion)/2, sm.GAMMA[0]), 'all six symmetric expansion directions enter the matter equation only through trace')
        # Compare quadratic-form coefficients, rather than selected matter states.
        for a, b in product(SPACE, repeat=2):
            source = m.scale(Q(-1, 4), m.add(m.mm(sm.HG[a+1], connection[b]),
                                                    m.mm(sm.HG[b+1], connection[a])))
            expected = sm.lin([expansion[a][k] for k in SPACE]+[-expansion[k][b] for k in SPACE],
                              [SPIN_FORMS[k][b] for k in SPACE]+[SPIN_FORMS[a][k] for k in SPACE])
            expected = m.scale(Q(1, 8), expected)
            require(m.add(source, m.transpose(source)) == m.scale(2, expected), 'universal anisotropic stress coefficient +Z/8 [K,Q]')
            count += 1
    return dict(conserved_spin_commutators=6, trace_only_connection_directions=6,
                universal_spatial_stress_forms=count, axial_time_identity=True)


def bianchi_fixture(seed, omit_spin_torque=False):
    norm, kappa, mass, rate = Q(4, 3), Q(1), Q(1, 2), Q(2)
    raw = [Q(i == 0) for i in range(8)]
    psi = m.mv(m.add(m.eye(8), sm.H), raw)
    if seed == 2:
        psi = m.mv(m.add(m.eye(8), sm.H), [Q((i*i+2*i) % 5-2, 3) for i in range(8)])
    if seed == 3:
        psi = [x+Q((i*i+3*i) % 7-3, 20) for i, x in enumerate(psi)]
        require(sm.bil(psi, ps.PSEUDO) != 0, 'full-component nonrest fixture')
    shear = [[Q(1, 3), Q(seed, 11), Q(-1, 7)],
             [Q(seed, 11), Q(-1, 5), Q(2, 13)],
             [Q(-1, 7), Q(2, 13), Q(-2, 15)]]
    expansion = m.add(m.scale(rate, m.eye(3)), shear)
    scalar, pseudo = sm.bil(psi, sm.H), sm.bil(psi, ps.PSEUDO)
    energy = norm*mass*scalar/2+3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
    pressure = 3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
    spin = spin_matrix(psi)
    anisotropic = m.scale(norm/8, m.comm(shear, spin))
    require(any(x for row in anisotropic for x in row), 'nonzero spin-shear stress')
    rate_dot = -kappa*(energy+pressure)/2-shear_square(shear)
    cosmological = 3*rate*rate-kappa*energy-shear_square(shear)
    require(cosmological > 0, 'positive-Lambda Bianchi fixture')
    shear_dot = m.scale(-3*rate, shear)
    if not omit_spin_torque:
        shear_dot = m.add(shear_dot, m.scale(kappa, anisotropic))
    expansion_dot = m.add(m.scale(rate_dot, m.eye(3)), shear_dot)
    dp = [x-3*rate*y/2 for x, y in zip(ps.nonlinear_flow(psi, mass, 3*kappa*norm/16), psi)]
    psi_j = [Jet(x, (dx, Q(0), Q(0), Q(0))) for x, dx in zip(psi, dp)]
    lengths = [[Q(1), Q(1, 7), Q(-1, 11)], [Q(0), Q(3, 2), Q(2, 9)], [Q(0), Q(0), Q(2)]]
    lengths_dot = m.mm(expansion, lengths)
    e = [[Jet.lift(0) for _ in IX] for _ in IX]
    e[0][0] = Jet.lift(1)
    for i, j in product(SPACE, repeat=2):
        e[i+1][j+1] = Jet(lengths[i][j], (lengths_dot[i][j], Q(0), Q(0), Q(0)))
    kj = [[Jet(expansion[i][j], (expansion_dot[i][j], Q(0), Q(0), Q(0))) for j in SPACE] for i in SPACE]
    omega_lc = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    omega = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    ax = sm.axial(psi_j)
    for p, (a, b) in enumerate(PAIRS):
        for mu in IX:
            if a == 0:
                omega_lc[p][mu] = sum(kj[b-1][j]*e[j+1][mu] for j in SPACE)
            omega[p][mu] = omega_lc[p][mu]+kappa*norm*sum(ETA[c]*ax[a, b, c]*e[c][mu] for c in IX)/8
    return dict(e=e, ef=mg.coframe_forms(e), omega=omega, omega_lc=omega_lc, psi=psi,
                dpsi=[dp]+[[Q(0)]*8 for _ in SPACE], mass=mass, kappa=kappa, norm=norm,
                rate=rate, rate_dot=rate_dot, cosmological=cosmological, energy=energy,
                pressure=pressure, n=norm*scalar/2, shear=shear, expansion=expansion,
                expansion_dot=expansion_dot, anisotropic=anisotropic)


def direct_stress(f):
    e, psi, dp = lr.matrix_values(f['e']), f['psi'], f['dpsi']
    inv, norm = m.inverse(e), f['norm']
    gam = [sm.lin(inv[mu], sm.GAMMA) for mu in IX]
    metric = mg.metric(e)
    lower = [sm.lin(metric[mu], gam) for mu in IX]
    connection = sm.connections(lr.matrix_values(f['omega_lc']), [Q(0)]*4, Q(0))
    cov = [[x+y for x, y in zip(dp[mu], m.mv(connection[mu], psi))] for mu in IX]
    kinetic = norm*sum(sm.bil(psi, m.mm(sm.H, gam[mu]), cov[mu]) for mu in IX)/2
    lag = kinetic-f['energy']
    require(lag == f['pressure'], 'anisotropy does not alter the on-shell scalar pressure')
    return [[-norm*(sm.bil(psi, m.mm(sm.H, lower[mu]), cov[nu])
                   +sm.bil(psi, m.mm(sm.H, lower[nu]), cov[mu]))/4+metric[mu][nu]*lag
             for nu in IX] for mu in IX]


def original_equation_checks():
    for seed in (1, 2, 3):
        f = bianchi_fixture(seed)
        for name, values in cs.first_order_residuals(f).items():
            require(all(x == 0 for x in values), 'general Bianchi original '+name+' equations')
        e = lr.matrix_values(f['e'])
        orthogonal = [[Q(0)]*4 for _ in IX]
        orthogonal[0][0] = f['energy']
        for i, j in product(SPACE, repeat=2):
            orthogonal[i+1][j+1] = Q(i == j)*f['pressure']+f['anisotropic'][i][j]
        require(direct_stress(f) == m.mm(m.mm(m.transpose(e), orthogonal), e), 'all coordinate stress components, including nondiagonal spatial source')
        require(m.trace(f['anisotropic']) == 0
                and m.trace(m.mm(f['shear'], f['anisotropic'])) == 0, 'spin anisotropic stress does no shear work')
    bad = cs.first_order_residuals(bianchi_fixture(1, omit_spin_torque=True))
    require(any(bad['coframe']), 'omitting spin torque violates the coframe equations')
    require(not any(bad['matter']) and not any(bad['cartan']), 'negative control still solves matter and Cartan equations')
    return dict(solution_jets=3, coframe_equations=48, connection_equations=72,
                matter_equations=24, effective_stress_components=48,
                nonrest_general_metric_jet=True,
                omitted_torque_nonzero_coframe_components=sum(x != 0 for x in bad['coframe']))


def volume_law(omega, mass, kappa, n0, q, z):
    cosmological = omega*omega/3
    d, b = kappa*mass*n0/(2*cosmological), q/omega
    ch, sh = (z+1/z)/2, (z-1/z)/2
    v = d*(ch-1)+b*sh
    dv = omega*(d*sh+b*ch)
    ddv = omega*omega*(d*ch+b*sh)
    return v, dv, ddv


def volume_and_spin_clock_checks():
    count = 0
    for omega, mass, kappa, n0, ratio, z in product((Q(1), Q(3, 2)), (Q(0), Q(1, 2)),
                                                  (Q(1),), (Q(4, 3),), (Q(0), Q(3, 5)), (Q(2), Q(3))):
        c = 3*kappa*n0/4
        q = c*(1+ratio*ratio)/(1-ratio*ratio)
        shear_constant = (q*q-c*c)/3
        v, dv, ddv = volume_law(omega, mass, kappa, n0, q, z)
        cosmological, n = omega*omega/3, n0/v
        rho, pressure = mass*n+3*kappa*n*n/16, 3*kappa*n*n/16
        rate = dv/(3*v)
        rate_dot = ddv/(3*v)-dv*dv/(3*v*v)
        require(dv*dv == 3*cosmological*v*v+3*kappa*mass*n0*v+c*c+3*shear_constant,
                'anisotropic volume first integral')
        require(3*rate*rate == cosmological+kappa*rho+shear_constant/(v*v), 'anisotropic Hamiltonian constraint')
        require(rate_dot == -kappa*(rho+pressure)/2-shear_constant/(v*v), 'anisotropic mean acceleration')
        d, b = kappa*mass*n0/(2*cosmological), q/omega
        u_dot = omega*z*(1/(z-1)-(d+b)/((d+b)*z-(d-b)))/q
        require(u_dot == 1/v, 'closed inverse-volume clock primitive')
        if ratio == 0:
            require((v, dv, ddv) == cs.positive_volume(omega/6, mass, kappa, n0, z), 'zero-shear limit exactly recovers CS')
        require((3*kappa*n0/8)/(kappa*n0/4) == Q(3, 2), 'contact phase to unwrapped spin-shear rotation ratio')
        count += 1
    return dict(volume_constraint_and_clock_checks=count, zero_shear_CS_limit=True,
                contact_phase_rotation_ratio='3/2', scope='Ratio within the declared rest-sector action, not quantum spin or measured calibration')


def series_product(a, b, order):
    result = []
    for n in range(order+1):
        coefficient = m.zero(3)
        for i in range(n+1):
            if i < len(a) and n-i < len(b):
                coefficient = m.add(coefficient, m.mm(a[i], b[n-i]))
        result.append(coefficient)
    return result


def exponential_series(matrix, order):
    answer, power = [], m.eye(3)
    for n in range(order+1):
        answer.append(m.scale(Q(1, factorial(n)), power))
        power = m.mm(power, matrix)
    return answer


def rotation_and_coframe_checks():
    order, count = 6, 0
    for shear in [m.scale(Q(1, 3), [[1, 2, 0], [2, -2, 1], [0, 1, 1]]),
                  [[Q(1), Q(0), Q(0)], [Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(-1)]]]:
        omega = [[Q(0), Q(1, 5), Q(-1, 3)], [Q(-1, 5), Q(0), Q(2, 7)], [Q(1, 3), Q(-2, 7), Q(0)]]
        require(m.transpose(shear) == shear and m.trace(shear) == 0, 'arbitrary symmetric traceless initial shear')
        require(m.trace(m.mm(shear, m.comm(shear, omega))) == 0, 'Lax shear invariant')
        rotation = exponential_series(m.scale(-1, omega), order+1)
        inverse_rotation = exponential_series(omega, order+1)
        rotated = series_product(series_product(rotation, [shear], order+1), inverse_rotation, order+1)
        coframe = series_product(rotation, exponential_series(m.add(shear, omega), order+1), order+1)
        rhs = series_product(rotated, coframe, order)
        require([m.scale(n+1, coframe[n+1]) for n in range(order+1)] == rhs,
                'ordered exponential coframe satisfies F_u=S(u) F through exact coefficients')
        require([m.scale(n+1, rotated[n+1]) for n in range(order+1)] ==
                [m.comm(rotated[n], omega) for n in range(order+1)], 'rotating shear solves the Lax equation')
        wrong = exponential_series(shear, order+1)
        require([m.scale(n+1, wrong[n+1]) for n in range(order+1)] != series_product(rotated, wrong, order),
                'dropping noncommuting spin rotation from coframe is detected')
        count += 1
    # At any instantaneous orientation the commutator acts skew on the five shear coordinates.
    basis = [[[Q(int(i == j == 0))-Q(int(i == j == 2)) for j in SPACE] for i in SPACE],
             [[Q(int(i == j == 1))-Q(int(i == j == 2)) for j in SPACE] for i in SPACE]]
    for i, j in combinations(SPACE, 2):
        unit = m.zero(3)
        unit[i][j] = unit[j][i] = Q(1)
        basis.append(unit)
    pairs = 0
    for a, b in product(basis, repeat=2):
        require(m.trace(m.mm(a, m.comm(b, omega)))+m.trace(m.mm(m.comm(a, omega), b)) == 0,
                'five-dimensional homogeneous shear generator is skew in Frobenius metric')
        pairs += 1
    return dict(noncommuting_coframe_series=2, exact_series_order=order,
                five_shear_mode_pairings=pairs, missing_rotation_control=True)


def diagonal_and_blind_readout_checks():
    xi = m.mv(m.add(m.eye(8), sm.H), [Q(i == 0) for i in range(8)])
    omega = m.scale(Q(1, 6), spin_matrix(xi))  # kappa=1, Z=4/3
    aligned = [[Q(1), Q(0), Q(0)], [Q(0), Q(-2), Q(0)], [Q(0), Q(0), Q(1)]]
    tilted = [[Q(1), Q(0), Q(0)], [Q(0), Q(1), Q(0)], [Q(0), Q(0), Q(-2)]]
    require(m.comm(aligned, omega) == m.zero(3) and m.comm(tilted, omega) != m.zero(3),
            'aligned locally rotationally symmetric sector versus rotating sector')
    require(shear_square(aligned) == shear_square(tilted) == 3, 'identical volume shear invariant')
    norms = [m.trace(m.mm(m.comm(s, omega), m.comm(s, omega))) for s in (aligned, tilted)]
    require(norms == [0, 2], 'equal scalar data do not fix anisotropic stress norm')
    benchmark = []
    for z in (Q(2), Q(3), Q(5), Q(9)):
        v, dv, _ = volume_law(Q(1), Q(1, 2), Q(1), Q(4, 3), Q(2), z)
        require(v == (z-1)*(3*z+1)/(2*z), 'explicit anisotropic benchmark volume')
        benchmark.append(dict(z=str(z), volume=str(v), mean_rate=str(dv/(3*v)),
                              shear_square=str(1/(v*v)), anisotropic_stress_square=str(8/(9*v**4))))
    return dict(diagonal_spin_alignment_required=True,
                equal_volume_invariant='D^2=3', distinct_scaled_stress_norms=[str(x) for x in norms],
                benchmark=benchmark, benchmark_future_clock='(1/2) log(7/3)',
                benchmark_future_axis_rotation='(1/6) log(7/3)',
                benchmark_future_contact_phase='(1/4) log(7/3)',
                full_spatially_inhomogeneous_stability='NOT_ESTABLISHED')


def run():
    return dict(algebra=universal_matrix_checks(), original=original_equation_checks(),
                volume=volume_and_spin_clock_checks(), coframe=rotation_and_coframe_checks(),
                benchmark=diagonal_and_blind_readout_checks())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
