"""CS exact controls for a sourced solution; native/SM engines remain upstream."""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, product

import model as m
import metric_dynamics as mg
import loop_curvature as lr
import spin_matter as sm
from thermo_gauge import Jet

IX, PAIRS, TOP, ETA = sm.IX, sm.PAIRS, sm.TOP, sm.ETA
require = mg.require


def cubic_identity_checks():
    scalar, pseudo = sm.quadratic_polynomial(sm.H), sm.quadratic_polynomial(m.mm(sm.H, sm.VOLUME))
    axial = {abc: sm.quadratic_polynomial(matrix) for abc, matrix in sm.HG3.items()}
    monomials = list(combinations_with_replacement(range(8), 3))
    for j in range(8):
        coefficients = {}
        for (a, b, c), polynomial in axial.items():
            for key, value in polynomial.items():
                for k, matrix_value in enumerate(sm.G3[a, b, c][j]):
                    index = tuple(sorted(key+(k,)))
                    coefficients[index] = coefficients.get(index, Q(0))+ETA[a]*ETA[b]*ETA[c]*value*matrix_value
        for polynomial, matrix in [(scalar, m.eye(8)), (pseudo, sm.VOLUME)]:
            for key, value in polynomial.items():
                for k, matrix_value in enumerate(matrix[j]):
                    index = tuple(sorted(key+(k,)))
                    coefficients[index] = coefficients.get(index, Q(0))+6*value*matrix_value
        require(all(coefficients.get(key, 0) == 0 for key in monomials), 'complete cubic Fierz identity')
    projector = m.scale(Q(1, 2), m.add(m.eye(8), sm.H))
    require(m.mm(projector, projector) == projector and m.rank(projector) == 4, 'positive scalar eigenspace')
    require(m.mm(sm.GAMMA[0], projector) == m.mm(sm.PHASE, projector), 'rest equation is native phase rotation')
    require(m.comm(projector, sm.PHASE) == m.zero(8), 'phase preserves the rest eigenspace')
    for matrix in sm.A+[m.mm(sm.H, sm.VOLUME), sm.HG3[1, 2, 3]]:
        require(m.mm(m.mm(projector, matrix), projector) == m.zero(8), 'rest spatial current, pseudoscalar and axial time vanish')
    return dict(cubic_components=8,monomials_per_component=len(monomials),
                rest_projector_rank=4, rest_sector_preserved=True)


def positive_volume(scale, mass, kappa, number, z):
    """z=exp(6 scale t); exact values and time derivatives, without transcendental rounding."""
    omega = 6*scale
    cosmological = 12*scale*scale
    d = kappa*mass*number/(2*cosmological)
    b = 3*kappa*number/(4*omega)
    ch, sh = (z+1/z)/2, (z-1/z)/2
    volume = d*(ch-1)+b*sh
    velocity = omega*(d*sh+b*ch)
    acceleration = omega*omega*(d*ch+b*sh)
    return volume, velocity, acceleration


def solution_fixture(scale, mass, kappa, z, a, seed):
    unit_volume = positive_volume(scale, mass, kappa, Q(1), z)[0]
    number = a**3/unit_volume
    volume, velocity, acceleration = positive_volume(scale, mass, kappa, number, z)
    rate = velocity/(3*volume)
    rate_dot = acceleration/(3*volume)-velocity*velocity/(3*volume*volume)
    n = number/volume
    energy, pressure = mass*n+3*kappa*n*n/16, 3*kappa*n*n/16
    cosmological = 12*scale*scale
    require(3*rate*rate == cosmological+kappa*energy, 'exact volume Hamiltonian constraint')
    require(rate_dot == -kappa*(energy+pressure)/2, 'exact volume acceleration equation')
    raw = [Q(i == 0) for i in range(8)] if seed == 1 else [Q(((i+1)*seed) % 5-2) for i in range(8)]
    psi = m.mv(m.add(m.eye(8), sm.H), raw)
    require(any(psi), 'nonzero rest-state fixture')
    scalar = sm.bil(psi, sm.H)
    norm = 2*n/scalar
    effective_mass = mass+3*kappa*n/8
    dpsi = [-3*rate*x/2-effective_mass*y for x, y in zip(psi, m.mv(sm.PHASE, psi))]
    psi_j = [Jet(x, (dx, Q(0), Q(0), Q(0))) for x, dx in zip(psi, dpsi)]
    lengths = [Q(1), Q(2), Q(3, 2)]
    e = [[Jet(Q(a0 == mu)*(Q(1) if a0 == 0 else a*lengths[a0-1]),
              (Q(a0 == mu and a0 != 0)*rate*a*(lengths[a0-1] if a0 else 1), Q(0), Q(0), Q(0)))
          for mu in IX] for a0 in IX]
    rate_j = Jet(rate, (rate_dot, Q(0), Q(0), Q(0)))
    ef, ax = mg.coframe_forms(e), sm.axial(psi_j)
    omega_lc = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    omega = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    for p, (b, c) in enumerate(PAIRS):
        for mu in IX:
            if b == 0:
                omega_lc[p][mu] = rate_j*e[c][mu]
            omega[p][mu] = omega_lc[p][mu]+kappa*norm*sum(ETA[d]*ax[b, c, d]*e[d][mu] for d in IX)/8
    return dict(scale=scale, mass=mass, kappa=kappa, cosmological=cosmological, norm=norm,
                number=number, volume=volume, rate=rate, rate_dot=rate_dot, energy=energy,
                pressure=pressure, n=n, effective_mass=effective_mass,
                e=e, ef=ef, psi=psi, psi_j=psi_j, dpsi=[dpsi]+[[Q(0)]*8 for _ in range(3)],
                omega=omega, omega_lc=omega_lc)


def curvature_and_torsion(ef, omega):
    om = [[{} for _ in IX] for _ in IX]
    for p, (a, b) in enumerate(PAIRS):
        om[a][b] = {(mu,): omega[p][mu] for mu in IX}
        om[b][a] = mg.fscale(-1, om[a][b])
    torsion = [mg.fadd(mg.exterior_d(ef[a]), *(mg.fscale(ETA[b], mg.wedge(om[a][b], ef[b])) for b in IX))
               for a in IX]
    curvature = [[mg.fadd(mg.exterior_d(om[a][b]),
                          *(mg.fscale(ETA[c], mg.wedge(om[a][c], om[c][b])) for c in IX))
                  for b in IX] for a in IX]
    return [[mg.values(x) for x in row] for row in curvature], [mg.values(x) for x in torsion]


def first_order_residuals(f, omit_contorsion=False):
    e, ef = lr.matrix_values(f['e']), [mg.values(x) for x in f['ef']]
    omega_j = f['omega_lc'] if omit_contorsion else f['omega']
    omega = lr.matrix_values(omega_j)
    rf, tf = curvature_and_torsion(f['ef'], omega_j)
    spin_curvature = mg.matrix_form([rf[a][b] for a, b in PAIRS],
                                   [m.scale(Q(1, 2), m.mm(sm.OLD_GAMMA[a], sm.OLD_GAMMA[b])) for a, b in PAIRS])
    psi, dp, norm, kappa, mass = f['psi'], f['dpsi'], f['norm'], f['kappa'], f['mass']
    w, inv, gauge = m.det(e), m.inverse(e), [Q(0)]*4
    gam = [sm.lin(inv[mu], sm.GAMMA) for mu in IX]
    conn = sm.connections(omega, gauge, Q(0))
    cov = [[x+y for x, y in zip(dp[mu], m.mv(conn[mu], psi))] for mu in IX]
    trace_t = sm.torsion_trace(e, tf)
    matter = [sum(m.mv(gam[mu], cov[mu])[j]+trace_t[mu]*m.mv(gam[mu], psi)[j]/2
                  for mu in IX)-mass*psi[j] for j in range(8)]
    coframe = []
    for a, mu in product(IX, repeat=2):
        ep, em = [row[:] for row in e], [row[:] for row in e]
        ep[a][mu], em[a][mu] = ep[a][mu]+Q(1, 17), em[a][mu]-Q(1, 17)
        gravitational = (mg.gravity_density(ep, spin_curvature, f['cosmological'], Q(0), kappa)
                         -mg.gravity_density(em, spin_curvature, f['cosmological'], Q(0), kappa))*Q(17, 2)
        matter_source = (sm.density(ep, omega, gauge, psi, dp, mass, norm, Q(0))
                         -sm.density(em, omega, gauge, psi, dp, mass, norm, Q(0)))*Q(17, 2)
        coframe.append(gravitational+matter_source)
    ax = sm.axial(psi)
    cartan = []
    for c, d in PAIRS:
        lhs = mg.fadd(*(mg.fscale(mg.sign((a, b, c, d)), mg.wedge(tf[a], ef[b]))
                       for a, b in product(IX, repeat=2)))
        spin = {tuple(i for i in IX if i != mu):
                norm*w*ETA[c]*ETA[d]*sum(inv[mu][a]*ax[a, c, d] for a in IX)*(-1)**mu/4 for mu in IX}
        residual = mg.fadd(lhs, mg.fscale(kappa, spin))
        cartan += [residual.get(abc, Q(0)) for abc in mg.TRIPLES]
    scalar_curvature = 2*kappa*mg.gravity_density(e, spin_curvature, Q(0), Q(0), kappa)/w
    scalar_expected = 4*f['cosmological']+kappa*mass*f['n']
    return dict(matter=matter, coframe=coframe, cartan=cartan, trace=trace_t,
                scalar_curvature=[scalar_curvature-scalar_expected])


def effective_stress_checks(f):
    e, psi, dp = lr.matrix_values(f['e']), f['psi'], f['dpsi']
    norm, mass, kappa = f['norm'], f['mass'], f['kappa']
    conn = sm.connections(lr.matrix_values(f['omega_lc']), [Q(0)]*4, Q(0))
    cov = [[x+y for x, y in zip(dp[mu], m.mv(conn[mu], psi))] for mu in IX]
    metric, inv = mg.metric(e), m.inverse(e)
    gam = [sm.lin(inv[mu], sm.GAMMA) for mu in IX]
    lower_gam = [sm.lin(metric[mu], gam) for mu in IX]
    scalar, pseudo = sm.bil(psi, sm.H), sm.bil(psi, m.mm(sm.H, sm.VOLUME))
    potential = norm*mass*scalar/2+3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
    kinetic = norm*sum(sm.bil(psi, m.mm(sm.H, gam[mu]), cov[mu]) for mu in IX)/2
    lag = kinetic-potential
    stress = [[-norm*(sm.bil(psi, m.mm(sm.H, lower_gam[mu]), cov[nu])
                           +sm.bil(psi, m.mm(sm.H, lower_gam[nu]), cov[mu]))/4
               +metric[mu][nu]*lag for nu in IX] for mu in IX]
    expected = [[(f['energy'] if mu == 0 else metric[mu][mu]*f['pressure'])
                 if mu == nu else Q(0) for nu in IX] for mu in IX]
    require(stress == expected, 'all Hilbert stress components from full spinor bilinears')
    require(potential == f['energy'] and lag == f['pressure'], 'on-shell energy and pressure signs')
    # Vary lapse BEFORE choosing cosmic time, and scale before using the field equation.
    a, adot = e[1][1], f['rate']*e[1][1]
    for vary_lapse in (True, False):
        lapse = Jet(Q(1), (Q(vary_lapse), Q(0), Q(0), Q(0)))
        scale = Jet(a, (Q(not vary_lapse), Q(0), Q(0), Q(0)))
        ej = [[(lapse if b == 0 else scale) if b == mu else Jet.lift(0) for mu in IX] for b in IX]
        om = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
        for p, (b, c) in enumerate(PAIRS):
            if b == 0:
                om[p][c] = adot/lapse
        matter = sm.density(ej, om, [Q(0)]*4, psi, dp, mass, norm, Q(0))
        contact = -lapse*scale*scale*scale*3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
        actual = (matter+contact).grad[0]
        target = -a**3*f['energy'] if vary_lapse else 3*a*a*f['pressure']
        require(actual == target, 'independent lapse/scale variation of the covariant eliminated density')


def volume_and_constraint_checks():
    count = zero_count = 0
    for scale, mass, kappa, number, z in product((Q(1, 6), Q(1, 4)), (Q(0), Q(1, 2), Q(2)),
                                               (Q(1), Q(3, 5)), (Q(2, 3), Q(5, 4)), (Q(3, 2), Q(2), Q(3))):
        v, dv, ddv = positive_volume(scale, mass, kappa, number, z)
        cosmological, omega = 12*scale*scale, 6*scale
        require(v > 0 and dv > 0, 'positive expanding branch')
        require(dv*dv == 3*cosmological*v*v+3*kappa*mass*number*v+9*kappa*kappa*number*number/16,
                'closed volume first integral')
        require(ddv == 3*cosmological*v+3*kappa*mass*number/2, 'closed volume acceleration')
        d, b = kappa*mass*number/(2*cosmological), 3*kappa*number/(4*omega)
        phase_rate = mass+omega*z*(1/(z-1)-(d+b)/((d+b)*z-(d-b)))/2
        require(phase_rate == mass+3*kappa*number/(8*v), 'derivative of explicit logarithmic spinor phase')
        n, rate = number/v, dv/(3*v)
        energy, pressure = mass*n+3*kappa*n*n/16, 3*kappa*n*n/16
        energy_dot = (mass+3*kappa*n/8)*(-3*rate*n)
        require(energy_dot+3*rate*(energy+pressure) == 0, 'matter continuity from exact solution')
        count += 1
    for mass, kappa, number, t in product((Q(0), Q(1, 2), Q(2)), (Q(1), Q(3, 5)),
                                         (Q(2, 3), Q(5, 4)), (Q(1, 3), Q(1), Q(2))):
        v = 3*kappa*number*t*(1+mass*t)/4
        dv = 3*kappa*number*(1+2*mass*t)/4
        require(dv*dv == 3*kappa*mass*number*v+9*kappa*kappa*number*number/16,
                'zero-Lambda contraction first integral')
        require(mass+1/(2*t*(1+mass*t)) == mass+3*kappa*number/(8*v), 'zero-Lambda exact phase')
        zero_count += 1
    # General off-constraint data: evolution by the spatial Einstein equation.
    constraint_count = 0
    for rate, n, mass, kappa, cosmological in product((Q(-1), Q(2, 3)), (Q(1), Q(3, 2)),
                                                    (Q(0), Q(2)), (Q(1),), (Q(0), Q(1, 3))):
        rho, p = mass*n+3*kappa*n*n/16, 3*kappa*n*n/16
        defect = 3*rate*rate-cosmological-kappa*rho
        rate_dot = (cosmological-3*rate*rate-kappa*p)/2
        rho_dot = (mass+3*kappa*n/8)*(-3*rate*n)
        defect_dot = 6*rate*rate_dot-kappa*rho_dot
        require(defect_dot == -3*rate*defect, 'Hamiltonian constraint propagation without assuming the constraint')
        constraint_count += 1
    return dict(positive_Lambda_volume_phase_checks=count, zero_Lambda_contraction_checks=zero_count,
                off_constraint_propagation_checks=constraint_count,
                lapse_constraint_retained=True)


def full_solution_checks():
    parameters = [(Q(1, 6), Q(1, 2), Q(1), Q(2), Q(1), 1),
                  (Q(1, 4), Q(2), Q(3, 5), Q(3), Q(3, 2), 2),
                  (Q(2, 5), Q(0), Q(7, 3), Q(5, 2), Q(2), 3)]
    for args in parameters:
        f = solution_fixture(*args)
        residuals = first_order_residuals(f)
        for name, components in residuals.items():
            require(all(x == 0 for x in components), 'full '+name+' equations, not only FLRW projection: '+str(components))
        effective_stress_checks(f)
        ax = sm.axial(f['psi'])
        require(any(ax.values()) and ax[1, 2, 3] == 0, 'nonzero polarized torsion is not a rotationally invariant spin average')
    f = solution_fixture(*parameters[0])
    bad = first_order_residuals(f, omit_contorsion=True)
    omitted = {name: sum(x != 0 for x in bad[name]) for name in ('matter', 'coframe', 'cartan')}
    require(all(omitted.values()), 'dropping the contorsion violates every equation family')
    e = lr.matrix_values(f['e'])
    source = -f['norm']*m.det(e)*sm.bil(f['psi'], sm.HG_PHASE[0])/2
    require(source != 0, 'unit-charge F=0 Gauss-law obstruction')
    return dict(solution_jets=len(parameters), coframe_equations=16*len(parameters),
                matter_equations=8*len(parameters), cartan_equations=24*len(parameters),
                effective_stress_components=16*len(parameters), lapse_scale_variations=2*len(parameters),
                omitted_contorsion_nonzero_residuals=omitted, charged_zero_field_rejected=True,
                full_scalar_curvature_checks=3,
                cosmological_sector='positive Lambda, finite native return-scale and readout')


def benchmark_checks():
    rows = []
    for z in (Q(3, 2), Q(2), Q(3), Q(5)):
        v, dv, ddv = positive_volume(Q(1, 6), Q(1, 2), Q(1), Q(4, 3), z)
        require(v == z-1 and dv == ddv == z, 'reported exponential-volume benchmark')
        n, rate = Q(4, 3)/v, dv/(3*v)
        rho, p = n/2+3*n*n/16, 3*n*n/16
        lc_scalar = Q(4, 3)+rho-3*p
        full_scalar = lc_scalar+3*n*n/8
        require(full_scalar == Q(4, 3)+n/2, 'contact cancels from full curvature trace')
        rows.append(dict(z=str(z), volume=str(v), rate=str(rate), number_density=str(n),
                         energy=str(rho), pressure=str(p), metric_scalar=str(lc_scalar),
                         connection_scalar=str(full_scalar), torsion_contraction=str(-3*n*n/2)))
    # For m=0 the full curvature trace is constant while this torsion invariant grows.
    massless = []
    for z in (Q(2), Q(3, 2), Q(5, 4)):
        v = positive_volume(Q(1, 6), Q(0), Q(1), Q(1), z)[0]
        n = 1/v
        lc_scalar, torsion_sq = Q(4, 3)-3*n*n/8, -3*n*n/2
        require(lc_scalar+3*n*n/8 == Q(4, 3), 'constant full scalar does not imply bounded torsion')
        massless.append(torsion_sq)
    require(massless[0] > massless[1] > massless[2], 'torsion contraction diverges toward zero volume')
    return dict(parameter_values=dict(kappa='1', Z='4/3', m='1/2', Lambda='1/3',
                                      loop_scale='1/6', readout_beta='18', conserved_number='4/3'),
                rows=rows, massless_constant_scalar_control=True,
                interpretation='Chosen dimensionless benchmark, not measured cosmological or particle constants')


def run():
    return dict(cubic=cubic_identity_checks(), volume=volume_and_constraint_checks(),
                coupled=full_solution_checks(), benchmark=benchmark_checks())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
