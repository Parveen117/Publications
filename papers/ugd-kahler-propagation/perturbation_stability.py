"""PS exact perturbation controls, consuming the unchanged SM and CS modules."""
from fractions import Fraction as Q
from itertools import product

import model as m
import metric_dynamics as mg
import loop_curvature as lr
import spin_matter as sm
import coupled_cosmology as cs
from thermo_gauge import Jet

IX, PAIRS, ETA = sm.IX, sm.PAIRS, sm.ETA
require = mg.require
SCALAR, PSEUDO, AXIAL = sm.H, m.mm(sm.H, sm.VOLUME), m.mm(sm.PHASE, sm.VOLUME)
PLUS = m.scale(Q(1, 2), m.add(m.eye(8), sm.H))
MINUS = m.add(m.eye(8), m.scale(-1, PLUS))


def outer(x, y=None):
    y = x if y is None else y
    return [[a*b for b in y] for a in x]


def nonlinear_flow(x, mass, coefficient):
    scalar, pseudo = sm.bil(x, SCALAR), sm.bil(x, PSEUDO)
    inside = [(mass+coefficient*scalar)*u+coefficient*pseudo*v
              for u, v in zip(x, m.mv(sm.VOLUME, x))]
    return m.mv(m.scale(-1, sm.GAMMA[0]), inside)


def bilinear_flow_checks():
    matrices = [SCALAR, PSEUDO, AXIAL]
    for i, x in enumerate(matrices):
        require(m.transpose(x) == x and m.mm(x, x) == m.eye(8), 'three real symmetric involutions')
        for y in matrices[i+1:]:
            require(m.anti(x, y) == m.zero(8), 'three-bilinear Clifford sphere bound')
    tables = [[m.zero(8), m.scale(-2, AXIAL), m.scale(2, PSEUDO)],
              [m.scale(2, AXIAL), m.zero(8), m.scale(-2, SCALAR)]]
    for generator, row in zip([sm.GAMMA[0], m.mm(sm.GAMMA[0], sm.VOLUME)], tables):
        require(m.transpose(generator) == m.scale(-1, generator), 'homogeneous nonlinear generator skew')
        for matrix, expected in zip(matrices, row):
            require(m.comm(matrix, generator) == expected, 'universal bilinear commutator')
    count = 0
    for seed, volume in product((1, 2, 3), (Q(1), Q(3, 2))):
        x = [Q((i*i+seed*i) % 7-3, 5) for i in range(8)]
        mass, kappa, norm, rate = Q(2, 3), Q(5, 7), Q(3, 2), Q(4, 5)
        coefficient = 3*kappa*norm/(16*volume)
        scalar, pseudo, axial = [sm.bil(x, a) for a in matrices]
        dx = nonlinear_flow(x, mass, coefficient)
        actual = [2*sm.bil(x, a, dx) for a in matrices]
        expected = [-2*coefficient*pseudo*axial,
                    2*(mass+coefficient*scalar)*axial, -2*mass*pseudo]
        require(actual == expected and m.dot(x, dx) == 0, 'direct native flow versus bilinear equations and conserved norm')
        require(sum(a*b for a, b in zip((scalar, pseudo, axial), actual)) == 0, 'conserved bilinear radius')
        require(scalar*scalar+pseudo*pseudo+axial*axial <= m.dot(x, x)**2, 'Clifford sphere inequality')
        sj, pj = [Jet(value, (derivative, Q(0), Q(0), Q(0))) for value, derivative in zip((scalar, pseudo), actual)]
        vj = Jet(volume, (3*rate*volume, Q(0), Q(0), Q(0)))
        rho = norm*mass*sj/(2*vj)+3*kappa*norm*norm*(sj*sj+pj*pj)/(64*vj*vj)
        pressure = 3*kappa*norm*norm*(sj*sj+pj*pj)/(64*vj*vj)
        require(rho.grad[0]+3*rate*(rho.value+pressure.value) == 0, 'general homogeneous stress continuity')
        count += 1
    return dict(universal_commutator_identities=6, involutions=3, general_bilinear_jets=count,
                conserved_number_and_bilinear_radius=True)


def nonrest_fixture(seed):
    kappa, norm, mass = Q(2, 3), Q(3, 2), Q(1, 2)
    reference = m.mv(m.add(m.eye(8), sm.H), [Q(i == 0) for i in range(8)])
    psi = [x+Q((i*seed+i*i) % 7-3, 30) for i, x in enumerate(reference)]
    scalar, pseudo = sm.bil(psi, SCALAR), sm.bil(psi, PSEUDO)
    require(scalar > 0 and pseudo != 0 and sm.bil(psi, AXIAL) != 0, 'fixture leaves the old rest sector')
    rho = norm*mass*scalar/2+3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
    pressure = 3*kappa*norm*norm*(scalar*scalar+pseudo*pseudo)/64
    rate = Q(1)+Q(seed, 5)
    rate_dot, cosmological = -kappa*(rho+pressure)/2, 3*rate*rate-kappa*rho
    require(cosmological > 0, 'positive cosmological fixture')
    dp = [x-3*rate*y/2 for x, y in zip(nonlinear_flow(psi, mass, 3*kappa*norm/16), psi)]
    psi_j = [Jet(x, (dx, Q(0), Q(0), Q(0))) for x, dx in zip(psi, dp)]
    lengths = [Q(1), Q(1), Q(2), Q(3, 2)]
    e = [[Jet(Q(a == mu)*lengths[a], (Q(a == mu and a != 0)*rate*lengths[a], Q(0), Q(0), Q(0)))
          for mu in IX] for a in IX]
    rate_j, ax = Jet(rate, (rate_dot, Q(0), Q(0), Q(0))), sm.axial(psi_j)
    omega_lc = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    omega = [[Jet.lift(0) for _ in IX] for _ in PAIRS]
    for p, (a, b) in enumerate(PAIRS):
        for mu in IX:
            if a == 0:
                omega_lc[p][mu] = rate_j*e[b][mu]
            omega[p][mu] = omega_lc[p][mu]+kappa*norm*sum(ETA[c]*ax[a, b, c]*e[c][mu] for c in IX)/8
    return dict(e=e, ef=mg.coframe_forms(e), omega=omega, omega_lc=omega_lc,
                psi=psi, dpsi=[dp]+[[Q(0)]*8 for _ in range(3)], mass=mass, kappa=kappa,
                norm=norm, rate=rate, rate_dot=rate_dot, cosmological=cosmological,
                energy=rho, pressure=pressure, n=norm*scalar/2)


def homogeneous_completion_checks():
    for seed in (1, 2):
        f = nonrest_fixture(seed)
        for name, values in cs.first_order_residuals(f).items():
            require(all(x == 0 for x in values), 'nonrest full '+name+' residual')
        cs.effective_stress_checks(f)
    return dict(nonrest_first_order_jets=2, coframe_components=32, cartan_components=48,
                matter_components=16, effective_stress_components=32,
                nonlinear_completion='Homogeneous FLRW metric with all eight matter components; not a solved general inhomogeneous PDE')


def rest_vectors():
    result = []
    for raw in [[Q(i == 0) for i in range(8)], [Q(i % 3-1, 3) for i in range(8)]]:
        x = m.mv(m.add(m.eye(8), sm.H), raw)
        result.append(x)
        phase = m.add(m.scale(Q(3, 5), m.eye(8)), m.scale(Q(4, 5), sm.PHASE))
        result.append(m.mv(phase, x))
    return result


def perturbation_matrix(x, mass, shift):
    number = m.dot(x, x)
    jx = m.mv(sm.VOLUME, x)
    projector = m.scale(1/number, m.add(outer(x), outer(jx)))
    return m.add(m.scale(2*(mass+shift), m.mm(sm.PHASE, MINUS)),
                 m.scale(-2*shift, m.mm(sm.PHASE, projector)))


def linearization_and_energy_checks():
    derivatives = frames = fourier = 0
    for x, mass, shift in product(rest_vectors(), (Q(0), Q(2, 3)), (Q(1, 5), Q(3, 4))):
        number, effective = m.dot(x, x), mass+shift
        coefficient = shift/number
        expected = perturbation_matrix(x, mass, shift)
        columns = []
        for j in range(8):
            xj = [Jet(value, (Q(i == j), Q(0), Q(0), Q(0))) for i, value in enumerate(x)]
            differentiated = nonlinear_flow(xj, mass, coefficient)
            columns.append([value.grad[0] for value in differentiated])
            derivatives += 1
        actual = m.add(m.transpose(columns), m.scale(effective, sm.PHASE))
        require(actual == expected, 'independent cubic Jacobian in the rotating native phase frame')
        ix, jx = m.mv(sm.PHASE, x), m.mv(sm.VOLUME, x)
        ijx = m.mv(sm.PHASE, jx)
        frame = [x, ix, jx, ijx]
        require([[m.dot(a, b) for b in frame] for a in frame] == m.scale(number, m.eye(4)),
                'orthogonal four-vector amplitude/phase/transverse frame')
        targets = [m.mv(m.scale(-2*shift, sm.PHASE), x), [Q(0)]*8,
                   [2*mass*y for y in ijx], [-2*effective*y for y in jx]]
        require([m.mv(actual, z) for z in frame] == targets, 'all nontrivial homogeneous perturbation blocks')
        pp = m.scale(1/number, m.add(outer(x), outer(ix)))
        pm = m.scale(1/number, m.add(outer(jx), outer(ijx)))
        plus_free, minus_free = m.add(PLUS, m.scale(-1, pp)), m.add(MINUS, m.scale(-1, pm))
        require(m.rank(plus_free) == m.rank(minus_free) == 2, 'two remaining components in each scalar sector')
        require(m.mm(actual, plus_free) == m.zero(8)
                and m.mm(actual, minus_free) == m.scale(2*effective, m.mm(sm.PHASE, minus_free)),
                'neutral orientation and skew-rotation complement')
        symmetric = m.scale(Q(1, 2), m.add(actual, m.transpose(actual)))
        support = m.add(pp, pm)
        require(m.mm(symmetric, symmetric) == m.scale(shift*shift, support), 'exact logarithmic norm support')
        grow, decay = [a-b for a, b in zip(x, ix)], [a+b for a, b in zip(x, ix)]
        require(m.mv(symmetric, grow) == [shift*a for a in grow]
                and m.mv(symmetric, decay) == [-shift*a for a in decay], 'instantaneous gain bound is sharp')
        require(symmetric != m.zero(8), 'linearized norm is not conserved merely because the nonlinear norm is')
        # Real cosine/sine Fourier block; test arbitrary directions without complex arithmetic.
        for k in [(Q(1), Q(0), Q(0)), (Q(2), Q(-1), Q(3))]:
            spatial = sm.lin(k, sm.A)
            block = [actual[i]+spatial[i] for i in range(8)]
            block += [[-z for z in spatial[i]]+actual[i] for i in range(8)]
            sym16 = m.scale(Q(1, 2), m.add(block, m.transpose(block)))
            require(sym16 == m.kron(m.eye(2), symmetric), 'spatial Fourier frequency drops out of the energy estimate')
            fourier += 1
        # A volume perturbation changes c=alpha/v and drives only the phase direction.
        delta_volume_ratio = Q(2, 7)
        cj = Jet(coefficient, (-coefficient*delta_volume_ratio, Q(0), Q(0), Q(0)))
        forcing = nonlinear_flow(x, mass, cj)
        require([v.grad[0] for v in forcing] == [shift*delta_volume_ratio*v for v in ix],
                'coupled homogeneous metric forcing retained')
        frames += 1
    return dict(independent_cubic_jacobian_columns=derivatives, full_eight_component_frames=frames,
                real_fourier_blocks=fourier, exact_symmetric_spectrum='(+s,+s,-s,-s,0,0,0,0)',
                positive_instantaneous_gain_witness=True, volume_phase_forcing_checked=True)


def coupled_metric_and_transverse_checks():
    metric_count = transverse_count = 0
    for scale, mass, n0, z in product((Q(1, 6), Q(1, 4)), (Q(0), Q(1, 2)),
                                     (Q(2, 3), Q(4, 3)), (Q(2), Q(3))):
        kappa, dn0, dt0 = Q(1), Q(2, 7), Q(-3, 11)
        volume, velocity, acceleration = cs.positive_volume(scale, mass, kappa, n0, z)
        third = (6*scale)**2*velocity
        rate = velocity/(3*volume)
        rate_dot = acceleration/(3*volume)-velocity*velocity/(3*volume*volume)
        rate_ddot = third/(3*volume)-velocity*acceleration/(volume*volume)+2*velocity**3/(3*volume**3)
        n = n0/volume
        dv, ddv = volume*dn0/n0-velocity*dt0, velocity*dn0/n0-acceleration*dt0
        dh = (ddv*volume-velocity*dv)/(3*volume*volume)
        dn = dn0/volume-n0*dv/(volume*volume)
        require(dh == -rate_dot*dt0 and dn == 3*rate*n*dt0, 'two constrained cosmological tangent modes')
        drho, dp = (mass+3*kappa*n/8)*dn, 3*kappa*n*dn/8
        require(6*rate*dh-kappa*drho == 0, 'linearized lapse constraint')
        require(-2*rate_ddot*dt0+6*rate*dh+kappa*dp == 0, 'linearized spatial Einstein equation')
        shift, shift_dot = 3*kappa*n/8, -9*kappa*rate*n/8
        amplitude = dn0/(2*n0)
        require(-2*shift*amplitude+shift*dv/volume == shift_dot*dt0,
                'metric backreaction cancels the spurious amplitude-induced phase shear')
        metric_count += 1
        for pseudo, axial in [(Q(2, 3), Q(-1, 5)), (Q(-1), Q(2))]:
            pseudo_dot, axial_dot = 2*(mass+shift)*axial, -2*mass*pseudo
            energy_dot = 2*mass*pseudo*pseudo_dot+shift_dot*axial*axial+2*(mass+shift)*axial*axial_dot
            require(energy_dot == shift_dot*axial*axial <= 0, 'transverse Lyapunov energy identity')
            transverse_count += 1
    return dict(constrained_metric_tangent_checks=metric_count,
                transverse_energy_checks=transverse_count, all_eight_homogeneous_components_accounted_for=True,
                scope='Linear homogeneous coupled FLRW sector; not anisotropic or inhomogeneous metric stability')


def gain_and_gradient_checks():
    gains = []
    for z in (Q(2), Q(3), Q(5), Q(9)):
        squared_gain = 2*(1-1/z)
        volume = z-1
        require(1 <= squared_gain < 2 and squared_gain/volume == 2/z,
                'benchmark integrated bound and physical field dilution')
        gains.append(dict(z=str(z), rescaled_gain_squared=str(squared_gain),
                          physical_gain_squared=str(2/z)))
    primitive_count = 0
    for scale, mass, n0, z0, z in product((Q(1, 6), Q(1, 4)), (Q(0), Q(2)),
                                         (Q(2, 3),), (Q(3, 2), Q(2)), (Q(3), Q(5))):
        kappa, omega, cosmological = Q(1), 6*scale, 12*scale*scale
        d, b = kappa*mass*n0/(2*cosmological), 3*kappa*n0/(4*omega)
        gain_sq = (z-1)*((d+b)*z0-(d-b))/((z0-1)*((d+b)*z-(d-b)))
        future_sq = ((d+b)*z0-(d-b))/((d+b)*(z0-1))
        v = cs.positive_volume(scale, mass, kappa, n0, z)[0]
        log_gain_derivative = omega*z*(1/(z-1)-(d+b)/((d+b)*z-(d-b)))/2
        require(log_gain_derivative == 3*kappa*n0/(8*v), 'exact gain primitive derivative')
        require(1 <= gain_sq <= future_sq, 'finite future gain bound for all tested positive-Lambda parameters')
        primitive_count += 1
    leakage_count = 0
    for k in [(Q(1), Q(0), Q(0)), (Q(0), Q(2), Q(0)), (Q(2), Q(-1), Q(3))]:
        spatial = sm.lin(k, sm.A)
        leak = m.mm(m.mm(MINUS, spatial), PLUS)
        require(m.mm(m.mm(PLUS, spatial), PLUS) == m.zero(8), 'gradients leave the positive rest sector')
        require(m.rank(leak) == 4 and m.mm(m.transpose(leak), leak) == m.scale(m.dot(k, k), PLUS),
                'every nonzero rest-sector plane wave has a transverse gradient component')
        leakage_count += 1
    # Massless zero-Lambda transverse system in logarithmic time L=log(t/t0).
    nilpotent = [[Q(0), Q(1)], [Q(0), Q(0)]]
    require(m.mm(nilpotent, nilpotent) == m.zero(2), 'zero-Lambda massless shear')
    shear_norms = []
    for log_time in (Q(0), Q(1), Q(3)):
        propagator = m.add(m.eye(2), m.scale(log_time, nilpotent))
        moved = m.mv(propagator, [Q(0), Q(1)])
        require(moved == [log_time, 1], 'unbounded relative transverse shear control')
        shear_norms.append(str(m.dot(moved, moved)))
    return dict(benchmark_start='t*=log(2), a(t*)=1', benchmark=gains,
                benchmark_future_rescaled_gain_squared='2', benchmark_physical_gain_squared='2 exp(-t)',
                integrated_gain_checks=primitive_count, gradient_leakage_ranks=[4]*leakage_count,
                zero_Lambda_massless_shear_norms=shear_norms,
                full_inhomogeneous_metric_stability='NOT_ESTABLISHED')


def run():
    return dict(bilinears=bilinear_flow_checks(), homogeneous=homogeneous_completion_checks(),
                linearization=linearization_and_energy_checks(),
                coupled=coupled_metric_and_transverse_checks(), bounds=gain_and_gradient_checks())


if __name__ == '__main__':
    import json
    print(json.dumps(run(), indent=2, sort_keys=True))
