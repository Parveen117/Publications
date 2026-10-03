"""Read-only by default. Exact controls accompany, but do not replace, proofs."""
import argparse
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

import model as m
import metric_dynamics
import loop_curvature
import spin_matter
import coupled_cosmology
import perturbation_stability
import anisotropic_cosmology
import speed_calibration
import thermo_gauge

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def flatten(a):
    return [x for row in a for x in row]


def linear_combination(coeffs, basis):
    answer = m.zero(len(basis[0]))
    for c, b in zip(coeffs, basis):
        answer = m.add(answer, m.scale(c, b))
    return answer


def algebra_checks():
    native, gamma, a, j = m.fixture()
    ident = m.eye(4)
    i, r, k, l = native
    require(m.mm(r, r) == m.scale(-1, i), 'native return square')
    require(m.mm(k, k) == i and m.mm(l, l) == i, 'native cut squares')
    require(m.anti(k, r) == m.zero(2), 'native anticommutation')
    require(m.comm(k, r) == m.scale(2, l), 'complex compatibility K coefficient')
    require(m.comm(l, r) == m.scale(-2, k), 'complex compatibility L coefficient')
    k_can = [[0, 1], [1, 0]]
    require(m.mm(k_can, r) == k and m.scale(-1, k_can) == l, 'canonical basis change')
    eta = [-1, 1, 1, 1]
    for p, q in product(range(4), repeat=2):
        require(m.anti(gamma[p], gamma[q]) == m.scale(2*eta[p] if p == q else 0, ident),
                'all Clifford anticommutators')
    for p, q in product(range(3), repeat=2):
        require(m.anti(a[p], a[q]) == m.scale(2 if p == q else 0, ident), 'spatial closure')
    require(all(m.transpose(x) == x and m.trace(x) == 0 for x in a), 'symmetric traceless A')
    require(j == m.kron(k, r) and m.mm(j, j) == m.scale(-1, ident), 'native volume square')
    require(m.transpose(j) == m.scale(-1, j), 'native volume skew')
    require(all(m.comm(j, x) == m.zero(4) for x in a), 'native volume commutant')
    tensor_basis = [m.kron(x, y) for x, y in product(native, repeat=2)]
    require(m.rank([flatten(x) for x in tensor_basis]) == 16, 'full coefficient-sector basis')
    gram = [[m.trace(m.mm(m.transpose(x), y))/4 for y in tensor_basis] for x in tensor_basis]
    require(gram == m.eye(16), 'derived trace Gram')

    # Solve the independent 16-unknown linear multiplier problem.
    units = []
    for p, q in product(range(4), repeat=2):
        e = m.zero(4)
        e[p][q] = Q(1)
        units.append(e)
    comm_columns = [[v for x in a for v in flatten(m.comm(s, x))] for s in units]
    skew_columns = [flatten(m.add(s, m.transpose(s))) for s in units]
    constraints = m.transpose(comm_columns)
    with_skew = constraints+m.transpose(skew_columns)
    with_mass = with_skew+m.transpose([flatten(m.comm(s, gamma[0])) for s in units])
    require(m.rank(constraints) == 14, 'commutant dimension two')
    require(m.rank(with_skew) == 15, 'skew commutant dimension one')
    require(m.rank(with_mass) == 16, 'mass leaves no constant skew multiplier')

    # A full native cut has two retained real directions, not one.
    cut = m.scale(Q(1, 2), m.add(ident, m.kron(l, i)))
    require(m.mm(cut, cut) == cut and m.rank(cut) == 2, 'native cut rank')
    split = m.add(ident, m.scale(-2, cut))
    require(m.rank(m.add(split, ident)) == 2 and m.rank(m.add(split, m.scale(-1, ident))) == 2,
            'whole cut flip is neutral')
    require(all(m.anti(gamma[0], x) == m.zero(4) for x in a), 'gapped generator cross terms')
    require(m.anti(j, gamma[0]) == m.zero(4), 'mass anti-compatibility')
    return dict(clifford_pairs=16, spatial_pairs=9, trace_gram_entries=256,
                commutant_dimension=2, skew_multiplier_dimension=1,
                massive_skew_multiplier_dimension=0, full_cut_rank=2)


def symbol_checks():
    _, _, a, _ = m.fixture()
    ident = m.eye(4)
    fixtures = [
        (m.eye(3), [Q(0)]*3),
        ([[2, 1, 0], [0, 3, 1], [1, 0, 2]], [Q(1, 3), Q(-2, 5), Q(1, 7)]),
        ([[Q(1, 2), 0, 0], [0, -2, 0], [1, 0, 3]], [Q(-1), Q(2), Q(0)]),
    ]
    checked = 0
    for v, beta in fixtures:
        require(m.det(v) != 0, 'nonsingular constitutive fixture')
        for tau, *k in product((-1, 0, 1), repeat=4):
            s = tau-m.dot(beta, k)
            vk = m.mv(v, k)
            symbol = m.add(m.scale(s, ident), m.scale(-1, linear_combination(vk, a)))
            q = -s*s+m.dot(vk, vk)
            require(m.det(symbol) == q*q, 'independent principal determinant')
            checked += 1
    singular = [[1, 0, 0], [0, 1, 0], [0, 0, 0]]
    require(m.det(singular) == 0 and m.dot(m.mv(singular, [0, 0, 1]), [0, 0, 1]) == 0,
            'singular coframe control')
    split_form = [[Q(s) if i == j else Q(0) for j in range(4)]
                  for i, s in enumerate([1, 1, -1, -1])]
    n, w = [1, 0, 0, 0], [0, 1, 0, 0]
    aa, bb, cc = m.dot(n, m.mv(split_form, n)), 2*m.dot(n, m.mv(split_form, w)), m.dot(w, m.mv(split_form, w))
    require(bb*bb-4*aa*cc < 0, 'split-signature nonreal-root control')
    speed_roots = []
    for speed in (Q(1, 2), Q(3)):
        symbol = m.add(m.scale(speed, ident), m.scale(-speed, a[0]))
        require(m.det(symbol) == 0, 'freely calibrated characteristic root')
        speed_roots.append(str(speed))
    return dict(constitutive_fixtures=len(fixtures), covectors_checked=checked,
                singular_coframe_rejected=True, split_null_cone_not_hyperbolic=True,
                freely_calibrated_speed_roots=speed_roots)


def geometry_checks():
    rows = []
    for ell, x1, x2 in product((Q(1), Q(2)), (Q(1), Q(2, 3)), (Q(1), Q(3, 2))):
        metric = [[{} for _ in range(4)] for _ in range(4)]
        for a, sign, powers in [(0, -1, (0, -2, 0, 0)), (1, 1, (0, -2, 0, 0)),
                                 (2, 1, (0, 0, -2, 0)), (3, 1, (0, 0, -2, 0))]:
            metric[a][a] = {powers: sign*2*ell*ell}
        g, ric, scalar = m.ricci(metric, [1, x1, x2, 1])
        chi = -1/(2*ell*ell)
        require(ric == m.scale(chi, g), 'direct product Ricci tensor')
        require(scalar == -2/(ell*ell), 'direct product scalar curvature')
        einstein = m.add(ric, m.scale(chi-scalar/2, g))
        require(einstein == m.zero(4), 'vacuum Einstein equation residual')
        rows.append(dict(ell=str(ell), x1=str(x1), x2=str(x2), scalar=str(scalar),
                         cosmological_coefficient=str(chi)))
    gradient = [[{} for _ in range(4)] for _ in range(4)]
    gradient[0][2] = gradient[2][0] = {(-1, 0, -1, 0): -2}
    gradient[1][1] = {(-2, 0, 0, 0): 2}
    gradient[3][3] = {(0, 0, -2, 0): 2}
    g, ric, scalar = m.ricci(gradient, [1, 1, 1, 1])
    require(ric == [[Q(s) if i == j else Q(0) for j in range(4)]
                    for i, s in enumerate([-1, 1, -1, 1])], 'gradient-reflection Ricci control')
    require(scalar == 1 and g[0][0] == 0 and ric[0][0] != 0, 'gradient flip not Einstein')

    # Test the cofactor divergence and Euler identity using independent jets.
    poly = {(2, 0): Q(1, 2), (0, 2): Q(1, 2), (2, 1): Q(1, 2), (0, 4): Q(1, 24)}
    ma_count = 0
    for point in product((Q(0), Q(1, 3), Q(2, 3)), repeat=2):
        d = lambda indices: m.derivative(poly, point, indices)
        cof = lambda i, j, indices=(): (1 if i == j else -1)*d((1-i, 1-j)+indices)
        h = [[d((i, j)) for j in range(2)] for i in range(2)]
        for j in range(2):
            require(sum(cof(i, j, (i,)) for i in range(2)) == 0, 'cofactor divergence')
        twice = sum(d((i, j))*cof(i, j)+d((i,))*cof(i, j, (j,))
                    +d((j,))*cof(i, j, (i,))+d(())*cof(i, j, (i, j))
                    for i, j in product(range(2), repeat=2))
        require((m.det(h)+twice)/3 == m.det(h), 'Monge-Ampere action first variation')
        ma_count += 1
    return dict(lorentzian_product_checks=rows,
                gradient_reflection=dict(ricci_diagonal=['-1', '1', '-1', '1'], scalar='1', einstein=False),
                monge_ampere_jet_checks=ma_count,
                existence_uniqueness='Explicit solution and comparison proof in GS-4; not inferred from samples')


def variable_balance_checks():
    _, _, a, j = m.fixture()
    ident = m.eye(4)
    count = 0
    failed_omission = 0
    for r, s in product(range(-2, 3), range(-1, 2)):
        rho = Q(r*r+s*s+1, 3)
        lt = Q(2*r+1, 5)
        lx = [Q(s+1, 3), Q(r-2, 7), Q(r+s+3, 11)]
        matrices = []
        derivatives = []
        for i in range(3):
            matrices.append(linear_combination([Q(r+i+1, 7)]+[Q((r+2)*(k+1)+s+i, 5)
                                                   for k in range(3)], [ident]+a))
            derivatives.append(linear_combination([Q(s-i, 3)]+[Q(r+k+i+1, 11)
                                                     for k in range(3)], [ident]+a))
        drift = m.scale(lt, ident)
        for i in range(3):
            drift = m.add(drift, m.scale(-1, m.add(derivatives[i], m.scale(lx[i], matrices[i]))))
        drift = m.scale(Q(1, 2), drift)
        c = m.mm(j, linear_combination([Q(r+1, 3), Q(s+2, 5)], [ident, a[0]]))
        require(m.transpose(c) == m.scale(-1, c) and m.comm(j, c) == m.zero(4), 'compatible skew C')
        require(m.comm(j, drift) == m.zero(4), 'variable drift commutes with J')
        psi = [Q(r*k+s+2, k+1) for k in range(4)]
        dx = [[Q((i+1)*(k+2)+s-r, i+k+2) for k in range(4)] for i in range(3)]
        transport = [sum(m.mv(matrices[i], dx[i])[k] for i in range(3)) for k in range(4)]
        dt = [x-y for x, y in zip(transport, m.mv(m.add(drift, c), psi))]
        flux = sum(lx[i]*m.dot(psi, m.mv(matrices[i], psi))/2
                   +m.dot(psi, m.mv(derivatives[i], psi))/2
                   +m.dot(psi, m.mv(matrices[i], dx[i])) for i in range(3))
        time_rate = rho*(lt*m.dot(psi, psi)/2+m.dot(psi, dt))
        require(time_rate-rho*flux == 0, 'independent local conservation jet')
        # Dropping Q keeps a skew C but in general breaks weighted conservation.
        dt_bad = [x-y for x, y in zip(transport, m.mv(c, psi))]
        if rho*(lt*m.dot(psi, psi)/2+m.dot(psi, dt_bad)-flux) != 0:
            failed_omission += 1
        formal_defect = m.scale(-lt, ident)
        for i in range(3):
            formal_defect = m.add(formal_defect, m.add(derivatives[i], m.scale(lx[i], matrices[i])))
        require(m.add(formal_defect, m.add(drift, m.transpose(drift))) == m.zero(4),
                'independent formal adjoint residual')
        count += 1
    require(failed_omission > 0, 'missing drift negative control detected')
    return dict(variable_coefficient_jets=count, missing_drift_failures=failed_omission)


def action_checks():
    _, gamma, a, j = m.fixture()
    d = m.centered_periodic(3)
    dt = m.kron(d, m.eye(3))
    dx = m.kron(m.eye(3), d)
    direction = linear_combination([Q(1), Q(2, 3), Q(-1, 5)], a)
    c = m.scale(Q(2, 7), j)
    b0 = m.add(m.kron(dt, m.eye(4)), m.scale(-1, m.kron(dx, direction)))
    b = m.add(b0, m.kron(m.eye(9), c))
    jbig = m.kron(m.eye(9), j)
    kinetic = m.mm(jbig, b)
    require(m.transpose(b) == m.scale(-1, b), 'periodic formal skew adjoint')
    require(m.transpose(kinetic) == kinetic, 'periodic action Hessian symmetry')
    state = [Q((i*i+3*i) % 17-8, (i % 3)+1) for i in range(36)]

    # Compute the Euler residual directly from neighboring field values,
    # separately from the full action matrix.
    expected = []
    field = [[state[4*(3*t+x):4*(3*t+x+1)] for x in range(3)] for t in range(3)]
    for t, x in product(range(3), repeat=2):
        td = [(u-v)/2 for u, v in zip(field[(t+1) % 3][x], field[(t-1) % 3][x])]
        xd = [(u-v)/2 for u, v in zip(field[t][(x+1) % 3], field[t][(x-1) % 3])]
        residual = [u-v+w for u, v, w in zip(td, m.mv(direction, xd), m.mv(c, field[t][x]))]
        expected.extend(m.mv(j, residual))

    def action(z, k=kinetic):
        return m.dot(z, m.mv(k, z))/2

    for i in range(36):
        plus, minus = list(state), list(state)
        plus[i] += Q(1, 3)
        minus[i] -= Q(1, 3)
        require((action(plus)-action(minus))/Q(2, 3) == expected[i],
                'exact finite action variation versus neighbor equation')

    # Bad multiplier and incompatible mass are detected, not silently fitted.
    bad = m.mm(m.kron(m.eye(9), gamma[0]), b0)
    require(m.transpose(bad) != bad, 'Gamma0 is not a variational multiplier')
    mass_term = m.kron(m.eye(9), m.scale(Q(3, 5), gamma[0]))
    bm = m.add(b0, mass_term)
    km = m.mm(jbig, bm)
    require(m.transpose(bm) == m.scale(-1, bm), 'mass conserves norm')
    require(m.scale(Q(1, 2), m.add(km, m.transpose(km))) == m.mm(jbig, b0),
            'quadratic action erases incompatible mass')
    require(m.mv(mass_term, state) != [0]*36, 'lost mass has nonzero equation residual')

    generator = m.add(m.kron(d, direction), m.scale(-1, m.kron(m.eye(3), c)))
    require(m.transpose(generator) == m.scale(-1, generator), 'spatial norm conservation')
    initial = state[:12]
    e0 = m.dot(initial, initial)/2
    current = initial
    steps = []
    for h in (Q(1, 7), Q(-2, 9), Q(3, 11)):
        next_state = m.cayley(generator, current, h)
        midpoint = [(x+y)/2 for x, y in zip(current, next_state)]
        require([(y-x)/h for x, y in zip(current, next_state)] == m.mv(generator, midpoint),
                'midpoint equation independently checked')
        require(m.dot(next_state, next_state)/2 == e0, 'exact norm conservation')
        steps.append(dict(step=str(h), conserved_norm=str(e0)))
        current = next_state
    gm = m.add(m.kron(d, direction), m.scale(-1, m.kron(m.eye(3), gamma[0])))
    require(m.mm(gm, gm) == m.add(m.kron(m.mm(d, d), m.mm(direction, direction)),
                                m.scale(-1, m.eye(12))), 'gapped generator square')

    # An independent two-state witness distinguishes Hamiltonian from norm.
    hmatrix = m.mm(m.kron(m.eye(3), j), m.kron(d, a[0]))
    p, q = next((p, q) for p, q in product(range(12), repeat=2) if p != q and hmatrix[p][q])
    plus, minus = [Q(0)]*12, [Q(0)]*12
    plus[p] = minus[p] = plus[q] = Q(1)
    minus[q] = Q(-1)
    hplus = m.dot(plus, m.mv(hmatrix, plus))/2
    hminus = m.dot(minus, m.mv(hmatrix, minus))/2
    require(m.dot(plus, plus) == m.dot(minus, minus) and hplus == -hminus and hplus != 0,
            'positive norm is not positive Hamiltonian')
    return dict(periodic_spacetime_sites=9, real_field_components=4,
                independent_action_variations=36, norm_preserving_steps=steps,
                incompatible_mass_erased=True, gapped_generator_square=True,
                hamiltonian_witness=[str(hplus), str(hminus)],
                wrong_multiplier_rejected=True)


def generate():
    require(__debug__, 'Do not use Python -O: recovered checks use assertions.')
    require(sys.version_info[:2] in ((3, 11), (3, 12)), 'Use the declared Python 3.11 or 3.12 environment.')
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for row in pins['local_upstream_sources']:
        data = (REPO/row['path']).read_bytes()
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'upstream SHA256 '+row['path'])
        git_hash = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        require(git_hash == row['git_blob_sha1'], 'upstream git blob '+row['path'])
    for row in pins['recovered_artifacts']:
        require(hashlib.sha256((HERE/row['path']).read_bytes()).hexdigest() == row['sha256'],
                'preserved original bytes '+row['path'])
    old = subprocess.run([sys.executable, str(HERE/'history/verify_bridge.py')],
                         check=True, capture_output=True, text=True)
    recovered = json.loads(old.stdout)
    require(recovered == json.loads((HERE/'history/verification.json').read_text()), 'recovered exact output')
    paths = sorted(str(p.relative_to(HERE)) for p in HERE.rglob('*')
                   if p.is_file() and '__pycache__' not in p.parts and p.suffix in ('.md', '.py', '.json')
                   and p.name != 'CERTIFICATE.json')
    return dict(protocol='UGD_KAHLER_PROPAGATION_R9', status='PASS_EXACT_CONDITIONAL_CONTROLS',
                algebra=algebra_checks(), symbols=symbol_checks(), geometry=geometry_checks(),
                balance=variable_balance_checks(), action=action_checks(), recovered=recovered,
                thermo_gauge=thermo_gauge.run(),
                metric_dynamics=metric_dynamics.run(),
                loop_curvature=loop_curvature.run(),
                spin_matter=spin_matter.run(),
                coupled_cosmology=coupled_cosmology.run(),
                perturbation_stability=perturbation_stability.run(),
                anisotropic_cosmology=anisotropic_cosmology.run(),
                speed_calibration=speed_calibration.run(),
                local_upstream_pins_checked=len(pins['local_upstream_sources']),
                external_pins_scope=pins['external_dependency_check_scope'],
                source_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
                boundaries=dict(universal_proofs='Written in GS/NP/CP/MG/LR/SM/CS/PS/BI/SC notes; coefficient ranks and complete polynomial coefficient checks exhaust their fixed algebraic problems; variable fixtures and truncated series are not universal proofs',
                                physical_signature_selection='NOT_DERIVED_FROM_NATIVE_AXIOMS',
                                full_seam_faithful_representation='NOT_CLAIMED',
                                lambda_to_propagation_coframe='NOT_SELECTED_BY_NATIVE_PROCESS',
                                thermo_to_metric_variation='EXPLICIT_LOCAL_CONSTITUTIVE_ADAPTER_MG5',
                                einstein_hilbert_action='CONDITIONAL_NATIVE_LOOP_LIMIT_LR1_LR4_NOT_PRIMITIVE_SELECTION',
                                loop_process='DECLARED_QUADRATIC_ORIENTED_READOUT_WITH_SUPPLIED_CONNECTION_AND_SCALE',
                                cosmological_relation='Lambda=-12 sigma u^2; kappa Lambda=6/beta; NOT_PHYSICALLY_CALIBRATED',
                                native_readout_and_module_selection='NOT_DERIVED_FROM_PRIMITIVE_LAWS',
                                curved_spin_matter_coupling='CONDITIONAL_DOUBLED_REAL_CLASSICAL_ACTION_SM1_SM6',
                                torsion_elimination='FULL_CONNECTION_SOLUTION_FOR_ZERO_HOLST_WITH_CONTACT_TERM_RETAINED',
                                massive_multiplier_minimality='TWO_COPIES_WITHIN_CONSTANT_SKEW_MULTIPLIERS_ON_FIXED_REAL_FOUR_MODULE_COPIES',
                                physical_particle_masses_and_quantum_statistics='NOT_DERIVED',
                                coupled_einstein_matter_PDE_existence='EXPLICIT_HOMOGENEOUS_NEUTRAL_REST_FLRW_AND_GENERAL_BIANCHI_I_FAMILIES_CS_BI; GENERAL_CAUCHY_PROBLEM_OPEN',
                                cosmological_perturbation_stability='LINEAR_COUPLED_HOMOGENEOUS_BIANCHI_I_PS2_PS3_BI6; SPATIALLY_INHOMOGENEOUS_GRAVITY_OPEN',
                                spatial_matter_stability='GLOBAL_FUTURE_LINEAR_Hr_BOUND_ON_PRESCRIBED_CS_GEOMETRY_PS5_PS6; NOT_FULL_EINSTEIN_MATTER_PERTURBATIONS',
                                homogeneous_matter_completion='ALL_EIGHT_COMPONENT_NONLINEAR_ODE_PS1_BI1_BI2; NONLINEAR_STABILITY_FOR_GENERAL_MATTER_NOT_ESTABLISHED',
                                nonlinear_anisotropic_result='EXACT_REST_BIANCHI_I_FAMILY_WITH_FUTURE_SHEAR_DECAY_BI3_BI6; NOT_GENERAL_NONLINEAR_STABILITY',
                                spin_shear_rotation='CONDITIONAL_NATIVE_STRESS_COMMUTATOR_AND_FIXED_PARALLEL_FRAME_PHASE_RELATION; NOT_QUANTUM_SPIN_OR_EMPIRICAL_CALIBRATION',
                                rest_sector_gradient_closure='FAILS_FOR_EVERY_NONZERO_WAVEVECTOR_PS4',
                                cosmological_benchmark='SUPPLIED_DIMENSIONLESS_PARAMETERS_NOT_PHYSICAL_CALIBRATION',
                                geometric_regularization='NOT_CLAIMED; SINGULAR_PAST_ENDPOINT_WITH_DIVERGENT_TORSION',
                                constants_G_c_Lambda_alpha_hbar='NOT_PREDICTED',
                                vacuum_speed_ratio='c_GW/c_gamma=1_CONDITIONAL_ON_SHARED_METRIC_ACTION_SC6; NOT_PRIMITIVE_ONLY_OR_DISTINCT_FROM_GR',
                                speed_calibration='NATIVE_WAVE_FAMILY_ADMITS_ALL_POSITIVE_SPEEDS; SI_c_IS_DEFINED_EXACT_NOT_A_FIT_TARGET',
                                observational_comparison='SC7_CONSISTENT_WITH_PUBLISHED_GW170817_EMISSION_DEPENDENT_BOUND; NO_RAW_DATA_REANALYSIS_OR_NEW_EMPIRICAL_DISCRIMINATOR',
                                empirical_or_independent_review='LITERATURE_CONSISTENCY_COMPARISON_SC7_ONLY; INDEPENDENT_EXPERIMENTAL_VALIDATION_AND_REVIEW_NOT_PERFORMED',
                                quantum_graviton='NOT_DERIVED; SC4_TWO_CLASSICAL_TENSOR_POLARIZATIONS_ONLY',
                                formal_proof_assistant='NOT_PERFORMED'))


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true', help='Read-only verification (default).')
    mode.add_argument('--write', action='store_true', help='Explicitly regenerate certificate and digest.')
    args = parser.parse_args()
    result = generate()
    data = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(data).hexdigest()
    if args.write:
        (HERE/'CERTIFICATE.json').write_bytes(data)
        (HERE/'EXPECTED.sha256').write_text(digest+'\n')
    else:
        require((HERE/'CERTIFICATE.json').read_bytes() == data, 'Certificate mismatch; nothing overwritten.')
        require((HERE/'EXPECTED.sha256').read_text().strip() == digest, 'Certificate digest mismatch.')
    print(json.dumps(dict(status=result['status'], mode='write' if args.write else 'check', sha256=digest,
                          source_files=len(result['source_sha256']),
                          action_variations=result['action']['independent_action_variations'],
                          thermo_gauge_variations=result['thermo_gauge']['euler_lagrange']['exact_lagrangian_variations'],
                          metric_action_variations=result['metric_dynamics']['einstein']['independent_coframe_variations'],
                          maxwell_metric_variations=result['metric_dynamics']['maxwell_stress']['independent_metric_action_variations'],
                          loop_seed_variations=result['loop_curvature']['action_variation']['independent_seed_coframe_variations'],
                          based_loop_controls=result['loop_curvature']['holonomy']['based_affine_connection_rectangles'],
                          massive_matter_euler_components=result['spin_matter']['matter']['independent_euler_components'],
                          matter_spin_variations=result['spin_matter']['sources']['independent_spin_variations'],
                          torsion_hessian_rank=result['spin_matter']['torsion']['full_connection_hessian_rank'],
                          coupled_solution_coframe_equations=result['coupled_cosmology']['coupled']['coframe_equations'],
                          coupled_solution_cartan_equations=result['coupled_cosmology']['coupled']['cartan_equations'],
                          coupled_solution_matter_equations=result['coupled_cosmology']['coupled']['matter_equations'],
                          perturbation_jacobian_columns=result['perturbation_stability']['linearization']['independent_cubic_jacobian_columns'],
                          homogeneous_metric_tangents=result['perturbation_stability']['coupled']['constrained_metric_tangent_checks'],
                          benchmark_future_rescaled_gain_squared=result['perturbation_stability']['bounds']['benchmark_future_rescaled_gain_squared'],
                          anisotropic_coframe_equations=result['anisotropic_cosmology']['original']['coframe_equations'],
                          anisotropic_stress_forms=result['anisotropic_cosmology']['algebra']['universal_spatial_stress_forms'],
                          homogeneous_shear_mode_pairings=result['anisotropic_cosmology']['coframe']['five_shear_mode_pairings'],
                          speed_ratio=result['speed_calibration']['comparison']['speed_ratio'],
                          speed_observational_comparison=result['speed_calibration']['comparison']['comparison'],
                          independent_principal_ricci_columns=result['speed_calibration']['characteristics']['independent_ricci_columns'],
                          symbol_checks=result['symbols']['covectors_checked']), sort_keys=True))


if __name__ == '__main__':
    main()
