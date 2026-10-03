"""Read-only by default. General results are written proofs, not sample inference."""
import argparse
from fractions import Fraction as Q
import hashlib
from itertools import product
import json
from pathlib import Path

import model as m

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
t = m.thermo


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def all_zero(a):
    return all(x == 0 for row in a for x in row)


def reject(fn, label):
    try:
        fn()
    except ValueError:
        return
    raise RuntimeError('Missing refusal: '+label)


def tensors(coeff, order):
    if order == 3:
        return tuple(t.matrix([[coeff[i+j+k] for k in range(2)] for j in range(2)])
                     for i in range(2))
    return tuple(tuple(t.matrix([[coeff[r+i+j+k] for k in range(2)] for j in range(2)])
                       for i in range(2)) for r in range(2))


def geometry_checks():
    counts = dict(cubic_jets=0, scalar_matrix_curvature_agreements=0,
                  coframe_connection_agreements=0, phase_cut_controls=0,
                  simple_pullback_controls=0)
    metrics = [t.I, t.matrix([[4, 2], [2, 5]]), t.matrix([[9, 3], [3, 2]])]
    fourths = [tensors((0, 0, 0, 0, 0), 4), tensors((1, -2, 3, 2, -1), 4)]
    for h in metrics:
        b = m.coframe(h)
        inv_b = t.inverse(b)
        require(t.multiply(t.transpose(b), b) == h, 'orthonormal coframe metric')
        require(t.multiply(t.multiply(b, t.quarter_turn(h)), inv_b) == t.R,
                'native phase in orthonormal frame')
        for coeff in product((-1, 0, 1), repeat=4):
            dh = tensors(coeff, 3)
            for derivative in dh:
                transformed = t.add(t.multiply(t.multiply(b, t.connection(h, derivative)), inv_b),
                                    t.scale(t.multiply(m.coframe_derivative(h, derivative), inv_b), -1))
                require(transformed == t.scale(t.R, m.scalar_connection(h, derivative)),
                        'scalar phase connection from moving orthonormal frame')
                counts['coframe_connection_agreements'] += 1
            for ddh in fourths:
                f = m.scalar_curvature_direct(h, dh, ddh)
                full = t.curvature_direct(h, dh, ddh)
                require(t.multiply(t.multiply(b, full), inv_b) == t.scale(t.R, f),
                        'scalar differentiated curvature versus full Christoffel route')
                counts['scalar_matrix_curvature_agreements'] += 1
            counts['cubic_jets'] += 1

    h = t.matrix([[2, 1], [1, 5]])
    dh = tensors((0, 1, 0, 0), 3)
    f = m.scalar_curvature_direct(h, dh, fourths[0])
    require([m.scalar_connection(h, x) for x in dh] == [Q(1, 6), Q(-1, 12)],
            'curved thermo scalar potential')
    require(f == Q(-5, 108), 'curved thermo scalar field value')

    reflection = t.matrix([[1, 0], [0, -1]])
    projector = t.matrix([[1, 0], [0, 0]])
    phases = [(Q(1), Q(0)), (Q(-1), Q(0)), (Q(0), Q(1)),
              (Q(3, 5), Q(4, 5)), (Q(-5, 13), Q(12, 13))]
    for c, s in phases:
        g = t.add(t.scale(t.I, c), t.scale(t.R, s))
        require(t.multiply(t.transpose(g), g) == t.I, 'unit phase')
        require((t.commutator(g, reflection) == t.ZERO) == (s == 0),
                'fixed reflection stabilizer is discrete')
        moved = t.multiply(t.multiply(g, reflection), t.transpose(g))
        require(t.multiply(moved, moved) == t.I and
                t.add(t.multiply(moved, t.R), t.multiply(t.R, moved)) == t.ZERO,
                'covariant reflection remains a cut')
        counts['phase_cut_controls'] += 1
    require(t.commutator(projector, t.R) != t.ZERO, 'rank-one readout loses phase')
    require(t.commutator(t.scale(t.R, f), reflection) != t.ZERO,
            'nonzero curvature cannot preserve a parallel reflection')

    vectors = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
               [1, 2, 3, 4], [2, -1, 0, 3], [Q(1, 3), 0, Q(2, 5), 1]]
    for u, v in product(vectors, repeat=2):
        form = m.wedge(u, v)
        require(m.rank(form) <= 2 and m.wedge_square_coefficient(form) == 0,
                'single-channel rank and wedge restriction')
        counts['simple_pullback_controls'] += 1
    f1, f2 = m.wedge(vectors[0], vectors[1]), m.wedge(vectors[2], vectors[3])
    combined = [[x+y for x, y in zip(r, s)] for r, s in zip(f1, f2)]
    require(m.rank(combined) == 4 and m.wedge_square_coefficient(combined) == 2,
            'two-channel full-rank witness')
    return dict(counts=counts, reference=dict(A_s='1/6', A_v='-1/12', f_sv=str(f)),
                rank_four_wedge_square='2', darbouxs_local_theorem='WRITTEN_PROOF_NOT_FINITE_INFERENCE')


def field_checks():
    d0, curl, d2 = m.periodic_complex(2)
    nv, ne, nf, nc = len(d0[0]), len(d0), len(curl), len(d2)
    require(all_zero(m.mm(curl, d0)) and all_zero(m.mm(d2, curl)), 'boundary squared zero')
    ranks = [m.rank(x) for x in (d0, curl, d2)]
    require(ranks == [7, 14, 7], 'periodic cube exact ranks including global modes')

    phases = [m.C.one(), m.C.iota(), m.C(Q(3, 5), Q(4, 5)), m.C(Q(-5, 13), Q(12, 13))]
    links = [phases[(i*i+3*i+1) % len(phases)] for i in range(ne)]
    vertex_phases = [phases[(2*i+1) % len(phases)] for i in range(nv)]
    moved = [u*m.phase_product(row, vertex_phases) for row, u in zip(d0, links)]
    loops = [m.phase_product(row, links) for row in curl]
    require(loops == [m.phase_product(row, moved) for row in curl], 'compact loop gauge invariance')
    require(all(m.phase_product(row, loops) == m.C.one() for row in d2), 'compact cell Bianchi')
    require(any(z != m.C.one() for z in loops), 'nonflat compact witness')

    a = [Q((7*i) % 11-5, 3) for i in range(ne)]
    adot = [Q((3*i) % 7-3, 2) for i in range(ne)]
    phi = [Q(i-3, 5) for i in range(nv)]
    chi = [Q(i*i-4, 7) for i in range(nv)]
    chidot = [Q(2*i-3, 9) for i in range(nv)]
    e = [-x-y for x, y in zip(adot, m.mv(d0, phi))]
    b = m.mv(curl, a)
    moved_a = [x+y for x, y in zip(a, m.mv(d0, chi))]
    moved_adot = [x+y for x, y in zip(adot, m.mv(d0, chidot))]
    moved_phi = [x-y for x, y in zip(phi, chidot)]
    require(m.mv(curl, moved_a) == b and
            [-x-y for x, y in zip(moved_adot, m.mv(d0, moved_phi))] == e,
            'time-dependent gauge transformation')

    # Independent positive catalogue: Y stacks I and a diagonal readout;
    # its recognized and discarded Gram matrices sum to the full weights.
    eye = m.ident(ne)
    diagonal = [[Q((i % 3)+1) if i == j else Q(0) for j in range(ne)] for i in range(ne)]
    catalogue = eye+diagonal
    full = m.mm(m.transpose(catalogue), catalogue)
    lost = m.mm(m.transpose(diagonal), diagonal)
    require(full == [[eye[i][j]+lost[i][j] for j in range(ne)] for i in range(ne)],
            'recognition information ledger supplies positive weights')
    eps = [full[i][i] for i in range(ne)]
    nu = [Q(1+(i % 2+1)**2) for i in range(nf)]
    generator = m.field_generator(curl, eps, nu)
    weights = [1/x for x in eps]+nu
    require(all(weights[i]*generator[i][j]+generator[j][i]*weights[j] == 0
                for i in range(ne+nf) for j in range(ne+nf)), 'energy-skew generator')

    displacement = m.mv(m.transpose(curl), [Q(i % 5-2) for i in range(nf)])
    state = displacement+b
    e0 = m.energy(state, eps, nu)
    require(e0 > 0, 'nonzero initial field energy')
    initial_gauss = m.mv(m.transpose(d0), state[:ne])
    initial_div_b = m.mv(d2, state[ne:])
    require(all(x == 0 for x in initial_gauss+initial_div_b), 'admissible initial constraints')
    trajectory = []
    for step in (Q(1, 7), Q(2, 9), Q(-1, 5)):
        old = state
        state = m.midpoint_step(generator, old, step)
        mid = [(x+y)/2 for x, y in zip(old, state)]
        require([(y-x)/step for x, y in zip(old, state)] == m.mv(generator, mid),
                'independent midpoint equation residual')
        require(m.energy(state, eps, nu) == e0, 'exact midpoint energy conservation')
        require(m.mv(m.transpose(d0), state[:ne]) == initial_gauss and
                m.mv(d2, state[ne:]) == initial_div_b, 'propagated field constraints')
        trajectory.append(dict(step=str(step), energy=str(e0)))

    # Source work/continuity at an independent state, including nonzero charge.
    d = [eps[i]*e[i] for i in range(ne)]
    state_source = d+b
    current = [Q(i % 5-2, 3) for i in range(ne)]
    dy = m.mv(generator, state_source)
    for i in range(ne):
        dy[i] -= current[i]
    energy_rate = m.dot([weights[i]*state_source[i] for i in range(ne+nf)], dy)
    require(energy_rate == -m.dot(e, current), 'source work balance')
    require([-x for x in m.mv(m.transpose(d0), dy[:ne])] == m.mv(m.transpose(d0), current),
            'source charge continuity')

    # Central differences differentiate the quadratic action exactly.
    def action(aa, vv, pp):
        ee = [-x-y for x, y in zip(vv, m.mv(d0, pp))]
        bb = m.mv(curl, aa)
        return (sum(w*x*x for w, x in zip(eps, ee))-
                sum(w*x*x for w, x in zip(nu, bb)))/2 + m.dot(current, aa)-m.dot(rho, pp)
    rho = [-x for x in m.mv(m.transpose(d0), d)]
    expected = [x+y for x, y in zip([-z for z in m.mv(m.transpose(curl),
                 [w*x for w, x in zip(nu, b)])], current)]
    variations = 0
    for slot, vector, gradient in ((0, a, expected), (1, adot, [-x for x in d]),
                                   (2, phi, [Q(0)]*nv)):
        for i in range(len(vector)):
            plus, minus = [list(a), list(adot), list(phi)], [list(a), list(adot), list(phi)]
            plus[slot][i] += 1
            minus[slot][i] -= 1
            require((action(*plus)-action(*minus))/2 == gradient[i], 'independent action variation')
            variations += 1

    # Deliberate wrong signs and a broken incidence detect real failures.
    explicit = [x+Q(1, 7)*y for x, y in zip(displacement+b, m.mv(generator, displacement+b))]
    require(m.energy(explicit, eps, nu) != e0, 'explicit Euler is not certified midpoint transport')
    wrong = [list(row) for row in generator]
    for i in range(ne, ne+nf):
        wrong[i] = [-x for x in wrong[i]]
    require(any(weights[i]*wrong[i][j]+wrong[j][i]*weights[j] != 0
                for i in range(ne+nf) for j in range(ne+nf)), 'wrong Faraday sign breaks energy law')
    broken = [list(row) for row in curl]
    broken[0][0] += 1
    require(not all_zero(m.mm(broken, d0)), 'broken incidence violates gauge invariance')
    return dict(cell_counts=dict(vertices=nv, links=ne, faces=nf, cells=nc), ranks=ranks,
                compact_face_checks=nf, compact_cell_checks=nc,
                independent_action_variations=variations, midpoint_steps=trajectory,
                source_work_rate=str(energy_rate), negative_controls=3)


def propagation_checks():
    phases = [m.C.one(), m.C(-1), m.C.iota(), m.C(Q(3, 5), Q(4, 5))]
    modes = 0
    for z in product(phases, repeat=3):
        delta, curl = m.fourier_curl(z)
        lam = sum(x.norm_square() for x in delta)
        gram = m.mm(m.native_adjoint(curl), curl)
        target = [[m.C(lam if i == j else 0)-delta[i]*delta[j].dagger()
                   for j in range(3)] for i in range(3)]
        require(gram == target, 'Fourier curl symbol Gram identity')
        require(all(x.is_zero() for x in m.mv(curl, delta)), 'longitudinal kernel')
        require(m.rank(m.real_representation(gram)) == (4 if lam else 0),
                'two complex transverse modes; separate constant mode')
        if lam:
            modes += 1
    # Only +/-1 phases are actual characters of the finite side-two witness.
    ranks = []
    for z in product((m.C.one(), m.C(-1)), repeat=3):
        d, curl = m.fourier_curl(z)
        ranks.append(m.rank(m.real_representation(curl))//2)
    require(sum(ranks) == 14, 'finite Fourier rank agrees with independently built incidence')
    eta, gap, speed, kval = Q(5), Q(4), Q(3), Q(1)
    omega = Q(5)
    require(omega*omega == gap*gap+speed*speed*kval*kval, 'declared gapped dispersion fixture')
    energy, momentum = eta*omega, eta*kval
    rest, mass = eta*gap, eta*gap/speed**2
    require(energy**2-speed**2*momentum**2 == mass**2*speed**4 and rest == mass*speed**2,
            'conditional mass-energy algebra')
    return dict(nonzero_native_fourier_characters=modes, zero_mode_checks=1,
                periodic_curl_rank=sum(ranks),
                conditional_mass_fixture=dict(eta=str(eta), gap=str(gap), c=str(speed),
                    E=str(energy), p=str(momentum), E0=str(rest), inertial_mass=str(mass)),
                physical_vacuum_or_mass_selection='NOT_DERIVED')


def generate():
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    for row in pins['local_sources']:
        data = (ROOT/row['path']).read_bytes()
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'source SHA256: '+row['path'])
        require(hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest() == row['git_blob_sha1'],
                'source Git blob: '+row['path'])
    reject(lambda: m.periodic_complex(1), 'collapsed periodic chart')
    reject(lambda: m.field_generator([[Q(1)]], [Q(0)], [Q(1)]), 'zero kinetic weight')
    reject(lambda: m.field_generator([[Q(1)]], [Q(1)], [Q(-1)]), 'negative curvature weight')
    reject(lambda: m.fourier_curl([m.C(2), m.C.one(), m.C.one()]), 'nonunit phase')
    reject(lambda: m.solve([[Q(0)]], [Q(1)]), 'singular solve')
    files = ['README.md', 'THEOREM.md', 'SOURCE_PINS.json', 'model.py', 'verify.py']
    return dict(protocol='THERMO_PHASE_FIELD_R1', status='PASS_DECLARED_EXACT_CONTROLS',
                geometry=geometry_checks(), field=field_checks(), propagation=propagation_checks(),
                refusal_controls=5, pinned_sources_checked=len(pins['local_sources']),
                source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
                scope=dict(proofs='THEOREM.md; finite checks do not replace the proofs',
                           field_dynamics='Declared quadratic phase action and selected cell complex',
                           two_channel_completeness='Local closed nondegenerate two-forms; not global selection',
                           c='Computed model speed; physical universal c not derived',
                           mass='Conditional dispersion/calibration corollary; no particle mass prediction',
                           empirical_validation='NOT_PERFORMED', proof_assistant='NOT_PERFORMED'))


def main():
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = generate()
    data = (json.dumps(result, sort_keys=True, indent=2)+'\n').encode()
    digest = hashlib.sha256(data).hexdigest()
    if args.write:
        (HERE/'CERTIFICATE.json').write_bytes(data)
        (HERE/'EXPECTED.sha256').write_text(digest+'\n')
    else:
        require((HERE/'CERTIFICATE.json').read_bytes() == data, 'Certificate differs; no file overwritten.')
        require((HERE/'EXPECTED.sha256').read_text().strip() == digest, 'Certificate digest differs.')
    print(json.dumps(dict(status=result['status'], sha256=digest,
                         mode='write' if args.write else 'check',
                         geometry=result['geometry']['counts'],
                         field_ranks=result['field']['ranks'],
                         action_variations=result['field']['independent_action_variations'],
                         transverse_mode_controls=result['propagation']['nonzero_native_fourier_characters']),
                     sort_keys=True))


if __name__ == '__main__':
    main()
