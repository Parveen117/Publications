"""Source-bound exact controls. Default/check mode does not rewrite evidence."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import subprocess
import sys

import model as m


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def require(value, label):
    if not value:
        raise RuntimeError(label)


def reject(call, label):
    try:
        call()
    except ValueError:
        return
    raise RuntimeError('Missing refusal: '+label)


def scalar_zero(a):
    return a == m.ZERO


def cubic(values):
    return tuple(m.matrix([[values[i+j+k] for k in range(2)]
                          for j in range(2)]) for i in range(2))


def quartic(values):
    return tuple(tuple(m.matrix([[values[r+i+j+k] for k in range(2)]
                                for j in range(2)]) for i in range(2))
                 for r in range(2))


def exact_checks():
    counts = dict(stable_hessians=0, normalized_phase_planes=0,
                  cubic_curvature_fixtures=0, direct_curvature_comparisons=0,
                  exact_transport_endpoints=0, refusal_controls=0,
                  negative_controls=0)
    change = m.matrix([[2, 1], [0, 3]])
    change_inv = m.inverse(change)
    dh = m.matrix([[1, 2], [2, 3]])
    for a, b, c in product(range(1, 9), range(-3, 4), range(1, 9)):
        if a*c-b*b <= 0:
            continue
        h = m.matrix([[a, b], [b, c]])
        d = m.det(h)
        r = m.multiply(m.R, h)
        require(m.multiply(r, r) == m.scale(m.I, -d), 'general determinant identity')
        require(scalar_zero(m.add(m.multiply(m.transpose(r), h), m.multiply(h, r))),
                'general H-skew identity')
        response = m.responses(h, 10, 3)
        require(response['C_P']/response['C_V'] == response['K_S']/response['K_T'],
                'thermal/mechanical closure')
        require(1-response['C_V']/response['C_P'] == Q(b*b, a*c), 'cross-response ratio')
        counts['stable_hessians'] += 1
        try:
            j = m.quarter_turn(h)
        except ValueError:
            continue
        require(m.multiply(j, j) == m.scale(m.I, -1), 'quarter-turn square')
        require(m.multiply(m.multiply(m.transpose(j), h), j) == h, 'metric isometry')
        h_new = m.multiply(m.multiply(m.transpose(change_inv), h), change_inv)
        require(m.quarter_turn(h_new) == m.multiply(m.multiply(change, j), change_inv),
                'oriented affine covariance')
        require(m.quarter_turn(m.scale(h, 7)) == j, 'constant energy-scale covariance')
        dj = m.quarter_turn_derivative(h, dh)
        require(scalar_zero(m.add(m.multiply(j, dj), m.multiply(dj, j))),
                'differentiated quarter-turn identity')
        ai = m.connection(h, dh)
        require(scalar_zero(m.add(dj, m.commutator(ai, j))), 'parallel iota')
        counts['normalized_phase_planes'] += 1

    # Independent ideal-gas response reference: C_V=3/2, R=1, T=2, V=3.
    gas = m.responses(m.matrix([[Q(4, 3), Q(-4, 9)],
                                [Q(-4, 9), Q(10, 27)]]), 2, 3)
    require(gas == dict(C_V=Q(3, 2), C_P=Q(5, 2), K_S=Q(10, 9), K_T=Q(2, 3)),
            'ideal-gas constrained-response reference')

    # Each symmetric cubic/quartic tensor is realized by a Taylor polynomial.
    # The full Christoffel derivative route retains fourth derivatives.
    hessians = [m.matrix([[2, 1], [1, 5]]), m.matrix([[3, 0], [0, 3]]),
                m.matrix([[1, Q(1, 3)], [Q(1, 3), 1]])]
    fourths = [quartic((0, 0, 0, 0, 0)), quartic((2, -1, 3, 1, -2)),
               quartic((1, 2, 1, 2, 1))]
    for h in hessians:
        for coefficients in product((-1, 0, 2), repeat=4):
            derivatives = cubic(coefficients)
            hs, hv = derivatives
            a_conn = (m.connection(h, hs), m.connection(h, hv))
            require(a_conn == m.christoffel_direct(h, derivatives), 'full Christoffel agreement')
            for i in range(2):
                require(m.add(m.multiply(m.transpose(a_conn[i]), h),
                              m.multiply(h, a_conn[i])) == derivatives[i], 'metric compatibility')
            f = m.curvature(h, hs, hv)
            require(scalar_zero(m.add(m.multiply(m.transpose(f), h), m.multiply(h, f))),
                    'curvature metric skew-adjointness')
            for fourth in fourths:
                require(f == m.curvature_direct(h, derivatives, fourth),
                        'curvature from direct differentiation, including quartic cancellation')
                counts['direct_curvature_comparisons'] += 1
            counts['cubic_curvature_fixtures'] += 1

    h = hessians[0]
    hs, hv = cubic((0, 1, 0, 0))
    f = m.curvature(h, hs, hv)
    expected_f = m.matrix([[Q(5, 324), Q(25, 324)], [Q(-5, 162), Q(-5, 324)]])
    require(f == expected_f and m.gaussian_curvature(h, hs, hv) == Q(5, 324),
            'curved polynomial exact value and sign')
    require(f == m.scale(m.quarter_turn(h), -Q(5, 108)), '2D curvature/quarter-turn orientation')
    require(m.curvature(h, m.ZERO, m.ZERO) == m.ZERO, 'quadratic zero-curvature control')
    require(f != m.ZERO and f != m.scale(expected_f, -1), 'detect dropped or reversed commutator')
    counts['negative_controls'] += 2
    wrong_j = m.scale(m.multiply(m.R, h), 1/m.det(h))
    require(m.multiply(wrong_j, wrong_j) != m.scale(m.I, -1), 'wrong determinant normalization')
    j = m.quarter_turn(h)
    wrong_dj = m.scale(m.multiply(m.R, dh), 1/m.rational_sqrt(m.det(h)))
    require(not scalar_zero(m.add(m.multiply(j, wrong_dj), m.multiply(wrong_dj, j))),
            'missing determinant derivative is detected')
    counts['negative_controls'] += 2

    # An actual nonconstant parallel transport, not a frozen-metric rotation.
    for s in (Q(-1, 2), Q(0), Q(1, 3), Q(2)):
        h = m.matrix([[(1+s)**2, 0], [0, 1]])
        hs = m.matrix([[2*(1+s), 0], [0, 0]])
        u = m.matrix([[1/(1+s), 0], [0, 1]])
        du = m.matrix([[-1/(1+s)**2, 0], [0, 0]])
        require(scalar_zero(m.add(du, m.multiply(m.connection(h, hs), u))), 'transport equation')
        require(m.multiply(m.multiply(m.transpose(u), h), u) == m.I, 'transported metric')
        require(m.multiply(u, m.R) == m.multiply(m.quarter_turn(h), u), 'transported iota')
        counts['exact_transport_endpoints'] += 1

    bad_h = [m.matrix([[1, 0], [0, 0]]), m.matrix([[1, 2], [2, 1]]),
             m.matrix([[1, 1], [0, 1]])]
    for h in bad_h:
        reject(lambda h=h: m.quarter_turn(h), 'invalid response metric')
        counts['refusal_controls'] += 1
    reject(lambda: m.responses(m.I, 0, 1), 'zero temperature')
    reject(lambda: m.responses(m.I, 1, 0), 'zero volume')
    reject(lambda: m.quarter_turn(m.matrix([[1, 0], [0, 2]])), 'unrepresented irrational root')
    counts['refusal_controls'] += 3
    wrong_hs, wrong_hv = m.ZERO, m.matrix([[1, 0], [0, 0]])
    reject(lambda: m.curvature(m.I, wrong_hs, wrong_hv), 'nonintegrable cubic jet')
    counts['refusal_controls'] += 1
    require(m.christoffel_direct(m.I, (wrong_hs, wrong_hv)) !=
            (m.connection(m.I, wrong_hs), m.connection(m.I, wrong_hv)),
            'Hessian simplification really needs cubic symmetry')
    counts['negative_controls'] += 1

    return dict(status='PASS_THERMO_COMPASS_EXACT_CONTROLS', counts=counts,
                curved_fixture=dict(H=[['2', '1'], ['1', '5']], K='5/324'),
                ideal_gas={k: str(v) for k, v in gas.items()},
                scope='Written TC-1--TC-6 proofs plus finite exact controls; no physical constant or gravity validation.')


def recovered_checks():
    helical_path = HERE.parent/'uncut-cut-measurement/certificates/native_helical_response.py'
    spec = importlib.util.spec_from_file_location('recovered_native_helical', helical_path)
    helical = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = helical
    spec.loader.exec_module(helical)
    h = m.matrix([[2, 1], [1, 5]])
    j = m.quarter_turn(h)
    def rho(z):
        return m.add(m.scale(m.I, z.radial), m.scale(j, z.turn))
    values = [helical.C(a, b) for a, b in product((Q(-1), Q(0), Q(2, 3)), repeat=2)]
    for z, w in product(values, repeat=2):
        require(rho(z*w) == m.multiply(rho(z), rho(w)), 'native scalar multiplication intertwiner')
    for z in values:
        adj = m.multiply(m.multiply(m.inverse(h), m.transpose(rho(z))), h)
        require(adj == rho(z.dagger()), 'native scalar dagger intertwiner')
    alpha = subprocess.run([sys.executable, str(HERE/'downstream/fine_structure_controls.py')],
                           check=True, capture_output=True, text=True)
    return dict(helical=helical.run(), native_iota_intertwiner_pairs=len(values)**2,
                native_iota_dagger_checks=len(values), downstream_alpha=json.loads(alpha.stdout))


def generate():
    require(__debug__, 'Run without Python -O: recovered assertion-based controls must remain active.')
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    checked = 0
    for row in pins['local_upstream_sources']:
        data = (REPO/row['path']).read_bytes()
        actual = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        require(actual == row['git_blob_sha1'], 'upstream source blob: '+row['path'])
        require(hashlib.sha256(data).hexdigest() == row['sha256'], 'upstream source SHA256')
        checked += 1
    sources = sorted(set([
        'README.md', 'THEOREM.md', 'RESULTS_INDEX.md', 'SOURCE_PINS.json',
        'RECOVERY_MANIFEST.json', 'model.py', 'verify.py',
        'downstream/FINE_STRUCTURE_R1.md', 'downstream/fine_structure_controls.py',
        '../uncut-cut-measurement/NATIVE_IOTA_HELICAL_RESPONSE.md',
        '../uncut-cut-measurement/NATIVE_HELICAL_SOURCE_PINS.json',
        '../uncut-cut-measurement/certificates/native_helical_response.py',
        '../uncut-cut-measurement/certificates/native_sources/native_summability.py',
    ] + [str(Path('..')/Path(row['path']).relative_to('papers'))
         for row in pins['local_upstream_sources']]))
    return dict(protocol='THERMO_COMPASS_FOUNDATIONS_R1',
                status='PASS_DECLARED_EXACT_AND_SOURCE_BOUND_CONTROLS',
                thermo=exact_checks(), recovered=recovered_checks(),
                upstream_local_pins_checked=checked,
                external_dependency_note=pins['external_dependency_check_scope'],
                source_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
                boundaries=dict(universal_proof='THEOREM.md, not finite sampling',
                                lean_or_other_proof_assistant='NOT_PERFORMED',
                                independent_review='NOT_CLAIMED',
                                empirical_validation='NOT_PERFORMED',
                                alpha_hbar_c_G='NOT_DERIVED',
                                geometry='thermodynamic response metric, not identified with spacetime'))


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true', help='Explicitly regenerate approved evidence.')
    mode.add_argument('--check', action='store_true', help='Read-only verification (default).')
    args = parser.parse_args()
    result = generate()
    data = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(data).hexdigest()
    if args.write:
        (HERE/'CERTIFICATE.json').write_bytes(data)
        (HERE/'EXPECTED.sha256').write_text(digest+'\n')
    else:
        require((HERE/'CERTIFICATE.json').read_bytes() == data, 'Certificate differs; no file was overwritten.')
        require((HERE/'EXPECTED.sha256').read_text().strip() == digest, 'Certificate digest differs.')
    print(json.dumps(dict(status=result['status'], counts=result['thermo']['counts'],
                          bound_sources=len(result['source_sha256']), sha256=digest,
                          mode='write' if args.write else 'check'), sort_keys=True))


if __name__ == '__main__':
    main()
