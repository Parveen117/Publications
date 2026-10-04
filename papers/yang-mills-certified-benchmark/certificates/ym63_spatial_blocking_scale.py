"""YM63: exact spatial blocking, returning memory and continuum scale gates.

Python 3.12 only. Reuses native polynomials and positive reference moments.
Leading finite-loop coefficients are not a uniform-area or continuum proof.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym62_gauge_lattice_gap as previous

native, p, y, a = previous.native, previous.p, previous.y, previous.a
source, exact = previous.source, previous.exact
RESULT = HERE / 'YM63_RESULT.json'
PIN = HERE / 'EXPECTED_YM63.sha256'
SOURCES = HERE / 'YM63_SOURCE_PINS.json'
MAX_RESPONSE_AREA = 12  # Resource guard, not a mathematical cutoff.


def positive(value):
    value = exact(value)
    if value <= 0:
        raise ValueError('A strictly positive exact scalar is required')
    return value


def scale_map(spacing, g_squared, electric=F(1), magnetic=F(4)):
    spacing, g_squared = positive(spacing), positive(g_squared)
    electric, magnetic = positive(electric), positive(magnetic)
    kappa = electric*g_squared/spacing
    theta = magnetic/(electric*g_squared*g_squared)
    return dict(kappa=kappa,theta=theta,magnetic_coefficient=magnetic/(spacing*g_squared),
                in_YM62_window=theta < F(1, 8000))


def scaled_gap(theta, kappa):
    theta, kappa = exact(theta), positive(kappa)
    gate = previous.block_gate(theta)
    if gate is None:
        return None
    return dict(lower=kappa*gate['gamma'].lo,upper=kappa*F(100, 19),
                exact_lower_expression='kappa*gamma(theta)')


@lru_cache(maxsize=None, typed=True)
def product_coordinates(length):
    native.width(length)
    out = (source.constant(1, length), {}, {}, {})
    for site in range(length):
        coords = tuple(native.lift(f, site, length) for f in p.COORD)
        out = previous.quaternion_product(out, coords)
    return out


def pullback(poly, length):
    coords = product_coordinates(length)
    out = {}
    for monomial, coefficient in poly.items():
        if len(monomial) != 4:
            raise ValueError('A one-link polynomial is required')
        term = source.constant(exact(coefficient), length)
        for coordinate, power in zip(coords, monomial):
            native.width(power, 0)
            for _ in range(power):
                term = p.mul(term, coordinate)
        out = p.add(out, term)
    return out


def blocked_geometry(coarse_side, factor):
    native.width(factor, 2)
    coarse = previous.lattice(coarse_side)
    fine_side = coarse_side*factor
    fine = previous.lattice(fine_side)
    lookup = {edge: i for i, edge in enumerate(fine['edges'])}
    paths = []
    for vertex, axis in coarse['edges']:
        start = tuple(factor*x for x in vertex)
        path = []
        for step in range(factor):
            x = list(start)
            x[axis] = (x[axis]+step)%fine_side
            path.append(lookup[tuple(x), axis])
        paths.append(tuple(path))
    used = [edge for path in paths for edge in path]
    assert len(set(used)) == len(used)
    retained = frozenset(used)
    missing = tuple(tuple(e for e, sign in face if e not in retained) for face in fine['faces'])
    assert all(missing)
    return dict(coarse=coarse,fine=fine,paths=tuple(paths),unretained_per_face=missing)


def rectangle(width, height):
    native.width(width)
    native.width(height)
    edges = tuple((0, i, j) for i in range(width) for j in range(height+1))
    edges += tuple((1, i, j) for i in range(width+1) for j in range(height))
    index = {e: k for k, e in enumerate(edges)}
    faces = tuple(((index[0, i, j], 1), (index[1, i+1, j], 1),
                   (index[0, i, j+1], -1), (index[1, i, j], -1))
                  for i in range(width) for j in range(height))
    outer = tuple((index[0, i, 0], 1) for i in range(width))
    outer += tuple((index[1, width, j], 1) for j in range(height))
    outer += tuple((index[0, i, height], -1) for i in range(width-1, -1, -1))
    outer += tuple((index[1, 0, j], -1) for j in range(height-1, -1, -1))
    return dict(edges=edges,faces=faces,outer=outer)


def boundary_masks(labels):
    labels = tuple(tuple(label) for label in labels)
    if not labels or any(not label for label in labels):
        raise ValueError('Nonempty labelled boundaries are required')
    flat = tuple(e for label in labels for e, sign in label)
    if any(isinstance(e, bool) or not isinstance(e, int) or e < 0 for e in flat):
        raise ValueError('Nonnegative integer links are required')
    previous.incidence(labels, max(flat)+1)
    masks = tuple(sum(1 << e for e, sign in label) for label in labels)
    result = [0]*(1 << len(labels))
    for mask in range(1, len(result)):
        bit = mask & -mask
        result[mask] = result[mask^bit] ^ masks[bit.bit_length()-1]
    return tuple(result)


def inverse_weights(labels, close_last=False):
    """Projected inverse factors; closed proper subsets are refused."""
    boundaries = boundary_masks(labels)
    full = len(boundaries)-1
    weights = [F(0)]*(full+1)
    weights[0] = F(1)
    for mask in range(1, full+1):
        incoming = sum((weights[mask^(1 << j)] for j in range(len(labels))
                        if mask >> j & 1), F(0))
        if close_last and mask == full:
            if boundaries[mask]:
                raise ValueError('The final source surface must close')
            weights[mask] = incoming
        else:
            count = boundaries[mask].bit_count()
            if not count:
                raise ValueError('A proper closed subset needs energy counterterms')
            weights[mask] = incoming/(F(3, 4)*count)
    return tuple(weights), boundaries


@lru_cache(maxsize=None, typed=True)
def loop_response(width, height):
    native.width(width)
    native.width(height)
    area = width*height
    if area > MAX_RESPONSE_AREA:
        raise ValueError('Exact enumeration exceeds the declared resource budget')
    patch = rectangle(width, height)
    faces, outer = patch['faces'], patch['outer']
    wave, boundaries = inverse_weights(faces)
    full = len(wave)-1
    moment = F(1, 4**area)
    bra_ket = moment*sum((wave[s]*wave[full^s] for s in range(full+1)), F(0))
    energy_weights, source_boundaries = inverse_weights(faces+(outer,), close_last=True)
    differentiated_energy = moment*energy_weights[-1]
    assert bra_ket == differentiated_energy > 0
    assert set(i for i, boundary in enumerate(source_boundaries) if not boundary) == {0, len(source_boundaries)-1}
    hist = sorted(Counter(b.bit_count() for b in source_boundaries[1:-1]).items())
    return dict(width=width,height=height,area=area,perimeter=2*(width+height),
                leading=bra_ket,energy_route=differentiated_energy,free_surface_moment=moment,
                fixed_volume_remainder_order=area+2,
                proper_inverse_witnesses=len(source_boundaries)-2,
                frontier_histogram=hist,matching_normalized_coarse_coefficient=6*bra_ket)


def memory_norm(plaquettes, theta):
    native.width(plaquettes, 0)
    theta = exact(theta)
    return theta*theta*plaquettes/4


def center_sign(vertex, axis):
    return (-1)**sum(vertex[:axis])


def parity_controls():
    count = 0
    for side in (2, 3, 4):
        cell = previous.lattice(side, False)
        signs = [center_sign(vertex, axis) for vertex, axis in cell['edges']]
        for word in cell['faces']:
            assert y.prod(signs[e] for e, direction in word) == -1
            count += 1
    for m in range(1, 5):
        for n in range(1, 5):
            patch = rectangle(m, n)
            signs = [center_sign((i, j, 0), axis) for axis, i, j in patch['edges']]
            assert y.prod(signs[e] for e, direction in patch['outer']) == (-1)**(m*n)
    odd = previous.lattice(5)
    odd_signs = [center_sign(vertex, axis) for vertex, axis in odd['edges']]
    assert any(y.prod(odd_signs[e] for e, direction in word) != -1 for word in odd['faces'])
    return dict(open_box_face_flips=count,rectangle_parity_checks=16,
                odd_periodic_extension_refused=True,remainder_order='area+2 at fixed volume and rectangle')


def multiply_average(left, right, sites, width):
    """Haar contractions before expansion, for independent two-face checks."""
    sites = frozenset(sites)
    out = {}
    for m, c in left.items():
        for n, d in right.items():
            exponent = [u+v for u, v in zip(m, n)]
            coefficient = c*d
            for site in sites:
                local = tuple(exponent[4*site:4*site+4])
                coefficient *= y.phi({local:F(1)})
                if not coefficient:
                    break
                exponent[4*site:4*site+4] = [0]*4
            if coefficient:
                term = tuple(exponent)
                out[term] = out.get(term, F(0))+coefficient
    return p.clean(out)


@lru_cache(maxsize=None)
def blocking_controls():
    checks = []
    samples = (p.COORD[0], p.mul(p.COORD[0], p.COORD[1]), p.mul(p.COORD[0], p.COORD[0]))
    for length in (2, 3):
        for j, f in enumerate(samples):
            actual = pullback(f, length)
            target = p.scale(pullback(y.lap(f), length), length)
            assert previous.electric(actual, length) == target
            assert native.phi(actual, length) == y.phi(f)
            checks.append(dict(length=length,sample=j,terms=len(actual)))
        coords = product_coordinates(length)
        assert native.phi(p.mul(coords[0], coords[1]), length) == 0
        assert native.phi(p.mul(coords[0], coords[0]), length) == F(1, 4)
    # Exact ordered endpoint cancellation with three distinct interior gauges.
    qs = y.turns()[:3]
    gs = (y.turns()[4], y.turns()[2], y.turns()[5], y.turns()[0])
    changed = previous.gauge_transform(qs, gs, ((0, 1), (1, 2), (2, 3)))
    word = lambda row: y.qmul(y.qmul(row[0], row[1]), row[2])
    assert word(changed) == y.qmul(y.qmul(gs[0], word(qs)), y.qdagger(gs[-1]))
    rows = []
    for side, factor in ((4, 2), (4, 3), (5, 2)):
        data = blocked_geometry(side, factor)
        rows.append(dict(coarse_side=side,factor=factor,coarse_links=len(data['paths']),
                         fine_links=len(data['fine']['edges']),
                         fine_plaquettes=len(data['fine']['faces']),
                         least_unretained_links_per_face=min(map(len, data['unretained_per_face']))))
    return dict(casimir_pullback_checks=checks,geometry=rows,
                product_reference_isometry=True,gauge_endpoint_cancellation=True,
                free_coarse_clock='factor times fine electric clock')


@lru_cache(maxsize=None)
def contraction_controls():
    patch = rectangle(2, 1)
    n = len(patch['edges'])
    left, right = (previous.loop_polynomial(word, n) for word in patch['faces'])
    outer = previous.loop_polynomial(patch['outer'], n)
    shared = set(e for e, sign in patch['faces'][0]) & set(e for e, sign in patch['faces'][1])
    assert len(shared) == 1
    joined = multiply_average(left, right, shared, n)
    assert joined == p.scale(outer, F(1, 4))
    assert previous.electric(outer, n) == p.scale(outer, F(9, 2))
    # Integrate all links in one product without constructing its huge square.
    all_averaged = multiply_average(outer, outer, range(n), n)
    assert all_averaged == source.constant(F(1, 4), n)
    moment = F(1, 16)
    direct = (2*F(2, 3)/F(9, 2)+F(2, 9))*moment
    assert direct == loop_response(2, 1)['leading'] == F(7, 216)
    # Universal single-link contraction on the native basis.
    for u, v in itertools.product(previous.UNITS, repeat=2):
        aq = [y.qmul(u, e)[0] for e in previous.UNITS]
        qb = [y.qmul(y.qdagger(e), v)[0] for e in previous.UNITS]
        assert sum((x*z for x, z in zip(aq, qb)), F(0))/4 == y.qmul(u, v)[0]/4
    return dict(shared_link_contraction=True,outer_terms=len(outer),outer_rate=F(9, 2),
                outer_free_norm2=F(1, 4),single_link_basis_checks=16,
                two_face_full_polynomial_response=direct)


def response_controls():
    cells = [loop_response(m, n) for m, n in ((1, 1), (1, 2), (1, 3), (1, 4), (2, 2), (2, 3), (3, 3))]
    expected = (F(1, 6), F(7, 216), F(11, 1728), F(73, 58320), F(35, 23328),
                F(57769, 806215680), F(66145860409, 77411255860592640))
    assert tuple(c['leading'] for c in cells) == expected
    assert cells[4]['area'] == cells[3]['area'] and cells[4]['leading'] != cells[3]['leading']
    assert cells[4]['leading'] != cells[0]['leading']**4
    assert cells[4]['leading']/cells[0]['leading']**4 == F(35, 18)
    assert cells[4]['matching_normalized_coarse_coefficient'] == F(35, 3888)
    assert 2*cells[4]['matching_normalized_coarse_coefficient'] == F(35, 1944)
    for m, n in ((1, 2), (1, 3), (2, 3)):
        assert loop_response(m, n)['leading'] == loop_response(n, m)['leading']
    # Boundary counts of the 2x2 patch range over different kinetic denominators.
    assert len(cells[4]['frontier_histogram']) > 1
    return dict(rectangles=cells,shape_information_survives=True,naive_factorization_ratio_2x2=F(35, 18),
                matching_2x2_potential_coefficient=F(35, 1944),
                matching_is_only_one_observable=True,uniform_in_area_remainder=False)


def memory_controls():
    # Genuine native plaquette: A1=-theta W, W is orthogonal to the retained
    # product-link algebra because it contains an unretained fine link.
    w = previous.loop_polynomial()
    assert not source.site_average(w, 1, 4)
    assert native.phi(p.mul(w, w), 4) == F(1, 4)
    rows = [dict(plaquettes=P,theta=theta,return_norm2=memory_norm(P, theta))
            for P in (6, 192, 1536) for theta in (F(0), F(1, 16000), F(-1, 16000))]
    assert all(r['return_norm2'] > 0 for r in rows if r['theta'])
    # An independent rational block fixture checks Schur signs and recurrence.
    C = [[F(2), F(1, 3)], [F(1, 3), F(3)]]
    determinant = C[0][0]*C[1][1]-C[0][1]**2
    assert C[1][1]/determinant == 1/(C[0][0]-C[0][1]**2/C[1][1])
    C2 = a.multiply(C, C)
    assert C2[0][0]-C[0][0]**2 == C[0][1]**2 == F(1, 9)
    return dict(native_unretained_link_projection_zero=True,cells=rows,
                compressed_first_derivative_interaction_zero=True,
                compressed_second_derivative_return_nonzero=True,
                positive_schur_subtraction_fixture=True,
                memory_is_finite_volume_all_content=True,volume_uniform_memory_bound=False)


def scale_controls():
    fixed, compensated = [], []
    for spacing in (F(1), F(1, 2), F(1, 4), F(1, 16)):
        cell = scale_map(spacing, 256)
        assert cell['theta'] == F(1, 16384) and cell['in_YM62_window']
        assert cell['kappa']*spacing == 256
        gap = scaled_gap(cell['theta'], cell['kappa'])
        assert gap['lower'] > 0 and gap['upper']*spacing == F(25600, 19)
        fixed.append(dict(spacing=spacing,**cell,gap_lower=gap['lower'],gap_upper=gap['upper']))
        changed = scale_map(spacing, 256*spacing)
        assert changed['kappa'] == 256 and changed['theta']*spacing*spacing == F(1, 16384)
        compensated.append(dict(spacing=spacing,**changed))
    assert compensated[0]['in_YM62_window'] and not any(c['in_YM62_window'] for c in compensated[1:])
    for electric, magnetic in ((F(1), F(4)), (F(1, 2), F(4)), (F(3, 2), F(7))):
        for spacing in (F(1), F(1, 10)):
            cell = scale_map(spacing, 3, electric, magnetic)
            assert cell['kappa']*cell['theta'] == cell['magnetic_coefficient']
            assert cell['in_YM62_window'] == (9 > 8000*magnetic/electric)
    assert scaled_gap(F(1, 8000), 1) is None
    return dict(convention_example=dict(electric=F(1),magnetic=F(4)),fixed_coupling=fixed,
                fixed_clock_leaves_window=compensated,physical_gap_upper_factor=F(100, 19),
                weak_g_to_zero_exits_certified_window=True,
                no_continuum_theory_exclusion_claim=True)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = native.digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = native.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = native.digest(SOURCES)
    return checked


def run():
    native.previous.runtime_check()
    return native.previous.encode(dict(
        certificate_type='YM63_NATIVE_SPATIAL_BLOCKING_MEMORY_LOOP_RESPONSE_AND_SCALE_GATES',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        blocking=blocking_controls(),contractions=contraction_controls(),
        response=response_controls(),parity=parity_controls(),memory=memory_controls(),scales=scale_controls(),
        claim_status='EXACT_BLOCKING_AND_FINITE_LOOP_RESPONSE__CONTINUUM_FLOW_OPEN',
        proof_scope='Written proofs plus exact finite controls; not a uniform-area or continuum construction',
        model=dict(YM62_operator_retained=True,free_reference_spatial_cut=True,
                   gauge_compatible_product_links=True,exact_free_heat_intertwining=True,
                   interacting_return_memory=True,leading_rectangle_response=True,
                   two_sided_scaled_gap_constraints=True,
                   exact_one_parameter_RG_closure=False,uniform_area_law=False,
                   interacting_reference_conditional_expectation=False,
                   continuum_R4_QFT=False,physical_mass_prediction=False),
        open='Uniform control of generated memory and operators under repeated spatial blocking; continuum trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM63 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM63 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],loop_2x2_leading='35/23328',
                          normalized_one_readout_matching='35/3888',continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
