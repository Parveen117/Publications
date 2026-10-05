"""NSB1 exact nonlinear response geometry controls, Python 3.12 stdlib.

The unchanged NSR1 Taylor engine and independent NT Christoffel calculation
are imported explicitly. Finite rational replays supplement the written
analytic identities; they are not experimental or continuum certificates.
"""

from fractions import Fraction as F
from functools import lru_cache
import importlib.util
import json
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SPEC = importlib.util.spec_from_file_location(
    'nsb1_nsr1_response', ROOT/'papers/native-symmetry-response/certificates/nsr1_response.py')
nsr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(nsr)
Jet, nt = nsr.Jet, nsr.nt
ZERO = [[F(0), F(0)], [F(0), F(0)]]
MODEL = dict(actual_convex_constitutive_potentials=True,
             entire_Hessian_spectral_field_determines_flatness=False,
             Maxwell_speed_field_determines_response_curvature=False,
             exponential_pair_connections_are_bundle_conjugate=True,
             halfplane_flat_and_curved_connections_are_gauge_equivalent=False,
             finite_response_jets_determine_neighborhood_flatness=False,
             mixed_state_protocol_is_response_parallel_transport=False,
             curvature_certifies_electromagnetic_or_Yang_Mills_mass=False,
             finite_replays_are_universal_analytic_proofs=False)


def exact(value):
    if isinstance(value, (float, bool)):
        raise ValueError('Supply exact rational data, not floating point or booleans')
    return F(value)


def positive(value):
    value = exact(value)
    if value <= 0:
        raise ValueError('A positive rational parameter is required')
    return value


def integer(value, minimum=2):
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError('An integer at or above the stated minimum is required')
    return value


def characteristic(h):
    return [F(1), -h[0][0]-h[1][1], nt.det(h)]


def log_normalized(jet):
    """Taylor log(jet/jet(0)); its omitted constant is explicit."""
    base = positive(jet.derivative())
    remainder = jet/base-1
    return sum((F((-1)**(k+1), k)*remainder**k for k in range(1, 5)), Jet(0))


def complex_log_potential_jet(x, y):
    """Degree 2--4 jet of Re(z log z-z), modulo affine terms.

    g^(n)(z)=(-1)^n (n-2)! z^(-(n-1)) for n>=2. The right-half-plane
    branch is single-valued. All higher coefficients are rational at x,y.
    """
    radius2 = x*x+y*y
    inverse = (x/radius2, -y/radius2)
    power, coefficients = (F(1), F(0)), {}
    for n in range(2, 5):
        power = (power[0]*inverse[0]-power[1]*inverse[1],
                 power[0]*inverse[1]+power[1]*inverse[0])
        factor = F((-1)**n*factorial(n-2), factorial(n))
        real, imag = factor*power[0], factor*power[1]
        for j in range(n+1):
            coefficients[n-j, j] = comb(n, j)*(real, -imag, -real, imag)[j % 4]
    return Jet(coefficients)


def halfplane_response(c, x, y, kind='curved'):
    c, x, y = positive(c), exact(x), exact(y)
    if x <= 1/c:
        raise ValueError('The common positive-response chart is x>1/c')
    if kind not in ('flat', 'curved'):
        raise ValueError('Choose flat or curved')
    X, Y = nsr.variables(x, y)
    energy = c*(X*X+Y*Y)/2-X*log_normalized(X)+X
    if kind == 'curved':
        energy += complex_log_potential_jet(x, y)
    row = nsr.analyze(energy)
    radius2 = x*x+y*y
    expected = ([[c-1/x, F(0)], [F(0), c]] if kind == 'flat' else
                [[c-y*y/(x*radius2), y/radius2], [y/radius2, c-x/radius2]])
    curvature = F(0) if kind == 'flat' else 1/(4*c*radius2*(c*x-1)**2)
    assert row['H'] == expected and row['kappa'] == curvature
    row.update(c=c, x=x, y=y, kind=kind, curvature_formula=curvature,
               H_characteristic=characteristic(row['H']),
               eigenvalues=[c, c-1/x], potential_jet_is_modulo_affine_terms=True)
    return row


def halfplane_pair(c, x, y):
    flat = halfplane_response(c, x, y, 'flat')
    curved = halfplane_response(c, x, y, 'curved')
    assert flat['H_characteristic'] == curved['H_characteristic']
    assert flat['F'] == ZERO and curved['F'] != ZERO
    return dict(flat=flat, curved=curved, identical_entire_spectral_field=True,
                selected_connections_gauge_equivalent=False,
                curved_oriented_small_loop_angle_coefficient_squared=
                    curved['kappa']**2*curved['delta'],
                curved_oriented_small_loop_angle_sign=1)


def exponential_response(A, B, sigma):
    A, B = positive(A), positive(B)
    if sigma not in (-1, 1) or isinstance(sigma, bool):
        raise ValueError('sigma must be -1 or +1')
    X, Y = nsr.variables(0, 0)
    def exp_jet(value):
        return sum((value**k/F(factorial(k)) for k in range(5)), Jet(0))
    energy = (1+F(sigma, 2))*X*X+(1-F(sigma, 2))*Y*Y
    energy += A*exp_jet(X+Y)+B*exp_jet(X-Y)
    row = nsr.analyze(energy)
    total, difference = A+B, A-B
    expected = [[2+total+sigma, difference], [difference, 2+total-sigma]]
    delta = 4*(1+A)*(1+B)-1
    curvature = 2*sigma*A*B/delta**2
    assert row['H'] == expected and row['delta'] == delta and row['kappa'] == curvature
    row.update(A=A, B=B, sigma=sigma, curvature_formula=curvature,
               H_characteristic=characteristic(row['H']),
               wave_characteristic=[F(1), -2*difference, -delta],
               impedance_squared=expected[0][0]/expected[1][1],
               bundle_reflection_conjugacy=True,
               global_response_lower_bound=F(1),
               A_and_B_are_positive_bias_exponentials=True,
               potential_jet_is_modulo_affine_terms=True)
    return row


def exponential_pair(A=1, B=1):
    plus, minus = exponential_response(A, B, 1), exponential_response(A, B, -1)
    reflection = [[F(0), F(1)], [F(1), F(0)]]
    assert plus['H_characteristic'] == minus['H_characteristic']
    assert plus['wave_characteristic'] == minus['wave_characteristic']
    assert plus['impedance_squared']*minus['impedance_squared'] == 1
    assert minus['F'] == nt.mul(nt.mul(reflection, plus['F']), reflection)
    return dict(plus=plus, minus=minus, reflection=reflection,
                same_entire_Hessian_and_Maxwell_speed_spectra=True,
                unordered_loop_spectra_separate_this_pair=False,
                calibrated_impedance_and_curvature_separate_this_pair=True,
                connections_bundle_conjugate=True,
                response_metrics_locally_Riemannian_isometric=False)


def scaled_connection_curvature(H, Hx, Hy, q):
    H = [[exact(value) for value in row] for row in H]
    derivatives = [[[exact(value) for value in row] for row in matrix] for matrix in (Hx, Hy)]
    nsr.check_positive(H)
    if any(matrix != nt.tr(matrix) for matrix in derivatives):
        raise ValueError('Derivatives of a symmetric response must be symmetric')
    q = exact(q)
    X = [nt.mul(nt.inv(H), matrix) for matrix in derivatives]
    commutator = nt.sub(nt.mul(X[0], X[1]), nt.mul(X[1], X[0]))
    exterior = nt.scale(commutator, -1)
    connection = [nt.scale(matrix, q) for matrix in X]
    curvature = nt.scale(commutator, q*(q-1))
    residuals = [nt.sub(derivatives[i], nt.add(nt.mul(nt.tr(connection[i]), H),
                    nt.mul(H, connection[i]))) for i in range(2)]
    assert residuals == [nt.scale(matrix, 1-2*q) for matrix in derivatives]
    return dict(q=q, X=X, commutator=commutator,
                Maurer_Cartan_exterior_derivative=exterior, connection=connection,
                F=curvature, metric_compatibility_residuals=residuals)


def dual_flat_control(response):
    H, derivatives = response['H'], response['dH']
    rows = [scaled_connection_curvature(H, *derivatives, q) for q in (0, F(1, 2), 1)]
    assert rows[0]['F'] == rows[2]['F'] == ZERO
    assert rows[1]['F'] == response['F']
    assert rows[1]['metric_compatibility_residuals'] == [ZERO, ZERO]
    duality = [nt.sub(derivatives[i], nt.add(nt.mul(nt.tr(rows[0]['connection'][i]), H),
                nt.mul(H, rows[2]['connection'][i]))) for i in range(2)]
    assert duality == [ZERO, ZERO]
    return dict(primal=rows[0], mean=rows[1], dual=rows[2], duality_residuals=duality,
                mean_of_endpoint_curvatures=ZERO,
                curvature_of_mean_connection=response['F'],
                mean_connection_is_the_Hessian_Levi_Civita_connection=True,
                both_affine_connections_flat=True)


def mixed_protocol(x, y, a, b):
    x, y, a, b = map(exact, (x, y, a, b))
    if 1+x <= 0 or 1+x+a <= 0:
        raise ValueError('The conjugate update requires positive 1+x throughout')
    high_y = y+b/(1+x+a)
    final_y = high_y-b/(1+x)
    points = [[x, y], [x+a, y], [x+a, high_y], [x, high_y], [x, final_y]]
    margins = [1+first-second*second for first, second in points]
    if min(margins) <= 0:
        raise ValueError('The complete protocol must stay inside 1+x-y^2>0')
    X, Y = nsr.variables(x, y)
    response = nsr.analyze((X*X+(1+X)*Y*Y)/2)
    shift = -a*b/((1+x)*(1+x+a))
    assert final_y-y == shift
    assert response['kappa'] == 1/(4*(1+x-y*y)**2)
    return dict(x=x, y=y, a=a, b=b, points=points,
                displacement=[F(0), shift], bracket=[F(0), -1/(1+x)**2],
                stable_path_margin=min(margins), response=response,
                protocol='+a X, +b Y, -a X, -b Y',
                X='partial_x at fixed y', Y='(1+x)^(-1) partial_y at fixed x',
                conjugate_variable='q=U_y=(1+x)y',
                endpoint_return_is_parallel_transport=False,
                all_continuous_segments_stable_by_endpoint_bounds=True)


def radial_jet_order(m):
    m = integer(m)
    return dict(m=m, first_energy_degree=2*m, first_response_jet_order=2*m-2,
                matched_response_jet_through=2*m-3, first_curvature_order=4*m-6)


def radial_response(m, eta, x, y):
    m, eta, x, y = integer(m), positive(eta), exact(x), exact(y)
    X, Y = nsr.variables(x, y)
    radius2 = X*X+Y*Y
    row = nsr.analyze(radius2/2+eta*radius2**m/F(2*m))
    row.update(m=m, eta=eta, x=x, y=y, global_response_lower_bound=F(1))
    return row


def radial_curvature(m, eta, t):
    m, eta, t = integer(m), positive(eta), exact(t)
    power = t**(2*m-2)
    return -2*(m-1)**3*eta**2*t**(4*m-6)/((1+(2*m-1)*eta*power)**2*(1+eta*power)**2)


def radial_blindness_control(m, eta, t):
    row = radial_response(m, eta, t, 0)
    formula = radial_curvature(m, eta, t)
    assert row['kappa'] == formula
    origin = radial_response(m, eta, 0, 0)
    assert origin['H'] == [[F(1), F(0)], [F(0), F(1)]] and origin['F'] == ZERO
    return dict(response=row, origin=origin, jet_orders=radial_jet_order(m),
                curvature_formula=formula, identical_finite_jet_requires_declared_order=True,
                degree_four_replay_alone_proves_arbitrary_order_vanishing=False)


@lru_cache(maxsize=None)
def all_controls():
    halfplanes = [halfplane_pair(c, x, y) for c, x, y in
                  ((1, 2, 0), (1, 2, 1), (1, F(3, 2), -2), (2, 1, F(1, 3)))]
    anchor = halfplanes[0]
    assert anchor['curved']['kappa'] == F(1, 16)
    assert anchor['curved']['F'] == [[F(0), F(1, 32)], [F(-1, 16), F(0)]]
    exponential = [exponential_pair(A, B) for A, B in
                   ((1, 1), (2, F(1, 2)), (F(3, 2), F(5, 3)), (3, 4))]
    assert exponential[0]['plus']['kappa'] == F(2, 225)
    return dict(halfplane_isospectral_flat_and_curved=halfplanes,
                exponential_same_wave_spectrum_opposite_curvature=exponential,
                dual_flat_connections_with_curved_mean=dual_flat_control(anchor['curved']),
                mixed_control_protocols=[mixed_protocol(*values) for values in
                    ((0, 0, F(1, 4), F(1, 4)), (1, F(1, 4), F(-1, 2), F(1, 5)),
                     (0, 0, 0, F(1, 4)), (0, 0, F(1, 4), 0))],
                radial_finite_jet_blindness=[radial_blindness_control(m, eta, t)
                    for m, eta, t in ((2, 1, 1), (3, F(1, 2), 1),
                                      (5, 2, F(1, 2)), (8, F(1, 3), F(2, 3)))],
                model=MODEL)


if __name__ == '__main__':
    nsr.require_runtime()
    print(json.dumps(nsr.encode(all_controls()), sort_keys=True))
