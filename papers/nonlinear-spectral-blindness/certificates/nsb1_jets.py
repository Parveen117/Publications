"""Exact finite Taylor-jet and nonlinear-response controls for NSB1.

Python 3.12, standard library only. Polynomial coefficients multiply monomials
directly: they are Taylor coefficients, not unnormalized derivatives. Matrix
columns are images of basis vectors. Characteristic coefficients are returned
in descending order. The finite replays do not prove the manuscript's general
theorems, an infinite Koopman-spectrum statement, or a quantum mass gap.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
import json
import sys


MODEL = dict(
    finite_polynomial_controls_use_exact_rationals=True,
    pullback_requires_a_fixed_base_point=True,
    all_jet_characteristic_data_determine_nonlinear_coefficients=False,
    all_similarity_invariants_are_blind=False,
    coordinate_images_recover_the_declared_finite_jet=True,
    minimum_probe_count_uses_an_unrestricted_coefficient_chart_and_C1_decoder=True,
    arbitrary_set_theoretic_encodings_obey_the_probe_count=False,
    smooth_germs_determined_by_all_Taylor_coefficients=False,
    infinite_Koopman_spectrum_identified=False,
    physical_Yang_Mills_mass_gap_proved=False,
    finite_controls_are_a_universal_proof=False,
)


def runtime_check():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError('NSB1 controls require Python 3.12')


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, F)):
        raise ValueError('An exact integer or Fraction is required')
    return F(value)


def width(value, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError('Invalid nonnegative integer dimension or order')
    return value


def _poly(poly, dimension=None):
    if not isinstance(poly, dict):
        raise ValueError('A polynomial must be an exponent-to-coefficient dictionary')
    result = {}
    for exponent, coefficient in poly.items():
        if not isinstance(exponent, tuple) or not exponent:
            raise ValueError('Monomial exponents must be nonempty tuples')
        for power in exponent:
            width(power)
        if dimension is None:
            dimension = len(exponent)
        if len(exponent) != dimension:
            raise ValueError('Polynomial dimensions do not match')
        coefficient = exact(coefficient)
        if coefficient:
            result[exponent] = coefficient
    return result


def _dimension(*polys):
    dimensions = {len(exponent) for poly in polys for exponent in poly}
    if len(dimensions) > 1:
        raise ValueError('Polynomial dimensions do not match')
    return next(iter(dimensions), None)


@lru_cache(maxsize=None, typed=True)
def monomials(dimension, order):
    """Degree-ascending, then lexicographic, monomial exponent tuples."""
    width(dimension, 1)
    width(order)
    def homogeneous(degree, size):
        if size == 1:
            yield (degree,)
        else:
            for first in range(degree+1):
                for tail in homogeneous(degree-first, size-1):
                    yield (first,)+tail
    return tuple(exponent for degree in range(order+1)
                 for exponent in homogeneous(degree, dimension))


def coordinate_polynomials(dimension):
    width(dimension, 1)
    return tuple({tuple(int(i == j) for j in range(dimension)): F(1)}
                 for i in range(dimension))


def poly_add(*polys):
    dimension = _dimension(*polys)
    result = {}
    for poly in polys:
        for exponent, coefficient in _poly(poly, dimension).items():
            result[exponent] = result.get(exponent, F(0))+coefficient
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def poly_scale(poly, coefficient):
    coefficient = exact(coefficient)
    return {exponent: coefficient*value for exponent, value in _poly(poly).items()
            if coefficient*value}


def poly_mul(first, second, order=None):
    if order is not None:
        width(order)
    dimension = _dimension(first, second)
    first, second = _poly(first, dimension), _poly(second, dimension)
    result = {}
    for left, a in first.items():
        for right, b in second.items():
            exponent = tuple(x+y for x, y in zip(left, right))
            if order is None or sum(exponent) <= order:
                result[exponent] = result.get(exponent, F(0))+a*b
    return {exponent: value for exponent, value in result.items() if value}


def poly_pow(poly, power, order=None, dimension=None):
    width(power)
    if order is not None:
        width(order)
    poly = _poly(poly, dimension)
    dimension = _dimension(poly) if dimension is None else dimension
    if dimension is None:
        if power:
            return {}
        raise ValueError('The zeroth power of a zero polynomial needs its dimension')
    width(dimension, 1)
    result = {(0,)*dimension: F(1)}
    while power:
        if power % 2:
            result = poly_mul(result, poly, order)
        power //= 2
        if power:
            poly = poly_mul(poly, poly, order)
    return result


def poly_derivative(poly, axis):
    width(axis)
    poly = _poly(poly)
    dimension = _dimension(poly)
    if dimension is not None and axis >= dimension:
        raise ValueError('Derivative axis outside the polynomial dimension')
    result = {}
    for exponent, coefficient in poly.items():
        if exponent[axis]:
            new = list(exponent)
            new[axis] -= 1
            result[tuple(new)] = coefficient*exponent[axis]
    return result


def poly_evaluate(poly, point):
    point = tuple(map(exact, point))
    if not point:
        raise ValueError('An evaluation point must have positive dimension')
    result = F(0)
    for exponent, coefficient in _poly(poly, len(point)).items():
        for value, power in zip(point, exponent):
            coefficient *= value**power
        result += coefficient
    return result


def poly_compose(poly, images, order=None):
    images = tuple(images)
    dimension = len(images)
    width(dimension, 1)
    if order is not None:
        width(order)
    images = tuple(_poly(image, dimension) for image in images)
    result = {}
    for exponent, coefficient in _poly(poly, dimension).items():
        term = {(0,)*dimension: coefficient}
        for image, power in zip(images, exponent):
            term = poly_mul(term, poly_pow(image, power, order, dimension), order)
        result = poly_add(result, term)
    return result


def _field(images, fixed_point=False):
    images = tuple(images)
    dimension = len(images)
    width(dimension, 1)
    images = tuple(_poly(poly, dimension) for poly in images)
    if fixed_point and any(poly.get((0,)*dimension, 0) for poly in images):
        raise ValueError('A jet endomorphism at zero requires a fixed base point')
    return images


def linear_part(images):
    images = _field(images)
    dimension = len(images)
    units = tuple(next(iter(poly)) for poly in coordinate_polynomials(dimension))
    return [[poly.get(unit, F(0)) for unit in units] for poly in images]


def linear_map_polynomials(matrix):
    matrix = _square(matrix)
    coordinates = coordinate_polynomials(len(matrix))
    return tuple(poly_add(*(poly_scale(poly, coefficient)
                            for poly, coefficient in zip(coordinates, row)))
                 for row in matrix)


def _matrix(matrix):
    matrix = [list(map(exact, row)) for row in matrix]
    if matrix and any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError('Matrix rows must have equal lengths')
    return matrix


def _square(matrix):
    matrix = _matrix(matrix)
    if any(len(row) != len(matrix) for row in matrix):
        raise ValueError('A square matrix is required')
    return matrix


def identity(size):
    width(size)
    return [[F(i == j) for j in range(size)] for i in range(size)]


def matrix_add(first, second, weight=1):
    first, second, weight = _matrix(first), _matrix(second), exact(weight)
    if len(first) != len(second) or any(len(a) != len(b) for a, b in zip(first, second)):
        raise ValueError('Matrix dimensions do not match')
    return [[x+weight*y for x, y in zip(a, b)] for a, b in zip(first, second)]


def matrix_mul(first, second):
    first, second = _matrix(first), _matrix(second)
    if not first or not second:
        if not first and not second:
            return []
        raise ValueError('Empty matrix dimensions do not match')
    if len(first[0]) != len(second):
        raise ValueError('Matrix product dimensions do not match')
    result = [[F(0) for _ in second[0]] for _ in first]
    for i, row in enumerate(first):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(second[k]):
                    if other:
                        result[i][j] += value*other
    return result


def matrix_power(matrix, power):
    matrix = _square(matrix)
    width(power)
    result = identity(len(matrix))
    while power:
        if power % 2:
            result = matrix_mul(result, matrix)
        power //= 2
        if power:
            matrix = matrix_mul(matrix, matrix)
    return result


def matrix_rank(matrix):
    matrix = _matrix(matrix)
    if not matrix:
        return 0
    pivot_row = 0
    for column in range(len(matrix[0])):
        pivot = next((row for row in range(pivot_row, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        value = matrix[pivot_row][column]
        matrix[pivot_row] = [entry/value for entry in matrix[pivot_row]]
        for row in range(pivot_row+1, len(matrix)):
            value = matrix[row][column]
            if value:
                matrix[row] = [entry-value*base for entry, base in zip(matrix[row], matrix[pivot_row])]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def matrix_charpoly(matrix):
    """Exact Newton identities; output coefficients of z^n,...,z^0."""
    matrix = _square(matrix)
    size = len(matrix)
    coefficients, traces = [F(1)], [F(0)]
    power = identity(size)
    for degree in range(1, size+1):
        power = matrix_mul(power, matrix)
        traces.append(sum((power[i][i] for i in range(size)), F(0)))
        coefficients.append(-sum((coefficients[degree-k]*traces[k]
                                  for k in range(1, degree+1)), F(0))/degree)
    return tuple(coefficients)


def word_matrix(matrices, word):
    matrices = tuple(_square(matrix) for matrix in matrices)
    if not matrices or any(len(matrix) != len(matrices[0]) for matrix in matrices):
        raise ValueError('A nonempty, equally sized operator family is required')
    result = identity(len(matrices[0]))
    for index in word:
        width(index)
        if index >= len(matrices):
            raise ValueError('Word index outside the operator family')
        result = matrix_mul(result, matrices[index])
    return result


def matrix_polynomial(matrices, terms):
    """Terms is a dictionary {operator-index word: exact scalar coefficient}."""
    matrices = tuple(matrices)
    unit = word_matrix(matrices, ())
    result = [[F(0) for _ in row] for row in unit]
    for word, coefficient in terms.items():
        if not isinstance(word, tuple):
            raise ValueError('Noncommutative monomials must be index tuples')
        result = matrix_add(result, word_matrix(matrices, word), exact(coefficient))
    return result


def jet_pullback(images, order):
    """The contravariant convention is C_F C_G = C_(G composed with F)."""
    images = _field(images, fixed_point=True)
    basis = monomials(len(images), order)
    index = {exponent: i for i, exponent in enumerate(basis)}
    matrix = [[F(0) for _ in basis] for _ in basis]
    for column, exponent in enumerate(basis):
        for output, coefficient in poly_compose({exponent: F(1)}, images, order).items():
            matrix[index[output]][column] = coefficient
    return matrix


def preserves_degree_filtration(matrix, dimension, order):
    matrix = _square(matrix)
    basis = monomials(dimension, order)
    if len(matrix) != len(basis):
        raise ValueError('Matrix dimension does not match the declared jet space')
    return all(not matrix[i][j] for i, left in enumerate(basis)
               for j, right in enumerate(basis) if sum(left) < sum(right))


def graded_blocks(matrix, dimension, order):
    matrix = _square(matrix)
    if not preserves_degree_filtration(matrix, dimension, order):
        raise ValueError('The operator does not preserve the degree filtration')
    basis = monomials(dimension, order)
    return tuple([[matrix[i][j] for j, right in enumerate(basis) if sum(right) == degree]
                  for i, left in enumerate(basis) if sum(left) == degree]
                 for degree in range(order+1))


def recover_coordinate_jet(matrix, dimension, order):
    width(order, 1)
    matrix = _square(matrix)
    basis = monomials(dimension, order)
    if len(matrix) != len(basis):
        raise ValueError('Matrix dimension does not match the declared jet space')
    index = {exponent: i for i, exponent in enumerate(basis)}
    return tuple({exponent: matrix[row][index[next(iter(coordinate))]]
                  for row, exponent in enumerate(basis)
                  if matrix[row][index[next(iter(coordinate))]]}
                 for coordinate in coordinate_polynomials(dimension))


def nonlinear_probe_count(dimension, order):
    """Free nonlinear coefficients; minimal probes assume a local C1 decoder."""
    width(dimension, 1)
    width(order)
    return dimension*sum(comb(dimension+degree-1, degree) for degree in range(2, order+1))


def lie_derivative(poly, vector_field, order=None):
    vector_field = _field(vector_field)
    poly = _poly(poly, len(vector_field))
    return poly_add(*(poly_mul(component, poly_derivative(poly, axis), order)
                      for axis, component in enumerate(vector_field)))


def jet_lie_generator(vector_field, order):
    vector_field = _field(vector_field, fixed_point=True)
    basis = monomials(len(vector_field), order)
    index = {exponent: i for i, exponent in enumerate(basis)}
    matrix = [[F(0) for _ in basis] for _ in basis]
    for column, exponent in enumerate(basis):
        for output, coefficient in lie_derivative({exponent: F(1)}, vector_field, order).items():
            matrix[index[output]][column] = coefficient
    return matrix


def nilpotent_matrix_exp(matrix, time=1):
    matrix, time = _square(matrix), exact(time)
    if not matrix:
        return []
    result, power = identity(len(matrix)), identity(len(matrix))
    for degree in range(1, len(matrix)+1):
        power = matrix_mul(power, matrix)
        if not any(any(row) for row in power):
            return result
        result = matrix_add(result, power, time**degree/factorial(degree))
    raise ValueError('The exact finite exponential helper requires a nilpotent matrix')


def quartic_response_controls(kappa=1):
    kappa = exact(kappa)
    q, momentum = coordinate_polynomials(2)
    field = (momentum, poly_scale(poly_pow(q, 3), -kappa))
    z = poly_pow(q, 2)
    responses = [z]
    for _ in range(6):
        responses.append(lie_derivative(responses[-1], field))
    energy = poly_add(poly_scale(poly_pow(momentum, 2), F(1, 2)),
                      poly_scale(poly_pow(q, 4), kappa/4))
    decoded = poly_scale(poly_add(responses[2], poly_scale(poly_mul(z, z), 3*kappa)), F(1, 4))
    derivative = lie_derivative(energy, field)
    crossing2 = {exponent: coefficient for exponent, coefficient in responses[2].items() if exponent[0] == 0}
    crossing6 = {exponent: coefficient for exponent, coefficient in responses[6].items() if exponent[0] == 0}
    assert decoded == energy and derivative == {}
    assert crossing2 == {(0, 2): F(2)}
    assert crossing6 == ({(0, 4): -72*kappa} if kappa else {})
    return dict(kappa=kappa, vector_field=field, z=z, z_first=responses[1],
                z_second=responses[2], z_sixth=responses[6], energy=energy,
                energy_from_response=decoded, energy_derivative=derivative,
                zero_crossing_z_second=crossing2, zero_crossing_z_sixth=crossing6,
                response_jets_are_time_derivatives=True,
                state_jet_spectral_values_and_time_response_entries_are_distinct=True,
                signed_state_recovered_from_even_observation=False)


def quartic_energy_from_response(z, z_second, kappa=1):
    z, z_second, kappa = exact(z), exact(z_second), exact(kappa)
    if z < 0:
        raise ValueError('The observed squared position must be nonnegative')
    return (z_second+3*kappa*z*z)/4


def quartic_parameter_from_response(z, z_first, z_second):
    z, z_first, z_second = exact(z), exact(z_first), exact(z_second)
    if z <= 0:
        raise ValueError('This parameter decoder requires strictly positive squared position')
    return (z_first*z_first-2*z*z_second)/(4*z**3)


def quartic_crossing_parameter(z_second, z_sixth):
    z_second, z_sixth = exact(z_second), exact(z_sixth)
    if z_second <= 0:
        raise ValueError('A nonzero crossing momentum is required')
    return -z_sixth/(18*z_second*z_second)


def quartic_energy_error_budget(kappa, radius, error_z, error_z_second):
    kappa = exact(kappa)
    radius, error_z, error_z_second = map(exact, (radius, error_z, error_z_second))
    if min(radius, error_z, error_z_second) < 0:
        raise ValueError('Radius and error budgets must be nonnegative')
    return error_z_second/4+3*abs(kappa)*radius*radius*error_z/2+3*abs(kappa)*error_z*error_z/4


def polynomial_records(poly):
    return [dict(exponent=list(exponent), coefficient=coefficient)
            for exponent, coefficient in sorted(_poly(poly).items(), key=lambda item: (sum(item[0]), item[0]))]


def _record_polynomials(row):
    result = {}
    for key, value in row.items():
        if key == 'vector_field':
            result[key] = [polynomial_records(poly) for poly in value]
        elif isinstance(value, dict):
            result[key] = polynomial_records(value)
        else:
            result[key] = value
    return result


def all_controls():
    x, y = coordinate_polynomials(2)
    family = (
        (poly_add(poly_scale(x, 2), y, poly_pow(x, 2)),
         poly_add(poly_scale(y, 3), poly_scale(poly_mul(x, y), -1))),
        (poly_add(x, poly_pow(y, 2)),
         poly_add(x, poly_scale(y, 2), poly_mul(poly_pow(x, 2), y))),
    )
    other_family = (
        (poly_add(poly_scale(x, 2), y, poly_scale(poly_pow(x, 2), -2), poly_scale(poly_mul(x, y), 3)),
         poly_add(poly_scale(y, 3), poly_pow(x, 3))),
        (poly_add(x, poly_scale(poly_mul(x, y), -2)),
         poly_add(x, poly_scale(y, 2), poly_pow(y, 3))),
    )
    terms = {(): F(2), (0,): F(3), (1, 0): F(-2), (0, 1, 1): F(1)}
    jet_rows = []
    for order in (2, 3):
        matrices = tuple(jet_pullback(images, order) for images in family)
        others = tuple(jet_pullback(images, order) for images in other_family)
        linear = tuple(jet_pullback(linear_map_polynomials(linear_part(images)), order) for images in family)
        for index, (matrix, other, graded) in enumerate(zip(matrices, others, linear)):
            assert preserves_degree_filtration(matrix, 2, order)
            assert graded_blocks(matrix, 2, order) == graded_blocks(graded, 2, order)
            assert matrix_charpoly(matrix) == matrix_charpoly(other) == matrix_charpoly(graded)
            expected = tuple({exponent: coefficient for exponent, coefficient in poly.items() if sum(exponent) <= order}
                             for poly in family[index])
            assert recover_coordinate_jet(matrix, 2, order) == expected
        word_rows = []
        for word in ((), (0,), (1,), (0, 1), (1, 0), (0, 1, 0, 1)):
            values = tuple(matrix_charpoly(word_matrix(bank, word)) for bank in (matrices, others, linear))
            assert values[0] == values[1] == values[2]
            word_rows.append(dict(word=list(word), characteristic_polynomial=values[0]))
        polynomial_values = tuple(matrix_charpoly(matrix_polynomial(bank, terms)) for bank in (matrices, others, linear))
        assert polynomial_values[0] == polynomial_values[1] == polynomial_values[2]
        jet_rows.append(dict(order=order, dimension=len(monomials(2, order)),
                             linear_parts=[linear_part(images) for images in family],
                             nonlinear_probe_count=nonlinear_probe_count(2, order),
                             recovered_coordinate_jets=[[polynomial_records(poly) for poly in recover_coordinate_jet(matrix, 2, order)]
                                                        for matrix in matrices],
                             word_spectra=word_rows,
                             noncommutative_polynomial_terms=[dict(word=list(word), coefficient=coefficient) for word, coefficient in terms.items()],
                             noncommutative_polynomial_characteristic=polynomial_values[0],
                             common_fixed_point_required=True))
    coordinate = coordinate_polynomials(1)[0]
    straight = jet_pullback((coordinate,), 2)
    curved = jet_pullback((poly_add(coordinate, poly_pow(coordinate, 2)),), 2)
    assert matrix_charpoly(straight) == matrix_charpoly(curved) == (1, -3, 3, -1)
    jordan = dict(order=2, characteristic_polynomial=matrix_charpoly(curved),
                  linear_map_rank_after_identity_subtraction=matrix_rank(matrix_add(straight, identity(3), -1)),
                  nonlinear_map_rank_after_identity_subtraction=matrix_rank(matrix_add(curved, identity(3), -1)),
                  nonlinear_square_after_identity_subtraction=matrix_power(matrix_add(curved, identity(3), -1), 2),
                  full_Jordan_information_claimed_blind=False)
    mechanics = []
    for kappa in (F(0), F(1), F(3, 2)):
        response = quartic_response_controls(kappa)
        generator = jet_lie_generator(response['vector_field'], 3)
        flow = nilpotent_matrix_exp(generator, F(1, 2))
        assert matrix_charpoly(generator) == (F(1),)+(F(0),)*10
        assert matrix_charpoly(flow) == tuple(F((-1)**k*comb(10, k)) for k in range(11))
        position, momentum = F(2, 3), F(5, 4)
        values = {key: poly_evaluate(response[key], (position, momentum))
                  for key in ('z', 'z_first', 'z_second', 'energy')}
        recovered_kappa = quartic_parameter_from_response(values['z'], values['z_first'], values['z_second'])
        crossing_kappa = quartic_crossing_parameter(poly_evaluate(response['z_second'], (0, momentum)),
                                                    poly_evaluate(response['z_sixth'], (0, momentum)))
        assert recovered_kappa == crossing_kappa == kappa
        assert quartic_energy_from_response(values['z'], values['z_second'], kappa) == values['energy']
        mechanics.append(dict(response=_record_polynomials(response), state_jet_order=3,
                              generator_characteristic_polynomial=matrix_charpoly(generator),
                              flow_time=F(1, 2), flow_characteristic_polynomial=matrix_charpoly(flow),
                              generator_power_ranks=[matrix_rank(matrix_power(generator, degree)) for degree in range(7)],
                              decoded_kappa_away_from_zero=recovered_kappa,
                              decoded_kappa_at_nonzero_momentum_crossing=crossing_kappa,
                              finite_response_identities_do_not_recover_global_signed_state=True))
    return dict(model=MODEL, jet_family_controls=jet_rows, Jordan_boundary=jordan,
                quartic_mechanics_controls=mechanics,
                basis_convention='Taylor-normalized monomials; matrix columns are images',
                dual_convention='algebraic transpose dual, not complex Hermitian adjoint')


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise ValueError('Serialize polynomial coefficient maps with polynomial_records first')
        return {key: encode(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(child) for child in value]
    if isinstance(value, float):
        raise ValueError('Inexact floating-point certificate value')
    return value


if __name__ == '__main__':
    runtime_check()
    print(json.dumps(encode(all_controls()), sort_keys=True, indent=2))
