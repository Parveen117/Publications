"""Independent exact controls for nonlinear pullback jets and flow responses.

Run with Python 3.12.  Determinants below use elimination independently of
the certificate's characteristic-polynomial routine.  Hand-expanded
polynomials and Lie-response formulas are intentionally kept explicit.
"""

from fractions import Fraction as F
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import sys
import unittest


CERTIFICATES = Path(__file__).resolve().parents[1] / "certificates"
sys.path.insert(0, str(CERTIFICATES))
import nsb1_jets as z


def canonical(poly):
    return {tuple(powers): F(value) for powers, value in poly.items() if value}


def matrix(rows):
    return tuple(tuple(F(value) for value in row) for row in rows)


def add_matrices(*terms):
    """Weighted matrix addition, separate from certificate matrix helpers."""
    size = len(terms[0][1])
    return tuple(
        tuple(sum(F(weight) * item[i][j] for weight, item in terms)
              for j in range(size))
        for i in range(size)
    )


def determinant(rows):
    """Exact pivoted elimination; used only as a characteristic oracle."""
    work = [list(map(F, row)) for row in rows]
    size = len(work)
    answer = F(1)
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if work[row][column]), None)
        if pivot is None:
            return F(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        diagonal = work[column][column]
        answer *= diagonal
        for row in range(column + 1, size):
            multiplier = work[row][column] / diagonal
            for entry in range(column + 1, size):
                work[row][entry] -= multiplier * work[column][entry]
            work[row][column] = F(0)
    return answer


def univariate_product(roots):
    """Descending coefficients of the hand-specified product of linear factors."""
    coefficients = [F(1)]
    for root in roots:
        following = [F(0)] * (len(coefficients) + 1)
        for index, value in enumerate(coefficients):
            following[index] += value
            following[index + 1] -= F(root) * value
        coefficients = following
    return tuple(coefficients)


def evaluate_univariate(coefficients, point):
    value = F(0)
    for coefficient in coefficients:
        value = value * point + coefficient
    return value


def independent_lie(poly, kappa):
    """Apply p*d/dq - kappa*q^3*d/dp by explicit monomial bookkeeping."""
    answer = {}
    for (q_degree, p_degree), coefficient in poly.items():
        if q_degree:
            powers = (q_degree - 1, p_degree + 1)
            answer[powers] = answer.get(powers, F(0)) + coefficient*q_degree
        if p_degree:
            powers = (q_degree + 3, p_degree - 1)
            answer[powers] = answer.get(powers, F(0)) - kappa*coefficient*p_degree
    return canonical(answer)


class PolynomialArithmeticTests(unittest.TestCase):
    def test_monomial_set_is_complete_without_assuming_its_order(self):
        for dimension, order in ((1, 0), (1, 5), (2, 3), (3, 3)):
            result = tuple(z.monomials(dimension, order))
            expected = {powers for powers in product(range(order + 1),
                                                     repeat=dimension)
                        if sum(powers) <= order}
            self.assertEqual(set(result), expected)
            self.assertEqual(len(result), comb(dimension + order, order))

    def test_coordinates_and_exact_arithmetic(self):
        x, y = z.coordinate_polynomials(2)
        self.assertEqual(canonical(x), {(1, 0): F(1)})
        self.assertEqual(canonical(y), {(0, 1): F(1)})
        quadratic = z.poly_pow(z.poly_add(x, z.poly_scale(y, F(2, 3))), 2)
        self.assertEqual(canonical(quadratic),
                         {(2, 0): F(1), (1, 1): F(4, 3), (0, 2): F(4, 9)})
        self.assertEqual(canonical(z.poly_add(quadratic,
                                             z.poly_scale(quadratic, -1))), {})
        self.assertEqual(canonical(z.poly_pow(x, 0)), {(0, 0): F(1)})

    def test_composition_matches_a_hand_expansion_and_truncation(self):
        polynomial = {(2, 0): 1, (1, 1): 2, (0, 1): -3, (0, 0): 5}
        images = ({(1, 0): 1, (0, 2): 1},
                  {(0, 1): 1, (2, 0): -1})
        expected = {(0, 0): 5, (0, 1): -3, (2, 0): 4, (1, 1): 2,
                    (3, 0): -2, (1, 2): 2, (0, 3): 2,
                    (0, 4): 1, (2, 2): -2}
        self.assertEqual(canonical(z.poly_compose(polynomial, images)),
                         canonical(expected))
        for order in (0, 1, 2, 3, 4):
            self.assertEqual(canonical(z.poly_compose(polynomial, images,
                                                      order=order)),
                             canonical({powers: value for powers, value in
                                        expected.items() if sum(powers) <= order}))

    def test_derivatives_and_substitution_at_rational_points(self):
        polynomial = {(3, 1): F(2, 5), (0, 2): -3, (0, 0): 7}
        self.assertEqual(canonical(z.poly_derivative(polynomial, 0)),
                         {(2, 1): F(6, 5)})
        self.assertEqual(canonical(z.poly_derivative(polynomial, 1)),
                         {(3, 0): F(2, 5), (0, 1): F(-6)})
        for q, p in ((F(2, 3), F(-5, 7)), (F(0), F(2)), (F(-2), F(1))):
            expected = F(2, 5)*q**3*p - 3*p**2 + 7
            self.assertEqual(z.poly_evaluate(polynomial, (q, p)), expected)

    def test_truncated_multiplication_really_discards_high_degree(self):
        p = {(1, 0): 2, (0, 2): 3}
        q = {(0, 1): 5, (2, 0): 7}
        self.assertEqual(canonical(z.poly_mul(p, q, order=2)),
                         {(1, 1): F(10)})
        self.assertEqual(canonical(z.poly_mul(p, q, order=3)),
                         {(1, 1): F(10), (3, 0): F(14), (0, 3): F(15)})

    def test_inexact_and_boolean_scalars_are_refused(self):
        x = {(1,): 1}
        for invalid in (0.5, True, False):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    z.poly_scale(x, invalid)
                with self.assertRaises(ValueError):
                    z.poly_evaluate(x, (invalid,))
                with self.assertRaises(ValueError):
                    z.jet_pullback(({(1,): invalid},), 2)


class PullbackAndSpectrumTests(unittest.TestCase):
    def assert_characteristic_by_determinants(self, item):
        coefficients = tuple(z.matrix_charpoly(item))
        size = len(item)
        self.assertEqual(len(coefficients), size + 1)
        self.assertEqual(coefficients[0], 1)
        # A degree-size polynomial is fixed by size+1 exact evaluations.
        for point in range(-2, size - 1):
            shifted = tuple(tuple(F(point)*(i == j) - item[i][j]
                                  for j in range(size)) for i in range(size))
            self.assertEqual(evaluate_univariate(coefficients, F(point)),
                             determinant(shifted))

    def test_characteristic_polynomial_against_independent_determinants(self):
        fixtures = (
            ((2,),),
            ((0, 1), (-2, 3)),
            ((2, 1, -1), (0, -3, 2), (4, 0, 1)),
            ((F(1, 2), 2, 0, 1), (-1, F(2, 3), 1, 0),
             (3, 0, -2, 1), (0, 1, 4, F(-1, 5))),
        )
        for fixture in fixtures:
            with self.subTest(size=len(fixture)):
                self.assert_characteristic_by_determinants(matrix(fixture))

    def test_columns_are_images_and_composition_reverses_map_order(self):
        f = ({(1, 0): 1, (0, 2): 1}, {(0, 1): 1})
        g = ({(1, 0): 1}, {(0, 1): 1, (1, 0): 1})
        basis = tuple(z.monomials(2, 3))
        index = {powers: i for i, powers in enumerate(basis)}
        cf = z.jet_pullback(f, 3)
        x_column = tuple(cf[row][index[(1, 0)]] for row in range(len(basis)))
        self.assertEqual(x_column,
                         tuple(F(powers in ((1, 0), (0, 2))) for powers in basis))
        # g after f: (x+y^2, x+y+y^2), distinct from f after g.
        g_after_f = ({(1, 0): 1, (0, 2): 1},
                     {(1, 0): 1, (0, 1): 1, (0, 2): 1})
        f_after_g = ({(1, 0): 1, (2, 0): 1, (1, 1): 2, (0, 2): 1},
                     {(1, 0): 1, (0, 1): 1})
        cg = z.jet_pullback(g, 3)
        self.assertEqual(matrix(z.matrix_mul(cf, cg)),
                         matrix(z.jet_pullback(g_after_f, 3)))
        self.assertEqual(matrix(z.matrix_mul(cg, cf)),
                         matrix(z.jet_pullback(f_after_g, 3)))
        self.assertNotEqual(matrix(z.matrix_mul(cf, cg)),
                            matrix(z.matrix_mul(cg, cf)))

    def test_filtration_and_graded_blocks_use_only_the_linear_part(self):
        nonlinear = ({(1, 0): 2, (0, 1): 1, (2, 0): 3, (1, 1): -2},
                     {(0, 1): -3, (0, 2): 4, (3, 0): 5})
        linear = ({(1, 0): 2, (0, 1): 1}, {(0, 1): -3})
        basis = tuple(z.monomials(2, 3))
        full = z.jet_pullback(nonlinear, 3)
        reference = z.jet_pullback(linear, 3)
        self.assertNotEqual(matrix(full), matrix(reference))
        for row, output_power in enumerate(basis):
            for column, input_power in enumerate(basis):
                if sum(output_power) < sum(input_power):
                    self.assertEqual(full[row][column], 0)
                if sum(output_power) == sum(input_power):
                    self.assertEqual(full[row][column], reference[row][column])
        eigenvalues = [F(2)**a * F(-3)**b for a, b in basis]
        self.assertEqual(tuple(z.matrix_charpoly(full)),
                         univariate_product(eigenvalues))
        self.assert_characteristic_by_determinants(full)

    def test_singular_linear_part_is_permitted(self):
        images = ({(2, 0): 1, (1, 1): -2}, {(0, 1): 2, (2, 0): 3})
        basis = tuple(z.monomials(2, 3))
        roots = [F(0) if a else F(2)**b for a, b in basis]
        self.assertEqual(tuple(z.matrix_charpoly(z.jet_pullback(images, 3))),
                         univariate_product(roots))

    def test_noncommutative_words_and_polynomial_preserve_blindness(self):
        f0 = ({(1, 0): 2, (0, 1): 1}, {(0, 1): 3})
        g0 = ({(1, 0): 1}, {(1, 0): 1, (0, 1): 2})
        f = (dict(f0[0]), dict(f0[1]))
        f[0].update({(2, 0): 1, (1, 1): -2})
        f[1].update({(0, 2): 3})
        g = (dict(g0[0]), dict(g0[1]))
        g[0].update({(0, 2): -1})
        g[1].update({(2, 0): 2, (1, 1): 1})
        pairs = [(z.jet_pullback(a, 3), z.jet_pullback(b, 3))
                 for a, b in ((f0, g0), (f, g))]
        outputs = []
        for first, second in pairs:
            ab = z.matrix_mul(first, second)
            ba = z.matrix_mul(second, first)
            aabb = z.matrix_mul(z.matrix_power(first, 2),
                               z.matrix_power(second, 2))
            abab = z.matrix_power(ab, 2)
            aab = z.matrix_mul(z.matrix_power(first, 2), second)
            polynomial = add_matrices((3, z.identity(len(first))), (2, ab),
                                      (-1, ba), (1, aab))
            self.assertNotEqual(matrix(ab), matrix(ba))
            outputs.append((ab, aabb, abab, polynomial))
        for linear_result, nonlinear_result in zip(*outputs):
            self.assertEqual(tuple(z.matrix_charpoly(linear_result)),
                             tuple(z.matrix_charpoly(nonlinear_result)))
        self.assertNotEqual(tuple(z.matrix_charpoly(outputs[0][1])),
                            tuple(z.matrix_charpoly(outputs[0][2])))
        self.assert_characteristic_by_determinants(outputs[1][-1])

    def test_public_word_and_polynomial_helpers_keep_the_declared_order(self):
        first = matrix(((2, 1), (0, 3)))
        second = matrix(((1, 0), (1, 2)))
        self.assertEqual(matrix(z.word_matrix((first, second), (0, 1))),
                         matrix(((3, 2), (3, 6))))
        self.assertEqual(matrix(z.word_matrix((first, second), (1, 0))),
                         matrix(((2, 1), (2, 7))))
        self.assertEqual(matrix(z.word_matrix((first, second), ())),
                         matrix(((1, 0), (0, 1))))
        # 3I + 2AB - BA has this direct two-by-two expansion.
        expression = {(): 3, (0, 1): 2, (1, 0): -1}
        self.assertEqual(matrix(z.matrix_polynomial((first, second), expression)),
                         matrix(((7, 3), (4, 8))))
        with self.assertRaises(ValueError):
            z.word_matrix((first, second), (2,))

    def test_coordinate_coefficient_recovery_and_probe_dimension(self):
        images = ({(1, 0): 2, (0, 1): -1, (2, 0): F(3, 7),
                   (0, 3): -4, (4, 0): 99},
                  {(1, 0): 1, (0, 1): 3, (1, 1): -2, (2, 1): F(5, 9)})
        recovered = z.recover_coordinate_jet(z.jet_pullback(images, 3), 2, 3)
        expected = tuple(canonical({power: value for power, value in item.items()
                                    if sum(power) <= 3}) for item in images)
        self.assertEqual(tuple(map(canonical, recovered)), expected)
        for dimension, order in ((1, 1), (1, 4), (2, 3), (3, 2)):
            self.assertEqual(z.nonlinear_probe_count(dimension, order),
                             dimension*(comb(dimension + order, order)-dimension-1))

    def test_each_nonlinear_coefficient_is_independently_spectrally_blind(self):
        dimension, order = 2, 2
        basis = tuple(z.monomials(dimension, order))
        nonlinear_powers = tuple(power for power in basis if sum(power) >= 2)
        labels = tuple((axis, power) for axis in range(dimension)
                       for power in nonlinear_powers)
        baseline = ({(1, 0): 2}, {(0, 1): 3})
        baseline_characteristic = z.matrix_charpoly(z.jet_pullback(baseline, order))
        response_columns = []
        for axis, power in labels:
            changed = [dict(item) for item in baseline]
            changed[axis][power] = F(1)
            operator = z.jet_pullback(tuple(changed), order)
            self.assertEqual(tuple(z.matrix_charpoly(operator)),
                             tuple(baseline_characteristic))
            recovered = z.recover_coordinate_jet(operator, dimension, order)
            response_columns.append(tuple(recovered[a].get(p, F(0))
                                          for a, p in labels))
        response = tuple(zip(*response_columns))
        self.assertEqual(matrix(response), matrix(z.identity(len(labels))))
        self.assertEqual(z.matrix_rank(response), z.nonlinear_probe_count(2, 2))
        self.assertEqual(z.matrix_rank(response[:-1]), len(labels)-1)

    def test_fixed_point_and_equilibrium_conditions_are_required(self):
        with self.assertRaises(ValueError):
            z.jet_pullback(({(0,): 1, (1,): 1},), 3)
        with self.assertRaises(ValueError):
            z.jet_lie_generator(({(0,): F(1, 3), (1,): -1},), 3)
        for invalid in (True, 1.5, -1):
            with self.subTest(order=invalid):
                with self.assertRaises(ValueError):
                    z.jet_pullback(({(1,): 1},), invalid)

    def test_jordan_structure_is_not_part_of_spectral_blindness(self):
        plain = z.jet_pullback(({(1,): 1},), 2)
        nonlinear = z.jet_pullback(({(1,): 1, (2,): F(7, 3)},), 2)
        expected = (F(1), F(-3), F(3), F(-1))
        self.assertEqual(tuple(z.matrix_charpoly(plain)), expected)
        self.assertEqual(tuple(z.matrix_charpoly(nonlinear)), expected)
        nilpotent = add_matrices((1, nonlinear), (-1, z.identity(3)))
        self.assertEqual(z.matrix_rank(nilpotent), 1)
        self.assertEqual(z.matrix_rank(z.matrix_power(nilpotent, 2)), 0)
        self.assertEqual(matrix(plain), matrix(z.identity(3)))


class QuarticFlowTests(unittest.TestCase):
    def test_response_polynomials_include_the_exact_sixth_derivative(self):
        for kappa in (F(0), F(1), F(3, 5), F(7)):
            result = z.quartic_response_controls(kappa)
            expected = {
                "z": {(2, 0): 1},
                "z_first": {(1, 1): 2},
                "z_second": {(0, 2): 2, (4, 0): -2*kappa},
                "z_sixth": {(8, 0): -144*kappa**3,
                             (4, 2): 936*kappa**2, (0, 4): -72*kappa},
                "energy": {(0, 2): F(1, 2), (4, 0): kappa/4},
            }
            for key, polynomial in expected.items():
                with self.subTest(kappa=kappa, key=key):
                    self.assertEqual(canonical(result[key]), canonical(polynomial))
            response = {(2, 0): F(1)}
            for _ in range(6):
                response = independent_lie(response, kappa)
            self.assertEqual(canonical(result["z_sixth"]), response)
            self.assertEqual(tuple(map(canonical, result["vector_field"])),
                             ({(0, 1): F(1)}, canonical({(3, 0): -kappa})))

    def test_energy_decoder_is_an_identity_and_energy_is_conserved(self):
        for kappa in (F(0), F(1), F(5, 2)):
            result = z.quartic_response_controls(kappa)
            self.assertEqual(canonical(result["energy_from_response"]),
                             canonical(result["energy"]))
            self.assertEqual(canonical(result["energy_derivative"]), {})
            self.assertEqual(independent_lie(result["energy"], kappa), {})
            for q, p in ((F(-2), F(3)), (F(0), F(-4, 3)), (F(1, 2), F(0))):
                observed = q*q
                second = 2*p*p-2*kappa*q**4
                self.assertEqual(z.quartic_energy_from_response(observed, second,
                                                                kappa=kappa),
                                 p*p/2+kappa*q**4/4)

    def test_energy_error_budget_is_attained_at_an_explicit_boundary(self):
        kappa, radius = F(5, 3), F(7, 4)
        error_z, error_second = F(1, 20), F(3, 40)
        squared_position, momentum = radius**2, F(2, 3)
        second = 2*momentum**2-2*kappa*radius**4
        exact_energy = momentum**2/2+kappa*radius**4/4
        observed_energy = z.quartic_energy_from_response(
            squared_position+error_z, second+error_second, kappa)
        expected_error = (error_second/4+3*kappa*radius**2*error_z/2
                          +3*kappa*error_z**2/4)
        budget = z.quartic_energy_error_budget(kappa, radius,
                                               error_z, error_second)
        self.assertEqual(budget, expected_error)
        self.assertEqual(observed_energy-exact_energy, budget)
        for q in (F(0), radius/2, radius):
            for sign_z, sign_second in product((-1, 1), repeat=2):
                # Test the polynomial error directly, including noisy negative z.
                difference = (sign_second*error_second/4
                              +3*kappa*((q*q+sign_z*error_z)**2-q**4)/4)
                self.assertLessEqual(abs(difference), budget)
        for parameters in ((-1, 0, 0), (1, -1, 0), (1, 0, -1)):
            with self.assertRaises(ValueError):
                z.quartic_energy_error_budget(kappa, *parameters)

    def test_parameter_recovery_in_configuration_and_crossing_charts(self):
        for kappa in (F(0), F(1, 7), F(3), F(11, 2)):
            for q, p in ((F(2), F(-3)), (F(-1, 3), F(2, 5)), (F(3), F(0))):
                observed, first = q*q, 2*q*p
                second = 2*p*p-2*kappa*q**4
                self.assertEqual(z.quartic_parameter_from_response(observed,
                                                                  first, second),
                                 kappa)
            for p in (F(-2), F(1, 3)):
                self.assertEqual(z.quartic_crossing_parameter(2*p*p,
                                                              -72*kappa*p**4),
                                 kappa)
            controls = z.quartic_response_controls(kappa)
            self.assertEqual(canonical(controls["zero_crossing_z_second"]),
                             {(0, 2): F(2)})
            self.assertEqual(canonical(controls["zero_crossing_z_sixth"]),
                             canonical({(0, 4): -72*kappa}))

    def test_invalid_parameter_chart_is_not_silently_divided_by_zero(self):
        for observed in (0, -1):
            with self.assertRaises(ValueError):
                z.quartic_parameter_from_response(observed, 0, 0)
        for second in (0, -1):
            with self.assertRaises(ValueError):
                z.quartic_crossing_parameter(second, 0)
        for inexact in (0.5, True):
            with self.assertRaises(ValueError):
                z.quartic_energy_from_response(1, 2, kappa=inexact)

    def test_time_jet_rank_repairs_crossings_but_not_the_equilibrium(self):
        controls = z.quartic_response_controls(F(2))
        observables = [controls[key] for key in ("z", "z_first", "z_second")]
        for p in (F(1), F(-2, 3)):
            jacobian = tuple(tuple(z.poly_evaluate(z.poly_derivative(item, axis),
                                                  (0, p)) for axis in range(2))
                             for item in observables)
            self.assertEqual(matrix(jacobian),
                             ((F(0), F(0)), (2*p, F(0)), (F(0), 4*p)))
            self.assertEqual(z.matrix_rank(jacobian[:1]), 0)
            self.assertEqual(z.matrix_rank(jacobian[:2]), 1)
            self.assertEqual(z.matrix_rank(jacobian), 2)
        at_equilibrium = tuple(tuple(z.poly_evaluate(z.poly_derivative(item, axis),
                                                     (0, 0)) for axis in range(2))
                                for item in observables)
        self.assertEqual(z.matrix_rank(at_equilibrium), 0)

    def test_central_symmetry_is_target_irrelevant_not_state_recovery(self):
        controls = z.quartic_response_controls(F(3, 2))
        for key in ("z", "z_first", "z_second", "z_sixth", "energy"):
            polynomial = controls[key]
            self.assertTrue(all(sum(power) % 2 == 0 for power in polynomial))
            self.assertEqual(z.poly_evaluate(polynomial, (F(2, 3), F(-5, 7))),
                             z.poly_evaluate(polynomial, (F(-2, 3), F(5, 7))))

    def test_lie_generator_spectrum_and_jordan_boundary_are_distinct(self):
        patterns = {0: (6, 3, 1, 0, 0, 0), 1: (7, 5, 3, 2, 1, 0)}
        for kappa, expected_ranks in patterns.items():
            field = ({(0, 1): 1}, {(3, 0): -kappa})
            generator = z.jet_lie_generator(field, 3)
            self.assertEqual(len(generator), 10)
            self.assertEqual(tuple(z.matrix_charpoly(generator)),
                             (F(1),)+(F(0),)*10)
            self.assertEqual(tuple(z.matrix_rank(z.matrix_power(generator, power))
                                   for power in range(1, 7)), expected_ranks)
            # A finite exact exponential suffices because the generator is nilpotent.
            time_one = add_matrices(*[(F(1, factorial(power)),
                                       z.matrix_power(generator, power))
                                      for power in range(6)])
            self.assertEqual(tuple(z.matrix_charpoly(time_one)),
                             univariate_product([F(1)]*10))

    def test_tangent_kernel_inequality_alone_does_not_give_a_decoder(self):
        # A(x)=x², Pi(x)=x³: differential domination holds, but fibers fail.
        for x in (F(-1), F(-1, 2), F(0), F(1, 2), F(1)):
            self.assertLessEqual(abs(3*x*x), F(3, 2)*abs(2*x))
        self.assertEqual(F(1, 2)**2, F(-1, 2)**2)
        self.assertNotEqual(F(1, 2)**3, F(-1, 2)**3)


class CertificateScopeTests(unittest.TestCase):
    def test_aggregate_records_exact_data_and_the_correct_scope(self):
        controls = z.all_controls()
        self.assertEqual(set(controls),
                         {"model", "jet_family_controls", "Jordan_boundary",
                          "quartic_mechanics_controls", "basis_convention",
                          "dual_convention"})
        self.assertTrue(controls["model"]["pullback_requires_a_fixed_base_point"])
        for key in ("all_similarity_invariants_are_blind",
                    "smooth_germs_determined_by_all_Taylor_coefficients",
                    "infinite_Koopman_spectrum_identified",
                    "physical_Yang_Mills_mass_gap_proved",
                    "arbitrary_set_theoretic_encodings_obey_the_probe_count"):
            self.assertIs(controls["model"][key], False)
        self.assertEqual([row["order"] for row in controls["jet_family_controls"]],
                         [2, 3])
        for row in controls["quartic_mechanics_controls"]:
            self.assertEqual(row["decoded_kappa_away_from_zero"],
                             row["decoded_kappa_at_nonzero_momentum_crossing"])
            self.assertFalse(row["response"]["signed_state_recovered_from_even_observation"])
        encoded = z.encode(controls)
        self.assertEqual(json.loads(json.dumps(encoded, sort_keys=True)), encoded)
        with self.assertRaises(ValueError):
            z.encode({"inexact": 0.1})
        with self.assertRaises(ValueError):
            z.encode({(1, 0): F(1)})


if __name__ == "__main__":
    unittest.main()
