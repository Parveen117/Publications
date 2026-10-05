import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym81_local_third_moment_vacuum_response as z
import ym80_bounded_phase_energy_certificates as stage80


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


def epsilon(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return 1 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def multiply(left, right):
    return [[sum((a*b for a, b in zip(row, column)), F(0))
             for column in zip(*right)] for row in left]


def subtract(left, right):
    return [[a-b for a, b in zip(row, other)] for row, other in zip(left, right)]


def determinant2(matrix):
    return matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]


def independent_bridges(side):
    """Tuple-labelled links and sets; independent of native integer edge IDs."""
    def shift(vertex, axis):
        moved = list(vertex)
        moved[axis] = (moved[axis]+1) % side
        return tuple(moved)
    def face(vertex, first, second):
        return frozenset(((vertex, first), (shift(vertex, first), second),
                          (shift(vertex, second), first), (vertex, second)))
    roots = tuple(product(range(0, side, 2), range(0, side, 2), range(side)))
    owners = {}
    for root in roots:
        for edge in face(root, 0, 1):
            if edge in owners:
                raise AssertionError('Selected link supports overlap')
            owners[edge] = root
    bridges, max_touched = {}, 0
    for vertex in product(range(side), repeat=3):
        for first, second in combinations(range(3), 2):
            edges = face(vertex, first, second)
            touched = sorted({owners[edge] for edge in edges if edge in owners})
            max_touched = max(max_touched, len(touched))
            if len(touched) == 2:
                self_counts = [sum(owners.get(edge) == block for edge in edges) for block in touched]
                if self_counts != [1, 1]:
                    raise AssertionError('A bridge has more than one link in one block')
                pair = tuple(touched)
                bridges[pair] = bridges.get(pair, 0)+1
    return roots, bridges, max_touched


class VacuumResponseMathematicsTests(unittest.TestCase):
    def test_native_Lie_connection_has_Ricci_one_half(self):
        # R(i,j)k = nabla_i nabla_j k - nabla_j nabla_i k - nabla_[i,j] k.
        ricci = [[sum((
            F(epsilon(j, k, m)*epsilon(i, m, i), 4)
            -F(epsilon(i, k, m)*epsilon(j, m, i), 4)
            -F(epsilon(i, j, m)*epsilon(m, k, i), 2)
            for i in range(3) for m in range(3)), F(0))
            for k in range(3)] for j in range(3)]
        self.assertEqual(ricci, [[F(int(i == j), 2) for j in range(3)] for i in range(3)])
        row = z.lie_ricci_control()
        self.assertEqual(row['Ricci'], ricci)
        self.assertEqual(row['fundamental_Casimir'], [[F(3*int(i == j), 4) for j in range(4)] for i in range(4)])
        for representation, action in zip(row['fundamental_representation_matrices'], row['native_coordinate_action_matrices']):
            self.assertEqual(action, [[-value for value in line] for line in representation])

    def test_connection_correction_is_half_the_score_energy(self):
        hessian = [[F(1), F(2), F(3)], [F(2), F(-4), F(5)], [F(3), F(5), F(6)]]
        score = (F(2), F(-3), F(5))
        ordered = [[hessian[i][j]+sum((F(epsilon(i, j, k), 2)*score[k]
                                      for k in range(3)), F(0)) for j in range(3)] for i in range(3)]
        defect = sum(value*value for row in ordered for value in row)-sum(value*value for row in hessian for value in row)
        self.assertEqual(defect, sum(value*value for value in score)/2)
        self.assertEqual(defect, 19)

    def test_native_Haar_derivatives_independently_fix_connection_sign_and_linear_response(self):
        native, poly, model = stage80.native, stage80.p, stage80.previous.model
        W = model.loop_polynomial(((0, 1), (1, 1), (2, -1), (3, -1)))
        first = {axis: native.derivative(W, 0, axis+1) for axis in range(3)}
        commutator = poly.add(native.derivative(first[1], 0, 1),
                              poly.scale(native.derivative(first[0], 0, 2), -1))
        self.assertEqual(commutator, first[2])
        electric = sum((native.phi(poly.mul(value, value), 4) for value in first.values()), F(0))
        ordered, covariant = F(0), F(0)
        for link in range(4):
            for left in range(3):
                for right in range(3):
                    derivative = native.derivative(first[right], link, left+1)
                    ordered += native.phi(poly.mul(derivative, derivative), 4)
                    corrected = derivative
                    if link == 0:
                        corrected = poly.add(corrected, *(poly.scale(first[k], -F(epsilon(left, right, k), 2))
                                                           for k in range(3) if epsilon(left, right, k)))
                    covariant += native.phi(poly.mul(corrected, corrected), 4)
        self.assertEqual(native.phi(poly.mul(W, W), 4), F(1, 4))
        self.assertEqual(electric, F(3, 16))
        self.assertEqual(ordered, F(9, 16))
        self.assertEqual(covariant, F(15, 32))
        self.assertEqual(ordered-covariant, electric/2)
        # h=1+(theta/3)V+O(theta²); four plaquettes meet a fixed link.
        m_linear = F(2, 3)*F(1, 4)
        tau_quadratic, R_quadratic = F(4, 9)*electric, F(4, 9)*covariant
        self.assertEqual((m_linear, tau_quadratic, R_quadratic), (F(1, 6), F(1, 12), F(5, 24)))
        self.assertEqual(R_quadratic+tau_quadratic/2, F(3, 2)*m_linear)

    def test_sparse_bridges_have_single_link_intersections_and_periodic_multiplicity(self):
        for side in (4, 6):
            roots, bridges, max_touched = independent_bridges(side)
            self.assertEqual(len(roots), side**3//4)
            self.assertEqual(max_touched, 2)
            self.assertEqual(sum(bridges.values()), 6*len(roots))
            for root in roots:
                neighbors = {pair[1] if pair[0] == root else pair[0]: count
                             for pair, count in bridges.items() if root in pair}
                expected = {}
                for axis, step in enumerate((2, 2, 1)):
                    for sign in (-1, 1):
                        moved = list(root)
                        moved[axis] = (moved[axis]+sign*step) % side
                        key = tuple(moved)
                        expected[key] = expected.get(key, 0)+(4 if axis == 2 else 1)
                self.assertEqual(neighbors, expected)
                self.assertEqual(sum(neighbors.values()), 12)
                self.assertEqual(len(neighbors), 4 if side == 4 else 6)
            replay = z.native_third_moment_geometry(side)
            self.assertEqual(replay['selected_plaquettes'], len(roots))
            self.assertEqual(replay['occupied_links'], 4*len(roots))
            self.assertEqual(replay['bridge_faces_by_axis'], dict(x=len(roots), y=len(roots), z=4*len(roots)))
            module_neighbors = {tuple((2*entry['cell'][0], 2*entry['cell'][1], entry['cell'][2])): entry['bridge_faces']
                                for entry in replay['origin_neighbor_bridges']}
            own_neighbors = {pair[1]: count for pair, count in bridges.items() if pair[0] == (0, 0, 0)}
            self.assertEqual(module_neighbors, own_neighbors)
            self.assertTrue(replay['every_bridge_shares_one_link_with_each_block'])
            self.assertEqual(replay['parallel_horizontal_bonds_counted'], side == 4)

    def test_mixed_bridge_normalization_is_sharp_for_positive_and_inverse_links(self):
        native, model = stage80.native, stage80.previous.model
        W = model.loop_polynomial(((0, 1), (1, 1), (2, -1), (3, -1)))
        flat = ((F(1), F(0), F(0), F(0)),)*4
        for other in (1, 2, 3):
            mixed = native.derivative(native.derivative(W, 0, 1), other, 1)
            # Each selected-block link gradient can have norm1/2, giving the extra1/4.
            contracted = model.evaluate(mixed, flat)/4
            self.assertEqual(abs(contracted), F(1, 16))
        geometry = z.native_third_moment_geometry(4)
        self.assertEqual(geometry['native_mixed_W_derivative_upper'], F(1, 16))
        self.assertEqual(geometry['coupling_over_theta_upper'], dict(x=F(1, 8), y=F(1, 8), z=F(1, 2)))

    def test_actual_Ward_row_budget_keeps_native_and_physical_scales_distinct(self):
        for theta in (F(25), F(256), F(65536)):
            row = z.ward_response_budget(theta, spacing=F(1, 8), clock_scale=3)
            self.assertEqual(row['universal_bare_Hessian_row'], dict(lower=F(3, 5)*theta, upper=F(3, 2)*theta))
            self.assertEqual(row['physical_Hessian_row'], dict(lower=F(3, 5)*9*64, upper=F(3, 2)*9*64))
            self.assertGreaterEqual(row['bare_Hessian_row']['lower'], F(3, 5)*theta)
            self.assertLessEqual(row['bare_Hessian_row']['lower'], row['bare_Hessian_row']['upper'])
            self.assertTrue(row['actual_gauge_invariant_vacuum'])
            self.assertFalse(row['score_vectors_are_gauge_invariant_states'])
            self.assertFalse(row['positive_local_gap_inferred'])
        self.assertEqual(z.ward_response_budget(25)['mean_plaquette_lower'], F(11, 20))

    def test_Ward_and_selected_geometry_reject_wrong_regimes(self):
        for args in ((24,), (25, 0), (25, 1, -1)):
            with self.assertRaises(ValueError): z.ward_response_budget(*args)
        for side in (2, 3, 5):
            with self.assertRaises((ValueError, TypeError)): z.native_third_moment_geometry(side)

    def test_noncommuting_moment_compression_has_a_strict_positive_defect(self):
        energies = (F(1), F(4), F(9), F(16))
        source = [[F(1, 2), F(7, 10)], [F(1, 4), F(-1, 20)],
                  [F(1, 6), F(1, 30)], [F(1, 8), F(-7, 40)]]
        moments = [multiply(transpose(source), [[energy**power*value for value in row]
                    for energy, row in zip(energies, source)]) for power in (1, 2, 3)]
        first, second, third = moments
        self.assertEqual(first, [[1, 0], [0, 1]])
        self.assertEqual(second, [[F(15, 2), F(-5)], [F(-5), F(423, 50)]])
        self.assertEqual(third, [[F(177, 2), F(-86)], [F(-86), F(1269, 10)]])
        defect = subtract(third, multiply(second, second))
        self.assertEqual(defect, [[F(29, 4), F(-31, 5)], [F(-31, 5), F(75821, 2500)]])
        self.assertEqual(determinant2(defect), F(1814409, 10000))
        self.assertGreater(defect[0][0], 0)
        self.assertGreater(determinant2(defect), 0)
        commutator = subtract(multiply(second, third), multiply(third, second))
        self.assertEqual(commutator[0][1], F(-2736, 25))
        # Equality, the reversed inequality, and treating these matrices as scalars all fail.
        self.assertNotEqual(defect, [[0, 0], [0, 0]])
        replay = z.matrix_moment_compression(energies, source)
        self.assertEqual(replay['M1'], first)
        self.assertEqual(replay['M2'], second)
        self.assertEqual(replay['M3'], third)
        self.assertEqual(replay['compression_deficit'], defect)
        self.assertEqual(replay['commutator'], commutator)

    def test_diagonal_finite_range_does_not_force_low_momentum_cancellation(self):
        for n in (2, 7, 32):
            energy, theta, source_scale = F(n*n), F(n**4), F(1, n)
            self.assertEqual(source_scale**2*energy, 1)
            self.assertEqual(source_scale**2*energy**3/theta, 1)
            self.assertEqual(source_scale**2*energy**2, n*n)
        self.assertFalse(z.MODEL.get('uniform_positive_continuum_survival', False))

    def test_circulant_identity_counts_parallel_bonds_and_keeps_zero_mean_intercept(self):
        for shape in ((2, 2, 4), (4, 4, 4)):
            cells = tuple(product(*(range(length) for length in shape)))
            coefficients = [F((-1)**cell[0]) for cell in cells]
            row = z.coefficient_gradient_identity(shape, coefficients, 3, (F(-1, 8), F(-1, 8), F(1, 2)))
            self.assertEqual(row['coefficient_sum'], 0)
            self.assertEqual(row['coefficient_squared_norm'], len(cells))
            self.assertEqual(row['gradient_squared_by_axis'], dict(x=4*len(cells), y=0, z=0))
            self.assertEqual(row['zero_momentum_coefficient'], F(7, 2))
            self.assertEqual(row['scalar_intercept_contribution'], F(7, 2)*len(cells))
            self.assertEqual(row['matrix_quadratic_form'], 4*len(cells))
            self.assertEqual(row['matrix_quadratic_form'], row['gradient_quadratic_form'])
            self.assertTrue(row['circulant_fixture_not_actual_Yang_Mills_vacuum'])

    def test_exact_rational_Fourier_mode_matches_circulant_quadratic_form(self):
        shape = (4, 2, 4)
        values = (F(1), F(0), F(-1), F(0))
        coefficients = [values[cell[0]] for cell in product(*(range(length) for length in shape))]
        couplings = (F(-1, 8), F(-1, 8), F(1, 2))
        row = z.coefficient_gradient_identity(shape, coefficients, 3, couplings)
        symbol = z.third_moment_symbol(3, couplings, (0, 1, 1))
        self.assertEqual(symbol['symbol'], F(15, 4))
        self.assertEqual(row['matrix_quadratic_form'], symbol['symbol']*row['coefficient_squared_norm'])
        self.assertNotEqual(row['scalar_intercept_contribution'], 0)
        with self.assertRaises(ValueError): z.third_moment_symbol(1, (1, 2, 3), (0, 0, 2))
        with self.assertRaises(ValueError): z.coefficient_gradient_identity((2, 2, 2), [1], 1, (0, 0, 0))

    def test_spatial_remainder_has_one_half_constant_and_a_conditional_intercept(self):
        rows = [z.physical_third_moment_budget(F(1, n), 2, 5, 7, clock_scale=3) for n in (4, 64)]
        for row in rows:
            self.assertEqual(row['exact_raw_physical_energy'], 6)
            self.assertEqual(row['spatial_gradient_remainder_upper'], F(135, 2))
            self.assertEqual(row['physical_third_moment_upper'], F(219, 2))
            self.assertFalse(row['actual_uniform_zero_mode_premise_certified'])
            self.assertFalse(row['zero_mean_coefficients_remove_scalar_intercept'])
        self.assertEqual(rows[0]['physical_third_moment_upper'], rows[1]['physical_third_moment_upper'])
        for args in ((0, 1, 1, 1), (1, -1, 1, 1), (1, 1, -1, 1), (1, 1, 1, -1)):
            with self.assertRaises(ValueError): z.physical_third_moment_budget(*args)

    def test_native_bridge_and_cell_normalizations_give_same_physical_half_in_all_axes(self):
        geometry = z.native_third_moment_geometry(4)
        for theta, root, spacing, clock in ((F(256), F(16), F(3, 10), F(5, 2)),
                                            (F(25), F(5), F(2, 7), F(3))):
            kappa = clock/(spacing*root)
            cell_volume = 4*spacing**3
            coefficient_squared_factor = 4*root*spacing**4
            self.assertEqual(kappa*coefficient_squared_factor/cell_volume, clock)
            axis_gradients = dict(x=F(2), y=F(3), z=F(5))
            physical_remainder = F(0)
            for axis, displacement in zip('xyz', (2*spacing, 2*spacing, spacing)):
                bridge_coupling = theta*geometry['coupling_over_theta_upper'][axis]
                # Sum |d_(b+step)-d_b|² <= 4sqrt(theta)a^4 * step²/(4a³) ||partial f||².
                coefficient_difference_per_gradient = coefficient_squared_factor*displacement**2/cell_volume
                physical_per_gradient = kappa**3*bridge_coupling*coefficient_difference_per_gradient
                self.assertEqual(physical_per_gradient, clock**3/2)
                # Fourier variation: 2|t|(1-cos(step*p)) <= |t|step²*p².
                self.assertEqual(kappa**2*bridge_coupling*displacement**2, clock**2/2)
                physical_remainder += physical_per_gradient*axis_gradients[axis]
            replay = z.physical_third_moment_budget(spacing, 2, sum(axis_gradients.values()), 7, clock_scale=clock)
            self.assertEqual(replay['spatial_gradient_remainder_upper'], physical_remainder)

    def test_a_fixed_zero_mode_premise_gives_only_conditional_uniform_phase_survival(self):
        rows = [z.conditional_phase_budget(F(1, n), 1, 2, 3) for n in (4, 64, 4096)]
        for row in rows:
            self.assertGreater(row['reflected_norm_squared_lower'], 0)
            self.assertTrue(row['conditional_on_actual_uniform_scaled_zero_mode_bound'])
            self.assertFalse(row['premise_disproved_or_proved_by_this_certificate'])
            self.assertFalse(row['actual_unconditional_continuum_survival'])
            self.assertEqual(row['conditional_phase_moment_scale'], 2)
        self.assertEqual(len({row['reflected_norm_squared_lower'] for row in rows}), 1)
        self.assertEqual(z.conditional_phase_budget(F(1, 8), 0, 2, 3)['reflected_norm_squared_lower'], 0)
        with self.assertRaises(ValueError): z.conditional_phase_budget(1, 1, 1, 1, beta=2)

    def test_rational_window_is_a_resolvent_difference_with_both_endpoint_budgets(self):
        for x in (F(1, 1000), F(1, 3), F(1), F(3), F(1000)):
            row = z.window_scalar(x)
            value = 2*(1/(1+x/3)-1/(1+3*x))
            self.assertEqual(row['filter'], value)
            self.assertEqual(row['resolvent_difference'], value)
            self.assertEqual(z.window_scalar(1/x)['filter'], value)
            self.assertLessEqual(x*value**2, 2)
            self.assertLessEqual(value**2/x, 2)
            b = x/(1+x)
            self.assertLessEqual(b**3/value, (3*x+b*b)/16)
            self.assertLessEqual(b**4/value, 3*x/16)
        self.assertEqual(z.window_scalar(0)['filter'], 0)
        self.assertEqual(z.window_scalar(1)['filter'], 1)
        with self.assertRaises(ValueError): z.window_scalar(-1)

    def test_rational_window_transfer_bounds_an_independent_three_atom_spectrum(self):
        energies, weights, time = (F(1, 4), F(1), F(4)), (F(1, 4), F(1, 3), F(1, 6)), F(3, 2)
        mass = sum(weights)
        energy = sum(weight*value for weight, value in zip(weights, energies))
        response = [time*value/(1+time*value) for value in energies]
        filters = [16*time*value/(3*(time*value)**2+10*time*value+3) for value in energies]
        R = sum(weight*b*b for weight, b in zip(weights, response))
        row = z.rational_window_transfer(time, R, energy, mass)
        actual_norm = sum(weight*g*g for weight, g in zip(weights, filters))
        self.assertLessEqual(row['filter_expectation_lower'], sum(weight*g for weight, g in zip(weights, filters)))
        self.assertLessEqual(row['filtered_norm_squared_lower'], actual_norm)
        self.assertLessEqual(row['basic_Cauchy_filtered_norm_squared_lower'], actual_norm)
        self.assertLessEqual(sum(weight*value*g*g for weight, value, g in zip(weights, energies, filters)), row['filtered_energy_upper'])
        self.assertLessEqual(sum(weight*g*g/value for weight, value, g in zip(weights, energies, filters)), row['filtered_inverse_energy_upper'])
        with self.assertRaises(ValueError): z.rational_window_transfer(1, F(1, 2), 1)
        with self.assertRaises(ValueError): z.rational_window_transfer(1, 2, 10)
        self.assertEqual(z.rational_window_transfer(1, 0, 0)['filtered_norm_squared_lower'], 0)

    def test_window_band_bound_refers_to_original_source_mass(self):
        energies, weights = (F(1, 128), F(1), F(128)), (F(1, 4), F(1, 2), F(1, 4))
        norm = sum(weight*z.window_scalar(energy)['filter']**2 for weight, energy in zip(weights, energies))
        row = z.window_band_budget(1, norm)
        self.assertEqual(row['spectral_window'], dict(lower=F(1, 32), upper=F(32)))
        self.assertEqual(row['filter_squared_outside_upper'], F(512, 3395)**2)
        self.assertGreater(row['original_source_spectral_mass_lower'], 0)
        self.assertLessEqual(row['original_source_spectral_mass_lower'], F(1, 2))
        self.assertTrue(row['bound_on_original_source_measure'])

    def test_actual_phase_window_has_a3_certificate_and_uniform_two_sided_energy_budgets(self):
        rows = [z.phase_window_budget(25, F(1, n), 1) for n in (4, 64, 4096)]
        for row in rows:
            self.assertGreater(row['filtered_norm_squared_lower'], 0)
            self.assertGreaterEqual(row['filtered_norm_squared_lower'], row['uncorrelated_transfer_lower'])
            self.assertLess(row['normalized_a3_lower_coefficient'], row['fixed_data_a3_asymptotic_coefficient_lower'])
            self.assertEqual(row['filtered_energy_upper'], 2)
            self.assertEqual(row['filtered_inverse_energy_upper'], 2)
            self.assertFalse(row['uniformly_positive_continuum_lower_bound'])
        self.assertTrue(all(first['normalized_a3_lower_coefficient'] < second['normalized_a3_lower_coefficient']
                            for first, second in zip(rows, rows[1:])))

    def test_actual_four_cube_window_reuses_actual_ground_state_anchor_in_bare_units(self):
        row = z.actual_cube_window()
        self.assertEqual(row['inherited_resolvent_lower'], F(32607, 250000))
        self.assertEqual(row['rational_window_transfer']['filtered_norm_squared_lower'], F(61632531115632, 1063322424390625))
        self.assertGreater(row['actual_original_source_band']['original_source_spectral_mass_lower'], F(1, 32))
        self.assertTrue(row['bare_energy_units'])
        self.assertTrue(row['hidden_space_untruncated'])
        self.assertFalse(row['uniform_continuum_survival_asserted'])


class VacuumResponseCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_manifest_and_parent_bind_the_expected_sources(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 302)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(len(self.cert['input_sha256']), 307)
        self.assertEqual(manifest['base_tree'], 'fa8e2cc81dbe343aaef0403012f36ef2ad3102cc')
        self.assertEqual(self.cert['parent_result_sha256'], '0f84160b93477194c7f18cec524f37d084b1cc7f46bbf67b1f58cc945a37f17f')
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_changed_upstream_parent_or_base_is_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'BASE_TREE', '0'*40):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_gap_zero_mode_and_charge_overclaims_are_rejected(self):
        for key in ('positive_continuum_mass_gap', 'uniform_positive_continuum_survival',
                    'actual_uniform_scaled_zero_mode_upper_bound_proved',
                    'vacuum_score_vectors_are_physical_gauge_invariant_states',
                    'zero_mean_smearing_removes_zero_mode_intercept'):
            wrong = copy.deepcopy(self.cert)
            wrong['model'][key] = True
            with tempfile.TemporaryDirectory() as directory:
                result, pin = Path(directory)/'result', Path(directory)/'pin'
                result.write_text(json.dumps(wrong))
                pin.write_text(z.p.canonical_sha(wrong))
                with self.assertRaises(ValueError): z.check(wrong, result, pin)
                with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_certificate_is_exact_and_JSON_roundtrippable(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                self.assertTrue(all(isinstance(key, str) for key in value))
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        self.assertEqual(json.loads(json.dumps(self.cert)), self.cert)

    def test_model_keeps_actual_identities_distinct_from_conditional_survival(self):
        model = self.cert['model']
        self.assertTrue(model['actual_vacuum_Hessian_row_Ward_identity'])
        self.assertTrue(model['actual_smooth_source_third_moment_has_nearest_neighbor_support'])
        self.assertTrue(model['actual_four_cube_rational_window_nonzero'])
        self.assertTrue(model['actual_smooth_phase_window_has_order_a3_lower_certificate'])
        self.assertFalse(model['actual_uniform_scaled_zero_mode_upper_bound_proved'])
        self.assertFalse(model['rational_filtering_alone_proves_uniform_continuum_survival'])


if __name__ == '__main__':
    unittest.main()
