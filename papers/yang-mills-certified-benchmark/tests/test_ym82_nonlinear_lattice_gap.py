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
import ym82_nonlinear_lattice_gap as z
import ym80_bounded_phase_energy_certificates as stage80


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


def independent_quadratic_extrema(theta, lam):
    points = [F(-1), F(1)]
    if lam and -1 <= -theta/(2*lam) <= 1:
        points.append(-theta/(2*lam))
    values = [theta*x+lam*x*x for x in points]
    return min(values), max(values)


def independent_partition(side):
    def shift(vertex, axis):
        out = list(vertex)
        out[axis] = (out[axis]+1) % side
        return tuple(out)
    def face(vertex, first, second):
        return frozenset(((vertex, first), (shift(vertex, first), second),
                          (shift(vertex, second), first), (vertex, second)))
    vertices = tuple(product(range(side), repeat=3))
    owners, selected_faces = {}, set()
    for root in product(range(0, side, 2), range(0, side, 2), range(side)):
        support = face(root, 0, 1)
        selected_faces.add(support)
        for edge in support:
            if edge in owners:
                raise AssertionError('Selected plaquettes overlap in links')
            owners[edge] = ('plaquette', root)
    for vertex in vertices:
        for axis in range(3):
            edge = (vertex, axis)
            owners.setdefault(edge, ('link', edge))
    degrees = {block: 0 for block in owners.values()}
    pair_multiplicities, cross_faces = {}, 0
    for vertex in vertices:
        for first, second in combinations(range(3), 2):
            support = face(vertex, first, second)
            if support in selected_faces:
                continue
            blocks = {owners[edge] for edge in support}
            if len(blocks) != 4:
                raise AssertionError('A weak plaquette meets fewer than four blocks')
            cross_faces += 1
            for block in blocks:
                degrees[block] += 1
            for block in blocks:
                for other in blocks-{block}:
                    pair_multiplicities[block, other] = pair_multiplicities.get((block, other), 0)+1
    return owners, degrees, pair_multiplicities, cross_faces


class NonlinearGapMathematicsTests(unittest.TestCase):
    def test_native_radial_chain_rule_verifies_quartic_ground_state_without_linearization(self):
        # For a loop class function, all four link Casimirs act identically.
        native, p, model = stage80.native, stage80.p, stage80.previous.model
        w = p.COORD[0]
        radius = p.add(*(p.mul(value, value) for value in p.COORD))
        constant = {(0, 0, 0, 0): F(1)}
        w2, w3, w4 = p.mul(w, w), p.mul(p.mul(w, w), w), p.mul(p.mul(w, w), p.mul(w, w))
        LE = lambda poly: p.scale(model.electric(poly, 1), 4)
        Gamma = lambda poly: p.scale(model.energy(poly, 1), 4)
        self.assertEqual(LE(w), p.scale(w, 3))
        self.assertEqual(LE(w2), p.add(p.scale(w2, 8), p.scale(radius, -2)))
        self.assertEqual(LE(w3), p.add(p.scale(w3, 15), p.scale(p.mul(radius, w), -6)))
        for t in (F(1, 3), F(2), F(7)):
            S = p.scale(w2, t)
            parent_potential = p.add(p.scale(LE(S), -F(1, 2)), p.scale(Gamma(S), F(1, 4)))
            nonlinear = p.add(p.scale(w4, t*t), p.scale(w2, 4*t-t*t))
            # Away from the unit sphere, the exact polynomial residual records radius-1.
            expected = p.mul(p.add(radius, p.scale(constant, -1)),
                             p.add(p.scale(constant, t), p.scale(w2, t*t)))
            residual = p.add(parent_potential, nonlinear, p.scale(constant, -t))
            self.assertEqual(residual, expected)
        self.assertEqual(native.phi(w2, 1), F(1, 4))
        self.assertEqual(native.phi(w4, 1), F(1, 8))
        replay = z.native_nonlinear_controls()
        self.assertEqual(replay['linear_mode_rate'], 3)
        self.assertEqual(replay['centered_quadratic_mode_rate'], 8)
        self.assertTrue(replay['native_holonomy_reduction_uses_four_link_Casimirs'])
        for row in replay['exact_doublewell_parent_identities']:
            self.assertEqual(row['ground_energy'], -row['t'])
            self.assertTrue(row['identity_holds_modulo_native_unit_radius'])

    def test_parent_gradient_square_retains_cross_plaquette_terms(self):
        dynamics = stage80.stage76
        model = stage80.previous.model
        configurations = ((F(3, 5), F(4, 5), 0, 0),)+(model.y.ONE,)*6
        words = (((0, 1), (1, 1), (2, -1), (3, -1)),
                 ((0, 1), (4, 1), (5, -1), (6, -1)))
        gradients = [dynamics.face_gradient(configurations, word)[1] for word in words]
        diagonal = sum(value*value for gradient in gradients for vector in gradient.values() for value in vector)
        cross = sum(a*b for a, b in zip(gradients[0][0], gradients[1][0]))
        self.assertEqual(diagonal, F(32, 25))
        self.assertEqual(cross, F(4, 25))
        self.assertEqual(diagonal+2*cross, F(8, 5))
        self.assertNotEqual(diagonal+2*cross, diagonal)

    def test_native_Haar_quadratic_response_separates_modes_and_does_not_identify_a_trial_vacuum(self):
        native, p, model = stage80.native, stage80.p, stage80.previous.model
        w = p.COORD[0]
        constant = {(0, 0, 0, 0): F(1)}
        w2 = p.mul(w, w)
        quadratic = p.add(w2, p.scale(constant, -F(1, 4)))
        self.assertEqual(native.phi(quadratic, 1), 0)
        self.assertEqual(native.phi(p.mul(w, quadratic), 1), 0)
        self.assertEqual(native.phi(p.mul(quadratic, quadratic), 1), F(1, 16))
        linear_inverse = p.scale(w, F(1, 3))
        quadratic_inverse = p.scale(quadratic, F(1, 8))
        self.assertEqual(native.phi(p.mul(w, linear_inverse), 1), F(1, 12))
        self.assertEqual(native.phi(p.mul(quadratic, quadratic_inverse), 1), F(1, 128))
        self.assertEqual(F(3, 2)*F(1, 6)-F(1, 2)*F(1, 12), F(5, 24))
        self.assertEqual(4*F(1, 64)-F(1, 2)*F(1, 128), F(15, 256))
        for theta, lam in ((F(1), F(0)), (F(0), F(1)), (F(-2, 3), F(5, 7))):
            row = z.nonlinear_ward_response(theta, lam)
            tau = theta**2/12+lam**2/128
            self.assertEqual(row['mean_potential_Haar'], lam/4)
            self.assertEqual(row['first_linear_response'], theta/6)
            self.assertEqual(row['first_centered_quadratic_response'], lam/64)
            self.assertEqual(row['ground_energy_per_face_through_order_two'], -lam/4-tau)
            self.assertEqual(row['electric_energy_per_link_quadratic'], tau)
            self.assertEqual(row['log_vacuum_Hessian_row_quadratic'], 5*theta**2/24+15*lam**2/256)
            self.assertEqual(row['mixed_theta_lambda_second_order_coefficient'], 0)
            self.assertTrue(row['response_coefficients_are_fixed_volume_perturbation'])
            self.assertFalse(row['perturbative_remainder_uniform_in_volume'])
            self.assertFalse(row['score_derivatives_are_physical_gauge_invariant_states'])

    def test_partition_counts_full_nonselected_faces_and_neighbor_multiplicity(self):
        for side in (4, 6):
            owners, degrees, pairs, face_count = independent_partition(side)
            selected = [block for block in degrees if block[0] == 'plaquette']
            singletons = [block for block in degrees if block[0] == 'link']
            self.assertEqual(len(owners), 3*side**3)
            self.assertEqual(len(selected), side**3//4)
            self.assertEqual(len(singletons), 2*side**3)
            self.assertEqual(face_count, 11*side**3//4)
            self.assertEqual({degrees[block] for block in selected}, {12})
            self.assertEqual({degrees[block] for block in singletons}, {4})
            for block, degree in degrees.items():
                row = sum(count for (first, _), count in pairs.items() if first == block)
                self.assertEqual(row, 3*degree)
                self.assertEqual(F(row, 4), 9 if block[0] == 'plaquette' else 3)
            replay = z.nonlinear_partition_geometry(side)
            self.assertEqual(replay['selected_four_link_blocks'], len(selected))
            self.assertEqual(replay['remaining_singleton_blocks'], len(singletons))
            self.assertEqual(replay['total_links'], len(owners))
            self.assertEqual(replay['nonselected_faces'], face_count)
            self.assertTrue(replay['every_nonselected_face_meets_four_distinct_blocks'])
            self.assertEqual(replay['selected_block_cross_degree'], 12)
            self.assertEqual(replay['singleton_cross_degree'], 4)
            self.assertEqual(replay['selected_block_mixed_row_per_abs_b'], 9)
            self.assertEqual(replay['singleton_mixed_row_per_abs_b'], 3)

    def test_doublewell_Haar_cap_and_Rayleigh_prefactor_are_consistent(self):
        # 4/pi * integral_0^(pi/6) (2r/pi)^2 dr = 16/(3*6^3).
        cap = F(16, 3*6**3)
        self.assertEqual(cap, F(2, 81))
        self.assertEqual(4/cap, 162)
        self.assertGreater(cap, F(1, 48))
        with localcontext() as context:
            context.prec = 70
            pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816')
            exact_cap = Decimal(1)/3-Decimal(3).sqrt()/(2*pi)
            self.assertLess(decimal(cap), exact_cap)
        row = z.quartic_doublewell_budget(4)
        rational_cap = 2*(F(15, 16)-F(7, 8))*F(5, 8)*F(11, 32)
        self.assertEqual(row['rational_native_Haar_cap_lower'], rational_cap)
        self.assertGreater(rational_cap, cap)
        self.assertEqual(row['simplified_cap_lower'], cap)

    def test_conditional_generator_cross_identity_matches_a_correlated_Gaussian(self):
        # Precision [[2,-1],[-1,5]], covariance [[5,1],[1,2]]/9.
        # For f=xy-E[xy], L1 f=2xy-y² and L2 f=5xy-x².
        xx, xy, yy = F(5, 9), F(1, 9), F(2, 9)
        x4, y4, x3y, xy3, x2y2 = 3*xx*xx, 3*yy*yy, 3*xx*xy, 3*yy*xy, xx*yy+2*xy*xy
        diagonal1 = 4*x2y2-4*xy3+y4
        diagonal2 = 25*x2y2-10*x3y+x4
        cross = 11*x2y2-2*x3y-5*xy3
        self.assertEqual(diagonal1, 2*yy)
        self.assertEqual(diagonal2, 5*xx)
        self.assertEqual(cross, 1-xy)
        self.assertEqual(cross, F(8, 9))
        self.assertEqual(diagonal1+diagonal2+2*cross, 5)
        # The negative cross-Hessian term, not a positive one, occurs in the identity.
        self.assertNotEqual(cross, 1+xy)
        self.assertGreater(diagonal1+diagonal2+2*cross, xx+yy)
        criterion = z.block_gap_matrix((2, 5), ((0, 1), (1, 0)), 1)
        self.assertEqual(criterion['comparison_matrix'], [[2, -1], [-1, 5]])
        self.assertTrue(criterion['accepted'])
        self.assertTrue(criterion['positive_gap_certified'])
        self.assertTrue(criterion['conditional_Poincare_and_cross_Hessian_premises_required'])
        self.assertFalse(z.block_gap_matrix((2, 5), ((0, 1), (1, 0)), 2)['accepted'])

    def test_parent_block_constants_do_not_depend_on_high_frequency_derivatives(self):
        A, bridge = F(1, 8), F(1, 192)
        conditional_oscillation = 2*A+24*bridge
        mixed_row_sum = 36*bridge/4
        self.assertEqual(conditional_oscillation, F(3, 8))
        self.assertEqual(mixed_row_sum, F(3, 64))
        self.assertEqual(F(3, 4)*(1-conditional_oscillation)-mixed_row_sum, F(27, 64))
        for n in (4, 17, 10000):
            # At W=pi/n, cos(nW)=-1 and its first derivative is zero.
            lower_Hessian = A*(n*n-10)-3*bridge
            self.assertGreaterEqual(lower_Hessian, F(47, 64))
            self.assertGreater(lower_Hessian, F(1, 2))
            row = z.parent_gap_budget(n)
            self.assertEqual(row['conditional_log_density_oscillation_upper'], conditional_oscillation)
            self.assertEqual(row['cross_Hessian_row_upper'], mixed_row_sum)
            self.assertEqual(row['elementary_parent_gap_lower'], F(27, 64))
            self.assertEqual(row['Hessian_S_direction_lower'], lower_Hessian)
            self.assertTrue(row['positive_gap_certified'])
            self.assertTrue(row['global_positive_Bakry_bound_fails'])
            self.assertTrue(row['gap_lower_independent_of_n'])
            self.assertFalse(row['identified_with_original_plaquette_Hamiltonian'])

    def test_quadratic_parameter_cell_and_directed_log_hand_floor_are_exact(self):
        endpoint = F(1, 16000)
        maximum = F(0)
        for theta_ratio, lam_ratio in product((F(-1), F(-1, 2), F(0), F(1, 2), F(1)), repeat=2):
            lo, hi = independent_quadratic_extrema(theta_ratio*endpoint, lam_ratio*endpoint)
            maximum = max(maximum, (hi-lo)/2)
        self.assertEqual(maximum, F(9, 128000))
        load = 4*maximum
        alpha = F(1, 2)+1000*load
        R = (1+1/alpha)/2
        self.assertEqual(alpha, F(25, 32))
        self.assertEqual(R, F(57, 50))
        self.assertEqual((2*(R-1)/(R+1))/50, F(7, 2675))
        with localcontext() as context:
            context.prec = 70
            self.assertGreater(decimal(R).ln()/50, decimal(F(7, 2675)))
        row = z.quadratic_gap_budget(endpoint, endpoint)
        self.assertEqual(row['interaction']['halfosc'], maximum)
        self.assertEqual(row['weighted_halfosc_load'], load)
        self.assertEqual(row['alpha_upper'], alpha)
        self.assertEqual(row['time_weight_ratio'], R)
        self.assertEqual(row['weighted_margin'], F(7, 64))
        self.assertGreater(row['gamma_interval']['lower'], F(7, 2675))
        self.assertTrue(row['full_quantum_Hilbert_space_gap'])
        self.assertTrue(row['finite_volume_uniform'])

    def test_quadratic_extrema_cover_sign_changes_and_interior_vertex(self):
        for theta, lam in product((F(-3), F(-1), F(0), F(1, 3), F(2), F(5)), repeat=2):
            low, high = independent_quadratic_extrema(theta, lam)
            row = z.quadratic_halfosc(theta, lam)
            self.assertEqual((row['minimum'], row['maximum']), (low, high))
            self.assertEqual(row['halfosc'], (high-low)/2)
            self.assertGreaterEqual(row['halfosc'], abs(theta))
            self.assertTrue(row['halfosc_dominates_linear_amplitude'])
            self.assertEqual(row['midpoint'], (high+low)/2)
            self.assertEqual(row['Lipschitz_upper'], abs(theta)+2*abs(lam))
            self.assertTrue(row['constant_shift_does_not_change_gap'])
            if row['critical_point'] is not None:
                self.assertEqual(theta+2*lam*row['critical_point'], 0)
                self.assertLessEqual(abs(row['critical_point']), 1)
        self.assertEqual(z.quadratic_halfosc(0, 0)['halfosc'], 0)

    def test_directed_log_intervals_enclose_independent_high_precision_log(self):
        with localcontext() as context:
            context.prec = 110
            for value in (F(1, 3), F(1), F(57, 50), F(3)):
                row = z.log_interval(value)
                actual = decimal(value).ln()
                self.assertLessEqual(decimal(row['lower']), actual)
                self.assertGreaterEqual(decimal(row['upper']), actual)
            for theta, lam in ((F(0), F(0)), (F(1, 16000), F(-1, 16000))):
                row = z.quadratic_gap_budget(theta, lam)
                actual = decimal(row['time_weight_ratio']).ln()/(10*row['J'])
                self.assertLessEqual(decimal(row['gamma_interval']['lower']), actual)
                self.assertGreaterEqual(decimal(row['gamma_interval']['upper']), actual)
        for bad in (0, -1):
            with self.assertRaises(ValueError): z.log_interval(bad)
        with self.assertRaises(ValueError): z.log_interval(2, terms=0)

    def test_threshold_failure_is_a_failed_sufficient_gate_and_not_gaplessness(self):
        threshold = F(1, 2000)
        for load in (threshold, 2*threshold):
            for J in range(1, 10):
                row = z.nonlinear_gap_budget(load, J)
                self.assertEqual(J*(row['alpha_upper']-1), F((J-5)**2, 10)+(200*load-F(1, 10))*J*J)
                self.assertFalse(row['accepted'])
                self.assertIsNone(row['gamma_interval'])
                self.assertIsNone(row['full_quantum_Hilbert_space_gap'])
                self.assertFalse(row['rejection_proves_gaplessness'])
        for load in (F(0), F(1, 4000)):
            row = z.nonlinear_gap_budget(load)
            self.assertTrue(row['accepted'])
            self.assertGreater(row['weighted_margin'], 0)
            self.assertTrue(row['actual_ground_measure_not_a_substituted_Gibbs_measure'])
        self.assertFalse(z.nonlinear_gap_budget(0, R=2)['accepted'])

    def test_parent_directed_interval_encloses_only_the_comparison_rate(self):
        with localcontext() as context:
            context.prec = 100
            rows = [z.parent_gap_budget(n) for n in (1, 4, 1000)]
            actual_rate = decimal(F(3, 4))*(-decimal(F(3, 8))).exp()-decimal(F(3, 64))
            for row in rows:
                self.assertLessEqual(decimal(row['parent_gap_interval']['lower']), actual_rate)
                self.assertGreaterEqual(decimal(row['parent_gap_interval']['upper']), actual_rate)
                self.assertTrue(row['interval_encloses_comparison_rate_not_actual_gap'])
            self.assertEqual(rows[0]['parent_gap_interval'], rows[-1]['parent_gap_interval'])
            self.assertIsNone(rows[0]['Hessian_S_direction_lower'])
            self.assertFalse(rows[0]['global_positive_Bakry_bound_fails'])
        self.assertFalse(z.parent_gap_budget(4, b=1)['positive_gap_certified'])

    def test_doublewell_exponential_upper_vanishes_without_claiming_finite_parameter_gap_zero(self):
        bounds = []
        with localcontext() as context:
            context.prec = 100
            for t in (F(0), F(1, 3), F(4), F(16), F(64)):
                row = z.quartic_doublewell_budget(t)
                actual = (-decimal(t)/2).exp()
                self.assertLessEqual(decimal(row['exponential_interval']['lower']), actual)
                self.assertGreaterEqual(decimal(row['exponential_interval']['upper']), actual)
                self.assertEqual(row['gauge_invariant_gap_upper'], 162*row['exponential_interval']['upper'])
                self.assertEqual(row['W4_coefficient'], t*t)
                self.assertEqual(row['W2_coefficient'], 4*t-t*t)
                self.assertEqual(row['exact_ground_energy'], -t)
                self.assertTrue(row['exact_vacuum_requires_one_or_edge_disjoint_active_plaquettes'])
                self.assertFalse(row['exact_vacuum_claimed_for_overlapping_all_plaquettes'])
                self.assertTrue(row['added_W_even_sector_would_exclude_trial'])
                self.assertFalse(row['zero_gap_at_any_fixed_finite_t_asserted'])
                bounds.append(row['gauge_invariant_gap_upper'])
        self.assertTrue(all(first > second for first, second in zip(bounds, bounds[1:])))
        self.assertLess(bounds[-1], F(1, 10**10))

    def test_exact_matrix_gate_handles_kernel_indefiniteness_and_invalid_premises(self):
        singular = z.block_gap_matrix((1, 1), ((0, 1), (1, 0)), 0)
        self.assertTrue(singular['accepted'])
        self.assertFalse(singular['positive_gap_certified'])
        self.assertFalse(z.block_gap_matrix((1, 1), ((0, 2), (2, 0)), 0)['accepted'])
        invalid = [((1,), ((1,),), 0), ((1, 1), ((0, 1), (0, 0)), 0),
                   ((1, 1), ((0,),), 0), ((0,), ((0,),), 0),
                   ((1, 1), ((0, -1), (-1, 0)), 0), ((1,), ((0,),), -1)]
        for gaps, cross, target in invalid:
            with self.assertRaises(ValueError): z.block_gap_matrix(gaps, cross, target)

    def test_native_quadratic_Ward_insertion_has_the_stated_constant_term(self):
        p, model = stage80.p, stage80.previous.model
        w = p.COORD[0]
        radius = p.add(*(p.mul(value, value) for value in p.COORD))
        w2 = p.mul(w, w)
        theta, lam = F(7, 9), F(-5, 11)
        potential = p.add(p.scale(w, theta), p.scale(w2, lam))
        half_LE = p.scale(model.electric(potential, 1), 2)
        expected = p.add(p.scale(w, F(3, 2)*theta), p.scale(w2, 4*lam), p.scale(radius, -lam))
        self.assertEqual(half_LE, expected)
        self.assertEqual(stage80.native.phi(half_LE, 1), 0)

    def test_nonphysical_parameter_inputs_are_rejected(self):
        for side in (2, 5):
            with self.assertRaises(ValueError): z.nonlinear_partition_geometry(side)
        for kwargs in ({'load': -1}, {'load': 0, 'J': 0}, {'load': 0, 'R': 1}):
            with self.assertRaises(ValueError): z.nonlinear_gap_budget(**kwargs)
        for kwargs in ({'n': 0}, {'n': 4, 'A': -1}, {'n': 4, 'b': -1}):
            with self.assertRaises(ValueError): z.parent_gap_budget(**kwargs)
        with self.assertRaises(ValueError): z.quartic_doublewell_budget(-1)

    def test_actual_strong_block_kernel_minorization_has_barrier_dependent_constants(self):
        rows = []
        with localcontext() as context:
            context.prec = 100
            for barrier in (F(0), F(1, 4), F(1, 2), F(8)):
                row = z.strong_block_gap_budget(barrier, 4, 12, 4, 0)
                if barrier == 0:
                    self.assertEqual(row['minorization_interval']['lower'], F(19, 21)**4)
                    self.assertEqual(row['minorization_interval']['upper'], F(19, 21)**4)
                else:
                    actual = (-10*decimal(barrier)).exp()*decimal(F(19, 21))**4
                    self.assertLessEqual(decimal(row['minorization_interval']['lower']), actual)
                    self.assertGreaterEqual(decimal(row['minorization_interval']['upper']), actual)
                self.assertGreater(row['epsilon_lower'], 0)
                self.assertLessEqual(row['epsilon_lower'], row['minorization_interval']['lower'])
                self.assertGreaterEqual(row['J']*row['epsilon_lower'], 4)
                self.assertEqual(row['alpha_upper'], 2/(row['epsilon_lower']*row['J']))
                self.assertLessEqual(row['alpha_upper'], F(1, 2))
                self.assertTrue(row['accepted'])
                self.assertTrue(row['actual_block_ground_transform_used'])
                self.assertFalse(row['internal_quantum_vacuum_assumed_Gibbs'])
                self.assertFalse(row['uniform_in_internal_amplitude'])
                self.assertTrue(row['smooth_potentials_and_finite_bridge_Gamma_required'])
                rows.append(row)
        self.assertTrue(all(first['epsilon_lower'] > second['epsilon_lower']
                            for first, second in zip(rows, rows[1:])))
        self.assertTrue(all(first['J'] < second['J'] for first, second in zip(rows, rows[1:])))

    def test_actual_strong_cosine_quantum_gap_is_uniform_in_frequency(self):
        epsilon, J, coupling = F(1, 32), 128, F(1, 1228800)
        temporal = 2/(epsilon*J)
        spatial = 40*12*(4+1)*coupling*J
        self.assertEqual(temporal, F(1, 2))
        self.assertEqual(spatial, F(1, 4))
        rows = [z.strong_cosine_gap_budget(n) for n in (1, 4, 64, 1024)]
        for row in rows:
            self.assertEqual(row['internal_oscillation_bound'], F(1, 4))
            self.assertEqual(row['max_block_links'], 4)
            self.assertEqual(row['bridge_degree'], 12)
            self.assertEqual(row['bridge_arity'], 4)
            self.assertEqual(row['epsilon_lower'], epsilon)
            self.assertEqual(row['J'], J)
            self.assertEqual(row['weak_coupling'], coupling)
            self.assertEqual(row['alpha_upper'], temporal+spatial)
            self.assertEqual(row['time_weight_ratio'], F(7, 6))
            self.assertEqual(row['weighted_margin'], F(1, 8))
            self.assertEqual(row['sufficient_weak_coupling_threshold'], F(1, 614400))
            self.assertGreater(row['gamma_interval']['lower'], F(1, 8320))
            self.assertTrue(row['full_quantum_Hilbert_space_gap'])
            self.assertTrue(row['finite_volume_uniform'])
            self.assertFalse(row['internal_quantum_vacuum_assumed_Gibbs'])
            self.assertFalse(row['interacting_ground_state_assumed_exp_S_over_2'])
            self.assertTrue(row['gap_bound_independent_of_fixed_frequency'])
        self.assertEqual(rows[0]['gamma_interval'], rows[-1]['gamma_interval'])
        with localcontext() as context:
            context.prec = 100
            rate = decimal(F(7, 6)).ln()/(10*J)
            self.assertLessEqual(decimal(rows[0]['gamma_interval']['lower']), rate)
            self.assertGreaterEqual(decimal(rows[0]['gamma_interval']['upper']), rate)

    def test_strong_block_gate_rejects_overstated_overlap_and_excessive_bridges(self):
        with self.assertRaises(ValueError):
            z.strong_block_gap_budget(F(1, 4), 4, 12, 4, 0, epsilon_lower=1)
        with self.assertRaises(ValueError): z.strong_block_gap_budget(-1, 4, 12, 4, 0)
        with self.assertRaises(ValueError): z.strong_block_gap_budget(0, 4, 12, 4, 0, J=0)
        row = z.strong_block_gap_budget(F(1, 4), 4, 12, 4, 1, J=128, epsilon_lower=F(1, 32))
        self.assertFalse(row['accepted'])
        self.assertIsNone(row['gamma_interval'])
        self.assertIsNone(row['full_quantum_Hilbert_space_gap'])
        self.assertFalse(row['uniform_in_internal_amplitude'])
        self.assertFalse(row['rejection_proves_gaplessness'])


class NonlinearGapCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_manifest_and_parent_bind_expected_sources(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 309)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(len(self.cert['input_sha256']), 314)
        self.assertEqual(manifest['base_tree'], '48e5eb6efd097809e649c57b5faa4db675392ae2')
        self.assertEqual(self.cert['parent_result_sha256'], 'e1720eeb6f2c6aed083e64bcbc8a9ad3bf712e4be5377d06bf0b9b08b13218e7')
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_changed_upstream_parent_or_base_is_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'BASE_TREE', '0'*40):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_continuum_parent_and_gaplessness_overclaims_are_rejected(self):
        for key in ('positive_continuum_mass_gap', 'spatial_continuum_Yang_Mills_constructed',
                    'unknown_quantum_vacuum_replaced_by_classical_Gibbs',
                    'parent_vacuum_identified_with_original_plaquette_vacuum',
                    'doublewell_exact_vacuum_applies_to_overlapping_all_plaquettes',
                    'all_nonlinear_models_have_one_uniform_positive_gap',
                    'strong_block_gap_uniform_in_internal_amplitude',
                    'rejected_sufficient_gate_proves_gaplessness'):
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

    def test_model_distinguishes_actual_nonlinear_quantum_gap_and_separate_parent(self):
        model = self.cert['model']
        self.assertTrue(model['actual_bounded_nonlinear_plaquette_gap_theorem'])
        self.assertTrue(model['actual_strong_internal_block_quantum_gap_theorem'])
        self.assertTrue(model['strong_block_gap_uniform_in_fixed_internal_frequency'])
        self.assertFalse(model['strong_block_gap_uniform_in_internal_amplitude'])
        self.assertTrue(model['actual_quantum_ground_measure_used'])
        self.assertTrue(model['weighted_half_oscillation_controls_volume_uniform_gap'])
        self.assertTrue(model['separately_defined_parent_has_exact_exp_S_over_2_vacuum'])
        self.assertTrue(model['parent_gap_is_uniform_in_oscillation_frequency'])
        self.assertFalse(model['global_Bakry_curvature_positivity_is_necessary'])
        self.assertFalse(model['upstream_linear_Yang_Mills_Hamiltonian_redefined'])
        self.assertFalse(model['perturbative_remainders_certified_uniform_in_volume'])
        self.assertFalse(model['finite_controls_are_universal_proof'])


if __name__ == '__main__':
    unittest.main()
