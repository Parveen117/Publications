import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym68_averaged_exterior_locality as z


class AveragedExteriorLocalityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.eta = F(1,16000)

    def test_fresh_certificate_and_hash_match(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_runtime_requires_python_312(self):
        with patch.object(z.native.previous.sys, 'version_info', (3,11,0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_211_upstream_files_and_four_local_inputs_are_pinned(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 211)
        self.assertEqual(len(manifest['local_inputs']), 4)
        for stage in ('YM65_SPATIAL_RESPONSE_LOCALITY', 'YM67_POINTWISE_COARSE_VACUUM'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+stage+'.md',
                          manifest['upstream_sha256'])

    def test_source_changes_are_detected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()

    def test_rehashed_worst_case_claim_does_not_match_fresh_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            wrong = copy.deepcopy(self.cert)
            wrong['model']['worst_case_exterior_configuration_locality'] = True
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_certificate_uses_exact_scalars_and_directed_intervals(self):
        def walk(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for item in value.values():
                    walk(item)
            elif isinstance(value, (list,tuple)):
                for item in value:
                    walk(item)
        walk(self.cert)

    def test_product_vacuum_has_no_exterior_information_or_score_effect(self):
        cell = z.exterior_budget(0, 1, 4, 2)
        for name in ('covariance', 'expected_L1', 'expected_TV',
                     'mutual_information', 'expected_relative_Fisher'):
            self.assertEqual(cell[name], 0)
        self.assertEqual(z.smoothing_budget(self.eta, 4, 0)['bound'], 0)

    def test_invalid_couplings_sizes_and_covariances_are_rejected(self):
        for eta in (-1, F(1,8000), 0.00001, True):
            with self.assertRaises(ValueError):
                z.exterior_budget(eta, 10, 4, 2)
        for args in ((self.eta,0,4,2), (self.eta,1,0,2), (self.eta,1,4,0),
                     (self.eta,1,True,2)):
            with self.assertRaises(ValueError):
                z.exterior_budget(*args)
        for delta in (-1, F(101,100), 0.1, True):
            with self.assertRaises(ValueError):
                z.smoothing_budget(self.eta, 4, delta)

    def test_unresolved_near_endpoint_gap_is_not_silently_used(self):
        with self.assertRaises(ArithmeticError):
            z.exterior_budget(F(1,8000)-F(1,10**40), 10, 4, 2)

    def test_native_moments_match_independent_catalan_formula(self):
        for k in range(15):
            monomial = {(2*k,0,0,0): F(1)}
            self.assertEqual(z.y.phi(monomial), z.coordinate_even_moment(k))

    def test_native_laplacian_and_energy_normalizations(self):
        X = z.p.COORD[0]
        self.assertEqual(z.laplace_poly(X), z.p.scale(X,F(-3,4)))
        self.assertEqual(z.y.phi(z.gamma_poly(X,X)), F(3,16))
        self.assertEqual(self.cert['native']['weighted_integration_by_parts_checks'],9)

    def test_symmetric_ball_error_is_quadratic_in_radius(self):
        cell = z.previous.parameters(self.eta)
        row = z.exterior_budget(self.eta, 512, 4, 2)
        r = row['smoothing']['radius']
        K = z.previous.power_upper(cell['kappa'].hi, 4)
        rhs = 40**4*row['covariance']/r**12 + 4*cell['density_second_ratio'].hi*K*r*r
        self.assertGreaterEqual(row['expected_L1'], rhs)
        self.assertLessEqual(row['expected_L1']-rhs, F(1,10**120))
        self.assertEqual(row['analytic_exponent'], F(1,7))

    def test_covariance_duality_has_no_hidden_factor_two(self):
        # Exact non-product 3-by-4 law, independent of any YM discretization.
        table = [[F(v,48) for v in line] for line in ((1,2,3,4),(5,2,4,1),(3,7,8,8))]
        py = [sum(line) for line in table]
        pz = [sum(table[i][j] for i in range(3)) for j in range(4)]
        self.assertEqual(sum(py), 1)
        for mask in range(1,7):
            rows = [i for i in range(3) if mask & (1<<i)]
            pa = sum(py[i] for i in rows)
            H = [sum(table[i][j] for i in rows)/pz[j]-pa for j in range(4)]
            g = [1 if h>=0 else -1 for h in H]
            cov = sum(table[i][j]*g[j] for i in rows for j in range(4))
            cov -= pa*sum(pz[j]*g[j] for j in range(4))
            self.assertEqual(cov, sum(pz[j]*abs(H[j]) for j in range(4)))

    def test_expected_TV_is_half_L1(self):
        for row in self.cert['exterior_budgets']:
            self.assertGreaterEqual(F(row['expected_TV']),F(row['expected_L1'])/2)
            self.assertLessEqual(F(row['expected_TV'])-F(row['expected_L1'])/2,F(1,10**120))

    def test_density_range_caps_close_exteriors(self):
        row = z.exterior_budget(self.eta, 1, 4, 2)
        self.assertEqual(row['smoothing']['method'], 'density_range')
        self.assertIsNone(row['smoothing']['radius'])
        self.assertLess(row['expected_L1'], F(15,100))
        self.assertLessEqual(row['expected_L1'], 2)

    def test_arbitrarily_many_pins_do_not_enter_target_dimension(self):
        # The inherited covariance input uses bm, an upper bound on the
        # smaller support for ANY exterior, and never b*(m+|Z|).
        row = z.exterior_budget(self.eta, 512, 4, 2)
        self.assertEqual(row['covariance'],z.previous.covariance_upper(self.eta,512,8))
        self.assertEqual(row['fine_target_links'],8)
        self.assertEqual(row['analytic_exponent'],F(1,7))
        larger = z.exterior_budget(self.eta,512,5,2)
        self.assertEqual(larger['analytic_exponent'],F(2,17))

    def test_long_distance_bounds_have_the_claimed_upward_orientation(self):
        row = z.exterior_budget(self.eta,512,4,2)
        self.assertLess(row['expected_TV'],F(134436,10**12))
        self.assertLess(row['mutual_information'],F(39644,10**12))
        self.assertLess(row['expected_relative_Fisher'],F(96216,10**13))
        far = z.exterior_budget(self.eta,1024,4,2)
        self.assertGreater(far['expected_L1'],0)
        self.assertLess(far['expected_L1'],F(123178,10**18))

    def test_averaged_variation_and_information_decrease_with_distance(self):
        for key in ('expected_L1','expected_TV','mutual_information'):
            values = [F(row[key]) for row in self.cert['exterior_budgets']]
            self.assertTrue(all(a>b for a,b in zip(values,values[1:])))
        values = [F(row['expected_relative_Fisher']) for row in self.cert['exterior_budgets'][1:]]
        self.assertTrue(all(a>b for a,b in zip(values,values[1:])))

    def test_Fisher_comparison_uses_three_native_directions_and_link_norms(self):
        row = z.exterior_budget(self.eta,512,4,2)
        cell = z.previous.parameters(self.eta)
        coefficient = 24*cell['log_second'].hi+8*cell['score'].hi**2
        self.assertEqual(row['Fisher_coefficient'],z.previous.up(coefficient))
        self.assertGreaterEqual(row['expected_relative_Fisher'],coefficient*row['expected_L1'])

    def test_nonpolynomial_Fisher_integrals_have_independent_exact_values(self):
        rows = self.cert['Fisher']['cases']
        self.assertEqual(F(rows[0]['exact_Fisher']),F(13,180))
        self.assertEqual(F(rows[1]['exact_Fisher']),F(299,40400))
        for row in rows:
            value = F(row['exact_Fisher'])
            self.assertLessEqual(F(row['series_enclosure']['lo']),value)
            self.assertGreaterEqual(F(row['series_enclosure']['hi']),value)
            self.assertLessEqual(value,F(row['comparison_upper']))

    def test_Fisher_series_tail_shrinks_without_float_quadrature(self):
        coarse = z.fisher_series(F(3,5),5)
        fine = z.fisher_series(F(3,5),20)
        self.assertGreater(fine.lo,coarse.lo)
        self.assertLess(fine.hi,coarse.hi)
        self.assertEqual(z.fisher_series(0,1).lo,0)
        with self.assertRaises(ValueError):
            z.fisher_series(1,20)

    def test_exceptional_event_bound_is_on_probability_not_every_pin(self):
        self.assertEqual(z.exceptional_probability(F(1,100),F(1,10)),F(1,10))
        self.assertEqual(z.exceptional_probability(2,F(1,10)),1)
        self.assertEqual(z.exceptional_probability(0,1),0)
        self.assertLess(F(self.cert['exceptional_L1_probability_at_512']),F(269,10**6))
        with self.assertRaises(ValueError):
            z.exceptional_probability(1,0)

    def test_rare_native_pin_has_small_average_but_nonvanishing_exact_score(self):
        rows = self.cert['limitations']['rare_pins']
        self.assertLess(F(rows[-1]['average_L1_upper']),F(1,10**38))
        self.assertGreater(F(rows[-1]['exact_pin_score']),F(499,100000))
        self.assertTrue(all(F(row['density_floor'])==F(99,100) for row in rows))
        self.assertTrue(all(F(row['hidden_variance'])>0 for row in rows))

    def test_rare_pin_variance_matches_direct_native_integration(self):
        # For two exterior links, expand g=(1+z0)(1+w0)/4-1/4.
        Z,W = (1,0,0,0,0,0,0,0),(0,0,0,0,1,0,0,0)
        ZW = tuple(a+b for a,b in zip(Z,W))
        g = {Z:F(1,4),W:F(1,4),ZW:F(1,4)}
        self.assertEqual(z.native.phi(g,2),0)
        self.assertEqual(z.native.phi(z.p.mul(g,g),2),z.rare_pin_case(2)['hidden_variance'])

    def test_nearby_pin_creates_dependence_in_exact_three_link_polynomial(self):
        X = (1,0,0,0)+(0,)*8
        Z = (0,)*8+(1,0,0,0)
        XYZ = (1,0,0,0)*3
        epsilon = F(1,100)
        w = {(0,)*12:F(1),XYZ:epsilon}
        xz = {tuple(a+b for a,b in zip(X,Z)):F(1)}
        self.assertEqual(z.native.phi(z.p.mul(w,xz),3),0)
        pinned = {}
        for monomial,c in w.items():
            self.assertFalse(any(monomial[5:8]))
            key = monomial[:4]+(0,)*4+monomial[8:]
            pinned[key] = pinned.get(key,0)+c
        self.assertEqual(z.native.phi(z.p.mul(pinned,xz),3),epsilon/16)

    def test_scope_distinguishes_averaged_infinite_pins_from_Schur_and_continuum(self):
        model = self.cert['model']
        self.assertTrue(model['averaged_bounds_independent_of_exterior_cardinality'])
        self.assertTrue(model['countably_infinite_exterior_conditionals_C1_almost_surely'])
        self.assertFalse(model['worst_case_exterior_configuration_locality'])
        self.assertFalse(model['uniform_mixing_after_nearby_conditioning'])
        self.assertFalse(model['full_Schur_memory_kernel_locality'])
        self.assertFalse(model['time_zero_Markov_identified'])
        self.assertFalse(model['continuum_R4_QFT'])


if __name__ == '__main__':
    unittest.main()
