import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym66_local_vacuum_limit as z


class LocalVacuumLimitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.eta = F(1,16000)

    def test_fresh_certificate_and_pin_match(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_python_312_is_enforced(self):
        with patch.object(z.native.previous.sys,'version_info',(3,11,0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_all_197_upstream_files_and_four_local_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']),197)
        self.assertEqual(len(manifest['local_inputs']),4)
        self.assertIn('papers/yang-mills-certified-benchmark/YM65_SPATIAL_RESPONSE_LOCALITY.md',
                      manifest['upstream_sha256'])

    def test_changed_upstream_and_rehashed_overclaims_are_rejected(self):
        with patch.object(z.native,'digest',return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            wrong = copy.deepcopy(self.cert)
            wrong['model']['pointwise_conditional_density_control'] = True
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError):
                z.check(self.cert,result,pin)

    def test_no_binary_floats_enter_the_serialized_evidence(self):
        def walk(value):
            self.assertNotIsInstance(value,float)
            if isinstance(value,dict):
                for item in value.values():
                    walk(item)
            elif isinstance(value,(list,tuple)):
                for item in value:
                    walk(item)
        walk(self.cert)

    def test_signed_and_deleted_face_paths_preserve_the_incidence_budget(self):
        controls = self.cert['profiles']
        self.assertEqual(F(controls['alpha']),F(3,4))
        self.assertEqual(F(controls['margin']),F(1,8))
        self.assertEqual({F(row['path_parameter']) for row in controls['profiles']},
                         {F(0),F(1,3),F(1)})
        for row in controls['profiles']:
            self.assertLessEqual(F(row['weighted_row_max']),16*self.eta)
            self.assertGreater(F(row['potential_norm_budget']),0)

    def test_strict_coupling_window_and_positive_resolution_are_required(self):
        for value in (F(-1,16000),F(1,8000),F(1),0.00001,True):
            with self.assertRaises(ValueError):
                z.parameters(value)
        with self.assertRaises(ArithmeticError):
            z.parameters(F(1,8000)-F(1,10**40))

    def test_profile_difference_cannot_exceed_the_cube_diameter(self):
        self.assertEqual(z.change_budget(self.eta,2*self.eta),2*self.eta)
        with self.assertRaises(ValueError):
            z.change_budget(self.eta,3*self.eta)
        with self.assertRaises(ValueError):
            z.change_budget(0,F(1,100000))

    def test_native_cut_matches_an_independently_constructed_inner_box(self):
        cut = self.cert['geometry']
        self.assertEqual(cut['interior_links'],300)
        self.assertEqual(cut['interior_faces'],240)
        self.assertGreater(cut['removed_faces'],cut['interior_faces'])
        self.assertGreaterEqual(cut['distance_to_removed_faces'],1)
        self.assertTrue(cut['discarded_links_remain_free'])

    def test_shell_count_covers_open_and_wrapped_distances(self):
        rows = self.cert['geometry']['shells']
        self.assertEqual({row['periodic'] for row in rows},{False,True})
        for row in rows:
            for d,count in row['shell_counts']:
                self.assertLessEqual(count,3*row['support_links']*(2*d+5)**3)

    def test_tail_includes_the_endpoint_and_has_the_exact_remainder(self):
        for q,R in ((F(1,2),1),(F(1,4),9),(F(3,4),5)):
            tail = z.tail_moments(q,R)
            farther = z.tail_moments(q,R+20)
            for k in range(5):
                direct = sum(F(d**k)*q**d for d in range(R,R+20))
                self.assertEqual(tail[k],direct+farther[k])
                self.assertGreater(tail[k],z.tail_moments(q,R+1)[k])

    def test_directed_tail_rounding_never_undershoots_exact_arithmetic(self):
        q,R = F(3,4),31
        for exact,upper in zip(z.tail_moments(q,R),z.tail_moments(q,R,rounded=True)):
            self.assertGreaterEqual(upper,exact)
        self.assertEqual(z.tail_moments(0,8),(0,0,0,0,0))
        self.assertEqual(z.tail_moments(F(1,2),1),(1,2,6,26,150))

    def test_invalid_tail_and_support_arguments_are_refused(self):
        for args in ((1,2),(F(-1,2),2),(0.5,2),(F(1,2),0),(F(1,2),True)):
            with self.assertRaises(ValueError):
                z.tail_moments(*args)
        for R,links in ((0,4),(1,0),(True,4),(1,False)):
            with self.assertRaises(ValueError):
                z.boundary_tail(self.eta,R,links)

    def test_trace_distance_boundary_budgets_are_nontrivial_at_large_buffer(self):
        rows = self.cert['boundary_budgets']
        values = [F(row['bound']) for row in rows]
        self.assertEqual(values[0],2)
        self.assertTrue(all(left > right for left,right in zip(values,values[1:])))
        for value,upper in zip(values[1:],(F(200366,10**6),F(4121,10**9),F(5490,10**14))):
            self.assertLess(value,upper)
        self.assertGreater(values[-1],F(5489,10**14))

    def test_boundary_budget_keeps_the_full_fine_support_and_change_size(self):
        one = z.boundary_tail(self.eta,256,4)['uncapped_upper']
        two = z.boundary_tail(self.eta,256,8)['uncapped_upper']
        doubled = z.boundary_tail(self.eta,256,4,2*self.eta)['uncapped_upper']
        self.assertLessEqual(abs(two-2*one),F(2,10**30))
        self.assertLessEqual(abs(doubled-2*one),F(2,10**30))

    def test_zero_change_and_free_model_give_zero_comparison(self):
        self.assertEqual(z.boundary_tail(self.eta,1,100,0)['bound'],0)
        self.assertEqual(z.boundary_tail(0,1,100)['bound'],0)
        self.assertEqual(z.dynamic_tail(0,1,100,1000),0)
        self.assertEqual(z.dynamic_tail(self.eta,1,100,0),0)

    def test_source_derivative_sign_matches_an_independent_energy_eigenvector(self):
        data = self.cert['source']
        self.assertEqual(F(data['derivative']),F(-1176,15625))
        self.assertEqual(F(data['integrated_readout_change']),F(-18,25))
        self.assertLess(F(data['derivative']),0)

    def test_quantum_trace_distance_and_classical_L1_are_distinct(self):
        data = self.cert['source']
        self.assertEqual(F(data['quantum_trace_distance']),F(6,5))
        self.assertEqual(F(data['classical_L1_distance']),F(18,25))
        self.assertLess(F(data['classical_L1_distance']),F(data['quantum_trace_distance']))

    def test_real_time_comparison_decays_in_buffer_and_grows_in_time(self):
        bounds = [z.dynamic_tail(self.eta,R,4,100) for R in (16,32,64)]
        self.assertTrue(all(x > y for x,y in zip(bounds,bounds[1:])))
        self.assertLess(bounds[-1],F(1,10**12))
        self.assertLess(z.dynamic_tail(self.eta,32,4,10),z.dynamic_tail(self.eta,32,4,100))
        self.assertEqual(z.dynamic_tail(self.eta,32,4,-100),z.dynamic_tail(self.eta,32,4,100))
        self.assertEqual(z.dynamic_tail(self.eta,1,4,10000),2)

    def test_real_time_word_requires_cumulative_times(self):
        data = self.cert['time']
        self.assertEqual(F(data['real_word']),-1)
        self.assertEqual(F(data['incorrect_reset_word']),1)
        self.assertTrue(data['negative_energy_substitution_rejected'])

    def test_poisson_tails_are_normalized_and_have_directed_bounds(self):
        for row in self.cert['time']['poisson_tails']:
            actual = row['tail']
            self.assertGreater(F(actual['lo']),0)
            self.assertLess(F(actual['hi']),F(row['upper']))
            self.assertLessEqual(F(row['upper']),1)
        self.assertEqual(z.poisson_tail_upper(0,1),0)
        with self.assertRaises(ValueError):
            z.poisson_tail_upper(1,0)

    def test_native_product_pushforward_contracts_total_variation(self):
        data = self.cert['marginals']
        self.assertEqual(F(data['original_tv']),F(1,6))
        self.assertEqual(data['product_link_tv'],data['original_tv'])
        self.assertEqual(F(data['forgotten_tv']),0)

    def test_conditioning_needs_a_nonvanishing_probability_budget(self):
        for row in self.cert['marginals']['rare_conditioning']:
            p0,tv,after = map(F,(row['probability'],row['unconditional_tv'],row['conditional_tv']))
            self.assertEqual(after,F(1,4))
            self.assertEqual(tv,p0/4)
            self.assertGreater(after,tv)
            self.assertLessEqual(after,z.conditioning_bound(tv,p0))
        self.assertEqual(z.conditioning_bound(F(1,10),F(1,100)),1)
        for epsilon,p0 in ((F(1,10),0),(2,F(1,2)),(F(1,10),2)):
            with self.assertRaises(ValueError):
                z.conditioning_bound(epsilon,p0)

    def test_two_atom_tower_recovers_energies_and_weights(self):
        data = self.cert['moments']
        self.assertEqual(list(map(F,data['inverse_energy_polynomial'])),[F(1,10),F(-7,10),F(1)])
        self.assertEqual(list(map(F,data['recovered_energies'])),[2,5])
        self.assertEqual(list(map(F,data['recovered_weights'])),[F(1,3),F(2,3)])

    def test_hankel_log_convexity_drives_monotone_visible_energy_bounds(self):
        data = self.cert['moments']
        m = list(map(F,data['inverse_moments']))
        ratios = list(map(F,data['energy_ratio_upper_bounds']))
        self.assertTrue(all(m[n]*m[n+2] > m[n+1]**2 for n in range(len(m)-2)))
        self.assertTrue(all(x > y > 2 for x,y in zip(ratios,ratios[1:])))
        self.assertLess(ratios[-1],F(10001,5000))

    def test_complete_probe_tower_does_not_identify_an_invisible_lower_level(self):
        data = self.cert['moments']
        self.assertEqual(F(data['invisible_full_gap']),1)
        self.assertEqual(F(data['visible_threshold']),2)
        self.assertEqual(z.spectral_moments((1,2,5),(0,F(1,3),F(2,3)),8),
                         z.spectral_moments((2,5),(F(1,3),F(2,3)),8))

    def test_native_free_plaquette_has_one_exact_excitation_energy(self):
        m = list(map(F,self.cert['moments']['free_plaquette_moments']))
        self.assertEqual(m[0],F(1,4))
        self.assertEqual({x/y for x,y in zip(m,m[1:])},{3})

    def test_a_finite_moment_prefix_does_not_generally_fix_the_spectrum(self):
        # Both positive measures have m0=1,m1=3/10; their second moments differ.
        mixed = z.spectral_moments((2,5),(F(1,3),F(2,3)),2)
        single = z.spectral_moments((F(10,3),),(1,),2)
        self.assertEqual(mixed[:2],single[:2])
        self.assertNotEqual(mixed[2],single[2])

    def test_spectral_measure_rejects_nonpositive_energies_and_signed_weights(self):
        for energies,weights in (((0,2),(1,1)),((1,2),(-1,2)),((),()),((1,),(0,)),
                                  ((1,2),(1,)),((1.0,),(1,))):
            with self.assertRaises(ValueError):
                z.spectral_moments(energies,weights,3)

    def test_clock_rescaling_changes_visible_energies_and_inverse_orders(self):
        base = z.spectral_moments((2,5),(F(1,3),F(2,3)),5)
        scaled = z.spectral_moments((6,15),(F(1,3),F(2,3)),5)
        for n,value in enumerate(base):
            self.assertEqual(scaled[n],value/F(3)**n)
        self.assertEqual(scaled[3]/scaled[4],3*base[3]/base[4])

    def test_scope_distinguishes_unique_admissible_limits_from_pointwise_conditioning(self):
        model = self.cert['model']
        for key in ('inhomogeneous_gap_cube','local_quantum_state_trace_norm_limit',
                    'admissible_vacuum_limit_unique','heat_history_limit_unique',
                    'fixed_coarse_marginal_L1_limit','fixed_inverse_moments_unique',
                    'probe_spectral_measure_determined_by_full_tower'):
            self.assertTrue(model[key])
        for key in ('all_algebraic_ground_states_classified','pointwise_conditional_density_control',
                    'Schur_memory_kernel_locality','time_zero_Markov_identified',
                    'interacting_spectral_measure_computed','continuum_R4_QFT','physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
