import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym53_anisotropic_interaction as z


class AnisotropicInteractionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.cert=z.run()

    def test_fresh_certificate_and_pin(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_all_content_tail_has_explicit_remainder(self):
        d=self.cert['all_content_heat']
        self.assertGreater(F(d['exp_9_over_2_lower']),84)
        self.assertLess(F(d['all_content_tail_upper']),F(1,20))
        self.assertEqual(d['tail_identity_checks'],25)
        self.assertGreater(z.geometric_tail(F(1,84)),F(4,84))

    def test_both_spin_parities_and_rank_two_are_checked(self):
        d=self.cert['spin_controls']
        self.assertEqual(d['finite_spin_inequalities'],108)
        self.assertEqual(d['half_integer_cases'],54)
        self.assertEqual(d['rank_two_cases'],24)

    def test_second_eigenvalue_not_the_smallest_is_the_gate(self):
        self.assertEqual(z.z.rates((0,F(3,2),F(3,2)))[0],F(3,4))
        self.assertIsNotNone(z.block_gate(F(3,2),F(1,4096),8,F(4,3)))
        with self.assertRaises(ValueError):z.block_gate(0,0,8,F(4,3))
        self.assertEqual(z.z.q.lap(z.z.EVEN,z.a.diag((0,0,3))),{})

    def test_mesh_bound_and_clock_covariant_ceil(self):
        self.assertEqual(self.cert['parameters']['fine_step_cases'],80)
        self.assertEqual(z.skeleton(F(3,2),F(2,3)),9)
        self.assertEqual(z.skeleton(F(3,2),F(1,7)),42)
        with self.assertRaises(ValueError):z.skeleton(F(3,2),1)

    def test_strict_window_and_budget_endpoint(self):
        self.assertEqual(self.cert['parameters']['endpoint_identity_checks'],150)
        for J in range(1,31):
            self.assertIsNone(z.block_gate(1,F(1,2400),J,F(10001,10000)))
        u=F(999,2400000);alpha=F(1,2)+1200*u
        self.assertIsNotNone(z.block_gate(1,u,5,(1+1/alpha)/2))

    def test_rank_two_interacting_cell_is_not_the_free_gap(self):
        d=self.cert['parameters']['rank_two_example']
        self.assertEqual(F(d['alpha']),F(5,8));self.assertEqual(F(d['margin']),F(1,6))
        self.assertGreater(F(d['gamma_lower']),F(539403,100000000))
        self.assertLess(F(d['gamma_upper']),F(3,4))
        self.assertIsNone(z.block_gate(F(3,2),F(1,16),8,F(4,3)))

    def test_zero_normalizer_is_rejected(self):
        with self.assertRaises(ValueError):z.short_bridge(0,3)
        self.assertEqual(sum(z.short_bridge(0,1)),1)
        self.assertEqual(sum(z.short_bridge(2,3)),1)
        self.assertTrue(self.cert['admissible_bridges']['zero_normalizer_rejected'])

    def test_grouped_endpoints_and_positive_tilt_support(self):
        d=self.cert['admissible_bridges']
        self.assertEqual(d['invalid_intermediate'],[0,3])
        self.assertLessEqual(F(d['grouped_endpoint_hamming_cost']),1)
        self.assertTrue(d['positive_tilt_preserves_support'])
        K=z.path_kernel()
        self.assertEqual(K[0][3],0)
        self.assertTrue(all(v>0 for row in z.a.multiply(z.a.multiply(K,K),K) for v in row))

    def test_bond_identity_includes_rotated_unequal_tensors(self):
        d=self.cert['bond_energy']
        self.assertEqual(d['coordinate_count'],8)
        self.assertEqual(d['site_gradient_identities'],2)
        self.assertEqual(d['inhomogeneous_tensor_pairs'],8)

    def test_refinement_reaches_the_normalized_gap_gate(self):
        d=self.cert['norm_refinement']
        self.assertEqual(d['moment_enclosures'],20)
        self.assertEqual(len(d['gap_transport_cells']),12)
        for row in d['gap_transport_cells']:
            self.assertLess(F(row['normalized_gap_ceiling_upper']),1)
            self.assertGreaterEqual(F(row['error_upper']),0)
        self.assertFalse(d['uniform_in_unbounded_trace'])

    def test_refinement_needs_trace_and_recovers_free_case(self):
        self.assertEqual(z.refinement_budget(1,F(1,4096),3,1,F(1,16)),0)
        self.assertEqual(z.refinement_budget(4,0,3,1,F(1,16)),0)
        low=z.refinement_budget(4,F(1,4096),3,1,F(1,16))
        high=z.refinement_budget(4,F(1,4096),12,1,F(1,16))
        self.assertGreater(high,low)
        self.assertIsNone(z.refine.gap_gate(F(1),F(1,2),F(1,4)))

    def test_ground_weight_and_cross_terms_are_retained(self):
        d=self.cert['ground_source_energy']
        self.assertEqual(d['weighted_form_matrix_entries'],16)
        self.assertEqual(d['ground_kernel_dimension'],1)
        self.assertTrue(d['omitting_ground_weight_rejected'])

    def test_existing_coupling_controls_are_reused(self):
        d=self.cert['reused_YM45_controls']
        self.assertEqual(d['bridges']['skeleton_bridge_cases'],18)
        self.assertEqual(d['tilts']['tilted_block_cases'],20)
        self.assertEqual(d['incidence']['clipped_layouts'],75)
        self.assertEqual(len(d['joint_updates']['cases']),2)

    def test_exact_domains_and_clock_scaling(self):
        self.assertEqual(self.cert['clock']['clock_covariant_cells'],9)
        for q,N in ((1,0),(-1,0),(F(1,2),-1),(0.5,2),(F(1,2),True)):
            with self.assertRaises(ValueError):z.geometric_tail(q,N)
        for beta,theta,J,R in ((0,0,5,2),(1,0,True,2),(1,0,5,1),(1,0.1,5,2)):
            with self.assertRaises(ValueError):z.block_gate(beta,theta,J,R)
        with self.assertRaises(ValueError):z.refinement_budget(2,1,0,1,F(1,4))

    def test_verdict_and_hash_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=Path(tmp)/'result';pin=Path(tmp)/'pin'
            changed=copy.deepcopy(self.cert);changed['parameters']['window']='unconditional'
            result.write_text(json.dumps(changed));pin.write_text(z.old.canonical_sha(changed))
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)
            result.write_text(json.dumps(self.cert));pin.write_text('0'*64)
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)

    def test_upstream_change_is_rejected(self):
        original=z.digest
        try:
            z.digest=lambda path:'0'*64
            with self.assertRaises(ValueError):z.source_checks()
        finally:z.digest=original

    def test_scope_does_not_promote_joint_or_physical_limit(self):
        self.assertIn('FIXED_WIDTH',self.cert['claim_status'])
        self.assertIn('JOINT_LIMIT',self.cert['claim_status'])
        self.assertIn('four-dimensional',self.cert['evidence_scope']['open'])
        self.assertFalse(self.cert['clock']['physical_clock_selected'])


if __name__=='__main__':unittest.main()
