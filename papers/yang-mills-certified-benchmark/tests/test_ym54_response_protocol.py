import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym54_response_protocol as z


class ResponseProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.cert=z.run()

    def test_fresh_certificate_and_pin(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_compact_embedding_and_dagger_have_independent_engine_replay(self):
        d=self.cert['native_embedding']
        self.assertEqual(d['quaternion_product_checks'],16)
        self.assertEqual(d['dagger_checks'],4)
        self.assertTrue(d['unlifted_shape_refused_as_skew_turn'])
        self.assertEqual(z.mm(z.K,z.K),z.I)
        self.assertEqual(z.mm(z.R,z.R),z.a.scale(z.I,-1))

    def test_metric_adjoint_is_not_ordinary_adjoint(self):
        X=z.mm(z.a.inverse(z.a.diag((2,1))),z.L)
        self.assertEqual(z.trace(z.mm(X,X)),1)
        self.assertEqual(z.trace(z.mm(z.a.transpose(X),X)),F(5,4))
        bad=z.mm(z.K,z.L)
        self.assertEqual(z.trace(z.mm(bad,bad)),-2)

    def test_scalar_trace_is_present_and_removed_from_shape(self):
        d=z.response_data(z.I,[z.I,z.K],z.I)
        self.assertEqual([z.trace(X) for X in d['X']],[2,0])
        full=[[z.trace(z.mm(x,w)) for w in d['X']] for x in d['X']]
        self.assertEqual(full,[[2,0],[0,2]])
        self.assertEqual(d['G'],[[0,0],[0,2]])
        self.assertEqual(d['marker'],0)

    def test_curvature_area_and_floor_cover_noncommuting_and_degenerate_jets(self):
        d=self.cert['response_geometry']
        self.assertEqual(d['response_jets'],48)
        self.assertEqual(d['rank_two_floor_checks'],24)
        self.assertEqual(d['degenerate_shape_cases'],24)
        f=z.fixture()
        self.assertEqual(z.det2(f['G']),16*f['marker']**2)
        self.assertEqual(z.det2(f['G']),8*f['E'])

    def test_shape_gram_and_protocol_have_distinct_size_and_normalization(self):
        d=z.from_shapes(((1,0),(0,2)))
        self.assertEqual(d['G'],z.a.diag((2,8)))
        self.assertEqual(d['C'],z.a.diag((F(1,2),2,0)))
        self.assertEqual(z.a.y37.inertia(d['C']),(2,0,1))
        self.assertEqual(d['m2'],d['tau']/4)

    def test_actual_signed_records_replay_the_generator(self):
        d=self.cert['counted_protocol']
        self.assertEqual(d['four_label_moment_checks'],4)
        self.assertEqual(d['independent_generator_checks'],16)
        self.assertEqual(d['finite_coefficient_positive_energy_checks'],16)
        self.assertEqual(d['degree_range'],[0,3])

    def test_factor_covariance_and_response_orientation(self):
        self.assertEqual(self.cert['covariance']['factor_and_direction_covariances'],6)
        d=z.fixture()
        reflected=z.response_data(d['H'],d['derivatives'],z.mm(z.K,d['factor']),2)
        self.assertEqual(reflected['marker'],d['marker'])
        reversed_pair=z.response_data(d['H'],d['derivatives'][::-1],d['factor'],2)
        self.assertEqual(reversed_pair['marker'],-d['marker'])
        self.assertEqual(reversed_pair['C'],d['C'])

    def test_same_heat_keeps_distinct_finite_event_moments(self):
        d=self.cert['covariance']['same_heat_different_finite_records']
        self.assertEqual(F(d['first_m4']),F(17,2))
        self.assertEqual(F(d['second_m4']),F(8033,1250))
        self.assertEqual(F(d['scalar_fourth_coefficient_difference']),F(27,5000))

    def test_actual_finite_turns_satisfy_refinement_budget(self):
        d=self.cert['refinement']
        self.assertEqual(d['independent_finite_turn_enclosures'],36)
        for row in d['cells']:
            self.assertLessEqual(F(row['linear_error_upper']),F(row['written_degree_one_budget']))
        self.assertEqual(z.cos_sqrt_interval(0).lo,1)
        for x,N in ((-1,20),(2,20),(0.5,20),(F(1,2),True)):
            with self.assertRaises(ValueError):z.cos_sqrt_interval(x,N)

    def test_original_potential_derivatives_reproduce_fixture_and_flat_control(self):
        d=self.cert['thermo_fixture']
        self.assertTrue(d['potential_derivative_replay'])
        self.assertTrue(d['point_H_does_not_select_protocol'])
        self.assertEqual(F(d['marker']),F(-5,108))
        self.assertEqual(F(d['tau']),F(65,162))
        self.assertEqual(F(d['beta_from_curvature']),F(5,234))

    def test_curvature_bridge_gives_conservative_interacting_rate(self):
        d=self.cert['thermo_fixture']
        self.assertEqual(d['sufficient_window'],'abs(theta) < 1/112320')
        self.assertGreaterEqual(F(d['gamma_lower']),F(79467076,10**12))
        gate=z.source_gate(z.fixture(),F(5,108),F(65,162),F(1,262144),6,F(5,4))
        self.assertGreater(gate['gamma'].lo,F(79467076,10**12))
        self.assertLess(F(d['gamma_upper']),F(d['free_full_rate']))
        self.assertEqual(F(d['free_full_rate']),F(65,2592))
        self.assertEqual(F(d['free_even_rate']),F(5,162))
        self.assertNotEqual(F(d['free_full_rate']),F(d['tau'])/2)

    def test_source_gate_checks_actual_bounds_and_strict_endpoint(self):
        d=z.fixture();f0=abs(d['marker']);T=d['tau']
        self.assertIsNotNone(z.source_gate(d,f0,T,F(1,262144),6,F(5,4)))
        self.assertIsNone(z.source_gate(d,f0,T,F(1,112320),5,F(10001,10000)))
        for f,budget in ((2*f0,2*T),(f0/2,T/2)):
            with self.assertRaises(ValueError):z.source_gate(d,f,budget,0,6,F(5,4))

    def test_independence_alone_has_no_uniform_floor(self):
        rows=self.cert['scale_and_failures']['independent_but_slow']
        self.assertEqual(len(rows),3)
        self.assertTrue(all(F(row['half_gram_trace'])>=1 for row in rows))
        self.assertEqual(F(rows[-1]['second_eigenvalue']),F(1,20000))

    def test_fixed_curvature_requires_an_upper_response_budget(self):
        rows=self.cert['scale_and_failures']['fixed_curvature_without_budget']
        self.assertTrue(all(F(row['marker'])==F(1,2) for row in rows))
        self.assertGreater(F(rows[-1]['tau']),F(rows[0]['tau']))
        self.assertLess(F(rows[-1]['second_eigenvalue']),F(rows[0]['second_eigenvalue']))

    def test_idle_counts_prevent_unique_physical_clock_selection(self):
        rows=self.cert['scale_and_failures']['fixed_response_idle_protocols']
        active=z.fixture()['m2']
        for row in rows:self.assertEqual(F(row['trace']),F(row['active_weight'])*active)
        self.assertFalse(self.cert['scale_and_failures']['physical_protocol_and_clock_selected'])

    def test_scalar_response_invariance_and_clock_covariance(self):
        self.assertEqual(self.cert['scale_and_failures']['scalar_and_clock_scale_checks'],3)
        self.assertTrue(self.cert['scale_and_failures']['central_response_does_not_restore_shape_rank'])

    def test_exact_domains_refuse_nonpositive_and_incompatible_response_data(self):
        for H in (z.ZERO,z.K,[[1,1],[0,1]],[[1.0,0],[0,1]],[[True,0],[0,1]]):
            with self.assertRaises(ValueError):z.response_data(H,[z.K,z.L],z.I)
        for ds,B,scale in (([z.K],z.I,1),([z.K,z.R],z.I,1),([z.K,z.L],z.ZERO,1),
                           ([z.K,z.L],z.I,-1),([z.K,z.L],z.I,2)):
            with self.assertRaises(ValueError):z.response_data(z.I,ds,B,scale)
        for floor,T in ((0,1),(-1,8),(1,7),(0.5,8),(True,8)):
            with self.assertRaises(ValueError):z.curvature_budget(floor,T)

    def test_certificate_and_pin_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=Path(tmp)/'result';pin=Path(tmp)/'pin'
            changed=copy.deepcopy(self.cert);changed['claim_status']='PHYSICAL_GAP_SOLVED'
            result.write_text(json.dumps(changed));pin.write_text(z.chain.old.canonical_sha(changed))
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)
            result.write_text(json.dumps(self.cert));pin.write_text('0'*64)
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)

    def test_upstream_change_rejected_and_scope_retained(self):
        original=z.digest
        try:
            z.digest=lambda path:'0'*64
            with self.assertRaises(ValueError):z.source_checks()
        finally:z.digest=original
        self.assertIn('PHYSICAL_SELECTION_AND_ANISOTROPIC_JOINT_LIMIT_OPEN',self.cert['claim_status'])
        self.assertIn('not mechanically formalized',self.cert['evidence_scope']['general_proof'])


if __name__=='__main__':unittest.main()
