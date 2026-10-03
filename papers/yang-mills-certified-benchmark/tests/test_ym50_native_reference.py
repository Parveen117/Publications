from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym50_native_reference as y


class NativeReferenceTests(unittest.TestCase):
    def test_native_quaternion_product_and_infinite_order_controls(self):
        result=y.quaternion_controls()
        self.assertEqual(result['native_matrix_product_checks'],64)
        self.assertEqual(result['infinite_order_congruence_controls'],16)
        self.assertNotEqual(y.qmul(y.turns()[0],y.turns()[2]),
                            y.qmul(y.turns()[2],y.turns()[0]))

    def test_finite_energy_fixed_spaces_and_poisson_invertibility(self):
        rows=y.finite_energy_controls()['finite_degree_spaces']
        self.assertEqual([r['dimension'] for r in rows],[1,4,10,20,35])
        self.assertEqual([r['stationary_dimension'] for r in rows],[1,0,1,0,1])

    def test_raw_word_counts_and_padded_native_trace_agree(self):
        result=y.counted_readout_controls()
        self.assertEqual(result['independent_word_transfer_checks'],24)
        self.assertEqual(result['padded_native_trace_checks'],18)
        self.assertEqual(result['complex_native_pairing_checks'],2)
        self.assertEqual(result['largest_raw_word_count'],1728)

    def test_finite_pairing_positivity_from_independent_records(self):
        x=y.y48.COORD
        p=y.y48.add(x[0],y.y48.scale(x[1],F(-7,3)),y.y48.mul(x[2],x[3]))
        for k in (0,1,2,3):
            mean=y.direct_count(p,k)
            square=y.direct_count(y.y48.mul(p,p),k)
            self.assertGreaterEqual(square,mean*mean)
            self.assertEqual(y.direct_count(y.y48.ONE,k),1)

    def test_constructive_telescoping_and_explicit_tails(self):
        result=y.tail_controls()
        self.assertEqual(result['exact_telescoping_cases'],20)
        self.assertEqual(result['pointwise_tail_checks'],160)
        for row in result['cases']:
            self.assertLessEqual(F(row['identity_readout_error']),
                                 F(row['uniform_error_upper']))

    def test_poisson_solution_has_expected_linear_and_radial_parts(self):
        x=y.y48.COORD[0]
        u,bound=y.poisson(y.y48.add(x,y.y48.scale(y.RADIUS,F(9,7))))
        self.assertEqual(u,y.y48.scale(x,5))
        self.assertEqual(bound,5)
        for N in (1,7,25):
            self.assertEqual(y.evaluate(y.cesaro(x,N),y.ONE),
                             5*(1-F(4,5)**N)/N)

    def test_projected_moments_match_existing_reference(self):
        result=y.moment_and_bridge_controls()
        self.assertEqual(result['independent_projected_moments'],70)
        self.assertEqual(result['legacy_moment_matches'],495)
        self.assertEqual(result['left_right_invariance_checks'],140)
        self.assertEqual(result['recognition_gram_inertia'],[30,0,5])
        self.assertEqual(result['bounded_multiplier_checks'],36)
        self.assertEqual(y.native_moment((2,2,0,0)),F(1,24))
        self.assertEqual(y.native_moment((4,0,0,0)),F(1,8))

    def test_unit_relation_and_observables_descend_through_null_space(self):
        zero=y.y48.add(y.RADIUS,y.y48.scale(y.y48.ONE,-1))
        factor=y.y48.add(y.y48.ONE,y.y48.scale(y.y48.COORD[1],F(5,4)))
        product=y.y48.mul(zero,factor)
        self.assertEqual(y.phi(y.y48.mul(product,product)),0)
        for N in (1,3):
            self.assertEqual(y.evaluate(y.cesaro(product,N),y.ONE),0)

    def test_raw_record_non_cauchy_energy_from_direct_extensions(self):
        x=y.y48.COORD[0]
        for n in (0,1,2):
            energy=sum(((y.evaluate(x,y.qmul(q,p))-y.evaluate(x,p))**2
                         for p in y.endpoints(n) for q in y.ALPHABET),F(0))/12**(n+1)
            self.assertEqual(energy,F(1,10)-F(1,50)*F(43,75)**n)
            self.assertGreaterEqual(energy,F(2,25))

    def test_single_axis_fails_the_full_reference(self):
        x=y.y48.COORD
        p=y.y48.add(y.y48.mul(x[2],x[2]),y.y48.mul(x[3],x[3]))
        for q in y.turns((1,)):
            self.assertEqual(y.substitute(p,y.qmatrix(q)),p)
        self.assertEqual(y.evaluate(p,y.ONE),0)
        self.assertEqual(y.phi(p),F(1,2))

    def test_native_heat_generator_and_outward_refinement_controls(self):
        result=y.heat_controls()
        self.assertEqual(result['skew_generator_checks'],15)
        self.assertEqual(result['native_casimir_polynomial_checks'],70)
        self.assertEqual(result['harmonic_eigenvalue_checks'],4)
        self.assertEqual(len(result['outward_heat_refinements']),18)
        for row in result['outward_heat_refinements']:
            self.assertLessEqual(F(row['independent_scalar_error_upper']),
                                 F(row['general_coefficient_bound']))

    def test_heat_respects_sphere_relation_and_has_correct_clock(self):
        x=y.y48.COORD[0]
        relation=y.y48.add(y.RADIUS,y.y48.scale(y.y48.ONE,-1))
        self.assertEqual(y.y48.sphere_reduce(y.lap(y.y48.mul(relation,x))),{})
        self.assertEqual(y.lap(x),y.y48.scale(x,F(3,4)))
        one_axis=y.y48.scale(y.native_generator(y.native_generator(x,1),1),-1)
        self.assertEqual(one_axis,y.y48.scale(x,F(1,4)))
        self.assertEqual(y.heat_budget(2,F(1),32),F(3,16))

    def test_refusals_cover_memory_normalization_and_clock(self):
        result=y.refusal_controls()
        self.assertEqual(len(result['controls']),12)
        self.assertTrue(all(result['controls'].values()))
        self.assertEqual(result['raw_record_energy_limit'],'1/10')

    def test_invalid_count_degree_and_heat_inputs_are_refused(self):
        calls=(lambda:y.basis(-1),lambda:y.vector(y.y48.ONE,1),
               lambda:y.endpoints(-1),lambda:y.cesaro(y.y48.ONE,0),
               lambda:y.heat_budget(1,F(-1),4),lambda:y.heat_budget(1,F(1),0),
               lambda:y.cosine_square(F(-1)),lambda:y.cosine_square(F(2)))
        for call in calls:
            with self.assertRaises(ValueError):
                call()

    def test_origin_ledger_separates_readout_from_raw_memory_and_physics(self):
        ledger=json.loads(y.ORIGIN.read_text())
        self.assertEqual(ledger['evidence']['runtime'],'Python 3.12 only')
        self.assertIn('not recognition-Cauchy',
                      ledger['refused_identifications']['raw_record_completion'])
        self.assertIn('specified',ledger['refused_identifications']['physical_uniqueness'])
        self.assertEqual(len(ledger['open_obligations']),7)

    def test_source_tampering_is_refused(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()

    def test_record_and_pin_tampering_are_refused(self):
        data={'verdict':'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'result.json',Path(folder)/'pin'
            record.write_text(json.dumps(data));pin.write_text(y.canonical_sha(data))
            y.check(data,record,pin)
            record.write_text(json.dumps({'verdict':'FAIL'}))
            with self.assertRaises(ValueError):
                y.check(data,record,pin)
            record.write_text(json.dumps(data));pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(data,record,pin)


if __name__=='__main__':
    unittest.main()
