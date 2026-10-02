from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym48_reflected_time as y


class ReflectedTimeTests(unittest.TestCase):
    def test_fold_matches_literal_reflected_paths(self):
        out=y.history_controls()
        self.assertEqual(out['independent_reflected_path_identities'],48)
        self.assertEqual(out['positive_gram_inertia'],[4,0,0])
        self.assertEqual(out['mixed_source_gap_and_composition_checks'],30)

    def test_null_histories_translate_but_do_not_form_an_ideal(self):
        f=y.SIGNS[1]
        null=[(F(1),[(2,f)]),(F(-1,2),[(1,f)])]
        self.assertEqual(y.direct_pair(null,null),0)
        for t in (0,1,7):
            self.assertEqual(y.fold(null,t),[F(0)]*4)
        multiplied=y.product_history(null,[(F(1),[(2,f)])])
        self.assertEqual(y.direct_pair(multiplied,multiplied),F(9,16))

    def test_literal_time_translation_composes_on_mixed_history(self):
        history=[(F(2,3),[(1,y.SIGNS[1]),(3,y.SIGNS[2])]),
                 (F(-1,5),[(2,y.SIGNS[3])])]
        source=y.fold(history)
        self.assertEqual(y.fold(history,5),y.apply(y.transfer(5),source))
        self.assertEqual(y.fold(history,5),y.apply(y.transfer(2),y.fold(history,3)))
        self.assertEqual(y.direct_pair(history,history,2),y.dot(source,y.fold(history,2)))

    def test_resolvent_inverse_domain_and_graph_reconstruction(self):
        out=y.resolvent_controls()
        self.assertEqual(out['two_sided_inverse_checks'],8)
        self.assertEqual(out['resolvent_identities'],6)
        self.assertEqual(out['generator_domain_and_gap_checks'],24)
        H=y.generator();lam=F(11,7);x=list(map(F,(1,-3,2,5)))
        u=y.apply(y.resolvent(lam),x)
        self.assertEqual(y.apply(H,u),[a-lam*b for a,b in zip(x,u)])

    def test_common_domain_change_of_resolvent_parameter(self):
        lam,mu=F(3,5),F(11,4);x=list(map(F,(2,1,-4,3)))
        R,S=y.resolvent(lam),y.resolvent(mu)
        u=y.apply(S,x)
        recovered=y.apply(R,[a+(lam-mu)*b for a,b in zip(x,u)])
        self.assertEqual(recovered,u)

    def test_time_cut_and_mesh_are_separate_errors(self):
        out=y.quadrature_controls()
        self.assertEqual(len(out['outward_quadratures']),12)
        self.assertTrue(out['tail_and_mesh_both_needed'])
        # At fixed time cut, fine mesh converges to a truncated integral,
        # so its distance to the full resolvent need not decrease monotonically.
        first=out['outward_quadratures'][:3]
        self.assertGreater(F(first[-1]['independent_resolvent_error_upper']),
                           F(first[0]['independent_resolvent_error_upper']))
        coarse=y.quadrature_budget(F(1),F(4),F(1,8),F(1,4))
        fine=y.quadrature_budget(F(1),F(4),F(1,64),F(1,4))
        self.assertLess(fine,coarse)
        self.assertGreater(fine,y.exp_negative(F(4)).lo)

    def test_fundamental_casimir_and_coordinate_factor(self):
        x0=y.COORD[0]
        self.assertEqual(y.lap(x0),y.scale(x0,F(3,4)))
        self.assertEqual(y.integral(y.mul(x0,x0)),F(1,4))
        self.assertEqual(y.integral(y.grad_square(x0)),F(3,16))
        for a in (1,2,3):
            self.assertEqual(y.deriv(x0,a),y.scale(y.COORD[a],F(1,2)))

    def test_weighted_nonconstant_vacua_preserve_energy_identity(self):
        out=y.coefficient_controls()
        self.assertEqual(out['leibniz_and_fundamental_checks'],35)
        self.assertEqual(out['integration_by_parts_checks'],108)
        self.assertEqual(len(out['nonconstant_vacuum_energy_controls']),12)
        for row in out['nonconstant_vacuum_energy_controls']:
            self.assertGreaterEqual(F(row['weighted_energy']),0)
            self.assertLessEqual(F(row['weighted_energy']),F(row['coefficient_gradient_ceiling']))

    def test_sphere_moment_recursion_not_only_selected_values(self):
        for m in itertools_product_even(4):
            raised=[]
            for j in range(4):
                n=list(m);n[j]+=2;raised.append(tuple(n))
            self.assertEqual(sum((y.moment(n) for n in raised),F(0)),y.moment(m))
        self.assertEqual(y.moment((1,0,0,0)),0)
        self.assertEqual(y.moment((4,0,0,0)),F(1,8))

    def test_nonzero_excitation_has_separate_lower_bound(self):
        out=y.coefficient_controls()
        variance=F(out['single_site_variance']);energy=F(out['single_site_energy'])
        time=F(out['excited_history_time'])
        self.assertEqual(variance-2*time*energy,F(5,32))
        self.assertGreater(variance-2*time*energy,0)
        self.assertEqual(len(y.parameter_controls()),4)

    def test_high_contents_do_not_give_norm_continuity_or_bounded_generator(self):
        rows=y.high_content_controls()
        self.assertEqual([F(r['generator_norm_lower']) for r in rows],[F(272),F(4160),F(65792)])
        for row in rows:
            self.assertGreater(F(row['time_identity_error_lower']),F(1,2))

    def test_invalid_integral_and_source_parameters_are_refused(self):
        for args in ((F(0),F(2),F(1,8),F(1,4)),(F(1),F(2),F(1,8),F(0)),
                     (F(1),F(2),F(3),F(1,4))):
            with self.assertRaises(ValueError):
                y.quadrature_budget(*args)
        with self.assertRaises(ValueError):
            y.resolvent(F(0))
        with self.assertRaises(ValueError):
            y.transfer(-1)
        with self.assertRaises(ValueError):
            y.exp_negative(F(-1))
        self.assertTrue(all(y.negative_controls().values()))

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


def itertools_product_even(dimension):
    import itertools
    return itertools.product((0,2,4),repeat=dimension)


if __name__=='__main__':
    unittest.main()
