from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1]/'certificates/nsr2_diagram_properties.py'
SPEC = importlib.util.spec_from_file_location('nsr2_under_test', SOURCE)
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


class DiagramPropertyTests(unittest.TestCase):
    def test_jacobian_derivative_respects_the_held_variable(self):
        x, y = m.previous.variables(2, 3)
        self.assertEqual(m.constrained(x*x+y, x, y), 4)
        self.assertEqual(m.constrained(x*x+y, x, x+y), 3)
        with self.assertRaises(ValueError):
            m.constrained(x*x+y, x, 2*x)

    def test_four_potential_differentials(self):
        u = m.polynomial_fixture(F(2), F(3), 7)
        d = m.chart(u, 2, 3)
        p, q = d['p'].derivative(), d['q'].derivative()
        for axis in (0,1):
            dx, dy = F(axis==0), F(axis==1)
            dp, dq = d['p'].derivative(axis), d['q'].derivative(axis)
            self.assertEqual(d['Fx'].derivative(axis), -2*dp+q*dy)
            self.assertEqual(d['Fy'].derivative(axis), p*dx-3*dq)
            self.assertEqual(d['G'].derivative(axis), -2*dp-3*dq)

    def test_nonlinear_derivative_definitions_match_all_six_formulas(self):
        for seed in range(12):
            x, y = F(seed+1, 3), F(seed+2, 5)
            u = m.polynomial_fixture(x,y,seed)
            self.assertEqual(m.derivative_readouts(u,x,y),m.formula_readouts(u,x,y))

    def test_original_pressure_orientation_is_recovered(self):
        s,v = m.previous.variables(1,1)
        ds,dv = s-1,v-1
        u = 10*ds-10*dv+ds*ds+ds*dv+F(5,2)*dv*dv
        thermo = m.thermal_readouts(u,1,1)
        generic = m.derivative_readouts(u,1,1)
        self.assertEqual(thermo['lambda_P'],F(-9,5))
        self.assertEqual(thermo['z_S'],2)
        self.assertEqual(thermo['lambda_T'],F(20,9))
        self.assertEqual(thermo['z_S'],-generic['z_x'])
        self.assertEqual(thermo['lambda_T'],-generic['lambda_p'])

    def test_Born_Infeld_six_property_values(self):
        r=m.derivative_readouts(m.previous.born_infeld(F(1,2),F(1,2)),F(1,2),F(1,2))
        self.assertEqual(tuple(r[k] for k in m.NAMES),
                         (F(-3,8),F(-2,5),F(-2,5),F(-3,8),F(5,8),F(2,3)))

    def test_electromagnetic_susceptibility_dictionary(self):
        u=m.previous.born_infeld(F(1,2),F(1,2))
        r=m.derivative_readouts(u,F(1,2),F(1,2))
        h=m.previous.metric_jets(u)[0]
        chi=m.previous.responses(h)
        self.assertEqual(r['lambda_q'],-F(1,2)/chi['epsilon_h'])
        self.assertEqual(r['lambda_y'],-F(1,2)/chi['epsilon_y'])
        self.assertEqual(r['lambda_x'],-F(1,2)/chi['mu_x'])
        self.assertEqual(r['z_p'],-F(1,2)/chi['mu_e'])
        self.assertEqual(r['z_x'],u.derivative(1)*chi['mu_x'])
        self.assertEqual(r['lambda_p'],u.derivative(1)*chi['mu_e'])

    def test_closure_and_conjugate_products(self):
        data=m.check_identities(m.previous.born_infeld(F(1,2),F(1,2)),F(1,2),F(1,2))
        r=data['values']
        self.assertEqual(r['lambda_x']*r['z_x'],F(-1,4))
        self.assertEqual(r['lambda_p']*r['z_p'],F(-1,4))
        self.assertEqual(data['gamma_1_extension'],F(15,16))
        self.assertEqual(data['gamma_2_extension'],F(16,15))

    def test_zero_coordinate_keeps_cross_product_but_not_zero_over_zero_ratio(self):
        data=m.check_identities(m.polynomial_fixture(0,1,2),F(0),F(1))
        self.assertFalse(data['literal_lambda_ratios_defined'])
        r=data['values']
        self.assertEqual(r['lambda_q'],0)
        self.assertEqual(r['lambda_y'],0)
        self.assertEqual(r['lambda_q']*r['lambda_x'],r['lambda_y']*r['z_p'])

    def test_linear_field_model_reduces_to_fields(self):
        x,y=m.previous.variables(2,3)
        u=x*x/10+y*y/14
        r=m.derivative_readouts(u,2,3)
        self.assertEqual(tuple(r[k] for k in m.NAMES),
                         (F(-2,5),F(-2,5),F(-3,7),F(-3,7),F(3),F(3)))

    def test_mixed_channel_separates_the_three_pairs(self):
        x,y=m.previous.variables(2,3)
        u=x*x/2+x*y/4+y*y/2
        r=m.check_identities(u,F(2),F(3))['values']
        self.assertNotEqual(r['lambda_q'],r['lambda_y'])
        self.assertNotEqual(r['lambda_x'],r['z_p'])
        self.assertNotEqual(r['z_x'],r['lambda_p'])

    def test_separable_nonlinearity_does_not_require_pair_splitting(self):
        x,y=m.previous.variables(2,3)
        u=x*x/2+x**4/4+y*y/2+y**4/4
        r=m.check_identities(u,F(2),F(3))['values']
        self.assertEqual(r['lambda_q'],r['lambda_y'])
        self.assertEqual(r['lambda_x'],r['z_p'])
        self.assertEqual(r['z_x'],r['lambda_p'])

    def test_same_six_values_can_hide_opposite_signed_coupling(self):
        minus,plus=m.signed_pair(-1),m.signed_pair(1)
        self.assertEqual(m.derivative_readouts(minus,1,1),m.derivative_readouts(plus,1,1))
        self.assertEqual(minus.derivative(0,1),-1)
        self.assertEqual(plus.derivative(0,1),1)

    def test_next_response_derivative_separates_the_signed_pair(self):
        for sign in (-1,1):
            r=m.formula_jets(m.signed_pair(sign),1,1)
            self.assertEqual(r['z_x'].derivative(0),F(sign,5))
            self.assertEqual(r['lambda_p'].derivative(0),F(2*sign,9))

    def test_reference_origin_is_not_silently_erased(self):
        u=m.signed_pair(1)
        physical=m.derivative_readouts(u,1,1)
        redefined=m.derivative_readouts(u,0,0)
        self.assertTrue(all(physical[k]!=0 for k in m.NAMES[:4]))
        self.assertTrue(all(redefined[k]==0 for k in m.NAMES[:4]))

    def test_constant_energy_shift_changes_no_property(self):
        u=m.signed_pair(1)
        self.assertEqual(m.derivative_readouts(u,1,1),m.derivative_readouts(u+37,1,1))

    def test_singular_or_indefinite_response_is_rejected(self):
        x,y=m.previous.variables(1,1)
        for u in (x*x/2,(x*x-y*y)/2):
            with self.assertRaises(ValueError):
                m.derivative_readouts(u,1,1)


if __name__=='__main__':
    unittest.main()
