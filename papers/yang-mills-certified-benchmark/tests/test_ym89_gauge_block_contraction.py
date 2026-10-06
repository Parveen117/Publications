import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym89_gauge_block_contraction as y

p,a = y.p,y.algebra


def dec(x):
    x = F(x)
    return Decimal(x.numerator)/Decimal(x.denominator)


def apply(matrix,vector):
    return [sum(v*x for v,x in zip(row,vector)) for row in matrix]


def dot(x,z):
    return sum(v*w for v,w in zip(x,z))


def weighted_norm(matrix):
    return max(sum(2**abs(i-j)*abs(v) for j,v in enumerate(row)) for i,row in enumerate(matrix))


def comparison(n=5):
    return [[F(3) if i==j else F(-1,4) if abs(i-j)==1 else F(0)
             for j in range(n)] for i in range(n)]


class MarginalSchurTests(unittest.TestCase):
    def test_actual_endpoint_constants_and_free_limit(self):
        out = y.marginal_budget(y.previous.THETA_MAX)
        self.assertEqual(out['marginal_mixed_Hessian_row_upper'],F(13,960))
        self.assertEqual(out['marginal_weighted_mixed_Hessian_row_upper'],F(13,440))
        self.assertEqual(out['marginal_weighted_fill_row_upper'],F(1,110))
        self.assertEqual(out['fine_gap'],F(3,10))
        free = y.marginal_budget(0)
        self.assertEqual(free['marginal_mixed_Hessian_row_upper'],0)
        self.assertEqual(free['marginal_weighted_mixed_Hessian_row_upper'],0)

    def test_Gaussian_marginal_covariance_matches_full_inverse_submatrix(self):
        C = comparison()
        covariance = a.inverse(C)
        for size in (1,2,3,5):
            for retained in combinations(range(5),size):
                out = y.schur_comparison(C,retained)
                self.assertEqual(a.inverse(out['effective']),y.submatrix(covariance,retained,retained))

    def test_completing_the_square_identifies_the_conditional_mean_and_fill(self):
        C = comparison(4)
        R,Z = [0,3],[1,2]
        out = y.schur_comparison(C,R)
        x = [F(2,3),F(-5,7)]
        z_star = apply(a.scale(p.multiply(out['hidden_inverse'],y.submatrix(C,Z,R)),-1),x)
        displacement = [F(1,5),F(-2,9)]
        z = [v+w for v,w in zip(z_star,displacement)]
        full = [x[0],z[0],z[1],x[1]]
        expected = dot(x,apply(out['effective'],x))+dot(displacement,apply(y.submatrix(C,Z,Z),displacement))
        self.assertEqual(dot(full,apply(C,full)),expected)

    def test_every_elimination_order_gives_the_same_bound(self):
        C = comparison()
        direct = y.schur_comparison(C,[0,4])['effective']
        for order in permutations([1,2,3]):
            current,labels = C,list(range(5))
            for label in order:
                keep = [i for i,v in enumerate(labels) if v!=label]
                current = y.schur_comparison(current,keep)['effective']
                labels = [labels[i] for i in keep]
            self.assertEqual(current,direct)

    def test_hidden_fill_creates_an_off_diagonal_coupling(self):
        C = comparison(3)
        out = y.schur_comparison(C,[0,2])
        self.assertEqual(C[0][2],0)
        self.assertEqual(out['fill'][0][1],F(1,48))
        self.assertEqual(out['effective'][0][1],F(-1,48))
        self.assertNotEqual(out['effective'],y.submatrix(C,[0,2],[0,2]))

    def test_pins_and_marginalization_preserve_a_uniform_positive_margin(self):
        C = comparison(8)
        margin = min(sum(row) for row in C)
        for retained in ([0],[0,3,7],[1,2,3,4,5,6],list(range(8))):
            out = y.schur_comparison(C,retained)
            self.assertGreaterEqual(out['effective_row_margin'],margin)
            self.assertTrue(p.psd(p.add(out['effective'],a.scale(p.identity(len(retained)),margin),-1)))
            pinned = y.submatrix(out['effective'],[0],[0])
            self.assertGreaterEqual(pinned[0][0],margin)

    def test_hidden_inverse_series_has_a_full_weighted_tail(self):
        J = [[F(1,40) if abs(i-j)==1 else F(0) for j in range(4)] for i in range(4)]
        rho = F(1)
        q = weighted_norm(J)/rho
        inverse = a.inverse(p.add(a.scale(p.identity(4),rho),J,-1))
        term,partial = p.identity(4),a.zeros(4)
        for degree in range(9):
            partial = p.add(partial,a.scale(term,1/rho))
            error = p.add(inverse,partial,-1)
            self.assertTrue(all(v>=0 for row in error for v in row))
            self.assertLessEqual(weighted_norm(error),q**(degree+1)/(rho*(1-q)))
            term = p.multiply(term,a.scale(J,1/rho))
        endpoint = y.marginal_budget(y.previous.THETA_MAX)
        for degree in (0,1,4,12):
            self.assertEqual(y.hidden_walk_tail(y.previous.THETA_MAX,degree),F(40,11)*F(2,13)**(degree+1))
            self.assertEqual(y.hidden_walk_tail(0,degree),0)
        self.assertEqual(endpoint['weighted_J_row']/endpoint['rho'],F(2,13))

    def test_truncated_hidden_inverse_underestimates_the_memory_fill(self):
        C = comparison(4)
        R,Z = [0,3],[1,2]
        out = y.schur_comparison(C,R)
        truncated = a.scale(p.identity(2),F(1,3))
        rz = y.submatrix(C,R,Z)
        finite_fill = p.multiply(rz,p.multiply(truncated,p.transpose(rz)))
        self.assertEqual(finite_fill[0][1],0)
        self.assertGreater(out['fill'][0][1],0)

    def test_matrix_class_and_partition_errors_are_rejected(self):
        for C in ([[1,1],[1,1]],[[1,-2],[-2,1]],[[1,0],[-1,2]],[[0]]):
            with self.assertRaises(ValueError):
                y.schur_comparison(C,[0])
        for keep in ([],[0,0],[True],[7],[-1]):
            with self.assertRaises(ValueError):
                y.schur_comparison(comparison(),keep)


class NativeGaugeTests(unittest.TestCase):
    def test_product_pullback_cross_derivatives_in_noncommuting_configurations(self):
        for factor in (2,3):
            out = y.native_fiber_control(factor)
            self.assertEqual(out['exact_derivative_components_checked'],18)
            self.assertEqual(out['coarse_gradient'],out['gradient_reconstructed_from_cross_Hessian'])
            self.assertTrue(any(out['coarse_gradient']))
            self.assertNotEqual(out['same_path_cross_Hessian'],p.transpose(out['same_path_cross_Hessian']))

    def test_skew_operator_norm_constant_is_not_a_Hilbert_Schmidt_substitute(self):
        for factor in (2,3):
            out = y.native_fiber_control(factor)
            skew,g = out['antisymmetric_part'],out['coarse_gradient']
            expected = a.scale(p.add(a.scale(p.identity(3),dot(g,g)),a.outer(g),-1),F(1,4))
            self.assertEqual(p.multiply(p.transpose(skew),skew),expected)
            self.assertEqual(apply(skew,g),[0,0,0])
            self.assertEqual(sum(v*v for row in skew for v in row),dot(g,g)/2)

    def test_prefix_adjoint_rotations_are_exactly_orthogonal(self):
        q = (F(1,2),)*4
        units = [tuple(F(i==j) for i in range(3)) for j in range(3)]
        rotations = [y.adjoint(q,u) for u in units]
        self.assertEqual(p.multiply(rotations,p.transpose(rotations)),p.identity(3))
        for u in units:
            self.assertEqual(y.adjoint(y.local.y.qdagger(q),y.adjoint(q,u)),u)

    def test_internal_gauges_flatten_fibers_without_changing_holonomy(self):
        for factor in (2,3,5,9):
            out = y.gauge_flatten_control(factor)
            self.assertEqual(out['internal_links_flattened'],factor-1)
            self.assertTrue(out['endpoint_gauges_fixed'])
            self.assertEqual(sum(v*v for v in out['holonomy']),1)

    def test_native_retained_network_has_degree_two_internal_vertices(self):
        for factor in (2,3,4):
            out = y.network_control(factor)
            self.assertTrue(out['every_internal_retained_degree_two'])
            self.assertEqual(out['retained_links'],factor*out['coarse_links'])
            self.assertEqual(out['internal_vertices'],(factor-1)*out['coarse_links'])
            self.assertGreater(out['ball_path_incidence_checks'],0)

    def test_adding_an_internal_branch_breaks_the_required_geometry(self):
        cell = y.blocking.blocked_geometry(4,2)
        path = cell['paths'][0]
        vertex = cell['fine']['ends'][path[0]][1]
        retained = {e for path in cell['paths'] for e in path}
        touching = [e for e,ends in enumerate(cell['fine']['ends']) if vertex in ends]
        self.assertEqual(sum(e in retained for e in touching),2)
        extra = next(e for e in touching if e not in retained)
        retained.add(extra)
        self.assertEqual(sum(e in retained for e in touching),3)


class BlockContractionTests(unittest.TestCase):
    def test_optimized_path_weight_matches_independent_closed_forms(self):
        for b in range(2,257):
            expected = 2 if b==2 else 3*2**(b//2-1) if b%2==0 else 2**((b+1)//2)
            out = y.path_weight(b)
            self.assertEqual(out['weight'],expected)
            self.assertEqual(out['weight'],out['guaranteed_far_links']*2**out['radius'])

    def test_score_decays_by_one_half_per_added_fine_link(self):
        theta = y.previous.THETA_MAX
        for b in (2,3,4,8,16,32):
            cell = y.block_budget(theta,b)
            next_cell = y.block_budget(theta,b+1)
            self.assertEqual(cell['coarse_log_density_score_upper'],F(13,110*2**(b-1)))
            self.assertEqual(next_cell['coarse_log_density_score_upper'],cell['coarse_log_density_score_upper']/2)
            self.assertEqual(cell['coarse_log_density_score_upper'],2*cell['coarse_log_wave_gradient_upper'])

    def test_coarse_Ricci_floor_includes_all_off_diagonal_blocks(self):
        for theta in (0,F(1,8000),y.previous.THETA_MAX):
            for b in (2,3,4,8,16):
                cell = y.block_budget(theta,b)
                self.assertEqual(cell['unit_clock_coarse_gap_lower'],F(1,2)-2*cell['coarse_full_Hessian_row_upper'])
                self.assertGreaterEqual(cell['unit_clock_coarse_gap_lower'],F(227,480))
                self.assertEqual(cell['compressed_diagonal_gap_lower'],b*cell['unit_clock_coarse_gap_lower'])
        self.assertEqual(y.block_budget(y.previous.THETA_MAX,2)['unit_clock_coarse_gap_lower'],F(227,480))

    def test_large_block_row_bound_uses_spatial_decay_beyond_path_counting(self):
        theta = y.previous.THETA_MAX
        t = y.marginal_budget(theta)['marginal_mixed_Hessian_row_upper']
        for b in (8,16,32):
            cell = y.block_budget(theta,b)
            self.assertLess(cell['coarse_full_Hessian_row_upper'],t/(b-1))
        cell,next_cell = y.block_budget(theta,16),y.block_budget(theta,18)
        self.assertLess(next_cell['coarse_full_Hessian_row_upper'],cell['coarse_full_Hessian_row_upper']/2)

    def test_conditional_TV_counts_targets_and_native_diameter(self):
        theta = y.previous.THETA_MAX
        for m in (1,4,100):
            out = y.block_budget(theta,8,m)
            self.assertEqual(out['conditional_log_oscillation_upper_pi_multiple'],F(13*m,7040))
            self.assertEqual(out['conditional_TV_upper_pi_multiple'],F(13*m,28160))
            self.assertEqual(out['conditional_TV_rational_upper'],min(F(1),F(13*m,8960)))
        self.assertLess(y.local.spatial.pi_interval().hi,F(22,7))

    def test_density_oscillation_bound_is_sharp_for_two_level_densities(self):
        with localcontext() as context:
            context.prec = 70
            for ratio_root in (F(3,2),F(2),F(5),F(20)):
                s = ratio_root
                high_weight = 1/(s+1)
                density_low,density_high = 1/s,s
                self.assertEqual(high_weight*density_high+(1-high_weight)*density_low,1)
                TV = high_weight*(density_high-1)
                self.assertEqual(TV,(s-1)/(s+1))
                self.assertLessEqual(dec(TV),dec(s).ln()/2)

    def test_free_blocks_have_zero_score_but_do_not_claim_exact_gap_one_half(self):
        for b in (2,8,32):
            out = y.block_budget(0,b)
            self.assertEqual(out['coarse_log_density_score_upper'],0)
            self.assertEqual(out['conditional_TV_rational_upper'],0)
            self.assertEqual(out['unit_clock_coarse_gap_lower'],F(1,2))
            self.assertLess(out['unit_clock_coarse_gap_lower'],F(3,4))

    def test_coupling_and_block_scope_are_enforced(self):
        for theta in (-1,F(1,399),1,25):
            with self.assertRaises(ValueError):
                y.block_budget(theta,2)
        for b in (1,0,-2,True,F(3,2)):
            with self.assertRaises(ValueError):
                y.block_budget(0,b)
        with self.assertRaises(ValueError):
            y.block_budget(0,2,0)
        with self.assertRaises(ValueError):
            y.hidden_walk_tail(0,-1)


class DynamicMemoryTests(unittest.TestCase):
    def test_markov_operator_and_projection_are_genuinely_reversible(self):
        out = y.memory_control(7,F(1,11))
        H = out['positive_generator']
        P = [[F(1,2),0,F(1,2)],[0,1,0],[F(1,2),0,F(1,2)]]
        self.assertEqual(H,p.transpose(H))
        self.assertEqual(apply(H,[1,1,1]),[0,0,0])
        self.assertTrue(p.psd(H))
        self.assertEqual(p.multiply(P,P),P)
        self.assertEqual(P,p.transpose(P))
        self.assertEqual(apply(P,[1,1,1]),[1,1,1])
        self.assertEqual(apply(P,[1,0,-1]),[0,0,0])

    def test_diagonal_and_coupling_are_computed_from_the_original_three_state_chain(self):
        v,u = [1,-2,1],[1,0,-1]
        for fast,slow in ((F(2),F(1,4)),(F(1),F(1)),(F(64),F(1,4096))):
            out = y.memory_control(fast,slow)
            H = out['positive_generator']
            self.assertEqual(dot(v,apply(H,v))/dot(v,v),out['compressed_centered_diagonal'])
            self.assertEqual(dot(u,apply(H,u))/dot(u,u),out['hidden_floor'])
            self.assertEqual(dot(v,apply(H,u))**2/(dot(v,v)*dot(u,u)),out['coupling_squared'])

    def test_full_gap_vanishes_even_while_both_diagonal_floors_grow(self):
        previous = None
        for b in (2,4,16,64):
            out = y.memory_control(b,F(1,b*b))
            self.assertLessEqual(out['full_gap_upper'],F(3,b*b))
            self.assertLess(out['full_gap_upper'],out['compressed_centered_diagonal'])
            self.assertGreater(out['hidden_floor'],F(b,2))
            if previous:
                self.assertGreater(out['relative_memory'],previous['relative_memory'])
                self.assertLess(out['full_gap_upper'],previous['full_gap_upper'])
            previous = out

    def test_eigenvalue_bounds_against_independent_characteristic_polynomial(self):
        Pi = a.scale([[1]*3 for _ in range(3)],F(1,3))
        centered = p.add(p.identity(3),Pi,-1)
        with localcontext() as context:
            context.prec = 100
            for b in (2,4,16,64):
                out = y.memory_control(b,F(1,b*b))
                H = out['positive_generator']
                trace = sum(H[i][i] for i in range(3))
                minors = sum(H[i][i]*H[j][j]-H[i][j]*H[j][i] for i,j in combinations(range(3),2))
                self.assertEqual(trace,out['nonzero_eigenvalue_trace'])
                self.assertEqual(minors,out['nonzero_eigenvalue_product'])
                gap = (dec(trace)-(dec(trace)**2-4*dec(minors)).sqrt())/2
                self.assertLessEqual(dec(out['full_gap_lower']),gap)
                self.assertLessEqual(gap,dec(out['full_gap_upper']))
                self.assertTrue(p.psd(p.add(H,a.scale(centered,out['full_gap_lower']),-1)))

    def test_relative_memory_closes_the_two_block_bound_only_when_supplied(self):
        for b in (2,4,16):
            out = y.memory_control(b,F(1,b*b))
            lower = y.two_block_floor(out['compressed_centered_diagonal'],out['hidden_floor'],out['relative_memory'])
            self.assertEqual(lower,out['full_gap_lower'])
        H = [[F(4),F(3)],[F(3),F(9)]]
        lower = y.two_block_floor(4,9,F(1,4))
        self.assertTrue(p.psd(p.add(H,a.scale(p.identity(2),lower),-1)))
        for eta in (-1,1,2):
            with self.assertRaises(ValueError):
                y.two_block_floor(1,1,eta)


class EvidenceTests(unittest.TestCase):
    def test_sources_reproduce_and_dynamic_and_continuum_gates_stay_open(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']),401)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        self.assertTrue(cert['model']['actual_gauge_product_block_score_decays_exponentially'])
        self.assertTrue(cert['model']['coarse_conditionals_uniform_over_full_coarse_exterior'])
        for name in ('fine_arbitrary_pins_preserve_the_gauge_fiber_argument',
                     'static_Schur_matrix_is_dynamic_Schur_memory',
                     'coarse_diagonal_generator_is_exact_reduced_dynamics',
                     'actual_weak_coupling_block_gate_discharged',
                     'positive_four_dimensional_continuum_mass_gap',
                     'analytic_proofs_formally_machine_verified'):
            self.assertFalse(cert['model'][name])

    def test_result_pin_and_scope_drift_are_rejected(self):
        cert = y.run()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            pin.write_text(p.canonical_sha(cert))
            changed = copy.deepcopy(cert)
            changed['controls']['marginal_endpoint']['fine_gap'] = '30'
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
            with self.assertRaisesRegex(ValueError,'scope mismatch'):
                y.check(changed,result,pin)

    def test_upstream_hash_drift_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        key = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][key] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'sources'
            path.write_text(json.dumps(manifest))
            with patch.object(y,'SOURCES',path):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    y.source_checks()

    def test_cli_is_read_only_and_rejects_optimization(self):
        inputs = (y.RESULT,y.PIN,y.SOURCES)
        before = [path.read_bytes() for path in inputs]
        script = str(y.HERE/'ym89_gauge_block_contraction.py')
        for args in ([],['--check']):
            run = subprocess.run([sys.executable,'-B',script,*args],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            self.assertIn('YM89 PASS',run.stdout)
            self.assertEqual(before,[path.read_bytes() for path in inputs])
        run = subprocess.run([sys.executable,'-B','-O',script,'--check'],capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('without optimization',run.stderr)


if __name__=='__main__':
    unittest.main()
