"""Consume the original UGD numeral and EMK geometry modules unchanged."""
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, REPO/relative)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def run():
    if sys.version_info[:2] not in ((3, 11), (3, 12)) or not __debug__:
        raise RuntimeError('Use Python 3.11/3.12 without -O.')
    ugd = module('native_correction_ugd', 'papers/emk-ugd-algebra/certificates/ugd1_numerals.py')
    geometry = module('native_correction_geometry',
                      'papers/emk-recognition-geometry/certificates/ugdg1_phase_scale_seam_geometry.py')
    g3 = geometry.g3
    g1 = module('native_correction_metric',
                'papers/emk-recognition-geometry/certificates/emkg1_rotational_seam_metric.py')
    a, b = ugd.digit(3, 0, 1), ugd.digit(3, 0, 0)
    assert a.classical_projection() == b.classical_projection()
    assert a.phase_content() == b.phase_content() and a.total_seam()-b.total_seam() == 1
    carried = ugd.add_digits(a, a)
    broken = ugd.add_digits_broken_carry(a, a)
    assert carried.total_seam() == 2 and carried.ledger == 1
    assert broken.total_seam() == 1
    # Numeral multiplication is its own supplied law, not arrow composition.
    x = ugd.digit(1, 0, 0)
    product = ugd.mul_digits(x, x, seam_rule='additive')
    assert product.digits[0][0] == 1
    assert (1+1) % x.K == 2
    start = (Q(0), Q(0), 0)
    alpha = Q(g3.KMOD, 4)  # quarter-turn subgroup in the geometric phase chart
    orbit = [g3.rho_power(start, n, alpha=alpha, beta=Q(0), q=1) for n in range(5)]
    assert orbit[4] == (0, 0, 4) and orbit[4] != start
    assert geometry.phase_closes(orbit[4][0]) and not geometry.sheet_closes(orbit[4][2])
    # EMK base (u,v) and the UGD fibre (phi,sigma,k) are distinct carriers.
    metric_rows = []
    for kappa in [Q(0), Q(3)]:
        w = g1.W_quadratic(Q(0), kappa)
        length_squared = g3.L**2*w
        assert w == 1 and g1.gaussian_curvature_closed(Q(0), kappa) == -kappa
        metric_rows.append(dict(kappa=str(kappa), seam_length_squared=str(length_squared),
                                curvature=str(-kappa)))
    assert metric_rows[0]['seam_length_squared'] == metric_rows[1]['seam_length_squared']
    # Clock-free pullback differential, using the actual EMK metric evaluator.
    kappa, v, beta, gamma = Q(3), Q(1, 4), Q(1, 3), Q(-2, 5)
    delta = g1.W_quadratic(v+beta, kappa)-g1.W_quadratic(v, kappa)
    total = g1.W_quadratic(v+beta+gamma, kappa)-g1.W_quadratic(v, kappa)
    later = g1.W_quadratic(v+beta+gamma, kappa)-g1.W_quadratic(v+beta, kappa)
    assert delta == kappa*(2*v*beta+beta**2) and total == delta+later
    clock_periods = [Q(1), Q(2)]
    protocol_speed_squared = [g3.L**2/time**2 for time in clock_periods]
    assert protocol_speed_squared[0] == 4*protocol_speed_squared[1]
    return dict(
        numeral=dict(classical_projection=str(a.classical_projection()),
                     equal_phase_content=True, distinct_seam_charges=[a.total_seam(), b.total_seam()],
                     correct_carry_charge=carried.total_seam(), dropped_carry_charge=broken.total_seam(),
                     numeral_product_exponent=product.digits[0][0], arrow_sum_exponent=2,
                     numeral_multiplication_is_not_transport_composition=True),
        helix=dict(phase_period=g3.KMOD, phase_step=str(alpha),
                   orbit=[[str(phi), str(sigma), k] for phi, sigma, k in orbit],
                   four_step_visible_return=True, four_step_full_return=False),
        metric=dict(base_coordinates=['u', 'v'], fibre_coordinates=['phi', 'sigma', 'integer k'],
                    rows=metric_rows, transition_differential=str(delta),
                    composed_differential=str(total), clock_free_chain_rule=True,
                    admitted_clock_periods=[str(x) for x in clock_periods],
                    protocol_speed_squared=[str(x) for x in protocol_speed_squared],
                    physical_clock_selected=False))


if __name__ == '__main__':
    print(json.dumps(run(), sort_keys=True))
