"""Retrospective external apparatus control; Newtonian gravity is explicit.

Requires NumPy. No phase fitting, no native-gravity verdict. Reference inputs
and observed summary are separated: forward() never reads the observation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent


def quad(lo, hi, n):
    x, w = leggauss(n)
    return lo + (x + 1) * (hi - lo) / 2, w * (hi - lo) / 2


def axial_field_control(order=28):
    """Independent closed-form cylinder field checks sign and normalization."""
    radius, length, height = 0.05, 0.15, 0.4
    r, wr = quad(0, radius, order)
    z, wz = quad(-length/2, length/2, order)
    dz = z[None, :]-height
    numerical = np.sum((2/ (radius**2*length)) * wr[:, None]*wz[None, :]
                       *r[:, None]*dz/(r[:, None]**2+dz**2)**1.5)
    # Unit G and total mass; exterior point on the cylinder's own axis.
    exact = -2/(radius**2*length) * (
        length + np.sqrt((height-length/2)**2+radius**2)
        - np.sqrt((height+length/2)**2+radius**2))
    if exact >= 0 or abs(float(numerical)-float(exact)) > 1e-10:
        raise AssertionError('Analytic axial-cylinder field control failed')


def forward(inputs, order):
    p = inputs['reported_nominal_geometry']
    q = inputs['external_reference_choices']
    radius, length = p['diameter_m']/2, p['height_m']
    time = p['pulse_separation_s']
    rho, wr = quad(0, radius, order)
    theta, wt = quad(0, 2*np.pi, 2*order)
    z, wz = quad(-length/2, length/2, order)
    r, th, z = rho[:, None, None], theta[None, :, None], z[None, None, :]
    weight = (wr[:, None, None]*wt[None, :, None]*wz[None, None, :]*r
              / (np.pi*radius**2*length))
    if abs(float(weight.sum())-1) > 1e-12:
        raise AssertionError('Cylinder mass normalization failed')
    t1, w1 = quad(0, time, order)
    t2, w2 = quad(time, 2*time, order)
    times = np.r_[t1, t2]
    sensitivity = np.r_[w1*t1, w2*(2*time-t2)]
    if abs(float(sensitivity.sum())-time**2) > 1e-14:
        raise AssertionError('Constant-acceleration response failed')
    k_eff = 4*np.pi/q['nominal_laser_wavelength_m']
    phase = {}
    for name in ('C', 'F'):
        centers = p[f'configuration_{name}_centers_relative_to_lower_apex_m']
        phase[name] = []
        for apex in (0, p['cloud_apex_separation_m']):
            path = apex - q['background_acceleration_m_s2']/2 * (
                times-time-p['apex_after_middle_pulse_s'])**2
            acceleration = np.zeros_like(times)
            for center in centers:
                # Six cylinders on each of two radii, on each platform.
                for offset in (2*radius, 2*radius*np.sqrt(3)):
                    xy2 = offset**2+r*r+2*offset*r*np.cos(th)
                    for j, point in enumerate(path):
                        dz = center+z-point
                        acceleration[j] += (q['G_m3_kg_s2'] *
                            p['total_cylinder_mass_kg']/p['cylinder_count'] * 6 *
                            np.sum(weight*dz/(xy2+dz*dz)**1.5))
            phase[name].append(float(k_eff*np.dot(sensitivity, acceleration)))
    # Signed convention fixed explicitly; no absolute-value matching.
    phase['CF_lower_minus_upper_rad'] = ((phase['C'][0]-phase['C'][1]) -
                                        (phase['F'][0]-phase['F'][1]))
    return phase


def build(inputs):
    axial_field_control()
    orders = inputs['numerics']['quadrature_orders']
    phases = [forward(inputs, n) for n in orders]
    values = [p['CF_lower_minus_upper_rad'] for p in phases]
    change = abs(values[-1]-values[-2])
    if change >= inputs['numerics']['successive_phase_tolerance_rad']:
        raise AssertionError('Quadrature stability check failed')
    p = inputs['published_summary_only']
    missing = [k for k, v in inputs['native_validation_readiness'].items() if not v]
    return {
        'classification': 'EXTERNAL_REFERENCE_CONTROL_ONLY',
        'native_validation': 'NOT_EVALUATED',
        'analytic_axial_cylinder_control': 'PASS_NUMERICAL_CONTROL',
        'missing_native_validation_inputs': missing,
        'quadrature_orders': orders,
        'CF_phase_by_order_rad': [round(v, 10) for v in values],
        'last_successive_difference_below_1e_minus_8_rad': True,
        'phase_components_rad': {k: [round(v, 10) for v in phases[-1][k]] for k in ('C', 'F')},
        'observed_aggregate_phase_rad': p['phase_CF_rad'],
        'unmodelled_observed_minus_cylinders_rad': round(p['phase_CF_rad']-values[-1], 10),
        'relative_shortfall_percent': round(100*(1-values[-1]/p['phase_CF_rad']), 6),
        'authors_full_simulation_rad': p['authors_full_simulation_rad'],
        'comparison_is_held_out': False,
        'native_field_used_in_reference': False,
        'scientific_z_score_or_pass_computed': False
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    inputs = json.loads((HERE/'atom_interferometer_2014_inputs.json').read_text())
    result = build(inputs)
    target = HERE/'ATOM_INTERFEROMETER_REFERENCE.json'
    if args.write:
        target.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    elif args.check:
        if json.loads(target.read_text()) != result:
            raise AssertionError('Reference result drift')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
