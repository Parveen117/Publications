"""Exact conditional source-to-phase adapter controls; no empirical verdict."""
from fractions import Fraction as Q
import json


def kernel_moment(T, n):
    """Integral of triangular acceleration weight times t**n, exactly."""
    T = Q(T)
    if T <= 0 or type(n) is not int or n < 0:
        raise ValueError('Positive pulse separation and nonnegative integer degree required')
    return ((2*T)**(n+2)-2*T**(n+2))/((n+1)*(n+2))


def phase(green, source, integrated_detector_row, scale=1):
    """scale * r * H^-1 * b with a supplied Green matrix and sensor row.

    r includes time sensitivity, source-independent trajectory sampling,
    and the calibrated field-gradient observation. It is not derived here.
    """
    n = len(source)
    if (not n or len(green) != n or any(len(row) != n for row in green)
            or len(integrated_detector_row) != n):
        raise ValueError('Incompatible Green/source/detector dimensions')
    field = [sum(Q(green[i][j])*Q(source[j]) for j in range(n)) for i in range(n)]
    return Q(scale)*sum(Q(r)*v for r, v in zip(integrated_detector_row, field))


def run():
    T = Q(4, 25)
    assert kernel_moment(T, 0) == T*T
    assert kernel_moment(T, 1) == T**3
    assert kernel_moment(T, 2) == Q(7, 6)*T**4
    green = ((Q(2, 3), Q(1, 3)), (Q(1, 3), Q(2, 3)))
    bC, bF, r = (1, 0), (0, 1), (1, -1)
    delta = phase(green, bC, r)-phase(green, bF, r)
    assert delta == Q(2, 3)
    assert delta == phase(green, (1, -1), r)
    assert phase(green, (1, 1), r) == 0  # common response cancels
    assert phase(green, (4, -4), r) == 4*delta
    # Changed paths must not be silently represented by one fixed row.
    changed_row_delta = phase(green, bC, r)-phase(green, bF, (2, -2))
    assert changed_row_delta == 1 and changed_row_delta != delta
    # Dimensions and illegal pulse spacing must fail loudly.
    for call in (lambda: phase(green, bC, (1,)), lambda: kernel_moment(0, 0)):
        try:
            call()
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid adapter input accepted')
    return {'triangular_kernel_moments_exact': True,
            'same_row_source_double_difference': str(delta),
            'changed_path_double_difference': str(changed_row_delta),
            'common_response_cancels': True,
            'native_force_to_phase_is_conditional': True,
            'experimental_native_prediction_computed': False}


if __name__ == '__main__':
    print(json.dumps(run(), sort_keys=True))
