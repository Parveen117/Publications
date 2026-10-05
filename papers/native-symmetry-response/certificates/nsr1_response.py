"""NSR1 exact finite controls. Python 3.12, standard library only.

Taylor coefficients are rational. The inherited independent Christoffel
calculation is compared with the represented native response commutator.
This is an evidence adapter, not a new native algebra implementation.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import factorial, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKET = HERE.parent
SPEC = importlib.util.spec_from_file_location(
    'nsr1_nt_controls', ROOT/'papers/native-thermodynamic-curvature/controls.py')
nt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(nt)


def require_runtime():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError('Use Python 3.12 only')
    if sys.flags.optimize:
        raise RuntimeError('Run without optimization: evidence assertions must execute')


def rational_sqrt(value):
    value = F(value)
    if value < 0:
        raise ValueError('Negative square-root argument')
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    if n*n != value.numerator or d*d != value.denominator:
        raise ValueError('This finite fixture requires a rational square root')
    return F(n, d)


class Jet:
    """Two-variable Taylor polynomial through total degree four.

    Coefficients multiply dx**i dy**j, including the Taylor factorials.
    Operations take place in the exact quotient by degree-five monomials.
    """
    degree = 4

    def __init__(self, value=0):
        if isinstance(value, Jet):
            value = value.c
        if not isinstance(value, dict):
            value = {(0, 0): F(value)}
        self.c = {k: F(v) for k, v in value.items()
                  if v and len(k) == 2 and min(k) >= 0 and sum(k) <= self.degree}

    @classmethod
    def variable(cls, value, axis):
        if axis not in (0, 1):
            raise ValueError('Two axes only')
        return cls({(0, 0): F(value), (1-axis, axis): F(1)})

    def __add__(self, other):
        other = Jet(other)
        return Jet({k: self.c.get(k, F(0))+other.c.get(k, F(0))
                    for k in self.c.keys() | other.c.keys()})

    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()})

    def __sub__(self, other):
        return self + (-Jet(other))

    def __rsub__(self, other):
        return Jet(other) + (-self)

    def __mul__(self, other):
        other = Jet(other)
        out = {}
        for (i, j), a in self.c.items():
            for (k, l), b in other.c.items():
                key = (i+k, j+l)
                if sum(key) <= self.degree:
                    out[key] = out.get(key, F(0)) + a*b
        return Jet(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, Jet):
            return self * other.power(F(-1))
        return self * (1/F(other))

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            return self.power(F(exponent))
        answer = Jet(1)
        for _ in range(exponent):
            answer = answer*self
        return answer

    def power(self, exponent):
        exponent = F(exponent)
        base = self.c.get((0, 0), F(0))
        if exponent == F(1, 2):
            scale = rational_sqrt(base)
        elif exponent.denominator == 1:
            scale = base**int(exponent)
        else:
            raise ValueError('Only square roots and integer powers are used')
        if not base:
            raise ValueError('Nonzero expansion value required')
        remainder = self/base-1
        total, term, coefficient = Jet(1), Jet(1), F(1)
        for k in range(1, self.degree+1):
            term = term*remainder
            coefficient *= (exponent-(k-1))/k
            total += coefficient*term
        return scale*total

    def derivative(self, *axes):
        if any(axis not in (0, 1) for axis in axes):
            raise ValueError('Two axes only')
        i, j = axes.count(0), axes.count(1)
        if i+j > self.degree:
            raise ValueError('Derivative exceeds the stored Taylor order')
        return self.c.get((i, j), F(0))*factorial(i)*factorial(j)


def variables(x, y):
    return Jet.variable(F(x), 0), Jet.variable(F(y), 1)


def born_infeld(x, y):
    x, y = variables(x, y)
    return ((1+x*x)*(1+y*y)).power(F(1, 2))-1


def radial(x, y, beta=1):
    x, y = variables(x, y)
    radius2 = x*x+y*y
    return radius2/2+F(beta)*radius2*radius2/4


def relaxed(x, y):
    x, y = variables(x, y)
    return (x*x+y*y-x*x*y*y)/2


def clamped(x, y, z):
    x, y = variables(x, y)
    return (x*x+y*y+F(z)**2)/2-F(z)*x*y


def quadratic_tangent(jet):
    return Jet({k: v for k, v in jet.c.items() if sum(k) <= 2})


def metric_jets(energy):
    h = [[energy.derivative(i, j) for j in range(2)] for i in range(2)]
    dh = [[[energy.derivative(p, i, j) for j in range(2)]
           for i in range(2)] for p in range(2)]
    ddh = [[[[energy.derivative(p, q, i, j) for j in range(2)]
             for i in range(2)] for q in range(2)] for p in range(2)]
    return h, dh, ddh


def check_positive(h):
    if h != nt.tr(h) or h[0][0] <= 0 or nt.det(h) <= 0:
        raise ValueError('A symmetric strictly positive response matrix is required')


def shape_kappa(h, dh):
    """Rational form of -f_xy/sqrt(det H), avoiding unnecessary radicals."""
    check_positive(h)
    a, b, c = h[0][0], h[0][1], h[1][1]
    rho, log_z = [], []
    for d in dh:
        da, db, dc = d[0][0], d[0][1], d[1][1]
        rho.append(db-b*(da/a+dc/c)/2)  # sqrt(ac) times derivative of r
        log_z.append((da/a-dc/c)/2)
    return a*c*(rho[0]*log_z[1]-rho[1]*log_z[0])/(2*nt.det(h)**2)


def analyze(energy):
    h, dh, ddh = metric_jets(energy)
    check_positive(h)
    hi = nt.inv(h)
    x = [nt.mul(hi, d) for d in dh]
    native_f = nt.scale(nt.sub(nt.mul(x[0], x[1]), nt.mul(x[1], x[0])), F(-1, 4))
    gamma, _, coordinate_f = nt.christoffel(h, dh, ddh)
    if native_f != coordinate_f:
        raise ArithmeticError('Native and independent Christoffel routes disagree')
    kappa = nt.mul(h, native_f)[0][1]/nt.det(h)
    determinant_kappa = -nt.det([h[0][:]+[h[1][1]],
                                dh[0][0][:]+[dh[0][1][1]],
                                dh[1][0][:]+[dh[1][1][1]]])/(4*nt.det(h)**2)
    if not kappa == determinant_kappa == shape_kappa(h, dh):
        raise ArithmeticError('Three curvature readouts disagree')
    for i in range(2):
        if nt.add(nt.mul(nt.tr(gamma[i]), h), nt.mul(h, gamma[i])) != dh[i]:
            raise ArithmeticError('Metric compatibility failed')
    return dict(H=h, dH=dh, ddH=ddh, delta=nt.det(h), F=native_f, kappa=kappa)


def responses(h):
    check_positive(h)
    a, b, c, delta = h[0][0], h[0][1], h[1][1], nt.det(h)
    return dict(epsilon_y=1/a, epsilon_h=c/delta, mu_x=1/c, mu_e=a/delta,
                mixed=-b/delta, gain=a*c/delta)


def wave_data(h):
    check_positive(h)
    a, b, c = h[0][0], h[0][1], h[1][1]
    s, impedance = rational_sqrt(a*c), rational_sqrt(a/c)
    return dict(s=s, r=b/s, Z=impedance, v_plus=b+s, v_minus=b-s)


def reconstruct(v_plus, v_minus, impedance):
    v_plus, v_minus, impedance = F(v_plus), F(v_minus), F(impedance)
    if not v_plus > 0 > v_minus or impedance <= 0:
        raise ValueError('Require the stated positive-response wave sector')
    s, b = (v_plus-v_minus)/2, (v_plus+v_minus)/2
    return nt.mat([[s*impedance, b], [b, s/impedance]])


def schur_one(hqq, hqz, hzz):
    hzz = F(hzz)
    if hzz <= 0:
        raise ValueError('The internal mode must have positive stiffness')
    return [[hqq[i][j]-F(hqz[i])*F(hqz[j])/hzz for j in range(2)] for i in range(2)]


def radial_ward(x, y, beta=1):
    h, dh, _ = metric_jets(radial(x, y, beta))
    rotation = nt.mat([[0, -1], [1, 0]])
    directional = nt.add(nt.scale(dh[0], -F(y)), nt.scale(dh[1], F(x)))
    return nt.add(directional, nt.add(nt.mul(nt.tr(rotation), h), nt.mul(h, rotation)))


def bi_radicand(d, b):
    cross = [d[1]*b[2]-d[2]*b[1], d[2]*b[0]-d[0]*b[2], d[0]*b[1]-d[1]*b[0]]
    return 1+sum(v*v for v in d)+sum(v*v for v in b)+sum(v*v for v in cross)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def verified_source_digest():
    path = HERE/'NSR1_SOURCE_PINS.json'
    raw = path.read_bytes()
    pins = json.loads(raw)
    for rel, expected in pins['sha256'].items():
        actual = hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError(f'Source mismatch: {rel}')
    return hashlib.sha256(raw).hexdigest()


def evidence():
    require_runtime()
    sources = verified_source_digest()
    zero = nt.mat([[0, 0], [0, 0]])
    bi_cases = []
    for t in (F(0), F(1, 4), F(1, 3), F(1, 2), F(3, 4)):
        for sign in (1, -1):
            x, y = t, sign*t
            u = born_infeld(x, y)
            info = analyze(u)
            expected = (x*x+y*y+2*x*x*y*y)/((1+t*t)*(1-x*x*y*y)**2)
            assert info['kappa'] == expected
            w = wave_data(info['H'])
            assert abs(w['v_plus']) <= 1 and abs(w['v_minus']) <= 1
            assert reconstruct(w['v_plus'], w['v_minus'], w['Z']) == info['H']
            matrix = [info['H'][1], info['H'][0]]
            for velocity, sign_mode in ((w['v_plus'], 1), (w['v_minus'], -1)):
                vector = [[F(1)], [sign_mode*w['Z']]]
                assert nt.mul(matrix, vector) == nt.scale(vector, velocity)
                fields = nt.mul(info['H'], vector)
                assert fields[0][0]/fields[1][0] == sign_mode*w['Z']
            resp = responses(info['H'])
            assert resp['epsilon_h']/resp['epsilon_y'] == resp['mu_e']/resp['mu_x'] == resp['gain']
            assert resp['gain'] == 1/(1-w['r']**2)
            bi_cases.append(dict(x=x, y=y, kappa=info['kappa'], wave=w, responses=resp))

    half = analyze(born_infeld(F(1, 2), F(1, 2)))
    linear = analyze(quadratic_tangent(born_infeld(F(1, 2), F(1, 2))))
    assert half['H'] == linear['H'] and half['kappa'] == F(128, 225)
    assert half['F'] != zero and linear['F'] == zero
    assert half['H'] == nt.mat([[F(4, 5), F(1, 5)], [F(1, 5), F(4, 5)]])

    # Variable common stiffness: this matrix family need not remain Hessian.
    w, dw = F(7, 3), [F(2, 5), F(-3, 7)]
    h = nt.scale(half['H'], w)
    dh = [nt.add(nt.scale(half['dH'][i], w), nt.scale(half['H'], dw[i])) for i in range(2)]
    hi = nt.inv(h)
    xs = [nt.mul(hi, d) for d in dh]
    rescaled_f = nt.scale(nt.sub(nt.mul(xs[0], xs[1]), nt.mul(xs[1], xs[0])), F(-1, 4))
    assert rescaled_f == half['F']
    assert shape_kappa(h, dh) == half['kappa']/w

    radial_cases = []
    for beta in (F(1, 2), F(1), F(3)):
        for t in (F(0), F(1, 2), F(1)):
            kappa = analyze(radial(t, 0, beta))['kappa']
            expected = -2*beta*beta*t*t/((1+3*beta*t*t)**2*(1+beta*t*t)**2)
            assert kappa == expected and radial_ward(t, F(1, 3), beta) == zero
            radial_cases.append(dict(beta=beta, t=t, kappa=kappa))
    origin = born_infeld(0, 0)
    assert all(i % 2 == j % 2 == 0 for i, j in origin.c)
    assert origin.derivative(0, 0, 1, 1) == 1
    assert analyze(origin)['F'] == zero

    # Full-vector duality is checked in all six channels, not on a false
    # rotation of only the transverse pair.
    d, b = [F(1, 3), F(2, 5), F(3, 7)], [F(2, 3), F(-1, 5), F(1, 7)]
    cosine, sine = F(3, 5), F(4, 5)
    dr = [cosine*x-sine*y for x, y in zip(d, b)]
    br = [sine*x+cosine*y for x, y in zip(d, b)]
    assert bi_radicand(dr, br) == bi_radicand(d, b)

    clamped_info = analyze(clamped(F(1, 2), F(1, 2), F(1, 4)))
    relaxed_info = analyze(relaxed(F(1, 2), F(1, 2)))
    assert schur_one(clamped_info['H'], [F(-1, 2), F(-1, 2)], 1) == relaxed_info['H']
    assert clamped_info['kappa'] == 0 and relaxed_info['kappa'] == F(128, 25)

    # Compare many Hessian jets rather than only named physical examples.
    nt_result = nt.run()
    for fixture in nt_result['fixtures']:
        h = nt.mat(fixture['H'])
        dh = [nt.mat(d) for d in fixture['dH']]
        assert shape_kappa(h, dh) == F(fixture['kappa'])

    rejected = []
    for name, fn in (
        ('degenerate_response', lambda: analyze(born_infeld(1, 1))),
        ('indefinite_response', lambda: analyze(born_infeld(2, 2))),
        ('unstable_internal_mode', lambda: schur_one(half['H'], [1, 1], -1)),
        ('invalid_wave_sector', lambda: reconstruct(2, 1, 1)),
    ):
        try:
            fn()
        except ValueError:
            rejected.append(name)
        else:
            raise AssertionError(f'Missing rejection: {name}')

    return encode(dict(
        status='PASS_NSR1_EXACT_RESPONSE_CONTROLS', source_manifest_sha256=sources,
        runtime='Python 3.12', born_infeld_cases=bi_cases,
        born_infeld_half_bias=half, point_matched_quadratic=linear,
        rotation_cases=radial_cases, common_scale_curvature_unchanged=True,
        clamped=clamped_info, relaxed=relaxed_info,
        independent_NT_fixtures_compared=len(nt_result['fixtures']),
        inherited_NT_rational_status=nt_result['status'],
        full_six_channel_duality_radicand_preserved=True,
        fourth_order_zero_bias_mixed_response=origin.derivative(0, 0, 1, 1),
        rejected_domains=rejected,
        not_established=['experimental_validation', 'literature_priority',
                         'automatic_optical_phase', 'photon_mass', 'Yang_Mills_continuum_gap']
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = evidence()
    raw = (json.dumps(result, indent=2, sort_keys=True)+'\n').encode()
    digest = hashlib.sha256(raw).hexdigest()
    output, expected = HERE/'NSR1_RESULT.json', HERE/'EXPECTED_NSR1.sha256'
    if args.write:
        output.write_bytes(raw)
        expected.write_text(digest+'\n')
    elif output.read_bytes() != raw or expected.read_text().strip() != digest:
        raise RuntimeError('Evidence differs; inspect changes before using --write')
    print(json.dumps(dict(status=result['status'], sha256=digest,
                          BI_cases=len(result['born_infeld_cases']),
                          rotation_cases=len(result['rotation_cases']),
                          NT_fixtures=result['independent_NT_fixtures_compared'])))


if __name__ == '__main__':
    main()
