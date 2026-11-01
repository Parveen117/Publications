"""Independent exact geometry checks for NSB1 (Python 3.12 stdlib).

The derivative oracle stores actual partial derivatives, uses the multivariate
Leibniz rule and set partitions for composition, and differentiates the full
Christoffel formula. It does not call the certificate's Jet or geometry helpers.
"""

from fractions import Fraction as Q
from functools import lru_cache
import importlib.util
from math import comb
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).resolve().parents[1] / "certificates/nsb1_geometry.py"
SPEC = importlib.util.spec_from_file_location("nsb1_geometry_under_test", MODULE_PATH)
geometry = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(geometry)


@lru_cache(None)
def partitions(n):
    """Set partitions of derivative positions, including multiplicities."""
    if n == 0:
        return ((),)
    result = []
    for old in partitions(n - 1):
        result.append(old + ((n - 1,),))
        for index in range(len(old)):
            result.append(old[:index] + (old[index] + (n - 1,),) + old[index + 1:])
    return tuple(result)


class Derivatives:
    order = 4

    def __init__(self, value=0, entries=None):
        self.d = {(0, 0): Q(value)} if entries is None else dict(entries)

    @classmethod
    def variable(cls, value, axis):
        return cls(entries={(0, 0): Q(value), (1 - axis, axis): Q(1)})

    def value(self, i=0, j=0):
        return self.d.get((i, j), Q(0))

    @staticmethod
    def coerce(other):
        return other if isinstance(other, Derivatives) else Derivatives(other)

    def __add__(self, other):
        other = self.coerce(other)
        return Derivatives(entries={key: self.d.get(key, 0) + other.d.get(key, 0)
                                    for key in self.d.keys() | other.d.keys()})

    __radd__ = __add__

    def __neg__(self):
        return Derivatives(entries={key: -value for key, value in self.d.items()})

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        other = self.coerce(other)
        result = {}
        for n in range(self.order + 1):
            for i in range(n + 1):
                j = n - i
                result[i, j] = sum(
                    Q(comb(i, a) * comb(j, b)) * self.value(a, b)
                    * other.value(i - a, j - b)
                    for a in range(i + 1) for b in range(j + 1)
                )
        return Derivatives(entries=result)

    __rmul__ = __mul__

    def compose(self, outer):
        """outer[k] is the kth univariate derivative at this value."""
        result = {(0, 0): Q(outer[0])}
        for n in range(1, self.order + 1):
            for i in range(n + 1):
                j = n - i
                axes = (0,) * i + (1,) * j
                total = Q(0)
                for partition in partitions(n):
                    term = Q(outer[len(partition)])
                    for block in partition:
                        dx = sum(axes[position] == 0 for position in block)
                        term *= self.value(dx, len(block) - dx)
                    total += term
                result[i, j] = total
        return Derivatives(entries=result)

    def reciprocal(self):
        a = self.value()
        return self.compose([1/a, -1/a**2, 2/a**3, -6/a**4, 24/a**5])

    def __truediv__(self, other):
        return self * self.coerce(other).reciprocal()

    def __pow__(self, exponent):
        if exponent < 0:
            return self.reciprocal() ** (-exponent)
        answer = Derivatives(1)
        for _ in range(exponent):
            answer *= self
        return answer

    def log_without_constant(self):
        a = self.value()
        return self.compose([0, 1/a, -1/a**2, 2/a**3, -6/a**4])

    def atan_without_constant(self):
        a = self.value()
        d = 1 + a*a
        return self.compose([0, 1/d, -2*a/d**2, (6*a*a - 2)/d**3,
                             24*a*(1-a*a)/d**4])

    def exp_at_zero(self):
        if self.value() != 0:
            raise ValueError("This rational oracle expands exponentials at zero")
        return self.compose([1, 1, 1, 1, 1])

    def partial(self, *axes):
        return self.value(axes.count(0), axes.count(1))


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def subtract(a, b):
    return [[a[i][j] - b[i][j] for j in range(2)] for i in range(2)]


def scale(a, coefficient):
    return [[coefficient * value for value in row] for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def determinant(a):
    return a[0][0]*a[1][1] - a[0][1]*a[1][0]


def metric_oracle(energy):
    """Full Christoffel derivative route, without the commutator shortcut."""
    h = [[energy.partial(i, j) for j in range(2)] for i in range(2)]
    dh = [[[energy.partial(p, i, j) for j in range(2)] for i in range(2)]
          for p in range(2)]
    ddh = [[[[energy.partial(p, q, i, j) for j in range(2)] for i in range(2)]
            for q in range(2)] for p in range(2)]
    delta = determinant(h)
    inverse = scale([[h[1][1], -h[0][1]], [-h[1][0], h[0][0]]], 1/delta)
    dinverse = [scale(matmul(matmul(inverse, p), inverse), -1) for p in dh]
    gamma = [[[sum(inverse[k][l] * (dh[i][l][j] + dh[j][i][l] - dh[l][i][j])
                      for l in range(2))/2 for j in range(2)]
              for k in range(2)] for i in range(2)]
    dgamma = [[[[sum(
        dinverse[p][k][l] * (dh[i][l][j] + dh[j][i][l] - dh[l][i][j])
        + inverse[k][l] * (ddh[p][i][l][j] + ddh[p][j][i][l] - ddh[p][l][i][j])
        for l in range(2))/2 for j in range(2)] for k in range(2)]
        for i in range(2)] for p in range(2)]
    product = subtract(matmul(gamma[0], gamma[1]), matmul(gamma[1], gamma[0]))
    curvature = [[dgamma[0][1][i][j] - dgamma[1][0][i][j] + product[i][j]
                  for j in range(2)] for i in range(2)]
    return dict(H=h, dH=dh, ddH=ddh, F=curvature, delta=delta,
                kappa=matmul(h, curvature)[0][1]/delta)


def halfplane_energy(c, x, y, curved):
    x, y = Derivatives.variable(x, 0), Derivatives.variable(y, 1)
    common = Q(c)*(x*x+y*y)/2
    # Discarded logarithm/atan constants multiply only affine coordinates.
    if curved:
        z = y/x
        return common + x*(1+z*z).log_without_constant()/2 - y*z.atan_without_constant()
    return common - x*x.log_without_constant() + x


def exponential_energy(a, b, sigma):
    x, y = Derivatives.variable(0, 0), Derivatives.variable(0, 1)
    return (2+sigma)*x*x/2 + (2-sigma)*y*y/2 + a*(x+y).exp_at_zero() + b*(x-y).exp_at_zero()


def radial_energy(m, eta, x, y):
    x, y = Derivatives.variable(x, 0), Derivatives.variable(y, 1)
    radius2 = x*x+y*y
    return radius2/2 + Q(eta)*radius2**m/(2*m)


class IndependentGeometryTests(unittest.TestCase):
    def assert_response(self, actual, expected):
        for key in ("H", "dH", "ddH", "F", "delta", "kappa"):
            self.assertEqual(actual[key], expected[key], key)

    def test_halfplane_potentials_have_equal_spectra_but_different_curvature(self):
        for c, x, y in [(Q(1), Q(2), Q(0)), (Q(3, 2), Q(2), Q(3, 2)),
                        (Q(2), Q(3, 4), Q(-2, 3))]:
            with self.subTest(c=c, x=x, y=y):
                flat = geometry.halfplane_response(c, x, y, kind="flat")
                curved = geometry.halfplane_response(c, x, y, kind="curved")
                for actual, kind in [(flat, False), (curved, True)]:
                    self.assert_response(actual, metric_oracle(halfplane_energy(c, x, y, kind)))
                    self.assertEqual(actual["H_characteristic"], [1, -2*c+1/x, c*(c-1/x)])
                self.assertEqual(flat["H_characteristic"], curved["H_characteristic"])
                self.assertEqual(flat["F"], [[0, 0], [0, 0]])
                self.assertGreater(curved["kappa"], 0)
                self.assertEqual(curved["kappa"], 1/(4*c*(x*x+y*y)*(c*x-1)**2))
        anchor = geometry.halfplane_response(1, 2, 0)
        self.assertEqual(anchor["F"], [[0, Q(1, 32)], [Q(-1, 16), 0]])

    def test_halfplane_rejects_degenerate_or_wrong_patch(self):
        for args in [(0, 2, 0), (-1, 2, 0), (1, 1, 0), (2, Q(1, 2), 0), (1, -2, 1)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                geometry.halfplane_response(*args)

    def test_exponential_spectra_and_waves_miss_signed_response_geometry(self):
        reflection = [[0, 1], [1, 0]]
        for a, b in [(Q(1), Q(1)), (Q(2), Q(1, 3)), (Q(1, 5), Q(4))]:
            plus, minus = [geometry.exponential_response(a, b, s) for s in (1, -1)]
            for s, response in [(1, plus), (-1, minus)]:
                self.assert_response(response, metric_oracle(exponential_energy(a, b, s)))
                h, delta = response["H"], response["delta"]
                self.assertEqual(response["H_characteristic"], [1, -h[0][0]-h[1][1], delta])
                self.assertEqual(response["wave_characteristic"], [1, -2*h[0][1], -delta])
                self.assertEqual(response["impedance_squared"], h[0][0]/h[1][1])
            self.assertEqual(plus["H_characteristic"], minus["H_characteristic"])
            self.assertEqual(plus["wave_characteristic"], minus["wave_characteristic"])
            self.assertEqual(plus["kappa"], -minus["kappa"])
            self.assertGreater(plus["kappa"], 0)
            self.assertEqual(plus["impedance_squared"]*minus["impedance_squared"], 1)
            # Equal spectra here do NOT establish gauge-inequivalent bundles.
            for key in ("H", "F"):
                self.assertEqual(minus[key], matmul(matmul(reflection, plus[key]), reflection))
            self.assertTrue(plus["bundle_reflection_conjugacy"])

    def test_exponential_parameter_domain(self):
        for args in [(0, 1, 1), (1, -1, 1), (1, 1, 0), (1, 1, 2)]:
            with self.subTest(args=args), self.assertRaises(ValueError):
                geometry.exponential_response(*args)

    def test_dual_flat_connections_and_curved_metric_mean(self):
        reference = metric_oracle(halfplane_energy(Q(1), Q(2), Q(1, 2), True))
        h, (hx, hy) = reference["H"], reference["dH"]
        for q in [Q(-1), Q(0), Q(1, 2), Q(1), Q(2)]:
            data = geometry.scaled_connection_curvature(h, hx, hy, q)
            self.assertEqual(data["F"], scale(reference["F"], 4*q*(1-q)))
            for index, derivative in enumerate([hx, hy]):
                connection = data["connection"][index]
                compatible_sum = [[u+v for u, v in zip(row1, row2)] for row1, row2 in zip(
                    matmul(transpose(connection), h), matmul(h, connection))]
                self.assertEqual(compatible_sum, scale(derivative, 2*q))
            if q in (0, 1):
                self.assertEqual(data["F"], [[0, 0], [0, 0]])
            elif q == Q(1, 2):
                self.assertNotEqual(data["F"], [[0, 0], [0, 0]])

    def test_mixed_protocol_reconstructs_each_constraint_step(self):
        for x, y, a, b in [(Q(0), Q(0), Q(1, 4), Q(1, 3)),
                            (Q(1), Q(1, 4), Q(-1, 3), Q(1, 5)),
                            (Q(0), Q(1, 4), Q(0), Q(1, 5))]:
            data = geometry.mixed_protocol(x, y, a, b)
            first = [x+a, y]
            second = [first[0], ((1+first[0])*first[1]+b)/(1+first[0])]
            third = [second[0]-a, second[1]]
            last = [third[0], ((1+third[0])*third[1]-b)/(1+third[0])]
            self.assertEqual(data["points"], [[x, y], first, second, third, last])
            self.assertEqual(data["displacement"], [last[0]-x, last[1]-y])
            self.assertEqual(data["bracket"], [0, -1/(1+x)**2])
            dx, dy = Derivatives.variable(x, 0), Derivatives.variable(y, 1)
            self.assert_response(data["response"], metric_oracle((dx*dx+(1+dx)*dy*dy)/2))
            self.assertFalse(data["endpoint_return_is_parallel_transport"])
        # A stable endpoint alone does not license a path through degeneracy.
        with self.assertRaises(ValueError):
            geometry.mixed_protocol(0, 0, Q(1, 10), Q(11, 10))

    def test_radial_curvature_is_computed_from_the_energy(self):
        for m in (2, 3, 5, 7):
            for eta, t in [(Q(1), Q(0)), (Q(2, 3), Q(1, 2)), (Q(3, 2), Q(-2, 3))]:
                expected = metric_oracle(radial_energy(m, eta, t, 0))
                self.assert_response(geometry.radial_response(m, eta, t, 0), expected)
                self.assertEqual(geometry.radial_curvature(m, eta, t), expected["kappa"])
                if t:
                    self.assertLess(expected["kappa"], 0)

    def test_rotational_symmetry_does_not_flatten_the_radial_response(self):
        rotation = [[Q(3, 5), Q(-4, 5)], [Q(4, 5), Q(3, 5)]]
        for m in (2, 4):
            axis = geometry.radial_response(m, Q(2, 3), 1, 0)
            rotated = geometry.radial_response(m, Q(2, 3), Q(3, 5), Q(4, 5))
            self.assert_response(rotated, metric_oracle(radial_energy(m, Q(2, 3), Q(3, 5), Q(4, 5))))
            self.assertEqual(rotated["H"], matmul(matmul(rotation, axis["H"]), transpose(rotation)))
            self.assertEqual(rotated["kappa"], axis["kappa"])
            self.assertNotEqual(rotated["kappa"], 0)

    def test_delayed_response_jet_order_from_polynomial_support(self):
        for m in range(2, 8):
            # Differentiate every monomial of (x²+y²)^m, independently of metadata.
            hessian_terms = []
            for k in range(m+1):
                powers = (2*k, 2*(m-k))
                for i in range(2):
                    for j in range(2):
                        remaining = list(powers)
                        coefficient = Q(comb(m, k), 2*m)
                        for axis in (i, j):
                            coefficient *= remaining[axis]
                            remaining[axis] -= 1
                        if coefficient:
                            hessian_terms.append((tuple(remaining), coefficient))
            first = min(sum(powers) for powers, _ in hessian_terms)
            metadata = geometry.radial_jet_order(m)
            self.assertEqual(metadata["first_response_jet_order"], first)
            self.assertEqual(metadata["matched_response_jet_through"], first-1)
            self.assertEqual(metadata["first_energy_degree"], first+2)
            # The leading derivatives of H have order first-1; their nonzero
            # commutator has twice that order because H(0)=I.
            self.assertEqual(metadata["first_curvature_order"], 2*(first-1))
            self.assertEqual(geometry.radial_response(m, 1, 0, 0)["H"], [[1, 0], [0, 1]])

    def test_radial_invalid_parameters(self):
        for m, eta in [(1, 1), (Q(5, 2), 1), (2, 0), (2, -1)]:
            with self.subTest(m=m, eta=eta), self.assertRaises(ValueError):
                geometry.radial_curvature(m, eta, Q(1, 2))


if __name__ == "__main__":
    unittest.main()
