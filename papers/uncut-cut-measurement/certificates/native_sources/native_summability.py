from __future__ import annotations

"""Exact native summability and recognition-energy dense core.

This module deliberately avoids ordinary complex numbers, NumPy, Hilbert
spaces, L-p spaces, operator norms, spectral theory, and Haar measure. Its
primitive data are rational oriented cuts and finite typed residue paths.

The companion theorem completes this dense core using two independently
declared Cauchy relations:

* cut-tail mass Cauchy presentations for the path algebra;
* closed-loop recognition-energy Cauchy presentations for recognized states.
"""

from dataclasses import dataclass
from fractions import Fraction
from numbers import Rational
from typing import Any, Iterable, Mapping


@dataclass(frozen=True)
class NativeCutScalar:
    """Exact radial plus quarter-turn cut scalar."""

    radial: Fraction = Fraction(0)
    turn: Fraction = Fraction(0)

    def __post_init__(self) -> None:
        object.__setattr__(self, "radial", Fraction(self.radial))
        object.__setattr__(self, "turn", Fraction(self.turn))

    @classmethod
    def zero(cls) -> "NativeCutScalar":
        return cls()

    @classmethod
    def one(cls) -> "NativeCutScalar":
        return cls(1, 0)

    @classmethod
    def iota(cls) -> "NativeCutScalar":
        return cls(0, 1)

    @classmethod
    def coerce(cls, value: Any) -> "NativeCutScalar":
        if isinstance(value, cls):
            return value
        if isinstance(value, Rational):
            return cls(Fraction(value), 0)
        raise TypeError(f"cannot coerce {type(value)!r} to NativeCutScalar")

    def __add__(self, other: Any) -> "NativeCutScalar":
        try:
            value = self.coerce(other)
        except TypeError:
            return NotImplemented
        return NativeCutScalar(
            self.radial + value.radial,
            self.turn + value.turn,
        )

    def __radd__(self, other: Any) -> "NativeCutScalar":
        return self + other

    def __neg__(self) -> "NativeCutScalar":
        return NativeCutScalar(-self.radial, -self.turn)

    def __sub__(self, other: Any) -> "NativeCutScalar":
        try:
            value = self.coerce(other)
        except TypeError:
            return NotImplemented
        return self + (-value)

    def __rsub__(self, other: Any) -> "NativeCutScalar":
        try:
            value = self.coerce(other)
        except TypeError:
            return NotImplemented
        return value - self

    def __mul__(self, other: Any) -> "NativeCutScalar":
        try:
            value = self.coerce(other)
        except TypeError:
            return NotImplemented
        return NativeCutScalar(
            self.radial * value.radial - self.turn * value.turn,
            self.radial * value.turn + self.turn * value.radial,
        )

    def __rmul__(self, other: Any) -> "NativeCutScalar":
        return self * other

    def dagger(self) -> "NativeCutScalar":
        """Reverse quarter-turn orientation."""

        return NativeCutScalar(self.radial, -self.turn)

    def norm_square(self) -> Fraction:
        """Return the radial remainder after orientation reversal."""

        return self.radial * self.radial + self.turn * self.turn

    def diamond_mass(self) -> Fraction:
        """Exact rational tail gauge used before any L1 notation."""

        return abs(self.radial) + abs(self.turn)

    def is_zero(self) -> bool:
        return self.radial == 0 and self.turn == 0


@dataclass(frozen=True, order=True)
class NativeArrow:
    """One typed path with a finite residue channel."""

    source: int
    target: int
    residue: int


class NativePathAlgebra:
    """Finite-support typed residue paths over exact cut scalars."""

    def __init__(self, residue_modulus: int = 2) -> None:
        modulus = int(residue_modulus)
        if modulus <= 0:
            raise ValueError("residue_modulus must be positive")
        self.residue_modulus = modulus

    def arrow(
        self,
        source: int,
        target: int,
        residue: int = 0,
    ) -> NativeArrow:
        return NativeArrow(
            int(source),
            int(target),
            int(residue) % self.residue_modulus,
        )

    def identity_arrow(self, state: int) -> NativeArrow:
        return self.arrow(state, state, 0)

    def compose(
        self,
        left: NativeArrow,
        right: NativeArrow,
    ) -> NativeArrow | None:
        """Compose as left after right."""

        if right.target != left.source:
            return None
        return self.arrow(
            right.source,
            left.target,
            right.residue + left.residue,
        )

    def dagger_arrow(self, arrow: NativeArrow) -> NativeArrow:
        return self.arrow(arrow.target, arrow.source, -arrow.residue)

    def zero(self) -> "FiniteNativePath":
        return FiniteNativePath(self)

    def basis(
        self,
        source: int,
        target: int,
        residue: int = 0,
    ) -> "FiniteNativePath":
        return FiniteNativePath(
            self,
            {self.arrow(source, target, residue): NativeCutScalar.one()},
        )

    def projection(self, state: int) -> "FiniteNativePath":
        return self.basis(state, state, 0)


class FiniteNativePath:
    """Finite linear combination of typed residue paths."""

    def __init__(
        self,
        algebra: NativePathAlgebra,
        coefficients: Mapping[NativeArrow, Any] | None = None,
    ) -> None:
        self.algebra = algebra
        cleaned: dict[NativeArrow, NativeCutScalar] = {}
        for arrow, raw in (coefficients or {}).items():
            normalized = algebra.arrow(
                arrow.source,
                arrow.target,
                arrow.residue,
            )
            value = NativeCutScalar.coerce(raw)
            if value.is_zero():
                continue
            cleaned[normalized] = (
                cleaned.get(normalized, NativeCutScalar.zero()) + value
            )
        self.coefficients = {
            arrow: value
            for arrow, value in cleaned.items()
            if not value.is_zero()
        }

    def _require_same(self, other: "FiniteNativePath") -> None:
        if self.algebra is not other.algebra:
            raise ValueError("paths belong to different native algebras")

    def __add__(self, other: "FiniteNativePath") -> "FiniteNativePath":
        self._require_same(other)
        result = dict(self.coefficients)
        for arrow, value in other.coefficients.items():
            result[arrow] = (
                result.get(arrow, NativeCutScalar.zero()) + value
            )
        return FiniteNativePath(self.algebra, result)

    def __sub__(self, other: "FiniteNativePath") -> "FiniteNativePath":
        return self + (-other)

    def __neg__(self) -> "FiniteNativePath":
        return FiniteNativePath(
            self.algebra,
            {
                arrow: -value
                for arrow, value in self.coefficients.items()
            },
        )

    def __rmul__(self, scalar: Any) -> "FiniteNativePath":
        coefficient = NativeCutScalar.coerce(scalar)
        return FiniteNativePath(
            self.algebra,
            {
                arrow: coefficient * value
                for arrow, value in self.coefficients.items()
            },
        )

    def __matmul__(self, other: "FiniteNativePath") -> "FiniteNativePath":
        """Convolve finite paths: left after right."""

        self._require_same(other)
        result: dict[NativeArrow, NativeCutScalar] = {}
        for left, left_value in self.coefficients.items():
            for right, right_value in other.coefficients.items():
                composed = self.algebra.compose(left, right)
                if composed is None:
                    continue
                result[composed] = (
                    result.get(composed, NativeCutScalar.zero())
                    + left_value * right_value
                )
        return FiniteNativePath(self.algebra, result)

    def dagger(self) -> "FiniteNativePath":
        return FiniteNativePath(
            self.algebra,
            {
                self.algebra.dagger_arrow(arrow): value.dagger()
                for arrow, value in self.coefficients.items()
            },
        )

    def canonical_trace(self) -> NativeCutScalar:
        return sum(
            (
                value
                for arrow, value in self.coefficients.items()
                if arrow.source == arrow.target and arrow.residue == 0
            ),
            NativeCutScalar.zero(),
        )

    def pairing(self, other: "FiniteNativePath") -> NativeCutScalar:
        self._require_same(other)
        return (self.dagger() @ other).canonical_trace()

    def coefficient_pairing(
        self,
        other: "FiniteNativePath",
    ) -> NativeCutScalar:
        self._require_same(other)
        arrows = set(self.coefficients) | set(other.coefficients)
        return sum(
            (
                self.coefficients.get(
                    arrow,
                    NativeCutScalar.zero(),
                ).dagger()
                * other.coefficients.get(arrow, NativeCutScalar.zero())
                for arrow in arrows
            ),
            NativeCutScalar.zero(),
        )

    def recognition_energy(self) -> Fraction:
        value = self.pairing(self)
        if value.turn != 0 or value.radial < 0:
            raise ArithmeticError(
                "recognition energy must be nonnegative radial"
            )
        return value.radial

    def cut_tail_mass(self) -> Fraction:
        return sum(
            (
                value.diamond_mass()
                for value in self.coefficients.values()
            ),
            Fraction(0),
        )

    def exact_equal(self, other: "FiniteNativePath") -> bool:
        self._require_same(other)
        return self.coefficients == other.coefficients

    def is_zero(self) -> bool:
        return not self.coefficients


def scalar_diamond_submultiplicative(
    left: NativeCutScalar,
    right: NativeCutScalar,
) -> bool:
    return (left * right).diamond_mass() <= (
        left.diamond_mass() * right.diamond_mass()
    )


def recognition_parallelogram_remainder(
    left: NativeCutScalar,
    right: NativeCutScalar,
) -> Fraction:
    """Return 2N(x)+2N(y)-N(x+y)=N(x-y)."""

    return (
        2 * left.norm_square()
        + 2 * right.norm_square()
        - (left + right).norm_square()
    )


def weighted_cut_cauchy_schwarz_identity(
    values: Iterable[NativeCutScalar],
    weights: Iterable[Fraction],
) -> tuple[Fraction, Fraction]:
    """Return both sides of an exact weighted cut-square identity.

    For positive rational weights d_i and cut scalars x_i,

        D sum N(x_i)/d_i - N(sum x_i)
        = sum_{i<j} N(d_j x_i-d_i x_j)/(d_i d_j).

    No square root or pre-existing inner-product theorem is used.
    """

    xs = tuple(values)
    ws = tuple(Fraction(value) for value in weights)
    if len(xs) != len(ws):
        raise ValueError("values and weights must have equal length")
    if any(weight <= 0 for weight in ws):
        raise ValueError("weights must be positive")

    total_weight = sum(ws, Fraction(0))
    total_value = sum(xs, NativeCutScalar.zero())
    left = total_weight * sum(
        (
            value.norm_square() / weight
            for value, weight in zip(xs, ws)
        ),
        Fraction(0),
    ) - total_value.norm_square()
    right = sum(
        (
            (
                ws[j] * xs[i] - ws[i] * xs[j]
            ).norm_square()
            / (ws[i] * ws[j])
            for i in range(len(xs))
            for j in range(i + 1, len(xs))
        ),
        Fraction(0),
    )
    return left, right


def native_action_bound(
    action: FiniteNativePath,
    state: FiniteNativePath,
) -> tuple[Fraction, Fraction]:
    """Return E(action*state) and M(action)^2 E(state)."""

    product = action @ state
    left = product.recognition_energy()
    right = (
        action.cut_tail_mass() ** 2
        * state.recognition_energy()
    )
    return left, right


def endpoint_blind_shadow(
    path: FiniteNativePath,
) -> dict[tuple[int, int], NativeCutScalar]:
    """Forget residue while retaining typed endpoints."""

    result: dict[tuple[int, int], NativeCutScalar] = {}
    for arrow, value in path.coefficients.items():
        key = (arrow.source, arrow.target)
        result[key] = (
            result.get(key, NativeCutScalar.zero()) + value
        )
    return {
        key: value
        for key, value in result.items()
        if not value.is_zero()
    }


def geometric_partial(
    algebra: NativePathAlgebra,
    maximum_index: int,
) -> FiniteNativePath:
    """Exact alternating-residue Cauchy presentation."""

    result = algebra.zero()
    for index in range(maximum_index + 1):
        coefficient = NativeCutScalar(
            Fraction(1, 2 ** (index + 1)),
            0,
        )
        result = (
            result
            + coefficient * algebra.basis(0, index, index % 2)
        )
    return result


def exact_geometric_mass_tail(
    after_index: int,
    through_index: int,
) -> Fraction:
    if through_index <= after_index:
        return Fraction(0)
    return sum(
        (
            Fraction(1, 2 ** (index + 1))
            for index in range(after_index + 1, through_index + 1)
        ),
        Fraction(0),
    )


def exact_geometric_energy_tail(
    after_index: int,
    through_index: int,
) -> Fraction:
    if through_index <= after_index:
        return Fraction(0)
    return sum(
        (
            Fraction(1, 4 ** (index + 1))
            for index in range(after_index + 1, through_index + 1)
        ),
        Fraction(0),
    )
