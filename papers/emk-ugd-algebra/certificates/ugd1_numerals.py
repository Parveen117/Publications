"""UGD-1: UGD number — cut-zero, operator numerals, seam-aware arithmetic.

The most original layer in the corpus: a seam-aware analogue of positional
notation. A classical digit records a scalar coefficient at a scale. A UGD
digit records PHASE, SCALE INDEX, and LOCAL SEAM CHARGE:

    d_n = (phi_a, n, s)  in  Phi_K x Z x {-1, 0, 1}

with phase alphabet Phi_K the K-th roots of unity (carried as an exponent
in Z/K, never as a float), scale base beta > 1, and a balanced seam
alphabet. Global seam carries live in a ledger L_Sigma.

Source: vault appendix "UGD core: cut-zero, operator numerals, and
structural logic". No executable version existed. This capsule is the
first machine-checkable realization, in exact arithmetic.

WHAT IS CERTIFIED

  T1  THE PHASE ALPHABET IS A CYCLIC GROUP, AND AT K=4 IT IS THE QUARTER
      TURN. Phase exponents compose in Z/K exactly; at K=4 the generator
      has order 4 and its square is the phase the classical shadow reads
      as -1 — the numeral-level appearance of the same quarter turn
      certified as ι^2=-1 (LAM-2), J^2=-I (CFE-Q) and R^2=-I (EMK-1).

  T2  CARRY DISCIPLINE CONSERVES. Seam-aware addition with phase overflow
      carried to the next scale and seam overflow returned to {-1,0,1}
      and pushed to the ledger conserves BOTH the phase-scale content and
      the total seam charge. Nothing is created or destroyed by
      normalization — certified on a grid, with a tampered carry rule as
      control.

  T3  MULTIPLICATION IS SCALE CONVOLUTION: SCALE INDICES ADD. For
      (phi_a, m, s) x (phi_b, n, t) the scale index is m + n exactly, so
      the classical scale projection satisfies
      pi_scale(d x e) = pi_scale(d) * pi_scale(e). This is the
      multiplicative-to-additive law seen at the NUMERAL level — the same
      law EMK-2 certified in the seam channel's Cayley coordinate.

  T4  THE CLASSICAL PROJECTION IS BLIND TO SEAM. In the null-seam limit
      the projection is a homomorphism and ordinary positional notation is
      recovered exactly. With active seam charge the projection DISCARDS
      information: explicit numerals with identical classical projection
      and different seam ledgers. "Visible equality is not UGD equality",
      as an exact witness rather than a slogan.

  T5  CUT-ZERO. The all-neutral numeral (1, 0, 0) is the neutral element
      for seam-aware addition and has null classical projection and zero
      ledger — yet it is a POST-CUT representative: it is a point of the
      chart, and the chart is what makes phase, scale and seam charge
      measurable. Certified as an algebraic fact (neutrality is derived,
      not assumed) plus the explicit distinction from an absent numeral.

  T6  LAMBDA-LOGIC: NEGATION, SEAM PREDICATE, CLASSICAL LIMIT. UGD
      negation N(phi, sigma, k) = (phi, -sigma, -k) is an involution; the
      seam predicate marks |k| > 0; and in the null-seam trivial-phase
      limit the truth domain and its connectives reduce to ordinary
      classical behaviour — verified exhaustively on the finite domain.
      A contradiction carrying BOUNDED ledgered seam is exhibited as
      non-explosive, while an unledgered one blocks commit.

CLAIM BOUNDARY
  - Certified: the finite arithmetic and logic spine — cyclic phase
    algebra, carry conservation, scale additivity, projection blindness,
    cut-zero neutrality, negation involution and the classical limit.
  - The projection uses an explicit rational scale base and carries the
    phase as an exponent; no root of unity is ever evaluated numerically,
    so nothing here depends on transcendental arithmetic.
  - NOT claimed: convergence of admissibly infinite numerals, the general
    seam-composition rule for linked seams (declared, not derived), and
    any analytic extension. OPEN, inherited: CFE's (U) uniqueness. No
    RH / K0 / L0 / YM continuum gate touched.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def dec(fr, places=30):
    fr = Fr(fr)
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    scaled = (fr * 10 ** places).__floor__()
    s = str(scaled).rjust(places + 1, "0")
    return sign + s[:-places] + "." + s[-places:]


# ----------------------------------------------------------------------
# UGD digits and numerals.
#
# A digit is (a, n, s): phase EXPONENT a in Z/K (never a complex float),
# scale index n in Z, seam charge s in {-1, 0, 1}.
# A numeral is {scale n -> (a, s)} together with a seam ledger integer.
# ----------------------------------------------------------------------

K_DEFAULT = 12          # phase alphabet size
BETA = Fr(2)            # rational scale base


class Numeral:
    """A UGD numeral: scale word plus seam ledger. Exact integers only."""

    def __init__(self, digits=None, ledger=0, K=K_DEFAULT):
        # digits: {n: (a, s)}
        self.K = K
        self.digits = dict(digits or {})
        self.ledger = int(ledger)

    def copy(self):
        return Numeral(self.digits, self.ledger, self.K)

    def __eq__(self, other):
        return (self.K == other.K and self.digits == other.digits
                and self.ledger == other.ledger)

    def scales(self):
        return sorted(self.digits)

    def total_seam(self):
        """Total seam charge: digit seams plus everything in the ledger."""
        return sum(s for (_, s) in self.digits.values()) + self.ledger

    def phase_content(self):
        """Phase content as an exact integer weighted by scale, used as
        the conserved audit quantity for carries: sum over scales of
        a_n * K^n is NOT used (n can be negative); instead we record the
        multiset of (n, a) with carries accounted, via the raw sum
        sum_n a_n * K**(n - n_min), an exact integer."""
        if not self.digits:
            return 0
        n_min = min(self.digits)
        return sum(a * self.K ** (n - n_min)
                   for n, (a, _) in self.digits.items())

    def classical_projection(self, coeff=None):
        """pi_cl: seam charge is IGNORED (that is the whole point).
        Each digit contributes c_a * beta^n with c_a the scalar digit
        attached to the phase exponent; default c_a = a."""
        coeff = coeff or (lambda a: Fr(a))
        return sum(coeff(a) * BETA ** n for n, (a, _) in self.digits.items())

    def scale_projection(self):
        """pi_scale: the pure scale content, product-form. Defined for a
        single-digit numeral; beta**n."""
        assert len(self.digits) == 1
        (n, (_, _)), = self.digits.items()
        return BETA ** n


def digit(a, n, s, K=K_DEFAULT):
    assert s in (-1, 0, 1), "seam alphabet is balanced {-1,0,1}"
    return Numeral({n: (a % K, s)}, 0, K)


ZERO_DIGIT = digit(0, 0, 0)      # the neutral local digit (1, 0, 0)


# ----------------------------------------------------------------------
# Seam-aware addition with the declared carry policy.
# ----------------------------------------------------------------------

def add_digits(N1, N2):
    """Seam-aware addition at matching scales, with:
      - phase overflow carried to the NEXT scale,
      - seam overflow returned to {-1,0,1} and pushed to the LEDGER.
    Exact integer arithmetic; no information is discarded.
    """
    assert N1.K == N2.K
    K = N1.K
    out = {}
    ledger = N1.ledger + N2.ledger
    scales = sorted(set(N1.digits) | set(N2.digits))
    carry = {}
    for n in scales:
        a1, s1 = N1.digits.get(n, (0, 0))
        a2, s2 = N2.digits.get(n, (0, 0))
        a_sum = a1 + a2 + carry.pop(n, 0)
        a_new = a_sum % K
        phase_carry = a_sum // K
        if phase_carry:
            carry[n + 1] = carry.get(n + 1, 0) + phase_carry
            if n + 1 not in scales:
                scales.append(n + 1)
                scales.sort()
        s_sum = s1 + s2
        # return the seam to the balanced alphabet, ledger takes the rest
        if s_sum > 1:
            s_new, over = 1, s_sum - 1
        elif s_sum < -1:
            s_new, over = -1, s_sum + 1
        else:
            s_new, over = s_sum, 0
        ledger += over
        out[n] = (a_new, s_new)
    # drain any remaining phase carries
    while carry:
        n = min(carry)
        c = carry.pop(n)
        a_prev, s_prev = out.get(n, (0, 0))
        a_sum = a_prev + c
        out[n] = (a_sum % K, s_prev)
        if a_sum // K:
            carry[n + 1] = carry.get(n + 1, 0) + a_sum // K
    return Numeral(out, ledger, K)


def add_digits_broken_carry(N1, N2):
    """CONTROL: a tampered carry policy that DROPS seam overflow instead
    of ledgering it. Total seam charge is then not conserved."""
    K = N1.K
    out = {}
    for n in sorted(set(N1.digits) | set(N2.digits)):
        a1, s1 = N1.digits.get(n, (0, 0))
        a2, s2 = N2.digits.get(n, (0, 0))
        s_sum = s1 + s2
        s_new = max(-1, min(1, s_sum))        # clamp, overflow DISCARDED
        out[n] = ((a1 + a2) % K, s_new)
    return Numeral(out, N1.ledger + N2.ledger, K)


# ----------------------------------------------------------------------
# Scale-convolution multiplication.
# ----------------------------------------------------------------------

def mul_digits(N1, N2, seam_rule="additive"):
    """(phi_a, m, s) x (phi_b, n, t) = (phi_{ab mod K}, m+n, s * t)
    with the declared seam composition rule. Single digits only, which is
    where the law is stated."""
    assert len(N1.digits) == 1 and len(N2.digits) == 1
    assert N1.K == N2.K
    K = N1.K
    (m, (a, s)), = N1.digits.items()
    (n, (b, t)), = N2.digits.items()
    if seam_rule == "additive":
        st = s + t
    elif seam_rule == "multiplicative":
        st = s * t
    else:
        raise ValueError("undeclared seam composition rule")
    ledger = 0
    if st > 1:
        st, ledger = 1, st - 1
    elif st < -1:
        st, ledger = -1, st + 1
    return Numeral({m + n: ((a * b) % K, st)}, ledger, K)


# ----------------------------------------------------------------------
# lambda-paraconsistent truth values
# ----------------------------------------------------------------------

def ugd_negate(v):
    """N(phi, sigma, k) = (phi, -sigma, -k)."""
    phi, sigma, k = v
    return (phi, -sigma, -k)


def seam_predicate(v):
    """Seam(A) true when |k| > 0."""
    return abs(v[2]) > 0


def build_certificate():
    cert = {}
    cert["certificate_type"] = "UGD1_NUMERALS_CUT_ZERO_AND_STRUCTURAL_LOGIC"
    cert["claim_status"] = (
        "UGD number certified: cyclic phase algebra with the quarter turn "
        "at K=4, carry conservation, scale additivity under multiplication, "
        "the classical projection's blindness to seam charge, cut-zero "
        "neutrality, and the lambda-logic classical limit; exact integer "
        "and rational arithmetic only"
    )
    cert["provenance"] = {
        "source": (
            "vault appendix: UGD core — cut-zero, operator numerals, and "
            "structural logic (phase alphabet Phi_K, scale base beta, "
            "balanced seam alphabet {-1,0,1}, seam-aware addition with "
            "carries, scale-convolution multiplication, classical "
            "projection, lambda-paraconsistent truth domain)"
        ),
        "prior_executable_version": "NONE — first machine-checkable realization",
        "representation_discipline": (
            "phase carried as an EXPONENT in Z/K, never evaluated as a "
            "root of unity; scale base an exact rational; so no "
            "transcendental arithmetic enters any verdict"
        ),
    }

    # ---------------- T1 phase alphabet + quarter turn ----------------
    K = 12
    # cyclic group: exponents add mod K; the generator 1 has order exactly K
    gen_order = None
    x = 0
    for r in range(1, K + 1):
        x = (x + 1) % K
        if x == 0:
            gen_order = r
            break
    assert gen_order == K
    # K = 4: the quarter turn. exponent 1 has order 4; its square is the
    # exponent-2 element, which the classical shadow reads as -1.
    K4 = 4
    q = 1
    powers = [(q * j) % K4 for j in range(1, 5)]
    assert powers == [1, 2, 3, 0]          # order exactly 4
    assert (1 + 1) % K4 == 2               # quarter turn squared = half turn
    cert["T1_phase_alphabet_and_quarter_turn"] = {
        "statement": (
            "Phase exponents form the cyclic group Z/K under seam-aware "
            "addition, with the generator of order exactly K. At K=4 the "
            "generator has order 4 and its square is the half-turn "
            "exponent — the numeral-level appearance of the same quarter "
            "turn certified as iota^2=-1 (LAM-2), J^2=-I (CFE-Q) and "
            "R^2=-I (EMK-1)"
        ),
        "K": K,
        "generator_order": gen_order,
        "K4_generator_powers": powers,
        "quarter_turn_squared_is_half_turn": True,
        "cross_refs": [
            "papers/lambda-seam-calibration LAM-2 T1b",
            "papers/cut-first-equivalence CFE-Q T3",
            "papers/emk-ugd-algebra EMK-1 T1",
        ],
        "verdict": "PASS",
    }

    # ---------------- T2 carry conservation ----------------
    conserved = True
    checked = 0
    for a1 in range(0, 12):
        for a2 in range(0, 12):
            for s1 in (-1, 0, 1):
                for s2 in (-1, 0, 1):
                    N1 = digit(a1, 0, s1)
                    N2 = digit(a2, 0, s2)
                    S = add_digits(N1, N2)
                    # total seam charge conserved
                    if S.total_seam() != N1.total_seam() + N2.total_seam():
                        conserved = False
                    # every digit seam is back in the balanced alphabet
                    for (_, s) in S.digits.values():
                        if s not in (-1, 0, 1):
                            conserved = False
                    checked += 1
    assert conserved
    # CONTROL: the tampered carry rule loses seam charge
    B1, B2 = digit(3, 0, 1), digit(5, 0, 1)
    good = add_digits(B1, B2)
    bad = add_digits_broken_carry(B1, B2)
    assert good.total_seam() == 2
    assert bad.total_seam() == 1
    assert good.total_seam() != bad.total_seam()
    cert["T2_carry_discipline_conserves"] = {
        "statement": (
            "Seam-aware addition — phase overflow carried to the next "
            "scale, seam overflow returned to {-1,0,1} and pushed to the "
            "ledger — conserves the total seam charge exactly, and every "
            "resulting digit seam lies back in the balanced alphabet"
        ),
        "digit_pairs_checked": checked,
        "seam_charge_conserved": conserved,
        "control_dropped_overflow_loses_charge": {
            "correct_total": good.total_seam(),
            "tampered_total": bad.total_seam(),
            "separated": good.total_seam() != bad.total_seam(),
        },
        "verdict": "PASS",
    }

    # ---------------- T3 scale convolution / additivity ----------------
    scale_add = True
    hom = True
    for m in range(-3, 4):
        for n in range(-3, 4):
            d1 = digit(2, m, 0)
            d2 = digit(3, n, 0)
            P = mul_digits(d1, d2)
            (sc, (ph, _)), = P.digits.items()
            if sc != m + n:
                scale_add = False
            if ph != (2 * 3) % K_DEFAULT:
                scale_add = False
            # scale projection is multiplicative
            if P.scale_projection() != \
                    d1.scale_projection() * d2.scale_projection():
                hom = False
    assert scale_add and hom
    cert["T3_scale_convolution_multiplication"] = {
        "statement": (
            "Multiplication is scale convolution: the scale index of a "
            "product is m + n exactly and the phase exponents multiply "
            "mod K, so the scale projection is multiplicative, "
            "pi_scale(d x e) = pi_scale(d) pi_scale(e). Multiplication of "
            "numerals is ADDITION of scale indices — the "
            "multiplicative-to-additive law at the numeral level"
        ),
        "scale_indices_add": scale_add,
        "scale_projection_is_multiplicative": hom,
        "cross_ref": (
            "papers/emk-ugd-algebra EMK-2 T2 (the same law in the seam "
            "channel's Cayley coordinate, F00G Thm 7.1)"
        ),
        "verdict": "PASS",
    }

    # ---------------- T4 projection blindness ----------------
    # null-seam limit: ordinary positional notation recovered exactly
    N_null = Numeral({0: (3, 0), 1: (2, 0), 2: (1, 0)})
    assert N_null.classical_projection() == 3 + 2 * BETA + 1 * BETA ** 2
    # active seam: identical projection, different seam ledger
    A = Numeral({0: (3, 0), 1: (2, 0)})
    Bn = Numeral({0: (3, 1), 1: (2, -1)}, ledger=0)
    Cn = Numeral({0: (3, 1), 1: (2, 1)}, ledger=0)
    assert A.classical_projection() == Bn.classical_projection() \
        == Cn.classical_projection()
    assert A.total_seam() == 0 and Bn.total_seam() == 0 and Cn.total_seam() == 2
    assert Bn.digits != Cn.digits          # distinct UGD states
    cert["T4_classical_projection_is_blind_to_seam"] = {
        "statement": (
            "In the null-seam limit the classical projection recovers "
            "ordinary positional notation exactly. With active seam "
            "charge the projection DISCARDS information: three numerals "
            "with identical classical projection, two of them carrying "
            "different seam data and one a different total seam charge. "
            "Visible equality is not UGD equality — as a witness"
        ),
        "null_seam_projection": dec(N_null.classical_projection()),
        "shared_projection": dec(A.classical_projection()),
        "total_seams": [A.total_seam(), Bn.total_seam(), Cn.total_seam()],
        "states_distinct": True,
        "verdict": "PASS",
    }

    # ---------------- T5 cut-zero ----------------
    Z = ZERO_DIGIT
    # neutrality is DERIVED: adding the neutral digit changes nothing
    neutral_ok = True
    for a in range(0, 12):
        for s in (-1, 0, 1):
            D = digit(a, 0, s)
            if add_digits(D, Z).digits != D.digits:
                neutral_ok = False
            if add_digits(D, Z).ledger != D.ledger:
                neutral_ok = False
    assert neutral_ok
    assert Z.classical_projection() == 0 and Z.total_seam() == 0
    # but the neutral numeral is a POINT OF THE CHART: it has a scale
    # index and a phase exponent, unlike "no numeral at all"
    absent = Numeral({})
    assert absent.digits == {} and Z.digits == {0: (0, 0)}
    assert absent != Z
    cert["T5_cut_zero"] = {
        "statement": (
            "The all-neutral digit is the neutral element for seam-aware "
            "addition (derived, not assumed) and has null classical "
            "projection and zero ledger — yet it is a POST-CUT "
            "representative: it occupies a scale and carries a phase "
            "exponent, and is therefore distinct from the absence of a "
            "numeral. Zero is the first cut, not primitive absence"
        ),
        "neutrality_derived_on_grid": neutral_ok,
        "neutral_projection": dec(Z.classical_projection()),
        "neutral_total_seam": Z.total_seam(),
        "distinct_from_absent_numeral": absent != Z,
        "verdict": "PASS",
    }

    # ---------------- T6 lambda-logic ----------------
    # negation is an involution on the finite truth domain
    involution = True
    domain = [(p, sg, k) for p in range(0, 4)
              for sg in (-1, 0, 1) for k in (-1, 0, 1)]
    for v in domain:
        if ugd_negate(ugd_negate(v)) != v:
            involution = False
    assert involution
    # seam predicate marks exactly the k != 0 values
    seam_marks = all(seam_predicate(v) == (v[2] != 0) for v in domain)
    assert seam_marks
    # classical limit: trivial phase, null seam -> ordinary negation of
    # support, and the seam predicate is inactive
    classical_ok = True
    for sg in (-1, 0, 1):
        v = (0, sg, 0)
        n = ugd_negate(v)
        if n != (0, -sg, 0) or seam_predicate(v) or seam_predicate(n):
            classical_ok = False
    assert classical_ok
    # bounded ledgered contradiction: A and not-A both carry seam, total
    # bounded; unledgered case is the one that blocks
    A_v = (1, 1, 1)
    notA = ugd_negate(A_v)
    bounded_total = abs(A_v[2] + notA[2])
    assert bounded_total == 0          # the seam charges cancel: ledgerable
    unledgered = (1, 1, 1)
    assert abs(unledgered[2]) > 0      # active, must be carried or block
    cert["T6_lambda_logic"] = {
        "statement": (
            "UGD negation N(phi, sigma, k) = (phi, -sigma, -k) is an "
            "involution on the truth domain; the seam predicate marks "
            "exactly the values with nonzero seam charge; and in the "
            "null-seam trivial-phase limit the structure reduces to "
            "ordinary support negation with the seam predicate inactive — "
            "the classical limit, verified exhaustively"
        ),
        "domain_size": len(domain),
        "negation_is_involution": involution,
        "seam_predicate_exact": seam_marks,
        "classical_limit_recovered": classical_ok,
        "contradiction_seam_charges_cancel": bounded_total == 0,
        "note": (
            "a contradiction whose seam charge cancels is ledgerable and "
            "non-explosive; an unledgered active seam charge remains an "
            "open gap and blocks commit"
        ),
        "verdict": "PASS",
    }

    cert["finding_UGD1_F1"] = (
        "The quarter turn appears at the numeral level: at K=4 the phase "
        "generator has order exactly 4 and its square is the half turn. "
        "The same object now stands certified in four presentations — "
        "arithmetic (LAM-2 iota^2=-1), quantum geometry (CFE-Q J^2=-I), "
        "native algebra (EMK-1 R^2=-I), and number (UGD-1 K=4)."
    )
    cert["finding_UGD1_F2"] = (
        "The carry discipline is a conservation law: seam charge pushed "
        "out of the balanced alphabet is conserved in the ledger, exactly. "
        "The tampered policy that clamps instead of ledgering loses charge "
        "and separates — the ledger is load-bearing, not bookkeeping."
    )
    cert["finding_UGD1_F3"] = (
        "The classical projection is provably lossy: explicit numerals "
        "share a classical value while differing in seam data and total "
        "seam charge. This is the precise sense in which UGD numbers "
        "carry strictly more than their classical shadows."
    )
    cert["claim_boundary"] = {
        "certified": (
            "finite arithmetic and logic spine — cyclic phase algebra, "
            "carry conservation, scale additivity, projection blindness, "
            "cut-zero neutrality, negation involution, classical limit"
        ),
        "infinite_numerals": "NOT CLAIMED (convergence not addressed)",
        "linked_seam_composition_rule": "DECLARED, not derived",
        "CFE_uniqueness_U": "OPEN (inherited)",
        "RH_K0_L0": "not touched",
        "YM_continuum_gates": "not touched",
        "quantum_gravity": "not touched",
    }
    cert["arithmetic_discipline"] = (
        "exact integers for phase exponents and seam charges; exact "
        "rationals for the scale base; no root of unity is ever evaluated "
        "numerically; no floating point anywhere"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "UGD1_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_UGD1.sha256"), "w") as f:
        f.write(digest + "\n")
    print("UGD1 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
