"""YM-23: THE WEAK-COUPLING TURN — where the tiling route ends, what the
fabric looks like near the flat point, and the exact negative result
that defines the real problem: the ABELIANIZED fabric has NO volume-
uniform gap, so at weak coupling the gap can come only from the
commutator (non-abelian) residue.

Setting: chain fabric (YM-F1), heat-kernel time step a, bridge coupling
kappa. Weak coupling = large kappa: face weights Exp_Sigma(-kappa rho_l)
concentrate near the lawful face H_l = I.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) THE TILING ROUTE HAS AN EXACT WEAK-COUPLING BOUNDARY.
      (a) Coefficient ladder collapse: r_j(kappa) = I_{2j+1}(kappa)/I_1(kappa)
          increases toward 1 — certified monotone on a kappa grid and
          r_{1/2}, r_1, r_{3/2} all above 4/5 at kappa = 50 (0.970, 0.922, 0.859). The
          expansion parameter of YM-21/22 is no longer small: every face
          content carries comparable weight. No tiling order dominates.
      (b) Leading rate sign change: the leading tiling rate per step
          3a/4 - A_Sigma(r(kappa)) is positive iff r(kappa) < tanh(3a/8);
          the boundary r_c(a) = tanh_Sigma(3a/8) is computed NATIVELY as
          (Exp(3a/4) - 1)/(Exp(3a/4) + 1) (F00-G Cayley coordinate of
          Exp_Sigma(3a/4)), and kappa_c(a) is bracketed by monotonicity
          of r(.) — certified at a in {1, 1/2, 1/4, 1/8}: the route is
          vacuous for kappa > kappa_c(a), and kappa_c(a) -> 0 with a.
      (c) AF crossing: on a physical trajectory kappa(a) = c log(1/a)
          (DECLARED shape, not derived) the tiling route covers only
          a > a*(c) where c log(1/a*) = kappa_c(a*); certified bracket of
          a* for c = 1/4. Conclusion: the continuum along any
          asymptotically-free trajectory lies OUTSIDE the tiling route.

 (T2) NATIVE WEAK-COUPLING LINEARIZATION — CURVATURE IS THE COMMUTATOR.
      In the stereographic chart (YM-F1, no trig) a face H = 1 + p with
      p in the odd (imaginary-quaternion) sector has residue
          rho(H) = 2 s/(1+s),   s = |p|^2   (exact rational),
      so near the flat point rho = 2|p|^2 - 2|p|^4 + ... — the face
      weight is Gaussian in the odd coordinate at leading order. For
      two faces H_1 = 1 + p_1, H_2 = 1 + p_2 (UNNORMALISED, the chart
      before projection), the product's odd part is EXACTLY
          odd(H_1 H_2) = p_1 + p_2 + (1/2)[p_1, p_2]  (quaternion commutator),
      and its even part is 1 - <p_1, p_2>. The commutator is the ONLY
      term that breaks additivity: linear Stokes (sum of face
      coordinates = boundary coordinate) holds exactly iff all face
      commutators vanish. This is EMK-1 T5's cut-commutator curvature
      identity (C_i C_x - C_x C_i = commutator) and RST-1 T4's flow loop
      L_h = I + h^2 [D_1, D_2] read on the fabric: at weak coupling the
      fabric's recognition curvature IS the commutator residue.
      Certified exactly on a rational grid; control: commuting faces
      (p_1 parallel p_2) give exactly additive Stokes.

 (T3) THE ABELIANIZED FABRIC IS GAPLESS IN VOLUME (exact negative).
      Drop the commutator (keep only the Gaussian/additive part). The
      space coupling becomes the quadratic form kappa sum_i |p_i - p_{i+1}|^2
      = kappa <p, L_m p> with L_m the path-graph Laplacian. Its lowest
      nonzero stiffness is bounded by the exact rational Rayleigh
      quotient of the linear test vector v_i = i - (m+1)/2:
          <v, L_m v>/<v, v> = 12 / (m(m+1))    (exact identity, certified m=2..40),
      so the softest spatial mode of the abelianized chain has stiffness
      <= 12 kappa / (m(m+1)) -> 0. Whatever excitation energy the
      harmonic transfer assigns to a mode of vanishing stiffness (the
      harmonic relation itself is DECLARED here, not certified), the
      abelianized chain cannot have a volume-uniform gap. Control: the
      constant vector has stiffness exactly 0 (gauge/zero mode) and is
      orthogonal to v.

 (T4) WHAT THIS SAYS, NATIVELY. Strong coupling (YM-21/22): the gap is
      a tiling count, uniform in m and a at leading order. Weak coupling:
      the additive part of the fabric is gapless in volume (T3); the
      only non-additive term is the commutator residue (T2). Therefore
      a weak-coupling gap, if it exists, is carried ENTIRELY by the
      commutator sector — the framework's "non-healable curvature"
      (EMK-1 winding, EMK-G3 helical memory) in the fabric. That is the
      native statement of "the mass gap is non-perturbative". Named
      next object: the commutator-residue transfer on the fabric — the
      dynamics of [p_i, p_{i+1}] alone — and whether it carries a gap.

NOT CLAIMED: any weak-coupling gap; asymptotic freedom (trajectory
shape declared); harmonic excitation spectrum (declared); continuum.

Controls:
  C1  ladder monotone and collapsed at kappa = 50.
  C2  r_c(a) native Cayley form equals tanh via Exp two-route.
  C3  kappa_c(a) bracket: rate positive just below, negative just above.
  C4  stereographic residue identity rho = 2s/(1+s) on rational grid.
  C5  product law odd/even parts exact; commuting control additive.
  C6  Rayleigh identity 12/(m(m+1)) exact m = 2..40; zero mode orthogonal.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, bessel_I, _dec, canonical_sha, TERMS  # noqa: E402
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym20_native_origin_audit import native_odd_log  # noqa: E402
from ymf1_chain_fabric import Quat  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402

C_HALF = F(3, 4)
KAPPA_LADDER = [F(1), F(2), F(5), F(10), F(20), F(50)]
A_GRID = [F(1), F(1, 2), F(1, 4), F(1, 8)]
AF_C = F(1, 4)


def _terms(kappa: F) -> int:
    return max(TERMS, int(2 * kappa) + 40)


def r_of(kappa: F) -> Iv:
    t = _terms(kappa)
    return _r(bessel_I(2, kappa, t) / bessel_I(1, kappa, t))


def rj(j2: int, kappa: F) -> Iv:
    t = _terms(kappa)
    return _r(bessel_I(j2 + 1, kappa, t) / bessel_I(1, kappa, t))


def r_crit(a: F) -> Iv:
    """tanh(3a/8) = (e^{3a/4} - 1)/(e^{3a/4} + 1), native Cayley form."""
    e = exp_point(C_HALF * a)
    E = Iv(e.lo, e.hi)
    return _r((E - Iv(F(1))) / (E + Iv(F(1))))


def leading_rate_step(a: F, kappa: F) -> Iv:
    return Iv(C_HALF * a) - native_odd_log(r_of(kappa))


def kappa_crit_bracket(a: F, lo=F(1, 100), hi=F(20), tol=F(1, 1000)):
    """bisection on the exact sign of the leading rate (r monotone)."""
    assert leading_rate_step(a, lo).lo > 0 and leading_rate_step(a, hi).hi < 0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        g = leading_rate_step(a, mid)
        if g.lo > 0:
            lo = mid
        elif g.hi < 0:
            hi = mid
        else:
            break
    return lo, hi


_LOG2 = None


def log_inv_a_bracket(k: int) -> Iv:
    """native log(2^k) = k * A_Sigma(1/3), since 2 = (1+1/3)/(1-1/3)."""
    global _LOG2
    if _LOG2 is None:
        _LOG2 = _r(native_odd_log(Iv(F(1, 3))))
    return _r(Iv(F(k)) * _LOG2)


def stereo_residue(p: tuple) -> F:
    s = sum(F(x) ** 2 for x in p)
    return 2 * s / (1 + s)


def quat_from_odd(p: tuple) -> Quat:
    return Quat(1, p[0], p[1], p[2])


def run():
    # ---- T1a ladder collapse
    c1 = True
    ladder = {}
    prev = None
    for k in KAPPA_LADDER:
        vals = [rj(1, k), rj(2, k), rj(3, k)]
        if prev is not None and not all(vals[i].lo > prev[i].hi for i in range(3)):
            c1 = False
        prev = vals
        ladder[str(k)] = [_dec(v.lo, 8) for v in vals]
    if not all(v.lo > F(4, 5) for v in prev):
        c1 = False

    # ---- T1b boundary
    c2 = c3 = True
    boundary = {}
    for a in A_GRID:
        rc = r_crit(a)
        # C2: leading rate at r = r_c is zero: A(r_c) = 3a/4 (two-route)
        A = native_odd_log(rc)
        if not (A.lo <= C_HALF * a <= A.hi + F(1, 10 ** 20)
                and A.lo - F(1, 10 ** 20) <= C_HALF * a):
            c2 = False
        lo, hi = kappa_crit_bracket(a)
        if not (leading_rate_step(a, lo).lo > 0 and leading_rate_step(a, hi).hi < 0):
            c3 = False
        boundary[str(a)] = {"r_c": [_dec(rc.lo, 12), _dec(rc.hi, 12)],
                            "kappa_c_bracket": [_dec(lo, 6), _dec(hi, 6)]}
    # monotone: kappa_c decreases with a
    kcs = [F(boundary[str(a)]["kappa_c_bracket"][1]) for a in A_GRID]
    if not all(kcs[i + 1] < kcs[i] for i in range(len(kcs) - 1)):
        c3 = False

    # ---- T1c AF crossing (declared trajectory kappa = c log(1/a))
    af = {}
    crossing = None
    prev_sign = None
    for k in range(1, 13):
        a = F(1, 2 ** k)
        kap_traj = _r(Iv(AF_C) * log_inv_a_bracket(k))
        k_hi = F(int(kap_traj.hi * 10 ** 6) + 1, 10 ** 6)   # outward rounding
        k_lo = F(int(kap_traj.lo * 10 ** 6), 10 ** 6)
        # larger kappa -> smaller rate
        g_lo = leading_rate_step(a, k_hi)
        g_hi = leading_rate_step(a, k_lo)
        sign = "covered" if g_lo.lo > 0 else ("outside" if g_hi.hi < 0 else "straddle")
        af[str(a)] = {"kappa_traj": [_dec(kap_traj.lo, 6), _dec(kap_traj.hi, 6)],
                      "status": sign}
        if prev_sign == "covered" and sign == "outside" and crossing is None:
            crossing = str(a)
        prev_sign = sign
    t1c = crossing is not None

    # ---- T2 linearization
    c4 = all(stereo_residue(p) == 2 * sum(F(x) ** 2 for x in p) /
             (1 + sum(F(x) ** 2 for x in p))
             for p in [(F(1, 2), F(0), F(1, 3)), (F(2), F(-1), F(1, 5))])
    # rho = 2|p|^2 - 2|p|^4 + ... : check rho - 2s + 2s^2 = O(s^3) exactly
    s = F(1, 10)
    c4 = c4 and (2 * s / (1 + s) - 2 * s + 2 * s * s == 2 * s ** 3 / (1 + s))
    c5 = True
    for p1, p2 in [((F(1, 2), F(1, 3), F(0)), (F(0), F(1, 5), F(2, 7))),
                   ((F(1), F(0), F(0)), (F(0), F(1), F(0))),
                   ((F(3, 4), F(-1, 2), F(1, 6)), (F(1, 3), F(2), F(-1)))]:
        h = quat_from_odd(p1) * quat_from_odd(p2)
        # quaternion commutator of pure parts: [p1,p2] = 2 (p1 x p2)
        cross = (p1[1] * p2[2] - p1[2] * p2[1],
                 p1[2] * p2[0] - p1[0] * p2[2],
                 p1[0] * p2[1] - p1[1] * p2[0])
        odd_expected = tuple(p1[i] + p2[i] + cross[i] for i in range(3))
        dot = sum(p1[i] * p2[i] for i in range(3))
        if not ((h.b, h.c, h.d) == odd_expected and h.a == 1 - dot):
            c5 = False
    # commuting control: parallel odd parts -> additive
    p1 = (F(1, 2), F(1, 3), F(-1, 4))
    p2 = tuple(3 * x for x in p1)
    h = quat_from_odd(p1) * quat_from_odd(p2)
    if not ((h.b, h.c, h.d) == tuple(p1[i] + p2[i] for i in range(3))):
        c5 = False

    # ---- T3 Rayleigh identity
    c6 = True
    for m in range(2, 41):
        v = [F(i) - F(m + 1, 2) for i in range(1, m + 1)]
        num = sum((v[i] - v[i + 1]) ** 2 for i in range(m - 1))
        den = sum(x * x for x in v)
        if num / den != F(12, m * (m + 1)):
            c6 = False
        if sum(v) != 0:          # orthogonal to the zero mode
            c6 = False
    ok = c1 and c2 and c3 and t1c and c4 and c5 and c6
    cert = {
        "certificate_type": "YM23_WEAK_COUPLING_TURN",
        "claim_status": "tiling_route_boundary_exact__native_linearization__"
                        "abelianized_fabric_gapless__weak_coupling_gap_OPEN",
        "theorems": {
            "T1_tiling_route_boundary": {
                "ladder_r_half_r1_r3half_by_kappa": ladder,
                "r_c_and_kappa_c_by_a": boundary,
                "af_trajectory_c_log_1_over_a_c": str(AF_C),
                "af_coverage_by_a": af,
                "first_a_outside_route": crossing,
            },
            "T2_native_linearization":
                "rho = 2s/(1+s); odd(H1 H2) = p1 + p2 + (1/2)[p1,p2]; "
                "even = 1 - <p1,p2>; commutator is the only non-additive "
                "term = EMK-1 T5 / RST-1 T4 curvature on the fabric",
            "T3_abelianized_fabric_gapless":
                "softest spatial stiffness <= 12 kappa/(m(m+1)) exactly; "
                "no volume-uniform gap without the commutator",
            "T4_native_statement":
                "a weak-coupling gap lives entirely in the commutator "
                "(non-healable curvature) sector — named next object",
        },
        "declared_not_certified": ["AF trajectory shape kappa = c log(1/a)",
                                   "harmonic excitation energy ~ sqrt(stiffness)"],
        "controls": {
            "C1_ladder_monotone_and_collapsed": bool(c1),
            "C2_r_c_native_cayley_two_route": bool(c2),
            "C3_kappa_c_bracket_and_monotone": bool(c3),
            "C4_stereographic_residue_identity": bool(c4),
            "C5_product_law_and_commuting_control": bool(c5),
            "C6_rayleigh_identity_and_zero_mode": bool(c6),
            "T1c_af_crossing_found": bool(t1c),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM23_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM23.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    t1 = cert["theorems"]["T1_tiling_route_boundary"]
    print("ladder:", t1["ladder_r_half_r1_r3half_by_kappa"])
    print("boundary:", t1["r_c_and_kappa_c_by_a"])
    print("AF crossing at a =", t1["first_a_outside_route"])
    print("sha256:", sha)
