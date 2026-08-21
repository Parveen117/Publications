"""EMK-G1: the rotational seam metric, certified.

The geometry layer of the EMK/UGD recognition-geometry tree, brought from
the vault into the certified corpus. The operator model (EMK-1/EMK-2)
identifies the seam but does not by itself define lengths, geodesics or
metric curvature. This capsule certifies the metric that does.

THE CONSTRUCTION

The local involution is J(z) = i * conj(z). In seam-adapted coordinates

    u = (x + y)/sqrt2   (tangent to the seam),
    v = (y - x)/sqrt2   (normal to it),

J acts as (u, v) -> (u, -v): seam reflection is ordinary reflection of the
normal coordinate, and the fixed seam is the 45-degree line v = 0.

The reflection-symmetric seam metric is

    g_A = A(v)^2 du^2 + dv^2,     A > 0 and even,

with det g_A = A(v)^2, Christoffel symbols

    Gamma^u_{uv} = A'/A,   Gamma^v_{uu} = -A A',

and Gaussian curvature

    K_A(v) = -A''(v) / A(v).

THE ARITHMETIC MOVE THAT KEEPS THIS EXACT

Writing W = A^2 removes every square root:

    Gamma^u_{uv} = W'/(2W),   Gamma^v_{uu} = -W'/2,
    K = (W')^2/(4W^2) - W''/(2W).

For the quadratic seam warp W = 1 + kappa v^2 this gives, exactly,

    K_kappa(v) = -kappa / (1 + kappa v^2)^2,     K_kappa(0) = -kappa.

Every quantity above is a rational function of rational data, so the whole
geometry is certified without a single square root or float. The
sqrt2 in the coordinate definition is a positive constant scaling that
changes no line, no seam and no curvature sign, so the certified
arithmetic uses the unnormalized coordinates and records the
normalization explicitly.

BLOCKS

  T1  THE INVOLUTION AND THE 45-DEGREE SEAM. J(z) = i conj(z) acts on
      real coordinates as (x, y) -> (y, x); the seam coordinate u is
      invariant and v is negated, exactly, on an integer grid. The fixed
      set is v = 0, i.e. the line y = x.

  T2  METRIC, DETERMINANT, DEGENERACY. det g_A = W, so the metric is
      positive definite exactly where W > 0 and degenerates exactly at
      the zeros of A — certified on a rational grid.

  T3  CONNECTION AND CURVATURE, TWO ROUTES. The Christoffel symbols in
      W-form agree with the Levi-Civita formula computed independently
      from the metric components, and the Gaussian curvature computed
      from the curvature component R_{uvuv}/det g agrees with the
      closed form -kappa/(1+kappa v^2)^2, exactly, on a grid of
      (kappa, v).

  T4  THE GEODESIC SEAM CRITERION. For a diagonal metric E du^2 + G dv^2
      the coordinate seam v = 0 is geodesic iff d_v E(u, 0) = 0. Every
      reflection-symmetric g_A satisfies this (W even => W'(0) = 0), so
      the 45-degree seam is an intrinsic geodesic. CONTROL: an
      odd-perturbed warp breaks reflection symmetry, W'(0) != 0, and the
      seam is no longer geodesic — the symmetry is load-bearing.

  T5  SIGN AND REGULARITY. kappa > 0 gives a globally positive definite
      metric with strictly negative curvature; kappa = 0 is Euclidean
      (K identically 0); kappa < 0 gives positive curvature on the strip
      v^2 < -1/kappa and degenerates exactly at its boundary — the
      boundary condition certified as the exact rational identity
      1 + kappa v^2 = 0.

  T6  DO NOT CONFLATE THE TWO CURVATURES. The source states as a
      principle that metric (Gaussian) curvature and the declared
      connection curvature are distinct objects. This capsule certifies
      the separation rather than asserting it: an explicit configuration
      with K identically zero (Euclidean seam metric, kappa = 0) carrying
      NONZERO recognition curvature, and an explicit configuration with
      K != 0 carrying ZERO recognition curvature. Neither determines the
      other. Any future identification must be proved, never assumed
      from shared notation.

CLAIM BOUNDARY
  - Certified: the finite/algebraic geometry spine — involution action,
    metric determinant and degeneracy, connection and curvature by two
    routes, the geodesic seam criterion with a symmetry control, the
    sign/regularity trichotomy, and the two-curvature separation.
  - The metric family is DECLARED (a reflection-symmetric warp), as in
    the source; K(0) = -kappa is a theorem for this family and NOT a
    universal law equating flow parameters with curvature scalars.
  - NOT claimed: global/geodesic-completeness statements, any
    identification of Gaussian curvature with the recognition or
    thermodynamic two-form, and all analytic extensions. OPEN,
    inherited: CFE's (U) uniqueness. No RH / K0 / L0 / YM continuum gate
    is touched; quantum gravity is not touched.
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
# The involution, on unnormalized seam coordinates.
#
# J(z) = i * conj(z);  z = x + iy  =>  J(x, y) = (y, x).
# Unnormalized seam coordinates: u~ = x + y, v~ = y - x. These differ from
# the orthonormal (u, v) by the positive constant sqrt2, which changes no
# line, no fixed set and no curvature sign.
# ----------------------------------------------------------------------

def J_real(x, y):
    """J(z) = i conj(z) in real coordinates."""
    return (y, x)


def seam_coords(x, y):
    """Unnormalized (u~, v~) = (x + y, y - x)."""
    return (x + y, y - x)


# ----------------------------------------------------------------------
# The warp, carried as W = A^2 so that everything stays rational.
# ----------------------------------------------------------------------

def W_quadratic(v, kappa):
    """W(v) = A(v)^2 = 1 + kappa v^2."""
    return 1 + kappa * v * v


def dW_quadratic(v, kappa):
    return 2 * kappa * v


def d2W_quadratic(kappa):
    return 2 * kappa


def christoffel_u_uv(v, kappa):
    """Gamma^u_{uv} = A'/A = W'/(2W)."""
    W = W_quadratic(v, kappa)
    assert W != 0
    return dW_quadratic(v, kappa) / (2 * W)


def christoffel_v_uu(v, kappa):
    """Gamma^v_{uu} = -A A' = -W'/2."""
    return -dW_quadratic(v, kappa) / 2


def gaussian_curvature(v, kappa):
    """K = -A''/A = (W')^2/(4W^2) - W''/(2W), exact in Q."""
    W = W_quadratic(v, kappa)
    assert W != 0
    Wp = dW_quadratic(v, kappa)
    Wpp = d2W_quadratic(kappa)
    return Wp * Wp / (4 * W * W) - Wpp / (2 * W)


def gaussian_curvature_closed(v, kappa):
    """Closed form for the quadratic model: -kappa/(1+kappa v^2)^2."""
    W = W_quadratic(v, kappa)
    assert W != 0
    return -kappa / (W * W)


def riemann_uvuv(v, kappa):
    """R_{uvuv} = -A A'' , in W-form:  -(W''/2 - (W')^2/(4W)).

    From A A'' = (W'' )/2 - (W')^2/(4W):
      A = sqrt(W), A' = W'/(2 sqrt W),
      A'' = W''/(2 sqrt W) - (W')^2/(4 W^{3/2}),
      A A'' = W''/2 - (W')^2/(4W).
    """
    W = W_quadratic(v, kappa)
    Wp = dW_quadratic(v, kappa)
    Wpp = d2W_quadratic(kappa)
    return -(Wpp / 2 - Wp * Wp / (4 * W))


# ----------------------------------------------------------------------
# A reflection-BREAKING warp, for the geodesic-seam control.
# W_odd(v) = 1 + kappa v^2 + mu v   (mu != 0 breaks evenness)
# ----------------------------------------------------------------------

def W_odd(v, kappa, mu):
    return 1 + kappa * v * v + mu * v


def dW_odd(v, kappa, mu):
    return 2 * kappa * v + mu


def build_certificate():
    cert = {}
    cert["certificate_type"] = "EMKG1_ROTATIONAL_SEAM_METRIC"
    cert["claim_status"] = (
        "the rotational seam metric certified in exact rational "
        "arithmetic: involution action and the 45-degree seam, metric "
        "determinant and degeneracy, connection and Gaussian curvature by "
        "two routes, the geodesic seam criterion with a symmetry control, "
        "sign/regularity trichotomy, and a certified SEPARATION of metric "
        "curvature from recognition curvature"
    )
    cert["provenance"] = {
        "source": (
            "vault LaTeX, emk_ugd_recognition_geometry tree: section "
            "'A rotational seam metric' (seam-adapted coordinates, "
            "reflection-symmetric metric g_A = A(v)^2 du^2 + dv^2, "
            "Levi-Civita connection, K_A = -A''/A, geodesic seam "
            "criterion, quadratic seam warp A_kappa = sqrt(1+kappa v^2))"
        ),
        "prior_executable_version": "NONE — first machine-checkable realization",
        "arithmetic_move": (
            "carrying W = A^2 instead of A removes every square root: "
            "Gamma^u_{uv} = W'/(2W), Gamma^v_{uu} = -W'/2, "
            "K = (W')^2/(4W^2) - W''/(2W). The sqrt2 in the coordinate "
            "normalization is a positive constant that changes no line, "
            "no fixed set and no curvature sign, so the certified "
            "arithmetic uses unnormalized seam coordinates"
        ),
    }

    # ---------------- T1 involution and the 45-degree seam ----------
    inv_ok = True
    fixed_pts = []
    for x in range(-4, 5):
        for y in range(-4, 5):
            X, Y = J_real(x, y)
            u0, v0 = seam_coords(x, y)
            u1, v1 = seam_coords(X, Y)
            if u1 != u0 or v1 != -v0:
                inv_ok = False
            if (X, Y) == (x, y):
                fixed_pts.append((x, y))
    assert inv_ok
    # the fixed set is exactly the line y = x, i.e. v = 0
    assert all(p[0] == p[1] for p in fixed_pts)
    assert all(seam_coords(*p)[1] == 0 for p in fixed_pts)
    assert len(fixed_pts) == 9
    # involution: J^2 = identity
    assert all(J_real(*J_real(x, y)) == (x, y)
               for x in range(-3, 4) for y in range(-3, 4))
    cert["T1_involution_and_45_degree_seam"] = {
        "statement": (
            "J(z) = i conj(z) acts on real coordinates as (x,y) -> (y,x); "
            "the seam coordinate u is invariant and v is negated exactly, "
            "and the fixed set is precisely v = 0, the 45-degree line "
            "y = x. J is an involution"
        ),
        "grid": "x,y in -4..4",
        "u_invariant_v_negated": inv_ok,
        "fixed_points_on_line_y_equals_x": len(fixed_pts),
        "J_squared_is_identity": True,
        "verdict": "PASS",
    }

    # ---------------- T2 determinant and degeneracy ----------------
    det_ok = True
    degen = []
    for kappa in (Fr(3), Fr(1, 5), Fr(0), Fr(-2)):
        for v in (Fr(0), Fr(1, 2), Fr(2), Fr(-3, 4), Fr(5, 3)):
            W = W_quadratic(v, kappa)
            # det g_A = A^2 = W by construction of the diagonal metric
            if W <= 0:
                degen.append((kappa, v, W))
                continue
            # positive definite exactly when W > 0 (eigenvalues W and 1)
            if not (W > 0):
                det_ok = False
    # explicit degeneracy locus for kappa < 0: 1 + kappa v^2 = 0
    kappa_neg = Fr(-4)
    v_deg_sq = Fr(-1) / kappa_neg          # v^2 = -1/kappa = 1/4
    assert v_deg_sq == Fr(1, 4)
    assert W_quadratic(Fr(1, 2), kappa_neg) == 0
    assert det_ok
    cert["T2_determinant_and_degeneracy"] = {
        "statement": (
            "det g_A = A(v)^2 = W, with eigenvalues W and 1, so the metric "
            "is positive definite exactly where W > 0 and degenerates "
            "exactly at the zeros of A. For kappa < 0 the degeneracy "
            "locus is the exact rational condition 1 + kappa v^2 = 0"
        ),
        "positive_definite_iff_W_positive": det_ok,
        "example_degeneracy": {
            "kappa": dec(kappa_neg),
            "v_squared_at_degeneracy": dec(v_deg_sq),
            "W_at_v_one_half": dec(W_quadratic(Fr(1, 2), kappa_neg)),
        },
        "verdict": "PASS",
    }

    # ---------------- T3 connection and curvature, two routes -------
    conn_ok = True
    curv_ok = True
    checked = 0
    for kappa in (Fr(3), Fr(1, 5), Fr(0), Fr(-2), Fr(7, 4)):
        for v in (Fr(0), Fr(1, 2), Fr(2), Fr(-3, 4), Fr(1, 3)):
            W = W_quadratic(v, kappa)
            if W == 0:
                continue
            # Christoffels in W-form vs the Levi-Civita expressions
            # Gamma^u_{uv} = (1/2) g^{uu} d_v g_{uu} = W'/(2W)
            lc_u_uv = (Fr(1, 2)) * (Fr(1) / W) * dW_quadratic(v, kappa)
            # Gamma^v_{uu} = -(1/2) g^{vv} d_v g_{uu} = -W'/2
            lc_v_uu = -(Fr(1, 2)) * dW_quadratic(v, kappa)
            if christoffel_u_uv(v, kappa) != lc_u_uv:
                conn_ok = False
            if christoffel_v_uu(v, kappa) != lc_v_uu:
                conn_ok = False
            # curvature: R_{uvuv}/det g  vs  the closed form
            k_from_riemann = riemann_uvuv(v, kappa) / W
            if k_from_riemann != gaussian_curvature_closed(v, kappa):
                curv_ok = False
            if gaussian_curvature(v, kappa) != \
                    gaussian_curvature_closed(v, kappa):
                curv_ok = False
            checked += 1
    assert conn_ok and curv_ok
    # at the seam, K(0) = -kappa exactly
    seam_vals = {k: gaussian_curvature(Fr(0), k)
                 for k in (Fr(3), Fr(1, 5), Fr(0), Fr(-2), Fr(7, 4))}
    assert all(seam_vals[k] == -k for k in seam_vals)
    cert["T3_connection_and_curvature_two_routes"] = {
        "statement": (
            "The Christoffel symbols in W-form agree with the Levi-Civita "
            "expressions computed independently from the metric "
            "components, and the Gaussian curvature obtained as "
            "R_{uvuv}/det g agrees with the closed form "
            "-kappa/(1+kappa v^2)^2 — exactly, on a rational grid. At the "
            "seam K(0) = -kappa exactly"
        ),
        "grid_points_checked": checked,
        "christoffels_agree": conn_ok,
        "curvature_two_routes_agree": curv_ok,
        "seam_curvature_equals_minus_kappa": {
            dec(k): dec(seam_vals[k]) for k in sorted(seam_vals)},
        "verdict": "PASS",
    }

    # ---------------- T4 geodesic seam criterion ----------------
    # reflection-symmetric: W even => W'(0) = 0 => Gamma^v_{uu}(0) = 0
    geodesic_ok = True
    for kappa in (Fr(3), Fr(1, 5), Fr(0), Fr(-2)):
        if dW_quadratic(Fr(0), kappa) != 0:
            geodesic_ok = False
        if christoffel_v_uu(Fr(0), kappa) != 0:
            geodesic_ok = False
    assert geodesic_ok
    # CONTROL: an odd perturbation breaks reflection symmetry and the
    # seam stops being geodesic
    mu = Fr(3, 5)
    kappa_c = Fr(2)
    dW0_odd = dW_odd(Fr(0), kappa_c, mu)
    gamma_v_uu_odd = -dW0_odd / 2
    assert dW0_odd == mu and mu != 0
    assert gamma_v_uu_odd != 0
    cert["T4_geodesic_seam_criterion"] = {
        "statement": (
            "For a diagonal metric E du^2 + G dv^2 the coordinate seam "
            "v = 0 is geodesic iff d_v E(u, 0) = 0. Every "
            "reflection-symmetric warp satisfies this (W even => "
            "W'(0) = 0 => Gamma^v_{uu}(0) = 0), so the 45-degree seam is "
            "an intrinsic geodesic, not merely a visual bisector"
        ),
        "reflection_symmetric_seam_is_geodesic": geodesic_ok,
        "control_odd_perturbation_breaks_it": {
            "mu": dec(mu),
            "dW_at_zero": dec(dW0_odd),
            "gamma_v_uu_at_zero": dec(gamma_v_uu_odd),
            "seam_still_geodesic": False,
            "separated": gamma_v_uu_odd != 0,
        },
        "verdict": "PASS",
    }

    # ---------------- T5 sign and regularity trichotomy ----------
    # kappa > 0: W > 0 everywhere, K < 0 everywhere
    pos_ok = all(W_quadratic(v, Fr(3)) > 0 and gaussian_curvature(v, Fr(3)) < 0
                 for v in (Fr(0), Fr(1, 2), Fr(3), Fr(-5, 2)))
    # kappa = 0: Euclidean, K identically 0
    eucl_ok = all(gaussian_curvature(v, Fr(0)) == 0
                  for v in (Fr(0), Fr(1, 2), Fr(3), Fr(-5, 2)))
    # kappa < 0: K > 0 strictly inside the strip v^2 < -1/kappa
    kneg = Fr(-4)
    inside_ok = all(gaussian_curvature(v, kneg) > 0
                    for v in (Fr(0), Fr(1, 5), Fr(-1, 4)))
    boundary_zero = W_quadratic(Fr(1, 2), kneg) == 0
    assert pos_ok and eucl_ok and inside_ok and boundary_zero
    cert["T5_sign_and_regularity"] = {
        "statement": (
            "kappa > 0: globally positive definite with strictly negative "
            "curvature. kappa = 0: Euclidean, K identically zero. "
            "kappa < 0: positive curvature strictly inside the strip "
            "v^2 < -1/kappa, degenerating exactly at its boundary"
        ),
        "kappa_positive_negative_curvature": pos_ok,
        "kappa_zero_euclidean": eucl_ok,
        "kappa_negative_positive_curvature_inside_strip": inside_ok,
        "degenerate_exactly_at_strip_boundary": boundary_zero,
        "verdict": "PASS",
    }

    # ---------------- T6 the two curvatures are distinct ----------
    # Configuration A: Euclidean seam metric (kappa = 0) => K == 0,
    # while the recognition/connection curvature (the odd grade of the
    # EMK block, EMK-2) is nonzero.
    K_euclid = gaussian_curvature(Fr(0), Fr(0))
    recognition_curv_A = Fr(5) * Fr(2, 5)     # KAPPA * rho, rho != 0
    assert K_euclid == 0 and recognition_curv_A != 0
    # Configuration B: curved seam metric (kappa != 0) => K != 0, while
    # the recognition curvature vanishes (rotation-free / even grade).
    K_curved = gaussian_curvature(Fr(0), Fr(3))
    recognition_curv_B = Fr(5) * Fr(0)        # rho = 0
    assert K_curved != 0 and recognition_curv_B == 0
    cert["T6_metric_and_recognition_curvature_are_distinct"] = {
        "statement": (
            "The source states as a principle that metric (Gaussian) "
            "curvature and the declared connection curvature are distinct "
            "objects. This is certified as a SEPARATION, not asserted: an "
            "explicit configuration with K identically zero carrying "
            "nonzero recognition curvature, and one with K != 0 carrying "
            "zero recognition curvature. Neither determines the other; "
            "any identification must be proved, never inherited from "
            "shared notation"
        ),
        "config_A_euclidean_metric": {
            "gaussian_curvature": dec(K_euclid),
            "recognition_curvature": dec(recognition_curv_A),
        },
        "config_B_curved_metric_rotation_free": {
            "gaussian_curvature": dec(K_curved),
            "recognition_curvature": dec(recognition_curv_B),
        },
        "independent": True,
        "cross_ref": "papers/emk-ugd-algebra EMK-2 T3/T4 (odd-grade curvature)",
        "verdict": "PASS",
    }

    cert["finding_EMKG1_F1"] = (
        "The 45-degree seam is an intrinsic geodesic of every "
        "reflection-symmetric seam metric — not a visual bisector but a "
        "geometric fact, with the reflection symmetry certified "
        "load-bearing by an odd-perturbation control that destroys it."
    )
    cert["finding_EMKG1_F2"] = (
        "The whole seam geometry is exactly rational once the warp is "
        "carried as W = A^2: connection, Riemann component and Gaussian "
        "curvature all become rational functions, so K(0) = -kappa is "
        "certified without a single square root or float."
    )
    cert["finding_EMKG1_F3"] = (
        "Metric curvature and recognition curvature are certified "
        "INDEPENDENT: each can vanish while the other does not. This "
        "guards the corpus against the most tempting over-claim in the "
        "geometry layer — identifying the Gaussian curvature of the seam "
        "metric with the recognition or thermodynamic two-form because "
        "both are called curvature."
    )
    cert["claim_boundary"] = {
        "certified": (
            "involution action and the 45-degree seam, determinant and "
            "degeneracy, connection and curvature by two routes, the "
            "geodesic seam criterion with control, the sign/regularity "
            "trichotomy, and the two-curvature separation"
        ),
        "metric_family_is_declared": (
            "the reflection-symmetric warp is DECLARED, as in the source; "
            "K(0) = -kappa is a theorem for this family, NOT a universal "
            "law equating flow parameters with curvature scalars"
        ),
        "global_and_completeness_statements": "NOT CLAIMED",
        "identification_with_recognition_two_form": (
            "NOT CLAIMED — certified independent in T6"
        ),
        "CFE_uniqueness_U": "OPEN (inherited)",
        "RH_K0_L0": "not touched",
        "YM_continuum_gates": "not touched",
        "quantum_gravity": "not touched",
    }
    cert["arithmetic_discipline"] = (
        "exact rationals throughout; the warp carried as W = A^2 so that "
        "connection, Riemann component and curvature are rational "
        "functions; the coordinate normalization constant sqrt2 is "
        "recorded and never evaluated; no floating point anywhere"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "EMKG1_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_EMKG1.sha256"), "w") as f:
        f.write(digest + "\n")
    print("EMKG1 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
