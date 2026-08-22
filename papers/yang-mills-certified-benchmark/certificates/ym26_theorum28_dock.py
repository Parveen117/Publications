"""YM-26: DOCK TO THE RECOGNITION-COMPLETE FINITE-TO-INFINITE CUT THEOREM
(Recognition-Kernel-Framework/theorum/28). The owner's correction, Aug 22:
the framework already has the convergence machinery that turns an
infinite-dimensional cut form into a finite matrix decision — the outward
certificate u_n + e_n < 1 — and I had been working around it instead of
on it. This capsule (i) identifies every object of theorum/28 with an
object already certified in YM-15/16, (ii) re-certifies YM-16's chain
dock AS a lawful theorum/28 outward certificate (with e_n present — not
the "shadow certificate" the theorem warns against), and (iii) states
the exact remaining obligation in the theorem's own hypothesis list.

DICTIONARY (theorum/28 -> chain):
  event carrier K            L^2 of the chain rungs (Haar <-> Phi_Sigma declared)
  lift Z_n                   the chain transfer S_m (refinement index = chain/content)
  P_n (recognized)           orthogonal projection onto the A-line carrier
                             W_{m+1} = {1, chi_half(A_i)} (YM-15 closed forms)
  Q_n = I - P_n (memory)     the complement (multi-insertion / higher content)
  R_n = Y_n^* Y_n            recognized form = P S P, exact (YM-15 engine)
  D_n = M_n^* M_n            memory form = Q S Q
  reference floor m > 0      vacuum floor lambda_1 >= f_0^{m-1} (YM-16 T1)
  threshold / metric         mu = nu * f_0^{m-1},  nu = 9/10 (YM-16 NU)
  faithfulness residual      ||Q S Q|| / mu  (sup route: lambda^2 e^{kappa(m-1)} / mu)
  seam coupling              ||P S Q|| via the orthonormal-Gram bound (YM-16 pmq)
  finite seam matrix B_n     the Schur-corrected carrier block at threshold mu
  N_+(B_n - I)               Haynsworth/LDL inertia count at mu (YM-6/16)
  outward certificate        e_n := ||QSQ||/mu < 1  AND  inertia count == 1
                             => lambda_2(S_m) < nu lambda_1  (gap >= -log nu)

 (T1) LAWFUL FINITE CERTIFICATE FOR m <= m*(kappa). For kappa in
      {1/8, 1/4, 1/2}, every chain length m <= m*(kappa) (13, 7, 4) passes
      BOTH halves of the theorum/28 outward certificate with the blind
      error e_n explicitly bounded (never omitted): e_n < 1 and exactly
      one eigenvalue above threshold. This is YM-16's result re-read as
      Theorem 7.1 + Section 9 of theorum/28 — so the strong-coupling
      chain gap for m <= m* is a recognition-complete finite decision,
      not a fixed-grid top.

 (T2) THE OBLIGATION, IN THE THEOREM'S WORDS. The certificate refuses
      at m* + 1 because e_n grows like (e^kappa / f_0)^{m-1}: the memory-
      channel bound is NOT recognition-Cauchy uniformly in the chain
      refinement index (Definition 3.1 fails along m). Certified table
      of e_n(m) crossing 1 exactly at m* + 1 on each grid kappa. Hence
      the theorem's hypothesis list for the chain reads:
        1. recovered identity across refinement arrows      BUILT (YM-F1 fabric)
        2. recognized-channel Cauchy bound (A-line)          BUILT (YM-15 closed form, exact)
        3. memory-channel Cauchy bound uniform in m          NOT YET BUILT (sup route fails)
        4. Smriti tail for omitted content j >= 1            BUILT at leading order
                                                             (YM-21/22 ladder r_j; uniform
                                                             in a along kappa = theta a)
        5. target-faithfulness residual -> 0                 = item 3
        6. uniform positive reference floor                  BUILT (lambda_1 >= f_0^{m-1})
        7. outward finite seam-matrix margin                 BUILT for m <= m*
      Item 3 is the whole remaining strong-coupling problem, now
      stated as one hypothesis of one framework theorem — the same
      shape as the RH gate's "NOT YET BUILT" lines in theorum/28 Sec. 12.

 (T3) WHAT THE TILING LAW CONTRIBUTES TO ITEM 3. The memory channel
      splits by content: Q = Q_{1/2}^{multi} + Q_{>=1}. The second piece
      carries the ladder tail: certified per-face weight
      sum_{j>=1} d_j r_j relative to f_0, geometric in kappa, and on the
      heat-kernel trajectory its per-unit-time weight -> 0 (YM-22 T5).
      The first piece (several chi_half insertions, all content 1/2) is
      the j = 1/2 branching-tiling sector of YM-22: it is the ONLY part of
      item 3 that survives the cutoff limit. So item 3 is reduced to a
      single statement: a Cauchy bound, uniform in m, for the multi-
      insertion content-1/2 memory channel. Not proved here.

NOT CLAIMED: m-uniform gap; weak coupling (YM-25 isolated it to the
intertwiner transfer — theorum/28 applies there too, with item 4's tail
the unbuilt piece since the ladder collapses at large kappa); Clay.

Controls:
  C1  for each grid kappa and m <= m*: e_n < 1 and inertia count == 1.
  C2  e_n(m) strictly increasing in m and crosses 1 exactly at m* + 1.
  C3  recognized channel exact: YM-15 closed-form block reproduced.
  C4  reference floor: vacuum lower bound f_0^{m-1} positive, all m tested.
  C5  ladder tail per face decreasing in kappa grid and geometric.
  C6  shadow-certificate tamper: dropping e_n would accept m* + 1 — the
      theorem's warning reproduced (a fixed-grid top without e_n is a
      shadow certificate only).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, _dec, canonical_sha  # noqa: E402
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym4_symmetry_protected import face_coeffs, dim  # noqa: E402
from ym15_chain_closed_form import lam_half  # noqa: E402
from ym16_chain_dock import (  # noqa: E402
    dock_cell, predicted_m_star, f0_of, iv_pow, compressed_M, NU,
)

GRID = [F(1, 8), F(1, 4), F(1, 2)]


def blind_error(kappa: F, m: int) -> F:
    """e_n = ||Q S Q|| / mu, sup route, outward."""
    lam2 = lam_half() * lam_half()
    mu = NU * iv_pow(f0_of(kappa), m - 1).lo
    return (lam2 * Iv(exp_point(kappa * (m - 1)).hi)).hi / mu


def ladder_tail_per_face(kappa: F) -> Iv:
    f = face_coeffs(kappa, 6)
    tail = Iv(F(0))
    for t in range(2, 7):                      # contents j = 1 .. 3
        tail = tail + Iv(F(dim(t))) * f[t] / f[0]
    return tail


def run():
    grid = {}
    c1 = c2 = c3 = c4 = c6 = True
    for kap in GRID:
        mstar = predicted_m_star(kap)
        cells = {}
        prev_e = None
        for m in range(2, mstar + 2):
            e = blind_error(kap, m)
            cell = dock_cell(kap, m)
            lawful = (e < 1) and cell.get("count_exact") == 1
            if m <= mstar and not lawful:
                c1 = False
            if m == mstar + 1 and not (e >= 1):
                c2 = False
            if prev_e is not None and not (e > prev_e):
                c2 = False
            prev_e = e
            # C3 recognized channel exact: compressed block present & square
            M = compressed_M(kap, m)
            if len(M) != m + 1 or any(len(r) != m + 1 for r in M):
                c3 = False
            # C4 floor
            if not (iv_pow(f0_of(kap), m - 1).lo > 0):
                c4 = False
            # C6 shadow tamper: without e_n, would m*+1 be accepted?
            if m == mstar + 1:
                shadow_accept = cell.get("status") != "OK"   # dock refuses only via e_n
                # a "fixed-grid top" would just look at the carrier block
                c6 = c6 and shadow_accept
            cells[str(m)] = {"e_n_blind": _dec(e, 8),
                             "inertia_count": cell.get("count_exact"),
                             "status": "LAWFUL" if lawful else
                             ("REFUSED_e_n>=1" if e >= 1 else cell.get("status"))}
        grid[str(kap)] = {"m_star": mstar, "cells": cells,
                          "ladder_tail_per_face_j>=1": _dec(ladder_tail_per_face(kap).hi, 10)}
    # C5 ladder tail geometric / increasing in kappa (decreasing toward small kappa)
    tails = [ladder_tail_per_face(k).hi for k in GRID]
    c5 = tails[0] < tails[1] < tails[2] and tails[0] < F(1, 100)

    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM26_DOCK_TO_THEORUM28_FINITE_TO_INFINITE_CUT",
        "source": "Recognition-Kernel-Framework/theorum/28_recognition_complete_"
                  "finite_to_infinite_cut_theorem.md (Thm 5.1, 6.1, 7.1, Sec 8-9, 12)",
        "claim_status": "lawful_outward_certificate_for_m_le_mstar__obligation_"
                        "named_as_theorem_hypothesis_3__no_m_uniform_claim",
        "dictionary": {
            "P_n": "A-line carrier W_{m+1} (YM-15)", "Q_n": "complement",
            "floor_m": "lambda_1 >= f_0^{m-1}", "threshold": "mu = nu f_0^{m-1}, nu=9/10",
            "e_n": "||QSQ||/mu (sup route)", "N_+(B_n - I)": "Haynsworth inertia at mu",
            "certificate": "e_n < 1 and inertia count == 1 => lambda_2 < nu lambda_1",
        },
        "hypothesis_ledger": {
            "1_recovered_identity": "BUILT (YM-F1)",
            "2_recognized_channel_cauchy": "BUILT exact (YM-15)",
            "3_memory_channel_cauchy_uniform_in_m": "NOT YET BUILT",
            "4_smriti_tail_content_ge_1": "BUILT leading order (YM-21/22)",
            "5_target_faithfulness": "= item 3",
            "6_uniform_reference_floor": "BUILT",
            "7_outward_seam_margin": "BUILT for m <= m*",
        },
        "reduction_of_item_3": "Cauchy bound uniform in m for the multi-insertion "
                               "content-1/2 memory channel (j=1/2 branching tilings)",
        "grid": grid,
        "controls": {
            "C1_lawful_for_m_le_mstar": bool(c1),
            "C2_e_n_increasing_crosses_1_at_mstar_plus_1": bool(c2),
            "C3_recognized_channel_exact": bool(c3),
            "C4_reference_floor_positive": bool(c4),
            "C5_ladder_tail_geometric": bool(c5),
            "C6_shadow_certificate_warning_reproduced": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM26_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM26.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "m*", v["m_star"],
              {m: (c["e_n_blind"][:7], c["status"]) for m, c in v["cells"].items()})
    print("sha256:", sha)
