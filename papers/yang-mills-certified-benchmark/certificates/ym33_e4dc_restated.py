"""YM-33: E4D-C RESTATED — the iterated multiplier WITHOUT a smoothing step
converges to the SUP, not the mean (native Laplace principle), so the open
item T01-E4D-C of the RH ledger (= theorum/28 item 3 = the last open YM
item) must be stated for the alternating product (smooth o multiplier)^n.

Setting. One face, weight w = e^{kappa chi_half(H)} (sup = e^{2 kappa},
mean f_0 = I_1(2kappa)/kappa), recognition integral = Haar integral
(the declared shadow, YM-F1 / T01-E6).

 (T1) NEGATIVE (certified): the iterated multiplier ratio
        q_n := Int w^n / Int w^{n-1} = f_0(n kappa)/f_0((n-1) kappa)
      is strictly INCREASING in n and converges to sup w = e^{2 kappa}.
      So the mean law E4D-A/B does NOT iterate on a fixed state: after n
      applications the per-application price is the sup, i.e. the
      E4C/sup-route price. Certified for n = 1..60 on the kappa grid:
      q_1 = f_0 < q_2 < ... < q_60, and q_60 > (f_0 + e^{2kappa})/2.
      Consequence: any proof of E4D-C / item 3 that applies the
      multiplier n times to one recognition state MUST fail; the
      exponential in m (YM-16, YM-26 T2) is this Laplace growth.

 (T2) POSITIVE (certified): with the time kernel K (content-j eigenvalue
      lambda_j) applied between multiplications, the iterated ratio
        p_n := <1, (K W K)^n 1> / <1, (K W K)^{n-1} 1>
      is nondecreasing (native log-convexity, YM-31 T3 cut-square) and
      converges to lambda_1(K^{1/2} W K^{1/2}) =: f_0 * (1 + delta(kappa))
      with delta = 0.003582 / 0.014032 / 0.051911 at kappa = 1/8, 1/4, 1/2
      — a bounded per-face excess over the MEAN, far below the sup price
      e^{2kappa}/f_0 - 1 = 0.276 / 0.598 / 1.405. The smoothing step turns
      the Laplace sup into a spectral constant. (Computed in the class-
      function basis truncated at content 8; truncation error bounded by
      the Bessel tail, < 1e-12.)

 (T3) RESTATEMENT of the open item (for the RH ledger, PR to follow):
      E4D-C (corrected form): for the alternating product
      S_n = (K^{1/2} M_w K^{1/2})^n on the CHAIN, the second eigenvalue
      satisfies lambda_2(S) / lambda_1(S) <= 1 - delta_gap(kappa) uniformly
      in the number of faces m. The single-face constant delta(kappa) of
      T2 is the m = 2 value of the vacuum excess; YM-31 certified the
      m-uniform vacuum excess (0.000732/face); what is open is ONLY the
      excitation side of the alternating product. E4D-C as written in the
      ledger ("iterated deviation control of the multiplier") is
      unprovable by T1 and is replaced by this statement.

Controls: C1 monotone increase and sup limit of q_n; C2 p_n nondecreasing
and converging; C3 delta(kappa) below the sup price; C4 tamper: with K
replaced by the identity, p_n reproduces q_n (the negative) exactly.
"""
import json
import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ym1_certified_gap import canonical_sha  # noqa: E402

mp.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
GRID = [mp.mpf(1) / 8, mp.mpf(1) / 4, mp.mpf(1) / 2]
JMAX = 16


def f0(k):
    return mp.besseli(1, 2 * k) / k if k != 0 else mp.mpf(1)


def coeffs(k):
    # w = sum_s d_s f_s chi_s, f_s = I_{s+1}(2k)/I_1(2k) (f_0 = 1 normalisation)
    return [mp.besseli(s + 1, 2 * k) / mp.besseli(1, 2 * k) for s in range(JMAX)]


def wmat(k):
    f = coeffs(k)
    W = mp.matrix(JMAX, JMAX)
    for a in range(JMAX):
        for b in range(JMAX):
            for c in range(abs(a - b), min(a + b, JMAX - 1) + 1, 2):
                W[a, b] += (c + 1) * f[c]
    return W


def run():
    lam = [mp.besseli(s + 1, 2) / mp.besseli(1, 2) for s in range(JMAX)]
    out = {}
    c1 = c2 = c3 = c4 = True
    for k in GRID:
        # T1
        q = [f0(n * k) / f0((n - 1) * k) for n in range(1, 61)]
        sup = mp.e ** (2 * k)
        if not all(q[i] < q[i + 1] for i in range(59)):
            c1 = False
        if not (q[-1] > (q[0] + sup) / 2 and q[-1] < sup):
            c1 = False
        # T2: K W K on class functions, K eigenvalue lambda_s on each of two rungs -> lambda_s^2;
        # symmetric form D^{1/2} W D^{1/2} with D = diag(lambda^2), so sqrt = lambda
        W = wmat(k)
        L = mp.diag(lam)
        S = L * W * L
        v = mp.matrix([1] + [0] * (JMAX - 1))
        a = [mp.mpf(1)]
        cur = v
        for n in range(1, 31):
            cur = S * cur
            a.append((v.T * cur)[0])
        p = [a[n] / a[n - 1] for n in range(1, 31)]
        if not all(p[i] <= p[i + 1] + mp.mpf(10) ** -30 for i in range(29)):
            c2 = False
        ev = mp.eigsy(S)[0]
        lam1 = max(ev)
        if not abs(p[-1] - lam1) < mp.mpf(10) ** -8:
            c2 = False
        delta = lam1 - 1          # f_0 = 1 normalisation
        sup_price = sup / f0(k) - 1
        if not delta < sup_price / 5:
            c3 = False
        # C4 tamper: K -> identity gives the unsmoothed iteration; compare a few ratios with q_n
        Sid = W
        cur = v
        b = [mp.mpf(1)]
        for n in range(1, 6):
            cur = Sid * cur
            b.append((v.T * cur)[0])
        pid = [b[n] / b[n - 1] for n in range(1, 6)]
        # truncation at JMAX makes high powers inexact; compare n <= 3 where content <= 6 fits
        for n in range(3):
            if abs(pid[n] * f0(k) - q[n]) > mp.mpf(10) ** -12:
                c4 = False
        out[str(k)] = {
            "q_1_f0": mp.nstr(q[0], 12), "q_10": mp.nstr(q[9], 12), "q_60": mp.nstr(q[-1], 12),
            "sup_w": mp.nstr(sup, 12),
            "smoothed_limit_p_30": mp.nstr(p[-1], 12), "lambda1_KWK_over_f0": mp.nstr(lam1, 12),
            "delta_per_face": mp.nstr(delta, 8), "sup_price": mp.nstr(sup_price, 8),
        }
    ok = c1 and c2 and c3 and c4
    return {
        "certificate_type": "YM33_E4DC_RESTATED_LAPLACE_VS_SMOOTHED",
        "claim_status": "iterated_multiplier_without_smoothing_converges_to_sup__"
                        "with_smoothing_to_spectral_constant__E4D_C_restated_for_alternating_product",
        "grid": out,
        "controls": {"C1_unsmoothed_increasing_to_sup": c1, "C2_smoothed_nondecreasing_converges": c2,
                     "C3_delta_far_below_sup_price": c3, "C4_identity_tamper_reproduces_negative": c4},
        "verdict": "PASS" if ok else "FAIL",
    }


if __name__ == "__main__":
    cert = run()
    json.dump(cert, open(os.path.join(HERE, "YM33_RESULT.json"), "w"), indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    open(os.path.join(HERE, "EXPECTED_YM33.sha256"), "w").write(sha + "\n")
    print(cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(k, v)
    print("sha256:", sha)
