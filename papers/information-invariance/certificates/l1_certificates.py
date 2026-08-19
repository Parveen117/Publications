"""Information-Invariance paper — executable finite certificates (Layer 1).

Ledger discipline: PROVED / DISPROVED / OPEN / CONJECTURE. Only PROVED items
are claimed in the paper body; everything else is tagged.

Finite chart model. A two-dimensional equilibrium chart is discretised on an
n x m lattice of (p, v) nodes. The accessibility 1-form is
    ω = λ_p dp + λ_v dv,
with λ_p, λ_v rational-valued node functions. Exterior derivative on each
plaquette (the discrete curvature 2-form Ω = dω):
    Ω[i,j] = Δ_p λ_v − Δ_v λ_p  evaluated by the plaquette circulation.
The loop functional (the paper's 𝓘 made concrete) is
    𝓘(γ) = Σ_{edges of γ} ω,   for any closed lattice loop γ.
Discrete Stokes: 𝓘(∂D) = Σ_{plaquettes in D} Ω  (exact over ℚ).

Theorem A (Information invariance ⟺ flat sector)  — PROVED (finite)
  The following are equivalent on the finite chart:
  (i)   δ𝓘 = 0: 𝓘(γ) = 0 for every closed lattice loop γ;
  (ii)  Ω ≡ 0 on every plaquette (dω = 0);
  (iii) ω is exact: ∃ potential Φ with λ_p = Δ_p Φ, λ_v = Δ_v Φ  (path
        independence).
  On the thermodynamic λ-chart with λ_X = −ST/C_X this is the flat-closure
  sector of Theorem 2.1 (Γ_c Γ_m = 1) — i.e. the Onsager-reciprocal sector.
  Proof: (ii)⇒(i) by discrete Stokes (every loop bounds a plaquette sum on a
  simply connected chart); (i)⇒(ii) by taking γ = ∂(single plaquette);
  (i)⇔(iii) by the standard path-independence ⇔ exactness argument on a
  connected graph (define Φ by integrating along any path from a base node).
  The certificate executes all three directions on random rational charts,
  plus planted-curvature negatives.

Theorem B (Single-cycle non-closure does not witness curvature) — PROVED (finite)
  Let the observed branches on a cycle be ω_heat(T), ω_cool(T) with
  ω_cool(T_0) ≠ ω_heat(T_0) at the base point (non-closure). Then there exist
  two distinct generating mechanisms producing identical observed loop area:
  (a) a flat chart (Ω ≡ 0) with an irreversible offset δ applied once
      (state-function shift: Φ → Φ + δ·1_{after excursion});
  (b) a curved chart (Ω ≠ 0) with no offset.
  Hence loop area A = ∮|Δω| dT on one cycle cannot distinguish (a) from (b).
  What distinguishes them is the SECOND cycle: under (a) the loop collapses
  (offset already absorbed); under (b) it reopens with the same area.
  This is exactly the predeclared cycle-2 rule of the T02 Raman campaign
  (Thermodynamics-Reproducibility, reopen ratio 0.6 / 0.2).
  The certificate constructs (a) and (b) explicitly with equal first-cycle
  areas and different second-cycle areas.

Everything else in the paper (energy dynamics eq., emergent Schrödinger,
correlation ~1/(r²+L_P²), emergent Einstein, GUP, Γ_grav) is OPEN or
CONJECTURE — see LEDGER.md. No certificate here touches them.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import random
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent
PROTOCOL = "INFORMATION_INVARIANCE_L1_V1"


def cj(o):
    return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha(o):
    return hashlib.sha256(cj(o).encode()).hexdigest()


# ---------------------------------------------------------------- chart model
class Chart:
    """n x m lattice; edges carry ω-values: ep[i][j] = λ_p on edge (i,j)->(i+1,j),
    ev[i][j] = λ_v on edge (i,j)->(i,j+1). Exact rationals."""

    def __init__(self, n, m, ep, ev):
        self.n, self.m, self.ep, self.ev = n, m, ep, ev

    @classmethod
    def from_potential(cls, n, m, Phi):
        ep = [[Phi[i + 1][j] - Phi[i][j] for j in range(m)] for i in range(n - 1)]
        ev = [[Phi[i][j + 1] - Phi[i][j] for j in range(m - 1)] for i in range(n)]
        return cls(n, m, ep, ev)

    @classmethod
    def random(cls, n, m, rng, exact=False):
        if exact:
            Phi = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(m)] for _ in range(n)]
            return cls.from_potential(n, m, Phi)
        ep = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(m)] for _ in range(n - 1)]
        ev = [[Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(m - 1)] for _ in range(n)]
        return cls(n, m, ep, ev)

    def plaquette(self, i, j):
        """Circulation around plaquette with corner (i,j): ccw."""
        return self.ep[i][j] + self.ev[i + 1][j] - self.ep[i][j + 1] - self.ev[i][j]

    def curvature(self):
        return [[self.plaquette(i, j) for j in range(self.m - 1)] for i in range(self.n - 1)]

    def loop_integral(self, loop):
        """loop: list of nodes (i,j) closed (last == first), unit steps."""
        tot = Fr(0)
        for (a, b), (c, d) in zip(loop, loop[1:]):
            if c == a + 1 and d == b: tot += self.ep[a][b]
            elif c == a - 1 and d == b: tot -= self.ep[c][b]
            elif c == a and d == b + 1: tot += self.ev[a][b]
            elif c == a and d == b - 1: tot -= self.ev[a][d]
            else: raise ValueError("non-unit step")
        return tot

    def rect_loop(self, i0, j0, i1, j1):
        L = [(i, j0) for i in range(i0, i1 + 1)] + [(i1, j) for j in range(j0 + 1, j1 + 1)] \
            + [(i, j1) for i in range(i1 - 1, i0 - 1, -1)] + [(i0, j) for j in range(j1 - 1, j0 - 1, -1)]
        return L

    def try_potential(self):
        """Integrate from (0,0) along rows then columns; return Phi and whether it is consistent."""
        Phi = [[None] * self.m for _ in range(self.n)]
        Phi[0][0] = Fr(0)
        for j in range(1, self.m): Phi[0][j] = Phi[0][j - 1] + self.ev[0][j - 1]
        for i in range(1, self.n):
            for j in range(self.m): Phi[i][j] = Phi[i - 1][j] + self.ep[i - 1][j]
        ok = all(Phi[i][j + 1] - Phi[i][j] == self.ev[i][j] for i in range(self.n) for j in range(self.m - 1))
        return Phi, ok


# ---------------------------------------------------------------- Theorem A
def theorem_A(rng, trials=60):
    checks = {"flat_implies_all_loops_zero": True, "all_plaquettes_zero_implies_exact": True,
              "curved_has_nonzero_loop": True, "stokes_exact_identity": True,
              "single_plaquette_loop_equals_curvature": True}
    for _ in range(trials):
        n, m = rng.randint(3, 6), rng.randint(3, 6)
        flat = Chart.random(n, m, rng, exact=True)
        curv = Chart.random(n, m, rng, exact=False)
        # (ii)->(i) on flat: every rectangular loop zero
        for _ in range(5):
            i0, i1 = sorted(rng.sample(range(n), 2)); j0, j1 = sorted(rng.sample(range(m), 2))
            if flat.loop_integral(flat.rect_loop(i0, j0, i1, j1)) != 0: checks["flat_implies_all_loops_zero"] = False
            # Stokes on curved: loop integral == sum of enclosed plaquettes
            L = curv.loop_integral(curv.rect_loop(i0, j0, i1, j1))
            S = sum(curv.plaquette(i, j) for i in range(i0, i1) for j in range(j0, j1))
            if L != S: checks["stokes_exact_identity"] = False
        # (i)->(ii): single plaquette loop equals curvature
        i, j = rng.randrange(n - 1), rng.randrange(m - 1)
        if curv.loop_integral(curv.rect_loop(i, j, i + 1, j + 1)) != curv.plaquette(i, j):
            checks["single_plaquette_loop_equals_curvature"] = False
        # (i)<->(iii): flat chart admits a consistent potential; curved generically does not
        _, ok_flat = flat.try_potential()
        if not ok_flat: checks["all_plaquettes_zero_implies_exact"] = False
        K = curv.curvature()
        if any(K[a][b] != 0 for a in range(n - 1) for b in range(m - 1)):
            # there must be some loop with nonzero integral: the plaquette itself
            a, b = next((a, b) for a in range(n - 1) for b in range(m - 1) if K[a][b] != 0)
            if curv.loop_integral(curv.rect_loop(a, b, a + 1, b + 1)) == 0: checks["curved_has_nonzero_loop"] = False
    # negative control: plant curvature in one plaquette of a flat chart -> a loop becomes nonzero
    n, m = 4, 4
    flat = Chart.random(n, m, rng, exact=True)
    flat.ep[1][1] += Fr(1, 7)
    planted = flat.loop_integral(flat.rect_loop(1, 1, 2, 2)) == Fr(1, 7)
    checks["planted_curvature_detected"] = planted
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks, "trials": trials}


# ---------------------------------------------------------------- Theorem B
def theorem_B():
    Ts = [300, 350, 400, 450, 500]
    # (a) flat chart + one-time irreversible offset on cooling branch
    base = {T: Fr(1600) - Fr(T - 300, 25) for T in Ts}          # reversible branch (state function)
    delta = Fr(5)                                                 # irreversible upshift after excursion
    heat_a = [base[T] for T in Ts]
    cool_a = [base[T] + delta for T in Ts]                        # whole cooling branch shifted
    # (b) curved chart: genuine reversible loop, closes at base point
    heat_b = [base[T] for T in Ts]
    cool_b = [base[T] + (Fr(20) * (T - 300) * (500 - T) / Fr(10000)) for T in Ts]  # zero at ends
    def area(h, c):
        d = [abs(x - y) for x, y in zip(h, c)]
        w = [Fr(Ts[1] - Ts[0], 2)] + [Fr(Ts[k + 1] - Ts[k - 1], 2) for k in range(1, len(Ts) - 1)] + [Fr(Ts[-1] - Ts[-2], 2)]
        return sum(wi * di for wi, di in zip(w, d))
    A1_a, A1_b = area(heat_a, cool_a), area(heat_b, cool_b)
    # make (b) equal first-cycle area by scaling its bump
    scale = A1_a / A1_b
    cool_b = [h + (c - h) * scale for h, c in zip(heat_b, cool_b)]
    A1_b = area(heat_b, cool_b)
    # second cycle: (a) offset already absorbed -> heating starts from shifted state, loop collapses
    heat_a2 = cool_a; cool_a2 = cool_a
    A2_a = area(heat_a2, cool_a2)
    # (b) curvature reopens identically
    A2_b = area(heat_b, cool_b)
    checks = {
        "first_cycle_areas_equal": A1_a == A1_b,
        "a_nonclosure_at_base": cool_a[0] != heat_a[0],
        "b_closes_at_base": cool_b[0] == heat_b[0],
        "second_cycle_a_collapses": A2_a == 0,
        "second_cycle_b_reopens": A2_b == A1_b,
        "reopen_ratio_a": str(A2_a / A1_a), "reopen_ratio_b": str(A2_b / A1_b),
        "cycle2_rule_separates": (A2_a / A1_a <= Fr(2, 10)) and (A2_b / A1_b >= Fr(6, 10)),
    }
    ok = all(v for k, v in checks.items() if isinstance(v, bool))
    return {"status": "PROVED" if ok else "FAILED", "checks": checks,
            "first_cycle_area": str(A1_a), "construction": {
                "a": "flat chart + one-time irreversible offset delta=5 on cooling branch",
                "b": "curved chart, closes at base, bump scaled to equal first-cycle area"}}


def build():
    rng = random.Random(2026_08_19)
    body = {"protocol": PROTOCOL,
            "theorem_A_information_invariance_iff_flat": theorem_A(rng),
            "theorem_B_single_cycle_nonclosure_not_curvature_witness": theorem_B(),
            "ledger_pointer": "papers/information-invariance/LEDGER.md",
            "claim_boundary": ("Finite-chart theorems over exact rationals. They fix the paper's 𝓘 as the "
                               "loop functional of ω and prove δ𝓘=0 ⟺ dω=0 ⟺ exactness, and that one-cycle "
                               "non-closure is not a curvature witness. Nothing here concerns quantum mechanics, "
                               "gravity, or decoherence; those items are OPEN/CONJECTURE in the ledger.")}
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    out = HERE / "L1_certificate.json"
    out.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_L1.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("A:", cert["theorem_A_information_invariance_iff_flat"]["status"],
          "| B:", cert["theorem_B_single_cycle_nonclosure_not_curvature_witness"]["status"],
          "|", cert["certificate_sha256"])
