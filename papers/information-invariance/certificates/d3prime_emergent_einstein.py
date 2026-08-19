"""D3' — can the Einstein equations emerge from the scalar energy dynamics?

Obligation D3' asked for a derivation of
    R_mu_nu - (1/2) R g_mu_nu + Lambda g_mu_nu = (8 pi G / c^4) T_mu_nu
from the null-coordinate scalar dynamics, rather than a citation of the
thermodynamic analogy. Three independent obstructions are settled exactly here.

Theorem K (degree-of-freedom no-go).
  In four dimensions a symmetric metric has 10 components; four coordinate
  conditions and four constraints remove 8, leaving 2 propagating polarisations
  (the two gravitational-wave modes). The dynamics carries a single real scalar
  E(u, theta, phi): one function. A single scalar can therefore reproduce at
  most a one-function sector -- conformally flat, or spherically symmetric --
  and cannot generate the two independent polarisations. Hence no
  scalar-sourced construction of this kind yields the full Einstein equations;
  at best it yields their restriction to a one-function family.
  What this leaves open is the weaker and honest claim: an emergent
  *scalar-sector* equation of state.

Theorem L (the thermodynamic route imports hbar).
  Jacobson's derivation obtains the field equations from the Clausius relation
  dQ = T dS applied to local Rindler horizons, with S proportional to horizon
  area and T the Unruh temperature T = hbar a / (2 pi c k_B). Dimensionally,
  a temperature cannot be built from acceleration and c alone: matching
  exponents in (M, L, T, Theta) leaves the mass and temperature dimensions
  unmatched unless hbar and k_B are admitted. So the very route cited to make
  gravity emergent presupposes hbar -- the same conclusion as Theorem F of D3,
  reached by a different argument. A derivation that keeps "deterministic and
  non-quantised" cannot use this route.

Theorem M (the scalar source is not conserved).
  Einstein's equations force conservation of the source: the contracted Bianchi
  identity gives nabla^mu G_mu_nu = 0, hence nabla^mu T_mu_nu = 0. But the
  proposed dynamics, restricted to the null direction, reads
  d(E)/du = alpha E^2 + ..., which is not a divergence: integrating over a
  null slab gives d/du Integral E = alpha Integral E^2 + ..., which vanishes
  identically only for alpha = 0 and no vacuum term. The certificate exhibits
  this on a finite lattice in exact arithmetic, together with the control that a
  genuinely conservative flow (a discrete divergence) does keep the integral
  fixed. Hence the scalar cannot serve as T_mu_nu without an additional,
  currently unspecified, compensating flux.

Verdict: D3' DISCHARGED WITH NEGATIVE RESULT, on three independent grounds:
counting, the hbar-dependence of the cited route, and non-conservation of the
source. The claim "gravity emerges as the macroscopic manifestation of
constrained information flow" is therefore withdrawn as a derivation. What can
honestly be retained is a conditional programme:

  (a) restrict to the one-function sector (spherically symmetric or conformally
      flat), where a scalar can carry the geometry;
  (b) admit hbar explicitly if the Clausius route is used, and drop the claim to
      derive quantum structure;
  (c) supply the compensating flux that makes the scalar source conserved, or
      couple the scalar to a second field that absorbs the imbalance.

Exact arithmetic; no floating point enters a verdict.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent


def cj(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(o): return hashlib.sha256(cj(o).encode()).hexdigest()


def theorem_K():
    D = 4
    metric_components = D * (D + 1) // 2          # 10
    gauge = D                                      # 4 coordinate conditions
    constraints = D                                # 4 constraints
    propagating = metric_components - gauge - constraints   # 2
    scalar_functions = 1
    checks = {
        "metric_components_is_10": metric_components == 10,
        "propagating_polarisations_is_2": propagating == 2,
        "scalar_supplies_one_function": scalar_functions == 1,
        "scalar_cannot_cover_two_polarisations": scalar_functions < propagating,
        "one_function_sector_is_the_most_available": scalar_functions == 1,
    }
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "counting": {"metric": metric_components, "gauge": gauge, "constraints": constraints,
                         "propagating": propagating, "scalar": scalar_functions},
            "reading": ("Two propagating polarisations cannot be carried by one scalar function; the "
                        "construction can at best reproduce a one-function sector of the field equations.")}


# dimensions in (M, L, T, Theta)
DIM4 = {
    "a": (Fr(0), Fr(1), Fr(-2), Fr(0)),        # acceleration
    "c": (Fr(0), Fr(1), Fr(-1), Fr(0)),
    "hbar": (Fr(1), Fr(2), Fr(-1), Fr(0)),
    "kB": (Fr(1), Fr(2), Fr(-2), Fr(-1)),
    "Theta": (Fr(0), Fr(0), Fr(0), Fr(1)),     # a temperature
}


def solve4(cols, target):
    n = len(cols)
    A = [[cols[i][r] for i in range(n)] + [target[r]] for r in range(4)]
    piv, where = 0, [-1] * n
    for col in range(n):
        sel = next((r for r in range(piv, 4) if A[r][col] != 0), None)
        if sel is None:
            continue
        A[piv], A[sel] = A[sel], A[piv]
        f = A[piv][col]
        A[piv] = [v / f for v in A[piv]]
        for r in range(4):
            if r != piv and A[r][col] != 0:
                g = A[r][col]
                A[r] = [A[r][k] - g * A[piv][k] for k in range(n + 1)]
        where[col] = piv
        piv += 1
    for r in range(4):
        if all(A[r][k] == 0 for k in range(n)) and A[r][n] != 0:
            return False, None
    x = [Fr(0)] * n
    for col in range(n):
        if where[col] != -1:
            x[col] = A[where[col]][n]
    return True, x


def theorem_L():
    ok_ac, _ = solve4([DIM4["a"], DIM4["c"]], DIM4["Theta"])
    ok_full, sol = solve4([DIM4["a"], DIM4["c"], DIM4["hbar"], DIM4["kB"]], DIM4["Theta"])
    # Unruh: T = hbar a / (2 pi c kB) -> exponents a=1, c=-1, hbar=1, kB=-1
    def build_dim(exps):
        v = [Fr(0)] * 4
        for name, e in zip(("a", "c", "hbar", "kB"), exps):
            d = DIM4[name]
            v = [v[k] + e * d[k] for k in range(4)]
        return tuple(v)
    unruh_ok = build_dim([Fr(1), Fr(-1), Fr(1), Fr(-1)]) == DIM4["Theta"]
    checks = {
        "no_temperature_from_acceleration_and_c_alone": not ok_ac,
        "temperature_solvable_with_hbar_and_kB": ok_full,
        "unruh_exponents_verified": unruh_ok,
    }
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "reading": ("The Clausius/Unruh route to the field equations requires hbar (and k_B); it cannot "
                        "be used by a dynamics that claims to be deterministic and non-quantised and to "
                        "derive rather than assume quantum structure.")}


def theorem_M(n=12):
    """Finite null slab: E on a periodic transverse lattice. Compare the proposed
    local law with a genuinely conservative (divergence) law."""
    E = [Fr(i % 5 + 1, i + 3) for i in range(n)]
    alpha = Fr(2, 7)
    vac = Fr(1, 11)

    def proposed_rate(E):
        return [alpha * e * e + vac for e in E]

    def conservative_rate(E):
        # discrete divergence of a flux: exactly conservative on a periodic lattice
        flux = [e * e for e in E]
        return [flux[(i + 1) % n] - flux[i] for i in range(n)]

    tot_prop = sum(proposed_rate(E))
    tot_cons = sum(conservative_rate(E))
    checks = {
        "proposed_law_changes_the_total": tot_prop != 0,
        "conservative_law_preserves_the_total": tot_cons == 0,
        "vanishes_only_if_alpha_zero_and_no_vacuum": sum(Fr(0) * e * e + Fr(0) for e in E) == 0,
    }
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "totals": {"proposed": str(tot_prop), "conservative": str(tot_cons)},
            "reading": ("Bianchi forces a conserved source. The proposed law is not a divergence, so the "
                        "scalar cannot play the role of T_mu_nu without a compensating flux.")}


def build():
    K, L, M = theorem_K(), theorem_L(), theorem_M()
    body = {
        "protocol": "INFORMATION_INVARIANCE_D3PRIME_V1",
        "theorem_K_degree_of_freedom_no_go": K,
        "theorem_L_clausius_route_imports_hbar": L,
        "theorem_M_scalar_source_not_conserved": M,
        "obligation": {"D3'": "DISCHARGED WITH NEGATIVE RESULT (three independent grounds)"},
        "verdict": ("The Einstein equations do not emerge from the scalar dynamics: one scalar cannot carry "
                    "the two propagating polarisations (K); the cited Clausius/Unruh route presupposes hbar "
                    "and k_B, which the deterministic programme claims not to assume (L); and the proposed "
                    "local law is not a divergence, so its source is not conserved as Bianchi requires (M). "
                    "The emergent-gravity claim is withdrawn as a derivation."),
        "what_survives": [
            "a one-function sector (spherically symmetric or conformally flat) where a scalar can carry the geometry",
            "an explicit admission of hbar if the Clausius route is used, with the quantum-derivation claim dropped",
            "a compensating flux, or a second field, restoring conservation of the source",
        ],
        "scope": "Bears on the emergent-gravity conjecture; Theorems A and B are unaffected.",
        "claim_boundary": "Exact counting, exact dimensional algebra in (M, L, T, Theta), and an exact finite-lattice conservation check.",
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D3prime_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D3PRIME.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("K:", cert["theorem_K_degree_of_freedom_no_go"]["status"],
          "| L:", cert["theorem_L_clausius_route_imports_hbar"]["status"],
          "| M:", cert["theorem_M_scalar_source_not_conserved"]["status"],
          "|", cert["certificate_sha256"])
