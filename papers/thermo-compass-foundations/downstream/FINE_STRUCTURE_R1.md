# RKF constants programme — R1: the fine-structure coupling target

Research note for Monty Dabas · 29 September 2026

**Outcome:** an explicit conditional reduction from the existing native phase and cut-square machinery to the coefficient that an electromagnetic adapter would identify with the fine-structure constant. **The numerical value of alpha is not derived in this note.** No fit to its measured value is performed.

## 1. The research target

The first target is the low-momentum electromagnetic coupling,

\[
\alpha(0)=\frac{e^2}{4\pi\epsilon_0\hbar c}.
\]

The 2022 CODATA recommended reference is

\[
\alpha(0)^{-1}=137.035\,999\,177(21).
\]

This is a comparison value, not an input to the native operator. The uncertainty in parentheses is 0.000000021 in the inverse. The low-momentum qualification matters: a coupling defined at a different scale or in a different renormalization convention requires an explicit matching calculation.

Sources: [NIST, 2022 CODATA table](https://physics.nist.gov/cuu/pdf/all.pdf); [PDG, Electroweak Model and Constraints on New Physics, 2025](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf), sections on electromagnetic coupling and the Thomson limit.

The reference value was known during this work. A future comparison must therefore be described as retrospective unless a separate prediction protocol genuinely reserves new data. Independence of a derivation is a dependency/provenance obligation, not a claim that the researcher has never seen 137.

## 2. Existing native sources actually used

| Source | What is inherited | Boundary relevant to constants |
|---|---|---|
| RKF F00-E | Native cut-complex scalar, dagger, factorial exponential, unit-norm phase orbit | Its stated theorem does not fix the fundamental period or an electromagnetic coupling |
| RKF T53 | Native weighted-square factorization and Schur elimination | Weights are input-dependent; a positive square is not automatically the physical vacuum action |
| Thermodynamic recognition-square manuscript | A conditional classification of the lowest-order local action | The coefficients remain constitutive data in the existing statement |
| Thermodynamic spacetime-action manuscript | Gauge term \(-\langle F,F\rangle/(4g_{\rm RK}^{2})\) | The group, representation, invariant inner product, coupling and spacetime contract must retain their types |

Immutable source links:

- [F00-E](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md)
- [T53](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/53_native_cut_square_factorization_theorem.md)
- [Action classification](https://github.com/Parveen117/Publications/blob/c1a35e0f509cd137d20f9504543f97b5ec7479fc/papers/thermodynamic-response-corrections/source/sections/06a0_thermodynamic_recognition_square.tex)
- [Recognition action](https://github.com/Parveen117/Publications/blob/c1a35e0f509cd137d20f9504543f97b5ec7479fc/papers/thermodynamic-response-corrections/source/sections/06a_spacetime_recognition_datum.tex)

The source bytes read locally were checked against the corresponding remote Git blob hashes. This note is a separate research artifact; it does not amend those theorems or the existing PR certificate.

## 3. From iota to a covariant charged phase family

F00-E gives, in its completed scalar field,

\[
U(\theta)=\operatorname{Exp}_{\Sigma}(\iota_\Sigma\theta),
\qquad U(\theta)^\dagger U(\theta)=1,
\qquad U(\theta+\eta)=U(\theta)U(\eta).
\]

For a declared radial charge weight \(q\), set \(U_q(\theta)=U(q\theta)\). On an oriented native link \(x\to y\), introduce a radial connection coordinate \(a_{xy}\) and the covariant difference

\[
(D_a\psi)_{xy}=\psi_y-U_q(a_{xy})\psi_x.
\]

Under the local change of phase

\[
\psi_x\mapsto U_q(\lambda_x)\psi_x,
\qquad a_{xy}\mapsto a_{xy}+\lambda_y-\lambda_x,
\]

the exponential addition law gives

\[
(D_a\psi)_{xy}\mapsto U_q(\lambda_y)(D_a\psi)_{xy}.
\]

Consequently its native norm square is unchanged. This is an explicit covariant construction inside the earned scalar calculus. It establishes a family of charge-carrying phase representations; it does not yet identify one member with the electron.

If a normalized compact phase with fundamental period \(2\pi\) is separately established and chosen, single-valued scalar characters have integer weights in that normalization. This conditional statement is not a derivation of the electron charge, nor does it fix interaction strength. The same normalization must be transported into the gauge action. The electron need not have unit weight in every convention for all particle charges.

## 4. The coefficient alpha actually measures

Here introduce an **explicit physical adapter**: one massless electromagnetic U(1) field in four-dimensional spacetime, a normalized electron representation, the dimensionless action \(S/\hbar\), and the low-momentum electromagnetic normalization. Using rationalized natural units at this interface only, write

\[
\mathscr L_{\rm EM}=-\frac{Z_\gamma}{4}f_{\mu\nu}f^{\mu\nu},
\qquad D_\mu=\partial_\mu+\iota_\Sigma q_e a_\mu,
\qquad f=da.
\]

The quantities \(q_e\) and \(Z_\gamma>0\) refer to the same generator, matter normalization and action normalization. For constant \(Z_\gamma\), let

\[
a^{\rm can}_\mu=\sqrt{Z_\gamma}\,a_\mu.
\]

Then the field term is canonical and the electron coefficient is \(q_e/\sqrt{Z_\gamma}\). Therefore

\[
\boxed{\alpha(0)=\frac{q_e^2}{4\pi Z_\gamma(0)}}.
\]

In the manuscript's convention \(Z_\gamma=g_{\rm RK}^{-2}\), this is \(\alpha=q_e^2g_{\rm RK}^{2}/(4\pi)\), **provided that its gauge sector is the normalized physical electromagnetic sector**. The factor \(4\pi\) is part of this stated physical convention; it has not been derived from F00-E.

The target is the ratio \(q_e^2/Z_\gamma\). Setting \(q_e=1\) by a choice of phase coordinate is permitted when applied consistently, but does not evaluate the ratio. Under \(a'=t a\),

\[
q'_e=q_e/t,\qquad Z'_\gamma=Z_\gamma/t^2,
\qquad (q'_e)^2/Z'_\gamma=q_e^2/Z_\gamma.
\]

This makes the proposed target invariant under this field-coordinate normalization. It is not independent of changes to the physical cost or the dimensionless action scale.

## 5. Exact cut-memory reduction of the stiffness

Consider a finite native cut-complex carrier with a positive self-dagger quadratic form, after removing any gauge-null directions. Split it into a declared visible channel \(x\) and retained hidden channels \(h\):

\[
\mathcal E(x,h;j)=\frac12
\begin{pmatrix}x\\h\end{pmatrix}^{\dagger}
\begin{pmatrix}Z&B\\B^\dagger&C\end{pmatrix}
\begin{pmatrix}x\\h\end{pmatrix}
-\operatorname{Re}(j^\dagger q_e x),
\qquad C>0.
\]

All entries and the visible/hidden split are declared inputs to this theorem. The source in this statement acts only on the visible channel. Native T53 supplies the relevant exact elimination calculus.

Completing the square gives

\[
2\mathcal E=x^\dagger(Z-BC^{-1}B^\dagger)x
+(h+C^{-1}B^\dagger x)^\dagger C(h+C^{-1}B^\dagger x)
-2\operatorname{Re}(j^\dagger q_e x).
\]

Hence the unique stationary hidden state and reduced coefficient are

\[
h_*=-C^{-1}B^\dagger x,
\qquad
\boxed{Z_{\rm eff}=Z-BC^{-1}B^\dagger}.
\]

Positivity of the full quadratic form implies \(Z_{\rm eff}>0\). If the visible channel is scalar, its stationary source response is \(x_*=q_ej/Z_{\rm eff}\). Thus the source-source coefficient is

\[
\boxed{\gamma_{\rm eff}=q_e^2/Z_{\rm eff}}.
\]

The same result follows from the visible block of the inverse of the full matrix. Hidden coordinates therefore affect the reduced coefficient through a computable term, rather than through a renamed observation or a discarded sheet label.

For an invertible hidden change of coordinates \(h=T h'\), the blocks become \(B'=BT\), \(C'=T^\dagger CT\). Direct substitution gives \(B'(C')^{-1}(B')^\dagger=BC^{-1}B^\dagger\). The reduction is independent of that choice of hidden basis.

**Physical interpretation requires another step.** For a field theory with momentum-dependent blocks, first compute

\[
\mathcal K_{\rm eff}(p)=\mathcal K_{vv}(p)
-\mathcal K_{vh}(p)\mathcal K_{hh}(p)^{-1}\mathcal K_{hv}(p).
\]

Then establish a massless transverse photon sector whose low-momentum kernel has the required form

\[
\mathcal K_{\rm eff}(p)=Z_\gamma(0)p^2P_T+o(p^2)
\]

in the declared Euclidean response representation, with a justified physical continuation and charge normalization. An arbitrary scalar Schur complement is not automatically this photon coefficient. If hidden channels couple directly to the source, the effective source also changes and must be reduced along with the quadratic form. Nonquadratic interactions and charged quantum fluctuations require their own effective-action calculation and matching; the finite theorem does not perform that calculation.

For infinitely many helical sheets, coercivity, a lawful inverse, convergence of the reduction and stability under refinement remain necessary. The finite result does not certify an infinite-sheet limit.

## 6. Exact arithmetic control

A deliberately small fixture, chosen for transparent arithmetic and unrelated to the measured alpha, is

\[
H=\begin{pmatrix}5&1&1\\1&2&0\\1&0&3\end{pmatrix},
\qquad q_e=2.
\]

It gives

\[
Z_{\rm eff}=5-\frac12-\frac13=\frac{25}{6},
\qquad\gamma_{\rm eff}=\frac{24}{25}.
\]

The independent full source solution is

\[
H^{-1}\begin{pmatrix}2\\0\\0\end{pmatrix}
=\begin{pmatrix}12/25\\-6/25\\-4/25\end{pmatrix}.
\]

The contraction with the source is again \(24/25\). Deleting both hidden channels instead gives \(4/5\), which the control detects. This fixture verifies the reduction, not a value of the fine-structure constant.

The attached reproducible code uses rational arithmetic and checks completion of squares at 125 points, agreement with full inversion, a nontrivial hidden-basis change, four visible normalizations, sequential elimination, and sensitivity to erased memory and changed stiffness. The general statement is established by the written square-completion proof; the finite checks are implementation controls.

## 7. The next mathematical obligation

The constants problem is now a specific construction problem: derive the physical normalized pair \((q_e,Z_\gamma(0))\) from a selected native law and state. Required inputs cannot be chosen by back-solving from 137, an equivalent electromagnetic datum, or a fitted action coefficient.

The immediate candidate route is to construct the gauge quadratic form from native weighted squares and lawful helical gluing, determine the relevant charged representation, and carry the full retained memory into the effective photon coefficient. A successful selection theorem must justify the graph/continuation law, weights, vacuum and boundary sector, source representation, action normalization, and continuum or refinement limit used. This list identifies the variables to derive; it is not a theorem that doing so is impossible.

Once a coefficient is calculated, match it at a stated scale and compare with the measured value and uncertainty. A second electromagnetic observable, with all independent inputs declared, should test the same fixed prediction. No numerical alpha prediction has passed these steps here.

## 8. Other constants in the programme

| Target | Required next native result |
|---|---|
| \(\alpha(0)\) | Normalized electron charge and the low-momentum photon stiffness |
| \(m_\mu/m_e\), later \(m_p/m_e\) | Identified particle sectors and their physical mass/pole ratios; the proton requires its bound-state dynamics |
| \(Gm_e^2/(\hbar c)\) | Normalized gravitational response tied to the same matter and action scales |
| Other gauge couplings and mixing ratios | Selected gauge sectors, representations and matching at specified scales |
| \(c,\hbar,k_B\) in laboratory units | Native propagation/action/entropy relations plus an explicit unit interface |

The last row is a different kind of numerical question: present SI definitions assign exact numerical values to \(c,h,e,k_B\). A derivation should explain their physical roles and dimensionless relations, while keeping the unit convention explicit. See [BIPM, SI defining constants](https://www.bipm.org/en/measurement-units/si-defining-constants).

## 9. Evidence status

| Layer | Status |
|---|---|
| Phase covariance, field normalization, finite Schur reduction | Written derivations under stated hypotheses |
| Exact arithmetic controls | PASS for the supplied fixture and controls |
| Native selection of physical electromagnetic operator | Open construction |
| Numerical fine-structure constant | NOT DERIVED |
| Empirical validation | Not performed |
| Proof-assistant verification | Not performed |
| Worldwide originality | Not claimed; covariance, canonical normalization and Schur elimination are established mathematics used here to formulate the RKF constants target |

This is a constructive first reduction for the constants programme. It identifies the quantity to compute without treating an arbitrary dimensionless expression as the electromagnetic constant.

## Appendix A. Source hashes

```json
[
  {
    "repository": "Recognition-Kernel-Framework",
    "path": "theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md",
    "git_blob_sha1": "9ea9c78aa12d8fa0db5b4c87b9463d64a4751f41",
    "sha256": "2654b5cab54e5f103b01471ffc06ca5e32eb56b11444744ffdf7edadb27daf8c"
  },
  {
    "repository": "Recognition-Kernel-Framework",
    "path": "theorum/53_native_cut_square_factorization_theorem.md",
    "git_blob_sha1": "0791e3fb107267a1224b6e1f3cb4d7ab7b489133",
    "sha256": "0b6276f8f25a66ae3e574c7936a8a75bca5fdca29a7c10e08363c42deade3323"
  },
  {
    "repository": "Publications",
    "path": "papers/thermodynamic-response-corrections/source/sections/06a0_thermodynamic_recognition_square.tex",
    "git_blob_sha1": "760c821f934a1ce4c18a05aee12746d93f8af8b0",
    "sha256": "f555e4326579a963a76345bf17bdd942f4f155ef4acb6946dcea6f5f762e5f21"
  },
  {
    "repository": "Publications",
    "path": "papers/thermodynamic-response-corrections/source/sections/06a_spacetime_recognition_datum.tex",
    "git_blob_sha1": "497bbba637fc212847cfbc5a8036d96c2f462799",
    "sha256": "249f33a2123bb4abeed930670ce21d887bd8930e933272d03cdac2256e0b4a03"
  }
]
```

## Appendix B. Executed result

```json
{
  "status": "PASS_EXACT_FINITE_COUPLING_CONTROLS",
  "Z_eff": "25/6",
  "gamma_q_squared_over_Z_eff": "24/25",
  "checks": {
    "schur_coefficient": "25/6",
    "full_source_solution_agrees": true,
    "completion_of_squares_cases": 125,
    "hidden_basis_invariance": true,
    "visible_normalization_cases": 4,
    "sequential_elimination_agrees": true,
    "memory_deletion_detected": true,
    "physical_stiffness_change_detected": true
  },
  "alpha_physical_status": "NOT_DERIVED",
  "matrix_status": "DECLARED_FIXTURE_NOT_NATIVE_VACUUM_SELECTION"
}
```

## Appendix C. Reproducible exact control

Save the following block as `exact_checks.py` and run `python exact_checks.py` with Python 3.11 or newer. It needs no third-party packages and contains no measured-alpha input.

```python
"""Exact finite controls for the conditional RKF coupling reduction.

This file has no experimental alpha input and predicts no physical constant.
The witness matrix is a declared arithmetic fixture, not a selected vacuum.
"""
from fractions import Fraction as Q
import json


def tr(a):
    return [list(r) for r in zip(*a)]


def mul(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def inv(a):
    n = len(a)
    m = [[Q(x) for x in row] + [Q(i == j) for j in range(n)]
         for i, row in enumerate(a)]
    for c in range(n):
        pivot = next(i for i in range(c, n) if m[i][c])
        m[c], m[pivot] = m[pivot], m[c]
        d = m[c][c]
        m[c] = [x/d for x in m[c]]
        for i in range(n):
            if i != c:
                d = m[i][c]
                m[i] = [x-d*y for x, y in zip(m[i], m[c])]
    return [r[n:] for r in m]


def schur(a):
    b = [a[0][1:]]
    c = [row[1:] for row in a[1:]]
    return Q(a[0][0]) - mul(mul(b, inv(c)), tr(b))[0][0]


H = [[Q(5), Q(1), Q(1)],
     [Q(1), Q(2), Q(0)],
     [Q(1), Q(0), Q(3)]]
q = Q(2)
Z = schur(H)
gamma = q*q/Z
checks = {}

assert Z == Q(25, 6) and gamma == Q(24, 25)
checks['schur_coefficient'] = str(Z)

# Independent solution of the full sourced system.
J = [[q], [Q(0)], [Q(0)]]
state = mul(inv(H), J)
assert state == [[Q(12, 25)], [Q(-6, 25)], [Q(-4, 25)]]
assert mul(tr(J), state)[0][0] == gamma
checks['full_source_solution_agrees'] = True

# Completion of squares over a finite exact audit lattice.
n = 0
for x in map(Q, range(-2, 3)):
    for h1 in map(Q, range(-2, 3)):
        for h2 in map(Q, range(-2, 3)):
            v = [[x], [h1], [h2]]
            full = mul(mul(tr(v), H), v)[0][0]
            reduced = Z*x*x + 2*(h1+x/2)**2 + 3*(h2+x/3)**2
            assert full == reduced
            n += 1
checks['completion_of_squares_cases'] = n

# Arbitrary invertible hidden-coordinate change h = T h_new.
T = [[Q(1), Q(1)], [Q(0), Q(2)]]
P = [[Q(1), Q(0), Q(0)],
     [Q(0), T[0][0], T[0][1]],
     [Q(0), T[1][0], T[1][1]]]
assert schur(mul(mul(tr(P), H), P)) == Z
checks['hidden_basis_invariance'] = True

# Visible coordinate change x_new = t*x transports both source and cost.
for t in (Q(-3), Q(1, 2), Q(2), Q(5, 3)):
    assert (q/t)**2 / (Z/t**2) == gamma
checks['visible_normalization_cases'] = 4

# Eliminate h1 then h2; compare with simultaneous elimination.
assert Q(5) - Q(1, 2) - Q(1, 3) == Z
checks['sequential_elimination_agrees'] = True

# Sensitivity controls: deleting memory and changing physical stiffness
# must change the coupling, unlike a mere coordinate change.
assert q*q/Q(5) != gamma
assert q*q/schur([[2*x for x in row] for row in H]) == gamma/2
checks['memory_deletion_detected'] = True
checks['physical_stiffness_change_detected'] = True

print(json.dumps({
    'status': 'PASS_EXACT_FINITE_COUPLING_CONTROLS',
    'Z_eff': str(Z),
    'gamma_q_squared_over_Z_eff': str(gamma),
    'checks': checks,
    'alpha_physical_status': 'NOT_DERIVED',
    'matrix_status': 'DECLARED_FIXTURE_NOT_NATIVE_VACUUM_SELECTION'
}, indent=2))
```
