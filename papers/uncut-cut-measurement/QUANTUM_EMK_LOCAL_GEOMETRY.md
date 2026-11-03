# Quantum-EMK local geometry bridge: represented iota, response-plane invariant, and helical-memory guard

**Monty Dabas — v1.8 research extension, 26 September 2026.**

## 1. Why this bridge exists

The gravity programme must not begin by importing a Riemann tensor and then
renaming its output "recognition". The repository already contains a more
primitive geometry route:

1. the cut-complex foundation derives the quarter-turn element
   \(\iota_\Sigma\), with \(\iota_\Sigma^2=-1\) and
   \(\iota_\Sigma^\dagger=-\iota_\Sigma\);
2. the represented response layer already uses a real quarter-turn \(Z\) and
   the exact tensor \(Q=G+iA\);
3. EMK-G1 already supplies the seam metric
   \(g_A=A(v)^2du^2+dv^2\), equivalently \(W=A^2\), and its exact seam
   curvature formula;
4. EMK-G3 proves that local flatness does not erase global sheet memory.

The missing obligation is therefore a typed bridge between (1), (2), and the
**local EMK seam geometry**, while keeping global helical memory separate.
This note constructs that bridge. It does not derive physical gravity.

## 2. Source pins and claim boundary

This extension consumes the following already certified or source-bound inputs.

- RKF F00-E, blob \`9ea9c78aa12d8fa0db5b4c87b9463d64a4751f41\`:
  oriented cut -> \(\iota_\Sigma\), with
  \(\iota_\Sigma^2=-1\) and
  \(\iota_\Sigma^\dagger=-\iota_\Sigma\).
- Publications EMK recognition geometry, README blob
  \`5aa32e20b1b8c74cccff7aab2b07cb682b7eb8fd\`:
  \(g_A=A(v)^2du^2+dv^2\),
  \[
  K={W'^2\over4W^2}-{W''\over2W},
  \qquad W=A^2,
  \]
  plus the G1/G2/G3 separation guards.
- The current branch response tensor, blob
  \`02c2b1389170b79023b3875095ea61fc8cf076c1\`, which explicitly says its
  represented \(Z\) is not automatically the native central scalar
  \(\iota_\Sigma\).

Nothing below equates recognition curvature with Gaussian curvature. Nothing
below identifies sheet memory with a physical double helix. A physical length
scale, probe law, spacetime interpretation and experimental gravity test remain
separate obligations.

## 3. QG-1 — represented-iota gate

Let \(V\) be a real response carrier with symmetric positive form \(H\). A
representation of the native cut-complex scalar on \(V\) is admitted only after
a real operator \(Z\) satisfies

\[
Z^2=-I,\qquad Z^T H=-HZ.
\]

Then

\[
\rho(a+\iota_\Sigma b)=aI+bZ
\]

respects multiplication and dagger on the quadratic cut-complex algebra.

**Proof.** Multiplication follows directly from \(Z^2=-I\). For the dagger,
\[
(aI+bZ)^{\dagger_H}=aI-bZ
=\rho(a-\iota_\Sigma b).
\]
Thus the represented quarter-turn may be typed as \(\rho(\iota_\Sigma)\)
only after these identities pass. An arbitrary matrix called \(i\) does not
pass merely by typography. **Proof complete.**

This closes the typing gap in the response-tensor module without claiming that
every representation of the native field is faithful, unique or physical.

## 4. QG-2 — a basis-invariant response-plane scalar

Take two calibrated response columns \(Y=[y_1,y_2]\). On an admitted
represented-iota carrier define

\[
G=Y^T H Y,\qquad A=-Y^T H ZY,\qquad Q=G+iA.
\]

Assume the real response plane is nondegenerate, so \(\det G>0\). Since a
two-dimensional antisymmetric matrix has the form

\[
A=\begin{pmatrix}0&a\\-a&0\end{pmatrix},
\]

define

\[
\boxed{\chi(Y)={\det A\over\det G}={a^2\over\det G}.}
\]

### Proposition QG-2.1 — label-basis invariance

For every \(B\in GL(2,\mathbb R)\), \(Y\mapsto YB\) gives

\[
G\mapsto B^TGB,\qquad A\mapsto B^TAB.
\]

Therefore both determinants acquire the same factor \(\det(B)^2\), hence

\[
\boxed{\chi(YB)=\chi(Y).}
\]

For the represented complex Gram tensor \(Q\ge0\),

\[
\det Q=\det G-a^2\ge0,
\]

so

\[
\boxed{0\le\chi\le1.}
\]

This scalar is deliberately orientation-blind: reversing the response-label
orientation flips \(A\) but leaves \(\chi\) unchanged. The signed orientation
record must therefore remain outside this scalar metric adapter.

## 5. QG-3 — local EMK seam two-jet, not an imported spacetime metric

Use the EMK-G1 reflection-symmetric seam metric in its exact \(W=A^2\)
presentation. Normalize one local seam chart by

\[
W(0)=1,\qquad W'(0)=0.
\]

The second equality is not an extra curvature ansatz: it is the EMK reflection
condition at the fixed seam. Supply a positive independently calibrated
length-squared \(\ell^2\). Define the candidate local coupling

\[
\boxed{\kappa_{\rm loc}={\chi\over\ell^2}}
\]

and constrain only the second jet,

\[
\boxed{W''(0)=2\kappa_{\rm loc}.}
\]

EMK-G1 then gives immediately

\[
K(0)
={W'(0)^2\over4W(0)^2}-{W''(0)\over2W(0)}
=-\kappa_{\rm loc},
\]

hence

\[
\boxed{K_{\rm EMK}(0)=-{\chi\over\ell^2}.}
\]

Equivalently,

\[
W(v)=1+{\chi\over\ell^2}v^2+O(v^4).
\]

The \(O(v^4)\) terms are **not selected** here. Thus this bridge does not
smuggle in the global quadratic warp, a Riemannian spacetime, Einstein's
equations, or even a global surface. It selects one local EMK metric two-jet
from one represented response invariant after a supplied physical scale.

The map \(\chi\mapsto W''(0)\) is a declared candidate constitutive adapter.
Its usefulness must be decided by a later physical calibration or rejection
test. EMK-G1's guard remains active: recognition curvature and metric curvature
are still distinct objects.

## 6. QG-4 — lawful complex cuts are not automatically geometry-faithful

A subtle obstruction appears immediately. Suppose a lossy readout projector
\(P\) commutes with \(Z\). Then the represented complex structure descends and
the U31 seam cross-term vanishes. Nevertheless \(\chi\) is a nonlinear ratio,
so in general

\[
\chi(PY)\ne\chi(Y).
\]

The exact certificate uses

\[
Y=
\begin{pmatrix}
1&0\\
0&1\\
0&1\\
0&0
\end{pmatrix},
\qquad
Z=R\oplus R,
\qquad
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix}.
\]

Uncut,

\[
G=\operatorname{diag}(1,2),\quad
A=\begin{pmatrix}0&1\\-1&0\end{pmatrix},\quad
\chi={1\over2}.
\]

Projecting lawfully onto the first \(Z\)-invariant real pair gives

\[
G_{\rm vis}=I,\quad A_{\rm vis}=A,\quad \chi_{\rm vis}=1.
\]

Therefore

\[
\boxed{\text{\(Z\)-compatibility is necessary for the complex ledger, but
not sufficient for the \(\chi\)-geometry target.}}
\]

A separate target-faithfulness gate is required. In the certificate, a
two-axis non-\(Z\)-invariant readout has rank two, while adjoining its
\(Z\)-partners raises the row rank to four. Exactly two additional real scalar
channels repair the represented-iota carrier on that fixture.

## 7. QG-5 — local geometry and helical sheet memory remain a pair

EMK-G3 already proves that a flat local connection can carry nontrivial global
monodromy. The present bridge therefore refuses to compress the geometry state
to \(K(0)\) alone. At minimum the research state must retain

\[
(\chi,\;K_{\rm EMK}(0),\;q_{\rm sheet})
\]

or an equivalent local/global ledger.

The exact control uses two records with identical \(\chi\) and identical local
EMK two-jet but sheet labels \(q=0\) and \(q=1\). A local metric observer cannot
distinguish them. This is precisely the information EMK-G3 warns us not to
erase.

This is a **helical-memory guard**, not a physical double-helix claim.

## 8. Exact rational witness

The executable certificate
\`certificates/quantum_emk_local_geometry.py\` uses exact
\`fractions.Fraction\` arithmetic.

For the fixture above,

\[
\det G=2,\qquad \det A=1,\qquad \chi={1\over2}.
\]

With supplied \(\ell^2=4\),

\[
\kappa_{\rm loc}={1\over8},\qquad
W''(0)={1\over4},\qquad
K_{\rm EMK}(0)=-{1\over8}.
\]

It additionally checks:

- \(Z^2=-I\) and \(Z^T=-Z\);
- rejection of \(I\) as a false representation of \(\iota_\Sigma\);
- exact \(GL(2,\mathbb R)\) label-basis invariance of \(\chi\);
- orientation reversal \(A\mapsto-A\) with unchanged \(\chi\);
- the \(Z\)-compatible-cut counterexample \(\chi:1/2\mapsto1\);
- two-channel paired repair on the non-\(Z\)-invariant fixture;
- same local two-jet with different sheet memory.

## 9. What is now earned, and what is not

Earned in this extension:

- an explicit gate for calling the response quarter-turn a represented native
  \(\iota_\Sigma\);
- a dimensionless two-response invariant \(\chi\) with exact label-basis
  covariance;
- a local EMK-G1 metric two-jet adapter using \(\chi\) and a supplied scale;
- a proof that complex-structure-safe cuts need not preserve the geometry
  target;
- an explicit local-curvature/global-sheet-memory separation.

Still open:

- deriving \(\ell\) from the native law rather than calibration;
- selecting this \(\chi\to W''(0)\) adapter uniquely against competitors;
- deriving spacetime dimension or a Lorentzian spacetime metric;
- identifying a universal gravitational source and probe coupling;
- deriving \(\hbar\), a physical clock, Born frequencies or quantum detector
  statistics;
- identifying EMK sheet memory with a physical double helix;
- deriving Einstein equations or any replacement field equation;
- empirical validation of a native quantum-gravity prediction.

The next discriminating step is no longer "assume curvature." It is to bind the
same represented-iota response invariant to a **physical source/probe scale**
already exposed by the atom-interferometer and source-response interfaces, and
ask whether one predeclared adapter survives more than one geometry and source
configuration.
