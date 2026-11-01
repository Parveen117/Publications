# YM86 — Smooth relational measurements and finite-energy response tails

Monty Dabas. 6 October 2026. Python 3.12 evidence packet.

## What is established

YM85 found that normalized exterior-frame coordinates can be physical
L2 functions with infinite electric energy. This continuation replaces
them by smooth local measurements and gives an energy-dependent bound
on response omitted by a complete native electric cutoff.

- Seven bounded polynomial measurements separate the common-conjugacy
  orbits of three SU(2) plaquette holonomies, including degenerate
  configurations. All have finite full electric energy. Their three
  single traces and three pair measurements alone miss an orientation
  sign; an antisymmetric ordered triple trace recovers it.
- For three edge-disjoint plaquettes meeting at a vertex, the complete
  electric energy matrix of the three relational measurements is explicit.
  In the actual homogeneous vacuum their energy upper bounds are uniform
  in total volume and tend to zero in the bare clock as theta increases.
  No variance lower bound or normalized continuum-state survival follows.
- For a normalized physical observable f supported on S, with full bare
  energy at most E, the actual-measure projection onto a **complete**
  native cutoff obeys
  \(\|(I-P_\Lambda)f\|_\mu^2\le D_SE/\Lambda\).
  The density ratio D_S has an explicit all-coupling bound from YM83,
  independent of total volume for fixed S.
- A finite matrix of full output responses then bounds the response of
  every such f by \(c_\Lambda+\sqrt{D_SE/\Lambda}\), capped at one.
  This is an energy-constrained statement, not an unrestricted operator
  norm or a completed Yang--Mills gap criterion.

The coordinate separation uses classical SU(2) trace/Gram algebra; the
cutoff uses established Peter--Weyl theory. The development here is their
explicit combination with the actual-vacuum energy and density estimates
to repair the measurement protocol identified in YM85. No new general
theory of trace invariants or spin networks is claimed.

## 1. Actual model and a local three-loop configuration

Keep the periodic cubic lattice of side L>=4, the native link Casimir
3/4, and

\[
 A_\theta=L_E-\theta\sum_pW_p,\quad H=A_\theta-E_0,
 \quad d\mu=h^2d\Phi,\quad
 \mathcal E(f,g)=\mu\sum_e\Gamma_e(f,g).
 \tag{1}
\]

h is its positive normalized ground state; the form of H on hf is
\(\mathcal E(f)\). No trial Gibbs density replaces this measure.
Use three plaquettes based at o, in the oriented quadrants

\[
 X:(+x,+y),\qquad A:(-x,+z),\qquad B:(-y,-z).
 \tag{2}
\]

Each notation means the four steps in the two displayed directions and
then their inverses. These loops have disjoint edges and meet only at o.
Their union S has twelve links. All three holonomies transform by the
same conjugation at o, while transformations at other vertices cancel.

Write unit quaternions
\(A=(a_0,\mathbf a)\), \(B=(b_0,\mathbf b)\),
\(X=(x_0,\mathbf x)\). Set

\[
 u=\mathbf a\cdot\mathbf b,\quad
 d_1=\mathbf a\cdot\mathbf x,\quad
 d_2=\mathbf b\cdot\mathbf x,\quad
 d_3=(\mathbf a\times\mathbf b)\cdot\mathbf x,
 \quad \mathcal I=(a_0,b_0,x_0,u,d_1,d_2,d_3).
 \tag{3}
\]

Every entry is a bounded smooth gauge-invariant polynomial of the
twelve link quaternions. At any finite regulator it is in the full
electric form domain, including at parallel or central holonomies.
There is no division by a reference-vector norm.

## 2. Theorem 86.1 — Complete three-loop information with finite energy

The seven measurements (3) separate all simultaneous SU(2) conjugacy
orbits of (A,B,X). With \(W(Y)=\tfrac12\operatorname{Tr}Y\), they
can be read from Wilson words using

\[
 u=a_0b_0-W(AB),\quad
 d_1=a_0x_0-W(AX),\quad d_2=b_0x_0-W(BX),
 \quad
 \boxed{d_3=\tfrac12[W(BAX)-W(ABX)].}
 \tag{4}
\]

Once the first six readings are known, one ordered triple trace is
sufficient: \(d_3=a_0b_0x_0-x_0u-a_0d_2-b_0d_1-W(ABX)\).
The antisymmetric expression in (4) displays the missing orientation
information directly.

They satisfy the Gram relation

\[
 d_3^2=\det\begin{pmatrix}
 1-a_0^2&u&d_1\\
 u&1-b_0^2&d_2\\
 d_1&d_2&1-x_0^2
 \end{pmatrix}.
 \tag{5}
\]

**Proof.** Quaternion multiplication gives (4). The scalar parts and
pair products give the Gram matrix of the three imaginary vectors.
Two vector triples with the same Gram matrix are related by an
orthogonal map: the map between their spans is an isometry and extends
to the ambient space. If the vectors have rank three, equality of the
signed triple product makes that map orientation preserving. If the
rank is less than three, its extension on the orthogonal complement
can be chosen to have positive determinant. Every such SO(3) map is
SU(2) conjugation on imaginary quaternions. This proves separation,
including every degenerate case. Equation (5) is the Gram determinant
identity. Smoothness gives finite energy on the compact lattice. QED.

This is completeness for the common-conjugacy orbit of these **three
holonomies**, not for the configuration of the entire lattice. On the
three-loop subgraph, internal degree-two vertex gauges reduce its
configuration precisely to that orbit. The seven readings separate
configurations; their seven-dimensional linear span is not the full
observable Hilbert space. Higher polynomial combinations are needed
when an electric cutoff is enlarged.

### Pairwise blindness and the unavoidable conditioning loss

For 0<t<1 put c=sqrt(1-t^2) and take

\[
 A=(c,t,0,0),\quad B=(c,0,t,0),\quad
 X_+=(c,0,0,t),\quad X_-=(c,0,0,-t).
 \tag{6}
\]

The first six readings in (3), equivalently all three single and three
specified pair traces, are identical. But \(d_3=\pm t^3\), so the
two triples are not gauge equivalent.

In the Euclidean quaternion distance minimized over common conjugation,
their orbit distance is exactly 2t, while the Euclidean difference of
their seven readings is 2t^3. To see the distance, compare tI_3 with
\(t\operatorname{diag}(1,1,-1)\) under left SO(3) rotations. An
orientation-reversing orthogonal 3x3 matrix has trace at most one,
so the minimum squared Frobenius distance is 4t^2, attained here.
Consequently no uniform inverse estimate with Holder exponent greater
than 1/3 can hold near the identity triple. This does not assert a
global Holder upper bound with exponent 1/3.

For a fixed nondegenerate reference pair put
\(\Delta=|\mathbf a\times\mathbf b|^2>0\) and let K have rows
\(\mathbf a^T,\mathbf b^T,(\mathbf a\times\mathbf b)^T\).
Then \(\mathbf d=K\mathbf x\), and

\[
 KK^T=\begin{pmatrix}
 |\mathbf a|^2&u&0\\u&|\mathbf b|^2&0\\0&0&\Delta
 \end{pmatrix}\ge\frac\Delta2I_3.
 \tag{7}
\]

The upper-left block has determinant Delta and trace at most two,
which proves the bound. In particular
\(\|\delta\mathbf x\|\le\sqrt{2/\Delta}\,
\|\delta\mathbf d\|\) when only the three d-readings have error
and the reference pair is held fixed. The polynomial readouts stay
smooth as Delta tends to zero; recovering normalized frame coordinates
does not stay uniformly conditioned. No error bound for uncertain
reference data is silently included in (7).

## 3. Theorem 86.2 — Exact full electric energy and actual-vacuum bounds

Write \(r_a^2=|\mathbf a|^2\), \(r_b^2=|\mathbf b|^2\),
\(r_x^2=|\mathbf x|^2\). For the geometry (2), the full native
gradient pairings of the three d-observables are

\[
 \begin{aligned}
 \Gamma(d_1)&=r_a^2+r_x^2-2d_1^2,\\
 \Gamma(d_2)&=r_b^2+r_x^2-2d_2^2,\\
 \Gamma(d_1,d_2)&=u-d_1d_2,\\
 \Gamma(d_1,d_3)&=-2d_1d_3,\qquad
 \Gamma(d_2,d_3)=-2d_2d_3,\\
 \Gamma(d_3)&=|\mathbf b\times\mathbf x|^2
                +|\mathbf x\times\mathbf a|^2
                +|\mathbf a\times\mathbf b|^2-3d_3^2.
 \end{aligned}
 \tag{8}
\]

These include the derivatives of **both reference loops**. They are
not the target-only partial energy matrix from YM85. Also

\[
 L_Ed_1=6d_1,\qquad L_Ed_2=6d_2,\qquad L_Ed_3=9d_3.
 \tag{9}
\]

Equation (9) concerns the Haar electric operator; it does not make
the d-observables eigenvectors of the interacting Hamiltonian.

For the actual homogeneous vacuum let
\(\delta=\mu(1-W_p)\) and \(q=\mu(1-W_p^2)\). Lattice symmetry
makes these common to all three loops. Then

\[
 \boxed{\mathcal E(d_1),\mathcal E(d_2)\le2q\le4\delta,
        \qquad\mathcal E(d_3)\le3q\le6\delta.}
 \tag{10}
\]

YM76 gives, for theta>=25,

\[
 \mathcal E(d_1),\mathcal E(d_2)\le9/\sqrt\theta,
 \qquad\mathcal E(d_3)\le27/(2\sqrt\theta).
 \tag{11}
\]

These full-energy upper bounds are independent of total volume.

**Proof.** A four-link loop contributes the unit-sphere quaternion
gradient: for a scalar function of its imaginary part with Euclidean
gradient v, the corresponding loop energy is
\(|v|^2-(\mathbf a\cdot v)^2\). This follows by summing the four
native link factors 1/4. The three supports are edge disjoint, so add
their energies. Applying this rule to dot and triple products gives
(8). Each imaginary coordinate has four-link electric eigenvalue 3;
multilinearity and disjoint supports give (9).

The first two upper bounds in (10) follow directly from (8). For the
third, discard the negative term and bound each squared cross product
by the product of squared radii. For 0<=v,w<=1, vw<=(v+w)/2, so the
sum of those three products is at most \(r_a^2+r_b^2+r_x^2\).
Take actual vacuum expectations and use q<=2delta. QED.

The norms also shrink: \(\mu d_k^2\le q\le2\delta\) for all
three k. Strict positivity of the finite-regulator measure gives
positive variance, but these estimates give no useful uniform lower
bound on it. Hence (11) does not control the quotient
\(\mathcal E(d_k)/\operatorname{Var}_\mu d_k\).
The physical energy is \(\kappa\mathcal E\), and
\(\kappa\sqrt\theta=c/a\). For example the upper bound for the
unnormalized d_1 energy becomes \(9c/(a\theta)\); neither this
expression nor a small bare energy proves continuum survival.

## 4. Theorem 86.3 — An actual-measure tail bound using energy

Let S be any finite link support, n=|S|, and w_S the density of its
actual marginal relative to product Haar. Write

\[
 m_S=\min w_S>0,\quad M_S=\max w_S,\quad D_S=M_S/m_S.
 \tag{12}
\]

YM83's ratio estimate implies, for \(\Omega_S>0\),

\[
 D_S\le\min\{(21/19)^{2n}e^{18\Omega_S},
               64^n e^{4\pi\sqrt{3n\Omega_S}}\},
 \qquad
 \Omega_S=2\theta\,\#\{p:p\cap S\ne\varnothing\}.
 \tag{13}
\]

If theta=0, D_S=1. At fixed S, (13) is independent of the total
volume, but deteriorates with coupling and support size.

Let V_Lambda be the **whole** gauge-invariant native electric spectral
space on S with Haar electric energy at most Lambda>0, including
constants. Let \(\Pi_\Lambda^\Phi\) be its Haar orthogonal projection,
and \(P_\Lambda\) the orthogonal projection onto the same function
space in the actual inner product \(L^2(\mu)\).
They need not be the same projection.

For every S-supported physical function in the full electric form domain,

\[
 \boxed{\|(I-P_\Lambda)f\|_\mu^2
       \le\frac{D_S}{\Lambda}\,\mathcal E(f).}
 \tag{14}
\]

**Proof.** Squaring and integrating YM83's pointwise bound
\(h(x,z)/h(y,z)\le R_S\) shows
\(w_S(x)/w_S(y)\le R_S^2\). Its two explicit R_S bounds give
(13); integrating out the exterior does not introduce a volume factor.

The gauge action commutes with the Haar electric Laplacian on S.
Its complete spectral projection therefore preserves physical
functions. The actual orthogonal projection is the best approximation
in its own norm, giving

\[
 \begin{aligned}
 \|f-P_\Lambda f\|_\mu^2
 &\le M_S\|f-\Pi_\Lambda^\Phi f\|_\Phi^2\\
 &\le\frac{M_S}{\Lambda}\mathcal E_\Phi(f)
 \le\frac{M_S}{m_S\Lambda}\mathcal E(f).
 \end{aligned}
 \tag{15}
\]

The middle inequality is the Haar spectral theorem on the omitted
eigenvalues, all greater than Lambda. The two density comparisons use
the marginal w_S and the full gradient of the S-supported function.
Form closure proves (14). QED.

The finite spaces consist of native harmonic polynomials followed by
gauge averaging, as in YM74. Their union is a form core: Haar spectral
truncations converge in Haar H1, and m_S,M_S make that norm equivalent
to the actual form norm on the fixed support. This assertion needs no
uniform equivalence as the regulator changes.

For the twelve-link geometry (2), exactly thirty plaquettes meet S.
To count them, the twelve links contribute 48 incidences. Each of the
three selected faces is counted four times, removing nine excess
incidences. The other nine faces meeting o each have two of the six
incident S-links, removing nine more. There are no other multiple
incidences for L>=4. Thus \(\Omega_S=60\theta\), and (13) becomes

\[
 D_S\le\min\{(21/19)^{24}e^{1080\theta},
                 64^{12}e^{48\pi\sqrt{15\theta}}\}.
 \tag{16}
\]

These deliberately conservative bounds certify finiteness. They are
not claimed computationally cheap at weak coupling.

## 5. Theorem 86.4 — Response control on a declared energy budget

Let T be a contraction on physical \(L^2(\mu)\) with T1=0. A concrete
example is \(T=P_BP_C-P_{B\cup C}\) for disjoint blocks, using the
actual conditional projections of YM83. Output functions can involve
links outside S.

Choose an independent centered basis \(\phi_1,\ldots,\phi_N\) for
\(V_\Lambda\cap\mathbf1^\perp\) in the actual inner product. Define

\[
 C_{ij}=\langle\phi_i,\phi_j\rangle_\mu,\qquad
 R_{ij}=\langle T\phi_i,T\phi_j\rangle_\mu,
 \quad c_\Lambda^2=\lambda_{\max}(C^{-1/2}RC^{-1/2}).
 \tag{17}
\]

For an empty centered space put c_Lambda=0. These are actual Gram
matrices; R contains the **full outputs**, not a truncated output
projection. If \(0\le c_\Lambda\le\bar c_\Lambda\le1\) is certified,
every centered S-supported physical f with norm one and
\(\mathcal E(f)\le E\) satisfies

\[
 \boxed{\|Tf\|_\mu\le
    \min\{1,\bar c_\Lambda+\sqrt{D_SE/\Lambda}\}.}
 \tag{18}
\]

**Proof.** Because constants belong to V_Lambda, actual orthogonal
projection preserves centering. Decompose
\(f=P_\Lambda f+(I-P_\Lambda)f\). The first response has norm at
most \(\bar c_\Lambda\|P_\Lambda f\|\le\bar c_\Lambda\).
Contraction bounds the second by its input norm, which (14) bounds.
Contraction also supplies the cap one. QED.

For the conditional-projection example its full output matrix has the
exact alternative expression

\[
 R_{ij}=\langle P_C\phi_i,P_BP_C\phi_j\rangle_\mu
        -\langle P_{B\cup C}\phi_i,P_{B\cup C}\phi_j\rangle_\mu.
 \tag{19}
\]

It follows by expanding T* T and using nesting of the projections.
No output Hilbert--Schmidt tail is dropped. If R is replaced by a
measured or truncated matrix, its uncertainty or omitted output must
be separately bounded. Finite positive semidefiniteness of
\(\bar c_\Lambda^2C-R\) certifies the finite part of (18).

To achieve an input-tail tolerance epsilon, it suffices that
\(\Lambda\ge D_SE/\epsilon^2\). For a physical energy budget E_phys,
use E=E_phys/kappa. At fixed finite regulator and support this proves
uniform approximation on the bounded-energy unit ball. It does not
bound the rank or cost uniformly in the continuum limit.

### Limits of the new estimate

- The complete native cutoff is essential. Seven orbit-separating
  readings do not constitute seven basis vectors capturing every
  bounded-energy function. An arbitrary small marker subspace cannot
  be substituted for V_Lambda in (14).
- A bound on a bounded-energy class cannot be substituted for the
  unrestricted coefficients in YM83/YM84's operator row inequality.
  Conditional projections of a state also require energy control
  before that argument can be localized to an energy class.
- The actual response matrices in (17) have not been bounded at weak
  coupling in this stage. Fixed-support density comparison and finite
  energy do not discharge that missing input.
- Whole-lattice supports grow under spatial refinement. The constants
  in (13), cutoff size, physical clock and nonzero variances all require
  additional scale control to construct continuum states.

Thus YM86 supplies a smooth local replacement and a valid response-tail
theorem with an energy budget. A uniform full-lattice gap and a nontrivial
four-dimensional continuum Yang--Mills theory remain open.

## 6. Verification and lineage

The Python packet checks quaternion trace identities, mirror pairs,
conditioning matrices, twelve-link incidence and full native energy
Grams. Independent weighted-Haar polynomial fixtures check the actual-
versus-Haar projection distinction and the energy-tail bound. Finite
matrix fixtures check full-output response bounds and reject unjustified
finite norm estimates. The fixtures are not interacting vacuum samples.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym86_finite_energy_response.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym86_finite_energy_response.py -v
```

Finite checks certify their declared arithmetic and source bytes. The
universal assertions are the written proofs, with no independent external
review or proof-assistant verification claimed. Prior evidence is retained.

- [YM85](YM85_EXTERIOR_FRAMES_AND_PHYSICAL_ENERGY.md): the normalized
  frame energy obstruction being repaired at the measurement level.
- [YM83](YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE.md): actual local
  density comparison and conditional projections.
- [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md): actual defect
  bound, native normalization and the physical clock.
- [YM74](YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md): construction of
  complete finite native electric spaces.
- W. M. Goldman, *Trace Coordinates on Fricke spaces of some simple
  hyperbolic surfaces*, Handbook of Teichmuller Theory II (2009),
  [arXiv:0901.1404](https://arxiv.org/abs/0901.1404), for established
  trace-coordinate invariant theory of matrix pairs and triples.
- J. C. Baez, *Spin Network States in Gauge Theory*, Adv. Math. 117
  (1996), 253--272, [arXiv:gr-qc/9411007](https://arxiv.org/abs/gr-qc/9411007),
  for the established physical basis on a fixed graph.
