# ASM-1 — Arithmetic source memory and finite cut reduction

17 October 2026. Written continuation at Publications commit
e2528feffbe28e1917eef35e60038af27dcffd31, following
[CLA-1](RH_CLOSURE_AUDIT.md). The mathematical inputs are the public
native cut calculus and the actual arithmetic source already derived
in [CWF](RH_CRITICAL_WAVE_FORMS.md), [PTE](RH_PRIME_TRANSPORT_ENERGY.md)
and [FPW](RH_FIRST_PRIME_WINDOW.md).

The new step is constructive. At every fixed finite support width,
the actual Gamma kernel makes all sufficiently high Legendre layers
strictly positive, with an explicit lower bound. The remaining sign
is carried by a finite matrix after all higher-layer memory is
retained through an exact Schur complement. Nested eliminations
give the same effective potential.

The matrix is finite in dimension; its exact memory term initially
contains an infinite-dimensional inverse. A second theorem supplies
finite trial matrices and an explicit residual bound for that term.
This is a source-preserving reduction, not a computation or sign
certificate for the resulting matrices.

**Status: RH remains unproved.** What is proved here is an explicit
positive tail, an exact finite cut reduction for every finite window,
a coherent refinement law, and outward residual inequalities for the
remaining matrix. Positivity of those actual matrices at every width
is still required. No numerical source evaluation, local scientific
computation, independent agent, replay or workflow dispatch is used.

## 1. Native dependency order and the actual source

The public [native foundation order](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/operator_foundation/theory/NATIVE_PRIMITIVE_ORDER.md)
places typed cuts, cut scalars, composition and retained residue
before a declared completion or representation.
The [clock-free cut calculus, T24](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/theorum/24_clock_free_recognition_seam_cut_calculus.md)
requires memory to survive composition.
The [finite native square theorem, T53](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/operator_foundation/sources/rkf_reference/theorum/53_native_cut_square_factorization_theorem.md)
factorizes an admitted square; it does not establish that every
signed arithmetic source is already one.

Here the positive Gamma difference energy and its completion have
already been constructed from the source in PTE. We use that declared
representation to build a particular lawful cut. No primitive
positive pairing for the signed Weil form is assumed.

Fix L>0 and I=(-L,L). Let H_L be the real Hilbert space of even
L2 functions on I, extended by zero outside I. Use
\[
T_u g(x)=g(x-u),\qquad
C_g(u)=\int_{\mathbb R}g(x)g(x-u)\,dx,\qquad
\ell(g)=\int_I g(x)\cosh(x/2)\,dx .
\tag{1}
\]
The actual arithmetic form is
\[
\boxed{
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}\|g\|_2^2
+\mathcal E_\Gamma[g]
-\sum_{n\in\mathcal N_L}a_n C_g(\log n),
}
\tag{2}
\]
where
\[
\begin{gathered}
\mathcal N_L=\{n\ge2:n\text{ is a prime power},\ \log n<2L\},\\
a_n=\frac{\Lambda(n)}{\sqrt n},\qquad
S_L=\sum_{n\in\mathcal N_L}a_n,\qquad
c_L=S_L-\frac{\Omega_0}{2},\\
\Omega_0=-\gamma_E-\frac{\pi}{2}-3\log2-\log\pi,\\
k(u)=\frac{e^{-u/2}}{1-e^{-2u}},\qquad
\mathcal E_\Gamma[g]=\frac12\int_0^\infty k(u)
\|g-T_ug\|_2^2\,du,\\
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}).
\end{gathered}
\tag{3}
\]
All prime powers meeting the support occur with their actual
coefficients. The scalar Omega_0/2 is retained.

Let D_L be PTE's closure of the even smooth compact tests in
the norm \((\|g\|_2^2+\mathcal E_\Gamma[g])^{1/2}\).
The form (2) is closed and bounded below on D_L: its other terms
are bounded forms on H_L. Denote its self-adjoint operator by
\(\mathcal T_L\). Positivity of this operator is not assumed.
Write \(q_L(g,h)\) for the symmetric bilinear polarization.

On the smooth core, \(q_L=\mathscr Q_\xi\) in CWF's normalization.
Thus
\[
\mathrm{RH}\quad\Longleftrightarrow\quad
q_L[g]\ge0\text{ for all }L>0,\ g\in D_L .
\tag{4}
\]
The extension from smooth tests in (4) is by form-core density,
not an additional zero-sum assertion on nonsmooth functions.

### Why the diagonal subtraction stays visible

PTE's coefficientwise identity is
\[
q_L[g]=\|\mathcal A_Lg\|^2-c_L\|g\|_2^2 .
\tag{5}
\]
The public [T36 source adapter](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/theorum/36_native_prime_gamma_debt_source_adapter.md)
similarly writes \(S=M_{\rm mis}-d_0I\) and uses
\((I-d_0M_{\rm mis}^{-1})^{1/2}\).
The positive square root requires the corresponding source
inequality on its declared support. Adapter algebra after that
root cannot supply the inequality retroactively.

The [T38 form-range capsule](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/theorum/38_native_compatibility_form_range_source_domination.md)
displays the positive symbol
\((g_0+m_{\rm p})I_2-\mathcal Q\), with \(g_0=G+m_\Gamma\).
Those are also the positive raw-mismatch contributions listed
in T36. To consume them as the post-debt source, an exact chart
and metric identity carrying the subtraction is needed; matching
the positive summands alone does not supply it. This note proves
its estimates directly from (2), and does not use the T38 scalar
decoder envelope as a sign certificate for (2).

## ASM-T1. Every finite window has an explicit positive unresolved tail

Let
\[
\phi_j(x)=\sqrt{\frac{4j+1}{2L}}P_{2j}(x/L),
\qquad j=0,1,\ldots,
\tag{6}
\]
be the even orthonormal Legendre basis of H_L. Define the aperture
and unresolved channel by
\[
V_N=\operatorname{span}\{\phi_0,\ldots,\phi_{N-1}\},
\qquad P_N:H_L\to V_N,\qquad Q_N=I-P_N .
\tag{7}
\]
The index N counts retained even modes. The unresolved fibre over
v is \(v+(Q_ND_L)\).

Put \(H_m=\sum_{j=1}^m1/j\), and
\[
\boxed{
d_{L,N}=F(L)-c_L+\frac{e^{-L}}2H_{2N}.
}
\tag{8}
\]
Then every \(h\in D_L\) with \(P_Nh=0\) satisfies
\[
\boxed{q_L[h]\ge d_{L,N}\|h\|_2^2.}
\tag{9}
\]
For every fixed finite L there is a finite N with \(d_{L,N}>0\).
For example, any integer N satisfying
\[
2N\ge
\exp\!\left(2e^L\bigl(\max\{c_L-F(L),0\}+1\bigr)\right)
\tag{10}
\]
has \(d_{L,N}\ge1\). This is an existence bound; it is not claimed
to be an efficient cut size.

**Proof.** The exact Gamma boundary decomposition is
\[
\begin{aligned}
\mathcal E_\Gamma[h]
={}&\frac14\int_I\int_I
k(|x-y|)|h(x)-h(y)|^2\,dy\,dx\\
&+\frac12\int_I |h(x)|^2
\bigl(F(L-x)+F(L+x)\bigr)\,dx .
\end{aligned}
\tag{11}
\]
Convexity of F bounds the second line below by \(F(L)\|h\|^2\).
For \(0<u\le2L\),
\[
1-e^{-2u}\le2u,\qquad e^{-u/2}\ge e^{-L},
\qquad k(u)\ge\frac{e^{-L}}{2u}.
\tag{12}
\]
Define
\[
J_L[h]=\frac18\int_I\int_I
\frac{|h(x)-h(y)|^2}{|x-y|}\,dy\,dx .
\]
FPW-T2 proves the Legendre energy identity by polynomial division
and symmetry:
\[
\int_I\frac{P_m(x/L)-P_m(y/L)}{|x-y|}\,dy
=2H_mP_m(x/L).
\tag{13}
\]
Finite projections and nonnegative remainder energy give
\[
J_L[h]\ge\frac12\sum_{j\ge0}
H_{2j}|\langle h,\phi_j\rangle|^2.
\tag{14}
\]
The proof of (13)--(14) is scale invariant and has no restriction
L<=1/2; that restriction in FPW concerned its sharper kernel bound.

If \(P_Nh=0\), (14) and completeness give
\(J_L[h]\ge H_{2N}\|h\|^2/2\). Therefore
\[
\mathcal E_\Gamma[h]\ge
\left(F(L)+\frac{e^{-L}}2H_{2N}\right)\|h\|^2.
\tag{15}
\]
Each translation is unitary on full-line L2, so
\(|C_h(u)|\le\|h\|^2\). Drop the nonnegative pole square in (2)
and bound the finite prime sum below by \(-S_L\|h\|^2\).
This proves (9).

Finally \(H_{2N}\ge\log(2N+1)\). Substitute (10) into (8).
If \(c_L\ge F(L)\), the bound is at least 1; if \(c_L<F(L)\),
it is greater than 1. QED.

This does not impose additional vanishing moments on the full
RH test family. Every test is decomposed into its observed part
and this tail; the observed part is retained below.

## ASM-T2. The cut and its coupling are defined by the actual source

The zero extensions of the polynomials (6) belong to D_L and
to the operator domain of \(\mathcal T_L\). For any polynomial v
on I, the actual operator column is, for almost every x in I,
\[
\boxed{
\begin{aligned}
(\mathcal T_Lv)(x)
={}&\cosh(x/2)\ell(v)+\frac{\Omega_0}{2}v(x)\\
&+\frac12\int_I k(|x-y|)\bigl(v(x)-v(y)\bigr)\,dy\\
&+\frac12\bigl(F(L-x)+F(L+x)\bigr)v(x)\\
&-\frac12\sum_{n\in\mathcal N_L}a_n
\bigl(v(x-\log n)+v(x+\log n)\bigr).
\end{aligned}
}
\tag{16}
\]
Values outside I in the last line are zero. In particular,
every column in (16) is in H_L.

The even polynomials are a form core for D_L.

**Proof.** On I a polynomial is bounded and Lipschitz. Since
\(k(u)=O(1/u)\) near zero, the difference in the integral in (16)
cancels that singularity. The integral is bounded on I.
The function F(v) has only a logarithmic singularity at v=0,
which is square integrable. The boundary multiplier times v
therefore belongs to L2. The remaining terms are bounded or finite
translations of bounded compactly supported functions.

For domain membership, approximate a zero-extended polynomial
by even smooth cutoffs inside I. The difference r_epsilon is
supported in boundary strips of total length O(epsilon), has
uniformly bounded supremum and total variation, and satisfies
\[
\|r_\epsilon\|_2^2=O(\epsilon),\qquad
\|r_\epsilon-T_ur_\epsilon\|_2^2
\le C\min\{u,\epsilon\}.
\]
The second inequality follows from the L1 translation bound
for bounded-variation functions and from the support bound.
Consequently
\[
\mathcal E_\Gamma[r_\epsilon]
=O\bigl(\epsilon(1+|\log\epsilon|)\bigr)\longrightarrow0.
\tag{17}
\]
This places the polynomial in D_L. Polarizing (11) and the
finite prime terms gives (16) against every smooth test.
Its right side is in L2, so form continuity and core density
give the operator-domain assertion.

Conversely, an even smooth test can be approximated by even
polynomials in C1([-L,L]). One construction uniformly approximates
its derivative by polynomials, integrates, adjusts the constant
and then symmetrizes. For the difference w, the interior part of
(11) is bounded by
\[
C_L\|w'\|_\infty^2
\]
because \(\int_0^{2L}u^2k(u)\,du<\infty\).
The boundary part is bounded by \(C_L\|w\|_\infty^2\), since
F is integrable at zero. These bounds give convergence in the
form norm. Density of the smooth tests proves the last claim. QED.

Thus \(P_ND_L\subset D_L\) and \(Q_ND_L\subset D_L\).
This is a lawful domain-preserving cut, rather than an arbitrary
discarding of source coordinates.

## ASM-T3. Exact finite memory matrix and square remainder

Fix N with \(d=d_{L,N}>0\). Restrict q_L to
\(D_L\cap Q_NH_L\). Its associated positive self-adjoint
operator on \(Q_NH_L\) is denoted \(C_N\). By ASM-T1,
\[
C_N\ge dI,\qquad \|C_N^{-1}\|\le d^{-1}.
\tag{18}
\]
Let \(U_N:\mathbb R^N\to V_N\) send the standard basis to (6).
Define
\[
\begin{gathered}
(A_N)_{ij}=q_L(\phi_i,\phi_j),\qquad 0\le i,j<N,\\
B_N=Q_N\mathcal T_LU_N:\mathbb R^N\to Q_NH_L,\\
\boxed{\mathsf S_{L,N}=A_N-B_N^*C_N^{-1}B_N.}
\end{gathered}
\tag{19}
\]
All coefficients and columns here come from (2) and (16).

For \(g=U_Na+h\), \(h\in Q_ND_L\), one has the exact identity
\[
\boxed{
q_L[g]=a^{\mathsf T}\mathsf S_{L,N}a+
\left\|C_N^{1/2}
\left(h+C_N^{-1}B_Na\right)\right\|_2^2.
}
\tag{20}
\]
Consequently
\[
\boxed{
q_L\ge0\text{ on }D_L
\quad\Longleftrightarrow\quad
\mathsf S_{L,N}\succeq0.
}
\tag{21}
\]
The number of negative directions, with multiplicity, is exactly
\[
n_-(q_L)=n_-(\mathsf S_{L,N})\le N.
\tag{22}
\]

**Proof.** The restricted form is closed and densely defined on
\(Q_NH_L\). This follows from domain preservation, continuity of
the finite projection in form norm, and polynomial core density.
Its positive lower bound gives (18).

ASM-T2 gives the bounded finite-column map B_N. Polarization yields
\[
q_L[U_Na+h]=a^{\mathsf T}A_Na
+2\langle B_Na,h\rangle+q_L[h].
\]
Since \(C_N^{-1}B_Na\in\operatorname{Dom}C_N\), completing the
square in the closed tail form proves (20). No inverse of a
possibly indefinite full source is used.

The change of variables
\[
(a,h)\longmapsto(a,h+C_N^{-1}B_Na)
\]
is an invertible linear map of the form domain. Its inverse
subtracts the same vector. The second summand in (20) is positive
definite. This proves (21) and the equality of negative inertias
in (22). QED.

The subtraction
\[
\mathsf M_{L,N}:=B_N^*C_N^{-1}B_N\succeq0
\tag{23}
\]
is the exact returned memory. Testing A_N alone is insufficient:
the source seen after minimization is \(A_N-\mathsf M_{L,N}\).

## ASM-T4. Potential on an observation fibre and coherent nested cuts

Define the unresolved fibre and its variational potential by
\[
\mathfrak Z_{L,N}(a)=
\{g\in D_L:P_Ng=U_Na\},
\qquad
\Phi_{L,N}(a)=\inf_{g\in\mathfrak Z_{L,N}(a)}q_L[g].
\tag{24}
\]
Then
\[
\boxed{
\Phi_{L,N}(a)=a^{\mathsf T}\mathsf S_{L,N}a,
\qquad
g_{\min}(a)=U_Na-C_N^{-1}B_Na.
}
\tag{25}
\]
The minimizing representative is unique.

This is an operational meaning of unresolved potential: a cut
specifies a fibre of states with the same reading. The source
and the declared minimization select its effective potential.
An unobserved component is retained through its minimizing return,
not identified with the scalar zero.

For M>N with \(d_{L,N}>0\), identify the first N coordinates of
\(\mathbb R^M\) with a. Then
\[
\boxed{
\Phi_{L,N}(a)
=\inf_{b\in\mathbb R^{M-N}}\Phi_{L,M}(a,b).
}
\tag{26}
\]
Equivalently, writing the exact memory matrix at M as
\[
\mathsf S_{L,M}=
\begin{pmatrix}E&F^*\\ F&D\end{pmatrix},
\]
one has \(D\ge d_{L,N}I\) and
\[
\boxed{\mathsf S_{L,N}=E-F^*D^{-1}F.}
\tag{27}
\]

**Proof.** Equation (20) makes the infimum in (24) attained
exactly when \(h=-C_N^{-1}B_Na\), proving (25).

The unresolved space at N decomposes into the finite modes
N through M-1 and the unresolved space at M. Taking the infimum
over this Cartesian product in either order gives (26).
For a=0, every such state lies in \(Q_ND_L\), so ASM-T1 gives
\[
q_L[U_{N:M}b+h]\ge
d_{L,N}(\|b\|^2+\|h\|^2)\ge d_{L,N}\|b\|^2.
\]
Infimizing over h proves \(D\ge d_{L,N}I\).
Finite-dimensional completion of the square gives (27). QED.

This is a precise refinement invariance: direct and successive
lawful eliminations preserve the same source potential. It is
not the numerical identification \(0=\infty\). It also does not
supply the missing sign of \(\Phi_{L,N}\).

## ASM-T5. Finite trial matrices and an outward memory residual

For M>N, form the actual M by M source matrix
\[
A_M=(q_L(\phi_i,\phi_j))_{0\le i,j<M}
=\begin{pmatrix}A_N&E_M^*\\ E_M&D_M\end{pmatrix}.
\tag{28}
\]
Here \(D_M\ge d_{L,N}I\) by ASM-T1. Define the finite trial matrix
\[
\mathsf S_{L,N}^{(M)}
=A_N-E_M^*D_M^{-1}E_M.
\tag{29}
\]
Let \(U_{N:M}\) synthesize the modes N through M-1, and put
\[
Z_Ma=U_Na-U_{N:M}D_M^{-1}E_Ma,\qquad
R_M=Q_M\mathcal T_LZ_M .
\tag{30}
\]
The columns of Z_M are polynomials; the columns of R_M are
well-defined L2 functions by (16). Their Gram \(R_M^*R_M\)
is a finite N by N matrix with exact arithmetic-source entries.

Then
\[
\boxed{
\mathsf S_{L,N}^{(M)}-\mathsf S_{L,N}
=R_M^*C_N^{-1}R_M\succeq0,
}
\tag{31}
\]
and therefore
\[
\boxed{
\mathsf S_{L,N}^{(M)}
-\frac1{d_{L,N}}R_M^*R_M
\ \preceq\ \mathsf S_{L,N}\
\preceq\ \mathsf S_{L,N}^{(M)}.
}
\tag{32}
\]

There are two usable conclusions:

- If the left matrix in (32) is positive semidefinite, the entire
  window form q_L is positive semidefinite, including all unresolved
  layers.
- If \(a^{\mathsf T}\mathsf S_{L,N}^{(M)}a<0\), the explicit trial
  polynomial Z_Ma has negative actual source form. Form-core
  approximation then gives a smooth compact negative test.

Neither conclusion has been evaluated for a new window in this note.

**Proof.** The finite normal equation gives
\[
q_L(Z_Ma,u)=0
\quad(u\in V_M\cap Q_NH_L).
\]
Thus the Q_N part of \(\mathcal T_LZ_Ma\) is precisely R_Ma.
Also \(q_L[Z_Ma]=a^{\mathsf T}\mathsf S_{L,N}^{(M)}a\).
For any \(h\in Q_ND_L\),
\[
q_L[Z_Ma+h]
=a^{\mathsf T}\mathsf S_{L,N}^{(M)}a
+2\langle R_Ma,h\rangle+q_L[h].
\]
Complete the positive C_N square and minimize over h.
Since \(P_N(Z_Ma+h)=U_Na\), its minimum is (25).
This proves (31). Apply (18) to obtain (32).
The two conclusions follow from (20) and the trial identity. QED.

Moreover,
\[
\mathsf S_{L,N}^{(M+1)}
\preceq\mathsf S_{L,N}^{(M)},\qquad
\mathsf S_{L,N}^{(M)}\longrightarrow\mathsf S_{L,N}
\quad(M\to\infty)
\tag{33}
\]
in finite-matrix norm. Indeed the trial minimization spaces
increase, and their union is a form core for \(Q_ND_L\), by
ASM-T2 followed by Q_N. Their scalar infima converge to (25).
Polarization in the fixed N-dimensional space proves norm
convergence.

Equation (33) gives convergence without a rate. The residual
bound (32) supplies a quantitative enclosure whenever its
columns are enclosed. We do not infer that the unweighted
residual norms themselves tend to zero merely from form-core
density, and no such unproved rate is used.

## ASM-T6. Actual square lift and the exact remaining RH condition

For a fixed L and an N satisfying (8) with \(d_{L,N}>0\), a
finite matrix factor \(R_N\) with
\[
R_N^*R_N=\mathsf S_{L,N}
\tag{34}
\]
exists exactly when \(\mathsf S_{L,N}\succeq0\). In that case
\[
\boxed{
\mathcal J_{L,N}(U_Na+h)=
\left(
R_Na,\
C_N^{1/2}(h+C_N^{-1}B_Na)
\right)
}
\tag{35}
\]
is an actual square lift of the arithmetic source:
\[
\boxed{
q_L(g,h)=
\langle\mathcal J_{L,N}g,\mathcal J_{L,N}h\rangle.
}
\tag{36}
\]
The high-layer square in (35) is constructed unconditionally.
The remaining finite square (34) is the sign obligation.

**Proof.** The finite positive-semidefinite matrix factorization
gives (34) under its stated condition. Substitute it into (20)
and polarize. Conversely, any positive square representation of
the whole form makes it nonnegative; ASM-T3 then forces
\(\mathsf S_{L,N}\succeq0\). QED.

Choose once and for all a finite N(m) satisfying (10) for each
positive integer m. Then the complete remaining arithmetic
statement can be written
\[
\boxed{
\mathrm{RH}
\quad\Longleftrightarrow\quad
\mathsf S_{m,N(m)}\succeq0
\quad\text{for every integer }m\ge1.
}
\tag{37}
\]
Every compact smooth profile fits in an integer-width window.
Hence (37) follows from (4) and (21), without a uniform bound
on N(m).

The dimension and memory inverse depend on the window.
A finite collection of successful windows is not (37).
A finite trial matrix without the residual bound is not the
exact memory matrix. These distinctions prevent an infinite
layer from being silently erased by a finite cut.

## Result and scope

The actual-source development is now
\[
\boxed{
\begin{gathered}
\text{prime--Gamma coefficients with their exact scalar subtraction}\\
\longrightarrow\text{explicit positive unresolved tail}\\
\longrightarrow
q_L=\text{finite memory form}+\text{positive tail square}\\
\longrightarrow\text{coherent nested-cut potential}.
\end{gathered}
}
\tag{38}
\]

| Proved here | Remaining obligation |
|---|---|
| Explicit positive high-mode bound at every finite L | Sign of the retained low-mode matrix |
| Source columns include every active prime power and the Gamma boundary | Actual outward entry and residual enclosures for a chosen new window |
| Exact memory subtraction and finite inertia | All integer-width matrices in (37) are positive semidefinite |
| Successive lawful cuts preserve the variational potential | No positivity follows from refinement invariance alone |
| Explicit residual inequality for a finite trial | Form convergence alone does not give a residual-norm rate |
| Conditional source-matching square lift | The finite square (34) must be proved, not assumed |

The restricted TPW margin \(1/500\) for \(L\le11/20\) and the
unrestricted PTE margin \(1/16\) for \(L\le1/32\) are unchanged.
No larger positive arithmetic window is claimed here.

The Schur complement, closed-form square completion, Legendre
basis and variational projection are standard tools. The
programme-specific work is their explicit source-preserving
assembly, the all-width tail estimate (8), and the retained-memory
cut law (27) for the stated arithmetic form. No literature-wide
novelty claim is made.

Only this note and its README entry are added. All cited
cross-repository mathematical inputs are public and pinned to
RKF commit e8e3089745bf7dc0532d01c62a5f660f5658c6f0.
Private review material is not transferred. These are written
proofs, not newly machine-checked certificates. The documentation
commit contains [skip ci], is stacked on draft PR #112, and is
not merged.
