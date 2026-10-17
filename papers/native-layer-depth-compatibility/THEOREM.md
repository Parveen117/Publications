# LDC-1 — Stable observer records across changing source depth

17 October 2026. Continues LAA-1 in the same source track. Physical
interpretation is deliberately open. The question here is precise:
which readings survive a change of finite carrier, and how does the
cut itself change? No additional actual derivatives are needed for
this development; the retained twenty-layer jets supply nine nested
transitions at each of three anchors.

Retain c_n=(-a)^n q^(n)(a)/(n!q(a)), c_0=1, and the ILC forms
H_N=(c_(i+j)), A_N=(c_(i+j+1)), B_N=H_N-A_N. Assume the two
forms at order N+1 are strictly positive. This implies positivity
at every smaller order. All adjoints below use the acquired source
pairings, not an independently chosen Euclidean metric.

Carrier dimension N, derivative order n, binary record length m,
wave iteration and displacement lambda=(t-a)/a remain distinct.

## LDC-T1. Exact compression and one new boundary direction

Let E map a coefficient vector to the same polynomial in the next
carrier, E x=(x,0), and define

\[
\begin{split}
Y_N&=H_N^{-1}A_N,&E^*&=H_N^{-1}E^T H_{N+1},\\
\boldsymbol h&=(c_N,\ldots,c_{2N-1})^T,&
v&=\binom{-H_N^{-1}\boldsymbol h}{1},\\
h_N&=c_{2N}-\boldsymbol h^T H_N^{-1}\boldsymbol h
=\frac{\det H_{N+1}}{\det H_N},&\ell&=e_{N-1}.
\end{split}
\tag{1}
\]

The vector v is the monic degree-N polynomial orthogonal to the old
carrier. Its squared norm is h_N>0. Then

\[
\boxed{E^*E=I,\qquad E^*Y_{N+1}E=Y_N,\qquad
Y_{N+1}E-EY_N=v\ell^T.}
\tag{2}
\]

Consequently the squared-operator compression has a nonzero rank-one
defect

\[
\boxed{R_N:=E^*Y_{N+1}^2E-Y_N^2
=h_N H_N^{-1}\ell\ell^T\succeq0.}
\tag{3}
\]

Its sole nonzero eigenvalue is

\[
\boxed{\beta_N=\frac{h_N}{h_{N-1}}\in(0,1/4),\qquad h_0=1.}
\tag{4}
\]

This scalar is invariant under an isometric change of source
coordinates; it is not assigned a physical mass or energy meaning.

**Proof.** The top-left block of H_(N+1) is H_N, proving the
isometry. The top-left block of A_(N+1) is A_N, proving the first
compression identity. For j<N-1, Y_N e_j=e_(j+1) and the same
identity holds in the next carrier. For the last column, Y_(N+1)
maps y^(N-1) to y^N exactly; its component outside the old carrier
is the monic polynomial v. Orthogonality determines its old component
and proves (2). The Schur norm in (1) is positive and equals the
determinant ratio by block elimination.

Put P=EE*. The defect is
E*Y_(N+1)(I-P)Y_(N+1)E, whose quadratic form at x is
h_N |ell^T x|^2. This proves (3). The last diagonal entry of
H_N^(-1) is 1/h_(N-1), again by a Schur complement, so the sole
nonzero eigenvalue is (4). Since the next operator is a strict
contraction, Y_(N+1)-Y_(N+1)^2>0. Compression gives
R_N<Y_N-Y_N^2<=I/4 in the source pairing, proving the scalar bound.
For N>1, R_N is positive semidefinite of rank one, not positive
definite. QED.

## LDC-T2. Exact positive innovation of the full response

Let D_N(lambda)=H_N+lambda A_N, u_N=H_N e_0 and
f_N(lambda)=u_N^T D_N(lambda)^(-1)u_N. For real lambda>-1,

\[
\boxed{f_{N+1}(\lambda)-f_N(\lambda)=
\frac{\lambda^{2N}(\det H_{N+1})^2}
{\det D_N(\lambda)\det D_{N+1}(\lambda)}
=\frac{h_N\lambda^{2N}}
{\det(I+\lambda Y_N)\det(I+\lambda Y_{N+1})}.}
\tag{5}
\]

It is strictly positive away from lambda=0. At the source anchor
both responses equal one; their first difference is h_N lambda^(2N).
Thus the first 2N Taylor readings are preserved exactly, and the
first newly resolved even reading is governed by the same h_N that
measures the boundary direction in (3).

**Proof.** Write k=(c_(N+1),...,c_(2N))^T and
alpha=H_N^(-1)h. Orthogonality of the monic polynomial gives
k-A_N alpha=h_N ell: the first N-1 entries vanish and the last is
<y^N,v>=h_N. The cross block of D_(N+1) is b=h+lambda k.
Block elimination yields

\[
f_{N+1}-f_N=\frac{r_N^2}{s_N},\quad
r_N=c_N-b^T D_N^{-1}u_N,\quad
s_N=\frac{\det D_{N+1}}{\det D_N}>0.
\tag{6}
\]

Because b=D_N alpha+lambda h_N ell and alpha^T u_N=c_N,
r_N=-lambda h_N ell^T(I+lambda Y_N)^(-1)e_0. The first N-1
columns of Y_N are e_1,...,e_(N-1), so the relevant inverse
cofactor is triangular, giving

\[
\ell^T(I+\lambda Y_N)^{-1}e_0
=\frac{(-\lambda)^{N-1}}{\det(I+\lambda Y_N)},\qquad
r_N=\frac{(-\lambda)^N\det H_{N+1}}{\det D_N}.
\tag{7}
\]

Substitution into (6) proves (5). All response pencils are positive
because D_N=B_N+(1+lambda)A_N. Their determinants are positive,
proving the strict sign. Expansion at zero gives the leading term.
QED.

This is a finite identity for nested positive data. It does not say
that the actual infinite response lies above every finite response
without the required source hypotheses. Nor does monotonicity of a
finite list prove positivity of its unobserved next source layer.

For outward execution, evaluate (5) as
h_N (det H_N/det D_N)(det H_(N+1)/det D_(N+1)) lambda^(2N).
Forming the two extremely small determinant products first would
unnecessarily lose precision in the inherited fixed-scale intervals.
The factored expression uses the same exact identity and engine.

## LDC-T3. Stable word readings and consistent record deletion

At carrier N use the inherited commuting branch operators
K_1=sqrt(Y_N), K_0=sqrt(I-Y_N) and seed e=1. An ordered binary
word of length m with j ones has squared norm
<e,Y_N^j(I-Y_N)^(m-j)e>. Whenever m<=2N-1, this is exactly

\[
D_{j,m-j}=\sum_{r=0}^{m-j}(-1)^r\binom{m-j}{r}c_{j+r}.
\tag{8}
\]

Therefore every such word probability is unchanged in ANY larger
certified carrier from the same source. A particular record length m
requires N>=ceil((m+1)/2). This invariance concerns the specified
probability readings, not all state amplitudes or all observers.

Let J_m count ones and omega_(m,j)=binom(m,j)D_(j,m-j). Deleting
one uniformly chosen bit gives the exact consistency law

\[
\boxed{\omega_{m,j}=
\frac{m+1-j}{m+1}\omega_{m+1,j}
+\frac{j+1}{m+1}\omega_{m+1,j+1}.}
\tag{9}
\]

Every factorial observer reading is retained:

\[
\boxed{\mathbb E\frac{(J_m)_r}{(m)_r}=c_r\quad(0\le r\le m).}
\tag{10}
\]

For m>=1, with both m and m+1 available, set Z_m=J_m/m in this
deletion coupling. Then

\[
\boxed{\mathbb E[Z_m\mid J_{m+1}]=Z_{m+1},\qquad
\mathbb E[(Z_m-Z_{m+1})^2]=\frac{c_1-c_2}{m(m+1)}.}
\tag{11}
\]

**Proof.** ILC-T4 matches all moments through degree 2N-1, so the
polynomial in (8) has the exact source expectation. All ordered words
with the same count have the same probability and their multiplicity
is binom(m,j). The identity D_(j,k)=D_(j+1,k)+D_(j,k+1), followed
by the two binomial coefficient ratios, proves (9). Its coefficients
are exactly the probabilities of deleting a zero or a one. The
falling-factorial binomial identity from ILC-T2 gives (10).

Conditional on J_(m+1)=j, deletion leaves j ones with probability
(m+1-j)/(m+1) and j-1 with probability j/(m+1). Direct calculation
gives conditional mean j/(m+1) and conditional variance
j(m+1-j)/[m^2(m+1)^2]. Equation (10) gives
E[J_(m+1)(m+1-J_(m+1))]=m(m+1)(c_1-c_2), proving (11). QED.

In particular, E Z_m=c_1 and
E Z_m^2=c_2+(c_1-c_2)/m. These finite identities require only the
available positive record lengths. Under the additional ALL-ORDER
source positivity premise, ILC supplies a probability measure nu on
[0,1]. Choose y with law nu and, conditionally on y, independent
Bernoulli(y) bits. Their finite word probabilities are exactly (8),
and

\[
\boxed{\mathbb E[(Z_m-y)^2]=(c_1-c_2)/m.}
\tag{12}
\]

Thus this specified infinite observation process converges in L2
with an explicit squared-error rate. Equation (12) supplies a
quantitative version of the conditional ILC moment recovery; the
finite execution does not establish its all-order premise. It is
record refinement, distinct from norm convergence of iterated W^m.

## LDC-T4. Stable records do not require identical compressed cuts

Let Pi_N be ILC's two-channel cut, V=E direct_sum E, and
C=V*Pi_(N+1)V. Then C is a positive contraction and

\[
\boxed{C-C^2=V^*\Pi_{N+1}(I-VV^*)\Pi_{N+1}V\succeq0.}
\tag{13}
\]

For x in the old carrier,

\[
\langle(x,0),(C-C^2)(x,0)\rangle
\ge h_N|\ell^T x|^2.
\tag{14}
\]

Since h_N>0, C is not a projection and C!=Pi_N. In particular,
neither Pi_(N+1)V=V Pi_N nor W_(N+1)V=V W_N is possible under
strict growth with this natural polynomial inclusion.

**Proof.** Expand C-C^2 using Pi_(N+1)^2=Pi_(N+1) and V*V=I.
The result is the squared norm of (I-VV*)Pi_(N+1)V. On input
(x,0) its first-channel contribution is (I-EE*)Y_(N+1)Ex, whose
squared norm is exactly h_N|ell^T x|^2 by LDC-T1. The second
channel adds a nonnegative norm. Thus the defect is nonzero.
Exact cut intertwining would force C=Pi_N. Exact wave intertwining
would force Y_(N+1)E=EY_N from the wave's top-left block, again
contradicting (2). QED.

An exact two-mode witness has atoms 1/5,4/5 and weights 1/2,1/2.
Here Y_1=1/2, S_1=1/2 while S_2=(2/5)I. The compressed deeper
cut is [[1/2,2/5],[2/5,1/2]], with defect (9/100)I. Its compressed
wave has W*W=(16/25)I. Both full finite cuts are exact projections;
discarding the new carrier direction changes that property.

## LDC-T5. Actual execution and the next condition

The pinned LAA jets at original s=29,31,33 each have twenty layers
through order 19. They certify transitions N=1,...,9 to N+1. The
execution checks positive boundary norms, beta_N in (0,1/4), monic
orthogonality, and the determinant/Schur innovation at the held-out
original coordinate s=35. The finite responses strictly increase
with carrier dimension. Independent source errors already retained
by LAA provide a separate value check at s=35.

At each anchor, record lengths through 19 are positive and consistent
under deletion. Every available factorial reading recovers the same
source moment, and all eighteen adjacent squared-error identities
are checked with outward intervals. Thus the ten-direction carrier
retains every ordered word probability of length at most nineteen.
The three anchors remain different source distributions; no equality
of their c_r readings is asserted.

The new boundary norms certify failure of exact natural cut and wave
intertwining, by LDC-T4, without approximating any matrix square root.
This coexists with exact stabilization of the specified record
probabilities. The actual infinite-depth construction still requires
compatible positive source forms at every N, equivalently the
retained ILC all-order condition. Physical labels are left open;
the next step is selected by this source obligation and its exact
boundary/innovation law.
