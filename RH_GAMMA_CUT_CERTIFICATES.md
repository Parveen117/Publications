# GCC-1 — Gamma cuts and finite certificate completeness

17 October 2026. Written continuation at Publications commit
595df50d206990e562383e0e3b5aabb0513479af, following
[FCP-1](RH_FINITE_CUT_POSITIVITY.md).

FCP supplies a finite-cut source-error bound and a sufficient
memory-residual enclosure. ASM proves convergence of variational
trial values but does not infer convergence of their unweighted
operator residuals. The present step resolves that particular
convergence issue through an auxiliary cut of the positive Gamma
energy.

At every positive separation cutoff, the source is a bounded
operator. Its finite-mode memory residual then converges to zero.
As that separation cutoff closes, its spectral bottom increases
to the actual source bottom. This gives the following conditional
completeness statements at every fixed finite support width:

- A strictly positive actual source admits a finite positive
  certificate with all prime-record and memory errors included.
- A source with a negative direction admits a finite polynomial
  witness with an outward upper bound.
- Any two strict bounds around the actual spectral bottom admit
  finite certificates for the corresponding enclosure.

The exact-zero case is handled by a vanishing tolerance in the
all-width criterion:
\[
\boxed{
\mathrm{RH}
\Longleftrightarrow
\forall j,r\ge1\
\exists\text{ a finite positive certificate for }
q_j+\frac1r\|\cdot\|^2 .
}
\]
Here j and r are integers, q_j is the actual source in the support
window (-j,j), and the certificate is the explicit finite matrix
constructed below. Its data contain no infinite-dimensional inverse.

**Status: RH remains unproved.** The criterion does not produce
the certificates at every j and r. It proves completeness of the
finite certificate mechanism, not the required all-width signs.
No new positive window is claimed: FCP's restricted range
\(L\le111/200\) and the earlier unrestricted range remain unchanged.
These are written proofs, with no scientific computation or replay.

## 1. Retained source, domain and the auxiliary Gamma cut

Fix L>0, let I=(-L,L), and let H_L be the real Hilbert space of
even L2 functions on I, extended by zero to the full line. Retain
\[
\begin{gathered}
N_g=\|g\|_2^2,\qquad
T_ug(x)=g(x-u),\qquad C_g(u)=\langle g,T_ug\rangle,\\
\ell(g)=\int_I g(x)\cosh(x/2)\,dx,\qquad
\mathcal N_L=\{n\ge2:n\text{ is a prime power},\ \log n<2L\},\\
a_n=\frac{\Lambda(n)}{\sqrt n},\qquad
c_L=\sum_{n\in\mathcal N_L}a_n-\frac{\Omega_0}{2},\\
\Omega_0=-\gamma_E-\frac{\pi}{2}-3\log2-\log\pi,\\
k(u)=\frac{e^{u/2}}{2\sinh u},\qquad
F(v)=\int_v^\infty k(u)\,du
=\operatorname{artanh}(e^{-v/2})+\arctan(e^{-v/2}),\\
\mathcal E_\Gamma[g]=\frac12\int_0^\infty
k(u)\|g-T_ug\|_2^2\,du,\\
q_L[g]=\ell(g)^2+\frac{\Omega_0}{2}N_g+
\mathcal E_\Gamma[g]-\sum_{n\in\mathcal N_L}a_nC_g(\log n).
\end{gathered}
\tag{1}
\]
D_L is the closure of the even smooth compact tests in the
Gamma form norm. PTE constructs its closed source operator
\(\mathcal T_L\), with compact resolvent. Its spectral bottom is
\[
\lambda_L=\inf_{\substack{g\in D_L\\N_g=1}}q_L[g].
\tag{2}
\]
This is the retained window profile, not an assumed positive number.

For \(0<\epsilon<2L\), define
\[
\mathcal E_{\Gamma,\epsilon}[g]
=\frac12\int_\epsilon^\infty k(u)\|g-T_ug\|_2^2\,du
\]
and let \(q_{L,\epsilon}\) be (1) with this energy in place of
\(\mathcal E_\Gamma\). The prime terms, pole square and negative
Gamma scalar are unchanged.

The parameter \(\epsilon\) cuts a separation integral in an
already constructed positive source channel. It is distinct from
the geometric prime-record depth K and from the number of profile
modes retained by a Legendre cut.

## GCC-T1. The Gamma cut is bounded, ordered and has the correct limiting domain

For every \(\epsilon>0\), \(q_{L,\epsilon}\) is a bounded
self-adjoint quadratic form on all of H_L. As \(\epsilon\)
decreases to zero,
\[
\boxed{
q_{L,\epsilon}[g]\uparrow q_L[g]\quad(g\in D_L),\qquad
q_{L,\epsilon}[g]\le q_L[g].
}
\tag{3}
\]
Moreover,
\[
\boxed{
D_L=\{g\in H_L:\mathcal E_\Gamma[g]<\infty\}.
}
\tag{4}
\]
Thus the monotone limit does not introduce a larger, unverified
source domain.

For later finite bounds, put
\[
h_\epsilon(t)=
\begin{cases}
k(|t|),&|t|\ge\epsilon,\\
0,&|t|<\epsilon,
\end{cases}
\]
and define the even kernel on I times I by
\[
\mathcal K_\epsilon(x,y)
=\frac14\bigl(h_\epsilon(x-y)+h_\epsilon(x+y)\bigr).
\tag{5}
\]
Its operator \(K_\epsilon\) is Hilbert--Schmidt, and on H_L
\[
\boxed{
\mathcal E_{\Gamma,\epsilon}[g]
=F(\epsilon)N_g-\langle g,K_\epsilon g\rangle.
}
\tag{6}
\]

**Proof.** Since \(\|g-T_ug\|\le2\|g\|\),
\[
0\le\mathcal E_{\Gamma,\epsilon}[g]
\le2F(\epsilon)N_g.
\]
All other terms in (1) are bounded on H_L. This proves boundedness.
The difference between the two Gamma energies is
\[
\frac12\int_0^\epsilon k(u)\|g-T_ug\|^2\,du\ge0,
\]
and monotone convergence proves (3) on its finite-energy domain.

For even g, the two kernel terms in (5) have the same action.
Changing variables in the correlation integral gives
\[
\langle g,K_\epsilon g\rangle
=\int_\epsilon^{2L}k(u)C_g(u)\,du.
\]
The correlations vanish beyond 2L. Since
\(\|g-T_ug\|^2=2(N_g-C_g(u))\), this proves (6).
The kernel is bounded on the finite square after the separation
cutoff, hence is Hilbert--Schmidt. It annihilates odd functions,
so its full-square Hilbert--Schmidt norm is also its norm on
the even subspace.

It remains to justify (4). PTE's retained Fourier multiplier is
\[
\mathcal E_\Gamma[g]
=\frac1{2\pi}\int_{\mathbb R}m_\Gamma(t)|\widehat g(t)|^2\,dt,
\qquad
m_\Gamma(t)=\int_0^\infty k(u)(1-\cos(tu))\,du.
\tag{7}
\]
The identity holds for arbitrary finite-energy profiles by
Tonelli and Plancherel. The multiplier is nonnegative and continuous,
and PTE's digamma asymptotic gives
\[
1+m_\Gamma(t)\asymp1+\log(1+|t|).
\tag{8}
\]
In particular, the maximal full-line energy space is a closed
weighted Fourier L2 space.

Let g have finite energy and support in [-L,L].
For \(a<1\) tending to one, the unitary dilations
\(g_a(x)=a^{-1/2}g(x/a)\) converge to g in this energy norm.
Here is a direct density justification. For \(1/2\le a\le1\),
(8) bounds the dilation operators uniformly on the weighted
Fourier space. On smooth compactly supported Fourier profiles
their strong continuity is immediate. Those profiles are dense
in weighted L2, which proves strong continuity for every g.
The supports of \(g_a\) lie in [-aL,aL].

Convolve \(g_a\) with real even smooth compact mollifiers whose
supports are smaller than the remaining boundary distance.
Their Fourier multipliers are bounded in absolute value by one
and converge pointwise to one. Dominated convergence in (7)
therefore gives energy-norm convergence. The resulting functions
are real even smooth compact tests inside I.
A diagonal choice of a and mollifier scale approximates g in
the defining norm of D_L. The reverse inclusion follows from
closedness of the maximal energy norm. This proves (4). QED.

## GCC-T2. The bounded-source spectral bottoms converge to the actual bottom

Define
\[
\lambda_{L,\epsilon}
=\inf_{\substack{g\in H_L\\N_g=1}}q_{L,\epsilon}[g].
\]
Then
\[
\boxed{
\lambda_{L,\epsilon}\uparrow\lambda_L
\quad(\epsilon\downarrow0).
}
\tag{9}
\]
Consequently \(\lambda_L>0\) implies
\(\lambda_{L,\epsilon}>0\) for some strictly positive
separation cutoff. No rate is asserted.

**Proof.** Monotonicity and
\(\lambda_{L,\epsilon}\le\lambda_L\) follow from (3).
Let \(\epsilon_j\downarrow0\), and choose unit vectors \(g_j\)
with
\[
q_{L,\epsilon_j}[g_j]
\le\lambda_{L,\epsilon_j}+1/j.
\]
Since \(q_{L,\epsilon}[g]\ge
\mathcal E_{\Gamma,\epsilon}[g]-c_LN_g\),
their cut energies are uniformly bounded by a finite constant C.

Pass to a weakly convergent subsequence \(g_j\rightharpoonup g\).
For each fixed \(\delta>0\), eventually
\(\epsilon_j<\delta\), and
\(\mathcal E_{\Gamma,\delta}[g_j]\le C\).
Compactness of \(K_\delta\) and (6) give
\[
F(\delta)-\langle g,K_\delta g\rangle\le C,
\]
or equivalently
\[
F(\delta)(1-\|g\|^2)+
\mathcal E_{\Gamma,\delta}[g]\le C.
\tag{10}
\]
Both terms on the left are nonnegative. Because
\(F(\delta)\to\infty\) as \(\delta\downarrow0\),
one has \(\|g\|=1\). Weak convergence of unit vectors to a unit
vector is strong convergence in H_L.

Equation (10) and monotone convergence also give
\(\mathcal E_\Gamma[g]\le C\), so (4) places g in D_L.
For fixed \(\delta\), the cut form is bounded and therefore
continuous under this strong convergence. Hence
\[
q_{L,\delta}[g]
=\lim_jq_{L,\delta}[g_j]
\le\lim_j\lambda_{L,\epsilon_j}.
\]
Let \(\delta\downarrow0\), use (3), and then (2):
\[
\lambda_L\le q_L[g]
\le\lim_j\lambda_{L,\epsilon_j}\le\lambda_L.
\]
This proves (9). The positive-bottom consequence follows. QED.

The compact-kernel argument in (10) rules out loss of unit norm
into unresolved modes along this cut limit. Monotone pointwise
convergence of forms alone would not justify that conclusion.

## GCC-T3. A positive high-mode bound for the bounded cut has finite data

Retain the even orthonormal Legendre basis
\[
\phi_j(x)=\sqrt{\frac{4j+1}{2L}}P_{2j}(x/L),
\qquad
U_Na=\sum_{j=0}^{N-1}a_j\phi_j,\qquad
P_N=U_NU_N^*,\quad Q_N=I-P_N.
\tag{11}
\]
Choose
\[
\epsilon_0=
L\exp\!\left(
-2e^L\bigl(\max\{c_L-F(L),0\}+2\bigr)
\right).
\tag{12}
\]
Then \(0<\epsilon_0<L\) and
\[
F(\epsilon_0)-c_L\ge2.
\tag{13}
\]

For this fixed cutoff define
\[
\begin{aligned}
t_N^2={}&
\int_I\int_I|\mathcal K_{\epsilon_0}(x,y)|^2\,dx\,dy\\
&-\sum_{i,j<N}
|\langle\phi_i,K_{\epsilon_0}\phi_j\rangle|^2 .
\end{aligned}
\tag{14}
\]
This is nonnegative, decreases to zero, and equals
\(\|K_{\epsilon_0}-P_NK_{\epsilon_0}P_N\|_{\rm HS}^2\).
Choose a finite N for which \(t_N<1\), and put
\[
d=F(\epsilon_0)-c_L-t_N>1.
\tag{15}
\]
For every \(0<\epsilon\le\epsilon_0\) and \(h\in Q_NH_L\),
\[
\boxed{q_{L,\epsilon}[h]\ge d\|h\|^2.}
\tag{16}
\]

**Proof.** The elementary kernel bound from ASM,
\(k(u)\ge e^{-L}/(2u)\) for \(0<u\le2L\), gives
\[
F(\epsilon_0)\ge
F(L)+\frac{e^{-L}}2\log(L/\epsilon_0)
\ge c_L+2.
\]
This proves (13).

Finite rank projections converge in Hilbert--Schmidt norm for
a Hilbert--Schmidt operator. More explicitly, the squared norm
is the sum of all squared matrix coefficients in the complete
even basis. The finite square sums in (14) increase to that sum.
This proves its stated identities and the finite existence of N.

Drop the nonnegative pole square and bound each prime
correlation in absolute value by \(\|h\|^2\). Equations (1) and
(6) then give
\[
q_{L,\epsilon_0}[h]\ge
(F(\epsilon_0)-c_L)\|h\|^2-\langle h,K_{\epsilon_0}h\rangle.
\]
For \(P_Nh=0\), the last quadratic term has absolute value at most
\(t_N\|h\|^2\). This proves (16) at \(\epsilon_0\).
Monotonicity proves it for every smaller cutoff. QED.

All data in (14) are finite integrals. An outward upper bound on
\(t_N\), below one, suffices; the exact Hilbert--Schmidt norm need
not be evaluated symbolically. This is an existence construction
and a finite acceptance condition, with no efficiency claim.

## GCC-T4. An ordered prime cut gives a finite sufficient source certificate

For each active prime \(p<e^{2L}\), retain
\[
r_p=p^{-1/2},\quad
\kappa_p=\frac{1+r_p}{1-r_p},\quad
V_{p,K_p}=\sqrt{1-r_p^2}
\sum_{k=0}^{K_p}r_p^kT_{k\log p},
\]
and define
\[
\eta_{\mathbf K}
=\sum_{p<e^{2L}}\log p\,\kappa_p r_p^{K_p+1}.
\tag{17}
\]
The sum is over primes. Let
\[
\begin{aligned}
\widetilde q_{\epsilon,\mathbf K}[g]
={}&\ell(g)^2+\frac{\Omega_0}{2}N_g+
\mathcal E_{\Gamma,\epsilon}[g]\\
&+\frac12\sum_{p<e^{2L}}\log p
\left(N_g-\|V_{p,K_p}g\|^2\right),\\
\widehat q_{\epsilon,\mathbf K}[g]
={}&\widetilde q_{\epsilon,\mathbf K}[g]
-\eta_{\mathbf K}N_g .
\end{aligned}
\tag{18}
\]
These are bounded forms on H_L. FCP-T1 gives the ordered bounds
\[
\boxed{
q_{L,\epsilon}-2\eta_{\mathbf K}\|\cdot\|^2
\le\widehat q_{\epsilon,\mathbf K}
\le q_{L,\epsilon}\le q_L.
}
\tag{19}
\]
The last comparison is on D_L. Thus the finite prime-record
error is compensated in a lower source before certification.

Fix \(\epsilon_0,N,d\) from GCC-T3. Choose
\(0<\epsilon\le\epsilon_0\) and finite depths with
\(\eta_{\mathbf K}<d/2\). Put
\[
\delta=d-2\eta_{\mathbf K}>0,\qquad
\widehat T=\text{the bounded operator of }
\widehat q_{\epsilon,\mathbf K}.
\tag{20}
\]
Its restriction to \(Q_NH_L\) is at least \(\delta I\).
For \(M>N\), partition the finite source matrix as
\[
U_M^*\widehat T U_M=
\begin{pmatrix}A&E_M^*\\ E_M&D_M\end{pmatrix},
\qquad D_M\ge\delta I .
\tag{21}
\]
Here \(U_{N:M}\) synthesizes the modes \(N,\ldots,M-1\).
Define
\[
\begin{gathered}
Z_M=U_N-U_{N:M}D_M^{-1}E_M,\qquad
X_M=A-E_M^*D_M^{-1}E_M,\\
R_M=Q_M\widehat T Z_M,\qquad
\boxed{\mathsf L_M=X_M-\delta^{-1}R_M^*R_M.}
\end{gathered}
\tag{22}
\]
Then
\[
\boxed{
\mathsf L_M\succeq0
\quad\Longrightarrow\quad
q_L[g]\ge0\quad(g\in D_L).
}
\tag{23}
\]

**Proof.** The record error in FCP is unchanged by the common
bounded Gamma cutoff:
\[
|\widetilde q_{\epsilon,\mathbf K}[g]
-q_{L,\epsilon}[g]|\le\eta_{\mathbf K}N_g.
\]
Subtracting \(\eta_{\mathbf K}N_g\) proves (19).
Equations (16) and (19) prove the tail bound in (20), hence
the finite inverse in (21) exists.

Let \(C=Q_N\widehat T|_{Q_NH_L}\),
\(B=Q_N\widehat T U_N\), and
\(S=A-B^*C^{-1}B\).
The finite normal equations imply
\(Q_N\widehat T Z_M=R_M\). Completing the C square on the
observation fibre, exactly as in ASM and FCP, gives
\[
X_M-S=R_M^*C^{-1}R_M,\qquad
\mathsf L_M\preceq S\preceq X_M.
\tag{24}
\]
If \(\mathsf L_M\succeq0\), the exact square decomposition
makes \(\widehat q_{\epsilon,\mathbf K}\) nonnegative.
Equation (19) transfers this to the actual source. QED.

The inverse of C occurs only in the proof of the enclosure.
It is absent from the certificate data (14), (17), (21)--(22).
For example, the residual Gram can be written entirely with
finite columns:
\[
R_M^*R_M
=(\widehat T Z_M)^*(\widehat T Z_M)
-(U_M^*\widehat T Z_M)^*(U_M^*\widehat T Z_M).
\tag{25}
\]
The columns use the bounded kernel (5), finitely many translated
polynomials, the pole function and finite matrix operations.

A strict certificate also yields a source margin. If
\(\mathsf L_M\ge sI\) with \(s>0\), and \(b\ge\|B\|\), then
\[
\boxed{
q_L[g]\ge
\frac{\min\{s,\delta\}}{(1+b/\delta)^2}\|g\|^2.
}
\tag{26}
\]
Indeed the square decomposition is
\[
\widehat q[U_Na+h]
=a^*Sa+\|C^{1/2}(h+C^{-1}Ba)\|^2.
\]
The inverse of the map \((a,h)\mapsto(a,h+C^{-1}Ba)\)
has norm at most \(1+b/\delta\). Apply \(S\ge sI\) and
\(C\ge\delta I\). The finite Gram \(B^*B\) supplies b.

## GCC-T5. The bounded-source residual converges, and strict positivity is finitely certifiable

Fix \(\epsilon,\mathbf K,N\) satisfying (20). Then
\[
\boxed{
\|R_M\|\longrightarrow0,\qquad
X_M\longrightarrow S,\qquad
\mathsf L_M\longrightarrow S
\quad(M\to\infty)
}
\tag{27}
\]
in finite-column or finite-matrix operator norm.

Consequently,
\[
\boxed{
\lambda_L>0
\quad\Longrightarrow\quad
\exists\ \epsilon>0,\ N<M,\ \mathbf K<\infty
\text{ satisfying the construction with }\mathsf L_M\succ0 .
}
\tag{28}
\]

**Proof of convergence.** The exact minimizing lift is
\[
Z_\infty=U_N-C^{-1}B,\qquad
Q_N\widehat T Z_\infty=0.
\]
On \(Q_NH_L\), the C energy norm is equivalent to the L2 norm:
C is bounded and \(C\ge\delta I\).
The finite intermediate Legendre spaces are dense in this
space. Their Galerkin minimizers therefore satisfy
\(Z_M\to Z_\infty\) in L2 operator norm. Explicitly, the
squared C error is minimized by the Galerkin equations and is
bounded by the squared C error of any finite-mode approximation
to each column of \(C^{-1}B\). Those errors tend to zero.
There are only N columns.

The finite normal equations and boundedness of \(\widehat T\)
now give
\[
\|R_M\|
=\|Q_N\widehat T(Z_M-Z_\infty)\|
\le\|\widehat T\|\,\|Z_M-Z_\infty\|\longrightarrow0.
\]
Equation (24) gives \(X_M\to S\), and then (22) gives
\(\mathsf L_M\to S\). This proves (27).

**Proof of finite positive certification.** Construct
\(\epsilon_0,N,d\) without assuming a sign, using GCC-T3.
If \(\lambda_L>0\), GCC-T2 permits
\(0<\epsilon\le\epsilon_0\) with \(\lambda_{L,\epsilon}>0\).
Because the active prime set is finite, finite record depths
can make \(\eta_{\mathbf K}\) smaller than both \(d/2\) and
\(\lambda_{L,\epsilon}/4\). Equation (19) then gives
\[
\widehat T\ge
(\lambda_{L,\epsilon}-2\eta_{\mathbf K})I
\ge\frac{\lambda_{L,\epsilon}}2I.
\]
Infimizing on each observation fibre shows
\(S\ge\lambda_{L,\epsilon}I/2\).
By (27), \(\mathsf L_M\succ0\) for sufficiently large finite M.
This proves (28). QED.

This convergence argument uses boundedness after the Gamma cut.
It does not retroactively assume that the unweighted residuals
for ASM's original unbounded operator converge merely from
form-core density.

The implication (28) is conditional on the actual positive gap.
It gives neither a uniform cut size nor a positive gap at a
previously unresolved width. Its finite parameters can depend
on L and on how close its source is to zero.

## GCC-T6. Negative signs have finite outward witnesses, and strict spectral brackets are complete

Let v be a real even polynomial on I, extended by zero.
For \(0<\epsilon<2L\), put
\[
M_0=\|v\|_{L^\infty(I)},\qquad
M_1=\|v'\|_{L^\infty(I)},\qquad
\omega_\epsilon(v)=
e^{\epsilon/2}
\left(\frac{\epsilon M_0^2}{2}
+\frac{L\epsilon^2M_1^2}{4}\right).
\tag{29}
\]
Then
\[
\boxed{
q_L[v]\le
\widetilde q_{\epsilon,\mathbf K}[v]
+\eta_{\mathbf K}\|v\|^2+\omega_\epsilon(v).
}
\tag{30}
\]
A strictly negative right side is therefore an actual negative
source witness, with both cutoff errors controlled.

Conversely,
\[
\boxed{
\lambda_L<0
\quad\Longrightarrow\quad
\text{some finite polynomial v and finite cut data
make the right side of (30) strictly negative.}
}
\tag{31}
\]

**Proof.** For \(0<u<2L\), the two boundary strips of the
translated supports have total length 2u. On the overlap, the
mean-value bound gives
\[
\|v-T_uv\|^2
\le2uM_0^2+2Lu^2M_1^2.
\]
For \(u\le\epsilon\),
\(k(u)\le e^{\epsilon/2}/(2u)\).
Consequently the omitted Gamma energy is at most
\[
\frac12\int_0^\epsilon
k(u)(2uM_0^2+2Lu^2M_1^2)\,du
\le\omega_\epsilon(v).
\]
Combine this with the prime-record error to obtain (30).

If \(\lambda_L<0\), ASM's polynomial form core supplies v
with \(q_L[v]<0\).
For this fixed polynomial,
\(\omega_\epsilon(v)\to0\) as \(\epsilon\downarrow0\), and
\(\eta_{\mathbf K}\to0\) as the finite prime depths increase.
The right side of (30) then tends to \(q_L[v]\), so it is
negative for some finite choices. This proves (31). QED.

The sup norms in (29) may be replaced by elementary upper
bounds from the finite polynomial coefficients. A nonzero
polynomial witness can also be chosen with rational coefficients:
strict negativity is open in each finite coefficient space.

These results give two-sided spectral enclosure completeness.
For rational a,b with
\[
a<\lambda_L<b,
\]
apply GCC-T5 to \(q_L-a\|\cdot\|^2\), and GCC-T6 to
\(q_L-b\|\cdot\|^2\). Their finite certificates establish
\[
\boxed{a<\lambda_L<b.}
\tag{32}
\]
Every argument above is unchanged by a scalar source shift;
in GCC-T3 replace \(c_L\) by \(c_L+a\) for the first shifted
source, or by \(c_L+b\) for the second.

Thus arbitrarily tight strict brackets have finite source data.
This does not provide an efficient algorithm or decide whether
\(\lambda_L\) equals zero exactly. No such bracket has been
numerically evaluated in this note.

## GCC-T7. Vanishing tolerance gives an exact all-width finite-certificate criterion for RH

For each pair of positive integers j,r, consider the shifted
actual source
\[
q_j^{(r)}[g]=q_j[g]+\frac1r\|g\|^2.
\tag{33}
\]
Apply the complete construction (12)--(22) to this source.
Its scalar in the tail estimate is
\(c_j^{(r)}=c_j-1/r\); the Gamma kernel and prime coefficients
remain unchanged. Denote its certificate matrix by
\(\mathsf L_M^{(j,r)}\).

Then
\[
\boxed{
\begin{aligned}
\mathrm{RH}\quad\Longleftrightarrow\quad&
\text{for every }j,r\in\mathbb N_{\ge1},\\
&\text{there exist finite admissible }
\epsilon,N,M,\mathbf K\\
&\text{for which }\mathsf L_M^{(j,r)}\succ0 .
\end{aligned}
}
\tag{34}
\]
Admissible means that the explicit tail and record-error
conditions (12)--(20) hold for the shifted source.
All scalar integrals, matrix entries and residual Grams needed
to check one certificate are finite.

**Proof.** If RH holds, the retained CWF/PTE criterion makes
\(q_j\ge0\) at every integer width. Therefore
\(q_j^{(r)}\ge r^{-1}I\), with a strictly positive gap.
GCC-T5 gives a finite strict certificate for every j,r.

Conversely, a certificate in (34) and GCC-T4 imply
\[
q_j[g]\ge-\frac1r\|g\|^2\qquad(g\in D_j).
\]
Fix j and g, and let r tend to infinity. Then \(q_j[g]\ge0\).
Every smooth compact even test fits in some integer-width
window. Its source value is unchanged by zero extension:
additional disjoint prime shifts have zero correlation.
Thus positivity at every integer width implies the full
CWF/PTE criterion, hence RH. QED.

The tolerance is not a positive probability floor.
Each positive \(1/r\) is retained during its finite certificate,
and the universal quantifier over r recovers exact nonnegativity.
It also covers a hypothetical zero source eigenvalue, so (34)
does not assume strict positivity of the unshifted source under RH.

Strict inequalities in the tail tests and in the final matrix
permit outward finite entry enclosures: sufficiently tight
bounds preserve their signs. The theorem leaves those bounds
as explicit finite source-integral obligations. It does not
assert that a new evaluator, formal verifier or implementation
has been written or run.

## What the three cuts do

| Cut | Operation | Retained control |
|---|---|---|
| Gamma separation \(\epsilon\) | Omits a nonnegative short-separation energy | Lower source, with spectral bottom converging upward |
| Prime record \(\mathbf K\) | Truncates the normalized geometric records | Subtracts \(\eta_{\mathbf K}\) to keep the source order |
| Profile modes M | Eliminates intermediate modes and bounds the remaining memory | Explicit residual Gram and a positive tail lower bound |

The first two cuts produce the lower bounded source
\(\widehat q_{\epsilon,\mathbf K}\le q_L\).
The third tests that source through (22).
Positivity is established only when the resulting finite
matrix passes its stated inequality.

For a negative test, the omitted positive Gamma energy cannot
be dropped from an upper bound; that is why (29)--(30) retain
it. A negative value of a lower source by itself is not a
negative witness for the actual source.

## Established result and remaining arithmetic obligation

FCP supplied a sufficient finite enclosure. GCC adds a
convergent route to such enclosures whenever a fixed window
has a strict sign, and the tolerance criterion (34) handles
nonnegative endpoints without assuming a gap.

The all-width obligation is still substantive:
\[
\text{prove the finite certificate inequalities for every }
(j,r).
\]
Completeness of a certificate mechanism does not prove that
all those inequalities have the positive sign. A finite set
of successful windows or tolerances would not establish (34).
No actual negative source witness is asserted either.

The existing positive ranges are retained:

- On the same real even family with zero constant and
  \(\cosh(x/2)\) moments, FCP gives margin \(1/6000\)
  for \(L\le111/200\).
- TPW and BPC retain their stronger margins on their earlier
  restricted ranges.
- The established unrestricted range remains PTE's
  \(L\le1/32\).

No larger positive arithmetic window, uniform finite-bank
size, physical constant or RH closure follows from this note.

## Dependencies and verification scope

- [FCP-1](RH_FINITE_CUT_POSITIVITY.md): finite prime-record
  source error, combined source-memory enclosure and the
  retained restricted positive window.
- [PCM-1](RH_PRIME_CUT_MEMORY.md): geometric records,
  coherent prime cuts and their exact source normalization.
- [ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md): polynomial form
  core, source columns, exact Schur memory and variational trials.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): closed source form,
  Fourier multiplier growth, compact resolvent and spectral gate.
- [CWF-1](RH_CRITICAL_WAVE_FORMS.md): actual arithmetic
  source matching and the all-width RH criterion.
- [FPW-1](RH_FIRST_PRIME_WINDOW.md),
  [BPC-1](RH_BOUNDARY_PRIME_CAPS.md) and
  [TPW-1](RH_TWO_PRIME_WINDOW.md): the retained restricted lineage.

Monotone energy limits, Hilbert--Schmidt approximation, Galerkin
projection and Schur complements are standard mathematical tools.
The contribution here is their source-preserving assembly with
explicit cutoff directions, finite error data and the quantified
RH certificate criterion. No literature-wide novelty claim is made.

Only this note and its README entry are changed. The document is
stacked on draft PR #115 with a [skip ci] documentation commit.
Verification concerns written argument review, document structure,
repository links, immutable-head file equality and PR metadata.
No scientific computation, source evaluation, numerical spectrum,
independent agent, replay or workflow dispatch is performed.
The PR remains an unmerged draft.
