# RLO-1 — Relational observers, coincident readings and the RH sign

17 October 2026. Written continuation at Publications commit
d15ef36b2087714b7eb728076d613726abadef59, following
[PRW-1](RH_POLE_RECOVERY_WINDOW.md).

The proposed change of viewpoint is made concrete: the object is a
relation between actual source readings, with its observer coordinates
and their joint limiting rule retained. An apparent 0/0 can then have
a well-defined relational continuation. Two divergent readings can
have a finite, cut-dependent difference whose exact correction is
known. These are specified constructions, not a single assigned value
for every expression called 0/0 or infinity minus infinity.

For the actual xi source, a divided-difference matrix is the central
object. Its determinant divided by the squared observer-spacing
Vandermonde has a regular coincident-observer limit. Newton changes
of observers preserve its inertia, and an explicit error bound
transfers a strict finite sign between coincident and separated
observers.

At any one finite source anchor, positivity on one fixed sequence of
observer positions is equivalent to RH. At infinite anchor its
normalized determinant is explicitly positive at every fixed finite
size. The missing step remains uniform control of all sizes at the
actual finite anchor.

**Status: RH remains unproved.** This note establishes the relational
identities, their finite error bound, the fixed-sequence equivalence,
the explicit limiting determinant and the Gamma cut correction.
It does not establish the all-size inequality in that equivalence.
PRW's full even positive window and FCP's restricted window are
unchanged. No scientific computation or numerical source evaluation
is performed.

## 1. Retained source and the roles of the coordinates

Keep CWR, OQD and NIF's actual source:
\[
\begin{gathered}
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
a\ge a_0=3249/4,\\
\mu=(\gamma-i\delta)^2,\quad \delta=\beta-\tfrac12,\qquad
u_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},\\
\mathcal L_a[f]=\sum_\mu w_\mu f(u_\mu),\qquad
c_n(a)=\mathcal L_a[y^n],\\
q(a)>0,\qquad q'(a)<0,\qquad
0<|u_\mu|<1,\qquad
\sum_\mu |w_\mu|\le2,\qquad \sum_\mu w_\mu=1.
\end{gathered}
\tag{1}
\]
The sum is over distinct source points with their multiplicities.
The support is infinite and conjugation symmetric, and its only
accumulation point is zero. Its values and weights are not assumed
real and positive. Products inside L_a remain bilinear products.

The retained identity is
\[
q(t)=\sum_\mu\frac{\nu_\mu}{t+\mu}.
\tag{2}
\]
On the real observer range used below this converges absolutely,
with its local derivatives, as proved in the source lineage.

| Symbol | Role in this note |
|---|---|
| a | A fixed finite arithmetic anchor |
| y and u_mu | Source variable and its actual spectral values |
| x,z | Dimensionless observer positions relative to a |
| s=a(1+x), t=a(1+z) | The corresponding absolute source arguments |
| N | Number of observers or moment directions |
| h | Upper bound on observer positions approaching zero |
| epsilon | Separation cutoff in the Gamma integral |
| rho | Ratio of two Gamma cutoffs in RLO-T7 |

These roles are not interchangeable because the same numeral can
label them. The new objects are relations on this declared data.
Their algebraic laws and exact source matching remain part of the
definition.

For 0<=x<1 set
\[
F_a(x)=\frac{q(a(1+x))}{q(a)},\qquad
r_x(y)=\frac1{1+xy}.
\tag{3}
\]
Equations (1)--(2) give F_a(x)=L_a[r_x]. All these rational
observers are bounded on the closed source disk when x is bounded
strictly below one.

## RLO-T1. The divided reading is an exact source pairing, including coincidence

Define
\[
\mathsf K_a(x,z)=
\begin{cases}
\dfrac{F_a(x)-F_a(z)}{z-x},&x\ne z,\\[6pt]
-F_a'(x),&x=z.
\end{cases}
\]
Then
\[
\boxed{
\mathsf K_a(x,z)
=\mathcal L_a\!\left[
\frac{y}{(1+xy)(1+zy)}
\right]
=\mathcal L_a[y\,r_x r_z].
}
\tag{4}
\]
This is a real symmetric kernel for real observer positions.
Every individual entry is strictly positive on the stated range.

**Proof.** Subtract the two rational functions in (3):
\[
\frac{r_x(y)-r_z(y)}{z-x}
=\frac{y}{(1+xy)(1+zy)}.
\]
The absolute budget in (1) and the denominator bound justify
source evaluation, the coincident limit and differentiation.
Conjugation symmetry gives reality. For distinct positions,
q'(t)<0 makes the divided difference positive; the same strict
inequality gives positivity on the diagonal. QED.

Entrywise positivity does not establish positivity of its finite
matrices. The latter retains relations among all selected observers.

There is also an anchor-free reading:
\[
\mathsf K^{\rm abs}(s,t)=
\begin{cases}
\dfrac{q(s)-q(t)}{t-s},&s\ne t,\\[6pt]
-q'(s),&s=t,
\end{cases}
\qquad
\boxed{\mathsf K^{\rm abs}(s,t)=\frac{q(a)}a\mathsf K_a(x,z).}
\tag{5}
\]
At the same absolute positions, another admissible anchor chart
gives this same kernel. The chart changes only the stated positive
normalization factor. Keeping numerical x labels fixed while
changing a instead changes the actual arguments of q.

The 0/0 at x=z in the first line of (4) retains a specific
numerator--denominator relation. Its continuation is -F_a'(x),
not an arbitrary value assigned to the symbols 0/0.

## RLO-T2. Coincident observer matrices retain a finite determinant relation

Take distinct real nodes x_0,...,x_(N-1) in [0,h], h<1, and write
\[
\begin{gathered}
\mathsf K_N=(\mathsf K_a(x_i,x_j))_{0\le i,j<N},\qquad
V(\mathbf x)=\prod_{0\le i<j<N}(x_j-x_i),\\
b_j(y)=\frac{y^j}{\prod_{k=0}^j(1+x_k y)},\qquad
\mathsf G_N=(\mathcal L_a[y\,b_i b_j])_{0\le i,j<N}.
\end{gathered}
\tag{6}
\]
There is an explicit real invertible Newton observer matrix D:
\[
D_{jk}=
\begin{cases}
\dfrac{(-1)^j}
{\displaystyle\prod_{\substack{0\le\ell\le j\\\ell\ne k}}
(x_k-x_\ell)},&0\le k\le j,\\[9pt]
0,&k>j.
\end{cases}
\tag{7}
\]
Empty products equal one. It obeys
\[
\boxed{
\mathsf G_N=D\mathsf K_ND^T,\qquad
\mathfrak D_N(a;\mathbf x):=
\frac{\det\mathsf K_N}{V(\mathbf x)^2}
=\det\mathsf G_N.
}
\tag{8}
\]
As every node tends to zero,
\[
\boxed{
\mathsf G_N\longrightarrow
\mathsf A_N(a):=(c_{i+j+1}(a))_{0\le i,j<N},\qquad
\mathfrak D_N(a;\mathbf x)\longrightarrow\det\mathsf A_N(a).
}
\tag{9}
\]
The limit is independent of the relative rates at which the nodes
coalesce, for fixed N.

**Proof.** The elementary divided-difference identity is
\[
r_{\,\cdot}(y)[x_0,\ldots,x_j]
=\frac{(-y)^j}{\prod_{k=0}^j(1+x_k y)}.
\]
It follows either by subtracting two consecutive lower-order
formulas or directly by partial fractions. Multiplying by (-1)^j
gives b_j=sum_k D_(jk) r_(x_k). Bilinearity of the actual source
pairing (4) proves the congruence in (8).

The diagonal of D gives
\[
\det D=\frac{(-1)^{N(N-1)/2}}{V(\mathbf x)}.
\]
Taking determinants proves the second identity in (8). As all
x_k tend to zero, b_j tends uniformly to y^j on |y|<=1.
The budget (1) proves (9). QED.

For N>=2 the numerator and denominator in (8) both vanish when
all nodes coincide. Their joint relation retains (9).
The kernel matrix at exact coincidence alone has repeated rows;
the Newton observer information is what recovers the higher bank.

This determinant is normalized for the stated Newton coordinates.
It is not a numerical invariant under every arbitrary basis change.
Under any real invertible observer mixing M,
\[
\det(M\mathsf K_NM^T)=(\det M)^2\det\mathsf K_N,
\tag{10}
\]
and inertia is unchanged. Hence sign and inertia are the intrinsic
congruence data; the determinant's magnitude also carries a basis
volume.

## RLO-T3. A finite confluence error transfers strict signs

For the same nodes,
\[
\boxed{
\|\mathsf G_N-\mathsf A_N(a)\|
\le
2N\bigl((1-h)^{-2N}-1\bigr).
}
\tag{11}
\]
In particular, if h<=1/(4N),
\[
\boxed{\|\mathsf G_N-\mathsf A_N(a)\|\le8N^2h.}
\tag{12}
\]
If a known finite bank satisfies A_N(a)>=delta I with delta>0,
then
\[
h\le\min\left\{\frac1{4N},\frac{\delta}{16N^2}\right\}
\quad\Longrightarrow\quad
\mathsf G_N\ge\frac\delta2 I,\qquad \mathsf K_N\succ0.
\tag{13}
\]
If instead a unit vector v has v^T A_N(a)v<=-delta, the same
choice gives v^T G_N v<=-delta/2 and a negative direction of K_N.

**Proof.** Products of the geometric expansions for the denominators
in (6) give, uniformly on |y|<=1,
\[
\left|
\frac{y^{i+j+1}}
{\prod_{k=0}^i(1+x_k y)\prod_{\ell=0}^j(1+x_\ell y)}
-y^{i+j+1}
\right|
\le(1-h)^{-(i+j+2)}-1.
\]
Applying the absolute budget two, then the row-sum bound for a
symmetric N by N matrix, proves (11).
Bernoulli's inequality gives
\[
(1-h)^{-2N}-1
\le\frac{2Nh}{1-2Nh}\le4Nh
\]
when h<=1/(4N), proving (12). The two strict-sign conclusions
follow from the quadratic-form norm bound and congruence (8).
QED.

This is an explicit transfer bound for genuine separated readings,
not only a limiting correspondence. A positive finite bank supplies
a nearby positive observer configuration. A negative direction, if
present, survives in sufficiently tight observer configurations.

Exact confluence is not a claim of numerical immunity. If an
entry-enclosed matrix is K_N+E, its Newton error is DED^T, so
\[
\|DED^T\|\le\|D\|^2\|E\|.
\tag{14}
\]
Small spacings can amplify raw reading errors. Certified use must
retain those errors; no such numerical evaluation is run here.

## RLO-T4. One fixed sequence of relational observers retains the entire RH condition

At one fixed finite a>=a_0, take
\[
x_j=2^{-(j+1)},\qquad j=0,1,2,\ldots,
\]
and use the first N nodes for K_N and the normalized determinant
in (8). Then
\[
\boxed{
\begin{aligned}
\mathrm{RH}
&\Longleftrightarrow
\mathsf A_N(a)\succeq0\quad\hbox{for every }N\\
&\Longleftrightarrow
\mathsf K_N\succ0\quad\hbox{for every }N\\
&\Longleftrightarrow
\mathfrak D_N(a;x_0,\ldots,x_{N-1})>0
\quad\hbox{for every }N.
\end{aligned}
}
\tag{15}
\]
Under RH all the A_N are also strictly positive.
If RH fails, some fixed finite observer count has a negative
direction in every sufficiently late consecutive block of this
sequence. Those blocks can lie arbitrarily close to x=0.

**Proof of the polynomial condition.** The matrices A_N are
the Gram banks of the form
\[
B_a^+(p,h)=\mathcal L_a[yph].
\tag{16}
\]
Under RH the actual u_mu belong to (0,1) and the weights
w_mu u_mu are positive. Every nonzero polynomial is nonzero on
some of the infinite distinct support, so B_a^+(p,p)>0.

Conversely suppose RH fails. Choose a nonreal conjugate support
pair z,conjugate(z). Let P be the real polynomial vanishing at all
other support values of modulus at least s=|z|. There are only
finitely many such points. The remaining support has modulus at
most s_0<s. On the two polynomials y^m P and y^(m+1) P, the
B_a^+ Gram matrix is
\[
\begin{pmatrix}
2\Re A_m&2\Re(A_m z)\\
2\Re(A_m z)&2\Re(A_m z^2)
\end{pmatrix}
+O(s_0^{2m}),\qquad
A_m=w_z z^{2m+1}P(z)^2\ne0.
\tag{17}
\]
The leading determinant is
\[
-4|A_m|^2(\Im z)^2.
\]
After division by |A_m|^2 the error tends to zero. Thus some
finite polynomial bank has a negative direction, contradicting
positivity of every A_N. This is OQD's finite detector with the
explicit additional nonzero weight u_mu; it does not replace
complex products by absolute squares.

**Proof for separated observers.** Under RH, (4) is a positive
Gram kernel. The functions r_(x_j) at distinct nonzero positions
are linearly independent on the actual support. Indeed a rational
combination vanishing there has infinitely many zeros accumulating
at zero, where its denominators do not vanish. Its numerator is
identically zero. Its distinct poles -1/x_j then make every
coefficient zero. Thus every K_N is strictly positive.

Conversely, positivity of every prefix makes every finite
consecutive block positive. Fix a size m and take nodes
x_n,...,x_(n+m-1), then let n tend to infinity. Equation (8)
preserves positivity, and (9) gives A_m(a)>=0. The polynomial
condition already proved now implies RH.

Finally the Vandermonde square is positive for every finite
distinct-node bank. Positivity of every normalized determinant
in (8) is positivity of every leading principal determinant
of the nested matrices K_N.
Sylvester's criterion gives their strict positive definiteness.
This proves all equivalences.

If RH fails, (17) supplies A_m with a negative unit-vector
reading -delta for some m and delta>0. Its late blocks have
h=x_n tending to zero. RLO-T3 preserves that negative sign
in every sufficiently late block. QED.

This theorem concerns the actual q, not a freely chosen positive
interpolating kernel. Its universal size condition is the remaining
arithmetic obligation. No finite prefix or coordinate identity
establishes that condition.

## RLO-T5. The infinite-anchor relational kernel has an explicit positive determinant

For every fixed finite collection of distinct nodes in [0,1/2],
NIF's actual-source limit gives
\[
F_a(x)\longrightarrow F_*(x)=(1+x)^{-1/2}
\quad(a\to\infty).
\]
Write alpha_i=sqrt(1+x_i). Then
\[
\boxed{
\mathsf K_a(x,z)\longrightarrow
\mathsf K_*(x,z)
=\frac1{\sqrt{1+x}\sqrt{1+z}
       (\sqrt{1+x}+\sqrt{1+z})}.
}
\tag{18}
\]
The limiting matrix is strictly positive and
\[
\boxed{
\lim_{a\to\infty}\mathfrak D_N(a;\mathbf x)
=
\frac1{
2^N\displaystyle\prod_i\alpha_i^3
\displaystyle\prod_{i<j}(\alpha_i+\alpha_j)^4
}>0.
}
\tag{19}
\]
Consequently, for every fixed finite prefix N in (15), there
exists a finite threshold a_N after which that prefix is
strictly positive. No threshold uniform in N is asserted.

For fixed N there is also the joint limit
\[
\boxed{
\lim_{\substack{a\to\infty\\ \max_i x_i\to0}}
\mathfrak D_N(a;\mathbf x)
=2^{-N(2N-1)}.
}
\tag{20}
\]
Nodes remain distinct before taking the extended determinant
limit. Both orders of these two limits give (20) for fixed N.

**Proof.** NIF-T1 gives locally uniform analytic convergence
F_a to F_* on |x|<1. Cauchy's formula gives the same for first
derivatives on smaller disks. The identity
\[
\mathsf K_a(x,z)
=-\int_0^1F_a'(x+t(z-x))\,dt
\]
then proves uniform convergence on the stated real node range,
including the diagonal. Subtracting F_* values yields (18).

The kernel has the positive representation
\[
\mathsf K_*(x_i,x_j)
=\int_0^\infty
\frac{e^{-\alpha_i t}}{\alpha_i}
\frac{e^{-\alpha_j t}}{\alpha_j}\,dt .
\tag{21}
\]
Distinct exponential rates are linearly independent, giving
strict positivity. The Cauchy determinant identity gives
\[
\det\mathsf K_*
=
\frac{\prod_{i<j}(\alpha_j-\alpha_i)^2}
{2^N\prod_i\alpha_i^3
 \prod_{i<j}(\alpha_i+\alpha_j)^2}.
\tag{22}
\]
For completeness, the identity
det(1/(u_i+v_j))=
prod_(i<j)(u_j-u_i)(v_j-v_i)/prod_(i,j)(u_i+v_j)
follows by multiplying through the common denominator:
the numerator is alternating in each variable family and has
the corresponding Vandermonde factors. Comparing a pole residue,
starting at size one, fixes its multiplicative constant to one.
Set u_i=v_i=alpha_i and include the diagonal factors alpha_i^-1
on both sides to obtain (22).

Since x_j-x_i=(alpha_j-alpha_i)(alpha_j+alpha_i),
division by V(x)^2 proves (19). Strict positivity persists
for each fixed finite matrix under sufficiently small perturbations,
proving the threshold statement.

To justify the joint limit, use the uniform-in-a bound (12):
G_N-A_N(a) tends to zero whenever h tends to zero at fixed N.
NIF's fixed-order moment convergence gives A_N(a) tending to
its shifted arcsine moment matrix. Letting all alpha_i tend to
one in (19), or taking the confluence of (22), gives its
determinant 2^(-N(2N-1)). This proves (20). QED.

The quantifiers matter:
\[
\forall N\ \exists a_N\ \forall a\ge a_N
\]
is established here. RH in (15) asks for every N at the same
chosen finite a. The positive limiting determinants become
small with N; their positivity alone supplies no uniform error
bound at growing bank size.

Thus the infinite-anchor reading and the zero-spacing reading
can be joined as a rigorous relation for each fixed bank.
They are not a proof of the missing all-size finite-anchor sign.

## RLO-T6. Angles are valid probability coordinates precisely on a positive pair

For two chosen real observers with Gram matrix
\[
G=\begin{pmatrix}A&B\\B&C\end{pmatrix},\qquad A,C>0,
\]
define the normalized cross-reading R=B/sqrt(AC). Then
\[
\boxed{
G\succeq0
\Longleftrightarrow |R|\le1
\Longleftrightarrow
R=\cos\theta\text{ for some }\theta\in[0,\pi].
}
\tag{23}
\]
When these equivalent conditions hold,
\[
p_{\parallel}=R^2=\cos^2\theta,\qquad
p_{\perp}=1-R^2=\sin^2\theta,\qquad
\frac{\det G}{AC}=\sin^2\theta .
\tag{24}
\]

**Proof.** With positive diagonal entries, the two-dimensional
positive-semidefinite condition is AC-B^2>=0. Division by AC
proves (23), and (24) follows. In a positive Gram realization,
these probabilities are the squared projection and complementary
norms of one normalized observer vector onto the other's line.
QED.

The signed R, rather than R^2 alone, retains the relative real
orientation. In a complex two-channel representation the relative
phase is additional data, as retained in PRW.

For the individual rational observers of RLO-T1, A,C and B are
strictly positive without assuming RH. Whether every bank admits
a common positive Gram realization is the condition in (15).
Pairwise angular consistency alone is not substituted for all
finite matrix sizes.

Angles change under a general observer mixing. The determinant
sign and inertia obey (10) even when a chosen ratio chart becomes
singular. At A=0, a positive-semidefinite two-by-two form requires
B=0; an indefinite form can instead have A=0 with B nonzero.
The full matrix retains this distinction when a normalized
scalar reading would encounter division by zero.

## RLO-T7. The Gamma infinity-minus-infinity reading has an exact cut correction

Return to PTE/GCC's real even wave profile g in its closed
finite-energy domain D_L. Let
\[
N_g=\|g\|^2,\qquad
C_g(u)=\langle g,T_u g\rangle,\qquad
k(u)=\frac{e^{u/2}}{2\sinh u},\qquad
F(\epsilon)=\int_\epsilon^\infty k(u)\,du.
\]
Define two cut readings
\[
A_\epsilon[g]=N_gF(\epsilon),\qquad
B_\eta[g]=\int_\eta^\infty k(u)C_g(u)\,du .
\]
For nonzero g, each tends to positive infinity as its cutoff
tends to zero. The same-cut relation has the finite limit
\[
\boxed{
\lim_{\epsilon\downarrow0}
(A_\epsilon[g]-B_\epsilon[g])
=\mathcal E_\Gamma[g]
=\frac12\int_0^\infty k(u)\|g-T_u g\|^2\,du .
}
\tag{25}
\]
For any fixed cutoff ratio rho>0,
\[
\boxed{
\lim_{\epsilon\downarrow0}
(A_\epsilon[g]-B_{\rho\epsilon}[g])
=\mathcal E_\Gamma[g]+\frac{\log\rho}{2}N_g .
}
\tag{26}
\]
Therefore the corrected relational reading
\[
\boxed{
\lim_{\epsilon\downarrow0}
(A_\epsilon[g]-B_{\rho\epsilon}[g])
-\frac{\log\rho}{2}N_g
=\mathcal E_\Gamma[g]
}
\tag{27}
\]
is independent of that observer-cut ratio.

**Proof.** The identity
N_g-C_g(u)=||g-T_u g||^2/2 gives (25), with the limiting
domain justified by GCC-T1. Translation continuity gives
C_g(u) tending to N_g, while
k(u)=1/(2u)+O(1) as u tends to zero. Hence
\[
A_\epsilon-B_{\rho\epsilon}
=(A_\epsilon-B_\epsilon)
+\int_\epsilon^{\rho\epsilon}k(u)C_g(u)\,du ,
\]
where an integral with reversed endpoints has its signed value.
The last term tends to N_g log(rho)/2. Indeed its difference
from that elementary logarithmic integral is bounded by a
constant times sup|C_g(u)-N_g| on the shrinking interval,
plus a term tending to zero from the bounded kernel remainder.
This proves (26)--(27). The same local behavior proves the
individual divergences when N_g>0. QED.

The correction composes additively:
log(rho sigma)/2=log(rho)/2+log(sigma)/2.
The ratio therefore has an explicit transformation law, and
(27) is its invariant finite reading.

This is an actual instance of the proposed relational treatment
of infinity minus infinity. If the correction is omitted,
changing the relative cut shifts the source by
N_g log(rho)/2 and can change a claimed sign. Retaining (27)
keeps the original RH source exactly. No value is assigned
to an unpaired infinity-minus-infinity expression.

## 2. The relational target uses explicit arithmetic readings

All positive observer arguments in (15) are in the Euler
half-plane of the retained arithmetic identity. NIF equation (9)
gives the exact source formula
\[
\boxed{
q(t)=\frac1{t-1/4}
+\frac{\psi(1/4+\sqrt t/2)-\log\pi}{4\sqrt t}
-\frac1{2\sqrt t}
\sum_{m\ge2}\Lambda_{\rm vM}(m)m^{-1/2-\sqrt t}.
}
\tag{28}
\]
The square root is positive for t>0. The prime-power series and
its local derivatives converge absolutely in this range.
Thus the matrices use q values and first responses at explicit
real positions. They are not built by fitting unknown zeros
or choosing a positive replacement measure.

Equation (15), with (28) retained in every entry, is the concrete
new target:
\[
\boxed{
\mathfrak D_N(a_0;1/2,1/4,\ldots,2^{-N})>0
\quad\text{for every }N\ge1.
}
\tag{29}
\]
RLO-T1 proves every entry is positive. RLO-T5 proves every fixed
finite bank is positive sufficiently far out in anchor. RLO-T3
controls coalescing observer errors. None of these statements
yet proves the simultaneous all-N assertion at a_0.

A next proof step would be a source-specific domination of the
actual finite-anchor kernel by a positive reference on every
finite observer bank, or another bound establishing (29).
An entrywise bound alone would not suffice: it must control
all signed combinations of observer readings. The prime and
Gamma contributions in (28) have to remain jointly accounted for.

## 3. What the change of viewpoint establishes

| Proposed idea | Exact retained object and result |
|---|---|
| Numbers have different roles in one equation | Source values, anchor labels, observer positions, cutoffs and bank size are separately typed in (1)--(3) |
| A 0/0 reading can be resolved relationally | The divided reading extends by (4), and the determinant-spacing relation extends by (8)--(9) |
| Observer changes should preserve an invariant | The absolute kernel (5), congruence inertia (10), and corrected Gamma reading (27) specify the transformations and invariants |
| Probability can be read as an angle | Positive pair relations have the exact coordinates (23)--(24), retaining orientation |
| Infinity minus infinity can contain finite information | Matched Gamma cuts give (25); a changed relative cut carries the explicit logarithmic correction (26) |
| Infinite anchor and zero spacing may meet | Their fixed-N joint limit is the explicit positive value (20) |
| These relations should prove RH | They give the exact actual-source target (29); its all-size sign remains open |

The viewpoint is structural: retain the coupled data before
collapsing them to one numerical reading. This is compatible
with the native shunya potential and cut language. It does not
identify an invariant conservation identity with a positive form.
CLA's comparison-source audit already shows why that extra
arithmetic step is necessary.

## Dependencies, attribution and verification scope

- [PRW-1](RH_POLE_RECOVERY_WINDOW.md): angular cut coordinates
  and the full even positive window L<=1/5.
- [CWR-1](RH_CUT_WAVE_RETURN.md): actual resolvent source,
  absolute budget, positivity of q and strict decrease of q.
- [OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md): faithful source
  and finite conjugate-pair detector.
- [NIF-1](RH_INFINITY_CUT_FIXED_POINT.md): actual arithmetic
  formula and locally uniform infinite-anchor source limit.
- [GCC-1](RH_GAMMA_CUT_CERTIFICATES.md): correct closed
  Gamma domain and finite cut identity.
- [PTE-1](RH_PRIME_TRANSPORT_ENERGY.md): Gamma difference
  energy and the complete arithmetic source.
- [SPT-1](RH_SHUNYA_POTENTIAL.md) and
  [CLA-1](RH_CLOSURE_AUDIT.md): native potential invariance
  and its distinction from source positivity.
- [FCP-1](RH_FINITE_CUT_POSITIVITY.md): unchanged restricted
  positive window.

Divided differences, Gram congruence, Cauchy determinants,
Sylvester's criterion and compensated divergent expressions
are standard mathematical tools. The contribution here is
their explicit combination on the retained arithmetic source,
with the observer confluence bound, fixed-sequence criterion,
limiting determinant and exact Gamma cut correction.
No literature-wide novelty or completed RH proof is claimed.

Only this note and its README entry change. The work is stacked
on draft PR #117 with a [skip ci] documentation commit.
Verification consists of written derivation review, document
structure and links, exact immutable-head content and PR metadata.
No scientific computation, source replay, independent agent,
workflow dispatch or merge is performed. The PR remains a draft.
