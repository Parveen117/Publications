# SFA-1 — Actual source flow, carried invariants and the acceleration gate

18 October 2026 (Asia/Kolkata). Author: Monty.
Written continuation of APD-1 at
28bb5bb1f292938d1e1fb79096aa8113710c8176.

## Result and scope

The question is whether information-preserving flow itself supplies the
missing positivity, with a second reading retaining what a single
invariant omits. This note uses the existing actual anchor hierarchy
of LAT, rather than imposing a new positive exponential ensemble.

The actual mean and fluctuation have coupled flow equations. Their
finite positive interpretation is valid on the established source
range; the third central reading remains in the fluctuation equation.
For both source forms, an exact transport generator carries the
quadratic pairing invariantly while the generator is regular.

A first loss of finite-form positivity at a regular finite anchor
requires a pole in a specific compressed boundary reading. Its leading
coefficient is determined by the order of the determinant zero.
This gives a precise noncollapse target, not a proof excluding the pole.

Finally the old forms and their first derivatives do not determine the
next-depth sign. The missing odd moment is recovered exactly from the
second anchor derivative of the last old layer. LCW's sharp extension
window becomes a sharp two-sided acceleration inequality.

No all-depth sign theorem, new RH closure, or physical position/momentum
identification is claimed. APD's finite positivity ranges are retained.

## 1. The existing source determines the flow

Let
\[
\mathcal X(w^2)=\xi(1/2+w),\quad q=\mathcal X'/\mathcal X,\qquad
a(\tau)=a_*e^\tau,
\]
\[
c_j(a)=\frac{(-a)^j q^{(j)}(a)}{j!\,q(a)},\qquad c_0=1.
\tag{1}
\]
Work on real anchor intervals where q is regular and positive.
LAT-T1 proves the actual hierarchy
\[
\boxed{\dot c_j=(j+c_1)c_j-(j+1)c_{j+1}.}
\tag{2}
\]
A dot denotes differentiation in \(\tau=\log(a/a_*)\). It is
neither differentiation in the wave coordinate nor a new physical time.

For the source coordinates and normalized weights
\[
y_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},
\]
the same motion reads
\[
\boxed{\dot y_\mu=y_\mu(1-y_\mu),\qquad
\dot w_\mu=(c_1-y_\mu)w_\mu.}
\tag{3}
\]
Thus
\[
w_\mu(\tau)=w_\mu(0)
\exp\!\int_0^\tau(c_1(s)-y_\mu(s))\,ds.
\tag{4}
\]
These are actual, possibly complex, source weights. The exponential
factor in (4) does not turn them into a positive measure.
In the separately acquired finite positive model, LAT proves the
corresponding probability transport with its poles and residues
unchanged. That finite-model statement is not assumed for all xi zeros.

Equation (2), including its dependence on the next layer, is retained
from LAT. The results below concern its finite quadratic forms and
the information needed to extend them.

## SFA-T1. Mean motion and fluctuation have distinct coupled readings

Set
\[
m=c_1,\qquad v=c_2-c_1^2,\qquad
k_3=c_3-3c_1c_2+2c_1^3.
\]
Then the actual hierarchy gives
\[
\boxed{\dot m=m(1-m)-2v,\qquad
\dot v=2(1-2m)v-3k_3.}
\tag{5}
\]
In particular
\[
\boxed{v=\tfrac12\{m(1-m)-\dot m\}.}
\tag{6}
\]
The second reading measures a departure from the single-rate logistic
law. It is not asserted to be a conserved scalar.

**Proof.** Use (2) for j=1,2 and differentiate \(c_2-c_1^2\).
Substitute \(c_2=m^2+v\) and \(c_3=m^3+3mv+k_3\).
This gives (5), and its first equation gives (6). QED.

If the actual two-dimensional forms \(A_2,B_2\) are strictly positive,
their acquired finite probability model has distinct interior nodes.
It follows that \(0<m<1\) and \(v>0\). The two scalar Schur tests give
\[
\frac{c_2^2}{m}<c_3<
c_2-\frac{(m-c_2)^2}{1-m}.
\]
Equivalently,
\[
-mv+\frac{v^2}{m}<k_3<
(1-m)v-\frac{v^2}{1-m}.
\]
Therefore
\[
\boxed{
-(1+m)v+\frac{3v^2}{1-m}
<\dot v<
(2-m)v-\frac{3v^2}{m}.
}
\tag{7}
\]
In particular \(-2v<\dot v<2v\).
On any interval where those actual finite signs have been proved,
\[
\boxed{
e^{-2|\tau-\tau_0|}
\le\frac{v(\tau)}{v(\tau_0)}
\le e^{2|\tau-\tau_0|}.
}
\tag{8}
\]
This follows by integrating the logarithmic derivative bound.

Thus a strictly positive finite fluctuation has a quantitative
noncollapse bound on its established source interval. Flow by itself
does not force v to be nonzero: a one-atom positive model has v=0
and follows the first equation of (5) with nonconstant m.
The strict assertion here uses the two-direction source signs.
No assertion is made that every event must have positive probability.

## SFA-T2. Both actual source forms have an exact carried invariant

For N>=1 retain
\[
G_A=A_N=(c_{i+j+1})_{i,j<N},\qquad
G_B=B_N=(c_{i+j}-c_{i+j+1})_{i,j<N}.
\tag{9}
\]
Put \(r_A=1,\ r_B=0\), and define their boundary columns
\[
h_A=(c_{N+1},\ldots,c_{2N})^T,\qquad
h_B=(c_N-c_{N+1},\ldots,c_{2N-1}-c_{2N})^T.
\tag{10}
\]
Let e be the last coordinate vector. Define the N-by-N matrix D by
\[
D_{ii}=i\quad(0\le i<N),\qquad
D_{i+1,i}=-(i+1)\quad(0\le i<N-1),
\]
with every other entry zero. For \(K=A,B\),
\[
\boxed{
\dot G_K=(c_1+r_K)G_K+D^TG_K+G_KD
-N(eh_K^T+h_Ke^T).
}
\tag{11}
\]
Where \(G_K\) is invertible, set
\[
\boxed{
\Gamma_K=D-NG_K^{-1}h_Ke^T
+\tfrac12(c_1+r_K)I.
}
\tag{12}
\]
Then
\[
\boxed{\dot G_K=\Gamma_K^TG_K+G_K\Gamma_K.}
\tag{13}
\]

**Proof.** For the A entry with k=i+j+1, equation (2) gives
\[
\dot c_k=(i+j+1+c_1)c_k-(i+j+2)c_{k+1}.
\]
For a B entry, subtraction of adjacent instances of (2) gives
\[
\frac d{d\tau}(c_{i+j}-c_{i+j+1})
=(i+j+c_1)(c_{i+j}-c_{i+j+1})
-(i+j+2)(c_{i+j+1}-c_{i+j+2}).
\]
The internal raised indices are exactly the subdiagonal terms in D.
Each last raised index contributes its boundary column with
coefficient -N. This proves (11). Substitution of (12) proves (13).
QED.

Let U solve
\[
\dot U=-\Gamma_K U,\qquad U(\tau_0)=I.
\]
On a regular finite interval,
\[
\boxed{U(\tau)^TG_K(\tau)U(\tau)=G_K(\tau_0).}
\tag{14}
\]
Indeed differentiation cancels all three terms using (13).
A continuous, or integrable, generator gives an invertible fundamental
matrix U. Hence the inertia of the carried form is unchanged there.

This is the precise preserved pairing. U is the ordered exponential
of the varying generator; it is not replaced by a constant-generator
exponential when the generators do not commute.

Each channel has its own generator and carried matrix U; the two
generators need not coincide.

The regularity requirement is essential. Equation (12) contains an
inverse of the actual form and the next boundary layer. Writing
(13) does not prove those quantities remain controlled through a
candidate singularity.

There is also no hidden nonnegative source injection in (11).
For N>=2 on a range represented by a positive (N+1)-node interior
model, the first component of each h_K is positive. Thus h_K is not
parallel to e. The two nonzero eigenvalues of its boundary term are
\[
-N\left(e^Th_K+\|h_K\|\right),\qquad
-N\left(e^Th_K-\|h_K\|\right).
\tag{15}
\]
They have opposite signs. This follows by restricting the symmetric
rank-two matrix to the span of e and h_K. The exact boundary flux
is therefore indefinite even on this proved positive range.
Absorbing it into the carried generator is an identity, not a
proof of a new positive direction.

## SFA-T3. Finite positivity loss requires a boundary pole

While \(G_K\succ0\), define
\[
\chi_K=e^TG_K^{-1}h_K.
\]
Then
\[
\boxed{
\frac d{d\tau}\log\det G_K
=N(c_1+r_K)+N(N-1)-2N\chi_K.
}
\tag{16}
\]
Consequently
\[
\frac{\det G_K(\tau)}{\det G_K(\tau_0)}
=\exp\!\int_{\tau_0}^{\tau}
\left[N(c_1+r_K)+N(N-1)-2N\chi_K\right]ds.
\tag{17}
\]

**Proof.** Take the trace of \(G_K^{-1}\dot G_K\) in (11).
The two D terms give \(2\operatorname{tr}D=N(N-1)\),
and the two boundary terms each give \(-N\chi_K\). QED.

Suppose q remains regular and positive at a finite endpoint T,
the form is positive for \(\tau<T\), and \(\det G_K(T)=0\).
The entries are real analytic there. If the determinant zero has
order d>=1, then
\[
\boxed{
\chi_K(\tau)=\frac{d}{2N(T-\tau)}+O(1)
\quad(\tau\uparrow T).
}
\tag{18}
\]
In particular its integral diverges to positive infinity.

**Proof.** Since the analytic determinant is positive immediately
before T and is not identically zero, it has the expansion
\[
\det G_K(\tau)=C(T-\tau)^d(1+O(T-\tau)),\qquad C>0.
\]
Its logarithmic derivative is \(-d/(T-\tau)+O(1)\).
The drift terms in (16) remain bounded. Solving (16) for chi
gives (18). QED.

When approaching a loss endpoint T from above by decreasing the
anchor parameter, the corresponding formula is
\[
\chi_K(\tau)=-\frac{d}{2N(\tau-T)}+O(1).
\tag{19}
\]
Thus a suitable integrable bound on the actual boundary reading
would exclude finite-anchor collapse in the relevant direction.
No such bound at arbitrary depth is asserted here.

The pole in (18) or (19) belongs to the compressed transport reading,
which uses \(G_K^{-1}\). The underlying source q is regular by
hypothesis. This is not an identification of that pole with a
zero of xi.

**A quantitative bound on the established finite range.**
Let the actual (N+1)-dimensional forms be positive throughout an
anchor interval. Its finite interior model represents the needed
moments. For the channel weight \(\omega_A(y)=y\) or
\(\omega_B(y)=1-y\), put
\[
(S_K)_{ij}=\mathcal L_a[\omega_K(y)y^{i+j+1}],
\qquad Y_K=G_K^{-1}S_K.
\]
The first N-1 columns of \(Y_K\) are \(e_1,\ldots,e_{N-1}\);
the last is \(G_K^{-1}h_K\). Hence
\[
\boxed{\chi_K=\operatorname{tr}Y_K.}
\tag{20}
\]
In the acquired G_K pairing, the finite model gives
\(0<Y_K<I\). Therefore \(0<\chi_K<N\), and also \(0<c_1<1\).
Equation (16) yields
\[
\boxed{
\left|\frac d{d\tau}\log\det G_K\right|<N(N+1).
}
\tag{21}
\]
Integrating gives
\[
\boxed{
e^{-N(N+1)|\tau-\tau_0|}
\le\frac{\det G_K(\tau)}{\det G_K(\tau_0)}
\le e^{N(N+1)|\tau-\tau_0|}.
}
\tag{22}
\]
The strict contraction follows from the positive model's interior
support and polynomial independence; all the forms involved use
moments within its matching budget.

RBP proves the premise, for example, along anchor paths satisfying
\(N+1\le N_*(a)\) throughout. This is an actual finite transport
bound on that known region. It does not use the next-depth positivity
premise to claim it at an uncovered anchor.

## SFA-T4. The first jet misses the next odd layer; acceleration recovers it

Fix an anchor with \(A_N,B_N\succ0\), and write
\[
k=2N-1,\qquad x=c_{2N},\qquad z=c_{2N+1}.
\]
The two old forms and their first anchor derivatives depend only
on \(c_0,\ldots,c_{2N}\). They contain no z. In contrast, differentiating
(2) once more gives
\[
\ddot c_k=
\left[(k+c_1)^2+\dot c_1\right]c_k
-(k+1)(2k+1+2c_1)c_{k+1}
+(k+1)(k+2)c_{k+2}.
\]
Define the retained drift
\[
\mathfrak d_N=\left[(k+c_1)^2+\dot c_1\right]c_k
-2N(4N-1+2c_1)x,
\qquad
\dot c_1=(1+c_1)c_1-2c_2.
\]
Then
\[
\boxed{
\ddot c_{2N-1}=\mathfrak d_N+2N(2N+1)c_{2N+1}.
}
\tag{23}
\]
Thus the missing odd layer is exactly recoverable from the
acceleration of an old reading:
\[
\boxed{
c_{2N+1}=\frac{\ddot c_{2N-1}-\mathfrak d_N}{2N(2N+1)}.
}
\tag{24}
\]
The even layer is already recoverable from its velocity:
\[
c_{2N}=\frac{(2N-1+c_1)c_{2N-1}-\dot c_{2N-1}}{2N}.
\tag{25}
\]

Use the boundary columns (10) evaluated with this actual x and set
\[
\ell_N=h_A^TA_N^{-1}h_A,\qquad
u_N=x-h_B^TB_N^{-1}h_B.
\tag{26}
\]
The two actual next-depth forms are strictly positive if and only if
\[
\boxed{
\mathfrak d_N+2N(2N+1)\ell_N
<\ddot c_{2N-1}<
\mathfrak d_N+2N(2N+1)u_N.
}
\tag{27}
\]

**Proof.** Differentiating (2) gives (23)--(25).
The enlarged matrices are
\[
A_{N+1}=\begin{pmatrix}A_N&h_A\\h_A^T&z\end{pmatrix},\qquad
B_{N+1}=\begin{pmatrix}B_N&h_B\\h_B^T&x-z\end{pmatrix}.
\]
Their Schur innovations are \(z-\ell_N\) and \(u_N-z\).
Substitute (24). This gives both directions of (27). QED.

Equation (27) is LCW's sharp odd-layer gate expressed as an
actual-flow acceleration condition. It is not an independently
proved sign for the unknown actual acceleration.

If \(u_N-\ell_N>0\), its normalized position is
\[
\boxed{
\theta_N=
\frac{\ddot c_{2N-1}-\mathfrak d_N-2N(2N+1)\ell_N}
 {2N(2N+1)(u_N-\ell_N)}.
}
\tag{28}
\]
The next signs hold exactly when \(0<\theta_N<1\); their innovations
are \((u_N-\ell_N)\theta_N\) and
\((u_N-\ell_N)(1-\theta_N)\).
This is the same canonical position already present in LCW,
now read from acceleration. Its probability-like split follows
after the positive budget and the interior position are proved.
It is not assumed to be an invariant or a probability in advance.

For data whose odd window has positive width, choosing z inside
or outside that window leaves all the old forms and their first
derivatives unchanged while changing the next signs. This is an
algebraic demonstration of insufficient first-jet information.
It does not alter the actual xi coefficients or assert an actual
counterexample to RH.

A complete analytic flow record includes second and higher
derivatives. The insufficiency concerns retaining only values
and velocities, not retaining the full analytic record.

## Consequence for the proposed closure route

There are now two precise dynamic targets:

1. Control the actual boundary readings \(\chi_A\) and \(\chi_B\) so that
   no compressed-form singularity occurs when carrying a known
   positive bank through the desired anchor range.
2. Control the actual acceleration in (27) at each new depth.

The conserved pairing (14) handles regular transport within a
fixed dimension. It does not manufacture the additional independent
direction tested by (27). Likewise the scalar fluctuation equation
retains its third central reading rather than closing on two
numbers by assumption.

Proving (27) for the actual source at every depth would establish
the existing all-order positivity criterion and hence RH.
This note identifies the exact dynamic inequality; it does not
prove that universal statement. APD's established finite ranges
and the existing source-identification requirements remain in force.

## Sources, attribution, and verification

- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md),
  T1--T2: actual anchor hierarchy and finite-model transport.
- [LCW-1](papers/native-layer-continuation-window/THEOREM.md),
  T1--T3: sharp two-layer windows and their canonical positions.
- [LOT-1](papers/native-layer-observer-transport/THEOREM.md):
  acquired channel pairings, carried observables and the
  distinct new-depth contribution.
- [RBP-1](RH_RELATIONAL_BANK_POSITIVITY.md), T2--T3:
  proved finite source forms and their interior probability carrier.
- [ORR-1](RH_OBSERVER_RESIDUAL_RECURSION.md):
  the corresponding mixed-observer and residual update problem.
- [APD-1](RH_OBSERVER_APERTURE_DEPTH.md):
  current finite observation-bank ranges.

Linear transport congruences, logarithmic determinant identities,
moment fluctuations, Schur complements and canonical moments are
standard tools. The additions here are their explicit actual-source
generator, boundary-pole condition, and acceleration reading of the
retained extension gate. No literature-wide novelty is claimed.

Verification consists of written algebra, differentiation and analytic
endpoint reasoning, plus text/readback and repository metadata checks.
No scientific computation, numerical source evaluation, independent
agent, certificate replay or workflow dispatch is performed.
Only this note and its README entry change. The commit uses
[skip ci]. Published as an unmerged draft stacked on #135.
