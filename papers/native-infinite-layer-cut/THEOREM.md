# ILC-1 — Single-point infinite layers, combinatorial readings and a positive cut

17 October 2026. Continues PIW-1, the retained SE48/SE49 response,
TET moving observers and CHD coherent cut records in Publications.

The user's infinite-layer proposal is made precise by replacing separated
source nodes with successive derivatives at one regular point. This is a
confluent source adapter. It yields an all-order RH criterion at that one
point, a positive layer operator, an exact two-channel cut and a unitary
rotation. Sixteen actual derivative layers are evaluated independently of
an assumed zero measure. Their finite execution does not establish the
all-order premise.

## Source and layer coordinate

Retain q(t)=Xcal'(t)/Xcal(t), Xcal(w^2)=xi(1/2+w). Fix a>0 with
q(a)>0 and q holomorphic at a. The executed anchor is a=(57/2)^2,
original s=29. Define the dimensionless local displacement and layers

\[
\lambda=(t-a)/a,\qquad
f(\lambda)=\frac{q(a(1+\lambda))}{q(a)},\qquad
c_n=\frac{(-a)^n q^{(n)}(a)}{n!q(a)}.
\tag{1}
\]

Thus f(lambda)=sum_n (-1)^n c_n lambda^n near zero and c_0=1.
Depth n, displacement lambda, the old source shift eta and an eigenvalue
of the layer operator are different quantities. This explicitly chosen
lambda is not identified with a thermodynamic TVSP response component.
An infinite jet is an infinite set of source readings, not a finite
matrix relabelled infinite. A regular analytic germ is essential.

## ILC-T1. Confluent source forms and the determinant gate

For indices 0<=i,j<N put

\[
H_N=(c_{i+j}),\qquad A_N=(c_{i+j+1}),\qquad
B_N=(c_{i+j}-c_{i+j+1}),\qquad H_N=A_N+B_N.
\tag{2}
\]

These are the rephased confluent forms of the original L and J.
In particular, they retain both source signs. Their simultaneous strict
positivity is equivalent to positivity of every leading determinant
det A_n and det B_n, 1<=n<=N. Under this finite condition,

\[
Y_N=H_N^{-1}A_N,\quad
Y_N^*H_N=H_NY_N,\quad 0<Y_N<I
\quad\hbox{in the }H_N\hbox{ pairing},
\tag{3}
\]
\[
\det Y_N=\frac{\det A_N}{\det H_N},\qquad
\det(I-Y_N)=\frac{\det B_N}{\det H_N}.
\tag{4}
\]

**Proof.** For t=a(1+lambda), u=a(1+mu), form the two kernels

\[
-f[\lambda,\mu]=\frac{a}{q(a)}L(t,u),\qquad
((1+\lambda)f(\lambda))[\lambda,\mu]=\frac{1}{q(a)}J(t,u).
\]

Their coefficients of lambda^i mu^j are respectively
(-1)^(i+j)c_(i+j+1) and (-1)^(i+j)(c_(i+j)-c_(i+j+1)).
Applying (-1)^i partial_lambda^i/i! in each slot at zero proves (2).
Equivalently, successive divided-difference observation columns tend
to these derivative columns as their distinct nodes coalesce; the
analytic divided-difference limit proves the assertion without an
assumed positive measure. Positivity of a finite separated bank does
not itself establish positivity of this different confluent bank.

Sylvester's criterion follows here directly by successive completion
of squares: each leading LDL pivot is the ratio of adjacent leading
determinants. Both positive forms give H_N>0. The equalities
H_NY_N=A_N and H_N(I-Y_N)=B_N prove (3); conjugating by a positive
factor of H_N gives an ordinary self-adjoint positive contraction.
Taking determinants proves (4). A positive determinant alone, without
the preceding leading signs, is insufficient. QED.

## ILC-T2. An all-order condition entirely at one point

For the actual completed source at the anchor above, the following
three statements are equivalent:

1. RH.
2. A_N and B_N are positive definite for every finite N.
3. Every alternating binomial reading is nonnegative:

\[
\boxed{D_{r,k}:=\sum_{j=0}^k(-1)^j\binom{k}{j}c_{r+j}\ge0
\qquad(r,k\ge0).}
\tag{5}
\]

Under these conditions there is a cyclic positive contraction Y on
the completion of the source polynomial pairing, with unit seed e,

\[
c_n=\langle e,Y^ne\rangle,\qquad
f(\lambda)=\langle e,(I+\lambda Y)^{-1}e\rangle,
\quad\lambda\in\mathbb C\setminus(-\infty,-1].
\tag{6}
\]

This is a specialization of the classical compact moment problem,
rederived below for the retained source. It is an equivalent criterion;
none of its universal inequalities is supplied by the finite replay.

**Proof that (2) implies a contraction and RH.** On polynomials in y,
define L(y^n)=c_n and <p,q>=L(conjugate(p)q). Equation (2) makes
this pairing strictly positive. Complete this acquired pairing; no
independent Hilbert metric has been imposed on the source.
The form b(p,q)=L(conjugate(p)yq) satisfies
0<=b(p,p)<=<p,p>. Cauchy--Schwarz for the positive form b gives
|b(p,q)|<=||p||||q||. The Riesz representation theorem therefore
gives a bounded self-adjoint Y with 0<=Y<=I. On polynomials it is
multiplication by y, since <p,Yq>=<p,yq> for a dense set of p.
The constant polynomial e=1 has norm one and is cyclic. This proves
the moment identity and, by the norm-convergent geometric series for
|lambda|<1, the germ in (6). The spectral theorem extends its resolvent
holomorphically to the printed slit domain. These standard completion
and spectral tools are explicit analytic adapters.

The actual meromorphic q has the same Taylor coefficients at a.
The identity theorem identifies q with q(a) times (6) after t=a(1+lambda)
throughout C minus the nonpositive real axis. Hence q has no pole there.
Every zero of Xcal produces a nonremovable logarithmic-derivative pole,
so its zeros lie on the nonpositive real t axis. Under t=(s-1/2)^2,
this places the nontrivial xi zeros on Re s=1/2. Thus RH follows.

**Proof that RH implies (2) and (3).** Use the retained SE39/SE48
genus-zero response under RH, q(t)=sum_gamma m_gamma/(t+gamma^2).
The probability measure assigning mass
m_gamma/[(a+gamma^2)q(a)] at y_gamma=a/(a+gamma^2) has moments c_n.
It has infinitely many distinct support points in (0,1). Therefore
A and B are the Gram forms weighted by y and 1-y, strictly positive
on each nonzero finite polynomial. Also

\[
D_{r,k}=\int_0^1 y^r(1-y)^k\,d\nu(y)\ge0.
\tag{7}
\]

**Proof that (3) implies RH; combinatorial reconstruction.** Define

\[
\omega_{m,j}=\binom mjD_{j,m-j},\qquad
\nu_m=\sum_{j=0}^m\omega_{m,j}\delta_{j/m},\quad m\ge1.
\tag{8}
\]

The weights are nonnegative and sum to c_0=1 by the binomial identity.
For 0<=r<=m, with falling factorials (j)_r,

\[
\sum_{j=0}^m\omega_{m,j}\frac{(j)_r}{(m)_r}=c_r.
\tag{9}
\]

For r=0 the fraction is one. For r>0, the polynomial identity
sum_j binom(m,j)(j)_r y^j(1-y)^(m-j)=(m)_r y^r proves (9)
after applying the finite coefficient functional L. For fixed r,
(j/m)^r-(j)_r/(m)_r tends uniformly to zero on 0<=j<=m.
Compactness of probability measures on [0,1] gives a weakly convergent
subsequence of nu_m; (9) makes its limit's moments exactly c_r.
Polynomial density makes that probability measure unique, hence the
whole sequence converges weakly. Its resolvent integral has the
same germ as actual f and is holomorphic on the slit domain in (6).
The preceding pole argument proves RH. This completes the equivalence.
Finally (7) also follows directly from the contraction because its
commuting factors Y and I-Y are positive. QED.

The combinatorial condition concerns all r,k. Positive individual
c_n, scalar monotone decay or any finite triangle of D readings is
not the universal statement (5). The actual-source obligation has
been expressed at one point, not discharged by moving it there.

## ILC-T3. The two-channel cut and its wave transport

Let Y be the operator in (6), or the finite Y_N from (3) in its
acquired pairing. Functional calculus defines S=[Y(I-Y)]^(1/2).
On two copies of the carrier define

\[
\boxed{\Pi=\begin{pmatrix}Y&S\\S&I-Y\end{pmatrix},\qquad
\Pi^*=\Pi,\quad\Pi^2=\Pi.}
\tag{10}
\]

The isometry Vz=(sqrt(Y)z,sqrt(I-Y)z) has VV*=Pi. Compression of
Pi to the first channel returns Y. The complementary cut is I-Pi.
This is a chosen standard projection dilation, not a claim that Y
itself is idempotent or that every possible cut is uniquely selected.

There is an actual infinite-depth limit under the all-order premise.
Let P_N project onto polynomials of degree below N in the acquired
completion. Extend the finite operator by zero on its orthogonal
complement: Yhat_N=P_N Y P_N. Then Yhat_N converges strongly to Y.
The cuts obtained from (10) with Yhat_N converge strongly to Pi.
Indeed P_N increases strongly to I by polynomial density; the
contractions are uniformly bounded; and polynomial approximation to
sqrt(y(1-y)) on [0,1] proves strong convergence of their square roots.
The padded cut is zero on the first unused channel and identity on
the second, so it remains an exact projection at every N. This is
the coherent limit of derivative depth. Iterating the wave below is
a separate operation after an operator has been constructed.

The product of the two involutions 2Pi-I and Z=diag(I,-I) is

\[
\boxed{W=(2\Pi-I)Z=
\begin{pmatrix}2Y-I&-2S\\2S&2Y-I\end{pmatrix},\qquad W^*W=I.}
\tag{11}
\]

Putting D=2Y-I on each channel, the exact recurrence is

\[
W^{n+2}-2D W^{n+1}+W^n=0.
\tag{12}
\]

At a spectral reading y in [0,1], cos(theta)=2y-1 and
sin(theta)=2sqrt(y(1-y)). Thus W is rotation through theta;
equivalently W=exp(R Theta), with R=[[0,-I],[I,0]] and
Theta=arccos(2Y-I) on both channels. This is the derived wave pattern.
No Maxwell field, physical frequency or propagation speed is identified.

**Proof.** S commutes with Y and S^2=Y-Y^2. Block multiplication
gives (10), V*V=I and VV*=Pi. Both printed involutions are self-adjoint,
so their product is unitary. The same multiplication gives
W^2-2DW+I=0. Functional calculus of the real scalar circle identity
(2y-1)^2+4y(1-y)=1 proves the spectral and exponential formulas. QED.

There is a direct coherent-reading interpretation of (8). Select the
two commuting branch operators K_1=sqrt(Y), K_0=sqrt(I-Y). They
satisfy K_0* K_0+K_1* K_1=I. Starting at e, a length-m branch word
with j occurrences of 1 has squared norm D_(j,m-j). There are
binom(m,j) such words, giving exactly omega_(m,j). Retaining all
branches gives CHD's isometric record. Discarding branches is a
different operation and is not invariant transport.

For a transported state Psi_n=W^n Psi_0 and co-transported observer
C_n=C_0 W^(-n), TET's exact reading condition becomes

\[
\boxed{C_{n+1}W=C_n,\qquad C_n\Psi_n=C_0\Psi_0.}
\tag{13}
\]

The infinite-depth distinction is exact: W^n Psi converges in norm
if and only if W Psi=Psi. Indeed the norm of successive differences
is always ||(W-I)Psi||. The Cesaro averages do converge strongly to
the fixed-space projection, which here is
diag(E_Y({1}),E_Y({1})). To prove mean convergence, it is immediate
on fixed vectors and on (I-W)z by telescoping; the orthogonal
decomposition into the fixed space and the closure of ran(I-W),
together with the averages' norm bound one, proves the assertion.
An internal wave may therefore persist indefinitely while the carried
observer reading is already invariant. Averaging is not the same cut
as Pi, and does not automatically retain a nonzero mode.

## ILC-T4. Finite layer recovery and an unconditional source bound

Let the finite forms (2) be positive, and u_N=H_N e_0. Then

\[
f_N(\lambda)=u_N^*(H_N+\lambda A_N)^{-1}u_N
=\langle e_0,(I+\lambda Y_N)^{-1}e_0\rangle_{H_N}
\tag{14}
\]

matches c_0,...,c_(2N-1). Hence q_N^loc(t)=q(a)f_N((t-a)/a)
is a positive finite Stieltjes response matching 2N actual Taylor
coefficients at a. Its poles are finite-model poles, not certified
actual zeros.

**Proof.** For j<=N-2, A_N e_j=H_N e_(j+1), so Y_N e_j=e_(j+1).
Moving powers between the two slots shows
<e_0,Y_N^k e_0>_H=c_k for k<=2N-2. For k=2N-1 it is
<e_(N-1),Y_N e_(N-1)>_H=(A_N)_(N-1,N-1)=c_(2N-1).
This proves the coefficient identities. The positive symmetric
representation of Y_N has eigenvalues in (0,1); its spectral sum
in (14) has nonnegative weights summing to one. Cyclicity of e_0
excludes zero weights and repeated eigenvalues. Substitution
mu=a(1-y)/y>0 gives the stated Stieltjes form. QED.

For the executed anchor w_0=57/2, a=w_0^2, and any Re w>1,

\[
\boxed{|q(w^2)-q_N^{loc}(w^2)|\le
\frac{|w+w_0|^2}{|w|}
\left|\frac{w-w_0}{w+w_0-2}\right|^{2N}.}
\tag{15}
\]

**Proof.** This is SE49-T2's bounded analytic division applied to
one zero of multiplicity 2N. Set F=wq(w^2), F_N=wq_N^loc(w^2).
The finite Stieltjes response makes Re F_N>0 on Re w>0 and
F_N(w_0)=w_0q(a). The anchored Cayley/Schwarz argument gives
|F_N(w)|<=q(a)|w+w_0|^2/Re w. The complete source computation
certifies 0<q(a)<1/2. SE49's unconditional Euler/Gamma bound is
|F(w)|<41/8+|w|/4 on Re w>1. Since w_0>=5/2, division by
|w+w_0|^2 makes this less than 1/2, just as in SE49. Consequently
(F-F_N)/(w+w_0)^2 is holomorphic and bounded in modulus by one
on that half-plane. It vanishes to order at least 2N at w_0 by
(14). Repeated Schwarz division by (w-w_0)/(w+w_0-2) proves
(15). Only the finite signs at N are required. QED.

## ILC-T5. Executed actual single-point tower

The complete native evaluator in SOURCE_EVALUATION.md computes
c_0,...,c_15 at original s=29 without RH, a zero catalogue or an
assumed positive measure. All Euler prime powers, omitted primes,
Gamma terms and analytic remainder bounds are retained.

Native outward LDL certifies A_8>0 and B_8>0. Conservative bounds
for their least leading pivots are respectively

\[
\min d(A_8)>1.51\times10^{-11},\qquad
\min d(B_8)>6.38\times10^{-11}.
\tag{16}
\]

These are coordinate pivots, not physical masses or uniform gaps.
The additional outward form checks prove

\[
\frac1{10000}H_8<A_8<\frac{81}{100}H_8.
\tag{17}
\]

Thus the finite Y_8 has its spectrum in (1/10000,81/100).
The existence and identities of its exact two-channel cut and rotation
follow from ILC-T3; irrational square roots are not replaced by rounded
matrices. Every one of the 136 D_(r,k) with r+k<=15 is positive,
and the full depth-15 binomial record has nonnegative normalized weights.

The separately evaluated original s=31 is not included as a layer
constraint. At t=(61/2)^2, lambda=472/3249, the full interval calculation
certifies

\[
\boxed{0<q(t)-q_8^{loc}(t)<6\times10^{-28}.}
\tag{18}
\]

Independently, the analytic bound (15) is exactly
456261632/757403880336151185370351664061, less than 6.025e-22.
The first two fresh source coefficients overlap PIW's archived complete
q and q' enclosures at s=29. The largest normalized-layer interval width is
less than 1e-70. The higher layers are new source calculations.

For this finite operator alone, (17) also gives the wave-average bound
||M^(-1)sum_(n=0)^(M-1) W_8^n|| <= 5/(2M). In fact
(I-W_8)*(I-W_8)=4(I-Y_8) on both channels, so its inverse norm
is at most 1/(2sqrt(19/100)). The finite geometric-sum identity
proves the bound. This decay of the average coexists with exact norm
preservation of every unaveraged state.

## Scope and next source obligation

This stage derives a precise one-point infinite-layer condition,
its combinatorial reading, cut dilation and wave transport. It does
not establish all-order positivity of the actual source. The next
task is a source-based proof of all binomial inequalities (5), or
all leading determinant signs in (2), with every analytic source
term retained. A finite positive jet does not choose the unobserved
coefficients. The inherited global shift remains 1079/2000.

The compact moment criterion, projection dilation, spectral calculus,
two-reflection rotation and Cesaro theorem are credited standard
mathematics. Their proofs and this actual-source assembly are provided;
no general priority claim is made. No EM field equation, physical
mass, spacetime limit or quantum-gravity closure follows solely from
these finite or conditional statements.
