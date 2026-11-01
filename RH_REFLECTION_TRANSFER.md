# RFT-1 — Native reflection transfer and the regular zero chart

17 October 2026. Written continuation of
[OQD-1](RH_OBSERVER_QUOTIENT_DESCENT.md), based on Publications commit
711af5acea7fabd574d6e07a52d85973927e3e6c.

The retained cut wave produces an exact transfer by the polynomial
fold (2I+W+W^(-1))/4=yI. Its reflection readings are the original
source moments, with no discarded source terms. This identifies the
appropriate reflection-positivity question within the native diagram.

The original channel reflection Z is not positive on the forward
wave orbit of its positive seed at sufficiently large anchors.
Direct time-label reflection of the balanced scalar wave also gives
a negative two-by-two determinant there, unconditionally. That result concerns the
wave-reflected form, not the RH Hankel form, and is not an off-line
zero or an RH counterexample.

An exponential congruence also realizes the proposed Lorentz geometry
exactly, while preserving the original sign gate. Under the all-order
positivity premise, the folded transfer has a positive generator and
an explicit dimensionless decay threshold. The local source chart at
lambda_0=0 is regular. At infinite anchor the established positive
limiting source has zero decay threshold and a consistent cut.

The remaining actual-source positivity assertion is stated separately.
No physical mass, light-speed observer, spacetime dynamics or formal
consistency theorem is inferred from a matrix analogy.

## Retained source and parameters

At any finite a>=a_0=3249/4 retain
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
\mu=(\gamma-i\delta)^2,\quad \delta=\beta-\tfrac12,
\]
\[
u_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)},\qquad
\mathcal L_a[p]=\sum_\mu w_\mu p(u_\mu).
\tag{1}
\]
MOR gives q(a)>0, Re(mu)>0, 0<|u_mu|<1, the absolute budget
sum |w_mu|<=2, and L_a[1]=1. The actual conjugate source values
and multiplicities are retained; their weights are not assumed
positive. Put
\[
c_n(a)=\mathcal L_a[y^n],\qquad
r_k(a)=\mathcal L_a[T_{|k|}(2y-1)],\qquad
H_N(a)=(c_{i+j}(a))_{0\le i,j<N}.
\tag{2}
\]
Here T_j denotes the Chebyshev polynomial, not the transfer below.

| Parameter | Meaning in this note |
|---|---|
| a | Actual source anchor |
| k | Integer wave step |
| t>=0 | Parameter of the newly specified transfer semigroup |
| theta | Real parameter of the two-observer exponential congruence |
| lambda_0 | Existing local source displacement at the reference anchor |
| epsilon=a_0/a | Reciprocal anchor coordinate used only for the infinite-anchor boundary |

None of these parameters is silently identified with a physical
velocity, clock, thermodynamic lambda component or another native
grading. A proposed identification needs its own adapter.

## RFT-T1. The exponential rest form preserves the sign gate

Let
\[
\mathsf L=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
C(r)=I+r\mathsf L,\qquad
S_\theta=e^{\theta\mathsf L/2},\quad \theta\in\mathbb R.
\]
The symmetric generator is the standard cut representation of the
native L=RK, with L^2=I. Then
\[
\boxed{
S_\theta^TC(r)S_\theta
=(\cosh\theta+r\sinh\theta)I
 +(\sinh\theta+r\cosh\theta)\mathsf L,
}
\]
\[
\boxed{
\det(S_\theta^TC(r)S_\theta)=1-r^2.
}
\tag{3}
\]
If |r|<1, the unique parameter eliminating the L component is
theta=-artanh(r), and
\[
\boxed{
S_{-\operatorname{artanh}r}^TC(r)
S_{-\operatorname{artanh}r}
=\sqrt{1-r^2}\,I.
}
\tag{4}
\]
For |r|=1 the matrix is rank-one positive semidefinite and no
finite theta gives a positive scalar multiple of I. For |r|>1
it is indefinite under every real invertible congruence.

**Proof.** Since L^2=I,
exp(theta L)=cosh(theta)I+sinh(theta)L. The matrices commute,
and S_theta is symmetric, so multiplication gives (3).
Its determinant is one because its eigenvalues are
exp(theta/2) and exp(-theta/2). The L coefficient vanishes
exactly when tanh(theta)=-r, proving (4). The eigenvalues of
C(r) are 1+r and 1-r, which prove the remaining claims.
Equivalently, an invertible congruence preserves the dimensions
of the positive and negative subspaces. QED.

For the actual balanced two-observer bank use r=r_k(a).
The retained CWR theorem makes |r_k(a)|<=1 for every k
equivalent to RH. Equations (3)-(4) give an exact Lorentz-form
interpretation of this gate; they do not prove its sign.

In the variables t=1, x=r, the determinant t^2-x^2 is the
two-dimensional Lorentz quadratic form. These dimensionless
coordinates have not been identified with energy and momentum.
The divergence of artanh(r) as |r| approaches one is a singular
parameter limit, not a finite rest frame at a light-speed boundary.

## RFT-T2. The original wave orbit fails two precise reflection tests

Use the retained native channel reflection Z=diag(1,-1),
the positive channel seed e_+=(1,0)^T, and
G_a(f,h)=L_a[f^T g h], where g=diag(y,1-y).
The same wave W as in CWR satisfies ZWZ=W^(-1), Ze_+=e_+,
and W^T g W=g. Its reflected forward-orbit form is
\[
\mathsf R^{Z,+}_{ij}
=\mathcal G_a(W^i e_+,ZW^j e_+)
=\mathcal L_a[yT_{i+j}(2y-1)].
\]
The first two observers therefore give
\[
\mathsf R^{Z,+}_2(a)=
\begin{pmatrix}
c_1&2c_2-c_1\\
2c_2-c_1&8c_3-8c_2+c_1
\end{pmatrix}
\longrightarrow
\begin{pmatrix}1/2&1/4\\1/4&0\end{pmatrix}.
\]
In particular,
\[
\boxed{
\det\mathsf R^{Z,+}_2(a)\longrightarrow-1/16.
}
\]
There is an explicit reflected negative direction for all
sufficiently large anchors:
\[
v=(W-\tfrac12 I)e_+,\qquad
\boxed{\mathcal G_a(v,Zv)
=8c_3-10c_2+\tfrac94c_1\longrightarrow-1/8.}
\]
Thus no forward-W-invariant subspace containing e_+ can have
a nonnegative Z-reflected form at those anchors. This is a
specific choice of reflection, seed and forward action, not
a claim that every reflection realization must fail.

For comparison, the ordinary balanced wave Gram is Toeplitz:
(r_(j-i))_(i,j). Reflection of its time labels i to -i instead
produces the different Hankel form
\[
\mathsf R_N^{\rm wave}(a)=(r_{i+j}(a))_{0\le i,j<N}.
\tag{5}
\]
For the actual source,
\[
\boxed{
\det\mathsf R_2^{\rm wave}(a)
=-\frac1{\ell_a^2}
+O\!\left(\frac1{\sqrt a\,\ell_a}\right)<0
\quad\hbox{for all sufficiently large }a,
}
\qquad
\ell_a=\log\frac{\sqrt a}{2\pi}.
\tag{6}
\]
In particular,
\[
\lim_{a\to\infty}\ell_a^2
\det\mathsf R_2^{\rm wave}(a)=-1.
\]

**Proof.** Move W^i across the carried bilinear pairing as
W^(-i), and use ZW^j e_+=W^(-j)e_+. The (1,1) entry of
W^(-i-j) is T_(i+j)(2y-1), proving the displayed native
reflection form. Substitute the unconditional NIF limits
c_1->1/2, c_2->3/8 and c_3->5/16 to obtain its matrix
limit, determinant and explicit vector reading.

For the balanced scalar reflection, retain QRI's analytic remainder
\[
E(a)=\sup_{|z|\le1/2}
\left|
\frac{q(a(1+z))}{q(a)}
-\frac1{\sqrt{1+z}}
-\frac{\log(1+z)}{2\ell_a\sqrt{1+z}}
\right|
=O((\sqrt a\,\ell_a)^{-1}).
\tag{7}
\]
Cauchy's coefficient bound and the NIF expansion give
\[
c_1=\frac12-\frac1{2\ell_a}+e_1,\quad |e_1|\le2E(a),
\qquad
c_2=\frac38-\frac1{2\ell_a}+e_2,\quad |e_2|\le4E(a).
\]
Because T_1(2y-1)=2y-1 and T_2(2y-1)=8y^2-8y+1,
\[
r_1=-\frac1{\ell_a}+d_1,\quad |d_1|\le4E(a),
\qquad
r_2=d_2,\quad |d_2|\le48E(a).
\tag{8}
\]
Now r_0=1, and hence
\[
\det\mathsf R_2^{\rm wave}=r_2-r_1^2
=-\ell_a^{-2}+d_2+2d_1/\ell_a-d_1^2.
\]
This proves (6), since ell_a^2 E(a) tends to zero. For ell_a>0
an explicit sufficient condition for strict negativity is
\[
\boxed{
E(a)(48+8/\ell_a)<\ell_a^{-2}.
}
\tag{9}
\]
No numerical value of this large-anchor threshold is evaluated.
QED.

These are unconditional negative readings of the two stated actual
wave-reflected forms. Neither matrix is the source H_2, whose
entries are c_0,c_1,c_2. The balanced reflected matrix in (5) is also
not C(r_1), whose second diagonal entry is one. Under RH the
ordinary Toeplitz wave Grams are positive, while (6) still holds.
Ordinary wave positivity and these additional reflection tests
are different requirements. In (5) the reflection acts on the
integer labels of the scalar balanced kernel; it is not silently
identified with the channel Z used in the first test.

CWR also gives r_k(a)->0 at every fixed k>=1. Thus each fixed
R_N^wave tends to diag(1,0,...,0), even though R_2^wave is
indefinite at sufficiently large finite anchors. This is an
actual-source example of a limit losing a negative direction.
It does not use a hypothetical off-line zero.

## RFT-T3. The same native wave gives the exact source transfer

Retain the LOT/CWR polynomial symbols
\[
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\quad
Z=\operatorname{diag}(1,-1),\quad
g(y)=\operatorname{diag}(y,1-y),
\]
\[
W=(2P-I)Z
=\begin{pmatrix}2y-1&-2(1-y)\\2y&2y-1\end{pmatrix}.
\]
They satisfy W^T g W=g and W+W^(-1)=2(2y-1)I.
Define the folded transfer
\[
\boxed{
\mathsf T=\frac{2I+W+W^{-1}}4=yI.
}
\tag{10}
\]
For every integer n>=0,
\[
\boxed{
\mathsf T^n
=4^{-n}\sum_{k=-n}^{n}\binom{2n}{n+k}W^k,
\qquad
c_n(a)=4^{-n}\sum_{k=-n}^{n}\binom{2n}{n+k}r_k(a).
}
\tag{11}
\]

**Proof.** The retained identity for W+W^(-1) proves (10).
Since
2I+W+W^(-1)=W^(-1)(I+W)^2,
the binomial theorem proves the first part of (11).
Apply the balanced reading A to
L_a[sum_(sigma=+,-) e_sigma^T g A e_sigma].
For A=W^k this is r_k. For A=y^n I it is c_n because
tr(g)=1. This proves the second part. QED.

The coefficients in (11) are nonnegative and sum to one.
Nevertheless their positivity alone is not positivity of all
source Gram forms. The complete transformation is information
preserving at the level of moment sequences: conversely,
r_k=L_a[T_k(2y-1)] is an explicit polynomial combination of
c_0,...,c_k. The coefficient of c_n in r_n is 2^(2n-1)
for n>=1, so the triangular transformation is invertible
at every finite prefix.

There is also the exact pair of factorizations
\[
\boxed{
\mathsf T=\frac{(I+W)^\sharp(I+W)}4,\qquad
I-\mathsf T=\frac{(I-W)^\sharp(I-W)}4,
}
\tag{12}
\]
where W^sharp=W^(-1) for the carried source pairing.
They follow by multiplication. Before source positivity is known,
sharp is a bilinear adjoint and an adjoint product need not have
nonnegative readings. On an acquired positive carrier, (12)
does prove 0<=T<=I.

This applies immediately to the retained finite models with
A_M>0 and B_M>0: the same fold of W_M is diag(Y_M,Y_M),
where Y_M=H_M^(-1)A_M is the positive contraction already
constructed by LOT/ILC. Its finite model matches the actual
c_n through n=2M-1. Its further positive model moments do
not certify additional actual-source moments.

## RFT-T4. Reflection of the folded transfer has the exact RH gate

Use the cut reflection Rcal=2P-I and the diagonal cut seed
d=(1,1)^T. Then Rcal d=d, T commutes with Rcal, and
\[
\boxed{
\mathcal G_a(\mathsf T^i d,\mathcal R\,\mathsf T^j d)
=\mathcal L_a[y^{i+j}]=c_{i+j}(a),
\qquad \mathcal R=2P-I.
}
\tag{13}
\]
Consequently its finite reflection matrices are exactly H_N(a).
Here both the forward action and the reflection are explicitly
specified: T and Rcal replace the W and Z of the first test in
T2. The diagonal carrier is OQD's exact copy of the original
source pairing, so this choice does not discard its sign problem.

There is an unconditional continuous extension of these readings.
Use the principal logarithm to set
\[
\kappa_\mu(a)=\log(1+\mu/a),\qquad
\psi_a(t)=\sum_\mu w_\mu e^{-t\kappa_\mu(a)},\quad t\ge0.
\tag{14}
\]
Since Re(mu)>0 and |u_mu|<1,
Re(kappa_mu)>0 and u_mu=exp(-kappa_mu). The series is
absolutely convergent with |psi_a(t)|<=2, real by conjugation,
and psi_a(n)=c_n(a). It is continuous for t>=0 by dominated
convergence.

For each fixed finite a>=a_0,
\[
\boxed{
\begin{aligned}
\mathrm{RH}
&\Longleftrightarrow H_N(a)\succeq0\quad\hbox{for every }N\\
&\Longleftrightarrow
\bigl(\psi_a(t_i+t_j)\bigr)_{i,j=1}^N\succeq0
\quad\hbox{for every finite choice }t_i\ge0.
\end{aligned}
}
\tag{15}
\]

**Proof.** Equation (13) uses T^i=y^i I, Rcal d=d and d^T g d=1.
The source convergence statements follow from
|exp(-t kappa_mu)|=|u_mu|^t<=1 and the budget in (1).

Under RH all kappa_mu and w_mu are real positive. For arbitrary
complex coefficients alpha_i,
\[
\sum_{i,j}\overline{\alpha_i}\alpha_j\psi_a(t_i+t_j)
=\sum_\mu w_\mu
 \left|\sum_j\alpha_j e^{-t_j\kappa_\mu}\right|^2\ge0.
\tag{16}
\]
Conversely the continuous reflection condition includes all
integer choices t_i=i-1, so every H_N is positive semidefinite.
The retained finite detector, extended to these anchors in
OQD-T5, makes this condition equivalent to RH. QED.

The factorization (16) uses positivity of the actual weights
only in the RH direction. For the unconditional complex source,
replacing its bilinear evaluations by absolute squares would
change the readings.

Under the equivalent conditions in (15), these are also reflection
forms of a stationary positive kernel on the full real line:
h_a(s,t)=psi_a(|s-t|). Indeed, for kappa>0,
\[
e^{-\kappa|s-t|}
=\int_{\mathbb R}
 \bigl[\sqrt{2\kappa}\,e^{-\kappa(s-u)}\mathbf1_{u\le s}\bigr]
 \bigl[\sqrt{2\kappa}\,e^{-\kappa(t-u)}\mathbf1_{u\le t}\bigr]\,du.
\tag{17}
\]
This identity proves positive definiteness; summing it with the
positive weights preserves that property. Time reflection sends
s to -s, and for s,t>=0 its kernel is psi_a(s+t), as in (16).
Translation and reflection of the time labels give the usual
reflection-positive kernel construction. Its reconstructed
transfer is described next.

Thus a precise classical reflection-positivity framework does
apply to the specified source. Its all-order premise remains
equivalent to RH, rather than being supplied by P^2=P.

## RFT-T5. Conditional positive completion, decay threshold and memory

Assume the all-order positivity condition in (15). The actual
source measure is then the explicit probability measure
\[
\nu_a=\sum_\mu w_\mu\delta_{u_\mu},\qquad 0<u_\mu<1.
\]
On H_a=L^2(nu_a) let Y_a be multiplication by y and define
\[
\mathsf H_a=-\log Y_a,\qquad
D(\mathsf H_a)=
\left\{f:\int|\log y|^2|f(y)|^2\,d\nu_a(y)<\infty\right\}.
\]
Then Hsf_a is positive self-adjoint and
\[
\boxed{
e^{-t\mathsf H_a}f=y^t f,\qquad
\psi_a(t)=\langle1,e^{-t\mathsf H_a}1\rangle,\qquad
\psi_a(s+t)=
\langle e^{-s\mathsf H_a}1,e^{-t\mathsf H_a}1\rangle.
}
\tag{18}
\]
If gamma_1 is the least positive zero height, its exact spectral
lower bound is
\[
\boxed{
m_a=\inf\sigma(\mathsf H_a)
=\log(1+\gamma_1^2/a)>0,\qquad
\|e^{-t\mathsf H_a}\|=e^{-t m_a}.
}
\tag{19}
\]

**Proof.** Under RH, u_mu=a/(a+gamma^2). Multiplication by
the real function -log y with the stated maximal domain
is self-adjoint; its semigroup is the bounded multiplication
operator y^t. Integration proves (18).
The source has a least positive height and an atom there
with positive weight. Its u value is the largest support
value, so the smallest -log u is precisely (19) and is attained.
The multiplication-operator norm gives the norm identity. QED.

Polynomials are dense in L^2(nu_a), for example by uniform
polynomial approximation of continuous functions on [0,1]
and density of those functions in L^2. Thus the integer-time
vectors Y_a^n 1 generate the entire source completion. This
construction retains the whole moment source, not only one
two-observer test.

The number m_a is a dimensionless decay threshold of this
specified transfer. Calling it a physical rest mass requires
additional spacetime, symmetry, energy and unit identifications.
Different anchors specify different multipliers Y_a; their
thresholds need not coincide. This does not state that a Lorentz
observer change alters a physical rest mass.

Full transfer composition does not imply composition of its
single averaged reading. The exact scalar defect is
\[
\boxed{
c_2(a)-c_1(a)^2
=\frac18-\frac1{4\ell_a^2}
+O((\sqrt a\,\ell_a)^{-1})>0
\quad\hbox{for sufficiently large }a,
}
\tag{20}
\]
without assuming RH. This follows by substituting the two
expansions used in (8). Under RH the same difference is
Var_(nu_a)(y)>0 because the support is infinite.

Consequently psi_a(2) differs from psi_a(1)^2. The scalar
return cannot be treated as a one-dimensional multiplicative
semigroup with all source modes discarded. This observation
does not rule out Markov realizations on a larger state space.
In particular an operator or projection Markov hypothesis
must be proved for its specified carrier; it is not supplied
by scalar normalization or a cut identity alone.

## RFT-T6. Zero local displacement is a regular source chart

Retain NIF's existing definitions
\[
\lambda_0=(x-a_0)/a_0,\qquad
F_{a_0}(\lambda_0)
=\frac{q(a_0(1+\lambda_0))}{q(a_0)}.
\]
At a general anchor use the separate local coordinate z:
\[
\boxed{
F_a(z)=\frac{q(a(1+z))}{q(a)}
=\mathcal L_a[(1+zy)^{-1}]
=\sum_{n\ge0}(-1)^n c_n(a)z^n,\qquad |z|<1.
}
\tag{21}
\]
In particular
\[
\boxed{
F_a(0)=1,\qquad
F_a^{(n)}(0)=(-1)^n n!\,c_n(a),\qquad
F_a'(0)=\frac{a q'(a)}{q(a)}<0.
}
\tag{22}
\]
Thus lambda_0=0, or z=0 in a moved local chart, is regular;
setting the coordinate to zero does not set its derivatives
to zero.

**Proof.** Since
a(1+z)+mu=(a+mu)(1+z u_mu),
the source product proves the first equality with L_a.
The geometric series is dominated on |z|<=r<1 by
sum |w_mu|/(1-r), proving the expansion and holomorphy.
The derivative formulas follow by coefficients. The strict
last sign is the retained CWR prelude's q'(a)<0. QED.

Evaluation of an analytic germ at zero is a legitimate but
non-injective readout: it discards its derivative data.
Here the nonzero first derivative in (22) makes that loss
explicit. It is not a contradiction in the source equations
or a failure of the cut relation P^2=P.

The source displacement lambda_0=0 is also not the boundary
a=infinity. The reciprocal anchor epsilon=a_0/a approaches
zero at that boundary; it is a different parameter. The next
theorem treats it explicitly rather than renaming it lambda.

## RFT-T7. The infinite-anchor boundary is positive with zero threshold

The unconditional NIF limit of every fixed polynomial moment is
the arcsine source
\[
d\nu_*(y)=\frac{dy}{\pi\sqrt{y(1-y)}},\qquad
c_n^*=\frac{(1/2)_n}{n!}.
\]
On H_*=L^2(nu_*) define Y_*=y and Hsf_*=-log Y_* with its
multiplication-operator domain. Then
\[
\boxed{
\inf\sigma(\mathsf H_*)=0,\qquad
\ker\mathsf H_*=\{0\},\qquad
\|e^{-t\mathsf H_*}\|=1\quad(t\ge0).
}
\tag{23}
\]
The continuous reflection-positive extension of the limiting
moments is
\[
\boxed{
\psi_*(t)=\int_0^1 y^t\,d\nu_*(y)
=\frac{\Gamma(t+1/2)}{\sqrt\pi\,\Gamma(t+1)}
\sim\frac1{\sqrt{\pi t}}\quad(t\to\infty).
}
\tag{24}
\]
It obeys e^(-t Hsf_*)f->0 strongly for every f in H_* as t
tends to infinity, despite the lack of operator-norm decay.

**Proof.** The pushforward variable x=-log y has density
\[
d\eta_*(x)=
\frac{e^{-x/2}}{\pi\sqrt{1-e^{-x}}}\,dx,\qquad x>0.
\tag{25}
\]
It assigns positive measure to every open subinterval of
(0,infinity), but no atom to zero. Thus the essential range
of multiplication by x has closure [0,infinity), proving the
spectral and kernel statements in (23). The essential supremum
of exp(-t x) is one. Pointwise convergence to zero for x>0,
dominated by |f|^2 in the squared norm, proves strong decay.
The beta integral gives (24), and the Gamma-ratio asymptotic
gives its last expression. As a positive mixture of exp(-t x),
psi_* has the reflection factorization of (16). QED.

If RH holds, the finite-anchor thresholds from (19) satisfy
\[
\boxed{
m_a=\log\left(1+\frac{\gamma_1^2}{a_0}\epsilon\right)
\sim\frac{\gamma_1^2}{a_0}\epsilon
\longrightarrow0,\qquad \epsilon=a_0/a\longrightarrow0.
}
\tag{26}
\]
The limiting source and (23)-(25) themselves require no RH
assumption. NIF supplies convergence of each fixed integer
moment and each fixed finite source bank. No unproved
interchange of all depths, noninteger times and anchor limits
is used here.

Every nonzero polynomial still has strictly positive squared
norm in nu_*. Its native weighted two-channel cut remains the
projection P(f,g)=(yf+(1-y)g,yf+(1-y)g), with P^2=P.
Hence zero decay threshold, loss of a bounded inverse for the
generator, and failure of uniform exponential decay coexist
with a consistent positive source and cut. There is no
logical inconsistency at this boundary.

This limit also does not repair a possible finite-anchor RH
obstruction. MOR's exact carried observer preserves each old
pairing, and its singularity approaches y=1 as the anchor
grows. CWR quantifies the increasing degree needed to retain
that information. Holding a finite observer bank fixed while
taking a to infinity does not implement that full return.

## Result and remaining source theorem

| Question | Result in this continuation |
|---|---|
| Can an exponential observer create a positive rest form from an indefinite two-observer bank? | No; the exact congruence preserves 1-r^2 |
| Are the tested forward wave reflections positive? | No; both the native Z-seed test and the balanced time-label test fail at large anchors |
| Does the native diagram supply a faithful transfer? | Yes; (2I+W+W^(-1))/4=yI and the moment transformation is invertible |
| What is the correct transfer-reflection condition? | Positivity of the actual H_N at every depth, equivalently RH |
| What happens at the existing local source lambda_0=0? | A regular normalized germ with nonzero first derivative |
| Can a zero decay threshold coexist with a consistent positive cut? | Yes; the actual infinite-anchor limiting source proves it |

The new unconditional results select the appropriate reflection
and rule out two incorrect shortcuts: taking the oscillating
wave's reflected form as automatically positive, and identifying
zero local displacement with a singular consistency failure.
The transfer construction also preserves the exact source rather
than choosing a replacement positive model.

The missing theorem remains positivity of the actual transfer
reflection kernel at every depth. Its conditional factorization
(16) must not be used as an unconditional proof of its own
positive weights. Likewise, a finite positive carrier supplies
its declared source horizon, not the all-order completion.

No formal arithmetic theory or proof predicate is constructed
here. The matrix identities, existence of a cut and regularity
of a source chart are not an internal consistency proof of such
a theory, so no conclusion about evading or specializing
Goedel's incompleteness theorem follows from this development.

## Sources, attribution and execution

Retained source lineage:
[LOT](papers/native-layer-observer-transport/THEOREM.md),
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[NIF](RH_INFINITY_CUT_FIXED_POINT.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md),
[QRI](RH_REFLECTION_RATIO_GATE.md),
[CWR](RH_CUT_WAVE_RETURN.md), and
[OQD](RH_OBSERVER_QUOTIENT_DESCENT.md).

The classical reflection-positivity context is Palle Jorgensen
and Feng Tian, *Reflection positivity and spectral theory*,
[arXiv:1705.05262](https://arxiv.org/abs/1705.05262), Section 2
and Section 3. The paper explicitly assumes a positive ambient
Hilbert space and the reflected positivity condition; its
projection and Markov results cannot be invoked solely from
our algebraic P^2=P. The kernel and transfer formulas required
here are proved directly above.

Equation (24) uses the
[NIST beta integral](https://dlmf.nist.gov/5.12#E1) and the
[Gamma-ratio asymptotic](https://dlmf.nist.gov/5.11#E12).
Congruence, the Lorentz quadratic form, polynomial folding
and multiplication-operator calculus are standard tools.
No general novelty claim is made for them. The actual-source
wave-reflection obstruction and the exact identification of
the retained transfer are the applications established here.

Only this note and its README entry change. Existing source
bounds and finite certified carriers are reused. No local
scientific computation, source evaluation, numerical replay,
independent agent or workflow change is performed. These are
written derivations, not new machine-checked certificates.
The commit contains [skip ci]. This remains a draft stacked
on PR #101 without a merge.
