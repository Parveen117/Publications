# NSB1 — Nonlinear spectral blindness and response flatness

Monty Dabas. 5 October 2026. Python 3.12 evidence.

## Result and source contract

The two spectral papers distinguish endpoint spectra, ordered information,
and target-faithful observations. This development applies that distinction
to genuinely nonlinear maps, nonlinear mechanics, and constitutive physics.

The central obstruction is stronger than failure of a first-order
linearization: every finite nonlinear-map jet-pullback **characteristic polynomial** can
forget the nonlinear coefficients. Mixed unmarked spectral words do not
necessarily repair the loss. Coordinate response markers do.

A second obstruction persists over an entire state region. Two strictly
convex energies can have exactly the same Hessian eigenvalue fields while
one response metric is flat and the other curved. A separate conservative
electromagnetic example has identical signed wave speeds at every bias
but different calibrated response curvature. Impedance supplies missing
information in that example.

The T–V–S–P diagram has a precise connection to flatness. Potential
integrability supplies two dual flat affine transports. Their
metric-compatible mean can be curved. Thus Maxwell identities, flat
conjugate transport, and flat response geometry are distinct statements.

The triangular jet mechanism has established Carleman/Koopman lineage.
Dual flatness is established Hessian geometry. The contribution here is
the explicit spectral-observation theorem, its calibrated physics
counterexamples and repairs, and a native compact-gauge consequence.
No exhaustive claim of priority is made.

The first spectral paper already states an arbitrary-map fibre
factorization theorem (T18.1); restating it would not be a new nonlinear
extension. The second paper's linear minimum-probe theorem is used only
on a declared finite jet-coefficient carrier. Its generic finite-section
bound in Section 16 requires the correction recorded in
[SOURCE_AUDIT.md](SOURCE_AUDIT.md). No claim here depends on the uncorrected
bound.

## 1. T1 — Nonlinear coefficients lie in a spectrally invisible ideal

Let F:(R^d,0)->(R^d,0) be a C^r germ fixing zero. Let

\[
 \mathcal J_r=\mathbb R[x_1,\ldots,x_d]/\mathfrak m^{r+1},
 \qquad \mathfrak m=(x_1,\ldots,x_d),
 \qquad D_r=\binom{d+r}{r}.
\]

The jet pullback is the linear operator on polynomial observables

\[
 C_F^{(r)}p=j_r(p\circ F).
 \tag{1}
\]

It represents nonlinear substitution on the finite local algebra. It is
not a claim that the physical state dynamics is linear.

**Theorem.** Put A=DF(0). The degree filtration is preserved, and its
degree-k quotient is acted on by the symmetric power of the algebraic
dual of A. In coordinates this is Sym^k(A^T), not the Hermitian adjoint.
Consequently

\[
 \boxed{\det(zI-C_F^{(r)})
  =\prod_{k=0}^r\det\!\big(zI-\operatorname{Sym}^k(A^T)\big).}
 \tag{2}
\]

Over the complex numbers, if lambda_1,...,lambda_d are the eigenvalues of
A with algebraic multiplicities, these diagonal-block eigenvalues are
lambda^alpha for |alpha|<=r. Diagonalizability is not required.

There is a stronger family statement. If F_1,...,F_m share the fixed
point zero, then for every noncommutative polynomial P with scalar
coefficients,

\[
 \boxed{
 \det\!\left[zI-P(C_{F_1}^{(r)},\ldots,C_{F_m}^{(r)})\right]
 =
 \prod_{k=0}^r
 \det\!\left[zI-P(\operatorname{Sym}^k(A_1^T),
                       \ldots,\operatorname{Sym}^k(A_m^T))\right].
 }
 \tag{3}
\]

All trace words, characteristic coefficients, and source-dependent
determinants built this way are independent of every nonlinear
coefficient of the germs.

**Proof.** If p has no term below degree k and F(0)=0, then p composed
with F also has no term below degree k. In its degree-k part, replacing
one linear factor of F by a nonlinear factor raises the degree.
Therefore the induced map on degree-k homogeneous polynomials is
p(x)->p(Ax), which is Sym^k(A^T).

Choose monomials in increasing total degree. Each C_F is block
triangular with the stated diagonal. Block triangular matrices form an
algebra, and taking their diagonal blocks is an algebra homomorphism.
Its kernel consists of the strictly degree-raising operators; a product
of r+1 such operators vanishes. Hence applying any P preserves the
diagonal formula. Multiplication of block determinants proves (2) and
(3). Triangularizing A over C gives the eigenvalue-product description.
QED.

This connects directly to the source papers: unmarked mixed spectral
products may be complete for a declared matrix trace-word model and
still omit the degree-raising part of a nonlinear jet representation.
The admissible target and representation must be stated.

### What the theorem does not erase

The full matrices retain nonlinear entries. Jordan structure can retain
some of them too. For F(x)=x+a x^2 at order two, all eigenvalues are one,
but C_F-I is zero when a=0 and has rank one when a is nonzero.
Thus “same spectrum” here means characteristic polynomial/eigenvalues
with algebraic multiplicity, not equal Jordan form or equal operators.

The shared fixed point is essential. If F(0) is nonzero, a discarded
high-order term can acquire a nonzero constant after composition;
(1) need not define an operator on this quotient.

No infinite-dimensional Koopman spectral equality is inferred by
letting r increase. Function spaces, domains, convergence and resolvent
control are separate requirements.

## 2. T2 — Exact nonlinear marker repair and an invisible ordered loop

Write the fixed-linear-part jet in Taylor-normalized coefficients:

\[
 F_i(x)=\sum_jA_{ij}x_j+
       \sum_{2\le|\alpha|\le r}f_{i,\alpha}x^\alpha
       \pmod{\mathfrak m^{r+1}}.
\]

There are

\[
 N_{d,r}=d\left[\binom{d+r}{r}-d-1\right]
 \tag{4}
\]

free nonlinear coefficients.

**Theorem.** The coefficient observations

\[
 \mathcal M_{i,\alpha}(C_F)
    =[x^\alpha](C_Fx_i)
    =\frac{\partial^\alpha F_i(0)}{\alpha!}
 \tag{5}
\]

recover the complete r-jet. On the unrestricted real coefficient model
with A fixed, N_(d,r) scalar observations are locally minimal for
recovery with a C^1 observation map and C^1 decoder. The coordinate bank
(5) attains this minimum and has reconstruction norm one in the
declared Euclidean coefficient norm.

**Proof.** C_F x_i=j_r F_i, so (5) reads every unknown coefficient
directly. If G:R^N->R^m and D:G(U)->R^N have a differentiable local
extension with D composed with G equal to the identity, the chain rule
gives DD DG=I_N; hence m>=N. The coordinate bank is the identity on
coefficient space, so noise eta produces exactly the coefficient error
eta. QED.

The normalization and norm matter: raw derivatives include factorials,
and function-value or trajectory error requires radius-dependent
bounds. Fewer probes may suffice for a specified lower-dimensional
physical model or a smaller target. No global minimum on an arbitrary
parameter manifold is asserted.

### A genuinely nonlinear loop with identity jet spectra

For lambda>0 and an integer n>=2, take the diffeomorphisms

\[
 F_\lambda(x,y)=(\lambda x,y),\qquad
 G_t(x,y)=(x,y+t x^n).
\]

The two endpoints obey

\[
 (G_t\circ F_\lambda-F_\lambda\circ G_t)_y
       =t(\lambda^n-1)x^n.
 \tag{6}
\]

They agree through order n-1. For the explicitly ordered point
composition

\[
 K=G_{-t}\circ F_{\lambda^{-1}}\circ G_t\circ F_\lambda
\]

one obtains

\[
 K(x,y)=(x,y+t(\lambda^n-1)x^n),\qquad DK(0)=I.
 \tag{7}
\]

For every r its jet characteristic polynomial is
(z-1)^binom(r+2,2), exactly that of the identity map. For r>=n the
marker [x^n](C_K y)=t(lambda^n-1) detects the surviving loop.
For lambda not equal to one and t nonzero, the loop is not the identity.

This is a nonlinear process-order defect. It is not automatically the
Levi–Civita curvature of a thermodynamic response metric; that requires
a specified transport law.

Even all coefficient jets fail to recover an arbitrary smooth germ:
x and x+exp(-1/x^2), with the perturbation extended by zero at zero,
have identical jets of every order. Analytic germ recovery requires
analyticity and remains a separate stability problem.

## 3. T3 — A nonlinear mechanical distinction invisible to every jet spectrum

Consider the actual classical equations

\[
 \dot q=p,\qquad \dot p=-kq^3,\qquad k\ge0.
 \tag{8}
\]

For every fixed time t, their flow Phi_t fixes the equilibrium (0,0)
and has derivative

\[
 D\Phi_t(0)=
 \begin{pmatrix}1&t\\0&1\end{pmatrix},
\]

independently of k. Therefore every finite state-jet pullback has

\[
 \boxed{\det(zI-C_{\Phi_t}^{(r)})
        =(z-1)^{\binom{r+2}{2}}}
 \tag{9}
\]

for every k, r and finite t. These spectral data do not determine even
Lyapunov stability: k=0 has an unstable free-particle equilibrium,
whereas each k>0 has a stable confining equilibrium.

**Proof.** Linearization at zero gives qdot=p, pdot=0. Equation (2)
then gives (9). For k=0, an arbitrarily small nonzero initial p leads
to unbounded q(t)=q(0)+tp. For k>0 the conserved, positive definite,
coercive energy

\[
 E_k=\frac{p^2}{2}+\frac{kq^4}{4}
\]

has arbitrarily small invariant sublevel neighborhoods of zero.
It gives Lyapunov stability and bounded trajectories. It does not
give asymptotic stability in this conservative system. QED.

### Marked nonlinear response recovers the physical target

Measure z=q^2 along the trajectory. Its Lie derivatives satisfy

\[
 \dot z=2qp,\qquad \ddot z=2p^2-2kq^4,\qquad
 \boxed{E_k=\frac{\ddot z+3kz^2}{4}.}
 \tag{10}
\]

The energy is recovered although the sign-related states
(q,p) and (-q,-p) remain indistinguishable by the entire z-history.

At a crossing state (q,p)=(0,p_0), p_0 not zero, dz vanishes but
D(z,zdot,zddot) has rank two: the last two rows are (2p_0,0) and
(0,4p_0). This is a local time-response observation statement away
from the equilibrium. At (0,0) that rank claim would be false.

If k is itself unknown, it can be recovered when z>0 by

\[
 k=\frac{\dot z^2-2z\ddot z}{4z^3}.
\]

At the crossing state the regular alternative is

\[
 z''(0)=2p_0^2,\qquad z^{(6)}(0)=-72kp_0^4,\qquad
 \boxed{k=-\frac{z^{(6)}(0)}{18[z''(0)]^2}.}
 \tag{11}
\]

These identities follow by repeated differentiation of (8). Their
denominators state the conditioning requirement. They are not
uniformly stable estimators as z or |p_0| tends to zero. Equations
(9) and (11) concern different data: eigenvalues of state-jet
operators versus marked time-response derivatives.

## 4. T4 — An entire nonlinear spectral field can be blind to flatness

Use the native response adapter from NT3/NT4 and NSR1. In a calibrated
affine chart, a C4 energy U with positive Hessian H=D^2U supplies

\[
 X_i=H^{-1}\partial_iH,\qquad
 A_i=\frac12X_i,\qquad
 F_{ij}=-\frac14[X_i,X_j].
 \tag{12}
\]

Because H is a Hessian, this selected connection is Levi–Civita for
the response metric. In two dimensions, with
H=[[a,b],[b,c_H]] and Delta=ac_H-b^2,

\[
 \kappa=-\frac{1}{4\Delta^2}
 \det\begin{pmatrix}
 a&b&c_H\\a_x&b_x&(c_H)_x\\a_y&b_y&(c_H)_y
 \end{pmatrix}.
 \tag{13}
\]

It is curvature over the space of constitutive states, not spacetime
gauge curvature.

### Identical Hessian spectra, flat versus curved

Fix c>0 and the simply connected stable domain x>1/c. Define

\[
 U=\frac c2(x^2+y^2)-x\log x+x,
\]
\[
 V=\frac c2(x^2+y^2)
   +\frac x2\log(1+y^2/x^2)-y\arctan(y/x).
 \tag{14}
\]

Logarithm arguments refer to the declared dimensionless coordinates.
Affine energy terms may be chosen to calibrate a common reference
energy and force at one chosen point without changing any statement
below.

Put r_2=x^2+y^2. Direct differentiation gives

\[
 H_U=\begin{pmatrix}c-1/x&0\\0&c\end{pmatrix},
\qquad
 H_V=\begin{pmatrix}
 c-y^2/(xr_2)&y/r_2\\
 y/r_2&c-x/r_2
 \end{pmatrix}.
 \tag{15}
\]

**Theorem.** At every point of the domain,

\[
 \operatorname{spec}H_U=\operatorname{spec}H_V
      =\{c,c-1/x\},
\qquad
 \boxed{\kappa_U=0,\quad
 \kappa_V=\frac1{4cr_2(cx-1)^2}>0.}
 \tag{16}
\]

Thus even the complete eigenvalue fields, including all their spatial
derivatives, do not determine whether the response connection is flat.

**Proof.** The first matrix is diagonal. The second is

\[
 H_V=cI-\frac{1}{xr_2}
       \begin{pmatrix}y\\-x\end{pmatrix}
       \begin{pmatrix}y&-x\end{pmatrix}.
\]

The rank-one term has eigenvalues 0 and 1/x, proving the shared
spectrum and positivity. H_U is a separated metric; reparametrizing
its first coordinate by integral sqrt(c-1/x) dx and its second by
sqrt(c)y proves flatness. Substituting the derivatives of (15) into
(13) gives the displayed strictly positive curvature for H_V.
Equivalently, its determinant is c(cx-1)/x and the numerator
determinant in (13) is -c/[x^2 r_2]. QED.

At c=1 and (x,y)=(2,0), both spectra are {1,1/2}, but

\[
 F_U=0,\qquad
 F_V=\begin{pmatrix}0&1/32\\-1/16&0\end{pmatrix},
 \qquad\kappa_V=\frac1{16}.
 \tag{17}
\]

For the selected response transport, a small positively oriented
rectangle therefore has return angle
area/(16 sqrt(2))+o(area) in the curved case, and zero in the flat
case. This is a mathematical response-transport protocol; a physical
apparatus must implement that protocol before its return is called a
measured phase.

Because vanishing curvature is invariant under bundle gauge changes,
these two selected connections are not gauge equivalent. Their
response metrics are not locally isometric.

A useful rigidity boundary prevents overgeneralization: if the
Hessian eigenvalues are themselves constant in space, then the
Hessian is constant. Indeed trace(H) and ||H||_F^2 are constant;
Delta_Euclidean H_ij=partial_ij trace(H)=0, so

\[
 0=\Delta_{\rm Euclidean}\|H\|_F^2=2\|D^3U\|_F^2.
\]

On a connected domain U is therefore quadratic plus affine.
The counterexample (14) has a varying eigenvalue c-1/x.

### Same electromagnetic wave speeds across every bias

There is a different, globally uniformly convex example. Set
A=exp(x+y), B=exp(x-y), and sigma in {+1,-1}:

\[
 U_\sigma=x^2+y^2+A+B+
        \frac{\sigma}{2}(x^2-y^2)-2-2x.
 \tag{18}
\]

Both energies and gradients vanish at zero. Their Hessians are

\[
 H_\sigma=
 \begin{pmatrix}
 2+A+B+\sigma&A-B\\A-B&2+A+B-\sigma
 \end{pmatrix}\ge I.
\]

The inequality follows from the positive diagonal base and the
two nonnegative rank-one exponential contributions. Their Hessian
eigenvalues are identical at every bias:

\[
 2+A+B\ \pm\sqrt{1+(A-B)^2}.
\]

In the closed conservative transverse Maxwell adapter of NSR1,
the signed dimensionless wave speeds are also identical:

\[
 v_\pm=A-B\pm\sqrt{(2+A+B)^2-1}.
 \tag{19}
\]

Nevertheless, with Delta=4(1+A)(1+B)-1,

\[
 \boxed{\kappa_\sigma=\frac{2\sigma AB}{\Delta^2}.}
 \tag{20}
\]

At zero the spectra are {5,3}, the speeds are +/-sqrt(15), and
the two curvatures are +/-2/225.

**Proof.** The trace and determinant of H_sigma do not depend on
sigma. The Maxwell characteristic polynomial is
(v-b)^2-ac_H, also independent of sigma. In (13), the derivative
rows are (A+B,A-B,A+B) and (A-B,A+B,A-B). The determinant numerator
is -8 sigma AB, proving (20). QED.

The missing calibrated marker is impedance:

\[
 Z_\sigma^2=\frac{2+A+B+\sigma}{2+A+B-\sigma},
 \qquad Z_- = Z_+^{-1}.
 \tag{21}
\]

Speeds give s=(v_+-v_-)/2 and b=(v_++v_-)/2. A positive impedance
measurement then gives a=sZ and c_H=s/Z. Measurements across two
bias directions recover the full H field and its curvature through
(12), or the existing coupling/impedance formula

\[
 f=-\frac{dr\wedge d\log Z}{2(1-r^2)^{3/2}},
 \qquad r=b/\sqrt{ac_H},\qquad \kappa=-f_{xy}/\sqrt\Delta.
 \tag{22}
\]

For this two-member model, one calibrated impedance measurement
already separates sigma. No universal one-probe theorem is asserted.

The two exponential connections are conjugate by a constant
reflection on the response fibres, so their unordered loop spectra
agree. The statement is about calibrated constitutive geometry;
they must not be called gauge-inequivalent connections. Their
Gaussian response curvatures have opposite signs because the
identification with coordinate tangent directions is retained.
The flat/curved pair (14) supplies the stronger gauge-inequivalent
example.

A full conservative medium can be declared by taking (18) for
(D_y,B_z) and adding positive quadratic energies for omitted
components. At zero omitted fields their conjugates vanish, so
this transverse sector is closed under the source-free Maxwell
equations. These are specified constitutive models, not a claim of
Lorentz-invariant vacuum electrodynamics, a measured material, or
a universal causal-speed bound.

## 5. T5 — The precise nonlinear flatness statement behind the diagram

The word “symmetry” must identify both the operation and the transport.
For a C1 response j(x), define the process one-form beta=j_i dx^i.
Then

\[
 \boxed{d\beta=0\iff \partial_i j_j=\partial_j j_i.}
 \tag{23}
\]

On a star-shaped domain about zero this is equivalent to

\[
 j=\nabla\Psi,\qquad
 \Psi(x)=\int_0^1j(tx)\cdot x\,dt.
 \tag{24}
\]

**Proof.** The coefficient of dx^i wedge dx^j in d beta is the
difference in (23). If it vanishes, differentiate (24) with respect
to x^k, use partial_k j_i=partial_i j_k, and integrate
d[t j_k(tx)]/dt. The result is j_k(x). The converse follows from
equality of mixed derivatives. QED.

For conservative conjugates, beta=dU has energy units. For a
dissipative force–flux law, an admitted Psi has the units of its
dissipation potential. Their exactness mathematics is shared;
their physical meanings are not interchangeable.

If j=L(x)x, the relevant differential response is

\[
 \partial_i j_j=L_{ji}
       +\sum_k(\partial_iL_{jk})x^k.
\]

Symmetry of the secant L(x) alone does not imply (23). Nor is
(23) an unrestricted nonlinear Onsager theorem. For example,
L(b)=ell I+bR with ell>0 and R=[[0,-1],[1,0]] obeys
L(b)=L(-b)^T and nonnegative production ell|x|^2, but

\[
 d[(L(b)x)\cdot dx]=2b\,dx\wedge dy.
\]

Field-reversal reciprocity need not be flatness at fixed nonzero
field. Microscopic nonlinear transport reciprocity has additional
fluctuation and time-reversal hypotheses.

### Two flat transports and their potentially curved mean

For a C4 strictly convex potential U on a stable affine chart, put
p=grad U and H=D^2U. Define

\[
 \nabla^0=d,\qquad
 \nabla^*=d+X,\qquad X=H^{-1}dH.
 \tag{25}
\]

**Theorem.** These are torsion-free flat connections, dual for H.
Conjugate coordinates p are affine for nabla-star. Their mean is
the selected Levi–Civita response connection,

\[
 \nabla^{LC}=\tfrac12(\nabla^0+\nabla^*)
            =d+\tfrac12X,\qquad
 \boxed{F^{LC}_{ij}=-\tfrac14[X_i,X_j].}
 \tag{26}
\]

Thus both conjugate affine descriptions may be flat while the
response metric is curved.

**Proof.** The Hessian condition makes
(X_i)^k_j=H^{ka}U_aij symmetric in i,j, proving torsion-freeness.
The identity dX+X wedge X=0 proves flatness of nabla-star;
nabla-zero is flat by definition. The compatibility relation

\[
 \partial_i H=H X_i
\]

is precisely the H-duality identity between the zero connection
and X. Since dp=H dx, the inverse coordinate fields

\[
 Y_i=\sum_k(H^{-1})_{ki}\partial_{x^k}
\]

are the coordinate fields partial/partial p_i and commute. Direct
differentiation gives nabla-star Y_i=0. The mean is torsion-free
and metric compatible, so is Levi–Civita; inserting X/2 into the
curvature formula gives (26). QED.

This is the established dual-flat Hessian structure expressed in
the programme's native response connection. For a general positive
matrix field H that is not a Hessian, the selected matrix connection
exists but the tangent Levi–Civita identification needs separate
justification.

### An exact TVSP-compatible control experiment

Take

\[
 U(x,y)=\frac{x^2}{2}+\frac{(1+x)y^2}{2},\quad
 H=\begin{pmatrix}1&y\\y&1+x\end{pmatrix},\quad
 \Delta=1+x-y^2>0.
 \tag{27}
\]

All potential/Maxwell identities hold. The paired conjugate fields

\[
 Y_1=\frac{(1+x)\partial_x-y\partial_y}{\Delta},\qquad
 Y_2=\frac{-y\partial_x+\partial_y}{\Delta}
\]

commute. At the origin the selected response metric nevertheless has
F_xy=-R/4 and kappa=1/4.

Now choose a different, explicitly mixed protocol:
X=(partial_x)_y and Y=(partial_q)_x, with q=U_y=(1+x)y.
Then

\[
 X=\partial_x,\qquad Y=(1+x)^{-1}\partial_y,\qquad
 [X,Y]=-(1+x)^{-2}\partial_y.
\]

The point path +aX,+bY,-aX,-bY has the exact displacement

\[
 \boxed{\Delta y=
  \frac b{1+x+a}-\frac b{1+x}
  =-\frac{ab}{(1+x)(1+x+a)}.}
 \tag{28}
\]

The loop must remain in the stated stable chart. It uses controls
from different coordinate charts and therefore does not contradict
commutation of the paired Y_1,Y_2.

| Structure | Exact flatness condition | What does not suffice |
| --- | --- | --- |
| Process one-form beta | d beta=0; locally a response potential | Symmetry of a variable secant L alone |
| Paired conjugate affine transport | Hessian integrability and the chosen conjugate coordinates | Arbitrarily mixing controls from different charts |
| Selected response-metric transport | [H^-1 H_i,H^-1 H_j]=0 for every pair | Maxwell relations or rotation symmetry alone |
| A nonlinear ordered flow | Its declared loop map is the identity | Identity eigenvalues of every finite jet loop |

The rotationally invariant nonlinear energies below give another
explicit check that symmetry can coexist with response curvature.
Symmetry restricts curvature and observability; it does not
universally set every curvature to zero.

## 6. Full finite jets and neighborhood flatness remain different targets

For any requested finite Hessian-jet depth N, choose an integer m>=2
with 2m-2>N and eta>0. The globally strictly convex energy

\[
 U_\eta(q)=\frac{|q|^2}{2}+\frac{\eta}{2m}|q|^{2m}
 \tag{29}
\]

has H>=I. Every derivative of H at zero through order N agrees
with the quadratic model. Yet in two dimensions, at (t,0),

\[
 \kappa_\eta(t,0)=
 -\frac{2(m-1)^3\eta^2t^{4m-6}}
 {[1+(2m-1)\eta t^{2m-2}]^2
  [1+\eta t^{2m-2}]^2}<0\quad(t\ne0).
 \tag{30}
\]

**Proof.** The Hessian correction has degree 2m-2, proving the jet
claim. Differentiation gives tangential stiffness
1+eta|q|^(2m-2) and radial stiffness
1+(2m-1)eta|q|^(2m-2), both at least one.
Inserting their third derivatives into (13) yields (30). QED.

Therefore no fixed finite list of local response derivatives
establishes flatness on an entire neighborhood without further
model-class assumptions. This differs from T1: even full finite
jets may miss still higher terms, whereas T1 says eigenvalue-only
compression forgets nonlinear coefficients already present in
the retained jet.

## 7. Return to Yang–Mills: a precise nonlinear gauge consequence

The physics results suggest retaining source-marked nonlinear
response and its quantitative conditioning. They do not imply
that an exact quantum Hamiltonian spectrum can hide its own gap.

A compact gauge example makes the actual limitation rigorous.
On one four-link SU(2) plaquette, use the established normalization

\[
 L_EW=3W,\qquad \Gamma(W)=1-W^2.
\]

For integer N>=0 put s=1-W^2 and f_N=1-s^(N+1). Define the
separate parent Hamiltonian

\[
 K_{t,N}=L_E-\frac t2L_Ef_N+\frac{t^2}{4}\Gamma(f_N),
 \qquad t\ge0.
 \tag{31}
\]

Its exact positive vacuum is proportional to exp(t f_N/2);
its ground energy is zero. Direct radial differentiation gives

\[
 L_Ef_N=2(N+1)s^N[(2N+4)W^2-1],
\]
\[
 \Gamma(f_N)=4(N+1)^2W^2s^{2N+1}.
 \tag{32}
\]

Near the trivial link configuration, native exponential coordinates
give s=O(|A|^2). For any N>=1, the multiplication potential in
(31) has no Taylor terms below degree 2N. All its local jets
through order 2N-1 agree exactly with the free electric operator,
for every value of t.

Nevertheless,

\[
 \boxed{\operatorname{gap}_{\rm Gauss} K_{t,N}
 \le162e^{-t c_N},\qquad
 c_N=(3/4)^{N+1}-(1/4)^{N+1}>0.}
 \tag{33}
\]

**Proof.** The ground-state identity follows by differentiating
exp(t f_N/2), or the exact positive-square factorization.
The odd gauge-invariant test g=clip(2W,-1,1) has zero mean in
the even vacuum density. Its Dirichlet numerator is at most

\[
 4\exp\!\left(t[1-(3/4)^{N+1}]\right)/Z.
\]

On |W|>=sqrt(3)/2 its squared value is one and f_N is at least
1-(1/4)^(N+1). The native Haar cap has mass at least 2/81, as
proved in YM82. The variance is therefore at least

\[
 (2/81)\exp\!\left(t[1-(1/4)^{N+1}]\right)/Z.
\]

Their Rayleigh quotient is (33). Smooth approximation makes the
clipped trial legitimate in the closed Dirichlet form. QED.

For N=2, the local potential jets through cubic order agree while
the gap is at most 162 exp(-13t/32). Choosing N larger defeats
any specified finite local jet depth.

This is an actual quantum compact-gauge counterfamily, but it is
not the original homogeneous Yang–Mills Hamiltonian. The exact
vacuum formula applies to one active plaquette or edge-disjoint
copies, not an overlapping all-plaquette sum. It shows why
local spectral or response jets cannot alone provide a uniform
mass-gap certificate.

[YANG_MILLS_DOCK.md](YANG_MILLS_DOCK.md) states the next admissible
target: marked estimates for the actual vacuum and source-resolved
spectral weight, with volume/cutoff-uniform margins and the
continuum field-identification requirements retained.

## 8. Evidence and scope

The Python 3.12 packet checks exact rational polynomial identities,
finite pullback and Lie-derivative matrices, characteristic
polynomials, observer coefficients, native curvature through
independent formulas, nonlinear control loops, gauge-parent
identities, outward exponential bounds and source integrity.

These finite checks do not machine-formalize the universal proofs.
No infinite Koopman completion, experimental material validation,
spacetime metric, photon mass, or four-dimensional Yang–Mills
continuum gap is asserted.

~~~sh
python3.12 -B papers/nonlinear-spectral-blindness/certificates/nsb1_nonlinear_blindness.py --check
python3.12 -B -m unittest discover -s papers/nonlinear-spectral-blindness/tests -v
~~~

## Sources and attribution

1. Monty Dabas, *When Spectra Forget Order*, SPECTRAL 1.pdf,
   T5.1/T5.4, T15–18, and T19.1. Its matrix-order blindness,
   marker stability and existing set-level fibre theorem are
   recorded in the source audit.
2. Monty Dabas, *Stable Recovery Beyond Spectral Blindness*,
   SPECTRAL_2.pdf, T7.2/T7.5 and Section 16. Use the explicit
   Section 16 correction in the source audit.
3. [Native thermodynamic curvature](../native-thermodynamic-curvature/THEOREM.md),
   NT3/4/6/7; [NSR1](../native-symmetry-response/THEOREMS.md),
   T2–T5; and [NSR2](../native-symmetry-response/NSR2_CONJUGATE_DIAGRAM_PROPERTIES.md).
   These define the response connection, orientation, physical
   work conventions and diagram transfer.
4. A. Amini, C. Zheng, Q. Sun and N. Motee, *Carleman
   Linearization of Nonlinear Systems and Its Finite-Section
   Approximations*, [arXiv:2207.07755](https://arxiv.org/abs/2207.07755).
   Established lifting lineage; convergence requires additional
   assumptions beyond the finite algebra used here.
5. H. Shima, *Hessian manifolds of constant Hessian sectional
   curvature*, J. Math. Soc. Japan 47 (1995),
   [primary paper](https://www.jstage.jst.go.jp/article/jmath1948/47/4/47_4_735/_pdf).
   Dual-flat Hessian geometry is established. Shima's separately
   named Hessian curvature tensor is not silently identified with
   the selected Levi–Civita curvature here.
6. L. Mezincescu, J. G. Russo and P. K. Townsend, *Hamiltonian
   birefringence and Born–Infeld limits*,
   [arXiv:2311.04278](https://arxiv.org/html/2311.04278v3).
   Source for the admitted Hamiltonian Maxwell convention, already
   used in NSR1.
7. K. Saito and Y. Utsumi, *Symmetry in full counting statistics,
   fluctuation theorem, and relations among nonlinear transport
   coefficients in the presence of a magnetic field*,
   [arXiv:0709.4128](https://arxiv.org/pdf/0709.4128).
   Context for the distinction between an assumed nonlinear
   response potential and microscopic nonlinear reciprocity.
8. [YM82](../yang-mills-certified-benchmark/YM82_NONLINEAR_LATTICE_SPECTRAL_GAPS.md),
   for the native plaquette identities, Haar cap, and the explicitly
   separate quantum parent construction extended in (31)–(33).
