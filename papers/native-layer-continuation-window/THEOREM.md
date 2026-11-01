# LCW-1 — Sharp continuation of the actual single-point layer tower

17 October 2026. Continues ILC-1 in the same retained xi-source framework.
The new result is an exact two-layer continuation gate and its shrinking
width law. Actual source derivatives 16 and 17 are evaluated, rather than
chosen to satisfy that gate. No all-depth source sign theorem is claimed.

Retain Xcal(w^2)=xi(1/2+w), q=Xcal'/Xcal, a=3249/4 and
c_j=(-a)^j q^(j)(a)/(j!q(a)). Thus c_0=1. All indices below are
derivative depth; the ILC displacement lambda=(t-a)/a is unchanged.

## LCW-T1. The next even layer has a sharp algebraic window

Suppose real c_0,...,c_(2n-1), c_0=1, have both strict forms

\[
A=(c_{i+j+1})_{i,j<n}>0,\quad
B=(c_{i+j}-c_{i+j+1})_{i,j<n}>0,\quad H=A+B.
\]

Set h=(c_n,...,c_(2n-1))^T. For n>1 define

\[
C=(c_{i+j+1}-c_{i+j+2})_{i,j<n-1},\quad
d=(c_n-c_{n+1},\ldots,c_{2n-2}-c_{2n-1})^T,
\]
\[
L=h^TH^{-1}h,\quad U=c_{2n-1}-d^TC^{-1}d,\quad W=U-L.
\tag{1}
\]

For n=1 the d/C contribution is zero. Then C>0 when present, W>0,
and, writing e for the last coordinate vector,

\[
\boxed{W^{-1}=e^T(A^{-1}+B^{-1})e.}
\tag{2}
\]

For a proposed x=c_(2n), the enlarged H_(n+1) is positive exactly when
x>L. The enlarged mixed form C_n is positive exactly when x<U.

**Proof.** Let E_0=[e_0,...,e_(n-2)] and
E_1=[e_1,...,e_(n-1)], with empty matrices for n=1. Put

\[
G=\operatorname{diag}(A,B),\quad
R=\binom II,\quad T=\binom{E_0-E_1}{-E_1},\quad k=\binom e0,
\quad Q=[R\ T\ k].
\]

Q is invertible: subtracting its lower output block from its upper
recovers the E_0 and e coordinates, then the common coordinate follows.
The Hankel identities HE_1=AE_0 and Ae=h give

\[
Q^TGQ=\begin{pmatrix}H&0&h\\0&C&d\\h^T&d^T&c_{2n-1}\end{pmatrix}>0.
\tag{3}
\]

Its last Schur complement is W, proving C>0 and W>0. The last
coordinate of Q^-1 applied to (u,v) is e^T(u-v). Computing the last
entry of the inverse in (3) both ways proves (2). Finally the two
enlarged forms are [[H,h],[h^T,x]] and
[[C,d],[d^T,c_(2n-1)-x]]; their Schur complements prove the claim.
For n=1 the latter is the scalar c_1-x. QED.

## LCW-T2. The odd layer splits a conserved positive budget

For proposed x=c_(2n), define

\[
v(x)=(c_{n+1},\ldots,c_{2n-1},x)^T,\quad w(x)=h-v(x),
\]
\[
\ell(x)=v^TA^{-1}v,\qquad u(x)=x-w^TB^{-1}w.
\tag{4}
\]

For n=1, v=(x). Then the exact identity, valid even outside the window, is

\[
\boxed{u(x)-\ell(x)=\frac{(x-L)(U-x)}{U-L}.}
\tag{5}
\]

For proposed z=c_(2n+1), the two new forms are strictly positive iff

\[
\boxed{L<x<U,\qquad \ell(x)<z<u(x).}
\tag{6}
\]

Thus these are sharp admissible intervals, not merely sufficient bounds.
Closed inequalities characterize semidefinite extensions of the strict
old forms. Endpoints do not give a new strict direction.

**Proof.** Put b=H^-1h. The first n-1 rows of A are the following rows
of H, while the last row of A is h^T. Hence
v-Ab=(x-L)e. Also w=Bb-(x-L)e. Completing the squares gives

\[
v^TA^{-1}v+w^TB^{-1}w=L+(x-L)^2/W.
\]

Subtract this from x to obtain (5). The enlarged forms are

\[
A_{n+1}=\begin{pmatrix}A&v\\v^T&z\end{pmatrix},\quad
B_{n+1}=\begin{pmatrix}B&w\\w^T&x-z\end{pmatrix}.
\]

Their Schur innovations are alpha=z-ell and beta=u-z. Both must
be positive. Their sum is (5), which is positive precisely for L<x<U.
This proves (6), and the same Schur argument proves the closed version.
QED.

The even innovation x-L is the squared norm of the residual of y^n
after projecting onto polynomials of degree below n in H_(n+1).
This follows by completing that Gram square. It measures a newly
independent polynomial observation. Neither positivity of individual
c_j nor a finite triangle of binomial readings substitutes for (6).

## LCW-T3. Exact width transport and the all-depth obligation

Within a strict extension define the dimensionless positions

\[
p=(x-L)/W,\qquad \theta=(z-\ell)/(u-\ell).
\tag{7}
\]

Both lie in (0,1). The odd budget is W p(1-p), and its two
innovations are that budget times theta and 1-theta. If W_next
denotes the next even window, then

\[
\boxed{W_{\rm next}=W p(1-p)\theta(1-\theta)\le W/16.}
\tag{8}
\]

**Proof.** Inverse block elimination gives last inverse entries
1/alpha and 1/beta for the enlarged A and B. Apply (2) at n+1:
W_next=alpha beta/(alpha+beta), proving (8). The inequalities
p(1-p)<=1/4 and theta(1-theta)<=1/4 are elementary squares. QED.

Starting with 0<c_1<1, the first even width is c_1(1-c_1).
It follows inductively that the admissible width for c_j is at most
4^(-(j-1)), for every j>=1 along any strict continuation. This is
an upper bound, not an asymptotic equality or a uniform positive
margin for the actual source. Certification precision must resolve
the actual distances to both endpoints, which may be much smaller.

Applied to the exact actual source at a, condition 0<c_1<1 followed
by (6) at every n is equivalent to all-order A_n,B_n>0 and hence,
by ILC-T2, to RH. This statement is a criterion. LCW proves (8)
for accepted data; it does not prove that the actual coefficients
remain accepted at every n. Choosing interior positions can construct
another positive moment sequence, but cannot identify that sequence
with the derivatives of xi.

The variables p and theta in (7) are classical canonical moment
positions. They are not primes, operator eigenvalues, the ILC
displacement lambda, or physical phase angles. The compact moment
coordinates and their theory are credited to the classical literature;
see Dette and Nagel, *Annals of Probability* 40 (2012), 2690–2704,
Section 2, equations (2.1)–(2.2),
<https://arxiv.org/pdf/1007.3369>. The finite matrix identities needed
here have been proved above; no new priority for canonical moments
or the Hausdorff moment problem is claimed.

## LCW-T4. Accepted layers extend the cut and the response

An accepted pair constructs H_(n+1)>0 and the acquired self-adjoint
contraction Y_(n+1)=H_(n+1)^-1 A_(n+1), with 0<Y_(n+1)<I.
ILC-T3 then supplies its exact two-channel projection and unitary
two-reflection wave. Thus (6) is an executable gate for extending
that construction by one polynomial direction. The top polynomial
compression of Y_(n+1) recovers Y_n, since both pairings on the old
subspace are unchanged. The dilated cuts themselves are not asserted
to be literal principal submatrices: compression need not commute
with the positive square root.

Every accepted finite pair has a finite positive spectral measure
reproducing c_0,...,c_(2n+1), by ILC-T4's moment matching proof.
Its support is in (0,1), so all binomial readings through that depth
are positive. This finite representing measure is not identified
with the actual zero measure. Under acceptance at all orders, ILC's
strong limit of padded cuts applies without any extra uniform gap.

The new response matches two additional actual Taylor coefficients.
Its unfitted-source error is bounded using ILC-T4's unconditional
finite-source Schwarz estimate with exponent 2(n+1). No RH premise
is used for that finite analytic estimate.

## LCW-T5. Actual source execution and next unresolved layer

The unchanged complete source evaluator is run at original s=29
through derivative order 17. All sixteen overlapping old moments
agree with their retained outward enclosures. The actual pair
c_16,c_17 satisfies (6), and independent direct LDL factorization
certifies both 9x9 forms. All 171 binomial readings r+k<=17 are
positive. Gamma, Euler-prime and omitted-tail terms are all retained.

The exact next even window for c_18 is computed from those accepted
data and checked against (8). This does not evaluate c_18 or assert
its membership. The positive nine-direction response is also checked
against the previously independently evaluated, unfitted s=31 value.
See SOURCE_EVALUATION.md and CERTIFICATE.json for bounds and controls.

The next source-specific theorem needed is control of the exact xi
coefficients within (6), uniformly through arbitrary depth. The width
law alone supplies no such theorem. Physical mass, Maxwell dynamics,
spacetime geometry, all-order RH and the old global eta are unchanged.
