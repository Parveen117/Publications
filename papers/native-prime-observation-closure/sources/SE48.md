# SE48: a spectral parameter in the five-source operator and actual diagonal zero certificates

13 October 2026 (Asia/Kolkata). Author: Monty.

SE47's neutral comparison state is not automatically a zero of zeta.
This stage adds the missing spectral parameter using the actual source
response, retains the existing 1+4 architecture, and separately certifies
five local zeros of the actual zeta function. An all-order diagonal
criterion is proved, but its universal positivity hypothesis remains open.

Typed Cut, derived scalars and iota, lawful composition and retained
readings precede the admitted source representation. SE21/SE39 already
identify that representation with the completed xi response. Complex
analysis is used explicitly through that adapter. Neither the classical
zeta definition nor its analytic identity is replaced by a new meaning
of zero. The formal marked unit Lambda, the source shift eta and the
spectral parameter mu below have different types and roles.

## 1. Theorem 1: acquire a spectral pencil from complete source readings

Write the SE39 entire function as Xcal(t), so that

\[
 \mathcal X(w^2)=\xi(1/2+w),\qquad
 q(t)=\frac{\mathcal X'(t)}{\mathcal X(t)}
     =\frac{Q(1/2+\sqrt t)}{2\sqrt t},\qquad
 \Phi(t)=2tq(t).                                      \tag{48.1}
\]

Xcal is entire; q is meromorphic and each zero a of Xcal is a pole
of q with positive integer residue equal to its multiplicity. The
positive-real observations use t_i=(i+3/2)^2, i=1,...,N, and the same
complete Gamma, prime and pole data as SE21. More general distinct
positive nodes in (9/4,infinity) are allowed in the written theorem.
Let u_i=q(t_i), T=diag(t_i), and define symmetric matrices

\[
 \begin{aligned}
 L_{ij}&=-q[t_i,t_j],& L_{ii}&=-q'(t_i),\\
 J_{ij}&=(tq(t))[t_i,t_j],& J_{ii}&=q(t_i)+t_iq'(t_i).
 \end{aligned}                                        \tag{48.2}
\]

The actual source matrix is exactly S_N=4J_N. L is a second acquired
response form, not an externally imposed coordinate metric. Its sign
must also be read or proved. Direct multiplication gives

\[
 \boxed{J+LT=u\mathbf1^T,\qquad J+TL=\mathbf1u^T.}      \tag{48.3}
\]

**Theorem.** If the actual finite L and J are positive definite, then

\[
 \mathfrak A_N=L^{-1}J,\qquad
 q_N(z)=u^T(J+zL)^{-1}u                              \tag{48.4}
\]

give an L-self-adjoint positive operator and a rational response with
N distinct simple negative poles. More precisely,

\[
 q_N(z)=\sum_{k=1}^N\frac{w_k}{z+\mu_k},\quad
 \mu_k>0,\quad w_k>0,\qquad
 q_N(t_i)=q(t_i),\quad q_N'(t_i)=q'(t_i).              \tag{48.5}
\]

The corresponding locations in the original s coordinate are

\[
 \boxed{s=\tfrac12\pm\iota\sqrt{\mu_k}.}              \tag{48.6}
\]

These are poles of the finite response approximant, not yet certified
zeros of the actual xi. In particular their residues w_k need not be
integers, and replacing q_N by a determinant's logarithmic derivative
would discard those readings.

**Proof.** L A=J=A^T L. The positive symmetric matrix
B=L^(-1/2) J L^(-1/2) is the admitted orthonormal chart of the same
acquired pairing. Its eigenvalues mu are positive. Spectral expansion
of (48.4) gives nonnegative residues.

For a generalized eigenvector Jv=mu Lv, (48.3) implies
(T+mu I)Lv=1(u^T v). If u^T v were zero, positivity of all t_i+mu
and invertibility of L would give v=0. Every eigenspace therefore has
dimension one: a larger one would contain a nonzero vector annihilated
by u^T. No pole has zero residue. At a probe, (J+t_i L)e_i=u, so the
value in (48.4) is u_i and its derivative is -L_ii=q'(t_i).
Finally z=(s-1/2)^2=-mu gives (48.6). QED.

Under an invertible source-coordinate change B_0, transport
L,J,u to B_0^T L B_0, B_0^T J B_0, B_0^T u. Then A transforms by
similarity and q_N is unchanged. In a complex chart replace transpose
by the admitted conjugate transpose. Changing the operator without
its probe is not this transport.

The finite pencil uses N internal source states. Its reduction below
has five output labels. These are different counts.

## 2. Theorem 2: retain the existing five-label architecture

On V_N retain exactly SE47's fixed accepted space
E=span(f_2,f_3,f_4,f_5), boundary functional and moment complement.
Transport both S_N and L_N by that same coordinate change and form

\[
 \mathcal P_N(\mu)=S_N-4\mu L_N
       =\begin{pmatrix}A_\mu&C_\mu^\dagger\\C_\mu&D_\mu\end{pmatrix},
 \quad R_\mu=D_\mu+bb^\dagger,\quad
 W_\mu=[b,C_\mu],\quad G_\mu=\operatorname{diag}(1,A_\mu).
                                                               \tag{48.7}
\]

Whenever R_mu is invertible, set

\[
 \boxed{K_{5,N}(\mu)=G_\mu-W_\mu^\dagger R_\mu^{-1}W_\mu.}        \tag{48.8}
\]

At mu=0 this is SE47's actual five-source operator, with no replacement
of its boundary column or accepted space. For real mu all displayed
forms are Hermitian; for complex mu use the polynomial block extensions
from those real coefficient matrices, not conjugation of the variable.

**Theorem.** In the SE47 coordinate change, whose determinant has
absolute value one, the meromorphic determinant identity is

\[
 \boxed{\det R_\mu\det K_{5,N}(\mu)
       =\det\mathcal P_N(\mu)
       =4^N\det(J_N-\mu L_N).}                       \tag{48.9}
\]

At a real mu with invertible R_mu, ker K_(5,N)(mu) and
ker P_N(mu) are isomorphic. A root of det R_mu alone is not a source
zero; the characteristic expression must retain the cancellation in
(48.9). Positivity of A_mu or R_mu is not a premise of this determinant
identity, although it remains necessary for the corresponding positive
Schur sign test.

**Proof.** Introduce the augmented form
P_N(mu)[e,y]+|t+b^dagger y|^2. Its unit-determinant triangular
change gives diag(P_N(mu),1). In the order (y,t,e), eliminate
R_mu to obtain diag(R_mu,K_(5,N)(mu)). This proves (48.9) where
the inverse exists, and clearing denominators extends the polynomial
identity. A five-label null vector z=(t,e) lifts by
y=-R_mu^(-1)W_mu z; the boundary equation is t=-b^dagger y.
These equations give mutually inverse kernel maps. QED.

For example P(mu)=diag(1-mu,...,5-mu), b=1 on its complementary
state, has R_mu=6-mu and K_00=(5-mu)/(6-mu). The apparent reference
pole at 6 cancels; the actual source root at 5 remains. A fixed 5 by 5
meromorphic readout can encode more than five roots as its internal
source space grows. This is not a fixed linear pencil of degree five.

### The actual endpoint cannot have a bounded positive reference inverse

The eta-zero full-complement reference in SE47 is necessarily
noncoercive. This follows from SE5 and does not assume RH. Take its
unit weakly escaping f_j with ||S_b f_j||->0. The fixed moment
projection Pi has finite rank, so Pi f_j->0. Thus y_j=(I-Pi)f_j has
norm tending to one, is weakly null in Y, and S_b y_j->0. Also
ell(y_j)->0. Consequently the restricted pole-restored R satisfies
R y_j->0. It cannot satisfy R>=c I_Y for any c>0.

SE47's finite/coercive theorem remains valid within its hypotheses.
Its ordinary bounded inverse cannot simply be extended to the global
eta-zero complement. The endpoint route must use SE1's energy-domain
inverse or a controlled limit of finite readings. The next theorem
provides such a limit contract without postulating a bounded R^(-1).

## 3. Theorem 3: an all-order observer criterion excludes every off-line zero

Choose distinct nodes t_j>9/4 accumulating at an interior positive
point t_*. For example t_j=4+1/j. At each N form (48.2) from the
**actual** response q, not from the previously fitted rational model.

**Theorem.** RH is equivalent to L_N>0 and J_N>0 for every N on this
one fixed nested node sequence. Under these conditions q_N converges
locally uniformly on Omega=C\(-infinity,0] to the actual meromorphic
response q, which is therefore holomorphic on Omega. No bounded global
Schur-reference inverse is required. The universal positivity premise
is not established by the finite certificates in this stage.

**Proof of the exclusion direction.** By Theorem 1, q_N has positive
weights and support -mu_k<0, and q_N(t_1)=q(t_1). For every compact
K in Omega, the finite constant

\[
 C_K=\sup_{z\in K,\,\mu\ge0}\frac{t_1+\mu}{|z+\mu|}
\]

exists: the ratio tends uniformly to one at infinity, and is continuous
on the remaining compact set with a nonzero denominator. Therefore

\[
 |q_N(z)|\le C_K\sum_k\frac{w_k}{t_1+\mu_k}
           =C_Kq(t_1).                                \tag{48.10}
\]

Cauchy estimates and compact exhaustion give a locally uniformly
convergent subsequence with a holomorphic limit. Each fixed probe
eventually belongs to every interpolant; the limit agrees with q at
all t_j. The interior accumulation point and the analytic identity
theorem identify the limit with q near t_* and then meromorphically
throughout the connected Omega. Thus q has no pole there. Every
subsequence has the same limit, proving full local uniform convergence.

An off-line xi zero rho=beta+iota gamma gives
a=(rho-1/2)^2 with imaginary part 2(beta-1/2)gamma !=0. The actual
q has residue equal to the positive multiplicity at a, so its pole
cannot cancel. This contradicts holomorphy in Omega. Xi has no real
zeros and Xcal(0)!=0 by SE39; the only possible zero support is the
negative real axis. Hence every nontrivial zeta zero has real part 1/2.

**Converse.** Under RH, SE39's genus-zero product gives
q(t)=sum_gamma m_gamma/(t+gamma^2). Hence L and J are the complete
Cauchy Grams with respective positive weights m_gamma and
m_gamma gamma^2. Their series converge. A zero quadratic reading
would make a rational function sum_i c_i/(t_i+x) vanish at infinitely
many distinct gamma^2, forcing every coefficient to vanish. Both
finite Grams are therefore strictly positive. QED.

This is a criterion and a proved convergence implication, not an
unconditional convergence proof. The executed integer-rate nodes do
not accumulate in the interior; their universal sign criterion instead
uses SE21's density and SE39. Either route still needs an all-order
arithmetic sign argument. Symmetry or a formal zero Cut does not supply
it. Standard analytic compactness and continuation are explicit inputs
of the admitted adapter; they are not new primitive RKF axioms.

## 4. Theorem 4: certify an actual zero from the pole-guided reading

The finite poles provide probe centres, whose validity is then decided
by independent actual-zeta readings. No known zero catalogue is input.
Use the Euler-transformed alternating series E(s)=(1-2^(1-s))zeta(s):

\[
 E_K(s)=\sum_{k=0}^{K-1}2^{-k-1}
       \sum_{j=0}^k(-1)^j\binom{k}{j}(j+1)^{-s}.        \tag{48.11}
\]

SE21's Mellin remainder is analytic for Re s>0 and gives
|E-E_K|<=2^(-K) Gamma(sigma)/|Gamma(sigma+iota tau)|.
For sigma>=1/4 and |tau|<=34, the Gamma product yields the explicit
rational bound

\[
 \frac{\Gamma(\sigma)}{|\Gamma(\sigma+\iota\tau)|}
 \le B:=\prod_{k=0}^{127}\left(1+\frac{136}{4k+1}\right)
       4^{\lceil c\rceil},\qquad
 c=\frac{34^2}{2}\left((128+1/4)^{-2}+(128+1/4)^{-1}\right).
                                                               \tag{48.12}
\]

Indeed each initial square-root product factor is at most
1+|tau|/(sigma+k). The logarithm of the remaining product is at
most tau^2/2 times the displayed integral-test bound, and exp(c)<=
4^ceil(c). For K=320 the full remainder T=2^(-320)B is less than
10^(-56), approximately 2.222e-57. On every outer disk of radius
delta=1/4 used here, Cauchy's estimate also gives |E'-E_K'|<=4T
at its centre.

On those outer disks, 1/4<=Re s<=3/4 and |Im s|<34. Each finite
Euler term has absolute value at most 1/2, while
|1-2^(1-s)|>=2^(1/4)-1>1/6, since (7/6)^4<2. Therefore
|zeta(s)|<=6(160+T)<1000=:M throughout the disk. None contains the
zeta pole at 1.

**Theorem.** At s_0=1/2+iota g, suppose outward readings give
|zeta(s_0)|<=U and |zeta'(s_0)|>=d. For 0<r<delta, if

\[
 \boxed{dr-U-\frac{Mr^2}{\delta(\delta-r)}>0,}           \tag{48.13}
\]

the disk |s-s_0|<r contains exactly one simple zeta zero, and that
zero has real part exactly 1/2.

**Proof.** Cauchy's coefficient bound puts the Taylor remainder beyond
the linear term below M r^2/(delta(delta-r)) on the circle. Inequality
(48.13) makes the homotopy from zeta to zeta'(s_0)(s-s_0) nonzero
on that circle. The argument principle, equivalently Rouche's theorem,
gives one zero counted with multiplicity. The functional equation and
conjugation symmetry send each zero rho to 1-conjugate(rho). This disk
is invariant under that reflection. Uniqueness forces the zero to be
fixed by it, proving Re rho=1/2 and simplicity. QED.

The phase was not forced to zero in the computation. Both complex
components of zeta and its derivative were enclosed. The centre lies
on the critical line, but the zero's exact real part follows from
existence, uniqueness and reflection, not from that choice alone.

Complex exponentials in (48.11) are computed by outward rectangle
arithmetic. Divide the argument by a power of two until its L1 bound
is <=1/2, sum through degree 96, add the rigorous tail
(1/2)^97/(97! (1-(1/2)/98)), then square with outward rounding.
No floating-point complex special function is used as evidence.

## 5. Executed results and the finite/global boundary

For N=5,8,12,16, both acquired L_N and actual J_N are interval-positive.
Their least L pivots are approximately 2.817e-19, 8.337e-28,
8.979e-38 and 4.396e-47. These are coordinate pivots, not eigenvalue
gap estimates. Signed interval LDL gives disjoint ordinate brackets
of radius 10^(-9) for all 5, 8 and 12 finite poles respectively, and
the lowest five poles of the N=16 pencil: 30 brackets total.

The narrow sixth N=16 bracket is inconclusive at the inherited
precision. Its entrywise interval pivot crosses zero. This does not
negate the proved positivity of L_16 and J_16 or turn a midpoint into
a certified pole.

The N=16 centres give the following independent actual-zeta certificates.
Each disk has radius 10^(-6); the margin is the left side of (48.13).

| Centre ordinate g | Proved margin greater than | Actual conclusion in this disk |
| --- | --- | --- |
| 14.134725141735 | 7.67e-7 | One simple zero, Re rho=1/2 |
| 21.022039638772 | 1.09e-6 | One simple zero, Re rho=1/2 |
| 25.010857580146 | 1.27e-6 | One simple zero, Re rho=1/2 |
| 30.424876126791 | 1.10e-6 | One simple zero, Re rho=1/2 |
| 32.935061722570 | 8.86e-7 | One simple zero, Re rho=1/2 |

This is a certified local source-to-zero pipeline. No count of all
zeros below these heights is performed, and no new zero-discovery or
priority claim is made. The pole's 10^(-9) bracket is not substituted
for the actual zero's 10^(-6) disk.

A supplied exact countermodel explains why this does not complete RH.
Take

\[
 \mathcal X_{\rm bad}(t)=
   \prod_{j=1}^5(t+j^2)^{10}\big((t-x)^2+y^2\big),\quad
 x=1/16-10000,\quad y=50.                              \tag{48.14}
\]

Its pullback Xcal_bad((s-1/2)^2) obeys both reflection symmetries,
yet has the off-line quartet 1/2 plus/minus (1/4 plus/minus 100 iota).
On the five nodes 3,4,5,6,7, its exact L and J are strictly positive.
The finite rational interpolant therefore has only negative poles,
while its full generating function has nonreal poles in t. This is a
polynomial control model, not a counterexample to actual RH. It proves
that finite positivity plus the two symmetries does not justify the
universal inference.

The [replay](../../../proof_lab/diagonal_spectral_pencil.py) records
36 controls, including 8 negative controls; [readings](../../../proof_lab/DIAGONAL_SPECTRAL_PENCIL_RESULT.json),
[scope](SE48_SCOPE.md), [sources](SE48_SOURCE_RECORD.md) and
[pins](SE48_SOURCE_PINS.json) retain the different evidence layers.
The global source inequality still needs proof. Eta remains 1079/2000
with floor 1/492000, and eta-zero source positivity remains on V_16.
The five local zero certificates do not change that global bound.
