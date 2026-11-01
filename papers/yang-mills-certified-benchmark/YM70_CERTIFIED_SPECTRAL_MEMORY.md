# YM-70: four-rate Wilson memory and certified spectral enclosures

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This continuation returns to YM63's scale obstruction and YM64's complete
energy-dependent memory. It supplies two distinct results:

1. An explicit four-rate formula for the **order-theta-squared coefficient**
   of the actual square Wilson loop's adapted memory, at every integer
   block factor r>=3.
2. Positive rational lower and upper bounds for the **full interacting**
   memory, using its inverse moments. A quantitative error bound holds
   at every positive Laplace parameter, including outside the Taylor disk.

The first is a computation in the declared gauge model. The second is a
theorem applicable to that model; numerical interacting moments have not
been computed here. Neither extends the bare-coupling window or constructs
the spatial continuum. Exact finite controls accompany the written proofs;
they are not proof-assistant verification, peer review or physical data.

## Model and notation

Retain YM62's native SU(2) links on a finite periodic three-dimensional
cubic lattice, with continuous heat time, electric fundamental rate 3/4,
and Hamiltonian A=L_E-theta V, V=sum_p W_p. Here
W_p=Re(product of the oriented plaquette links), half the fundamental
trace. The word "unnormalized" below means not rescaled to unit L2 norm;
the trace convention remains half-trace.

For |theta|<1/8000, use the actual vacuum h and dmu=h^2 dPhi from YM64.
Its product-link cut retains the vacuum, has diagonal form C_r, hidden
operator D_r>=gamma(theta)>0, and energy-form map T_r. At fixed finite
volume these objects have the domains established in YM64-T1/T2. Set

\[
 X=D_r^{-1},\qquad 0<X\le\gamma^{-1}I,
 \qquad M_n[f,g]=\langle T_rf,X^nT_rg\rangle .             \tag{1}
\]

The inequality 0<X means injective and positive; it does not assert a
strict positive lower bound for X. An unbounded D_r is allowed. The
memory and exact reduced resolvent are

\[
 \Sigma_r(s)[f,g]=\langle T_rf,(I+sX)^{-1}T_rg\rangle,
 \quad F_r(s)=C_r-\Sigma_r(s),\quad
 (F_r(s)+s)^{-1}=I_r^*(\mathcal L+s)^{-1}I_r .             \tag{2}
\]

Throughout, s>=0 is a real Laplace parameter. No claim about uniform
analytic continuation to the real-time or complex-frequency boundary is
made. A finite-volume coefficient remainder below is not a bound uniform
in volume, r or a continuum trajectory.

Positive moment approximation and Schur elimination are established
tools: see [Golub and Meurant, *Matrices, Moments and Quadrature II*](https://www.gerard-meurant.fr/cgestbit_1997.pdf)
for inverse quadratic forms and quadrature bounds, and
[Griesemer and Hasler, *On the Smooth Feshbach--Schur Map*](https://arxiv.org/abs/0704.3244).
We give the specific domain, matrix, error and scaling arguments needed
here. We make no novelty claim for the general Galerkin, Stieltjes or
Chebyshev methods. The response properties in NSR1/NSR2 do not identify
their thermodynamic lambda coefficients with the moments (1).

## YM70-T1: the actual loop's four leading rates

Let f=W_C be one coarse elementary Wilson plaquette, pulled back to its
r-by-r fine square boundary, of length 4r. Assume r>=3 and coarse side
at least four. At theta=0, Phi(f)=0 and Phi(f^2)=1/4. Define
R_f=Gamma(V,I_r f). There are exactly

\[
 n_1=16r-8,\qquad n_2=4                              \tag{3}
\]

elementary plaquettes sharing respectively one link or a contiguous
two-link path with this boundary; all other plaquettes give zero.
For k=1,2 put

\[
 \lambda_0(k)=3r+3-\tfrac32k,\qquad
 \lambda_1(k)=3r+3+\tfrac12k .                         \tag{4}
\]

At each fixed volume, r and s>=0,

\[
 \boxed{\Sigma_r(s)[f,f]=\theta^2\sigma_r(s)
                       +O_{\Lambda,r,s}(\theta^3),}
\]
\[
 \boxed{\sigma_r(s)=\sum_{k=1}^2 n_k k^2
   \left\{\frac{1}{256(s+\lambda_0(k))}
          +\frac{1}{768(s+\lambda_1(k))}\right\}.}       \tag{5}
\]

Correspondingly the order-theta-squared coefficient of M_n[f,f] is

\[
 m_{r,n}=\sum_{k=1}^2 n_k k^2
  \left\{\frac{1}{256\lambda_0(k)^{n+1}}
         +\frac{1}{768\lambda_1(k)^{n+1}}\right\}.       \tag{6}
\]

**Proof.** Each fine boundary link meets four elementary plaquettes.
Only the four interior corner plaquettes meet the boundary twice.
Thus n_1+2n_2=16r and n_2=4, proving (3).

Fix such a plaquette. By orienting and cyclically starting both loops
along their common path, write W_p=Re(UV), f=Re(UW). U is the product
of the k common links, V the remaining plaquette path and W the remaining
large-loop path. The three paths are disjoint. Under product Haar
integration their products are independent native unit quaternions.
For P=Re(UV)Re(UW), its singlet and triplet parts on the shared path are

\[
 P_0=\tfrac14\operatorname{Re}(VW^{-1}),\qquad
 P_1=P-P_0 .                                           \tag{7}
\]

The native second moments Phi(U_i U_j)=delta_ij/4 give (7) directly.
Integrating first in V and W, or first in U, gives

\[
 \|P\|^2=\tfrac1{16},\quad
 \|P_0\|^2=\tfrac1{64},\quad
 \|P_1\|^2=\tfrac3{64},\quad
 \langle P_0,P_1\rangle=0 .                            \tag{8}
\]

On each shared link the singlet has Casimir zero and the triplet has
Casimir two. Every nonshared link has fundamental Casimir 3/4. There
are 4r+4-2k nonshared links, so L_E P_j=lambda_j(k)P_j. The product
identity L_E(ab)=a L_E b+b L_E a-2Gamma(a,b), with electric rates 3
and 3r for the two loops, now gives

\[
 R_p:=\Gamma(W_p,f)=\tfrac{k}{4}(3P_0-P_1).              \tag{9}
\]

These identities can also be obtained entirely by native quaternion
derivatives. The executable checks both k values in that algebra,
including the norms and the shared/nonshared Casimir equations.

Cross terms for distinct p,p' vanish, including with any free inverse
resolvent inserted. Indeed p has at least two links outside C, whereas
two distinct elementary plaquettes share at most one link. Some link
of p outside C is absent from p'. Both components (7) are fundamental
on that link and have zero Haar mean there. The free resolvent preserves
these Casimir components. This is orthogonality from a unique external
link; it does **not** require the plaquette supports to be disjoint.

YM64-T5 puts R_f in the free hidden space and gives the leakage
Q_theta L_theta I_r f=-(2theta/3)R_f+O(theta^2). Thus the leading
self-energy is (4theta^2/9)<R_f,(D_0+s)^(-1)R_f>. Combining (8),
(9) and the vanished cross terms yields (5) and (6). Analytic
perturbation is only used at fixed volume and fixed readout. Subtracting
the actual mean from f leaves the memory unchanged because T_r1=0. QED.

For r=3 the answer is particularly simple:

| Hidden rate lambda | Residue in sigma_3(s) |
| --- | --- |
| 9 | 1/16 |
| 21/2 | 5/32 |
| 25/2 | 5/96 |
| 13 | 1/48 |

Consequently

\[
 \sigma_3(0)=\frac{113}{4095},\qquad
 \lim_{s\to\infty}s\sigma_3(s)=\frac7{24},\qquad
 \|R_f\|^2=\frac{21}{32}.                              \tag{10}
\]

Replacing f by its unit-norm free version 2f multiplies the residues
and moments by four. These are four relaxation rates of the **leading
memory coefficient**, not four measured particles or the full
interacting spectrum. At nonzero theta, higher-order contributions
need not have only four rates. In particular, a matched scalar coarse
plaquette coupling misses a nonzero memory already at order theta^2.

## YM70-T2: finite moment pencils bound the full response

Take any finite family f_1,...,f_m in the retained form domain, and
let Z:c->sum_alpha c_alpha T_r f_alpha. In this section M_n denotes
the m-by-m matrix Z*X^nZ and Sigma(s)=Z*(I+sX)^(-1)Z. The readouts
need not be independent or orthonormal. Matrix inequalities are in
Loewner order and allow noncommuting M_n.

For N>=0 define the two block pencils, indexed by i,j=0,...,N,

\[
 K_N(s)=(M_{i+j}+sM_{i+j+1})_{ij},\qquad
 b_N=(M_0,M_1,\ldots,M_N)^\mathsf T,                    \tag{11}
\]
\[
 K'_N(s)=(M_{i+j+1}+sM_{i+j+2})_{ij},\qquad
 b'_N=(M_1,M_2,\ldots,M_{N+1})^\mathsf T.               \tag{12}
\]

Use the inverse on the common Gram-kernel quotient, equivalently the
Moore--Penrose inverse, and put

\[
 L_N(s)=b_N^*K_N(s)^\dagger b_N,\qquad
 U_N(s)=M_0-s\,b_N'^*K'_N(s)^\dagger b'_N .              \tag{13}
\]

Then for every s>=0,

\[
 \boxed{0\le L_N\le L_{N+1}\le\Sigma
                   \le U_{N+1}\le U_N\le M_0.}         \tag{14}
\]

Only M_0 through M_(2N+2) are needed for the pair at level N. Both
bounds equal M_0 at s=0. Gram singularity, including an exactly retained
vacuum, is allowed and must not be repaired by arbitrary regularization.

**Proof.** Put A_s=I+sX and let V_N map block vectors to
sum_(i=0)^N X^i Zc_i. Then K_N=V_N* A_s V_N and b_N=V_N*Z. Since
I<=A_s<=(1+s/gamma)I, ker K_N(s)=ker V_N, independent of s>=0,
and b_N takes values in the range of K_N. On a finite Gram quotient,
the variational identity

\[
 \langle z,A_s^{-1}z\rangle
 =\sup_y\{2\operatorname{Re}\langle z,y\rangle
                         -\langle y,A_sy\rangle\}     \tag{15}
\]

restricted to ran V_N gives L_N<=Sigma. The trial spaces increase
with N. Apply the same argument to Z'=X^(1/2)Z, which is bounded
because X is bounded. Its moments are M_(n+1). The exact identity

\[
 M_0-\Sigma(s)=s Z'^*(I+sX)^{-1}Z'                     \tag{16}
\]

gives the upper bound and its monotonicity. At s=0 the vector being
approximated already lies in the trial space. This proves (14) for
every coefficient vector c, hence as a matrix inequality. Nothing
requires simultaneous diagonalization of the M_n. QED.

The lower pencil is a positive rational Stieltjes response: in an
orthonormal basis of ran V_N it is Z_N*(I+s X_N)^(-1)Z_N, with
0<X_N<=gamma^(-1)I. Its pole residues are positive and its poles lie
at s<=-gamma. The upper bound is also a positive rational response
plus a positive constant term. To see this, diagonalize the finite
compression of X in the primed trial space and use
s/(1+s x)=x^(-1)-x^(-1)/(1+s x). Its constant term is the s->infinity
limit of U_N, positive by (14). This finite-dimensional explanation
of the pole signs is not a discretization of D_r.

For T1's scalar leading model there are four distinct atoms x=1/lambda.
The vectors 1,X,X^2,X^3 span its cyclic space by the Vandermonde
identity. Hence L_3(s)=sigma_r(s)=U_3(s) for that coefficient model,
at every s>=0. This termination does not extend to the full interacting
memory without an additional finite-rank theorem.

## YM70-T3: a quantitative bound beyond the Taylor disk

Let T_j denote the Chebyshev polynomial T_0(t)=1, T_1(t)=t,
T_(j+1)(t)=2tT_j(t)-T_(j-1)(t). For s>0 define

\[
 \epsilon_N(s;\gamma)
   =\frac{1}{T_{N+1}(1+2\gamma/s)^2},\qquad
 \epsilon_N(0;\gamma)=0 .                             \tag{17}
\]

Then

\[
 0\le\Sigma-L_N\le\epsilon_N\Sigma,\qquad
 0\le U_N-\Sigma\le\epsilon_N(M_0-\Sigma),              \tag{18}
\]
\[
 \boxed{0\le U_N-L_N\le\epsilon_N M_0\le\epsilon_N C.} \tag{19}
\]

Here C is the diagonal energy matrix of the same readouts. The bound
is valid for every positive s/gamma, with no convergence-disk condition.
For 0<=s<=S one may replace epsilon_N(s;gamma) by epsilon_N(S;gamma).

**Proof.** The spectrum of A_s lies in [1,b], b=1+s/gamma. The
degree-N+1 residual polynomial

\[
 p_{N+1}(t)=
 \frac{T_{N+1}((b+1-2t)/(b-1))}
      {T_{N+1}((b+1)/(b-1))}                           \tag{20}
\]

has p_(N+1)(0)=1 and sup_[1,b]|p_(N+1)|^2<=epsilon_N. Therefore
q_N(t)=(1-p_(N+1)(t))/t is a degree-N polynomial. For z=Zc,
q_N(A_s)z lies in ran V_N. Galerkin minimization of the squared
A_s-energy error, followed by spectral calculus, gives

\[
 \langle z,A_s^{-1}z\rangle-c^*L_Nc
 \le\|p_{N+1}(A_s)A_s^{-1}z\|_{A_s}^2
 \le\epsilon_N\langle z,A_s^{-1}z\rangle .             \tag{21}
\]

Repeat with z=Z'c and use (16) for the second inequality (18).
Adding the two errors proves (19), with **one** epsilon_N M_0, not
two. YM64 supplies M_0<=C. For t>=1, T_j(t) is positive and increasing,
so epsilon_N increases with s. The s=0 identities were already proved.
QED.

Equivalently, if q=(sqrt(1+s/gamma)-1)/(sqrt(1+s/gamma)+1), then
epsilon_N=4q^(2N+2)/(1+q^(2N+2))^2. Thus the width tends to zero on
every bounded positive spectral interval. Formula (17) is preferable
for the exact certificate because rational s and gamma give rational
bounds with no floating-point square root. If the actual gap is known
only through a certified positive rational lower bound, substituting
that lower bound weakens (19) safely.

This is an approximation theorem from **exact moments**. An algorithm
using noisy or interval-valued moments needs an additional perturbation
analysis, particularly near a Gram-rank change. Necessary positive
moment checks alone do not reconstruct the interacting vacuum or
certify an arbitrary proposed infinite moment sequence.

## YM70-T4: physical-clock covariance and a conditional scale criterion

Rescale the full generator by kappa>0. The cut is unchanged, and

\[
 D^{(\kappa)}=\kappa D,\quad
 T^{(\kappa)}=\sqrt\kappa T,\quad
 M_n^{(\kappa)}=\kappa^{1-n}M_n .                       \tag{22}
\]

Diagonal rescaling of the Krylov columns in (11)--(12) proves

\[
 L_N^{(\kappa)}(s)=\kappa L_N(s/\kappa),\qquad
 U_N^{(\kappa)}(s)=\kappa U_N(s/\kappa).                 \tag{23}
\]

The physical-band error bound is therefore

\[
 0\le U_N^{(\kappa)}(s)-L_N^{(\kappa)}(s)
 \le\frac{C^{(\kappa)}}
 {T_{N+1}(1+2\kappa\gamma/S)^2},\quad 0\le s\le S.     \tag{24}
\]

For a proposed sequence of scales, **if** kappa gamma>=g_*>0 and
the renormalized readouts have C^(kappa)<=E_*I uniformly, one common
N controls the memory width on [0,S] by
E_*/T_(N+1)(1+2g_*/S)^2. This follows directly from (24), independent
of the hidden dimension. It is a conditional criterion for preserving
the spectral family across scales, not a proof of those two hypotheses.

The known model still has theta=c_B/(c_E g^4) and
kappa=c_E g^2/a. A trajectory with g(a)->0 exits |theta|<1/8000.
This stage supplies no estimate outside that window, no continuum
flow of the complete interaction and no bounded renormalized-energy
family meeting the criterion. Keeping a full memory approximation is
one necessary piece of that work, not a replacement for it.

## YM70-T5: enclosures survive further exact elimination

There is a distinction between sampling a response and defining a
complete retained system. T2/T3 hold on any finite test family in the
actual YM64 cut. **Inverting that sampled matrix does not in general
give the compressed full resolvent:** other retained channels may
still couple to it.

For an exact finite retained system, instead take the orthogonal
projection onto the vacuum and a finite orthonormal family of smooth
centered gauge-invariant readouts. Place **all other modes** in the
hidden space. This projection preserves the form domain, its hidden
operator obeys D>=gamma, and the YM64 energy-form construction and
T2/T3 apply afresh. The moments for this cut need not equal the sampled
moments of the original spatial cut. Here C and F are genuinely the
complete finite retained matrices, and (2) is exact.

Write F_-=C-U_N and F_+=C-L_N. At s>0,

\[
 (F_++s)^{-1}\le(F+s)^{-1}\le(F_-+s)^{-1}.              \tag{25}
\]

The inverses exist: F_->=C-M_0, which is at least gamma on centered
retained vectors and zero on the vacuum. This bound follows by
minimizing the full energy over hidden vectors while keeping the
retained vector fixed; their norms add orthogonally. Adding s makes
the full retained matrix strictly positive.

For any further isometry J into that retained system, compress (25)
and reverse order on inversion. The resulting next effective matrix
obeys

\[
 \boxed{
 [J^*(F_-+s)^{-1}J]^{-1}-s
 \quad\le F_{\rm next}(s)\le\quad
 [J^*(F_++s)^{-1}J]^{-1}-s .}                           \tag{26}
\]

This is an order-preserving enclosure at the **same** spectral
parameter, with no static reset. It iterates for any finite chain of
compatible retained systems. We do not claim the numerical width in
(19) remains unchanged through these inversions. A direct final-cut
application has its own bound (19); a propagated enclosure can widen.
The certificate verifies (26) and gives an exact counterexample to
inverting only a sampled diagonal compression.

## Evidence and the next unresolved development

`certificates/ym70_certified_spectral_memory.py` checks native singlet/
triplet decomposition and exact Haar norms, actual lattice incidence,
four-pole reconstruction, singular and noncommuting moment pencils,
Chebyshev error widths, clock covariance, and nested order. Finite
matrix examples are labelled algebraic fixtures; no truncated fixture
is substituted for the interacting Yang--Mills vacuum. Source pins,
tamper controls, strict rational inputs and Python 3.12 are enforced.

The next substantial bridge is to **compute or enclose the interacting
moments with spatial and scale-uniform errors**, while extending control
of the generated interactions toward a continuum trajectory. This
stage solves neither full spatial Schur-memory locality nor the weak
bare-coupling limit. The 4D continuum existence/mass-gap problem remains
open in this programme. The concrete gain is an actual leading Wilson
memory spectrum and a complete positive-energy approximation theorem
that no longer relies on truncating a small-s Taylor series.
