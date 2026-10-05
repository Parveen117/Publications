# YM-72: growing blocks, matrix memory and the visible relaxation law

Monty Dabas. 5 October 2026. **Python 3.12 only.**

YM71 controlled one fixed Wilson-loop readout uniformly in lattice
volume. Here both the block size and the number of observed loops may
grow. The interacting memory is bounded as a matrix, with the single
dimensionless budget

\[
 \boxed{\rho=|\theta|\,rN.}                              \tag{1}
\]

For N square-loop readouts with disjoint one-plaquette neighborhoods,
r>=3, and 0<rho<=10^(-8), the actual matrix memory obeys

\[
 \boxed{\frac{99}{100}\theta^2\sigma_r(s)I_N
 \le\boldsymbol\Sigma_{r,N,\theta}(s)
 \le\frac{51}{50}\theta^2\sigma_r(s)I_N,
 \qquad s\ge0.}                                        \tag{2}
\]

This holds in every admissible finite volume, with **no upper bound
on r or N**. The coupling must decrease as their product grows; (2)
does not assert uniform accuracy for fixed nonzero theta and r->infinity.

Under rho->0 and r->infinity the rescaled memory has a universal
one-pole profile. A separate hierarchy controls each higher inverse
moment. The selected loop correlations approach a rate-three law on
the block clock, although the full gap on that clock tends to zero.
This last distinction is essential: a persistent observed rate is
not by itself a positive gap for the complete theory.

## Model and the family of readouts

On finite periodic cubic lattices, keep the actual YM62 Hamiltonian and the YM64 vacuum-adapted spatial
cut P_r, with Q_r=1-P_r, hidden generator D_r and vacuum representation
L. No vacuum or link-representation truncation is made. The native
fundamental electric rate is 3/4 and |theta|<1/8000 remains the inherited
coupling window. Write eta=|theta|.

Let f_1,...,f_N be congruent half-trace Wilson loops, each the pullback
of one coarse plaquette to an r-by-r fine square. Let S_i contain its
boundary and every elementary plaquette meeting that boundary. Assume
the fine **link sets** S_i are pairwise disjoint. This is a geometric
separation condition, not an assumption of independence in the
interacting vacuum.

Such families exist for arbitrary r,N: take xy-plane squares with
fine origins (3ri,0,0), i=0,...,N-1, on a periodic fine cube of side
(3N+1)r. Each star lies in a coordinate box expanded by one link from
its square. Adjacent origins differ by 3r and the wraparound gap is
larger, so their stars are disjoint for r>=3. The corresponding coarse
side is 3N+1>=4. A larger admissible volume also works.

Set q_i=Q_r L f_i and define the readout matrices

\[
 \boldsymbol\Sigma(s)_{ij}
   =\langle q_i,(D_r+s)^{-1}q_j\rangle_\mu,
 \qquad \mathbf M_n{}_{ij}=\langle q_i,D_r^{-(n+1)}q_j\rangle_\mu.
                                                               \tag{3}
\]

Matrix order refers to the Euclidean coefficient space of these
specified readouts. The f_i need not be orthogonal in the interacting
vacuum. These matrices sample the full spatial-cut forms. Their
inverses are not automatically compressed full resolvents.

The proof uses YM70's exact native four-rate coefficients and YM71's
hidden-lift estimates. Positive moment and residual arguments have
established precedents, including [Golub and Meurant](https://www.gerard-meurant.fr/cgestbit_1997.pdf)
and the [Feshbach--Schur framework](https://arxiv.org/abs/0704.3244).
The matrix extension, explicit growing-family budget and scope below
are proved here. General norm inequalities and approximate-eigenvector
arguments are not claimed as new methods.

## YM72-T1: a full interacting matrix enclosure

Take YM71's scalar constants at r and eta: G^2=(2r+1)/24,
lambda_-=3r, lambda_+=3r+4, Delta, E_*, d_*, and N_*=52r-16.
For N readouts replace them by

\[
 K_N=e^{4\pi BNN_*},\quad B=24\eta,
 \quad\Delta_N=\sqrt N\Delta,
 \quad E_N=\sqrt N E_*,\quad d_N=\sqrt N d_* ,           \tag{4}
\]
\[
 G_N=\sqrt{K_NG^2},\quad H_N=G_N+\Delta_N,
 \quad A_N=(K_N-1)G^2+2G_N\Delta_N+\Delta_N^2 .          \tag{5}
\]

The error constants

\[
 e_-={\lambda_+\over\lambda_-G^2}
       \left(A_N+2d_NH_N+{H_NE_N\over\lambda_-}\right),
\]
\[
 e_+=e_-+{\max(1,\lambda_+/\gamma)\over G^2}
                   \left(d_N+{E_N\over\lambda_-}\right)^2           \tag{6}
\]

give, for every s>=0,

\[
 \max(0,1-e_-)\theta^2\sigma_r(s)I_N
 \le\boldsymbol\Sigma(s)
 \le(1+e_+)\theta^2\sigma_r(s)I_N .                     \tag{7}
\]

**Proof.** For a coefficient vector c, combine the YM71 forcing
polynomials and their exact interacting hidden lifts linearly over
the readouts. The free forcing and its free resolvent have squared
norm G^2||c||^2 and quadratic response sigma_r(s)||c||^2. Their
cross terms vanish because disjoint stars are independent under
product Haar and each forcing has zero reference mean. Interacting
cross terms are not set to zero.

The triangle inequality and sum_i|c_i|<=sqrt(N)||c|| give the three
residual bounds in (4). The union of free polynomial supports has
at most NN_* links; the local density comparison gives K_N. The
lifted interacting vectors themselves may be nonlocal, just as in
YM71, and are estimated in full L2(mu). Apply the scalar completion
of the energy square to this arbitrary linear combination. Its
constants are (5)--(6), proving (7) in Loewner order. QED.

For theta nonzero, put S(s)=boldsymbolSigma(s)/(theta^2 sigma_r(s)).
Because e_+>=e_-,

\[
 \|S(s)-I_N\|\le e_+,\qquad
 \|[S(s),S(t)]\|\le2e_+^2,\quad s,t\ge0 .              \tag{8}
\]

The second estimate follows by subtracting the identity from both
factors before taking the commutator. No simultaneous diagonalization
of the interacting response matrices is presumed. Off-diagonal entries
are bounded by e_+ theta^2 sigma_r(s). This is not a distance-decay
theorem; separation is used only to obtain free orthogonality.

## YM72-T2: one explicit budget for arbitrarily large r and N

On |theta|<=10^(-6), the inherited gap formula gives
gamma>=gamma_*=1/125. Its endpoint is checked by directed rational
logarithm intervals. Hence this gap floor applies whenever
rho<=3*10^(-6), since rN>=3.

Define the following quantities depending only on rho:

\[
 \begin{split}
 K&=e^{4992\pi\rho},\qquad t_0=\sqrt{24/7},\qquad
 d_0={512\over3}\rho,\\
 d_1&=d_0+{768\rho\over\gamma_*},\\
 e_0&={16\over9}
       \left[(240+8/9+8\sqrt{80})\rho+7680\rho^2\right],\\
 d&=t_0d_0,\quad b=t_0d_1,\quad e=t_0e_0,
 \qquad h_0=\sqrt K+d .
 \end{split}                                           \tag{9}
\]
\[
 \begin{split}
 \epsilon_-(\rho)&={13\over9}
       [K-1+2\sqrt Kd+d^2+(2b+e)h_0],\\
 \epsilon_+(\rho)&=\epsilon_-(\rho)
                +{312\over7\gamma_*}(d_1+e_0)^2 .
 \end{split}                                           \tag{10}
\]

Then (7)--(8) hold with e_-,e_+ replaced by epsilon_-,epsilon_+,
for every integer r>=3 and N>=1 with eta rN=rho<=3*10^(-6).

**Proof.** The inequalities
NN_*<=52rN, G^2>=7/24, lambda_+/lambda_-<=13/9, and
lambda_+/G^2<=312/7 give the exponential and scalar factors above.
Moreover

\[
 \Delta_N\le d_0,\qquad d_N\le d_1,\qquad
 {E_N\over\lambda_-}\le e_0 .                          \tag{11}
\]

For the last bound, YM71 has ell=4r+3<=5r and
E_*=(16r/3)(2B ell+8eta/3+8B^2 ell+8sqrt(2eta B ell)).
After division by 3r and multiplication by sqrt(N), its four bracket
terms are at most 240rho, 8rho/9, 7680rho^2, and 8sqrt(80)rho.
Here eta=rho/(rN), r>=3 and N>=1 were used separately; no limiting
argument enters the finite inequality. Substitute (11) into (6) and
use gamma_* to obtain (9)--(10). QED.

All quantities in (9)--(10) increase with rho. At rho=10^(-8),
the directed certificate gives

\[
 \epsilon_-<0.005396<0.01,\qquad
 \epsilon_+<0.010608<0.02 .                             \tag{12}
\]

Thus (2) holds throughout that entire product-budget interval, not
merely at sampled block sizes. Also epsilon_+(rho)=O(rho) as rho->0.
The condition is sufficient; it is not claimed to be the largest
possible growing-block region.

## YM72-T3: the matrix memory has a one-pole scaling profile

Let z>=0 and rescale the spectral parameter by s=rz. The four native
rates are 3r plus offsets 0, 3/2, 7/2 and 4. Their residues sum to
(2r+1)/24. Therefore, for

\[
 \phi(z)={1\over12(z+3)},
\]

the exact free coefficient satisfies the uniform relative enclosure

\[
 {6r+3\over6r+8}\phi(z)
 \le\sigma_r(rz)\le\left(1+{1\over2r}\right)\phi(z).
                                                               \tag{13}
\]

For example, this follows by bounding every rate between 3r and 3r+4.
Combining (13) with T2 proves, for any sequence of admissible finite
volumes and disjoint readout families with

\[
 r\longrightarrow\infty,\qquad
 0<\rho=|\theta|rN\longrightarrow0,
\]

that

\[
 \boxed{\sup_{z\ge0}
   \left\|{\boldsymbol\Sigma(rz)\over\theta^2\phi(z)}-I_N\right\|
       \longrightarrow0 .}                            \tag{14}
\]

N may grow with r. Equation (14) compares matrices in their own
coefficient spaces and is a norm estimate; it does not assert a
fixed-Hilbert-space operator limit for an increasing family.

There is also an explicit first finite-block correction:

\[
 \sigma_r(rz)=\phi(z)+{\psi(z)\over r}+R_r(z),\qquad
 \psi(z)={z-1\over24(z+3)^2},
\]
\[
 { |R_r(z)|\over\phi(z)}\le {43\over36r^2}.             \tag{15}
\]

To verify it directly, put c=z+3. Write each residue as ra_j+b_j
and its rate as 3r+delta_j. In increasing-rate order the triples
(a_j,b_j,delta_j) are
(0,1/16,0), (1/16,-1/32,3/2), (1/48,-1/96,7/2),
and (0,1/48,4). Their sums give the displayed zeroth and first
orders. The exact remainder of one term is
delta_j(a_j delta_j-b_j c)/(r^2 c^2(c+delta_j/r)). Since
sum a_j delta_j^2=19/48 and sum |b_j|delta_j=1/6,
its total magnitude is at most
[19/(48c^3)+1/(6c^2)]/r^2, implying (15) for c>=3.

The full interacting matrix differs from the corrected profile by

\[
 \left\|{\boldsymbol\Sigma(rz)\over\theta^2\phi(z)}
       -\left(1+{z-1\over2r(z+3)}\right)I_N\right\|
 \le\epsilon_+(\rho)\left(1+{1\over2r}\right)+{43\over36r^2}.
                                                               \tag{16}
\]

If eta N r^2->0, then r times the right-hand error tends to zero.
This isolates the first correction in the actual interacting family.
It is a joint small-theta/large-block statement, not a claim at fixed
nonzero coupling.

## YM72-T4: each higher inverse order has its own block budget

Set Gamma_r(z)=boldsymbolSigma(rz)/theta^2 for theta nonzero.
It is a matrix Stieltjes response with hidden spectral support at
least gamma_*/r. Its nth inverse moment is

\[
 \mathbf A_{r,n}={r^n\over\theta^2}\mathbf M_n .         \tag{17}
\]

For h>0, the positive finite-difference identity from YM71 holds in
matrix order:

\[
 \mathbf F_{n,h}:={1\over n!h^n}
     \sum_{k=0}^n(-1)^k\binom nk\Gamma_r(kh),\qquad
 \mathbf F_{n,h}\le\mathbf A_{r,n}
 \le\prod_{k=1}^n\left(1+{krh\over\gamma_*}\right)\mathbf F_{n,h}.
                                                               \tag{18}
\]

The proof integrates the scalar partial-fraction identity against
the positive matrix spectral measure. The measured readouts and
moment matrices need not commute. Signed endpoint arithmetic in (18),
using T2, gives computable scalar lower and upper bounds multiplying
I_N for the actual interacting matrix moments. No uncertain Gram
matrix is inverted.

For each fixed n, a sufficient convergence condition is

\[
 \boxed{|\theta|N r^{n+1}\longrightarrow0
 \quad\Longrightarrow\quad
 \left\|\mathbf A_{r,n}-{I_N\over12\,3^{n+1}}\right\|
             \longrightarrow0,\qquad r\longrightarrow\infty.}     \tag{19}
\]

**Proof.** Let epsilon=epsilon_+(rho). The error in the finite
difference relative to its known free counterpart is at most
2^n epsilon sup_z sigma_r(rz)/(n!h^n). The supremum is bounded by
7/216 for r>=3. Choose h=(gamma_*/r)d with
d=(epsilon r^n)^(1/(n+1)). Under (19), epsilon r^n->0 and d->0.
The response-input error tends to zero, being bounded by a constant
depending on n and gamma_* times (epsilon r^n)^(1/(n+1)). The factor
on the right of (18) is product_(k=1)^n(1+kd), which tends to one.

The free scaled rates all lie in [3,3+4/r], with total residue
1/12+1/(24r). Its finite difference therefore converges to its
nth inverse moment, and that moment converges to 1/(12*3^(n+1)).
The lower and upper matrix enclosures squeeze A_(r,n) to this value.
For n=0 the result follows directly from (14) at z=0. QED.

The condition through order n controls the entire finite prefix of
the tower. Constants need not be uniform as n itself grows. In
particular, rho->0 alone does not justify differentiating (14).

The distinction is sharp for the available positive-measure information.
Let delta_r=gamma_*/r and w_r=r^(-n), n>=1. The positive scalar response

\[
 \widetilde\Gamma_r(z)
 ={1-w_r\over12(z+3)}+{w_r\delta_r\over36(z+\delta_r)}   \tag{20}
\]

lies between (1-w_r)phi(z) and phi(z) at every z>=0 and has the
required spectral floor. All inverse orders below n converge to
the profile values, but its nth moment has the additional limiting
term 1/(36 gamma_*^n). This is an algebraic counterexample to a
stronger inference from those bounds, not a counterexample claimed
for the actual Yang--Mills vacuum.

## YM72-T5: selected loop dynamics and the full gap separate

Use the block clock A_r=L/r, which fixes the free coarse-loop rate
at three. Let chi_i=f_i-mu(f_i), let G be their actual covariance
matrix, and let J:c->sum_i chi_i(G^(-1/2)c)_i. Then J is an isometry
into the gauge-invariant, vacuum-orthogonal space.

The union of loop boundaries has 4rN fine links. With
K_b=exp(384pi rho), the local-density comparison and the variational
definition of variance give

\[
 {1\over4K_b}I_N\le G\le{K_b\over4}I_N .               \tag{21}
\]

For the lower bound, minimize the weighted squared norm of
sum c_i f_i minus a constant, compare with product Haar, then
minimize there. Under Haar the centered loops have Gram matrix
I_N/4. No independence in the actual vacuum is used.

The exact generator equation L f_i=3r f_i-2Gamma(u,f_i) gives

\[
 (A_r-3)\chi_i=-{2\over r}
           [\Gamma(u,f_i)-\mu\Gamma(u,f_i)].
\]

Each loop has 4r links and each one-link gradient has norm at most
1/2. The score bound B=24eta, centering as an L2 contraction, and
(21) therefore imply

\[
 \boxed{\|(A_r-3)J\|\le R_{\rm obs}:=
                192\eta\sqrt{NK_b}.}                  \tag{22}
\]

Duhamel's identity on these smooth readouts and positivity of A_r give

\[
 \|e^{-tA_r}J-Je^{-3t}\|
 \le{R_{\rm obs}\over3}(1-e^{-3t}),
\]
\[
 \boxed{\sup_{t\ge0}
       \|J^*e^{-tA_r}J-e^{-3t}I_N\|
       \le64\eta\sqrt{NK_b}.}                         \tag{23}
\]

This is an actual interacting time-correlation statement. It is
derived directly from L and the isometry J, with no inversion of a
sampled spatial Schur matrix. In the product-budget region of (2),
the bound is below 2.134*10^(-7); as r grows with rho bounded it
tends to zero, since eta sqrt(N)=rho/(r sqrt(N)).

The spectral theorem also bounds unseen low-energy weight:

\[
 \|J^*1_{[0,3/2]}(A_r)J\|
             \le{4R_{\rm obs}^2\over9}.                \tag{24}
\]

Indeed |lambda-3|>=3/2 on that interval, so (22) bounds the norm
of its spectral projection applied to J. The vacuum contribution
is exactly zero because the readouts are centered.

Nevertheless YM63's physical-sector upper bound gives

\[
 0<\Delta_\Lambda(A_r)\le{100\over19r}
           \longrightarrow0 .                         \tag{25}
\]

Thus an observed rate tending to three coexists with a full gap
tending to zero on this particular clock. Slow modes have diminishing
weight in the selected large-loop observations. A finite matrix
counterexample in the certificate illustrates the distinction; (25)
itself is an inherited bound for the actual model.

## What the growing-block result does and does not construct

The four-rate memory profile in T3 has a nonzero limit only after
division by theta^2. On the block clock,
Sigma^(1/r)(z)=(theta^2/r)Gamma_r(z); its unnormalized amplitude
vanishes along the joint limit. The rate-three profile and correlation
law concern this chosen family of observations.

These results control arbitrarily large blocks and increasing finite
readout families, with explicit coupling-size budgets. They do not
prove a distance-decaying full Schur kernel, a uniform infinite-order
moment tower or a fixed-coupling large-block limit. The physical
dictionary theta=c_B/(c_E g^4) sends theta->0 toward large bare g,
whereas the usual weak-bare-coupling trajectory has g->0. The block
clock 1/r is a declared coarse-time normalization and has not been
identified with c_E g^2/a on a continuum trajectory. The 4D continuum
construction and its positive physical mass gap remain open.

The next scale problem is to control a family with persistent
interaction while retaining its slow modes. The present work supplies
both a growing-family estimate and an explicit test for when the
chosen observations lose sensitivity to those modes.

## Reproduction

The certificate checks native star geometry, matrix order and
noncommutation, exact one-pole asymptotics and finite-block correction,
directed universal budgets, matrix-moment stencils, and the separation
between visible relaxation and the full gap. It preserves its pinned
predecessors and uses Python 3.12 only. The numerical budgets follow
from the written proofs; finite fixtures do not compute the actual
interacting vacuum or constitute proof-assistant verification.
