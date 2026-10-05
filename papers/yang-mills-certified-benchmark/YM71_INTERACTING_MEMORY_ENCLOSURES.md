# YM-71: actual interacting memory with a volume-uniform remainder

Monty Dabas. 5 October 2026. **Python 3.12 only.**

YM70 computed a four-rate leading Wilson-memory coefficient but left its
perturbative remainder dependent on volume. Its full-response moment
pencils also needed moment inputs that had not been enclosed. This stage
supplies explicit two-sided intervals for the **actual interacting**
response and its scalar inverse moments on the same spatial cut.

The central construction uses the weighted one-link electric operator
to turn each free return component into an **exactly hidden interacting
trial vector**. The differentiated vacuum equation controls its equation
residual without a volume factor. Completing the energy square then
turns that residual into a response enclosure.

For a block factor of three, the resulting concrete consequence is

\[
 \boxed{\frac{49}{50}\theta^2\sigma_3(s)
 \le \Sigma_{3,\theta}(s)[f,f]
 \le \frac{107}{100}\theta^2\sigma_3(s)}                 \tag{1}
\]

for **every admissible finite lattice volume, every s>=0 and every
0<|theta|<=10^(-8)**. Here f is the half-trace coarse Wilson plaquette,
and sigma_3 is the four-rate function computed in YM70. These rates
describe a controlled approximation to the memory; they are not a
claim that the full interacting generator has only four eigenvalues.

This is a fixed-block, strong-bare-coupling result. It does not enlarge
YM62's coupling window, supply bounds uniform as r increases, construct
the spatial continuum or solve the four-dimensional mass-gap problem.
The calculations below certify constants in written proofs. They are
not proof-assistant verification, peer review or experimental evidence.

## Setup and the quantities being enclosed

Keep A=L_E-theta V_0, V_0=sum_p W_p, the actual positive vacuum h,
dmu=h^2 dPhi, and the vacuum-adapted product-link projection P=P_r of
YM64. Let Q=1-P and u=log h. In the vacuum representation,

\[
 \mathcal L=L_E-2\Gamma(u,\cdot),\qquad
 D=\mathcal L|_{Q}\text{ in the form sense},\qquad
 D\ge\gamma>0 .                                       \tag{2}
\]

Here gamma may be any certified positive lower bound on YM62's gap.
Write eta=|theta|<1/8000. YM67 gives, on every fine link e,
||nabla_e u||_infinity<=B, and we use the simple upper bound

\[
 B=24\eta .                                            \tag{3}
\]

All finite-volume functions below are smooth. The native derivatives
commute with L_E; no interchange of noncommuting ordered derivatives
is used. The gap is the inherited full-space gap, also applicable to
the derivative vectors that need not be gauge invariant.

For f=W_C on the r-by-r fine boundary, r>=3 and coarse side at least
four, put q=Q L f. The actual response and moments are

\[
 \Sigma(s)=\langle q,(D+s)^{-1}q\rangle_\mu,
 \qquad M_n=\langle q,D^{-(n+1)}q\rangle_\mu .           \tag{4}
\]

These are precisely YM64's Sigma_r(s)[f,f] and M_(r,n)[f,f]. Changing
f by its actual mean changes neither quantity. The cut and the full
hidden Hilbert space are retained; there is no finite spin truncation.

The residual-energy identity used below is an established variational
tool; see [Golub and Meurant, *Matrices, Moments and Quadrature II*](https://www.gerard-meurant.fr/cgestbit_1997.pdf).
YM67 provides the native vacuum derivative estimates, with its stated
[Feynman--Kac background](https://arxiv.org/abs/1811.03471).
The additional work here is the explicit hidden lift, its native
residual estimates and their application to the actual lattice memory.

## YM71-T1: an exact score equation and a uniform first-order error

Let w_(e,a)=D_(e,a)u. Differentiating Ah=Eh and then dividing by h gives

\[
 \mathcal L w_{e,a}=\theta D_{e,a}V_0,\qquad
 \mu(w_{e,a})=0 .                                      \tag{5}
\]

The second identity is Haar integration by parts:
int h D_(e,a)h dPhi=0. Integrating the first also gives
mu(D_(e,a)V_0)=0 when theta is nonzero. At theta=0 the score is zero.
On summing its Dirichlet energies over a=1,2,3,

\[
 \sum_a\mathcal E(w_{e,a})
 =\theta\mu\!\left(\sum_a w_{e,a}D_{e,a}V_0\right)
 \le 2\eta B .                                        \tag{6}
\]

Indeed at most four plaquettes contain e and each has one-link gradient
norm at most 1/2. The energy on the left sums derivatives over **all**
fine links. Thus (6) controls a possibly nonlocal score with no count
of the ambient lattice links.

Define the vector remainder on one link by

\[
 z_e=\nabla_e u-\frac{\theta}{3}\nabla_eV_0 .
\]

It is centered. Since L_E D_(e,a)V_0=3D_(e,a)V_0,

\[
 \mathcal L z_e=\frac{2\theta}{3}\Gamma(u,\nabla_eV_0),
 \qquad \boxed{\|z_e\|_{L^2(\mu;\mathbb R^3)}
                 \le\frac{8\eta B}{\gamma}.}           \tag{7}
\]

For the bound, every ordered second derivative of one plaquette is
at most 1/4 in absolute value. Its three-by-three derivative matrix
therefore has Frobenius norm at most 3/4. There are four incident
plaquettes and four differentiated link choices per plaquette. Hence
||Gamma(u,nabla_e V_0)||<=12B pointwise. Apply the centered inverse
of L, whose norm is at most 1/gamma, to obtain (7). This argument is
an exact equation estimate, not an assumed analytic perturbation series.

## YM71-T2: free return components have controlled exact hidden lifts

Use YM70's decomposition R_f=Gamma(V_0,f). For each incident plaquette
sharing k=1 or 2 boundary links, write its singlet/triplet functions as
P_0=(1/4)Re(VW^(-1)), P_1=W_p f-P_0. Decompose the free forcing
g=-(2/3)R_f into individual terms

\[
 g_{p,0}=-\frac{k}{2}P_0,\qquad
 g_{p,1}=\frac{k}{6}P_1 .                              \tag{8}
\]

Their electric rates are lambda_0(k)=3r+3-3k/2 and
lambda_1(k)=3r+3+k/2. Each p has a fine link e_p absent from **all**
retained product paths, not merely absent from C. Both components are
fundamental on that link, with L_(e_p)g_(p,j)=nu g_(p,j), nu=3/4.
Define the weighted electric operator and lift by

\[
 L_e^\mu=L_e-2\Gamma_e(u,\cdot),\qquad
 \widehat g_{p,j}=\nu^{-1}L_{e_p}^\mu g_{p,j}
 =g_{p,j}-\frac2\nu\sum_a w_{e_p,a}D_{e_p,a}g_{p,j} .   \tag{9}
\]

**Exact hiddenness.** For any retained smooth F, D_(e_p)F=0.
Weighted integration by parts gives <F,L_(e_p)^mu g>=0. By density,
P widehat g=0. The lift is gauge invariant because h, g and the
one-link Casimir form are gauge invariant. It is smooth at each finite
volume, belongs to Dom D, and D widehat g=Q L widehat g. This last
statement follows directly from the restricted form, not from assuming
an operator-domain block decomposition.

Here are explicit constants for all its errors. Each term is a sum of
trace words/products with total absolute scalar coefficient

\[
 A_{p,0}=k/8,\quad A_{p,1}=5k/24,\quad
 A_*:=\sum_{p,j}A_{p,j}=16r/3 .                         \tag{10}
\]

At most ell=4r+3 fine links occur in each term, with degree at most two
per link and degree one on e_p. Consequently
||nabla_l g||<=A_(p,j), ||nabla_(e_p)g||<=A_(p,j)/2, and each ordered
entry D_(l,b)D_(e_p,a)g has magnitude at most A_(p,j)/2.
These are trace-word insertion bounds: a native derivative contributes
one generator of norm 1/2 per occurrence.

It follows from (9) that

\[
 \sum_{p,j}\|\widehat g_{p,j}-g_{p,j}\|_\mu
 \le \Delta:=\frac43 B A_* .                           \tag{11}
\]

For a term with free rate lambda, use (5), the product rule for L,
and L_E D_(e,a)g=lambda D_(e,a)g. The exact equation residual before
the contractive projection Q is

\[
 \begin{split}
 (\mathcal L-\lambda)\widehat g={}&-2\Gamma(u,g)
 -\frac{2\theta}{\nu}\sum_a(D_{e,a}V_0)(D_{e,a}g)\\
 &+\frac4\nu\sum_a w_{e,a}\Gamma(u,D_{e,a}g)
 +\frac4\nu\sum_a\Gamma(w_{e,a},D_{e,a}g).
 \end{split}                                           \tag{12}
\]

The four terms have L2 bounds respectively
2B ell A, (8/3)eta A, 8B^2 ell A, and
8 sqrt(2eta B ell) A. For the last, apply Cauchy--Schwarz to (6)
and to the second-derivative matrices of g. Therefore

\[
 \sum_{p,j}\|(D-\lambda_{p,j})\widehat g_{p,j}\|_\mu
 \le E_*:=A_*\left(2B\ell+\frac83\eta+8B^2\ell
                           +8\sqrt{2\eta B\ell}\right).\tag{13}
\]

Finally let H=sum widehat g_(p,j). Since L_E f=3r f is retained,
q=-2Q Gamma(u,f). Apply (7), ||nabla_e f||<=1/2 on its m=4r links,
and QH=H to obtain

\[
 \boxed{\|q-\theta H\|_\mu\le\eta d_*,\qquad
 d_*:=\Delta+\frac{8mB}{\gamma}.}                       \tag{14}
\]

For example, the score-remainder part contributes at most
2 sum_(e in C)(1/2)||z_e||<=8m eta B/gamma. Equations (11)--(14)
have no volume factor. The **trial functions** involve the unknown
score, but their bounds do not require evaluating it; all constants
are known. Thus this is a computable enclosure, not a computed vacuum.

## YM71-T3: a uniform bound on the complete positive-parameter response

The union of the supports of the g_(p,j) contains at most
N_*=52r-16 fine links. This counts the 4r boundary links and at most
three or two new links for each one-link or two-link incident plaquette.
YM67's marginal-density argument gives

\[
 K_*^{-1}\le\frac{d\mu_{\rm support}}{d\Phi_{\rm support}}
       \le K_*,\qquad K_*:=e^{4\pi B N_*}.              \tag{15}
\]

This local comparison is used only for the free polynomials g and
their free resolvents. We do **not** assume the nonlocal lifted vectors
have that finite support. Their errors are controlled in full L2(mu)
by T2.

From YM70's exact free orthogonality, set

\[
 G^2=\|g\|_\Phi^2=\frac{2r+1}{24},\quad
 \lambda_-=3r,\quad\lambda_+=3r+4,\quad
 G_*=\sqrt{K_*G^2},\quad H_*=G_*+\Delta,
\]
\[
 A_0=(K_*-1)G^2+2G_*\Delta+\Delta^2 .                  \tag{16}
\]

Define explicitly

\[
 \varepsilon_-=
 \frac{\lambda_+}{\lambda_-G^2}
       \left(A_0+2d_*H_*+\frac{H_*E_*}{\lambda_-}\right),
\]
\[
 \varepsilon_+=\varepsilon_-+
 \frac{\max(1,\lambda_+/\gamma)}{G^2}
            \left(d_*+\frac{E_*}{\lambda_-}\right)^2 .  \tag{17}
\]

Then for every s>=0, at the full interacting theta,

\[
 \boxed{\max(0,1-\varepsilon_-)\theta^2\sigma_r(s)
  \le\Sigma(s)\le(1+\varepsilon_+)\theta^2\sigma_r(s).} \tag{18}
\]

**Proof.** Let v_s=sum widehat g_(p,j)/(lambda_(p,j)+s) and
u_s=sum g_(p,j)/(lambda_(p,j)+s). Then
||v_s||_mu<=H_*/(lambda_-+s), ||v_s-u_s||<=Delta/(lambda_-+s),
and ||(D+s)v_s-H||<=E_*/(lambda_-+s). With x=theta v_s and
e_s=q-(D+s)x, (14) implies

\[
 \|e_s\|\le\eta\left(d_*+\frac{E_*}{\lambda_-+s}\right).
\]

Completing the square gives the exact identity

\[
 \Sigma(s)=J_s+\langle e_s,(D+s)^{-1}e_s\rangle,
 \quad J_s=2\operatorname{Re}\langle q,x\rangle
                           -\langle x,(D+s)x\rangle . \tag{19}
\]

The last term is between zero and ||e_s||^2/(gamma+s).
Expanding J_s using q=theta H+(q-theta H), the norm estimates above,
and the marginal comparison (15) yields

\[
 |J_s-\theta^2\sigma_r(s)|
 \le\frac{\eta^2}{\lambda_-+s}
       \left(A_0+2d_*H_*+\frac{H_*E_*}{\lambda_-+s}\right).\tag{20}
\]

For clarity, the reference-to-mu term is at most
(K_*-1)||g||_Phi||u_s||_Phi. The two lift errors and their product
are bounded by (2G_*Delta+Delta^2)/(lambda_-+s). These are the terms
in A_0. No free Gram inner product is silently used for the lifted
interacting vectors.

Use sigma_r(s)>=G^2/(lambda_++s),
(lambda_++s)/(lambda_-+s)<=lambda_+/lambda_-, and
(lambda_++s)/(gamma+s)<=max(1,lambda_+/gamma). Equations (19)--(20)
give (17)--(18). QED.

At fixed r, gamma has a positive limit as eta decreases to zero,
K_*-1=O_r(eta), Delta,E_*,d_*=O_r(eta). Thus epsilon_-=O_r(eta)
and epsilon_+-epsilon_-=O_r(eta^2). In particular,

\[
 \Sigma_{r,\theta}(s)[f,f]
 =\theta^2\sigma_r(s)
   +O_r\!\left(\frac{|\theta|^3}{\lambda_-+s}\right)     \tag{21}
\]

uniformly over admissible finite volumes and s>=0 as theta->0. This
is the specific volume-uniform remainder that YM70 did not establish.
It is not an all-order expansion uniform in r or in the continuum.

## YM71-T4: certified inputs for the inverse-moment tower

There are two independent ways to get moments from T2/T3. First, for
j=n+1 define

\[
 T_j(\lambda_-,\gamma)=
  \sum_{k=1}^j\lambda_-^{-(j-k+1)}\gamma^{-k} .           \tag{22}
\]

This symbol denotes an inverse-power sum, not YM70's Chebyshev
polynomial. The exact identity for an approximate eigenvector is

\[
 D^{-j}\widehat g-\lambda^{-j}\widehat g
 =-\sum_{k=1}^j\lambda^{-(j-k+1)}D^{-k}(D-\lambda)\widehat g.
\]

Using (14), (16), and this identity gives the explicit moment radius

\[
 \boxed{|M_n-\theta^2m_{r,n}|\le R_n,\quad
 R_n=\eta^2\left[
  \frac{A_0}{\lambda_-^{n+1}}+H_*E_*T_{n+1}
    +(2H_*d_*+d_*^2)\gamma^{-(n+1)}\right].}            \tag{23}
\]

Indeed the change from q to theta H costs at most
eta^2(2H_*d_*+d_*^2)gamma^(-(n+1)); the inverse-power residual costs
eta^2 H_*E_*T_(n+1). The free-to-weighted and lift comparison costs
eta^2 A_0/lambda_-^(n+1), exactly as in T3. Here m_(r,n) is the
known positive coefficient in YM70 equation (6). For every fixed n
and r, (23) is O_(r,n)(|theta|^3) uniformly in volume. Its constants
can grow rapidly with n, so this is not one uniform bound for the
whole tower.

A second, often tighter method uses only the scalar intervals (18).
For a step h_s>0 and n>=0 define

\[
 F_{n,h_s}=\frac{1}{n!h_s^n}
     \sum_{k=0}^n(-1)^k\binom nk\Sigma(kh_s).
\]

The partial-fraction identity
sum_k (-1)^k binom(n,k)/(lambda+k h_s)
=n!h_s^n/product_(k=0)^n(lambda+k h_s) and the positive spectral
measure in (4) prove

\[
 \boxed{F_{n,h_s}\le M_n
       \le\prod_{k=1}^n(1+kh_s/\gamma)F_{n,h_s}.}        \tag{24}
\]

Evaluate the signed sum by interval arithmetic, using lower endpoints
for positive coefficients and upper endpoints for negative coefficients
to obtain its lower bound, and conversely for its upper bound. Clip
the lower bound at zero and apply the positive factor in (24).
Intersect these intervals over any declared steps with (23) and
0<=M_n<=gamma^(-n)Sigma(0). The result is a rigorous enclosure of
the **interacting** scalar moment, without evaluating the vacuum.

It is invalid to differentiate a response inequality as if it were an
equality. For example, with 0<delta,epsilon<1, the response
(1-delta)/(s+1)+delta epsilon/(s+epsilon) lies between
(1-delta)/(s+1) and 1/(s+1) for every s>=0, yet its first inverse
moment is 1-delta+delta/epsilon and can be arbitrarily large. Both
(23) and (24) use the missing spectral information explicitly.

These intervals are usable moment information, not permission to put
their midpoints into YM70's exact Gram inversion. A robust rational
pencil method with interval inputs still needs its own conditioning
analysis. The present interval method avoids that inversion.

## YM71-T5: a nonzero uniform window and physical-clock covariance

For r=3, certified upward bounds at eta=10^(-8) give

\[
 \varepsilon_-<0.016<1/50,\qquad
 \varepsilon_+<0.061<7/100.                             \tag{25}
\]

Every factor in (17) increases with eta when gamma is taken as
YM62's decreasing gap formula (or its endpoint lower bound on an
interval). Thus the same estimates hold for all 0<=eta<=10^(-8).
Equation (1) follows. Since sigma_3(s)>0 at finite s, the full adapted
cut has **strictly nonzero memory**, uniformly over all these volumes
and nonzero couplings. This is an actual-coupling result, not only a
formal nonzero coefficient.

The following outward decimal display describes ratios to the known
leading response theta^2 sigma_3(s) or moment theta^2 m_(3,n), as
appropriate. The exact rational endpoints are saved
in `YM71_RESULT.json`; different rows intentionally use different
couplings because high inverse orders amplify small-energy uncertainty.

| Absolute coupling eta | Quantity | Certified ratio interval |
| --- | --- | --- |
| 10^(-8) | Sigma(s), every s>=0 | [0.984, 1.061] |
| 10^(-12) | M_1 | [0.996, 1.133] |
| 10^(-16) | M_2 | [0.996, 1.226] |
| 10^(-20) | M_3 | [0.997, 1.029] |

The interval widths near the outer edge of the original coupling
window can be very large. Equation (25) is not a claim of seven-percent
accuracy throughout |theta|<1/8000.

Under a physical-clock change L->kappa L, kappa>0,

\[
 \Sigma^{(\kappa)}(s)=\kappa\Sigma(s/\kappa),\qquad
 M_n^{(\kappa)}=\kappa^{1-n}M_n .                       \tag{26}
\]

Scale the endpoints by the same factors. The relative enclosure (18)
then holds on the entire positive physical spectral axis with the
same epsilons. This does not provide a continuum coupling trajectory:
theta=c_B/(c_E g^4) still leaves the established window when g->0.
The constants N_*, Delta and E_* also depend on r. Controlling the
generated interactions as spatial blocks grow, and going beyond the
strong-bare-coupling region, remain separate requirements.

## Verification and next mathematical target

The executable checks exact native weighted integration by parts and
hiddenness, electric derivative identities, the actual unretained-link
geometry, residual-energy and inverse-power algebra, directed budgets,
signed Stieltjes stencils, and clock scaling. A weighted polynomial
fixture is explicitly distinguished from the actual Yang--Mills vacuum.
The intervals for the actual model follow from T1--T5 and the inherited
gap/score bounds, not from fitting a finite matrix to that vacuum.

This stage encloses one fixed square Wilson readout on the actual full
spatial cut. It does not establish decay of the full memory kernel
between separated regions or the existence of an infinite-volume
operator-valued memory limit. The next structural target is to extend
these residual estimates to coupled readouts while controlling their
spatial dependence and the growth with block size. The 4D continuum
construction and physical mass prediction remain open.
