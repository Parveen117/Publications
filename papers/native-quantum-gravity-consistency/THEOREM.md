# QGC-1 — Quantum consistency, outcome retention and geometric identification

14 October 2026.

## 1. Premises and relation to the existing programme

Use the existing finite native dagger algebra over the rational cut field
`Q + iota Q`, with `iota^2=-1`, its positive carrier and trace. Put

\[
 Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
 X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
 R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\quad Y=\iota R.
\]

Thus `Z X=-R`, `X Z=R`, and `Z,X,Y` are Hermitian involutions.
A positive source is

\[
 \rho=\tfrac12(NI+zZ+xX+yY),\qquad
 N\ge0,\qquad z^2+x^2+y^2\le N^2.                 \tag{1}
\]

For `N>0`, the normalized state is `sigma=rho/N`. Its rank-one cuts obey
`P^2=P`, `Tr P=1`. Zero amount has no normalized conditional state.

**Operational adapter.** Interpret a positive effect `E` as a reading
`Tr(E rho)`, and divide by `N` for its probability. Classical randomized
preparations with the same density operator are operationally equivalent
when the preparation label is not supplied to the operation. Independent
registers compose by the native tensor product. These are declared quantum
statistical premises; they are not derived from idempotence alone. The
matrix realization is already part of the native framework. No physical
Hilbert-space geometry is inferred from the fact that it has a finite
matrix representation.

EGR-1 gives the two-mode stationary response

\[
 \mathbf y=\chi K^{-1}\mathbf r,\quad
 K=\operatorname{diag}(10,47),\quad
 \chi=\alpha N/\gamma,\quad
 \mathbf r=(\operatorname{Tr}Z\sigma,\operatorname{Tr}X\sigma).
                                                               \tag{2}
\]

Its constants depend on the specified nine-cell Euclidean Regge geometry
and source readout. They are not universal gravitational constants.
EGR-2 constructs the local nonlinear stationary branch. CRC-1 shows that
a quadratic source coefficient `kappa` changes this branch without
changing its first curvature response. CID-D3 already distinguishes
calibration from held-out prediction; we reuse that principle here.

## 2. Exact obstruction to the fundamental-quantum interpretation of the old relaxation

EGR-1 explicitly selected a gradient relaxation. At fixed `beta=1` its
pure-cut vector field is

\[
 H(P)=-\left(\frac{r_1}{10}Z+\frac{r_2}{47}X\right),\qquad
 F(P)=[[P,H(P)],P].                                             \tag{3}
\]

For real cuts `P=(I+r1 Z+r2 X)/2`, direct multiplication gives

\[
 F(P)=g r_1r_2^2 Z-g r_1^2r_2 X,\qquad g=37/470.                \tag{4}
\]

**Theorem 1 (ensemble obstruction).** No affine, decomposition-independent
initial-velocity map on density operators can agree with (3) on every
pure cut while acting componentwise on classical randomizations.

**Proof.** Define `P(z,x)=(I+z Z+x X)/2`. The ensembles

\[
 \mathcal E_A=\tfrac45 P(1,0)+\tfrac15P(-1,0),\qquad
 \mathcal E_B=\tfrac12P(3/5,4/5)+\tfrac12P(3/5,-4/5)
\]

have the same density `diag(4/5,1/5)`. Equation (4) yields

\[
 \overline F_A=0,\qquad
 \overline F_B=\frac{888}{29375}Z\ne0.                          \tag{5}
\]

An affine velocity map assigns the same velocity to equal inputs, which
contradicts (5). In particular the average `Z` reading has initial rates
`0` and `1776/29375`. This is an exact finite contradiction, independent
of numerical tolerances. QED.

This theorem does not invalidate EGR-1's gradient-flow mathematics.
Nor does it prohibit all nonlinear effective equations, preparation-aware
feedback or conditional evolution. Those have different input contracts.
If an external controller knows a source preparation, that classical
record must be included, and its operation is not an autonomous channel
of the unknown source alone. The theorem also does not establish an actual
superluminal-signalling protocol without additional bipartite premises.

## 3. A fixed filter with a complete quantum outcome record

Let `Q=P+ = (I+Z)/2`, `Qbar=P- = I-Q`, and choose real nonnegative `a,b`
with `a^2+b^2=1`. Define

\[
 K_s=Q+a\bar Q,\quad K_f=b\bar Q,\qquad
 \mathcal E_a(\rho)=K_s\rho K_s^\dagger+K_f\rho K_f^\dagger.
                                                               \tag{6}
\]

**Theorem 2 (instrument and unitary realization).** Equation (6) is an
affine, completely positive, trace-preserving operation. It has the
following exact unitary realization with a two-component outcome register:

\[
 W=Q\otimes I+\bar Q\otimes R_a,\qquad
 R_a=\begin{pmatrix}a&-b\\b&a\end{pmatrix}.                     \tag{7}
\]

**Proof.** Orthogonality of the two cuts gives
`Ks^dagger Ks+Kf^dagger Kf=Q+(a^2+b^2)Qbar=I`.
Each amplified branch maps a positive factor `B B^dagger` to
`(K_j tensor I) B B^dagger (K_j^dagger tensor I)`, which is positive
for every ancillary dimension. Trace preservation follows by cyclicity.
Linearity is explicit in (6). Since `Ra^dagger Ra=I`, (7) satisfies
`W^dagger W=W W^dagger=I`. With the register initially in its first
basis cut, its two output blocks are exactly `Ks` and `Kf`. The full
joint operation is reversible; reduction to one register need not be.
Finally, for a remote observable `B`,

\[
 \sum_j K_j^\dagger K_j\otimes B=I\otimes B,
\]

so the remote marginal is unchanged by the nonselective local operation.
Conditional marginals need not be unchanged. QED.

The construction is covariant under a common change of carrier frame:
conjugate `rho,Q,Ks,Kf` by the same unitary. All trace readings agree.
It is a **selected alternative operation**, not an affine extension
retaining the state-dependent vector field (3), which Theorem 1 excludes.

### 3.1 Retain amounts before normalizing branches

The success and failure sources, before normalization, have moments

\[
 \begin{aligned}
 N_s&=\tfrac12[N+z+a^2(N-z)],&
 z_s&=\tfrac12[N+z-a^2(N-z)],& x_s&=ax,&y_s&=ay,\\
 N_f&=\tfrac12 b^2(N-z),&z_f&=-N_f,&x_f&=0,&y_f&=0.
 \end{aligned}                                                  \tag{8}
\]

Consequently

\[
 (N,z,x,y)\longmapsto(N,z,ax,ay).                               \tag{9}
\]

The branch probability is `N_j/N`, not `N_j` when the input amount is
not one. For any nonzero branch, `sigma_j=rho_j/N_j`; then
`rho_out=N_s sigma_s+N_f sigma_f`. No amount is recreated by normalization.
The identity

\[
 N_s^2-z_s^2-x_s^2-y_s^2
 =a^2(N^2-z^2-x^2-y^2)                                         \tag{10}
\]

also shows directly that a successful pure source remains pure.

For `a=3/5`, `b=4/5`, input `N=1`, `r=(3/5,4/5,0)`,

\[
 N_s=109/125,\quad N_f=16/125,\quad
 \mathbf r_s=(91/109,60/109,0),\quad
 \mathbf r_{\rm out}=(3/5,12/25,0).                            \tag{11}
\]

The successful source looks more aligned with `Z`, but the unconditional
`Z` moment is unchanged. Using its normalized direction while assigning
it the original total amount would give a false gravitational source
increment. For the linear readout (2), keeping `alpha/gamma` fixed,

\[
 y_1=(\alpha/\gamma)z/10,\quad y_2=(\alpha/\gamma)x/47,
 \qquad y_{1,\rm out}=y_{1,\rm in},\quad
 y_{2,\rm out}=a\,y_{2,\rm in}.                               \tag{12}
\]

Equation (12) is an ensemble readout identity. It does not specify how an
independent spacetime field responds to individual measurement outcomes.

### 3.2 Exponential composition and the meaning of conditional relaxation

Choose `a=exp(-s)` with `s>=0`, a dimensionless operation parameter.
Then `Ks=exp(-s Qbar)`. The unconditional maps obey

\[
 \mathcal E_s\mathcal E_t=\mathcal E_{s+t},\qquad
 \frac{d\rho}{ds}=\tfrac12(Z\rho Z-\rho)
                =-\tfrac14[Z,[Z,\rho]].                       \tag{13}
\]

This follows from (9), not from an assumed physical clock. Repeated
reduced operations use fresh outcome registers; the family of unitary
dilations (7) is not itself asserted to be a one-parameter group on a
single fixed register.

After `m` successive successes, `Ks^m=Q+a^m Qbar`. In a protocol that
stops at its first failure, the failure branches are `Kf Ks^j`,
`j=0,...,m-1`. Their completeness is the geometric-series identity

\[
 a^{2m}+b^2\sum_{j=0}^{m-1}a^{2j}=1.                          \tag{14}
\]

For normalized input with `p_+=Tr(Q sigma)` and `p_-=1-p_+`, success
probability is `p_+ + a^(2m)p_-`. If `0<a<1` and `p_+>0`, the normalized
successful state tends to `Q`, with limiting probability `p_+`.
If `p_+=0`, its probability tends to zero and the surviving normalized
state does not tend to `Q`. Thus concentration cannot be asserted at
unit probability for arbitrary input.

For a fixed positive Hermitian `H`, normalized filtering by `exp(-s H)`
satisfies

\[
 \sigma'=-\{H,\sigma\}+2\operatorname{Tr}(H\sigma)\sigma.
                                                               \tag{15}
\]

For a rank-one cut, `P H P=Tr(H P) P`, so (15) equals `[[P,H],P]`.
This derives a legitimate **conditional** double-commutator flow with
a fixed operator. It does not permit replacing `H` by the unknown
state-dependent `H(P)` of (3) while keeping a single fixed instrument.

## 4. Response observables and the cost of a simultaneous classical reading

Hold `N`, hence `chi`, fixed. Equation (2) for every source state defines
the pulled-back response observables

\[
 A_1=\frac{\chi}{10}Z,\qquad A_2=\frac{\chi}{47}X,\qquad
 [A_1,A_2]=\frac{\iota\chi^2}{235}Y.                          \tag{16}
\]

These are observables on the source algebra. Their noncommutation does
not by itself construct independent quantum gravitational degrees of
freedom. Mean response alone also does not specify their joint output
statistics in an actual experiment.

**Theorem 3 (added-noise bound).** Suppose a joint readout is represented
by an isometry `V` into a larger finite positive carrier, with commuting
Hermitian pointers `B1,B2`, and is exactly unbiased for all source states:

\[
 V^\dagger V=I,\qquad V^\dagger B_i V=A_i.
\]

Define its added second-moment noise operators

\[
 \mathsf N_i=V^\dagger B_i^2V-A_i^2.
\]

Then `N_i>=0` and, in the induced operator norm,

\[
 \boxed{\|\mathsf N_1\|\,\|\mathsf N_2\|
             \ge\frac{\chi^4}{220900}.}                       \tag{17}
\]

In particular both responses cannot have an exactly unbiased, noiseless,
commuting joint reading when `chi != 0`.

**Proof.** Set `P=V V^dagger`, `Li=(I-P) Bi V`. Direct multiplication gives

\[
 \mathsf N_i=L_i^\dagger L_i,\qquad
 [A_1,A_2]=L_2^\dagger L_1-L_1^\dagger L_2.                    \tag{18}
\]

Positivity follows from the first identity. The triangle and operator-norm
product inequalities give `||[A1,A2]|| <= 2 ||L1|| ||L2||`.
Since `||Li||^2=||N_i||` and `||Y||=1`, (16) gives (17). QED.

This is a bound on universal unbiased measurement noise, not the claim
that every source has a nonzero commutator expectation. Intrinsic state
variance and added pointer noise are distinct. At `chi=0` the response
vanishes and (17) is vacuous; this does not make the entire source theory
classical or prove a classical-gravity limit.

### 4.1 Exact saturating four-outcome construction

For signs `s,t` in `{+1,-1}`, define

\[
 E_{st}=\tfrac14(I+s\tfrac35 Z+t\tfrac45 X),\qquad
 b_{1,st}=s\chi/6,\quad b_{2,st}=5t\chi/188.                  \tag{19}
\]

Each `P_st=2E_st` is a rank-one cut, so every effect is positive and
the four effects sum to `I`. Their weighted first moments are precisely
`A1,A2`. Their noise operators are

\[
 \mathsf N_1=\frac{4\chi^2}{225}I,\quad
 \mathsf N_2=\frac{9\chi^2}{35344}I,\quad
 \|\mathsf N_1\|\|\mathsf N_2\|=\frac{\chi^4}{220900}.        \tag{20}
\]

An exact native dilation avoids irrational coefficients: take
`K_st=((1+iota)/2) P_st`, whose squared modulus factor is `1/2`,
and stack the four `K_st` into an `8 x 2` isometry. On its four outcome
blocks the diagonal pointers have the values (19). The certificate
verifies all identities in (18)–(20) directly.

For the source in (11), ordered outcomes `++,+-,-+,--` have probabilities
`1/2,9/50,8/25,0`. The naive replacement of both coefficients `3/5,4/5`
by `1` gives a negative probability `-1/10` on the opposite source,
and is not a valid effect model. General joint-measurement limitations
are established quantum theory, not a new gravitational uncertainty
principle claimed here.

## 5. Identify the nonlinear source law and retain a held-out prediction

Retain precisely the CRC-1 family, its labelled nine-cell complex and
source direction `r=(3/5,4/5)`. In its notation,

\[
 F_\kappa(h)=r^T A h+\frac\kappa2(\Gamma h)^T M(r)(\Gamma h),
 \quad d_\kappa=d_0+\frac{G\eta_0}{1-\kappa},\quad
 d_0=SK^{-1}r.                                                 \tag{21}
\]

For `kappa != 1`, the predecessor's local analytic branch has
`q(epsilon)=q0+epsilon d_kappa+O(epsilon^2)`. No uniform neighbourhood
in `kappa` or extension through `kappa=1` is assumed. Let an experimental
or simulated control `t` obey `epsilon(0)=0`, `epsilon'(0)=c != 0`.
For each cell class define

\[
 v_i=\left.\frac{d\log V_i}{dt}\right|_{t=0},\qquad
 (i=o,c\!a,i\!n)=(\text{outer},\text{cap},\text{inner}).
\]

There are four outer cells, one cap and four inner cells. To avoid using
`c` for both the scale and the cap, write their slopes as `v_o,v_cap,v_in`.

For a cell Gram matrix `Gcell`, `d log V=(1/2)Tr(Gcell^-1 dGcell)`.
Applying this directly to the existing certified geometry gives

\[
 \mathbf v=c(\mathbf a+\tau\mathbf b),\qquad \tau=(1-\kappa)^{-1},
\]
\[
 \mathbf a=\left(\frac3{10},0,\frac{64}{47}\right),\qquad
 \mathbf b=\left(-\frac{456869}{4241280},
                 \frac{4410869}{1060320},
                -\frac{531631}{1060320}\right).                \tag{22}
\]

**Theorem 4 (two-calibration identification and held-out slope).** Under
the preceding assumptions, the outer and cap slopes identify `c,kappa`:

\[
 c=\frac{10}{3}\left(v_o+\frac{456869}{17643476}v_{\rm cap}\right),
 \qquad
 \kappa=1-\frac{4410869\,c}{1060320\,v_{\rm cap}}.              \tag{23}
\]

The third, unused slope must satisfy

\[
 \boxed{v_{\rm in}=\frac{640}{141}v_o
               -\frac{1860931}{621932529}v_{\rm cap}.}         \tag{24}
\]

**Proof.** The first two rows of (22), considered as equations for
`c` and `c tau`, have determinant `a_o b_cap != 0`.
Eliminating `c tau` gives the first formula (23). The cap equation gives
`tau=v_cap/(b_cap c)`, whence the second formula. Expressing the last row
of (22) as a linear combination of its first two gives (24). The exact
coefficient identities are also serialized in the native proof compiler.
QED.

For finite `kappa != 1` and nonzero `c`, the cap slope is nonzero.
Zero total signal cannot identify the coupling. Nonzero fitted `c` with
zero cap slope is outside this finite-parameter family, not evidence
that a particular finite `kappa` has been measured.

The certificate rederives (22) from all nine Gram matrices, reproduces
the old `kappa=0,2` results, and checks inverse reconstruction also for
`kappa=3,-1` and a second source scale. Altering only the inner slope by
`1/100` keeps the two fitted parameters exactly unchanged but violates
(24) by `1/100`. Thus calibration cannot manufacture this validation.

Equation (24) tests this particular source-law/triangulation contract.
It is not universal under changes of source direction, triangulation,
higher-order source terms affecting this derivative, or observable
definitions. The cell volumes are registered discrete-model quantities;
their status as measurable, relational gravitational observables has not
been derived. Neither an experimentally usable finite-step error bound
nor a mapping of `epsilon` to mass, length or laboratory time follows
from the infinitesimal identity alone.

## 6. Verification boundary and prior art

`verify.cjs` uses the unchanged native field and matrix engine, with no
floating-point input. Polynomial identities are reduced in presentations
whose critical diamonds are checked; serialized witnesses are replayed.
All source files and the consumed predecessor certificates are hash pinned.
The separate `--check` process recomputes the certificate bytes.

The general arguments about positive amplification, trace preservation,
operator norms, conditional normalization and limits are written proofs
in this document. They are not represented as machine-checked analytic
proofs. The certificate does not establish experimental observations,
global existence, Lorentzian dynamics or a completed theory of gravity.

The mixture-linearity issue has a long history. Simon, Buzek and Gisin,
[The no-signaling condition and quantum dynamics](https://arxiv.org/abs/quant-ph/0102125),
derive linear completely positive dynamics under explicit quantum-state,
probability and no-signalling assumptions. Our finite counterexample
uses only the mixture contract stated here and does not rely on a
stronger or interpretation-independent version of their conclusion.

Busch and Heinosaari,
[Approximate joint measurements of qubit observables](https://arxiv.org/abs/0706.1415),
provide established joint-measurement context. The local contribution is
the exact specialization to the existing response modes and the native
saturating witness, not a claim to have invented quantum measurement
uncertainty, instruments or dilation.
