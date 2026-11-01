# YM-75: using the RH programme's source-resolved spectral tools

Monty Dabas. 5 October 2026. **Python 3.12 only.**

The owner's instruction was to check the Riemann-hypothesis programme
and reuse its proved tools in the actual Yang--Mills problem. This stage
reads the current source repositories, matches the hypotheses, and
applies three tools: source-resolved Schur inertia, certified source
solves, and stable recovery from retained observations.

The first application improves YM74's actual side-four SU(2) gap bound
at theta=1/48 to

\[
 \boxed{0.49069036\le\Delta_{4^3}(1/48)\le3.00712209.}       \tag{1}
\]

These are bare energy units for the full physical link Hilbert space.
The bound uses the same complete 241-dimensional electric cutoff as
YM74. It is sharper because the mixed source matrix is retained instead
of replaced by its largest eigenvalue. Neither a larger cutoff nor a
computed interacting vacuum is needed.

The second application gives a quadratic-residual, one-sided matrix
enclosure for the exact hidden memory. The third turns YM74's whole-band
tail bound into an exact stable-decoder theorem. The actual finite
positive continuum gap and nontrivial local 4D fields remain unproved.

## 1. Source audit and attribution

The audit covered the theorem indexes and target-relevant statements
in MP, RH-Framework, RH-project, and Recognition-Kernel-Framework (RKF).
The collected volume in this repository has 128 source renderings; that
is a navigation aid, not a count of newly verified theorems. This audit
does not independently re-prove the entire RH programme.

Current source pins used below:

- MP main: b8dff4583952723a84b639fd91758d9e41c4f628.
- RKF main: 3cc5a33b05c16d59c90994ddda69dedc0d392424.
- RH-Framework T09--T11 development snapshot:
  c8d7a1eb24fa4f22eeaff06d9d8ebdfd9106820a; this is a branch snapshot,
  not a claim that those files are on its default branch.

| Source result | Exact useful content | Yang--Mills application or remaining hypothesis |
|---|---|---|
| MP [finite-rank inertia and residual enclosure](https://github.com/Parveen117/MP/blob/b8dff4583952723a84b639fd91758d9e41c4f628/papers/rh_initial_journey/sections/07j_virial_birman_schwinger_final.tex) | A positive source minus a finite source Gram is decided by its source-resolved Birman--Schwinger matrix; solve residuals give outward error | T1 below retains the actual YM source matrix; T3 encloses approximate hidden solves |
| MP [target-faithful reconstruction](https://github.com/Parveen117/MP/blob/b8dff4583952723a84b639fd91758d9e41c4f628/papers/rh_initial_journey/sections/07c_spectral_blindness_and_freeze.tex) | Missing directions require a compatible operator-tail or hidden lower bound | YM74 supplies an independent electric hidden bound; changing a chart alone cannot supply it |
| MP [common-matrix covariance](https://github.com/Parveen117/MP/blob/b8dff4583952723a84b639fd91758d9e41c4f628/papers/rh_initial_journey/sections/07k2_common_five_matrix_convergence.tex) | Faithful unitary transport preserves the source Gram and threshold count | Supports the symmetry decomposition in T2; it does not change source rank |
| RH-Framework [T11 scaled-Gram solve modulus](https://github.com/Parveen117/RH-Framework/blob/c8d7a1eb24fa4f22eeaff06d9d8ebdfd9106820a/docs/T11_UNIFORM_SCALED_GRAM_SOLVE_MODULUS.md) | Residuals in the correct dual metric control the compressed inverse | T3 uses the actual YM hidden floor and residual Gram, with no borrowed RH constants |
| RH-Framework [T09 diagonal completion](https://github.com/Parveen117/RH-Framework/blob/c8d7a1eb24fa4f22eeaff06d9d8ebdfd9106820a/docs/T09_DIAGONAL_T03_COMPLETION_SCHEDULE.md), [T10 adaptive refinement](https://github.com/Parveen117/RH-Framework/blob/c8d7a1eb24fa4f22eeaff06d9d8ebdfd9106820a/docs/T10_ADAPTIVE_E5_NUMERICAL_REFINEMENT.md) | Choose joint refinement only after uniform error moduli and margins are established | Supplies an error-allocation pattern; does not establish the YM physical spectral margin |
| RKF [N10 uniform observer control](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/theory/02_COMPLETION_AND_MEMORY.md) | Bounded recovery exists iff target size is uniformly controlled by readout size | T4 discharges that hypothesis for each actual bounded-energy band using YM74 |
| RKF [T28 finite-to-infinite cut theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/28_recognition_complete_finite_to_infinite_cut_theorem.md) | Common-carrier norm convergence, faithful memory, reference floor and outward spectral margin | A framework for the continuum comparison; common transports and the uniform physical margin are still needed |
| RKF [T49 local-to-uniform gap](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/49_local_to_uniform_seam_gap_theorem.md) | Exact invariant product channels give a size-independent contraction; interacting corrections need a uniform reserve | Actual plaquettes share links. Independent product channels and a volume-uniform correction bound cannot be assumed |
| RKF [T39 no-blindness and strict positivity](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/39_native_source_kernel_no_blindness_strict_odd.md) | Source injectivity and strict relative positivity are distinct from uniform ambient coercivity | Prevents confusing recovery/injectivity with a positive physical mass gap |

Two source-version findings matter.

First, Publications' RH-paper pointer was stale. MP's current
[fixed-node section](https://github.com/Parveen117/MP/blob/b8dff4583952723a84b639fd91758d9e41c4f628/papers/rh_initial_journey/sections/07j1_fixed_shift_outward_final.tex)
reports the native outward upper bound 0.8388901637561339<1 at
eta=0.01, while its
[endpoint criterion](https://github.com/Parveen117/MP/blob/b8dff4583952723a84b639fd91758d9e41c4f628/papers/rh_initial_journey/sections/07l_endpoint_criterion_final.tex)
still requires the sign of K_0. The pointer is updated to describe
that source status. This YM audit has not replayed the entire private
RH numerical certification chain.

Second, use RKF's
[corrected archived T50](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/operator_foundation/sources/rkf_reference/theorum/50_native_seam_gap_odd_sector_covariance_theorem.md),
which withdraws exact tensor-mass multiplication and retains its
submultiplicative bound. For D(a+i b)=|a|+|b|,
D((1+i)^2)=2 while D(1+i)^2=4. The historical top-level copy is
not a replacement for the corrected source.

The RH source rank at most five comes from its specific boundary and
accepted-source structure. No five-dimensional or rank-one YM source
is inferred from that fact. At YM74's actual R=3 cutoff, the mixed
source has rank 240 when theta is nonzero.

These are applications of established variational, Schur,
Birman--Schwinger and bounded-recovery principles already developed in
the programme. No priority claim is made for those general principles.

## YM75-T1. Preserve the whole source matrix in the spectral enclosure

Use the unshifted actual finite-lattice Hamiltonian and the complete
physical electric cut from YM74:

\[
 A=\begin{pmatrix}C&B\\B^*&D_A\end{pmatrix},\qquad
 D_A\ge dI,\qquad G\ge BB^*,\qquad G\ge0.                 \tag{2}
\]

The retained space has dimension n>=2; C is finite self-adjoint,
B is bounded, and the hidden operator is self-adjoint on its full
domain. YM has compact resolvent at each finite regulator.
Let c_0<=c_1 be the first two Ritz eigenvalues of C and define the
auxiliary comparison matrix

\[
 {\cal L}_G=
 \begin{pmatrix}C&G^{1/2}\\G^{1/2}&dI_n\end{pmatrix}.
                                                               \tag{3}
\]

If ell_0<=ell_1 are its first two eigenvalues, then

\[
 \boxed{\ell_j\le E_j(A)\le c_j,\qquad j=0,1.}            \tag{4}
\]

Consequently

\[
 \boxed{\max(0,\ell_1-c_0)\le\Delta\le c_1-\ell_0.}       \tag{5}
\]

**Proof.** Interlacing with the auxiliary dI_n block gives
ell_j<=d for j=0,1. For z<d the actual Schur complement satisfies

\[
 S_A(z)=C-z-B(D_A-z)^{-1}B^*
 \ \ge\ C-z-\frac{G}{d-z}=:S_G(z).                     \tag{6}
\]

The right side is the Schur complement of (3) at z. Below ell_j,
S_G(z) has at most j negative eigenvalues, so the same holds for
S_A(z). The hidden block is positive; Schur inertia and a limit in
z prove the lower bounds. Rayleigh--Ritz proves the upper bounds.
QED.

If q=||G||, the same argument with G<=qI gives

\[
 F(c_j,d,q)\le\ell_j,\qquad
 F(c,d,q)=\frac{c+d-\sqrt{(c-d)^2+4q}}2.                \tag{7}
\]

Thus the source-resolved bound is at least as strong as the scalar
bound using the same d and q. No commutation of C and G is required.
It preserves directions in which the source is absent.

The comparison matrix is not an exact finite replacement for D_A.
For singular G, its uncoupled auxiliary levels at d must not be
discarded. With G=0 and c_0>d, the lower levels must still be capped
at d. Keeping (3), or checking (6) at rational z<d, avoids this error.

In a nonorthogonal native basis with positive Gram W, use

\[
 C_{\rm raw}-zW-\frac{G_{\rm raw}}{d-z}.                \tag{8}
\]

Congruence by W^{-1/2} preserves inertia. Exact rational elimination
therefore implements the test without computing G^{1/2}.
Here C_raw and G_raw are the matrices of the corresponding quadratic
forms in the chosen basis. Products of raw coefficient matrices
cannot omit the intervening Gram factors.

## YM75-T2. Actual side-four lattice improvement

YM74 proves completeness of the R=3 physical cut and computes

\[
 BB^*=\theta^2\operatorname{diag}
   \left(0,\frac{191I+J}{4},48I\right).                 \tag{9}
\]

Its sectors are the constant, 192 elementary plaquettes and 48 winding
four-cycles. Split the plaquette sector into its uniform vector and
191 transverse directions. Set

\[
 k^2=48\theta^2,\quad q_u=\frac{383}{4}\theta^2,\quad
 q_t=\frac{191}{4}\theta^2,\quad q_w=48\theta^2,\quad
 d=\frac92-192|\theta|.                                \tag{10}
\]

The lower Schur comparison consists of a constant/uniform block

\[
 \begin{pmatrix}
  -z&-k\\
  -k&3-z-q_u/(d-z)
 \end{pmatrix},                                       \tag{11}
\]

plus 191 copies of 3-z-q_t/(d-z) and 48 copies of
3-z-q_w/(d-z). This is a decomposition of the comparison pencil,
not a claim that the interacting hidden Hamiltonian decouples.

For theta=1/48, d=1/2. At the exact rational points

\[
 z_0=-0.00712209,\qquad z_1=0.48376192,                 \tag{12}
\]

all pivots are positive at z_0. At z_1 the first uniform pivot -z_1
is negative, the second

\[
 3-z_1-\frac{q_u}{d-z_1}+\frac{k^2}{z_1}
\]

is positive, and every transverse/winding entry is positive.
Thus the full lower comparison has respectively zero and one
negative direction. Exact rational signs, not rounded eigenvalues,
certify

\[
 E_0(A)\ge z_0,\qquad E_1(A)\ge z_1.                   \tag{13}
\]

The retained eigenvalues are
c_0=(3-sqrt(9+192theta^2))/2 and c_1=3. Directed arithmetic gives
c_0<=-0.00692844. Substituting in (5) proves (1).

For an independent algebraic check, the active comparison block
has rational characteristic polynomial

\[
 9216z^3-32256z^2+13249z+96.                            \tag{14}
\]

Its first two roots are approximately -0.00712208573309 and
0.483761924917664. These decimal roots explain the chosen rational
thresholds; the certificate uses the exact signs in (11)--(13).

| Same actual lattice, coupling and retained space | Bare gap enclosure |
|---|---|
| YM74, scalar mixed-norm enclosure | [0.4904, 3.07874] |
| YM75, source-resolved enclosure | [0.49069036, 3.00712209] |

The source matrix annihilates the retained constant but not the
uniform plaquette vector. C mixes those vectors, so C and BB* do
not commute. Replacing BB* by qI charges the constant for a direct
hidden coupling that is absent. Keeping the matrix removes that
loss of information. The actual cube result remains a finite-volume
bound, with no extrapolation to a physical continuum mass.

## YM75-T3. Quadratic-residual enclosure of the hidden memory

The RH residual method controls solves in the correct native norm.
For the YM application let S=D_A-z>=sI>0, where s=d-z and z<d.
Put F=B^*. For any finite-column trial X whose range is in D(S),
define

\[
 R=F-SX,\qquad
 T_X=F^*X+X^*F-X^*SX,\qquad M=F^*S^{-1}F.             \tag{15}
\]

Then the exact identity and Loewner bounds are

\[
 \boxed{M=T_X+R^*S^{-1}R,}
\]
\[
 \boxed{T_X\le M\le T_X+s^{-1}R^*R.}                  \tag{16}
\]

**Proof.** Expand
(F-SX)^*S^{-1}(F-SX). The cross terms are -F^*X-X^*F,
and the last term is X^*SX. Positivity and S^{-1}<=s^{-1}I
give (16). QED.

In particular, the enclosing matrix width is at most ||R||^2/s.
The simpler Hermitian centre sym(F^*X) has the familiar linear
error bound ||F|| ||R||/s from the RH source-solve theorem.
The new use here is the variational centre T_X, which keeps the
error positive and quadratic rather than replacing its orientation
by a scalar norm.

The actual Schur matrix is therefore enclosed by

\[
 C-z-T_X-s^{-1}R^*R\ \le S_A(z)\le C-z-T_X.             \tag{17}
\]

The two finite endpoints support exact inertia tests. At X=0,
(16) reduces to 0<=M<=BB*/s, recovering T1. Improving X can
reduce the interval width, although it need not monotonically
improve every chosen lower endpoint.

For the actual electric cut, F has a finite native electric shell
by YM74-T4. Choosing X with finite polynomial columns in the
hidden physical space makes SX, R and their Gram matrices finite
native calculations too. No tail is dropped: the finite-shell
statement concerns application to these columns, while S^{-1}
remains the full hidden inverse.

Along a chosen physical trajectory, a sufficient source-solve
allocation is

\[
 \kappa_a\,\|R_a\|^2/s_a\le\varepsilon_a,\qquad
 \varepsilon_a\longrightarrow0.                        \tag{18}
\]

Here s_a must already be a proved floor, and the residual must be
enclosed in the actual Hilbert metric. A Euclidean coefficient
residual in a nonorthogonal basis does not certify (18).
This is an error-control theorem; it does not prove a positive
retained spectral margin or efficient complexity.

## YM75-T4. Stable recovery of every retained bounded-energy state

Let E_M=1_[0,M](kappa H) be the actual physical spectral projection,
P=P_R the known electric cut, and use YM74's explicit budget

\[
 \zeta=\frac{\delta+M/\kappa}{\lambda_>(R)}<1.            \tag{19}
\]

For every psi in ran E_M, YM74 gives ||(I-P)psi||^2<=zeta||psi||^2.
Orthogonality then yields

\[
 (1-\zeta)\|\psi\|^2\le\|P\psi\|^2\le\|\psi\|^2.       \tag{20}
\]

Apply RKF N10 with observation C=P|ran E_M and target L=inclusion.
There is an exact bounded decoder D_R on the band image with

\[
 \boxed{D_RP\psi=\psi,\qquad
       \|D_R\|\le(1-\zeta)^{-1/2}.}                    \tag{21}
\]

At finite volume an extension to the whole retained space is

\[
 D_R=(E_MPE_M|_{\operatorname{ran}E_M})^{-1}E_M
       \big|_{\operatorname{ran}P}.                    \tag{22}
\]

Indeed T=E_MPE_M>=1-zeta on the band, D_RPE_M=E_M,
and D_RD_R^*=T^{-1}. This proves (21). It also gives
dim ran E_M<=rank P_R and amplifies retained readout error by
at most (1-zeta)^(-1/2).

Under YM74's growing-cutoff recipe zeta tends to zero for every
fixed M. The observation becomes an asymptotic isometry and the
recovery condition tends to one.

Formula (22) contains the unknown interacting projection E_M.
It proves existence and conditioning, not an executable
vacuum-free recovery algorithm. Nor does it produce an excited
state in that band: the band may contain only the vacuum.

## 2. What the RH tools do and do not close

The sources genuinely help: they tighten an actual lattice bound,
give a matrix-valued solve-error contract, and convert whole-band
coverage into stable exact recovery. T28 and T09 then state the
remaining uniformity requirements explicitly.

A finite positive physical mass still needs control of

\[
 \kappa_a\bigl(E_1(a)-E_0(a)\bigr)
\]

on a physical weak-coupling trajectory as spatial volume grows
and spacing shrinks. A nontrivial local field limit, with the
required four-dimensional structure and a surviving finite-energy
local excitation, is also still needed.

RH's source positivity, its rank-five structure, or a change of
chart cannot be assigned to this different interacting gauge
operator without a theorem identifying the corresponding sources.
Strict positivity for every nonzero vector is also insufficient:
the diagonal operator diag(1,1/2,1/3,...) is injective and strictly
positive as a quadratic form but has no positive spectral gap.

The work stays in Publications. The source repositories and their
shared engine remain unchanged; their exact paths and commit pins
are recorded for attribution.

## Reproduction

The Python certificate verifies the exact cube pivots, source-aware
inertia comparisons, independent finite-matrix spectra, quadratic
residual identities and stable-recovery examples. Finite matrix
fixtures are labelled as such and are not Yang--Mills vacua.
Source hashes bind the written argument and its computational
controls. This is not proof-assistant verification.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym75_rh_source_resolved_transfer.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym75_rh_source_resolved_transfer.py -v
~~~
