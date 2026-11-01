# YM-52: native energy, discarded records and observer-dependent relaxation

Monty Dabas. 3 October 2026. Verification: **Python 3.12 only**.

## Result and lineage

YM51 left the symmetric heat tensor C as protocol data. This chapter derives
what that data does: its derivative energy is the rate at which a mean
readout loses squared distinction; the complementary record variance keeps
an exact account of that loss. A positive density has a corresponding
relative-entropy dissipation law. The weighted native order brackets
construct the cofactor tensor, and two active independent directions suffice
for relaxation on this compact carrier.

For eigenvalues 0 <= a <= b <= c of C, the exact optimal centered relaxation
constants are

\[
\boxed{\gamma_{\rm full}=\min\{(a+b+c)/4,\ a+b\},\qquad
       \gamma_{\rm even}=a+b,\qquad
       \gamma_{\rm odd}=(a+b+c)/4.}
\tag{1}
\]

Consequently a sign-insensitive observer and the full observer have
different fixed-budget optimization problems. Maximizing the even-sector
gap selects isotropy; maximizing the full-sector gap does not.

These are results for YM50--51's specified compact coefficient carrier,
reference and free symmetric turn protocols. They are not an extension
of YM45's interacting-chain gap, a physical energy-unit calibration or a
four-dimensional field-theory construction.

The positive-definite full/even spectral formulas are known mathematics:
Emilio A. Lauret, *The smallest Laplace eigenvalue of homogeneous 3-spheres*,
Theorems 1.1--1.2,
[arXiv:1801.04259v3](https://arxiv.org/abs/1801.04259v3),
[doi:10.1112/blms.12213](https://doi.org/10.1112/blms.12213).
His derivative normalization is twice ours; his metric parameters squared,
in reversed order, are (c,b,a). Dividing his eigenvalues by four gives
(1). No priority claim is made for that formula, diffusion energy,
entropy dissipation or bracket generation. The proof below also covers
singular C directly and connects these tools to the native record protocol,
the existing observer and the explicit selection question.

Existing energy is not discarded or renamed. R18 and R31 have conserved
source-pairing energies; CF-3 specifies a thermodynamic formation energy;
PF-5 has a conserved field energy with source work; YM48-T6 already proves
the isotropic derivative energy and interacting ground-source identity.
This chapter supplies a different missing bridge: the anisotropic free
protocol's retained/noise/entropy balances and their exact relaxation
constants. The source manifest records these reads and preserves prior
certificates. General written proofs and finite exact controls have
different scopes; neither is a proof-assistant or external expert certificate.

## Contract

Use Q_Sigma, A_coeff, Phi_Q, the uniform completion and the recognition
completion from YM50. Write Phi=Phi_Q and
\[
D_v f(x)=\left.\partial_u f(\operatorname{Exp}_\Sigma(-u v/2)x)
              \right|_{u=0},\qquad [D_u,D_v]=D_{u\times v}.
\]
YM51 supplies C=sum_r w_r v_r v_r^T >= 0,
L_C=-sum_ab C_ab D_a D_b and the positive, unital, Phi-preserving
contraction E_t=Exp_Sigma(-t L_C). The symmetric finite protocols converge
to E_t with the explicit YM51 refinement bound.

Coordinates and coefficients are the completed native radial field;
complexification, when useful for finite spin coefficients, is the
previously derived iota representation. The Hilbert description is
YM50's proved completion representation, not a new primitive.

All time variables here are declared heat ticks. C is a constant internal
tensor, not a spacetime metric. Connection curvature still means
[nabla_u,nabla_v]-nabla_[D_u,D_v]. The raw brackets below are order defects
of the derivative frame; a flat connection may have this noncommuting
frame. No identification of those two curvature targets is made.

## YM52-T1: energy and the exact cost of averaging records

For real coefficient polynomials define
\[
\Gamma_C(f,g)=\sum_{a,b}C_{ab}(D_af)(D_bg),\quad
\Gamma_C(f)=\Gamma_C(f,f),\quad
\mathcal E_C(f,g)=\Phi(\Gamma_C(f,g)).
\tag{2}
\]
The protocol factorization makes Gamma_C(f) nonnegative pointwise.
Leibniz and Phi's turn invariance give
\[
L_C(fg)=fL_Cg+gL_Cf-2\Gamma_C(f,g),\qquad
\mathcal E_C(f,g)=\Phi(fL_Cg).
\tag{3}
\]
For f_t=E_t f, differentiation in its finite coefficient space proves
\[
\boxed{\partial_t\Phi(f_t^2)=-2\mathcal E_C(f_t),\qquad
       \partial_t\mathcal E_C(f_t)=-2\Phi((L_C f_t)^2).}
\tag{4}
\]
The first quantity is retained squared distinction. The second is a
nonnegative derivative energy, or dissipation form; it is not a second
conserved total energy.

To account for the complementary records set
\[
N_t(f)=E_t(f^2)-(E_tf)^2.
\]
Then
\[
\boxed{N_t(f)=2\int_0^t E_{t-s}\Gamma_C(E_s f)\,ds\ge0,}
\tag{5}
\]
\[
\boxed{\Phi(f^2)=\Phi((E_t f)^2)+\Phi(N_t(f)),\qquad
\Phi(N_t(f))=2\int_0^t\mathcal E_C(E_s f)\,ds.}
\tag{6}
\]
Proof: differentiate N_t, apply (3), and solve
N'_t=-L_C N_t+2 Gamma_C(f_t), N_0=0, by the finite coefficient series.
Positivity and invariance of E_t give (5)--(6).

There is an exact finite-record precursor. If p_l are normalized labelled
word counts and z_l=f(q_l x), then
sum_l p_l(z_l-mean z)^2=sum_l p_l z_l^2-(mean z)^2.
Every word preserves Phi. Thus the full averaged squared record retains
Phi(f^2), while the mean and its discarded record variance split it
exactly. Taking YM51's refinement on f and f^2 gives (5).
This noise is derived from a specified averaging operation. It is not a
claim that every physical noise source has been selected by the framework.

## YM52-T2: a scoped thermodynamic dissipation bridge

Let rho be a real coefficient polynomial with Phi(rho)=1 and
0<m<=rho<=M on Q_Sigma. Then rho_t=E_t rho remains between m and M.
Define the dimensionless relative entropy and its production by
\[
\mathcal H(\rho)=\Phi(\rho\log\rho),\qquad
\mathcal I_C(\rho)=\Phi\!\left(\frac{\Gamma_C(\rho)}{\rho}\right).
\]
Here log and reciprocal are continuous scalar functions on [m,M].
Polynomial approximation constructs their readouts from Phi. Uniform
approximation of the first two derivatives, or their convergent local
scalar series on a finite cover, extends the derivative chain rule to
these smooth functions.

Differentiating Phi(rho_t) and using (3) gives
\[
\boxed{\partial_t\mathcal H(\rho_t)=-\mathcal I_C(\rho_t)\le0.}
\tag{7}
\]
In detail, Phi(L_C rho_t)=0, and integration by parts against log rho_t
gives Phi((log rho_t)L_C rho_t)=Phi(Gamma_C(rho_t)/rho_t).
The pointwise convex inequality r log r-r+1>=0 gives H>=0.

Let Var(rho)=Phi((rho-1)^2). Taylor's integral remainder, using
(r log r-r+1)''=1/r, proves
\[
\frac{\operatorname{Var}(\rho)}{2M}\le \mathcal H(\rho)
\le\frac{\operatorname{Var}(\rho)}{2m},\qquad
\frac{\mathcal E_C(\rho)}{M}\le\mathcal I_C(\rho)
\le\frac{\mathcal E_C(\rho)}{m}.
\tag{8}
\]
If the centered energy inequality has constant gamma, then
\[
\mathcal H(\rho_t)\le
e^{-2m\gamma t/M}\mathcal H(\rho).
\tag{9}
\]
Indeed I>=gamma Var/M>=2m gamma H/M, and integrate (7).
For even densities T4 below permits gamma_even instead of gamma_full.
This is a bounded-positive-density estimate, not a claimed global
logarithmic Sobolev constant. For rho=1+epsilon f with Phi(f)=0,
H=epsilon^2 Phi(f^2)/2+O(epsilon^3) and
I=epsilon^2 E_C(f)+O(epsilon^3).

Equations (7)--(9) bridge record relaxation to a defined entropy
functional. Multiplication by temperature or an action/energy scale,
identification with material heat/work, and the system-plus-environment
first law require an additional physical adapter. They are not supplied
by calling the dissipation form energy.

## YM52-T3: the weighted order brackets construct the cofactor tensor

Write b_r=sqrt(w_r)v_r, so C=sum_r b_r b_r^T. Define
\[
\mathcal B(C)=\sum_{r<s}(b_r\times b_s)(b_r\times b_s)^T.
\]
Direct expansion of the two-by-two minors (Cauchy--Binet) yields
\[
\boxed{\mathcal B(C)=\operatorname{cof}(C),\qquad
\sum_{r<s}[D_{b_r},D_{b_s}]f\,[D_{b_r},D_{b_s}]g
=\Gamma_{\operatorname{cof}(C)}(f,g).}
\tag{10}
\]
The formula is polynomial and remains true at singular C. It is independent
of the chosen finite factorization, although the full record protocol is not.
For a rotated diagonal tensor the eigenvalues of the cofactor are
(bc,ac,ab). Thus:

* rank 3: the active directions already span all three;
* rank 2: one pair bracket supplies the missing direction;
* rank 1 or 0: every active pair bracket vanishes, and no new direction is supplied.

This is an algebraic span statement and a quantitative bracket-energy
identity; rank alone gives no uniform numerical lower bound. Cof(C) scales
quadratically when the clock rate C scales linearly. We do not add it
to C as if their time dimensions were identical.

## YM52-T4: exact relaxation for full and sign-insensitive observers

The antipodal operation Zf(x)=f(-x) preserves Phi and commutes with L_C
and E_t. P_even=(I+Z)/2 and P_odd=(I-Z)/2 are orthogonal projections.
The even functions form the sign-insensitive observation algebra.
The odd sector is an orthogonal contrast sector, not a unital algebra.

For the even sector remove its constant. On each of the full centered,
even centered and odd spaces, respectively, the optimal inequality is
\[
\mathcal E_C(f)\ge\gamma\,\Phi(f^2),\qquad
\|E_tf\|_\Phi\le e^{-\gamma t}\|f\|_\Phi,
\tag{11}
\]
with the constants in (1).

**Proof, including arbitrary coefficient content.** Use the existing
finite tensor/coefficient decomposition from YM48--50: a coefficient
polynomial splits into finite spin-j blocks, j=0,1/2,1,...; Phi makes
inequivalent blocks orthogonal. Finite tensor powers exhaust the coefficient
algebra. On a block, -sum D_i^2=j(j+1)I and a unit-axis generator has
weights -j,-j+1,...,j. These follow by differentiating the native two-role
action on its symmetric tensor powers. In particular -D_i^2<=j^2 I.
Z acts as (-1)^(2j).

YM51's native conjugations diagonalize C without changing the pairing.
On an integer block j>=1,
\[
L_C\ge-aD_1^2-b(D_2^2+D_3^2)
=b\,j(j+1)I+(b-a)D_1^2
\ge(aj^2+bj)I\ge(a+b)I.
\tag{12}
\]
On every half-integer block each squared unit-axis generator is at least
I/4: its weights are half-integers and the axes are conjugate.
Hence L_C>=tr(C)I/4 there.
Both bounds are attained. All four linear coefficients have rate tr(C)/4.
For a diagonal C with its largest eigenvalue along axis 3 the polynomial
\[
f_3=x_1^2+x_2^2-\tfrac12\sum_{\mu=0}^3x_\mu^2
\]
is even, Phi(f_3)=0, Phi(f_3^2)=1/12, D_3 f_3=0 and
D_1^2 f_3=D_2^2 f_3=-f_3. Its rate is a+b.
Rotate this witness along with C.

Sum the block inequalities. Differentiating the norm gives (11) on
polynomials. Density and contraction extend it to the corresponding
recognition completions. This establishes a statement on all completed
sources, not merely on the degrees sampled in the executable certificate.
Equivalently the block direct sum defines the closed energy form with
domain sum_blocks <f_j,L_C f_j><infinity.

For rank at least two, a+b>0 and the constants are the sole stationary
sector. For rank one a nonzero centered even witness is stationary.
The zero value in (1) is the gap above the constant sector, not a denial
that a rank-one generator has positive eigenvalues above its larger kernel.
For C=diag(0,3/2,3/2), gamma_full=3/4 despite the missing direct direction.
For YM51's family, gamma_full=min(3/4,2 eta): its earlier upper witness
is now an exact all-source rate when eta<=3/8.

## YM52-T5: observer choice changes the variational selector

Fix the clock budget s=tr(C)>0, and optimize over positive semidefinite C.
This is an explicit performance criterion, not a postulate that nature
optimizes it. Formula (1) proves
\[
\max\gamma_{\rm full}=s/4,\qquad
\gamma_{\rm full}=s/4\ \Longleftrightarrow\ c\le3s/4.
\tag{13}
\]
Thus a continuum of anisotropic protocols, including rank-two examples,
is as fast as isotropy for the full observer's slowest mode.

For the sign-insensitive observer,
\[
\boxed{\max\gamma_{\rm even}=2s/3,\qquad
\gamma_{\rm even}=2s/3\ \Longleftrightarrow\ C=(s/3)I.}
\tag{14}
\]
Proof: gamma_even=s-c and c>=s/3. Equality forces a=b=c.
At fixed s, a deficit epsilon below the optimum determines
c=s/3+epsilon, and
\[
\|C-(s/3)I\|_{\rm op}\le2\epsilon.
\tag{15}
\]
Indeed a=s-b-c>=s-2c=s/3-2epsilon and every eigenvalue is at most c.
The constant two is sharp for eigenvalues
(s/3-2epsilon,s/3+epsilon,s/3+epsilon) while they remain nonnegative.

This supplies an observer-conditioned alternative to YM51's fixed-law
symmetry condition. The additional assumption is now the even-observer
maximin objective, not an assumed isotropic C. Neither that objective nor
the budget s is physically selected here. Maximizing a restricted readout
does not certify that the excluded odd information is absent.
YM51's six quadratic channels recover C and therefore both objective
values; optimizing only a linear decay cannot distinguish any fixed-s C.

## YM52-T6: units, scale and the actual Yang--Mills gate

Scaling C to kC scales E_C(f), both rates and entropy production by k,
while Cof(C) scales by k^2 and E_t^{kC}=E_{kt}^C.
The dimensionless choice s=3 makes comparison with YM45--51 convenient;
it is not a derivation of seconds, joules, hbar or a particle mass.
YM51's interaction ratio theta/k remains necessary when changing clocks.

This chapter strengthens the free compact input: it replaces a slow-mode
upper bound with an exact full-source rate, includes rank-two protocols,
derives the discarded-record and entropy balances, and supplies a precise
observer-dependent selector with a stability bound. It does not transport
these anisotropic rates through an interacting chain automatically.

The next interacting question is whether a declared anisotropic family
admits the positive block-overlap and weighted estimates used by YM45,
with constants controlled by its tensor and normalization. The physical
question is whether the process actually selects the even-observer
criterion, its clock budget and its relative coupling. Those are distinct
tasks; the actual YM49 infinite-chain row defect remains separate.

The [official Jaffe--Witten problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf),
section 4, additionally requires nontrivial four-dimensional quantum
Yang--Mills theory, local fields and the stated axiomatic/short-distance
properties. A compact relaxation gap, or an interacting chain gap, does
not by itself establish those properties. The native construction may
provide ingredients for that bridge; this chapter does not prove it.

## Reproduction and evidence limits

Run the chapter's certificate in check mode, then its unittest file.
The source manifest pins earlier proof/code/evidence and the external
energy reads. The certificate checks exact polynomial energy identities,
independent labelled-record variances, the cofactor identity for different
factorizations, finite spin inequalities, all-source polynomial energy
matrices, sharp eigenfunction witnesses, parity separation, optimization
and its stability bound, and rational entropy-series enclosures.

Negative controls reject rank-one relaxation, an isotropy conclusion from
the full-observer objective, missing factors of two, a universal entropy
claim without positivity bounds, and certificate/source tampering.
Finite checks support the displayed general proofs; they do not replace
the argument across all degrees or certify physical interpretation.
