# YM-53: rank-two native heat survives a weak interaction

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

For the declared open chain, full isotropy is unnecessary for a
width-uniform interacting relaxation bound. At site i let C_i be a fixed
positive semidefinite symmetric protocol tensor with ordered eigenvalues
0<=a_i<=b_i<=c_i. Suppose b_i>=beta>0 uniformly. Rank two is allowed;
the tensors and their orientations may differ between sites. Then

\[
\boxed{|\theta|<\beta/2400}
\]

is a sufficient window for a positive rate gamma, independent of finite
width m and fine heat step 0<h<=1/beta. It bounds the actual square-sourced
transfer on its entire vacuum complement and survives time refinement at
each fixed finite width. The interacting ground-source energy is the
native derivative energy weighted by that ground source, and satisfies
the same Poincare lower bound. The constants are conservative, not a
phase diagram or a physical mass prediction.

## Carrier and dependencies

Use YM50's positive native coefficient functional Phi_Q, uniform
coefficient completion and positive recognition completion. The compact
quaternion/SU(2) coordinates and their coefficient Hilbert space are the
proved representation in YM50, not new primitive postulates. YM51's
symmetric-turn protocol supplies L_C=-sum C_ab D_a D_b and E_t^C.
YM52 supplies its product rule, integration by parts and coefficient
energy. All exponentials and logarithms are the existing completed cut
scalars. We use their represented notation below.

For m sites put L=sum_i L_{C_i}, K_t=tensor_i E_t^{C_i}, and

\[
v_i(x)=\tfrac12\chi_{1/2}(x_ix_{i+1}^{-1})=x_i\cdot x_{i+1},\quad
V=\sum_{i=1}^{m-1}v_i,\quad B=\theta M_V,\quad b_* = |\theta|(m-1).
\]

Thus |v_i|<=1 and ||B||<=b_*. The declared trajectory is kappa(h)=theta h:

\[
S_h=K_{h/2}e^{hB}K_{h/2},\qquad T_h=e^{hB/2}K_he^{hB/2}.
\]

The tensors, beta and theta stay fixed during refinement. There are at
most two spatial neighbours per site. This is the earlier chain carrier,
not the four-dimensional gauge lattice. No common diagonalization of
different sites' tensors or of interaction and heat is assumed.

**Lineage.** The temporal-block construction and full-source extraction
are YM43--45's. Their method belongs to the Dobrushin--Shlosman block
comparison lineage; see Rebeschini and van Handel,
[Comparison Theorems for Gibbs Measures](https://web.math.princeton.edu/~rvan/dobrushin130819.pdf).
No theorem from that paper is substituted for the estimates below.
YM44 supplies ordered insertion packets and norm-gap transport, and
YM48 supplies the ground-source energy argument. The new work verifies
their hypotheses for noncentral, possibly rank-two native heat and
handles zero-normalizer fine bridges without assuming ellipticity.

## YM53-T1: an all-content coarse kernel floor

Rotate one tensor in its native frame. On spin j, with j integer or
half-integer, the coefficient derivative matrices obey the Casimir
identity sum J_k^2=j(j+1)I and 0<=J_k^2<=j^2 I. Therefore

\[
L_C\geq aJ_1^2+b(J_2^2+J_3^2)
\geq (a j^2+b j)I\geq\beta jI.
\]

Indeed b*j(j+1)I-(b-a)J_1^2 has the stated lower bound. This does not
require a>0. By the finite coefficient decomposition inherited from
YM50--52, the nonconstant kernel block at n=2j has absolute contribution
at most (n+1)^2 exp(-beta*t*n/2). This follows from
|tr(AU)|<=dim(A)||A|| for a unitary coefficient matrix U and
A=exp(-t L_C|_j); centrality of A is unnecessary. Hence, with
q=exp(-beta*t/2),

\[
|k_t^C(x,y)-1|\leq G(q):=\sum_{n\geq1}(n+1)^2q^n
=\frac{1+q}{(1-q)^3}-1.                                      \tag{1}
\]

This converges uniformly for every t>0 and constructs a continuous
kernel representing the native positive semigroup. Its row integral is
one, and it is symmetric by the coefficient adjoint identity. It is
nonnegative: a negative value and continuity would give a nonnegative
continuous test function with negative image. Approximate that test
function by squares of coefficient polynomials and use positivity of
E_t^C and the full-support compact reference proved in YM50. This is a
representation argument, not a new choice of measure.

An explicit all-content remainder after degree N>=0 is

\[
q^{N+1}\left(\frac{(N+2)^2}{1-q}
+\frac{2(N+2)q}{(1-q)^2}+\frac{q(1+q)}{(1-q)^3}\right).       \tag{2}
\]

The exact rational inequalities sum_{k=0}^8 (9/2)^k/k!>84 and
G(1/84)<1/20 imply, for every t>=9/beta,

\[
19/20\leq k_t^C(x,y)\leq21/20,\qquad
\epsilon:=4/5\leq(19/21)^2.                                \tag{3}
\]

For any common nonnegative future weight H with positive integral,
the normalized laws k_t(x,y)H(y) and k_t(x',y)H(y) have overlap at
least epsilon. Both normalizers are included in the ratio bound.
No upper bound on tr C is needed for this coarse mixing statement.

## YM53-T2: admissible bridges, including singular fine kernels

Set r=ceil(9/(beta*h)); then 9/beta<=hr<10/beta. In a free bridge with
k interior sites, sample a skeleton every r fine steps. Its next law
is proportional to k_{hr}(x,y)H(y), where H includes the remaining
bridge to the fixed right endpoint. Its integral is positive. Couple
the skeleton by (3), keep it together after meeting, and fill the fine
segments with their exact conditional bridges. Such fillings are needed
only almost surely under their skeleton marginals. The expected Hamming
cost for one changed temporal endpoint is at most

\[
r\sum_{j\geq0}(1-\epsilon)^j=r/\epsilon.                    \tag{4}
\]

For k<r, any two admissible bridge laws instead have cost <=k<r/epsilon.
Crucially, a short bridge can have zero normalizer for some endpoints.
We do not divide by it or assert that every fine kernel is strictly
positive. Conditional laws are only needed at admissible exteriors,
almost surely under each strip law; their definitions on null exteriors
do not affect its marginal. For k>=r all endpoint pairs have positive
normalizer by (3) and the semigroup identity.

For a short block with both temporal endpoints changed, change them
**as a group**, directly coupling the two admissible laws with cost <=k.
Do not insert an intermediate endpoint pair which may be inadmissible.
Then change spatial exterior coordinates one at a time: their strictly
positive interaction tilts preserve the same free-bridge support.
The group cost is bounded by r/epsilon times the number of changed
temporal endpoints. Long blocks permit separate endpoint changes.

The interacting conditional tilt has log magnitude <=2h|theta|k.
Including its normalizer, the tilted/free half-L1 distance is at most
4h|theta|k. Gluing tilted-to-free, free-to-free and free-to-tilted
couplings therefore gives the valid per-endpoint upper cost

\[
D(k)=r/\epsilon+8h|\theta|k^2.                              \tag{5}
\]

One changed spatial exterior coordinate at one time costs at most
4h|theta|k. These are precisely YM45-T2's two different costs, now
justified on admissible supports. If both temporal endpoints change,
(5) charged twice overestimates the grouped cost and remains valid.
Spatial interpolation and the Hamming triangle inequality add all
costs. There is no independence assumption on exterior differences.

## YM53-T3: width-uniform full-source discrete gap

Take ell=Jr, J a positive integer, and use the same clipped block labels
as YM45-T3: starts 2-ell,...,N in each column on a strip with N interior
time rows. Duplicated clipped blocks retain their distinct labels.
Every coordinate is removed ell times, influences at most two temporal
endpoint labels and at most 2ell spatial labels. Consequently

\[
\alpha=\frac{2r}{\epsilon\ell}+24h|\theta|\ell
\leq\bar\alpha:=\frac{2}{\epsilon J}
+240\frac{|\theta|}{\beta}J.                              \tag{6}
\]

Choose R>1 with R*baralpha<1 and put

\[
\gamma=\frac{\beta\log R}{10J}.                            \tag{7}
\]

Distance-to-time-boundary weights w_u=exp(gamma*h*min(u,L-u)) change
by at most R over a block and its exterior endpoints since h*ell<10J/beta.
Let M=m(N+ell-1) and D=r/epsilon+8h|theta|ell^2. Coupled conditional
replacement preserves each strip marginal, including its admissible
support. The expected weighted Hamming distance obeys exactly

\[
W_{n+1}\leq(1-d)W_n+F/M,\quad
d=\ell(1-R\bar\alpha)/M,\quad F=2m\ell RD.
\]

Here 0<d<=1, and iteration leaves the boundary influence constant
2mRD/(1-R*baralpha). Thus the oscillation bound for a whole-row
observable a distance u from the temporal boundary decays as
exp(-gamma*h*u). The prefactor may depend on m,h; the rate does not.
It suffices to use strips long enough that all boundary pairs are
admissible; arbitrary long strips are available by (3).

For clarity, the vacuum extraction hypotheses of YM43 also hold:

\[
e^{-hr b_*}K_{hr}^{\rm row}(x,y)
\leq T_h^r(x,y)\leq e^{hr b_*}K_{hr}^{\rm row}(x,y).
\]

The strictly positive coarse bounds construct a unique bounded positive
vacuum by YM43's normalized iteration. T_h commutes with T_h^r, so its
positive vacuum is the same one; its eigenvalue lambda_h is its norm.
Conditioning long strips and then averaging an endpoint against this
vacuum gives decay of ground-stationary whole-row correlations, with
the above finite prefactor. If h_h is the normalized vacuum, take a
bounded row function f with Phi(h_h^2 f)=0 and the source g=f h_h.
Its even moments satisfy
<g,(T_h/lambda_h)^{2n}g><=C_g exp(-2n*gamma*h).
The even-moment norm extraction of YM43-T4 removes C_g. Such sources
are dense in the vacuum complement, so the bound extends to every complementary
source in the recognition completion, including real and iota parts.

Finally X=K_{h/2}e^{hB/2} gives T_h=X^*X and S_h=XX^*. The same
intertwining argument as YM45-T6 transfers the positive nonzero spectrum
and its complement ceiling. With Q_h the complement of S_h's vacuum,

\[
\boxed{\|Q_h(S_h/\lambda_h)^nQ_h\|\leq e^{-\gamma nh},
\qquad n\geq1.}                                           \tag{8}
\]

No coefficient cutoff enters (8). Independent free B-factors with the
same second-eigenvalue floor may also be included: T1 gives their gap
at least beta/2, whereas R<epsilon*J/2 implies gamma<beta/25<beta/2.

For every |theta|/beta<1/2400, J=5 gives
baralpha=1/2+1200|theta|/beta<1; choosing
R=(1+baralpha^{-1})/2 passes the strict gate. At or beyond the endpoint,

\[
J(\bar\alpha-1)=(J-5)^2/10
+(240|\theta|/\beta-1/10)J^2\geq0.
\]

This is a limit of this sufficient budget, not proof of gap loss.

## YM53-T4: norm refinement with noncentral heat

At fixed finite width let M_C=max_i tr C_i. A bond is linear in each
quaternion variable, so

\[
L v_i=\lambda_i v_i,\quad
\lambda_i=(\operatorname{tr}C_i+\operatorname{tr}C_{i+1})/4.
\]

For one endpoint, sum_a(D_a v_i)^2=(1-v_i^2)/4. Since
C_i<=tr(C_i)I, it follows that Gamma(v_i)<=lambda_i. Using YM52's
variance identity and K_s v_i=e^{-lambda_i s}v_i,

\[
K_s[(v_i-v_i(x))^2](x)
\leq (1-e^{-2\lambda_i s})+(1-e^{-\lambda_i s})^2
=2(1-e^{-\lambda_i s})\leq2\lambda_i s\leq M_C s.           \tag{9}
\]

To obtain the first term integrate
2 K_{s-u} Gamma(K_u v_i) from 0 to s, with
Gamma(K_u v_i)<=lambda_i exp(-2lambda_i u). For a symmetric Markov
kernel, Cauchy--Schwarz in the kernel and invariance turn (9) into
||[K_s,M_{v_i}]||<=sqrt(M_C*s). Summing gives

\[
\|[K_s,B]\|\leq b_*\sqrt{M_Cs}.                            \tag{10}
\]

Apply YM44-T2--T4's ordered insertion packets with (10), not a Taylor
series in the unbounded L. Moving one bounded insertion through heat
costs (10); absolute packet sums cost exp(b_*t). The same local
comparison constant is (1+1/sqrt(2))*sqrt(M_C)<sqrt(3M_C).
Telescoping over a partition with mesh delta therefore yields

\[
\boxed{\|S_{h_n}\cdots S_{h_1}-U(t)\|
\leq b_*t e^{b_*t}\sqrt{3M_C\delta},\quad \sum h_i=t.}       \tag{11}
\]

The ordered packet limit is independent of the partition; on coefficient
polynomials its infinitesimal generator is -L+B, with energy operator
A=L-B. For C_i=I, M_C=3, (11) recovers
YM44's 3*b_*t*exp(b_*t)*sqrt(delta). For b_*=0 the error is zero.
This estimate is at fixed finite profile; it does not claim uniform
refinement if tensor traces diverge with an additional cutoff.

## YM53-T5: the rate survives the fixed-width time limit

The packet construction gives a positive self-adjoint semigroup U(t).
Kernel comparison gives e^{-b_*t}K_t<=U(t)<=e^{b_*t}K_t. For
t>=9/beta the strictly positive bounded coarse kernel constructs a
unique positive vacuum. Commutation of the semigroup makes this a
common vacuum for every t>0; write lambda(t)=||U(t)||. In particular
lambda(t)>=exp(-b_*t).

For h=t/n, (11) gives delta_n->0 in operator norm, and (8) gives
the complement ratio q=exp(-gamma*t) for S_h^n. YM44-T5's normalized
gap gate applies once

\[
2\delta_n<e^{-b_*t}(1-q).
\]

The isolated rank-one vacuum projections converge, and the transferred
ceiling tends to q. Thus for Q the common limiting vacuum complement,

\[
\boxed{\|Q U(t)Q\|/\|U(t)\|\leq e^{-\gamma t},\quad t>0.} \tag{12}
\]

The rate is independent of finite m; the refinement prefactor need not
be. This statement constructs each finite-width time limit, not a new
anisotropic infinite-volume or joint-cutoff state. The old YM46--49
conclusions must be audited again before transferring them to this family.

## YM53-T6: interacting energy uses the ground source

Let A=L-B, Ah_0=e_0 h_0, h_0>0, Phi(h_0^2)=1, and H=A-e_0>=0.
Define pi(f)=Phi(h_0^2 f). For every real coefficient polynomial F,

\[
\boxed{\langle Fh_0,H Fh_0\rangle
=\pi\!\left(\sum_i\Gamma_{C_i}(F,F)\right)
\geq\gamma\,[\pi(F^2)-\pi(F)^2].}                          \tag{13}
\]

Here is the domain justification also for rank two. L is the
self-adjoint closure of its finite coefficient blocks, and B is bounded;
h_0 lies in D(L). Factor C_i as sum of b*b^T. For the corresponding
active native derivatives D_b, coefficient integration by parts controls
sum ||D_b h_0||^2=<h_0,Lh_0>. Polynomial multiplication preserves D(L):
apply the product rule to finite coefficient truncations of h_0, and
use bounded polynomial derivatives and this derivative norm bound to
pass in graph norm. This proof never uses C_i^{-1} or control of an
inactive derivative.

Expand <Fh_0,L(Fh_0)> by native integration by parts and subtract
<F^2h_0,Lh_0>. Cross terms cancel, leaving
Phi(h_0^2 sum Gamma_{C_i}(F,F)); the multiplication terms B and e_0
cancel using the ground equation. This proves the equality in (13).
Apply (12) to (F-pi(F))h_0 and differentiate the normalized semigroup
quadratic form at zero to get the inequality. For complex coefficients
replace squares by modulus squares and use the sesquilinear energy.

Thus the interaction changes the ground-source weight; it is incorrect
to substitute the free Phi energy or its exact free gap into (13).
No identification of this energy with physical joules, thermodynamic
heat/work, or a measured clock is made.

## YM53-T7: clock covariance, a rank-two cell and failure gates

For k>0, rescale C_i,beta,theta by k and h by 1/k. The actual transfers,
|theta|/beta, skeleton size and comparison margins are unchanged;
gamma scales by k. Equation (11) is unchanged when t and mesh scale by
1/k. Equation (13) scales linearly with k. These are clock-covariant
relations, not a physical calibration.

For C=diag(0,3/2,3/2), beta=3/2, the sufficient window is
|theta|<1/1600. In particular theta=1/4096, J=8, R=4/3 give

\[
\bar\alpha=5/8,\quad1-R\bar\alpha=1/6,\quad
\gamma=3\log(4/3)/160>0.00539403.
\]

This holds for every finite chain width and 0<h<=2/3, and in its
fixed-width time limit. Sitewise rotations of this tensor are permitted.
For rank one beta=0 is refused; YM52 has an explicit nonconstant
stationary even source. For C_eta=diag(eta,eta,3-2eta), the available
window beta/2400 and rate collapse as eta->0. Outside the sufficient
window the certificate abstains; it does not assert a gapless phase.

## Evidence and remaining obligations

The written all-content arguments above are separate from finite exact
controls in `ym53_anisotropic_interaction.py`. The controls include a
closed infinite-tail bound, integer and half-integer spin inequalities,
8-variable bond identities, outward rate/refinement enclosures and
independent finite fixtures for zero bridge normalizers and weighted
ground-source energy. YM45's bridge, tilt, incidence and joint-update
controls are reused without changing their engines. The finite fixtures
are not discretizations of the native compact carrier.

Source SHA-256 pins, fresh certificate comparison and tamper tests bind
the proof/code/tests/workflow. This is a written proof with reproducible
checks, not mechanical proof formalization or external expert review.
Physical protocol, clock, action/state and observer selection; actual
row-memory closure; the anisotropic joint/volume limit; NCG quantum
measure; four-dimensional QFT, asymptotic freedom, Clay and quantum
gravity remain open. The old isotropic certificates are unchanged.
