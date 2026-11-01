# YM-57: neighbour turns restore control at zero local curvature

Monty Dabas. 4 October 2026. **Python 3.12 only.**

YM56's second local response direction disappears at lambda=0. Can
neighbours restore a gap? The answer depends on the interaction. The
existing scalar overlap potential preserves the missing-direction
invariants. A specified correlated-turn kinetic interaction instead
recovers every missing local rotation through a finite native word.

For the new finite open chain of N>=2 sites, put
\[
 A_i=D_{1,i},\qquad B_i=D_{2,i}-D_{2,i+1},\qquad
 H_{N,k}=-\frac12\sum_{i=1}^N A_i^2
          -\frac{k}{2}\sum_{i=1}^{N-1}B_i^2,\qquad k>0.
\]
On the product of the existing compact native carriers, its constants
are the unique zero-energy sector and its centered gap obeys
\[
 \boxed{\operatorname{gap}(H_{N,k})\ge
 \min\{1/(63\pi^2),\ k/(42\pi^2)\}
 \ge \min\{1/630,\ k/420\}>0.}                            \tag{1}
\]
This is uniform in N. All *on-site YM56 response curvatures* may be zero.
The control is supplied by new mixed turn words across sites. This is a
new selected kinetic model, not an assertion that YM53's multiplicative
interaction has changed its principal derivative directions.

## Sources, representation and scope

Reuse YM50's compact quaternion carrier, positive reference Phi_Q and
coefficient/recognition completions, YM51's native turn derivatives,
YM52's square energies, YM53-T6's domain argument for ground transforms,
and YM56's exact stationary witness. Work at fixed finite N with product
reference Phi_N=Phi_Q tensor ... tensor Phi_Q. Different sites commute;
[D_{1,i},D_{2,i}]=D_{3,i}. Every real D is skew. No new one-site algebra,
scalar field, quaternion multiplication or single-site operator engine
is introduced. Tensoring and labelled correlated turns are explicit choices.

The proofs use established unitary telescoping, compact group averaging
and product-variance tensorization, derived below for this carrier. They
do not claim those general methods as inventions. The contribution is
their quantitative application to this native neighbour-turn model and
the explicit obstruction for the old potential-only model. Finite exact
checks support the written proof; they are not a formal proof-assistant
verification or independent external mathematical review.

## YM57-T1: scalar bonds do not restore a lost local derivative

In the original model let one site i have lambda_i=0 and set
\[
 w_i=q_{2,i}^2+q_{3,i}^2-\tfrac12\sum_{a=0}^3q_{a,i}^2.
\]
YM56 proves D_{1,i}w_i=0, Phi_Q(w_i)=0 and Phi_Q(w_i^2)=1/12.
All other sites' derivatives also annihilate w_i. Thus multiplication
by w_i commutes with the free heat. It commutes with every scalar
potential multiplier, including theta sum(q_i dot q_{i+1}). Consequently
it commutes with each square-sourced transfer and the bounded-potential
semigroup. This argument does not assume the rank-two mixing hypotheses
that fail at this site.

If that original operator had a strictly positive normalized ground
vector h, then (w_i-c)h, c=Phi_N(h^2 w_i), would be a nonzero orthogonal
ground vector. Indeed w_i is bounded, smooth and constant along every
active direction, so it preserves the heat-generator domain and the
commutation extends to the bounded potential. Its weighted variance is
positive because h>0 almost everywhere and w_i is not constant. Thus
the model cannot have a *unique strictly positive* ground vector there.
This is not a claim that every gap above a larger kernel vanishes.

The same mechanism is visible in brackets: [D,M_V]=M_{DV}; it produces
another multiplier, not a missing transverse derivative. Conditional
on the existence of a positive ground h at small nonzero lambda_i, the
ground-transform identity gives the explicit upper witness
\[
 \Delta\leq
 \frac{\lambda_i^2\,\pi_h((D_{2,i}w_i)^2)}{2\operatorname{Var}_{\pi_h}(w_i)}
 \leq\frac{\lambda_i^2}{8\operatorname{Var}_{\pi_h}(w_i)}.     \tag{2}
\]
Here pi_h(F)=Phi_N(h^2F). The last inequality follows from
|D_2 w|<=1/2 on the unit quaternion: its expression is a signed sum of
two products of disjoint coordinates, bounded by half the squared norm.
One cannot replace the weighted variance by 1/12 without a measure
comparison. Equation (2) rules out a uniform positive gap as lambda_i
vanishes only if that weighted variance has a positive uniform floor.

## YM57-T2: construct the correlated-turn square source

The derivative B_i generates opposite second-axis turns at two adjacent
sites. It is skew and preserves Phi_N. Let W=N+k(N-1). Give each
on-site label A_i weight 1/W and each bond label B_i weight k/W, and
pair each label with its negative with equal sign counts. At small
epsilon, use the actual operator turns Exp(+/-epsilon sqrt(W) D).
The average S_epsilon is positive, unital, symmetric and Phi_N-preserving.
On each fixed finite tensor coefficient space,
\[
 S_\epsilon=I-\epsilon^2 H_{N,k}+O(\epsilon^4),\qquad
 (S_{\sqrt{t/n}})^n\longrightarrow e^{-tH_{N,k}}.            \tag{3}
\]
Proof: the odd powers cancel pairwise; the second moment is
epsilon^2[sum A_i^2+k sum B_i^2]/2. At fixed N and coefficient degrees,
every generator is a bounded finite matrix and the Taylor remainder is
bounded there, so telescoping proves (3). Positivity and contraction
extend the compatible limit to the inherited product completion.
Equivalently use its orthogonal finite spin blocks; no infinite-degree
operator-norm Taylor estimate is claimed. This is the declared clock
normalization; it will matter in T6.

Its closed energy form is
\[
 \mathcal E(\psi)=\frac12\sum_i\|A_i\psi\|^2
                         +\frac{k}{2}\sum_i\|B_i\psi\|^2. \tag{4}
\]
On each finite coefficient block the operator is self-adjoint and
nonnegative in the native product pairing. The block closures supply
the positive self-adjoint generator on the completion, consistently
with (3).

## YM57-T3: the neighbour echo generates a missing local turn

Fix i and an adjacent site j. Orient B=D_{2,i}-D_{2,j} toward i and put
U=Exp(pi D_{1,j}). Quaternion multiplication, or the rotation of the
two native bracket directions, gives
\[
 U B U^{-1}=D_{2,i}+D_{2,j}.
\]
These two second-axis combinations commute. Hence, for every real s,
\[
 \boxed{e^{sD_{2,i}}
       =e^{sB/2}Ue^{sB/2}U^{-1}.}                          \tag{5}
\]
This is an exact finite-word identity, with no small-commutator limit.
The finite half-turn transfers access to a direction; it does not
change the operator algebra or identify frame brackets with NT curvature.

For a skew generator D and psi in its domain,
||(e^{tD}-I)psi||<=|t| ||Dpsi||, by integrating its unitary orbit.
For unitary factors V_r, expand their product minus I with left prefixes
to obtain ||(product V_r-I)psi||<=sum ||(V_r-I)psi||. Applying this to
(5) proves
\[
 ||(e^{sD_{2,i}}-I)\psi||
       \leq |s|\,||B\psi||+2\pi\,||D_{1,j}\psi||.         \tag{6}
\]
All factors act on psi in the bound; no uncontrolled derivative of a
transported state has been inserted.

## YM57-T4: compact averaging turns the echo into a uniform gap

Every unit quaternion has a first/second/first-axis Euler factorization
with |alpha|,|gamma|<=2pi and |beta|<=pi. To see this directly, write
q=z_1+z_2 e_2 in the e_1 complex plane. Choose phases phi_1,phi_2 in
[-pi,pi] and |z_1|=cos(beta/2), |z_2|=sin(beta/2) with beta in [0,pi].
Then alpha=phi_1+phi_2, gamma=phi_1-phi_2 give the factorization; when
one component vanishes choose its phase arbitrarily. The signs of the
point-action convention only reverse these bounded angles.

Use (6) for the middle factor. The action T_g on site i obeys
\[
 ||(T_g-I)\psi||\leq
 4\pi||A_i\psi||+\pi||B\psi||+2\pi||A_j\psi||.
\]
Cauchy--Schwarz gives the upper bound
21pi^2(||A_i psi||^2+||B psi||^2+||A_j psi||^2) for its square.
The constant 21 is conservative: (4,1,2) has squared norm 21.

Let P_i be averaging over site i using the existing invariant reference.
Invariance and unitarity give
\[
 ||(I-P_i)\psi||^2
  =\tfrac12\int ||(T_g-I)\psi||^2\,d\Phi_Q(g).             \tag{7}
\]
The commuting orthogonal projections P_i satisfy
I-product P_i <= sum(I-P_i): decompose into their joint zero/one
eigenspaces. Consequently total variance is at most the sum in (7).

Choose j(i)=i+1 for i<N and j(N)=N-1. Each site is selected as a
neighbour at most twice; each undirected edge is used at most twice.
Summing (7) gives
\[
 \operatorname{Var}_{\Phi_N}(\psi)
 \leq\frac{63\pi^2}{2}\sum_i||A_i\psi||^2
                   +21\pi^2\sum_i||B_i\psi||^2.           \tag{8}
\]
Comparison with (4) proves (1) on the polynomial core. Closure extends
it to the full form domain. It forces the kernel to be constants and
gives ||e^{-tH}(I-P_0)||<=e^{-gamma t}, where P_0=product P_i.
The inequality pi^2<10 gives the stated rational lower bound. No
finite spin truncation is the proof of this uniform statement.

## YM57-T5: local rank loss and the constitutive family

Add any fixed on-site YM56 second directions:
\[
 H_{N,k,\lambda}=H_{N,k}
                       -\tfrac12\sum_i\lambda_i^2D_{2,i}^2.
\]
The additional form is nonnegative, fixes constants and preserves the
same finite coefficient sectors. Therefore (1) still holds, even if
every lambda_i is zero. One may choose the cubic potentials of YM56
at all sites and freeze their jets at the origin. The on-site selected
curvatures are lambda_i R/2; the cross-site turn protocol is an explicit
additional selection, not derived from those zero local curvatures.

The simultaneous first-axis half-turn conjugation changes all B_i to
-B_i and preserves all squares. Antipodal parity is also preserved.
Thus the gap does not require breaking those symmetries.

The original scalar bond model at lambda_i=0 and this kinetic bond model
must remain distinct. Their equilibrium reference, ground-source weight
and dynamics are not interchangeable. Adding the original extensive
potential requires a further uniform interacting argument. YM53--55's
site-product heat hypotheses do not apply to this correlated heat merely
because both models use the same compact variables.

## YM57-T6: failure controls and the next physical gate

At N=1 with lambda=0, or at k=0 with every lambda=0, YM56's w_i remains
stationary. An isolated vertex likewise defeats the stated proof. The
gap above constants in those cases is zero, despite positive excited
eigenvalues above the enlarged kernel.

The clock is essential. If turns have fixed amplitude and weights
normalized by W, the limiting generator is H/W. For the centered local
w_i, exact integration by parts gives the Rayleigh quotient
\[
 \frac{\mathcal E(w_i)}{\operatorname{Var}(w_i)}
                         =\frac{k\deg(i)}2
\]
when all lambda=0. Thus the endpoint witness in the normalized clock
has quotient k/(2W), which tends to zero as N grows. A uniform statement
cannot ignore the distinction between an extensive-rate clock and a
single fixed-rate update for the entire chain.

For k=1, (1) supplies the dimensionless bound 1/630 for every N>=2.
The normalization and edge strength are specified, not physically
calibrated. This chapter establishes a finite-chain uniform energy gap
with correlated native transport. It does not construct its infinite-volume
history limit, add YM55's scalar interaction, remove spatial lattice
spacing, implement gauge constraints or identify a four-dimensional
Yang--Mills mass. The next mathematical task is the correlated model's
consistent local-history limit and its physical action/clock adapter.

## Reproduction

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym57_neighbour_turn_gap.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym57_neighbour_turn_gap.py -v
```

The certificate checks exact quaternion echoes, two-site polynomial
commutation obstructions, the actual correlated square generator,
finite tensor energy inequalities, sharp rank-loss controls, graph
incidence budgets and the clock-dependent Rayleigh witnesses. All run
on Python 3.12 with rational arithmetic. Source pins preserve YM56 and
its dependencies. The all-N and all-content statements rest on T1--T6.
