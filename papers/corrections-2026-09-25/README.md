# Corrections to the native-theory prerequisites

Date: 2026-09-25. These corrections supersede the identified statements at the
source revisions below. They repair proofs and evidence scope; they do not
certify the complete theory, RH, the Yang--Mills mass gap, or quantum gravity.

| Source examined | Revision |
|---|---|
| Recognition-Kernel-Framework main | 86198d29cbf30390059079f38675c952e506ea9c |
| RH-Framework main | fd104f46cc8c99c05f0e05b10f18c9332e8219de |
| Publications main | a3ff8e9f8e6afdea44c55557346e8b5fea0c5d22 |
| MP paper-source-floor-five-by-five-theorem | 686b5506c7c2fedd53641dbf550218b60c996a13 |
| Thermodynamic manuscript | arXiv:2603.20773v3 |

## 1. Tensor cut-tail mass and function completion

For z=r+iota t, D(z)=|r|+|t| is submultiplicative, not multiplicative.
Expanding zw and applying the real triangle inequality proves
D(zw)<=D(z)D(w). Summing over tensor indices proves
M(A tensor B)<=M(A)M(B).

T50 now uses that inequality with positive recognized floors and 0<=rho<1.
Each memory word is bounded by rho^k times the product floor. A word with
one memory factor attains its local ratio, so the worst sheet ratio still
is the maximum local ratio. Applying T01's action bound on disjoint sheets
preserves the energy contraction theorem.

Exact negative control: D=[[0,1+iota],[0,0]], B=(I+D)/8, f0=1, rho=1/2.
M(B)=1/2, but M(B tensor B)=7/32<1/4. The earlier fixtures included only one
mixed face, so they missed this failure.

T54's memory increment is bounded by Pi_n mu_(n+1), with equality for its
pure-axis calibration. Summable local masses still imply a mass-Cauchy
limit: split off a finite prefix and choose a tail sum s<=1/2; the tail
product is at most 1/(1-s)<=1+2s. The constant-mass pure-axis example still
separates a uniform finite gap from existence of the infinite product.

RH-Framework T01-E4 uses I1(zf)<=D(z)I1(f). Multiplication by a fixed scalar
is uniformly continuous, which suffices for completion. The unit-cell
example z=f=1+iota has left side 2 and product bound 4.

## 2. SPECTRAL 2, Corollary 16.2 and Remark 16.3

This subsection is the replacement statement for the archived SPECTRAL_2.pdf
at the Publications revision above, Git blob
549a48bd8e536fc85abfb0bf5f3f5772d1d415d0. The original PDF is retained as a
historical publication; its uncorrected Corollary 16.2 must not be used.

**Corrected upper spectral enclosure.** Let K be bounded and self-adjoint,
P an orthogonal projection with nonzero finite-dimensional range, and
mu=sup spectrum(PKP restricted to PX). If e>=||K-PKP||, then

\[
\sup\sigma(K)\le \max\{\mu,0\}+e.
\]

For P=0, use the bound e alone. Thus max(mu,0)+e<=1 is sufficient for
I-K>=0. It is sufficient, not necessary. The reference/congruence hypotheses
in the original corollary remain required when translating this to W.

**Proof.** Write x=Px+(I-P)x. Then
<x,PKPx>=<Px,KPx><=mu ||Px||^2<=max(mu,0)||x||^2.
Also <x,(K-PKP)x><=e||x||^2. Add the bounds and take the supremum over unit
vectors. This proves the claim without identifying the compressed spectrum
with the spectrum of its zero extension. When the complement is nontrivial,
the zero extension includes zero; when P=I, the displayed max is conservative.

**Counterexample to the old bound.** K=diag(-2,2), P=diag(1,0), R=I.
mu=-2 and e=2 give the old upper bound 0, but sup spectrum(K)=2 and
W=I-K=diag(3,-1) is not positive. The corrected bound is 2 and refuses
acceptance. P_1=P and P_n=I for n>=2 satisfy strong convergence, so this
is also a counterexample within a convergent projection family.

## 3. Information invariance: full response versus coarse observation

With W>0 and an involution Theta, the full antisymmetric correlation A_W
vanishes for every pair of observables iff W(Theta a)=W(a) for every state.
Reindexing the sum proves one direction; the indicator pair at a and Theta a
isolates W(a)-W(Theta a), proving the other.

For a declared linear observer K, the correct chain is
chi=1 iff A_W=0 => K(A_W)=Omega=0. Under Theorem A's chart hypotheses this
last condition is equivalent to delta I=0. Reversing the coarse arrow requires
ker(K) intersect admissible antisymmetric responses = {0}.
The two-orbit weight (3,1,1,3) produces (+2,-2); its sum observation is zero
but the microscopic response is nonzero. The source, abstract, ledger,
summary table and MICRO certificate now agree. The certificate's executable
status is PASS_FINITE_CHECKS; the general proof is written in the paper.

## 4. Twisted trace and Morphic Calculus

For W_k=O_k composed with W_(k-1), twisted additivity gives

\[
\operatorname{tr}_\omega(W_n)
=\sum_{k=1}^n\operatorname{tr}_\omega(O_k)
+\sum_{k=2}^n\omega(O_k,W_{k-1}).
\]

This follows by telescoping; the cocycle law makes the result independent
of parenthesization. The previous sum over every prefix/suffix cut double
counted. With additive arrows, t(a)=a^2/2, omega(b,a)=ba and three unit
arrows, the old formula gives 11/2 instead of 9/2. This counterexample tests
twisted additivity and the cocycle law, not every other overlay axiom.

For bounded represented generators and small h, Euler factors C_i=I+h D_i
give h^-1 log(C_1 C_2)=D_1+D_2+h([D_1,D_2]-D_1^2-D_2^2)/2+O(h^2).
The symmetric square-root Euler composition cancels the commutator term,
but retains -h(D_1^2+D_2^2)/2. Exponential symmetric splitting has the
second-order property claimed previously for the Euler factors.

A loop H_h=I+h^2[D_1,D_2]+O(h^3) satisfies (H_(T/n))^n -> I.
The nontrivial area scaling is (H_(T/sqrt(n)))^n -> exp(T^2[D_1,D_2]),
as follows by multiplying the local logarithm by n. Nonzero commutator
gives a nonconstant limiting family, not a guarantee at every finite T.

The torus equation (1+i k.Omega) U_k=0 kills every mode, including k=0.
For finite Markov chains, irreducibility alone does not make Shannon entropy
monotone: the positive matrix with both rows (3/4,1/4) decreases entropy from
the uniform distribution. Use negative relative entropy to a positive
stationary distribution; log-sum contraction proves monotonicity. Shannon
entropy works when the stationary distribution is uniform. A clocked rate
requires nonzero clock derivative. Spectral expansions include Jordan terms
unless diagonalizability is assumed.

Morphic Geometry's abstract now states the admitted-diamond holonomy result,
matching its proof; it does not assert confluence iff flatness. F00-I's
Taylor majorant starts with n^(-(1+2 delta)); multiplication by n^|h| with
|h|<=delta gives the asserted summable n^(-(1+delta)) bound.

## 5. Thermodynamic Pluecker test

The complete corrected source edition is in
[thermodynamic-response-corrections](../thermodynamic-response-corrections/).
For any estimated vectors u~,v~, their reconstructed bracket
B~=u~ v~^T-v~ u~^T has Pluecker residual zero identically. No violation can
be observed using that construction alone.

Let C instead contain six independently observed skew entries. With
b=beta(B), c=beta(C), ||c-b||<=eta from an independently justified error
contract, and P(B)=0, write e=c-b. The signed coordinate-pairing matrix J
has ||J||=1 and P(B)=b^T J b/2. Therefore

\[
|P(C)|\le\|b\|\|e\|+\tfrac12\|e\|^2
\le\|c\|\eta+\tfrac32\eta^2.
\]

A violation rejects the conjunction of model and measurement assumptions.
The error budget must not be fitted from the residual being tested. The
original gradient reconstruction estimate remains valid, but is a separate
positive control. No experiment or nonlinear gravity claim is certified here.

## 6. MP normalization and current status

The active-band derivative gap tends to zero at the origin. A positive
unweighted infimum on (0,7) was therefore wrong. The actual driver separately
bounds ((1-eta)g'-|rho'|)/x on (0,h] and (1-eta)g'-|rho'| on [h,7].
The theorem and certificate metadata now state those domains. Since g'>=0,
replacing eta by smaller theta improves the conservative derivative lower
bound (1-theta)g'. The source-floor shift argument therefore survives with
the two separately normalized bounds; no fictitious global constant is used.

Earlier MP route ledgers are explicitly historical. The latest source reports
fixed-shift five-matrix closure through its companion certificate; the actual
unshifted K0 sign remains open. This correction does not independently
revalidate that companion's full numerical certificate. Publications now
states that provenance and does not claim RH.

## Evidence layers and reproduction

- Written proofs: the corrected arguments above and the corresponding source.
- Exact finite checks: certificate.py and its pinned report; RKF's correction
  regressions and T50/T54 tests; RH's strict multiplier-bound control.
- Formal proof assistant: no new Lean certification is claimed.
- Provenance: source revisions and file hashes identify the checked objects.
- External peer review: not performed by this correction pass.

Run `python certificate.py --check` and `python -m unittest test_certificate -v`
from this directory. Verification details and affected source hashes are in
VERIFICATION.json. Original manuscripts, archived book PDFs and old source
snapshots remain versioned historical records; this correction record controls
the identified statements in those editions.
