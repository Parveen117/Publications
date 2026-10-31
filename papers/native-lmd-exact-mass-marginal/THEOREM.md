# Exact total-mass law for the remaining uniform LMD pulse marginal

Author: Monty. 31 October 2026.
Frozen source parent: `ace56a0a3faa21708e4a6779c8e4a8f02dab70e0`.

The remaining probability-accuracy case can be approached through an exact
law of the total mass instead of a many-mode conditional partition. This
packet adds one explicit external analytic dependency: Rémy's proved
Fyodorov--Bouchaud law. It audits that law against the *existing* LMD
covariance, Wick subtraction and circle normalization before applying it.
It does not derive a new universal source law or replace the quantum model.

## 1. External input and its source match

Guillaume Rémy, *The Fyodorov-Bouchaud formula and Liouville conformal field
theory*, arXiv:1710.06897v3, Theorem 1.1 and equation (1.9), proves a full
distribution identity for normalized circular GMC. In his convention the
field X has boundary covariance -2 log|exp(i theta)-exp(i theta')|, and the
Wick exponential has coefficient beta/2. For 0<beta<2 its normalized mass is
distributed as E^(-beta^2/4)/Gamma(1-beta^2/4), with E~Exp(1).
This theorem, not a fit to finitely many moments, is the imported result.
The precise statement is on PDF page 3; conventions are on page 2:
https://arxiv.org/pdf/1710.06897v3.

Here is the source-specific applicability calculation. The existing LMD
field splits into its independent zero mode and nonzero modes,

    s(z)=q+eta(z),  q~N(0,h),
    Cov(eta(z),eta(w))=-h log|exp(2 pi i z)-exp(2 pi i w)|.

Set beta=sqrt(2h), theta=2 pi z and eta=(beta/2)X. Then

    Cov((beta/2)X)=(beta^2/4)(-2 log|chord|)
                 =-h log|chord|.

Normalized source length dz becomes dtheta/(2 pi). The nonzero-mode Wick
subtraction is (h/2) sum_(n<=N)1/n, while the independent zero-mode factor
is exp(q-h/2). The LMD-T3 Fourier/mollifier identification supplies the same
Wick measure; no extra unknown smooth covariance term is inserted. Thus
beta^2/4=h/2 and, for the existing total mass J=M_h(S1),

    J =_law exp(q-h/2) E^(-h/2) / Gamma(1-h/2),
    q independent of E,  E~Exp(1).                         (1)

The model's actual measurement uses 0<h<1/2, comfortably inside both this
subcritical range and LMD's L2 range. Formula (1) is an auxiliary probability
representation of this one total-mass law. It does not identify E with an
individual Fourier amplitude, create a field-state preparation, or provide
an operator/dynamics intertwiner. The independent zero mode is essential.

For p<2/h, direct integration of the auxiliary law gives

    E J^p = exp(h(p^2-p)/2)
            Gamma(1-ph/2)/Gamma(1-h/2)^p.                 (2)

In particular E J=1 and

    E J^2=exp(h) Gamma(1-h)/Gamma(1-h/2)^2=I(h),

agreeing with the independently frozen source moment. These checks verify
normalization, but do not themselves prove equality in law. For example,
a two-point variable taking values 0 and I(h), with probability 1/I(h) at
I(h), shares both moments and generally has another characteristic function.

## 2. The actual pulse target is unchanged

Retain the frozen common-parent characteristic formula

    Phi_t(k)=d_h(k) E exp(i a_t(k) J),
    d_h(k)=exp(-h k^2/4),
    a_t(k)=(2t/h) sinh(hk/2),
    p_cos,+=(1+Re Phi_t(k))/2,
    p_sin,+=(1+Im Phi_t(k))/2.                            (3)

Since the expectation in (3) depends only on the distribution of J, (1)
can be used without constructing a joint law of the field and its conjugate
detector. The inherited derivation of (3) still carries that operator claim.
The auxiliary variable is not substituted for the source inside an arbitrary
nonuniform pulse or joint coframe observable.

The all-amplitudes-only obstruction from the predecessor remains true for
its specified conditional sigma algebra. Formula (1) uses the complete
distribution of the collective mass, including its phase correlations. It
does not independently erase the Fourier phases or contradict that bound.
This changes the expectation algorithm and its theorem dependencies, not the
physical source or target. No 8-mode solver is claimed here.

## 3. Positive conditional cells, including the singular endpoint

Put alpha=h/2, s=1-alpha, C=Gamma(s), and

    Q(q)=exp(q-h/2),  W(E)=E^(-alpha)/C.

Both Q and W have mean one. Gaussian q bins retain the inherited tilted
first masses a_i and probabilities p_i. For an E cell (u_j,v_j],

    p_j=exp(-u_j)-exp(-v_j),
    b_j=E[W(E) 1_(u_j<E<=v_j)]
       = [gamma(s,v_j)-gamma(s,u_j)]/Gamma(s),            (4)

where gamma(s,x)=integral_0^x exp(-v) v^(s-1) dv is the lower incomplete
Gamma integral. The first cell starts exactly at zero; the last reaches
infinity. The divergence of W at zero is integrable, and its squared moment
is finite because h<1. Neither endpoint is discarded or replaced by a
finite cutoff followed by renormalization.

For the product partition, define

    P_ij=p_i p_j,  A_ij=a_i b_j,  m_ij=A_ij/P_ij,
    S_q=sum_i a_i^2/p_i,  S_E=sum_j b_j^2/p_j.

These are conditional means of the auxiliary mass, whose law is exactly
J's. The complete residual is

    delta_G^2=I(h)-S_q S_E
       = [exp(h)-S_q] E[W^2] + S_q [E[W^2]-S_E],        (5)
    E[W^2]=Gamma(1-h)/Gamma(1-h/2)^2.

There is no omitted Fourier-mode term in (5), conditional on the external
law and its source match. The two terms are q and E quantization, including
variation inside their retained tail cells. Fixed undivided tails still
leave a nonzero budget; simply shrinking interior cells is not claimed to
give arbitrary accuracy. Expanding the domains and refining a generating
partition gives L2 convergence by the conditional-expectation martingale.

The old Taylor cancellation applies unchanged:

    |E exp(i a J)-sum_(i,j) P_ij exp(i a m_ij)|
          <= min(2, (a^2/2) delta_G^2).                  (6)

Its proof needs only the full auxiliary second moment and exact conditional
cell means, not a physical measurement of q and E. Combining (6), the
existing phase interpolation remainder, its calculator-node enclosures,
damping in (3) and the factor 1/2 certifies each output probability.

## 4. Native Gamma enclosures without floating special functions

For s>0 and x>0, expanding exp(x-v) under the integral and integrating
the positive terms gives

    gamma(s,x)=x^s exp(-x)/s
       * sum_(n>=0) x^n/[(s+1)...(s+n)].                 (7)

The empty product is one. If t_N is the last retained series term, its
remaining sum is bounded by

    t_N * x/(s+N+1) / [1-x/(s+N+2)]                     (8)

when the denominator is positive. The implementation uses N=384 and x<=64.
All terms, powers, exponentials and remainders are directed rational
intervals. Dyadic s values use repeated certified square roots followed by
integer powers, so x^s needs no floating logarithm.

The complete Gamma normalizer is enclosed using its integral to 64 and
the positive tail bound, for 0<s<=1,

    0 < integral_64^infinity exp(-v)v^(s-1)dv
          <= exp(-64)64^(s-1).                          (9)

For the independent recurrence checks the implementation also accepts
1<s<=2. Writing v=64+u and using
(1+u/64)^(s-1)<=exp((s-1)u/64) gives the bound
exp(-64)64^(s-1)/[1-(s-1)/64]. Thus those checks integrate Gamma(s+1)
independently; they do not compute it by the recurrence being tested.

The s=1 identity and Gamma(s+1)=s Gamma(s), tested by separately computed
integrals, provide analytic checks. The first and second auxiliary mass
moments are compared with their old independent source enclosures. Gamma
CDF endpoints 0 and infinity are assigned their analytic values exactly.
Known CDF/bin ranges [0,1] may be intersected with interval bounds, but no
lost tail is silently normalized away.

## 5. Frozen execution scope

The new request fixes the remaining case h=3/8,t=1/2,k=2, its time reverse,
a previously certified h=1/8,t=1/4,k=2 cross-check and a zero-pulse baseline.
The original requested probability halfwidth remains 1/200. Two nested
partitions use q bins 32 and 64 on [-2,2] plus both tails. E boundaries cover
each octave [2^j,2^(j+1)] for -16<=j<4 using 2 then 4 equal subdivisions,
plus the cell from zero to 2^-16 and the tail from 16 to infinity.

All nine original NRQC phase results are reused with their pinned hashes;
there are no fresh calculator or hardware requests. Node uses the existing
native rational engine. Python 3.12 orchestrates and hashes only. The exact
law, source audit, algorithms and request are committed before evaluation.
All achieved widths and failures are retained. A passing certificate closes
the selected accuracy gate under the stated dependencies, not full quantum
gravity or an independently implemented whole-field circuit.

Nonuniform sources, joint geometry and current laws, detector calibration,
fixed-h removal of geometric smoothing, length domains, autonomous gravity
and physical 3+1 constraints remain separate obligations.
