# A certified executable solver for a uniform LMD pulse marginal

Author: Monty. 30 October 2026.
Frozen parent: `ce60e932bb46142eafeef304c099c15ea89ad088`.

The previous packet derived an exact characteristic law but did not evaluate
its oscillatory source expectation. This packet implements that evaluation
with a deterministic, certified interval. It uses native rational arithmetic
for conditional mass averages and error bounds, and fresh NRQC phase-readout
queries for the oscillatory kernel. The certification includes all omitted
source modes; it is conditional on the same LMD/common-parent construction.
This is one vacuum marginal under a prescribed uniform pulse, not a full
field-state simulation, joint coframe solver or autonomous gravity model.

## 1. Retained source and operational output

For 0<h<1/2, let J=M_h(S1) in the LMD vacuum and chi_h(a)=E exp(i a J).
The [frozen characteristic theorem](../native-lmd-pulse-characteristic/THEOREM.md)
gives, for the existing common-parent conjugate reading v_0,

    Phi_t(k) = E_t exp(i k v_0)
             = d_h(k) chi_h(a_t(k)),
    d_h(k)   = exp(-h k^2/4),
    a_t(k)   = (2t/h) sinh(hk/2).                            (1)

The solver returns intervals for the real/imaginary parts and for

    p_cos,+ = (1+Re Phi_t(k))/2,
    p_sin,+ = (1+Im Phi_t(k))/2.                             (2)

These are legitimate binary measurement probabilities: after observing v_0,
output '+' with probability (1+cos(kv_0))/2 or (1+sin(kv_0))/2. Equivalently,
integrate these bounded positive functions against the existing parent POVM.
They are not P(v_0>0), nor samples or frequency estimates from an experiment.

For comparison, the rejected additive vacuum-label rule has

    Phi_naive,t(k) = d_h(k) chi_(2h)(tk).                    (3)

Both source expectations in (1) and (3) are evaluated with their own mass
law and the same explicitly audited algorithm. The interval comparison,
rather than just a difference between point estimates, decides separation.

## 2. Conditional averaging improves the continuum error bound

Let J_gamma be a positive source mass with E J_gamma=1 and finite second
moment I(gamma), 0<gamma<1. For any finite measurable partition G, set

    Jbar = E[J_gamma | G],
    delta_G^2 = E[(J_gamma-Jbar)^2]
              = I(gamma)-E[Jbar^2].                        (4)

Then, for every real a,

    |E exp(i a J_gamma)-E exp(i a Jbar)|
       <= min(2, (a^2/2) delta_G^2).                       (5)

Proof. The scalar Taylor formula with integral remainder gives

    exp(i a x) = exp(i a y) + i a exp(i a y)(x-y) + R,
    |R| <= (a^2/2)(x-y)^2.

Take y=Jbar. The linear term has zero expectation by the defining conditional
mean identity; it is integrable by L2. Taking expectations proves the
quadratic bound. The modulus-one bound supplies the alternative constant 2.
Orthogonality of conditional expectation proves (4). This is a direct
application of standard conditional averaging and Taylor remainder, not a
new general positivity or probability theorem.

In particular, the same improvement applies to the previous Fourier mass
cutoff J_N=E[J_gamma|F_N]: its old O(|a| delta_N) error can be replaced by
O(a^2 delta_N^2). The source's omitted modes are accounted for by an exact
second-moment remainder, without assuming that a finite state is the field.
The improvement is especially useful near zero pulse.

## 3. An exactly integrable finite source compression

Use the zero-mode Gaussian q~N(0,gamma). The nonzero-mode Wick mass is
independent of q and has mean one, so

    E[J_gamma | q] = exp(q-gamma/2).                        (6)

Partition the entire real q-axis into finitely many bins B_b=(l_b,r_b],
including both unbounded tails. For F_gamma the N(0,gamma) CDF,

    P_b = F_gamma(r_b)-F_gamma(l_b),
    A_b = E[J_gamma 1_(q in B_b)]
        = F_gamma(r_b-gamma)-F_gamma(l_b-gamma),
    m_b = A_b/P_b.                                         (7)

The shift identity follows by completing the square in
exp(q-gamma/2) times the Gaussian density. Every P_b is strictly positive.
No tail is discarded or silently renormalized. Also sum P_b=sum A_b=1.
Therefore the implementable approximation and its complete source budget are

    chi_G(a) = sum_b P_b exp(i a m_b),
    delta_G^2 = I(gamma)-sum_b A_b^2/P_b.                   (8)

The I(gamma) intervals are the frozen native log-Gamma calculations in the
predecessor's source_intervals.json, with its source hash pinned. They include
gamma=1/8,1/4,3/8,1/2,3/4. Reusing certified source values is explicit; the
NRQC calculator is not being credited with deriving continuum gamma moments.

On nested partitions, E[Jbar^2] increases and delta_G^2 decreases. Refining
only q cannot recover the nonzero modes: even a complete q observation leaves

    delta_zero_mode^2 = I(gamma)-exp(gamma).                 (9)

The implemented central interval is [-2,2], with 32,64 or 128 equal bins and
two tail bins. Refinement reduces its extra quantization error. Neither
central refinement nor changing a requested tolerance can remove (9).
To reduce the source part of this certificate, additional Fourier modes or
another justified mass approximation must be included. A bound exceeding a
requested tolerance is reported explicitly; it is not relabeled a success.

## 4. Exact calculator nodes and certified phase interpolation

The unmodified, pinned NRQC cyclotomic64 backend executes H|0>, followed
by j copies of p64, for j=0,...,16. X and Y expectations are respectively
cos(j pi/32) and sin(j pi/32), evaluated exactly in its cyclotomic field.
Nine requests use at most two qubits each. The 17 first-quadrant nodes extend
to a complete 64-node circle by exact signs and reflection.

An arbitrary phase a m_b is evaluated by cubic Lagrange interpolation at
four adjacent nodes. Write Delta=pi/32 and phase/Delta=j+u. The weights for
nodes j-1,j,j+1,j+2 are

    -u(u-1)(u-2)/6,
    (u+1)(u-1)(u-2)/2,
    -(u+1)u(u-2)/2,
    (u+1)u(u-1)/6.                                        (10)

The fourth derivatives of sine and cosine have absolute value at most one.
The interpolation error for each component is therefore at most

    Delta^4 |(u+1)u(u-1)(u-2)|/24 <= Delta^4/24.           (11)

The implementation verifies -0.01<=u<=1.01. On this slightly enlarged cell
the absolute polynomial product is less than one: with v=u-1/2, the product
is v^4-(5/2)v^2+9/16; inspection on |v|<=0.51 gives the bound. This handles
tiny interval crossings of a grid node without an unjustified branch choice.
Periodicity supplies nodes across 0 and 2pi. NRQC's 24-decimal outward
enclosures and all coefficient/phase uncertainties remain intervals.

These are classical interpolation weights, which can be negative. They are
not amplitudes of a prepared source state or probabilities of a quantum
mixture. NRQC supplies actual phase-readout computations to a native solver;
no unitary implementing the full density pulse has been compiled. A separate
native sine/cosine Taylor implementation checks every node and every assembled
discrete characteristic value independently of the NRQC phase engine.

## 5. Elementary-function arithmetic and total error

The new interval adapter uses the existing native F rational field and
192-bit outward rounding. Signed floor and ceiling are explicit. Square roots
use integer square-root bounds. Pi is enclosed by Machin's identity
16 atan(1/5)-4 atan(1/239), with convergent alternating remainder bounds.
Exponentials use argument halving, a positive Taylor series and a geometric
tail; negative arguments use interval reciprocals. Gaussian CDFs use the
integrated exponential series with a bounded alternating remainder, and
Gaussian normalization uses the certified pi/square root. Sine/cosine use
period reduction and absolute Taylor-tail bounds. The supported input range
is asserted before these routines run; arbitrary unbounded inputs are not
silently accepted.

Combining interval (8), phase error (11) and source error (5) gives a certified
enclosure of chi_gamma(a). Multiplication by the positive damping in (1) and
the affine map in (2) gives the final probability enclosure. It is intersected
with [0,1], using the analytic probability property. No floating-point number
enters a scientific bound. Floating time measurements are execution metadata.

The central value is not presented without its interval. The output reports
the actual probability halfwidth and separately flags whether it meets the
requested halfwidth. Narrower numerical arithmetic cannot hide a coarse
source-compression budget. The naive comparison uses gamma=2h, so it may have
a substantially larger source error than the actual model.

The discrete expression chi_G is a probability characteristic function of
the compressed mass as a function of a. Its composition with nonlinear
a_t(k), or the cubic approximation to that composition, is not asserted to
be the characteristic function of an independently prepared quantum state.
In particular, replacing J by a piecewise constant mass destroys the exact
zero-mode scaling used to derive (1). The solver approximates the expectation
in the already-proved formula; it does not re-derive that formula for a
different pulse. Each final interval encloses the original source probability.

## 6. Scope and next refinement

This closes an executable expectation gate for the declared uniform-vacuum
measurement: the code returns bounded-reading probabilities with a rigorous
continuum error budget, conditional on the frozen source theorem. The final
request and all intervals, phase requests, calculator outputs and checks are
retained in results/; EXECUTION.md reports which requested accuracies and
separations were achieved.

A general next refinement uses a finite sigma algebra G involving additional
Gaussian Fourier modes and still sets Jbar=E[J|G]. Formula (5) and the second
moment budget remain valid. Conditional source masses can be obtained by
Gaussian tilts followed by a certified spatial integral. The size and cost
of that representation must be measured; no efficient high-mode algorithm
or quantum advantage follows from this first one-mode partition.

Nonuniform pulses, non-vacuum/joint coframe laws, all-material-flow detector
covariance, fixed-h removal of geometric smoothing, pulse-evolved length
domains, autonomous gravitational dynamics and physical 3+1 recovery remain
separate obligations. The solver does not change the chosen physical source,
the prior infinite-mean-length result, or the external residual-gauge work.
