# First spatial-mode conditioning for the certified LMD pulse solver

Author: Monty. 31 October 2026.
Frozen parent: `13389614767542792097c189dcf6e11747806b42`.

This packet advances the four unresolved accuracy cases of the preceding
uniform-vacuum pulse solver. It conditions the same source on its first
nonzero sine/cosine pair, integrates the pair's orientation exactly, and
constructs a finite positive source partition including every tail. The
physical source, pulse, detector and probability target are unchanged.
All continuum statements remain conditional on the frozen LMD/common-parent
construction. This is a source-expectation refinement, not gravitational
Hamiltonian dynamics or a 3+1 theory.

## 1. Reuse the actual characteristic law

For the total LMD Wick mass J_gamma, 0<gamma<1, retain

    E J_gamma = 1,
    I(gamma) = E J_gamma^2
             = exp(gamma) Gamma(1-gamma)/Gamma(1-gamma/2)^2.

For the actual uniform-vacuum pulse, gamma=h, 0<h<1/2, the measured
characteristic and binary postprocessing probabilities are

    Phi_t(k) = exp(-h k^2/4) E exp(i a J_h),
    a = (2t/h) sinh(hk/2),
    p_cos,+ = (1+Re Phi_t(k))/2,
    p_sin,+ = (1+Im Phi_t(k))/2.                         (1)

The frozen solver proved the conditional-mean bound: for a finite sigma
algebra G, Jbar=E[J|G] and delta_G^2=I(gamma)-E[Jbar^2],

    |E exp(i a J)-E exp(i a Jbar)|
        <= min(2, (a^2/2) delta_G^2).                   (2)

Its q-only approximation leaves the budget I(gamma)-exp(gamma), even
before finite-bin error. We reuse (2); we do not claim a new general
conditional-expectation theorem or a lower bound on the actual characteristic
error. A large bound means this certificate cannot resolve the requested
tolerance, not that the true error must be large.

## 2. A first-mode source with no orientation quadrature

In LMD's orthonormal real Fourier convention, the zero coefficient q has
variance gamma and each first sine/cosine coefficient has variance gamma/2.
Writing the first-mode field instead in the unnormalized cosine/sine basis,

    s_1(z)=q+X cos(2 pi z)+Y sin(2 pi z),
    q,X,Y independent N(0,gamma),

is essential: confusing these two basis normalizations gives the wrong
source law. The diagonal variance of s_1 is 2 gamma. The cutoff Wick mass is

    J_1 = exp(q-gamma) (1/(2 pi))
          integral_0^(2 pi) exp(X cos(theta)+Y sin(theta)) dtheta
        = exp(q-gamma) B_0(R),
    R=sqrt(X^2+Y^2),
    B_0(r)=sum_(n>=0) (r^2/4)^n/(n!)^2.                (3)

B_0 denotes the standard modified Bessel function I_0 here, to distinguish
it from the already named mass moment I(gamma). Its series and circle
integral are standard identities, credited in SOURCE_MAP.md. The last
equality follows by translating theta by the polar angle of (X,Y).

The frozen Wick-martingale identity gives E[J|q,X,Y]=J_1. Formula (3)
is measurable in (q,R); hence the tower property also gives

    E[J|q,R]=J_1.                                       (4)

Thus removing the first pair's orientation loses no conditional mass
information for this uniform total-mass target. This does not say the full
field or its nonuniform measurements are independent of that orientation.
No separate observer measurement of three exact Gaussian coordinates is
being granted: conditioning is the classical expectation solver's internal
representation, not additional physical detector records.

Set S=R^2/(2 gamma). Polar integration of the two-dimensional Gaussian gives
S~Exp(1), independent of q. Therefore

    J_1=Q(q) W_gamma(S),
    Q(q)=exp(q-gamma/2),
    W_gamma(s)=exp(-gamma/2) B_0(sqrt(2 gamma s)).         (5)

Both factors have mean one. For W, introduce two independent uniform unit
vectors u,v in the plane and the Gaussian vector Z=(X,Y). Positive Fubini
and the Gaussian exponential identity yield

    E W^2 = exp(-gamma) E_(u,v) exp((gamma/2)|u+v|^2)
           = E_(u,v) exp(gamma u.v) = B_0(gamma).

Consequently

    E J_1^2 = exp(gamma) B_0(gamma),
    delta_1^2 = I(gamma)-exp(gamma) B_0(gamma).           (6)

The gain over q-only conditioning is exactly

    delta_0^2-delta_1^2
       = exp(gamma)[B_0(gamma)-1] > 0, gamma>0.          (7)

This agrees with the N=1 constant Laurent coefficient already retained by
the predecessor. The advance is an executable conditional representation of
the oscillatory expectation, rather than a new formula for that moment.

## 3. Exact positive cell masses from a Poisson/Gamma series

Partition the full q-axis into Gaussian bins B_i=(l_i,r_i], including both
unbounded tails. For the N(0,gamma) CDF F_gamma, use the inherited identities

    p_i=F_gamma(r_i)-F_gamma(l_i),
    a_i=F_gamma(r_i-gamma)-F_gamma(l_i-gamma).

Partition S in [0,infinity) into C_j=(u_j,v_j], including the last unbounded
tail. Boundary choices have probability zero. Define

    p_j=exp(-u_j)-exp(-v_j),
    b_j=integral_(u_j)^(v_j) exp(-s) W_gamma(s) ds.       (8)

With lambda=gamma/2 and

    T_n(x)=exp(-x) sum_(m=0)^n x^m/m!,
    T_n(infinity)=0,

positive termwise integration of (5) gives

    b_j=exp(-lambda) sum_(n>=0) lambda^n/n!
                      [T_n(u_j)-T_n(v_j)].             (9)

The bracket is the probability of a Gamma(n+1,1) random variable in C_j;
it lies in [0,1]. Thus this is a positive Poisson-weighted mixture of Gamma
bin probabilities. Formula (9) follows directly from the series and the
integer-shape incomplete-gamma integral; no noncentral special-function
library or floating quadrature is needed.

After retaining n=0,...,N, a uniform omitted positive mass bound is

    epsilon_N <= exp(-lambda) lambda^(N+1)/(N+1)!
                  /[1-lambda/(N+2)],                   (10)

provided lambda/(N+2)<1. Each b_j is enclosed by its truncated interval sum
and this one-sided remainder. Gaussian, radial and source-mass totals satisfy

    sum_i p_i=sum_i a_i=sum_j p_j=sum_j b_j=1.

No tail is dropped, and there is no renormalization of retained cells.
The upper endpoint infinity is treated analytically, not as a large finite
substitute. The integral representation also proves b_j>0 for every cell.

On the product partition G, independence and (4) imply

    P_ij=p_i p_j,
    A_ij=E[J 1_(B_i x C_j)]=a_i b_j,
    m_ij=A_ij/P_ij,
    E[Jbar^2]=(sum_i a_i^2/p_i)(sum_j b_j^2/p_j),
    delta_G^2=I(gamma)-E[Jbar^2].                        (11)

The two sums in (11) converge respectively to exp(gamma) and B_0(gamma)
under an exhausting nested partition. Finite Gaussian refinement and radial
refinement both reduce delta_G^2 by conditional-expectation orthogonality.
Retaining a fixed central window and one undivided tail leaves extra tail
quantization error; the implementation does not identify its finite
partition with the ideal limit (6).

The squared-L2 budget factors exactly into distinct obligations:

    delta_G^2 = [I-exp(gamma) B_0(gamma)]
                + [exp(gamma)-S_q] B_0(gamma)
                + S_q [B_0(gamma)-S_S],                 (12)

where S_q=sum a_i^2/p_i and S_S=sum b_j^2/p_j. Each bracket is nonnegative.
The first term is the omitted higher spatial modes; the other terms are
finite q and radial bin errors. This is a budget decomposition, not three
independent measurements of the unknown source.

## 4. Error transfer and execution contract

The computed expectation is

    chi_G(a)=sum_(i,j) P_ij exp(i a m_ij).

Insert the complete (11) into (2), then multiply by exp(-h k^2/4) and
divide by two to obtain an error bound for each probability in (1).
Retain the previous exact NRQC phase nodes, their decimal enclosures and
their certified cubic interpolation remainder. They are reused data from
nine earlier calculator requests, not fresh calls or a prepared field state.
The phase interpolation is classical postprocessing, not a quantum mixture.

All scientific arithmetic uses the frozen native rational interval adapter.
The new Poisson/Gamma recurrence and Bessel series have explicit positive
remainders. Every cell is retained. A sparse deterministic selection of
cell phases is checked against the separately implemented native Taylor
trig routine; this does not replace the analytic interpolation bound.
The full assembled expectation is not advertised as independently replayed
by a second algorithm in this run.

The request freezes the four previously wide cases, k=2, a requested
halfwidth 1/200, and two nested product partitions: q has 64 then 128 bins
on [-2,2] plus both tails; S has 16 then 32 bins on [0,16] plus its tail.
It is committed before these new probability evaluations. A result is
certified only if its complete interval meets the tolerance. An unmet
tolerance is retained explicitly and is not a failed arithmetic check.
The earlier eight passing cases remain predecessor results, with their
original probabilities and provenance. They are not counted as new runs.

## 5. Scope of the next gate

If the ideal first-mode budget in (6), after the transfer factor from (2),
still exceeds 1/200, finer cells alone cannot make this *variance-based
certificate* uniformly small enough. More source information or a sharper
error theorem is then required. This is not a lower bound on the unknown
actual oscillatory error and not a physical uncertainty or no-go theorem.

For two nonzero modes, simultaneous spatial translation removes one common
orientation, but leaves a relative angle alpha_2-2 alpha_1. Treating both
orientations as disposable would generally change the total mass. A next
solver must retain that relative phase, or justify another conditional
compression. The attractive rotational simplification here is specific to
the first-pair uniform-mass target.

Nonuniform/joint coframe distributions, general states, detector calibration
covariance, fixed-h removal of geometric smoothing, pulse-evolved length
domains, autonomous gravitational dynamics and physical 3+1 recovery remain
separate obligations. No speed, memory, quantum or hardware advantage is
claimed by this refinement.
