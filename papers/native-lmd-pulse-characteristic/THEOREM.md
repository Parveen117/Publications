# Uniform LMD pulses: an exact characteristic law and a bounded solver test

Author: Monty. 30 October 2026.
Frozen research parent: `70aef2535a1c124a410e4a4bf96780236ca25c23`.

This packet executes the existing NRQC calculator on the new common-parent
pulse-moment contract, then advances the analytic comparison from a variance
to the complete one-dimensional conjugate-label characteristic function.
The finite computation and continuum derivation are separate evidence.
There is no claim of a continuum pulse simulation or solved physical gravity.

## 1. Source, state and detector

Reuse [LMD](../native-logarithmic-material-density/THEOREM.md) and the
[common noisy parent](../native-lmd-common-parent/THEOREM.md), without changing
their covariance or geometry prescription. Work in the positive sector,
in the vacuum psi=1, with 0<h<1/2. All variables here are dimensionless.
Set

    J = M_h(S1),             D = D(1) = i h partial_q - i q/2,
    q = s_0 ~ N(0,h),       U_t = exp(-i t J/h).

The other Fourier modes are independent of q. Constant Cameron--Martin
translation scales the total source mass exactly:

    J(q+a,eta) = exp(a) J(q,eta).

D is the total affine conjugate, not the material-reference generator R[1].
The measured v_0 marginal is the spectral law of D in U_t 1 convolved with
an independent N(0,h/4) detector noise. This is the existing parent marginal;
it does not assert that all quantum observables have simultaneous sharp values.

Write chi_h(a)=E_mu exp(i a J). Its argument is real, so it exists without
an exponential moment of J. E J=1 and E J^2=I(h), where

    I(g) = integral_(S1 x S1) exp(g K(z-y)) dz dy
         = exp(g) Gamma(1-g) / Gamma(1-g/2)^2,      0<g<1.       (1)

To obtain (1), use K(r)=1-log(2|sin(pi r)|), substitute theta=pi r,
evaluate the sine power integral with Euler's beta integral, and apply the
gamma duplication formula. These are standard identities, credited in
SOURCE_MAP.md; the gamma ratio is not a new special-function result.

## 2. Exact characteristic function of the actual pulse

For every real t,k, the actual measured characteristic function is

    Phi_t(k) = E_t exp(i k v_0)
             = exp(-h k^2/4) chi_h(a_t(k)),
    a_t(k)   = (2t/h) sinh(h k/2).                             (2)

Proof. Multiply a q-wavefunction by the square root of its N(0,h) density.
This unitary map takes D to i h partial_q on Lebesgue L2, hence exp(i k D)
translates the flattened wavefunction by -h k. The vacuum matrix element
after U_t is therefore the overlap of the two Gaussian amplitudes, with
phase exp(i t[J(q,eta)-J(q-hk,eta)]/h). Center at r=q-hk/2. The product of
Gaussian square roots becomes the original N(0,h) density in r times
exp(-h k^2/8), while the phase difference becomes
2 sinh(hk/2) J(r,eta). Integrate over r and eta. The independent marginal
detector contributes the second factor exp(-h k^2/8). This proves (2).

All translations and matrix elements in this proof are bounded unitaries;
the integrand's phase has modulus one. No D^n J^m product-domain assumption,
Taylor expansion in an unbounded J, or samplewise classical label update
is needed. The cutoff identity below gives the same result by an L2 limit.

Equation (2) determines this one marginal, not the full joint (u,v) outcome
law or the nonlinear coframe law. Although chi_h is the characteristic
function of a positive random variable, k -> a_t(k) is nonlinear. Equation
(2) must not be read as the characteristic function of independent Gaussian
noise plus tJ; that would use the argument tk.

The known moments follow by two derivatives at k=0:

    E_t v_0 = t,
    Var_t(v_0) = h/2 + t^2 [I(h)-1].                           (3)

The L2 mass moment suffices for these derivatives. Their interpretation
as spectral moments also follows directly from the frozen Stone-domain
calculation D U_t 1 = U_t(tJ-iq/2).

## 3. A bounded discriminator, already at linear response

The tempting vacuum-label rule v_0 -> v_0+t Mhat(u,S1) has characteristic
function

    Phi_naive,t(k) = exp(-h k^2/4) chi_(2h)(tk),               (4)

because the noisy Wick measure has parameter 2h and is independent of v_0
in the vacuum parent law. It gives the correct mean and the wrong variance,
h/2+t^2[I(2h)-1], as proved in the predecessor.

There is a more sensitive test that uses only a bounded reading. Since
E J = E Mhat(S1)=1, dominated difference quotients in t give

    d/dt E_t sin(k v_0)|_(t=0)
       = exp(-h k^2/4) (2/h) sinh(h k/2),
    d/dt E_naive,t sin(k v_0)|_(t=0)
       = exp(-h k^2/4) k.                                    (5)

The slopes differ for every nonzero k. For k>0 their difference is positive;
it is of order h^2 k^3 near k=0. Thus even replacing Mhat in the naive rule
by a mass variable with the correct source distribution would not fix the
full law. The nonlinear characteristic argument is an additional quantum
translation effect, distinct from the doubled density covariance.

For fixed h,k!=0 the nonzero derivative difference implies different sine
expectations for all sufficiently small nonzero t. This is a local existence
statement; no finite-t threshold is inferred from the derivative alone.
The observed slope can be estimated by repeated finite-difference experiments,
but no experimental or shot-based data are claimed here.

The native calculation encloses the positive slope difference at k=2 for
h=1/8,1/4,3/8 using positive sinh/exponential series with remainder bounds.
Those values are evaluations of (5), not NRQC evaluations of chi_h.

## 4. A controlled cutoff route for an actual marginal solver

Let J_N be LMD's total Fourier-cutoff Wick mass, including the zero mode,
with J_0=exp(q-h/2). The positive martingale property and L2 convergence give

    delta_N^2 = ||J-J_N||_2^2 = I(h)-I_N(h),
    I_N(h) = exp(h) [z^0] exp(sum_(n=1)^N h(z^n+z^-n)/(2n)).   (6)

Here [z^0] denotes the constant Laurent coefficient, equivalently the circle
integral. Although K_N is not pointwise monotone, I_N is monotone as the
second moment of an L2 martingale. Do not replace this argument with a
pointwise comparison of truncated logarithmic kernels.

The cutoff pulse has the exact marginal characteristic function

    Phi_(t,N)(k) = exp(-h k^2/4) E exp(i a_t(k) J_N).           (7)

By |exp(ix)-exp(iy)|<=min(2,|x-y|),

    |Phi_t(k)-Phi_(t,N)(k)|
       <= exp(-h k^2/4) min(2, |a_t(k)| delta_N).              (8)

This is uniform for t,k in compact sets and also gives a frequency-weighted
bound on the whole real k-axis. It turns a certified mass-cutoff error into
an error budget for bounded sine and cosine measurements. A future solver
must additionally control its quadrature/state truncation, gate approximation
and calculator error; (8) supplies only the continuum-to-mode-cutoff part.

The native code encloses I_N for N=0,1,2,4,8,16. Its 48-term exponential
Laurent expansion has nonnegative coefficients; each rational convolution
is rounded outwards. The unretained constant coefficient is bounded by the
scalar exponential tail at x=h sum_(n<=N) 1/n. The N=1 calculation is also
compared with the independent positive Bessel series

    I_1(h)=exp(h) sum_(m>=0) (h^2/4)^m/(m!)^2.

These executed source calculations show strict second-moment progress at
the declared cutoffs. They do not evaluate E exp(i a J_N), and do not certify
that N=16 meets an arbitrary requested continuum tolerance.

Equation (2) also yields |Phi_t(k)|<=exp(-h k^2/4), uniformly in t. Fourier
inversion gives a smooth marginal density p_t with

    ||p_t||_infinity <= 1/sqrt(pi h).

Every derivative exists since |k|^m exp(-h k^2/4) is integrable. The cutoff
densities converge uniformly, with the explicit bound

    ||p_t-p_(t,N)||_infinity <= C(h,T) delta_N,  |t|<=T,
    C(h,T) = (1/(2pi)) integral_R exp(-h k^2/4)
                                  (2T/h)|sinh(hk/2)| dk < infinity.   (9)

This regularity concerns a measured conjugate marginal. It gives no
exponential moment of the geometric length and does not repair the previous
ungated infinite-mean-length result or prove pulse-invariance of gated domains.

## 5. What the fresh calculator run does

In the vacuum, the source-centered covariance of D and J is

    Q = [[h/4, i h/2],[-i h/2, I(h)-1]].                     (10)

For the sign, D1=-iq/2 and Gaussian integration by parts gives E(qJ)=h;
therefore <D1,(J-1)1>=+ih/2. Also I(h)>exp(h)>1+h, so Q is strictly positive.

The existing NRQC signed-Pauli interface is used exactly as documented.
Set A=h/4, B=I(h)-1, c=h/2 and choose 28-bit dyadic approximations

    alpha ~= sqrt(A), beta ~= c/alpha,
    gamma ~= sqrt(B-beta^2).

On three levels embedded in two qubits, use finite controls

    D_f = alpha (|1><0|+|0><1|),
    J_f = I + i beta |1><0| - i beta |0><1|
              + gamma (|2><0|+|0><2|).                       (11)

The vacuum means are 0,1 and the Gram entries are alpha^2,
beta^2+gamma^2,+i alpha beta. The coefficients come from native rational
intervals for the source; NRQC does not derive those continuum coefficients.
The error of each supplied Gram entry is certified at most 10^-6 for the
three h values. J_f is a moment-matching control, not a positivity-preserving
representation of the unbounded positive source J.

NRQC computes means and covariance G for six Pauli words

    X0, Z1 X0, Y0, Z1 Y0, X1, X1 Z0.

Native JS then forms W G W^dagger using the explicit Pauli decomposition
of (11). It applies the source Heisenberg moment shear D -> D+tJ for
t=-1,-1/2,0,1/2,1, and adds the h/4 detector variance. Independent explicit
4x4 matrices check both the means and centered covariances at every point.

Eight fresh NRQC executions use observer and clifford_t engines in four
states: vacuum, Bell, controlled-S on |++>, and controlled-S-dagger on |++>.
The controlled phase circuits execute non-Clifford T gates. Their explicit
independent state vectors are (1,1,1,+i)/2 and (1,1,1,-i)/2. Only the vacuum
is a source-moment encoding; the other states test the finite adapter and
its complex phase conventions. These phase circuits are not LMD pulse gates.

The finite control pair cannot satisfy [D_f,J_f]=ih J_f: a finite commutator
has trace zero, whereas ih J_f has trace 4ih. This is an executed negative
control and an elementary analytic obstruction. Consequently the supplied
Heisenberg shear is not being implemented by exp(-it J_f/h) on that finite
carrier. It is essential not to report these covariance checks as quantum
time evolution, a full affine-algebra representation, or a solution for
the common-parent outcome law.

## 6. Closure and next implementation gate

The executed calculation checks the stated finite contract, verifies the
variance separation against a native exact baseline, and supplies certified
source cutoff errors. The new written result (2) gives an exact all-t,k
marginal law, the bounded response witness (5), and a controlled cutoff
route (8)-(9). These are conditional results within the frozen LMD source
and chosen common measurement.

The next computational task is E exp(i a J_N) with certified approximation
and an honest circuit/representation map, followed by nonuniform and joint
readouts. Detector/calibration covariance under general material flows,
fixed-h removal of spatial smoothing, pulse-evolved length domains, local
apparatus, autonomous gravitational dynamics and physical 3+1 recovery remain
open. No speed, memory or quantum advantage is inferred from these tiny tests.
