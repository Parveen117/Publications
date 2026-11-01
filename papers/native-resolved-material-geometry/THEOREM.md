# RMG-1 — Resolved material geometry and its continuum obstruction

16 October 2026. Frozen parent: RGR-1 at
d9c70208520f53f31d93bd94c3190f48b34acb4c. This packet acts on the
selected positive physical space already constructed by RGR. It neither
classifies further gauge channels nor changes the source, gauge measure,
MPQ covariance, sign sectors or frozen pointer coupling.

The additional contract is a finite-resolution coherent-state joint
measurement of MPQ's adapted physical coordinates. Its nonlinear output
is the original reduced material geometry. This is a POVM prescription,
not an assertion of mutually sharp canonical field operators. We prove
its normalization, useful operator domains, information content and a
controlled classical limit. We also prove that its raw pointwise geometry
has no continuum limit at fixed nonzero hbar, even in the vacuum. Thus
the prescription does not close the continuum composite-field problem.

## 1. Source geometry and the declared resolution

MRF-T3 gives the physical profiles r,m and the representative
Z=x, phi=tau, q=r, p=m, k=0, rho=m'-m r'. LNC's q=log e therefore
identifies exp(r) dz as the spatial coframe on this material slice.
The slice metric is exp(2r) dz^2 and its circumference is integral exp(r).
These statements use the inherited independent-coframe model; they do
not impose a new global circle-cut winding condition. A spatial profile
alone does not fix spacetime curvature, lapse or normal evolution.

Keep the exact MPQ inverse chart, with epsilon=+1 or -1:

    a(theta)=1+(1/2)cos theta, Y(theta)=integral_0^theta dt/a(t),
    j=epsilon exp(s), y=d/j, theta=Y^{-1}(y),
    m=j/a(theta), P=integral m, r=theta+log|m/P|.           (1)

In particular P is computed from the reconstructed profile; it is not
held fixed during inversion. theta and y are real lifts. Let
omega=sqrt(3)/2. The elementary tangent-half-angle integral gives
Y(theta+2pi)=Y(theta)+2pi/omega. Consequently

    Y^{-1}(y)=omega y+b(y),                                (2)

where b is smooth, bounded and periodic with period 2pi/omega.
The map Y^{-1} is odd. Also 1/2<=a<=3/2 everywhere.

Use MPQ's real Fourier modes through n=N, dimension k_N=2N+1,
and diagonal covariance C_N with lambda_0=1 and
lambda_n=(1+n^2)^(-4), twice for n>=1. At physical hbar=h>0,
H_epsilon=L^2(mu_h), mu_h=N(0,h C); h=1 is the existing endpoint.
Only the MPQ semiclassical family reinstates h. RGR's residual gauge
covariance is not rescaled here. Product factorization gives
H_epsilon=H_N tensor H_>N with H_N=L^2(N(0,h C_N)).

Measurement labels u,v are arbitrary real k_N vectors. Their Fourier
profiles are u_N(z),v_N(z). Substituting s=u_N,d=v_N in (1) defines
smooth outputs r_N,m_N,e_N=exp(r_N), and L_N=integral e_N.
They are generally not band limited. Each output has definite sign m_N
and finite positive L_N. No cutoff of the old constraints is performed.
The reconstructed reference momentum also has the exact simplified form

    rho_N=m_N'-m_N r_N'=-j_N y_N'=-v_N'+v_N u_N'.         (1a)

Indeed r'=theta'+m'/m, theta'=a y', and d=j y with j'=j s'.
This retains the reference constraint on each outcome; it does not
define a sharp continuum rho operator.

## 2. A normalized measurement on the existing physical Hilbert space

For s in the finite coordinate space define normalized vectors

    chi_(u,v)(s)=exp[u.C_N^{-1}s/(2h)-u.C_N^{-1}u/(4h)]
                 exp[-i v.s/h],
    dgamma_N=du dv/(2pi h)^{k_N}.                          (3)

**RMG-T1 (physical geometry POVM).** In the weak operator sense,

    integral |chi_(u,v)><chi_(u,v)| dgamma_N=I_HN.         (4)

Tensor with I_H>N and include epsilon as a recorded central label.
Push this POVM forward by (1). This defines a countably additive,
normalized positive operator-valued measure E_N on the smooth reduced
profile outcomes (r_N,m_N). For every normal physical state rho,
Tr[rho E_N(B)] is a probability measure. It applies to general continuum
states; the state need not lie in MPQ's finite-mode image.

Proof. Multiplication by the square root of the Gaussian density maps
H_N unitarily to L^2(ds) and maps (3) to

    (det(2pi h C_N))^{-1/4}
    exp[-(s-u).C_N^{-1}(s-u)/(4h)] exp[-i v.s/h].          (5)

Each vector has norm one. Integrating v in its rank-one kernel gives
delta(s-t); the remaining u Gaussian integral is one. This proves (4)
first on Schwartz functions and then by positivity and completion.
The smooth inverse (1) is a continuous map from finite-dimensional
labels into smooth profiles, so its pushforward is a Borel POVM.

More explicitly the analysis map

    (V_N psi)(u,v)=(<chi_(u,v)| tensor I)psi               (6)

is an isometry into L^2(dgamma_N;H_>N). Then
E_N(B)=V_N^* 1_{(r_N,m_N) in B} V_N. This supplies a concrete dilation,
not an inference from kinematic positivity. The inner product is the
already derived RGR/MPQ physical norm. The bounded effects act on its
entire domain and can be lifted on RGR's gauge-test tensor physical
space as I tensor E_N(B); the quotient map intertwines these effects.
No extra classification of invariant distributions is needed.

An explicit instrument exists: with respect to dgamma_N use Kraus
density K_(u,v)=|chi_(u,v)><chi_(u,v)| tensor I. Since K^*K=K,
(4) proves normalization. This is a declared finite-mode measurement
instrument, not a derivation of a new source apparatus. For a fixed
clock slice or frozen pulse, conjugating its effects/instrument by the
existing MPQ unitary gives the corresponding relational measurement.

## 3. Nonlinear observables, domains and retained phase information

For any bounded real Borel outcome function F(r,m), the operator
A_N(F)=V_N^* F(r_N,m_N)V_N is bounded self-adjoint on all H_red, with
norm at most sup|F|. This includes bounded readings of coframe,
circumference, profile values and material-window integrals.

**RMG-T2 (unbounded moment domains and positive momentum readings).**
For a real output f=f(u,v), define the second-moment domain and operator

    D_2(f)={psi: integral |f|^2 ||V_N psi||^2 dgamma_N<infinity},
    T_f psi=V_N^* [f V_N psi], psi in D_2(f).             (7)

This operator is symmetric on that domain. For f=r_N(z), m_N(z),
or their smooth material smearings, D_2(f) is dense. All finite powers
have dense moment domains. No self-adjointness of T_r or field CCR is
asserted. For f=|m_N(z)|, or f=integral w|m_N| with w>=0 smooth,

    q_f(psi)=integral f ||V_N psi||^2 dgamma_N,
    D(q_f)={psi:q_f(psi)<infinity}                        (8)

is a densely defined closed positive form. Its associated positive
self-adjoint operator M_f is uniquely specified by this form. Explicitly
its domain consists of psi in D(q_f) for which there is an eta in H
such that q_f(phi,psi)=<phi,eta> for every phi in D(q_f), and M_f psi=eta.
Multiplication by epsilon gives the signed momentum mean observable.
The spectral measurement of M_f need not equal the original POVM
marginal; agreement of their first moments is not agreement of all
outcome probabilities.

Proof. Put U=||u_N||_infinity. The bounds on a give

    (2/3)exp(u_N)<=|m_N|<=2exp(u_N),
    |log|m_N/P_N||<=2U+log 3,
    |r_N(z)|<=omega |v_N(z)|exp(U)+||b||_infinity+2U+log 3. (9)

At fixed N these bounds are polynomial in v times an exponential in
the Euclidean norm of u. The coherent transforms of finite Fourier
Hermite polynomials are polynomials in u,v times the vacuum transform.
Their squared moduli have Gaussian tails in both labels. Thus all finite
powers in (9) are integrable on their tensor products with arbitrary
tail vectors. That subspace is dense. Smearing preserves these estimates.
Symmetry in (7) follows from real multiplication in the dilation space.
For (8), M_sqrt(f)V_N is closed: convergence of psi and of its image
reduces the assertion to the closed multiplication operator on the
dilation space. Its norm square is therefore a closed positive form.
The closed-form representation theorem supplies exactly the stated
self-adjoint domain, without an essential-self-adjointness guess.

The joint measurement is not a joint projection-valued measure. For
example the u marginal is the sharp S_N measurement convolved with
N(0,h C_N); its second moment is S_j^2+h lambda_j, rather than S_j^2.
Likewise the v marginal adds variance h/(4lambda_j) to D_j's spectral
variance. The extra variances have product h^2/4. In particular

    A_N(epsilon exp[u_N(z)])
       =exp[h K_N(z,z)/2] epsilon exp[s_N(z)],             (10)

where K_N=sum lambda_j e_j e_j. This finite-resolution reading is not
the unchanged sharp J probe. The old G_w pointer is not relabelled.

**RMG-T3 (information beyond the frozen probe).** Already at N=0,
the coherent states with the same s center and opposite nonzero d
centers have identical sharp J(f) distributions but different geometry
POVM distributions. In fact r_0=Y^{-1}(epsilon v exp(-u)), so
the probability of r_0>0 is the probability of epsilon v>0. The v laws
are Gaussians of opposite means and identical nonzero variance. Their
probabilities differ strictly. Thus this measurement distinguishes
physical conjugate information that the frozen probe alone misses.
This is not a claim that one simultaneous outcome determines the state,
nor that the original pointer apparatus implements this added instrument.

## 4. Exact noise law and a fixed-hbar continuum obstruction

In a smooth MPQ coherent state centered at s_*,d_*, the label law is
exactly a product of independent Gaussian pairs:

    u=P_N s_*+xi,   Cov(xi)=2h C_N,
    v=P_N d_*+eta,  Cov(eta)=(h/2) C_N^{-1},
    xi independent of eta.                              (11)

Indeed the squared coherent overlap is
exp[-(u-u_*).C_N^{-1}(u-u_*)/(4h)
    -(v-v_*).C_N(v-v_*)/h]. With dgamma_N it integrates to one.
This law follows from (3), including both state and measurement noise.

**RMG-T4 (no raw continuum point geometry for this prescription).**
Fix h>0 and the normalized physical vacuum psi=1. At each material
point z, for every finite R,

    Pr(|r_N(z)|<=R) -> 0 as N -> infinity.                (12)

Hence these resolved geometry POVMs do not converge, on this state, to
a finite real point reading r(z), or to a probability law on continuous
r profiles with the corresponding evaluation marginals. This is a
failure of the specified inverse-chart measurement, not a theorem that
no other renormalized or distributional quantum geometry can exist.

Proof. The independent-mode laws (11) permit a common Gaussian coupling.
The u_N converge almost surely in H^3, since

    E||u_N||_H3^2=2h[1+2 sum_(n=1)^N (1+n^2)^(-1)]

is uniformly bounded and the Gaussian series converges in H^3. Thus
sup_N ||u_N||_infinity is finite almost surely. At every z, the identity
(e_n^c(z))^2+(e_n^s(z))^2=2 gives the exact variance

    sigma_N^2=Var(v_N(z))
       =(h/2)[1+2 sum_(n=1)^N (1+n^2)^4]
       ~ (h/9) N^9.                                    (13)

It is independent of u_N. On the event ||u_N||_infinity<=M, (2),(9)
imply that |r_N(z)|<=R requires
|v_N(z)|<=exp(M)(R+||b||_infinity+2M+log 3)/omega.
The centered Gaussian probability of that fixed interval is at most
its length divided by sqrt(2pi) sigma_N, hence tends to zero. The
complement probability can be made uniformly small by increasing M.
This proves (12), not merely a divergence of second moments. A heavy
tail with vanishing weight would not have sufficed for this argument.

There is a second, distinct warning at fixed resolution. In the vacuum,
even N=0 has

    E[L_0]=infinity, although 0<L_0<infinity per outcome.  (14)

Here m_0=P_0 and L_0=exp[Y^{-1}(epsilon v exp(-u))]. By (2), this
is at least exp(-||b||_infinity) exp[omega epsilon v exp(-u)].
Integrating v first using (11) gives
exp[(omega^2 h/4) exp(-2u)]. Its integral against N(0,2h) diverges,
since its positive double-exponential tail dominates the Gaussian
quadratic exponent as u tends to minus infinity. This proves (14)
by Tonelli. In particular the vacuum is not in the first-length-moment
form domain. We claim neither density of that domain nor an unbounded
length mean operator. Bounded detector readings such as arctan L_N and
length-bin effects remain well-defined on every physical state.

## 5. Controlled recovery of smooth classical geometry

**RMG-T5 (simultaneous classical/resolution limit).** Fix any smooth
regular MRF classical profile r_*,m_*. Use the exact MPQ chart to obtain
smooth s_*,d_*. In its h-coherent states, choose N(h)->infinity such that

    h N(h)^15 -> 0 as h -> 0.                            (15)

Then the geometry measurement outcomes converge in probability to
(r_*,m_*) in H^3 and hence in C^2. The coframe e_N, slice metric
exp(2r_N) and circumference L_N converge in probability to their
classical values. Every bounded continuous test of these outcomes has
the corresponding classical limiting expectation. Unbounded moment
convergence is not claimed; (14) explicitly prevents it for length in
the zero-center case. Condition (15) is sufficient, not asserted optimal.

Proof. The two exact noise budgets in H^3 are

    E||xi||_H3^2=2h[1+2 sum_(n=1)^N (1+n^2)^(-1)]=O(h),
    E||eta||_H3^2=(h/2)[1+2 sum_(n=1)^N (1+n^2)^7]
                  =O(h N^15).                           (16)

Markov's inequality and smooth Fourier convergence imply (u_N,v_N)
converge in probability in H^3 to (s_*,d_*). In one dimension H^3 is
an algebra embedded in C^2. On a bounded neighborhood of a fixed smooth
pair, exponentiation, division by epsilon exp(u), composition with
Y^{-1}, and division by a(theta) are continuous H^3 maps. The integral
P stays nonzero and log|m/P| is continuous there. Thus (1) is continuous
in H^3 at that pair. The continuous mapping theorem proves the profile
claim and then the coframe, metric and circumference claims. Bounded
continuous tests follow from convergence in probability to a constant.
At fixed N the same proof gives the ordinary h->0 classical limit of
the truncated adapted profiles. The h=constant, N->infinity order of
limits is instead excluded by T4.

## 6. What this closes and what it leaves to the two research tracks

Constructed on the existing selected positive physical space: normalized
finite-resolution nonlinear geometry probabilities, a concrete instrument,
bounded geometry observables, dense moment domains for r and m, positive
self-adjoint m mean observables, a phase-sensitive geometry reading, and
a controlled C^2 classical recovery with a quantified resolution budget.
The fixed-hbar point limit and the naive mean circumference fail under
this particular added measurement contract. These are analytic theorems;
finite exact witnesses below do not certify infinite-dimensional analysis.

The geometry track should next address a justified coarse/renormalized
composite prescription with controlled fixed-hbar behavior and a stated
length domain, or prove a broader obstruction with explicit assumptions.
The present results do not produce sharp continuum r,m fields, a
spacetime curvature operator or an experimentally calibrated apparatus.
Do not identify one-dimensional spatial coframe data with spacetime
curvature. Full physical 3+1 quantum gravity remains open.

The parallel external report at 60fd188 on
claude/residual-gauge-invariant-dual concerns gauge distributions. Its
next independent task is actual residual-Diff completeness/uniqueness
and transfer of its Hida result beyond the abelianized action. None of
that classification is duplicated or assumed here, and that branch is
not merged by this packet. RGR's selected existence is the dependency.

Coherent-state/anti-Wick quantization is established mathematics, credited
in SOURCE_MAP.md. The source-specific results are the exact nonlinear
profile reconstruction, domain bounds, physical phase reading, N^9
non-tightness proof, length-tail counterexample and h N^15 recovery
budget in the frozen covariance. Predecessor proofs remain unchanged.
