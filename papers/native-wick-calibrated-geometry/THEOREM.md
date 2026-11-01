# WCG-1 — Wick-calibrated joint readings give resolved material geometry

16 October 2026. Frozen parent: `82b362b94f460cc4679bdc51ebc9cee5af9a4beb`.

## Source contract and result

PCR recovered the phase from original current means but left a nonlinear
quantum material reading to construct. Work in the same selected LMD
sector H_epsilon=L2(mu_h), with lambda_0=1 and lambda_n=1/(2n) for the
paired real Fourier modes n>=1. Keep its positive Wick measure M, original
S(f), D(k), corrected local currents and classical material inverse.
Here epsilon is the fixed orientation sign, and the circle has length one.

For this measurement construction assume 0<h<1/2. The system itself
retains LMD's wider 0<h<1 range. The stricter bound is the elementary L2
range for a detector density with doubled Gaussian covariance; it is
neither a measured Planck constant nor a necessity theorem outside this
particular L2 argument.

Adapt the already constructed CMG/MDD coherent joint detector to this
LMD covariance. Its auxiliary Gaussian probe and complementary readout
remain declared resources. This does not identify the LMD representation
with the old MPQ one, change the source engine, or derive a local physical
detector from the old density-only interaction.

The new density output is a measured Wick measure N_u, rather than the
exponential of a smoothed configuration label. Its first-moment form is
exactly the existing M(w). The original current D(k) is the other first
moment, with its full moment domain. A fixed positive spatial response
scale ell gives a normalized nonlinear material-profile POVM, bounded
metric observables and a thresholded positive self-adjoint length mean.
PCR's phase-pulsed coherent states recover the smooth classical geometry
when h ell^(-8)->0. At fixed h, removing ell still fails for raw point r.

## WCG-T1. Joint detector with exact Wick-density and current calibration

Let nu_h be the countable product Gaussian probability law with independent
real coordinates

    u_j ~ N(0,2h lambda_j), v_j ~ N(0,h/(2lambda_j)).      (1)

For normalized finite Hermite multi-indices H_alpha in H_epsilon, put

    zeta_j=u_j/(2 sqrt(h lambda_j))+i sqrt(lambda_j/h)v_j,
    W H_alpha=product_j zeta_j^alpha_j/sqrt(alpha_j!).     (2)

Then W extends to an isometry, and E(B)=W*1_B W is a normalized countably
additive POVM. Every normal-state output law is absolutely continuous
with respect to nu_h. Typical infinite sequences do not have to label
individual normalizable coherent vectors.

Put c_N=1+sum_(n=1)^N 1/n and u_N=P_Nu. The positive measures

    N_(u,N)(dz)=exp[u_N(z)-h c_N] dz                     (3)

converge weakly almost surely to a finite nonzero measure N_u with full
support and no atoms. For bounded w they converge in L2(nu_h). For
nonnegative bounded w and every psi in H_epsilon,

    integral N_u(w)|W psi|² dnu_h
       =integral M_s(w)|psi(s)|² dmu_h(s),               (4)

as an extended nonnegative identity. Thus the positive first-moment
form has exactly domain Dom M(w)^(1/2) and associated self-adjoint
operator M(w). Multiplication by epsilon gives the signed density.
Equation (4) does not assert equality of sharp and noisy spectral laws
or identify the second-moment domain of the POVM with Dom M(w).

For smooth real k, v(k)=sum k_j v_j exists as a distributional pairing.
Its marginal is the spectral distribution of D(k) convolved with a
centered Gaussian of variance sigma_k²=(h/4)||k||_C². In particular

    Dom(M_(v(k)) W)=Dom D(k),
    W* M_(v(k)) W psi=D(k)psi,
    ||v(k)W psi||²=||D(k)psi||²+sigma_k²||psi||².          (5)

The last identity is extended-valued on arbitrary vectors and finite
exactly on Dom D(k). Density and current are jointly unsharp outputs;
their calibration does not give a joint sharp system PVM.

**Proof.** The circular Gaussian identity
E[zeta^a conjugate(zeta)^b]=delta_ab a! makes (2) an isometry on the
dense finite Hermite span. This is CMG-T1's covariance-independent
construction with the explicitly changed lambda_j. Compression of
multiplication projections proves the POVM statements. A normal state
sum p_a |psi_a><psi_a| has density sum p_a |W psi_a|² with respect to nu_h.

MDD's balanced mixing proof applies to two identical copies of mu_h:
rotating each system/probe Gaussian pair preserves their product law.
The countable rotations are strongly continuous by their action on
dense cylinder vectors. Complementary output quadratures pull back to
S_j+A_j and D_j-B_j, where A,B belong to the independent probe. Their
strong commutation follows from unitary conjugation of the separate
output PVMs. Preparing the probe vacuum and compressing reproduces (2),
first on finite marginals and then on all cylinder-generated events.
Thus this is an actual ideal dilation with the already declared extra
probe resource, not an inference from a formal vanishing commutator.

Under nu_h, u has covariance 2h C. Apply LMD-T3 at parameter 2h<1 to
obtain (3) and its L2 convergence. The full zero-mode subtraction is
included. The configuration marginal of the joint detector is the
kernel u=s+a, where a is an independent sample with law mu_h. At each
cutoff, Gaussian integration gives

    E_a exp[s_N+a_N-h c_N]=exp[s_N-h c_N/2].             (6)

In L2(mu_h tensor mu_h), conditional-expectation contraction passes
(6) to the limit:

    E_a[N_(s+a)(w)|s]=M_s(w).                            (7)

The noisy-configuration kernel follows from the finite coherent
marginals, and its extension is fixed by cylinder sets. Tonelli applied
to |psi(s)|² and (7) proves (4), even if its value is infinite. The
multiplication positive form then specifies exactly its operator/domain.

The pulled-back current is D(k)-B(k). Its two terms strongly commute
on different factors, and the probe vacuum B(k) law is the stated
Gaussian. Its characteristic function proves the convolution claim.
Tonelli applied to the squared output proves (5)'s exact domain. On
that domain, integrating the zero-mean probe term gives the compression
D(k); polarization gives the corresponding operator identity.

For example, the vacuum has

    E N_u(w)=integral w,
    E[N_u(w)N_u(v)]=integral integral w(z)v(y)e^{2hK(z-y)}, (8)

whereas the sharp system M second moment has e^{hK}. The detector
fluctuations are different despite the exact first moment. Bare exp(u),
or using only the system half-variance counterterm in (3), would fail
the calibration. No arbitrary high positive moments of N_u are asserted.

## WCG-T2. Nonlinear material observables at fixed spatial resolution

Let Lambda e_n=n² e_n and use the positive circle heat response

    H_ell=exp(-ell² Lambda), t_n=exp(-ell²n²), ell>0.

This is a declared material resolution, not physical time evolution.
Its convolution kernel p_ell is smooth, strictly positive, and integrates
to one. Form the actual outcome profiles

    q_ell(z)=integral p_ell(z-y)N_u(dy),
    j_ell=epsilon q_ell, d_ell=H_ell v,
    y_ell=d_ell/j_ell, vartheta_ell=Y^(-1)(y_ell),
    a(vartheta)=1+cos(vartheta)/2, Y(vartheta)=integral_0^vartheta dt/a(t),
    m_ell=j_ell/a(vartheta_ell), P_ell=integral m_ell,
    r_ell=vartheta_ell+log|m_ell/P_ell|,
    g_ell(z)=exp(2r_ell(z)), L_ell=integral exp(r_ell).     (9)

At every fixed ell these are smooth almost surely in every normal state;
q_ell is strictly positive, m_ell has sign epsilon, and P_ell never
vanishes. Pushforward of E by (9) defines a normalized geometry POVM.
Every bounded real Borel function F of its outputs defines a bounded
self-adjoint operator

    A_ell(F)=W* F(r_ell,m_ell) W, ||A_ell(F)||<=||F||inf. (10)

An explicit material-metric reading on all H_epsilon is

    G_ell(z)=W* [g_ell(z)/(1+g_ell(z))] W.               (11)

It is a positive contraction, with strictly positive expectation of
both G_ell(z) and I-G_ell(z) in every nonzero vector. The bounded
monotone response in (11) retains the outcome metric value; its operator
expectation is not the metric of a reconstructed mean state.

The outcome satisfies the original inverse/source identity

    rho_ell=m_ell'-m_ell r_ell'
           =-j_ell y_ell'
           =-d_ell'+d_ell (log q_ell)'.                 (12)

This is an identity of the reconstructed smooth outcomes, not a sharp
composite quantum constraint or a claim about spacetime curvature.

Positive finite-resolution approximants are obtained by

    q_(ell,N)=p_ell*N_(u,N), d_(ell,N)=H_ell P_Nv.        (13)

They and their inverse profiles converge almost surely in every C^k,
with the C^infinity x C^infinity topology on the pair (r,m).
Consequently their bounded continuous output operators converge strongly
to (10), and all normal-state output laws converge weakly. For finite
detector partitions with nu_h-null inverse-image boundaries, the bin
effects converge strongly and their square-root instruments converge
in trace norm on each trace-class input.

Here the chosen instrument is I_b(rho)=sqrt(E(B_b)) rho sqrt(E(B_b)).
It need not coincide with MDD's apparatus-induced instrument; the same
effects do not uniquely determine the postmeasurement state.

**Proof.** Gaussian Fourier moments give u in H^(-eta) and v in
H^(-1-eta) almost surely for every eta>0. Heat smoothing sends the latter
to C^infinity. The finite positive measure N_u gives
q_ell>=min(p_ell)N_u(S1)>0 and a smooth convolution. The same events
hold for all normal states by T1's absolute continuity; this transfers
almost-sure regularity, not vacuum moment estimates for arbitrary states.
The inverse chart has smooth compositions and denominators bounded away
from zero on each realized compact circle. It therefore gives (9).
Positivity, bounded domains and (11) follow by compression.

Differentiate r=vartheta+log|m/P| spatially, retaining spatially constant
P. Since vartheta'=a y' and m a=j, the first equality in (12) is
-m vartheta'=-j y'. Using d=j y proves the last one.

Weak convergence of positive measures N_(u,N), including convergence
and boundedness of their total masses, gives uniform convergence of
their integrals against the compact family of translates of each
derivative of p_ell. Hence q_(ell,N)->q_ell in every C^k. Heat-smoothed
Fourier v converges likewise. Continuity of the inverse at strictly
positive q proves the profile statement. In particular (13) smooths
the positive measure itself; Fourier truncation of q_ell could lose
positivity and is not the prescription. Bounded dominated convergence
on the dilation space proves strong convergence. Continuous functional
calculus on [0,1] gives convergence of square-root effects; finite-rank
approximation gives the trace-norm instrument conclusion.

## WCG-T3. Accepted circumference has a positive self-adjoint mean

Fix delta>0, accept on A_(ell,delta)={min q_ell>=delta}, and record a
separate rejection symbol on its complement. No probability is discarded
or silently renormalized. Define f_L=1_A L_ell, assigning the number zero
to the rejection flag only for bookkeeping. The positive form

    Q_L(psi)=integral f_L |W psi|² dnu_h,
    Dom Q_L={psi:Q_L(psi)<infinity}                       (14)

is closed and densely defined. Its associated positive self-adjoint
operator mathsf L_(ell,delta) has the exact domain

    {psi in Dom Q_L: exists eta in H such that
       Q_L(phi,psi)=<phi,eta> for every phi in Dom Q_L},
    mathsf L_(ell,delta)psi=eta.                         (15)

Every finite Fourier Hermite vector, including the vacuum, has all
finite accepted-length POVM moments and belongs to (15). On the
sufficient second-moment domain it acts as W*(f_L W psi). Its sharp
spectral law need not be the original length POVM law; higher POVM
moments are not asserted to equal powers of this compressed operator.

The rejection-inclusive bounded readings have the same strong cutoff
limit as T2. Accepted-length cutoff moments converge on the finite
Hermite span. This is not an assertion of strong-resolvent convergence
for all unbounded cutoff length operators.

**Proof.** RMG's inverse obeys
Y^(-1)(y)=omega y+b(y), omega=sqrt(3)/2, B=||b||inf<infinity.
The definite sign and global normalization in (9) give

    L_ell=integral exp(vartheta_ell) |m_ell|/|P_ell|,
    integral |m_ell|/|P_ell|=1,
    1_A L_ell<=exp[B+omega ||d_ell||inf/delta].           (16)

Under the reference nu_h, d_ell is a smooth Gaussian whose H3 covariance
is trace class. Its linear norm exponential moments are finite, as in
CMG-T4's displayed Gaussian product proof. Finite Hermite images under
W are polynomials in finitely many Gaussian labels. Holder's inequality
and (16) give every finite accepted POVM moment on this dense span.
The closed multiplication operator by sqrt(f_L), composed with the
isometry W, is closed; its squared norm is (14). The closed positive
form representation theorem yields exactly (15).

Write u=u_0+u_perp. From the full counterterm in (3),
N_u=exp(u_0-h)N_perp, with u_0~N(0,2h) independent of N_perp and v.
The random min q_ell has a continuous conditional logarithmic law.
Thus nu_h(min q_ell=delta)=0. Cutoff convergence from T2 then preserves
the accepted/rejected outcome almost surely. For all N,
||H_ell P_Nv||inf<=C_emb ||H_ell v||H3; use (16) to dominate the
accepted lengths and their finite moments, including each Hermite weight.
This proves the cutoff statements.

These exponential moment estimates are for the specified dense domain
under the reference measure. They do not give finite length means in
every normal state. A conditional mean among accepted runs also requires
the acceptance probability and is generally a state-dependent ratio.

## WCG-T4. PCR phases are visible and recover classical material geometry

Use exactly PCR-T5's separate phase-pulsed coherent family

    phi_h=exp[beta_(s_*)/(2h)-||s_*||_C²/(4h)],
    Phi_(h,b)=exp[-iM(b)/h]phi_h,                         (17)

where s_*,b are fixed smooth real profiles. Set rho_*=exp(s_*),
j_*=epsilon rho_* and d_*=b rho_*. The joint outcome has

    E N_u(w)=integral w rho_*, E v(k)=integral k d_*,
    Var N_u(w)=integral integral w(z)w(y)rho_*(z)rho_*(y)
                                    [exp(2hK(z-y))-1],
    Var v(k)=(h/2)||k||_C²+Var_(phi_h) M(kb).            (18)

The configuration marginal u is Gaussian with mean s_* and covariance
2hC. The full joint law and the conjugate marginal are generally not
Gaussian after this nonlinear phase pulse.

For every fixed ell, the complete geometry POVM distinguishes the two
states Phi_(h,b),Phi_(h,-b) whenever b is not identically zero, although
all their density statistics coincide. This is a statement about laws
of repeated readings, not single-shot state reconstruction.

Fix h0<1/2. There are constants C_rho,C_d depending on s_*,b,h0 but
not on h<=h0 or ell, such that, with

    A(ell)=1+2 sum_(n>=1)(1+n²)^3 exp(-2ell²n²),
    B(ell)=1+4 sum_(n>=1)n(1+n²)^3 exp(-2ell²n²),

the exact error bounds are

    E||q_ell-H_ell rho_*||H3² <= C_rho h A(ell),
    E||d_ell-H_ell d_*||H3² <= (h/2)B(ell)+C_d h A(ell). (19)

For 0<ell<=1, explicit sufficient bounds are

    A(ell)<=2085 ell^-7, B(ell)<=8283 ell^-8.             (19a)

Thus the density error is at most 2085 C_rho h ell^-7 and the
conjugate error at most (8283/2+2085 C_d)h ell^-8. If

    ell(h)->0, h ell(h)^(-8)->0,                        (20)

the outcome profiles r_ell,m_ell converge in probability in H3, hence
C2, to the original classical inverse of j_*,d_*. The coframe, spatial
metric and circumference converge in probability as well. Every bounded
continuous geometry reading has its corresponding classical limiting
expectation. For fixed delta<min rho_*, acceptance probability tends to
one and accepted length converges in probability to classical length.
No unbounded mean-length convergence is asserted. One explicit admissible
choice is ell=h^(1/9); the displayed exponents are sufficient, not optimal.

**Proof.** The density phase in (17) changes no configuration law.
The noisy configuration marginal therefore has the stated Gaussian
law. Cameron--Martin shifting in (3) gives N_u with law rho_*N_0,
where N_0 has parameter 2h. Equation (8) proves the first variance.
PCR-T5 gives the current mean and variance; T1's independent probe
noise adds h||k||_C²/4, giving the second variance in (18).

The geometry outputs contain q_ell and d_ell: from r,m compute P,
vartheta=r-log|m/P|, j=m a(vartheta), d=j Y(vartheta). Since H_ell
has no zero Fourier multipliers, H_ell(b rho_*) is not identically zero
for b!=0. Some smeared d_ell mean therefore changes sign between the
two states, and is finite by (18). Their geometry laws must differ,
so some bounded Borel geometry reading distinguishes them.

For each real normalized Fourier mode e_j, ||e_j||inf<=sqrt(2).
The kernel K(r)=1-log(2|sin(pi r)|) is positive and

    exp(2hK)-1 <= 2h K exp(2h0 K).

The right side divided by h is integrable since 2h0<1. For example take

    C_rho=4 exp(2||s_*||inf) integral K(r)exp(2h0 K(r))dr,
    C_d=2||b||inf² exp(2||s_*||inf)
                                integral K(r)exp(h0 K(r))dr.

Equation (18) gives uniformly Var N_u(e_j)<=C_rho h and
Var v(e_j)<=h/(2lambda_j)+C_d h. Summing the nonnegative expected
squared Fourier coefficients with their H3 weights proves (19);
cross-mode independence is not required. For n>=1, (1+n²)^3<=8n^6.
For explicit constants put L=ceil(1/ell), so 1/ell<=L<=2/ell. Below
L, sum n^6<=L^7 and sum n^7<=L^8. Above L, the exponential-series
inequality exp(-x)<=6!/x^6 gives

    (1+n²)^3 exp(-2ell²n²)<=90 ell^-12 n^-6,
    n(1+n²)^3 exp(-2ell²n²)<=90 ell^-12 n^-5.

Integral tail bounds are sum_(n>L)n^-6<=1/(5L^5) and
sum_(n>L)n^-5<=1/(4L^4). Including the respective paired factors
2 and 4 and the zero mode gives
1+2048+36=2085 and 1+8192+90=8283 in (19a). These deliberately
conservative constants prove the rates without numerical heat evaluation.

Markov's inequality, (19),(20) and smooth heat approximation of rho_*,d_*
give H3 convergence in probability. The inverse (9) is continuous in
H3 near these smooth profiles with min rho_*>0 and nonzero P_*. Sobolev
algebra/composition and H3 embedding into C2 give the claimed geometry
limit. This reuses PCR's actual density phase states, without assigning
them DCD's separate positive-curvature preparation or a Gaussian current
law. Variance control supplies probability and bounded-test limits; it
does not provide uniform integrability for unbounded length.

## WCG-T5. Fixed-h point and ungated-length obstructions remain explicit

At each fixed 0<h<1/2, in the normalized LMD vacuum, for every ell>0
the ungated circumference is finite almost surely but has infinite mean:

    E L_ell=infinity.                                   (21)

As ell decreases to zero, the raw point r reading is not tight. More
precisely, at any point z and finite R>=0,

    Pr(|r_ell(z)|<=R)
      <=2(R+B+log3)/(omega sqrt(2pi) sigma_ell) ->0,
    sigma_ell²=Var d_ell(z)
      =(h/2)[1+4 sum_(n>=1)n exp(-2ell²n²)].              (22)

These are obstructions to the declared raw point/ungated prescription,
not to every possible renormalized or distributional quantum geometry.

**Proof of (21).** In the decomposition used in T3, put
Q_perp=p_ell*N_perp and d_ell=v_0+d_perp. The latter is smooth Gaussian,
with v_0 independent of d_perp and all u. Some finite M>1 gives positive
probability to max Q_perp<=M and ||d_perp||inf<=M. Intersect with the
independent positive-probability event epsilon v_0>=M+1. Uniformly in z,

    y_ell(z)>=exp(h-u_0)/M,
    L_ell>=exp[-B+(omega exp(h)/M)exp(-u_0)].

Its expectation over u_0~N(0,2h) is infinite. Tonelli gives (21).
As delta decreases to zero, the accepted length increases pointwise
to L_ell, so its vacuum mean diverges by monotone convergence.

**Proof of (22).** Write T=N_u(S1). Since 1/2<=a<=3/2 and the heat
kernel preserves mass,

    2T/3<=|P_ell|<=2T,
    r_ell(z)=omega epsilon d_ell(z)/q_ell(z)
                +log(q_ell(z)/T)+E_ell(z),
    |E_ell(z)|<=B+log3.                                 (23)

Condition on the entire N_u. In the vacuum d_ell(z) is independent of
N_u, centered Gaussian with variance in (22), obtained from the paired
real Fourier modes. The event |r_ell(z)|<=R requires this Gaussian to
fall into an interval centered at
-epsilon q_ell(z)log(q_ell(z)/T)/omega, of length
2q_ell(z)(R+B+log3)/omega. Although E_ell can depend on d, its uniform
bound makes this interval a necessary condition. Bound the Gaussian
density by 1/(sqrt(2pi)sigma_ell), then use E q_ell(z)=1. This proves
(22). The variance diverges, in fact sigma_ell²~h/(2ell²), by comparison
with integral_0^infinity x exp(-2ell²x²)dx=1/(4ell²).

## Exact closure and next obligation

The PCR phase/current bridge now has an actual nonlinear material
measurement in the same LMD sector: a normalized continuum POVM with
exact Wick-density first moment, exact original-current moment domain,
smooth resolved inverse, bounded metric operators and an accepted
self-adjoint circumference mean. Its cutoff and classical-resolution
limits have explicit proofs. The generic Gaussian detector, positive
form method and classical inverse are inherited, not newly discovered.

The measured measure N_u has Gaussian parameter 2h. First-moment
calibration does not make its entire outcome law covariant under the
system's corrected local flow: the sample Wick transformation at 2h
would use correction 1+h, whereas the system retains 1+h/2. The chosen
filter and detector remain stated measurement data, not a canonical
map or a new material-flow symmetry theorem.

The next source obligation is a compatible material-flow transformation
law for the calibrated detector and resolution, or a justified alternative
composite prescription addressing the fixed-h point obstruction. Neither
this material-slice measurement nor its classical limit supplies a
Lorentzian spacetime, gravitational field equations, or the full 3+1
quantum-gravity closure.
