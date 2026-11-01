# MDD-1 — Source lift, joint detector dilation and calibrated readings

16 October 2026. Frozen parent: CMG-1 at
3bec523b8b89a7c98b204da012920d0f8b2dea47. Retain the selected RGR/MPQ
physical space, covariance, source chart and sign superselection. The
old J-pointer action is unchanged. This packet constructs a separate,
explicit apparatus extension and states exactly what it adds.

The finite-mode mixing and joint readout are the established Gaussian
heterodyne construction, rederived below with this source's signs and
normalizations. The contribution is its source lift, countable physical
realization, calibration and compatibility boundaries, not discovery of
heterodyne measurement. See SOURCE_MAP.md for attribution.

## 1. Contract and the old apparatus's information boundary

Write h=hbar>0 and H_epsilon=L^2(mu_h), mu_h=N(0,h C), with MPQ's real
Fourier basis and lambda_n=(1+n^2)^(-4), including lambda_0=1 and the two
real modes at each n>=1. In this section j labels a real Fourier mode.
The source coordinates and physical operators are

    s=log|j_density|, d=j_density Y(theta),
    S_j=s_j, D_j=ih partial_(s_j)-i s_j/(2 lambda_j),
    [D_j,S_k]=ih delta_jk, J(f)=epsilon integral f exp(s). (1)

D_j has the self-adjoint domain specified by MPQ's translation group,
with h restored. The common polynomial core suffices for all displayed
commutator calculations; unitary groups specify the closed operators.

Define the J-only apparatus class to have independent prepared probes,
arbitrary probe-only controls/readouts, and system-probe interactions
whose unitaries are decomposable in s. In particular this includes any
finite sequence of frozen exp[-i alpha J(w) tensor Pi/h] pulses, and
their outcome-conditioned choices, with classical processing of the
complete outcome record. It excludes a new system D-control and excludes
initial system-probe correlations encoding the unknown state.

**MDD-T1 (phase blindness of the frozen interaction class).** Every
effect of such a protocol is multiplication by a bounded function of s
in each sign sector. Therefore two states with the same |psi(s)|^2
have identical outcome laws for every protocol in this class. The
RMG/CMG joint adapted-coordinate POVM cannot be obtained from it by
classical postprocessing alone.

Proof. On H_epsilon tensor K_probe a decomposable unitary acts as U(s)
on the probe fiber. Taking matrix elements in the fixed probe state
after any probe effect produces a scalar measurable multiplier in s.
Probe-only operations, products, conditioning on earlier records, sums
and integrals over records preserve decomposability. The resulting
effects therefore lie in the multiplication algebra. This also covers
limits in which the effects converge weakly, since that algebra is
weakly closed. It does not classify protocols with other interactions.

For a concrete separation use psi_0=1 and
psi_b=exp[-i b s_0/h], b!=0. Both are normalized and have identical
configuration laws. RMG's overlap formula, or T3 below, gives v_0 laws
N(0,h/2) and N(b,h/2), respectively. Their characteristic functions
differ by exp(itb). Thus some bounded v_0 reading distinguishes them,
while no J-only record can do so. This is an information obstruction,
not an assertion that the original J instrument is inconsistent.

## 2. A canonical probe extension and its actual source lift

For finitely many real Fourier modes through N introduce independent
canonical apparatus pairs (a_j,b_j) with {b_j,a_k}=delta_jk. Their
one-form is -sum b_j da_j. Let s_j,d_j be the exact reduced physical
coordinates (1), without truncating the old constraints, and set

    G_N=sum_(j through N) (s_j b_j-d_j a_j),
    C_T=E+g(T) G_N, {T,E}=1.                             (2)

Here g is a real smooth compactly supported pulse. Its area is an
angle theta, with no identification as a physical energy or light
travel time. These are global mode probes, not a local apparatus field.

**MDD-T2 (finite first-class source lift with feedback).** Pull s,d in
(2) back through the exact MPQ canonical chart. On its regular source
patch, adding (2) and the apparatus one-form to the MRF action yields
a classical first-class extension. Reduction gives exactly (2), and
the coupled physical/probe flow at pulse area theta is

    (s_j,a_j) -> (c s_j+t a_j, -t s_j+c a_j),
    (d_j,b_j) -> (c d_j+t b_j, -t d_j+c b_j),
    c=cos(theta), t=sin(theta).                          (3)

Modes beyond N are unchanged. The system backreaction is part of this
flow; it is not a measurement on a fixed profile. All finite pulse
areas preserve the sign sector of the reconstructed smooth reduced
geometry. Off-shell patch preservation for every gauge representative
is not claimed.

Proof. MPQ-T1 gives canonical pairs (c_m,m), (phi,K_phi), (Z,R_Z) with

    H_tilde=p K_phi+u Z' R_Z, D_tilde=v K_phi+Z' R_Z,
    v=phi', u=v/p.                                      (4)

The reduced s,d depend only on c_m,m and satisfy {d_j,s_k}=delta_jk.
The lift of G_N is independent of phi,Z,K_phi,R_Z, so it strongly
commutes with the two gauge momenta. Its brackets with (4) are exactly

    {G_N,H_tilde[N_0]}=integral N_0({G_N,p}K_phi
                                      +{G_N,u Z'}R_Z),
    {G_N,D_tilde[f]}=integral f({G_N,v}K_phi
                                      +{G_N,Z'}R_Z).     (5)

The coefficient functionals are smooth on the same regular chart.
K_phi,R_Z are invertible smooth linear combinations of the old local
constraints there. Thus (5) is in their constraint ideal for arbitrary
smearings. The old brackets are MRF-T1, {C_T,C_T}=0, and the new clock
and apparatus pairs commute with them. This proves first-class closure
and conditional preservation of the common zero set within the chart.

The extended action is the MRF old-field canonical action plus
E Tdot-sum b_j adot_j and multiplier term -nu C_T, replacing the old
pointer extension for this experiment. Reduce the old local families
using MPQ, set T as clock and eliminate E. The physical one-form is
integral m delta r= -sum d_j delta s_j plus an exact differential;
all physical modes remain, while only the displayed finite subset
couples. Hamilton's equations give sdot_j=g a_j, adot_j=-g s_j,
ddot_j=g b_j, bdot_j=-g d_j, and integrate to (3).

Use j_density=epsilon exp(s), y=d/j_density,
theta_source=Y^{-1}(y), m=j_density/a(theta_source),
P=integral m and r=theta_source+log|m/P| to reconstruct the system.
These are exactly the frozen inverse formulas, with the global P
feedback included. Smooth finite-mode changes in s,d remain smooth;
j_density,m,P keep their sign. The representative phi'=0,Z=x stays
in the regular patch. General off-shell representatives need not do so.

No pointwise infinite sum G_infinity is required as a classical source
functional on rough quantum samples. T3 constructs its reduced quantum
generator by a unitary group. This does not specify quantum orderings
for every nonlinear unreduced field in (4).

## 3. Exact physical dilation and countable limit

Prepare a probe copy H_a=L^2(mu_h) in its vacuum Omega_a=1. It has
A_j=a_j and B_j=ih partial_(a_j)-i a_j/(2lambda_j). The probe covariance
is matched to the selected system covariance; this is preparation data.
Let V psi=psi tensor Omega_a. The sign label is unchanged throughout.

For a finite set of modes define on H_epsilon tensor H_a

    (U_(theta,N) Psi)(s,a)
       =Psi(c s-t a, t s+c a) on those modes,             (6)

with the other coordinates unchanged. The two equal Gaussian factors
make this pullback unitary. Its generator is the self-adjoint closure
of Ghat_N=sum(S_j B_j-D_j A_j) on finite Hermite polynomials, and
U_(theta,N)=exp[-i theta Ghat_N/h]. The Gaussian score terms cancel
in each S_j B_j-D_j A_j. Its Heisenberg action is exactly (3).

**MDD-T3 (the RMG/CMG POVM is realized exactly).** At theta=pi/4 measure
the commuting output quadratures sqrt(2) S_j on the system and
-sqrt(2) B_j on the probe. Before interaction their pulled-back readings
are the commuting self-adjoint families

    U_j=S_j+A_j, V_j=D_j-B_j.                            (7)

The resulting system POVM through N is exactly RMG-T1, including its
du dv/(2pi h) normalization and the identity on unused system modes.
There is a strongly continuous countable rotation U_theta with
U_(theta,N) -> U_theta strongly, uniformly on compact theta intervals.
Its full output POVM is exactly CMG-T1, for every normal physical state.
Ghat_N -> Ghat in strong resolvent sense, where Ghat is the Stone
generator of U_theta; this is not a formal assertion of convergence
of unbounded sums on arbitrary vectors.

Proof. Strong commutativity in (7) follows by conjugating the output
PVM, whose system S_j and probe B_k strongly commute. In particular
[S_j,D_k]=-ih delta_jk cancels [A_j,-B_k]=+ih delta_jk.
Commutators alone would not have established the joint PVM.

To identify it, pass in one mode to Lebesgue coordinates by multiplying
by the square roots of the Gaussian densities. The probe wavefunction is
eta(a)=(2pi h lambda)^(-1/4) exp[-a^2/(4h lambda)]. At balance the
postinteraction wavefunction is psi((s-a)/sqrt(2)) eta((s+a)/sqrt(2)).
The B eigenfunctions are exp[-ika/h]/sqrt(2pi h). Set u=sqrt(2)s,
v=-sqrt(2)k. Changing integration variable to x=(s-a)/sqrt(2) shows
that its joint density, including ds dk=(du dv)/2, is

    |integral eta(x-u) exp[ivx/h] psi(x) dx|^2
                                               du dv/(2pi h). (8)

The vector eta(x-u) exp[-ivx/h] is precisely RMG (5); the harmless
overall phase drops out. Tensor products establish the finite identity.

For all modes, (6) defines a measure-preserving transformation of the
countable product. Cylinder functions are dense. If a function depends
on modes through M, U_(theta,N) equals U_theta on that function for all
N>=M and every theta. Unitarity gives the claimed uniform strong limit
on all vectors and strong continuity. The Stone generator has domain
{Psi: lim_(theta->0) ih(U_theta Psi-Psi)/theta exists in norm}.
Projections onto finitely many modes and bounded total Hermite degree
commute with the rotations and increase strongly to I. Their finite
dimensional ranges give a common graph core and the stated polynomial
formula for the generator. Convergence of the groups gives strong
resolvent convergence (also by integrating the resolvent kernel).

The countable joint output PVM exists: configuration multiplication on
one factor and the commuting momentum PVM on the other. In the latter,
Hermite phase rotation diagonalizes each B_j as a scaled Gaussian
coordinate; its vacuum has independent variances h/(4lambda_j) and is
cyclic. Thus the scaled joint output vacuum law is precisely kappa_h
from CMG (1), and the output vacuum is cyclic. All its normal-state
spectral laws are therefore absolutely continuous with respect to
kappa_h. Compressing by U_(pi/4)V gives the desired system POVM.
Its finite cylinder marginals are RMG's by (8). Cylinder sets determine
all scalar matrix-element measures, so it equals CMG-T1 everywhere.

The ancillary probe and output quadratures are extra resources. A
countable ideal dilation is not a finite hardware device, a proof of
spacetime locality, or a realization by the old J-only interaction.

## 4. Clock norm, geometry outputs and induced instruments

**MDD-T4 (relational physical apparatus and cutoff instruments).** Let
alpha(T)=integral_(-infinity)^T g, W(T)=U_(alpha(T)) and let mathcal W
multiply H_ext=H_red tensor H_a-valued histories by W(T). Then

    C_quant=mathcal W(-ih partial_T)mathcal W^*,
    Dom C_quant=mathcal W H^1(R_T;H_ext)                  (9)

is self-adjoint and on its transported spectral core equals
-ih partial_T+g(T)Ghat. On the test space
mathcal W[S(R_T) tensor_alg D_G], where D_G is the union of bounded
spectral subspaces of Ghat, group averaging exp[it C_quant/h] gives

    Rf=integral W(T)^*f(T)dT,
    integral dt <f,exp[it C_quant/h]k>=<Rf,Rk>_Hext.      (10)

The physical completion is H_ext. Preparation V is an isometry of the
existing physical system into it. The same statements hold at finite N,
where (9) is the specified quantization after the classical reduction
in T2. RGR's selected gauge quotient is unaffected: these operations
act on its physical factor, with the identity on its gauge-test factor.

After balance, apply the same CMG filter R_ell=(I+ell^2 Lambda)^(-4),
inverse chart and rejection-inclusive threshold map F_(ell,delta) to
the commuting outputs. This gives exactly CMG's geometry POVM and its
accepted-length form/operator. Its finite-bin readout also has an
explicit apparatus instrument, not just an unspecified CP extension.

Proof. Unitary conjugation proves (9) with its domain. Removing W from
the test space reduces (10) to integration of ordinary clock translations.
Choose a Schwartz function of integral one to see dense range D_G.
The null quotient and completion give H_ext. This is the already proved
RPM/MPQ clock construction with the present Ghat, not an independent
assumption that kinematic and physical norms agree.

Let P_out be the joint PVM in T3, and let {B_b} be a finite partition of
the coarse output space, including rejection. Define projections
Q_b=P_out(F_(ell,delta)^(-1)(B_b)), K_b=Q_b U_(pi/4)V, and

    I_b(rho)=Tr_a[K_b rho K_b^*].                        (11)

These are normal completely positive operations, sum_b K_b^*K_b=I,
and their effects are exactly E_(ell,delta)(B_b). They retain the
disturbed system and discard the probe after the output projection.
They need not equal CMG's separately selected square-root instruments:
an effect does not uniquely specify the postmeasurement state.

For finite N use U_(pi/4,N), replace filtered output profiles by their
Fourier projections through N, and use the corresponding Q_(b,N).
Assume every bin boundary has kappa_h-null inverse image under F.
CMG-T2/T3 prove almost-sure convergence of F_N to F, including the
null threshold boundary. Cyclicity and absolute continuity in T3 imply
Q_(b,N) -> Q_b strongly on the whole extended Hilbert space. Hence
K_(b,N) -> K_b strongly; all are contractions. Rank-one approximation,
the trace norm inequality for outer products, and contraction of
partial trace prove I_(b,N)(rho) -> I_b(rho) in trace norm for every
trace-class input. This proves the cutoff instrument statement.

The ideal quadrature PVM is a declared terminal readout. Neither (9)
nor (11) derives a microscopic local detector or resolves every
unreduced quantum constraint ordering. The finite source lift, the
reduced quantum limit and their precise distinction are retained.

## 5. Calibration, resolution-independent readings and their limits

For real smooth f,g write C(f,f)=sum lambda_j f_j^2 and
||g||_C^2=sum g_j^2/lambda_j. Define u(f),v(g) by distributional
pairing with the raw CMG outputs u in H^3, v in H^-5. The smeared
momentum below extends to g in H^4 by the Gaussian/Segal limit; the
pointwise pathwise pairing claims here use smooth g.

**MDD-T5 (noise and source-density calibration).** The joint POVM has
the following bounded operator characteristic function:

    E[exp(i[u(f)+v(g)])]
       =exp[-h C(f,f)/2-h ||g||_C^2/8]
                                      exp(i[S(f)+D(g)]). (12)

The last exponential is the Weyl operator of the self-adjoint Segal
field sum, not a product of two commuting system observables. On the
common polynomial domain, and as the corresponding quadratic moments,

    E[u(f)]=S(f), E[v(g)]=D(g),
    E[u(f)^2]=S(f)^2+h C(f,f) I,
    E[v(g)^2]=D(g)^2+(h/4)||g||_C^2 I.                  (13)

The square-integrability domains of the two recorded readings are
exactly Dom S(f) and Dom D(g). Their added noise variances obey

    n_u(f)n_v(g)=(h^2/4)C(f,f)||g||_C^2
                      >=(h^2/4)|integral fg|^2.         (14)

Equality holds for one matched mode. For any fixed 0<theta<pi/2,
rescaling the two outputs by 1/c and -1/t instead gives U=S+(t/c)A,
V=D-(c/t)B. The mode noises become h lambda tan^2(theta) and
h cot^2(theta)/(4lambda); their product remains h^2/4. Balanced mixing
is the particular calibration giving exactly CMG, not the only possible
joint detector. This is the added-noise product in this construction,
not a claimed optimality theorem for all measurement schemes.

There are also exact resolution-independent records. Before the geometry
threshold is applied, retain the adapted readings (or their complete
filtered profiles). For smooth f,g and any ell>0,

    u_ell((I+ell^2 Lambda)^4 f)=u(f),
    v_ell((I+ell^2 Lambda)^4 g)=v(g).                    (15)

Thus these probability laws and calibrated moments do not depend on
ell. Without deconvolution, readings with the same fixed smearing instead
measure R_ell f,R_ell g. Discarding rejected records invalidates the
unconditional equalities; the linear records must be retained even
when the nonlinear geometry result is rejection.

Finally put c_C=C(z,z)=1+2 sum_(n>=1)lambda_n, independent of z. For
smooth f>=0, the positive recorded density

    j_cal(f)=integral f(z) exp[u(z)-h c_C/2] dz           (16)

has exactly the positive first-moment form of epsilon J(f). Equivalently,
epsilon j_cal(f) is unbiased for the original signed J(f). Its form
domain and associated positive self-adjoint mean operator are precisely
those of epsilon J(f); equality includes possible infinite expectations
outside the form domain. It is a noisy measurement of that mean,
not its spectral measurement. It can be recovered from complete filtered
configuration data by (15), or directly from the retained raw u record.

Proof. In the dilation (7), u(f)+v(g) is the sum of a system Segal field
S(f)+D(g) and its commuting probe field A(f)-B(g). Probe vacuum has
zero means, variances h C(f,f), h||g||_C^2/4, and zero symmetrized cross
covariance. Its Weyl expectation is the Gaussian scalar in (12).
Finite-mode identities pass to smooth smearings in Gaussian variance
and strongly in the unitary groups, proving (12). The one-variable
marginals are therefore the spectral laws of S(f) and D(g) convolved
with their stated independent Gaussian noises. Tonelli and the scalar
identity E[(x+Z)^2]=x^2+Var Z give (13), including exact second-moment
domains. Cauchy--Schwarz applied to sqrt(lambda_j)f_j and
g_j/sqrt(lambda_j) gives (14); (3) gives the unbalanced formulas.

The filter is self-adjoint, with nonzero multiplier
t_n=(1+ell^2 n^2)^(-4). Its inverse maps smooth test functions to smooth
test functions, so (15) holds as distributional duality, pathwise on
CMG's full-measure regularity set and hence for every normal state.
This is not a bounded inversion theorem for arbitrary detector errors:
on the n-th mode the inverse gain is (1+ell^2 n^2)^4, unbounded in n.
No fixed spatial resolution or finite precision yields arbitrary point
recovery from this algebraic identity.

For (16) the u marginal is exactly the system configuration s plus an
independent Gaussian probe field a of covariance h C. Both fields have
continuous representatives. Conditioning on s, Gaussian integration
gives E_a exp[s(z)+a(z)-h c_C/2]=exp[s(z)]. Positivity and Tonelli
justify integration against f and arbitrary |psi(s)|^2, even when the
answer is infinite. Polarization on the finite form domain identifies
the forms; the multiplication operator epsilon J(f) is their associated
operator. Gaussian exponential integrability gives a dense domain
containing the finite Hermite polynomials. For a configuration-point
reading, the conditional second moment is exp[2s(z)+h c_C], displaying
the residual noise rather than claiming a sharp measurement.

These results calibrate a known observable family and smooth adapted
smearings. They do not reconstruct a sharp continuum r or remove the
inverse-density tail in the length. CMG-T5 still proves non-tight raw
point r at fixed h and infinite ungated vacuum mean circumference.
CMG's ell,delta remain declared nonlinear-geometry detector settings;
the present construction selects neither them nor C or h from gravity.

## 6. Closure and next obligation

Closed under the stated contracts: the J-only apparatus information
boundary; a finite first-class source lift of the additional joint
probe; its exact RMG/CMG quantum measurement, countable rotation and
positive relational clock norm; induced finite-bin instruments with
controlled cutoffs; exact added-noise and original J-density calibration;
and resolution-independent smooth adapted records on the full record.

Still open: a local material detector model selecting these resources;
calibration from independently specified dynamics or experiment; a
resolution-independent nonlinear metric/curvature observable with a
justified domain, or a larger obstruction theorem; quantization of all
original nonlinear source constraints in this extended apparatus;
physical matter/potentials and physical 3+1 gravity. No residual-Diff
uniqueness result is assumed or developed here.

The next geometry target is the nonlinear observable algebra: determine
whether source-geometric quantities beyond J admit a calibrated
distributional or smeared continuum definition, with physical domains,
compatible dynamics and controlled classical recovery. Exact linear
deconvolution alone does not establish it.
