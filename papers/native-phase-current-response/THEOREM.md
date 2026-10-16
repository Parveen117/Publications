# PCR-1 — Phase and current complete the density-response coordinates

16 October 2026. Frozen parent: `7664681c901f4d063b09e0d021d1d74b0a615ca1`.

## Source contract and result

Keep the selected LMD Hilbert sector L²(mu_h), 0<h<1, its density M,
the QRT windows X_i=M(w_i), and a fixed real smooth bounded-mass
preparation psi_* from DCD-T1. Write T=M(S1) for total mass, not time.
Let k_i be smooth real profiles equal to one on supp(w_i) and zero
on the other window support, and retain the original currents

    P_i=D(k_i),  i=1,2.

They are actual Cameron--Martin Stone generators. In particular they
are not finite matrices assigned to a finite spectral comparison model.

For real theta,beta in R² define

    b_beta(z)=sum_i beta_i w_i(z),
    V_beta=exp[-i M(b_beta)/h]=exp[-i beta.X/h],
    psi_(theta,beta)=V_beta Z(theta)^(-1/2)exp(-theta.X/2)psi_*.
                                                               (1)

The phase unitary is an existing density pulse. The modulus is QRT's
same normalized family. This adds two phase controls without adding a
Hilbert-space factor, a representation or a primitive spacetime dimension.
There are four real controls and two complex tangent directions.

The exact observer data are

    m_i=<X_i>, p_i=<P_i>,
    boxed: p_i=beta_i m_i,   beta_i=p_i/m_i.                (2)

These four means distinguish the entire declared four-control family.
In q_i=log(m_i),p_i, its Berry curvature has the exact Darboux form

    boxed: F_Berry=(1/h)sum_i dq_i wedge dp_i.              (3)

q_i means log of a density expectation, not expectation of log density
and not a new sharp operator. A separate phase-pulsed coherent family
gives controlled smeared classical material recovery. The positive
curvature preparation and that semiclassical family have distinct roles.

## PCR-T1. Actual current readout and closed operator domains

For any smooth real k,v, the exact conjugations are

    V_beta* D(k) V_beta = closure[D(k)+M(k b_beta)],
    V_beta* R_(a_h)[v] V_beta
                       = closure[R_(a_h)[v]-M(v b_beta')],
    a_h=1+h/2.                                             (4)

The closures in (4) start on the bounded-mass graph cores defined below.
Their domains are respectively V_beta*Dom D(k) and
V_beta*Dom R_(a_h)[v]. No equality with an unrestricted algebraic
sum domain is assumed. In particular (4) gives

    V_beta* P_i V_beta=closure[P_i+beta_i X_i].              (5)

Every state (1), and its required parameter derivatives, lies in these
current domains. Its density statistics are beta-independent. Since
its unphased amplitude is real, <D(k)>=<R_(a_h)[v]>=0 before the pulse.
Consequently (2) holds, and also

    <R_(a_h)[v]>_(theta,beta)=-<M(v b_beta')>_theta.         (6)

The P_i commute strongly with each other, but P_i and its own X_i do
not. Their means can be measured on separately prepared ensembles;
these formulas do not assert a joint sharp density/current measurement.

**Proof and domains.** LMD's Cameron--Martin density law differentiates
M(b) along k to M(kb); its corrected material flow differentiates the
window by v b'. Acting on exp(-iM(b)/h) gives the positive sign for D(k)
and the negative sign for the material current in (4), agreeing with
the inherited density-feedback identity.

For A=D(k) or R_(a_h)[v], put

    C_A=union_(L<infinity) [Dom A intersect 1_{T<=L}H].      (7)

This is a graph core. For R, the mass projections reduce its unitary
group. For D(k), approximate f in Dom D(k) by chi(T/L)f, with smooth
chi equal to one below one and zero above two. The commutator is

    [D(k),chi(T/L)]=i h chi'(T/L) M(k)/L.

Its norm is uniformly bounded, using |M(k)|<=||k||_infinity T, and
it tends strongly to zero since it is supported in the tail T>=L.
The product rule therefore gives graph convergence to f.

On any bounded-mass support, V_beta times a slightly larger smooth
mass cutoff is a multiplier with norm-differentiable current orbit.
For D(k), the common support enlargement under translation is bounded
using exp(-|t|||k||)T<=M(e^(tk))<=exp(|t|||k||)T. For the material
current total mass is conserved. The first derivatives of the phase
are bounded on these supports. The Stone product rule proves that
V_beta and its inverse map each bounded-mass current domain onto
itself and gives (4) there. Hence C_A is also a core for V_beta* A V_beta;
the displayed sums close to exactly that unitary conjugate.

DCD's preparations and derivatives have this bounded-mass smoothness.
The real unphased vector has a purely imaginary current action; its
self-adjoint current expectation must also be real, so it is zero.
Equation (5) now gives (2). Strong commutation of the P_i follows from
the genuine abelian Cameron--Martin translation representation. The
affine relation [P_i,X_j]=i h delta_ij X_j explains why the own-channel
joint sharp reading is not supplied by this result.

## PCR-T2. Quantum response tensor and its exact normalization

Let Y_i=X_i-m_i and H_ij=<Y_iY_j>, the unchanged QRT covariance.
With inner product conjugate-linear in its first argument, the horizontal
state derivatives are

    nabla_(theta_i)psi=-Y_i psi/2,
    nabla_(beta_i)psi=-iY_i psi/h.

In the control order (theta_1,theta_2,beta_1,beta_2), their quantum
geometric tensor and its real metric are

    Qgeom = [[H/4,       iH/(2h)],
             [-iH/(2h),  H/h²]],
    g=Re Qgeom=diag(H/4,H/h²).                              (8)

g is positive definite on the four real controls. Qgeom itself has
complex rank two, not four; the two beta derivatives are complex
multiples of the theta derivatives. The complex structure

    J=[[0,-2I/h],[hI/2,0]]                                 (9)

satisfies J²=-I and J^T g J=g. The manifold is Kähler, with convention

    omega(u,v)=g(Ju,v),
    omega=(1/(2h))sum_ij H_ij dtheta_i wedge dbeta_j.

For the Berry convention A_Berry=i<psi,dpsi>,

    A_Berry=(1/h)sum_i m_i dbeta_i,
    F_Berry=dA_Berry=-(1/h)sum_ij H_ij dtheta_i wedge dbeta_j
                   =-2 omega.                             (10)

The factor and sign in (10) are part of the statement. The positive
four-control metric is not a Lorentzian 3+1 metric.

**Proof.** Differentiate (1) and subtract the component along psi.
Taking the four blocks of the resulting Gram matrix proves (8).
H>0 follows from DCD's actual spectral law. Formula (9) is direct
matrix multiplication. Put z_i=theta_i/2+i beta_i/h. The unnormalized
section exp(-z.X)psi_* is holomorphic, since X is bounded on its support.
Its logarithmic squared norm is u(z+conj(z))=log Z(theta), with
partial_(z_i)partial_(conj z_j)u=H_ij. This gives the stated Kähler
metric; alternatively the Hessian symmetry directly proves d omega=0
and (9) is the integrable constant complex structure. Differentiating
the normalized section gives (10).

The covariance H, its original TVSP response connection and its QRT
curvature are unchanged by beta. In particular the selected DCD state's
1/8<kappa<3/8 statement still refers to H. It is not a claim about the
scalar curvature of the whole four-control metric. The coefficients
of (8) alone remain insensitive to beta's value; the current readings
in T1 supply that missing coordinate information.

## PCR-T3. Four observed means identify the family and give Darboux coordinates

The map (theta,beta)->(m_1,m_2,p_1,p_2) is globally one-to-one and a
smooth diffeomorphism onto m(R²) x R². Its Jacobian determinant is

    det(H) m_1 m_2 > 0.                                   (11)

The mean image m(R²) is open; it is not asserted to be the whole
positive quadrant. In q=log m,p coordinates, (3) holds exactly.
Thus the phase/current completion turns the declared observer data
into canonical coordinates for the explicitly signed Berry form.

**Proof.** Since partial_(theta_j)m_i=-H_ij, for theta!=eta,

    (theta-eta).(m(theta)-m(eta))
      =-integral_0^1 (theta-eta)^T
           H(eta+t(theta-eta))(theta-eta) dt <0.

The mean map is injective; its nonsingular derivative gives a smooth
local inverse, and these inverses agree globally on its open image.
The measured positive m fixes theta and then beta_i=p_i/m_i fixes
beta. In blocks the readout Jacobian is

    [[-H,              0],
     [-diag(beta) H, diag(m)]],

which proves (11). Finally

    F_Berry=(1/h)sum_i dm_i wedge dbeta_i
      =(1/h)sum_i (dm_i/m_i) wedge dp_i
      =(1/h)sum_i dq_i wedge dp_i.                         (12)

This is identification within the prescribed family with a known
preparation and windows. Four means do not tomographically determine
an arbitrary LMD state or all original material profiles. These are exact
ensemble means, not one-shot state recovery. Positive H gives local
inverse conditioning on compact parameter patches; no uniform stable
inverse over all theta or finite-sample precision bound is asserted.

For comparison with the source sign convention, using i_X Omega=dE
for Hamiltonian evolution requires Omega=-h F_Berry=sum_i dp_i wedge
dq_i. The density Hamiltonian expectation E=m_i=exp(q_i) then gives
dot p_i=m_i, exactly the positive momentum shift in (5). This check
fixes the sign; it does not identify every source Hamiltonian with a
flow confined to this finite control manifold. No finite-dimensional
matrix realization of the full current commutation algebra is implied.

## PCR-T4. Source lift, mixed-loop phase and what it measures

Let K_theta(t,r) be DCD's unitary along the theta path. Define

    L(t,r)=V_beta(t) K_theta(t,r) V_beta(r)*.                (13)

This is an exact unitary propagator carrying the base state family (1).
If E_beta=V_beta E_theta is its density-response frame, then

    L(t,r)E_beta(r)=E_beta(t)U_TVSP(t,r).                    (14)

Its instantaneous self-adjoint generator is

    H_control(t)=dot(beta).X+i h V_beta B_theta V_beta*,
    Dom H_control(t)=Dom(dot(beta).X),                      (15)

where the density combination uses its joint-spectral realization.
The second term is bounded self-adjoint. The differential equation
i h partial_t L=H_control L is proved on the invariant common domain
Dom X_1 intersect Dom X_2, and hence on the density polynomial domain
intersection. It is not asserted on every changing maximal domain
Dom(dot(beta).X) without transverse density regularity.

For each fixed original current A=D(k) or R_(a_h)[v], L and its inverse
preserve Dom A intersect Dom T. This includes the prepared states.
An unrestricted phase pulse is not claimed to preserve every full
current domain independently of total-mass integrability.

The Berry phase for a control rectangle with order
+eta theta_i,+zeta beta_j,-eta theta_i,-zeta beta_j is exactly

    gamma_B=(zeta/h)[m_j(theta+eta e_i)-m_j(theta)]
       =-(zeta/h)integral_0^eta H_ji(theta+t e_i)dt.         (16)

For i=j and eta,zeta>0 this is strictly negative. Sufficiently small
such loops have nontrivial line holonomy. Its small-area coefficient
is -H_ii eta zeta/h. This is a mixed amplitude/phase loop; the Berry
curvature on the theta-only slice is zero even when the response-plane
holonomy of QRT is nonzero.

**Proof.** Differentiate the unitary product (13) on the common density
domain and use DCD's bounded generator. Bounded self-adjoint perturbation
proves (15)'s instantaneous domain. Equation (14) follows from the QRT
intertwiner by conjugation. V_beta commutes with T. Its current orbit
difference quotients are bounded pointwise by a constant times T, so
the Stone product rule applies on Dom A intersect Dom T. DCD's K
preserves Dom A, is identity off a fixed bounded-total-mass support,
and preserves Dom T. This proves the current-domain assertion.
Integrating (10) along the rectangle gives (16); positivity of H_ii
proves the small-loop claim without an unknown third-moment calculation.

The explicit L follows the chosen normalized section and returns its
base state with zero total phase on a closed control loop. Its expected
generator is dot(beta).m, so the dynamical phase -integral A_Berry
cancels the geometric phase +integral A_Berry. A horizontal line lift is

    psi_horizontal(t)=exp[i integral_0^t A_Berry]psi(t).

It has zero inner product with its derivative and returns the phase
exp(i gamma_B). Distinguishing that phase experimentally requires a
coherent reference and dynamical-phase calibration; the pulse product
alone does not automatically produce a measurable net Berry phase.

In particular (14) is a conjugated source transport identity. It does
not claim that L is Kato parallel transport for the full moving complex
response projection under arbitrary beta motion.

## PCR-T5. Controlled material information from the phase pulse

For this theorem use LMD's separate real coherent amplitude centered at
a fixed smooth classical profile s_*, with zero initial conjugate phase:

    phi_h(s)=exp[beta_(s_*)(s)/(2h)-||s_*||_C²/(4h)],
    Phi_(h,b)=exp[-i M(b)/h]phi_h,

where b is any fixed smooth real function (in the two-channel family,
b=sum beta_i w_i). For smooth real k,v and bounded real w,

    <M(w)>=integral w exp(s_*),
    <D(k)>=integral k b exp(s_*),
    <R_(a_h)[v]>=-integral v b' exp(s_*).                   (17)

The actual current variances are

    Var D(k)=h||k||_C²/4+Var_phi M(kb),
    Var R_(a_h)[v]=h||v s_*'+a_h v'||_C²/4
                    +h²||S_v||_HS²/2+Var_phi M(vb').       (18)

Here ||.||_C and S_v are LMD's Cameron--Martin norm and whitened
current's symmetric Hilbert--Schmidt part. Every variance in (18)
tends to zero as h->0. The density variance does likewise. These are
separate spectral concentration statements for the smeared observables;
no joint sharp law of noncommuting density and momentum is postulated.

Their classical fields are

    j_*=epsilon exp(s_*), d_*=b exp(s_*), y_*=epsilon b.

The inherited inverse chart therefore gives the smooth classical labels

    vartheta=Y^(-1)(epsilon b), a(vartheta)=1+cos(vartheta)/2,
    m_mat=j_*/a(vartheta), P_tot=integral m_mat,
    r=vartheta+log|m_mat/P_tot|,
    ds_material²=exp(2r) dz².                              (19)

The sign epsilon matters because the pulse in (1) uses positive X=M(w),
not the signed J=epsilon M(w).

**Proof.** Under |phi_h|²mu_h the field has the shifted Gaussian law
of LMD, so its total density has a finite second moment. Current-orbit
phase difference quotients are dominated by a constant times T.
The Stone product rule and dominated convergence put Phi_(h,b) in
the actual current domains and give (4) acting on phi_h. The unphased
coherent current action is purely imaginary, while every density
multiplier is real. The real cross terms in the squared norms vanish.
LMD's coherent variances and density covariance then give (17),(18).
Explicitly, for a real test f,

    Var_phi M(f)=integral integral f(z)f(y)exp(s_*(z)+s_*(y))
                         [exp(h K(z-y))-1] dz dy,

which tends to zero by the inherited integrable logarithmic bound.
Chebyshev proves the claimed spectral concentration. Finally
R[v]=integral d(v'+vs_*')=-integral v b'exp(s_*) for d=b exp(s_*),
and substitution into the source inverse chart proves (19).

An explicit one-channel homogeneous specialization takes positive
orientation, w=k=1, s_*=log j_0 and b=beta. For every 0<h<1,

    <T>=j_0, <D(1)>=beta j_0,
    beta=<D(1)>/<T>,
    r=Y^(-1)(beta), length=exp(r).                          (20)

This is a calibrated reconstruction of the classical family label from
exact quantum means. All material-current means are zero in this case
because b'=0; the Cameron--Martin momentum supplies the missing signal.

Equations (19),(20) do not define a finite-h sharp operator D(1)/T,
a pointwise quantum r or a quantum material metric. The semiclassical
coherent family in this theorem is not the curvature-enclosed DCD cell
preparation, and its TVSP curvature has not been assigned that enclosure.
The general classical source recovery in LMD is inherited; the new
statement is its realization by these density phase pulses, their
phase/current identification and their connection to the response model.

## Closure and next obligation

The density-response family now has phase-sensitive current readings,
global identification by four means, an explicit Kähler response tensor
and observed Darboux coordinates. Its exact source lift and mixed-loop
Berry phase retain their distinct dynamical and geometric meanings.
A controlled coherent comparison reaches the missing conjugate material
information and gives an explicit homogeneous reconstruction.

The next geometric task is an actual conjugate-sensitive quantum
observable for the nonlinear material inverse, with a stated resolution
and continuum/domain limit. The displayed response manifold supplies
neither a Lorentzian signature nor gravitational field equations, and
the original pointwise LMD inverse and full 3+1 source remain open.
