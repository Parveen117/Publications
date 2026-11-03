# Vacuum propagation speed, calibration and an experimental consistency check

Author: Monty Dabas. Development R9, 30 September 2026.

**Result.** The existing common-coframe Einstein–Maxwell–matter action has one
local vacuum characteristic cone. After removing gauge directions it carries
two photon and two classical metric-wave polarizations. Its geometric-optics
speed ratio is exactly c_GW/c_gamma=1. This is the familiar Einstein–Maxwell
result, now checked through the native coefficient realization, and is
consistent with the specified GW170817 comparison below. It is conditional on
the action and common metric already supplied in CP/MG/SM; it is not a new
discriminator against general relativity or a primitive-only numerical
prediction. The native equations do not select a dimensional speed calibration.

This development addresses propagation only. It neither completes general
classical Einstein–matter existence/stability nor derives a quantum graviton,
hbar, alpha, a cosmological constant, or a particle mass.

## Inputs and conventions

We consume unchanged [NP](NATIVE_PROPAGATION_ACTION.md),
[CP](THERMO_GAUGE_COMPLETION.md), [MG](NATIVE_METRIC_DYNAMICS.md),
[LR](NATIVE_LOOP_CURVATURE.md), [SM](NATIVE_MASSIVE_MATTER.md), and
[PF-6](../thermo-phase-field/THEOREM.md). [Source pins](SOURCE_PINS.json) bind
the consumed bytes. R8's anisotropic results remain unchanged.

The spacetime coframe e is nondegenerate, with g=e^T eta e and
eta=diag(-1,1,1,1). This is the declared Lorentzian propagation metric, not a
claim that a four-real-dimensional pseudo-Kähler metric has Lorentzian index.
Use the nonzero Palatini, zero-Holst sector of LR/SM, a nonzero Maxwell kinetic
coefficient, and a nonzero matter normalization. For the coupled linearization
take a torsion-free Einstein vacuum with background Psi=0 and F=0. Lambda may
be nonzero: the explicit LR curved vacua are allowed. No flat global vacuum or
zero-Lambda limit is needed for the local characteristic argument.

At a point let p be a nonzero real covector, p^mu=g^{mu nu}p_nu, and
q=g^{mu nu}p_mu p_nu. A principal symbol retains the highest derivative terms.
We use derivative substitution d_mu -> p_mu; an overall Fourier sign has no
effect on its kernel. Characteristics describe wavefronts and the
geometric-optics limit. They do not determine source emission delays,
finite-wavelength curvature scattering, or quantum/discrete corrections.

## SC-1. The native wave law admits a continuous speed family

NP supplies real matrices Gamma^a with

\[
\{\Gamma^a,\Gamma^b\}=2\eta^{ab}I_4,
\quad A_i=\Gamma^0\Gamma^i=A_i^T,
\quad\{A_i,A_j\}=2\delta_{ij}I_4.
\]

Its native volume multiplier J is invertible, J^T=-J, and [J,A_i]=0.
In dimensionless coordinates (tau,y), **every** v>0 admits

\[
B_v=\partial_\tau-vA_i\partial_{y_i},\qquad
D_v=\Gamma^0\partial_\tau+v\Gamma^i\partial_{y_i},
\]
\[
D_v^2=(-\partial_\tau^2+v^2\Delta_y)I_4,
\qquad
\det(\omega I_4-vA_i k_i)=(\omega^2-v^2|k|^2)^2.
\tag{1}
\]

The Clifford anticommutators prove (1); Gamma^0 D_v=-B_v connects the two
first-order forms. The evolution obeys

\[
\partial_\tau(\psi^T\psi/2)
 -\partial_{y_i}(v\psi^TA_i\psi/2)=0.
\]

For constant nonzero Z, the action Z/2 integral psi^T J B_v psi has Euler
equation Z J B_v psi=0: both J and JA_i are skew, so their first-derivative
operators are formally self-adjoint after integration by parts. These facts
hold for all v>0. The family is symmetric hyperbolic; positivity of its norm,
the variational multiplier and the native phase return do not select v.

With a supplied length/time adapter x=L_0 y, t=T_0 tau,

\[
\boxed{c_{\rm prop}=v\,L_0/T_0.}                         \tag{2}
\]

Setting v=1 chooses normalized native coordinates. It does not determine
L_0/T_0. The phase identity R^2=-I likewise fixes an angular period, while
exp(R omega t) admits any physical rate omega until a clock law is supplied.
Changing units and transforming the standards changes no physics. Holding
external rods/clocks fixed while varying an admitted constitutive coefficient
does change a reported speed. The present construction supplies neither a
selection of that coefficient nor an independent identification of those
standards. This is an underdetermination result, not a proof that all coordinate
representatives of a generally covariant metric describe different physics.

## SC-2. The lattice and SI audits reach the same boundary

The unchanged PF-6 model has

\[
\omega_\tau^2=4\frac{\nu_0}{\epsilon_0}
                   \sum_i\sin^2(q_i/2),\qquad
c_{\rm eff}^2=\frac{\nu_0}{\epsilon_0}\frac{\ell^2}{T_0^2}
\quad(q_i=\ell k_i).
\tag{3}
\]

Here epsilon_0, nu_0, ell, T_0 are supplied positive data. Multiplying nu_0 by
four at fixed phase geometry doubles the long-wavelength speed. The finite
lattice is dispersive; it is not automatically identical to the continuum
action used below, and its phase/group speed is not identically (3) at every
wavevector. A native phase period does not remove these parameters.

The [BIPM definition of the metre](https://www.bipm.org/en/si-base-units/metre)
fixes c=299792458 m/s exactly, with the second supplied by its own definition.
Consequently, choosing L_0/T_0 to reproduce that numerical value is a
calibration, not an independent experimental prediction. A physical theory
can predict a universal causal cone and dimensionless relations among speeds;
there is no demand that pure algebra select the human SI numeral. Equations
(2)–(3) identify the missing constitutive and operational input precisely.

## SC-3. Maxwell characteristics after removing gauge freedom

The CP-4 common-cone constitutive hypothesis and MG-7 action use the same g
as gravity. Form the principal field strength f_mu nu=p_mu a_nu-p_nu a_mu.
Contracting the Maxwell equation gives

\[
P_\gamma(a)_\mu=q a_\mu-p_\mu(p^\nu a_\nu).
\tag{4}
\]

The nonzero scalar kinetic normalization cancels. For every p, a=p chi lies
in the kernel, so the ungauged determinant is identically zero; it cannot
serve as a characteristic test by itself.

If q is nonzero, (4)=0 implies a=p(p.a)/q: the entire kernel is gauge.
If q=0, then (4)=0 is the single independent condition p.a=0. Its kernel has
dimension three, containing the one-dimensional gauge subspace. The physical
quotient therefore has **two** dimensions precisely on the null cone.
In an orthonormal frame with p=(1,n), |n|=1, representatives are
a=(0,u), (0,v), where u,v are independent and perpendicular to n.

## SC-4. Einstein characteristics and the two tensor polarizations

Eliminating the algebraic connection in the stated vacuum gives the ordinary
Einstein principal operator. For a symmetric metric variation h, write
h_tr=g^{mu nu}h_mu nu, s_mu=p^alpha h_alpha mu, and pph=p^mu s_mu. Directly
linearizing the Levi-Civita connection and contracting its curvature gives

\[
P_T(h)_{\mu\nu}=\frac12\left[
p_\mu s_\nu+p_\nu s_\mu-qh_{\mu\nu}
-p_\mu p_\nu h_{\rm tr}
-g_{\mu\nu}(pph-qh_{\rm tr})\right].
\tag{5}
\]

Background curvature, Lambda and the nonzero overall gravitational coupling
do not change this principal kernel. The four independent gauge variations
h_mu nu=p_mu xi_nu+p_nu xi_mu lie in it for every nonzero p. Contracting (5)
with p^mu gives zero, the principal linearized Bianchi identity.

For q nonzero, define C_nu=s_nu-p_nu h_tr/2. Under a gauge variation,
C_nu -> C_nu+q xi_nu. Set C=0 by choosing xi=-C/q. Equation (5) then reads
-q(h_mu nu-g_mu nu h_tr/2)/2=0, so h=0 in that gauge. The original kernel
is exactly the four gauge directions: the ten-dimensional symbol has rank six.

For q=0, all nonzero null covectors are related by a Lorentz transformation
and nonzero scaling. It suffices to take p=(1,0,0,1). Equation (5) becomes
the four independent conditions

\[
h_{00}-2h_{03}+h_{33}=0,\qquad
h_{13}-h_{01}=0,\qquad h_{23}-h_{02}=0,\qquad
h_{11}+h_{22}=0.
\tag{6}
\]

Thus rank P_T=4, dim ker P_T=6, and quotienting the four-dimensional gauge
image leaves **two** polarizations. Representatives have h_0mu=0, spatial
trace zero, and h_ij n_j=0. For transverse equal-length orthogonal u,v they are
u_i u_j-v_i v_j and u_i v_j+v_i u_j. These span the physical quotient. This
orbit argument proves the statement for all nonzero null p, rather than
inferring it from randomly sampled covectors. The modes are classical tensor
waves; their existence is not a quantization or a graviton construction.

## SC-5. Native matter has the same wavefront but a massive group delay

On the doubled real SM module,

\[
P_\psi(p)=\gamma^a e_a{}^\mu p_\mu,\qquad
P_\psi(p)^2=q I_8,\qquad \det P_\psi=q^4.
\tag{7}
\]

It is invertible off the null cone and has rank four on it. The latter follows
by the same Lorentz/scaling reduction and the rank-two null symbol on each
real-four copy. Mass, the background spin connection and the algebraic contact
term are lower derivative order. At background Psi=0 the contact term also
vanishes in the linearization. None changes (7).

For clarity, a local constant-coefficient massive mode can use the existing
real native phase Ical, with Ical^2=-I and [Ical,gamma^a]=0, instead of positing
a new complex scalar. Substitution of exp[Ical(-omega tau+k.y)] into
(gamma^0 d_tau+v gamma^i d_yi-m)Psi=0 gives

\[
\left[\mathcal I(-\omega\gamma^0+vk_i\gamma^i)-mI\right]
\left[\mathcal I(-\omega\gamma^0+vk_i\gamma^i)+mI\right]
=(\omega^2-v^2|k|^2-m^2)I_8.
\]

The positive-frequency branch has omega^2=v^2|k|^2+m^2 and

\[
v_{\rm group}=\frac{v^2|k|}{\sqrt{v^2|k|^2+m^2}}\leq v.
\tag{8}
\]

For m>0 and finite momentum the inequality is strict. A massive mode's phase
velocity may exceed v; this does not move its characteristic wavefront. The
parameter m here is the uncalibrated native equation parameter, not a measured
particle mass in kilograms and not a derivation of hbar.

## SC-6. Conditional common-speed theorem and its countermodel

In the stated Einstein vacuum the matter stress, current and spin source
start at quadratic order in Psi, and Maxwell stress is quadratic in F.
They have no linear sources about Psi=F=0. Conversely, variations of e and
the connection multiply the zero background matter field in its linearized
equation. The algebraic connection has no additional propagating torsion mode
in this sector. Thus SC-3–SC-5 are the decoupled vacuum principal blocks of
the existing coupled equations, not an assumption that arbitrary nonzero
matter backgrounds have already been analyzed.

Their physical characteristic set is the same q=0. In the orthonormal frame
of any one timelike observer the photon and metric wavefronts therefore obey

\[
\boxed{r=\frac{c_{\rm GW}}{c_\gamma}=1,\qquad
\delta=\frac{c_{\rm GW}-c_\gamma}{c_\gamma}=0.}
\tag{9}
\]

Any common positive time/length calibration cancels from r. A nondiagonal or
anisotropic coframe can change coordinate speeds without separating the two
local cones. The Maxwell, matter and Palatini overall coefficients do not
affect (9) provided they remain nonzero. This theorem concerns the continuum
vacuum and its geometric-optics interpretation; no quantum, medium or
finite-lattice dispersion prediction is being appended to it.

**Why the shared metric must remain a stated hypothesis.** On the same
native coefficient algebra take two independent Maxwell sectors with positive
kinetic energy and inverse metrics diag(-1,1,1,1) and diag(-1,4,4,4). The covector
(1,1,0,0) is null for the first and non-null for the second; their coordinate
speeds relative to a fixed common clock are 1 and 2. Each coframe also admits
the NP native propagation action. The native
matrix identities alone do not forbid this pair. CP's same-metric choice
does forbid it. The pair is a counterexample to inferring universal speed
from the algebra alone, not a counterexample to (9) with its hypotheses.

## SC-7. A published observation checks the ratio, not the SI numeral

[Abbott et al., ApJL 848 L13 (2017), section 4.1, equation (1)](https://doi.org/10.3847/2041-8213/aa920c)
report for GW170817/GRB170817A

\[
-3\times10^{-15}\leq\delta\leq 7\times10^{-16}.
\tag{10}
\]

Their interval uses a conservative 26 Mpc distance, a 1.74 +/- 0.05 s arrival
lag, and intrinsic gamma-ray emission endpoints of 0 and 10 s after the GW
signal. These emission assumptions matter; (10) is not a model-independent
measurement of exact equality. The paper already identifies equal speeds as
the standard minimally coupled Einstein–Maxwell prediction.

Substituting the model result delta=0 into the published interval gives
**consistent with the stated bound**. No speed parameter was adjusted in this
comparison. It is a retrospective consistency check of a recovered known
prediction, not an independent new prediction distinguishing this construction
from GR, a fit to raw event data, or a measurement of exact equality to fifteen
decimal places. It is not presented as the latest or strongest available bound.

For an elementary fixed-distance, nondispersive propagation comparison,

\[
\Delta t_{\rm arrival}=\tau_{\rm emission}
 +D\left(\frac1{c_\gamma}-\frac1{c_{\rm GW}}\right)
=\tau_{\rm emission}+\frac{D}{c_\gamma}\frac{\delta}{1+\delta}.
\tag{11}
\]

At delta=0 the arrival lag is the intrinsic source lag. Its value is not
calculated by this vacuum principal-symbol theorem. Equation (11) explains
the distinction; it is not a new cosmological propagation or merger-emission
model. [Reference data](SPEED_REFERENCE_DATA.json) keep external definitions
and the published comparison separate from model parameters. The
[primary paper](https://arxiv.org/abs/1710.05834) and its
[institutional PDF](https://ntrs.nasa.gov/api/citations/20170011356/downloads/20170011356.pdf)
provide the provenance.

## SC-8. Verification, outcome and the remaining physical selection

[Exact controls](speed_calibration.py) use rational arithmetic only. They
construct Maxwell's symbol from its field strength and independently compare
the Einstein symbol with a linearized-Christoffel/Ricci calculation on every
symmetric metric basis direction. Three coframes, including an off-diagonal
one, cover twelve null and nine non-null covectors; gauge ranks, Bianchi
identities, transverse representatives, matter squares and determinants are
checked. The family v=1/2,1,3 retains action/norm compatibility and common
characteristics. Massive on-shell controls separate group and front speeds.
The countermodel checks failure of a common cone when its metric hypothesis
is removed. These finite controls accompany the general proofs above; their
quantity does not substitute for independent review or experimental evidence.

| Question | R9 result |
|---|---|
| Native wave speed in a selected chart | v, with c_prop=v L_0/T_0 after a supplied adapter |
| Unique primitive-only dimensional c | Not derived; the current inputs leave the constitutive and calibration family open |
| Local vacuum metric/light speed ratio | Exactly 1 under the existing common-coframe action |
| Specified observational comparison | Consistent with the published, emission-dependent interval |
| New empirical distinction from GR | Not established |
| Quantum graviton or hbar | Not derived |

A stronger next result needs an independently specified native process that
selects the constitutive coefficients and relates them operationally to clocks
and rods, or a distinct dimensionless speed/dispersion correction fixed before
comparison. In the lattice route that would require selecting its weights and
spacing/time relation; equations (3) alone leave them free. Using the target c
to set those inputs would repeat calibration. The undeveloped action-selection
and quantum questions remain separate. A successful consistency check is
useful, but it alone establishes neither a unique fundamental action nor a
solution of quantum gravity.
