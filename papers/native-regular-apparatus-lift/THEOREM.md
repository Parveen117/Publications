# RAL-1 — Regular-source apparatus lift and fixed-window obstruction

16 October 2026. Frozen parent: RPM-1 at
5c75d47d51262a1382b8893f9b818bc4e57adb43. Retain the periodic, real,
V=0 matter-free LNC fields and LPS/GOA conventions, with hbar=1.

RPM added a clock and pointer after reduction. This packet gives an explicit
canonical action **before the old local constraint reduction**, proves its
classical first-class algebra, includes its geometry backreaction, and
recovers RPM exactly on the regular P branch. The apparatus still consists
of global canonical pairs. Its coupling contains a global momentum; it is
not a spatially localized apparatus field or a new local matter action.

## 1. Regular currents and an exact off-shell observable lift

Use x in R/Z, {q(x),p(y)}={phi(x),k(y)}=delta(x-y), and

    A_+=p+phi', A_-=p-phi', B_+=k+q', B_-=k-q',
    T_+=(A_+ B_+-A_+')/2, T_-=(A_- B_-+A_-')/2,
    P=integral p, c(x)=q-(1/2)log|A_+ A_-|.                  (1)

Work on the open smooth periodic off-shell patch A_+ A_->0 everywhere,
with both currents nowhere zero. Their integrals coincide, so their
constant signs agree and P!=0. The inequalities and inverses are actual
domain conditions. CCS proves that every smooth on-shell P!=0 solution
is in this patch; no such implication is asserted off shell.

Let a(c,P) be any real smooth function for P!=0, and define

    G_a=integral p(x) a(c(x),P) dx.                          (2)

The coordinate P in the second argument is the integral of p, not a fixed
external parameter. Subscripts c,P denote partial derivatives in these
two arguments. G_a has smooth Hamiltonian gradients on this patch.

**RAL-T1 (explicit local constraint-ideal bracket).** Put

    Z=(A_-/A_+)T_+ +(A_+/A_-)T_-.

For every smooth periodic smearing f, the following are exact off-shell
identities, with delta_epsilon F={F,T_epsilon[f]}:

    c'=T_+/A_+ - T_-/A_-,
    delta_epsilon c=f T_epsilon/A_epsilon, delta_epsilon P=0,
    {G_a,T_+[f]}={G_a,T_-[f]}=(1/2)integral f a_c Z.         (3)

In particular G_a strongly commutes with D[f]=T_+[f]-T_-[f]
and is weak Dirac under both local families. On shell c is spatially
constant and G_a=P a(c,P).

Proof. The frozen LPS current transformations are

    delta_epsilon A_epsilon=epsilon(f A_epsilon)',
    delta_epsilon A_(-epsilon)=0,
    delta_epsilon q=(f B_epsilon+epsilon f')/2.

Differentiating (1) gives the first line of (3), with no constraint imposed.
Integration of the current derivative gives delta P=0. Since
delta_epsilon p=epsilon(f A_epsilon)'/2, integration by parts gives

    delta_epsilon G_a
      =integral f a_c[-epsilon A_epsilon c'/2
                                +p T_epsilon/A_epsilon].    (4)

There is no a_P contribution because delta P=0. Substitute c' from (3)
and p=(A_++A_-)/2. Both signs give (a_c/2)Z. The boundary term is a
derivative of a smooth periodic function and integrates to zero. This is
an actual linear combination of the local constraints with smooth
field-dependent coefficients on the patch, stronger than just checking
that the bracket vanishes at homogeneous solutions.

For the GOA/RPM probe choose

    a(c,P)=1+(1/2)cos(c+log|P|),
    a_P=a_c/P.                                             (5)

LPS's homogeneous coordinate is Q=c+log|P|. Thus on every regular
constraint orbit, not just the homogeneous representative,

    G_a=P[1+(1/2)cos Q]=F_class.                             (6)

Equation (2) is a new off-shell representative of the same regular
physical observable. GOA's globally smooth Q_sharp representative is
unchanged. No equality of their off-shell functionals is claimed.
Neither log|P| nor a nonlinear function of Q is evaluated numerically by
the certificate; the analytic chain rule supplies (5).

## 2. The unreduced canonical extension and its backreaction

Introduce global canonical pairs (T,E), (X,Pi), with {T,E}={X,Pi}=1,
commuting with the old fields. Choose the same real compact smooth pulse
g(T) as RPM and add the single constraint

    C_T=E+g(T) Pi G_a.                                      (7)

The proposed action, now on the unreduced old field phase space, is

    S_ext=integral ds [integral(p dot q+k dot phi)dx
                 +E dot T+Pi dot X
                 -T_+[n_+]-T_-[n_-]-nu C_T].               (8)

All three multipliers are independent; n_+,n_- are fields and nu is global.
Their primary momenta can be added as the usual vanishing constraints.
The old local generators themselves are not altered.

**RAL-T2 (first-class lift with retained feedback).** On the stated patch
the secondary constraints of (8) obey

    {T_epsilon[f],T_epsilon[h]}=epsilon T_epsilon[fh'-hf'],
    {T_+[f],T_-[h]}=0,
    {C_T,T_epsilon[f]}=(g Pi/2)integral f a_c Z,
    {C_T,C_T}=0.                                           (9)

They are first class for arbitrary admissible local smearings, not only
in a zero-mode truncation. Their Hamilton equations preserve their common
zero set while a smooth solution remains in the patch. The apparatus
contributions to the old field equations are, writing D=A_+ A_-=p^2-phi'^2
and J_P=integral p a_P dx,

    dot q|_app=nu g Pi [a-p^2 a_c/D+J_P],
    dot p|_app=-nu g Pi p a_c,
    dot phi|_app=0,
    dot k|_app=nu g Pi partial_x[p a_c phi'/D].              (10)

The apparatus equations are

    dot T=nu, dot Pi=0, dot X=nu g G_a,
    dot E=-nu g' Pi G_a.                                  (11)

G_a is constant along the extended flow on the constraint surface; P
generally is not. In particular dot P|_app=-nu g Pi integral p a_c.

Proof. The first two brackets are the frozen LNC/LQA classical algebra.
The global variables commute with those constraints, so (3) proves the
third bracket exactly; the fourth is antisymmetry for one scalar
constraint. Thus time derivatives of every constraint are linear
combinations of constraints. This proves conditional preservation, not
global PDE existence or preservation of the off-shell patch for every
unreduced representative.

For (10), the first variation is obtained from

    delta c=delta q-p delta p/D+phi' delta phi'/D,
    delta G_a=integral [p a_c delta q
                   +(a-p^2 a_c/D+J_P)delta p
                   +p a_c phi' delta phi'/D]dx.             (12)

Integrate the last term by parts and use the canonical brackets. The
J_P term is the response of the global P inside a. Omitting it changes
the generator and the reduced q equation. Equation (11) follows by varying
the new canonical pairs. Since {G_a,C_T}=0 and its old gauge variations
are given by (3), G_a is conserved on shell. The apparatus really changes
geometry; it is not a reading imposed on a frozen source.

The action is a selected extension with an exact consistency check. Its
new global clock, pointer and g remain inputs. First-class closure neither
selects g nor makes the apparatus spatially local. This result concerns
classical constraints; it does not define inverse quantum current fields.

## 3. Exact regular reduction and recovery of RPM

**RAL-T3 (reduction and quantization-after-reduction).** Reduce (8) first
by its old local gauge families, using the LPS current straightening and
the old scalar-clock condition. On the homogeneous slice

    p=P, q=Q, k=0, phi'=0, phi=constant,
    c=Q-log|P|,

the canonical one-form and remaining constraint are exactly

    P dQ+E dT+Pi dX,
    C_T,red=E+g(T) Pi P[1+(1/2)cos Q].                      (13)

Gauge fixing the new clock T=s and eliminating E yields

    S_red=integral dT [P Qdot+Pi Xdot-g(T)Pi P a_0(Q)],
    a_0(Q)=1+(1/2)cos Q.                                   (14)

The corresponding equations are

    Qdot=g Pi a_0(Q), Pdot=-g Pi P a_0'(Q),
    Xdot=g P a_0(Q), Pidot=0.                              (15)

For all finite pulse area and finite initial Pi, the homogeneous classical
solution exists globally in T and remains in its initial P sign. In terms
of y(Q)=integral_0^Q du/a_0(u), it is specified exactly by

    y(Q(T))=y(Q(T_0))+Pi integral_(T_0)^T g,
    P(T)=P(T_0) a_0(Q(T_0))/a_0(Q(T)),
    X(T)-X(T_0)=P(T_0)a_0(Q(T_0)) integral_(T_0)^T g.         (16)

Thus the classical pointer shift and its geometry feedback come from the
same source action. They are not separately assigned evolution rules.

Proof. RAL-T1 makes G_a a function on the old gauge quotient; the global
canonical pairs are invariant under those old gauge transformations.
LPS gives the pullback P dQ. This proves (13) without truncating before
checking the local algebra. Since {T,C_T}=1, the new clock is a regular
global gauge parameter; solving for E gives (14). Ordinary differentiation
proves (15). On the slice J_P=a_c by (5), so (10) gives Qdot=g Pi a_0,
consistent with (14). Without that feedback term it would incorrectly
give g Pi(a_0-a_0'). The homogeneous slice is preserved by the apparatus
part of (10); old local multipliers may be set to zero there.

Because 1/2<=a_0<=3/2, y is an increasing diffeomorphism of R onto R.
Differentiating (16) gives (15); P a_0 is constant and the sign cannot
change. These statements concern the homogeneous reduced representative;
they do not assert that every inhomogeneous representative stays inside
the inverse-current chart under zero compensating local multipliers.

Finally choose the **same** reduced ordering and full-line matching as
GOA, on H_P=L^2(dP/|P|). Then the quantization of P a_0(Q) is exactly its
self-adjoint F on Dom P. Quantize the pointer by Pi=-i partial_X and the
clock by E=-i partial_T. RPM's joint-spectral product F tensor Pi, unitary
W(T)=exp[-i(integral g)F tensor Pi], and

    C_T,quantum=W Pi_T W^*, Dom C_T,quantum=W H^1            (17)

are precisely the quantization of the reduced extension (13). RPM-T1
then supplies its positive group-averaged norm and completed physical
space H_P tensor L^2(pointer); RPM-T2 through T5 supply the readout,
instrument, phase contrast, variance and coherent-record interface.
The identification is the identity in these inherited coordinates and
norms; no small-pulse or finite-mode approximation is involved.

This is a demonstrated classical reduction followed by the specified
quantization, not a theorem that quantizing all unreduced extended
constraints gives the same theory. The original LQA anomaly cancellation
is frozen and valid for its old generators; quantum commutators involving
an unreduced C_T still require an operator construction. The classical
action does not choose the inherited full-line quantum matching, whose
orientation mixing differs from the classical sign-preserving flow.

This source lift covers H_P. RPM's zero-action extension to H_K and the
origin is still an inherited quantum choice, not derived from this patch.
GOA's global smooth observable values and COS's two scales are unchanged.

## 4. A precise obstruction to fixed-window localization

A natural attempted localization replaces (2) by

    G_(a,w)=integral w(x) p(x) a(c(x),P) dx,                 (18)

where w is a prescribed real smooth periodic function, not a dynamical
reference field and not varied by the old local gauge generators.

**RAL-T4 (fixed-window obstruction in this ansatz).** Its exact bracket is

    {G_(a,w),T_epsilon[f]}
      =(1/2)integral f w a_c Z
           -(epsilon/2)integral f A_epsilon a w'.           (19)

For the chosen a in (5), G_(a,w) is weak Dirac under all old local
constraints if and only if w is constant. Thus a nonzero fixed window
supported in a proper spatial interval fails this coupling test.
With C_(T,w)=E+g Pi G_(a,w), the unchanged constraint set is not first
class whenever g Pi!=0 and w is nonconstant.

Proof. Repeat (4), now integrating the derivative of w f A_epsilon a.
The additional w' term gives (19). On the old constraint surface c and P
are constant, Z=0, and 1/2<=a<=3/2. If w' is nonzero, select at a chosen
regular constrained configuration the admissible smearing
f=A_epsilon w'. The residual is

    -(epsilon a/2) integral A_epsilon^2 (w')^2 !=0.          (20)

Choose E=-g Pi G_(a,w) to satisfy also the new constraint. The bracket
still does not vanish there when g Pi!=0, proving failure of first-class
closure on the full common zero set. Conversely w constant reduces to
RAL-T1. A smooth nonzero function with support in a proper interval cannot
be constant on the whole circle.

For an explicit periodic witness take the homogeneous solution Q=0,P=2,
so a=3/2, with w=1+(1/2)cos(2 pi x) and f=sin(2 pi x). The plus residual
is 3 pi/4 and the minus residual -3 pi/4. The native Fourier calculation
records the coefficients after division by the common 2 pi, namely
+3/8 and -3/8. This normalization is not a change of the source period.

This is not a no-go theorem for dynamical material rods, extra compensating
fields, gauge fixing or all local matter actions. A covariantly transported
smearing alone also does not specify a canonical apparatus: its conjugate
variables and constraints must be supplied. The obstruction identifies the
missing spatial-reference dynamics in the fixed-background-window attempt.

## 5. Scope and next obligation

New closure: an explicit regular unreduced canonical apparatus action,
its arbitrary-smearing classical first-class algebra, source feedback,
and exact recovery of RPM after reduction with inherited quantization.
New obstruction: a fixed nonconstant spatial window breaks the unchanged
constraint set. Neither statement is inferred from finite test counts.

The coupling remains global through its apparatus variables and P inside
a(c,P). The regular current logarithms and bracket coefficients are not
defined at current zeros. GOA already gives another globally smooth weak
observable representative; this chart failure does not refute that result.

Next: introduce a dynamical spatial reference/apparatus field and derive
its own constraints so that its variation cancels the window residual,
then determine whether a controlled reduction recovers (8) or (14).
Quantum anomaly, domains and physical norms for that unreduced extended
field theory remain a new obligation. Matter, potentials, critical-branch
source coupling, physical 3+1 recovery and calibration are still open.
