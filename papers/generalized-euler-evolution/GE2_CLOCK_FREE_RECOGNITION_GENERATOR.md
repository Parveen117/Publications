# GE2: clock-free records, recognition memory and the native heat generator

Monty Dabas. 3 October 2026. **Python 3.12 only.**

This continuation starts with finite Recognition–Seam arrows, not a
primitive time parameter. For the declared symmetric, independent record
protocol below, it proves on the native polynomial core

\[
\boxed{\lim_{\gamma\to0}\frac{S_\gamma-I}{h_\gamma}
 =\frac12\lim_{\gamma\to0}\frac{\Delta_\gamma P}{h_\gamma}
 =-L_C.}                                                     \tag{1}
\]

Here the record mean is S, the genuine typed Recognition form difference
is Delta P, and h is computed from the arrow's native turn coefficients.
The equality is a limit theorem under stated hypotheses; Delta is not
renamed as a single-space derivation. Irregular refinement produces the
existing native heat law. Deterministic refinement instead produces the
existing phase derivation. Their exponentials have different algebraic roles.

**Scope.** Written proofs plus exact finite controls, not formal
proof-assistant or external expert certification. The native scalar,
quaternion carrier, coefficient pairings, reference Phi and heat completion
are the unchanged YM50–YM52/GE1 constructions. No primitive Hilbert space,
abstract generation theorem, stochastic limit theorem or physical clock
is a premise. Completed positive-pairing spaces are used only after their
native construction. This is a free compact record model, not the
interacting row process or four-dimensional Yang–Mills.

## Primitives and conventions

The source is [T24](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/24_clock_free_recognition_seam_cut_calculus.md):
an arrow A maps one source core to another and transports forms by
Delta_gamma G=A^dagger G_target A-G_source. For composable arrows,

\[
\Delta_{\delta\gamma}G=\Delta_\gamma G+
 A_\gamma^\dagger(\Delta_\delta G)A_\gamma.                   \tag{2}
\]

Daggers below use the declared native positive pairing. Work first on
degree-d homogeneous coefficient spaces with YM50's positive tensor
pairing G_d; the direct sum of degrees is the polynomial core. The native
reference pairing Phi(f^dagger g), after its null quotient, gives the
function completion. Rotations are isometries in both. A native rotation
also preserves the uniform norm. Operator products act right to left;
pullbacks reverse the order of quaternion point maps.

The already derived imaginary words obey
e_a e_b=-delta_ab+epsilon_abc e_c. Write v=sum v_a e_a,
r^2=sum v_a^2. The normalization D_a corresponds to left point motion
x -> Exp_Sigma(-s e_a/2)x, as in YM50–YM51. The native bracket is retained;
no Riemannian curvature is substituted for the framework's commutator
curvature. Constants used below are

\[
B_d=\frac{d^2(d^2+2)}{384},\qquad
K_d=\frac{d^2(2d^2+1)}{96}.                                  \tag{3}
\]

## GE2-T1: rational arrows recover the phase derivation

Define, using native arithmetic alone,

\[
q(v)=\frac{1-r^2/16-v/2}{1+r^2/16},\qquad
U_v f(x)=f(q(v)x).                                          \tag{4}
\]

Then q(v)^dagger q(v)=1 and q(-v)=q(v)^dagger. Thus U_v is an
invertible isometry. On each polynomial degree,

\[
\lim_{s\to0}\frac{U_{sv}-I}{s}=D_v=\sum v_aD_a,
\quad \|D_v\|_d\le d|v|/2.                                 \tag{5}
\]

It is a derivation on products of polynomials. No exponential is needed
to define (4). For a comparison after that definition, set
b(r)=integral_0^r (1+u^2/16)^(-1) du in the earned native scalar calculus.
Direct differentiation of (4), with initial value 1, gives
q(v)=Exp_Sigma(-b(r)v/(2r)) for r>0. At zero use continuity. Moreover
0<=r-b(r)<=r^3/48, by integrating
u^2/(16+u^2)<=u^2/16. Therefore

\[
\|U_{Tv/n}^{\,n}-\operatorname{Exp}_\Sigma(TD_v)\|_d
 \le\frac{d|T|^3|v|^3}{96n^2}.                              \tag{6}
\]

**Proof.** Quaternion multiplication proves the unit and inverse
identities. The linear term of (4) is -sv/2. Substitution into monomials
gives (5), the existing generator, and the ordinary product expansion
gives its Leibniz rule. The degree-d tensor action is the sum of d
degree-one skew actions of norm |v|/2, proving the bound. Same-axis
actions commute; their parameter discrepancy after n steps is at most
|T|^3|v|^3/(48n^2). Integrating the skew evolution gives (6).
These are finite-core arguments with native scalar integration.

We will also use the even-action estimates

\[
\left\|\frac{U_v+U_{-v}}2-I-\frac{D_v^2}2\right\|_d
 \le B_d|v|^4,\qquad
\left\|I-\frac{U_v+U_{-v}}2\right\|_d\le d^2|v|^2/8.          \tag{7}
\]

For completeness, the fourth-order integral remainder of the even
skew exponential costs d^4 r^4/384. Replacing r by b(r) costs at most
(r^2-b(r)^2)d^2/8<=d^2 r^4/192: integrate the even action twice,
using isometry to bound its second derivative. Their sum is B_d r^4.
The second inequality follows by the same double integration and b<=r.
This proves (7) for all r, without a smallness assumption.

## GE2-T2: an actual typed record lift and its memory law

Let records s have positive rational weights p_s with sum one and
isometries U_s. Zero-weight labels are removed. Form the weighted direct
sum K=direct_sum_s H, with pairing sum_s p_s <f_s,g_s>. Define

\[
Jf=(f)_s,\quad J^\dagger(f_s)=\sum_s p_s f_s,\quad
Af=(U_s f)_s,\quad P=JJ^\dagger,\quad Q=I-P.
\]

Here J and A are isometries, P is an orthogonal cut, and S=J^dagger A
is a contraction. Choosing the initial cut P_0=I gives exactly

\[
M=A^\dagger QA=I-S^\dagger S\succeq0,\qquad
\Delta_\gamma P=A^\dagger PA-I=-M.                            \tag{8}
\]

Thus <f,Mf>=sum_s p_s ||U_s f-Sf||^2. This is YM52's record
variance expressed as T24's typed form difference, not a newly discovered
variance identity. Weighted direct sums avoid square roots of weights.

Retain all old labels when appending an independent new record. The
next arrow acts by the same U_s on every old-label component. Write P_n
for the cut to constant record functions at stage n. Then

\[
P_{n+1}A_{n+1}(I-P_n)=0,\qquad
S_{1:n}=S_n\cdots S_1,
\]
\[
I-S_{1:n}^\dagger S_{1:n}
 =\sum_{j=1}^n S_{1:j-1}^\dagger M_j S_{1:j-1}.              \tag{9}
\]

**Proof.** All assertions through (8) follow by expanding finite sums.
For no return, first sum over old labels: its weighted mean is zero on
I-P_n, and the new U_s does not depend on those labels. At one append,
Delta P=-J_n M_{n+1}J_n^dagger. Pull this identity back through the old
history, using J_n^dagger A_history=S_{1:n}, and apply (2) to get (9).
Alternatively telescope the contractions directly. Record independence
and absence of feedback are hypotheses; (9) is not a proof that a
general native observer has no returning memory.

## GE2-T3: an additive history budget and a recognition generator

A symmetric step is a finite list (w_r,v_r), w_r>=0, sum w_r=1,
representing two labelled turns +v_r and -v_r, each of weight w_r/2.
Define from the actual arrow records

\[
Q_\gamma=\sum_r w_r v_rv_r^T,\quad
h_\gamma=\tfrac12\operatorname{tr}Q_\gamma,\quad
m_{4,\gamma}=\sum_r w_r|v_r|^4,\quad
C_\gamma=Q_\gamma/(2h_\gamma).                              \tag{10}
\]

For h>0, C is positive semidefinite of trace one. If h=0 all positive
weight turns vanish and S=I; no normalized shape C is assigned. The
history ledger tau=sum_gamma h_gamma is additive under concatenation.
The factor 1/2 and native coefficient norm are explicit conventions.
This is a clock on retained record histories, not a scalar state function,
not a unique clock, and not a physical-seconds calibration. A path and
its reversed point motion can close at the endpoint while retaining a
positive history ledger; the history arrow has not become the empty word.

For h>0 use the already defined native operator
L_C=-sum_ab C_ab D_aD_b. It satisfies 0<=L_C and ||L_C||_d<=d^2/4.
Equations (7) give

\[
\|S-I+hL_C\|_d\le B_dm_4,\quad
\|I-S\|_d\le hd^2/4,
\]
\[
\left\|\frac{M}{2h}-L_C\right\|_d
 \le B_dm_4/h+d^4h/32.                                    \tag{11}
\]

If all active turns obey |v|<=epsilon, then m_4<=2h epsilon^2,
h<=epsilon^2/2, and the last bound is at most K_d epsilon^2.
For fixed shape C and epsilon->0, equations (8) and (11) prove (1).

**Proof.** Sum (7) with weights w_r and use
sum w_r D_v^2=-2hL_C. Symmetry gives S=S^dagger. If delta=S-I,
M=-2delta-delta^2. Substitute the first estimate and
||delta||<=hd^2/4 to obtain (11). The stated constants reduce to (3).
The generator in (1) is the density of recognition-form loss with
respect to a specified internal ledger on this protocol class. On the
native completion its closure is the existing GE1/YM51 heat generator;
the limit initially asserts a core identity, not convergence on every
vector in an unspecified domain.

## GE2-T4: unequal-step refinement recovers heat evolution

For any finite ordered symmetric history, allow different normalized
shapes C_j and step budgets h_j. Let tau=sum h_j, and epsilon bound
every active turn. Let E_C(h) be the already proved native heat on the
completion, used here as a comparison target, not as an arrow primitive.
Then

\[
\left\|S_n\cdots S_1-E_{C_n}(h_n)\cdots E_{C_1}(h_1)\right\|_d
 \le\sum_j\left(B_dm_{4,j}+\frac{d^4h_j^2}{32}\right)
 \le K_d\epsilon^2\tau.                                    \tag{12}
\]

**Proof.** Twice integrate the finite-core heat equation to obtain
||E_C(h)-I+hL_C||<=h^2||L_C||^2/2<=h^2d^4/32. Combine with
(11) and telescope the two products; every exact step and reference
factor is a contraction. The moment inequalities give the last bound.

For fixed C the reference word is E_C(tau). Therefore histories with
epsilon->0 and tau->T give E_C(T), without equal steps or a primitive
external time. A changing final budget adds at most d^2|tau-T|/4 to
the degree bound. Approximate shapes are also admitted if
sum_j h_j ||C_j-C||_entry1 ->0: finite-core Duhamel integration adds
(d^2/4) sum_j h_j ||C_j-C||_entry1 to (12).

For a fixed polynomial, the finite coefficient norm controls uniform
evaluation. Uniform contraction of all word averages and native heat
then extends convergence to the uniform polynomial completion by
approximating once by a fixed polynomial. The same argument applies to
the Phi-pairing completion, where convergence on a fixed coefficient
core implies convergence in Phi norm. Approximants whose degrees grow
arbitrarily are allowed if they converge in the chosen norm: insert the
fixed polynomial between them and use the two contractions. No uniform
degree-free rate is inferred from (12).

For variable C_j, (12) certifies comparison with the *ordered* reference
word. Refinements of a specified finite piecewise-constant shape profile
converge block by block. Arbitrary evolving protocols do not acquire a
continuum limit from the scalar ledger alone.

## GE2-T5: three necessary distinctions and exact counterexamples

**Phase versus heat.** T1's D_v is a derivation, whereas the generator
G=-L_C reconstructed in (1) obeys the inherited YM52 rule

\[
G(fg)=G(f)g+fG(g)+2\Gamma_C(f,g),\quad
\Gamma_C(f,g)=\sum C_{ab}(D_af)(D_bg).                        \tag{13}
\]

Consequently the Morphic manuscript's derivation-based proof that
Exp(tD) preserves multiplication applies to phase flow, not this heat.
For f=x_0 and C=diag(0,1/2,1/2), Gamma_C(f,f)=(x_2^2+x_3^2)/8
is nonzero. The sign in (13) is positive for G=-L_C.

**Memory does not determine arbitrary motion.** A deterministic unitary
record has M=0 while its phase generator can be nonzero. Also S and -S
have identical I-S^dagger S. Equation (1) requires the symmetric,
near-identity family and its quadratic ledger; it is not a universal
inversion of lost information. Outside the no-feedback append model,
let U be the two-mode rotation with cosine 3/5 and sine 4/5, and P the
first-coordinate cut. Then

\[
PU^{-1}UP=P,\qquad PU^{-1}PUP=\tfrac9{25}P,
\quad PU^{-1}(I-P)UP=\tfrac{16}{25}P.                         \tag{14}
\]

This is the existing T24/N06 returning-memory obstruction in an exact
control, not a resolution of the interacting YM row defect.

**The scalar clock does not erase order.** Put C_0=diag(0,1/2,1/2)
and rotate it by native conjugation with quaternion (3/5,0,0,4/5) to C_1.
On degree two, [L_C0,L_C1] is nonzero. Thus the mixed coefficient of
E_C1(b)E_C0(a)-E_C0(a)E_C1(b) is [L_C1,L_C0], nonzero. Both words
have the same scalar budget a+b and the same integrated tensor aC_0+bC_1.
Native order information survives clock recovery. Exact rational record
words also distinguish these two orders. This uses the existing EMK-T2
ordered-transport principle; it does not define a new curvature invariant.

## GE2-T6: response curvature expressed in the record clock

Take the *existing* YM54 two-direction shape Gram matrix G_response,
with T=tr(G_response)>0 and det(G_response)=16f^2. Avoid confusing
this response trace with the history ledger tau. YM54 gives a distinct
three-dimensional protocol tensor C_YM with trace m=T/4 and nonzero
eigenvalues g_min/4, g_max/4. Normalize C=C_YM/m as in (10). Then

\[
\kappa=\frac{64f^2}{T^2}\in[0,1],\qquad
\operatorname{spec}(C)=
 \left\{0,\frac{1-\sqrt{1-\kappa}}2,
             \frac{1+\sqrt{1-\kappa}}2\right\}.             \tag{15}
\]

The unchanged YM52 all-content centered compact rate in this clock is

\[
\gamma_{\rm rec}=
 \min\left\{\frac14,\frac{1-\sqrt{1-\kappa}}2\right\}.
                                                                    \tag{16}
\]

**Proof.** Solve the response matrix's quadratic characteristic equation
using T and 16f^2, divide its eigenvalues by T, and apply YM52's
already proved rate min(tr C/4,a+b). No new spectral theorem is needed.
For fixed m, tau=m t_old and tau L_C=t_old L_CYM; the change in the
numerical rate is clock normalization, not a new mass gap. Balanced
rank-two shape has gamma_rec=1/4, while diag(0,1,1) has old rate 1/2.

If every block of an ordered free heat word satisfies kappa>=kappa_0>0,
then its centered Phi norm contracts by at most exp(-gamma_* tau),
where gamma_*=min(1/4,(1-sqrt(1-kappa_0))/2). Each block fixes the same
Phi reference and has this bound, so multiply the bounds. Any admitted
strong limit inherits it. Finite rational words have the additional
degree error (12). This reuses the existing compact/ordered estimates;
it neither proves convergence for an arbitrary profile nor replaces the
interacting ground-source weight with the free reference.

## Certification boundary and next obligation

The new bridge is **rational recognition arrows -> retained record
memory -> intrinsic quadratic history ledger -> the established native
phase/heat generators**, with explicit unequal-step errors. Algebraic
curvature, operator calculus, reference pairing, energy/variance and
compact spectral rates were already present and are credited as inputs.

Outstanding choices include the record law, independent fresh labels,
observer cut, symmetric turns, response directions and clock calibration.
General Morphic arrows with feedback require their returning-memory
terms; a singular state clock cannot be repaired by deleting them.
Actual interacting row closure, physical state/protocol selection,
four-dimensional continuum gauge theory, Clay mass gap, RH and quantum
gravity remain open. See [source audit](GE2_SOURCE_AUDIT.md) and the
read-only executable certificate for the exact finite evidence.
