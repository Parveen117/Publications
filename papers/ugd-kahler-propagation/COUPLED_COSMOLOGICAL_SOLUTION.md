# An explicit coupled native matter–gravity solution

Author: Monty Dabas. Research development R6, 30 September 2026.

This note constructs an analytic, self-consistent solution of the classical
action in [SM-1–SM-8](NATIVE_MASSIVE_MATTER.md). The solution includes a massive
eight-real-component matter field, a spatially flat FLRW metric, and a nonzero
independent Lorentz connection with algebraic torsion. It lies in the positive-
Lambda native return sector of [LR](NATIVE_LOOP_CURVATURE.md), at finite
supplied loop scale and readout normalization. Every first-order coframe,
connection and matter equation is retained.

The new result is a solved family within the declared model. It does not
select that action, determine physical constants, fit the observed universe,
or establish general existence for arbitrary Einstein–matter initial data.
Homogeneous classical spinor cosmology and torsion-induced self-interactions
have established precedents, credited below. Here their coefficients and the
complete first-order equations are tied to the native module and the source-
bound R5 action, without changing its signs or introducing a new potential.

## CS-1. The cubic equation and an invariant rest sector

Retain SM's real matrices, orientation, signature (-,+,+,+), gamma_H=0,
commuting field Psi, H=-Ical gamma^0, Z>0 and kappa>0. Write

\[
S=\bar\Psi\Psi,\qquad P=\bar\Psi\mathcal J\Psi,\qquad
A^{abc}=\bar\Psi\gamma^{abc}\Psi,\qquad \bar\Psi=\Psi^TH.
\]

Differentiating SM's universal quartic identity
A_abc A^{abc}=-6(S^2+P^2), using the symmetry of H gamma^{abc}, yields

\[
\boxed{A_{abc}\gamma^{abc}\Psi=-6(S+P\mathcal J)\Psi.}
\tag{CS.1}
\]

Indeed, variation of the two sides gives respectively
4 delta Psi^T H A_abc gamma^{abc} Psi and
-24 delta Psi^T H(S+P Jcal)Psi. Since H is invertible, equality for every
delta Psi proves (CS.1). `coupled_cosmology.py` also exhausts all 120 cubic
monomials in each of its eight components, independently of sampled states.

The eliminated action is consequently

\[
S_{\rm eff}=S_g[e,\omega_{\rm LC}]+S_{\rm EM}
 +\int w\left[\frac Z2\bar\Psi\gamma^\mu D_\mu^{\rm LC}\Psi-U(S,P)\right]d^4x,
\]
\[
U(S,P)=\frac{Zm}{2}S+\frac{3\kappa Z^2}{64}(S^2+P^2),\qquad m\geq0,
\tag{CS.2}
\]

and its matter equation reads

\[
\left[\gamma^\mu D_\mu^{\rm LC}-m
       -\frac{3\kappa Z}{16}(S+P\mathcal J)\right]\Psi=0.
\tag{CS.3}
\]

The positive sign of the quadratic term in U follows from SM's actual
commuting-field Fierz identity and contact-action sign. It is not an
independently adjustable cosmological coupling in this model.

Choose the rank-four real eigenspace

\[
P_+=\tfrac12(I+H),\qquad H\Psi=\Psi.
\tag{CS.4}
\]

Direct native matrix identities give

\[
S=\Psi^T\Psi>0,\quad P=0,\quad
V^a=(S,0,0,0),\quad A^{123}=0,\quad
\gamma^0\Psi=\mathcal I\Psi,\quad [P_+,\mathcal I]=0.
\tag{CS.5}
\]

For example, P_+ Acal_i P_+=P_+ H Jcal P_+=0. Thus native phase rotation
and scalar dilution preserve this rest sector. It selects a positive scalar
branch of the classical homogeneous field, not a positive-energy quantum
Hilbert space. The R5 Hamiltonian counterexample outside this sector remains
unchanged.

## CS-2. Homogeneous field, neutral Maxwell sector and stress

Use cosmic time and a spatially flat coframe

\[
e^0=dt,\qquad e^i=a(t)dx^i,\qquad
h=\frac{\dot a}{a},\quad a>0,\qquad
\omega_{\rm LC}^{0i}=h e^i,\quad \omega_{\rm LC}^{ij}=0.
\tag{CS.6}
\]

The symbol h denotes the expansion rate and is unrelated to Planck's constant.
For a spatially homogeneous spinor, SM's lowered-index spin representation
gives

\[
\gamma^\mu D_\mu^{\rm LC}\Psi
=\gamma^0\left(\dot\Psi+\tfrac32h\Psi\right).
\]

In the rest sector define the classical number density and effective phase rate

\[
n=\frac Z2S>0,\qquad M=m+\frac{3\kappa}{8}n.
\]

Equation (CS.3) becomes

\[
\boxed{\dot\Psi=-\tfrac32h\Psi-M\mathcal I\Psi,\qquad
       \dot n=-3hn,\qquad n=\frac{n_0}{a^3},\quad n_0>0.}
\tag{CS.7}
\]

Here n_0 is an integration constant; it is not a prediction of a particle
number or measured cosmological density.

**Maxwell consistency.** The explicit solution uses the neutral sector
g=0, a_mu=0, F=0, so all Maxwell equations hold. If instead one left a nonzero
charge coupling g in this one-field ansatz, the current would be j^0=gn and
j^i=0. A homogeneous electromagnetic field has
partial_i(w F^{i0})=0, whereas Gauss's law would demand wgn!=0. In particular,
F=0 with a nonzero charged field is not a solution. Adding compensating
charged species would require varying a different coupled matter system; it
has not been silently assumed here. The global native phase symmetry and
its conserved number remain meaningful at g=0.

**Stress from the action.** On shell, the effective symmetric metric source is

\[
T_{\mu\nu}^{\rm eff}
=-\frac Z4\left[\bar\Psi\gamma_{(\mu}D_{\nu)}^{\rm LC}\Psi
 -(D_{(\mu}^{\rm LC}\bar\Psi)\gamma_{\nu)}\Psi\right]
 +g_{\mu\nu}\mathcal L_{\rm eff}/w,
\tag{CS.8}
\]

with normalized symmetrization. This follows by varying the coframe and the
induced Levi–Civita spin connection in the eliminated action; the connection
variation supplies the symmetric spin contribution. Its sign also follows
from SM's convention delta S_m/int delta e^a wedge T_a. On (CS.5)–(CS.7),
the spatial off-diagonal spin bilinears cancel in the symmetric combination,
and V^i=0 removes the momentum density. The result is

\[
T_{00}=\rho,\qquad T_{0i}=0,\qquad T_{ij}=p\,g_{ij},
\]
\[
\boxed{\rho=mn+\frac{3\kappa}{16}n^2,\qquad
       p=\frac{3\kappa}{16}n^2.}
\tag{CS.9}
\]

An independent lapse check fixes the energy sign. Before imposing cosmic
time, set e^0=N(t)dt. Per unit coordinate spatial volume, the homogeneous
matter action is

\[
L_{\rm hom}=\frac Z2a^3\Psi^T\mathcal I\dot\Psi-Na^3U(S,0).
\tag{CS.10}
\]

The expansion part of the spin connection contributes the skew bilinear
Psi^T Ical Psi=0. Varying N at fixed field gives -a^3 rho. Varying a and only
then using its field equation gives 3Na^2 p, since the on-shell kinetic
scalar is ZMS/2 and ZMS/2-U=3 kappa n^2/16. Fixing N before varying would
discard this independent energy constraint. The verifier checks both
variations directly in the covariant density, and all sixteen components of
(CS.8), rather than inserting a dust or fluid stress by assumption.

The two contributions scale as a^{-3} and a^{-6}. The latter has p=rho for
its own contribution, conventionally called a stiff equation of state. This
is a coherent classical spinor solution, not a thermal fermion or radiation
gas.

Although the metric and effective stress are isotropic, the unaveraged spin
and torsion generally select a spatial axis. In particular A^{123}=0 while
A^{0ij} is nonzero. The full matter–connection configuration is homogeneous
but is **not** asserted to be rotationally invariant. No spin averaging is
used to obtain (CS.9).

## CS-3. The gravitational constraint and local existence

The remaining metric equations are

\[
\boxed{3h^2=\Lambda+\kappa\rho,\qquad
       2\dot h+3h^2=\Lambda-\kappa p.}
\tag{CS.11}
\]

Together with dot n=-3hn, these imply
dot rho+3h(rho+p)=0. Their Hamiltonian defect
C=3h^2-Lambda-kappa rho satisfies

\[
\dot C=-3hC
\tag{CS.12}
\]

when evolving n and h by the continuity and spatial Einstein equations,
without first setting C=0. Thus constrained initial data remain constrained;
the energy equation is not an optional normalization fitted afterwards.

In the expanding branch let v=a^3. Equations (CS.7), (CS.9) and (CS.11) give

\[
\boxed{\dot v^{\,2}=3\Lambda v^2+3\kappa m n_0v
                          +\frac9{16}\kappa^2n_0^2,\qquad
       \ddot v=3\Lambda v+\frac32\kappa m n_0.}
\tag{CS.13}
\]

The coefficient 9 kappa^2 n_0^2/16 is fixed by the torsion interaction. It
cannot be chosen freely when integrating the second-order volume equation.

For v>0, n>0 and fixed parameters, the equations for (v,n,h,Psi) have smooth
right-hand sides in a fixed coframe/spin gauge. Ordinary differential equation
uniqueness therefore applies to their homogeneous constrained initial-value
problem. The explicit family below proves existence and smoothness on t>0;
it does not supply a theorem for arbitrary spatially dependent Einstein–matter
data or for passage through v=0.

## CS-4. Closed solution in the positive-Lambda native sector

Take LR's external-return representation, sigma=-1, with supplied inverse
length u_l>0 and readout beta>0:

\[
\Lambda=12u_\ell^2,\qquad
\kappa=\frac1{2\beta u_\ell^2},\qquad
\kappa\Lambda=\frac6\beta.
\tag{CS.14}
\]

Thus no zero-scale or infinite-readout limit is needed for the main family.
Define

\[
\omega=\sqrt{3\Lambda}=6u_\ell,\qquad
d=\frac{\kappa m n_0}{2\Lambda},\qquad
b=\frac{3\kappa n_0}{4\omega}>0.
\]

Choosing the singular endpoint as t=0, the expanding solution is

\[
\boxed{v(t)=d[\cosh(\omega t)-1]+b\sinh(\omega t),\quad
       a(t)=v(t)^{1/3},\quad n(t)=n_0/v(t),\qquad t>0.}
\tag{CS.15}
\]

This has v>0 and dot v>0 because d>=0 and b>0. Its derivative at the endpoint
is omega b=3 kappa n_0/4; direct substitution proves both equations (CS.13).
An arbitrary translation of t sets a different time origin.

For completeness the native matter phase also has a closed form. Put
z=exp(omega t), choose any t_*>0 with z_*=exp(omega t_*), and define

\[
\boxed{\theta(t)=\theta_*+m(t-t_*)
 +\frac12\log\!\left[
 \frac{z-1}{z_*-1}\,
 \frac{(d+b)z_*-(d-b)}{(d+b)z-(d-b)}\right].}
\tag{CS.16}
\]

Every logarithm argument here is positive and dimensionless. To derive it,
factor v=[(z-1)((d+b)z-(d-b))]/(2z) and integrate
dot theta=m+3 kappa n_0/(8v), using omega b=3 kappa n_0/4. It follows that

\[
\boxed{\Psi(t)=v(t)^{-1/2}
       \exp[-\mathcal I\theta(t)]\Xi,\qquad
       H\Xi=\Xi,\qquad \frac Z2\Xi^T\Xi=n_0.}
\tag{CS.17}
\]

Since Ical^2=-I, the exponential is the real matrix
I cos(theta)-Ical sin(theta). Equations (CS.15)–(CS.17) specify real fields
without a numerical time integrator. Finally reconstruct the full connection:

\[
\boxed{\omega^{ab}=\omega_{\rm LC}^{ab}+K^{ab},\qquad
K^{ab}=\frac{\kappa Z}{8}\eta_{cc}A^{abc}e^c,\qquad
T_{abc}=-\frac{\kappa Z}{4}A_{abc}.}
\tag{CS.18}
\]

There is a sum on c in K^{ab}. These fields are analytic at every t>0.
The nonzero spinor and torsion decay as t tends to infinity; h tends to
sqrt(Lambda/3). Their past endpoint is treated in CS-6, not removed by
calling the full connection scalar finite.

## CS-5. Why the full first-order equations hold

SM establishes that the independent connection equation has the unique
solution (CS.18), and that substituting it into the action retains the
contact term (CS.2). At that stationary connection, varying the eliminated
action with respect to the other fields is equivalent to varying the full
action: the additional chain-rule term is the connection Euler residual
times delta omega, which vanishes. The same reasoning applies to the matter
equation. Hence (CS.3), (CS.8), (CS.11) and (CS.18) prove that the ansatz
solves the complete first-order system, including the off-diagonal equations.

For independent computational controls, `coupled_cosmology.py` constructs
coframe and spinor first jets, then constructs omega_LC+K and differentiates
it directly. Its curvature is d omega+omega wedge omega; its torsion is
de+omega wedge e. It does not replace that curvature with a Friedmann
expression. At three rational solution jets, including nontrivial spin
orientations, a massless case and unequal constant spatial coordinate
scalings, it checks:

- all 16 coframe variations of the original native gravitational trace
  action plus the original SM matter action, with curvature and connection
  held independent during each variation;
- all 24 connection equations epsilon_abcd T^a wedge e^b=-kappa S_cd;
- all 8 original first-order matter components, including the torsion-trace
  convention, and the vanishing torsion trace;
- the effective stress, lapse/scale variations and scalar-curvature readout
  by separate computations.

Taking z rational makes cosh(omega t) and sinh(omega t) rational. Constants
and normalization are chosen so the checked local field jets are rational;
the analytic solution is not claimed to be rational at every time. A free
phase constant can set the phase at each checked point. These fixtures
support the written general deduction; sampling them alone would not prove
the solution for arbitrary parameters.

A negative control deletes K from the original action while keeping the
claimed fields. Matter, coframe and connection residuals all become nonzero.
Thus the solution has not been obtained by silently discarding the very
interaction responsible for its a^{-6} term.

## CS-6. What the interaction changes, and what a scalar readout misses

For kappa>0, m>=0, n>0 and Lambda>=0, (CS.11) gives

\[
3h^2=\Lambda+\kappa mn+\frac{3\kappa^2}{16}n^2>0.
\tag{CS.19}
\]

Consequently this flat rest branch has no bounce at positive volume. Its
torsion contribution is positive stiff energy. Changing it to the negative
stiff term of another coupling or spin-fluid model would change the R5
action; no such change has been made. The effect of torsion on a bounce is
coupling-dependent in the established literature, not automatic.

As t approaches zero from above,

\[
v\sim\frac{3\kappa n_0}{4}t,\qquad
a\sim\left(\frac{3\kappa n_0}{4}t\right)^{1/3},\qquad
n\sim\frac4{3\kappa t}.
\]

For m>0 there can be an intermediate dust-dominated interval when both the
stiff and Lambda contributions are small compared with mn. Its existence
and duration depend on the chosen parameters; it is not guaranteed for every
member of the positive-Lambda family. The future is Lambda dominated.

The scalar curvature of the metric and that of the independent connection
are different readouts. With SM's torsion convention and the axial solution,
the divergence terms vanish and

\[
R(\omega)=R(g)-K_{abc}K^{abc},\qquad
K_{abc}K^{abc}=-\frac{3\kappa^2}{8}n^2.
\]

Taking the metric trace of (CS.11) gives

\[
\boxed{R(g)=4\Lambda+\kappa mn-\frac{3\kappa^2}{8}n^2,\qquad
       R(\omega)=4\Lambda+\kappa mn,\qquad
       T_{abc}T^{abc}=-\frac{3\kappa^2}{2}n^2.}
\tag{CS.20}
\]

The last contraction is Lorentzian, so its negative sign is not a positivity
violation. In the massless member, **R(omega)=4 Lambda is constant for the
entire solution while T_abc T^{abc} diverges at t=0**. R(g) diverges there
as well. The cancellation in one curvature scalar cannot establish a smooth
metric–connection extension or nonsingularity. This gives an explicit
coupled example of a limited geometric readout hiding a singular invariant.

The zero-Lambda contraction is also elementary:

\[
v(t)=\frac{3\kappa n_0}{4}t(1+mt),\qquad
\theta(t)-\theta_* =m(t-t_*)+
\frac12\log\!\left[\frac{t(1+mt_*)}{t_*(1+mt)}\right].
\tag{CS.21}
\]

For m>0 it interpolates from a proportional to t^{1/3} to a proportional to
t^{2/3}. At fixed kappa this is the limit u_l->0, beta->infinity of
(CS.14), or a Lambda=0 member of the more general MG action. It is not a
finite-parameter Lambda=0 point of LR's nonzero loop-scale action.

## CS-7. A complete benchmark with exact numbers

Choose dimensionless units and parameter values

\[
\kappa=1,\quad Z=\frac43,\quad m=\frac12,\quad
u_\ell=\frac16,\quad \beta=18,\quad \Lambda=\frac13,
\qquad \Xi=(1,0,0,0,0,0,-1,0)^T.
\]

Then H Xi=Xi, Xi^T Xi=2, n_0=4/3, d=b=1 and omega=1. The full solution is

\[
\boxed{v=e^t-1,\quad a=(e^t-1)^{1/3},\quad
\Psi=v^{-1/2}\left[I\cos\!\left(\tfrac12\log v\right)
 -\mathcal I\sin\!\left(\tfrac12\log v\right)\right]\Xi,\quad t>0.}
\tag{CS.22}
\]

It has n=4/(3v), h=e^t/(3v), theta-dot=e^t/(2v), and

\[
\rho=\frac2{3v}+\frac1{3v^2},\quad p=\frac1{3v^2},\qquad
A^{013}=-\frac2v,\quad K_{013}=\frac1{3v},\quad T_{013}=-\frac2{3v}.
\]

The other nonzero axial components are fixed by antisymmetry. The connection
is completely specified by (CS.18), and F=0 with g=0. For z=e^t:

| z | v | h | rho | p | R(g) | R(omega) | T_abc T^{abc} |
|---|---|---|---|---|---|---|---|
| 3/2 | 1/2 | 1 | 8/3 | 4/3 | 0 | 8/3 | -32/3 |
| 2 | 1 | 2/3 | 1 | 1/3 | 4/3 | 2 | -8/3 |
| 3 | 2 | 1/2 | 5/12 | 1/12 | 3/2 | 5/3 | -2/3 |
| 5 | 4 | 5/12 | 3/16 | 1/48 | 35/24 | 3/2 | -1/6 |

For instance, at t=log 2, 3h^2=4/3=Lambda+kappa rho and dot h=-2/3.
The code verifies the table with exact fractions. These are values of a
chosen solved model, not estimates of G, particle masses or cosmological
observables. This example's relation m=3u_l is a convenient parameter choice;
the general solution (CS.15) permits arbitrary m>=0 and does not derive that
relation.

## CS-8. Reproduction, attribution and the next open problem

Run the unchanged entry point with Python 3.11 or 3.12:

```bash
python papers/ugd-kahler-propagation/verify.py --check
```

The R6 certificate adds the coupled solution controls, pins the unchanged SM
proof and code to their R5 commit, and preserves all ten prior result groups.
The complete cubic polynomial check, 72 positive-Lambda volume/phase cases,
36 zero-Lambda contraction cases, 16 off-constraint propagation checks and
the full first-order solution jets are recorded separately. The verifier is
read-only by default; these exact computational controls are not a
proof-assistant certificate or independent review.

The former open existence question now has an explicit answer for this
homogeneous neutral rest sector: real analytic coupled fields exist for t>0,
with local uniqueness after constrained initial data and gauges are fixed.
The generic spatially varying Cauchy problem, stability against perturbations
outside this sector, charged multi-field solutions and a native selection of
the action and scales remain open. No empirical cosmology or quantum
statistics is inferred from the benchmark.

Primary attribution:

1. C. Armendáriz-Picón and P. B. Greene, *Spinors, Inflation, and Non-Singular
   Cyclic Cosmologies*, General Relativity and Gravitation **35**, 1637–1658
   (2003), [doi:10.1023/A:1025783118888](https://doi.org/10.1023/A:1025783118888),
   [arXiv:hep-th/0301129](https://arxiv.org/abs/hep-th/0301129). Homogeneous
   classical spinors as self-consistent cosmological sources are established;
   this note derives the particular potential fixed by SM's native action.
2. J. Magueijo, T. G. Zlosnik and T. W. B. Kibble, *Cosmology with a spin*,
   [arXiv:1212.0585](https://arxiv.org/abs/1212.0585). This studies cosmological
   spinor sources in Einstein–Cartan–Holst theory and explains the dependence
   of the effective interaction's energy sign on the coupling choices. Our
   no-bounce conclusion is specific to (CS.2) and its stated branch, not a
   claim about every Einstein–Cartan cosmology.
