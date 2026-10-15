# NAC-1 — From a calibrated native angle to finite transport and curvature

15 October 2026, India.

## 1. Contract and what is reused

Use the existing real native completion, central quarter-turn extension,
positive dagger pairing and finite operator representation. In the native
two-mode sector take `R^2=J^2=-I`, `RJ=-JR`, and `L=RJ`, as in QG-R1.
Here `J` is the anti-dagger direction `iota K`, not NJC-1's real reflection.
For real vector words put

\[
 \operatorname{sc}(Q)=\tfrac12\operatorname{Tr}Q,\qquad
 \langle A,B\rangle=-\operatorname{sc}(AB).
\]

Real unit combinations of `R,J,L` are anti-dagger, have square `-I`,
and have unit length for this pairing. Two such combinations satisfy

\[
 A^2=B^2=-I,\qquad AB+BA=-2\mu I,\qquad
 \mu=\langle A,B\rangle\in[-1,1].                         \tag{1}
\]

The angular integral and its lift come from
[phase observation, PM1/AC1/AC2](../native-phase-observation-certification/CLAIM_COVERAGE.md).
The circular exponential is NJC-1 T4. The four-edge curvature convention,
finite relative quotient and variable-field terms are
[QG-R1 R1–R7](../native-relative-curvature/THEOREM.md).
N4H-1 H3/H6 already establish loss of the geometric angle lift and the
integer-character criterion. Those results are inputs, not rediscoveries.

This packet organizes their angular calibration, proves the exact finite
loop coefficient for general direction overlap, and gives the precise
winding residue required by a fractional rescaling. It does not claim
priority for exponential composition, circle characters, quaternion
commutators or the classical Chebyshev identities.

## 2. Exponential transport and the native turn unit

For a fixed unit circular generator `A`, the native Cayley curve is

\[
 C_A(t)=(I-tA)^{-1}(I+tA)
       =\frac{1-t^2}{1+t^2}I+\frac{2t}{1+t^2}A,
 \quad
 \Phi(t)=2\int_0^t\frac{du}{1+u^2},\quad
 \pi=\lim_{t\to+\infty}\Phi(t).                            \tag{2}
\]

This is the inherited native half-turn convention. Differentiation of
the rational coefficients gives

\[
 \frac{dC_A}{dt}=\frac{2}{1+t^2}A C_A,
 \qquad \frac{dC_A}{d\Phi}=A C_A.
\]

Consequently `C_A(t)=exp(Phi(t) A)`. Indeed multiply the last equation
by `exp(-Phi A)` and differentiate; the result is constant and equals
identity at zero. Existence of the exponential and termwise
differentiation use the inherited bounded finite-generator completion.
Continuation with the retained angle lift gives

\[
 U(\theta)=e^{\theta A}=\cos\theta I+\sin\theta A,
 \quad U(\theta+\eta)=U(\theta)U(\eta),
 \quad U(\pi)=-I,\quad U(2\pi)=I.                         \tag{3}
\]

The cosine and sine here are the even and odd coefficients of the native
exponential. No decimal trigonometric evaluation is needed by the certificate.

**Calibration theorem.** Retain a nonzero reference angle `theta0` and
its generator increment `L0=theta0 A`, including the lift. Then for every
real `theta` in the completed sector,

\[
 \boxed{U(\theta)=\exp\!\left(\frac{\theta}{\theta_0}L_0\right).}
                                                               \tag{4}
\]

**Proof.** Substitute `L0=theta0 A` in the exponent. Under an angular-unit
change `a=k theta`, `a0=k theta0`, `A_a=A/k`, both the ratio and the
generator increment remain unchanged. This proves unit independence.
Equivalently `nu=theta/(2 pi)` is a turn fraction and
`U=exp(nu G_turn)`, `G_turn=2 pi A`. QED.

Equation (4) requires `L0`, not merely its exponential. For example,
`exp(0 A)=exp(2 pi A)=I`, but their half-angle transports are `I` and `-I`.
The scalar constants `e=exp(1)` and `pi` are fixed in this completion.
They are not position-dependent fields. The local data can be `A(x)`,
the angular increment and the ordered path. Their relationship is carried
by the exponential argument and angle ratios, not a bare number `e/pi`.

## 3. Exact angular rescaling of one loop

In a fixed oriented circular plane let `H=c I+s A`, `c^2+s^2=1`, and
`kappa=1-c`. This is the compact defect; it does not retain the entire
geometric deficit. Multiplication gives

\[
 H^2-2cH+I=0,
 \qquad H^{n+1}=2cH^n-H^{n-1}\quad(n\geq1).                \tag{5}
\]

Define scalar polynomials `P0=1`, `P1=c`,
`P_(n+1)=2 c P_n-P_(n-1)`. Taking scalar parts in (5), with the two
initial values, proves for all integers `n>=0`

\[
 \boxed{\operatorname{sc}(H^n)=P_n(c),\qquad
 \kappa_n=1-P_n(1-\kappa).}                                \tag{6}
\]

In particular

\[
 \boxed{\kappa_2=4\kappa-2\kappa^2.}                       \tag{7}
\]

The proof is induction from the certified quadratic relation. Negative
integer powers have the same compact defect as their positives. These
are the usual first-kind Chebyshev polynomials; NIST DLMF
[18.5.1](https://dlmf.nist.gov/18.5.E1) supplies the classical identification
`P_n(cos theta)=cos(n theta)`. The recurrence is derived natively here,
not imported as a computational engine.

For `H=(3 I+4 A)/5`, `kappa=2/5`, so `kappa2=32/25` and
`kappa3=242/125`. The doubled phase is the geometric phase of MGH-1;
this checks angular consistency without selecting its matter interaction.
The law concerns repeating the same holonomy, or rescaling a fixed
generator. It does not assert that an arbitrary larger geometric cell
has holonomy `H^n`.

## 4. Exact finite-angle loop law

For fixed directions satisfying (1), set

\[
 U=\cos\theta I+\sin\theta A,\quad
 V=\cos\phi I+\sin\phi B,\quad
 H=V^\dagger U^\dagger VU.
\]

The direction commutator obeys

\[
 [A,B]^\dagger[A,B]=4(1-\mu^2)I.                            \tag{8}
\]

Define the normalized squared operator length
`||Q||_sc^2=sc(Q^dagger Q)` and the existing return defect

\[
 \chi(H)=\tfrac12(1-\operatorname{sc}H)
        =\tfrac14\|I-H\|_{\rm sc}^2.
\]

**Finite-loop theorem.** For all real angles in this fixed-direction
protocol,

\[
 \boxed{\chi(H)=(1-\mu^2)\sin^2\theta\sin^2\phi.}           \tag{9}
\]

**Proof.** Write `a=cos theta`, `b=sin theta`, `c=cos phi`,
`d=sin phi`. Reduce the products using only (1), `a^2+b^2=1` and
`c^2+d^2=1`. Direct noncommutative multiplication yields

\[
 H+H^\dagger=2I-4b^2d^2(1-\mu^2)I.
\]

The scalar part of a real native unit word equals that of its dagger,
so this proves (9). Unitarity also gives
`(I-H)^dagger(I-H)=2I-H-H^dagger`, proving the defect expression.
For (8), use `[A,B]=2(AB+mu I)` and
`(AB)^2=-I-2mu AB`. These identities have native rewrite witnesses.
QED.

If a reference loop has `sin(theta0) sin(phi0) != 0`, its observed defect
therefore determines

\[
 \boxed{\mathcal K=
 \frac{\chi(\theta_0,\phi_0)}{\sin^2\theta_0\sin^2\phi_0}
 =1-\mu^2=\frac14\|[A,B]\|_{\rm sc}^2,}
                                                               \tag{10}
\]
\[
 \boxed{\chi(\theta,\phi)=\mathcal K\sin^2\theta\sin^2\phi.}
                                                               \tag{11}
\]

This is the precise sense in which an angular calibration predicts any
other angle. It determines a scalar return defect, not the full oriented
holonomy matrix. Opposite signs of `mu` can give the same scalar defect.
Prediction (11) remains defined at zero sine, although such a reference
loop cannot be used for the division in (10).

For nonzero `theta phi` and `|mu|<1`, the Gram area of the retained
generator increments is
`A_ang^2=(theta phi)^2(1-mu^2)`. Thus

\[
 \boxed{\frac{\chi(H)}{\mathcal A_{\rm ang}^2}
 =\left(\frac{\sin\theta}{\theta}\right)^2
  \left(\frac{\sin\phi}{\phi}\right)^2.}                    \tag{12}
\]

It approaches one as both angles vanish. This is an exact finite-angle
correction to the fixed-direction small-loop limit. It uses generator
increments, so changing an angular unit together with its generator
leaves the area unchanged. It is a different finite readout from QG-R1's
Cayley-area quotient; their shared infinitesimal normalization does not
make them equal at finite angles.

An exact example uses Cayley parameters `t=1/2`, `u=1/3`, so
`sin theta=4/5`, `sin phi=3/5`. With perpendicular directions,
`chi=144/625`, reproducing QG-R1 R6/R7. With
`B=(3R+4J)/5`, it becomes `2304/15625`. With `B=R`, it is zero.
All three have the same angular calibration and the same constants.
Consequently angle or `e/pi` alone does not specify the loop response.

## 5. Exactly what a fractional angle needs

Fix an oriented axis and an endpoint phase with a canonical angle `phi`.
Its possible lifts are `theta=phi+2 pi n`, `n` an integer. For a specified
rescaling `r=p/q` in lowest terms, `q>0`, define

\[
 U_r(n)=\exp\!\left(r(\phi+2\pi n)A\right).
\]

**Winding-residue theorem.** Under the faithful circular convention of
(3),

\[
 \boxed{U_r(n_1)=U_r(n_2)
 \quad\Longleftrightarrow\quad n_1\equiv n_2\pmod q.}       \tag{13}
\]

**Proof.** Cancel the common phase. Equality means
`r(n1-n2)` is an integer. Coprimality of `p,q` makes this equivalent to
`q` dividing `n1-n2`. QED.

Thus an integer rescaling uses the endpoint alone; a half-rescaling
needs winding parity; a `3/4` rescaling needs winding modulo four.
For an irrational `r`, equality forces `n1=n2`, since a nonzero integer
multiple of an irrational number is not an integer. The full lift is
therefore sufficient for every real rescaling and is necessary for
exact unrestricted irrational rescaling. This is a mathematical exactness
statement; no robustness claim for noisy observations is made.

The result refines N4H-1's integer-character condition and CRC-1's
operation-dependent record rule. An axis and comparator contract are
still needed when the full operator, rather than just its scalar defect,
is requested. A cut-conjugation channel can itself lose an additional
sign, as already recorded in PM1/PM2.

## 6. Scope, controls and use in the gravity programme

The calibration has material rejection controls:

1. Equal angles with parallel and perpendicular directions give different
   defects. The relative-direction record cannot be omitted.
2. Noncommuting generators can have `H=I` at a finite half-turn because
   `U(pi)=-I` commutes with every `V`. Zero finite return alone does not
   establish a flat connection on a neighborhood.
3. Equal endpoint rotations with different winding parity give different
   half-angle transports. Endpoint matrix fractional powers are ambiguous.
4. `kappa2` differs from `4 kappa` at finite angle. The missing term in
   a quadratic extrapolation is exactly `-2 kappa^2`.
5. If directions or connection coefficients vary along the controls,
   the general curvature includes the already certified derivative terms
   `F_xy=partial_x A_y-partial_y A_x+[A_x,A_y]`. QG-R1 R3/R4 gives
   identical instantaneous directions with different curvature. The
   fixed-direction calibration is not applied to that unrestricted case.

For the geometric hinge already constructed in N4H-1, retain its oriented
normal generator and lifted deficit. Equations (4)–(7) then describe
rescaled normal transport and its compact defect. For the native two-cut
loop, (10) isolates the direction-commutator strength from finite angular
effects. This improves the transport/curvature dictionary without adding
a new matter coupling or claiming an Einstein field equation.

The executable checks exact native matrix fixtures and serialized universal
rewrite witnesses. All-angle completion, induction and lift arguments are
written proofs under the contracts above; they are not inferred from a finite
sample grid and are not advertised as proof-assistant formalizations.
