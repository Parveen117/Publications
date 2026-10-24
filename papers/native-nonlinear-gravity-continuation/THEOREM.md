# EGR-2 — Nonlinear cut-source response, local continuation and focusing

14 October 2026. All lengths, areas, volumes and coefficients below are
dimensionless. The operation parameter is not a physical clock.

**Result.** For the explicit source functional below, the nine-cell model
has an analytic local solution of all ten nonlinear metric equations for
each fixed source direction with both channels nonzero. The source
`r=(3/5,4/5)` gives an expanding saddle branch. Its complementary cut,
`r=(-3/5,-4/5)`, gives a strict local minimum with all nine volumes initially
decreasing. Both branches have the same source amount and cut spectrum.
The first two response coefficients, the continuation Jacobian and the
initial volume derivatives are obtained exactly in the native arithmetic.

The stationary branch also supplies a nonlinear cut-feedback equation.
Its cut transport remains exponential; its scalar channel-ratio evolution
has a nonzero correction to EGR-1's simple exponential law.

This is a finite Euclidean source-coupled geometry result. The analytic
theorem is proved below using the analytic implicit-function theorem;
the executable checks its finite algebraic hypotheses. There is no
claimed finite convergence radius, global collapse trajectory,
Lorentzian quantum measure or full quantum-gravity certification.

## 1. Inherited geometry and the additional source contract

Use exactly the nine cells, boundary lengths, flat point q0, interior edge
order and lifted deficits of [EGR-1](../native-exponential-gravity-response/THEOREM.md).
The edge order is

`AO, BO, BP, CO, CP, DO, DP, EO, EP, OP`.

P in those edge labels is the geometric vertex P_geo. Below, P in an
operator equation is the observation cut. Boundary vertices are A=0 and
B,C,D,E the four coordinate unit vectors. At the flat point,
O=(1/5,1/5,1/5,1/5) and P_geo=(6/25,6/25,6/25,6/25).

Let R(q) be the Regge curvature functional including the boundary term
and sign convention of N4H-1/EGR-1. For fixed boundary lengths its
interior derivative is

\[
f_i(q)=\partial_i\mathcal R(q)
=\sum_{h\text{ interior}}\partial_i a_h(q)\,\delta_h(q).
\tag{1}
\]

Near q0 all cell Grams are positive definite. Use positive area and
volume roots, interior angles in (0,pi), and the continuous deficit lift
delta_h=2pi-sum(theta_sigma,h). These choices make f real analytic there.
The flat subdivision gives f(q0)=0.

The pinned predecessor gives matrices

\[
\mathsf H=f'(q_0)=A^TKA,\quad K=\operatorname{diag}(10,47),\quad
AS=I_2,\quad AG=0,\quad \operatorname{rank}G=8.
\tag{2}
\]

G is the 10-by-8 differential of the squared lengths under displacing
O and P_geo. The columns of B=[S G] form a basis of the ten edge variables.
S selects edges AO and OP. The readout is

\[
\begin{aligned}
y_1&=dq_{AO}+dq_{BO}+dq_{CO}+dq_{DO}+dq_{EO},\\
y_2&=dq_{OP}+dq_{BP}+dq_{CP}+dq_{DP}+dq_{EP}
-\tfrac15(dq_{BO}+dq_{CO}+dq_{DO}+dq_{EO}),\qquad y=A\,dq.
\end{aligned}\tag{3}
\]

Take the native real rank-one cut P=(I+r1 Z+r2 X)/2, with r1^2+r2^2=1,
and calibrated amount N>=0. Z=diag(1,-1) and X has off-diagonal entries 1.
No primitive Hilbert-space premise is added. The finite positive carrier,
dagger, trace and scalar completion are those already declared in EGR-1.

**Selected nonlinear extension.** For gamma>0, alpha>=0 and a fixed cut,
set epsilon=alpha*N/gamma and choose

\[
\boxed{\mathcal S_\epsilon(q)/\gamma
=\mathcal R(q)-\epsilon\,r^TA(q-q_0).}
\qquad
\boxed{f(q)=\epsilon A^Tr.}
\tag{4}
\]

This extends the previous source term affinely in squared lengths while
retaining the full nonlinear geometric functional R. It is a specified
model choice. A is invariant under the old displacement directions only
to first order at q0. Equation (4) is not an exact nonlinear gauge-invariant
matter action, nor an identification of an arbitrary information measure
with stress-energy. A nonlinear completion with different second
derivatives can change the conclusions below.

## 2. The complete third derivative from native geometric identities

Write T[i,j,k]=partial_i partial_j partial_k R(q0). Thus

\[
f(q_0+t d)=t\mathsf H d+\tfrac12t^2T[d,d]+O(t^3).
\tag{5}
\]

For a cell Gram g along a straight squared-length variation, g''=0. Its
inverse J satisfies

\[
J'=-Jg'J,\qquad J''=2Jg'Jg'J,\qquad
\frac{v'}v=\tfrac12\operatorname{tr}(Jg').
\tag{6}
\]

Obtain the 5-by-5 barycentric Gram B from J as in EGR-1. For a hinge
opposite vertices i,j define

\[
\ell=\frac{B'_{ii}}{B_{ii}}+\frac{B'_{jj}}{B_{jj}},\quad
k=B'_{ij}-\tfrac12B_{ij}\ell,
\]
\[
\ell'=\frac{B''_{ii}}{B_{ii}}-\left(\frac{B'_{ii}}{B_{ii}}\right)^2
+\frac{B''_{jj}}{B_{jj}}-\left(\frac{B'_{jj}}{B_{jj}}\right)^2,
\quad k'=B''_{ij}-\tfrac12B'_{ij}\ell-\tfrac12B_{ij}\ell'.
\tag{7}
\]

The identity a_h theta'_ij=12 v k gives

\[
C_h=a_h\delta'_h=-\sum_{\sigma\supset h}12v_\sigma k_{\sigma h},
\quad
C'_h=-\sum_{\sigma\supset h}12(v'_\sigma k_{\sigma h}+v_\sigma k'_{\sigma h}).
\tag{8}
\]

C'_h differentiates the entire product a_h delta'_h. For s_h=a_h^2,
s_i=partial_i s_h, s_d=partial_d s_h and s_id=partial_d partial_i s_h,
differentiating (1) twice at the flat point gives

\[
\boxed{T[d,d]_i=\sum_h
\left(\frac{s_{id}}{s_h}-\frac{3s_i s_d}{4s_h^2}\right)C_h
+\frac{s_i}{2s_h}C'_h.}
\tag{9}
\]

For completeness: at delta=0, f''_i=2 a_id delta'_h+a_i delta''_h.
Insert a_i=s_i/(2a), a_id=s_id/(2a)-s_i s_d/(4a^3), and
delta''=(C'-a' delta')/a to obtain (9). No numerical differentiation or
angle approximation is involved.

The executable evaluates (9) on all ten coordinate vectors and their
45 pairwise sums, polarizes to recover the entire T, and checks all its
index symmetries. For each direction it checks both inverse-derivative
identities, Schlaefli cancellation and its derivative, and the inherited
first derivative H d. It also compares the final source contraction
against a fresh direct geometric evaluation.

## 3. Why fixing eight old null directions fails

Let z0=K^{-1}r and d0=S z0. Every first-order solution of (4) is
d=d0+G eta. A possible second-order edge jet is

\[
q=q_0+\epsilon d+\epsilon^2 e+O(\epsilon^3).
\]

The second-order equation is

\[
\mathsf H e+\tfrac12T[d,d]=0,
\qquad G^TT[d,d]=0.
\tag{10}
\]

The second condition is necessary because G^T H=0. For the source
r=(3/5,4/5), keeping eta=0 violates it. For example, the first four
entries of b=G^T T[d0,d0] are -245989/1767200; the last four are
-1318/2209. The certificate retains the full residual. Solving only the
two selected edge equations cannot remove this eight-component obstruction.

The needed correction is already part of q'(0). Adding a second-order
term along G cannot repair a wrong first derivative, since H G=0.

## 4. Exact solution of all eight compatibility conditions

**Flat-family identity.** For every eta,

\[
G^TT[G\eta,G\eta]=0.
\tag{11}
\]

Proof: small coordinate displacements of the two interior vertices
preserve a flat subdivision with the same external boundary. Let
q_flat(t)=q0+t G eta+t^2 b_flat(eta). Its interior gradient vanishes.
Twice differentiating f(q_flat(t))=0 gives
T[G eta,G eta]+2 H b_flat(eta)=0. Multiply by G^T. Polarization gives
the mixed version. The executable verifies every pure-null cubic
contraction of the recovered T, as well as this written geometric reason.

Consequently (10) is an affine linear equation in eta:

\[
b+D\eta=0,\quad b=G^TT[d_0,d_0],\quad
D_{ab}=2g_a^TT[d_0,g_b].
\tag{12}
\]

Set M=D/2. In the order O's four displacements followed by P_geo's four,

\[
\boxed{M(r)=
\begin{pmatrix}
-10r_1-\tfrac25r_2&2r_2\\
2r_2&-10r_2
\end{pmatrix}\otimes I_4,
\qquad \det M_{2\times2}=100r_1r_2.}
\tag{13}
\]

The native calculation checks both source coefficients of this linear
matrix identity. There is also a direct interpretation. For one ambient
component of the flat displacement, call the vertex motions o and p.
The quadratic flat-path changes of the two readouts (3) are

\[
5o^2,\qquad \tfrac15o^2+5p^2-2op.
\]

Contracting the twice-differentiated flat equation with d0 yields

\[
(o,p)M(r)(o,p)^T
=-2\left[5r_1o^2+r_2(\tfrac15o^2+5p^2-2op)\right],
\tag{14}
\]

which proves (13). It also displays the nonlinear displacement dependence
of the selected source functional.

For r1*r2!=0, eta=-D^{-1}b is unique. At r=(3/5,4/5),

\[
\eta_O=-\frac{456869}{21206400}\quad\text{in each of four components},
\qquad
\eta_{P_{geo}}=-\frac{4410869}{106032000}\quad\text{in each component}.
\tag{15}
\]

These are coefficients in a decomposition of a curved edge perturbation
into the flat tangent G and S. They are not an embedding of the curved
geometry into the original four-dimensional Euclidean coordinate frame.

## 5. A genuine local analytic solution, beyond a truncated jet

**Local continuation theorem.** Fix r with r1*r2!=0, the geometry, source
contract (4) and local angle branches. There is a real analytic solution
q(epsilon) of all ten equations (4), with q(0)=q0 and
q'(0)=d0+G eta selected by (12). It is locally unique in the scaled
coordinate neighborhood used below, and every cell remains nondegenerate
for sufficiently small |epsilon|.

Proof: introduce ten unknowns z in R^2 and eta in R^8 by

\[
q=q_0+\epsilon(Sz+G\eta).
\]

For nonzero epsilon define

\[
F_1=\frac{S^T[f(q)-\epsilon A^Tr]}\epsilon,
\qquad
F_2=\frac{G^T[f(q)-\epsilon A^Tr]}{\epsilon^2}.
\tag{16}
\]

Analytic Taylor expansion and (2) show that both have analytic extensions
to epsilon=0, with

\[
F_1(0,z,\eta)=Kz-r,
\quad
F_2(0,z,\eta)=\tfrac12G^TT[Sz+G\eta,Sz+G\eta].
\tag{17}
\]

At z0=K^{-1}r and the eta selected by (12), both vanish. The Jacobian
with respect to (z,eta) has block form

\[
\begin{pmatrix}K&0\\ *&M(r)\end{pmatrix}.
\tag{18}
\]

K is invertible and det M=(100r1r2)^4 is nonzero. The analytic
implicit-function theorem therefore gives analytic z(epsilon),
eta(epsilon). Since [S G] is invertible, vanishing of both projected
residuals at nonzero epsilon implies the full ten-component residual
vanishes. Continuity gives the solution at zero as well. Positive Gram
matrices form an open set, proving local admissibility. This is a proof
of an actual local solution, not an assertion that a finite Taylor
polynomial solves the exact equations. QED.

For each of the two displayed complementary sources, the determinant
in (18) is exactly 470*48^4=2494955520. Local uniqueness means uniqueness
near the specified (epsilon,z,eta) point; it excludes no remote branch
and gives no global continuation range. Source directions on a channel
axis have singular M and require a separate, higher-order analysis.

## 6. Exact response coefficients

For r+=(3/5,4/5), symmetry under permuting B,C,D,E gives four groups of
equal edge derivatives:

| Squared edge(s) | Base value | d=q'(0), positive source orientation |
|---|---:|---:|
| AO | 4/25 | 338371/13254000 |
| OB, OC, OD, OE | 19/25 | 456869/53016000 |
| PB, PC, PD, PE | 469/625 | 4410869/1325400000 |
| OP | 4/625 | 878369/82837500 |

Let

\[
z_2=-\tfrac12K^{-1}S^TT[d,d].
\]

The readout of the actual analytic solution has expansion

\[
\boxed{y_+(\epsilon)=
\epsilon\binom{3/50}{4/235}
+\epsilon^2
\binom{521098203253/22485570048000}{336288173/11711234400}
+O(\epsilon^3).}
\tag{19}
\]

The actual full edge coefficient at epsilon^2 is S z2+G eta1. The
coefficient eta1 requires the next order and is not computed here.
Nevertheless A G=0 makes (19) the exact second-order readout coefficient.
The certificate checks H S z2+T[d,d]/2=0 in all ten components, so S z2
is one valid edge jet through that order.

For r-=-r+, the action at positive epsilon is the plus action at negative
epsilon. Local uniqueness gives

\[
q_-(\epsilon)=q_+(-\epsilon),\qquad
y_-(\epsilon)=-\epsilon z_0+\epsilon^2z_2+O(\epsilon^3).
\tag{20}
\]

The first curvature variation remains EGR-1's L z0, now with a nonlinear
continuation. Because L has rank two and z0!=0, these small nonzero
solutions are curved. A(q-q0) is still a selected length readout; the
identity delta-linearization=L y is not asserted as an all-order formula
for every deficit.

## 7. Stability selects a contracting complementary branch

The affine source in (4) has zero metric Hessian. At the branch,

\[
G^T\mathcal R''(q(\epsilon))G
=\epsilon M(r)+O(\epsilon^2).
\tag{21}
\]

Indeed the leading coefficient is G^T T[d,G]=G^T T[d0,G]=M; the
pure-null cubic term vanishes by (11). In coordinates [S G] the top block
is K+O(epsilon), and the cross block is O(epsilon). Its Schur complement
is epsilon M+O(epsilon^2). Sylvester inertia and continuity therefore give
the following for sufficiently small positive epsilon:

| Source signs | Sign of M | Full metric Hessian inertia (positive, negative, zero) |
|---|---|---|
| r1>0, r2>0 | Negative definite | (2,8,0): saddle |
| r1<0, r2<0 | Positive definite | (10,0,0): strict local minimum |
| r1*r2<0 | Four positive, four negative | (6,4,0): saddle |

For r+, each ambient component has the block

\[
M_+=\begin{pmatrix}-158/25&8/5\\8/5&-8\end{pmatrix},
\qquad \det M_+=48.
\tag{22}
\]

The complementary cut has M-=-M+. Explicitly,

\[
P_+=\begin{pmatrix}4/5&2/5\\2/5&1/5\end{pmatrix},\qquad
P_-=I-P_+=\begin{pmatrix}1/5&-2/5\\-2/5&4/5\end{pmatrix}.
\tag{23}
\]

Both are positive rank-one cuts, with trace P=trace(P^2)=1. Both give
rho=N P with the same source amount and spectral invariants. The native
quarter-turn R0=[[0,-1],[1,0]] transports P+ to P- by conjugation.
Here the readout operators Z,X are held fixed. This is a source-orientation
change, not a simultaneous observer-frame transformation of P,Z,X.

The eight former null directions have been lifted by this source-coupled
nonlinear problem. They must not be silently divided out as exact gauge
directions at nonzero epsilon. The Hessian inertia is a property of the
specified Euclidean variational model. It does not count physical
graviton polarizations or prove stability of a quantum or time evolution.

## 8. Cell volumes, focusing and a precise degeneracy condition

The local branch is permutation symmetric: the action and source are
symmetric, and the locally unique solution must preserve that symmetry.
Use the exact finite squared lengths

\[
x=q_{AO},\quad u=q_{OB}=q_{OC}=q_{OD}=q_{OE},\quad
v=q_{PB}=q_{PC}=q_{PD}=q_{PE},\quad w=q_{OP}.
\]

The nine cell volumes fall into three classes. Their Gram determinants are

\[
\begin{aligned}
\Delta_{out}&=x-3\left(\frac{x+1-u}{2}\right)^2,\\
\Delta_{cap}&=4v-3,\\
\Delta_{in}&=(3u-2)w-3\left(\frac{w+u-v}{2}\right)^2,
\qquad v_\sigma^2=\Delta_\sigma/576.
\end{aligned}\tag{24}
\]

There are four outer cells O A with three of B,C,D,E, one cap P_geo BCDE,
and four inner cells O P_geo with three of B,C,D,E. The bases used to
derive these determinants are respectively A, P_geo and O. Expanding
their Gram determinants gives polynomial native word certificates.

**Exact symmetric admissibility.** All nine cell Grams are positive
definite if and only if

\[
\boxed{\Delta_{out}>0,\quad v>3/4,\quad u>2/3,\quad\Delta_{in}>0.}
\tag{25}
\]

Proof: the outer Gram has an identity 3-by-3 block, whose Schur complement
is Delta_out. The cap Gram is I4+(v-1)11^T, with eigenvalues 1,1,1,4v-3.
The inner Gram has a three-vertex block I3+(u-1)11^T with eigenvalues
1,1,3u-2; its remaining Schur complement is Delta_in/(3u-2).
These conditions are necessary and sufficient. Positive edge squares
alone are insufficient; the certificate includes the degenerate cap
at v=3/4 with all edge squares positive.

At q0 the determinants are 1/25,1/625,1/625. The exact initial logarithmic
volume derivatives are:

| Cell class | Number | Base volume | v_sigma'(0)/v_sigma(0) for r+ | For stable r- |
|---|---:|---:|---:|---:|
| Outer | 4 | 1/120 | 163103/848256 | -163103/848256 |
| Cap | 1 | 1/600 | 4410869/1060320 | -4410869/1060320 |
| Inner | 4 | 1/600 | 912209/1060320 | -912209/1060320 |

The native calculation obtains these rates by the Gram trace (6) and
independently by differentiating (24). They are strictly negative on the
stable complementary branch. Analyticity makes every volume decrease
on a sufficiently small positive source interval, while all volumes
remain positive. This is initial geometric focusing in this model.
The plus branch instead has nonzero curvature with all volumes increasing.

A density additionally requires an assigned receiving amount. If fixed
weights omega_sigma>0 with sum omega_sigma=1 allocate n_sigma=N omega_sigma,
define I_sigma=n_sigma/v_sigma. This allocation is an extra measurement
contract; it is not supplied by (4). At fixed positive n_sigma,

\[
\frac{d\log I_\sigma}{d\epsilon}
=-\frac{d\log v_\sigma}{d\epsilon},\qquad
\boxed{I_\sigma^2\Delta_\sigma=576n_\sigma^2.}
\tag{26}
\]

Thus the stable branch initially increases the assigned density when
epsilon varies through the coupling with N fixed. If instead alpha/gamma
is fixed and N varies, the numerator also varies: n_sigma is proportional
to epsilon and d log I/d epsilon=1/epsilon-d log v/d epsilon.
Neither operation is a spatial conservation or transport equation.

At fixed nonzero n_sigma a reachable limit Delta_sigma->0+ would force
I_sigma to diverge as 24 n_sigma/sqrt(Delta_sigma). The local theorem
does not show such a limit is reached. Extrapolating these initial slopes
to a collapse threshold is unjustified. Density blowup also must not be
identified with every curvature-singularity diagnostic without a separate
comparison theorem; QG-R1/R2's return-defect/area criteria remain distinct.

## 9. Nonlinear geometry backreacts on the exponential cut flow

For r near either displayed direction, the implicit-function argument
also applies with r as a parameter. Keep alpha,N,gamma fixed and small
enough for this local branch. Let q(P) denote its stationary metric and
y(P)=A(q(P)-q0). Define its on-shell functional

\[
W(P)=\gamma\left[\mathcal R(q(P))
-\epsilon r(P)^TA(q(P)-q_0)\right].
\tag{27}
\]

The metric derivative vanishes by stationarity, so the envelope derivative
along a cut variation is

\[
dW=\operatorname{Tr}(H_{eff}\,dP),\qquad
\boxed{H_{eff}(P)=-\alpha N\,[y_1(P)Z+y_2(P)X].}
\tag{28}
\]

Retain EGR-1's explicitly chosen cut gradient relaxation:

\[
\boxed{P'=[[P,H_{eff}],P],\qquad
U'=[P,H_{eff}]U,\qquad P=UP(0)U^\dagger.}
\tag{29}
\]

This preserves the cut constraint, source amount and spectrum, and
W'=-||[P,H_eff]||_F^2<=0 by the existing native cut identities. It is a
local gradient relaxation, not a derived physical clock law. On the
negative-quadrant branch the metric elimination uses a local minimum;
the positive-quadrant branch remains stationary but is a metric saddle.

Let R0=[[0,-1],[1,0]]. In the real sector,

\[
[P,H_{eff}]=\omega R_0,\qquad
\omega=\alpha N(r_1 y_2-r_2 y_1),\quad
r'_1=-2\omega r_2,\quad r'_2=2\omega r_1.
\tag{30}
\]

For r1!=0, the exact local ratio equation is

\[
\boxed{m'=\frac{2\alpha N}{r_1}(y_2-m y_1),\qquad m=r_2/r_1.}
\tag{31}
\]

Substitute y=epsilon K^{-1}r+epsilon^2 z2(r)+O(epsilon^3):

\[
m'=-\frac{37}{235}\frac{\alpha^2N^2}{\gamma}m
+\frac{2\alpha N\epsilon^2}{r_1}
\left(z_{2,2}-m z_{2,1}\right)+O(\alpha N\epsilon^3).
\tag{32}
\]

For both displayed cuts the second term has a nonzero coefficient,
checked exactly. Therefore EGR-1's simple scalar exponential ratio law
is the leading metric approximation, not an all-order law for (4).

The exponential transport itself remains exact. With
r=(cos(2phi),sin(2phi)) on a continuous local angular branch, phi'=omega,
and all generators commute as multiples of R0. Hence

\[
\boxed{U(\lambda)=\exp[(\phi(\lambda)-\phi(0))R_0].}
\tag{33}
\]

The nonlinear solution changes phi(lambda); it does not remove the
continuous unitary transport of a discrete-spectrum cut. Equations
(4), (28) and (29) are the local coupled geometry/cut system supplied by
this development. They do not include a derived spatial source-transport
law. At r1*r2=0 the continuation Jacobian degenerates; no global metric
branch through those orientations is inferred from this theorem.

## 10. Scope of the result and certificate

| Claim | Status |
|---|---|
| Full nonlinear metric equations have a local analytic source branch | Proved under (4), away from the channel axes |
| All ten second-order metric equations | Exact native check |
| All old flat-null directions can stay fixed | Rejected by a nonzero second-order residual |
| Stable complementary cut with initial volume focusing | Proved for the specified Euclidean action and sufficiently small source |
| Equal amount and spectrum force equal geometry or stability | Rejected by the two complementary cuts |
| Nonlinear geometric response feeds the cut force | Envelope identity (28) and local flow (29) |
| Simple ratio exponential remains exact to all metric orders | Rejected by a nonzero second-order coefficient |
| Exponential cut transport survives the nonlinear correction | Established locally by (30)-(33) |
| Finite-source singular collapse or transported singularity | Not established |
| Exact nonlinear diffeomorphism symmetry or a quantum measure | Not established |

The executable records the full T tensor, the rank/basis checks, general
null-lift coefficients, all second-order compatibility equations, both
source responses, inertia inputs, volume conditions and local cut-feedback
coefficients. Symbolic determinant and polynomial identities use the
unchanged Paaninian compiler, with serialized witnesses and replay.
Python 3.12 handles only extraction, hashing and process launching.

The proof of analyticity, the geometric flat family, analytic
implicit-function theorem, inertia continuity, local minima and
arbitrary-parameter cut-flow preservation are written mathematical
arguments. The exact arithmetic checks do not make this a generic
formal-proof-assistant certificate. Neither the new packet nor its
replay modifies EGR-1, earlier certificates or either native engine.

The dependency on nonlinear source choice is consequential: replacing
the affine source F=r^TA(q-q0) by a nonlinear F with the same first
derivative adds -epsilon F'' to the metric Hessian and changes the
null lift to M-G^T F''(q0)G. One must specify that completion before
promoting this branch's focusing or stability to a universal law.

For context on the difference between flat discrete displacement
symmetries and curved Regge equations, see B. Bahr and B. Dittrich,
[(Broken) Gauge Symmetries and Constraints in Regge Calculus](https://arxiv.org/abs/0905.1670),
Class. Quantum Grav. 26 (2009) 225011. That work is background context;
none of its coefficients are imported as a proof of this source model.
This packet claims a new step in this development sequence, not a
literature-wide novelty result.
