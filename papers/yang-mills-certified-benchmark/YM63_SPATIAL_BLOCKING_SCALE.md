# YM-63: spatial blocking retains memory and loop-shape response

Monty Dabas. 4 October 2026. **Python 3.12 only.**

We keep YM62's actual native link Hamiltonian
A(theta)=L_E-theta V. This continuation supplies an exact spatial cut,
its interacting return term, leading larger-loop responses and the scale
conditions needed before a continuum claim can be made. The principal
new coefficient is

\[
 \boxed{\langle W_{2\times2}\rangle_\theta
       =\frac{35}{23328}\theta^4+O_\Lambda(\theta^6).}          \tag{1}
\]

The same-area 1x4 loop has coefficient 73/58320. Shape therefore survives
the higher response tower. The error here is at fixed finite volume and
fixed loop, not uniform in loop area or a continuum limit.

Under an exact block of b fine links per coarse link, the free electric
heat closes with clock factor b. The interacting heat does not close:
its compressed semigroup has a nonzero second-order return term. This
term is retained by a well-defined Schur operator on the full compact
coefficient Hilbert space. It is not replaced by a finite-spin surrogate.

## Carrier, source order and scope

Use YM62's oriented SU(2)/native quaternion links, normalized real trace
W=Re(q)=Tr(q)/2, positive product reference Phi and Casimir
L_E=-sum_(e,a)D_(a,e)^2. Its one-link fundamental rate is 3/4, with no
factor 1/2. The plaquettes are counted once, with mu<nu. The magnetic
potential is V=sum_p W_p. The actual ground state is denoted h_Lambda,
its energy E_Lambda, and H=A-E_Lambda. None of these is replaced by
a truncated source expansion when a gap is invoked.

T1 uses YM62's periodic cubes and their subsequential history limits.
T2--T3 use finite periodic cubes divisible by the block factor, with
coarse side at least four. T4--T5 use finite open three-dimensional
boxes containing a fixed planar rectangle. These distinctions avoid
unproved boundary or winding identifications.

The ordered native algebra and reference calculus are inherited from
YM50--62. The source insertion method follows YM61, with a different
geometry and the independent Casimir of YM62. The elimination principle
already appears in YM49 and
[GE3's returning-memory theorem](../generalized-euler-evolution/GE3_RETURNING_MEMORY.md).
Here the free spatial projection reduces the unbounded electric
operator, while its interaction coupling is bounded. This verifies the
domain conditions for an all-content finite-volume extension of that
principle. No general moving-cut theorem is assumed.

The Hamiltonian lattice-gauge setting is established; see
[Kogut and Susskind, Phys. Rev. D 11 (1975), 395](https://doi.org/10.1103/PhysRevD.11.395).
For an explicit electric g^2 / magnetic 1/g^2 convention with a freely
displayed magnetic normalization constant, see equation (42) of
[Mathur and Sreeraj, arXiv:1509.04033](https://arxiv.org/abs/1509.04033).
Only that coefficient structure is used here; none of its weak-coupling
approximations is imported. Our trace and clock constants remain explicit.

## YM63-T1: a scaled gap has both lower and upper constraints

For the YM62 periodic vacuum, put m=mu(W_p), q=mu(W_p^2). The physical
trial vector (W_p-m)h has squared norm q-m^2>19/100. The exact ground
source identity gives energy mu(Gamma(W_p))=1-q<=1. The Rayleigh
principle and YM62's lower gap therefore yield

\[
 \boxed{\gamma(\theta)\le\Delta_\Lambda(\theta)\le100/19,
       \qquad |\theta|<1/8000.}                             \tag{2}
\]

Here Delta is the physical-sector gap above the vacuum, not an exact
plaquette eigenvalue. The same upper bound holds in each constructed
physical history representation: the time-zero vector has the same
variance lower bound, and its difference quotient in the heat semigroup
has energy at most one by YM62's local time estimate. This also verifies
the quadratic-form domain required by the infinite-volume trial bound.

Now introduce an explicit positive energy/clock scale kappa(a) at spatial
spacing a. Subtract the actual vacuum energy after scaling. Then

\[
 H_a=\kappa(a)H(\theta(a)),\quad \tau=\kappa(a)t,\quad
 \kappa(a)\gamma(\theta(a))\le\Delta_a\le100\kappa(a)/19.    \tag{3}
\]

An additive scalar in the uncentered Hamiltonian does not change (3).
For |theta(a)|<=theta_*<1/8000, gamma(theta(a))>=gamma(theta_*)>0.
Consequently kappa->0 forces the gap to zero; kappa->infinity forces
it to infinity. A finite positive limiting gap within this uniformly
controlled interior requires kappa to stay bounded above and away from
zero. This is a necessary scale condition, not a continuum existence proof.
No lower bound uniform in a follows by setting theta_* exactly at the
certificate endpoint, where the available gamma tends to zero.

The stronger time-correlation diagnostic is, for centered physical
bounded readouts F_a,G_a,

\[
 |\langle F_a h,e^{-tH_a}G_a h\rangle|
 \le e^{-t\kappa(a)\gamma(\theta(a))}
       \|F_a\|_\infty\|G_a\|_\infty.                      \tag{4}
\]

This follows from the complement norm bound and Cauchy--Schwarz and
also holds in the constructed history representation. If kappa>=c/a,
theta stays in the above interior, and the readout norms grow at most
as powers of 1/a, (4) tends to zero for each fixed positive physical t.
Indeed a^(-r)exp(-c gamma_* t/a)->0. This statement includes any chosen
bounded Wilson loops, even if their link lengths grow. It does not
exclude theories on other coupling trajectories or with other operator
normalizations; it identifies what this particular scaling can preserve.

For coefficient calibration, explicitly declare constants c_E,c_B>0
and the family, in natural energy-length units,

\[
 \mathcal H_{a,g}=\frac{c_Eg^2}{a}L_E
       +\frac{c_B}{ag^2}\sum_p(1-W_p),\quad
 \kappa=\frac{c_Eg^2}{a},\quad
 \boxed{\theta=\frac{c_B}{c_Eg^4}.}                         \tag{5}
\]

This is an explicit calibration choice; the native diagram does not
select c_E,c_B or a running g. Its constant shift is kappa theta |P|.
The certified window becomes g^4>8000 c_B/c_E. Any trajectory g->0
eventually leaves it, regardless of the positive fixed convention
constants. At fixed g inside the window, (3) diverges as 1/a.
Keeping kappa finite by g^2 proportional to a instead makes theta
grow as a^(-2), again leaving the window. None of these conclusions
asserts a gapless phase outside the certificate.

## YM63-T2: exact gauge-compatible product-link blocking

On a fine periodic cube of side bL, b>=2 and L>=4, retain the products

\[
 Q_{X,\mu}=q_{bX,\mu}q_{bX+\mu,\mu}\cdots
                  q_{bX+(b-1)\mu,\mu}.                    \tag{6}
\]

Distinct coarse links use disjoint fine links. Intermediate vertex
gauge factors cancel in order, so Q_(X,mu) transforms by the gauge
values at bX and b(X+mu). Thus pullback I_b f=f(Q) embeds coarse
gauge-invariant readouts into the fine physical subspace.

Independent Haar/native-reference links have independent Haar products.
Consequently I_b is an isometry of the completed coefficient Hilbert
spaces. Put P_b=I_b I_b* and Q_b=1-P_b. P_b is the free-reference
conditional expectation onto the retained product-link algebra.
It is not the interacting-vacuum conditional expectation. The fine
gauge group preserves this algebra and the reference, so P_b commutes
with its action.

For a polynomial f on one coarse link, differentiation at the k-th
fine factor rotates its three native derivative directions by the
adjoint of the fixed preceding product. The adjoint is orthogonal and
that preceding product is independent of the differentiated factor.
Summing the three squares gives exactly one coarse Casimir. Summing
the b fine factors proves

\[
 \boxed{L_{\rm fine}I_b=bI_bL_{\rm coarse},\qquad
 e^{-tL_{\rm fine}}I_b=I_b e^{-btL_{\rm coarse}}.}           \tag{7}
\]

The identity holds on each finite coefficient block. Their dense
orthogonal decomposition and the self-adjoint spectral closures give
the semigroup identity on the full completion. In particular P_b reduces
L_fine and preserves its operator and form domains. Blocking factors
compose associatively, I_(bc)=I_b I_c with consistent lattice origins,
and the free clock factors multiply. The b in (7) is not b^2 or 1/b.

## YM63-T3: the interaction returns through an exact memory operator

The union of retained fine paths is the subdivision of the coarse
cubic graph. Its shortest cycles have length at least 4b. Therefore
every fine elementary plaquette contains an unretained fine link.
That link is independent of all retained products and occurs linearly
in W_p. Averaging it annihilates the plaquette. Hence

\[
 \boxed{P_b V P_b=0,\qquad I_b^* A I_b=bL_{\rm coarse}.}     \tag{8}
\]

Equation (8) is a compression identity, not an interacting evolution
law. To see the missing term without a formal unbounded Taylor series,
define F_b(t)=I_b* exp(-tA)I_b. Since A is self-adjoint and bounded below,

\[
 F_b(2t)-F_b(t)^2
       =I_b^* e^{-tA}Q_b e^{-tA}I_b\succeq0.                \tag{9}
\]

The constant vector 1 is in D(A). Strong differentiation there gives
Q_b exp(-tA)1=-t Q_b A1+o(t)=t theta V+o(t), using (8).
Distinct plaquettes have zero free cross moment and Phi(W_p^2)=1/4.
It follows that

\[
 \boxed{\lim_{t\downarrow0}\frac{\langle1,
                [F_b(2t)-F_b(t)^2]1\rangle}{t^2}
              =\theta^2|P|/4>0\quad(\theta\ne0).}           \tag{10}
\]

Thus F_b is not a semigroup at nonzero coupling, although its first
derivative on the coarse core is -bL_coarse. The degree-two return
coefficient retains information erased by first-order compression.
It is a cut-coupling quantity; it is not identified with spacetime
curvature or with an independently defined thermodynamic lambda.

The missing evolution is concrete. In the I_b/Q_b decomposition write

\[
 A=\begin{pmatrix}C&B\\B^*&D\end{pmatrix},\quad
 C=bL_{\rm coarse},\quad B=-\theta I_b^*VQ_b,\quad
 D=Q_b A Q_b\big|_{Q_b\mathcal H}.                         \tag{11}
\]

Because Q_b reduces the free Casimir and V is bounded at finite volume,
D is self-adjoint on Q_b D(L_E); B is bounded with norm at most
|theta||P|. Also D>=3/4-|theta||P|, since Q_b is orthogonal to the
free constant. For z below both spectra, for example z<-|theta||P|,
ordinary block elimination on these domains gives

\[
 I_b^*(A-z)^{-1}I_b
 =\left[C-z-B(D-z)^{-1}B^*\right]^{-1}.                    \tag{12}
\]

The term subtracted in (12) is positive for such real z. Its time-domain
counterpart is K_b(t)=B exp(-tD)B*, with nonnegative quadratic form.
For initial hidden component zero, the retained component satisfies

\[
 x'(t)=-Cx(t)+\int_0^t K_b(t-s)x(s)\,ds.                   \tag{13}
\]

For general vectors interpret (13) in its mild semigroup form; for
domain vectors it is the differentiated equation. It follows from the
variation-of-constants formula for the two diagonal semigroups and
bounded B. D need not be positive, so no unshifted decay of K_b is
assumed. Formula (10) is exactly <1,BB*1>.

There is also an exact equation for the interacting coarse vacuum.
The full finite-volume gap Delta_Lambda>0 and P_b h!=0 imply

\[
 D-E_\Lambda\ge\Delta_\Lambda\|P_bh\|^2 Q_b>0.             \tag{14}
\]

Indeed for a unit vector y in ran Q_b, |<h,y>|^2<=||Q_bh||^2;
apply the full-space gap to y. Positivity of h makes its projection
onto the constant, and hence P_b h, nonzero. Here Delta_Lambda in
(14) is the full-space gap (at least YM62's gamma), distinct from the
physical-sector notation in (2). At z=E_Lambda block elimination gives

\[
 \boxed{[bL_{\rm coarse}-E_\Lambda
       -B(D-E_\Lambda)^{-1}B^*]I_b^*h=0.}                  \tag{15}
\]

This is the actual vacuum equation, with an energy-dependent generated
operator. Neither the overlap ||P_b h|| nor the extensive bound on B
has been controlled uniformly in volume and repeated blocks. Those
are explicit obligations before (15) can support a continuum argument.

## YM63-T4: surfaces select the first possible loop response

Let R be an m by n planar rectangle in a finite open three-dimensional
box, with area A_R=mn and boundary C. W_C is the normalized real trace
of its ordered boundary word. Under independent sign changes q_e->-q_e
the free reference and L_E are invariant. A nonzero source coefficient
for W_C therefore requires the odd-incidence plaquette set to have
boundary C modulo two. If a link has odd total multiplicity, its native
reference average vanishes; Casimir inverses preserve this parity.

Project that finite plaquette set onto the plane of R. Vertical faces
project to zero and parallel faces project to their planar unit cells.
The projected finite two-chain has boundary C. In the planar square
complex a finite two-chain with this boundary is exactly R: their
difference would be a nonempty finite two-cycle, which has an exposed
edge and thus cannot be closed. There must therefore be at least A_R
plaquette insertions. At equality all are parallel faces, one over
each cell of R; cancellation of shared edges and the prescribed
boundary force every face to be at R's height. The minimal surface is
the single planar sheet, with each face inserted once.

This proves that all orders below A_R vanish and that its coefficient
can be computed on that sheet alone, without exterior plaquettes.
It is the same coefficient in every finite open box containing R.
No finite-volume Taylor remainder bound uniform in the box is inferred.

The stronger sign involution on the open cubic lattice is
q_(x,mu)->eta_(x,mu)q_(x,mu), where
eta_(x,mu)=(-1)^(sum_(nu<mu)x_nu). Every elementary plaquette changes
sign. A rectangular boundary changes by (-1)^A_R, since its interior
edges cancel in the product of face signs. Unitary equivalence of
A(theta) and A(-theta), and uniqueness of the positive vacuum, give

\[
 \langle W_C\rangle_{-\theta}=(-1)^{A_R}
                                 \langle W_C\rangle_\theta. \tag{16}
\]

Thus the remainder after the leading term starts at A_R+2. This
particular center-sign assignment is not periodic on odd tori; no
such boundary extension is used in the response theorem.

The native integration constant is fixed, not fitted:

\[
 \Phi_q[\operatorname{Re}(Uq)\operatorname{Re}(q^{-1}V)]
                  =\tfrac14\operatorname{Re}(UV).           \tag{17}
\]

It follows directly from Phi(q_i q_j)=delta_ij/4 and the signed native
product. Peel faces successively from an outside corner of the sheet.
Contract a shared boundary link using (17); other common boundary
segments cancel in the ordered merged word. Repeat along successive
rows. This gives

\[
 \boxed{\Phi\!\left(W_C\prod_{p\in R}W_p\right)=4^{-A_R}.} \tag{18}
\]

The exact two-face polynomial check verifies that integration of their
single shared link gives one quarter of the 1x2 boundary loop. The
argument (17)--(18) handles arbitrary finite rectangles, without
expanding their exponentially many coordinate monomials.

## YM63-T5: two independent exact response recursions

Introduce a distinct scalar j_p for each face in R. In the normalized
vacuum series with Phi(h)=1, let h_S be the coefficient of the monomial
product_(p in S)j_p, with each label occurring once. Every nonempty
subset of this sheet has nonempty boundary, so its scalar vacuum-energy
coefficient vanishes by link parity. The coefficient equation is
L_E h_S=sum_(p in S) W_p h_(S minus p).

For the final loop response, a link whose two incident labels have
both been inserted may now be averaged. It appears in no remaining
insertion, and reference averaging on that link commutes with L_E
and its centered inverse. Averaging an unfinished link is forbidden.
After these legitimate projections, each active boundary link is
linear and all completed links are constant. The projected source is
therefore an exact Casimir eigenfunction with rate

\[
 \lambda(S)=\tfrac34|\partial S|.                           \tag{19}
\]

This is a projected coefficient calculation, not the assertion that
the full unprojected h_S has only this eigenvalue. Define rational
weights d_empty=1 and

\[
 d_S=\frac{\sum_{p\in S}d_{S\setminus p}}
                    {(3/4)|\partial S|},\qquad S\ne\varnothing.
                                                               \tag{20}
\]

The completed-link averaging is compatible with every subsequent
insertion in this sheet, so the bra/ket expansion and (18) give

\[
 \boxed{c_{m,n}=4^{-mn}\sum_{S\subseteq R}d_S d_{R\setminus S},
 \quad\langle W_C\rangle_\theta
           =c_{m,n}\theta^{mn}+O_\Lambda(\theta^{mn+2}).}     \tag{21}
\]

Normalization corrections cannot affect this leading term because
all lower numerator orders vanish. Each d is strictly positive.

There is a separate differentiated-energy calculation. Adjoin one
source label C for the boundary loop. Every fine link now occurs in
exactly two labels. The only closed subsets of these labels are empty
and full: without C a proper sheet subset has a boundary, and with C
closure forces inclusion of the entire sheet. For every proper nonempty
label set T use the recursion (20), with its own odd-link boundary.
At the full set omit the inverse and sum the predecessor weights.
Its energy coefficient is minus 4^(-mn) times that sum. Differentiating
with respect to the boundary source gives the same c_(m,n). The tests
compare these two recursions, and a full native two-face contraction
independently checks c_(1,2)=7/216.

| Rectangle | First order | Exact leading coefficient |
|---|---:|---:|
| 1x1 | 1 | 1/6 |
| 1x2 | 2 | 7/216 |
| 1x3 | 3 | 11/1728 |
| 1x4 | 4 | 73/58320 |
| 2x2 | 4 | 35/23328 |
| 2x3 | 6 | 57769/806215680 |
| 3x3 | 9 | 66145860409/77411255860592640 |

In particular c_(2,2)/(c_(1,1))^4=35/18. Treating plaquettes as
independent scalar responses loses the changing boundary lengths in
the inverse denominators. Equal area does not imply equal coefficient.
The code's area-12 enumeration guard is only a resource policy;
(18)--(21) are finite-rectangle statements for arbitrary m,n.

## YM63-T6: an observable matching coefficient is not a closed flow

Suppose, as a deliberately restricted comparison, one fits the mean
of a b by b fine loop with a coarse model bL_coarse-eta sum W_coarse.
Its normalized coupling is eta/b. YM62's elementary source response
then starts at eta/(6b). Matching only the first nonzero loop mean
in (21) forces

\[
 \eta_{\rm match}=6b\,c_{b,b}\theta^{b^2}
                       +O_\Lambda(\theta^{b^2+2}).           \tag{22}
\]

For b=2,

\[
 \boxed{\eta_{\rm match}=\frac{35}{1944}\theta^4+O_\Lambda(\theta^6),
 \quad (\eta/b)_{\rm match}=\frac{35}{3888}\theta^4
                                         +O_\Lambda(\theta^6).} \tag{23}
\]

The fit exists locally for fixed finite coarse volume since the
derivative 1/(6b) is nonzero. These are one-observable matching
coefficients, not a beta function and not an exact renormalization
group map. Equations (9)--(15) show why keeping only eta discards
interacting dynamics: the retained cut has a nonzero memory term and
an energy-dependent generated operator. Matching a static mean does
not match either of those.

The concrete next object is the Schur term in (15) under repeated
spatial blocks. A usable continuum route needs control of its range,
source derivatives and low-energy action, together with a nontrivial
spacing/coupling trajectory and renormalized observables. Repeating
the scalar coefficient (23) without that operator control is not such
a route. The existing weak-magnetic-ratio gap remains valid within
its stated window; no extrapolation to weak bare gauge coupling is made.

## Verification and physical gate

The certificate uses Python 3.12 exact fractions and the existing native
polynomial, multiplication, Casimir and reference engines. It checks
gauge endpoint cancellation, disjoint path geometry, the free clock
factor, a genuine unretained-link projection, the return coefficient,
the Schur sign, shared-link contraction, both surface recursions,
center signs, strict scale gates and claim/source tampering. The full
written proofs are distinct from these finite controls and have not
been mechanically formalized or independently peer reviewed.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym63_spatial_blocking_scale.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym63_spatial_blocking_scale.py -v
```

This chapter constructs neither a four-dimensional continuum field
theory nor a physical mass value. It also does not prove a Wilson area
law uniform in arbitrarily large loops. The
[Clay formulation](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf)
requires a nontrivial quantum Yang--Mills theory on R4 with a positive
gap. Exact lattice blocking and finite-loop response provide concrete
inputs toward that question, with the unresolved scale and operator
requirements now specified.
