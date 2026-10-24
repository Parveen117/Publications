# EGR-1 — Invariant source strength, exponential cut transport and curvature response

14 October 2026. Continuation of N4H-1, responding to the proposed bridge
between quantised cuts, continuous exponential transport and information
concentration. All variables below are dimensionless.

The result is a concrete nine-simplex model. Ten interior edge variations
reduce to two curvature coordinates after removing eight flat-background
vertex-displacement directions. Their exact stiffnesses are 10 and 47.
A selected invariant source coupling gives a linear curvature response in
its strength N. A selected relaxation of the observation cut then gives
an exact exponential law for the ratio of those two source channels.

The geometry calculation is new to this development sequence. Exponential
transport, cut identities, gradient relaxation and linearized Regge gauge
directions have established precedents. No literature-wide novelty is claimed.

## 1. Precise meaning of the proposed bridge

Take a finite positive native carrier and a self-adjoint rank-one cut
P^2=P, with trace P=1. Encode a calibrated source weight as

\[
\rho=NP,\qquad N\geq0.
\]

N=Tr(rho) is invariant under unit-frame changes and unitary transport.
It is a declared source strength, not an assertion that trace equals every
notion of quantum information or entropy. P has the discrete spectrum
{0,1}, while its orientation can vary continuously.

For an anti-dagger generator V, the transport equation is

\[
U'=VU,\qquad P(\lambda)=U(\lambda)P(0)U(\lambda)^\dagger.
\tag{1}
\]

For a constant V, U=Exp_Sigma(lambda V); for a variable V, use the ordered
exponential. The derivative is P'=[V,P]. Dagger and product differentiation
prove U^dagger U=I and preservation of P^2=P, trace, rank and the spectrum
of rho. The parameter lambda orders operations; no physical clock is used.
The scalar completion and exponential are inherited from GE1 and NJC-1.

Thus the spectrum of a cut can be discrete while its transport is continuous.
This is the precise bridge used here. A continuum limit of spacetime is a
different requirement and is not obtained by merely writing an exponential.

Two qualifications follow directly from the algebra:

1. Exp(lambda P)=I+(exp(lambda)-1)P commutes with P. Exponentiating the cut
   itself leaves that cut fixed under conjugation. Motion requires a
   generator with a nonzero commutator with P.
2. Continuous exponential transport does not imply flatness. Transports
   in two directions can fail to commute and have nontrivial loop return.
   The exponential is a frame-covariant transport operator; its preserved
   spectral quantities and conjugacy-class readings are invariants.

The native certificate checks both points. Sampled return data also leave
generator ambiguity: exp(0)=exp(2pi R)=I for R^2=-I, although the paths
between the samples differ. N4H-1's curvature-lift record remains necessary.

## 2. A four-dimensional complex with ten interior edges

Begin with a Euclidean four-simplex with vertices

\[
A=0,\quad B=e_1,\quad C=e_2,\quad D=e_3,\quad E=e_4.
\]

Insert its barycentre O=(1,1,1,1)/5 and subdivide into five four-simplices.
Then insert P_geo=(6,6,6,6)/25, the barycentre of OBCDE, and subdivide
that cell into five. Here P_geo is a **vertex label**, distinct from the
observation cut P. The resulting nine cells are

\[
OACDE,\ OABDE,\ OABCE,\ OABCD,
\ P_{geo}BCDE,\ P_{geo}OCDE,\ P_{geo}OBDE,
\ P_{geo}OBCE,\ P_{geo}OBCD.
\]

The external boundary is the original five tetrahedra. There are ten
interior edges and twenty interior triangular hinges. All cells are
nondegenerate. Four have volume 1/120, five have volume 1/600, and their
sum is 1/24, the original volume. This is an actual flat subdivision,
so every interior deficit angle vanishes at the base geometry q_0.

Use squared edge lengths q and the Regge bulk-plus-boundary curvature
functional R from N4H-1. Its boundary convention and positive overall
sign are unchanged. Boundary edge lengths remain fixed. Write

\[
\mathcal R(q_0+dq)=\mathcal R(q_0)
+\tfrac12dq^T\mathsf H\,dq+O(\|dq\|^3).
\tag{2}
\]

There is no linear interior term because the base deficits vanish.

### Exact Hessian and curvature derivatives

For each cell use its four-edge Gram matrix G, its inverse, and the
5-by-5 Gram matrix B of the affine barycentric gradients. For an edge
variation, differentiate the inverse exactly:

\[
(G^{-1})'=-G^{-1}G'G^{-1}.
\]

N4H-1's angle-derivative identity gives, for the hinge opposite vertices
i,j in that cell,

\[
a_h\theta'_{ij}=12v_\sigma
\left[B'_{ij}-\frac{B_{ij}}2
\left(\frac{B'_{ii}}{B_{ii}}+\frac{B'_{jj}}{B_{jj}}\right)\right].
\tag{3}
\]

All entries on the right are rational at these coordinates. Assemble

\[
\mathsf C_{he}=a_h\frac{\partial\delta_h}{\partial q_e},
\qquad
\mathsf H_{ef}=\sum_h
\frac{\partial_e(a_h^2)}{2a_h^2}\mathsf C_{hf}.
\tag{4}
\]

This computes the area-weighted deficit derivative and the Hessian without
approximating angles or square roots. The executable derives (4) from all
nine cells, verifies each cell's Schläfli cancellation, and records all
individual angle contributions. Neither H nor its rank is inserted as an
assumption.

Let G_v be the 10-by-8 differential of edge squares under moving the two
interior vertices in the four ambient directions:

\[
(G_v)_{(ij),(v,a)}
=2(x_i^a-x_j^a)(1_{i=v}-1_{j=v}).
\]

The exact results are

\[
\operatorname{rank}G_v=8,\qquad
\mathsf H G_v=0,\qquad \mathsf C G_v=0,\qquad
\operatorname{rank}\mathsf H=\operatorname{rank}\mathsf C=2.
\tag{5}
\]

These are linearized vertex-displacement directions at the specified flat
background. They do not prove nonlinear diffeomorphism constraints or an
anomaly-free quantum constraint algebra. The two remaining directions are
the curvature coordinates of this finite complex, not a calculation of
continuum graviton polarizations.

## 3. The two diagonal curvature coordinates

Define a 2-by-10 map A through

\[
\begin{aligned}
y_1&=\sum_{v=A,B,C,D,E}dq_{Ov},\\
y_2&=dq_{OP_{geo}}+\sum_{v=B,C,D,E}dq_{P_{geo}v}
-\tfrac15\sum_{v=B,C,D,E}dq_{Ov},\qquad y=A\,dq.
\end{aligned}\tag{6}
\]

The barycentre identities give A G_v=0. These are invariant **to first
order at q_0** under the stated vertex displacements; they are not asserted
to be exact nonlinear observables at every metric.

The full native calculation gives

\[
\boxed{\mathsf H=A^TKA,\qquad K=\begin{pmatrix}10&0\\0&47\end{pmatrix}.}
\tag{7}
\]

Let S inject y into dq by setting dq_AO=y_1, dq_OPgeo=y_2 and the other
eight edge variations to zero. Then AS=I and K=S^T H S. Also

\[
\mathsf C=LA,\qquad L=\mathsf C S,\qquad \operatorname{rank}L=2.
\tag{8}
\]

Thus a nonzero y gives nonzero linearized curvature, and every curvature
variation factors through these two coordinates. The ten-edge Hessian has
nonzero mixed entries; the coordinates (6) diagonalize its physical
quadratic part. The numbers 10 and 47 depend on the geometry and this
normalization. They are not universal constants or a quantized area spectrum.

## 4. Invariant amount produces proportional curvature in a specified coupling

Use the real two-component cut sector and fixed self-adjoint readout
operators

\[
Z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
X=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
r(P)=\binom{\operatorname{Tr}(PZ)}{\operatorname{Tr}(PX)}.
\]

A real rank-one cut is

\[
P=\tfrac12(I+r_1 Z+r_2 X),\qquad r_1^2+r_2^2=1.
\tag{9}
\]

An observer-frame change transforms P, Z and X together, preserving r.
The chosen readouts calibrate how the source couples to the two geometric
coordinates. N alone does not specify this directional information.

For positive gamma and nonnegative alpha, select the quadratic coupling

\[
\boxed{\Phi(y,P)=\frac\gamma2 y^TKy-\alpha N\,y^Tr(P).}\tag{10}
\]

This bilinear source term is an explicit model choice. It is not obtained
from invariance alone, and it is not an identification of an arbitrary
information statistic with the stress-energy tensor. In particular it
differs from the positive mismatch functional considered in N4H-1.

The stationary metric equation and its unique solution are

\[
\boxed{\gamma Ky=\alpha N r(P),\qquad
y=\frac{\alpha N}{\gamma}
\begin{pmatrix}r_1/10\\r_2/47\end{pmatrix}.}\tag{11}
\]

In the original variables the source is j=alpha N A^T r, so G_v^T j=0.
All ten linearized metric equations H dq=j/gamma hold for dq=Sy and for
every representative differing by G_v xi. An arbitrary edge force need
not satisfy this compatibility condition; the certificate provides a
counterexample.

**Strength-response theorem.** At fixed P, alpha, gamma and background,
multiplying N by t>=0 multiplies y and every component of the first-order
area-weighted deficit L y by t. Define the positive curvature response score

\[
\boxed{\mathcal B^2=y^TKy
=\left(\frac{\alpha N}{\gamma}\right)^2
\left(\frac{r_1^2}{10}+\frac{r_2^2}{47}\right).}\tag{12}
\]

It scales by t^2. This proves a precise version of "more invariant source
strength gives more bending" for this coupling and fixed directional data.
It is exact in the quadratic model. Extending a large-N extrapolation to
the full nonlinear geometry is not justified by this local calculation.

Conversely, at the same N a source aligned with the first channel has
B^2=(alpha N/gamma)^2/10; alignment with the second gives the same factor
divided by 47. Source magnitude alone therefore does not determine bending.
Adding a constant to a source potential also leaves its derivative and
force unchanged. These prevent a universal rule based on the numerical
value of any chosen invariant without its coupling law.

## 5. Exponential alignment of the cut and increasing curvature response

Eliminate y from (10). The effective cut functional is

\[
\Phi_{\rm eff}(P)=-\frac\beta2
\left(\frac{r_1^2}{10}+\frac{r_2^2}{47}\right),\qquad
\beta=\frac{\alpha^2N^2}{\gamma}.
\tag{13}
\]

Its ambient first variation is Tr(H_eff dP), where

\[
H_{\rm eff}=-\beta(r_1Z/10+r_2X/47).
\]

Choose the cut-tangent **gradient relaxation**

\[
\boxed{P'=[[P,H_{\rm eff}],P],\qquad
V=[P,H_{\rm eff}],\qquad U'=VU.}\tag{14}
\]

This choice supplies attraction. It replaces the conservative second-order
choice in CID-D2 for this model; attraction does not follow from an
exponential alone. With a general conserved kinetic term the cut could
oscillate instead. The state-dependent flow (14) is also not asserted to
be a derived linear Schrodinger law.

**Preservation and monotonicity.** V is anti-dagger and P' is self-adjoint.
Using P^2=P gives PP'+P'P=P', hence the cut constraint is preserved.
Trace cyclicity gives

\[
\frac{d\Phi_{\rm eff}}{d\lambda}
=-\|[P,H_{\rm eff}]\|_F^2\leq0.
\tag{15}
\]

Equation (1) supplies the unitary lift, preserving N and the spectrum of
rho. In the real two-mode sector set

\[
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
g=\frac1{10}-\frac1{47}=\frac{37}{470}.
\]

Direct native multiplication gives

\[
V=-\beta g r_1r_2 R,\qquad
r_1'=2\beta g r_1r_2^2,\qquad
r_2'=-2\beta g r_1^2r_2.
\tag{16}
\]

On a branch with r_1!=0 put m=r_2/r_1. Combining (16) and the cut circle
constraint yields

\[
\boxed{m'=-\frac{37\beta}{235}m,
\qquad m(\lambda)=m(0)
\exp\left[-\frac{37}{235}\frac{\alpha^2N^2}{\gamma}\lambda\right].}
\tag{17}
\]

This is an exact exponential ratio law, not a numerical fit. The cases
r_1=0 are fixed points, and must not be evaluated by division. If beta=0
the whole flow is stationary. For beta>0 and r_1(0)!=0, the sign of r_1 is
preserved and the cut approaches the corresponding first-channel cut.

Write r_1=cos(2phi), r_2=sin(2phi) with a continuous branch. All generators
in (16) are proportional to the same R, so the ordered exponential reduces to

\[
\boxed{U(\lambda)
=\exp\big[(\phi(\lambda)-\phi(0))R\big].}\tag{18}
\]

The preserved cut spectrum remains {0,1}. The continuous change is its
orientation and the metric response y(P). Finite-parameter transport is
invertible; its limiting alignment need not retain the initial ratio if
the transport record is discarded.

Let chi=alpha N/gamma. From (16),

\[
\boxed{\frac{d\mathcal B^2}{d\lambda}
=4\chi^2\beta\left(\frac{37}{470}\right)^2r_1^2r_2^2\geq0.}\tag{19}
\]

Thus this selected feedback simultaneously preserves the source strength,
moves the cut continuously and increases its curvature response score.
The role of the exponential is now explicit and connected to the
calculated geometric response operator.

### An exact finite transport witness

Choose source orientations

\[
(r_1,r_2)_{initial}=(7/25,24/25),\qquad
(r_1,r_2)_{later}=(119/169,120/169).
\]

Their ratios satisfy m_later/m_initial=5/17. For beta>0 they lie on the
same solution (17), separated by
Delta lambda=(235/(37 beta)) log(17/5). The exact transport is

\[
U=\frac{63}{65}I-\frac{16}{65}R
=\operatorname{cayley}(R,-1/8).
\]

The native engine verifies U P_initial U^dagger=P_later, unitarity, the
ratio factor, both instantaneous derivatives, and strict increase of the
response score. This checks the bridge with exact rational entries rather
than relying on an exponential approximation.

## 6. Concentration and singularity are different conclusions

For finite N, alpha, gamma>0, the exact regular model has

\[
\frac{\chi^2}{47}\leq\mathcal B^2\leq\frac{\chi^2}{10}.
\tag{20}
\]

The cut stays bounded, V stays bounded, and the exponential flow is smooth
at every finite parameter. The limit in (17) is a finite projector, not a
singularity. In the tested small-response geometries every cell Gram
matrix remains positive; those finite admissibility checks are not a
solution of the full nonlinear Regge equations.

A density statement additionally needs a receiving region and a measure.
For example N/v can diverge at fixed positive N only if the specified
receiving volume v tends to zero. The earlier QG-R1/R2 criterion instead
uses a return-defect numerator divided by squared comparison area and
tests its asymptotic rate. These observables must not be identified without
an explicit comparison theorem.

To prove that a singularity moves along a cut direction one must also
map that direction to spatial cells and derive their motion or collapse.
Equation (14) moves an observation cut in its carrier; it does not yet
provide that spatial focusing map. The present result is exponential
alignment with a bounded curvature response. Singular focusing remains
a nonlinear geometric question.

## 7. What the certificate establishes

The executable includes the complete rational nine-cell geometry, all
cell angle derivatives, the 10-by-10 Hessian, the 10-by-8 displacement
matrix, the twenty-hinge curvature derivative, exact quotient
factorizations, source covariance, force equations, rational flow
snapshots and negative controls. Symbolic witnesses and critical-pair
audits use the unchanged native Pāṇinian compiler. Python 3.12 handles
only file I/O and process launching.

The exponential completion, real-angle interpretation, geometric
subdivision proof, arbitrary-parameter flow and asymptotic statements
retain their written hypotheses. Finite symbolic coefficients are not
used as a substitute for the all-order analytic proof. This packet is
not a general-purpose proof-assistant formalization or a certificate of
complete quantum gravity.

| Claim | Result |
|---|---|
| Discrete cut spectrum can have continuous transport | Established under the native exponential and anti-dagger generator |
| Larger invariant strength increases curvature | Established for fixed source direction and the selected linearized coupling |
| Exponential alignment increases response at fixed strength | Established for the explicit relaxation (14) |
| Exponential of P itself moves P | False: it commutes with P |
| Exponential continuity guarantees flatness | False: loop transports can fail to commute |
| Any larger information statistic determines more gravity | Not established; coupling, directional data and derivatives matter |
| This bounded cut flow transports a physical singularity | Not established; the required nonlinear focusing law is absent |

Lineage: N4H-1 supplies the Regge geometry and Schläfli adapter; NJC-1
supplies typed exponentials and the normal-transport distinction; CID-D2
supplies the cut-tangent variational language; GE1 supplies the native
completion and continuity; QG-R1/R2 supply the concentration/regularity
distinction. No previous engine or certificate is modified.

For the established linearized Regge context, see P. A. Hohn,
*Canonical linearized Regge calculus: Counting lattice gravitons with
Pachner moves*, [Phys. Rev. D 91, 124034](https://doi.org/10.1103/PhysRevD.91.124034).
The paper's linearized gauge framework is background context; the nine-cell
matrices above are independently computed in the native engine. Its results
are not imported as a proof of our nonlinear or quantum constraint closure.
