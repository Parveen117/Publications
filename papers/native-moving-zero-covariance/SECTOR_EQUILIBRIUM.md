# MZC-3 — Sector admissibility and the common equilibrium carrier

15 October 2026. Uses CRC-1's existing native row-closure theorem and MZC-1/2.

## 1. Define the sector first

Let V be the admitted unitary circular comparison operator, with spectral
measure E_V. Declare an allowed Borel arc/set Sigma, its boundary convention,
orientation and reference mark before testing a state. Set

    P_Sigma=E_V(Sigma),  Q_Sigma=I-P_Sigma.

A state is sector-admissible exactly when P_Sigma x=x, equivalently Q_Sigma x=0.
For a positive normalized density matrix rho in finite dimension the equivalent
condition is rho=P_Sigma rho P_Sigma (or Tr(Q_Sigma rho)=0). Positivity of rho
is an explicit density-state contract; it is not a proof of the RH source sign.

Sector admissibility is the primitive acceptance predicate in this adapter.
It need not be defined by a scalar plus/minus label. This does not create a
canonical global order on a circle. An unspecified or retrospectively chosen
arc gives no test. Full phase and state information inside the arc is retained.

For the actual RH adapter the arc is already fixed by SE26/SE29: take the
upper oriented semicircle and marked balance point 1 of V=U_Delta. The chart
excludes -1. Under that source identification this admissibility corresponds
to the nonnegative comparison sector. The scalar sign criterion is recovered
as a coordinate theorem; renaming or enlarging Sigma cannot settle RH.

## 2. Put the cut and sector in the same record

Let X be a finite-dimensional complex inner-product carrier, C:X->Y a linear
cut reading, G a fixed skew-adjoint generator and W(t)=exp(tG). No quantization
of G or its outcomes follows merely from naming C an operator.
Define the stacked defect record

    D x = (C x, Q_Sigma x).

Its zero space E0=ker C intersect ran P_Sigma consists of sector-admissible
equilibria at the initial observation. The second entry retains sector leakage;
checking only the scalar cut reading would omit it.

The co-moving cut and sector are C_t=C W(t)* and P_t=W(t)P_Sigma W(t)*.
For x in E0, the transported state W(t)x always obeys C_t W(t)x=0 and
(I-P_t)W(t)x=0 by cancellation. This is frame covariance, not a dynamical
prediction. The fixed cut and fixed sector agree for that trajectory iff

    D W(t)x=0 for every real t.

## 3. Maximal common-equilibrium theorem

In dimension n, with this declared G, the exact common carrier is

    E_infinity = intersection_(k=0,...,n-1) ker(D G^k).

It is the largest G-invariant subspace contained in E0. In particular
E_infinity=E0 iff G E0 is contained in E0. Equivalently there is a well-defined
operator B on im(D) with DG=BD. This is CRC-1's quotient/intertwining criterion
applied to the combined cut-and-sector record, not a new sufficiency theorem.

Proof. If D exp(tG)x vanishes for all t, differentiating at zero gives DG^k x=0
for every k. Conversely Cayley-Hamilton expresses every G^k for k>=n in the
span of lower powers. The first n vanishing conditions therefore annihilate
the entire convergent matrix exponential. They also show that G preserves
E_infinity. Every invariant subspace of ker D satisfies all these conditions,
which proves maximality. The quotient B(Dx)=DGx is well-defined exactly when
ker D is invariant. This proof needs neither a scalar sign choice nor G=0.

The existing native rowClosure(D,[G]) computes the same space by its annihilator:

    O = basis span{rows(D G^k): k>=0},
    E_infinity=ker O,  dim E_infinity=n-rank(O).

At saturation OG=B_O O can be checked with the existing descendedAction.
For several declared generators use rowClosure(D,generators), which retains
all finite words, not just each generator's powers. Invariance then covers
finite compositions of their flows. This does not prove an infinite-dimensional
closure, an arbitrary unbounded-generator result or a continuum limit.

## 4. Sector preservation alone and genuine internal dynamics

For skew-adjoint G, the full allowed sector is preserved for all t iff
[G,P_Sigma]=0. Indeed invariance of ran P under the unitary group also preserves
its orthogonal complement, so W commutes with P; differentiating proves the
generator condition, and exponentiating proves the converse.

This condition does not ensure cut equilibrium: C may change inside the allowed
sector. Conversely a fixed cut may read zero while an unrecorded sector leaks.
The stacked D checks both. E_infinity may carry nontrivial internal rotations;
zero defect means zero in the specified measurement, not an absent state or
motionless source. Whether those internal distinctions matter for another
target is CRC-1's separate record-completeness question.

## 5. Exact controls with circular operators

Use V=diag(i,i,-i) and the declared upper arc. Thus P=diag(1,1,0).
Let e1,e2,e3 be coordinate columns. Except in case 4 use C=(0,0,1).
All generators are skew-adjoint, independent of the sector classifier V.

| Case | Generator action | Initial dim E0 | Persistent dim E_infinity |
| --- | --- | --- | --- |
| Internal circle | G e1=e2, G e2=-e1, G e3=0 | 2 | 2 |
| One leaking direction | G e1=0, G e2=e3, G e3=-e2 | 2 | 1 |
| Coupled circle planes | G e1=e2, G e2=-e1+e3, G e3=-e2 | 2 | 0 |
| Cut changes inside allowed arc | Internal-circle G, C=(1,0,0) | 1 | 0 |

In the coupled case x=e1 satisfies Dx=DGx=0 but DG^2x is nonzero.
Thus a first-order equilibrium check can miss later escape. The next retained
response level supplies the missing constraint; calling it a quantum jump or
identifying its order with a physical lambda requires a separate adapter.

Exact finite transport controls use W=(I-sG)^(-1)(I+sG), a unitary rational
Cayley step for real s. This step is not identified with exp(sG). Since it is
a rational function of G with invertible denominator, it preserves every
G-invariant carrier. General exponential persistence follows from section 3.
For these fixtures the certificate verifies closure ranks, kernel bases,
descended actions, internal motion, leakage and co-moving cancellation.

## 6. What this gives the RH application

For the actual bounded source comparison Delta, P is its already defined
allowed-sector projection. Choosing G=i Delta automatically commutes with P;
with C=Delta it freezes the exact kernel and preserves every spectral sector.
That is an available compatible construction, not independent evidence that
the arithmetic source has no forbidden sector. With a genuinely separate
source-derived generator, the combined record D supplies the concrete test.

We have not extracted such an independent RH generator, or a certified finite
restriction retaining its full source and tails, from the pinned sources.
The four fixtures are therefore exact structural controls, not new RH data.
The next source-specific gate is to specify that generator/domain and prove
its cut-and-sector intertwining, with the original arithmetic comparison intact.

Even a nontrivial E_infinity establishes compatibility only for those states.
It does not show that every source direction is admissible. RH's universal
sector-support obligation remains unchanged. Likewise an essential spectral
threshold at zero need not supply a nonzero exact kernel: approximate balances
require norm estimates and cannot be substituted into this exact-zero theorem.
