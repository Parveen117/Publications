# MZC-2 — Circular operators and sector-relative positivity

15 October 2026. Continuation of MZC-1, using existing UGD/EMK and RH sources.

## 1. The carrier is not a two-state object

UGD-1 already retains cyclic phase, scale and seam ledger, with a finite
structural-logic realization. We do not construct a replacement number system.
Its strict numeral carry operation and an associative operator algebra have
different equality contracts; UGD-A0 states the gate explicitly. We therefore
use the established associative native operator carrier for the calculation.

EMK-2 splits that algebra into even span{I,K} and odd span{R,RK}, where K^2=I,
R^2=-I and KR=-RK. These are algebra grades. They are not two possible values
of the full operator, and they are not automatically its positive/negative
spectral sectors. Likewise a projector's eigenvalues 0 and 1 specify membership
in a subspace, not the complete information or phase inside that subspace.
Neither circular logic nor this grading alone supplies an order on all operators.

## 2. Reuse the existing lossless circular operator

Let A=A* be the actual bounded comparison on its admitted positive pairing.
Fix a calibration a>0. SE26 Theorem 3 already proves

    U=(aI+iA)(aI-iA)^(-1),    U*U=I,
    A=-ia(U-I)(U+I)^(-1),    ker A=ker(U-I),
    F=(U-U*)/(2i)=2a R* A R,    R=(aI-iA)^(-1).

This is the Cayley representation, not a newly selected transport law. For
bounded A the spectrum avoids -1 and the chart is faithful. Circular phase
retains a continuum of readings, not only the signs. For a scalar eigenvalue t,
the circle coordinate is

    u_a(t)=(a^2-t^2+2iat)/(a^2+t^2).

The marked balance point is u=1. The positive spectrum maps to the open upper
semicircle, the negative spectrum to the open lower semicircle, and zero to 1.
The second real point -1 is the chart boundary, not a second finite zero.

The positive pairing, orientation, calibration and marked balance point are
part of the comparison contract. Reversing orientation or silently moving the
mark changes the sign convention/target; it does not prove a source inequality.

## 3. Positivity stated as a sector-support condition

Using the admitted spectral calculus, define Eplus, Ezero and Eminus as the
spectral projectors of A onto (0,infinity), {0} and (-infinity,0). Equivalently
use U's upper arc, marked point 1 and lower arc. These orthogonal projectors
sum to I, commute with A and U, and retain the full operator within each sector.
They need not be finite rank, and zero may be a spectral limit without a zero
eigenvector. Do not infer a positive spectral gap from Ezero=0.

Then

    A>=0  iff  F>=0  iff  Eminus=0.

Proof: SE26's invertible congruence proves the first equivalence. The spectral
representation writes <x,Ax> as the integral of t against its positive spectral
measure. If Eminus=0, this integral is nonnegative. Conversely, a nonzero lower
sector contains a nonzero spectral subspace with t<=-1/n for some positive
integer n, since these subsets exhaust the negative half-line. A vector there
has strictly negative quadratic reading. This proves the second equivalence.
The sector formulation is the desired positivity condition, not a proof that
the actual source obeys it. Even when Eminus=0, all phases and weights within
the allowed sector remain part of the operator record.

For a different proposed sector decomposition {P_j}, testing P_j A P_j>=0
is sufficient only when all sectors are included and cross blocks vanish
(equivalently, these are reducing sectors for A). In a finite orthogonal
decomposition the proof is <x,Ax>=sum_j <P_j x,A P_j x>. Without reduction,
A=[[1,2],[2,1]] has positive diagonal-sector compressions but the vector
(1,-1) has quadratic reading -2. Its eigenvalues are 3 and -1. The mixed
operator information cannot be discarded merely by naming sectors.

## 4. Circular composition does not automatically preserve the positive sector

At a=1 let A=2I. Its U=(-3+4i)I/5 is in the upper semicircle. But

    U^2=(-7-24i)I/25

lies in the lower semicircle; its inverse Cayley reading is -4I/3. By contrast,
the positive sum A+A=4I maps to (-15+8i)I/17 in the upper semicircle.
Thus circular composition is not source addition, and unrestricted phase
rotation does not preserve the chosen positive cone. This is an exact control,
not a counterexample to positivity of sums of positive operators.

Indeed no proper positive arc is invariant under every phase rotation: the
rotation orbit of any point is the entire circle. Either the reference moves
with the observer, or the allowed dynamics must preserve the required sector.
For observer relabeling U'=VUV* with unitary V, sectors transform by
E_j'=VE_jV*, so Eminus=0 is preserved. Nonunitary frame changes require MZC-1's
retained metric. A unit-modulus scalar alone drops phase; a phase alone drops
winding whenever the target depends on the lifted path.

## 5. Actual RH application and result

Use the actual SE29 normalized comparison Delta=(I-T)(I+T)^(-1), not an
arbitrarily chosen positive model. Let U_Delta be its SE26 Cayley image.
Under the pinned source hypotheses,

    Sb>=0  iff  Delta>=0  iff  E_lower(U_Delta)=0
             iff  1_(1,infinity)(T)=0.

For any unitary relabeling V, retaining all sectors gives
E_lower(V U_Delta V*)=V E_lower(U_Delta) V*. The absence of a negative sector
is therefore frame-independent under this stated equivalence.

This is the requested sector formulation of the existing RH positivity gate.
SE29 already identified this missing sector; we reuse that result. The new
integration here makes its circular realization, grading distinction, mixed
block gate and composition controls explicit in the moving-zero packet.

The result after imposing these structural conditions is an exact equivalent
criterion, not an RH proof. We have not established that the actual lower arc
has zero spectral support. Proving that from the arithmetic source, or proving
SE49's universal actual innovation positivity with its analytic bridge, remains
the substantive task. Discarding the lower sector would alter the source.

## 6. Evidence and scope

The native engine checks exact Cayley inversion/unitarity, oriented sign
readings, transported sectors, mixed-block failure and composition failure.
General spectral statements use the written proof and the already admitted
source spectral calculus; finite controls do not certify that analytic bridge.
No replacement UGD arithmetic, new universal logic, physical dynamics or
solution of RH/quantum gravity is asserted.

Sources: pinned SE26/SE29 in `sources/`; unchanged UGD-A0 contract,
UGD-1 numeral implementation and EMK-2 carrier implementation in
`../emk-ugd-algebra/`, all bound in SOURCE_PINS.json.
