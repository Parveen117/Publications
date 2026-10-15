# HCM-2 — Phase-aware response transport and the UGD interface boundary

15 October 2026. Extends the evaluated HCM-1 source without changing UGD arithmetic.

## 1. What the existing UGD implementation actually supplies

UGD-1 has a finite cyclic phase alphabet (K=12 by default, K=4 admitted),
integer scale indices, seam digits and an integer seam ledger. Phase is stored
as an exponent; the implementation does not evaluate a continuous root-of-unity
wave amplitude. Its classical projection defaults to sum a_n 2^n, not
sum exp(2 pi i a_n/K) 2^n. Its neutral cut-zero digit is (a,n,s)=(0,0,0).
Sending that whole digit naively to 2^n i^a would send its zero to 1.

There is another operation mismatch. In the pinned source, mul_digits uses
phase exponent (a*b) mod K, whereas cyclic phase composition uses (a+b) mod K.
For K=4 and a=b=1, these give 1 and 2 respectively. Thus the phase-label map
a -> i^a respects cyclic composition, but is not a multiplicative homomorphism
of that existing numeral operation. We do not change or silently substitute
either operation. UGD-A0 separately records the strict carry-associativity gate.

Consequently a faithful map from the entire existing UGD numeral algebra to
the two heat-response amplitudes has NOT been established. The construction
below is a new typed response-record adapter using its existing quarter-turn
labels and scale base, with an explicitly retained residual. It is faithful
for the two real HCM channels only. An independent legacy seam ledger, if part
of the source, must be retained separately; it cannot be inferred from them.

## 2. Exact source evolution in the native EMK operator carrier

Set z=q+i h for HCM-1's visible and hidden readings. Its source-derived law is

    z' = alpha z + beta conjugate(z),
    alpha=-65/1296,    beta=(7+24i)/1296.

Proof: substitute q'=-(29/648)q+h/54 and h'=q/54-h/18, collect the real and
imaginary parts, and then write q=(z+conjugate(z))/2 and
h=(z-conjugate(z))/(2i). No new source parameters are selected.

In the native real two-mode representation let Rz=iz and Kz=conjugate(z).
They obey R^2=-I, K^2=I, KR=-RK. Then

    G=-(65/1296)I+(7K+24RK)/1296,
    S=(7K+24RK)/25,     S^2=I,
    Pslow=(I+S)/2,      Pfast=(I-S)/2,
    exp(tG)=exp(-5t/162)Pslow+exp(-5t/72)Pfast.

The native word proof checks (7K+24RK)^2=625I. The sector formula follows by
the complementary idempotents and gives the same rates as the actual YM54
coefficient block. These are spectral sectors with retained amplitudes, not
two possible states of the whole operator.

The operator is real-linear and is not complex-linear in this chosen record:
the coefficient beta is nonzero. Thus a single scalar complex exponential
cannot transport every initial z. A quarter-turn rephasing does not leave this
fixed anisotropic law unchanged; its conjugate coefficient must transform.

## 3. Amplitude and phase transport are coupled

For z=r exp(i phi), r>0, in a local unwrapped angular chart, the exact equations
are

    r'/r = [-65+7 cos(2phi)+24 sin(2phi)]/1296,
    phi' = [24 cos(2phi)-7 sin(2phi)]/1296.

They follow from z'/z=alpha+beta exp(-2i phi), taking real and imaginary parts.
At z=1 the phase derivative is 1/54, so a phase initially in the finite K=4
alphabet immediately moves between its symbols. A fixed finite alphabet alone
cannot retain this continuous, nonconstant phase exactly. An integer winding
or seam count does not contain the missing within-sector angle.

The equation is for the actual two-channel classical heat response encoding.
It is not a Born-rule identification, a qubit state construction or a physical
quantum evolution. The complex ratio of amplitudes and IRT's real ratio of
expectation readings have different declared targets.

## 4. A reversible quarter-turn/scale/residual record

For nonzero z choose the unique integer n satisfying

    4^n <= |z|^2 < 4^(n+1).

Choose a in {0,1,2,3} so z lies in the half-open quadrant
[a pi/2,(a+1)pi/2), with angles read modulo 2pi. Let

    w=2^(-n) i^(-a) z.

Then Re w>0, Im w>=0 and 1<=|w|^2<4. The exact record (a,n,w) reconstructs
z=2^n i^a w. The half-open conventions make it unique. Zero has a separate
zero tag and no phase. Rational q,h give a rational pair w without computing
an angle, root or logarithm; comparisons of |z|^2 and component signs suffice.
For nonrational completed values this is a mathematical chart construction,
not an effective exact-decision algorithm for arbitrary real inputs.

This record uses UGD's quarter-turn labels and binary scale but adds the full
residual w. It does not assert that w is encoded by the historical finite
numeral operations. Equal labels a,n can hide different w and different future
readings; both residual components are necessary on the unrestricted source.
The residual is measured/retained information, not information manufactured by
renaming a hidden channel as phase.

Within a chart interval where a,n are constant,

    w'=alpha w+(-1)^a beta conjugate(w).

Indeed i^(-a) conjugate(i^a)=i^(-2a)=(-1)^a. A change of chart uses

    w_new=2^(n_old-n_new) i^(a_old-a_new) w_old.

This preserves the decoded z exactly. Crossing a quadrant or scale boundary
requires this transition, not an assertion that the discontinuous labels have
a smooth derivative. The native replay verifies the transition and the local
operator intertwining for every quarter-turn label and positive/negative scales.

## 5. Ledger and winding cannot be identified by name

The existing phase-observation chapter requires an anchor and additional path
information to recover a lifted angle; its doubled-channel convention is not
automatically the single complex phase used here. UGD's seam ledger is a
separate algebraic record, not a derived angular-current measurement.

For this particular autonomous heat law there is a stronger restriction.
For any nonzero initial vector x0=(q0,h0), Hobs is positive definite and

    x0^T exp(-t Hobs) x0 > 0 for every finite t>=0.

Therefore the angle between x(t) and x0 stays in the open interval
(-pi/2,pi/2). The trajectory cannot hit zero at finite t or complete a full
winding. Its phase may change, and a chosen chart label may change, but an
unwrapped lift anchored at x0 has no full-turn accumulation. An arbitrary
winding ledger cannot be extracted from this heat trajectory.

Proof: positive eigenvalues of exp(-t Hobs) give the strict quadratic inequality;
invertibility excludes a zero trajectory. The positive dot product keeps the
direction in the initial open half-plane, which admits one continuous anchored
angular chart. The claim concerns this fixed autonomous source, not arbitrary
moving protocols or all UGD paths.

## 6. What is established and what remains

Established: an exact EMK decomposition of the published heat source, its
conjugate-coupled amplitude/phase law, a reversible phase/scale/residual record,
and correct chart transport. The known two-channel memory repair is now tied
to a precise phase-aware representation.

Not established: a homomorphism of the full historical UGD numeral algebra,
identification of seam ledger with phase winding, or recovery of an unknown h
from q alone. The direct finite-label bridge fails on this source without
additional residual data. A future full UGD interface must explicitly select
an equality/product contract and prove operator descent, as UGD-A0 requires.
