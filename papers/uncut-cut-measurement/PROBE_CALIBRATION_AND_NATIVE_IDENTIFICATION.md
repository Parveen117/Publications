# Calibrated probes and identification of the native interaction

Monty Dabas. Research continuation v0.9, 25 September 2026.

U32 supplied a signed native order response and an exact information ledger.
To test whether different probes obey one continuation law, we must first ask
whether that response identifies the law. It does not: an exact reciprocal
ambiguity survives the complete signed response and its tensor. This extension
repairs that ambiguity, constructs a calibrated cross-probe test and gives
finite-error rejection/enclosure rules. It does not identify gravity by name.

## 1. Source audit and observation contract

The [six pinned sources](IDENTIFICATION_SOURCE_PINS.json) establish:

| Source | Result used and boundary |
|---|---|
| ID1, thermo residual non-selection | The Recognition square and declared action class do not fix numerical couplings or the potential; identifying a model is different from selecting it |
| ID2, thermo identification section | Complete transport plus independent scale/gauge calibration identifies its declared geometric target; a finite response subset is not automatically complete transport |
| ID3, thermo fixed-decoder theorem | Quantitative stability needs a uniform inverse bound and an independently justified error budget |
| ID4, LTB-1 | A reference/ledger fitted after seeing the outcome can erase the very mismatch being tested |
| ID5, U27–U29 | Native-compatible mixing and local R controls, explicit real carrier and transcript contract |
| ID6, U30–U32 | Signed order target, Gram tensor, target-faithful cut and exact response scale |

ID1–ID3 are from the corrected thermodynamic manuscript in open Publications
PR #4 at its immutable commit, not asserted merged. The native derivation
below uses finite algebra, not its spacetime metric, action or field equations.
The earlier Spectral I/II target-factorization and catalogue principles remain
the lineage for blindness and repair; no priority claim is made for them.

Use one **declared oriented pair** of channels, W=R^4. Let

    U_s = mixing Cayley step on channels (0,1),
    V_t = local R Cayley step on channel 0,
    t,s real, finite and nonzero;              L(t,s)=U_s V_t-V_t U_s.

These are U27's operators. Channel order, signed coordinate orientations,
normalization and the calibration frame are fixed before collecting data.
In this inverse problem t,s are fixed but unknown; the coordinate calibration
is known. They were supplied parameters in the earlier forward construction.
No physical clock or distance is used. The U27 convention is q=t or s;
RKF 51's corresponding parameter is h=2q.

Let e1 and e3 be the first coordinate of each channel (one-based indices).
On the same prepared e3, record the first two signed coordinates after each
word V,U and U,V. Their differences are d and e. Independently record the
first coordinate of U_s e1, called c. Thus

    d = (L e3)_1,       e = (L e3)_2,       c = (U_s e1)_1.

There are **three reported statistics from five raw scalar readings**: two
differences and one direct reading. This count is not an optimum over all
possible instruments. The extra direct reading is minimal relative to the
already available complete order-response target, as U33 proves.

The same-input preparation and signed-response interface is a model contract,
not an assertion that an unknown quantum state can be copied. Real physical
preparation, calibration uncertainty and raw-data acquisition remain to be
specified for an instrument.

## U33. Reciprocal blindness and its exact one-channel repair

Put k=2s/(1+s^2), a=(1-t^2)/(1+t^2), b=2t/(1+t^2). Then

    L(t,s) = [[0, k(T_t-I)], [k(T_t-I), 0]],
    T_t = [[a,-b],[b,a]],
    d=k(a-1),             e=kb,
    c=(1-s^2)/(1+s^2).

Consequently, the complete signed matrix L determines

    t=-d/e,              k=-(d^2+e^2)/(2d),

but generally does not determine s. Precisely,

    L(t,s)=L(t',s') iff t'=t and (s'=s or s'=1/s).

At s=1 or -1 the reciprocal alternatives coincide. For every other allowed s,
U_s and U_(1/s) are different transports with the same entire L and hence the
same response tensor Q(L), all its input responses, and its covariance ledgers
for the same ensemble. Their direct readings are c and -c.

Adding c removes the ambiguity:

    s=k/(1+c).

An exact real triple (d,e,c) belongs to this family iff

    d != 0, e != 0,       -1<c<1,
    c^2 + [(d^2+e^2)/(2d)]^2 = 1.

When it belongs, the inverse formulas give the unique pair (t,s).

**Proof.** Multiply the displayed Cayley blocks; the mixing cosine cancels
from UV-VU. The local rotation gives a-1=-2t^2/(1+t^2), b=2t/(1+t^2), so
-d/e=t. Also (a-1)^2+b^2=2(1-a), giving the formula for k.
Solving 2s/(1+s^2)=2s'/(1+s'^2) gives
(s-s')(1-ss')=0. The stated classification follows; the diagonal mixing
cosine changes sign under s->1/s. Its value cannot be recovered by any function
of L on a nontrivial reciprocal fibre. Thus at least one extra distinguishing
readout is necessary there, and c supplies one.

For sufficiency of the triple criterion, set t=-d/e and k as above. The circle
condition and strict c range give s=k/(1+c) != 0 and imply
2s/(1+s^2)=k, (1-s^2)/(1+s^2)=c. Substitution into the local rotation recovers
d and e. Necessity follows from the forward formulas. The argument is valid
over the reals; the executable fixtures use exact rationals.

**W39.** t=1/2 with s=1/2 or s=2 gives exactly

    (d,e)=(-8/25,16/25),        c=3/5 or -3/5.

The raw mixing matrices differ. Their complete order response and full Gram
tensor agree. Even a perfect measurement of L therefore needs the extra c
channel to choose between these continuation laws.

## U34. Calibrated cross-probe agreement and its scope

For each probe label p, let an independently fixed invertible frame F_p map
canonical preparations into its four-dimensional carrier. The canonical
readout uses F_p^-1; if it changes units, its metric is transported as in U30.
Assume the candidate class is exactly

    U_p=F_p U_(s_p) F_p^-1,      V_p=F_p V_(t_p) F_p^-1.

Run the protocol in those fixed frames. Then equality of the three reported
statistics for two probes is equivalent to equality of their two parameters,
and hence to equality of their full canonical U,V pair and every generated
word. Equality of order-response data alone does not imply this conclusion.

**Proof.** Pulling back the preparation and readout cancels F_p, leaving the
canonical three statistics. U33 is injective on this augmented observation,
so equal triples give equal t,s; substituting gives equal operators and words.
The converse is immediate. W39 proves the weaker observation fails.

This is a **proposed native probe-agreement gate**, named NP0. It checks whether
the selected probes, calibration and model class share one control pair.
It is not a derivation of the physical equivalence principle. Repeating any
chosen nonzero (t,s) across probe labels passes NP0; another parameter pair
also passes when repeated. Thus NP0 alone selects neither the numerical law,
an interaction graph, a physical source coupling nor gravity.

**W40: predeclared calibration can reveal disagreement.** The reciprocal
W39 probes pass every order-response comparison and fail the augmented test
because 3/5 != -3/5. The failure has a specified new readout, not a posterior
redefinition of the target.

**W41: fitted decoding can fake agreement.** Since U32 gives invertible L on
this block, for any two responses L_0,L_p an output-only decoder fitted as
D_p=L_0 L_p^-1 makes D_p L_p=L_0. This decoder is not the independently fixed
two-sided frame required above. Permitting this fit as a calibration would make
the order-response agreement diagnostic vacuous, exactly LTB-1's warning.

**W42: the declared catalogue matters.** A third unobserved native channel
may have a different local R step while all five raw readings on channels 0,1
agree. U34 identifies its four-dimensional family, not arbitrary hidden
extensions. U28 supplies the relevant full-catalogue closure obligations.

After fitting the three statistics, the inverse predicts the remaining U,V
matrix responses. Basis preparations not used in the fit are held-out model
checks. A mismatch rejects at least one calibration, carrier, preparation,
error-budget or candidate-family assumption. The finite certificate checks
these predictions on synthetic exact fixtures, not on empirical probe data.

## U35. Finite-error reconstruction and honest rejection

Let the reported statistics be (d_hat,e_hat,c_hat) with independent justified
component error bounds (epsilon_d,epsilon_e,epsilon_c). A difference of two raw
readings with error bounds alpha,beta has bound alpha+beta; calibration and
preparation errors must also be included. Independence here means the budget
is justified separately from the desired answer, not statistical independence
or a Gaussian-noise assumption.

Form exact real intervals D,E,C centred at the reported values with these
radii. If the divisions are separated from zero, interval evaluation of

    T=-D/E,
    K=-(D^2+E^2)/(2D),
    S=K/(1+C)

encloses every parameter pair in the declared family compatible with the data.
The interval for C^2+K^2-1 must contain zero; otherwise the family is rejected.
Containment of zero is only a necessary consistency test, not an existence
proof. Interval dependency can make these enclosures wider than necessary.

For example, if |e_hat|>epsilon_e, the phase estimate t_hat=-d_hat/e_hat obeys

    |t-t_hat| <=
    (epsilon_d |e_hat| + |d_hat| epsilon_e)
      / (|e_hat| (|e_hat|-epsilon_e)).

The same ratio bound applies to s=k/(1+c) once a bound on k has been obtained.
The margin condition is essential: exact identifiability does not guarantee
uniform precision as the response or inverse denominators approach zero.

**Proof.** Interval addition, multiplication, squaring and division away from
zero include the corresponding operations on every member, so composition
with the U33 formulas gives the enclosure. The circle condition is necessary
for every compatible model and therefore proves rejection if its interval
excludes zero. For the ratio bound, subtract a/b from a_hat/b_hat, bound the
numerator by epsilon_a |b_hat|+|a_hat|epsilon_b, and use
|b|>=|b_hat|-epsilon_b>0 in the denominator. No stochastic law is assumed.

For cross-probe comparison, disjoint reported-statistic intervals or disjoint
decoded parameter enclosures reject a shared pair. Overlap returns
**UNRESOLVED**, not proof of universality. Zero-error data satisfying the exact
gate may return **EXACT_SHARED_PAIR** within the stated family.

**W43.** In W39 the direct predictions differ by 6/5. If each probe's direct
error radius is epsilon, disjoint prediction-centred intervals require
2epsilon<6/5. At epsilon=3/5 they touch: the bound cannot reject a shared
reading solely from that coordinate. Smaller budgets give a robust separation.
For arbitrary measurement centres the code tests actual interval separation.

**W44.** A weak phase response can leave D or E straddling zero, giving
INSUFFICIENT_RESOLUTION. At very large |s|, c approaches -1 and the inverse
mixing denominator becomes ill-conditioned. No bounded-error guarantee over
all nonzero finite parameters is claimed. Increasing precision or changing the
predeclared protocol may repair the observation; declaring a PASS cannot.

## 2. Boundary and next experiment

This extension makes native coupling identification and cross-probe agreement
testable before choosing a spacetime representation. It exposes a blindness
that survived even the full U32 response tensor and repairs it with an exact
additional measurement. The source theorem on noise-stable identification
motivates the decoder margin; no spacetime assumptions are used to prove it.

The remaining physical step is to select actual preparations, readouts, sources
and probe species, then justify their relation to these native controls. A
gravitational interpretation needs discriminating physical response criteria
beyond NP0. A quantum reading still needs a probability/instrument adapter.
No such data or laws are invented by the synthetic fixtures. The general
proofs are under stated hypotheses, finite PASS is separately scoped, and
proof-assistant verification, external review and empirical validation are
not claimed. RH and Yang–Mills closure are untouched.
