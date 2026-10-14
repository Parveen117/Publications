# QG-R2 — Many centres, one observation base

**Executed development 1.0, 14 October 2026.** Native loops at different
sites can be compared at one base while retaining their source, reference
and connector records. Their combined loop can cancel completely even
when their retained information density diverges.

The companion [CID-D1 cut–information–curvature equation](CUT_INFORMATION_DYNAMICS.md)
now derives a covariant tensor from a moving observation cut. Its symmetric
part is a cut information metric and its antisymmetric part is induced
curvature. One transport equation evolves both, with an explicit mapping
to the earlier QTH-1 score tensor in the trace-one cut sector.

The [theorem and proofs](THEOREM.md) continue
[QG-R1](../native-relative-curvature/README.md). The motivating Kennedy
three-centre theorem concerns collinear instantaneous centres in planar
kinematics; it supplies no density-collapse rule. The native development
derives its own composition and concentration conditions.

## Main results

| Result | Output |
| --- | --- |
| F1 | Reversible transport to a common observation base; reference prefixes give the exact relative product of source and reference loop sequences. |
| F2 | Three-loop fusion is reconstructed by three scalar parts, three cross comparisons and an oriented triple reading. Six pairwise scalars alone can miss the answer. |
| F3 | For `m` unit loops, combined defect satisfies `chi(product) <= m sum chi_i`. Exact positive native squares prove it; complete cancellation remains possible. |
| F4 | A labelled block record retains every loop. Numerator/budget pairs merge associatively, and the resulting density is a positive weighted mean of child densities. |
| F5 | Receiving density separates into the budget ratio times the mean density. A measured composition factor connects retained and combined readings; a conditional regularity obstruction follows. |
| F6 | Three exact shrinking cells have diverging retained density and an identically zero combined defect, with independent comparison directions throughout. |

For the orientation witness, the triples `(R,J,L)` and `(R,J,-L)` have
identical individual scalar parts and pairwise comparisons. Their products
in observation order are respectively `I` and `-I`, with return defects
`0` and `1`. The missing oriented triple reading is respectively `+1`
and `-1`. The general orientation-loss principle already existed in
Recognition-Seam Topology; this is its specific native loop reconstruction.

## Concentration with a retained record

For positive squared native areas `b_i` and a receiving area budget `b_*`,

\[
S=\sum_i\chi_i,\qquad B=\sum_i b_i,\qquad
\overline I=S/B,\qquad I_{\rm ret,*}=(B/b_*)\overline I.
\]

The additive budget `B` is a record budget, not an assumed geometric union
area. Keeping it prevents an increase in record count from being silently
interpreted as a local singularity. The receiving budget must come from
an admitted comparison pair.

The exact three-loop family

\[
Q_1=\operatorname{cayley}(R,-t),\quad
Q_2=\operatorname{cayley}(J,-t),\quad
Q_3=(Q_2Q_1)^\dagger
\]

has, for `t=2^(-n)` and `b_i=b_*=t^4`,

\[
I_{\rm ret,*}=\frac{2(2+t^2)}{t^2(1+t^2)^2}\to\infty,
\qquad I_{\rm fused,*}=0.
\]

All three loops approach the identity, and every comparison shape is `1`.
The budget ratio is fixed at `3`. The document proves the all-refinements
bound `I_ret,* >= 64/(25t^2)` and the associated conditional obstruction.
Six rational refinements are executed independently of that analytic claim.

This supplies an observation and retention mechanism. Physical focusing,
gravitational geometry and a source law selecting the connectors and
receiving area remain open. A common observation base alone does not
physically bring distinct sites together or make them singular.

## Executed evidence

**22 native algebraic certificates, 22 serialized replays, 172 named exact
checks and 19 negative controls.** The new check count includes presentation
audits and serialized replays; negative controls are separate.

- 12 three-loop reconstruction cases and 34 frame comparisons;
- exact positive-square bounds and labelled block records for 1–7 loops;
- six refinements each of cancelling, aligned and increasing-record-count families;
- four rational moving cuts and five joint information/curvature evolution cases;
- controls for lost orientation, wrong ordering, omitted reference prefixes,
  changed connectors, invalid budgets, nonunit inputs and mutated witnesses.
- cut controls distinguish a motion generator from the cut itself, reject
  non-tangent directions and detect omitted connection derivatives.

QG-R1, the earlier phase-observation certificate and all six source packets
also pass read-only replay. The native engine and predecessor files are
unchanged. General induction, positivity consequences and C2/limit arguments
have written proofs; finite examples do not formalize those quantifiers.

## Reproduce

From the Publications root, with **Python 3.12** and **Node v24.19.0**:

```sh
python3.12 papers/native-loop-concentration/verify.py --check
```

The default is read-only `--check`. It validates 26 direct source pins,
reruns the predecessor chain and reproduces this certificate byte for byte.
No third-party mathematical packages or network access are required.
`--write` explicitly generates only this stage's certificate/checksum;
it never updates source pins or predecessor results.

- [Written proofs and source audit](THEOREM.md)
- [Covariant cut–information–curvature equation](CUT_INFORMATION_DYNAMICS.md)
- [Native certificate](CERTIFICATE.json) and [checksum](EXPECTED.sha256)
- [Source pins](SOURCE_PINS.json), [native implementation](verify.cjs),
  [Python orchestration](verify.py)
- [Dedicated CI workflow](../../.github/workflows/native-loop-concentration.yml)

The next physical obligation is a source/transport law that determines
receiving geometry and composition while retaining enough records to
distinguish cancellation, increasing multiplicity and local concentration.
