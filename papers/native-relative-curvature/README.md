# QG-R1 — Native relative curvature without a clock

**Executed development 1.0, 14 October 2026.** A nondegenerate native
four-vertex cell now has a dimensionless source/reference curvature
quotient. Four scalar comparisons reconstruct it exactly. Independent
unitary frames at all four vertices leave it unchanged.
The return-defect information density also distinguishes divergence with
nondegenerate comparisons from a loop-chart pole or collapsing directions.

The [theorem and proofs](THEOREM.md) connect this finite observable to the
full variable-field curvature, derive the missing third-channel curl and
give a conditional obstruction to smooth continuation from concentration.
The result advances the quantum-gravity programme's relational geometry
layer. It supplies an observable and its transport interpretation; it does
not yet select a gravitational field law, quantum measure or physical metric.

## Main results

| Result | What changed |
| --- | --- |
| R1 | Source and reference links transform together. Every intermediate vertex frame cancels in the relative loop. |
| R2 | Relative-edge Gram area normalizes the relative-loop Cayley readout. Four scalar readings suffice; external step/clock/length labels cancel. |
| R3 | Four varying edges recover both derivative terms and the ordered commutator. The limit measures `F_source-F_reference`. |
| R4 | The inherited third word `L=RJ` carries a derivative contribution. An oriented loop reading recovers `d gamma = F_L-2 alpha wedge beta` in the stated reference frame. |
| R5 | The smooth nondegenerate limit is invariant under control-coordinate changes. Closed two-axis phase currents have calibration 1; a variable nonclosed example has 3/2. |
| R6 | The finite constant-cell benchmark is exactly **324/481**, with its finite-size correction derived explicitly. |
| R7 | The bounded return defect divided by squared native area defines an information-density observable. Its divergence excludes a C2 connection realization under the stated nondegenerate calibration. |

The benchmark's four readings are

```text
a = 3/5       relative first edge scalar part
b = 4/5       relative second edge scalar part
c = 12/25     cross-edge scalar comparison
d = 337/625   relative loop scalar part
```

They reconstruct `324/481` without a generating `q,p`, floating trigonometry,
a square root, or an external metric. This finite value is distinct from
its small-cell calibration of 1. The reference experiment remains part of
the data; replacing it need not preserve the result.

## Information concentration along a path

Attach a declared family of shrinking comparison cells to the path. Define

\[
\chi_{\rm ret}=\frac{1-\operatorname{sc}P_\square}{2},\qquad
\mathcal I_{\rm ret}=\frac{\chi_{\rm ret}}{\mathcal A^2}.
\]

The concentration condition is `limsup I_ret = infinity`, equivalently
`chi_ret` is not `O(A^2)`. Under a common C2 source/reference realization
with nondegenerate normalized comparisons, this cannot happen. The density
is defined even at the loop antipode, where the Cayley quotient fails.

An exact family with `D_x=tR`, `D_y=tJ`, `t=2^(-n)` and
`P=cayley(L,-t)` has `A^2=t^4`, shape `1`, and

\[
\chi_{\rm ret}=\frac{t^2}{1+t^2}\to0,\qquad
\mathcal I_{\rm ret}=\frac1{t^2(1+t^2)}\to\infty.
\]

Thus vanishing loop mismatch can still concentrate relative to area.
The first six densities run from `16/5` to `16777216/4097`; the theorem
proves the all-refinements lower bound `I_ret >= 4/(5t^2)`.
A second family diverges only because comparison directions collapse,
despite bounded normalized curvature. These mechanisms must be distinguished.

This is an operational return-defect density, separate from the signed
loop functional in the earlier information-invariance paper. Its physical
information interpretation and identification with gravitational singularity
remain open. Bounded density does not prove a regular extension.

## Executed evidence

**42 native algebraic certificates, 42 serialized replays, 327 exact checks,
21 negative controls.** The run includes:

- 35 admitted constant cells and 86 independent vertex-frame comparisons;
- 10 exact rational variable-field refinements, with targets `3/2` and `1`;
- six refinements each of bounded-density, concentrating and collapsing-direction families;
- a curved reference, an exactly flat endpoint-generated transport, and
  two full cells showing why the cross-edge scalar is necessary;
- rejection of degenerate area, antipodal Cayley readouts, nonunit links,
  inconsistent scalar data and mutated proof witnesses;
- a finite information-density reading at the loop antipode, with frame invariance;
- a control separating exact finite link transformation from pointwise
  resampling of a transformed continuum connection.

The prior phase-observation certificate and all six uploaded source packets
also pass read-only replay. No predecessor file or native engine is changed.
The new check count includes confluence/dagger audits and serialized replays;
negative controls are counted separately.

General frame and coordinate arguments, the C2 continuum limit, the
concentration bounds and the continuation obstruction have written proofs with explicit
domains in [THEOREM.md](THEOREM.md). Executed polynomial identities and finite
cases are recorded as such; they do not stand in for formalized real analysis.

## Run

From the Publications root, using **Python 3.12** and **Node v24.19.0**:

```sh
python3.12 papers/native-relative-curvature/verify.py --check
```

The default is also `--check`. It validates 17 direct source pins, reruns the
pinned predecessor (including all six original packets), extracts the same
unchanged native engine and reproduces the new certificate byte for byte.
No external packages or network access are needed by the run.

`--write` explicitly generates this stage's certificate and checksum only.
It never refreshes source pins or upstream results. The source engine remains
RKF `e8e3089745bf7dc0532d01c62a5f660f5658c6f0`.

- [Native certificate](CERTIFICATE.json) and [checksum](EXPECTED.sha256)
- [Source pins](SOURCE_PINS.json), [native implementation](verify.cjs),
  [Python orchestration](verify.py)
- [Dedicated CI workflow](../../.github/workflows/native-relative-curvature.yml)

## Quantum-gravity development route

| Stage | Retained output | Next issue it exposes |
| --- | --- | --- |
| [Clock-free native lattice](../Clock_Free_Native_Lattice_Source.zip) | Minimum declared cell and noncommuting transport | Which dimensionless parameters are observable? |
| [Cut quantisation](../Clock_Free_Cut_Quantisation_Source.zip) | Exact admitted cut/response reconstruction | Which source scales realize the cell? |
| [Physics cubits](../Clock_Free_Physics_Cubits_Source.zip) | Relations among missing scales | Separate source calibration from invariance |
| [Scale compass](../Clock_Free_Scale_Compass_Source.zip) | Controlled scale reconstruction | Recover repeated transport and branch data |
| [Clock relations](../Clock_Relation_Reconstruction_Source.zip) | Optional duration-comparison adapter | Retain phase aliases and coupling dependence |
| [Phase words](../Native_Phase_Word_Reconstruction_Source.zip) | Native word and orientation reconstruction | Recover accumulated winding |
| [Phase Memory / Angular Current](../native-phase-observation-certification/README.md) | Anchored phase record, gain/frame current identities | Specify full varying-edge transport |
| **QG-R1, this stage** | Relative-loop/area observable, third-channel curl and concentration obstruction | Select and realize a gravitational source/transport law |

The full connection identities and Gram-area method already existed in
EMK-C1 and YM54; the new work is the source/reference lattice readout,
four-scalar reconstruction, its variable-field matching and the carefully
admitted concentration criterion. The lineage
table in the theorem identifies reused results explicitly.
