# Thermo-compass foundations for the physics programme

**Current development entry point · 29 September 2026**

Start from the existing T–V–S–P compass, constrained thermodynamic responses and covariant Jacobian tower. Construct transport and its geometry before attempting numerical fundamental constants. The [corrected thermo manuscript](../thermodynamic-response-corrections/README.md) is preserved on this research branch at its exact PR #4 source edition.

## Read first

- [TC-1–TC-6: derivations, proofs and exact examples](THEOREM.md)
- [Important earlier results and their repository locations](RESULTS_INDEX.md)
- [Source provenance](SOURCE_PINS.json) and [recovered-draft provenance](RECOVERY_MANIFEST.json)
- [Generated finite evidence](CERTIFICATE.json) and [read-only verifier](verify.py)

## This step's result

For a supplied stable U(S,V), the response Hessian yields C_V, C_P, K_S and K_T and their closure relation. After constant unit normalization, H>0 and a declared orientation give

\[
J_H=\frac{RH}{\sqrt{\det H}},\qquad J_H^2=-I,
\qquad J_H^TH=-HJ_H.
\]

This explicitly represents the already-derived native iota on the thermodynamic response plane. Once the Hessian is declared the response metric, metric compatibility and zero torsion select its connection. In its natural affine chart,

\[
A_i=\tfrac12H^{-1}\partial_iH,\qquad
F_{ij}=-\tfrac14[H^{-1}\partial_iH,H^{-1}\partial_jH].
\]

The higher response tower thus gives a computable local transport geometry. The curved polynomial fixture in the proof has K=5/324 at its reference state; a separate nonconstant metric has an exact path transport preserving both the response norm and the represented iota.

These are thermodynamic state-space results under the declared metric choice. They do not identify a physical clock, spacetime gravity, photon sector, or numerical alpha. Standard geometric identities are consumed explicitly; worldwide novelty and independent empirical validation are not claimed.

## Reproduce

Python 3.11 or 3.12; standard library only. From the repository root:

```bash
python papers/thermo-compass-foundations/verify.py --check
```

This performs exact controls for 388 positive Hessians, 54 rational normalized phase planes, 243 symmetric cubic-response fixtures, 729 comparisons with differentiated full Christoffels, four nonconstant transport endpoints, 81 multiplication and nine dagger intertwiners against the pinned native scalar implementation, seven input refusals and five negative controls. It also runs the recovered helical witness and the earlier conditional coupling controls. Finite counts state the computation's scope; the general arguments are in THEOREM.md.

Default mode is also read-only. Intentional evidence regeneration after reviewed source changes is:

```bash
python papers/thermo-compass-foundations/verify.py --write
python papers/thermo-compass-foundations/verify.py --check
```

The new workflow runs this entry point on Python 3.11 and 3.12. A local PASS is separate from the hosted workflow's status.

## What was recovered

The previous helical source calculation was present as an unpublished local draft. Its proof, source pins, exact implementation and native scalar dependency now live under [Physics from Cuts](../uncut-cut-measurement/NATIVE_IOTA_HELICAL_RESPONSE.md). Its earlier “v1.8” heading is retained as draft provenance; it is separate from the previously published Quantum–EMK v1.8 note. The exact sheet-local value 2/3, sheet-summed value 6/7 and difference 4/21 are preserved.

The [fine-structure R1 working note](downstream/FINE_STRUCTURE_R1.md) and [its exact controls](downstream/fine_structure_controls.py) are retained downstream. It is a conditional coupling reduction and does not derive numerical alpha. It is not the starting assumption of this thermo-first programme.

The earlier Physics-from-Cuts master certificate remains byte-for-byte unchanged; the new certificate separately binds and checks the recovered helical files and the thermo development.

## Next construction

Develop the response-fibre connection and holonomic jet carrier from this tangent/cotangent geometry, declare a native process law, and compute its transport/raiser compatibility defect. Reversible transport, irreversible response and retained memory must keep their distinct meanings. Physical propagation, matter identification and constants follow only through further proved interfaces.
