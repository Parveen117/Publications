# Native cut–transport feedback — CID-D2

The cut and its transport generator now evolve together:

\[
D_\lambda P=[V,P],\qquad
D_\lambda V=[P,H_I(P)+H_C(P)].
\]

The selected positive information/curvature functional is explicit, and
its first variation derives both force terms. The conserved total is
transport information plus directional information plus squared induced
curvature. [THEOREM.md](THEOREM.md) gives the complete derivation, action,
gradient, covariance, conserved quantity and scope.

An exact rank-one two-mode solution has a rotating cut and a rotating
generator with a fixed quarter-turn offset. It supplies a precise
two-strand helical interpretation in operation/record coordinates. Both
return after one cycle; their accumulated transport returns after **five**.
The coefficients selecting this example are model data, not predicted
physical constants.

## What is new relative to the predecessors

| Predecessor | Consumed result | New step here |
| --- | --- | --- |
| CID-D1 / QTH-1 | One tensor has symmetric information and antisymmetric induced-curvature parts | A selected action makes those quantities feed back into an evolving generator |
| QG-R2 | Retained records and fused return are different observations | Even the pair `(cut,current generator)` can return before the transport record |
| EMK-G3 | Visible return and lifted memory can differ | An explicit coupled native operator evolution produces a nontrivial return phase |
| GE1 / EMK-T2 | Native exponential and ordered transport | Exact solution of the varying generator, with its order retained |

## Reproduce

```bash
python3.12 papers/native-cut-transport-feedback/verify.py --check
```

Python **3.12 only** orchestrates. Node **24.19.0** runs the unchanged
pinned RKF native engine. The command checks every direct source pin,
replays CID-D1/QG-R2 and its full QG-R1/phase/six-packet chain, safely
extracts the same engine, and reproduces the new certificate byte for byte.
It changes no committed result file. To deliberately regenerate the new
certificate, use `--write` instead.

New local evidence: **39 algebraic certificates and serialized replays,
150 exact checks, 11 negative controls**, six exact orbit cases, 24
independent first-variation checks and six frame cases. The four native
presentations pass critical-pair and applicable dagger checks.

| Result | Evidence |
| --- | --- |
| Cut and horizontal-generator constraints; acceleration equation | Native noncommutative word certificates plus differentiated-constraint proof |
| Energy cancellation and variation identities | Native commutator certificates; their traces vanish by cyclicity |
| Explicit information and curvature gradient | Written arbitrary-carrier variation, independent exact rational variations |
| Two-mode information/curvature force and orbit | Native symbolic sphere quotient, with checked confluence and scalar dagger |
| Five-cycle full transport return | Native fifth-root quotient, unit-transport certificate and four nonreturn controls |
| Variational/ordered-exponential and helical statements | Written analytic proofs; not represented as numerical simulations |

`CERTIFICATE.json`, `EXPECTED.sha256` and `SOURCE_PINS.json` hold the
deterministic evidence and its lineage. Historical EMK-G3 is audited and
pinned as a predecessor; its Python arithmetic is not used as a new engine.

## Scope

This is a closed, reversible **selected finite model**. The action, control
operators, weights, operation parameter and initial data are declared.
It does not show that every solution is helical, that any single
information statistic is preserved, or that cut return erases transport.
It supplies no universal gravitational action, physical clock, continuum
quantum-gravity limit, singularity result or measured constants.

The control curvature penalty is precisely the induced increment from
CID-D1. Ambient curvature is kept distinct. Lambda-dependent controls or
weights would add explicit source terms to the conserved-total equation.
