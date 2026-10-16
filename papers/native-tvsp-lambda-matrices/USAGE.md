# Reading a lambda-layered entry

Think of the outer number as a displayed reading and its context as the local
TVSP compass from which that reading came. Opening a layer reveals response
derivatives generated from that compass, not arbitrary added information.

For the exact source family in TLM-T5:

    visible entry:         -lambdaS(0,0)=5
    deformation control:   mu
    retained source layer: partial_v lambdaV(0,0)=-mu
    curvature readout:     kappa(0,0)=5mu²/324

At mu=0 and mu=5 the outer entry is identical. The first retained context
layer distinguishes the two. At mu=5 and mu=-5 the curvature agrees, but the
oriented response layer still distinguishes them.

The displayed first-layer source matrix is

\[
\operatorname{Reg}_1(H_\mu)=
\begin{pmatrix}
2&1&0&0&0&0\\
1&5&0&0&0&0\\
0&\mu&2&1&0&0\\
\mu&0&1&5&0&0\\
\mu&0&0&0&2&1\\
0&0&0&0&1&5
\end{pmatrix}.
\]

The lower blocks are the two directional first layers. The usual point
matrix occupies every diagonal block. Reading only eigenvalues forgets
the new layers; reading the constant-input columns recovers them exactly.
The matrix is an encoding of multiplication on local records. Its ordinary
transpose is not the encoded formal adjoint, and it is not presented as a
six-dimensional physical space or one self-adjoint quantum measurement.

Use the adapter only after choosing a context and depth. Independent contexts
or depths are rejected rather than silently combined. Project explicitly to
a common lower depth. Differentiation consumes one depth; an unknown next
layer is never filled in with zero and called a proved derivative.

The six-response source remains visible in the record. Its covariant tower
uses the actual connection formula of THEOREM.md; scalar partial derivatives
are storage coordinates for the tower and are not silently substituted for
every covariant tensor slot.

For example, from this folder:

```javascript
const L = require('./tvsp_adapter.cjs');
const source = L.buildFixture('5', 3);
const reading = L.outer(source.r);                 // 5
const firstLayer = L.regularMatrix(L.matrixProject(source.H, 1));
const tower = L.responseTower(source, 3);           // L1 through L4
const orientedResponse = L.outer(tower[1].components.lambdaV['1']); // -5
```

In tower component keys, 0 is the s direction and 1 is the v direction;
the leftmost index is the newest covariant derivative. Layer k stores scalar
jets at remaining depth N-k. The six-response fibre connection is flat;
the lower tangent slots use the actual Hessian connection. An arbitrary
caller-supplied matrix still needs the source positivity/integrability
hypotheses; an algebraic shape check alone does not certify a new material.

`serialize` and `serializeMatrix` export coefficient-only snapshots. A full
contextual record must additionally retain its source specification, H/Lambda
context layers, calibration and frame; those functions do not claim to infer
or serialize the missing source data. `buildFixture` returns that source and
the actual response records together for the demonstrated family.
