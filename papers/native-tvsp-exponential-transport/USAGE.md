# Following the observer through the source layers

The displayed number is a reading. Its TVSP context supplies the geometry
and response derivatives needed to compare that reading with another point.
The new adapter constructs the transport coefficients from that retained
source; the movement path remains an explicit input.

## Two observer questions

| Question | Comparison | Result |
|---|---|---|
| Does a transported vector keep its reading for a co-moving observer? | O(t)=O_0 U(t)^(-1), w(t)=U(t)w_0 | O(t)w(t)=O_0w_0 |
| Did a loop rotate the frame relative to its starting reference? | Fixed base O_0 reads the returned U_box w_0 | Leading signal is -area times O_0 F_sv w_0 |
| Does the actual scalar response remain constant along a chosen path? | Read f(gamma(t)) from the source | Depends on the source and path; transport alone does not force it |

For the original TLM family, the s-axis keeps f=-lambdaS=5, while the
ambient curvature changes for nonzero mu. The v-axis gives f=5(1+t).
Both have exact closed transport U=(I+t H_0^(-1)D)^(-1/2).
On the curved route (s,v)=(t,t²), the coordinate-frame transport needs
ordering: replacing it with the ordinary exponential of the integrated
generator changes its cubic coefficient by F_sv(0)/6.

In the oriented orthonormal frame, NT's rank-two connection has one
rotation generator, so its exponential can be written as an ordinary
phase exponential. That is an explicit frame conversion; it is not
permission to omit ordering in the original coordinates.

## Application API

From this folder:

```javascript
const L = require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const T = require('./transport_adapter.cjs');
const N = require('../native-common-cone-limit/native_engine/native_operator.cjs');

const source = L.buildFixture('1', 4);
const path = T.sourcePath(source, {a: 1, b: 1, p: 1, q: 2});
// s(t)=a*t^p, v(t)=b*t^q. Coefficients are exact native Cut matrices.
const U3 = T.coefficient(path.U, 3);
const observer = T.observer(path.U, N.matrix([[1, 0]]));
const loop = T.smallLoop(source);
const minusCurvature = T.coefficient(loop.holonomy, 2);
```

`path.H` and `path.U` retain depth 4; `path.A` retains depth 3.
`coefficient` rejects requests beyond the retained depth. `project` can
only lower depth. `derivative` consumes one degree. `transport` determines
one additional coefficient by the specified ODE and its initial value;
it does not invent a higher source derivative.

`restrictScalar` and `restrictMatrix` support zero-origin monomial paths
s=a*t^p, v=b*t^q with positive integer p,q and real rational a,b.
The written theorem allows arbitrary supplied smooth paths; this finite
adapter deliberately implements the stated path class. Contexts record
the source and path identity, and different contexts or depths cannot be
silently multiplied. The full constitutive source and path must remain
with serialized coefficient records; `serialize` alone is not that context.

`connectionFromMetric` checks algebraic Hermitian symmetry and invertibility
of the constant term. Positivity on a path, Hessian integrability and the
physical source interpretation remain caller hypotheses. `sourcePath` uses
the supplied TLM source with the stated stable-patch restriction.

`affineBinomialTransport` computes the finite coefficients of the affine
metric formula. Zero higher retained coefficients alone do not establish
that an unspecified full source is affine. In the original TLM straight
paths that fact follows from the full supplied potential.

`ordinaryExponential` accepts a square exponent with zero constant term;
its truncated series is then finite. It does not evaluate arbitrary
constant matrix exponentials. `smallLoop` constructs the order-two return
for equal side lengths epsilon, including actual connection derivatives;
it is not an exact finite-loop integrator.

## Exact endpoint and finite accuracy

For mu=5 on the v-axis, t=369/10000 is inside the source patch. With
Q=[[1,0],[-1/5,0]], the full analytic formula gives exactly U=I-Q/21.
Its transported metric and cut can therefore be checked with rational
arithmetic, without estimating a transcendental or truncation error.

For a general supplied source, use the residual bound in TET-T6 when
interval bounds are available. Equal finite jets alone do not guarantee
equal endpoint transports, even for positive analytic TVSP Hessians.
The finite regular matrix of a jet is an information encoding; its formal
adjoint is not an ordinary Hilbert-space adjoint or a physical unitary.
