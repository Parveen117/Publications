# Following the density and coherent-cut controls

The new generator is B=[W,D], where D=dot(theta).X uses the original
density channels and W=Pi+Q/2 assigns different phases to the base state,
its response plane and the remaining Hilbert space. The weights are
fixed by the required response motion. A base-state-only check cannot
determine them: Pi/2 moves that state correctly while missing the higher
response block.

| Operation | What is supplied | What it establishes |
|---|---|---|
| Density phase exp(i delta_p D) | Existing strongly commuting density operators | A spectral unitary on the selected quantum sector |
| Cut phase exp(-i delta_p W) | Declared coherent phases of Pi and Q | Relative phase control between the three subspaces |
| Four-pulse cycle | Rightmost-first product in DCD-T4 | Generator B at order delta_p² with a bounded remainder |
| Refined cycles along the path | Specified controls and uniform norm bounds | Operator-norm convergence to QRT's K |
| Smooth bounded-mass preparation | Smooth spectral bumps and finite total-mass cap | Nonzero curvature and preservation of original current domains |

Measuring Pi or Q does not automatically implement its coherent phase.
The pulse construction assumes those controls and does not invoke
measurement collapse or postselection. The old density/current flows
are proved insufficient for the stated small-loop base/response action.
This identifies the required control extension rather than silently
assigning it to the original material current.

## Exact finite comparison

From this packet's directory:

```javascript
const QRT = require('../native-quantum-response-transport/response_adapter.cjs');
const T = require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const DCD = require('./density_cut_adapter.cjs');

const source = QRT.finiteSource(
  ['1/2', '1/2', '1/2', '1/2'],
  [[1, 1], [2, 1], [1, 2], [3, 4]],
  3,
  'four-atom response-leakage comparison'
);
const path = DCD.densityCutPath(source, {a: 1, b: 1, p: 1, q: 2});
const pulse = DCD.pulseJet(
  T.coefficient(path.W, 0),
  T.coefficient(path.D, 0),
  3
);
const secondCoefficient = T.coefficient(pulse.product, 2);
const wrongGenerator = DCD.weightedGenerator(path, '1/2', 0);
const budget = DCD.pulseResourceBound({delta: '1/1000', w: 1, d: 4, steps: 1000000});
```

`densityCutPath` uses QRT's state and response geometry with TET's path
calculus. For source depth N, B retains depth N-1 and its transport K
retains depth N. `Bprojection` and `Kprojection` retain the independently
constructed QRT expressions for exact comparison. Different contexts and
depths must be projected explicitly, as in the inherited adapters.

`pulseJet` creates a separate formal pulse-parameter context. Its
chronological order at orientation +1 is -iD,+iW,+iD,-iW. Thus the
displayed product is exp(-i delta_p W) exp(i delta_p D)
exp(i delta_p W) exp(-i delta_p D). Its first coefficient vanishes,
its second is B, and its cubic coefficient is (i/2)[D-W,B]. Reversing
the orientation gives -B at second order. These exact coefficients
are not an evaluation of the exponential at finite delta_p.

The three-atom saturated comparison has no remaining response direction
and hides the incorrect Pi/2 weight. The four- and five-atom comparisons
retain a nonzero higher response component and reject it. They are
finite algebraic checks; continuum injectivity uses the actual joint
spectral absolute-continuity theorem.

## Reading the error budget

`pulseResourceBound` computes the rational inequality

    local_error <= delta_p³(w+d)b + delta_p^4 b²/2,
    repeated_error <= steps * local_error,
    effective_duration = steps * delta_p².

The caller supplies valid upper bounds w,d,b for the corresponding
operator norms; b defaults to 2wd. The adapter cannot infer such bounds
for an unspecified continuum source. In the example above the effective
duration is one and the bound is 1251/31250, below 1/20. This is a
conservative mathematical budget, not a laboratory gate estimate.
The moving-path error additionally includes the freezing term in DCD-T4.

The analytic theorem uses bounded D_B on the preparation support. The
same four-pulse product can use the original unbounded D globally,
because its two phases cancel outside that support. Norm convergence
alone is not used to infer current-domain preservation: DCD-T3 proves
that independently using graph norms.

## Scope and replay

All six theorems retain the selected LMD sector and its actual density
law. The material-metric result is a nonidentifiability statement:
density geometry misses conjugate phase information needed by the
original inverse chart. No pointwise quantum metric is manufactured
by the pulse adapter.

From the repository root:

```bash
python3.12 papers/native-density-cut-dynamics/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The packet contains 143 exact checks, 13 symbolic certificates and
23 negative controls. Python 3.12 handles I/O and hashes; all scientific
arithmetic and symbolic certificate replay use the frozen native engine.
