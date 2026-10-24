# EGR-1 — Exponential cut and gravity response

14 October 2026. A concrete realization of the proposed
invariant-source / exponential-transport bridge.

Read [THEOREM.md](THEOREM.md) for the equations, proofs, assumptions and
the distinction between cut alignment and spatial singularity.

## Results

* Nine actual four-simplices, ten interior squared edges, eight
  flat-background vertex-displacement null directions, and two curvature
  response coordinates.
* Exact reduced Regge Hessian `K=diag(10,47)`; the full Hessian factors as
  `H=A^T K A`, and every linearized deficit variation factors through A.
* For the selected invariant source coupling,
  `y=(alpha*N/gamma)*(r_1/10,r_2/47)`. At fixed source direction, response
  scales with N and its quadratic score scales with N squared.
* Selected cut relaxation preserves the source strength and cut spectrum
  while giving
  `m(lambda)=m(0)*exp[-(37/235)*(alpha^2*N^2/gamma)*lambda]`.
  Its curvature response score increases monotonically.
* A fully rational transport connects two points on this exact flow,
  with ratio-decay factor `5/17`.

The numbers 10, 47 and 37/235 depend on the specified geometry, source
normalization and relaxation rule. They are not universal constants.
The geometry is linearized about a flat Euclidean subdivision. The
coupling and first-order relaxation are explicit choices. No physical
clock, universal information-to-gravity identity, singular collapse,
nonlinear constraint closure or Lorentzian quantum measure is claimed.

## Reproduce

From the repository root:

```sh
python3.12 papers/native-exponential-gravity-response/run.py --check
```

Node.js runs the exact mathematics. Python only verifies the archive and
extracts the unchanged native engines into a temporary directory. The
source is the existing `Native_Phase_Word_Reconstruction_Source.zip`.

[CERTIFICATE.json](CERTIFICATE.json) stores the full geometry, response
matrices, word proofs and negative controls. `--check` validates the pins,
reconstructs the payload and compares it byte for byte. To regenerate
after deliberately updating source pins, use `--write`.

Exact-check counts include symbolic proof/replay checks. They must not be
added together as independent theorems. General analytic statements use
the written hypotheses and proofs in the theorem document.

The next geometric target is nonlinear continuation of these two curvature
coordinates with a specified source and retained angle lifts. That is where
loss of regularity, changing cell volumes and actual focusing can be tested.
