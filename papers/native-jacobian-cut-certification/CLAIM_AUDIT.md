# NJC-1 source and claim audit

Eight supplied files contain five TeX manuscripts and three Git reflog entries.
The two circular–hyperbolic entries point to the **same commit**, so they are
one mathematical source, not two independent certificates.

| Supplied item | Actual content / provenance | Decision |
|---|---|---|
| `sine correction jacobian.tex` | Full manuscript | Exact reflection/sine identity accepted; lag, phase-product and physical identifications require correction |
| `emk_tensor_jacobian.tex` | Full manuscript | Involution core accepted; seam, tensors and connection rebuilt with existing CID-D1; proposed Einstein derivation does not pass |
| `JACOBIAN SERIES.tex` | Full synthesis manuscript | Reflection and typed recurrence retained; convergence equivalences and Millennium conclusions not certified |
| `cut-graded-lambda-jacobian-tower` | Reflog pointer to RKF commit `9f5792ee62ce3a1a71d9ed242ff9cb10e6745f3e` | Source recovered; ten legacy tests pass; native port passes after making hypotheses explicit |
| Both `circular-hyperbolic-jacobian-semigroup-bridge` files | Reflog pointers to MP commit `b4f14ef4e41c6e2a18918d7322fd1b4968230a90` | Full theorem, implementation and tests recovered; algebraic core ported to native engine |
| `extended exponential.tex` | Full late-added manuscript | Scalar identities retained with normalization/completion hypotheses; foundational priority is an interpretation |
| `expo jacobian series.tex` | Full late-added manuscript | Real-linear exponential accepted; relaxation, imaginary-unit typing and memory derivative corrected; thermodynamic/topological claims not derived |

The original attachments are unchanged. Original byte hashes, recovered Git
blob hashes and predecessor pins are in `SOURCE_PINS.json`.

## Sine correction manuscript

1. `J(z)-z=-2 iota r exp(iota*pi/4) sin(delta)` is correct for all `delta`.
   The fixed set is the real line through `1+iota`, including zero and both
   directions. `Jz=z` does not impose `|z|=1`.
2. For a differentiable curve `f(t)=J(z(t))`, the limit of the centered
   difference equals `f'(t)`. The draft's limiting lag is therefore zero.
   It is not the generally nonzero reflection residual. Finite-resolution
   error could be defined separately, but would require a resolution and
   regularity/error theorem.
3. A product of complex phase factors contains both cosine and sine cross
   terms. It is not the product of their sine components: phases `1,iota`
   give product `iota`, while the sine-component product is zero.
4. The zeta multiplier has root `q=2^(-s)=(1-iota)/2`. Under the usual
   exponential interpretation its roots are
   `s=1/2+iota*(pi/4+2*pi*k)/log(2)`, `k` integral. The formula in the
   synthesis draft for the imaginary part is incorrect. Multiplication
   retains every existing zero of zeta in the critical strip; it cannot
   exclude off-line zeta zeros. No RH conclusion is promoted.
5. Historical priority concerning the Surya Siddhanta, the measurement
   interpretation, and universal physical correction claims are outside
   this algebraic certificate.

## Tensor manuscript

| Claim | Finding and repair |
|---|---|
| Opposite eigenspaces identified along a nontrivial seam | Their intersection is `{0}`. Use a projector seam or two copies with a comparison map (T2). |
| Proposed parallel/perpendicular pair formulas | Components do not have the declared eigenspace type; their reconstruction and distance factors fail explicit two-component examples. |
| Eigenvalues of anti-linear `J` | The `+/-` decomposition is over the real module. It is not the usual complex-linear spectral theorem. |
| Arbitrary `Phi(e)=J(e)exp(i theta(e))` maps `E` to `M` | A general angle does not map one specified eigenline into the other. A comparison map must be supplied and checked. |
| Displayed tensor transformation law | Measurer/measured labels do not make tensor slots upper/lower indices. Use actual vector and dual actions. The extra phase factor is not derived and can double-count transformations. |
| Cross-fibre pairing `omega(v)` | Needs an identification between the fibres. A scalar Gaussian weight does not supply it. |
| `omega_J=tr(J nabla J)/2` generates nonzero transport | It vanishes for a differentiable involution and a connection compatible with composition. Retain the operator-valued term and reuse CID-D1 (T3). |
| Sum of two connections is automatically a connection | Types and Leibniz coefficient must be checked. On one bundle, summing two ordinary connections doubles the derivative of a scalar. The draft supplies neither a corrected product-bundle definition nor the claimed coefficients' derivation. |
| Mixed Riemann term from scalar gradient commutator | `(X kappa)(Y kappa)-(Y kappa)(X kappa)=0` for scalar kappa. It cannot yield the proposed nonzero mixed Ricci tensor. |
| Ricci and scalar formulas agree | The trace of `(d kappa tensor d kappa-g||d kappa||^2)/2` in dimension n is `(1-n)||d kappa||^2/2`, not `+||d kappa||^2/2`. |
| All Riemann symmetries/Bianchi contractions | Need an appropriate tangent connection, metric compatibility and torsion hypotheses. A general bundle connection has a covariant Bianchi identity, not automatically all Levi-Civita curvature symmetries. |
| Einstein equation derived from the split | No action or source-identification theorem is supplied. The proposed source cannot be inferred merely by naming part of a curvature expression stress-energy. |
| Ordinary seam delta gives a nonzero `kappa^2 delta_seam` source | If the delta is a measure supported where smooth kappa vanishes, this product is zero. Derivatives of delta would be a different distributional contract. |
| `kappa=0` removes derivative corrections | Value zero does not in general imply derivative zero. A norm-distance may instead be nonsmooth on its zero set; a regularity statement is needed. |
| One-dimensional ordinary real example `J(x)=x` | Original `J(x)=iota x`, so this changes the realization without stating it. |
| Rotational-surface example with arbitrary f(r) | Curvature cannot be independent of the unspecified f. For f(r)=4r, u=2 sqrt(r), w=sin(2 phi)/2 gives `du^2+u^2 dw^2`, locally flat where cos(2 phi) is nonzero, contradicting the claimed nonzero component. At cos(2 phi)=0 the proposed metric is degenerate. |

Thus the manuscript as a whole is **not certified**. T1–T3 and T8 give a
consistent, useful replacement for its algebraic and connection core.

## Jacobian synthesis manuscript

* `J^2=I` and the real fixed-line description are valid. Bare powers alternate
  and do not yield arbitrary accumulating phases.
* The source does not fully define `J_alpha`. Even for the natural completion
  `P_++alpha P_-`, transverse decay alone is insufficient for convergence of
  a forced series. T5 supplies the exact recurrence and a divergent on-seam
  counterexample.
* `lim (1+i*pi/n)^n=-1`, not Euler's number e, under the usual completed
  exponential. A nonterminating sequence alone also does not prove
  transcendence; rational sequences can be nonterminating too.
* RH, Yang–Mills, Navier–Stokes, BSD, Hodge, Collatz and force-identification
  statements do not come with the requisite concrete operators, faithful
  representations and implication proofs. They are not certified theorems.
  No claim about the truth or falsity of the open problems follows from
  rejecting an unsupported derivation.
* A fixed involution on a gauge-labelled space does not define a physical
  Hamiltonian gap. A contraction parameter becomes an energy only after a
  justified transfer/Hamiltonian and scale relation. `Delta=1-alpha` is not
  such a derivation.
* Fixedness on cohomology or a statement about a fundamental group cannot
  replace the missing algebraicity or topological hypotheses. For example,
  simply connectedness alone without closedness does not characterize S^3.
* Perception, universal constants and universal-object assertions remain
  interpretations or proposals; the certificate does not turn them into
  mathematical consequences of a two-dimensional reflection.

## Recovered tower and bridge

The recovered tower is substantially more precise than the synthesis draft.
Its finite diagonal-cut examples pass. The amended use contract requires:

* a tangent connection as well as a response-bundle connection;
* affine constant base cut for the displayed ordinary higher-derivative
  formula, or properly covariant cut-compatible differentiation;
* `J^T=J` for the stated cut/transpose joint grading;
* blind target images `Pi(ker A_N)`, rather than kernel intersection with an
  unspecified target support;
* explicit topology/domain hypotheses for an analytic tower exponential.

The MP bridge explicitly limits itself to bounded quadratic generators and
does not claim a universal quarter threshold or a Millennium solution. Its
old executable uses SymPy. This packet rederives its operative algebra with
the unchanged native engine. The original MP test suite was inspected, not
run or relabelled as a native test. Circular evolution and algebraic sector
projectors are not automatically positive-unitary or orthogonal.

## Extended exponential manuscript

The identity `sqrt(2) exp(iota*pi/4)=1+iota` is valid in the completed angular
representation, and the native algebra verifies `(1+iota)^2=2 iota` and
`|1+iota|^2=2`. Fixedness alone does not fix radius or select a primordial
length. Repeated application of J has an orbit of at most two points and
does not by itself generate the full Gaussian-integer lattice.

The Gaussian formula is a valid rescaling of a supplied Gaussian integral:
use `x=sqrt(2)u` directly in `integral exp(-x^2/2) dx`. It is an analytic
comparison input, not a new cut-only derivation of its integration measure.
The staircase issue is the failure of path-length continuity under uniform
convergence of these approximations, not inability of a grid to express a
diagonal endpoint.

On the full line `s=1/2+it`,
`2^s pi^(s-1)=sqrt(2/pi) exp(it log(2*pi))`. The draft's real value is its
value at `t=0`; its magnitude is constant along the line. Rewriting this
factor does not restrict the zeros. This source correctly acknowledges
that it does not prove RH.

`D_v N_lambda=lambda N_lambda` implies the advertised exponential evolution
under a defined domain/completion; the eigenvector equation is a hypothesis
to establish. Fractional dimension flow and preservation of an unnamed
invariant are not implied by scalar exponential evolution.

## Exponential–Jacobian thermodynamics manuscript

| Claim | Decision |
|---|---|
| `exp(tJ)=cosh(t)I+sinh(t)J` for real t | Accepted, T4/T5. |
| Distance to fixed line decays | Accepted; distance equals `exp(-t)||P_-z||`. |
| Full flow is a contraction and noninvertible | Rejected: plus-sector growth and explicit inverse. Use T5 for bounded relaxation. |
| `exp(i tau J)(z)=cos(tau)z-sin(tau)conjugate(z)` is unitary | Rejected for the original anti-linear J. Internal quarter-turn anticommutes with J; scalar-complexified J gives a different, valid formula (T5). |
| `sigma=R exp(-lambda t)` and its integral | Mathematically valid assumed positive rate for positive R, lambda. Its identification as thermodynamic entropy production is an axiom in the source, not derived from J. |
| Time-dependent cost exponential derivative | Incorrect; T6 supplies both the correct derivative and actual initial-value propagator. |
| Memory limit implies Markovianity/reversibility | A scalar accumulated cost does not specify a stochastic process, memory kernel or invertible physical dynamics. |
| `Omega=dJ` is a curvature two-form | Type error for an endomorphism-valued zero-form J: dJ is a one-form. Use a declared connection and its curvature. |
| Onsager reciprocity iff that curvature is zero | No constitutive bridge is supplied. Existing EMK-G1/LTB1 already record independence of response symmetry and independently chosen connection curvature. |
| Curvature integral over a transition curve is an integer jump | A two-form does not integrate over a one-dimensional curve as written. Integrality also needs a specified global bundle/cycle/normalization. A curvature-zero locus need not be a topological transition. |
| Claimed chromium/copper jumps | No data or calculation provided; not tested or certified. |
| `R/lambda=i*pi(2k+1)` for positive real R, lambda | Incompatible real/imaginary types. It does not follow from existence of a limiting retained state. |
| Maximum entropy and fluctuation–dissipation theorem | No statistical state space, entropy functional, noise covariance or response theorem supplied. Retained as proposals only. |

## Sources consumed without duplicating their claims

* RKF T42 at the supplied commit: `theorum/42_cut_graded_lambda_jacobian_tower_theorem.md`
  and its original proof-lab source/test/pin files. The local copies' Git
  blob hashes match that exact remote commit.
* MP at the supplied commit:
  `programmes/circular_hyperbolic_jacobian_semigroup_bridge/CIRCULAR_HYPERBOLIC_JACOBIAN_SEMIGROUP_BRIDGE_THEOREM.md`
  and its source/test files.
* Publications CID-D1 D2e, EMK-C1 and EMK-G1/LTB1, cited by repository paths.
* Standard zeta definitions/comparison:
  [NIST DLMF 25.2](https://dlmf.nist.gov/25.2) and
  [NIST DLMF 25.10](https://dlmf.nist.gov/25.10). No external theorem proves
  any new arithmetic or physical claim in this packet.
