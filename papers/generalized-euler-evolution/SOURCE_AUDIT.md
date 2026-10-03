# GE1 source audit: Generalized Euler and exponential change drafts

Both uploaded LaTeX drafts were read completely before choosing this work.
They are recorded as **Draft A** (Generalized Euler flow-space proposal) and
**Draft B** (exponential change on the torus). Their byte hashes are in
`certificates/GE1_SOURCE_PINS.json`; originals and personal filenames are
not republished. A source hash identifies evidence, not mathematical truth.

## What was already present

| Existing source | Relevant result | Consequence for this development |
| --- | --- | --- |
| RKF F00-E/F00-G | Native factorial exponential, radial completion and logarithmic chart | Do not introduce an unexplained new scalar exponential |
| RKF T55 | Generalized Euler/EMK phase–seam exchange, mixed bracket, determinant-channel exchange, lawful-residue test | The algebraic Euler/EMK connection is already present |
| RKF weighted operator completion WC1–WC4 | Weighted normal-word completion and bounded-element functional calculus; unbounded domains explicitly open | Does not license a Taylor series for an arbitrary unbounded derivative |
| Extra Ideas R42 | Finite native phase generator, reversal cut and controlled factorial interpolation | Finite phase reconstruction is not a new GE1 result |
| Existing Generalized Euler phase–ratio text | Dimensionless phase lift, declared ratios, winding and branch memory; typed curvature and projection boundary | A pure torus sector does not certify the whole lift |
| RH paper's MP source-of-record, native EMK Euler-scale carrier and audit | Multiplicative scale line, inversion cut, dilation flow, weighted energy completion, logarithmic/Mellin chart and spectral-blindness boundary | Infinite Euler flow and its native/chart distinction already exist; GE1 does not rediscover them or identify this carrier with a torus |
| RH-Framework T01-E4 and current native ledger | Native simple-function completion; actual UGD adapters and further targets separately gated | GE1 does not close the general RH carrier or RH itself |
| RH-Framework T03 action-domain adapter | Bounded multiplier/convolution and maximal multiplication classes, with actual source membership still a separate task | Mentioning an Euler action does not prove it belongs to a required domain |
| YM50–52 | Counted compact reference, symmetric turn heat, protocol tensor, energy and exact compact relaxation constants | A drift-to-heat claim must retain the protocol; known compact gaps are reused |
| YM55 | Fixed-profile anisotropic joint local-history limit | No change to its step law, window, evidence or remaining four-dimensional gates |

This audit covers these actual sources and both drafts; it is not a claim
to have exhaustively certified every file in every related repository.

The Publications RH landing page points to MP PR 244. Its source commit
`686b5506c7c2fedd53641dbf550218b60c996a13` was inspected, including the
native Euler-scale carrier, its audit addendum and the Euler–virial
spectral-zeta proposal. This is a source-of-record comparison, not a fresh
certification of the RH endpoint or of every status in historical source
files. The separate RH-Framework snapshot retains its own domain ledger.

There is a consequential norm distinction in that source. The dilation
`U_t F(r)=F(e^-t r)` is unitary in the unweighted dilation pairing. Its
weighted spatial energy uses `w(r)=max(r^3,r^-3)`. Substitution gives
`w(e^t r)<=e^(3|t|)w(r)`; the Mellin part of the stated energy is phase
invariant. Consequently `||U_t F||_(times,3)<=e^(3|t|/2)||F||_(times,3)`.
It is not generally an isometry in this weighted energy: for `t>0` and
nonzero core data supported in `r>1`, the spatial term is multiplied by
`e^(3t)` while the Mellin term is unchanged. GE1's isometric phase theorem
must therefore not be applied to that weighted carrier without a separate
growth/domain argument. This comparison uses the source's displayed energy;
it does not import its classical comparison results as premises for GE1.

## Draft A: claim-by-claim disposition

| Draft claim or construction | Assessment | Certified treatment or remaining requirement |
| --- | --- | --- |
| Exponential flow equals pullback by a flow | Valid for a complete real flow on an admitted carrier; Taylor representation needs an additional domain | GE1-T1 gives the constant phase sector and GE1-T2 its boundary |
| Orbit of one state equals the union of images of the entire state space | False in general; for a complete flow the latter union is the original space | No orbit quotient or enlarged state space is inferred |
| Countably many phase derivatives define an operator on all observables | Missing topology and domain | Finite-support character core and an explicit completed domain in T1 |
| Constant phase and additive ratio drift generate non-Abelian curvature | Constant coordinate drifts commute; connection curvature needs more data | Keep the certified EMK bracket and EMK-C1 connection separate |
| Listed invariants are eigenmodes of the displayed drift | Several fail: an angle differentiates to a constant; additive ratio drift sends `log Gamma` to `gamma/Gamma`; a constant seam phase has derivative zero | Do not promote the table; actual character eigenvalues are derived in T1 |
| Reciprocal-ratio product is invariant and a Casimir | Invariance requires a tangent evolution preserving that constraint; a Casimir additionally requires the relevant bracket structure | Neither structure is supplied by the displayed additive ratio drift |
| A flow commutator is curvature | Leading generator bracket can be checked, but it need not equal a connection curvature | Existing T55 and EMK-C1 supply the correct typed comparison |
| A globally defined symplectic form and quantized holonomy follow | New conjugate coordinates, convergence, global charts, bundle data and integrality are missing; the exponential also loses a factor of `iota` in the displayed formulas | Not certified by GE1 |
| Measurement is projection to the zero mode | Needs a specified observer/limiting procedure and pairing | T6 derives a particular time-average projection, without identifying it with all measurement |
| Winding derivative gives probability and Born's rule | Positivity/normalization and the link to amplitude are unproved; equating density with amplitude square assumes the conclusion | No probability or Born-law derivation is claimed |
| Lax form, complete integrability and a Jacobian torus tower follow | Declared equations do not establish their compatibility, conserved family or completeness | Remain proposals |
| All flow eigenvalues are roots of unity | False for generic frequencies and durations | T1 retains the actual exponential multiplier |
| A minimum absolute spectral value changes sign | An absolute value is nonnegative; it cannot diagnose a sign change this way | Use a specified dissipative generator and a proved uniform bound |
| `ker D` equals a zero vector field | Type mismatch; constants are annihilated by nonzero derivations | Corrected in T6 |
| A constant exponential is a general Magnus integrator with seam jumps | Time-dependent ordering, jump maps and domains are additional data | Existing EMK-T2 only certifies its own declared sector; no unrestricted Magnus claim |

## Draft B: claim-by-claim disposition

| Draft claim or construction | Assessment | Certified treatment or remaining requirement |
| --- | --- | --- |
| Character multiplier of `I+D_omega` is `1+iota k.omega` | Correct on the character core | T2 states the actual maximal-domain kernel |
| Kernel condition `k.omega=iota` | Correct; real frequencies give a zero kernel | Retained with completion, rather than just an algebraic finite span |
| For `omega=iota`, the kernel contains `exp(-iota theta)` | Sign error | Correct mode is `exp(+iota theta)` |
| Complex coefficients are an ordinary real torus flow | Generally false | Complex exponential has the weighted maximal domain (3), unbounded on the bilateral completion when its imaginary part is nonzero |
| Taylor exponential acts on every smooth function | False | One fixed smooth vector lies in every power domain while its Taylor terms diverge, T2 |
| `(I+tD/n)^n` always gives the exponential limit | True on each fixed finite core; not an unrestricted growing-cutoff assertion | T2 supplies a fixed-vector counterexample; T3 proves a stable alternative |
| Kernel eigenvectors decay as `exp(-t)` | Correct on the admitted eigenmodes | Retained without promoting it to full-space diffusion |
| Replacing the first derivative by a Laplacian derives heat | A new generator has been selected | YM51 supplies an actual symmetric counted protocol; T4–T5 close its domain/resolvent dock |
| `I+Delta` is automatically screened Poisson | With the displayed negative Fourier symbol of Delta, this is a resonant Helmholtz-type sign; a positive screened operator uses `m^2-Delta` | No elliptic solvability statement is imported |
| A minimal-coupling square equals a scalar potential plus one directional derivative | Missing quadratic-potential, divergence and ordering terms in general | Not used as a derivation of gauge dynamics |
| Schrödinger or gravitational equations emerge simply by choosing those generators | The displayed Hamiltonian/source laws are inputs; dropping a term also requires dimensionally controlled approximation | Not certified as derived physics |

## What GE1 adds and what it does not

GE1 adds explicit unbounded domains, graph-core closure, two stable evolution
approximations, a derived phase time-average projection, and an exact
near-resonance obstruction to promoting finite gaps. Its native YM
application is the **existing free compact heat generator**, including
rank-two and rotated anisotropic protocols.

The proofs are written mathematical arguments. Exact certificates provide
finite controls and mutation checks; a PASS is neither an automated proof of
every infinite statement nor a claim that either draft is certified wholesale.
No new Born law, physical clock, universal curvature, RH proof, interacting
row closure or four-dimensional Yang–Mills mass gap is asserted.
