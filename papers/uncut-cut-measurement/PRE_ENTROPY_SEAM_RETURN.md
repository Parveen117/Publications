# Pre-entropy seam return: local flatness, global memory and the gravity interface

**Version 1.5; status: exact conditional theorem and countermodels, not an empirical gravity result.** No spacetime is needed for the finite results. The literal uncut ground remains a proposed primitive; every vector, loop and observable below belongs to a declared *post-cut representation*.

## Source audit and type ledger

The private thermodynamic seed, the public corrected [thermodynamic response edition](https://github.com/Parveen117/Publications/tree/c1a35e0f509cd137d20f9504543f97b5ec7479fc/papers/thermodynamic-response-corrections), the public [RKF cut grading](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/927cdb6ca98221c0b4285da953a2c8b689fb202c/theorum/41_cut_graded_universal_generator_theorem.md), [lifted memory](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/927cdb6ca98221c0b4285da953a2c8b689fb202c/theorum/recognition_topology/03_principal_holonomy_blindness_and_lifted_memory.md), and sectoral atomic implementation were compared on 2026-09-25. The private sources are not reproduced or assumed as public dependencies. The corrected edition's TVSP compass uses overlapping charts of a two-dimensional equilibrium state manifold; its spacetime Recognition quotient begins with a **supplied** spacetime local-equilibrium map, and its gravitational action and Lorentzian propagation require further declared input. The sectoral atomic implementation describes itself as a pedagogical, internal-unit model. Its energy-sector relations and coefficients are prescribed, not a calibrated correction of physical relativity.

| Object | Meaning and input | Logical role |
| --- | --- | --- |
| `I_th = Gamma_c Gamma_m` | Positive product of caloric and mechanical response ratios on the admitted stable thermodynamic equilibrium chart | `I_th=1` is its closure condition; `I_th=-1` is not an alternative stable equilibrium |
| `sigma=+1/-1` | Chosen sector of a cut-graded return operator `U_gamma`; `D_sigma=U_gamma-sigma 1` | Both sectors exist in the RKF cut-graded algebra; this sign is **not** `I_th` |
| `F=0` | Local curvature of a **given** transport on a chart or cell complex | Contractible-loop flatness; global return needs topology and the full connection |
| `k` | A separately retained integer lift of an abelian phase path | Terminal principal holonomy cannot recover `k` |
| `S` | Thermodynamic/statistical entropy on an additional calibrated state and dynamics | A scalar derived reading; neither an uncut primitive nor automatically a holonomy norm |

## SR-1. Target-relative return defect precedes an entropy reading

Choose a finite complex inner-product carrier `H`, a closed path `gamma`, a unitary continuation `U_gamma`, a sector `sigma in {+1,-1}`, an input `x`, and a **fixed target observation** `P:H -> Z`. Define the full sector return `D_{gamma,sigma}=U_gamma-sigma 1`, and the target-visible defect `R_{gamma,sigma,P}(x)=P D_{gamma,sigma}x`. For unitary `U_gamma`,

`D*D=2 1-sigma(U_gamma+U_gamma*) >= 0`, and `||D x||²=0` iff `U_gamma x=sigma x`.

**Proof.** Expand `(U*-sigma 1)(U-sigma 1)`, use `U*U=1` and `sigma²=1`; a squared norm vanishes precisely for a zero vector. For an eigenphase `e^{i theta}`, the eigenvalue is `|e^{i theta}-sigma|²`. Applying `P` loses a difference exactly when that difference lies in `ker P`. Thus the defect is explicitly target-relative. □

This is a *cut-side compatibility property* that can be assessed before assigning probabilities, temperature or entropy. Calling it a property of literal uncut ground would exceed the construction. Neither the operator norm of holonomy nor the endpoint alone is a faithful general memory measure: **every** unitary `U` has `log ||U||=0`, even when `U=-1` and `D_{+,gamma}=-2 1`. A nontrivial unitary loop need not produce positive `log ||U||`, and curvature alone does not establish monotone entropy growth. The pre-entropy record must retain the target, sector and enough path/lift data for the question being asked.

## SR-2. Flat plaquettes with a nontrivial global return

Consider the `3 x 3` periodic square grid, with vertices `(a,b) mod 3`, one-dimensional unitary edge carriers, vertical edge transports `+1`, and horizontal edge transports `-1` only on edges from `a=2` to `a=0` (otherwise `+1`). Opposite edges use inverse transport. Each elementary square has transport `+1`: the same horizontal seam factor occurs twice and cancels. Every local plaquette is flat. Yet the horizontal loop around the torus has return `U_gamma=-1`. Hence `D_{+,gamma}=-2`, `D*D=4`, whereas `D_{-,gamma}=0`. A flat local connection can therefore have a nontrivial global `-1` sector. Thermodynamic `I_th=1` can hold independently: no equation connects that positive ratio to this return sign.

**Proof.** The plaquette product is `h(a,b) v(a+1,b) h(a,b+1)^(-1) v(a,b)^(-1)=h(a,b)^2=1`; the horizontal winding product is `(1)(1)(-1)=-1`. The defect identities follow by arithmetic. □

On a simply connected chart with sufficiently regular flat connection and contractible loops, flat transport has trivial contractible holonomy. That qualified statement cannot erase noncontractible global memory. Nor is `U_gamma=1` sufficient to identify a *lifted* phase history: `U(t)=exp(2 pi i t)`, `0<=t<=1`, ends at `+1`, so `D_+=0`, while its retained winding is `k=1`. This latter blindness of principal holonomy is already established in RKF's lifted-memory theorem; its appearance here is a **dependency and countercontrol**, not a novelty claim.

## SR-3. What a spacetime/thermodynamic bridge must supply

An observation of `I_th=1`, vanishing local response curvature, or a scalar entropy cannot by itself determine the global transport sector, its winding, or target-hidden return. To connect this pre-entropy record to the corrected thermodynamic paper, specify a map from native cut/continuation data to its state-response jets, an observation/target pair, and compatible transport; the corrected edition constructs a target-relevant quotient *after those data are supplied*. To call its local recognition curvature physical gravity, additionally derive or calibrate a Lorentzian propagation metric, universal source coupling and action/dynamics, and compare a gravitational observable with data. No flatness-equilibrium equivalence beyond a stated faithful chart and adapter is proved here.

**Research question now sharpened:** can the native continuation law select a target-faithful return record and a physical propagation law from the same admissible cuts, such that thermodynamic closure and measured gravity become compatible projections? A finite answer to the return-record problem is given above; the physical selection remains open.

Run `python certificates/pre_entropy_seam_return.py` and `python certificate.py --check` for exact integer checks and a source-bound pin. The examples have no adjustable floating-point tolerance.
