# v1.6: a typed holonomy-to-geometry gate

**Status:** conditional mathematical bridge and exact nonselection theorem. This does not identify the v1.5 native return with a physical spin connection, derive a spacetime metric or derive gravity. The finite phase controls are exact; the smooth arguments are written proofs.

## What the earlier repositories already provide

[U20](GRAVITY_BEFORE_CURVATURE.md) proves that a **supplied** nondegenerate coframe `e`, a metric-compatible internal connection `A`, and the torsion equation `de+A wedge e=0` turn the connection into the Levi-Civita connection of `g=e^T H e`, with `R^g=e^(-1) F e` and `F=dA+A wedge A`. [RKF's lifted phase theorem](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/927cdb6ca98221c0b4285da953a2c8b689fb202c/theorum/recognition_topology/03_principal_holonomy_blindness_and_lifted_memory.md) already shows that an endpoint phase forgets integer turns. [v1.5](PRE_ENTROPY_SEAM_RETURN.md) constructs flat local plaquettes with a global `-1` return. These are dependencies, not discoveries of the present note. Cone angle and spin double-cover formulas are standard geometry; this note tests **which representation and extra inputs** let a native sign use them.

For the standard conical relation between vector holonomy and deficit, see the [original cosmic-string holonomy calculation](https://arxiv.org/abs/1211.4365). That physical realization is an external comparison, not an axiom for the native carrier.

## SH-1: one sign admits distinct geometric interpretations

Fix the orientation convention in which a vector holonomy around an oriented planar cone is `R(delta)` for deficit `delta=2 pi(1-alpha)`, with `ds²=dr²+alpha² r² dphi²`, `0<alpha<=1`, and `0<=delta<2 pi`. Its `Spin(2)` lift, once a spin structure and sign convention are selected, is `exp(i delta/2)`; switching the spin structure around the puncture multiplies this spin return by `-1`. In particular:

| Supplied identification of the v1.5 `-1` | Exact result | What follows |
| --- | --- | --- |
| Vector holonomy `-I_2` of a cone | `R(pi)=-I_2`, `alpha=1/2`, `delta=pi` within the chosen deficit range | A conditional conical realization; spin lift is `+i` or `-i`, not `-1` |
| Spin holonomy `-1` on a flat torus with an antiperiodic cycle | Vector holonomy `+I_2`; Levi-Civita curvature zero | A flat spin-structure realization with no cone angle |
| Spin holonomy `-1` with untwisted lift `exp(i delta/2)` | `delta=2 pi (mod 4 pi)` | No nondegenerate cone with `0<=delta<2 pi` realizes it under these hypotheses |

**Proof.** A frame rotates by the cone deficit, and the double covering `Spin(2)->SO(2)` maps the half-angle spin phase to the vector angle. Evaluating at `pi` gives vector `-I_2` and spin `i` (with inverse/orientation convention `-i`). The spin phase equals `-1` exactly when `delta/2=pi (mod 2 pi)`. The bounded cone deficit rules this out. For the flat torus, take its Euclidean metric and a spin bundle whose horizontal deck translation acts as `-1` on spinors while acting trivially on tangent vectors. Its lifted Levi-Civita transport around that horizontal cycle is `-1`, and every local curvature component vanishes. All three calculations are independent of entropy. □

Consequently a scalar `-1` does **not** select spin versus vector representation, a spin structure, an angle branch, or local curvature. A scalar line carrying sign holonomy can serve as a *candidate* flat spin local system; identifying it with a spin bundle requires a coframe, double-cover map and compatible path transport. The v1.5 `3 x 3` sign lattice supplies the flat holonomy pattern, but no such geometric map by itself.

## SH-2: how curvature can be recovered after the missing inputs are supplied

Given a smooth connection on a fixed chart and oriented *contractible* shrinking coordinate rectangles of side lengths `a,b`, the holonomy satisfies `U_rectangle=1+F_uv(p) ab+o(ab)` in a fixed representation (up to the declared orientation/sign convention). Therefore `lim_(a,b->0) (U_rectangle-1)/(ab)=F_uv(p)`. A **single** sign on a nonshrinking or noncontractible path supplies no such limit. If the coframe and torsion/metric compatibility hypotheses of U20 hold, the recovered `F` becomes metric curvature `R^g=e^(-1) F e`. Without them it remains the curvature of the declared internal connection.

**Exact compatible example.** On a post-cut two-dimensional Riemannian chart let `e¹=du`, `e²=f(u) dv`, `H=I_2`, `J=[[0,-1],[1,0]]`, and `A=f'(u)J dv`. Then `de+A wedge e=0`, `g=du²+f(u)²dv²`, and `F=f''(u) J du wedge dv`. For `f=1+u²`, curvature is `2J du wedge dv`; for `f=1`, it vanishes. Both coframe profiles are **selected as inputs**. This gives a full mathematical holonomy → internal curvature → metric curvature route, with each adapter visible; it does not infer `f` or a gravitational source equation from a native `-1`.

## SH-3: the action-scale gate

All the signs, winding numbers and angles above are dimensionless. Relabeling an action unit multiplies the numerical value of physical `hbar` without altering any of those return records. Thus a sign-only/angle-only return structure cannot select a numerical dimensionful `hbar`; an independent action calibration or a primitive with action dimension and a law relating it to the return would be required. This is a dimensional nonselection statement, not a no-go for future richer native primitives.

**Next empirical obligation:** specify a target-faithful map from native carriers to *spin and vector* probes, compare both returns on contractible loops at several areas and on noncontractible loops, independently calibrate the metric/coframe and physical sources, and look for a shared fit. The comparison between spin `-1` with vector `+I_2` versus spin `±i` with vector `-I_2` separates the two conditional examples before claiming gravity. No experiment is reported here.

The exact controls live in [certificates/spin_holonomy_geometry_gate.py](certificates/spin_holonomy_geometry_gate.py); `certificate.py --check` binds this note and those controls to the versioned evidence pin.
