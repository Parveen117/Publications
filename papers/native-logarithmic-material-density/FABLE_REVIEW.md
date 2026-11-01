# Review of Fable's infinitesimal-to-finite bridge

Inspected commit:
[97bb7fe4ddccff1a5371689ecc2055b409606048](https://github.com/Parveen117/Publications/commit/97bb7fe4ddccff1a5371689ecc2055b409606048)
on claude/residual-gauge-invariant-dual.
Read BRIDGE.md and its native verifier against the frozen RGR/CGD
domains. Replayed bridge/run.py --check in a detached worktree:
12 exact checks, 5 negative controls, certificate SHA-256
9053cb7e491150fa200e868abd5e1e1ba90b3e2df50240b4004690c422d326aa.

## What is useful

The report separates the two outstanding notions correctly. On
Dbar=closure(D) in the full Psi topology, continuity and orbit
differentiability integrate the infinitesimal constraints to the finite
group action. It does not assume every circle diffeomorphism is one
autonomous flow. It uses generation by flows and the correct residual
circle-diffeomorphism product from CGD's central quotient.

Given its imported RGR and external ergodicity/uniqueness hypotheses,
the structural conclusion is

    K_inf = C ell direct_sum Ann(Dbar).

This explicitly locates any additional continuous solutions in the
annihilator of the untested quotient Psi/Dbar. The positive-form
classification is correspondingly a core-closure classification.
This makes the remaining completeness obligation concrete and avoids
duplicating our physical-observable track.

Its two-point counterexample also correctly shows why Hilbert density
of cylinders does not imply density for uniform-on-Hilbert-balls
seminorms. The declared revised contract Psi*=Dbar is explicit.
Neither that revision nor a full-Psi uniqueness theorem is silently
adopted in LMD. The report remains an unmerged research dependency
candidate, not a formal continuum certification.

## Corrections needed before integration

1. Section 5(e) and the dependency ledger call the A-approximation
   problem *equivalent* to density of D. The same paragraph says
   failure for A need not annihilate the larger Dbar. The equivalence
   has not been proved. Retain it only as a sufficient route unless
   closure(A)=closure(D) is established. A plausible repair is to prove
   every fixed U_g a has all Frechet derivatives bounded on each
   H^3 ball, then show theta_R U_g a belongs to A and converges to
   U_g a in every Psi seminorm. The cocycle and its derivative bounds,
   not a Hilbert-density assertion, must justify this step.

2. Section 5(d)'s C^{1,1} approximation does not establish the claimed
   absence of all first/second-order obstructions, nor force a hypothetical
   obstruction to third order. C^{1,1} is not C^2, and no theorem relating
   derivative order to the Psi-continuous annihilator is supplied.
   Keep the valid approximation result and remove that stronger
   classification unless a precise functional/domain theorem is proved.

3. Section 9 describes W11 on a space with an extra complement but
   claims only two scalar weights. The actual code instead tests a
   full two-dimensional core, with no untested complement. The theorem
   in Section 7 correctly restricts to Dbar, so align W11's prose/code
   with that restriction. A positive form supported entirely on an
   untested complement annihilates all core constraint images and
   is an immediate counterexample to the unrestricted finite claim.
   LMD's native audit fixture checks this explicit complement form.

These are scoped review findings. They do not invalidate the
infinitesimal integration argument merely because full Psi density
is unresolved. They do prevent treating every advertised equivalence
or finite-model description as already established.

## Next prompt for Fable

Continue only on claude/residual-gauge-invariant-dual from 97bb7fe.
Preserve frozen RGR and all Codex geometry packets. Repair the three
items above, then focus on the exact density bottleneck.

First prove or refute closure(A)=closure(D) in topology (19), with
explicit cocycle derivative estimates and cutoff convergence in
p_R, L2(nu_1) and L1(nu_2). If proved, state the remaining smooth
uniform approximation problem with its exact global growth/local
derivative hypotheses. Establish it with a precisely applicable
theorem or construction, or supply a continuous annihilator/counterexample.
Do not infer it from cylinder density, C^{1,1} regularization or
uncontrolled C^infinity partitions.

Retain K_inf=C ell direct_sum Ann(Dbar) as the correctly scoped
structural endpoint under its declared dependencies. State which
positive-form classification holds on the core, which on Psi*,
and which on all of frozen Psi. Correct W11 and add the complement
counterexample as a negative control. Keep finite witnesses separate
from continuum proofs. Use Python 3.12 only for I/O/orchestration and
the native JS/Paninian engines for scientific calculations. Commit,
push and report the exact SHA and surviving gap.

Codex owns local material-flow representations, renormalized density,
physical observables and the 3D/3+1 source bridge. Do not duplicate
those tasks.
