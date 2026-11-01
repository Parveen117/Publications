# Fable follow-up: repair and audit residual-Diff ergodicity

Continue on `claude/residual-gauge-invariant-dual` in
`Parveen117/Publications`. Own only residual-Diff ergodicity,
invariant-dual classification and rigging uniqueness. Codex owns
nonlinear geometry, detector calibration and apparatus compatibility.
Preserve the shared branch, PR #71 and frozen predecessor proofs.

Audit and repair `ERGODICITY.md` before calling this closure:

1. Section 6 incorrectly identifies `nu_2` as an invariant probability
   measure. RGR gives quasi-invariance and a weighted invariant functional.
   Remove that inference. State separately what, if anything, is proved
   about ordinary invariant measures.

2. Full two-sector unitary equivalence under
   `W = sqrt(c_+)P_+ + sqrt(c_-)P_-` requires both weights strictly
   positive. Classify zero-weight cases as quotients retaining only the
   supported sectors; include the zero form.

3. In Lemma 2, reconcile the finite expression for `u_K` with its
   definition: `G_k` requires `e'_k`, or use `Phi_k` with `e_k`.
   Check every normalization and derivative factor.

4. Make the spectator compatibility contract precise: use BC1's
   continuous extension when J maps outside D_red. BC3 may replace
   the J/sign condition only; retain S and D compatibility. State
   whether invariance is required separately in each gauge slot:
   simultaneous sesquilinear invariance alone does not imply
   factorization through the invariant dual.

5. Audit the main continuum proof in full: coefficient membership in A;
   H^4 composition/inversion bounds; Gaussian D^{1,2} membership;
   Hilbert--Schmidt summability; integrable dominating functions
   independent of K; and convergence of the complete divergence
   expression. Give precise external theorem references and verify
   their hypotheses. Separate proved estimates from assumed RGR
   dependencies.

Deliver corrected theorem statements, a concise audit note listing
changes and remaining gaps, and targeted native JS regression controls
for the actual corrections. Python 3.12 is I/O only. Finite checks do
not certify ergodicity or any continuum theorem. Keep sharp Condition
(L) undecided unless independently resolved.

Commit and push on your separate branch. Report the exact commit, proof
scope, remaining assumptions and whether the main ergodicity argument
survives. Do not expand into the geometry track.

Coordination update: Codex's new MDD-1 packet constructs the additional
joint detector, its source lift, calibration and exact RMG/CMG POVM.
It does not depend on or certify your ergodicity result. Codex next owns
nonlinear continuum geometric observables beyond the J-density family.
