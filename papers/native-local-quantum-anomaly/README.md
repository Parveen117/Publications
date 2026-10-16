# LQA-1 — Local quantum constraints: anomaly and a positive kinematic alternative

The next working representation closes the selected inhomogeneous V=0,
matter-free constraint algebra **without a central anomaly**, on a common
dense algebraic domain with positive kinematic norm. Reversing the
negative-sign oscillator vacuum gives central contributions +1 and -1;
the inherited null improvement leaves their cancellation intact. No new
counterterm or zero-mode shift is required. This is a selected standard
free-current representation, now derived for the actual LNC source.

Physical constraint states and their inner product remain open. In particular,
the algebraic finite-polynomial test domain has no nonzero common all-mode
kernel. Its mode grading is unbounded in both directions. The positive
kinematic construction is not yet a positive physical quantum theory.

The comparison with the conventional vacuum is also proved: with both null
currents using positive-index annihilators, the normal-ordered generators
have central charge **2 per chirality**.
Their derivative improvement is null and cannot cancel it. Consequently the
unchanged generators have no nonzero common all-mode constraint kernel, even
after a constant zero-mode shift. This is a representation-specific obstruction,
not a no-go theorem for quantizing the classical model.

Within the explicit linear-improvement family

    L_m(r,s,gamma) = S_m + i m (r a_m + s b_m) + gamma delta_m0,

the algebra closes without a central term exactly when

    r s = -1/12,    gamma = -1/12                 (hbar = 1).

Keeping the inherited coefficient r=1 selects s=-1/12 in this family.
Restoring the action scale gives s=gamma=-hbar/12: a possible quantum
counterterm, not an automatic consequence of the cut. Its zero-mode equation
differs from HQS, and the oscillator inner product is still indefinite.
This candidate is documented, **not adopted as a completed physical theory**.

Read [THEOREM.md](THEOREM.md) for normalization, all-mode written proofs,
common invariant domain, the mixed H/D anomaly, and the homogeneous comparison.
[SOURCE_MAP.md](SOURCE_MAP.md) credits the established Virasoro/free-field
methods and explains exactly what is reused. Frozen predecessors are unchanged.

Run with Python 3.12 and the existing native engines:

```sh
python3 papers/native-local-quantum-anomaly/run.py --check
```

The certificate passes **407 exact checks, 9 symbolic proofs/replays and
25 negative controls**. It uses native rational/cut arithmetic for finite oscillator
actions, central contractions, counterterm controls and serialized Paninian
word proofs. These are finite executable witnesses; the infinite-mode theorem
and domain statements have written proofs. No finite matrix cutoff is used to
infer a continuum commutator. No positive local physical Hilbert space,
coframe exponential operator, full quantum gravity or FTL is certified here.

Next: construct a justified physical extension of the reversed-vacuum test
domain, solve the constraints, derive a positive physical inner product and
compare the reduced sector with HQS. The matching zero-mode equation PK=0
alone is not a proof that the two physical quantizations agree.
