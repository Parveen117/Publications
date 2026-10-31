# Source and certificate review

Author: Monty. 31 October 2026. Internal review, not independent peer review.

1. **Same target.** The calculation evaluates the original uniform-vacuum
   characteristic formula. It adds source conditioning without changing
   LMD, its h range, the pulse, detector or '+' event definition.
2. **Fourier normalization.** The orthonormal first-pair coefficients have
   variance h/2; the coefficients of raw cos/sin have variance h. The Wick
   subtraction through the first pair is h, not h/2. The Gaussian polar
   variable R^2/(2h) is therefore Exp(1).
3. **Conditional identity.** The Wick martingale gives E[J|q,X,Y]=J_1.
   Orientation integration makes J_1 measurable in q,R, justifying the
   tower identity E[J|q,R]=J_1. Merely choosing an approximate radial
   value would not justify cancellation of the first Taylor term.
4. **Positive bins and all tails.** Each Gaussian and radial cell has a
   positive probability and exact conditional first moment, enclosed with
   directed rationals. Both infinite q tails and the infinite radial tail
   are retained. Removing the radial tail fails normalization.
5. **Series remainder.** The positive Poisson/Gamma series is truncated at
   index 32 with its geometric Poisson-tail majorant. Gamma bin probabilities
   are intersected with [0,1] using their known analytic probability property.
   There is no heuristic cancellation or zeroing of a scientific error.
6. **Second moment.** The product partition's captured second moment factors
   into the q and radial moments. The ideal first-pair moment agrees with
   the independently computed frozen Laurent bound. The complete residual
   is decomposed into higher-mode, q-bin and radial-bin budgets.
7. **Nonnegativity and refinement.** Every budget component is positive on
   the declared cases. Both nested product partitions are evaluated, and
   their total variance budgets strictly decrease. Existing and refined
   probability intervals overlap as two enclosures of the same target.
8. **Phase provenance.** All nine inherited NRQC phase records are pinned.
   No fresh calculator execution is claimed. Interpolation coefficients are
   classical postprocessing. Independent native trig checks cover every
   retained node and seven selected cells per case/refinement, not the
   entire assembled sum.
9. **Tolerances.** PASS describes successful arithmetic/consistency checks.
   The four separate tolerance statuses retain one unresolved case. The
   combined 11/12 statement explicitly includes eight earlier certifications;
   it does not imply twelve new executions or hardware validation.
10. **Post-run diagnosis.** MODE_BUDGET.json and NEXT_STEP.md are a later
    planning diagnostic using frozen moments. Their conclusions did not
    change the case selection, tolerance, partition or source in the run.
    N=8 is the first tested candidate, not a proved minimal or completed
    finite-source algorithm.
11. **Amplitude-only bound.** The all-radii conditional martingale is L2
    bounded. The positive-series bound B_0(x)<=exp(x^2/4) bounds its captured
    second moment even with infinitely many radii. The transferred bound
    exceeds the target for this certificate; it is not a physical no-go or
    a lower bound on actual probability error. Relative spatial phases
    are one concrete next source resource, exhibited by the two-mode mass.
12. **Open physics.** A source expectation solver is not a gravitational
    Hamiltonian, apparatus experiment, calibration-free geometry or 3+1
    constraint solution. Frozen predecessor obstructions are preserved.
