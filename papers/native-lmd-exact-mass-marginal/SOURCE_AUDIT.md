# Source applicability and evidence map

Author: Monty. 31 October 2026. Source parent: `ace56a0a3faa21708e4a6779c8e4a8f02dab70e0`.

The new dependency is Guillaume Rémy, *The Fyodorov-Bouchaud formula and
Liouville conformal field theory*, [arXiv:1710.06897v3](https://arxiv.org/pdf/1710.06897v3),
dated 4 November 2019. Theorem 1.1, equation (1.9), on page 3 supplies the
normalized circular GMC total-mass law. Equations (1.2), (1.3), (1.5), on
page 2 supply covariance, Wick convention and normalized circle measure.
The unversioned PDF was retrieved on 31 October 2026 and its header identifies
v3; the pinned version above identifies the dependency. This external theorem
is credited, assumed as proved, and not numerically re-proved by this packet.

| Item | Existing source | Applicability calculation |
|---|---|---|
| Covariance | [LMD section 2, equation (3)](../native-logarithmic-material-density/THEOREM.md) | K=1-log|chord|; subtract the independent constant mode of variance h, then scale the external field by beta/2 with beta^2=2h. |
| Zero mode | Same source, paired Fourier construction | q is independent of nonzero modes; retain exp(q-h/2). |
| Measure convention | LMD-T3 | z in R/Z gives dz=dtheta/(2 pi); no extra circumference factor. |
| Wick and cutoff convention | LMD-T3 and its Fourier/mollifier identification | Nonzero-mode subtraction is (h/2) sum 1/n. No additional smooth covariance is allowed. |
| Actual measured target | [Pulse characteristic section 2](../native-lmd-pulse-characteristic/THEOREM.md) | Only chi_h(a), the distribution of total mass, enters the selected vacuum uniform marginal. |
| Independent second moment | Pulse characteristic equation (1) and pinned source intervals | Euler sine integral provides I(h), independently of the imported law's application. |
| Conditional-mean error and phase nodes | [Certified solver](../native-lmd-certified-pulse-solver/THEOREM.md) | Reuse the quadratic Taylor remainder and nine exact phase outputs. |
| Relative-phase obstruction | [First-mode diagnosis](../native-lmd-first-mode-pulse-solver/NEXT_STEP.md) | Its amplitude-only projection is unchanged. A full collective mass law is a different auxiliary representation. |

`SOURCE_PINS.json` hashes all inherited numerical inputs and their mathematical
contracts. The new implementation integrates Gamma functions with rational
enclosures, records both unbounded source tails, and verifies independent
recurrences. Matching two moments alone is explicitly insufficient; the
frozen campaign includes a moment-matched two-point counterexample.

The finite output does not verify the external theorem, an operator
intertwiner, a joint field distribution, nonuniform dynamics or a physical
apparatus. This packet changes the scalar expectation algorithm while
retaining the existing physical model. Broader fields or covariance kernels
require a new applicability argument; the law cannot be carried over by name.
