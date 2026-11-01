# Retained actual-source evaluation

Input is the `actual_depth_compatibility` fixture in the byte-pinned
LDC certificate, with its LAA source lineage: s=29,31,33, a=(s-1/2)^2,
q(a) and c_0,...,c_19. LAA-T1 supplies the universal coefficient bound.
No source function, higher derivative, zero or spectral root is evaluated.

At each anchor, native arithmetic evaluates tau=1/2,1,2. Each reading
retains the coarse S_10 enclosure and the refined P_19 enclosure with
tail cutoff J=24. All finite missing terms through j=23 are bounded
explicitly; the remaining infinite tail is geometric. The larger cutoff
is for bounding unknown terms and does not claim 47 observed derivatives.

The certificate records exact rational endpoints, the certified q0
interval, M, C, A, both partial values, missing-tail budgets, and the
normalized interval B/c_1. Tests check radius improvement, consistency
between enclosures, matching of S_10 with the retained LWP prefix,
initial value and retention of input uncertainty. Compatibility with a
positive-model interval is not evidence of actual all-order positivity.

The finite certificate additionally checks signed majorant saturation,
two exact coefficient regroupings, positive rational-atom comparison,
and an entire nonunitary counterexample. Infinite convergence is proved
in THEOREM.md, not inferred from a finite test count.

Execution and certificate emission occur on GitHub Actions. The imported
provenance identifies the tested source tree and generated output hashes.
No local scientific replay or independent audit agent is used.
