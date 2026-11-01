# MRM-1 — Physical material-reference momentum

This packet constructs the total source reference momentum as a
self-adjoint physical observable on the selected RGR/MPQ space. It
separates that operator from a divergent noisy joint-readout estimator
and from a nonconstant local-current implementability obstruction.

Read [THEOREM.md](THEOREM.md) for the analytic proofs and domains,
[SOURCE_MAP.md](SOURCE_MAP.md) for dependencies and scope, and
[FABLE_NEXT_PROMPT.md](FABLE_NEXT_PROMPT.md) for the disjoint residual-gauge
task following the corrected external report at a9f44d5.

    python3.12 papers/native-material-reference-momentum/run.py --check

Python performs I/O and orchestration only. Scientific fixtures use the
unchanged native exact arithmetic and Paninian engines. Certificates
check finite identities and negative controls; they do not mechanically
prove the analytic continuum theorems.

The total charge has spectrum 2 pi hbar Z and is not positive energy.
Nonconstant material-label flows fail in this physical Gaussian measure
class. This is distinct from the source gauge action in inverse
embeddings; the earlier gauge representation and selected rigging remain
unchanged. Local metric/curvature observables and full physical quantum
gravity are still open.
