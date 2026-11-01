# Fable next task — from source constraints to the invariant dual

Continue Parveen117/Publications on your separate branch
claude/residual-gauge-invariant-dual from corrected commit
a9f44d5645448cc97dbc78418fd60cae09a90e96. Keep the shared Codex branch
and PR #71 untouched. Read the existing RGR/CGD common-domain
contracts and your corrected ERGODICITY.md, UNIQUENESS.md and AUDIT.md.
Preserve frozen predecessor packets.

Your next target is the actual continuous physical constraint kernel,
not another proof of nonsingular ergodicity:

    K_inf={F in Psi' : F(Q_epsilon(f) psi)=0
           for every smooth real f, both chiralities and psi in D}.

Prove or refute K_inf=C ell for the exact declared source test space
Psi, topology and common source domain D. ell is the existing
doubled-covariance invariant functional. Distinguish this target from
the finite-group invariant dual already studied.

1. Audit D=span{U_g a : a in A}. Establish or refute its density in
   the FULL Psi topology, including the ball-sup and weighted
   integrability seminorms. Hilbert-space density alone is insufficient.
   The frozen RGR statement claims stronger compatibility/density;
   its proof must be made explicit here, without silently changing Psi
   or editing the frozen predecessor.

2. Prove the orbit differentiability needed in Psi, identify which
   domain vectors stay in D, and justify the step from annihilation
   of every Q_epsilon(f) to invariance under its one-parameter group.
   Then justify passage to the complete residual group with its
   actual topology. Do not assume every diffeomorphism is the time-one
   map of a single autonomous vector field.

3. If D is not dense in Psi, characterize its closure and the extra
   annihilators of Psi/closure(D), or state and analyze an explicitly
   revised test-space contract. A counterexample or precisely located
   obstruction is a valid result; do not declare completeness by
   replacing the topology implicitly.

4. Only after that bridge, classify continuous positive slotwise
   constraint-annihilating forms under the corrected spectator
   contract. State each continuity/boundedness hypothesis and allowed
   sector weights. Separate a two-sector unitary equivalence from
   one-sector and zero quotients. Simultaneous invariance alone does
   not give slotwise factorization.

Ownership: Codex is working on physical material-reference momentum,
nonlinear geometry, detector laws and observable domains. Do not work
on those. Do not divert into sharp Condition L, sigma-finite invariant
measures or unrestricted unbounded forms in this task.

Deliver a written theorem or counterexample with every topology,
domain and integration step explicit; native JS/Paninian witnesses and
negative controls only where they test a concrete algebraic obligation;
source pins and a Python 3.12 I/O-only replay. Finite certificates are
not continuum proofs. Update your report with the exact surviving
scope, commit and push on your separate branch, and return the commit
SHA, proof boundary and next genuine gap.
