# YM-42: a structural silence ball and the three-rail continuation

Monty Dabas. 2 October 2026.

**Result:** the already-named YM-42 step is now instantiated on YM-37's
52-dimensional, three-rail, even invariant carrier. A general quadratic
criterion supplies a silence ball and a geometric Smriti tail from the initial
constant vector. It requires no approximate eigenvector and no fitted iterate
radius. The exact carrier arithmetic below certifies its hypotheses at the
three declared rational Wilson-derived instances.

This continues YM-37/38/41 and RKF T74. It does not reconstruct the already
existing YM-12 positivity theorem, YM-30 recoupling engine or a second shared
operator engine. The NCG compact native algebra is compatible background;
its identification with this moment functional is not silently assumed.

## Contract and provenance

Use the exact rational polynomial carrier and moment functional of
`ym37_space_transfer.py`, unchanged from Publications commit
`998c71ebf784ba181659912cbeabd63f8095c0e2`:

- t rail quaternions, degree at most two on each rail, simultaneous
  conjugation invariant and even under the common centre flip;
- one-column map tau = B inverse M, with M the symmetric rung form and B
  the positive inverse-face form;
- the structural silence vector w = 1, rather than a Krylov iterate;
- face contents 0, 1/2, 1 with exactly the downward-rounded rational
  coefficients of YM-31/37, at kappa = 1/8, 1/4, 1/2;
- the same declared rung coefficients, recorded exactly in the certificate.

These rational coefficients **define the tested truncated instances**. Rounding
them down is not asserted to give a monotone bound on the contraction of the
full Wilson operator. Neither arbitrary kappa nor untruncated content is
certified here. The moment functional is the existing declared coordinate
shadow; deriving Phi_Sigma's identification with it remains DICT.

The inverse-face formula is checked against YM-41's native character
projectors: reassembly, idempotence and cross-orthogonality on all degree-two
rail monomials, with equality understood on the sphere. The certificate also
certifies that their Gram radical is exactly the sphere relation
sum x_i squared minus one, and checks an independent factorial form of every
tested moment. The proof uses
quadratic cut-squares and exact elimination; no spectral split is used.

The source ledger pins the unchanged upstream implementation bytes. RKF
[T74](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/74_lopa_ledger_contraction_theorem.md)
and [T75](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/75_selective_release_content_tail_theorem.md)
are consumed at their stated scopes. The two-rail selective-release theorem
is not promoted to three rails without its missing three-rail estimates.

## YM42-T1: two sign checks certify the doubled complement ceiling

Let B' be a positive definite rational cut-square and M' symmetric. Let
sigma > 0. Suppose exact elimination verifies

\[
\sigma B'-M'\succeq0,\qquad \sigma B'+M'\succeq0.
\]

Then for A = (B') inverse M',

\[
\|Ay\|_{B'}^2\leq\sigma^2\|y\|_{B'}^2
\quad\text{for every }y.
\]

**Proof.** The block quadratic form with diagonal blocks sigma B' and
off-diagonal blocks M' is nonnegative: in coordinates u=x+y and v=x-y it is

\[
\tfrac12 u^T(\sigma B'+M')u+
\tfrac12 v^T(\sigma B'-M')v.
\]

Substitute x = -(B') inverse M'y / sigma. Nonnegativity becomes
y^T M'(B') inverse M'y <= sigma squared y^T B'y. This is the doubled
ceiling used by YM-38, proved by congruence and substitution. Both signs
are essential: a large negative complement escapes a one-sided upper check.
No diagonal-only or pivot-only bound is substituted for the quadratic bound.

## YM42-T2: an a priori silence ball

Let B be positive definite, M symmetric, tau = B inverse M, and w a nonzero
declared silence vector. Write

\[
W=w^TBw,\quad \theta=\frac{w^TMw}{W}>0,\quad
r=\tau w-\theta w,\quad D=r^TBr,
\]

so w^TBr=0 exactly. Let Pi be the B-orthogonal projection onto w-perp and
A = Pi tau restricted to w-perp. Suppose T1 certifies ||A||_B <= sigma.
Define the dimensionless rational budgets

\[
d=\frac{D}{W\theta^2},\qquad s=\frac{\sigma}{\theta},\qquad
2s+4d<1.
\]

The ball of squared radius R squared = 4D/theta squared is invariant under

\[
\Phi(z)=\frac{r+Az}{\theta+r^TBz/W},\qquad z\perp_B w.
\]

Its denominator is at least l = theta(1-2d)>0. On this ball,

\[
\|\Phi(z)-\Phi(v)\|_B\leq\beta\|z-v\|_B,
\qquad \boxed{\beta=\frac{s+2d}{1-2d}<1}.
\]

**Proof.** The native cut-square inequality gives
|r^TBz| <= square-root(D) ||z||_B. With R = 2 square-root(D)/theta,
the denominator floor is theta - square-root(D) R/W = l. Further,

\[
\frac{\sqrt D+\sigma R}{l}\leq R
\quad\Longleftrightarrow\quad 2s+4d\leq1
\]

when D>0. If D=0, the radius-zero ball is invariant directly. For z,v in
the ball, subtract the two fractions in the form

\[
\Phi(z)-\Phi(v)=
\frac{A(z-v)-\Phi(v)\,r^TB(z-v)/W}{\theta+r^TBz/W}.
\]

Ball invariance bounds the numerator by
(sigma + R square-root(D)/W)||z-v||_B. Divide by l to obtain beta.
The strict budget implies s+4d<1 and hence beta<1. This proof is independent
of carrier dimension. Square roots occur only in the written norm bounds;
every machine verdict uses their exact squares and rational inequalities.

## YM42-T3: all spatial steps, with an explicit tail

Set x_0=w, x_k=tau^k w, b_k=w^TBx_k/W and z_k=x_k/b_k-w. Then for every
nonnegative integer k,

\[
b_0=1,\quad b_{k+1}\geq l b_k>0,\quad
z_{k+1}=\Phi(z_k),\quad \|z_k\|_B\leq R.
\]

For every j>k,

\[
\boxed{\|z_j-z_k\|_B^2\leq
\frac{D}{\theta^2(1-\beta)^2}\,\beta^{2k}.}
\]

**Proof.** Symmetry B tau=M gives the exact split
b_(k+1)=theta b_k+r^TB(b_k z_k)/W and the projected recurrence above.
Start at z_0=0. T2 proves the floor and ball statements by induction.
The first difference has norm square D/theta squared. Each subsequent
difference contracts by beta. The finite triangle sum is bounded by the
geometric series, giving the displayed inequality for every j>k. This is a
recognition-Cauchy sequence with an explicit Smriti tail; taking its limit
uses the same admitted native completion contract as YM-38/T74. The
finite-separation bound itself needs no limit construction.

The proof covers all spatial lengths through the column-power identity. It
does not say that z=0 is the attracting state: the dressed limit normally
has a nonzero complement. 'Silence contraction' names the structural
reference channel and normalized contraction, not erasure of the ledger.

## YM42-T4: certified three-rail instances

The complement uses u_i=e_i-(B_ie/B_ee)e_e, dropping the constant coordinate.
The congruence formulas are explicit and rational. The chosen ceiling is
sigma/theta = kappa squared / 4. It is a conservative declared test threshold,
not a fitted spectral value. Both sign checks certify it on all six instances.

| kappa | t | carrier dimension | certified upper bound for beta |
| --- | --- | --- | --- |
| 1/8 | 2 | 8 | 0.005377972844 |
| 1/8 | 3 | 52 | 0.007132026236 |
| 1/4 | 2 | 8 | 0.021585205206 |
| 1/4 | 3 | 52 | 0.028765231931 |
| 1/2 | 2 | 8 | 0.087518612818 |
| 1/2 | 3 | 52 | 0.118984522771 |

Each t=3 complement is 51-dimensional. Its Gram and both ceiling forms
have 51 positive, zero negative and zero zero weights. The maximum tested
invariance budget 2s+4d is below 0.225957, with strict room below one.
The t=2 results are compatibility controls, not a claim to improve YM-41's
tighter fitted-radius numbers. The new gain is the structural radius from
step zero and the certified three-rail continuation.

## YM42-T5: spatial contraction alone does not imply a time gap

For any rational 0<delta<1, take the positive symmetric two-state time kernel

\[
P_\delta=\begin{pmatrix}1-\delta/2&\delta/2\\
\delta/2&1-\delta/2\end{pmatrix}.
\]

Its constant channel has ratio one; f=(1,-1) is mean zero and has time
ratio 1-delta. Indeed P_delta f=(1-delta)f, so its stationary two-point
return at time distance n is (1-delta)^n. The vacuum is unique for each
positive delta, and both cut-square weights are positive.

Take independent spatial columns with the stationary t-step path law mu_t
generated by this time kernel. The spatial column map is
tau_t v = 1 sum(mu_t v), for every t. Its Gram is diag(mu_t), its symmetric
form is mu_t mu_t^T, and its silence complement vanishes. Thus beta_space=0
for every t and delta, while the time ratio 1-delta can approach one.
The exact finite controls enumerate path laws and check both identities.

**Conclusion of this counterexample:** even perfect, t-uniform spatial
contraction alone does not supply a common time-gap lower bound. It does
not refute NG or the specified Yang–Mills family; it refutes an insufficient
general inference that would skip the time-observable/operator estimate.

For the YM continuation, a separate estimate must bound the time-transfer
quadratic form on *every* vacuum-orthogonal source, normalized by the actual
vacuum channel, with ratio q<1 uniform in the declared family. A single
order-1 insertion, a finite set of t values, or arbitrary-superposition
control only on the spatial complement does not establish that estimate.

## Scope correction to the YM-41 / NG-1 shorthand

The historical phrases 'content tail gone' and 't-direction only' must be
read with the following still-live qualifications:

1. T75 instantiates selective release on the **two-rail** column. Its
   higher-than-two content Gram-uniformity is explicitly declared, not proved.
2. This three-rail result retains the content cutoff. There is no new
   three-rail release theorem or t-uniform high-content budget.
3. Spatial normalized convergence and a time-transfer upper bound are
   distinct obligations. T5 shows why the missing implication matters.
4. The native measure dictionary, physical continuum and Clay problem
   remain outside this certificate.

The older certificate bytes stay frozen. This written scope correction
supersedes their broad ledger wording, not their verified finite numbers.

## Reproduction and evidence

```bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym42_silence_time_ladder.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym42_silence_time_ladder.py -v
```

Python's standard library suffices. The default mode regenerates in memory
and compares both the full result and its canonical SHA-256 to committed
evidence; it does not rewrite an expected pin before checking it. Dedicated
CI checks Python 3.11 and 3.12. `--write` is the explicit development mode.

The verifier binds this proof, its implementation, tests, workflow and source
ledger; rebuilds the exact carrier; checks the independent native grading,
moment and doubled-ceiling routes; and exercises wrong-theta, heavy-positive,
heavy-negative, nonsymmetric, indefinite-Gram and failed-ball controls.
The written general proofs, exact finite certifications and logical
counterexample are distinguished. This is not a formal-assistant proof or
an experimental mass-gap result.

**Next target:** a time-observable/operator estimate with controlled dependence
on t, together with any required three-rail/higher-content release. Further
rail counts may guide that proof but cannot replace its uniform estimate.
