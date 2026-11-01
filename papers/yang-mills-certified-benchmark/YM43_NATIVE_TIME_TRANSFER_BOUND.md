# YM-43: a time-transfer bound from local overlap and cut-squares

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result under the declared positive-functional carrier:** on the full
heat-kernel bounded-overlap chain used by YM-19, a local, anisotropic
overlap criterion bounds the normalized time transfer on **every**
vacuum-orthogonal source, uniformly in chain length. The three certified
parameter cells admit ratios q = 1/2, 1/4 and 1/64 respectively.

YM-42 controlled spatial normalized convergence at fixed rail count. This
chapter proves the additional time-direction estimate, using the full
positive kernel in a coarse heat-kernel regime. It does not extrapolate
YM-42's truncated Wilson numbers or T75's two-rail release to all time rails.

## Source, cut and adapter contract

Begin with a declared positive normalized recognition functional on a
single-site coefficient carrier, its finite tensor products, positive bond
weights, and lawful conditional normalization. A finite cut is a finite
positive weighted sum. The compact-coordinate adapter uses the already
admitted SU(2) reference integral and bounded measurable coefficient
functions. Positivity, product integration and uniform-limit continuity are
part of this **declared functional contract**, not newly derived from EMK.
The identification of this reference functional with Phi_Sigma remains DICT.

The source is the finite strip of coefficient configurations. The observer
is an arbitrary bounded function of one complete spatial row. The retained
channel is the normalized positive vacuum; its complementary cut contains
every source orthogonal to that vacuum. Local conditional replacement is a
proof device for controlling unresolved boundary memory, not a postulate
about physical stochastic events. Time below means transfer-step count;
physical units and the continuum limit remain separate.

All norms follow the positive cut form, null-seam quotient and declared
completion. Cauchy–Schwarz follows by nonnegativity of the cut-square of
u-zv for every scalar z. No Hilbert space, eigenvector decomposition,
Perron/Jentzsch theorem or Dobrushin comparison theorem is a proof premise.
Where infinite sequences are used, explicit geometric Smriti tails are
given, in the sense of RKF T28. T53 elimination supplies finite controls.

The local-comparison method has established Dobrushin lineage, already cited
by YM-19. See the authors' primary paper
[Rebeschini–van Handel, *Comparison Theorems for Gibbs Measures*](https://web.math.princeton.edu/~rvan/dobrushin130819.pdf)
for that lineage. This is a self-contained derivation and scoped integration,
not a claim of priority for the general comparison method.

## YM43-T1: local overlap bounds are constructive

At site i write the conditional normalized readout as

\[
p_i(dx\mid\eta)=\frac{\prod_{j\sim i}w_{ij}(x,\eta_j)}
 {\Phi_i(\prod_{j\sim i}w_{ij}(\cdot,\eta_j))}\,\Phi_i(dx).
\]

Suppose a bond has pointwise bounds 0<L<=w<=U and delta=U/L. If only
neighbour j changes, the ratio of unnormalized weights lies in
[delta inverse, delta]. The ratio of normalizers has the same bounds;
hence p/q>=delta inverse squared. Their common positive part has mass
at least delta inverse squared. Consequently their mismatch can be bounded by

\[
c_{ij}=1-\delta^{-2}.
\]

**Proof of the comparison operation.** Let a=min(p,q), pointwise relative to
their common reference, and r=1-Phi(a). Draw a common value with weight a;
with the remaining weight use the product of the two normalized nonnegative
remainders. Both marginals are exactly p and q, and mismatch mass is at most
r. This is explicit positive-sum algebra; in the compact adapter it is the
same construction with the declared integral. For several changed neighbours,
insert intermediate configurations one neighbour at a time. The half-L1
triangle inequality, obtained pointwise, gives the sum of their c_ij bounds.
No independence of the changes is assumed.

## YM43-T2: finite conditional updates give a boundary-memory comparison

For two boundary configurations let C contain the nonnegative interior
influences and b the sum of influences from changed boundary sites. Suppose
max row sum(C)<=alpha<1. Construct the two finite-strip positive laws directly
by normalizing the product of bonds. Start any joint law with these marginals.
At each proof step choose one of the N interior sites with weight 1/N and
apply T1's common-part replacement to that site's two conditionals.

Each marginal remains its original strip law: summing over the replaced
coordinate verifies this identity exactly. If p_i(k) is mismatch mass at
site i, then

\[
p(k+1)\leq A p(k)+b/N,\qquad
A=(1-1/N)I+C/N.
\]

If v>=0 satisfies v>=Cv+b, induction gives

\[
p(k)\leq v+A^k\mathbf1,
\qquad \|A^k\mathbf1\|_\infty\leq
\left(1-\frac{1-\alpha}{N}\right)^k.
\]

For any bounded observable F, let osc_i(F) be its change when only site i
changes. A coordinate-by-coordinate telescoping gives

\[
|\Phi(F)-\widetilde\Phi(F)|
\leq\sum_i\operatorname{osc}_i(F)v_i.
\]

**Proof.** The same expectation difference is bounded at every proof step
by sum osc_i(F)p_i(k). Substitute the finite recurrence and its geometric
tail, then send k to infinity. No stationary joint law or infinite-volume
existence theorem is needed. The construction and bound apply directly to
every finite positive cut and to the declared compact integration adapter.

## YM43-T3: an anisotropic time barrier, all widths and strip lengths

On the open strip of spatial width m and time rows 0,...,L, suppose each
site has at most two space neighbours with influence <=c_s and two time
neighbours with influence <=c_t. Choose a rational 0<q<1 such that

\[
\boxed{\eta(q)=2c_s+c_t(q+q^{-1})<1.}
\]

This also implies alpha=2c_s+2c_t<1. The comparison barrier

\[
v_{i,r}=q^r+q^{L-r},\qquad 1\leq r\leq L-1,
\]

satisfies v>=Cv+b for **arbitrary** differences on both end rows.

**Proof.** Extend the displayed formula to r=0,L, where v>=1 bounds every
boundary mismatch. At an interior site the sum of all four neighbour
contributions is at most eta(q)v_(i,r), because
v_(r-1)+v_(r+1)=(q+q inverse)v_r. Missing spatial neighbours only improve
the inequality. T2 now gives for an arbitrary whole-row F at row n,

\[
|\Phi^{x,z}(F)-\Phi^{x',z'}(F)|
\leq m\,\operatorname{osc}(F)(q^n+q^{L-n}).
\]

No claim that this is uniform in m as a *prefactor* is made. The rate q is
uniform in m; T5 removes the fixed-m prefactor in the operator conclusion.

## YM43-T4: construct the positive vacuum with a declared tail

For each fixed m let T be the symmetric row kernel

\[
T(x,y)=W(x)^{1/2}\,\prod_{i=1}^m K(x_i,y_i)\,W(y)^{1/2},
\quad W(x)=\prod_{i=1}^{m-1}w_s(x_i,x_{i+1}).
\]

The positive scalar square root is in the completed cut-scalar carrier.
Assume 0<ell<=T(x,y)<=U<infinity. These bounds may depend on m.
Let f_k=T^k1 and r_k=f_(k+1)/f_k. If a_k=inf r_k, b_k=sup r_k, then

\[
[a_{k+1},b_{k+1}]\subseteq[a_k,b_k],\qquad
b_{k+1}-a_{k+1}\leq\zeta(b_k-a_k),\quad
\zeta=1-(\ell/U)^2<1.
\]

**Proof.** r_(k+1)(x) is the positive normalized average of r_k under
weights proportional to T(x,y)f_k(y). The weights for two x values have
ratio at least (ell/U) squared. Subtract their common part as in T1;
the range of their averages is bounded by zeta times osc(r_k). In particular
a_0>=ell>0 and the nested scalar intervals have one common limit lambda.

Set H_k=f_k/Phi(f_k). For k>=1, ell/U<=H_k<=U/ell, and

\[
\|H_{k+1}-H_k\|_\infty\leq
\frac{U}{\ell}\,\frac{b_0-a_0}{a_0}\,\zeta^k.
\]

This follows by writing H_(k+1)/H_k=r_k/Phi(H_k r_k), whose denominator
lies in [a_k,b_k]. The geometric tail constructs a uniform limit h,
with Phi(h)=1 and ell/U<=h<=U/ell. Uniform continuity of the bounded
kernel map gives Th=lambda h. This proves the positive vacuum rather than
assuming a spectral theorem. Uniqueness follows also from T5.

Define

\[
Pf=\frac{T(hf)}{\lambda h},\qquad
\pi(dx)=\frac{h(x)^2}{\Phi(h^2)}\Phi(dx).
\]

Direct substitution proves P1=1, invariance of pi and symmetry for its
positive cut form. Conditional cut-square Cauchy–Schwarz gives
||Pf||_pi<=||f||_pi. All finite path laws pi(x_0) product P(x_r,x_(r+1))
are normalized and consistent by finite summation. No infinite path law is
required in the next argument.

## YM43-T5: from the time barrier to the full complementary operator bound

The finite P-path weight telescopes to

\[
\frac{h(x_0)h(x_L)}{\Phi(h^2)\lambda^L}
\prod_{r=0}^{L-1}T(x_r,x_{r+1}).
\]

Conditioning on its two end rows removes the endpoint factors. Every
interior row has exactly the space weight W and the declared time bonds K.
Thus T3 applies. Mixing over the end row x_L, and then letting L grow at
fixed m,n, gives

\[
\operatorname{osc}(P^nF)\leq m\operatorname{osc}(F)q^n
\]

for every bounded whole-row observable F. The discarded far-end memory
is at most m osc(F) q^(L-n), an explicit Smriti tail.

Let pi(F)=0 and u_n=||P^nF||_pi squared. Symmetry and cut-square
Cauchy–Schwarz give

\[
u_n^2\leq u_{n-1}u_{n+1},\qquad
u_n\geq u_0(u_1/u_0)^n
\]

when u_0,u_1>0. Zero cases are immediate. On the other hand stationarity,
mean zero and the oscillation bound give

\[
0\leq u_n=\langle F,P^{2n}F\rangle_\pi
\leq m\,\pi(|F|)\operatorname{osc}(F)q^{2n}.
\]

If u_1/u_0>q squared, these inequalities contradict growth of a positive
geometric power. Therefore

\[
\boxed{\|PF\|_\pi^2\leq q^2\|F\|_\pi^2
\quad\text{for every mean-zero source}.}
\]

The bounded core is dense by the declared positive-form completion, and P
is contractive, so the same inequality extends to the completed core.
For derived complex cut-scalars, apply the real estimate to the real and
iota components and add their squares; the kernel preserves both components.
Multiplication by h transports it isometrically, up to the fixed factor
Phi(h squared), to T/lambda on the h-orthogonal source cut Q:

\[
\boxed{\|Q(T/\lambda)Q\|\leq q<1
\quad\text{for every finite spatial width }m.}
\]

This is an operator bound on arbitrary superpositions, not a one-probe
correlation claim. It also implies uniqueness of the vacuum channel. The
proof does not identify individual elimination weights with eigenvalues.
Symmetry and Th=lambda h give QT=TQ, so the same argument yields
||Q(T/lambda)^nQ||<=q^n for every transfer-step count n.

## YM43-T6: transfer the bound to the existing square-sourced chain

The existing chain uses S=K_half M_W K_half, whereas the preceding
row kernel is T=M_sqrt(W) K M_sqrt(W). The full coefficient semigroup gives
K=K_half squared, already supplied by YM-9/12 on the declared carrier.
Set X=K_half M_sqrt(W). Then S=XX dagger and T=X dagger X.

If Th=lambda h, put g=Xh/square-root(lambda). This has the same norm
as h and Sg=lambda g. If v is orthogonal to g, then X dagger v is
orthogonal to h. The positive-form bound T<=q lambda on that complement
implies

\[
\|Sv\|^2
=\langle X^\dagger v,T X^\dagger v\rangle
\leq q\lambda\langle v,Sv\rangle
\leq q\lambda\|v\|\|Sv\|.
\]

Dividing when Sv is nonzero proves ||Sv||<=q lambda ||v||. Thus the
bound is for the actual YM-19 A-chain transfer, not a differently ordered
proxy. The independent B-channels of YM-16 have free ratio exp(-3a/4).
Where that ratio is <=q, expand the tensor product into its mutually
orthogonal vacuum/complement cuts: every nonvacuum block contracts by at
most q. This supplies the same bound for the factorized A/B chain.

## YM43-T7: exact certified parameter cells

Use YM-19's unchanged all-content heat-kernel series and its explicit tail:

\[
S_a=\sum_{n\geq1}(n+1)^2\operatorname{Exp}_\Sigma
       [-a n(n+2)/4],\quad 1-S_a\leq K_a\leq1+S_a.
\]

For the space bond Exp_Sigma(kappa Tr(xy inverse)/2), its ratio is
delta_s=Exp_Sigma(2 kappa). With outward coefficient bounds,

\[
c_t\geq1-\left(\frac{1-S_a}{1+S_a}\right)^2,\qquad
c_s\geq1-\operatorname{Exp}_\Sigma(-4\kappa).
\]

The certificate rounds these **up**, then checks eta(q)<1 exactly. It
does not drop the character tail or replace the full heat kernel by a
finite-content kernel. The reference functional/character-convolution
dictionary remains the same admitted coordinate adapter as YM-19.

| a | kappa | certified q | upper bound for eta(q) | gap per transfer step |
| --- | --- | --- | --- | --- |
| 6 | 1/16 | 1/2 | 0.850216 | >= Log_Sigma(2) |
| 8 | 1/8 | 1/4 | 0.952217 | >= 2 Log_Sigma(2) |
| 12 | 1/8 | 1/64 | 0.913217 | >= 6 Log_Sigma(2) |

For all three cells, the free B-ratio is below q. The inequalities hold
for every m and every transfer power in this declared chain family.
The displayed decimal eta bounds are upward rounded; exact rationals and
directed logarithm bounds are stored in the certificate. Division by the
declared a gives the corresponding rate in that parameter's units.

The previous refusal cells (6,1/8), (4,1/16), (10,1/4) remain refused by
this coefficient budget: alpha>=1 implies eta(q)>=alpha for every 0<q<1.
This refuses the route, not the physical existence of a gap at those cells.

## What is now closed, and what is not

| Obligation | Current status |
| --- | --- |
| YM-19's cited comparison/extraction step | Replaced by written T1–T6 proofs within the declared positive-functional/heat-kernel carrier |
| Arbitrary-source time-transfer bound, all finite m, in the three coarse cells | Proved under that carrier contract; exact parameter hypotheses certified |
| Perron simplicity used as an external premise | Unneeded on this route: T4 constructs h and T5 proves uniqueness |
| All-content kernel estimate on this route | Explicit convergent coefficient-tail bound; no T75 higher-content Gram assumption is used |
| YM-42's original truncated Wilson grid at a different kernel/parameter regime | Not identified with these cells; its pending operator bound is unchanged |
| Uniform cutoff a tending to zero, AF trajectory, generic E4D-C, NG for all a | Open; the coarse-kernel overlap criterion deteriorates as the kernel sharpens |
| Native Phi_Sigma/reference-measure identification and NCG observable dictionary | Open |
| Four-dimensional interacting continuum, Clay problem, quantum gravity | Not established |

The historical YM-19 certificate remains byte-preserved and labelled ANCHORED.
The new theorem supersedes its missing step only in the stated carrier and
certified coarse cells. Re-deriving a known comparison mechanism does not
make the admitted measure, action or physical gauge-group choice primitive.

## Executable evidence and runtime scope

```bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym43_native_time_transfer.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym43_native_time_transfer.py -v
```

The standard-library verifier uses exact Fraction arithmetic, the existing
outward exponential/logarithm arithmetic and native T53 elimination. It
checks the parameter cells, temporal barriers on varied strips, explicit
common-part couplings and exact whole-row marginal comparisons. Independent
finite rational transfer matrices check the complementary threshold and the
vacuum iteration. Refusal controls catch an omitted return-time term, wrong
normalization, hidden slow source, omitted vacuum projection, nonpositive
kernel and noncontracting budget. The general proof is written, not
mechanically formalized; finite fixtures are not labelled the universal proof.

Source pins bind the unchanged dependencies, this proof, the verifier,
tests and workflow. The default is a read-only fresh certificate comparison.
Only Python 3.12 is scheduled; `--write` explicitly regenerates evidence.
The accompanying YM-42 runtime-only rebind changes its workflow/proof input
hashes while checking that every mathematical output is unchanged.
