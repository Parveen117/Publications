# YM-49: time-zero observables and an exact row-memory closure test

Monty Dabas. 3 October 2026. Runtime: **Python 3.12 only**.

**Result on the existing adapter:** bounded local coefficient observables
act as a unital dagger representation on YM-48's reflected history
completion. The construction proves the time-zero collision limit, its
norm bound and its product law. It also reproduces ordered coefficient
histories. This supplies an observable action alongside the already
constructed continuous heat-time dynamics and self-adjoint generator.

For the instantaneous coefficient row, a positive two-time defect measures
exactly the norm lost into its orthogonal complement. It gives a necessary
and sufficient closure criterion, and the existing RKF memory equation
now has a specified infinite carrier and geometric remainder bound.
The criterion is **not evaluated as zero or nonzero for the interacting
infinite chain**. Its retained-memory channel must not be discarded.

## Inputs and scope

Keep the full SU(2) coefficient heat/reference-functional adapter,
abs(theta)<1/1680, kappa(a)=theta a, and the spatial chain. Reuse YM-47's
functional Omega, YM-48's H_ref, vacuum e, semigroup T(t), generator H,
coefficient-row isometry J and rate gamma>0. Times are normalized heat
parameters. Neither a physical clock nor unitary real time is selected.

Write A_coeff for local finite coefficient polynomials, closed under
complex conjugation and multiplication. Its uniform completion is
A_unif. No density claim for *all bounded measurable* row functions is
built into that definition. Let R=Ran J, and distinguish it from
K, the closed span in H_ref of coefficient-polynomial future histories.
Then R subset K subset H_ref; these equalities are development questions.

RKF N11's bounded-action obligation, N09's exact block elimination,
N08's bounded inverse algebra and T28's explicit completion tails are
reused. N09 is not rediscovered as a new mathematical law. The new work is
the time-zero multiplication bound and its collision construction on this
carrier, followed by a carrier-specific positive closure diagnostic and
gap-controlled memory. The canonical RKF engine and prior evidence are
unchanged. YM-48 records the reflection-positive reconstruction lineage.

## YM49-T1: a uniform form multiplier bound

At finite spatial interval I, retain YM-48's
H_I=L_I-B_I+E_I, positive unit vacuum h_I and V_I(t). On the free
coefficient form domain

\[
\mathcal Q_I=\{v:\sum_C C|v_C|^2<\infty\},
\qquad q_I(v)=\sum_{i,\alpha}\|D_{i,\alpha}v\|^2
                 +\langle v,(E_I-B_I)v\rangle ,
\tag{1}
\]

the form is nonnegative and closed. Here C is the total coefficient
Casimir. The bounded potential makes its form norm equivalent, at this
fixed I, to the free form norm after adding a sufficiently large multiple
of ||v||^2. Coordinate truncation is a core. Nonnegativity extends from
D(L_I) by that approximation. Lower semicontinuity follows from the
nonnegative coefficient sum and continuity of the bounded-potential term.

We need a domain argument before using the vacuum as a change of source.
At each fixed I, h_I is smooth and bounded above and away from zero:

* Multiplication by the finite coefficient polynomial B_I is bounded on
  every coefficient Sobolev norm. Tensor fusion shifts each spin by a
  bounded amount, so the weights (1+C)^p on coupled source/target
  blocks are comparable; only finitely many blocks meet each block.
  The equation L_I h_I=(B_I-E_I)h_I therefore bootstraps h_I to every
  power domain of L_I.
* The number of coefficient modes below a Casimir cutoff and the
  sup norms of their derivatives have polynomial growth at fixed I.
  Cauchy--Schwarz with a sufficiently high power of (1+C) proves
  uniform convergence of each differentiated coefficient series.
  Thus the preceding power-domain statement gives smoothness.
* YM-44's lower kernel comparison U_I(t)>=exp(-b_I t)K_I(t), at any
  fixed t>0, and the strictly positive continuous finite heat kernel
  give h_I>=c_I>0 from h_I=exp(-E_I t)U_I(t)h_I.

These constants may depend on I; they are used only for this domain/core
argument. For any coefficient polynomial v, v/h_I is smooth. Its
coefficient truncations converge with first derivatives, by the same
coefficient estimate. Multiplying back by h_I proves that h_I A_coeff(I)
is a form core for (1).

Expansion and the vacuum equation, as in YM48-T6, now give on this core

\[
q_I(h_I f)=\Phi_I\!\left(h_I^2\sum_{i,\alpha}|D_{i,\alpha}f|^2\right).
\tag{2}
\]

For F in A_coeff set f_F=||F||_infinity and
C_F=||sum |D_(i,alpha)F|^2||_infinity. The product rule and
|a+b|^2<=2|a|^2+2|b|^2 in (2) give

\[
\boxed{q_I(M_Fv)\leq
  2f_F^2q_I(v)+2C_F\|v\|^2,\qquad v\in\mathcal Q_I.}
\tag{3}
\]

First prove this for v=h_I f, then apply it to differences of core
approximants. Closedness extends the multiplier to the entire form
domain and identifies it with ordinary bounded multiplication. The
constants in (3) depend only on F, not I. This is not obtained by
dropping the interaction or the vacuum subtraction from H_I.

## YM49-T2: smoothing and a uniform time-zero collision estimate

The form version of YM48-T6 is

\[
\|(V_I(t)-I)v\|^2\leq2t q_I(v),\qquad v\in\mathcal Q_I.
\tag{4}
\]

It extends from D(H_I) by the form core and boundedness of V_I.
In addition,

\[
\boxed{q_I(V_I(s)w)\leq \|w\|^2/(2s),\qquad s>0.}
\tag{5}
\]

Here is a proof that checks the domain. Regularize w by
w_lambda=lambda(H_I+lambda)^(-1)w. YM-48's resolvent construction gives
||w_lambda||<=||w||, w_lambda->w, and commutation with V_I.
For v_lambda=V_I(s)w_lambda, the generator derivative and YM47-T3 imply

\[
q_I(v_\lambda)=\lim_{a\downarrow0}
 \frac{\langle w_\lambda,V_I(2s)(I-V_I(a))w_\lambda\rangle}{a}
 \leq\|w\|^2/(2s).
\]

Lower semicontinuity of (1) gives both membership and (5).
No norm-convergent series in the unbounded H_I is used.

Let G be a strict-future history sum with all nonconstant times at least
epsilon>0, and let M_G be its product-sum sup-norm envelope. Its finite
fold has the form Psi_I(G)=V_I(epsilon)w_I, ||w_I||<=M_G; constants can
also be included because V_I(epsilon)h_I=h_I.
For 0<delta<epsilon/2,

\[
\Psi_I(F(\delta)G)
   =V_I(\delta)M_F V_I(\epsilon-\delta)w_I .
\]

Subtract M_F Psi_I(G). Apply (4), (3) and (5) to the first difference,
and YM47-T3 to the remaining time displacement. This proves

\[
\boxed{\|\Psi_I(F(\delta)G)-M_F\Psi_I(G)\|
\leq E_F(\delta;\epsilon,M_G)
:=2M_G\sqrt{\delta(f_F^2/\epsilon+C_F)}
       +2M_G f_F\delta/\epsilon .}
\tag{6}
\]

This is uniform in I and tends to zero. The target on the left lives
in the finite I space; it is not declared to converge in a common
finite-volume norm.

## YM49-T3: bounded time-zero action and its null-space contract

Apply (6) twice, pass squared norms of common histories through YM-47,
and obtain a Cauchy sequence in H_ref. Define

\[
\boxed{M_0(F)[G]=\lim_{\delta\downarrow0}[F(\delta)G].}
\tag{7}
\]

The finite inequality ||M_F Psi_I(G)||<=f_F||Psi_I(G)|| and (6) give

\[
\|M_0(F)[G]\|\leq f_F\|[G]\|.
\tag{8}
\]

Thus (7) preserves the reflected null space, is independent of the
history representative, and extends uniquely to every H_ref.
The convergence in (7) has error at most (6). Limits of finite pairings
also give, for arbitrary strict-future histories G,K,

\[
\langle[K],M_0(F)[G]\rangle
 =\lim_I\langle\Psi_I(K),M_F\Psi_I(G)\rangle
 =\Omega(\overline{rK}\,F(0)G).
\tag{9}
\]

The last expression has separated nonzero history times and the extra
time-zero row, so YM-47 applies. Equations (8)--(9), followed by density,
prove M_0(1)=I and M_0(F)^dagger=M_0(conjugate F).
On the vacuum, (7) is exactly YM-48's embedding: M_0(F)e=JF.

This constructs time-zero multiplication. It does not represent
multiplication by an arbitrary positive-time history. YM-48's null-history
counterexample still forbids that shortcut.

## YM49-T4: multiplication, uniform completion and ordered readouts

The product identity requires a collision argument, not only (9).
Fix G and first eta>0 below its time margin. By (6), finite surrogates
for M_0(F)[A(eta)G] and M_0(FA)[G] differ, after their own boundary
limits, by at most f_F E_A(eta;epsilon,M_G). To make this statement
solely in H_ref, approximate both vectors by their defining positive-time
insertions, pass the finite squared norm to the common-history form,
and then remove the two insertion errors using (6).
Now let eta tend to zero and use (8):

\[
\boxed{M_0(F)M_0(A)=M_0(FA).}
\tag{10}
\]

Linearity, the dagger identity and (8) therefore give a bounded unital
dagger representation. It extends uniquely, in operator norm, from
A_coeff to A_unif; products and daggers pass through uniform limits.
No faithfulness of the entire represented algebra is inferred.

For 0<t_1<...<t_k and F_j in A_coeff, the same finite fold and insertion
bound prove

\[
[F_1(t_1)\cdots F_k(t_k)]
 =T(t_1)M_0(F_1)T(t_2-t_1)\cdots M_0(F_k)e.
\tag{11}
\]

One may prove this by induction from the last insertion. At each step
the target finite vector is the ordered fold, and (6) controls its
replacement before passing history forms. Consequently,

\[
\Omega(F_0(0)F_1(t_1)\cdots F_k(t_k))
 =\langle e,M_0(F_0)T(t_1)M_0(F_1)\cdots M_0(F_k)e\rangle .
\tag{12}
\]

Coincident coefficient times combine using (10). Uniform limits extend
the readout identity to A_unif. The closed coefficient-history sector K
is the cyclic sector generated from e by T and this algebra; it reduces
both, since their adjoints have the same form.

For A,F in A_coeff, (10) gives M_0(F)JA=J(FA). Thus R reduces the
observable algebra and P=JJ^dagger commutes with every M_0(F), including
its uniform completion. Reduction by the time action is a separate test.

## YM49-T5: a positive and exact instantaneous-row closure diagnostic

On H_row^coeff put

\[
C(t)=J^\dagger T(t)J,\quad
Q_R=I-P,\quad L(t)=Q_R T(t)J.
\]

C(t) is a strongly continuous positive self-adjoint contraction, fixes
the row vacuum, and obeys the inherited centered bound exp(-gamma t).
It is not yet called a semigroup. Actual composition of T gives

\[
C(s+t)-C(s)C(t)=J^\dagger T(s)Q_R T(t)J,
\]
\[
\boxed{\Delta(t):=C(2t)-C(t)^2=L(t)^\dagger L(t)\geq0,\qquad
 \|\Delta(t)\|=\|L(t)\|^2.}
\tag{13}
\]

For any row source f, its defect readout equals ||L(t)f||^2. Hence the
following are equivalent:

1. Delta(t)=0 for every t>0.
2. R is invariant, hence reducing, for every T(t).
3. C is a semigroup on H_row^coeff.
4. R=K, the whole coefficient-history sector.

For 1 -> 2 use (13); self-adjointness gives reduction. Reduction proves
3, and 3 implies 1. For 2 -> 4 use (11), the observable invariance of
R and its vacuum; for 4 -> 2 use invariance of K. This does **not**
add the separate assertion K=H_ref for arbitrary bounded histories.
Zero defects on finitely many probes alone do not establish item 1.
The interacting-chain value of (13) remains open here.

## YM49-T6: the exact memory equation with a controlled tail

For any fixed s>0 use the bounded U=T(s), and split H_ref=R direct-sum
Q_R H_ref. Write A=PUP, B=PUQ_R, C=Q_RUP=B^dagger, D=Q_RUQ_R,
as maps between their indicated summands. For z_n=U^n z_0 set
x_n=Pz_n and y_n=Q_Rz_n. RKF N09 applies without unbounded domains:

\[
x_{n+1}=Ax_n+BD^ny_0+
      \sum_{j=0}^{n-1}BD^{\,n-1-j}Cx_j .
\tag{14}
\]

Here e belongs to R, so Q_R H_ref is vacuum-orthogonal. With
r=exp(-gamma s)<1 and ell=||C||=||B||<=r, the inherited gap yields
||D||<=r; also D>=0. Thus each memory coefficient

\[
K_n=BD^nC=C^\dagger D^nC\geq0,\qquad
\|K_n\|\leq\ell^2r^n.
\tag{15}
\]

The initial hidden term has norm at most ell r^n||y_0||. If only memory
ages 0,...,N-1 are retained, the omitted part of (14) is bounded by

\[
\boxed{\frac{\ell^2r^N}{1-r}
 \sup_{0\leq j<n}\|x_j-\langle e,x_j\rangle e\|.}
\tag{16}
\]

The finite sum is bounded by the full geometric tail; B,C kill the
vacuum. Equations (13)--(16) derive the return term from the existing
time action. They do not select independent stochastic noise or identify
it with a physical heat bath. The word "positive" in (15) means a
positive Hilbert quadratic form, not pointwise transition probabilities.

The visible step closes for every hidden initial state exactly when
B=0, equivalently C=0, equivalently Delta(s)=0 at that step.
This is a statement about sampling at s; the continuous closure criterion
above explicitly checks all s. No unproved blocks PHQ_R or Q_RHQ_R of
the unbounded generator are used.

For real z>1, bounded Neumann inversion and block substitution give

\[
P(z-U)^{-1}|_R
 =[z-A-\Sigma(z)]^{-1},\qquad
\Sigma(z)=B(z-D)^{-1}C\geq0 .
\tag{17}
\]

Both inverses exist: z-U is invertible, and solving its two block
equations proves bijectivity of the Schur operator. Moreover
||Sigma(z)||<=ell^2/(z-r). Expanding only the bounded D gives
Sigma(z)=sum_(n>=0) z^(-n-1)K_n, with remainder after N terms at most

\[
\frac{\ell^2}{z}\frac{(r/z)^N}{1-r/z}.
\tag{18}
\]

This is N08/N09's lawful bounded resolvent calculus on the actual carrier.

## YM49-T7: observable and time action fail to commute

Take the same one-site F=x_0=chi_(1/2)/2 as YM48-T8. On this carrier,
||JF||=1/2 and JF is vacuum-orthogonal. Since M_0(F)e=JF,

\[
\boxed{\|[T(t),M_0(F)]e\|
 =\|(T(t)-I)JF\|\geq(1-e^{-\gamma t})/2>0,\quad t>0.}
\tag{19}
\]

The triangle inequality and gap prove the lower bound. The exact energy
3/16 from YM48-T8 also gives the upper bound sqrt(3t/8), by (4) and
the row limit. Thus the commuting time-zero coordinate algebra has a
noncommuting time action, with a nonzero witness on the same carrier.
This persists at theta=0 and is not an interaction/non-Gaussianity test.
Nor is the commutator automatically the native curvature tensor or a
gauge-invariant physical observable. Those identifications require their
own native operator/observation dictionary.

## YM49-T8: evidence, refusal controls and the next gate

The written infinite proof is not mechanically formalized. Exact finite
controls independently check time-zero folds against full reflected-path
sums, complex dagger/product laws, positive form multiplier bounds,
collision errors at rational exponential samples, and both a closing and
a nonclosing observation cut. The latter checks (13), the full memory
recursion, positive memory coefficients, Schur inverses and geometric
tails. These fixtures do not measure the infinite chain's Delta(t).

The four-state fixture reuses YM-48's rational heat samples and uniform
vacuum. Its nonclosing observation is the partition {0,1,2}|{3};
the closing control is {0,1}|{2,3}. These are deliberately different
observers of that finite fixture, not truncations of the chain's entire
instantaneous row. Positive nonconstant S^3 weights h=1+/-x_0/3 test
the weighted energy identity with potential (Lh)/h; they are not claimed
to be chain vacua.

Refusals include future multiplication of a null history, a naive
compressed semigroup, dropping the returned term or its initial source,
a vacuum-only closure test, a false zero memory tail, an omitted Schur
term, a wrong dagger, a missing product cross term and an illegal
time-collision margin. Coefficient/high-content or generator claims
remain at the domains proved by YM-48 and T1--T2.

**Next gate:** evaluate the row-closure defect of the actual interacting
chain with a uniform spatial/content bound. If it vanishes, identify
R=K and its closed row dynamics. If it does not, certify a nonzero
memory source and retain (14). Density of coefficient histories in the
larger bounded-history completion, native Phi_Sigma/NCG intertwining,
physical gauge-observable selection, clock, spatial/relativistic
continuum, AF, Clay and QG remain separate obligations.

~~~bash
python3.12 papers/yang-mills-certified-benchmark/certificates/ym49_time_zero.py --check
python3.12 -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym49_time_zero.py -v
~~~

Default/--check is read-only; --write regenerates only this chapter's
certificate. CI uses a single Python 3.12 job. No prior proof or
certificate is rewritten.
