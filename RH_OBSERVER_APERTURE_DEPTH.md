# APD-1 — Observer aperture and enlarged positive banks

18 October 2026 (Asia/Kolkata). Author: Monty.
Written continuation of ORR-1 at
b0ace198e2af9419f3b59a699403a2b2fc6e411e.

## Result and scope

The complete-source estimate now admits more actual observers on
the same interval. For
\[
N=N_*(a)=\left\lfloor\frac{e^{-4}\sqrt a}{4\pi}\right\rfloor
\ge2048,
\]
every distinct-position bank in \([a/2,3a/2]\) with at most
\(\lfloor N/9\rfloor\) positions has both actual source pencils
strictly positive. The preceding uniform schedule was
\(\lfloor N/12\rfloor\). Thus the guaranteed size coefficient
increases by a factor \(4/3\), with all mixed readings included.

A variable-aperture theorem gives a second schedule:
\(\lfloor N/6\rfloor\) arbitrary positions in
\([2a/3,4a/3]\), for \(N\ge256\).
A direct bound also admits all N known polynomial directions as
distinct observers in an explicitly narrower interval.

ORR's residual recursion and curvature bounds carry to these
enlarged finite ranges. The extension follows from the complete
source estimates below; it is not inferred from the signs of three
scalar derivatives alone.

| Source interval | Arbitrary distinct positions admitted | Range |
|---|---:|---|
| \([a/2,3a/2]\), previous UCB schedule | \(\lfloor N/12\rfloor\) | \(N\ge128\) |
| \([a/2,3a/2]\), this note | \(\lfloor N/9\rfloor\) | \(N\ge2048\) |
| \([2a/3,4a/3]\), this note | \(\lfloor N/6\rfloor\) | \(N\ge256\) |
| \([a(1-h_N),a(1+h_N)]\), this note | \(N\) | \(N\ge128\) |

Here
\[
h_N=\frac{3}{16\ell N^2 49^N},\qquad
\ell=\log(\sqrt a/(2\pi)).
\]
The word aperture denotes only this relative source-position
half-width. It is not an identification of h with native lambda,
cut depth, or the wave coordinate.

Every result is finite at a fixed a. An unlimited compatible bank
that identifies the actual source, and hence RH closure, remains open.
No numerical source evaluation or scientific computation is used.

## 1. Retained actual source and polynomial budget

Use the unchanged source and anchor functional:
\[
\mathcal X(w^2)=\xi(1/2+w),\quad q=\mathcal X'/\mathcal X,
\]
\[
y_a(\mu)=\frac a{a+\mu},\qquad
\mathcal L_a[f]=\frac1{q(a)}
 \sum_\mu\frac{\nu_\mu}{a+\mu}f(y_a(\mu)).
\tag{1}
\]
At the anchors covered here, RBP/UCB retain
\[
q(a)>0,\quad |y_a(\mu)|<1,\quad
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{|a+\mu|}\le2.
\tag{2}
\]
Actual source values \(\mu\) may be nonreal. Products at these
points remain bilinear.

Write \(c_j=\mathcal L_a[y^j]\). For every \(1\le m\le N=N_*(a)\),
the proved polynomial forms obey
\[
A_m=(c_{i+j+1})_{i,j<m},\quad
B_m=(c_{i+j}-c_{i+j+1})_{i,j<m},
\]
\[
\boxed{A_m,B_m\succeq\frac{12}{\ell49^m}I.}
\tag{3}
\]
RBP also supplies an N-node positive probability measure \(\nu_N\)
on \((0,1)\) matching \(c_0,\ldots,c_{2N-1}\).
Every subsequent actual-minus-model moment has absolute value
at most three, by the actual budget two and model budget one.
These are retained finite results, not an all-order positive
representation of the actual source.

Let
\[
I_{a,h}=[a(1-h),a(1+h)],\qquad 0<h<1,
\]
and take distinct \(t_i=a(1+x_i)\) with \(|x_i|\le h\), \(i<m\).
Define their common denominator
\[
D_{\mathbf x}(y)=\prod_{i<m}(1+x_i y)
\]
and the actual coefficient matrices
\[
Q_y=\left(\mathcal L_a\left[
 \frac{y^{i+j+1}}{D_{\mathbf x}(y)^2}\right]\right)_{i,j<m},
\]
\[
Q_c=\left(\mathcal L_a\left[
 \frac{(1-y)y^{i+j}}{D_{\mathbf x}(y)^2}\right]\right)_{i,j<m}.
\tag{4}
\]
Their model versions integrate against \(\nu_N\).

The actual position pencils remain
\[
L_{ij}=-q[t_i,t_j],\qquad J_{ij}=(tq(t))[t_i,t_j].
\]
If C has as its i-th column the coefficients of
\(D_{\mathbf x}(y)/(1+x_i y)\), UCB's algebra gives
\[
\boxed{
L=\frac{q(a)}a C^TQ_yC,\qquad
J=q(a)C^TQ_cC,\qquad
|\det C|=\left|\prod_{i<j}(x_j-x_i)\right|.
}
\tag{5}
\]
Thus coefficient positivity proves both actual mixed signs for
every distinct configuration in the interval. No minimum spacing
is imposed. Duplicate raw positions still represent duplicate rows.

## APD-T1. A variable-aperture, variable-radius source gate

Suppose \(m\le N\), \(0<h<1\), and \(1<R<1/h\).
Set
\[
B(h,R)=\frac{49(1+h)^2R^2}{(1-hR)^2},
\]
\[
\boxed{
G_{m,N}(h,R)=
\frac{\ell m(1+R)}{2R(R-1)}
 B(h,R)^m R^{-2N}.
}
\tag{6}
\]
If \(G_{m,N}(h,R)\le1\), then
\[
\boxed{
Q_y,Q_c\succeq\delta_m(h)I,\qquad
\delta_m(h)=\frac6{\ell[49(1+h)^2]^m}>0.
}
\tag{7}
\]
Consequently \(L,J\succ0\) for every distinct m-position bank
in \(I_{a,h}\). All estimates are uniform over those configurations.

**Proof of the model margin.** On the real model support,
\(0<D_{\mathbf x}(y)\le(1+h)^m\). For
\(p(y)=\sum_{i<m}v_i y^i\),
\[
\int\frac{yp(y)^2}{D_{\mathbf x}(y)^2}\,d\nu_N
\ge(1+h)^{-2m}\int yp(y)^2\,d\nu_N
=(1+h)^{-2m}v^TA_mv.
\]
The equality uses the matched moments through degree \(2m-1\).
The same reasoning with \(1-y\) uses \(B_m\). Hence (3) gives
\[
Q_y^{\rm mod},Q_c^{\rm mod}\succeq2\delta_m(h)I.
\tag{8}
\]

**Proof of the complete-source remainder.** On \(|y|=R\),
\[
|D_{\mathbf x}(y)|\ge(1-hR)^m.
\]
The two entry functions in (4) have maximum moduli bounded by
\[
(1-hR)^{-2m}R^{i+j+1},\qquad
(1+R)(1-hR)^{-2m}R^{i+j}.
\]
Their Taylor series converge on the actual source and model
supports. Cauchy's coefficient estimate, moment agreement below
order \(2N\), and the remaining moment budget three give
\[
\begin{aligned}
\|Q_y-Q_y^{\rm mod}\|
&\le\frac{3m}{R-1}(1-hR)^{-2m}R^{2m-2N},\\
\|Q_c-Q_c^{\rm mod}\|
&\le\frac{3m(1+R)}{R-1}(1-hR)^{-2m}R^{2m-1-2N}.
\end{aligned}
\tag{9}
\]
For example, the first entry tail is bounded by its circle
maximum times
\[
3\sum_{k\ge2N}R^{-k}=\frac{3R^{1-2N}}{R-1}.
\]
Use \(i+j\le2m-2\) and the symmetric row-sum bound to get (9).

The second bound in (9) dominates the first. Its ratio to
\(\delta_m(h)\) is precisely (6).
Thus under the stated gate the entire error has norm at most
half the model margin. Equations (7) and (5) follow. QED.

At \(h=1/2,\ R=3/2\), (6) is exactly UCB's previous
\[
(5/3)\ell m\,3969^m(3/2)^{-2N}\le1.
\]
The radius is an auxiliary coefficient-estimate choice.
Changing it does not change the source, the observer positions,
the normalization anchor, or the measured signs.

## APD-T2. The same interval admits the larger schedule N/9

If \(N\ge2048\), every
\[
\boxed{1\le m\le\lfloor N/9\rfloor}
\tag{10}
\]
satisfies the gate at \(h=1/2,\ R=7/4\).
Therefore every such bank in \([a/2,3a/2]\) is positive, with
the same coefficient margin as UCB:
\[
\boxed{
Q_y,Q_c\succeq
\frac6{\ell(441/4)^m}I.
}
\tag{11}
\]

**Proof.** Substitution in (6) gives
\[
B(1/2,7/4)=147^2,\qquad
G_{m,N}=\frac{22}{21}\ell m\,147^{2m}(7/4)^{-2N}.
\tag{12}
\]
The exact integer inequality
\[
25\cdot147\cdot4^9<24\cdot7^9
\]
implies \(147<(24/25)(7/4)^9\).
The expression \(m147^{2m}\) increases with positive m.
For \(m\le N/9\), it follows that
\[
G_{m,N}
\le\frac{22}{189}\ell N
 \left(\frac{147}{(7/4)^9}\right)^{2N/9}
<\frac{22}{189}\ell N(24/25)^{2N/9}.
\]
Use the retained floor bound \(\ell<N/8\) for \(N\ge128\),
and \(\log(1-1/25)<-1/25\), to obtain
\[
\boxed{
G_{m,N}<\frac{11}{756}N^2e^{-2N/225}.
}
\tag{13}
\]
The right side decreases for \(N>225\).
At \(N=2048\), its polynomial factor is below \(2^{16}\),
and its exponential denominator exceeds \(e^{18}>2^{24}\).
Thus it is below \(2^{-8}<1\).
This proves the gate for every \(N\ge2048\). Equation (7)
then gives (11). QED.

The previous \(N/12\) schedule is retained for \(N\ge128\).
The improvement here is an actual enlargement on the same interval,
not a new interpretation of the old scalar curvature signs.

## APD-T3. A narrower fixed aperture admits arbitrary banks of size N/6

If \(N\ge256\), every
\[
\boxed{m\le\lfloor N/6\rfloor}
\]
distinct-position bank in \([2a/3,4a/3]\) is positive, with
\[
\boxed{
Q_y,Q_c\succeq
\frac6{\ell(784/9)^m}I.
}
\tag{14}
\]

**Proof.** Choose \(h=1/3,\ R=2\). Then
\[
B(1/3,2)=56^2,\qquad
G_{m,N}=\frac34\ell m\,56^{2m}2^{-2N}.
\]
Monotonicity in m gives
\[
G_{m,N}\le\frac{\ell N}{8}(7/8)^{N/3}
<\frac{N^2}{64}e^{-N/24}.
\tag{15}
\]
This bound decreases for \(N>48\).
At \(N=256\), the polynomial factor is \(2^{10}\);
the exponential denominator exceeds \(e^9>2^{12}\).
The bound is below \(1/4\), so T1 applies. QED.

This schedule includes arbitrary configurations and both sides of
the anchor. MOP's separate dyadic \(N/6\) result remains available
on its own observation set; the two domains are not identified.

## APD-T4. A direct source estimate recovers the full finite polynomial budget

For a given \(m\le N\), define
\[
\boxed{h_m=\frac{3}{16\ell m^2 49^m}.}
\tag{16}
\]
For every \(0<h\le h_m\) and every distinct m-position bank
in \(I_{a,h}\),
\[
\boxed{
Q_y,Q_c\succeq\widetilde\delta_m I,\qquad
\widetilde\delta_m=\frac6{\ell49^m}>0.
}
\tag{17}
\]
In particular a single interval \(I_{a,h_N}\) admits every
distinct bank of size at most N.

**Proof.** This estimate compares directly with the proved actual
polynomial forms (3). It does not require an omitted-moment
model comparison.

For \(|y|\le1\), the absolutely convergent product series gives
\[
\left|D_{\mathbf x}(y)^{-2}-1\right|
\le(1-h)^{-2m}-1=:E_D.
\tag{18}
\]
Indeed the absolute coefficients are majorized by those of
\((1-h|y|)^{-2m}\); their constant coefficient is one.
Using the actual absolute source budget two, \(|y_a|<1\),
and \(|1-y_a|\le2\), each entry difference satisfies
\[
|(Q_y-A_m)_{ij}|\le2E_D,\qquad
|(Q_c-B_m)_{ij}|\le4E_D.
\]
Both matrices are real symmetric by conjugation symmetry.
Thus
\[
\boxed{
\|Q_y-A_m\|\le2mE_D,\qquad
\|Q_c-B_m\|\le4mE_D.
}
\tag{19}
\]
No order comparison at the complex source points was used.

Here \(\ell>1\), so (16) ensures \(h\le1/(4m)\le1/2\).
The elementary estimates
\[
-\log(1-h)\le2h,\qquad e^u-1\le2u\quad(0\le u\le1)
\]
give
\[
E_D\le e^{4mh}-1\le8mh.
\]
Therefore
\[
4mE_D\le32m^2h
\le\frac6{\ell49^m}.
\tag{20}
\]
Subtracting (20) from the minimum eigenvalue in (3) proves (17).
Equation (5) then proves the two actual position signs.

Finally \(h_N\le h_m\) when \(m\le N\), so the same interval
\(I_{a,h_N}\) works simultaneously for every covered size. QED.

The full N in this theorem is the already known finite polynomial
budget at that anchor. The explicit width can be extremely small.
The result admits distinct positions without a minimum spacing;
it does not assert a spacing-independent raw matrix eigenvalue gap.

## APD-T5. The enlarged ranges retain the residual and mixed-curvature margins

Use any one of the schedules above. For total bank size n let
\(\Delta_n\) be its proved coefficient margin:
\[
\begin{array}{c|c}
\text{Schedule}&\Delta_n\\ \hline
\text{T2, }h=1/2&6/[\ell(441/4)^n]\\
\text{T3, }h=1/3&6/[\ell(784/9)^n]\\
\text{T4, fixed }h\le h_N&6/[\ell49^n].
\end{array}
\tag{21}
\]
Let \(A=\{t_1,\ldots,t_m\}\) lie in the selected interval.
If m+1 is covered, a new position \(\tau\notin A\) has the
actual Schur margins
\[
\boxed{
\begin{aligned}
\alpha_A(\tau)&\ge
\frac{q(a)}a\Delta_{m+1}
\frac{\prod_i((\tau-t_i)/a)^2}
     {S_{m+1}(\tau/a-1)},\\
\beta_A(\tau)&\ge
q(a)\Delta_{m+1}
\frac{\prod_i((\tau-t_i)/a)^2}
     {S_{m+1}(\tau/a-1)},
\end{aligned}
}
\tag{22}
\]
where \(S_n(v)=\sum_{j=0}^{n-1}v^{2j}\).
These hold at every new position in the whole selected interval.

**Proof.** The common numerator argument in UCB-T3 uses only
the lower coefficient bound and the polynomial identity in (5).
Replacing its coefficient margin by the already proved
\(\Delta_{m+1}\) gives (22) without any other change. QED.

In particular all mixed Schur losses in ORR-T3 are controlled
whenever the complete union is within the chosen size budget.
This is a statement about the whole union, not merely its
separately positive subbanks.

For explicit curvature bounds suppose m+2 is covered.
Retain ORR's functions for this bank:
\[
q-q_A=H_A^2\mathcal T_A,\qquad
V=1/\mathcal T_A,\quad F=-V/\tau,
\]
\[
b_i(\tau)=\frac{\prod_{k\ne i}(\tau-t_k)}{d_A(\tau)},\qquad
c_A(\tau)=b(\tau)^T(J+\tau L)^{-1}b(\tau).
\]
Define
\[
\boxed{
\mathcal E_A(\tau)=
\frac{q(a)\Delta_{m+2}d_A(\tau)^2}
 {a^{2m+3}S_{m+2}(\tau/a-1)\mathcal T_A(\tau)^2}>0.
}
\tag{23}
\]
Then throughout the chosen interval, including fitted positions,
\[
\boxed{-\tfrac12V''-c_A\ge(a+\tau)\mathcal E_A>0,}
\tag{24}
\]
\[
\boxed{
\frac16V'''-(b')^TL^{-1}b'
-\frac{(\frac12V''-b^TL^{-1}b')^2}
 {V'-b^TL^{-1}b}
\ge\mathcal E_A>0,
}
\tag{25}
\]
\[
\boxed{
\frac16F'''-(b')^TJ^{-1}b'
-\frac{(\frac12F''-b^TJ^{-1}b')^2}
 {F'-b^TJ^{-1}b}
\ge\frac a{\tau^2}\mathcal E_A>0.
}
\tag{26}
\]
Primes here are source-position derivatives with A fixed.
In particular
\[
\boxed{V'>0,\qquad V''<0,\qquad V'''>0}
\tag{27}
\]
now hold on these expanded declared ranges.

**Proof.** Divide (22) by the known interpolation factor \(H_A^2\).
DSG's exact identities then give the positive normalized margins,
positivity of \(\mathcal T_A\), and the inverse-residual slope
window. Apply the same argument to \(B=A\cup\{\tau\}\), leaving
the second extra position available.

ORR-T1 supplies the unchanged exact values
\[
d_B(\tau)=d_A(\tau)\mathcal T_A(\tau)/\mathfrak a_A(\tau),
\quad
H_B'(\tau)=\mathfrak a_A(\tau)H_A(\tau)/\mathcal T_A(\tau),
\]
and its fitted residual identity. Insert the new margin
\[
\varepsilon_L^B(\tau)=
\frac{q(a)\Delta_{m+2}d_B(\tau)^2}
 {a^{2m+3}S_{m+2}(\tau/a-1)},\qquad
\varepsilon_J^B=a\varepsilon_L^B
\]
into ORR-T4--T5. The cancellations are identical and give
(23)--(26). Strict diagonal margins make all denominators positive.
Analytic continuation includes the old fitted positions.
The sign conclusions follow as in ORR. QED.

For the same-interval result, the former restriction
\(m+2\le\lfloor N/12\rfloor\) is thus replaced by
\(m+2\le\lfloor N/9\rfloor\) when \(N\ge2048\).

## What the enlargement does and does not settle

The extra finite range is proved from the unchanged complete
source, with all omitted moments bounded. The positive finite
model is used only where its exact moment matching is available.
The narrow-interval theorem instead uses the actual polynomial
forms and the actual absolute source bound directly.

The proved polynomial budget at a fixed anchor remains finite, N.
Shrinking the interval
does not turn (3) into an all-order assertion. Increasing a gives
larger capacities, but banks chosen at different anchors must still
have their mixed signs proved if they are to form one compatible
sequence retaining the old observations. Separately covered banks
do not automatically supply that union.

The sufficient NOE/DSG route to actual-source recovery remains:
one nested positive sequence, with positions tending to infinity
and divergent sum of reciprocal square roots. Neither the improved
constant nor the shrinking-window construction establishes that
sequence. RH is not proved.

## Sources, attribution, and verification

- [RBP-1](RH_RELATIONAL_BANK_POSITIVITY.md), T1--T3:
  complete wave estimate, finite polynomial positivity, and the
  positive carrier matching the first 2N source moments.
- [UCB-1](RH_OBSERVER_BANK_UNION.md), T1--T3:
  common denominator, actual position congruence, previous
  N/12 schedule, and spacing-normalized source innovations.
- [ORR-1](RH_OBSERVER_RESIDUAL_RECURSION.md), T1--T5:
  exact residual recursion, mixed signs and confluent curvature.
- [DSG-1](RH_DENSE_OBSERVER_SLOPE_GATE.md), T2--T5:
  inverse-residual slope window and source-identification coverage.
- [NOE-1](RH_NESTED_OBSERVER_EXTENSION.md):
  compatible infinite model limits and the remaining source-recovery
  condition.

Cauchy coefficient estimates, perturbation bounds, rational
interpolation and Schur complements are standard tools.
The additions here are the variable-aperture source bound, the
larger explicit finite schedules, and the transfer of the retained
curvature margins to those ranges. No literature-wide novelty
or optimality of the constants is claimed.

Verification consists of written algebra and inequalities, plus
text/readback and repository metadata checks. No scientific
computation, numerical source evaluation, independent agent,
certificate replay or workflow dispatch is performed.
Only this note and its README entry change. The commit uses
[skip ci]. Published as an unmerged draft stacked on #134.
