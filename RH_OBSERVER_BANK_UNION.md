# UCB-1 — Compatible unions of actual observer banks

17 October 2026. Written continuation of MOP-1 at
22833267656c2df2560968bd9e280ee9ba406b00.

MOP supplies one common positive finite spectrum on a growing
dyadic bank. Its final qualification is that separately positive
banks do not establish the signs of their mixed readings. This
note supplies those signs on an explicit finite range: every set
of distinct observer positions in a declared interval is admitted,
provided the total number of positions is within the bound below.
The set may be a union of banks chosen at different local anchors.

The new estimates use the whole actual source in a common rational
basis. They give both actual Schur innovations, an exact finite
compression system, strict interlacing under one-position extension,
and a positive error between a finite response and the actual
response throughout the interval when one extra position is allowed.

The established interval and total-size restriction remain essential.
An unlimited compatible positive family with sufficient sampling
density, and hence RH, is not proved.
No scientific computation or source replay is used.

## 1. The retained source and the new scope

Keep the actual completed source
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,
\]
\[
y_a(\mu)=\frac a{a+\mu},\qquad
\mathcal L_a[f]=
\frac1{q(a)}\sum_\mu\frac{\nu_\mu}{a+\mu}f(y_a(\mu)).
\tag{1}
\]
The source spectral coordinate is
\(\mu=(\gamma-i(\beta-1/2))^2\). No reality of every actual
\(\mu\), or positive infinite source measure, is assumed.

For \(a\ge a_0=3249/4\), the retained source has
\(q(a)>0,\ |y_a(\mu)|<1\), total absolute normalized weight at
most two, and real moments \(c_k=\mathcal L_a[y^k]\), \(c_0=1\).
Products in this functional remain bilinear at the actual,
possibly complex, source points.

Write
\[
b=\sqrt a,\qquad \ell=\log\frac b{2\pi},\qquad
N_*(a)=\left\lfloor\frac{e^{-4}b}{4\pi}\right\rfloor,\qquad
I_a=[a/2,3a/2].
\tag{2}
\]
For \(b\ge64\) and \(m\le N\le N_*(a)\), RBP and MOP retain
\[
A_m=(c_{i+j+1}),\qquad B_m=(c_{i+j}-c_{i+j+1}),
\]
\[
A_m\succeq\frac1{2\ell}A_m^*,\qquad
B_m\succeq\frac1{2\ell}B_m^*,\qquad
A_m^*,B_m^*\succeq24\,49^{-m}I.
\tag{3}
\]
The stars use the arcsine reference
\(d\nu_*=dy/(\pi\sqrt{y(1-y)})\).
RBP also constructs an \(N\)-node positive probability measure
\(\nu_N\) on \((0,1)\) matching \(c_0,\ldots,c_{2N-1}\).
Its later moments can differ from the actual source moments.

Choose any distinct real positions
\[
t_i=a(1+x_i)\in I_a,\qquad -1/2\le x_i\le1/2,\qquad i<m.
\]
No dyadic pattern, minimum spacing, integer spacing, or prime
selection is imposed. Define the actual kernels
\[
K_y(x,z)=\mathcal L_a\!\left[
\frac{y}{(1+xy)(1+zy)}\right],\qquad
K_c(x,z)=\mathcal L_a\!\left[
\frac{1-y}{(1+xy)(1+zy)}\right].
\tag{4}
\]
With the actual pencils
\[
L_{ij}=-q[t_i,t_j],\qquad J_{ij}=(tq(t))[t_i,t_j],
\]
\[
L_{ii}=-q'(t_i),\qquad J_{ii}=q(t_i)+t_iq'(t_i),
\]
the exact identities are
\[
L=\frac{q(a)}a K_y[\mathbf x],\qquad J=q(a)K_c[\mathbf x].
\tag{5}
\]
MOP proves these identities for nonnegative displacements.
The same algebra applies here because \(|x_i y_a|<1\);
all denominators remain nonzero.

The value \(a\) in this estimate is a common comparison anchor.
The matrices \(L,J\) themselves use actual \(q,q'\) at the
selected physical positions and do not depend on a label
previously assigned to a smaller bank's local anchor.

## UCB-T1. A common rational basis controls every mixed reading

Put \(R=3/2\) and
\[
D_{\mathbf x}(y)=\prod_{i<m}(1+x_i y),\qquad
Q_y=\left(\mathcal L_a\!\left[
\frac{y^{i+j+1}}{D_{\mathbf x}(y)^2}\right]\right)_{i,j<m},
\]
\[
Q_c=\left(\mathcal L_a\!\left[
\frac{(1-y)y^{i+j}}{D_{\mathbf x}(y)^2}\right]\right)_{i,j<m}.
\tag{6}
\]
Let \(Q_y^{\rm mod},Q_c^{\rm mod}\) replace \(\mathcal L_a\)
by integration against the acquired \(\nu_N\).
Then, uniformly over the positions just specified,
\[
\boxed{
Q_y^{\rm mod},Q_c^{\rm mod}\succeq
\frac{12}{\ell(441/4)^m}I,
}
\tag{7}
\]
\[
\boxed{
\begin{aligned}
\|Q_y-Q_y^{\rm mod}\|
&\le9m\,16^m R^{2m-1-2N},\\
\|Q_c-Q_c^{\rm mod}\|
&\le\frac{45}{2}m\,16^m R^{2m-2-2N}.
\end{aligned}}
\tag{8}
\]
In particular the sufficient gate
\[
\boxed{\frac53\,\ell m\,3969^mR^{-2N}\le1}
\tag{9}
\]
implies
\[
\boxed{
Q_y,Q_c\succeq\delta_m(a)I,\qquad
\delta_m(a)=\frac6{\ell(441/4)^m}>0.
}
\tag{10}
\]

**Model margin.** On the real model support, every factor of \(D\)
is positive and \(D(y)\le R^m\). If \(p(y)=\sum_{i<m}v_i y^i\),
the model \(y\)-form is bounded below by
\[
R^{-2m}\int y p(y)^2\,d\nu_N
=R^{-2m}v^TA_m v
\ge\frac{12}{\ell49^mR^{2m}}\|v\|^2.
\]
The equality uses degree at most \(2m-1\le2N-1\).
Replace \(y\) by \(1-y\) and use \(B_m\) for the second form.
Since \(49R^2=441/4\), this proves (7).
These order comparisons occur only on the positive real model
and on the already proved polynomial forms.

**Complete-source remainder.** On \(|y|=R\),
\[
|D(y)|\ge(1-R/2)^m=4^{-m}.
\]
The two entry functions in (6) therefore have modulus at most
\(16^mR^{i+j+1}\) and
\((1+R)16^mR^{i+j}\), respectively.
Cauchy's coefficient estimate applies because their poles lie
outside this circle. The actual and model moments agree below
order \(2N\); every later moment difference has modulus at most
three. Consequently the first entry error is at most
\[
3\,16^m\sum_{k\ge2N}R^{i+j+1-k}
=9\,16^mR^{i+j+1-2N}.
\]
The second is at most
\(9(1+R)16^mR^{i+j-2N}\).
Use \(i+j\le2m-2\) and the symmetric row-sum norm bound to
obtain (8). The second displayed matrix bound dominates the first.

For it to be at most half the model margin, it suffices that
\[
\frac{15}{4}\ell m(49\cdot16)^mR^{4m-2-2N}\le1.
\]
Here \(49\cdot16\cdot R^4=3969\) and
\((15/4)R^{-2}=5/3\), giving (9) and (10). QED.

The common denominator retains all cross terms. The argument
does not assemble a positive matrix from its separate diagonal
blocks or substitute absolute squares at complex source points.

## UCB-T2. Every finite union in the interval is positive within an explicit size budget

Take \(N=N_*(a)\ge128\) and set \(M(a)=\lfloor N/12\rfloor\).
Then every \(1\le m\le M(a)\) satisfies (9). Consequently,
for every \(m\)-position bank in \(I_a\),
\[
\boxed{L_m\succ0,\qquad J_m\succ0.}
\tag{11}
\]
This applies to an arbitrary union of banks after repeated
positions have been identified and its total distinct size counted.

**Proof of the schedule.** The retained floor estimate gives
\[
\ell<4+\log(2(N+1))<N/8\qquad(N\ge128).
\]
Also \(3969<R^{21}\): \(R^8>25,\ R^4>5\), so
\(R^{21}>25^2\cdot5\cdot(3/2)=9375/2>3969\).
For \(m\le N/12\), the left side of (9) is therefore less than
\[
\frac53\ell mR^{21m-2N}
\le\frac5{288}N^2R^{-N/4}.
\]
This decreases for \(N\ge128\), since \(\log R>1/3\).
At \(N=128\), the numerator is below \(2^9\) and
\(R^{32}>2^{16}\), so the expression is below \(2^{-7}<1\).

**Passage to the actual position matrices.** Let \(C_{\mathbf x}\)
have as its \(i\)-th column the coefficients of
\[
p_i(y)=D_{\mathbf x}(y)/(1+x_i y).
\]
The polynomials have degree at most \(m-1\), and
\[
K_y=C_{\mathbf x}^TQ_y C_{\mathbf x},\qquad
K_c=C_{\mathbf x}^TQ_c C_{\mathbf x}.
\tag{12}
\]
For distinct positions the coefficient matrix is invertible, with
\[
|\det C_{\mathbf x}|=
|V(\mathbf x)|,\qquad
V(\mathbf x)=\prod_{i<j}(x_j-x_i).
\tag{13}
\]
For nonzero \(x_i\), evaluate a linear combination of the \(p_i\)
at \(y=-1/x_j\) to establish independence. Comparing the resulting
Vandermonde determinants gives (13); polynomial continuity includes
a position \(x_i=0\).
Thus (10), (12) and (5) prove (11). They also give
\[
\boxed{
\frac{\det K_y}{V(\mathbf x)^2},
\ \frac{\det K_c}{V(\mathbf x)^2}
\ge\delta_m(a)^m.
}
\tag{14}
\]
QED.

Equations (10) and (14) use declared coefficient and spacing
normalizations. They do not supply a spacing-independent lower
bound on the eigenvalues of the raw position matrices.

The new size coefficient is deliberately compared on its own
scope. MOP retains the larger \(N/6\) schedule for its specified
dyadic bank at \(N\ge2048\). Here \(N/12\) admits EVERY distinct
configuration in the wider interval \(I_a\), including mixed
banks and observers on either side of \(a\). It is not an
improvement of MOP's numerical dyadic coefficient.

## UCB-T3. Both actual innovations have explicit spacing-normalized margins

Suppose an old bank contains \(x_0,\ldots,x_{m-2}\), and add
a distinct \(x=x_{m-1}\), with the full union covered by T1
or T2. Put
\[
\tau=a(1+x),\qquad r=q(\tau),\qquad d=q'(\tau),
\]
\[
u_i=q(t_i),\qquad
l_i=\frac{u_i-r}{\tau-t_i},\qquad
j_i=\frac{\tau r-t_i u_i}{\tau-t_i}.
\]
The actual SE49 innovations obey
\[
\boxed{
\begin{aligned}
\alpha&=-d-l^TL_{\rm old}^{-1}l
 \ge\frac{q(a)}a\,\delta_m(a)
 \frac{\prod_{i<m-1}(x-x_i)^2}{S_m(x)},\\
\beta&=r+\tau d-j^TJ_{\rm old}^{-1}j
 \ge q(a)\,\delta_m(a)
 \frac{\prod_{i<m-1}(x-x_i)^2}{S_m(x)},
\end{aligned}}
\tag{15}
\]
where \(S_m(x)=\sum_{k=0}^{m-1}x^{2k}\le4/3\).
Both signs are strictly positive. The formulas for the first
position use empty products equal to one and omit the old inverses.

**Proof.** The Schur complement of a real positive matrix is the
minimum of its quadratic form over old coefficients, with the new
coefficient fixed to one. For either kernel, (12) makes the
corresponding common numerator
\[
p(y)=p_{m-1}(y)+\sum_{i<m-1}v_i p_i(y).
\]
If \(x\ne0\), every old \(p_i\) vanishes at \(y=-1/x\), so
\[
x^{m-1}p(-1/x)=\prod_{i<m-1}(x-x_i).
\]
Writing \(p(y)=\sum_{k<m}p_k y^k\), Cauchy--Schwarz gives
\[
\left|\prod_{i<m-1}(x-x_i)\right|
\le\left(\sum_{k<m}|p_k|^2\right)^{1/2}
   \left(\sum_{k<m}|x|^{2(m-1-k)}\right)^{1/2}.
\]
At \(x=0\) the same inequality follows from the leading
coefficient of \(p\): the old numerators have degree at most
\(m-2\), while the new one has leading coefficient
\(\prod_{i<m-1}x_i\).
Apply (10), minimize over \(v\), and insert the two scale
factors from (5). This proves (15). QED.

The derivative window from SE49 now has actual proved margins:
\[
\frac{j^TJ_{\rm old}^{-1}j-r}{\tau}
<q'(\tau)<-l^TL_{\rm old}^{-1}l.
\tag{16}
\]
The established premise is the whole union's size and location
gate; positivity of an old bank alone is still insufficient.
Small spacing makes raw innovations small, as (15) records.
It does not change their established signs.

## UCB-T4. Finite unions form an exact compression system

Let \(U\) be any covered union. For each nonempty subbank
\(A\subseteq U\), use its actual principal pencils and source
column \(r_A=(q(t_i))_{i\in A}\). SE48 gives
\[
q_A(z)=r_A^T(J_A+zL_A)^{-1}r_A
      =\sum_{j=1}^{|A|}\frac{R_j^A}{z+\mu_j^A},
\qquad R_j^A,\mu_j^A>0,
\tag{17}
\]
with distinct model poles and exact value/derivative matching
at every position of \(A\).

Define the symmetric carrier and its source vector
\[
\widehat A_A=L_A^{-1/2}J_A L_A^{-1/2},\qquad
v_A=L_A^{-1/2}r_A.
\]
For \(A\subseteq B\subseteq U\), let \(I_{B,A}\) insert the old
position coordinates into the larger bank and set
\[
E_{B,A}=L_B^{1/2}I_{B,A}L_A^{-1/2}.
\]
Then
\[
\boxed{
E_{B,A}^TE_{B,A}=I,\qquad
E_{B,A}^T\widehat A_B E_{B,A}=\widehat A_A,\qquad
E_{B,A}^Tv_B=v_A.
}
\tag{18}
\]
For \(A\subseteq B\subseteq C\subseteq U\),
\[
\boxed{E_{C,B}E_{B,A}=E_{C,A}.}
\tag{19}
\]

**Proof.** The actual pencils and source columns restrict exactly:
\(I_{B,A}^TL_B I_{B,A}=L_A\), and likewise for \(J,r\).
Substitute these identities into (18). The middle square roots
cancel in (19), leaving the composite coordinate insertion.
SE48 applies because T2 proves both signs at every subbank.
This proves the assertions. QED.

There is also one explicit positive spectral realization of ALL
these subbank pairings. Write
\(q_U(z)=\sum_{k=1}^{|U|}R_k/(z+\mu_k)\), put
\(\Lambda_U=\operatorname{diag}(\mu_k)\), and define
\[
(\Phi_A)_{ki}=\frac{\sqrt{R_k}}{t_i+\mu_k}\qquad(i\in A).
\]
Exact value and derivative interpolation gives
\[
\boxed{
L_A=\Phi_A^T\Phi_A,\qquad
J_A=\Phi_A^T\Lambda_U\Phi_A,\qquad
r_A=\Phi_A^T(\sqrt{R_k})_k.
}
\tag{20}
\]
The off-diagonal identities follow by divided differences,
and the diagonal identities by the derivative matches.

Thus one finite positive carrier realizes the actual mixed
readings of a covered union. Its subbank operators are
compressions; their spectra need not coincide with its spectrum.
Compressing a resolvent is not generally the resolvent of the
compressed operator. Equation (18) does not assert exact operator
intertwining or equality of the full responses outside their
interpolation positions.

Within the fixed union model, MOP's observer probabilities are
\[
\pi_k^U(t)=\frac{R_k}{(t+\mu_k)q_U(t)},\qquad \sum_k\pi_k^U(t)=1.
\tag{21}
\]
They obey its exact anchor-transport law. Equation (20) supplies
the additional compatibility of finite subbanks; it does not
identify their individual probability atoms across changing depth.

## UCB-T5. One-position growth gives strict interlacing

Suppose \(B=A\cup\{\tau\}\) is covered, with \(|A|=m\).
Order both model spectra increasingly. Then
\[
\boxed{
\mu_j^B<\mu_j^A<\mu_{j+1}^B,\qquad 1\le j\le m.
}
\tag{22}
\]
In particular, adding an actual independent position changes the
finite spectrum in a controlled way; its old eigenvalues are not
required to remain fixed.

**Proof.** Realize the old pairing inside the \((m+1)\)-mode
model of \(B\) by (20). Every \(m\)-by-\(m\) row minor of
\(\Phi_A\) is nonzero: after the positive row factors
\(\sqrt{R_k}\) are removed, it is a Cauchy determinant with
distinct positive \(\mu_k^B\) and distinct positive \(t_i\).
Its determinant is a nonzero product of the two Vandermonde
factors divided by \(\prod_{k,i}(t_i+\mu_k^B)\).

Consequently the one-dimensional orthogonal complement of
\(\operatorname{range}\Phi_A\) has a unit vector \(n\) with
every component nonzero, as follows by the cofactor formula.
The old symmetric carrier is the compression of
\(\Lambda_B=\operatorname{diag}(\mu_k^B)\) to \(n^\perp\).
Its eigenvalues solve
\[
\sum_{k=1}^{m+1}\frac{n_k^2}{\mu_k^B-\lambda}=0.
\tag{23}
\]
No diagonal value \(\mu_k^B\) is a compression eigenvalue,
because \(n_k\ne0\). Between consecutive diagonal values the
left side of (23) is strictly increasing from minus infinity
to plus infinity. It has exactly one root in each interval.
These are all \(m\) compression eigenvalues, proving (22). QED.

The same source inclusion also gives the retained SE49 response
update. With
\[
P(z)=J_A+zL_A,\qquad k(z)=j+zl,\qquad
s(z)=r+(\tau-z)d-k(z)^TP(z)^{-1}k(z),
\]
\[
e(z)=r-k(z)^TP(z)^{-1}r_A,
\]
one has, for real \(z>0\),
\[
\boxed{
q_B(z)-q_A(z)=\frac{e(z)^2}{s(z)}\ge0,\qquad
s(z)\ge\beta+z\alpha>0.
}
\tag{24}
\]
This is SE49's established identity, now applied where T3
supplies the actual innovations. It is not introduced as a
new Schur-complement theorem.

## UCB-T6. One extra admissible position proves an actual lower envelope

Suppose \(N_*(a)\ge128\), let \(A\subset I_a\) have \(m\ge1\)
distinct positions, and require
\[
m+1\le M(a)=\lfloor N_*(a)/12\rfloor.
\]
Then the same finite model is below the actual source throughout
the interval:
\[
\boxed{q_A(t)\le q(t)\qquad(t\in I_a),}
\tag{25}
\]
with equality at the selected positions and strict inequality
at every other position in \(I_a\). For \(t=a(1+x)\notin A\),
the explicit lower error is
\[
\boxed{
q(t)-q_A(t)\ge
q(a)(1+t/a)\,\delta_{m+1}(a)\,
\frac{\prod_{i\in A}(x-x_i)^2}{S_{m+1}(x)}>0.
}
\tag{26}
\]

**Proof.** Add the actual observation at this \(t\). The enlarged
union is covered by T2. At the response test point \(z=t\),
the mixed column of its pencil \(J+zL\) is exactly \(r_A\),
and its new diagonal entry is \(q(t)\); the derivative terms
cancel. Thus its Schur complement is
\[
s(t)=q(t)-r_A^T(J_A+tL_A)^{-1}r_A=q(t)-q_A(t).
\]
Equation (24) bounds it below by \(\beta+t\alpha\).
Apply the two bounds (15) with full size \(m+1\) to obtain
(26). At an old position, SE48 gives exact equality.
This proves (25) and the strict statement. QED.

This proves the sign of the actual error, not only monotonicity
of finite approximants. The extra position is used universally
in the written theorem; no fresh numerical observation is
evaluated at an interval of points.

For a fully explicit upper error bound, choose any \(t_*\in A\)
and write
\[
w_*=\sqrt{t_*},\qquad
C_A=q(t_*)+\frac{41}{8w_*^2}+\frac1{8w_*}.
\tag{27}
\]
For every covered bank, even without the extra-position budget,
the SE49/MOP analytic argument gives
\[
\boxed{
|q(w^2)-q_A(w^2)|
\le C_A\frac{|w+w_*|^2}{|w|}
\prod_{i\in A}
\left|\frac{w-\sqrt{t_i}}{w+\sqrt{t_i}-2}\right|^2,
\qquad \Re w>1.
}
\tag{28}
\]
On the extra-position range, combine this with (25) to get
\[
\boxed{
0\le q(t)-q_A(t)\le9\sqrt a\, C_A\,4^{-m}
\qquad(t\in I_a).
}
\tag{29}
\]

**Proof of the upper bound.** The full actual source satisfies
\(|wq(w^2)|<41/8+|w|/4\) on \(\Re w>1\), as proved in SE49.
The finite model \(F_A(w)=wq_A(w^2)\) has positive real part
on \(\Re w>0\), and \(F_A(w_*)=w_*q(t_*)>0\).
The half-plane Schwarz argument gives
\[
|F_A(w)|\le q(t_*)\frac{|w+w_*|^2}{\Re w}.
\]
Indeed its initial bound is \(c(1+\rho)/(1-\rho)\), with
\(c=w_*q(t_*)\), \(\rho=|w-w_*|/|w+w_*|\); use
\(|w-w_*|\le|w+w_*|\) to obtain the displayed expression.

Also \(|w+w_*|^2\ge|w|^2+w_*^2\), so
\[
\left|\frac{wq(w^2)}{(w+w_*)^2}\right|
\le\frac{41}{8w_*^2}+\frac1{8w_*}.
\]
Thus \((wq(w^2)-F_A(w))/(C_A(w+w_*)^2)\) is holomorphic
and bounded by one on \(\Re w>1\). It has a double zero at
every \(\sqrt{t_i}\). Successive Schwarz division proves (28),
with the full actual source and its tails retained.

For real \(w\) and every \(\sqrt{t_i}\) in
\([b/\sqrt2,\sqrt{3/2}\,b]\), use
\[
2b/3<b/\sqrt2,\qquad \sqrt{3/2}\,b<5b/4.
\]
The numerator of each unsquared factor is less than \(7b/12\),
while its denominator is at least \(4b/3-2\ge7b/6\),
since \(b\ge64\). Every such factor is at most \(1/2\).
The prefactor obeys
\[
\frac{(w+w_*)^2}{w}
\le \frac{4(3a/2)}{b/\sqrt2}
=6\sqrt2\,b<9b.
\]
This proves (29), using (25) for its lower side. QED.

For \(A,B\subseteq U\) with \(|U|+1\le M(a)\), the same
conclusions give
\[
\boxed{
q_A(t),q_B(t)\le q_U(t)\le q(t)\qquad(t\in I_a).
}
\tag{30}
\]
The first two inequalities follow by repeated finite updates
(24); the last uses one extra admissible position. This is a
common refinement of the actual finite response envelopes.
It does not equate individual eigenvalues in different models.

## 2. What the finite compatibility establishes

| Question | Established statement |
|---|---|
| Can independently selected finite banks be joined? | Yes, when their distinct union lies in \(I_a\) and has total size at most \(M(a)\) |
| Are the mixed actual source signs proved? | Both pencils and every sequential innovation have bounds (10), (15) |
| Is there one positive carrier for a covered union? | Yes; (20) realizes every subbank's actual pairings |
| Does compatibility compose through intermediate subbanks? | Yes; exact isometries and compressions satisfy (18), (19) |
| Must the eigenvalues remain the same as depth grows? | No; a one-position extension gives strict interlacing (22) |
| Is the finite response below the actual source? | Throughout \(I_a\), when one extra position remains in the budget |
| Does this identify all model poles with xi zeros? | No |
| Does this prove unlimited compatible positivity or RH? | No |

The construction supplies mixed signs directly, rather than
assuming that overlapping positive banks already imply a positive
union. A previous bank's local anchor is just a label when its
actual \(L,J,r\) data are carried into the common comparison.

The asymptotic size is
\[
M(a)\sim\frac{e^{-4}}{48\pi}\sqrt a.
\tag{31}
\]
At fixed \(a\), this is a finite budget. As \(a\) grows, the
whole interval \(I_a\) moves outward. A fixed old positive
position eventually lies outside it. Repeatedly joining covered
finite unions therefore does not establish every prefix of one
fixed infinite observation bank. Neither overlapping intervals
nor the exact finite compression identities remove this issue.

The next global obligation remains actual positivity on an
unlimited compatible bank with sufficient analytic sampling
density, or an equivalent all-order source condition. POC/SE49
already prove the implication to RH once those universal signs
hold. UCB supplies a new unconditional finite region of those
mixed signs; it does not assume their infinite extension.

All eigenvalues in (17), (20)--(23) belong to finite interpolating
models. The actual source coordinate \(\mu\) in (1) remains
distinct. In particular, positivity of model poles or their
interlacing cannot replace actual-zero certification.

CCI's full source negative-index count has not been shown to be
zero. BRR's unrestricted width requirement and the previously
declared full/restricted support-window results are not enlarged
by this note.

## References, attribution, and verification scope

- [MOP-1](RH_MULTI_OBSERVER_PROBABILITY.md): both dyadic pencils, common finite spectral response, normalized probability transport, and moving-interval limitations.
- [RBP-1](RH_RELATIONAL_BANK_POSITIVITY.md), T2--T3: actual reference-form comparisons and the positive moment-matching carrier.
- [RLO-1](RH_RELATIONAL_OBSERVER_CONFLUENCE.md): divided readings, lawful observer changes, and spacing-normalized determinants.
- [SE48](papers/native-prime-observation-closure/sources/SE48.md), theorem 1: positive finite spectral realization and value/derivative interpolation.
- [SE49](papers/native-moving-zero-covariance/sources/SE49.md), theorems 1--3: both innovation signs, finite response update, full-source analytic bound, and the conditional infinite limit.
- [LDC-1](papers/native-layer-depth-compatibility/THEOREM.md): prior polynomial-carrier compression and the distinction between stable readings and identical cuts.
- [LAT-1](papers/native-layer-anchor-transport/THEOREM.md): probability transport within a fixed model.
- [POC-1](papers/native-prime-observation-closure/THEOREM.md): prime-bank sampling density and its unresolved universal source signs.
- [CCI-1](RH_CONJUGATE_CANCELLATION_INDEX.md): actual-source negative index.
- [BRR-1](RH_BOUNDARY_OBSERVER_RETURN.md): the remaining all-width return condition.

Finite spectral realization, Schur updates, Cauchy determinants,
and Schwarz division are retained or standard tools. The new
written estimates are the position-uniform complete-source
bounds, their mixed-sign range, and the resulting actual
finite-union and lower-envelope conclusions. No literature-wide
novelty, optimality, or machine-checked analytic proof is claimed.

Only this note and its README entry change. No scientific
computation, numerical source evaluation, independent agent,
source replay, or workflow dispatch is performed. Frozen data,
engines, certificates, and workflows remain unchanged. The
documentation commit uses [skip ci]. The PR is stacked on
draft #122 and remains an unmerged draft.
