# YM-59: scalar interaction, a local Poisson tower, and nonlinear vacuum energy

Monty Dabas. 4 October 2026. **Python 3.12 only.**

Add the original scalar overlap interaction to YM57's correlated kinetic
chain, with the same extensive-rate heat clock:
\[
 A_{N,k}(\theta)=H_{N,k}-\theta V_N,\qquad
 V_N=\sum_{i=1}^{N-1}v_i,\quad v_i=q_i\cdot q_{i+1},\quad k>0.
                                                               \tag{1}
\]
No equilibrium weight or ground transform is substituted for this
operator. Let E_N(theta)=inf spec A_{N,k}(theta). Its thermodynamic
energy density exists. At unit neighbour rate the new quantitative result is
\[
 \boxed{e(\theta):=\lim_{N\to\infty}{E_N(\theta)\over N}
       =-\frac{\theta^2}{4}-\frac{1123}{47376}\theta^4
                                      +O(\theta^6).}           \tag{2}
\]
The remainder follows from a local logarithmic Poisson recursion and
is uniform in volume before taking the limit. This is more than an
exchange of finite-volume Taylor coefficients and limits. At every
fixed order the same construction gives a uniform asymptotic expansion.
It does not prove convergence of the infinite series or a positive
interacting spectral gap at fixed nonzero theta.

## Carrier, lineage and scope

Use the existing native compact quaternion carrier Q, positive reference
Phi_Q and its coefficient/recognition completions, as in YM50--58.
Write A_i=D_{1,i}, Z_i=D_{2,i}, B_i=Z_i-Z_{i+1}; then
\[
 H_{N,k}=-\tfrac12\sum_i A_i^2-\tfrac{k}{2}\sum_i B_i^2,
 \qquad \gamma_k=\min\{1/630,k/420\}.
\]
All on-site YM56 lambda parameters may still be zero. Phi_N is the
product native reference, not the interacting vacuum distribution.
Norms and pairings below are in its completed representation. Work
with N>=2; the isolated one-site rank-one model has a larger kernel.

The derivative product rule, Rayleigh--Schrodinger perturbation and
logarithmic ground-source factorization are established methods. The
model-specific contribution is the local inverse supplied by native
support preservation, the volume-uniform fixed-order recursion, and
the exact connected scalar response calculated below. No new one-site
operator engine is introduced. Written proofs and exact finite controls
are distinct; this is not proof-assistant formalization or external review.

Known general gap-stability results require their own hypotheses. For
example, [Yarotsky's relatively bounded perturbation theorem](https://arxiv.org/abs/math-ph/0412040)
specifies a classical reference Hamiltonian and relative form bounds.
We do not identify this correlated, unbounded native kinetic operator
with that reference or import an interacting gap from that theorem.
YM53--55's site-product bridge hypotheses also do not directly cover it.

## YM59-T1: a cost for each excited site

Let P_i average over site i, Q_i=I-P_i, and let
\(\mathcal N_N=\sum_iQ_i\). This counts sites excited away from their
constant sector; it is not a physical particle-number postulate.
YM57's echo and compact-averaging argument actually proves
\[
 \sum_i\|Q_i\psi\|^2
 \le {63\pi^2\over2}\sum_i\|A_i\psi\|^2
                     +21\pi^2\sum_i\|B_i\psi\|^2.
\]
Its use of total variance was only a final weakening. Consequently
\[
                 \boxed{H_{N,k}\ge\gamma_k\mathcal N_N}         \tag{3}
\]
as closed forms. Each P_i commutes with H: invariance annihilates the
integral of every native derivative, derivatives annihilate constants,
and different sites commute. The same is true of sitewise antipodal
parity. Excitation-support and parity sectors are therefore reducing.

For a polynomial supported on S, H preserves that support. Its restriction
includes all bonds meeting S, exactly as in YM58, with only the actual
physical endpoints having degree one. On mean-zero functions of S its
inverse R_S exists, has L2 norm at most 1/gamma_k, and preserves every
finite invariant coefficient space. The inverse is genuinely local as
an operation on functions: it does not add coordinates outside S.
Its coefficients can depend on a physical boundary touching S.

## YM59-T2: a thermodynamic energy and uniform interaction bounds

The scalar multiplier V_N is real and bounded by M=N-1. Thus (1) is a
self-adjoint bounded perturbation, with the same form domain. Splitting
an open chain into two pieces leaves one positive kinetic square and
one scalar bond bounded below by -abs(theta). Hence
\[
 E_{N+L}\ge E_N+E_L-|\theta|,\qquad
            -|\theta|(N-1)\le E_N\le0.                        \tag{4}
\]
The upper bound tests the constant vector. For this energy tiling only,
E_1=0. The shifted sequence
E_N-abs(theta) is superadditive. Tiling into blocks of any fixed length,
and using (4) on the bounded remainder, proves that E_N/N converges to
the supremum of (E_N-abs(theta))/N. This proves existence of e(theta)
for every real theta. It is even, concave and 1-Lipschitz. Evenness
follows by changing q_i to -q_i on alternating sites: H is unchanged
and every scalar bond changes sign.

There is a useful stronger small-coupling bound. On a bond b={i,i+1}
put P_b=P_iP_{i+1}, R_b=Q_iQ_{i+1}, C_b=I-P_b. Native moments give
\[
 P_bM_{v_b}P_b=0,\qquad
 M_{v_b}P_b=R_bM_{v_b}P_b,\qquad \|M_{v_b}P_b\|=1/2.
\]
For a normalized vector and any delta>0, split its vacuum/active
blocks and apply Young's inequality to get
\[
 |\theta\langle v_b\rangle|
 \le\delta\langle R_b\rangle
       +{\theta^2\over4\delta}\langle P_b\rangle
       +|\theta|\langle C_b\rangle .                         \tag{5}
\]
The chain incidences imply \(\sum_b R_b\le\mathcal N_N\) and
\(\sum_b C_b\le2\mathcal N_N\). Set a=gamma_k-2 abs(theta)>0. From (3)--(5),
\[
 E_N(\theta)\ge-{(N-1)\theta^2\over4a},                     \tag{6}
\]
and the alternative choice delta=a/2 gives
\[
 A_{N,k}(\theta)\ge {a\over2}\mathcal N_N
                            -{(N-1)\theta^2\over2a}.         \tag{7}
\]
Every normalized state with energy <=0 therefore has
\(\langle\mathcal N_N\rangle\le(N-1)\theta^2/a^2\).
These bounds are uniform in width. They apply to minimizing sequences;
no assertion about attainment of the bottom at arbitrary coupling is needed.

For a strict upper bound take disjoint neighbouring dimers and at each
dimer use the real trial vector 1+c v. Its norm squared is 1+c^2/4.
Cross-dimer scalar means vanish by simultaneous antipodal parity, and
cross-dimer kinetic mixed expectations vanish by skewness. Each dimer's
energy is at most its full two-boundary marginal value. Write
\(\omega_k=(1+3k)/4\). Optimizing in the span of 1 and 2v gives
\[
 E_N(\theta)\le\left\lfloor N/2\right\rfloor
                  {\omega_k-\sqrt{\omega_k^2+\theta^2}\over2}.
                                                               \tag{8}
\]
This is strictly negative for theta!=0. At endpoints the omitted
boundary noise only lowers this trial energy. Consequently, for
abs(theta)<gamma_k/2,
\[
 -{\theta^2\over4a}\le e(\theta)
                \le {\omega_k-\sqrt{\omega_k^2+\theta^2}\over4}<0
\]
when theta!=0. These are energy and excitation-density bounds, not a
gap above the displaced interacting vacuum.

## YM59-T3: the local logarithmic Poisson tower

For real smooth native polynomials define
\[
 \Gamma(f,g)=\tfrac12\sum_i(A_if)(A_ig)
                 +\tfrac{k}{2}\sum_i(B_if)(B_ig).
\]
The exact product identity is
\(H(fg)=fHg+gHf-2\Gamma(f,g)\). If F is real and h=exp F, then
\[
 {A(\theta)h\over h}=HF-\Gamma(F,F)-\theta V,\qquad
 A(\theta)=\tfrac12\sum_X w_X(D_X-M_{D_XF})^*(D_X-M_{D_XF})
                         +M_{HF-\Gamma(F,F)-\theta V}.        \tag{9}
\]
The weighted labels are the existing site and bond directions. Read
(9) first on the polynomial form core and then by closure; smooth exp F
is admissible as a trial vector. The
square part annihilates h. Thus the spectral bottom lies between the
infimum and supremum of the displayed scalar residual: the lower bound
uses the squares and the upper bound tests h. This does not assume h
is the true ground vector.

Construct mean-zero f_n recursively:
\[
 \begin{split}
 &Hf_1=V,\qquad c_{1,N}=0,\\
 &g_n=\sum_{a+b=n}\Gamma(f_a,f_b),\quad
 c_{n,N}=-\Phi_N(g_n),\quad
 Hf_n=g_n-\Phi_N(g_n),\quad\Phi_N(f_n)=0\quad(n\ge2).
 \end{split}                                                \tag{10}
\]
All inversions are on centered finite coefficient spaces. The coupling
theta has not been frozen inside the recursion. Each f_n is a sum of
local polynomials, of total degree at most 2n and interval length at
most 2n, with only finitely many coefficient patterns at fixed n,k.

Here is the locality proof. It holds for f_1 because each v_i is a local
eigenfunction, computed below. In Gamma(f_a,f_b), two local terms can
contribute only if their intervals overlap or are adjacent across an
active kinetic bond. Their union then has length at most 2(a+b), and
their product degree is at most 2(a+b). Center each such term separately
and apply the local inverse in T1; neither operation enlarges support.
At a fixed order there are only finitely many relative offsets, degrees
and physical-endpoint patterns. The positive finite-dimensional inverse
has finite coefficients in each case. This proves both the induction
and a volume-independent bound on local coefficients. It also treats
short chains, where only finitely many additional boundary patterns occur.

For explicit remainder bookkeeping choose polynomial representatives
in a fixed basis modulo the native unit-radius relations, and let
L_n(k) bound
\(\sup_i\sum_{X\ni i}\|(f_n)_X\|_{\mathrm{coeff},1}\)
uniformly in N. Such a finite bound follows from the preceding finite
list of inverses. Coefficient norm bounds pointwise norm on Q. An axis
derivative on degree d has coefficient norm at most d/2 times the
input norm; a bond uses the two endpoint sums. It follows that
\[
 \|\Gamma(f_a,f_b)\|_\infty
       \le N\,{1+4k\over8}\,(2a)(2b)L_a(k)L_b(k).          \tag{11}
\]
The derivative estimate may be applied to the chosen representatives
before reducing their output modulo the radius relations.

For F_p=sum_{n=1}^p theta^n f_n, equations (9)--(10) cancel every
nonconstant coefficient through degree p. For abs(theta)<=1,
\[
 \boxed{\left|E_N(\theta)-\sum_{n=1}^p c_{n,N}\theta^n\right|
        \le NC_p(k)|\theta|^{p+1},}                        \tag{12}
\]
where one valid finite constant is
\[
 C_p(k)={1+4k\over8}
       \sum_{\substack{1\le a,b\le p\\a+b>p}}(2a)(2b)L_a(k)L_b(k).
\]
This is constructive finite-order control; this chapter does not list
numerical L_n or C_p for arbitrarily high n. No bound on their growth
with p, and hence no nonzero convergence radius, is inferred.
In particular C_5 is not numerically evaluated here: (2) is an asymptotic
statement, not a numerical error bar at a selected finite coupling.

At each fixed n, only a bounded number of terms in c_{n,N} touch a
physical endpoint. The remaining local coefficient patterns are
translates, so c_{n,N}/N converges to a bulk c_n. Taking the already
proved energy limit in (12) justifies the asymptotic expansion of e,
with the same C_p. Alternating antipodal parity commutes with H and
Gamma and reverses V. Uniqueness of each centered inverse gives
f_n parity (-1)^n; hence c_{n,N}=0 at odd n. Taking p=5 therefore
justifies the O(theta^6) remainder in (2).

At any fixed finite N, bounded perturbation of the isolated constant
eigenvalue also supplies the usual analytic eigenbranch near zero.
On the circle abs(z)=gamma_k/2 the free resolvent has norm at most
2/gamma_k. Its Neumann series applies when 2 abs(theta)(N-1)<gamma_k;
the enclosed projection remains rank one.
Matching coefficients in (9) shows that (10) gives those Taylor
coefficients. Its volume-uniform conclusion (12) comes from local
residual bounds, not from that finite-volume analytic radius.

## YM59-T4: exact linear vacuum response and quadratic susceptibility

Let d_i be the degree of site i in the finite open chain, equal to one
at endpoints and two in its interior. The native dot product obeys
\[
 H v_i=\omega_i v_i,\qquad
 \omega_i={2+k(d_i+d_{i+1}+2)\over8}.                       \tag{13}
\]
An individual endpoint derivative squared gives -v_i/4. The mixed
Z_i Z_{i+1} derivative gives +v_i/4: simultaneous equal-axis rotations
preserve the dot product. This fixes the cross-term sign in (13).
Thus bulk omega_i=omega_k=(1+3k)/4; an endpoint bond in N>=3 has
rate (2+5k)/8, and the N=2 rate is (1+2k)/4.

Native moments give Phi(v_i)=0 and Phi(v_i v_j)=delta_ij/4, even when
the bonds are adjacent. Therefore
\[
 f_1=\sum_i{v_i\over\omega_i},\qquad
 c_{2,N}=-\sum_i{1\over4\omega_i},\qquad
 \boxed{c_2=-{1\over1+3k}.}                                \tag{14}
\]
The unit-normalized finite-volume ground branch has first derivative
f_1, since its derivative is orthogonal to constants. Its bond response
at theta=0 is
\[
 {d\over d\theta}\langle v_i\rangle_{h_\theta^2}\big|_0
                  ={1\over2\omega_i}.
\]
The bulk quadratic energy susceptibility is 2/(1+3k), equal to 1/2
at k=1. This is the centered second-difference response of e at zero,
with the uniform expansion supplied by T3. It is a model response in
the declared heat clock, not a newly calibrated physical constant.

## YM59-T5: the connected fourth-order coefficient

For an independent check and efficient exact calculation, let
R=H^{-1} on the centered subspace and Q=I-|1><1|. Fixed-volume
perturbation with intermediate normalization gives
\[
 c_{4,N}=-\langle QV R V, R QV R V\rangle
                   +\langle V,RV\rangle\langle V,R^2V\rangle.   \tag{15}
\]
This follows by solving the order-one, order-two and order-three
eigenvector equations; the cubic energy vanishes by parity. The second
term is essential: it cancels contributions of spatially disconnected
bond pairs. It is not optional vacuum normalization.

First compute the bulk local inverses using YM58's marginal generators,
so every involved site's total lattice degree is two. Put
\[
 u_i=v_i^2-1/4,\quad
 a_0=\langle u_0,Ru_0\rangle,\quad
 a_1=\langle v_0v_1,R(v_0v_1)\rangle,\quad
 a_2=\langle v_0v_2,R(v_0v_2)\rangle.
\]
For k=1, exact native coefficient inversion gives
\[
 \boxed{\omega_k=1,\qquad a_0=1/30,\quad a_1=31/705,\quad a_2=2/63.} \tag{16}
\]
These are complete inversions of the invariant spaces generated by the
specified polynomial sources, of dimensions 3, 4 and 2. They are not
cutoff estimates for these sources. The certificate constructs their
rational polynomial inverse witnesses and verifies H R f=f in the
native pairing, using the unchanged polynomial engine and exact solver.

Different u_i are orthogonal even after R: u_i is centered separately
at both sites, and R commutes with every P_i. Bond-pair products have
distinct sitewise antipodal parity patterns unless the unordered edge
pair is identical. They are also orthogonal to every u_i. Thus (15)
splits into individual bonds and unordered bond pairs.

If two bonds are separated by at least one unused lattice site, their
supports have no connecting kinetic edge. The marginal generator
splits, so their product has rate omega_i+omega_j and squared norm 1/16.
Its contribution to the first term of (15) cancels its contribution
to the normalization term exactly. Pairs two edge-indices apart are
different: a kinetic edge joins their supports. This is the a_2 term.

At k=1 the surviving contributions per translated configuration are

| Scalar-bond configuration | Fourth-order connected contribution |
|---|---:|
| One repeated bond | 1/16-a_0 = 7/240 |
| Two adjacent bonds | 1/8-4a_1 = -287/5640 |
| Two bonds joined by one kinetic edge | 1/8-4a_2 = -1/504 |
| More distant bonds | 0 |

There is one of each of the first three configurations per bulk site.
Their sum is -1123/47376, proving the coefficient in (2). For a marginal
interval with M>=3 scalar bonds, the exact coefficient is
\[
 c_{4}^{m}(M)=-{1123\over47376}M+{6497\over118440}.           \tag{17}
\]
This auxiliary marginal boundary convention makes the calculation
uniform; it has not replaced the physical open-chain operator (1).
For the latter, use its actual omega_i and local boundary inverses.
The general self contribution is
\(1/(16\omega_i^3)-a_{0,i}/\omega_i^2\); a close pair contributes
\[
 {\omega_i+\omega_j\over16\omega_i^2\omega_j^2}
       -(1/\omega_i+1/\omega_j)^2a_{ij}.                    \tag{18}
\]
All farther pairs still cancel, with their actual possibly unequal
rates. Only finitely many configurations touch an endpoint, so they
change the boundary intercept and not the bulk coefficient. The
certificate independently verifies this stabilization for open chains.

## YM59-T6: what the new energy result does and does not settle

The scalar interaction now has a thermodynamic energy, uniform
quadratic bounds, a controlled local response tower, and an exact
nonlinear bulk coefficient. Equations (3), (9) and (10) are the mechanism:
local native transport makes the centered Poisson inverses local, while
the logarithmic square identity converts finite-order cancellation into
uniform energy bounds.

Energy control alone does not give a gap above the interacting vacuum.
For fixed N, the elementary bounded-perturbation estimate gives at least
gamma_k-2 abs(theta)(N-1) as long as this is positive. Its coupling
window shrinks with N. Neither that estimate nor (6)--(12) proves a
fixed nonzero theta window for the interacting gap or its history limit.
The residual bound in (12) is extensive; it must not be relabeled as
a small operator-norm error uniform in N after ground-energy subtraction.

Required negative controls include deleting the normalization term in
(15), treating a_2 as the disconnected value 1/32, omitting boundary
noise, and using the shrinking finite-volume gap window as a uniform
one. The old scalar interaction has been retained throughout this chapter.
Physical gauge action and clock selection, gauge constraints, spatial
continuum and identification with four-dimensional Yang--Mills remain open.

## Reproduction

The implementation reuses YM58/YM57's native derivatives, product
moments, coefficient pairing and inherited exact linear solver.
Its source manifest preserves 146 preceding upstream files. The written
locality proof establishes all fixed orders; numerical certificates
check specified finite inverses, coefficients and bounds rather than
pretending to enumerate an infinite tower.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym59_scalar_vacuum_response.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym59_scalar_vacuum_response.py -v
```

Default and `--check` are read-only. `--write` changes only YM59's own
result and hash. Python 3.12 with assertions enabled is required.
