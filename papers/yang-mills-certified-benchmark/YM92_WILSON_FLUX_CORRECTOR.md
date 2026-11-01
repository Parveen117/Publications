# YM92 — A native Wilson flux with a controlled nonlinear residual

Monty Dabas. 6 October 2026. Python 3.12 evidence packet.

YM91 supplied a conditional-flux criterion but its explicit native family
was an engineered comparison Hamiltonian. This continuation constructs a
first-order plaquette flux for the **actual Wilson model**, bounds its
complete nonlinear residual, and corrects it to an exact feasible flux.
The remainder is controlled on the full compact configuration space;
it is not discarded as a formal higher-order term.

For 0 <= theta <= 1/400 the actual, volume-uniform gap floor improves from
1/2-80 theta to

\[
\boxed{d_{92}(\theta)=\tfrac12-72\theta-1000\theta^2.}
\tag{1}
\]

For the retained straight b-link products of YM90, b>=2, the relative
zero-frequency memory satisfies

\[
\boxed{\Sigma_b(0)\le\eta_{92}(\theta)C_b,\qquad
\eta_{92}(\theta)=
\left(\frac{4\theta/9+1112\theta^2}{d_{92}(\theta)}\right)^2.}
\tag{2}
\]

At the endpoint,

\[
d_{92}=\frac{251}{800},\qquad
\eta_{92}=\left(\frac{1451}{56475}\right)^2
          =\frac{2105401}{3189425625}<\frac1{1500}.
\tag{3}
\]

The memory upper coefficient is more than 10.5 times smaller than YM90's
1/144. This compares proved upper bounds, not measured values of the
actual memory. The bare gap floor improves from 3/10 to 251/800.
The coupling window is unchanged. It is still a strong-coupling lattice
window, not the weak-coupling continuum trajectory.

## 1. Actual model and the inherited estimates

On the periodic spatial cubic lattice of side Lb, L>=4, use

\[
A_\theta=-\Delta-\theta\sum_pW_p,\quad
A_\theta h=E_0h,\quad h>0,\quad\chi=\log h,\quad d\mu=h^2d\Phi.
\tag{4}
\]

The native SU(2) factor is the radius-two round three-sphere, with
fundamental Casimir 3/4 and Ricci tensor (1/2)g. Every link meets four
plaquettes. Write S_ef=Hess_ef chi, using covariant Hessians throughout.
YM88 proves, on the entire configuration space and uniformly in volume,

\[
G:=\max_e\sup|\nabla_e\chi|\le42\theta,\quad
\delta:=\max_e\sup\|S_{ee}\|\le35\theta,\quad
\beta:=\max_e\sum_{f\ne e}\sup\|S_{ef}\|\le5\theta.
\tag{5}
\]

It also proves for e!=f the tensor equation

\[
\mathcal T_\chi S_{ef}
 =\operatorname{Hess}_{ef}V+2\sum_g S_{eg}S_{gf},\quad
\mathcal T_\chi=-\Delta^{\rm rough}-2\nabla_{\nabla\chi}+1,
\quad V=\theta\sum_pW_p.
\tag{6}
\]

The connection diffusion is contractive in the supremum fiber operator
norm, so the inverse in (6) has that norm at most one. The factor 1 is
the sum of the two factor Ricci contributions. Replacing the rough
Laplacian by a scalar component Laplacian would change (6).

## 2. T1 — Resolve the leading native response before bounding the rest

Define the explicit first-order field

\[
\chi_1=\frac\theta3\sum_pW_p,\qquad
S^1_{ef}=\operatorname{Hess}_{ef}\chi_1,\qquad E_{ef}=S_{ef}-S^1_{ef}.
\tag{7}
\]

This definition does not identify chi with chi_1. Each Wilson plaquette
has electric eigenvalue 3, whence

\[
(-\Delta^{\rm rough}+1)S^1_{ef}=\operatorname{Hess}_{ef}V,
\quad
\mathcal T_\chi E_{ef}
 =2\sum_gS_{eg}S_{gf}+2\nabla_{\nabla\chi}S^1_{ef}.
\tag{8}
\]

**Native derivative lemma.** For distinct links e,f of a plaquette,
the operator norm of Hess_ef W_p is at most 1/4. For every link g of
that plaquette,

\[
\|\nabla_v\operatorname{Hess}_{ef}W_p\|_{\rm op}
 \le |v|/8,\qquad v\in T_gSU(2).
\tag{9}
\]

For three distinct factors this follows from three unit quaternion
insertions, each with norm 1/2. If g=e, multilinearity of W_p in the
e quaternion gives Hess_ee W_p=-(W_p/4)g_e. Differentiation in the
different f factor gives
nabla_e Hess_ef W_p=-(d_fW_p/4) tensor g_e in its appropriate slots.
Its norm is at most 1/8. The case g=f is the same. Covariant connection
terms are essential for these repeated-factor identities. A link
outside the plaquette gives zero. This proves (9) for arbitrary unit
directions, not merely the native coordinate axes.

Let N_ef be the number of plaquettes containing e and f. Its off-diagonal
row sum is twelve. Therefore

\[
\max_e\sum_{f\ne e}\sup\|S^1_{ef}\|\le\theta,
\qquad
\max_e\sum_{f\ne e}\sup\|2\nabla_{\nabla\chi}S^1_{ef}\|
 \le4\theta G.
\tag{10}
\]

For the second inequality, one incidence (e,f,p) contributes at most
2(theta/3) times 4G/8 = theta G/3. There are twelve such incidences
in a row. No factor for the total number of links appears.

Set eps=max_e sum_(f!=e) sup||E_ef||. Taking the tensor inverse bound
in (8), separating the two diagonal factors in S squared, gives

\[
\varepsilon\le4\delta\beta+2\beta^2+4\theta G
                \le918\theta^2.
\tag{11}
\]

Also beta<=theta+eps by (10). Substituting this and (5) back into
the first inequality yields

\[
\varepsilon\le310\theta^2+144\theta\varepsilon+2\varepsilon^2.
\tag{12}
\]

The inherited bound (11), rather than an unjustified choice of a
quadratic root, now closes the estimate:

\[
\varepsilon\le
 \frac{310\theta^2}{1-144\theta-1836\theta^2}
 \le\frac{12400000}{25141}\theta^2<500\theta^2.
\tag{13}
\]

The denominator is positive on the full declared window. Thus

\[
\boxed{\varepsilon\le500\theta^2,\qquad
\beta\le q:=\theta+500\theta^2,\qquad
\|\operatorname{Hess}\chi\|_{\rm op}\le
 h_*:=36\theta+500\theta^2.}
\tag{14}
\]

These are all-orders estimates for the actual logarithmic vacuum. The
block Schur estimate and weighted Bochner formula give
Ric-2 Hess chi >= (1/2-2h_*)I=d_92 I. This proves (1), both for the
full actual ground-state generator and for every coordinate conditional
with arbitrary fixed exterior. Multiplication by h transfers the full
bound to A_theta-E_0. The Gauss-invariant subspace inherits it.

## 3. T2 — A local plaquette flux and its exact nonlinear correction

Use the retained links R and unretained links Z of YM90. Gauge fixing
the internal path vertices gives the actual conditional measure

\[
d\nu_Q(z)\propto h(\operatorname{section}(Q),z)^2d\Phi_Z(z),\qquad
A_Q=-\nu_Q^{-1}\operatorname{div}_Z(\nu_Q\nabla_Z).
\tag{15}
\]

The preceding theorem gives A_Q>=d=d_92 on its centered space,
uniformly in Q. Let u_E be arbitrary coarse tangent covectors, fixed
at Q, and copy them to the fine retained links using YM90's prefix
adjoint frames. Denote the copied vector by u_R; then
|u_R|^2=b sum_E|u_E|^2. Prefix frames do not depend on z.

Write S_E for the actual horizontal logarithmic gradient in YM90.
The scalar conditional source is

\[
r_u=2\sum_Eu_E\cdot(S_E-\nu_Q S_E).
\tag{16}
\]

For each plaquette meeting R, put k_p=|p intersection Z| and

\[
r^1_{p,u}=\frac{2\theta}{3}
       \sum_{e\in p\cap R}\nabla_eW_p\cdot(u_R)_e,
\quad \phi^0_u=\sum_{p:p\cap R\ne\varnothing}
       \frac{r^1_{p,u}}{3k_p/4},\qquad J^0_u=\nabla_Z\phi^0_u.
\tag{17}
\]

An elementary fine plaquette meets at most two retained links. To see
this, its opposite links have transverse coordinates differing by one;
both cannot belong to the retained spacing-b grid when b>=2. There is
at most one retained edge of either orientation. Hence k_p is 2 or 3.
Each unretained link in r^1_p occurs in the fundamental representation,
so

\[
-\Delta_Z r^1_{p,u}=(3k_p/4)r^1_{p,u},\qquad
-\operatorname{div}_{\Phi_Z}J^0_u=r^1_u:=\sum_p r^1_{p,u}.
\tag{18}
\]

This is an exact Haar flux for the leading response, not yet a flux
for the actual conditional source. All Wilson words, adjoint frames,
source components and their cross terms are retained.

### Pointwise flux and derivative bounds

The block mapping a retained e vector to a flux on z is bounded by
2theta/(9k_p) per incident plaquette. Its column sum is at most 8theta/9:
the k_p unretained outputs cancel its denominator and e meets four faces.
Its row sum is at most 8theta/9 as well, since |p intersection R|/k_p<=1.
The rectangular Schur estimate gives

\[
|J^0_u|\le\sqrt b\,j\,|u|,\qquad j=8\theta/9.
\tag{19}
\]

The covariant derivative term (nabla_Z J^0_u)^* grad_Z chi has,
for each (w,e,p) with w in Z and e in R, norm at most theta G/9.
Indeed the sum over its k_p flux links cancels k_p in (17), and (9)
supplies 1/8. Its row sum is at most 8theta G/9 and column sum at
most 4theta G/3. In particular,

\[
|(\nabla_ZJ^0_u)^*\nabla_Z\chi|
 \le\sqrt b\,(4\theta G/3)|u|.
\tag{20}
\]

The bound uses the largest of those two sums; no unproved spectral
norm is inserted. It holds pointwise and hence under the actual nu_Q.

### Full residual, including the actual density drift

Define

\[
e_u=r_u+\operatorname{div}_{\nu_Q}J^0_u
    =r_u-r^1_u+2\nabla_Z\chi\cdot J^0_u.
\tag{21}
\]

Its nu_Q mean is exactly zero, by weighted integration by parts and
centering of r. Haar centering of r^1 cannot be substituted for this
actual centering. Differentiating (21) in Z and using (14), (19), (20),

\[
|\nabla_Ze_u|\le\sqrt b\,m\,|u|,\qquad
m=2\varepsilon_*+2h_*j+\tfrac83\theta G_*,
\quad\varepsilon_*=500\theta^2,\quad G_*=42\theta.
\tag{22}
\]

The three terms respectively bound the true mixed-Hessian remainder,
Hess_ZZ chi acting on the trial flux, and the covariant derivative of
that flux contracted with the full logarithmic gradient. The subtracted
conditional mean in r is constant in Z, so introduces no Z derivative.

Conditional Poincare and the positive inverse now give

\[
\nu_Q e_u^2\le b m^2|u|^2/d,
\quad K_u=\nabla_Z A_Q^{-1}e_u,\quad
\nu_Q|K_u|^2=\langle e_u,A_Q^{-1}e_u\rangle\le b m^2|u|^2/d^2.
\tag{23}
\]

Consequently J_u=J^0_u+K_u is an **exact actual feasible flux**:
-div_nu J_u=r_u. Its L2 norm is at most sqrt(b)(j+m/d)|u|.
The inverse in K is the full conditional inverse; it is not numerically
truncated or replaced by the free inverse used only in (17).
Equations (22)–(23) bound every mode of this correction.

YM91's minimum-flux principle therefore proves

\[
T_0(Q)\le b(j+m/d)^2 I.
\tag{24}
\]

Since d+2h_*=1/2, the bound simplifies exactly:

\[
j+m/d=\frac{j/2+2\varepsilon_*+(8/3)\theta G_*}{d}
       =\frac{4\theta/9+1112\theta^2}{d}.
\tag{25}
\]

Insertion in the actual conditional-to-dynamic comparison of YM91
proves (2). This is the desired actual Wilson flux estimate. The
factor b cancels against C_b[f]=b integral |grad f|^2 w_b dPhi.

## 4. Frequencies, projected resolvents and physical interpretation

The improved mixed-Hessian bound also gives the actual source estimate

\[
\|Bf\|^2\le\sigma^2 C_b[f],\qquad\sigma^2=4q^2/d.
\tag{26}
\]

Retain both valid bounds, for s>=0:

\[
0\le\Sigma_b(s)\le\alpha(s)C_b,\qquad
\alpha(s)=\min\{\eta_{92},\sigma^2/(d+s)\}.
\tag{27}
\]

For s>0 the projected resolvent R_b(s)=I_b^*(L+s)^-1 I_b obeys

\[
0\le R_b(s)-(C_b+s)^{-1}
 \le\frac{\alpha(s)}{1-\alpha(s)}(C_b+s)^{-1}.
\tag{28}
\]

At s=0 this holds on centered coarse functions, using the inherited
positive coarse floor. At the endpoint the relative coefficient is
2105401/3187320224 < 1/1499, compared with the old 1/143.
It is a Loewner form bound, not an entrywise or semigroup error.
One cannot replace (27) by eta_92 d/(d+s) without an additional argument;
the corresponding inherited counterexample is replayed.

The gap bound (1) holds on the full finite configuration space and hence
on the physical Gauss-invariant sector. A physical clock factor kappa
multiplies it. The estimate is uniform in finite volume; an infinite-volume
or continuum representation still needs its own construction.

The connection with NSB2 is methodological and concrete: the correction
requires marked mixed response and covariant third derivatives, together
with a norm for the entire omitted response. A same-matrix power or a
list of curvature eigenvalues supplies none of these estimates. NSB2's
Euclidean Hessian rigidity theorem is not applied as a theorem about the
curved SU(2) configuration manifold or the quantum Hamiltonian.

## 5. Scope, evidence and the next remaining gate

The following are proved here for the actual Wilson vacuum, conditional
on the unchanged YM88 estimates that this packet explicitly inherits:

- the covariant first-response remainder in (13)–(14);
- the improved full and conditional bare gap floor in (1);
- the actual, fully corrected conditional flux in (17)–(25);
- the volume- and block-factor-uniform memory/resolvent bounds (2), (27), (28).

The proof has no interacting vacuum ansatz and no discarded perturbation
tail. Exact finite checks cover the native third derivatives and rough
resolvent normalization, Haar plaquette flux identities, actual-density
drift signs on polynomial comparison controls, geometric incidence,
residual Gram enclosures, and rational endpoint arithmetic. The polynomial
comparison controls and finite matrices are not numerical Wilson vacua.
Finite checks support the analytic proof; they are not formal verification
or independent peer review.

The interval remains 0<=theta<=1/400. Under the existing physical scaling,
theta=c_B/(c_E g^4) grows when g tends to zero. No bound in this packet is
claimed on that trajectory. The remaining target is a conditional flux
and a physical hidden-sector floor controlled in that regime, together
with nonzero continuum states and the continuum axioms. The 4D continuum
mass gap remains open in this programme.

The general conditional-expectation and Poisson methods are established;
see F. Legoll and T. Lelievre, *Effective dynamics using conditional
expectations*, Nonlinearity 23 (2010), 2131–2163,
<https://arxiv.org/abs/0906.4865> (read 6 October 2026).
The numerical constants and actual Wilson combination above are proved
in the stated native normalization; no general priority claim is made.

Parent commit: bec73061e92b56f4772d5a31063dde4da0bcaa65 (NSB2).
The source manifest preserves that stage and the YM91/YM88 chain.
Default and --check are read-only. Only --write regenerates this stage's
result and pin.

```sh
python -B papers/yang-mills-certified-benchmark/certificates/ym92_wilson_flux_corrector.py --check
python -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym92_wilson_flux_corrector.py -v
python -B papers/yang-mills-certified-benchmark/certificates/ym91_conditional_memory_flux.py --check
python -B papers/yang-mills-certified-benchmark/certificates/ym88_spatial_score_contraction.py --check
```
