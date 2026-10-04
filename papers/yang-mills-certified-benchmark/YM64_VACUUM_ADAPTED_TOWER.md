# YM-64: vacuum-adapted blocks and the positive inverse-response tower

Monty Dabas. 4 October 2026. **Python 3.12 only.**

YM63's free-reference cut left a hidden-sector estimate containing
the overlap with the interacting vacuum. That overlap need not stay
positive uniformly in volume. Here we use the actual vacuum measure
to define the cut. The retained space then contains the whole vacuum,
and the hidden-sector lower bound becomes

\[
 \boxed{D_r\ge\gamma(\theta)I\quad\text{at every finite block depth}.}
                                                               \tag{1}
\]

This gives an exact composition rule for repeated spatial blocks and
uniform **energy-relative** bounds on every inverse-response coefficient.
It does not give an absolute operator-norm bound or all-order spatial
locality. A native gauge-invariant loop has nonzero first-order leakage;
the adapted cut therefore does not erase the interaction by definition.

## Model, dependencies and the change of cut

Keep YM62's Hamiltonian A=L_E-theta V on a finite periodic cubic lattice,
with native SU(2)/quaternion links, fundamental Casimir rate 3/4 and
V=sum_p W_p. Let E be its lowest eigenvalue and h>0 its normalized
vacuum. Put H=A-E. For |theta|<1/8000, YM62 supplies
H>=gamma(theta)(I-|h><h|), with gamma independent of volume.
All results below concern finite volumes with coarse side at least four,
and any finite sequence of integer blocks that fits inside that volume.
Uniform constants alone are not a construction of their joint limit.

Define the probability measure dmu=h^2 dPhi and the unitary map
U:L2(mu)->L2(Phi), Uf=hf. In L2(mu), write

\[
 \mathcal L=U^*HU,\qquad \mathcal L1=0,\qquad
 \mathcal E(f,g)=\sum_{e,a}\int
       \overline{D_{a,e}f}\,D_{a,e}g\,d\mu.                 \tag{2}
\]

The identity follows by integration by parts and Ah=Eh. In particular
the gap is equivalent to E(f,f)>=gamma Var_mu(f). The same positive
vacuum is used at all block levels; it is not solved anew in a truncated
coarse model.

Let pi_r multiply r successive fine links along each retained coarse
link, as in YM63. Let mu_r=(pi_r)_*mu, I_r f=f o pi_r, and
P_r=I_r I_r*, Q_r=I-P_r. Now I_r* is conditional expectation for **mu**,
not for the free product reference Phi. In the original Hilbert space,
the retained vectors are h f(pi_r), with f in L2(mu_r). Thus their
projection retains h exactly. YM63's free projection is not silently
redefined or its certificate changed.

Schur elimination and conditional-expectation reduction are established
mathematical tools. See [Griesemer and Hasler, arXiv:0704.3244,
Introduction](https://arxiv.org/abs/0704.3244) for the Feshbach--Schur
construction and [Nueske, Koltai, Boninsegna and Clementi,
Entropy 23 (2021), 134](https://doi.org/10.3390/e23020134) for reversible
conditional-expectation forms. The proofs here specify the native
product-link cut, its domains and the bounds inherited from YM62.
[GE4](../generalized-euler-evolution/GE4_MINIMAL_OBSERVER_MEMORY.md)
already uses positive moment matrices on bounded native sectors. The
tower below instead uses inverse powers of the gapped, generally
unbounded hidden operator. No priority claim is made for Schur
associativity or positive moment algebra.

## YM64-T1: the adapted projection and its exact coarse energy

At each fixed finite volume h is smooth and strictly positive. This
follows from the elliptic equation Ah=Eh on the compact product of
native link manifolds and positivity of the vacuum. Thus h and h^(-1)
are bounded there, though no volume-independent bounds are assumed.
The form domain in (2) is the ordinary first Sobolev space with its
equivalent weighted norm.

There is a global smooth product coordinate change for pi_r: replace
the last link on each disjoint retained path by the ordered product of
the path, and keep its other r-1 links and all unretained links. The
inverse solves for that last link by an ordered inverse product. This
change preserves product reference integration. In these coordinates
the conditional density is smooth and positive on a compact product.
Integration over the hidden coordinates, followed by division by the
positive marginal density, maps H1 continuously to coarse H1. Pullback
is also continuous in H1. Consequently P_r and Q_r preserve the form
domain. The constants in this domain argument may depend on volume
and r; they are not the uniform bounds asserted below.

The native adjoint rotation at each factor of a product link preserves
the sum over the three derivative directions. Disjoint paths give no
cross-link term. The pointwise polarized identity is therefore

\[
 \Gamma_{\rm fine}(I_r f,I_r g)
       =r I_r\Gamma_{\rm coarse}(f,g).                    \tag{3}
\]

Integrating against the actual mu gives the closed diagonal coarse form

\[
 \boxed{\mathcal C_r[f,g]=\mathcal E(I_r f,I_r g)
       =r\int\Gamma_{\rm coarse}(f,g)\,d\mu_r.}           \tag{4}
\]

If w_r=dmu_r/dPhi_coarse, its generator on smooth functions is
C_r f=-r w_r^(-1) sum_(e,a) D_(a,e)(w_r D_(a,e)f).
The principal diffusion coefficient is exactly r. The drift depends
on the actual coarse vacuum density, which can be spatially nonlocal.
Equation (4) is the diagonal form, not the exact interacting reduced
time evolution.

Gauge invariance of h and equivariance of pi_r make P_r commute with
the fine gauge action. Restricting to gauge-invariant readouts preserves
all subsequent forms and bounds. The argument is not restricted to
gauge-variant probes.

## YM64-T2: a uniform hidden gap and a Schur form without bounded B

Because P_r1=1, every hidden vector is orthogonal to the vacuum 1.
The restriction of E to Q_r H1 is a densely defined closed form on
Q_r L2(mu). Its self-adjoint operator D_r obeys (1) directly. There
is no factor ||P_free h||^2. This is a change of observation space,
not a stronger estimate for the old free-reference projection.

Write C=C_r, D=D_r and omit r temporarily. The off-diagonal coupling
need not be a bounded operator on L2; we define it through the form.
For each coarse f in the form domain, energy Cauchy--Schwarz gives

\[
 |\mathcal E(I f,g)|^2\le \mathcal C[f]\,\|D^{1/2}g\|^2,
                    \qquad g\in QH^1.                    \tag{5}
\]

Since D>=gamma, D^(1/2) maps its form domain onto the hidden Hilbert
space. Riesz representation defines a unique vector T f by
E(g,I f)=<D^(1/2)g,T f>, with ||T f||^2<=C[f]. On smooth coarse
readouts T f=D^(-1/2)Q L I f. This is a form-energy bounded map;
no claim that T or the formal coupling B is L2-bounded is needed.

For s>=0 define the positive correction form

\[
 \Sigma_r(s)[f,g]
   =\langle T_rf,D_r(D_r+s)^{-1}T_rg\rangle,
 \qquad \mathcal F_r(s)=\mathcal C_r-\Sigma_r(s).           \tag{6}
\]

Minimizing the hidden component by completing the square proves

\[
 (\mathcal F_r(s)+s)[f]
 =\inf_{g\in Q_rH^1}\{
     \mathcal E(I_rf+g)+s(\|f\|^2+\|g\|^2)\}.             \tag{7}
\]

In (7), '+s' means adding s||f||^2. The minimizer is
g_s=-D^(-1/2)(I+sD^(-1))^(-1)T f, which lies in the hidden
form domain. Form continuity of P_r and the gap on the centered
space make the minimized norm equivalent to the coarse form norm
at each finite volume. Hence the minimized forms are closed,
including at s=0 after separating the constant. Variation of (7)
also proves, for s>0,

\[
 \boxed{(F_r(s)+s)^{-1}=I_r^*(\mathcal L+s)^{-1}I_r.}       \tag{8}
\]

The inverse of the compressed resolvent is an unbounded positive
operator on its dense range; (7) specifies its form and domain. This
avoids asserting a block-operator domain decomposition that has not
been proved. For centered f, every I_r f+g in (7) is centered, so

\[
 \boxed{F_r(s)\ge\gamma(I-\Pi_r),\quad F_r(s)1=0,
 \quad 0\le\Sigma_r(0)\le C_r-\gamma(I-\Pi_r).}            \tag{9}
\]

Here Pi_r projects onto the constant. Inequalities are quadratic-form
inequalities. They hold with the same gamma at every admissible volume
and finite block depth. They depend on the already proved YM62 gap;
they are not an independent proof of that gap.

## YM64-T3: higher responses form a positive operator moment tower

For n>=0 set

\[
 M_{r,n}[f,g]=\langle T_rf,D_r^{-n}T_rg\rangle.
                                                               \tag{10}
\]

M_(r,0)=Sigma_r(0), and functional calculus gives

\[
 \boxed{0\le M_{r,n}\le\gamma^{-n}\Sigma_r(0)
          \le\gamma^{-n}C_r.}                             \tag{11}
\]

For |s|<gamma the expansion
Sigma_r(s)=sum_(n>=0)(-s)^n M_(r,n) converges in the energy-relative
form norm. For real s it is the Taylor expansion with
Sigma_r^(n)(0)=(-1)^n n! M_(r,n). For complex s the same power
series defines the analytic resolvent form in this disk.
For q=|s|/gamma<1 and a fixed coarse f,

\[
 \left|\Sigma_r(s)[f]-\sum_{n=0}^N(-s)^nM_{r,n}[f]\right|
 \le\frac{q^{N+1}}{1-q}\,\Sigma_r(0)[f].                  \tag{12}
\]

Polarization gives the corresponding bilinear bound with
||T_rf|| ||T_rg||. The disk and the relative constants are uniform;
C_r[f] itself need not be bounded for a sequence of growing readouts.
For real s>=0 a sharper signed finite-order bound holds, even when
s>=gamma and the infinite Taylor series is not certified to converge:

\[
 0\le(-1)^{N+1}\left(\Sigma_r(s)
             -\sum_{n=0}^N(-s)^nM_{r,n}\right)
 \le\frac{(s/\gamma)^{N+1}}{1+s/\gamma}\Sigma_r(0).        \tag{13}
\]

To prove it use the exact scalar geometric remainder on the spectrum
x of D^(-1), where 0<=x<=1/gamma. The function
(sx)^(N+1)/(1+sx) is increasing for x>=0. At N=0 this yields

\[
 0\le F_r(s)-F_r(0)
       \le\frac{s}{\gamma+s}\Sigma_r(0)
       \le\frac{s}{\gamma+s}C_r.                          \tag{14}
\]

More generally (-1)^n Sigma_r^(n)(s) is positive for s>=0.
This is a positive resolvent, or Stieltjes, response structure.
The operators M_n need not commute: the certificate includes exact
noncommuting M_0,M_1.

The constraints connect different orders, rather than merely assigning
each a sign. For any finite sequence of coarse form vectors f_i,

\[
 \sum_{i,j}M_{r,i+j}[f_i,f_j]
       =\left\|\sum_iD_r^{-i}T_rf_i\right\|^2\ge0.        \tag{15}
\]

The shifted matrix (M_(i+j+1)) is also positive. So is the matrix
(gamma^(-1)M_(i+j)-M_(i+j+1)), because
gamma^(-1)I-D^(-1)>=0. These are block Hankel moment constraints
with a spectral support bound, valid for noncommuting response
operators. Positivity and the hidden gap are the mathematical
requirements forcing this particular tower.

In time, the return kernel is the form
K_r(t)[f,g]=<T_rf,D_r exp(-tD_r)T_rg> for t>0. Spectral calculus gives

\[
 \int_0^\infty\frac{t^n}{n!}K_r(t)\,dt=M_{r,n},\qquad
 0\le\int_T^\infty K_r(t)\,dt
       \le e^{-\gamma T}\Sigma_r(0).                     \tag{16}
\]

For t>0, K_r(t)<=m_gamma(t)Sigma_r(0), where m_gamma(t)=1/(e t)
for t<=1/gamma and gamma exp(-gamma t) otherwise. The bound may
diverge as t decreases to zero; bounded coupling at t=0 is not
assumed. Its integrated form (16) remains finite. Zero M_0 is
equivalent to zero cross form T, hence to a reducing cut and closed
retained heat dynamics. It does not characterize symmetry in general.

## YM64-T4: repeated cuts compose at the same spectral parameter

Let r=b_1...b_k. The ordered product maps obey pi_r=pi_(b_k) o ... o
pi_(b_1) with consistent origins. Their marginal measures are actual
pushforwards at each level. Therefore their isometries compose,
I_r=I_(b_1)...I_(b_k), and their conditional expectations obey the
usual tower identity. This uses no independence of the interacting
coarse links.

For two steps b,c, (8) immediately gives, at the same s>0,

\[
 \boxed{F_{bc}(s)+s
    =\left[I_c^*(F_b(s)+s)^{-1}I_c\right]^{-1}.}           \tag{17}
\]

The identity iterates any finite number of times. A second derivation
comes from (7): minimizing over the first hidden space and then the
second minimizes over their orthogonal sum. Coarse form domains are
compatible by T1. Direct total blocking obeys the same hidden-gap,
moment and integrated-memory bounds (1), (11)--(16); their constants
do not acquire a factor per block. Its diagonal energy is exactly
r times the weighted coarse gradient energy in (4).

Replacing F_b(s) by F_b(0), C_b, or a fitted single plaquette coupling
changes (17). In particular a static zero-energy Schur matrix is not
the whole spectral family. The certificate gives a rational positive
operator for which this replacement disagrees with direct elimination
at each of three positive spectral parameters. Nested cuts with the
full parameter dependence are exact; scalar coupling closure has not
been established. Equation (17) also applies to weighted finite
conditional expectations, checked separately from the orthogonal
matrix fixture.

## YM64-T5: the actual gauge model has a local first return coefficient

At fixed finite volume, bounded analytic perturbation and simplicity
of the free vacuum give a vacuum expansion near theta=0. Normalize
it by Phi(h)=1 temporarily; this does not change its logarithmic
derivatives. Every elementary W_p has free electric rate 3 and
Phi(V)=0. Hence

\[
 h=1+\tfrac\theta3 V+O_\Lambda(\theta^2),\qquad
 \mathcal L f=L_E f-\tfrac{2\theta}{3}\Gamma(V,f)
                         +O_\Lambda(\theta^2).            \tag{18}
\]

The expansion is on a fixed smooth readout. Its radius and remainder
are not asserted uniform in volume or readout complexity. Let P_0
be the free-reference product-link projection and r>=2. YM63 gives
P_0 V=0 and reduction of L_E. For a retained f, the product rule
L_E(VI_r f)=V L_E I_r f+(L_E V)I_r f-2Gamma(V,I_r f)
then implies P_0 Gamma(V,I_r f)=0. Also L_E I_r f is retained for
every theta. Thus the actual adapted hidden leakage has expansion

\[
 \boxed{Q_{r,\theta}\mathcal L_\theta I_r f
     =-\tfrac{2\theta}{3}\Gamma(V,I_r f)
                         +O_\Lambda(\theta^2).}           \tag{19}
\]

This response is spatially local at its first order: only plaquettes
sharing a link with I_r f contribute. Write R_f=Gamma(V,I_r f) and
let L_(0,Q) be the free Casimir on the hidden space. At fixed volume,
after identifying the smooth positive weighted spaces, (10) has
leading coefficient

\[
 M_{r,n}[f,g]
   =\frac{4\theta^2}{9}
       \langle R_f,L_{(0,Q)}^{-(n+1)}R_g\rangle_\Phi
                         +O_\Lambda(\theta^3).            \tag{20}
\]

The leading term is independent of first-order changes of inner
product, since each coupling already contributes one theta. The
free inverse preserves the set of link variables on which a centered
readout depends. Consequently the leading bilinear coefficient in
(20) is zero whenever the link supports of R_f and R_g are disjoint:
product-reference independence and their zero means prove this.
This is a local leading coefficient, not an all-order locality bound.

There is an explicit physical witness that (19) is nonzero. Use r=3
and let f be a coarse elementary Wilson plaquette. Its fine boundary
has twelve links. Put every fine link at identity except the middle
link of one of its sides, which is q=(3/5,4/5,0,0). Only the four
fine plaquettes incident on that link have nonzero gradients of W_p
at this configuration. Each shares just that link with the coarse
boundary. Its contribution is (4/5)^2/4=4/25. Hence

\[
 \Gamma(V,I_3 f)=16/25,\qquad
 [\theta]\,Q_{3,\theta}\mathcal L_\theta I_3 f=-32/75
                                                        \tag{21}
\]

at this rational configuration. Ordered word derivatives verify both
orientations directly in the native algebra. The response is smooth
and nonzero on an open set, so its reference squared norm is positive.
Since L_(0,Q)^(-1) is strictly positive, M_(3,0)[f] has a strictly
positive theta^2 coefficient. Thus at each fixed such volume the
adapted physical cut has nonzero memory for all sufficiently small
nonzero theta. No common punctured neighborhood over all volumes is
claimed from this coefficient alone. Gauge symmetry and nonzero
memory coexist in this very model.

## YM64-T6: clock covariance and the remaining continuum obligations

For a positive physical clock factor kappa, H becomes kappa H and
the vacuum measure and cut stay the same. Direct substitution yields

\[
 F_r^{(\kappa)}(s)=\kappa F_r(s/\kappa),\quad
 \Sigma_r^{(\kappa)}(s)=\kappa\Sigma_r(s/\kappa),\quad
 M_{r,n}^{(\kappa)}=\kappa^{1-n}M_{r,n}.                   \tag{22}
\]

The inherited gap and the guaranteed Taylor disk scale to kappa gamma.
These covariance laws do not choose kappa or a coupling trajectory.
YM63's two-sided physical gap constraints and its explicit
theta=c_B/(c_E g^4), kappa=c_E g^2/a calibration still apply. In
particular this construction cannot extend the proven coupling window
by a change of projection.

What has been closed here is the overlap-dependent hidden inverse,
finite repeated-block composition, and the bounds on the complete
spectral response tower relative to coarse energy. These are more
specific than a claim to control every generated operator in space.
The following still require proof or construction:

1. Uniform spatial decay for the full coarse drift and response forms,
   beyond the local leading coefficient (20).
2. Constructive control of the actual marginal densities mu_r and
   renormalized local readouts along an increasing sequence of blocks.
3. A nontrivial spacing/coupling trajectory and the continuum field
   theory limits and symmetries needed for the four-dimensional target.

A positive uniform gap alone cannot prove item 1. On N sites arranged
on a ring, the reversible matrix
H_N=I-Pi_1+(e_0-e_(N/2))(e_0-e_(N/2))^T has vacuum 1 and gap at
least one, but its coupling between sites distance N/2 apart has
magnitude 1+1/N. This exact counterexample is not the YM Hamiltonian;
it shows why spectral bounds do not supply the missing spatial proof.

Similarly the adapted measure is defined using the actual h. It is an
exact mathematical construction, not an algorithm that has computed
the interacting vacuum in arbitrary volume. The moment positivity
structure does not by itself identify these energy responses with a
thermodynamic lambda or a spacetime curvature tensor.

## Verification

The Python 3.12 certificate checks pointwise native gradient pullbacks,
both ordered derivative signs, the physical coefficient (21), weighted
conditional towers, two routes to Schur compression, noncommuting
positive moment matrices, signed remainder bounds, clock covariance,
static-reset failures and the two counterexamples. Rational operator
fixtures are explicitly labelled and are not spin truncations of YM62.
The written all-content form arguments and finite arithmetic controls
have different roles; they have not been mechanically formalized or
independently peer reviewed.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym64_vacuum_adapted_tower.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym64_vacuum_adapted_tower.py -v
```

This chapter establishes neither a continuum theory on R4 nor a
physical mass prediction. It supplies exact nested cuts and uniform
relative memory estimates inside the existing fixed-spacing gauge model.
