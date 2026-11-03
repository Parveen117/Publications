# Perturbations of the coupled native cosmological solution

Author: Monty Dabas. Research development R7, 30 September 2026.

The exact solution in [CS](COUPLED_COSMOLOGICAL_SOLUTION.md) now admits two
distinct stability results: bounded linear perturbations of the **coupled
homogeneous FLRW sector**, and a wavelength-independent energy bound for
**linear matter perturbations on the prescribed CS geometry**. The latter is
not a stability theorem for spatially varying Einstein–matter perturbations.
For the CS benchmark, starting at t*=log(2), the future amplification of the
rescaled matter norm on that fixed geometry is at most sqrt(2).

All statements use the unchanged classical commuting-field action of
[SM](NATIVE_MASSIVE_MATTER.md), its zero-Holst connection solution, the neutral
sector g=0, and fixed kappa>0, Z>0, m>=0. The principal positive-Lambda results
assume the expanding CS branch and an initial time t*>0. No quantum statistics,
empirical calibration, nonlinear stability, or singularity-crossing result is
inferred. Sources, exact controls and limitations are recorded in PS-8.

## PS-1. Homogeneous completion with all eight matter components

Retain the SM matrices H, Ical, Jcal and gamma^a, with signature (-,+,+,+).
Here Ical and Jcal square to -I, H=-Ical gamma^0, H^2=I, and
gamma^0=Ical H. Ical commutes with H, Jcal and the spatial matrices
Acal_i=gamma^0 gamma^i; H anticommutes with Jcal and Acal_i. Use

\[
v=a^3,\qquad \chi=v^{1/2}\Psi,\qquad
\alpha=\frac{3\kappa Z}{16},\qquad \nu=\frac\alpha v,
\]
\[
\Sigma=\chi^TH\chi,\quad
\Pi=\chi^TH\mathcal J\chi,\quad
B=\chi^T\mathcal I\mathcal J\chi,\quad N=\chi^T\chi.
\]

The symbol alpha here is the stated contact coefficient, not a prediction of
the electromagnetic fine-structure constant. Without imposing H chi=chi,
the homogeneous CS matter equation becomes

\[
\boxed{\dot\chi=-\gamma^0\big[(m+\nu\Sigma)\chi
                         +\nu\Pi\mathcal J\chi\big].}
\tag{PS.1}
\]

Both gamma^0 and gamma^0 Jcal are skew real matrices. Hence N is conserved.
For the ordered symmetric matrices (H,H Jcal,Ical Jcal), the commutators are

\[
[\,(H,H\mathcal J,\mathcal I\mathcal J),\gamma^0\,]
 =(0,-2\mathcal I\mathcal J,2H\mathcal J),
\]
\[
[\,(H,H\mathcal J,\mathcal I\mathcal J),\gamma^0\mathcal J\,]
 =(2\mathcal I\mathcal J,0,-2H).
\]

Differentiating a quadratic form along (PS.1) gives the universal equations

\[
\boxed{\dot\Sigma=-2\nu\Pi B,\qquad
       \dot\Pi=2(m+\nu\Sigma)B,\qquad
       \dot B=-2m\Pi.}
\tag{PS.2}
\]

Consequently Sigma^2+Pi^2+B^2 is also conserved. Each of the three quadratic
form matrices is a symmetric involution, and they anticommute pairwise. For
any real vector q, their linear combination has operator norm |q|. Applying
this with q=(Sigma,Pi,B) gives

\[
\Sigma^2+\Pi^2+B^2\leq N^2.
\tag{PS.3}
\]

This proves a bound on each homogeneous solution's bilinears. It does not
prove that the distance between two solutions, or a linearized norm, is
conserved: the skew generator in (PS.1) itself depends on the solution.

The effective stress remains isotropic for every homogeneous chi:

\[
\boxed{\rho=\frac{Zm\Sigma}{2v}
           +\frac{3\kappa Z^2}{64v^2}(\Sigma^2+\Pi^2),\qquad
       p=\frac{3\kappa Z^2}{64v^2}(\Sigma^2+\Pi^2).}
\tag{PS.4}
\]

To see this directly, substitute the homogeneous Levi-Civita derivative and
(PS.1) into the Hilbert stress in CS.8. The matrices H gamma^i, H Acal_i and
H Acal_i Jcal are skew, so the candidate mixed components vanish for real
commuting chi. Symmetrized off-diagonal spatial Clifford products vanish;
the diagonal spatial kinetic terms give the same pressure. On shell the
kinetic density is Z[m S+alpha(S^2+P^2)]/2, with S=Sigma/v, P=Pi/v.
Subtracting U in CS.2 gives p in (PS.4), and T_00=U gives rho.

Equations (PS.2) also give
dot(rho)+3h(rho+p)=0 exactly. Thus the full homogeneous completion is

\[
\dot v=3hv,\qquad
3h^2=\Lambda+\kappa\rho,\qquad
2\dot h+3h^2=\Lambda-\kappa p,
\tag{PS.5}
\]

together with (PS.1) and SM's algebraic contorsion
K_abc=kappa Z A_abc/8. The usual constraint defect obeys dot(C)=-3h C
when the spatial equation and continuity hold. No additional homogeneous
rest constraint is required. Isotropy of the metric and effective stress
does not make the polarized torsion rotationally invariant.

This is a nonlinear homogeneous ODE completion, not an explicit solution of
all its trajectories or a theorem for the general inhomogeneous PDE. The
verifier checks two nonrest jets with both Pi and B nonzero against every
original first-order coframe, connection and matter equation, as well as the
effective stress; the universal reduction is the argument above.

## PS-2. Coupled homogeneous bilinear and metric perturbations

Linearize about the CS rest solution

\[
\chi_b=e^{-\mathcal I\theta}\Xi,\qquad H\Xi=\Xi,\qquad
\Xi^T\Xi=N=\frac{2n_0}{Z},\qquad
s(t)=\nu N=\frac{3\kappa n_0}{8v},\qquad M=m+s,\quad \dot\theta=M.
\]

Here Sigma=N, Pi=B=0, and dot(s)=-3hs<0. Equation (PS.2) gives

\[
\delta\dot\Sigma=0,\qquad
\delta\dot\Pi=2M\delta B,\qquad
\delta\dot B=-2m\delta\Pi.
\tag{PS.6}
\]

There is no first-order metric forcing in the last two equations because
their background Pi and B vanish. For m>0, the transverse quadratic energy

\[
E_\perp=m(\delta\Pi)^2+M(\delta B)^2,\qquad
\boxed{\dot E_\perp=\dot s(\delta B)^2\leq0}
\tag{PS.7}
\]

is coercive: E_perp>=m[(delta Pi)^2+(delta B)^2]. Both transverse bilinears
are therefore bounded to the future. For m=0, delta B is constant and
delta Pi=delta Pi*+2 delta B* integral(s dt); PS-6 proves this integral
finite on [t*,infinity) for positive Lambda.

At this order delta Sigma=delta N, since H chi_b=chi_b. Define
delta n_0=Z delta N/2. The scalar stress variation is exactly the tangent
variation of CS's rho=m n+3 kappa n^2/16 and p=3 kappa n^2/16.
At fixed couplings, the expanding constrained metric solutions have
v=n_0 f(t-t_0), so all their linear tangents are

\[
\boxed{\frac{\delta v}{v}=\frac{\delta n_0}{n_0}-3h\delta t_0,
\qquad \delta h=-\dot h\delta t_0,
\qquad \delta n=3hn\delta t_0.}
\tag{PS.8}
\]

The first-order expanding constraint supplies one integration constant t_0
for a given n_0; these two parameter variations exhaust the constrained
scalar tangents. Direct differentiation verifies both
6h delta h=kappa delta rho and
2 delta(dot h)+6h delta h=-kappa delta p.
For positive Lambda, h tends to sqrt(Lambda/3), while delta h and delta rho
decay at least as O(exp(-sqrt(3 Lambda)t)). The relative variation delta v/v
is bounded; the absolute delta v need not tend to zero. A constant spatial
scale normalization can be changed locally by rescaling coordinates; on a
fixed compact spatial manifold it can instead change physical volume.

## PS-3. All eight homogeneous components and metric forcing

The three bilinear equations alone would not account for all matter
perturbations. Write the variation of the volume-rescaled field as
delta chi=exp(-Ical theta) zeta, and introduce the orthogonal projectors

\[
P_\pm=\tfrac12(I\pm H),\qquad
\mathcal P=\frac{\Xi\Xi^T+(\mathcal J\Xi)(\mathcal J\Xi)^T}{N}.
\]

Differentiating the cubic field law, including delta nu=-nu delta v/v,
gives the exact homogeneous linear equation

\[
\boxed{\dot\zeta=C(t)\zeta+s\frac{\delta v}{v}\mathcal I\Xi,\qquad
C=2M\mathcal I P_- -2s\mathcal I\mathcal P.}
\tag{PS.9}
\]

The four vectors Xi, Ical Xi, Jcal Xi, Ical Jcal Xi are orthogonal, each
with squared norm N, and

\[
C\Xi=-2s\mathcal I\Xi,\quad C\mathcal I\Xi=0,\quad
C\mathcal J\Xi=2m\mathcal I\mathcal J\Xi,\quad
C\mathcal I\mathcal J\Xi=-2M\mathcal J\Xi.
\tag{PS.10}
\]

The remaining two-dimensional P_+ complement is stationary; the remaining
two-dimensional P_- complement evolves by the skew rotation 2M Ical.
Both complements are invariant. The Jcal Xi/Ical Jcal Xi block is precisely
(PS.6): if its coefficients are c,d, then delta Pi=-2Nc and delta B=2Nd.

For the amplitude/phase part zeta=A Xi+F Ical Xi, A is constant and
dot(F)=-2s A+s delta v/v. Since delta n_0/n_0=2A, the metric tangent gives

\[
\dot F=\dot s\delta t_0,\qquad
F(t)=F(t_*)+[s(t)-s(t_*)]\delta t_0.
\tag{PS.11}
\]

Thus the coupled metric variation cancels the amplitude-induced phase shear
one would obtain by freezing v in this homogeneous calculation. Together,
(PS.7), (PS.8), (PS.10) and (PS.11) account for all eight matter directions
and prove future boundedness of zeta and delta v/v on [t*,infinity).
The physical field variation is

\[
\delta\Psi=v^{-1/2}e^{-\mathcal I\theta}
 \left(\zeta-\tfrac12\frac{\delta v}{v}\Xi\right).
\]

It has bounded size relative to the background amplitude and decays in
absolute component size with v^(-1/2). Contorsion follows by differentiating
SM's algebraic solution; it has no independent propagating connection mode
in this zero-Holst sector. This is a constrained **linear homogeneous FLRW
stability statement**, with neutral parameter/orientation modes retained.
It does not include Bianchi anisotropy, inhomogeneous metric perturbations,
nonlinear orbital stability, or a positive quantum Hamiltonian.

## PS-4. Spatial gradients do not preserve the rest sector

The same native matrices give {H,Acal_i}=0. For a spatial covector k,
Acal(k)=sum_i k_i Acal_i satisfies Acal(k)^2=|k|^2 I. Therefore

\[
P_+\mathcal A(k)P_+=0,\qquad
L(k)=P_-\mathcal A(k)P_+,\qquad
L(k)^TL(k)=|k|^2P_+.
\tag{PS.12}
\]

For nonzero k the leakage map has rank four. Every nonzero rest-sector
plane-wave amplitude acquires a P_- component from its gradient. The
homogeneous rest truncation is consequently not a closed spatial matter
perturbation system.

For example, the background derivative

\[
\frac{dp}{d\rho}=\frac{3\kappa n/8}{m+3\kappa n/8}=\frac{s}{M}
\]

does not, by itself, derive the sound speed of a closed cosmological
perturbation system. Such an identification would require a justified
closure including the retained spinor components and gravitational
constraints. The symmetric principal matrices still give the characteristic
cone established in NP; that principal-symbol fact is distinct from a
fluid reduction or a coupled cosmological stability theorem.

## PS-5. Spatial matter energy theorem on prescribed geometry

Now prescribe the CS metric and its background field, and linearize the
effective matter equation only. Set delta g=0 in this separate problem.
Retain every component of delta chi=exp(-Ical theta) zeta. The equation is

\[
\boxed{\partial_t\zeta=\frac1a\sum_i\mathcal A_i\partial_i\zeta
                                 +C(t)\zeta,}
\tag{PS.13}
\]

with C from (PS.9), without its metric forcing. Work on a periodic flat
three-torus, or with square-integrable perturbations on R^3. Initially one
can take smooth periodic/decaying data and extend by density. All norms
below use the Euclidean eight-component norm and comoving spatial measure.

The spatial term integrates to zero in the real L^2 energy pairing because
the Acal_i are symmetric and a depends only on time. Also 2M Ical P_- is
skew. The symmetric part of C is

\[
C_s=\tfrac12(C+C^T)
    =-s(\mathcal I\mathcal P-\mathcal P\mathcal I),\qquad
\operatorname{spec}(C_s)=(s,s,-s,-s,0,0,0,0).
\tag{PS.14}
\]

For a direct spectral proof, restrict to the four-vector frame of PS-3.
On each amplitude/phase or transverse pair the symmetric block has
eigenvalues +/-s; on its orthogonal complement it vanishes. Equivalently,
C_s^2=s^2 times the projector onto that four-plane. Xi-Ical Xi and
Xi+Ical Xi give explicit eigenvectors with eigenvalues +s and -s. Thus the
instantaneous norm bound is sharp and is not a norm-conservation identity.

Integration by parts and Gronwall's inequality give

\[
\frac{d}{dt}\|\zeta\|_{L^2}^2\leq2s\|\zeta\|_{L^2}^2,\qquad
\boxed{\|\delta\chi(t)\|_{L^2}
 \leq G(t,t_*)\|\delta\chi(t_*)\|_{L^2},\quad
G=\exp\!\left(\int_{t_*}^t s(\tau)d\tau\right).}
\tag{PS.15}
\]

The estimate is uniform in spatial frequency and direction. It also holds
in every spatial Sobolev norm H^r: the coefficients depend only on time, so
the Fourier weight commutes with the equation. Fourier transformation gives
a finite-dimensional linear ODE at each wavevector; smooth coefficients on
every [t*,T] and the uniform estimate construct a unique global future
linear evolution for H^r data, with continuous dependence on the initial
data. This supplies existence as well as an upper bound for this **linear
prescribed-geometry problem**.

Holding delta g=0 is generally inconsistent with the full linearized
Einstein equations when delta T is nonzero. Equation (PS.15) must therefore
not be presented as full inhomogeneous Einstein–matter stability. Nor does
the sharp instantaneous logarithmic norm prove that the integrated bound
is attained by a single solution throughout its history.

## PS-6. Exact gain and a numerical benchmark without numerical integration

For the positive-Lambda CS solution, define

\[
\omega=\sqrt{3\Lambda},\quad z=e^{\omega t},\quad
d=\frac{\kappa m n_0}{2\Lambda},\quad
b=\frac{3\kappa n_0}{4\omega},\quad
v=\frac{(z-1)[(d+b)z-(d-b)]}{2z}.
\]

Here d>=0, b>0 and t>0; a shifted origin is handled by t-t_0. Direct
differentiation of the following logarithm gives d(log G)/dt=s:

\[
\boxed{G(t,t_*)^2=
\frac{(z-1)[(d+b)z_*-(d-b)]}
     {(z_*-1)[(d+b)z-(d-b)]},\qquad
G(\infty,t_*)^2=
\frac{(d+b)z_*-(d-b)}{(d+b)(z_*-1)}<\infty.}
\tag{PS.16}
\]

The finite limit proves the integrability used in PS-2 and global future
boundedness in PS-5 for every fixed t*>0. It does not control the singular
limit t* down to zero.

Use the exact CS benchmark kappa=1, Z=4/3, m=1/2, Lambda=1/3,
n_0=4/3, with d=b=omega=1 and v=exp(t)-1. Choose t*=log(2), so a(t*)=1.
Then

\[
\boxed{\|\delta\chi(t)\|_{H^r}
 \leq\sqrt{2(1-e^{-t})}\|\delta\chi(t_*)\|_{H^r}
 \leq\sqrt2\|\delta\chi(t_*)\|_{H^r}.}
\tag{PS.17}
\]

Since the geometry is fixed in PS-5, delta Psi=v^(-1/2) delta chi, giving

\[
\boxed{\|\delta\Psi(t)\|_{H^r}
 \leq\sqrt2 e^{-t/2}\|\delta\Psi(t_*)\|_{H^r}.}
\tag{PS.18}
\]

| z=exp(t) | Upper bound on rescaled norm gain squared | Upper bound on original field norm gain squared |
|---|---|---|
| 2 | 1 | 1 |
| 3 | 4/3 | 2/3 |
| 5 | 8/5 | 2/5 |
| 9 | 16/9 | 2/9 |
| infinity | 2 | 0 |

These are bounds, not simulated growth factors. The rescaled L^2 norm
equals the physical-volume-weighted L^2 norm of delta Psi. Equation (PS.18)
instead uses comoving measure and includes cosmological dilution. Neither
norm is being interpreted as a quantum transition probability or a
particle-creation rate. The constants are supplied dimensionless benchmark
parameters, not measured physical predictions.

## PS-7. A limit where the relative bound fails

The role of the expanding positive-Lambda background has a concrete
countercontrol. In CS's zero-Lambda, massless contraction,

\[
v=\frac{3\kappa n_0}{4}t,\qquad s=\frac1{2t}.
\]

The coupled homogeneous transverse equations become

\[
\delta B(t)=\delta B_*,\qquad
\delta\Pi(t)=\delta\Pi_*+\delta B_*\log(t/t_*).
\tag{PS.19}
\]

For delta B* nonzero the rescaled transverse perturbation is unbounded.
In logarithmic time L=log(t/t*) its propagator is
I+L [[0,1],[0,0]], acting on (delta Pi,delta B). The initial vector (0,1)
has squared norm 1+L^2; the exact checks use L=0,1,3 and return 1,2,10.
The physical field may still decay like log(t)/sqrt(t), so this is growth
relative to the background amplitude, not absolute field blow-up. A
linear secular term alone is not a proof of nonlinear instability.

Positive Lambda makes this contact integral finite. It is not uniquely
selected by stability: in the massive zero-Lambda CS member, v grows
quadratically at late times and s is also integrable. This countercontrol
tests the stated hypotheses, without promoting them to a physical
selection principle.

## PS-8. Sources, reproducible controls and remaining equations

The written arguments above consume the unchanged SM/CS action and solution.
Their source bytes and Git identities are pinned in [SOURCE_PINS.json](SOURCE_PINS.json).
The new standard-library module [perturbation_stability.py](perturbation_stability.py)
checks exact rational identities and independently differentiates the cubic
law using first-order jets. Its controls include:

- six universal bilinear commutators and the three Clifford involutions;
- two nonrest jets, retaining 32 coframe, 48 connection and 16 matter
  equation components, plus 32 effective stress components;
- 128 independent cubic Jacobian columns, 16 complete eight-component
  decompositions, and 32 real cosine/sine Fourier blocks;
- 16 constrained metric tangents and 32 transverse energy identities;
- 16 exact gain-primitive checks, three rank-four gradient leakage controls,
  the benchmark bounds, and the massless zero-Lambda secular control.

The finite parameter fixtures corroborate the universal written identities;
they are not substitutes for those proofs, a proof assistant, independent
peer review or experimental evidence. All eleven R6 numerical result groups
are retained unchanged when the R7 certificate is regenerated.

Homogeneous classical spinor cosmology is established prior work, for example
C. Armendariz-Picon and P. B. Greene, *Spinors, Inflation, and Non-Singular
Cyclic Cosmologies*, General Relativity and Gravitation 35 (2003), 1637–1658,
[arXiv:hep-th/0301129](https://arxiv.org/abs/hep-th/0301129),
[DOI](https://doi.org/10.1023/A:1025783118888). Expansion-dependent Dirac
existence and lifespan estimates also have established literature; see
K. Yagdjian, *Global in time self-interacting Dirac fields in the de Sitter
space*, Journal of Evolution Equations 22, 22 (2022),
[primary institutional record](https://scholarworks.utrgv.edu/mss_fac/259/),
[DOI](https://doi.org/10.1007/s00028-022-00769-8).
That work's nonlinear theorems have additional hypotheses and are not
imported into this different, self-consistent CS background. Here the
linear estimates and their coefficients are derived explicitly from the
source-bound native action.

The next unresolved stability problem is to include spatial metric and
coframe perturbations, fix their gauge, solve the linearized gravitational
constraints, and estimate the resulting coupled evolution. Anisotropic
homogeneous metric modes are also outside this note. No observable
cosmological spectrum, nonlinear Einstein–matter Cauchy theorem, primitive
action selection, or prediction of G, c, Lambda, fine-structure alpha,
hbar or particle masses is claimed.
