# YM-82 — Nonlinear lattice spectral gaps

## Claim and scope

The proposed nonlinear extension is useful, provided its hypotheses control
the **entire interacting Hamiltonian**, uniformly in the number of lattice
sites. This stage supplies three explicit constructions:

1. A weighted version of the existing YM62 argument gives a full spectral
   gap for nonlinear plaquette potentials in the original compact quantum
   lattice family. The quadratic potential is treated exactly, without
   linearizing its dynamics.
2. A strong-block extension treats actual nonlinear quantum interactions
   of arbitrary fixed internal amplitude, with sufficiently weak bridges
   between disjoint blocks. An explicit cosine model has a volume- and
   frequency-independent gap greater than 1/8320.
3. A conditional block Poincare theorem gives a uniform gap for a separately
   defined nonlinear, gauge-invariant parent Hamiltonian even when its
   ground-state measure has negative Bakry curvature somewhere. This
   illustrates a method for handling strong internal block nonlinearity.

An exact quartic example also has gaps tending to zero. Nonlinearity alone
therefore supplies neither a gap nor a continuum solution.

The first construction is an application and sharpening of YM62's bounded
hyperedge method. The general conditional-gap principle in the third
construction belongs to the established functional-inequality literature;
we give its proof and an explicit native lattice application. No priority
claim is made for these general principles. The new packet records the
specific nonlinear parameter region, strong-block quantum extension,
native examples, Ward extension, failure example, and independently
replayable rational controls.

All quantum lattice claims use the untruncated SU(2) link Hilbert space.
No finite representation cut is substituted for that space. A lattice gap
at fixed spatial spacing is not a proof of the four-dimensional continuum
Yang–Mills existence and mass-gap problem.

## 1. Native normalization and the actual nonlinear Hamiltonian

Let Lambda be a finite periodic three-dimensional cubic lattice with side
at least four. Each link carries SU(2) with its bi-invariant metric for which
the spin-one-half Casimir is 3/4. Let dPhi be product normalized Haar measure,

\[
 L_E=-\Delta=-\sum_e\Delta_e,\qquad
 W_p=\tfrac12\operatorname{Tr}(U_{\partial p})\in[-1,1].
\]

The native identities are

\[
 L_EW_p=3W_p,\qquad
 \Gamma(W_p)=1-W_p^2,\qquad
 \Gamma_e(W_p)=\tfrac14(1-W_p^2)\quad(e\in p).
 \tag{1}
\]

Every link meets four plaquettes. For real smooth potentials Phi_p on
[-1,1], define

\[
 A_\Lambda=L_E-\sum_p\Phi_p(W_p),\qquad
 H_\Lambda=A_\Lambda-E_0,\qquad A_\Lambda h=E_0h.
 \tag{2}
\]

Positivity improvement gives a unique normalized positive ground state h.
Gauge symmetry makes h gauge invariant. Its measure is the **actual**
measure dmu=h^2 dPhi; it is not assumed equal to exp(sum Phi_p) dPhi/Z.
For physical Hamiltonian kappa H_Lambda, all gaps below are multiplied by
kappa.

For each plaquette put

\[
 c_p=\frac{\max\Phi_p+\min\Phi_p}{2},\qquad
 a_p=\frac{\operatorname{osc}_{[-1,1]}\Phi_p}{2},\qquad
 s=\max_e\sum_{p\ni e}a_p.
 \tag{3}
\]

Subtracting c_p from Phi_p shifts E_0 but changes neither H_Lambda nor h.

## 2. T1 — A gap for actual nonlinear plaquette interactions

**Theorem.** If s<1/2000, set

\[
 \bar\alpha=\tfrac12+1000s,\qquad
 R=\frac{1+\bar\alpha^{-1}}2,\qquad
 \gamma=\frac{\log R}{50}.
 \tag{4}
\]

Then gamma>0 and, on the full link Hilbert space,

\[
 H_\Lambda\ge\gamma(I-P_0),\qquad
 \|e^{-tH_\Lambda}(I-P_0)\|\le e^{-\gamma t}.
 \tag{5}
\]

The same inequality holds after restriction to the gauge-invariant Hilbert
space. The constants are independent of volume and require no
representation cutoff.

**Proof.** We expose the weighted change to YM62 T2–T3, keeping its
endpoint-conditioned heat bridges and all-source spectral extraction.
Write the centered interaction as B=sum_p(Phi_p(W_p)-c_p).
Its factor on a temporal row of width h has log-amplitude at most h a_p.
Use 0<h<=1, r=ceil(9/h), and blocks of length ell=Jr. The free heat
minorization is epsilon=4/5 and 9<=hr<10.

For a target link e, the total log-tilt on k rows is bounded by
h k sum_{p containing e} a_p. The normalization comparison gives
total-variation error at most 2hks. The free endpoint coupling plus this
comparison gives the endpoint disagreement budget

\[
 D(k)\le r/\epsilon+4hsk^2.
 \tag{6}
\]

Changing an exterior link f on one time row costs at most

\[
 4hk\sum_{p\ni e,f}a_p.
 \tag{7}
\]

The temporal endpoint contributions to the outgoing block influence are
at most 2r/(epsilon ell)+8sh ell. Summing the spatial contributions uses

\[
 \sum_{p\ni e}(|p|-1)a_p=3\sum_{p\ni e}a_p\le3s,
\]

so those contributions are at most 12sh ell. Consequently

\[
 \alpha\le\frac{2r}{\epsilon\ell}+20sh\ell
 \le\frac{2}{\epsilon J}+200sJ.
 \tag{8}
\]

J=5 gives alpha<=baralpha<1. The admissible temporal weight R in (4)
satisfies R baralpha=(1+baralpha)/2<1. The same bridge contraction
therefore controls every bounded, centered, finite-volume readout with
rate at least log(R)/50 in continuous time.

Here the smooth nonlinear potential requires a refinement check. Define
D=sum_p Lip(Phi_p) and b=sum_p a_p, both finite at fixed volume.
The free increment bound

\[
 \sup_x\mathbb E_x|W_p(X_u)-W_p(x)|^2\le6u
 \tag{9}
\]

follows from K_uW_p=e^{-3u}W_p, the variance identity, and Gamma(W_p)<=3.
For example, the variance is at most 1-e^{-6u}; adding the squared drift
gives at most 2(1-e^{-3u})<=6u. The stronger identity (1) is consistent
with this bound. Minkowski and Lipschitz continuity give the symmetric
kernel Schur estimate

\[
 \|[K_u,M_B]\|\le D\sqrt{6u}.
 \tag{10}
\]

Thus the YM62 insertion-packet argument converges in operator norm at
fixed volume, with error bounded by D t e^{bt} sqrt(18 delta) for the
declared mesh refinement. The derivative budget D replaces the linear
potential coefficient only in this refinement estimate; it does not
enter (8). Volume-dependent refinement prefactors do not affect the
uniform gap rate, since the refinement is first performed at each
fixed volume.

For completeness, the all-source step is essential. At each discrete
time step let h_h denote its positive transfer vacuum. The stationary
bridge estimate applies to bounded centered f in L2(h_h^2 dPhi).
Such f h_h are dense in that vacuum complement. Even stationary moments
remove source-dependent correlation prefactors; the positive spectra
of X*X and XX* agree. Norm refinement then transports the isolated
projection to the continuous-time vacuum h. The spectral theorem
gives (5), as in YM62 T3.
A bound on only one Wilson correlation would not justify this step.
Finally the full-space bound restricts to the gauge-invariant subspace
because the vacuum and Hamiltonian are gauge invariant. QED.

### Exact quadratic region

For Phi(w)=theta w+lambda w^2, let t=|theta| and l=|lambda|.
Reflection of w and an overall sign reduce its oscillation to that of
t w+l w^2. Checking the endpoints and the possible vertex gives

\[
 a(\theta,\lambda)=
 \begin{cases}
 t,&t\ge2l,\\[2mm]
 (t+2l)^2/(8l),&0\le t<2l.
 \end{cases}
 \tag{11}
\]

The first branch includes l=0. On the homogeneous cubic lattice s=4a,
so the strict sufficient region is a<1/8000. Pure quadratic coupling
therefore permits |lambda|<1/4000.

A genuinely mixed cell is

\[
 \theta=\lambda=\frac1{16000},\quad
 a=\frac9{128000},\quad s=\frac9{32000},\quad
 \bar\alpha=\frac{25}{32},\quad R=\frac{57}{50}.
\]

It has

\[
 \boxed{\operatorname{gap}H_\Lambda
 \ge\frac{\log(57/50)}{50}>\frac7{2675}.}
 \tag{12}
\]

Indeed log R>2(R-1)/(R+1)=14/107, and R baralpha=57/64.
The cruder triangle estimate |theta|+|lambda|=1/8000 lies exactly on
the excluded boundary; centering at the actual extrema proves this cell.
The result concerns the nonlinear operator (2) itself.

There is also a useful exact limit to this improvement:
a(theta,lambda)>=|theta| for every lambda, since the two endpoint
values differ by 2theta. An even quadratic term cannot cancel that
oscillation. In the programme's physical scaling from YM76,

\[
 \theta(a)=\frac{c_B}{c_Eg(a)^4},\qquad
 \kappa(a)=\frac{c_Eg(a)^2}{a}.
\]

The weak-coupling trajectory g(a)->0 sends theta(a) to infinity.
It therefore leaves this sufficient small-oscillation region,
whatever quadratic coefficient is added. That is a limitation of
this criterion, not a proof that the physical gap closes.

## 3. T2 — Strong nonlinear quantum blocks with weak bridges

The small-oscillation condition in T1 is not necessary when strong
interactions can be placed inside disjoint blocks. This extension
keeps the **actual quantum Hamiltonian** and its actual vacuum.

Partition the links into disjoint blocks B, each containing at most
r_0 links. Let U_B be smooth real internal potentials with
osc(U_B)<=B_0. Let v_X be smooth real multiplication operators with
|v_X|<=1, each depending on at most q blocks, with every block meeting
at most d such terms. Here r_0,q,d are positive integers; unused
bounds may be enlarged. Set

\[
 A=\sum_B(L_{E,B}+U_B)-u\sum_Xv_X,\qquad
 \epsilon=e^{-10B_0}(19/21)^{r_0}.
 \tag{13}
\]

Choose any integer J>=4/epsilon. If

\[
 |u|<\frac1{80d(q+1)J},\qquad
 \bar\alpha_B=\frac2{\epsilon J}
               +40d(q+1)|u|J<1,
\]

then, with E_0 the actual ground energy and P_0 its projection,

\[
 A-E_0\ge\gamma_B(I-P_0),\qquad
 \gamma_B=\frac1{10J}
  \log\left(\frac{1+\bar\alpha_B^{-1}}2\right)>0.
 \tag{14}
\]

The constants are independent of the number of blocks. They depend
on the bound B_0; no uniformity in arbitrarily increasing internal
barrier amplitude is asserted. If all the potentials are gauge
invariant, the full-space bound also holds in the Gauss sector.
More generally the proof gives the same bound for any integer J>=1
for which the displayed alpha_B is below one. The conditions
J>=4/epsilon and the stated weak-coupling threshold are a convenient
sufficient gate. Any 1<R<1/alpha_B yields the rate log(R)/(10J).

**Proof.** For A_B=L_{E,B}+U_B, let h_B>0 and E_B be the
finite-block ground state and energy. Their existence follows from
elliptic positivity on the compact block, without assuming a
quantitative block gap. At coarse times 9<=t<=10 the Feynman–Kac
kernel obeys

\[
 e^{-t\max U_B}(19/20)^{|B|}
 \le K^B_t(x,y)\le
 e^{-t\min U_B}(21/20)^{|B|}.
\]

Consequently

\[
 \frac{\inf K^B_t}{\sup K^B_t}
 \ge e^{-10\operatorname{osc}U_B}(19/21)^{|B|}
 \ge\epsilon.
 \tag{15}
\]

Let h_0=prod_B h_B and dmu_0=h_0^2 dPhi. The exact reference
transformation is

\[
 h_0^{-1}\left[\sum_B(A_B-E_B)\right]h_0
       =\mathcal L_{\mu_0}.
 \tag{16}
\]

Multiplication by each v_X is unchanged by this unitary
identification. The reference measure is a product of actual
block-vacuum measures, not a guessed Gibbs measure.

The ground-state factors telescope in endpoint-conditioned
reference bridges. To see the required overlap directly, if
m<=K(x,y)<=M and F(y)>=0 is any common future weight with
0<integral F dPhi<infinity, then

\[
 \frac{K(x,y)F(y)}{\int K(x,z)F(z)\,d\Phi(z)}
 \ge\frac mM\,\frac{F(y)}{\int F(z)\,d\Phi(z)}.
 \tag{17}
\]

Thus any two endpoint-conditioned laws share a common minorant of
mass at least epsilon. The h_B factors can be absorbed into F, so
no bound on the individual ground states is needed.

Apply the same temporal-block construction as YM62, now with the
disjoint blocks as sites. Its reference endpoint loss is
2/(epsilon J); the weak interaction loss is bounded by
40d(q+1)|u|J. This proves (14) for the discrete approximants by
the all-source extraction used in T1.

For completeness the reference is no longer the free heat, but
the needed norm refinement still holds. Its Dirichlet form has
carre du champ Gamma(f)=|grad f|^2, independent of the unknown
block scores. With P_t=exp(-t L_mu0), the form commutator and
spectral calculus give

\[
 |\langle g,[\mathcal L_{\mu_0},f]k\rangle|
 \le\|\sqrt{\Gamma f}\|_\infty
  \big(\|\mathcal L_{\mu_0}^{1/2}g\|\,\|k\|
       +\|g\|\,\|\mathcal L_{\mu_0}^{1/2}k\|\big),
\]
\[
 \|[P_t,M_f]\|
 \le4\sqrt{\frac{t}{2e}}\,
                 \|\sqrt{\Gamma f}\|_\infty.
 \tag{18}
\]

Indeed insert P_s and P_(t-s), use
||L^(1/2)P_s||<=(2es)^(-1/2), and integrate from zero to t.
The smooth weak potentials have finite gradient bounds at fixed
volume. Equation (18) therefore supplies the same O(sqrt(mesh))
insertion-packet norm convergence. It transports the isolated
vacuum projection and gap to the actual continuous-time Hamiltonian.
Constants in this refinement need not be uniform in the number of
blocks. QED.

### Explicit strong nonlinear plaquette model

On the even cubic lattice use selected xy plaquettes rooted at
(2i,2j,k) as four-link blocks; all unused links are singleton blocks.
Every nonselected face meets four distinct blocks. A selected block
meets twelve nonselected faces, and a singleton meets four. Thus
r_0=4, q=4, and d=12 are valid.

For every positive integer n take the actual Hamiltonian

\[
 A^{\rm strong}_n
 =L_E-\frac18\sum_{b\ {\rm selected}}\cos(nW_b)
       -\frac1{1228800}\sum_{p\ {\rm nonselected}}W_p.
 \tag{19}
\]

Its internal oscillation is at most B_0=1/4. Since
e^{5/2}<16 and (19/21)^4>1/2, (15) gives epsilon>1/32.
Choose J=128 and use the certified floor 1/32. Then

\[
 \bar\alpha_B\le\frac12+\frac14=\frac34,\qquad
 R=\frac76,\qquad
 \boxed{\operatorname{gap}(A^{\rm strong}_n-E_0)
       \ge\frac{\log(7/6)}{1280}>\frac1{8320}.}
 \tag{20}
\]

The inequality log(7/6)>2/13 proves the rational floor.
It holds for every even side L>=4 and every n>=1, on both the
full and gauge-invariant Hilbert spaces.

This is genuinely outside the small-amplitude argument T1.
For n>=4 a selected cosine ranges over [-1,1], giving a link load
s>=1/8, far above T1's threshold 1/2000. The unknown actual
ground state has not been replaced by an explicit trial state.
No uniform derivative bound on the internal cosine is required;
the commutator (18) uses only the weak potentials' gradients.

The price is explicit: the bridge coefficient is very small, and
the interactions distinguish selected plaquettes. Homogeneous
large coupling on every plaquette, and the original continuum
trajectory, are not covered by this strong-block theorem.

## 4. T3 — Conditional block gaps beyond positive curvature

This is a spectral-gap version of the classical conditional
functional-inequality principle. Otto and Reznikoff's preprint,
Remark 4, discusses the spectral-gap result of Ledoux; its Theorem 1 is
the logarithmic-Sobolev analogue. Their discussion also covers
Riemannian configuration spaces. We include a direct proof to make
the assumptions used here explicit.

Let the configuration space be a finite product of compact connected
Riemannian manifolds without boundary, with blocks B grouping its
factors. Let dPhi be product normalized Riemannian volume and
dmu=exp(S)dPhi/Z with S smooth. Define the positive reversible generator

\[
 \mathcal L=-\Delta-\nabla S\cdot\nabla=\sum_B\mathcal L_B.
 \tag{21}
\]

Suppose every conditional distribution of block B, uniformly in its
exterior configuration, has Poincare gap rho_B>0. Suppose also

\[
 \|\operatorname{Hess}_{BC}S\|_{\rm op}\le c_{BC}=c_{CB}
 \quad(B\ne C).
\]

Define the symmetric matrix M by M_BB=rho_B and M_BC=-c_BC.
If M>=delta I with delta>0, then

\[
 \operatorname{gap}\mathcal L\ge\delta.
 \tag{22}
\]

**Proof.** Conditional integration by parts and commuting derivatives
on distinct product blocks give

\[
 \|\mathcal Lf\|_2^2
 =\sum_B\|\mathcal L_Bf\|_2^2+
 2\sum_{B<C}\left[
 \mu\|\nabla_B\nabla_Cf\|_{\rm HS}^2
 -\mu\,\operatorname{Hess}_{BC}S(\nabla_Bf,\nabla_Cf)\right].
 \tag{23}
\]

Conditional spectral calculus gives
||L_B f||_2^2>=rho_B mu|grad_B f|^2.
Drop the nonnegative mixed derivative terms, and put
v_B=(mu|grad_B f|^2)^(1/2). Cauchy–Schwarz yields

\[
 \|\mathcal Lf\|_2^2\ge v^TMv
 \ge\delta\sum_Bv_B^2=\delta\langle f,\mathcal Lf\rangle.
 \tag{24}
\]

Smooth functions form a core. The spectral theorem excludes spectrum
in (0,delta); connectedness makes the kernel the constants.
This proves (22). No derivative of a conditional expectation was
assumed bounded. QED.

For the actual Hamiltonian (2), the appropriate S is **2 log h**.
Establishing the conditional rho_B and mixed Hessian bounds for that
actual S would directly prove its gap. Taking S=sum Phi_p without
proving it equals 2 log h would change the model.

### A nonlinear gauge-invariant parent with a uniform gap

Let the side L be even and at least four. Select the xy plaquettes
rooted at (2i,2j,k); their four-link boundaries are disjoint.
Use those four-link sets as blocks and all remaining links as
singleton blocks. Every nonselected plaquette uses one link in each
of four distinct blocks. A selected block meets twelve nonselected
plaquettes; a singleton meets four.

For alpha,j>=0 and a positive integer n, define

\[
 S=\alpha\sum_{b\ {\rm selected}}\cos(nW_b)
       +j\sum_{p\ {\rm nonselected}}W_p,
 \qquad h_S=\frac{e^{S/2}}{\|e^{S/2}\|_2},
\]
\[
 K_S=L_E-\tfrac12L_ES+\tfrac14\Gamma(S).
 \tag{25}
\]

Direct differentiation proves K_Sh_S=0, and the ground-state transform
identifies K_S with (21). Thus it is nonnegative, local, gauge
invariant, and has exactly the indicated vacuum. It is a separately
defined parent Hamiltonian, not a replacement for (2).

The conditional oscillation of S is at most 2alpha+24j on a selected
block and 8j on a singleton. Product Haar on any block has gap 3/4.
The elementary density comparison therefore gives the common
conditional lower bound

\[
 \rho\ge\tfrac34e^{-2\alpha-24j}.
 \tag{26}
\]

To see the comparison, center a test function by its Haar mean:
Var_mu(f)<=e^{max S} Z^{-1} Var_Haar(f), while the Dirichlet
energy under mu is at least e^{min S} Z^{-1} times Haar energy.

For two different links of one face the native mixed Hessian of W has
operator norm at most 1/4. Since each bridge face meets four distinct
blocks, the sum of mixed block norms in a selected-block row is at
most 12*3*j/4=9j; singleton rows are smaller. Gershgorin and (22)
give

\[
 \boxed{\operatorname{gap}K_S
 \ge\tfrac34e^{-2\alpha-24j}-9j.}
 \tag{27}
\]

At alpha=1/8 and j=1/192, use e^{-3/8}>=1-3/8=5/8 to obtain

\[
 \boxed{\operatorname{gap}K_S\ge\frac{27}{64}}
 \quad\hbox{for every even }L\ge4\hbox{ and every }n\ge1.
 \tag{28}
\]

The internal nonlinear frequency n may grow without weakening this
bound: conditional oscillation, not the size of internal derivatives,
controls (26).

### Why a curvature-only test fails on this example

For n>=4 choose a selected plaquette with W=pi/n and a unit vector
X parallel to its gradient, supported on its four-link block.
Equation (1) and sin(pi)=0 give

\[
 \operatorname{Hess}[\alpha\cos(nW)](X,X)
   =\alpha(n^2-\pi^2).
 \tag{29}
\]

The native face estimate

\[
 |\operatorname{Hess}W_p(X,X)|
 \le\tfrac14\Big(\sum_{e\in p}|X_e|\Big)^2
 \le\sum_{e\in p}|X_e|^2
\]

bounds the nonselected contribution on this selected-block support
below by -3j: each of its links belongs to three nonselected faces.
Other selected blocks have disjoint support. Since Ric(X,X)=1/2,

\[
 \operatorname{Hess}S(X,X)
 \ge\alpha(n^2-\pi^2)-3j
 >\frac{47}{64}>\frac12
 \tag{30}
\]

at the parameters of (28), using pi^2<10. Thus Ric-Hess S has a
negative direction, while the global gap remains at least 27/64.
Local lack of convexity does not prevent this block method from working.

### The remaining bridge to the original model

For V=sum_p W_p and S=tV, equation (25) reads

\[
 K_{tV}=L_E-\tfrac32tV+\tfrac14t^2\Gamma(V).
\]

Putting t=2theta/3 gives

\[
 K_{tV}=(L_E-\theta V)+\frac{\theta^2}{9}\Gamma(V).
 \tag{31}
\]

The extra potential cannot be discarded in a gap proof. Operator
ordering controls eigenvalues separately, not their difference:
diag(0,1) minus the nonnegative matrix diag(0,1-epsilon) has gap
epsilon. A valid comparison of gaps needs additional estimates.

Likewise, YM81's actual vacuum-Hessian row identity is an averaged
L2 identity. It does not supply the uniform mixed multiplication-form
bound in (23) for arbitrary gradients, which can concentrate on small
sets. It identifies relevant quantities but does not close (22) for
the weak-coupling actual Yang–Mills vacuum. Full-space conditional
bounds can be restricted to the Gauss sector after proving the
global inequality; a conditional bound on an unspecified physical
subclass cannot be inserted without an invariance argument.

## 5. Actual nonlinear Ward identity and perturbative normalization

The nonlinear potential also has an exact response identity. Put
chi=log h for the actual ground state in (2), and

\[
 \tau_e=\mu|\nabla_e\chi|^2,\qquad
 R_e=\sum_f\mu\|\operatorname{Hess}_{f,e}\chi\|_{\rm HS}^2.
\]

The link chain rule is

\[
 \Delta_e\Phi_p(W_p)=\tfrac14\big[
 (1-W_p^2)\Phi_p''(W_p)-3W_p\Phi_p'(W_p)\big].
\]

Applying the YM81 invariant-derivative Ward calculation to
B=sum_p Phi_p(W_p), followed by the Ricci correction 1/2, gives

\[
 R_e+\tfrac12\tau_e
 =\tfrac18\sum_{p\ni e}\mu\big[
 3W_p\Phi_p'(W_p)-(1-W_p^2)\Phi_p''(W_p)\big].
 \tag{32}
\]

For a homogeneous quadratic potential on the periodic cubic lattice,
set m=mu W and Y=W^2-1/4. Then

\[
 R+\tfrac12\tau=\tfrac32\theta m+4\lambda\,\mu Y.
 \tag{33}
\]

These are actual-vacuum identities; no Gibbs ansatz is used.

For a fixed finite lattice, native Haar moments and electric
eigenfunctions give

\[
 L_EY_p=8Y_p,\quad \Phi(W_p^2)=\tfrac14,\quad
 \Phi(Y_p^2)=\tfrac1{16},\quad \Phi(W_pY_p)=0.
\]

Here Phi without subscript denotes Haar expectation only.
Different plaquettes are orthogonal for these W and Y insertions,
by integration over a private edge. Analytic perturbation near
(theta,lambda)=(0,0) yields

\[
 h=1+\tfrac\theta3\sum_pW_p+\tfrac\lambda8\sum_pY_p
     +O_\Lambda((|\theta|+|\lambda|)^2),
\]
\[
 \frac{E_0}{|P|}=-\tfrac\lambda4-\tfrac{\theta^2}{12}
                    -\tfrac{\lambda^2}{128}
                    +O_\Lambda((|\theta|+|\lambda|)^3),
\]
\[
 \mu W=\tfrac\theta6+O_\Lambda(2),\qquad
 \mu Y=\tfrac\lambda{64}+O_\Lambda(2),
\]
\[
 \tau=\tfrac{\theta^2}{12}+\tfrac{\lambda^2}{128}
       +O_\Lambda(3),\qquad
 R=\tfrac{5\theta^2}{24}+\tfrac{15\lambda^2}{256}
       +O_\Lambda(3).
 \tag{34}
\]

The absence of a mixed theta-lambda quadratic term follows from Haar
orthogonality. These fixed-volume expansions check native
normalizations. Their remainders are not asserted uniform in volume,
and the expansions are not used as the proof of (5).

## 6. T4 — Exact nonlinear gaps can collapse

On a single four-link closed plaquette define, for t>=0,

\[
 F_t(w)=t^2w^4+(4t-t^2)w^2,\qquad A_t=L_E-F_t(W).
 \tag{35}
\]

Equations (1) show directly that

\[
 h_t\propto e^{tW^2/2},\qquad A_th_t=-t h_t,\qquad
 A_t+t=K_{S=tW^2}.
 \tag{36}
\]

The positive ground state is exact, and its measure is
dmu_t=Z_t^{-1}e^{tW^2}dPhi.

The gauge-invariant odd form-domain test
f(W)=max(-1,min(1,2W)) has mu_t f=0. Its derivative is supported
on |W|<1/2, so

\[
 \mu_t\Gamma(f)\le4e^{t/4}/Z_t.
\]

On |W|>=sqrt(3)/2 its squared value is one. The Haar cap obeys

\[
 \Phi(|W|\ge\sqrt3/2)
 =\frac4\pi\int_0^{\pi/6}\sin^2r\,dr
 \ge\frac4\pi\int_0^{\pi/6}(2r/\pi)^2\,dr
 =\frac2{81}.
 \tag{37}
\]

Thus its variance is at least (2/81)e^{3t/4}/Z_t. Rayleigh's
principle in the physical gauge-invariant Hilbert space proves

\[
 \boxed{\Delta_t\le162e^{-t/2}\longrightarrow0.}
 \tag{38}
\]

Smooth approximation justifies the clipped test in the closed
Dirichlet form. The usual Gauss constraint does not exclude it;
an additional W-even restriction would be a different spectral
question. Formula (36) is a one-active-plaquette construction,
not an identity for an overlapping all-plaquette sum.

For an independent entirely rational cap check, both intervals
[7/8,15/16] and its negative lie in the cap. The Haar density
(2/pi)sqrt(1-w^2) on them exceeds (5/8)(11/32).
Their total length is 1/8, giving cap>55/2048>2/81.
This is the certificate's finite rational control of (37).

This obstruction does not disprove a continuum Yang–Mills gap.
It shows why a theorem for arbitrary nonlinear lattice systems must
retain quantitative hypotheses, including control of conditional
barriers or the actual vacuum coupling.

## 7. What this advances and what remains

| Result | Established here | Remaining obligation |
| --- | --- | --- |
| Actual nonlinear plaquette family | Uniform-volume quantum gap for s<1/2000; exact quadratic region | Reach the original weak-coupling continuum trajectory |
| Actual strong nonlinear blocks | Uniform-volume and frequency gap greater than 1/8320 in an explicit cosine model | Handle homogeneous large coupling instead of weak bridges |
| Conditional block method | Full gap from uniform conditional gaps and a positive coupling matrix | Prove those estimates for S=2 log h_YM |
| Explicit nonconvex parent | Gap at least 27/64 for all volumes and internal frequencies n | Compare to the original Hamiltonian without losing the gap |
| Actual nonlinear response | Ward identity and fixed-volume quadratic coefficients | Uniform estimates for the nonlinear continuum vacuum |
| Quartic native model | Exact collapsing-gap upper bound | No universal gap conclusion follows merely from nonlinearity |
| Four-dimensional continuum | No new claim of construction or mass gap | Nontrivial continuum fields, physical scale control, covariance and positive physical gap |

For the present programme the useful next analytic target is precise:
choose blocks for the **actual** vacuum and prove a uniform lower
conditional gap together with cross-block coupling small enough for
(22). Internal nonlinearity can then be handled inside each block.
A small-amplitude nonlinear lattice theorem alone does not establish
that target along the desired continuum scaling.

## 8. Reproduction and evidence limits

The companion Python 3.12 certificate uses exact rational operations
and outward transcendental brackets. It checks the half-oscillation
region, weighted contraction constants, strong-block kernel minorization,
native block geometry,
conditional comparison budgets, loss of positive curvature,
nonlinear Ward/Haar algebra, exact parent factorization, and the
double-well upper bound. Independent tests include invalid-parameter,
scope-tamper and source-tamper controls.

The manifest preserves 309 upstream inputs and records four local
sources: this manuscript, the certificate module, the independent
tests, and the dedicated workflow. Result hashes certify the
reproduced finite arithmetic and declared source bytes. They are
not a machine formalization of the analytic infinite-dimensional
proofs above.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym82_nonlinear_lattice_gap.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym82_nonlinear_lattice_gap.py -v
~~~

## References

- [YM53: anisotropic interacting gap](YM53_ANISOTROPIC_INTERACTING_GAP.md),
  for endpoint-conditioned temporal block comparison.
- [YM62: gauge lattice gap](YM62_GAUGE_LATTICE_GAP.md), especially T2–T3,
  for the full SU(2) heat, bounded hyperedge argument, all-source
  extraction, and continuous-time refinement.
- [YM81: local third moments and vacuum response](YM81_LOCAL_THIRD_MOMENT_AND_VACUUM_RESPONSE.md),
  for the actual-vacuum Ward identity and the unresolved continuum
  zero-momentum coefficient.
- [YM76: local defects and weak-coupling dynamics](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md),
  for the original physical coupling and clock used in the continuum
  comparison.
- Felix Otto and Maria G. Reznikoff, *A new criterion for the
  logarithmic Sobolev inequality and two applications*, J. Funct.
  Anal. 243 (2007), 121–157.
  [Primary preprint](https://webdoc.sub.gwdg.de/ebook/serien/e/sfb611/274.pdf),
  especially Remark 4 for the earlier spectral-gap criterion.
