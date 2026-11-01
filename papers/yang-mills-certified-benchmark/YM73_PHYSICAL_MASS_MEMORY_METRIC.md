# YM-73: physical mass from static memory and its normalization

Monty Dabas. 5 October 2026. **Python 3.12 only.**

The next continuum calculation must retain low-energy physical states
and the normalization carried by their hidden components. For a
finite physical observation cut, the two required memory moments are

\[
 M_0=BD^{-1}B^*,\qquad M_1=BD^{-2}B^*,\qquad Z=I+M_1 .
\]

The resulting effective energy operator in physical units is

\[
 \boxed{K=\kappa Z^{-1/2}(C-M_0)Z^{-1/2}.}                 \tag{1}
\]

If the hidden operator has an independently established lower bound
D>=d>0, put h=κd and m_eff=min spec K. The full physical gap satisfies

\[
 \boxed{\frac{m_{\rm eff}h}{m_{\rm eff}+h}
       \le\Delta_{\rm phys}\le m_{\rm eff}.}              \tag{2}
\]

For positive m_eff the lower bound is optimal given only these two
numbers. Thus m_eff->m>0 and h->infinity force the physical gaps to
converge to m. A uniform resolvent estimate below also controls the
reduced dynamics. These implications do not assume a full gap in
order to prove one.

We prove an actual slow-mode enrichment estimate for the YM62 model.
We also prove that the existing YM72 scaling sends the physical gap
to infinity under the declared Hamiltonian calibration. It cannot
itself produce the required finite-mass continuum theory. The
independent estimates needed to apply (1)--(2) on a weak-bare-coupling
trajectory, and the construction of local fields on R4, remain open.

## Setting, cuts and established tools

The actual model remains YM62's SU(2) link Hamiltonian
A=L_E-theta V, its positive vacuum h_vac, and H=A-E>=0.
The vacuum representation is L=U*HU on L2(mu), dmu=h_vac^2 dPhi.
We use the actual vacuum throughout; finite matrix examples below
are labelled algebraic examples and do not compute that vacuum.

T1 retains YM72's full spatial cut and its selected Wilson-loop
observations. T2--T3 introduce a different, explicitly declared cut:
inside the **full physical vacuum-orthogonal Hilbert space**, retain
the span of finitely many smooth gauge-invariant centered states.
Let P be that finite-rank orthogonal projection and Q=I-P on this
space. In the full space the vacuum is retained separately.

Every retained vector lies in D(H). Thus P preserves D(H), the
off-diagonal coupling is bounded, and subtracting it from H gives
a self-adjoint block diagonal operator with hidden part D. Write

\[
 H\big|_{\Omega^\perp_{\rm phys}}=
 \begin{pmatrix}C&B\\B^*&D\end{pmatrix}.                  \tag{3}
\]

The hidden Hilbert space is not truncated. The retained dimension
may increase along a volume/cutoff sequence. The formulae do not
identify this finite physical observation cut with the full
conditional-expectation cut P_r of YM64. In particular, one cannot
invert YM72's sampled memory matrix and thereby obtain (3).

We use the established Schur/Feshbach elimination principle; see
[Griesemer and Hasler, *On the Smooth Feshbach--Schur Map*](https://web.ma.utexas.edu/mp_arc/c/07/07-102.pdf).
The normalization, gap bracket and explicit error estimates here are
derived directly below. No priority claim is made for Schur
elimination, the variational principle or spectral calculus.
The continuum target is the
[Jaffe--Witten problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf):
a nontrivial local four-dimensional quantum Yang--Mills theory with
a finite positive gap for every compact simple gauge group. The
current concrete model is SU(2). A gapped finite matrix or lattice
is not that construction.

## YM73-T1: a microscopic probe restores a quantified slow component

For any elementary plaquette p, set

\[
 w_p=\frac{W_p-\mu(W_p)}{\sqrt{\operatorname{Var}_\mu(W_p)}}.
\]

YM63 proves Var_mu(W_p)>19/100 and
E(W_p,W_p)=mu(1-W_p^2)<=1. Hence

\[
 \|w_p\|=1,\qquad
 \langle w_p,Lw_p\rangle\le C_0:=100/19.                 \tag{4}
\]

For the block clock A_r=L/r and any b>0, positivity of the spectral
measure and Markov's inequality give

\[
 \boxed{\|1_{[0,b]}(A_r)w_p\|^2
       \ge\max\{0,1-C_0/(rb)\}.}                        \tag{5}
\]

The vacuum weight is zero. At b=3/2 the bound is
max(0,1-200/(57r)); it tends to one.

Let J be YM72's isometry for N separated, normalized large-loop
observations. Its generator residual is bounded by

\[
 R=192|\theta|\sqrt{N\exp(384\pi|\theta|rN)},\qquad
 \|(A_r-3)J\|\le R .
\]

YM72 gives ||1_[0,3/2](A_r)J||<=2R/3. Split the inner product at
that spectral projection and use (5). This proves

\[
 \boxed{\|J^*w_p\|\le
       \min\{1,\;2R/3+\sqrt{200/(57r)}\}.}               \tag{6}
\]

No geometric disjointness between p and the large loops is needed.
Their different spectral weights give the estimate.

The joint Gram matrix of J's columns and w_p is
[[I,J*w_p],[(J*w_p)*,1]]. Its eigenvalues lie between 1-t_r and
1+t_r, where t_r is the unclipped expression in (6). For r>=16
and |theta|rN<=10^(-8), directed constants give t_r<0.47, so

\[
 0.53I\le G_{\rm joint}\le1.47I .                        \tag{7}
\]

At r=16 the added microscopic state has at least 89/114 of its
weight below block-clock energy 3/2. It contributes an independent,
well-conditioned direction. Increasing r improves that lower bound.

This is an actual interacting-model enrichment result. One added
probe retains one demonstrated slow component; it does not prove
that every low-energy state has been retained. Unobserved sectors
remain relevant to any full-gap claim.

As a normalization check, the free physical electric gap on a
periodic cubic lattice of side at least four is exactly three.
In the Peter--Weyl edge decomposition, each nontrivial edge label
j>=1/2 contributes j(j+1)>=3/4. A gauge-invariant nonconstant
sector cannot have a vertex incident to only one nontrivial edge:
that representation has no invariant vector at that vertex.
Its nonempty support therefore contains a cycle, of length at
least four. Its electric energy is at least three. A fundamental
plaquette trace attains three. An r-by-r square trace instead
has rate 3r. This argument concerns the free electric model,
not an exact interacting mass value.

## YM73-T2: the first inverse moment fixes the physical normalization

For (3), assume **independently** that D>=dI with d>0. Define
X=D^(-1)B*, S_0=C-M_0, M_0=BD^(-1)B*, M_1=X*X and Z=I+M_1.
The finite matrices are well defined and Z>=I. Completing the
energy square gives S_0>=0 because H>=0.

The harmonic extension of a retained coefficient f is

\[
 \mathcal W f=(f,-Xf),\qquad
 \|\mathcal W f\|^2=\langle f,Zf\rangle,\qquad
 \langle\mathcal W f,H\mathcal W f\rangle
       =\langle f,S_0f\rangle.                           \tag{8}
\]

Thus W=mathcalW Z^(-1/2) is an isometry and
W*(κH)W=K from (1). The Rayleigh principle gives
Delta_phys<=m_eff. Ignoring M_1 would use the wrong norm for these
extended states.

For 0<e<h=κd, the exact Schur complement of κH-e is

\[
 \kappa S_0-eZ-\frac{e^2}{\kappa}
                  BD^{-2}(D-e/\kappa)^{-1}B^* .         \tag{9}
\]

After congruence by Z^(-1/2), its final term R_-(e) obeys

\[
 0\le R_-(e)\le\frac{e^2}{h-e}I,                        \tag{10}
\]

because (D-e/κ)^(-1)<=(d-e/κ)^(-1)I and M_1<=Z.
Consequently the Schur complement is positive whenever
m_eff>=e+e^2/(h-e). Choosing
e=m_eff h/(m_eff+h) establishes (2). If m_eff=0, the harmonic
extension supplies a zero-energy vector and the infimum of the
centered spectrum is zero. No positive-gap conclusion is made.

The lower estimate is optimal without further information about Z.
For scalars take physical hidden energy h and

\[
 \kappa H_x=
 \begin{pmatrix}
 hx^2+m_{\rm eff}(1+x^2)&hx\\
 hx&h
 \end{pmatrix}.                                        \tag{11}
\]

Here Z=1+x^2 and the effective energy is exactly m_eff. The trace
is (h+m_eff)(1+x^2), and the determinant is hm_eff(1+x^2).
As x->infinity, its smaller eigenvalue tends to
hm_eff/(h+m_eff). If x=0 and h>=m_eff, the upper bound is attained.
These are algebraic sharpness examples, not Yang--Mills vacua.

Under a change of clock H->κH, M_0->κM_0 but M_1 and Z remain
unchanged. Static memory has energy units; its first inverse
moment corrects the dimensionless state normalization.

## YM73-T3: two moments control the low-energy resolvent

For s>0, let R_P(s)=P(κH+s)^(-1)P restricted to ran P.
The exact identity

\[
 BD^{-1}B^*-B(D+z)^{-1}B^*
   =zM_1-z^2BD^{-2}(D+z)^{-1}B^*,\qquad z=s/\kappa,
\]

gives

\[
 Z^{1/2}R_P(s)Z^{1/2}
     =[K+sI-R_+(s)]^{-1},\qquad
 0\le R_+(s)\le\frac{s^2}{h+s}I .                       \tag{12}
\]

Thus, uniformly for every positive s,

\[
 \boxed{0\le Z^{1/2}R_P(s)Z^{1/2}-(K+sI)^{-1}
             \le\frac1h I.}                            \tag{13}
\]

**Proof.** Let A=K+sI and r_s=s^2/(h+s). Then A>=sI,
A-R_+>=(s-r_s)I, and inversion reverses matrix order.
The resolvent identity bounds the norm of the positive difference
by r_s/[s(s-r_s)]=1/h. This proves (13), without assuming that
K, Z or the memory matrices commute. QED.

Higher moments have not been discarded without a bound: the exact
remainder in (12) controls all of them at once. When h->infinity,
the first two memory moments suffice for this resolvent estimate.
It does not by itself identify local continuum field operators.

The criterion can accept certified matrix intervals. Suppose
C_L<=C<=C_U and M_j,L<=M_j<=M_j,U, j=0,1, with
Z_L=I+M_1,L>0 and Z_U=I+M_1,U. Put
S_L=C_L-M_0,U and S_U=C_U-M_0,L. If

\[
 \kappa S_L-\ell Z_U\ge0,\qquad \ell\ge0,                \tag{14}
\]

then m_eff>=ell. An explicit nonzero f gives the upper bound

\[
 m_{\rm eff}\le
   \kappa\frac{\langle f,S_Uf\rangle}{\langle f,Z_Lf\rangle}.
                                                               \tag{15}
\]

Combining (14)--(15) with (2) certifies a physical gap interval.
The inequalities require enclosure of the **whole retained pencil**
and a hidden-sector bound; isolated sampled entries are insufficient.
The executable uses exact rational positive-semidefinite checks.

## YM73-T4: the present growing-block trajectory loses finite energy

Use the declared positive calibration constants of YM63:

\[
 \kappa(a)=\frac{c_Eg(a)^2}{a}
     =\frac{\sqrt{c_Ec_B/\theta(a)}}{a},\qquad
 \theta(a)=\frac{c_B}{c_Eg(a)^4}>0 .                     \tag{16}
\]

In YM72's product-budget region, 0<theta rN<=10^(-8),
the actual bare lower gap is at least 1/125. Therefore

\[
 \boxed{\Delta_{\rm phys}(a)\ge
         \frac{80\sqrt{c_Ec_BrN}}a\ \longrightarrow\infty
         \quad(a\downarrow0).}                         \tag{17}
\]

This is stronger than observing that our lower bound is not yet
the desired finite number: (17) excludes a finite physical mass
along this particular calibrated trajectory. It leaves other
coupling trajectories unrestricted. In particular g->0 corresponds
to theta->infinity, outside the present small-theta estimates.

The operator consequence also closes a normalization loophole.
Let u_a be any centered physical states with ||u_a||^2->v>0
and bounded norms. They need not come from uniformly bounded
unrenormalized observables. By (17),

\[
 0\le\langle u_a,e^{-t\kappa H}u_a\rangle
       \le\|u_a\|^2e^{-t\Delta_{\rm phys}(a)}
       \longrightarrow0,\qquad t>0.                   \tag{18}
\]

Together with the nonzero t=0 limit, these correlations cannot
be those of a nonzero vector under a strongly continuous
nonnegative Hamiltonian. Multiplying a probe by a scalar does
not lower its spectral support. Rescaling time changes the
declared physical clock and therefore changes the proposed limit.
Distributional contact terms alone do not evade the stated
Hilbert-space conclusion.

Equation (17) pinpoints why further precision entirely within
YM72's region cannot complete the physical continuum problem.
The needed improvement is a new coupling trajectory with
controlled retained cancellation and hidden spectral separation.

## YM73-T5: a precise conditional continuum passage

Consider a sequence of finite physical lattice models, with their
volume/cutoff limits chosen and observation cuts as in (3).
Assume the following three inputs are established independently:

1. Hidden bounds d_a>0, with h_a=κ_a d_a->infinity, and certified
   effective energies m_eff,a->m in (0,infinity), from (1).
2. A nontrivial candidate continuum theory and a dense centered
   physical state core whose inner products and heat correlations
   are limits of those lattice states. Its positive heat semigroup
   is strongly continuous. The local fields, four-dimensional
   covariance, locality and Yang--Mills identification must be
   constructed, not inferred from a finite matrix.
3. At least one centered local state has approximants u_a with
   ||u_a||^2->v>0 and physical energy forms bounded by E_0<infinity.

Then the continuum Hamiltonian has a finite positive spectral gap.
More precisely, its gap is at least m and at most E_0/v.
Equality to m in the continuum additionally needs low-energy
states attaining that edge to survive the chosen limit.

**Proof.** T2 implies lattice gaps Delta_a->m. For centered
approximating vectors, the full complement bound gives

\[
 |\langle u_a,e^{-tH_a}v_a\rangle|
       \le e^{-t\Delta_a}\|u_a\|\|v_a\|.
\]

Pass to the limit on the stated dense core, then extend by
continuity to the whole centered physical Hilbert space. The
spectral theorem gives H_cont>=m on that space.
This argument needs the full centered core, not a few observed
Wilson loops or their complete moment tower.

For the state in input 3, 1-e^(-tE)<=tE gives
C_a(t)>=||u_a||^2-tE_0. Its limiting correlation is at least
v-tE_0 and is nonzero at sufficiently small positive t.
The increasing difference quotient
(v-C(t))/t as t decreases to zero is bounded by E_0;
spectral monotone convergence places the state in the form
domain with energy at most E_0. Its Rayleigh quotient bounds
the continuum gap above by E_0/v. QED.

This is a passage theorem, not a claim that its open inputs
have been established for four-dimensional Yang--Mills.
The distinction between equality and a lower bound matters:
the lowest lattice states could fail to survive a chosen limit.

## What is established and what the next calculation must supply

| Object | Result in this stage | Remaining continuum input |
|---|---|---|
| Large-loop blindness | An actual microscopic probe retains quantified slow weight | A family that retains every relevant low-energy sector |
| Static memory | Exact subtraction S_0=C-M_0 | Uniform evaluation on a physical weak-coupling trajectory |
| State normalization | Exact Z=I+M_1 and the generalized energy (1) | Certified whole-pencil bounds along that trajectory |
| Hidden dynamics | Sharp gap bracket and uniform error 1/(κd) | A hidden bound independent of the desired full gap |
| Physical scaling | The YM72 calibrated route has escaping spectrum | A different controlled coupling trajectory |
| Four-dimensional theory | Conditional gap passage on a dense core | Local field limits, spacetime covariance and Yang--Mills ultraviolet structure |

In particular, substituting YM64's hidden bound, which was deduced
from an already-known full gap, cannot establish a new weak-coupling
gap. A dense set of observations or an independent hidden estimate
must rule out an additional low-energy sector invisible to all
selected probes. The certificate contains explicit counterexamples
to ignoring either requirement.

The concrete target is now a critical cancellation calculation:
show that κ Z^(-1/2)(C-M_0)Z^(-1/2) has a controlled finite positive
lower spectral edge while κd diverges, and retain a nonzero local
state of bounded physical energy in a genuine continuum construction.
The actual small-theta trajectory fails that target by (17).
This stage supplies the spectral bridge and the slow-probe repair;
it does not claim completion of the 4D continuum mass-gap problem.

## Reproduction

The Python 3.12 certificate checks exact block elimination, the
harmonic-extension norm, noncommuting resolvent remainders, interval
pencils, the optimal scalar boundary, physical-clock covariance,
slow-probe budgets and invisible-sector counterexamples. It also
checks the calibration algebra and preserves all upstream pins.
Finite fixtures check these identities and inequalities; they are
not numerical evaluations of the interacting vacuum, formal
proof-assistant verification or a numerical physical mass prediction.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym73_physical_mass_memory_metric.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym73_physical_mass_memory_metric.py -v
~~~
