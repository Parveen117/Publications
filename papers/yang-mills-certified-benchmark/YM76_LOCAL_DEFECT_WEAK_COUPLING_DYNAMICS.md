# YM-76: local gradient defects and actual weak-coupling dynamics

Monty Dabas. 5 October 2026. **Python 3.12 only.**

**Status:** five written actual-model results, with exact finite replay controls.
This stage proves two-sided vacuum scales and a nonzero local spectral band
uniformly over finite periodic lattices. The band is microscopic: in physical
units its energies scale as \(1/a\). A positive finite continuum mass gap and
a nontrivial local four-dimensional quantum field theory are **not proved**.

**Continuation:** YM74's all-coupling vacuum bound and YM75's source-resolved
spectral transfer. This stage works directly with the actual SU(2) Hamiltonian;
it does not assume that a variational trial state is the vacuum.

## 1. Model, normalization and the new mechanism

Use the unchanged periodic cubic spatial lattice of side \(L\ge4\).
There are \(p=3L^3\) positive oriented links and \(p\) elementary plaquettes.
On the product of SU(2) link groups, with normalized product Haar measure,
write
\[
 A_\theta=L_E-\theta V,\qquad V=\sum_b W_b,\qquad
 X=p-V,\qquad H=A_\theta-E_0,\qquad \theta>0.
\]
Here \(W_b\) is the real half-trace of the oriented plaquette holonomy;
the fundamental one-link electric eigenvalue is \(3/4\).
The positive normalized finite-volume ground state is \(h\), and
\(\mu=h^2\,d\Phi\). Ellipticity on the connected compact configuration
manifold gives a unique smooth positive ground state. Its invariance under
gauge transformations and lattice symmetries places it in the physical
sector and makes the following quantities independent of link or plaquette:
\[
 \tau=\langle h,L_eh\rangle,\qquad
 \delta=\mu(1-W_b),\qquad
 \mathscr D=\frac{E_0+\theta p}{p}=\tau+\theta\delta.
\]
The symbol \(\tau\) denotes electric energy per link, not physical time.
The carré du champ is \(\Gamma(f,g)=\nabla f\cdot\nabla g\);
\(\Gamma(f)=\Gamma(f,f)\). Native normalizations give
\[
 L_EW_b=3W_b,\quad \Gamma(W_b)=1-W_b^2,\quad
 L_eW_b=\tfrac34W_b,\quad
 \Gamma_e(W_b)=\tfrac14(1-W_b^2)\quad(e\in b).
 \tag{1}
\]
Every link belongs to four plaquettes.

YM74 proved, without identifying \(h\),
\[
 \mathscr D\le D(\theta):=
 \min\left\{\theta,\frac3{\sqrt2}\sqrt\theta+\frac3{16}\right\}.
 \tag{2}
\]
The exact product-state bound retained by the certificate is also available:
for \(n=\lceil\sqrt{32\theta}\rceil\),
\[
 m_n=\frac{2n}{2n+3},\quad e_n=\frac{3n^2}{4(4n+1)},\qquad
 \mathscr D\le \min\{\theta,e_n+\theta(1-m_n^4)\}.
 \tag{3}
\]

The new step is to expose a positive local remainder in the Hamiltonian,
then use the **plaquette angle**, whose Dirichlet energy is exactly one,
to recover actual vacuum variance and spectral moments. The last step
turns energy estimates into a nonvanishing dynamical state.

## 2. Theorem 76.1 — exact local defect factorization

Set
\[
 t=t_\theta=\frac{\sqrt{9+32\theta}-3}{16},
 \qquad 3t+8t^2=\theta,\qquad Q_t=\nabla-t\nabla V.
\]
Then, as quadratic forms on the full configuration space and hence on
its gauge-invariant subspace,
\[
 \boxed{L_E+\theta X=Q_t^*Q_t+3tp+t^2\mathcal R,}
 \tag{4}
\]
where the multiplication operator \(\mathcal R\) is explicitly local:
\[
 \begin{split}
 \mathcal R&=8X-\Gamma(V)\\
 &=4\sum_b(1-W_b)^2+
 \sum_e\sum_{\substack{b<c\\b,c\ni e}}
       |\nabla_eW_b-\nabla_eW_c|^2\ \ge0.
 \end{split}
 \tag{5}
\]
The inner sum has six comparisons per link: two coplanar and four
orthogonal plaquette pairs. The derivatives are expressed in the same
orthonormal tangent frame at their common link, so their difference is
well-defined and its squared norm is gauge invariant.

**Proof.** At each link apply
\[
 4\sum_{i=1}^4|u_i|^2-\left|\sum_{i=1}^4u_i\right|^2
 =\sum_{i<j}|u_i-u_j|^2.
\]
Sum over links and use (1). Subtracting \(4\sum_b(1-W_b^2)\)
from \(8\sum_b(1-W_b)\) leaves \(4\sum_b(1-W_b)^2\), proving (5).
Integration by parts, using \(L_EV=3V\), gives
\[
 Q_t^*Q_t=L_E-3tV+t^2\Gamma(V).
\]
Insert \(3t+8t^2=\theta\) to obtain (4). These identities hold first
on smooth functions and then on the closed form domain. \(\square\)

Consequently the actual vacuum has the matching-order bounds
\[
 \boxed{\frac3{16}\bigl(\sqrt{9+32\theta}-3\bigr)
 \le\mathscr D\le D(\theta).}
 \tag{6}
\]
Moreover,
\[
 \mathscr D-3t=\frac{\|Q_th\|^2}{p}
                +t^2\frac{\mu\mathcal R}{p},\qquad
 \frac{\mu\mathcal R}{p}\le\frac{D(\theta)-3t}{t^2}.
 \tag{7}
\]
The last upper bound is \((6\sqrt2+o(1))/\sqrt\theta\).
Thus the actual vacuum suppresses both plaquette amplitude defects and
the displayed neighboring gradient mismatches in spatial average.

Although \(Q_t e^{tV}=0\), the remainder \(\mathcal R\) does not vanish
identically. Formula (4) does not identify \(e^{tV}\) as the vacuum.
After subtracting the actual vacuum energy, (4) subtracts the unknown
positive excess in (7); it does not give a positive excitation gap.
The value \(t=0\) is allowed in (4), but division by \(t^2\) in (7)
requires \(\theta>0\).

## 3. Theorem 76.2 — two-sided electric and magnetic scales

Let \(m=\mu W_b=1-\delta\) and \(q=\mu(1-W_b^2)\).
For any link \(e\) of that plaquette, integration by parts and
Cauchy–Schwarz give the actual-vacuum uncertainty inequality
\[
 \tfrac34m=2\langle\nabla_eh,h\nabla_eW_b\rangle,\qquad
 \boxed{9m^2\le16\tau q\le32\tau\delta.}
 \tag{8}
\]
One may sharpen \(q\le2\delta\) to \(q\le\delta(2-\delta)\)
by \(\mu W_b^2\ge(\mu W_b)^2\).

For every \(\theta\ge25\), (2) gives
\(D(\theta)\le(9/4)\sqrt\theta\) and
\(\delta\le9/(4\sqrt\theta)<1/2\). Inserting \(m\ge1/2\)
in (8), and then bounding each factor by (2), proves
\[
 \boxed{\frac{\sqrt\theta}{32}\le\tau\le\frac94\sqrt\theta,\qquad
 \frac1{32\sqrt\theta}\le\delta\le\frac9{4\sqrt\theta}.}
 \tag{9}
\]
In particular,
\[
 \frac1{16\sqrt\theta}\le q\le\frac9{2\sqrt\theta}.
 \tag{10}
\]
All constants are independent of \(L\). This adds lower bounds for
the actual electric energy and magnetic defect; an upper concentration
bound alone did not supply them.

## 4. Theorem 76.3 — radial uncertainty and an actual local excitation

For one plaquette define its angle \(r=\arccos W_b\in[0,\pi]\).
Almost everywhere,
\[
 \Gamma(r)=1,\qquad L_Er=-2\cot r.
 \tag{11}
\]
Define the first-order expression
\[
 J=\nabla r\cdot\nabla+\cot r.
 \tag{12}
\]
The integrations below are legitimate despite the two radial poles.
To see this, change one plaquette link to its holonomy, leaving all
other links fixed. This is a smooth Haar-preserving change of variables.
The radial Haar density is \((2/\pi)\sin^2r\,dr\).
The functions \(h\) and its derivatives are bounded at every fixed
finite volume and coupling. Thus \(\cot r\,h\in L^2\);
flux terms at \(r=\epsilon,\pi-\epsilon\) vanish as \(O(\epsilon^2)\)
for the skew identity and as \(O(\epsilon)\) for the norm identity.
These statements need no volume-uniform bound on pointwise derivatives.

In distributions, (11) has no pole-supported delta: the radial flux
through a shrinking pole surface is \(O(\epsilon^2)\).
Its right-hand side is \(L^2\). Elliptic domain characterization on
the compact configuration manifold gives \(r\in D(L_E)=H^2\).
Multiplication by smooth \(h\) preserves this domain, and the bounded
potential preserves \(D(H)=D(L_E)\). In particular \(rh\in D(H)\),
which will justify the second spectral moment in the next theorem.
No self-adjoint radial momentum operator is asserted.

Cutoff integration gives \(J^*=-J\) for the pairings used here,
\([J,r]=1\), and
\[
 \|Jh\|^2=\|\nabla r\cdot\nabla h\|^2+1\le4\tau+1.
 \tag{13}
\]
For completeness, the cross term in the squared norm is
\[
 \int\cot r\,\nabla r\cdot\nabla(h^2)
 =-\int h^2(-\csc^2r+2\cot^2r);
\]
adding the \(\cot^2r\) term leaves exactly \(\|h\|^2=1\).
The upper bound in (13) uses \(\Gamma(r)=1\) and only the four
link derivatives belonging to the plaquette.

Writing \(w=(r-\mu r)h\), the commutator gives
\[
 2\langle w,Jh\rangle=-1,\qquad
 \boxed{\operatorname{Var}_\mu(r)\ge\frac1{16\tau+4}.}
 \tag{14}
\]
On the other hand \(r^2\le(\pi^2/2)(1-\cos r)\), so
\[
 \operatorname{Var}_\mu(r)\le\frac{\pi^2\delta}{2}.
 \tag{15}
\]
The ground-state form identity gives the exact numerator
\[
 \langle w,Hw\rangle=\mu\Gamma(r)=1.
 \tag{16}
\]
Thus \(w\) is an explicitly specified, nonzero, gauge-invariant state
orthogonal to the actual vacuum. For all finite \(\theta>0\),
\[
 \Delta_\Lambda(\theta)\le16\tau+4\le16D(\theta)+4.
 \tag{17}
\]
For \(\theta\ge25\), simpler uniform constants are
\[
 \boxed{\frac1{40\sqrt\theta}\le\operatorname{Var}_\mu(r)
       \le\frac{9\pi^2}{8\sqrt\theta}.}
 \tag{18}
\]

There is also a trace-variance consequence. For \(r,s\in[0,\pi]\),
\(|\cos r-\cos s|\ge2|r-s|^2/\pi^2\). Apply this to independent
copies, then Jensen's inequality, to obtain
\[
 \operatorname{Var}_\mu(W_b)
 \ge\frac8{\pi^4}\operatorname{Var}_\mu(r)^2
 \ge\frac1{2\pi^4(1+4\tau)^2}.
 \tag{19}
\]
This lower bound does not assume an unproved vacuum susceptibility.

## 5. Theorem 76.4 — positive spectral weight and nonconstant local dynamics

The operator product rule, (11), and \(Hh=0\) give
\[
 \boxed{Hw=-2Jh,\qquad \|Hw\|^2\le16\tau+4.}
 \tag{20}
\]
Normalize \(u=w/\|w\|\), and let \(\nu_{\theta,L}\) be its spectral
probability measure for the microscopic-clock operator \(H/\sqrt\theta\).
For \(\theta\ge25\), its first two moments obey
\[
 m_1=\frac1{\sqrt\theta\,\operatorname{Var}_\mu(r)},\qquad
 \boxed{\frac4{45}\le m_1\le40,\qquad m_2\le40m_1\le1600.}
 \tag{21}
\]
Indeed (18) and \(\pi^2<10\) give the lower bound, while
\((16\tau+4)/\sqrt\theta\le36+4/5<40\) gives both upper bounds.

Set \(\alpha=4/45\), \(\beta=40\). The first moment below
\(\alpha/4\) is at most \(\alpha/4\). Above \(2\beta\), it is
at most \(m_2/(2\beta)\le m_1/2\). The first moment in the
closed middle band is therefore at least
\(m_1/2-\alpha/4\ge\alpha/4\). Dividing by the band's upper endpoint
proves the explicit, volume-independent statement
\[
 \boxed{\nu_{\theta,L}([1/45,80])\ge\frac1{3600}.}
 \tag{22}
\]
This is spectral weight of an actual interacting physical state,
not a spectral claim about an auxiliary finite matrix.

For microscopic time \(s\ge0\), its heat correlation satisfies
\[
 \boxed{
 e^{-40s}\le
 C_{\theta,L}(s):=\langle u,e^{-sH/\sqrt\theta}u\rangle
 \le1-\frac1{3600}(1-e^{-s/45}).
 }
 \tag{23}
\]
The lower bound is Jensen's inequality applied to (21);
the upper bound uses the mass in (22).

Along any sequence of finite lattices and couplings \(\theta\ge25\),
the measures are tight, since their second moments are bounded.
A weakly convergent subsequence is a probability measure on \([0,\infty)\).
Its heat transform is continuous at zero and is nonconstant by (23).
First moments also pass to the limit by uniform integrability from
the second-moment bound. The compact closed band in (22) retains at
least the displayed weight by the closed-set Portmanteau inequality.
This is a **scalar local spectral subsequence**, not a joint family of
fields, a common continuum Hilbert space, or a four-dimensional QFT.
The limiting measure may have an atom at zero.

An equivalent unnormalized choice makes the surviving dynamics explicit:
\[
 v_\theta=\theta^{1/4}w,\qquad
 \frac1{40}\le\|v_\theta\|^2\le\frac{9\pi^2}{8},\qquad
 \langle v_\theta,(H/\sqrt\theta)v_\theta\rangle=1,\qquad
 \|(H/\sqrt\theta)v_\theta\|^2\le40.
 \tag{24}
\]

## 6. Physical scaling: the result and its obstruction

Keep the physical clock used in YM73–75:
\[
 \theta(a)=\frac{c_B}{c_Eg(a)^4},\qquad
 \kappa(a)=\frac{c_Eg(a)^2}{a},\qquad
 \kappa(a)\sqrt{\theta(a)}=\frac{\sqrt{c_Ec_B}}a.
 \tag{25}
\]
The band (22) becomes
\[
 \left[\frac{\sqrt{c_Ec_B}}{45a},
       \frac{80\sqrt{c_Ec_B}}a\right].
 \tag{26}
\]
Thus a fixed positive fraction of this local state's spectral weight
survives at ultraviolet energies. Its normalized physical mean energy
has order \(1/a\); the physical form energy of \(v_\theta\) is
exactly \(\sqrt{c_Ec_B}/a\).

Neither statement says that the *lowest* excitation diverges or that
infrared states are absent. Bounds on one state's first two moments
cannot exclude other states or an additional spectral component near zero.
For example, a probability measure with most of its mass at zero and
a positive mass at a nonzero point can satisfy (21). This is a logical
control, not a proposed Yang–Mills vacuum.

The band proves more than a divergent mean: it forces a fixed positive
amount of ultraviolet spectral weight to escape on the physical clock.
For the same normalized states and any physical time \(t_{\rm phys}>0\),
\[
 C_{{\rm phys},a}(t_{\rm phys})
 \le1-\frac{1-\exp[-t_{\rm phys}\sqrt{c_Ec_B}/(45a)]}{3600}.
 \tag{26a}
\]
As \(a\to0\), the limiting upper bound is \(3599/3600\) for every
positive time, whereas \(C_{{\rm phys},a}(0)=1\).
Therefore this unfiltered, normalized microscopic probe cannot have
a norm-preserving physical Hilbert-state limit with a strongly continuous
heat semigroup whose correlations are the limits of these lattice heat
correlations. State convergence without convergence of the dynamics is
not excluded by this argument. A surviving low-energy projection may lose norm; spatial
smearing or another choice of probe may change the conclusion.
This obstruction concerns this state sequence, not continuum Yang–Mills
or local operator-valued distributions, which need not create Hilbert
vectors without smearing.

## 7. Theorem 76.5 — a concrete collective covariance target

For all plaquettes let
\[
 G_{bc}=\operatorname{Cov}_\mu(r_b,r_c),\qquad
 K_{bc}=\mu\Gamma(r_b,r_c).
 \tag{27}
\]
Both matrices are positive semidefinite, \(K_{bb}=1\), and
\[
 \boxed{0\le K\le4I.}
 \tag{28}
\]
To prove the upper bound, use Cauchy–Schwarz on the four plaquette
gradient terms incident at each link:
\[
 \Gamma\!\left(\sum_b c_b r_b\right)
 \le4\sum_b c_b^2\Gamma(r_b)=4\|c\|^2.
\]
For real \(c\) with \(c^TGc>0\), the actual centered collective state
\(F_ch=\sum_b c_b(r_b-\mu r_b)h\) has
\[
 \|F_ch\|^2=c^TGc,\qquad
 \langle F_ch,HF_ch\rangle=c^TKc,\qquad
 \boxed{\Delta_{\rm phys}\le\kappa\frac{c^TKc}{c^TGc}
                       \le\frac{4\kappa\|c\|^2}{c^TGc}.}
 \tag{29}
\]
This is an exact finite-volume variational statement.

It specifies a useful next sufficient target: for coefficients supported
in a declared fixed physical region, prove for a selected family \(c_a\)
an actual covariance bound \(c_a^TGc_a\ge\chi(a)\|c_a\|^2\).
For a target energy \(E_*>0\), the sufficient condition
\(\chi(a)\ge4\kappa(a)/E_*\) would give normalized states with physical
mean energy at most \(E_*\). No all-directions bound \(G\ge\chi I\)
is asserted. Relative to the individual
variance scale \(1/\sqrt\theta\), this sufficient route asks for covariance
amplification of order \(1/a\), with the same clock constants.
It is **not necessary**: small \(c^TKc\) through gradient cancellation
could give a different route. No such covariance amplification is
proved here, and bounded mean energy alone would still require
state/field convergence and a separate positive lower spectral edge.

For \(\theta\ge25\), there is a quantitative size requirement for this particular route.
If the chosen coefficients use \(N\) plaquettes, (18) implies
\(\lambda_{\max}(G_{\rm support})\le\operatorname{tr}G_{\rm support}
\le9\pi^2N/(8\sqrt\theta)\).
Consequently the sufficient covariance target can hold only if
\[
 N\ge\frac{32\kappa\sqrt\theta}{9\pi^2E_*}
   =\frac{32\sqrt{c_Ec_B}}{9\pi^2E_*a}.
 \tag{30}
\]
This is a constraint on certifying energy by that covariance target;
it is not a necessary plaquette count for every possible low-energy state.

For response-language compatibility, \(G\) is exactly the Hessian at zero of
\(\log\mu\exp(\sum_b s_b r_b)\). This finite-volume identity gives a
precise susceptibility object; it does not identify an unverified
thermodynamic or electromagnetic symmetry with gauge covariance.

## 8. Evidence, reproducibility and remaining work

The executable evidence has three distinct roles:

1. Exact native rational-quaternion configurations replay the local
   gradient identity on periodic graphs and check the incidence constants.
   Finite configurations test the implementation, while the proof of
   the universal identity is Section 2.
2. Exact rational or directed bounds replay the weak-coupling budgets,
   radial integral controls and spectral-band constants. Radial trial
   densities are labeled controls and never substituted for the actual vacuum.
3. Scope controls retain the possibility of a zero-energy spectral atom
   and distinguish microscopic energy from physical energy.

The new source manifest pins 267 unchanged upstream files and four local
inputs. The result and SHA-256 pin are generated under **Python 3.12 only**.
The targeted suite contains 32 new tests and 68 unchanged YM75/YM74 tests.
This is conventional written mathematics with reproducible finite checks,
not a proof-assistant formalization or a numerical proof of the continuum gap.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym76_local_defect_weak_coupling_dynamics.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym76_local_defect_weak_coupling_dynamics.py -v
~~~

The next infrared work must control actual collective covariance or a
different low-energy construction, propagate a positive physical spectral
lower bound through the source-resolved hidden memory, and construct
nontrivial local continuum fields. The factorization and microscopic
dynamic state above supply new input, but do not discharge those obligations.

## 9. Sources and attribution

- [YM74: all-coupling electric enclosures](YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md):
  native SU(2) conventions and the upper bounds (2)–(3).
- [YM75: RH source audit and spectral transfer](YM75_RH_SOURCE_RESOLVED_SPECTRAL_TRANSFER.md):
  the pinned RH/RKF source records, matrix-valued hidden-memory bounds,
  and stable band reconstruction. Those tools remain conditional on their
  actual Yang–Mills hypotheses; no RH gap is imported.
- Kogut and Susskind, *Hamiltonian formulation of Wilson's lattice gauge
  theories*, Physical Review D 11 (1975), 395,
  [DOI](https://doi.org/10.1103/PhysRevD.11.395): Hamiltonian lattice context.
- Bronzan, *Analytic approach to weak-coupling SU(2) Hamiltonian lattice
  gauge theory*, Physical Review D 37 (1988), 1621,
  [publisher abstract](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.37.1621):
  historical variational weak-coupling context. Only the abstract was
  consulted for this reference; the estimates above are derived here.

“New” means a development in this research programme. No claim of external
priority for square completion, radial uncertainty, or moment inequalities
is made.
