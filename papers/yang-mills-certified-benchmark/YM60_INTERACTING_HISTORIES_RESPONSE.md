# YM-60: interacting vacuum histories and a response-to-spectrum bound

Monty Dabas. 4 October 2026. **Python 3.12 only.**

Keep exactly YM59's scalar interaction and extensive-rate clock:
\[
 A_{N,k}(\theta)=H_{N,k}-\theta V_N,\quad
 H_{N,k}=-\tfrac12\sum_iD_{1,i}^2
          -\tfrac k2\sum_{i=1}^{N-1}(D_{2,i}-D_{2,i+1})^2,
 \quad V_N=\sum_i q_i\cdot q_{i+1},\quad k>0.                 \tag{1}
\]
This chapter constructs subsequential, spatially stationary, interacting
infinite-volume histories for every fixed real theta. They are positive,
stationary in continuous heat time, continuous in mean square for smooth
local readouts, and reflection positive. Their reconstructed heat
Hamiltonian is nonnegative. Their uniqueness, time-zero Markov property,
and a positive uniform interacting gap are not proved.

A second result connects the response tower to an observable spectral
constraint. Write c=<q_0 dot q_1> in any constructed limit and, for
nonzero theta, s=sign(theta), x=s c. Then 0<x<=1 and the centered readout
F=q_0^a+s q_1^a, for any real quaternion coordinate a, satisfies
\[
 \boxed{\frac{\langle F(0)F(t)\rangle}{\langle F^2\rangle}
       \ \ge\ \exp[-t\,m_{\rm trial}(x)],\qquad
 m_{\rm trial}(x)=\frac{1+2k-kx}{8(1+x)}.}                  \tag{2}
\]
Thus its spectral threshold is **at most** m_trial. This is an upper
bound on a possible decay gap, not a positive lower mass bound or an
identification of the true mass. At differentiability points of YM59's
energy density, c=-e'(theta). At k=1, all the constructed limits obey
\[
 \begin{split}
 c&=\theta/2+(1123/11844)\theta^3+O(|\theta|^5),\\
 m_{\rm trial}&=3/8-|\theta|/4+\theta^2/8
                    -(5207/47376)|\theta|^3+O(\theta^4).
 \end{split}                                               \tag{3}
\]
These are asymptotics, without a numerically certified coupling interval.
The absolute value selects the plus readout for positive coupling and
the minus readout for negative coupling.

## Lineage and analytic input

Q is the existing compact native unit-quaternion carrier with reference
Phi_Q, and Phi_N is its product reference. Use the completed
representation from YM50--59; no new site algebra or finite-matrix
surrogate replaces it. Put A_i=D_{1,i}, Z_i=D_{2,i}, B_i=Z_i-Z_{i+1}.
All statements about (1) use N>=2. The isolated rank-one N=1 operator
is excluded from the finite-vacuum claim.

The finite-dimensional analytic input is the established theory of
bracket-generating sums of squares on compact connected manifolds:
subelliptic compactness, smooth heat kernels and positivity from
controllability. The bracket hypothesis is verified explicitly below;
no constant uniform in N is imported from that theory. See
[L. Hormander, *Hypoelliptic second order differential equations*, Acta Math. 119 (1967), 147--171](https://doi.org/10.1007/BF02392081).
For the reflection-positive Hilbert-space quotient and its continuous
contraction semigroup, see Remark 2.8(b) and Definition 2.2 of
[Jorgensen, Neeb and Olafsson, SIGMA 12 (2016), 058](https://sigma-journal.com/2016/058/sigma16-058.pdf).
The latter also distinguishes reflection positivity from the additional
Markov condition. These are established methods, not new general laws.

The contribution here is their application with volume-independent
local time control to this actual interacting native model, and the
exact relation of its scalar response to the trial decay rate. Written
infinite-volume proofs and exact finite controls are distinguished.
This is not proof-assistant verification, independent peer review, or
a construction of physical four-dimensional Yang--Mills theory.

## YM60-T1: the actual finite-volume vacuum and its ground transform

Choose any incident edge b at site i, and let epsilon=+1 when i is its
left endpoint and -1 when it is its right endpoint. The native brackets give
\[
 [A_i,B_b]=\epsilon D_{3,i},\qquad
 [A_i,[A_i,B_b]]=-\epsilon Z_i.                              \tag{4}
\]
Every site has an incident edge when N>=2. Thus the vector fields in H
and their brackets span all three tangent directions independently at
each site. Compact subellipticity gives compact resolvent; bounded
smooth multiplication by -theta V preserves the form domain, compact
resolvent and regularity of eigenfunctions. The connected carrier and
bracket-generating control give a positivity-improving heat semigroup.
The bounded, strictly positive Feynman--Kac weight preserves that
property for exp(-t A). Consequently its lowest eigenvalue E_N is
simple, with a smooth strictly positive real normalized vector h_N.
There is a positive gap at each fixed finite N, without a uniform
lower bound from this argument.

Let K_N=A-E_N and d mu_N=h_N^2 d Phi_N. The unitary ground transform is
\[
 P_{N,t}f=h_N^{-1}e^{-tK_N}(h_N f),\qquad
 L_N^h=h_N^{-1}K_N h_N=H-2\Gamma(\log h_N,\cdot),             \tag{5}
\]
where
\[
 \Gamma(f,g)=\tfrac12\sum_i(A_if)(A_ig)
                   +\tfrac k2\sum_i(B_if)(B_ig).
\]
P is a positive, unital, reversible Markov semigroup in L2(mu_N).
The diffusion generator is -L_N^h; this fixes the sign in (5).
This uses the true h_N of (1), not exp of a truncated Poisson tower.

For real smooth f,g the cancellation identity is
\[
 \langle h f,A(hg)\rangle_{\Phi_N}
       -\langle hfg,A h\rangle_{\Phi_N}
          =\int h^2\Gamma(f,g)\,d\Phi_N.                   \tag{6}
\]
It holds for arbitrary real smooth h. Taking the actual ground vector gives
\[
 \langle h_N f,K_N(h_Ng)\rangle=\int\Gamma(f,g)\,d\mu_N.     \tag{7}
\]
The scalar potential cancels because multiplication commutes. An
independent pointwise check is
\([M_f,[A,M_g]]=2M_{\Gamma(f,g)}\).
The exact tests use arbitrary positive polynomial h only to check (6);
they do not pretend that those polynomials are the interacting vacuum.

## YM60-T2: a volume-independent time estimate and interacting histories

The ground-state semigroup supplies positive stationary finite-volume
histories through its transition kernels. For a future cylinder
functional G, reversibility and the Markov property give
\[
 \mathbb E_{\mu_N}[\overline{G\circ\mathcal R}\,G]
       =\|\mathbb E[G\mid q(0)]\|_{L^2(\mu_N)}^2\ge0,       \tag{8}
\]
where R reflects heat time about zero. Thus these histories are
reflection positive before any volume limit.

If f is a real smooth local readout, the spectral theorem and (7) imply
\[
 \mathbb E[(f(q(t))-f(q(0)))^2]
 =2\langle h_N f,(1-e^{-tK_N})h_N f\rangle
 \le 2t\int\Gamma_N(f,f)\,d\mu_N
 \le 2t\|\Gamma_\infty(f,f)\|_\infty.                     \tag{9}
\]
Here Gamma_infinity uses all site fields on the support and all bonds
meeting it. It has only finitely many terms. Actual open endpoints
only delete nonnegative squares. In particular the last bound is
independent of N and does not contain the extensive potential norm.
For a polynomial, its coefficient 1-norm supplies a rational upper
bound on the final supremum. No uniform bound on grad(log h_N) is used.

Here is an explicit compactness construction which does not assume
that a unique interacting vacuum already exists.

1. Extend the N-site stationary process to Z by independent Haar sites
   outside [1,N], each frozen in time. This is just a compact auxiliary
   embedding; it is not claimed to have a finite-volume global gap.
2. Average its N translates in which the observation origin is at each
   site r=1,...,N. Call the resulting history law nu_N. For every bounded
   cylinder G, its spatial-translation defect is at most
   2||G||_infinity/N, by telescoping the two end translates.
3. First restrict times to the rationals. The history space
   Q^(Z x rationals) is compact and metrizable. Probability compactness
   gives a common subsequence for all finite-dimensional distributions.
   Normalization, consistency, stationarity, spatial invariance and
   reflection positivity pass to the limit; the last assertion is a
   closed quadratic inequality for each continuous cylinder G.
4. Estimate (9) passes to that limit. Applying it to coordinates, and
   telescoping bounded products of local readouts, makes rational time
   translations strongly continuous on the dense cylinder subspace of
   the path L2 space. They extend uniquely to real time. Equivalently,
   coordinate variables at real t are the L2 limits along rational
   times tending to t; their values remain in the closed compact Q.
   Consistent real-time laws follow, independently of the approximating
   rational sequences. Positivity, stationarity and reflection
   positivity extend by this same continuity argument.

This proves existence, at every fixed real theta, of spatially invariant
interacting histories in continuous heat time, with local mean-square
continuity. It does not assert sample-path continuity. The positive-time
path algebra, reflected inner product, and its null-space quotient give
a Hilbert space with a strongly continuous self-adjoint contraction
semigroup exp(-t K_infinity), K_infinity>=0, and invariant unit vector 1.
There are no spatial Euclidean-invariance or relativistic claims here.

Convex mixtures and limits preserve (8), but the time-zero Markov
factorization need not follow. For example, mixing equally two stationary
symmetric two-state chains with one-step correlation eigenvalues 1/2
and 1/4 gives C(1)=3/8 and C(2)=5/32, whereas C(1)^2=9/64. Each component
is Markov and reflection positive, and the mixture is reflection
positive, but its observed two-state process is not Markov. Accordingly,
K_infinity lives on the reconstructed history space; identification
with a closed generator on time-zero configurations remains open.

## YM60-T3: exact symmetry identities for the interacting vacuum

Common right multiplication q_i -> q_i r commutes with all native left
fields and preserves each overlap. Common left multiplication by any
unit in {+/-1,+/-i,+/-j,+/-k} conjugates each first and second field to
itself or its negative, preserving the kinetic squares and the scalar
potential. Uniqueness and positivity of h_N force invariance under
these symmetries. In particular each one-site marginal of mu_N is Haar,
even though neighbouring sites are correlated.

Writing c_ij=<q_i dot q_j>, these symmetries give exactly
\[
 \langle q_i^a\rangle=0,\qquad
 \langle q_i^a q_j^b\rangle=\tfrac14\delta_{ab}c_{ij},
 \qquad c_{ii}=1.                                         \tag{10}
\]
One elementary proof averages over common left and right quaternion
units: its action on a bilinear has average delta_ab(q_i dot q_j)/4.
Common right multiplication alone would allow an additional quaternionic
commutant for cross-site bilinears; the left symmetries are used too.
These identities survive extension, spatial averaging and the limit.

Let d_i be the physical degree, one at an endpoint and two otherwise.
Since D_a acts on the four coordinates by a skew matrix with
D_a^T D_a=I/4, (10) yields the exact form matrix
\[
 \begin{split}
 \langle\Gamma(q_i^a,q_i^b)\rangle
     &=\delta_{ab}(1+kd_i)/32,\\
 \langle\Gamma(q_i^a,q_{i+1}^b)\rangle
     &=-\delta_{ab}k c_{i,i+1}/32,\\
 \Gamma(q_i^a,q_j^b)&=0\quad (|i-j|>1).
 \end{split}                                               \tag{11}
\]
The neighbouring entry is negative because the bond field is Z_i-Z_j.
Neither the sign nor the factor 1/32 can be dropped. Thus for
f=sum_i b_i q_i^a the finite-vacuum Rayleigh quotient is
\[
 \frac{\sum_i(1+kd_i)b_i^2-2k\sum_i c_{i,i+1}b_i b_{i+1}}
      {8\sum_{i,j}c_{ij}b_i b_j},                          \tag{12}
\]
whenever the denominator is positive. This is a variational upper bound,
not the spectrum of an assumed invariant coordinate subspace. The
interaction couples that subspace to higher polynomial degrees.

For comparison, a single physical endpoint already gives
Delta_N(theta)<= (1+k)/8. Such boundary readouts do not determine the
bulk spectrum obtained by the spatial average in T2.

## YM60-T4: the energy derivative fixes a surviving bulk correlation

Finite-volume simplicity and bounded analytic dependence on theta give
Hellmann--Feynman:
\[
 -E_N'(\theta)=\sum_{i=1}^{N-1}\langle v_i\rangle_{\mu_N}.
                                                               \tag{13}
\]
In nu_N, <q_0 dot q_1> is this sum divided by N: the one term crossing
into the independent frozen exterior has zero mean. Along any history
subsequence it converges to c. Concavity and the pointwise convergence
E_N/N -> e from YM59 imply, by the two secant inequalities,
\[
              -e'_-(\theta)\le c\le-e'_+(\theta).          \tag{14}
\]
One-sided derivatives exist because e is concave. At differentiability
points, c=-e'. No unproved exchange of differentiation and volume limits
is required. At theta=0 both one-sided derivatives vanish, so c=0.

For theta!=0, evenness, concavity and the strict dimer energy bound give
\[
 1\ge x:=\operatorname{sign}(\theta)c
    \ge-\frac{e(\theta)}{|\theta|}
    \ge\frac{\sqrt{\omega_k^2+\theta^2}-\omega_k}{4|\theta|}>0,
 \qquad\omega_k=(1+3k)/4.                                 \tag{15}
\]
Thus every constructed history limit retains a nonzero scalar bond
correlation. It is not accidentally the theta=0 product reference.
This remains true even if different subsequences select different states.

At k=1, YM59 supplies uniform energy asymptotics at every fixed order.
To obtain the derivative in (3), use its ninth-order construction:
e=P_8(theta)+O(theta^10), with P_8 even. Secants at
h=|theta|^5 sandwich all supergradients between
P_8'(theta)+O(h)+O(|theta|^10/h). Hence both sides of (14) equal
\[
 c=\theta/2+(1123/11844)\theta^3+O(|\theta|^5).              \tag{16}
\]
The unknown sixth- and eighth-order energy coefficients only enter the
stated error. Using just e=-theta^2/4-b theta^4+O(theta^6) would not
justify differentiating its remainder to O(theta^5). For example,
a sufficiently small multiple of theta^6 sin(theta^-2), added to the
negative quadratic on a small interval, preserves concavity but has an
oscillating derivative of order theta^3. This is why the higher-order
uniform tower, rather than a formal derivative of one remainder, matters.

## YM60-T5: from response to a rigorous spectral upper constraint

Take F=q_0^a+s q_1^a with s=sign(theta), or s=1 at theta=0. Equations
(10)--(11), after the spatial average removes the endpoint density, give
\[
 C(0)=\langle F^2\rangle=(1+x)/2,
 \qquad \lim_N\int\Gamma_N(F,F)\,d\nu_N
                       =(1+2k-kx)/16.                    \tag{17}
\]
Here the Gamma integral means the average of the appropriate component
forms; the mixed history law itself need not be a Markov process.
For each component its centered autocorrelation is a positive spectral
mixture of exp(-t lambda), lambda>=0. The same holds after averaging.
The first moments are bounded by (9), so no spectral mass escapes to
infinity. Every limiting autocorrelation has a positive spectral
measure rho with total mass C(0), and first moment at most the limit
in (17). Equality of first moments is not assumed: a first-moment
bound alone does not prove uniform integrability of those moments.

Jensen's inequality, or its finite-volume version followed by the
limit, now gives precisely (2). The infimum of the support of rho,
including a possible zero atom, is at most m_trial. If a proposed
uniform gap entails decay of every centered vacuum readout at rate
Delta, it must satisfy Delta<=m_trial. A zero atom is not ruled out;
in a vacuum with additional stationary sectors, centered need not
mean orthogonal to the entire zero-energy space.

At k=1, inserting (16) into (2) gives the expansion (3). In particular
m_trial decreases at first order by |theta|/4. This is a change of an
upper variational constraint, not proof of a shifted eigenvalue. For
any k>0 its derivative with respect to x is
\[
 m_{\rm trial}'(x)=-(1+3k)/[8(1+x)^2]<0.                    \tag{18}
\]
The full four real coordinate choices share this bound by symmetry.

An exact finite-chain check follows by adding the Rayleigh inequalities
for F_i=q_i^a+s q_{i+1}^a. Set
bar c=(N-1)^-1 sum_i c_(i,i+1). Then
\[
 \Delta_N\le
  \frac{1+2k-k/(N-1)-ks\bar c}{8(1+s\bar c)},             \tag{19}
\]
since sum_i(d_i+d_(i+1))=4N-6. The boundary correction and the separate
endpoint bound must be retained when discussing open-chain gaps.

## YM60-T6: what would turn this into a positive mass-gap result?

For the actual finite vacuum, a width-uniform lower gap is equivalent to
\[
 \int\Gamma_N(f,f)\,d\mu_N
      \ge\Delta\,\operatorname{Var}_{\mu_N}(f)
 \quad\hbox{for every real smooth f and every N,}           \tag{20}
\]
with a single Delta>0 on a specified fixed nonzero coupling interval.
Equation (7) proves this equivalence by the unitary ground transform.
YM57 proves the analogous inequality for Phi_N, not for mu_N. The
finite-order extensive residual in YM59 does not transfer it.

Neither reflection positivity nor an exact first spectral moment
supplies (20). For example, for any m>0 and 0<epsilon<m the positive
spectral measure (delta_epsilon+delta_(2m-epsilon))/2 has first moment
m and threshold epsilon. Letting epsilon tend to zero preserves that
moment. Adding a zero atom can also retain a positive first moment.
These are explicit refusal controls against relabelling (2) as a lower gap.

The next lower-bound route must control the actual interacting measure
or the full centered heat evolution uniformly in width. A convergent
local vacuum construction with an applicable weighted Poincare theorem,
or a direct uniform correlation estimate, would meet this obligation.
A chosen Gibbs weight, a finite polynomial projection, and a finite-N
positive gap do not by themselves do so. Physical gauge constraints,
clock/action selection and a four-dimensional continuum limit remain
additional obligations.

## Reproducible evidence

The exact certificate checks the native brackets, symmetry intertwinings,
left/right bilinear average, scalar cancellation and double commutator,
weighted first-moment identities, physical-boundary incidence, response
coefficients, and rational counterexamples for moment and Markov claims.
The controls do not numerically approximate the true h_N, prove compactness
by sampling, or replace the written limiting arguments.

```sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym60_interacting_histories_response.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym60_interacting_histories_response.py -v
```

Default and --check are read-only. Only --write updates this chapter's
result and canonical hash. The source manifest preserves all 153
upstream inputs. The dedicated workflow also replays YM59 and YM58.
