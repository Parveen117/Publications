# OQD-1 — Observer removal and faithful cut descent

17 October 2026. Written continuation of the actual source and native
cut in [CWR-1](RH_CUT_WAVE_RETURN.md), based on Publications commit
53f6454444c809f4ae4f42921b26a43ac061378b.

The observer-removal question has an exact algebraic answer on the
retained polynomial carrier. A direction can be removed while preserving
every source pairing precisely when it lies in the pairing's radical.
For this actual source that radical is zero, without assuming RH.

A finite observer bank can nevertheless miss a nonzero direction.
If the full source had a negative direction, every finite positive bank
would leave a negative residual invisible to that particular bank.
This is finite observational blindness, with an explicit Schur residual;
it is not a direction invisible to every lawful layer observation.

The native two-channel cut does have a nonzero kernel. Compressing to
its diagonal range preserves an exact copy of the entire original
source pairing. Consequently this cut quotient retains, rather than
resolves, the all-order RH positivity question. Replacing the source
by an exactly matching probability law requires the same condition.

All statements concern the specified algebra, source functional and
family of readings. An observer here is a test polynomial or declared
readout; no claim about consciousness, intrinsic essence or a physical
measurement process is needed.

## Retained source and admissible layer readings

Fix any finite anchor a>=a_0=3249/4. Keep
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad q=\mathcal X'/\mathcal X,\qquad
\mu=(\gamma-i\delta)^2,\quad \delta=\beta-\tfrac12,
\]
\[
u_\mu=\frac a{a+\mu},\qquad
w_\mu=\frac{\nu_\mu}{q(a)(a+\mu)}
      =\frac{\nu_\mu u_\mu}{a q(a)}.
\tag{1}
\]
The sum is over distinct source values, with multiplicities nu_mu.
The retained product and MOR/CWR bounds give
\[
q(a)>0,\qquad 0<|u_\mu|<1,\qquad
\sum_\mu |w_\mu|\le2,\qquad
\sum_\mu w_\mu=1.
\tag{2}
\]
The support S_a={u_mu} is infinite, conjugation symmetric, and has
only 0 as an accumulation point. Its elements and all weights are
nonzero. No source weight is assumed positive.

On the real vector space V=R[y], set
\[
\mathcal L_a[p]=\sum_\mu w_\mu p(u_\mu),\qquad
B_a(p,h)=\mathcal L_a[ph],\qquad Yp=yp,
\]
\[
O_n(p)=\mathcal L_a[y^n p]=B_a(p,y^n),\qquad n\ge0.
\tag{3}
\]
These are real by conjugation symmetry. Products remain bilinear
products, including when their evaluations are complex. They are
not replaced by absolute squares. The source matrix on degree-below-N
polynomials is H_N(a)=(c_{i+j}(a)), where c_j(a)=L_a[y^j].

The observable family is explicit: every finite layer n is allowed.
Equivalently one may use CWR's wave polynomials tau_n=T_n(2y-1).
They have nonzero leading coefficients and span the same polynomial
space, so testing against all tau_n and all y^n has the same kernel.

## OQD-T1. Exact pairing descent and layer-compatible removal

For any symmetric bilinear form B on V and any linear subspace K,
the rule
\[
\overline B([p],[h])=B(p,h)\quad\hbox{on }V/K
\]
is independent of representatives if and only if
\[
\boxed{
K\subseteq\operatorname{rad}B,\qquad
\operatorname{rad}B=\{k:B(k,h)=0\ \hbox{for every }h\in V\}.
}
\tag{4}
\]

**Proof.** Since [k]=[0] for k in K, a well-defined rule requires
B(k,h)=0 for every h. Conversely that condition, symmetry and
bilinearity make both changes of representative vanish. QED.

For the source in (3), retaining only the seed reading O_0 still
forces the same condition if multiplication by y must descend:
\[
YK\subseteq K,\quad O_0(K)=0
\quad\Longrightarrow\quad K\subseteq\operatorname{rad}B_a.
\tag{5}
\]
Indeed Y^n k belongs to K, so O_n(k)=O_0(Y^n k)=0 for every n.
Linearity then gives B_a(k,h)=0 for every polynomial h.

Thus keeping the seed and every lawful layer step already retains
the whole pairing. Dropping the layer-stability requirement changes
this conclusion. The distinction is the algebraic counterpart of
[NT-6](papers/native-thermodynamic-curvature/THEOREM.md): a quotient
must respect its declared operations. Connection descent there
requires preservation of a kernel; unchanged pairing descent here
requires the additional, explicitly proved radical condition (4).

## OQD-T2. The actual source has no nonzero fully invisible polynomial

At every finite a>=a_0,
\[
\boxed{
\operatorname{rad}B_a
=\bigcap_{n\ge0}\ker O_n
=\{0\}.
}
\tag{6}
\]

**Proof.** For a real polynomial p form the source resolvent
\[
\Phi_{a,p}(z)
=\sum_{n\ge0}O_n(p)z^n
=\sum_\mu\frac{w_\mu p(u_\mu)}{1-z u_\mu}
\quad (|z|<1).
\tag{7}
\]
Absolute convergence follows from (2) and boundedness of p on the
closed unit disk, so the geometric series may be interchanged.

The sum on the right extends meromorphically to the whole z-plane.
Its possible poles z_mu=1/u_mu=1+mu/a have no finite accumulation.
On a compact set avoiding those points, all but finitely many
denominators have modulus at least 1/2 because u_mu tends to zero.
The remaining tail is uniformly bounded by a constant times
sum |w_mu|. This proves local normal convergence away from the
possible poles and the corresponding statement after omitting any
one potentially singular term.

At each distinct source point,
\[
\operatorname*{Res}_{z=1/u_\mu}\Phi_{a,p}(z)
=-\frac{w_\mu p(u_\mu)}{u_\mu}
=-\frac{\nu_\mu}{a q(a)}p(u_\mu).
\tag{8}
\]
If p belongs to the radical, every coefficient in (7) is zero.
The meromorphic identity theorem makes the continued function zero,
so every residue in (8) vanishes. The prefactor is nonzero, hence
p vanishes on the infinite set S_a. A nonzero polynomial has only
finitely many roots. Therefore p=0. The equality with the
intersection of kernels follows from the monomial basis. QED.

This proof uses no positive completion, probability measure or
all-order sign assumption. A nondegenerate bilinear form need not
be positive; (6) by itself does not settle RH.

Several precise consequences follow.

- Every nonzero polynomial is detected by some finite layer
  O_n, and by some finite wave-polynomial test.
- No nonzero subspace can be removed while preserving every B_a
  reading. By (5), there is also no nontrivial quotient retaining
  both the seed reading and the full layer-multiplication rule.
- In particular, an exact finite-dimensional quotient with both
  of those properties cannot represent this infinite carrier.
  Finite models can still match a specified finite set of readings.
- For each fixed degree bound d, some finite prefix O_0,...,O_M
  separates all polynomials of degree at most d. To see this, the
  decreasing kernels restricted to that finite-dimensional space
  eventually stabilize; their intersection is zero by (6).
  This supplies existence of M, not an explicit bound or a new
  finite certification.
- A zero self-reading B_a(p,p)=0 is not the definition of invisibility.
  Invisibility requires every cross-reading B_a(p,h) to vanish.

## OQD-T3. A finite positive bank can leave an invisible negative residual

Let U=span(u_1,...,u_N) be any finite subspace whose actual Gram
matrix G_U=(B_a(u_i,u_j)) is positive definite. For p in V define
\[
b(p)_j=B_a(u_j,p),\qquad
\Pi_U p=\sum_{i=1}^N u_i\bigl(G_U^{-1}b(p)\bigr)_i,\qquad
k=p-\Pi_U p.
\tag{9}
\]
Then
\[
B_a(k,U)=0,\qquad
V=U\oplus U^{\perp_{B_a}},\qquad
\boxed{B_a(p,p)=b(p)^T G_U^{-1}b(p)+B_a(k,k).}
\tag{10}
\]
If B_a(p,p)<0, the resulting residual obeys
\[
\boxed{
O_U(k)=0,\qquad
B_a(k,k)\le B_a(p,p)<0,
}
\tag{11}
\]
where O_U(k) is the vector of its pairings against the chosen u_j.

**Proof.** Multiplying the coefficient vector in (9) by G_U
gives b(p), which proves orthogonality. The intersection of U with
its orthogonal complement is zero because G_U is invertible;
(9) gives the decomposition for every p. Expand B_a(p,p), remove
the orthogonal cross terms and use positivity of G_U to obtain
(10)-(11). QED.

For the monomial bank U=span(1,y,...,y^{N-1}) with H_N(a)>0,
take p=y^N and h=(c_N,...,c_{2N-1})^T. Its residual is precisely
the already retained Schur direction:
\[
k_N=y^N-\sum_{j=0}^{N-1}y^j(H_N^{-1}h)_j,\qquad
\boxed{B_a(k_N,k_N)=c_{2N}-h^T H_N^{-1}h.}
\tag{12}
\]
This scalar is the next sign gate, not a sign fixed by projection.

Equations (10)-(11) formalize one version of the moving-shadow
intuition. Conditional on a negative full-source direction, every
accepted finite positive bank leaves a negative residual invisible
to that bank. The residual depends on the bank. This does not say
that every finite monomial bank is positive, or that one fixed
nonzero direction is invisible to all observers. The latter is
excluded by (6). No actual negative direction is exhibited here.

Quotienting by U's invisible subspace preserves its selected
readout vector but does not preserve the full B_a form. The omitted
term B_a(k,k) in (10) records exactly what was discarded.

## OQD-T4. The native cut contains the entire original pairing

Use the exact polynomial version of the retained LOT/CWR cut:
\[
\mathcal E=V\oplus V,\qquad
\mathcal G_a((p,q),(h,k))
=\mathcal L_a[yph+(1-y)qk],
\]
\[
g(y)=\operatorname{diag}(y,1-y),\qquad
P(y)=\begin{pmatrix}y&1-y\\y&1-y\end{pmatrix},\qquad
Jf=(f,f).
\tag{13}
\]
Polynomial multiplication gives P^2=P and P^T g=gP, independently
of source positivity. Its range and kernel are
\[
\operatorname{ran}P=J V,\qquad
K_P=\ker P=\{((1-y)h,-yh):h\in V\}.
\tag{14}
\]
For x=(p,q), put r=yp+(1-y)q and h=p-q. The exact decomposition is
\[
x=Jr+((1-y)h,-yh),
\]
\[
\boxed{
\mathcal G_a(x,x)
=\mathcal L_a[r^2]+\mathcal L_a[y(1-y)h^2].
}
\tag{15}
\]
Its two summands are orthogonal for the bilinear form. On the
retained diagonal,
\[
\boxed{
PJ=J,\qquad
\mathcal G_a(Jf,Jh)=B_a(f,h).
}
\tag{16}
\]

**Proof.** Substitution proves the decomposition and (14), since
r=0 exactly characterizes the kernel and h=p-q recovers its
coefficient. The cross term of (15) is
L_a[y(1-y)rh-(1-y)yrh]=0. Its complementary square is
L_a[(y(1-y)^2+(1-y)y^2)h^2]=L_a[y(1-y)h^2].
The two weights add to one on the diagonal, proving (16). QED.

The full two-channel form also has zero radical. If (p,q) pairs
to zero with every (h,0), then B_a(yp,h)=0 for every h. By (6),
yp=0 and p=0. Testing (0,k) similarly gives (1-y)q=0 and q=0.
Hence G_a itself cannot descend unchanged through the nonzero K_P.

There is a legitimate compressed quotient form:
\[
\boxed{
\widehat{\mathcal G}_a([x],[z])
=\mathcal G_a(Px,Pz)
\quad\hbox{on }\mathcal E/K_P.
}
\tag{17}
\]
It is well defined because P kills K_P. The map [x] to r in (15)
identifies this quotient and its form exactly with (V,B_a).
Equivalently the form (x,z) to G_a(Px,Pz) on E has radical
exactly K_P, by (6) and (16).

Thus the cut can reduce the two channels to one while retaining
the complete original polynomial source. It does not reduce
that source to finitely many readings. If a polynomial f had
B_a(f,f)<0, its negative direction Jf would already lie on the
cut and satisfy PJf=Jf. Removing the complementary channel
would not remove it.

The full polynomial construction above requires no Hilbert-space
completion. LOT-T6's orthogonal Hilbert completion remains conditional
on its explicitly stated all-order positivity premise.

## OQD-T5. Exact probability replacement retains the RH condition

For each fixed finite a>=a_0, the following statements are equivalent:
\[
\boxed{
\begin{array}{ll}
\text{(i)}&\mathrm{RH},\\
\text{(ii)}&B_a(p,p)\ge0\quad\hbox{for every real polynomial }p,\\
\text{(iii)}&\widehat{\mathcal G}_a\text{ is positive semidefinite},\\
\text{(iv)}&\text{There is a probability measure }\nu\text{ on }[0,1]\\
&\text{with }\mathcal L_a[p]=\int p\,d\nu
       \quad\hbox{for every real polynomial }p.
\end{array}
}
\tag{18}
\]

**Proof.** Under RH the actual u_mu lie in (0,1) and the actual
w_mu are positive and sum to one. Therefore the explicit discrete
measure sum w_mu delta_{u_mu} proves (iv), and (iv) proves (ii)
by integrating p^2. Equation (16)-(17) makes (ii) and (iii)
equivalent.

For completeness, the retained [finite detector D1-D3](RH_FINITE_OBSERVER_DETECTOR.md)
applies at every anchor used here, not only its original a_0.
Its argument needs exactly the discrete conjugate support, the
nonzero weights in (1), the support bound and absolute budget
in (2), all retained at this a. If RH failed, take a nonreal
z in S_a and let Q be the real polynomial annihilating all
other support values of modulus at least s=|z|. The remaining
support has maximum modulus s_0<s. The two observers
y^n Q(y), y^{n+1} Q(y) have a Gram matrix
\[
P(A_n,z)+R_n,\qquad
A_n=w_z z^{2n}Q(z)^2,\qquad
|(R_n)_{ij}|=O(s_0^{2n}),
\]
\[
P(A,z)=
\begin{pmatrix}2\Re A&2\Re(Az)\\2\Re(Az)&2\Re(Az^2)\end{pmatrix},
\qquad
\det P(A,z)=-4|A|^2(\Im z)^2.
\tag{19}
\]
Here |A_n|=|w_z Q(z)^2|s^{2n}, with nonzero constant factor.
Divide the determinant by |A_n|^2. Every perturbation term
tends to zero, so the result tends to -4(Im z)^2<0.
Some finite Gram compression is therefore indefinite, contradicting
(ii). This proves (ii) implies RH, and completes (18). QED.

The same equivalence holds for positivity of the full G_a:
under RH both weighted channels are positive; conversely its
diagonal restriction is B_a.

A finite positive probability model for finitely many readings
does not assert (iv). An exact infinite replacement must match
every polynomial reading of the actual source. Likewise, any
positive bilinear representation that preserves every original
pairing B_a(p,h) would already imply (ii). Choosing positivity
in the representing space does not establish that matching map.

Neither (6) nor the projection identity P^2=P implies any of the
equivalent positivity assertions in (18). No actual negative
witness or proof of all-order positivity is supplied.

## OQD-T6. Common normalization does not identify distinct prime anchors

Use the existing prime-source sampling convention from
[POC-T3](papers/native-prime-observation-closure/THEOREM.md):
a_p=(p-1/2)^2 for a prime p. For this note's anchor range take
p>=29. CWR's source prelude already proves q'(a)<0 for a>=a_0,
so distinct such primes have distinct unnormalized q(a_p).

Even after dividing out q(a), no two distinct finite anchors
have the same complete fixed-coordinate moment sequence:
\[
\boxed{
a\ne b,\quad a,b\ge a_0
\quad\Longrightarrow\quad
(c_n(a))_{n\ge0}\ne(c_n(b))_{n\ge0}.
}
\tag{20}
\]
Equivalently, their normalized source germs differ:
\[
F_a(z)=\frac{q(a(1+z))}{q(a)}
      =\sum_{n\ge0}(-1)^n c_n(a)z^n.
\tag{21}
\]

**Proof.** Suppose F_a=F_b as germs and, after exchanging the
anchors if necessary, put rho=b/a>1. Setting t=a(1+z) yields
\[
q(\rho t)=\frac{q(b)}{q(a)}q(t)
\tag{22}
\]
near t=a and hence as a meromorphic identity. The actual logarithmic
derivative q has a nonempty discrete pole set {-mu}; there is no
pole at zero, and its smallest pole modulus is some d>0. Every
such pole has its nonzero multiplicity as residue. The poles of
q(rho t) therefore have smallest modulus d/rho, whereas multiplying
q(t) by the nonzero scalar q(b)/q(a) leaves smallest modulus d.
Equation (22) would force d/rho=d, impossible for rho>1. QED.

Thus some finite derivative layer distinguishes any specified
pair of these normalized anchors. No uniform separating depth is
claimed. Their common c_0=1 is one normalization identity, not
equality of their whole source responses.

This is a comparison at the same polynomial labels. MOR's exact
co-transport changes the observer as well as the anchor and preserves
the appropriate carried pairing; that covariance is compatible with
(20). It does not identify different prime integers. Describing all
primes as one structural class or assigning them one selected reading
would be a different declared equivalence, not a consequence that
their complete actual-source readings agree.

## Decision and remaining obligation

The proposed observer removal can now be assessed without assuming
a positive source norm.

| Proposed removal or replacement | Exact result |
|---|---|
| Forget directions invisible to one finite bank | Possible for that readout; generally discards part of the original form |
| Remove directions invisible to every polynomial layer observer | Only the zero direction exists, by OQD-T2 |
| Retain the seed and all layer steps on a quotient | No nontrivial quotient, by OQD-T1-T2 |
| Compress the native two-channel cut to its diagonal | Exact reduction to the original full pairing B_a |
| Replace the entire actual source by a matching probability law | Equivalent to RH, by OQD-T5 |
| Identify distinct normalized prime-anchor responses | Excluded for the complete fixed-coordinate germs, by OQD-T6 |

The completed result is a source-specific obstruction to removing
the RH sign problem by this quotient route. It does not prove that
an observer-independent physical object or essence exists. Nor is
it a formal consistency proof or an independence result: no theory,
proof predicate or metamathematical independence argument is
constructed here.

The open mathematical target is unchanged: establish positivity
of every actual-source B_a0 polynomial square, equivalently CWR's
all-step bound abs(r_k(a_0))<=1. The present result shows why an
exact removal of unseen directions cannot supply that missing
inequality. A finite source check or transport consistency identity
must not be substituted for it.

## Evidence and execution scope

The proofs use the retained actual source, the identity theorem,
simple-pole residues, polynomial algebra, orthogonal decomposition
for a declared finite nondegenerate form, and the existing finite
off-line-pair detector. These are standard tools; no general novelty
claim is made for quotient descent or these algebraic operations.

Source lineage:
[NT-6](papers/native-thermodynamic-curvature/THEOREM.md),
[POC](papers/native-prime-observation-closure/THEOREM.md),
[LOT](papers/native-layer-observer-transport/THEOREM.md),
[finite detector](RH_FINITE_OBSERVER_DETECTOR.md),
[MOR](RH_MOVING_OBSERVER_RETURN.md), and
[CWR](RH_CUT_WAVE_RETURN.md).

Only this note and the README entry change. No local scientific
computation, independent agent, source evaluation or scientific
replay is performed; no workflow is changed. These are written
derivations, not new machine-checked certificates. The commit
contains [skip ci], and the change remains a draft stacked on
PR #100 without a merge.
