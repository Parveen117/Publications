# PCM-1 — Prime cuts and arithmetic memory response

17 October 2026. Written continuation at Publications commit
3d3151823a9eddd330abd98cec5048762a9f03e6, following
[ASM-1](RH_ARITHMETIC_SEAM_REDUCTION.md).

This note connects an exponential prime-power probability packet
to the actual arithmetic source and then to ASM's Schur memory.
It constructs a normalized event lift, an explicit coherent cut,
and the exact signed source contribution of each primitive prime.
A prime insertion has a precise returned-memory correction.
Its infinitesimal version gives symmetric mixed prime responses
and a negative-semidefinite response Hessian.

The arithmetic source is unchanged. Every active prime power,
the Gamma boundary, the pole square and the scalar subtraction
are retained. The event probabilities are positive; the completed
Weil form is a signed balance of those events.

**Status: RH remains unproved.** The source-matching event
construction and memory-response identities below have written
proofs. No new positive window, all-width matrix sign, numerical
certificate or physical probability law is asserted. No scientific
computation, numerical source evaluation, independent agent,
replay or workflow dispatch is performed.

## Source order and notation

The inputs are the actual wave form in
[CWF](RH_CRITICAL_WAVE_FORMS.md), the root-normalized transport
and prime resolvent in [PTE-T1](RH_PRIME_TRANSPORT_ENERGY.md),
and the domain-preserving cut and positive tail in ASM.
The public native
[Euler-product construction, F00-I](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/theorems/foundation/F00I_NATIVE_DIRICHLET_ZETA_EULER_PRODUCT.md)
supplies prime powers and their logarithmic weights.
The [clock-free cut calculus, T24](https://github.com/Parveen117/Recognition-Kernel-Framework/blob/e8e3089745bf7dc0532d01c62a5f660f5658c6f0/theorum/24_clock_free_recognition_seam_cut_calculus.md)
supplies the recognition/memory split and composition residue.

The positive norm and completion used here are those already
derived for PTE's transport and difference energies. The new
record space is the square-summable completion of the declared
prime-power labels. Positivity of the signed arithmetic source
is not used to construct this event carrier.

Write
\[
(U_pg)(x)=g(x-\log p),\qquad
r_p=p^{-1/2},\qquad
\kappa_p=\frac{1+r_p}{1-r_p}.
\tag{1}
\]
The operator U_p is unitary on full-line L2. It acts before
restriction to the even finite window.

## PCM-T1. Exponential prime-power records have a lawful coherent cut

Fix a prime p and suppress its subscript temporarily. Declare
the potential and its normalized weights by
\[
\Phi_p(k)=\frac{k}{2}\log p,\qquad
w_{p,k}=\frac{e^{-\Phi_p(k)}}{\sum_{j\ge0}e^{-\Phi_p(j)}}
=(1-r)r^k,\qquad k\ge0.
\tag{2}
\]
On the record space with orthonormal basis e_k, define
\[
b_p=\sum_{k\ge0}\sqrt{(1-r)r^k}\,e_k,\qquad
J_pg=\sum_{k\ge0}\sqrt{(1-r)r^k}\,
e_k\otimes U_p^kg .
\tag{3}
\]
Then b_p has norm one and J_p is an isometry.

There are two specified readings of this same packet. The
orthogonal record cuts e_k resolve the exponent k, with weight
w_(p,k) for a normalized input. The coherent cut uses
\[
\begin{gathered}
\Pi_p=|b_p\rangle\langle b_p|\otimes I,\\
R_p=(b_p^*\otimes I)J_p
=(1-r_p)(I-r_pU_p)^{-1},\\
W_p=(I-\Pi_p)J_p.
\end{gathered}
\tag{4}
\]
Its recognized event is \(b_p\otimes R_pg\), and
\[
\boxed{
J_p^*J_p=I,\qquad
R_p^*R_p+W_p^*W_p=I.
}
\tag{5}
\]
Thus the cut retains an exact memory complement.

For a finite record depth K, let J_(p,K) be (3) with k<=K,
without renormalization, and put
\(R_{p,K}=(b_p^*\otimes I)J_{p,K}\). Then
\[
\boxed{
\|(J_p-J_{p,K})g\|^2=r_p^{K+1}\|g\|^2,\qquad
\|R_p-R_{p,K}\|\le r_p^{K+1}.
}
\tag{6}
\]
These are explicit completion tails for this prime-power tower.

**Proof.** The geometric sum of w_(p,k) equals one.
Orthogonality of the record coordinates and unitarity of U_p
therefore give \(\|J_pg\|^2=\|g\|^2\).
Pairing with b_p gives
\[
R_pg=(1-r_p)\sum_{k\ge0}r_p^kU_p^kg.
\]
This series converges in operator norm because r_p<1, and is
the displayed resolvent. The operator Pi_p is an orthogonal
projection. Pythagoras on its two event channels gives (5).

The omitted record squared norm is the geometric tail
\(\sum_{k>K}(1-r)r^k=r^{K+1}\).
The omitted return series has operator norm at most that same
sum. This proves (6). QED.

The probabilities in (2) belong to the specified exponent
reading. The coherent projection in (4) is a different reading
of the same packet; its probability for normalized g is
\(\|R_pg\|^2\). Exponent-only probabilities do not retain the
cross-correlations needed for this coherent return.

The labels p^k are supplied by arithmetic factorization.
This construction defines how a cut reads them; it does not
derive primality from measurement.

## PCM-T2. The coherent cut reproduces the actual prime source

PTE's prime observer is
\[
V_p=(1-r_p^2)^{1/2}(I-r_pU_p)^{-1}
=\sqrt{\kappa_p}\,R_p.
\tag{7}
\]
For real compact g, with
\(C_g(u)=\langle g,T_ug\rangle\), its exact norm is
\[
\boxed{
\|V_pg\|^2
=\|g\|^2+2\sum_{k\ge1}r_p^kC_g(k\log p).
}
\tag{8}
\]
Consequently the complete contribution of prime p is
\[
\boxed{
\begin{aligned}
q_p[g]
&=-\log p\sum_{k\ge1}p^{-k/2}C_g(k\log p)\\
&=\frac{\log p}{2}\bigl(\|g\|^2-\|V_pg\|^2\bigr)\\
&=\frac{\log p}{2}
\left(\|W_pg\|^2-(\kappa_p-1)\|R_pg\|^2\right).
\end{aligned}
}
\tag{9}
\]
For g supported in (-L,L), only k log p<2L contributes to
the correlation sum. Thus only finitely many primitive primes
are active, even though each normalized event packet has all
nonnegative powers.

Let \(\mathcal P_L=\{p\text{ prime}:\log p<2L\}\).
The full actual source has the positive-channel decomposition
\[
\boxed{
\begin{aligned}
q_L[g]&=\mathscr R_L[g]-\mathscr D_L[g],\\
\mathscr R_L[g]
&=\ell(g)^2+\mathcal E_\Gamma[g]
+\frac12\sum_{p\in\mathcal P_L}\log p\,\|W_pg\|^2,\\
\mathscr D_L[g]
&=-\frac{\Omega_0}{2}\|g\|^2
+\frac12\sum_{p\in\mathcal P_L}
\log p\,(\kappa_p-1)\|R_pg\|^2 .
\end{aligned}
}
\tag{10}
\]
Here \(\ell\) and \(\mathcal E_\Gamma\) have ASM's definitions,
and \(\Omega_0<0\) is its fixed Gamma scalar. Both displayed
channel forms are nonnegative.

**Proof.** Expanding the squared norm of the geometric
resolvent and grouping pairs of indices with difference k gives
(8): the diagonal sum is \(1/(1-r_p^2)\), and the two
off-diagonal sums give the paired correlations. The double
series is absolutely convergent by unitarity. This is also
PTE-T1's coefficientwise proof.

Equation (7) and the isometry identity (5) give the last line
of (9). Sum (9) over the active primes and add ASM's exact
Gamma, pole and scalar terms. This yields (10). QED.

The direct sum of the square channels in (10) constructs one
explicit event lift Z_L with a fixed sign involution j_L:
\[
q_L(g,h)=\langle Z_Lg,j_LZ_Lh\rangle .
\]
The positive event Gram is \(Z_L^*Z_L\); the arithmetic source
is \(Z_L^*j_LZ_L\). The normalization factor kappa_p is essential.
Contractivity of R_p does not make \(I-V_p^*V_p\) nonnegative.

For the finite record return, put \(V_{p,K}=\sqrt{\kappa_p}R_{p,K}\).
Equations (6)--(7) imply the outward form bound
\[
\left|
\frac{\log p}{2}
\bigl(\|V_pg\|^2-\|V_{p,K}g\|^2\bigr)
\right|
\le \log p\,\kappa_p r_p^{K+1}\|g\|^2 .
\tag{11}
\]
Indeed \(\|V_p\|,\|V_{p,K}\|\le\sqrt{\kappa_p}\) and
\(\|V_p-V_{p,K}\|\le\sqrt{\kappa_p}r_p^{K+1}\).
No finite record is silently renormalized into a different source.

### Finite prime probabilities and the infinite normalization boundary

For a finite prime set F, joint exponent readings have weights
\[
\boxed{
\mathbb P_F(n)=Z_F^{-1}n^{-1/2},\qquad
Z_F=\prod_{p\in F}(1-p^{-1/2})^{-1},
}
\tag{12}
\]
on naturals whose prime factors lie in F. The joint event lift
uses the product record amplitudes and the unitary
\(T_{\sum_{p\in F}k_p\log p}\). Orthogonality gives the product
probabilities, and unique factorization gives the number n.
Thus retaining every exponent gives exact reconstruction.

This finite-prime law does not converge to a probability law
on the naturals as F exhausts all primes. Once F contains every
prime <=M, the Euler expansion gives
\[
Z_F\ge\sum_{n=1}^M n^{-1/2}\ge\sqrt M.
\]
Hence each fixed n, and each fixed finite set of naturals,
has probability tending to zero. A probability limit on the
discrete naturals would therefore give zero mass to every
point, contradicting total mass one.

For exponent sigma>1, the analogous global law
\(n^{-\sigma}/\zeta(\sigma)\) is normalized by F00-I.
The source exponent 1/2 in (7)--(12) is handled here through
local prime packets and a finite active prime set. No global
critical-line probability normalization is used.

## PCM-T3. Prime events feed the actual Schur memory coefficientwise

Fix L and ASM's real even space H_L and form domain D_L.
Let \(\iota_L:H_L\to L^2(\mathbb R)\) be zero extension and put
\[
\begin{gathered}
\mathcal C_u=\iota_L^*\frac{T_u+T_{-u}}2\iota_L,\\
\mathcal P_p
=\log p\sum_{\substack{k\ge1\\k\log p<2L}}
p^{-k/2}\mathcal C_{k\log p},\qquad
s_p=\log p\sum_{\substack{k\ge1\\k\log p<2L}}p^{-k/2}.
\end{gathered}
\tag{13}
\]
These are bounded self-adjoint operators, with
\[
\boxed{
\|\mathcal P_p\|\le s_p,\qquad
\mathcal P_p=
\frac{\log p}{2}
\iota_L^*(\kappa_pR_p^*R_p-I)\iota_L .
}
\tag{14}
\]
Thus \(-\langle g,\mathcal P_pg\rangle=q_p[g]\).
The identity in (14) follows from the event construction,
not from a fitted matrix factor.

Introduce independent arithmetic coupling coordinates
\(t=(t_p)_{p\in\mathcal P_L}\in[0,1]^{\mathcal P_L}\):
\[
\boxed{
\mathcal T(t)=\mathcal T_\circ-\sum_p t_p\mathcal P_p,\qquad
q_t[g]=\ell(g)^2+\frac{\Omega_0}{2}\|g\|^2
+\mathcal E_\Gamma[g]-\sum_p t_p\langle g,\mathcal P_pg\rangle .
}
\tag{15}
\]
The operator \(\mathcal T_\circ\) is the Gamma--pole--scalar
source with no active prime term. All forms have domain D_L.
At t=1 this is the actual source; the other t are declared
diagnostic couplings.

Choose ASM's N so that
\[
\boxed{
d=F(L)+\frac{\Omega_0}{2}-\sum_p s_p
+\frac{e^{-L}}2H_{2N}>0.
}
\tag{16}
\]
Let \(P=P_N,\ Q=I-P\), and let U synthesize the first N even
Legendre modes. For every t in the whole coupling cube,
\[
q_t[h]\ge d\|h\|^2\qquad(h\in QD_L).
\]
Let C(t) be the operator of this restricted tail form, and define
\[
\begin{gathered}
A(t)=U^*\mathcal T(t)U,\qquad
B(t)=Q\mathcal T(t)U,\qquad C(t)\ge dI,\\
\mathsf S(t)=A(t)-B(t)^*C(t)^{-1}B(t),\qquad
E(t)=U-C(t)^{-1}B(t).
\end{gathered}
\tag{17}
\]
Here E(t)a is the unique minimizing source state over the
observation fibre \(Pg=Ua\). In particular,
\[
q_t[E(t)a]=a^{\mathsf T}\mathsf S(t)a,\qquad
\mathcal T(t)E(t)=U\mathsf S(t),\qquad PE(t)=U.
\tag{18}
\]

The prime-resolved coupling columns are
\[
B_\circ=Q\mathcal T_\circ U,\qquad
B_p=Q\mathcal P_pU,\qquad
B(t)=B_\circ-\sum_p t_pB_p .
\tag{19}
\]
Set \(G_i(t)=C(t)^{-1/2}B_i\), for i equal to circle or an
active prime. Then the complete returned memory is
\[
\boxed{
B(t)^*C(t)^{-1}B(t)
=\left(G_\circ-\sum_p t_pG_p\right)^*
 \left(G_\circ-\sum_p t_pG_p\right).
}
\tag{20}
\]
Every Gamma--prime and prime--prime cross-term is retained.

**Proof.** Each paired translation in (13) is a self-adjoint
contraction after extension and restriction, proving the norm
bound. Equation (8) proves (14) by polarization; the translated
supports are disjoint when k log p>=2L.

The bounded prime operators preserve the common operator
domain under the perturbation (15). ASM-T1 applies with
\(-t_p\langle h,\mathcal P_ph\rangle\ge-t_ps_p\|h\|^2\).
Since \(0\le t_p\le1\), its lower bound is at least (16).
ASM-T2--T4 give (17) and the minimizing state. Stationarity
against every tail variation gives \(Q\mathcal T(t)E(t)=0\)
in the form sense. Pairing with the finite observed component
then yields \(\mathcal T(t)E(t)=U\mathsf S(t)\); its right side
is in H_L, so the minimizing columns belong to the full operator
domain. This proves (18), including domain preservation.

Linearity of the source gives (19). Apply C(t)^(-1/2) to
that identity and take its Gram to get (20). QED.

In particular, adding the individual positive Grams
\(G_p^*G_p\) is not the total memory: (20) also contains
\(-t_p(G_\circ^*G_p+G_p^*G_\circ)\) and
\(t_pt_q(G_p^*G_q+G_q^*G_p)\).
Moreover, C(t)^(-1/2) itself depends on all active prime
couplings. Independence of the arithmetic exponent records
does not assert independence of these returned source channels.

The event memory W_p lives in the prime record carrier.
The source-tail memory lives in QH_L. Equations (14), (19)
and (20) provide the explicit connection between them; the
two spaces are not identified by notation.

## PCM-T4. Exact prime insertion includes a positive memory correction

Let t and t+Delta both lie in the coupling cube. Write
\[
\mathcal V_\Delta=\sum_p\Delta_p\mathcal P_p,\qquad
D_\Delta(t)=Q\mathcal V_\Delta E(t),\qquad
C_+=C(t+\Delta).
\tag{21}
\]
The exact update is
\[
\boxed{
\mathsf S(t+\Delta)
=\mathsf S(t)-E(t)^*\mathcal V_\Delta E(t)
-D_\Delta(t)^*C_+^{-1}D_\Delta(t).
}
\tag{22}
\]
In particular, adding one prime by an amount delta gives
\[
\boxed{
\mathsf S(t+\delta e_p)
=\mathsf S(t)-\delta E(t)^*\mathcal P_pE(t)
-\delta^2D_p(t)^*C(t+\delta e_p)^{-1}D_p(t),
}
\tag{23}
\]
where \(D_p(t)=Q\mathcal P_pE(t)\).

The last term is positive semidefinite before subtraction.
It obeys the explicit enclosure
\[
\boxed{
\begin{aligned}
\mathsf S(t)-E^*\mathcal V_\Delta E
-\frac1dD_\Delta^*D_\Delta
&\preceq\mathsf S(t+\Delta)\\
&\preceq\mathsf S(t)-E^*\mathcal V_\Delta E .
\end{aligned}
}
\tag{24}
\]
For a nonzero single-prime step, equality with the update
omitting memory holds exactly when \(D_p(t)=0\). On one
observation vector a, it holds exactly when \(D_p(t)a=0\).

The minimizing source state updates by
\[
\boxed{
E(t+\Delta)=E(t)+C_+^{-1}D_\Delta(t).
}
\tag{25}
\]

**Proof.** Every state in the fibre Pg=Ua has the form
\(E(t)a+h\), with h in QD_L. By (18), its source at t+Delta is
\[
\begin{aligned}
q_{t+\Delta}[E(t)a+h]
={}&a^{\mathsf T}
\bigl(\mathsf S(t)-E(t)^*\mathcal V_\Delta E(t)\bigr)a\\
&-2\langle D_\Delta(t)a,h\rangle
+q_{t+\Delta}[h].
\end{aligned}
\]
Complete the positive C_+ square. Its unique minimum occurs
at \(h=C_+^{-1}D_\Delta(t)a\). The minimum is the right side
of (22), and its minimizing state gives (25).

The bounds \(0<C_+^{-1}\le d^{-1}I\) give (24).
The inverse is injective and positive, so
\(D^*C_+^{-1}D=0\) if and only if D=0. This proves the equality
criteria and the single-prime statement. QED.

This answers the prime-cut matching question for the actual
source: the two sides match with the returned-memory term.
Prime labeling alone does not force that term to vanish.

The first-order term \(E^*\mathcal P_pE\) can have either
sign. Thus (23) does not say that adding a prime always lowers
the reduced source. It identifies a definite negative correction
relative to the frozen-state first-order prediction.

## PCM-T5. Prime response has reciprocity and an exact memory Hessian

The finite matrix S(t) and the finite-column lift E(t) are real
analytic on a neighborhood of the coupling cube. Their derivatives
satisfy
\[
\boxed{
\partial_p E(t)=C(t)^{-1}D_p(t),\qquad
\partial_p\mathsf S(t)=-E(t)^*\mathcal P_pE(t).
}
\tag{26}
\]
The mixed response is
\[
\boxed{
\partial_q\partial_p\mathsf S(t)
=-D_q(t)^*C(t)^{-1}D_p(t)
-D_p(t)^*C(t)^{-1}D_q(t).
}
\tag{27}
\]
In particular,
\(\partial_q\partial_p\mathsf S=\partial_p\partial_q\mathsf S\).

For a real direction v in coupling space, put
\(\mathcal V_v=\sum_pv_p\mathcal P_p\) and
\(D_v=Q\mathcal V_vE(t)\). Then
\[
\boxed{
\frac{d^2}{ds^2}\mathsf S(t+sv)\bigg|_{s=0}
=-2D_v^*C(t)^{-1}D_v\preceq0.
}
\tag{28}
\]
Equivalently, for the effective potential
\(\Phi_t(a)=a^{\mathsf T}\mathsf S(t)a\),
\[
\boxed{
\frac{d^2}{ds^2}\Phi_{t+sv}(a)\bigg|_{s=0}
=-2\|C(t)^{-1/2}D_v(t)a\|^2 .
}
\tag{29}
\]

**Proof.** The operators C(t) are bounded perturbations of
one self-adjoint tail operator with a common domain. The
uniform positive bound (16) allows a local resolvent
Neumann expansion at every point of the cube. There are
finitely many prime coordinates; a sufficiently small open
neighborhood still has a positive tail bound. Therefore
C(t)^(-1), E(t), and S(t) are real analytic in their displayed
bounded-map norms.

Differentiating the inverse gives
\[
\partial_p C(t)^{-1}
=C(t)^{-1}(Q\mathcal P_pQ)C(t)^{-1}.
\]
Differentiate E in (17) and use (19) to obtain
\(\partial_pE=C^{-1}Q\mathcal P_pE=C^{-1}D_p\).
The first derivative of (23), or direct differentiation of
the Schur complement, gives the derivative of S in (26).
Differentiating that formula once more gives (27).
Taking the same direction in both slots gives (28), and
evaluating on a proves (29). QED.

Consequently S is matrix concave on the coupling cube:
\[
\mathsf S((1-\theta)t+\theta u)
\succeq(1-\theta)\mathsf S(t)+\theta\mathsf S(u),
\qquad 0\le\theta\le1.
\tag{30}
\]
This also follows because \(\Phi_t(a)\) is the infimum,
over one fixed observation fibre, of functions affine in t.

Equation (27) is a precise Maxwell-type reciprocity for the
declared prime-coupling responses. Equation (29) defines the
response Hessian used here. No identification with a physical
thermodynamic coefficient or with the framework's separate
connection-curvature tensor is required.

Concavity and reciprocity are structural constraints. They do
not supply the endpoint sign \(\mathsf S(\mathbf1)\succeq0\).

## PCM-T6. Complete prime ledgers are independent of insertion order

Choose a sequence
\(t^{(0)},t^{(1)},\ldots,t^{(m)}\) in the cube, with
\(\Delta^{(j)}=t^{(j+1)}-t^{(j)}\).
Let \(E_j=E(t^{(j)})\),
\(D_j=Q\mathcal V_{\Delta^{(j)}}E_j\), and
\(C_{j+1}=C(t^{(j+1)})\). Then
\[
\boxed{
\begin{aligned}
\mathsf S(t^{(m)})-\mathsf S(t^{(0)})
={}&-\sum_{j=0}^{m-1}
 E_j^*\mathcal V_{\Delta^{(j)}}E_j\\
&-\sum_{j=0}^{m-1}
 D_j^*C_{j+1}^{-1}D_j .
\end{aligned}
}
\tag{31}
\]
The complete sum depends only on the endpoints. Individual
prime contributions and individual memory corrections can
depend on the insertion order.

In particular, primes may be inserted one at a time from
t=0 to t=1. Equation (31) then gives the same actual ASM
matrix for every ordering, provided E_j and the returned
memory are updated at every step.

**Proof.** Apply (22) at each step and telescope.
There is one source operator (15) at each endpoint and one
unique fibre minimum (17), so the total is independent of
the chosen ordering. No commutation of the compressed prime
operators is needed. QED.

Dropping the memory sum, or using the same frozen E for every
step, is not the exact identity (31). Conversely, exact
endpoint agreement is an algebraic consistency test; it does
not by itself prove that the common endpoint is positive.

The sufficient finite-step inequality furnished by (24)
is concrete: if
\[
\mathsf S(t)-E^*\mathcal V_\Delta E-d^{-1}D_\Delta^*D_\Delta
\succeq0,
\]
then the next exact matrix is nonnegative. Such a bound may be
used where its entries and the current state lift are enclosed.
No sign of all intermediate matrices is assumed or required
by the exact transport formulas.

## Connection to the RH sign obligation

At the endpoint t=1, this note uses exactly ASM's matrix:
\[
\boxed{
\mathsf S(\mathbf1)=\mathsf S_{L,N},\qquad
q_L\ge0\Longleftrightarrow\mathsf S(\mathbf1)\succeq0.
}
\tag{32}
\]
Thus the all-width obligation is unchanged. Choosing one
admissible N(m) for each integer width m gives
\[
\mathrm{RH}
\quad\Longleftrightarrow\quad
\mathsf S_{m,N(m)}(\mathbf1)\succeq0
\quad\text{for every integer }m\ge1.
\]
It is only the actual endpoint, not every diagnostic point in
the coupling cube, that enters this criterion.

The new constructive information is:

| Object | Exact source-matching statement |
|---|---|
| Exponential prime-power packet | Normalized weights, unitary branch transport and explicit geometric completion tail |
| Coherent prime cut | Its return reconstructs PTE's resolvent with the exact factor kappa_p |
| Actual prime source | Signed recognized/memory balance (9), with every prime power retained |
| Completed arithmetic event lift | Positive and negative channels in (10), including the Gamma scalar subtraction |
| ASM memory | Prime-resolved columns and all cross-terms in (20) |
| One prime insertion | Frozen-state contribution plus the exact returned-memory correction (23) |
| Prime response | Symmetric mixed derivatives and the negative-square Hessian (27)--(29) |
| Multiple insertions | Exact endpoint agreement after the complete memory ledger (31) |

The local event probabilities at exponent 1/2 are normalized.
The infinite product over all primes is not a probability
normalization on the naturals at that exponent. Neither local
normalization nor exact prime-order agreement replaces the
actual endpoint inequality.

The event memory and the source-tail inverse have different
roles. Equations (14)--(20) connect them explicitly; no
probability scalar is substituted for \(B^*C^{-1}B\).
The latter remains a finite matrix with an infinite-dimensional
tail inverse until it is enclosed. ASM-T5's finite trial and
residual bounds remain available for that task.

The known TPW restricted margin \(1/500\) for \(L\le11/20\)
and PTE unrestricted margin \(1/16\) for \(L\le1/32\) remain
unchanged. No new numerical matrix or larger positive window
is claimed here.

Geometric resolvents, orthogonal event projections, Schur
complements and variational response differentiation are
standard tools. This continuation assembles them with the
actual arithmetic coefficients and declared native cuts.
It makes no literature-wide novelty claim.

Only this note and its README entry are added. Sources from
other repositories are public and pinned to RKF commit
e8e3089745bf7dc0532d01c62a5f660f5658c6f0. No private review
material is transferred. These are written proofs, not newly
machine-checked certificates. The documentation commit
contains [skip ci], is stacked on draft PR #113, and is
not merged.
