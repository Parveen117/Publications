# LCT-1 — Layer-tower cut transfer with a pre-supplied correction factor

17 October 2026. Written continuation on branch `layer-cut-transfer-20261017`,
based on the draft-chain head 167eee544a2c5cfb5cef5b8032c82ac26a4b2e0d (PR #125).
It reuses ILC-1, LCW-1, LAT-1, LDC-1, LOT-1, the finite observer detector,
CLA-1, RLO-1, RCD-1 and FWT-1 without changing them.

The request: when a cut is placed on the layer tower at a reference layer,
all information should pass through it. If the cut forces a correction
factor, that factor should be supplied in advance.

This note proves that this is possible in an exact sense, and states its price.

- The correction factor is K_N(λ)=λ^{2N}/Q_N(λ)^2. It is built only from the
  retained readings c_0,…,c_{2N−1} at the plane.
- The whole tower above the cut is carried, without loss, by one function E_N.
  E_N is a bilinear source pairing and starts with the boundary norm h_N.
- Moving the cut by one layer has its own pre-supplied factor
  g_N=Q_N/Q_{N+1}. It needs only the next pair (c_{2N},c_{2N+1}).
- The sign of E_N is not supplied. For the actual xi source, positivity of
  h_N at every depth is the existing RH criterion. No new sign is claimed.

**Status.** LCT-T1 to LCT-T6 have written proofs. Finite exact regression
witnesses pass: 389 checks and 65 negative controls, in [verify.cjs](verify.cjs).
The witnesses use three finite sources over Q+iQ: positive, indefinite real,
and a conjugate-pair comparison source. They are not evidence about the xi
source and do not certify any all-depth statement. RH remains unproved.
No physical time, speed, mass, field or curvature identification is made.

## Source contract and typed roles

Retain RLO-1 (1) and CWR-1's source class. A source is a normalized bilinear
functional

\[
\mathcal L[g]=\sum_\mu w_\mu\,g(u_\mu),\qquad
\sum_\mu w_\mu=1,\qquad \sum_\mu|w_\mu|<\infty,\qquad 0<|u_\mu|<1 ,
\tag{1}
\]

with conjugation-symmetric points and weights. Weights and points may be
complex. Every product inside L stays a bilinear product; nothing is replaced
by an absolute square. The layers and response are

\[
c_n=\mathcal L[y^n],\qquad
f(\lambda)=\mathcal L\!\left[\frac1{1+\lambda y}\right]
=\sum_{n\ge0}(-\lambda)^n c_n\quad(|\lambda|<1),
\tag{2}
\]

continued meromorphically wherever the sum converges. For the actual source,
u_μ=a/(a+μ) and w_μ=ν_μ/(q(a)(a+μ)). Then f(λ)=q(a(1+λ))/q(a), and the c_n
are exactly ILC-1's layers c_n=(−a)^n q^{(n)}(a)/(n! q(a)).

| Symbol | Role | Not identified with |
|---|---|---|
| a | the reference anchor: the plane where the cut is placed | a time or a zero |
| λ=(t−a)/a | displacement from the plane; λ>0 above, λ<0 below | a TVSP response component (ILC-1) |
| t | absolute source position | a spatial coordinate |
| N | cut depth: the retained readings are c_0,…,c_{2N−1} | a bank size at another anchor |
| τ=log ρ | anchor-flow label, b=ρa | a physical clock |
| H_N=(c_{i+j}) | unshifted Hankel bank (finite detector, D3) | a positive Gram matrix unless proved |

## LCT-T1. The universal generator of the tower

For every source of the form above, write f_a(λ)=q(a(1+λ))/q(a). Then

\[
\boxed{\;\partial_\tau f_a=\mathcal G f_a,\qquad
\mathcal G=(1+\lambda)\partial_\lambda+c_1(a),\qquad \tau=\log a.\;}
\tag{3}
\]

On coefficients this is LAT-1's hierarchy
dc_n/dτ=(n+c_1)c_n−(n+1)c_{n+1}. In the source chart it is the vector field
Du=u(1−u), Dw=(c_1−u)w, which keeps Σw=1. The finite flow is the Möbius chart
change

\[
u\mapsto\frac{\rho u}{1+(\rho-1)u},\qquad
w\mapsto\frac{w}{1+(\rho-1)u}\Big/\frac{q(\rho a)}{q(a)},
\tag{4}
\]

with normalization cocycle q(ρa)/q(a)=exp(−∫_0^{log ρ}c_1 dτ).

**Proof.** Differentiate f_a at fixed λ:
a(1+λ)q'(a(1+λ))/q(a)=(1+λ)f_a'(λ), and a q'(a)/q(a)=−c_1.
Comparing coefficients of λ^n gives the hierarchy. For the source chart,
D log w_μ=−D log q(a)−a/(a+μ)=c_1−u_μ and Du_μ=u_μ(1−u_μ). Then
Σ_μ Dw_μ=c_1−c_1=0. Substituting b=ρa gives
(b+μ)/(a+μ)=1+(ρ−1)u_μ, which proves (4). QED.

In the absolute coordinate t, (3) is the Euler field t∂_t acting on the
observer position. The source points μ and multiplicities ν do not move.
The tower is the jet of the potential Φ=log X along this flow: q=Φ', and
c_n are its normalized higher derivatives.

The generator is universal in a precise sense: it is the same for every
analytic q. The source enters only through the scalar c_1 and the initial
tower. By CLA-1, such source-independent structure cannot supply the sign of
any source form.

The analogy with the TVSP potential tower (potential, responses, layers) is
structural only. SPT-1 proves that the TVSP corners share only their contact
one-form. ILC-1 does not identify this λ with a TVSP response.

## LCT-T2. The cut at depth N and its pre-supplied factor

Assume det H_N≠0. Let P_N(y)=Σ_k p_k y^k be the monic polynomial with
L[P_N y^k]=0 for k<N. Define

\[
Q_N(\lambda)=(-\lambda)^N P_N(-1/\lambda),\qquad
R_N=\bigl[Q_N f\bigr]_{\deg<N},\qquad
E_N(\lambda)=\mathcal L\!\left[\frac{P_N(y)^2}{1+\lambda y}\right].
\tag{5}
\]

Q_N(0)=1 and deg Q_N=N. If P_N has roots y_j, then Q_N=∏(1+λy_j), but no
root is needed. Then, as meromorphic functions on the common domain,

\[
\boxed{\;f(\lambda)=\frac{R_N(\lambda)}{Q_N(\lambda)}
+K_N(\lambda)\,E_N(\lambda),\qquad
K_N(\lambda)=\frac{\lambda^{2N}}{Q_N(\lambda)^2},\;}
\tag{6}
\]
\[
\boxed{\;E_N(0)=h_N:=\mathcal L[P_N^2]=\frac{\det H_{N+1}}{\det H_N}.\;}
\tag{7}
\]

f_N:=R_N/Q_N is the [N−1/N] Padé reading. When ILC-1's finite premise
A_N,B_N≻0 holds, it is ILC-T4's positive local model, by uniqueness of the
Padé approximant of that type.

**Proof.** For fixed λ≠0, divide in the variable y:
P_N(−1/λ)−P_N(y)=(−1/λ−y)S(y,λ). Here S has degree N−1 in y, and
λ^{N−1}S is polynomial in λ of degree at most N−1. Multiply by (−λ)^N and
divide by 1+λy:

\[
\frac{Q_N(\lambda)}{1+\lambda y}
=\frac{(-\lambda)^N P_N(y)}{1+\lambda y}-(-1)^N\lambda^{N-1}S(y,\lambda).
\]

Apply L. The last term is a polynomial in λ of degree below N, so

\[
Q_Nf-\widetilde R_N=(-\lambda)^N\mathcal L\!\left[\frac{P_N}{1+\lambda y}\right].
\]

Next multiply the same identity by P_N(y) before applying L. Orthogonality
removes the S term, since it has degree N−1 in y:

\[
Q_N(\lambda)\,\mathcal L\!\left[\frac{P_N}{1+\lambda y}\right]
=(-\lambda)^N E_N(\lambda).
\]

Combine the two displays: Q_Nf−R̃_N=λ^{2N}E_N/Q_N. The right side is
O(λ^{2N}), so R̃_N is the truncation of Q_Nf below degree N, which is R_N.
Dividing by Q_N gives (6). Absolute convergence (1) justifies applying L to
the bounded rational functions on |λ|<1. The meromorphic continuation
extends the identity. For (7), L[P_N^2]=L[P_N y^N]. Cramer's rule on the
bordered Hankel system gives det H_{N+1}/det H_N. QED.

Relational reading. K_N depends only on the retained plane readings: its
coefficients are those of P_N, which solve the H_N system. It has the same
universal form λ^{2N}/Q_N^2 for every source. All the unretained information
sits in one function E_N, which is still a bilinear pairing of the actual
source. In absolute coordinates, (6) reads

\[
q(t)-q(a)f_N\!\left(\tfrac{t-a}a\right)
=q(a)\left(\frac{(t-a)^N}{a^N Q_N((t-a)/a)}\right)^{2}E_N\!\left(\tfrac{t-a}a\right).
\]

This is the confluent, single-plane form of RCD-1 (7). It is proved here
directly, without nodes, and it holds for indefinite and complex sources.

## LCT-T3. The tower passes the cut without loss

Put e_{N,k}=L[P_N^2y^k]. Then, for |λ|<1,

\[
E_N(\lambda)=\sum_{k\ge0}(-\lambda)^k e_{N,k},\qquad
e_{N,k}=\sum_{i,j=0}^N p_ip_j\,c_{i+j+k},
\tag{8}
\]
\[
c_{2N+k}=e_{N,k}-\sum_{(i,j)\ne(N,N)}p_ip_j\,c_{i+j+k}.
\tag{9}
\]

Hence the map

\[
(c_n)_{n\ge0}\;\longleftrightarrow\;
\bigl(c_0,\dots,c_{2N-1}\,;\,e_{N,0},e_{N,1},\dots\bigr)
\]

is a bijection. It is triangular with unit diagonal (p_N=1). Every layer above
the cut is carried by the single function E_N, and no layer is discarded.

**Proof.** Expand 1/(1+λy) under |u_μ|<1 and the absolute budget, which gives
(8). In (8) the coefficient of c_{2N+k} is p_N^2=1, and every other term uses
a lower layer. This gives (9) by induction on k. QED.

Two boundaries are part of the statement.

- **Completeness of the tower itself.** For analytic q, the full tower at a
  regular plane determines q by the identity theorem. This holds for the
  actual xi source, which is meromorphic (ILC-T2 proof, FWT-T2). For merely
  smooth potentials the tower is not complete. The flat deformation
  u+ηe^{−1/s²}v²/2 has the same jets and different curvature (JED-T5 on
  branch claude/jet-entry-development at 72d6530). The cut cannot transfer
  information that the tower does not hold.
- **Saturation.** If L has exactly m support points, then h_m=0, E_m≡0 and
  f=f_m. A cut deeper than m does not exist, because det H_{m+1}=0.

## LCT-T4. Moving the cut: the depth cocycle

Assume also det H_{N+1}≠0. Then

\[
\boxed{\;f_{N+1}-f_N=\frac{h_N\,\lambda^{2N}}{Q_NQ_{N+1}},\qquad
R_{N+1}Q_N-R_NQ_{N+1}=h_N\lambda^{2N},\;}
\tag{10}
\]
\[
\boxed{\;E_N=g_N\bigl(h_N+\lambda^2g_N\,E_{N+1}\bigr),\qquad g_N=\frac{Q_N}{Q_{N+1}}.\;}
\tag{11}
\]

For M>N with det H_M≠0,

\[
E_N=\frac{Q_N^2\,(f_M-f_N)}{\lambda^{2N}}+\Gamma_{M\leftarrow N}\,E_M,\qquad
\Gamma_{M\leftarrow N}=\lambda^{2(M-N)}\frac{Q_N^2}{Q_M^2}.
\tag{12}
\]

The factors compose: Γ_{K←M}Γ_{M←N}=Γ_{K←N}. The composition of the one-step
maps (11) is the direct map (12).

**Proof.** R_{N+1}Q_N−R_NQ_{N+1} is a polynomial of degree at most 2N. By (6)
at depths N and N+1 it equals Q_NQ_{N+1}(K_NE_N−K_{N+1}E_{N+1}), which is
h_Nλ^{2N}+O(λ^{2N+1}), because Q_N(0)=Q_{N+1}(0)=1 and E_N(0)=h_N. This proves
(10). Then (6) at both depths gives K_NE_N=(f_{N+1}−f_N)+K_{N+1}E_{N+1}.
Dividing by K_N gives (11). The same argument with M in place of N+1 gives
(12). The Γ are ratios, so they compose. QED.

Iterating (11) is the Stieltjes/Jacobi continued fraction, written as a
relational transfer law. The classical theory is credited, not claimed.

The data the one-step factor needs are exactly the next pair
(c_{2N},c_{2N+1}), the same pair that LCW-1's continuation windows test. The
ratio h_N/h_{N−1} is LDC-1's sole nonzero boundary eigenvalue. LOT-1's
continuation position is p_N=h_N/W_N. In words: the correction for moving the
cut one layer is known in advance as soon as the next pair is read. The
correction for the whole remaining tower is the single function E_{N+1}.

## LCT-T5. The reference plane, flatness and the two sides

Let the source be conjugation-symmetric. Then the c_n, P_N, Q_N and R_N are
real. For real λ in the domain with Q_N(λ)≠0:

1. K_N(λ)≥0, with equality only at λ=0, where the zero has exact order 2N.
   The cut reading is exact at the plane. The difference f−f_N is flat there
   to order 2N−1, and its first non-flat coefficient is h_N.
2. On both sides of the plane, above (λ>0) and below (λ<0),
   sign(f−f_N)=sign E_N(λ). Near the plane both sides carry the sign of h_N.
   The factor has even order, so it cannot tell the sides apart. Any asymmetry
   between the sides lives in E_N.
3. If ILC-1's premise A_N,B_N≻0 holds, then Q_N(λ)>0 for every λ>−1, because
   the roots y_j lie in (0,1).
4. Under RH the actual E_N(λ)=∫P_N^2/(1+λy)dν is positive for all λ>−1 and
   all N. So the actual correction is positive on both sides at every depth.
   Conversely, RH ⇔ h_N>0 for every N. This is the finite observer detector's
   D3, det H_N>0 for every N, restated through h_N=det H_{N+1}/det H_N. It is
   **not** a new criterion.

**Proof.** Real coefficients follow from conjugation symmetry, as the
retained notes prove. Part 1 holds because Q_N(0)=1 and K_N is a square
divided by a square. Part 2 follows from (6) and part 1. Part 3: the positive
finite model has nodes y_j in (0,1) (ILC-T4), so 1+λy_j>0. Part 4: under RH
the weights form a probability measure on (0,1) with infinite support
(ILC-T2), so P_N^2/(1+λy) is a nonnegative function that is not almost
everywhere zero. The converse is D3. QED.

The witnesses show why the plane alone does not decide the sides. The
conjugate-pair comparison source has every moment real. At N=1 it has
h_1>0 at the plane, but E_1<0 on the side λ<0. At N=2 it has h_2<0. The
indefinite real source has h_1<0. In both cases the universal factor
K_N≥0 is unchanged. The sign information is entirely in the pairing E_N, as
CLA-1 and CCI-1 predict.

The two cut channels have the same structure. Applying T2 to the channel
functionals L^+[g]=L[yg]/c_1 and L^−[g]=L[(1−y)g]/(1−c_1), which are the A_N and
B_N banks, gives the same exact factorization in each channel. No new proof is
needed; only the source changes.

## LCT-T6. Observer transport across a fixed cut

Fix the cut (a,N). In absolute coordinates, E_N(λ)=Σ_μ w_μP_N(u_μ)^2(a+μ)/(t+μ)
is a function of the position t alone. Hence an observer at b=ρa reads its own
layers as Taylor coefficients at t=b of

\[
\frac{q(t)}{q(b)}=\frac{q(a)}{q(b)}\Bigl[f_N+K_NE_N\Bigr]\!\left(\tfrac{t-a}a\right).
\tag{13}
\]

- The model part f_N transports exactly: its poles and residues are invariant
  under (4), by LAT-T2.
- The factor K_N is explicit.
- The unknown is the same single function E_N.

Nothing is refitted and nothing is lost. If the observer at b instead places
its own cut at depth M, it obtains a different model and a different E^b_M.
The two cuts are related exactly by the cross-anchor law

\[
q(a)K^a_NE^a_N=\bigl(q^b_M-q^a_N\bigr)+q(b)K^b_ME^b_M,
\tag{14}
\]

where q^a_N=q(a)f_N((t−a)/a), and likewise at b. The models fitted at two
anchors differ (LAA-1), and (14) keeps their difference as data rather than
identifying them.

**Proof.** Substitute t=a(1+λ) in (6) and multiply by q(a). Write (6) at b
in the same way. Subtracting the two expressions gives (14). QED.

## What the cut does and does not supply

| Question | Answer | Where |
|---|---|---|
| Does information above the cut survive? | Yes, completely, as the single function E_N; for analytic sources the tower is complete | T3 |
| Is the correction factor known in advance? | Yes: K_N=λ^{2N}/Q_N^2 from the retained readings | T2 |
| Is the factor for moving the cut known in advance? | Yes: g_N=Q_N/Q_{N+1} and h_N, from the next pair only | T4 |
| Does every observer get the information? | Yes: model + known factor + the same E_N | T6 |
| Is the value of E_N determined by finitely many readings? | No: CLA-T3 matches any finite prefix with a positive and with an indefinite source | CLA-1 |
| Is the sign of E_N supplied? | No; for the actual xi source, h_N>0 for all N is RH | T5, D3 |
| Is the flow a physical time evolution? | No: it is the observer's Euler flow, with the source fixed | T1 |

The cut is therefore lossless and its correction is pre-supplied. The remaining
open content is exactly one sign: positivity of the actual pairing E_N at the
plane, h_N>0, at every depth. No part of this note proves it.

## Finite exact witnesses

[verify.cjs](verify.cjs) uses only the unchanged native engine
(papers/native-common-cone-limit/native_engine/native_operator.cjs, Q+iQ).
Run it through Python 3.12 orchestration: `python3.12 run.py --check`.

Sources:

- S1: six positive rational points, the positive control.
- S2: an indefinite real source with one negative weight.
- S3: a conjugate pair (1±i)/2 with conjugate weights plus three real points.
  All of its moments are real.

| Group | What is checked exactly |
|---|---|
| T1 | generator hierarchy for n≤8 on all sources; Σ Dw=0; reality of S3 moments; Möbius chart and weight tilt on a positive pole source; cocycle q(b)/q(a) |
| T2 | the polynomial form of (6) multiplied by ∏(1+λu_μ), at every tested depth; Q_N(0)=1; degrees; reality; E_N(0)=h_N=det H_{N+1}/det H_N |
| T3 | (8) and (9) for k≤3; saturation h_6=0, E_6≡0 for S1 |
| T4 | (10), (11), the two-step composition against the direct (12), and order-2N vanishing of f_{N+2}−f_N |
| T5 | f−f_N=K_NE_N with K_N≥0 at λ=±1/2, ±1/3, 2, −3/4; sign(f−f_N)=sign E_N; all S1 corrections positive on both sides |

Negative controls (65) include:

- order 2N−1 instead of 2N;
- the unsquared factor λ^{2N}/Q_N;
- a non-orthogonal cut polynomial;
- a wrong innovation h_N+1/1000;
- the cut without its correction (f≠f_N);
- a cut deeper than the support;
- untilted transported weights;
- w|P(u)|^2 in place of the bilinear wP(u)^2;
- the negative corrections of S2 and S3.

Counts and the certificate digest are in [CERTIFICATE.json](CERTIFICATE.json)
and [EXPECTED.sha256](EXPECTED.sha256).

## Attribution and scope

Padé approximants, the Gauss–Christoffel remainder, orthogonal polynomials,
Hankel determinants and Stieltjes continued fractions are classical. See
Stieltjes (1894) and Baker & Graves-Morris, *Padé Approximants* (2nd ed., 1996).
The contribution here is their node-free, source-bilinear assembly as an
exact cut-transfer law on the retained layer tower, with the typed roles, the
pre-supplied factors, the observer-transport form, and the explicit statement
of which part is open. No literature-wide novelty claim is made.

No predecessor packet, engine, certificate, workflow or shared branch is
changed. This branch is not a pull request.
