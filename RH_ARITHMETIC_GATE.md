# RH proof attempt: exact arithmetic gate and remaining gap

17 October 2026. Research continuation from Publications commit
`6a5c4a9e780b69f04fb3e2f5740d458eb44fd403`.

**Status: the RH proof is incomplete.** This note derives an exact
all-order arithmetic expression for the existing layer criterion and
identifies a failed primewise positivity argument. It proves neither the
required universal signs nor their negation for the actual xi source.
The finite certificates in the predecessor packets remain unchanged.

## 1. Retain the actual source

Write
\[
\mathcal X(w^2)=\xi(1/2+w),\qquad
q(t)=\mathcal X'(t)/\mathcal X(t),\qquad
w_0=57/2,\quad a=w_0^2,\quad s_0=29.
\]
The retained source enclosure proves \(q(a)>0\). Set
\[
c_n=\frac{(-a)^n q^{(n)}(a)}{n!\,q(a)},\qquad c_0=1.
\]
The displacement \(\lambda=(t-a)/a\), derivative depth \(n\), and the
signed thermodynamic deformation \(\mu\) in VSR are different variables.

[ILC-T2](papers/native-infinite-layer-cut/THEOREM.md) proves
\[
\mathrm{RH}\quad\Longleftrightarrow\quad
D_{r,k}:=\sum_{j=0}^{k}(-1)^j\binom{k}{j}c_{r+j}\ge0
\quad\hbox{for every }r,k\ge0.
\tag{1}
\]
Equivalently both \(A_N=(c_{i+j+1})_{i,j<N}\) and
\(B_N=(c_{i+j}-c_{i+j+1})_{i,j<N}\) are strictly positive for every
\(N\ge1\). These are proved equivalences; their universal premises
have not been established.

[SE49](papers/native-moving-zero-covariance/sources/SE49.md) separately
certifies the integer observation bank through \(V_{24}\).
[LAT](papers/native-layer-anchor-transport/THEOREM.md) and
[LAA](papers/native-layer-anchor-atlas/THEOREM.md) retain actual
\(c_0,\ldots,c_{19}\) and ten-dimensional positive layer forms at three
anchors. These two finite bank sizes count different directions.

## 2. An unconditional arithmetic formula

Let \(\Lambda(m)\) be the von Mangoldt function: \(\log p\) when
\(m=p^j\), and zero otherwise. For real \(t>1/4\), write \(w=\sqrt t\).
Logarithmic differentiation of the completed xi formula and the
absolutely convergent Euler product gives
\[
\boxed{
q(t)=\frac1{t-1/4}
+\frac{\psi(1/4+\sqrt t/2)-\log\pi}{4\sqrt t}
-\frac1{2\sqrt t}\sum_{m\ge2}\Lambda(m)m^{-1/2-\sqrt t}.
}
\tag{2}
\]
Here \(\psi=\Gamma'/\Gamma\).

**Proof.** For \(s=1/2+w>1\),
\[
\frac{\xi'}{\xi}(s)=
\frac1s+\frac1{s-1}-\frac12\log\pi
+\frac12\psi(s/2)-\sum_{m\ge2}\Lambda(m)m^{-s}.
\]
Divide by \(2w\); the first two rational terms combine to
\(1/(t-1/4)\). All subsequent differentiations at \(a\) take place
inside the Euler-product half-plane. QED.

The definition of xi is [NIST DLMF 25.4.4](https://dlmf.nist.gov/25.4#E4);
the Mangoldt series is [DLMF 27.4.12](https://dlmf.nist.gov/27.4#E12).
Equation (2) uses no critical-strip zero assumption.

## 3. Exact prime-power layers at arbitrary finite depth

Define polynomials by
\[
P_0(z)=1,\qquad
P_{n+1}(z)=(z+2n+1)P_n(z)-zP'_n(z).
\tag{3}
\]
Their first values are
\[
P_1=z+1,\quad P_2=z^2+3z+3,\quad
P_3=z^3+6z^2+15z+15.
\]
For \(b\ge0\),
\[
\boxed{
\frac{(-1)^n}{n!}\frac{d^n}{dt^n}
\left(\frac{e^{-b\sqrt t}}{\sqrt t}\right)
=
\frac{e^{-b\sqrt t}}{2^n n!\,t^{n+1/2}}P_n(b\sqrt t).
}
\tag{4}
\]

**Proof.** The case \(n=0\) is immediate. Apply
\(-\frac1{n+1}\frac{d}{dt}
=-\frac1{2(n+1)w}\frac{d}{dw}\) to the right side and collect
powers of \(w\). This gives exactly (3) and the next denominator.
If \(P_n=\sum_{j=0}^n p_{n,j}z^j\), then
\[
p_{n+1,j}=(2n+1-j)p_{n,j}+p_{n,j-1},
\]
with coefficients outside their ranges set to zero. Induction proves
that every coefficient in the polynomial's range is positive. QED.

Set
\[
\begin{aligned}
R_n&=\frac{a^n}{(a-1/4)^{n+1}},\\
G_n&=\left.\frac{(-a)^n}{n!}\frac{d^n}{dt^n}
  \frac{\psi(1/4+\sqrt t/2)-\log\pi}{4\sqrt t}\right|_{t=a},\\
E_n&=\frac1{2^{n+1}n!\,w_0}
 \sum_{m\ge2}\Lambda(m)m^{-s_0}P_n(w_0\log m).
\end{aligned}
\tag{5}
\]
Then for every fixed \(n\ge0\),
\[
\boxed{q(a)c_n=R_n+G_n-E_n.}
\tag{6}
\]
The entire prime-power sum, including all omitted tails of any numerical
truncation, is present. For each fixed \(n\), its absolute convergence
follows from \(\Lambda(m)\le\log m\) and convergence of
\(\sum m^{-s_0}(1+\log^{n+1}m)\). On a sufficiently small complex
neighbourhood of \(a\), the same estimates with a fixed real exponent
greater than one justify termwise differentiation. No interchange with
the limit \(n\to\infty\) is asserted.

## 4. Positive scalar prime terms do not give a positive moment matrix

For one fixed prime-power term let \(z=w_0\log m\). Divide its layer
sequence by its positive zeroth layer. Equations (4)-(5) give
\[
k_0=1,\qquad k_1=(z+1)/2,\qquad k_2=(z^2+3z+3)/8.
\]
Consequently
\[
\boxed{
\det\begin{pmatrix}k_0&k_1\\k_1&k_2\end{pmatrix}
=\frac{1-z-z^2}{8}.
}
\tag{7}
\]
At the actual anchor, every \(m\ge2\) has
\(z\ge(57/2)\log2>57/4>1\), using
\(\log2=\int_1^2 dx/x>1/2\). Thus this determinant is strictly
negative. Its symmetric matrix is indefinite. Already
\(k_0-k_1=(1-z)/2<0\) fails the first compact-moment difference.

This is an exact symbolic calculation, not a floating-point experiment.
It refutes the attempted inference
“positive prime-power weights imply positive layer moment matrices.”
It does **not** prove any negative sign for the full actual source:
the Euler contribution is subtracted in (6), and the Gamma and rational
terms remain essential. Nor does it assert that a sum of these
indefinite individual matrices must be indefinite.

The rational component also fails compact-moment positivity by itself.
With \(\rho=a/(a-1/4)>1\),
\[
R_n=R_0\rho^n,\qquad R_0-R_1=R_0(1-\rho)<0.
\tag{8}
\]
The full source therefore cannot be justified by separately declaring
each displayed contribution to be a positive compact-moment measure.

## 5. The precise arithmetic inequality still required

For a sequence \(h\), let
\[
\Delta_{r,k}h=\sum_{j=0}^{k}(-1)^j\binom{k}{j}h_{r+j}.
\]
Equations (6) and (8) give the unconditional identity
\[
\boxed{
q(a)D_{r,k}
=\frac{\rho^r(1-\rho)^k}{a-1/4}
+\Delta_{r,k}G-\Delta_{r,k}E.
}
\tag{9}
\]
The prime expression in the last term is explicitly
\[
\Delta_{r,k}E
=\frac1{2w_0}\sum_{m\ge2}\Lambda(m)m^{-s_0}
\sum_{j=0}^{k}(-1)^j\binom{k}{j}
\frac{P_{r+j}(w_0\log m)}{2^{r+j}(r+j)!}.
\tag{10}
\]
Every sum here is legitimate for fixed \(r,k\).

**Unproved arithmetic gate.** Establish, for the actual Gamma and
Mangoldt terms in (5), that
\[
\boxed{
\Delta_{r,k}G-\Delta_{r,k}E
\ \ge\
-\frac{\rho^r(1-\rho)^k}{a-1/4}
\qquad\hbox{for all }r,k\ge0.
}
\tag{11}
\]
Because \(q(a)>0\), (11) is exactly (1), expressed in the retained
arithmetic source. It is an equivalent RH obligation, not a weaker
lemma already proved by this note.

The standard completion after (11) is available: Hausdorff moment
reconstruction supplies a probability measure \(\nu\) on \([0,1]\),
so the actual normalized germ equals
\[
\frac{q(a(1+\lambda))}{q(a)}
=\int_{[0,1]}\frac{d\nu(y)}{1+\lambda y}.
\]
The right side is holomorphic off \((-\infty,-1]\). Equality of the
germs and meromorphic continuation identify it with the actual
logarithmic derivative. Every zero of \(\mathcal X\) produces a
nonremovable pole of that derivative, so all its zeros must lie on
the nonpositive real \(t\)-axis. This implies RH. This paragraph is
a conditional proof, with (11) as its still-missing hypothesis; the
full details and converse are retained in ILC-T2.

The shrinking window law in
[LCW-T3](papers/native-layer-continuation-window/THEOREM.md) does not
establish (11). It states that already accepted successive windows
have \(W_{\mathrm{next}}\le W/16\), without a lower bound on the
actual source's distance to either endpoint. An upper error bound
or a finite positive prefix does not supply that distance at all orders.

## 6. What reflection does and does not transport

The functional equation gives
\[
\xi(1/2+w)=\xi(1/2-w).
\]
It makes \(\mathcal X(w^2)\) well defined but leaves \(t=w^2\) unchanged.
It therefore does not implement \(\lambda\mapsto-\lambda\) about the
positive anchor \(a\).

Likewise VSR's proved
\(\mathcal G_{-\mu}=D\mathcal G_\mu D\) concerns its specified
thermodynamic cubic source. No identity intertwining that source with
the xi forms \(A_N,B_N\) has been established. Even when a real
symmetric form has a reflected congruence \(M_-=D^TM_+D\) with
invertible \(D\), its positive and negative direction counts agree.
A reflection can preserve positivity already present; that identity
alone cannot establish positivity of either side.

A finite-data example makes the missing source identity explicit.
Take any finitely many distinct real regular anchors \(a_j\), with
\(q(a_j)\ne0\), and
derivative orders \(m_j\ge0\). Put
\[
R(t)=\prod_j(t-a_j)^{m_j+2},\quad
P(t)=1+\epsilon R(t)^2,\quad \epsilon>0,\quad
\widetilde{\mathcal X}(t)=\mathcal X(t)P(t).
\]
Then \(P>0\) on the real axis, and
\[
\widetilde q-q=P'/P
=O((t-a_j)^{2m_j+3})
\]
at each anchor. Thus all prescribed finite response jets, their
normalized layers, and every matrix built from those jets agree
**exactly**, including their normalization \(q(a_j)\).
Yet \(P\) is a nonconstant real polynomial without real roots, so it
has nonreal roots; they are zeros of \(\widetilde{\mathcal X}\).
The entire function
\(\widetilde\xi(s)=\widetilde{\mathcal X}((s-1/2)^2)\) retains reflection
and conjugation symmetry while having zeros off the critical line.

This constructed function is **not** the Riemann xi function and
does not preserve its Euler arithmetic. It is not a counterexample to
RH. It proves only that the finite jets and these symmetries alone
do not force the desired zero location. Keeping the actual arithmetic
in (11) is essential.

## 7. Proof and execution record

| Statement | Status |
|---|---|
| Actual source decomposition (2), derivative recurrence (3)-(6) | Proved here by differentiation and absolute convergence |
| Primewise matrix determinant (7) and rational-component failure (8) | Proved exactly here |
| Full arithmetic identity (9)-(10) | Proved here |
| Universal arithmetic inequality (11) | **Open** |
| (11) implies RH | Conditional implication, using retained ILC-T2 |
| Finite-jet reflection alone forces RH | Refuted as an inference by Section 6 |
| New actual-source derivative values or larger positive bank | None computed |
| RH proof completion or global eta improvement | **Not established** |

This is a written derivation added to the existing repository, not a
new numerical certificate or replay stage. No scientific computation
has been executed locally or on GitHub for this note. Previous engines,
source data, certificate files and workflow definitions are unchanged.
Existing passing CI applies to the previously tested source, and is not
presented as machine verification of the new written arguments.
No independent agent was used.

Classical source formulas, compact moment theory and the identity
theorem retain their existing attribution. No new priority claim is
made for those tools or for the polynomial differentiation identity.
The [Clay Mathematics Institute status page](https://www.claymath.org/millennium/riemann-hypothesis/)
continues to list RH as unsolved (checked 17 October 2026).
