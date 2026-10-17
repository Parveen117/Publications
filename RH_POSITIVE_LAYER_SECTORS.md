# RH continuation: positive layer sectors at unbounded depth

17 October 2026. Continues [the arithmetic gate](RH_ARITHMETIC_GATE.md)
at Publications commit `f4aad12c93c204a50bc2f4d4823667c86dbc3ed9`.

The full Gamma, rational and prime expression is retained through the
unconditional zero product of the same actual xi source. This yields
positive regions containing infinitely many layer indices. No global RH
hypothesis is used, but one published finite-height verification theorem
is an explicit input. The uncovered indices still prevent an RH conclusion.

## Inputs and notation

Retain
\[
\mathcal X(w^2)=\xi(1/2+w),\quad q=\mathcal X'/\mathcal X,\quad
a=(57/2)^2=3249/4,
\]
\[
c_n=\frac{(-a)^nq^{(n)}(a)}{n!\,q(a)},\qquad
D_{r,k}=\sum_{j=0}^{k}(-1)^j\binom{k}{j}c_{r+j}.
\tag{1}
\]
All indices are nonnegative integers. The following inputs are reused.

1. The standard entire order-one xi function, functional equation,
   conjugation symmetry, absence of real xi zeros, and the nontrivial
   zero strip \(0<\beta<1\). These are the analytic source already
   used in SE39/SE48; the canonical product is also recorded in
   [NIST DLMF 25.2.12](https://dlmf.nist.gov/25.2#E12).
2. Platt and Trudgian, *The Riemann hypothesis is true up to*
   \(3\cdot10^{12}\), **Bulletin of the London Mathematical Society**
   53 (2021), 792–797,
   [DOI 10.1112/blms.12460](https://doi.org/10.1112/blms.12460),
   [arXiv:2004.09765, Theorem 1](https://arxiv.org/pdf/2004.09765).
   We use only the conservative statement that every zero with
   \(0<\gamma\le T:=3\cdot10^{12}\) has \(\beta=1/2\).
   This theorem is cited, not independently rerun or replaced by an
   assumed global RH.
3. The retained [SE48 local zero certificate](papers/native-prime-observation-closure/sources/SE48.md),
   Section 5, supplies an actual simple critical-line zero with
   \(32<\gamma_\star<33\): its certified disk is centered at
   \(1/2+32.935061722570\,i\), with radius \(10^{-6}\).
   Only the coarse inequalities \(32<\gamma_\star<33\) are used.
4. [SE49](papers/native-moving-zero-covariance/sources/SE49.md), Section 2,
   proves the unconditional majorant
   \(\lvert \xi'(s)/\xi(s)\rvert<10+|s|/2\) for \(\Re s>3/2\).
   At \(s=29\) it gives
   \[
   0<q(a)<\frac{49}{114}<\frac12.
   \tag{2}
   \]
   Positivity is retained from the existing source certificate and
   also follows from the zero sum below.

The use of input 2 is essential to the numerical regions stated here.
These results are consequences of published mathematics plus retained
source results, not a new independent finite-height zero verification.

## Theorem S1. Exact combined source and a high-zero remainder budget

Index the nontrivial zeros with positive imaginary part, including
multiplicity. For \(\rho=\beta+i\gamma\), put
\[
\delta_\rho=\beta-\tfrac12,\quad
\mu_\rho=(\gamma-i\delta_\rho)^2
=\gamma^2-\delta_\rho^2-2i\delta_\rho\gamma,\quad
y_\rho=\frac{a}{a+\mu_\rho}.
\tag{3}
\]
The multiset of these \(\mu_\rho\) is closed under conjugation; a
critical-line zero contributes a positive real \(\mu_\rho=\gamma^2\).
Then
\[
q(t)=\sum_{\gamma>0}\frac1{t+\mu_\rho},
\tag{4}
\]
and for every finite \(r,k\),
\[
\boxed{
q(a)D_{r,k}
=\sum_{\gamma>0}\frac{y_\rho^r(1-y_\rho)^k}{a+\mu_\rho}
=\frac1a\sum_{\gamma>0}y_\rho^{r+1}(1-y_\rho)^k.
}
\tag{5}
\]
The sums are absolutely convergent. Nonreal \(\mu_\rho\) are retained:
(4) is not an assumed positive Stieltjes representation.

Split (5) at \(T\). Define
\[
L_{r,k}(T)=\sum_{0<\gamma\le T}
  \frac1{a+\gamma^2}
  \left(\frac{a}{a+\gamma^2}\right)^r
  \left(\frac{\gamma^2}{a+\gamma^2}\right)^k,\quad
u_T=\frac{a}{a+T^2-1/4}.
\tag{6}
\]
The remainder \(\mathcal R_{r,k}(T)=q(a)D_{r,k}-L_{r,k}(T)\)
is real and satisfies
\[
\boxed{
|\mathcal R_{r,k}(T)|\le2q(a)\,u_T^r<10^{-22r}.
}
\tag{7}
\]
For \(r=0\), the last display means the strict bound \(<1\).

**Proof of (4)-(5).** The even entire function
\(w\mapsto\xi(1/2+w)\) has order one; hence \(\mathcal X(t)\)
has order \(1/2\). Its genus-zero Hadamard product is
\[
\mathcal X(t)=\mathcal X(0)\prod_{\gamma>0}(1+t/\mu_\rho).
\]
There is no nonconstant exponential factor at this order, and
\(\sum_{\gamma>0}|\mu_\rho|^{-1}<\infty\). The upper-half-plane
index selects exactly one representative of each pair \(w,-w\);
an off-line quartet gives the two conjugate \(t\)-zeros. Logarithmic
differentiation gives (4), with locally uniform convergence away from
poles. Taking any fixed number of derivatives at \(a\) and then the
finite binomial sum gives (5).

This is the same q as in the Gamma/prime identity of the preceding
note: both follow from the actual completed xi function. No term of
that arithmetic source has been discarded.

**Proof of (7).** At low zeros \(\mu_\rho=\gamma^2>0\).
At high zeros \(\gamma>T\) and \(|\delta_\rho|<1/2\), so
\(\Re\mu_\rho>\gamma^2-1/4>0\). Thus every term of (4) at \(a\)
has strictly positive real part. For a high zero,
\[
|y_\rho|<u_T,\qquad |1-y_\rho|<1.
\tag{8}
\]
The second inequality follows from
\(|a+\mu|^2-|\mu|^2=a^2+2a\Re\mu>0\).
Write \(d=a+\mu_\rho\). Since \(a>1/4\),
\[
\frac{|\Im d|}{\Re d}
\le\frac{\gamma}{a+\gamma^2-1/4}<\frac1\gamma<1.
\]
Consequently
\[
\left|\frac1d\right|
=\sqrt{1+\left(\frac{\Im d}{\Re d}\right)^2}
  \Re\frac1d
<2\Re\frac1d.
\]
Absolute convergence and the positivity of all real parts give
\[
\sum_{\gamma>T}\left|\frac1{a+\mu_\rho}\right|
\le2\sum_{\gamma>T}\Re\frac1{a+\mu_\rho}
\le2q(a).
\]
Combine this with (8). Finally \(a<900\), \(T^2=9\cdot10^{24}\),
and \(a>1/4\) imply \(u_T<10^{-22}\); (2) gives \(2q(a)<1\).
For \(r=0\), use the last inequality directly. QED.

The bound is uniform in the finite index \(k\). It does not neglect
the high zeros or assert that each high-zero summand is positive.

## Theorem S2. An unbounded positive sector

For the actual source at \(a=3249/4\),
\[
\boxed{
D_{r,k}>\frac23\,10^{-22r}>0
\qquad(r\ge1,\quad 0\le k\le60r).
}
\tag{9}
\]
In particular,
\[
\boxed{c_n>0\quad\hbox{for every }n\ge0.}
\tag{10}
\]
For each fixed \(k\), (9) proves positivity for all integers
\(r\ge\max(1,\lceil k/60\rceil)\).

**Proof.** The single retained zero \(\gamma_\star\in(32,33)\)
in input 3 has, using \(800<a<900\),
\[
\frac1{a+\gamma_\star^2}>\frac1{2000},\qquad
\frac{a}{a+\gamma_\star^2}>\frac25,\qquad
\frac{\gamma_\star^2}{a+\gamma_\star^2}>\frac12.
\tag{11}
\]
For the middle bound,
\(a/(a+\gamma_\star^2)>800/(800+1089)>2/5\).
For the last bound \(\gamma_\star^2>1024>a\).
All other low-zero terms are nonnegative. Therefore
\[
q(a)D_{r,k}>
\frac1{2000}\left(\frac25\right)^r2^{-k}-10^{-22r}.
\tag{12}
\]
When \(k\le60r\), the first term is at least
\[
\frac1{2000}\left(\frac25\,2^{-60}\right)^r.
\]
The exact integer comparison
\[
2^{60}=1152921504606846976
<1500000000000000000
\]
gives
\[
\frac25\,2^{-60}>\frac{8000}{3}\,10^{-22}.
\]
For every integer \(r\ge1\),
\[
\frac1{2000}\left(\frac{8000}{3}\right)^r-1
\ge\frac13.
\]
Substitution in (12) proves
\(q(a)D_{r,k}>\frac13\,10^{-22r}\).
Since \(q(a)<1/2\), (9) follows. For \(n\ge1\), take \(r=n,k=0\);
the remaining case is \(c_0=1\). QED.

The constants are deliberately coarse. Equation (12) and the exact
retained pole interval permit larger regions without changing the
argument. No optimal sector boundary is claimed.

## Theorem S3. A second positive region from phase control

For the same source,
\[
\boxed{
D_{r,k}>0
\quad\text{if}\quad
0\le r\le3\cdot10^{12}-1,\qquad 0\le k\le10^{34}.
}
\tag{13}
\]
This rectangle is finite, while (9) has no upper bound on depth.

**Proof.** Every low-zero summand in (5) is strictly positive.
For a high zero, conjugation allows us to consider \(\delta\ge0\).
Write
\[
u=\gamma^2-\delta^2>0,\quad v=2\delta\gamma\ge0,\quad
\theta=\arctan\frac{v}{a+u},\quad
\phi=\arctan\frac vu,\quad \tau=\phi-\theta.
\]
Then
\[
\arg y=\theta,\qquad \arg(1-y)=-\tau,\qquad
0\le\theta<\frac1T.
\tag{14}
\]
If \(v=0\), both angles vanish. Otherwise, the tangent subtraction
formula, with \(0<\tau<\pi/2\), gives
\[
\tau=\arctan\frac{av}{u(a+u)+v^2}
\le\frac{a\gamma}{(\gamma^2-1/4)(a+\gamma^2-1/4)}
<\frac{2a}{\gamma^3}<\frac{1800}{T^3}.
\tag{15}
\]
For the strict middle bound use \(a>1/4\), \(\gamma>T>1\),
and \(\gamma^2-1/4>\gamma^2/2\).

The phase of the last summand in (5) is
\((r+1)\theta-k\tau\). In the range (13),
\[
0\le(r+1)\theta<1,\qquad
0\le k\tau<
\frac{10^{34}\cdot1800}{27\cdot10^{36}}=\frac23.
\]
Its phase lies between \(-2/3\) and \(1\), strictly inside
\((-\pi/2,\pi/2)\). Its real part is therefore positive. The conjugate
summand has the same positive real part. Absolute convergence now
proves strict positivity of the real sum in (5), and (2) completes
the proof. QED.

The phase is bounded at each finite pair \(r,k\); it is not replaced
by zero. Thus possible off-line zeros above T are included in this
argument rather than excluded by assumption.

## Theorem S4. The remaining all-order condition

Let
\[
\mathcal P=
\{(r,k):r\ge1,\ 0\le k\le60r\}
\ \cup\
\{(r,k):0\le r<T,\ 0\le k\le10^{34}\},
\]
and let \(\mathcal U=\mathbb Z_{\ge0}^2\setminus\mathcal P\).
Then, for this unchanged actual source,
\[
\boxed{\mathrm{RH}\quad\Longleftrightarrow\quad
D_{r,k}\ge0\ \hbox{for every }(r,k)\in\mathcal U.}
\tag{16}
\]
**Proof.** ILC-T2 establishes equivalence with nonnegativity at every
pair. Theorems S2 and S3 settle every pair in \(\mathcal P\), so only
its complement remains. The converse restriction is immediate. QED.

In elementary inequalities, the unsettled complement is
\[
k>60r
\quad\hbox{and}\quad
\bigl(r\ge3\cdot10^{12}\ \hbox{or}\ k>10^{34}\bigr).
\tag{17}
\]
At \(r=0\) this formula has the same meaning, because \((0,0)\)
is already in the positive rectangle.

Neither (9) nor (13) proves all finite Hankel forms positive. The
Hausdorff criterion needs every \(r,k\); positivity of all individual
\(c_n=D_{n,0}\) is weaker. No new all-order positive contraction,
unitary xi dynamics, global eta reduction, or RH conclusion is claimed.

## Provenance and validation

| Item | Evidence |
|---|---|
| Same complete arithmetic source | Actual xi identity, unconditional product, equations (4)-(5) |
| Unknown high-zero contribution | Absolute bound (7), with all multiplicities retained |
| Unbounded positive sector | Written proof S2; uses one retained local zero and published finite-height theorem |
| Finite positive rectangle | Written phase proof S3 |
| Remaining RH implication | Existing ILC-T2 and the exact set complement in S4 |
| All signs on the remaining set | Open |
| New numerical experiment or zero verification | None |
| New machine-checked certificate | None |

This is a written theorem continuation within Publications. No local
scientific computation, independent agent, new scientific CI stage or
rerun of the Platt–Trudgian computation is used. The documentation commit
uses skip-CI to avoid replaying unchanged repository checks. Existing
source files, numerical certificates, native engines and workflows are
unchanged. The native finite certificates are cited as predecessor
evidence, not as machine verification of these new arguments.

The canonical product, finite-height verification, and analytic source
majorant retain their original attribution. The explicit consequences
in (9) and (13) are derived here; no literature-wide novelty claim is made.
