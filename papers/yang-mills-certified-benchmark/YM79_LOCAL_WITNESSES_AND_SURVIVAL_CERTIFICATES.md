# YM-79: local witnesses, nonlinear phase rescue and survival certificates

Monty Dabas. 5 October 2026. **Python 3.12 only.**

This stage improves the actual YM78 reflected-norm lower certificate
from order \(a^4\) to order \(a^2\), for the same source and Hamiltonian.
A smooth local witness supplies the better energy-moment estimate.
Bounded phase observables are then shown to rescue nonzero reflected
states when the original linear field loses energy at the infrared
endpoint. An explicit rational witness also certifies a reflected
local Wilson state of norm squared greater than \(13/100\) in the
actual untruncated side-four model at \(\theta=1/48\).

**This does not prove a positive continuum mass gap.** The \(a^2\)
certificate still vanishes. Nonlinear rescue requires excluding
complete ultraviolet escape for a fixed source. A pure infrared
limit, in fact, produces a gapless enlarged phase algebra.
The side-four certificate is at a fixed regulator.

## 1. Fixed model and the witness principle

Keep the unchanged periodic SU(2) Hamiltonian and actual vacuum:
\[
 A_\theta=L_E-\theta V,\quad V=\sum_bW_b,\quad
 H=A_\theta-E_0,\quad \mathsf H_a=\kappa_aH,\quad
 \mu_a=h_a^2d\Phi,\quad \kappa_a\sqrt{\theta_a}=c/a .
 \tag{1}
\]
The bare Laplacian is \(L_E=-\Delta\). The constant \(c>0\) retains
the original physical calibration. Take even \(L\ge4\), \(a\to0\),
and \(aL\to\infty\).

Use YM78's selected \(xy\) plaquettes based at \((2i,2j,k)\), with
disjoint four-link supports. Their spatial blocks have sides
\(2a,2a,a\). For real \(f\in C_c^\infty(\mathbb R^3)\), let
\(\bar f_b\) be its block average and \(P_af\) its block-constant
orthogonal projection. The original centered source is
\[
 F_a=\Psi_a(f)=\sum_bd_b(r_b-\mu_ar_b),\quad
 d_b=2\theta_a^{1/4}a^2\bar f_b,\quad r_b=\arccos W_b,\quad
 v_a=F_ah_a .
 \tag{2}
\]
YM78 proves the deterministic identity
\[
 \Gamma_{\rm phys}(F_a)=q_a
 =\langle v_a,\mathsf H_av_a\rangle
 =c\|P_af\|_2^2\longrightarrow q_0=c\|f\|_2^2 .
 \tag{3}
\]
For nonzero fixed \(f\), \(q_a>0\) at sufficiently small spacing.
The source belongs to \(D(\mathsf H_a)\).

For nonnegative self-adjoint \(K\), define
\[
 D_T(v)=\|(1-e^{-TK})v\|^2,\qquad
 R_T(v)=\|TK(I+TK)^{-1}v\|^2 .
 \tag{4}
\]
For \(v\in D(K^{1/2})\), YM78 gives
\[
 D_T(v)\ge R_T(v)=\sup_{z\in D(K)}
 \{2T\operatorname{Re}\mathcal E_K(v,z)-\|(I+TK)z\|^2\}.
 \tag{5}
\]
Here \(\mathcal E_K\) is the closed form of \(K\).
If one actual witness \(z\) satisfies
\[
 \operatorname{Re}\mathcal E_K(v,z)\ge b>0,\quad
 \|z\|^2\le N,\quad \mathcal E_K(z,z)\le Q,\quad \|Kz\|^2\le S,
\]
then scaling it by \(Tb/(N+2TQ+T^2S)\) proves
\[
 \boxed{D_T(v)\ge R_T(v)\ge\frac{T^2b^2}{N+2TQ+T^2S}.}
 \tag{6}
\]
This optimizes an admissible direction, not an assumed exact solve.

## 2. A smooth local witness improves the cutoff power

In this section \(\theta_a\ge25\). Set
\[
 \eta=\frac1{16},\quad \epsilon=\eta\theta_a^{-1/4},\quad
 \rho_\epsilon(r)=
 \arccos\!\left(\frac{\cos r}{\sqrt{1+\epsilon^2}}\right),
\]
\[
 G_a=\sum_bd_b[\rho_\epsilon(r_b)-\mu_a\rho_\epsilon(r_b)],
 \quad z_a=G_ah_a,\quad K_f=\frac{9\pi^2}{32}\|f\|_1^2 .
 \tag{7}
\]
As a function of \(W_b\in[-1,1]\), the regularizer is analytic on a
neighborhood of the interval. Thus \(G_a\) is smooth at both radial
poles, and \(z_a\) belongs to all finite-regulator operator domains.

### Theorem 79.1 — an actual \(a^2\) certificate

Uniformly in volume and \(\theta_a\ge25\),
\[
 \mathcal E_{\mathsf H_a}(v_a,z_a)\ge\tfrac9{10}q_a,\quad
 \mathcal E_{\mathsf H_a}(z_a,z_a)\le q_a,\quad
 \|z_a\|^2\le K_fa^{-2},
 \tag{8}
\]
\[
 \boxed{\|\mathsf H_az_a\|^2
       \le2^{18}c^2a^{-1}\|P_af\|_2^2.}
 \tag{9}
\]
The improved graph bound belongs to the smooth witness, not the raw
angle source. With
\[
 M_a=K_fa^{-2}+2Tq_a+2^{18}T^2c^2a^{-1}\|P_af\|_2^2,
\]
the explicit local witness \((9Tq_a/10M_a)z_a\) gives
\[
 \boxed{D_T(v_a)\ge R_T(v_a)\ge
 B_a(T,f):=\frac{81T^2q_a^2}{100M_a}.}
 \tag{10}
\]
For fixed \(T>0,f\ne0\),
\[
 B_a(T,f)\sim
 \frac{81T^2c^2\|f\|_2^4}{100K_f}\,a^2 .
 \tag{11}
\]
This improves YM78's asymptotic power. Conservative constants mean
the new numerical value need not exceed the old one at every
spacing; taking their maximum is valid. Equation (11) describes
the certificate, not the actual norm's asymptotic behavior.

### 2.1. Mixed energy survives smoothing

YM76's radial operator \(J=\nabla r\cdot\nabla+\cot r\) satisfies
\(\|Jh_a\|^2\le4\tau+1\) and
\(\|\nabla r\cdot\nabla h_a\|^2\le4\tau\).
Hence
\[
 \mu_a\cot^2r\le16\tau+2,\quad
 \mu_a\csc^2r\le16\tau+3
       \le\frac{183}{5}\sqrt{\theta_a}.
 \tag{12}
\]
We used the actual per-link estimate \(\tau\le(9/4)\sqrt{\theta_a}\).
Direct differentiation gives
\[
 \rho_\epsilon'=\frac{\sin r}{\sqrt{\sin^2r+\epsilon^2}}\in[0,1],
 \quad |\rho_\epsilon''|\le\epsilon^{-1},\quad
 1-\rho_\epsilon'\le\tfrac12\epsilon^2\csc^2r .
\]
It follows that
\[
 \mu_a\rho_\epsilon'\ge1-\frac{183}{2560}>\frac9{10}.
 \tag{13}
\]
Disjoint links turn the mixed form into a diagonal sum
\(\kappa_a\sum_bd_b^2\mu_a\rho_\epsilon'\), proving its lower bound.
The derivative bound proves the energy upper bound in (8).
Lipschitz contraction gives
\(\operatorname{Var}_\mu(\rho_\epsilon(r))\le\operatorname{Var}_\mu(r)\);
combine this with YM78's \(L^1\) sum estimate to get the norm bound.
No plaquette independence is assumed.

The \(9/10\) lower bound is for \(\mu\rho_\epsilon'\), not for
\(\mu(\rho_\epsilon')^2\).

### 2.2. An exact local third-moment identity

For smooth real \(F\) on the compact link manifold, write
\[
 X=\nabla F,\quad S=\operatorname{Hess}F,\quad g=\Delta F,\quad
 U=-\theta V-E_0,\quad H=-\Delta+U.
\]
Then
\[
 \begin{aligned}
 \langle Fh,H^3Fh\rangle={}&
 \int(8S^2-4\nabla_XS)(\nabla h,\nabla h)\,d\Phi\\
 &+\mu[X(\Delta g)-2\operatorname{div}(S\nabla g)]
       +2\mu X^2U .
 \end{aligned}
 \tag{14}
\]
Here \(S^2\) denotes the square of the Hessian endomorphism and
\(X^2U=X(XU)\).

To prove (14), let \(J_F=X\cdot\nabla+g/2\). It is skew on smooth
functions, \(H(Fh)=-2J_Fh\), and
\[
 [-\Delta,J_F]=-2\operatorname{div}(S\nabla)-\tfrac12\Delta g .
 \tag{15}
\]
This Riemannian divergence identity includes the Ricci terms.
Using \(Hh=0\), the third moment is
\(-2\langle h,[J_F,[H,J_F]]h\rangle\).
Expand
\[
 \nabla(J_Fh)=S\nabla h+\operatorname{Hess}(h)X
                       +(g/2)\nabla h+(h/2)\nabla g .
\]
Integrating \(X[S(\nabla h,\nabla h)]\) cancels the \(gS\) term.
The mixed term integrates to
\(-2\mu\operatorname{div}(S\nabla g)\); the scalar term leaves
\(\mu X(\Delta g)\). Finally, \([U,J_F]=-XU\) supplies \(2\mu X^2U\).
This proves (14) without boundary or pole terms.

For \(F=\sum_bd_b\rho_\epsilon(r_b)\), each Hessian block is supported
on its own four links. Its derivative along \(X\) involves only that
block. Thus the kinetic terms have no cross-block contributions.

### 2.3. Explicit constants and physical rescaling

A four-link block has dimension 12. Ordered native invariant-frame
derivatives of \(W\) of order \(k\) have absolute value at most
\(2^{-k}\). Connection coefficients are constant, zero or
\(\pm1/2\), with at most one nonzero output index. A component of
\(\nabla^kW\) is a sum of ordered frame derivatives. Weighting a
derivative of order \(r\) by \(2^{-r}\), the absolute coefficient
sum obeys \(B_1=1/2\), \(B_{k+1}\le(k+1)B_k/2\).
Consequently, for \(k\ge2\),
\[
 \|\nabla^kW\|_{\rm HS}\le(\sqrt3)^kk!\le M_k:=2^kk!,
 \qquad \|\nabla W\|=\sin r .
 \tag{16}
\]
The inherited order of tensor indices is retained; higher
covariant tensors are not assumed symmetric.

Put \(D=(\epsilon^2+\sin^2r)^{1/2}\) and
\(f_\epsilon(w)=\arccos(w/\sqrt{1+\epsilon^2})\). For \(1\le j\le5\),
\[
 |f_\epsilon^{(j)}(w)|\le C_jD^{1-2j},\quad
 (C_1,C_2,C_3,C_4,C_5)=(1,1,4,24,204).
\]
For example, the fifth-derivative numerator is
\(9(1+\epsilon^2)^2+72(1+\epsilon^2)w^2+24w^4\).
The bound uses \(\epsilon\le1\).
For a chain-rule partition with \(j\) blocks and \(n_1\) singleton
blocks, the singleton gradients supply
\(\sin^{n_1}r\le D^{n_1}\), and \(2j-n_1\le k\).
The Bell-partition sums give
\[
 \|\nabla^k\rho_\epsilon(r)\|_{\rm HS}\le A_k\epsilon^{1-k},
 \quad (A_2,A_3,A_4,A_5)=(9,76,984,17484).
 \tag{17}
\]
For example \(A_3=M_3+3M_2+4=76\).
All four sums are replayed exactly by the certificate.

Set
\[
 K=8A_2^2+4A_3=952,\qquad
 Z=12A_5+24A_3^2+4\sqrt3 A_2A_4\le419280 .
\]
The first line of (14) is bounded above by
\(4\tau K\epsilon^{-2}\sum_bd_b^2\).
Tracing the derivative tensors in the second line gives
\(Z\epsilon^{-4}\sum_bd_b^2\).
Here the trace factors 12 and \(\sqrt{12}\) are the four-link
block's dimension factors.

For the potential, native holonomy differentiation gives
\(\|\operatorname{Hess}W_b\|_{\rm op}\le1\).
Each link belongs to four plaquettes, so
\(\|\operatorname{Hess}V\|_{\rm op}\le4\).
Also, for each selected block,
\[
 \mu|\nabla_{\rm block}V|^2
 \le16\mu(1-W_b^2)\le72\theta^{-1/2}.
 \tag{18}
\]
This uses incidence, Cauchy–Schwarz and YM76's actual vacuum bound.
Since \(\Gamma(r_b)=1\), one has
\(\nabla_{\nabla r_b}\nabla r_b=0\) away from the poles.
The acceleration of \(\nabla\rho_\epsilon(r_b)\) is therefore
\(\rho_\epsilon'\rho_\epsilon''\nabla r_b\).
The potential contribution in (14) is at most
\[
 (8+12\sqrt2/\eta)\,\theta\sum_bd_b^2 .
\]
Combining the bounds with \(\epsilon=\eta\theta^{-1/4}\) gives
\[
 \langle G_ah_a,H^3G_ah_a\rangle
 \le C_\eta\theta\sum_bd_b^2,\quad
 C_\eta\le9K\eta^{-2}+419280\eta^{-4}+8+24/\eta
 =27\,480\,127\,880<2^{36}.
 \tag{19}
\]
The rational estimate used \(\sqrt3,\sqrt2\le2\).
Multiplying by \(\kappa_a^3\) and using (1)–(3) yields
\[
 m_3(z_a)\le2^{36}c^3a^{-2}\|P_af\|_2^2 .
\]
Spectral Cauchy–Schwarz, \(m_2^2\le m_1m_3\), with
\(m_1(z_a)\le q_a\), proves (9).
Equation (6) now proves Theorem 79.1. No clustering or uniform
spectral gap was assumed.

## 3. Bounded phases rescue infrared loss

Let \(X_a(t)\) be the actual stationary reversible vacuum diffusion.
For fixed \(f\ne0,T>0\), set
\[
 Y_a(T)=F_a(X_a(T))-F_a(X_a(0)),\quad
 V_a(T)=\mathbb E|Y_a(T)|^2 .
\]
For real \(\lambda\), define the bounded local phase observable
\[
 U_{a,\lambda}=e^{i\lambda F_a},\quad
 w_{a,\lambda}=(U_{a,\lambda}-\mu_aU_{a,\lambda})h_a .
\]
It is gauge invariant and center neutral, with the same support as
\(F_a\), and
\[
 \|w_{a,\lambda}\|^2=1-|\mu_aU_{a,\lambda}|^2\le1,\qquad
 \mathcal E_{\mathsf H_a}(w_{a,\lambda},w_{a,\lambda})
       =\lambda^2q_a .
 \tag{20}
\]
Write
\[
 I_{a,\lambda}(T)=
 \mathbb E|U_{a,\lambda}(X_a(T))-U_{a,\lambda}(X_a(0))|^2,\qquad
 D^W_{a,\lambda}(T)=D_T(w_{a,\lambda}).
\]

### Theorem 79.2 — nonlinear conversion of increment survival

At every regulator,
\[
 \mathbb E|Y_a(T)|^4\le12q_a^2T^2,\quad
 I_{a,\lambda}(T)\ge\lambda^2V_a(T)-\lambda^4q_a^2T^2,\quad
 \boxed{D^W_{a,\lambda}(T)\ge I_{a,\lambda}(T)^2/4.}
 \tag{21}
\]
Thus, if \(V_a(T)\ge d>0\) along a family and \(q_a\le\bar q\),
the fixed amplitude \(\lambda^2=d/(2\bar q^2T^2)\) gives
\[
 \boxed{D^W_{a,\lambda}(T)\ge d^4/(64\bar q^4T^4)>0.}
 \tag{22}
\]
This stage does not establish the premise \(V_a(T)\ge d\) on a
continuum trajectory.

**Proof.** The forward-backward decomposition is
\(Y_a(T)=(M_+-M_-)/2\). The deterministic carré du champ (3)
makes each martingale increment Gaussian with variance \(2q_aT\).
They need not be independent. Convexity gives
\[
 \mathbb E|(M_+-M_-)/2|^4
 \le\tfrac12(\mathbb E|M_+|^4+\mathbb E|M_-|^4)=12q_a^2T^2.
\]
The YM76 radial domain argument justifies this decomposition.
The global inequality \(2(1-\cos x)\ge x^2-x^4/12\) then gives
the lower bound on \(I\).
For the spectral measure \(\nu^W_a\) of the centered phase source,
\[
 I=2\int(1-e^{-TE})d\nu^W_a,\qquad
 D^W=\int(1-e^{-TE})^2d\nu^W_a .
\]
Its total mass is at most one. Cauchy–Schwarz proves \(D^W\ge I^2/4\).
The stated choice of \(\lambda\) proves (22). \(\square\)

### Exactly which endpoint is removed

Let \(\eta_a(dE)=E\|P_{\mathsf H_a}(dE)v_a\|^2\).
Along a subsequence with
\(\eta_a\Rightarrow\bar\eta\) on \([0,\infty]\),
\[
 V_a(T)=2\int\frac{1-e^{-TE}}E\,d\eta_a(E)
 \longrightarrow2\int\frac{1-e^{-TE}}E\,d\bar\eta(E).
 \tag{23}
\]
The kernel equals \(T\) at zero, is positive at every finite energy,
and vanishes at infinity. Hence:

- If \(\bar\eta([0,\infty))>0\), (22) supplies a nonzero reflected
  phase state for some fixed \(\lambda\).
- If \(\bar\eta=q_0\delta_\infty\), then \(V_a(T)\to0\).
  Since \(I_{a,\lambda}\le\lambda^2V_a(T)\) and \(D^W\le I/2\),
  every fixed-amplitude filtered phase state collapses.

This equivalence is for the selected source family along that
subsequence. A subsequential result is not a full-trajectory lower
bound. The obstruction to this rescue criterion is **complete
ultraviolet energy escape**, rather than infrared escape as well.

For a countable additive group of phase tests, circle-valuedness
and the uniform time-increment moments give joint subsequential
stationary path laws by diagonal Kolmogorov tightness.
Multiplication identities \(U(f+g)=U(f)U(g)\) and reflection positivity
pass to these laws. Phases use ordinary even time reflection, as
functions of the undifferentiated diffusion configuration.
Continuity in test parameters at time zero, spatial field regularity
of this algebra, and full Euclidean covariance are not inferred.

### Infrared rescue does not provide a mass gap

For a stationary Ornstein–Uhlenbeck coordinate with covariance
\((q/E)e^{-E|t|}\),
\[
 C_{E,\lambda}(t)=
 \mathbb E[\overline{U_\lambda(0)}U_\lambda(t)]
 =\exp[-(\lambda^2q/E)(1-e^{-E|t|})].
 \tag{24}
\]
Its spectral law is \(E\) times a Poisson variable of mean
\(\lambda^2q/E\), including the vacuum atom.
As \(E\downarrow0\), it tends to \(\delta_{\lambda^2q}\), and
\[
 D^W_{E,\lambda}(T)\to(1-e^{-\lambda^2qT})^2>0
 \quad(\lambda\ne0).
 \tag{25}
\]
The raw linear reflected norm tends to zero in this same example.
For \(E\to\infty\), the centered phase norm and dynamics collapse.
These are logical controls, not interacting Yang–Mills vacua.

An actual-model conditional statement is also available. If
\(\eta_a\Rightarrow q_0\delta_0\), write the forward decomposition
\(Y_a=M_a+A_a\). Its martingale is Brownian with bracket \(2q_at\),
and stationarity gives
\[
 \mathbb E|A_a(t)|^2=
 2\int\left[t-\frac{1-e^{-Et}}E\right]d\eta_a(E)\longrightarrow0 .
 \tag{26}
\]
The integrand is bounded on the compactified half-line and vanishes
at zero. Finite-dimensional convergence and increment tightness
yield Brownian increments. The deterministic-bracket characterization
makes those increments independent of the initial configuration.
Stationarity forces the limiting circle phase to have Haar law.
Each nonzero amplitude therefore gives circle Brownian motion with
rate \(\lambda^2q_0\).

If the limiting algebra includes all rational amplitudes, its phase
states have spectral energies \(\lambda^2q_0\); choosing \(\lambda=1/n\)
forces accumulation at zero. **The enlarged phase algebra is gapless
in this pure infrared regime.** Its time-zero character is nonregular
at \(\lambda=0\). Selecting a single circle subalgebra must not be
mistaken for constructing a complete Yang–Mills theory.

## 4. A computable witness when the vacuum is only enclosed

### Theorem 79.3 — ground-state error and resolvent residual

Let \(A\) have normalized simple ground state \(h\), ground energy
\(E_0\), and next energy at least \(\Lambda\). Suppose
\[
 L\le E_0\le U,\quad \|g\|=1,\quad
 L\le\rho=\langle g,Ag\rangle<\Lambda .
\]
With the phase of \(g\) aligned with \(h\),
\[
 1-|\langle h,g\rangle|^2\le s=\frac{\rho-L}{\Lambda-L},\qquad
 \|h-g\|^2\le2(1-\sqrt{1-s}).
 \tag{27}
\]
Indeed \(\rho\ge E_0+(\Lambda-E_0)(1-|\langle h,g\rangle|^2)\),
and the resulting ratio decreases with \(E_0\).
The trial's phase alignment must be established.

Let \(O\) be bounded, \(K=A-L\ge0\), \(H=A-E_0\ge0\), \(w=Og\).
For any \(y\in D(A)\), compute
\[
 r=w-(I+TK)y,\quad
 e=\|O\|\epsilon+\frac{T(U-L)}{1+T(U-L)}\|w\|+\|r\|,
 \qquad \|h-g\|\le\epsilon .
\]
Then
\[
 \boxed{D_T((O-\langle h,Oh\rangle)h)
 \ge R_T(Oh)\ge(\|w-y\|-e)_+^2.}
 \tag{28}
\]
This also applies to a non-self-adjoint \(O\), with complex centering.

**Proof.** Put \(B_T(H)=I-(I+TH)^{-1}\).
Since \(K=H+(E_0-L)\), spectral calculus gives
\[
 \|B_T(K)-B_T(H)\|\le T(U-L)/(1+T(U-L)).
\]
Also \(w-y=B_T(K)w+(I+TK)^{-1}r\).
Subtract it from \(B_T(H)Oh\), use the operator norm bounds by one,
and apply the triangle inequality. Centering drops out because
\(B_T(H)h=0\). \(\square\)

For gauge-invariant finite electric polynomials, the residual
arithmetic is finite and native. YM74 supplies convergent
finite-regulator energy enclosures. The electric core approximates
the vacuum in form norm and a resolvent solve in graph norm.
Thus the method is arbitrarily refinable at a fixed regulator,
with phase control. It supplies no uniform complexity or nonzero
limiting value as the regulator is removed.

## 5. An actual local Wilson witness exceeding \(13/100\)

### Theorem 79.4 — a rational untruncated-model result

Take side four, \(p=192\), \(\theta=1/48\), and one elementary
plaquette \(W\). In bare time \(T=1\), for \(v=(W-\mu W)h\),
\[
 \boxed{\|(1-e^{-H})v\|^2\ge R_1(v)
       \ge\frac{32607}{250000}>\frac{13}{100}.}
 \tag{29}
\]
Physical time is \(1/\kappa\). This is a fixed-regulator statement.

**Proof.** Normalized Haar expectation on the actual graph gives
\[
 \langle V^2\rangle_0=48,\quad \langle V^3\rangle_0=0,\quad
 \langle W^2\rangle_0=\tfrac14,\quad
 \langle W^4\rangle_0=\tfrac18,\quad
 \langle W^6\rangle_0=\tfrac5{64},
\]
\[
 \langle VW^2\rangle_0=\langle VW^4\rangle_0=0,\quad
 \langle V^2W^2\rangle_0=\tfrac{193}{16},\quad
 \langle V^2W^4\rangle_0=\tfrac{387}{64}.
 \tag{30}
\]
Odd link multiplicity kills the relevant cross terms by a Haar
center flip. A different plaquette has a private edge, giving
conditional second moment \(1/4\). The remaining even moments are
single SU(2) Haar moments. Native graph enumeration verifies these
parity and private-edge premises. No vacuum factorization is used.

The normalized trial is
\[
 g=\frac{1+V/144}{\sqrt{433/432}},\qquad
 \rho=\langle g,A_\theta g\rangle=-3/433 .
 \tag{31}
\]
The unchanged YM75/78 enclosures give
\[
 E_0\ge L=-712209/10^8,\quad E_1\ge\Lambda=48376192/10^8,\quad
 E_0\le\rho<0 .
\]
For this trial,
\[
 \frac{\rho-L}{\Lambda-L}<\frac{9999}{25000000}
       =1-(1-1/5000)^2 ,
\]
so (27) gives \(\|h-g\|\le1/50\).
The phase is known because \(A_\theta\Omega=-\theta V\), \(\Omega=1\),
implies
\[
 \langle h,1+V/144\rangle
       =(1-E_0/(144\theta))\langle h,\Omega\rangle>0 .
 \tag{32}
\]
The trial itself need not be pointwise positive.

Let \(a_0=3-E_0\), so \(3\le a_0\le301/100\). Since \(L_EW=3W\),
\[
 HW=(a_0-\theta V)W,\qquad
 \mathcal E_H(v,W)=\langle h,(a_0-\theta V)W^2\rangle .
\]
The exact moments give
\[
 \langle g,(a_0-\theta V)W^2\rangle
   =\frac{a_0/4-193/110592}{\sqrt{433/432}},\qquad
 \|(a_0-\theta V)W^2\|^2=a_0^2/8+\theta^2\,387/64 .
\]
Thus
\[
 \|(a_0-\theta V)W^2\|\le107/100,\quad
 \sqrt{433/432}\le501/500,
\]
\[
 \mathcal E_H(v,W)\ge
 (3/4-193/110592)\frac{500}{501}-\frac{107}{5000}\ge29/40,
\]
\[
 \|(I+H)W\|^2=(1+a_0)^2/4+\theta^2\,193/16\le403/100 .
 \tag{33}
\]
Insert the actual local witness \(z=(9/50)W\) in (5):
\[
 2\frac9{50}\frac{29}{40}
 -\left(\frac9{50}\right)^2\frac{403}{100}
       =\frac{32607}{250000}.
 \tag{34}
\]
This proves (29). The full neutral Hilbert space is untruncated;
finite Haar arithmetic evaluates the trial and witness, while
inherited spectral bounds control the actual vacuum. \(\square\)

## 6. Scope and reproducibility

The new weak-coupling witness improves the original source's
certificate to \(a^2\). The explicit side-four calculation also
turns a resolvent-survival condition into a positive actual-model
certificate at one regulator.

The phase theorem changes the nontriviality question: surviving
ordinary angle increments can yield nonzero reflected phase states,
even if the linear reflected state disappears. For this source
family, excluding complete ultraviolet energy escape suffices for
subsequential survival. Pure infrared survival nevertheless gives
a gapless enlarged phase algebra.

The immediate missing uniform estimate is one nonvanishing local
witness value, or one fixed \(f,T\) with
\(\liminf_aV_a(T)>0\). The \(a^2\) certificate does not establish either.
Full Euclidean covariance, regular local Yang–Mills field
identification, the required short-distance behavior, and a uniform
positive local gap remain unproved.

The exact certificate replays derivative-partition constants,
third-moment budgets, physical scaling, cutoff comparisons, nonlinear
conversion bounds, endpoint examples, native graph moments and the
rational Wilson witness. It pins 288 upstream inputs and four local
proof/code/test/workflow inputs. The shared engine and upstream
sources are unchanged.

~~~sh
python3.12 -B papers/yang-mills-certified-benchmark/certificates/ym79_local_witness_survival_certificates.py --check
python3.12 -B -m unittest discover -s papers/yang-mills-certified-benchmark/tests -p test_ym79_local_witness_survival_certificates.py -v
~~~

The dedicated workflow also replays YM78–YM75. Finite checks do not
prove the open continuum assertions.

## References

1. A. Jaffe and E. Witten, *Quantum Yang–Mills Theory*, official Clay
   Mathematics Institute problem statement, §4.
   [Primary text](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf).
   The target includes nontrivial local fields, the prescribed
   axiomatic and short-distance properties, and a positive gap.
2. L. Wu, *Forward-backward martingale decomposition and compactness
   results for additive functionals of stationary ergodic Markov
   processes*, Annales de l'Institut Henri Poincaré, Probabilités
   et Statistiques **35** (1999), 121–141.
   [Primary archive](https://www.numdam.org/item/AIHPB_1999__35_2_121_0/).
   The programme derives its specific decomposition and bounds in
   YM77–YM78.
3. [YM74](YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md),
   [YM75](YM75_RH_SOURCE_RESOLVED_SPECTRAL_TRANSFER.md),
   [YM76](YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS.md),
   [YM77](YM77_PHYSICAL_SMEARING_AND_CONTINUUM_SURVIVAL.md),
   [YM78](YM78_STATIONARY_DERIVATIVE_AND_ENERGY_SURVIVAL.md).
   Their Hamiltonian, domains, source identities and certificate
   hypotheses are retained.
