# PR1 — Sector speed, mass and curvature from one coin turn; each sector has its own relativity

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
First stage of the propagation line.

Owner's thesis for this line: what matters is not the physical quantity
but the curvature; after a cut the reading is flat; curvature propagates,
the speed belongs to the sector, and mass is the coupling of the response.

Sources read before building: **R28** (source roles H² = K² = I,
HK = −KH, R = KH, curvature [H,K] = −2R; walk U = DC with D = SP + S⁻¹Q),
**R29** (V = U²), **R38** (sharp cone c_Σ = 1/√2), **RMG1** (odd ratio in
the Klein disc), **RMG9-T1** (square law of uK + vS + τR), **PH2**,
**QC3-T1** (E² − O² = I).

## 1. One coin turn fixes speed, mass and curvature

Let the coin be the circular turn Exp(θR) = cos θ + sin θ·R and
U = D·Exp(θR): shift by the cut H, then turn.

**T1 (exact).**

```text
U + U⁻¹ = cos θ · (S + S⁻¹)                     every component: future + past = cos θ × (left + right)
cos ω = cos θ · cos k                           event turn ω, address turn k
1 − E_ω = (1 − E_θ) + E_θ (1 − E_k)             event defect = coin defect + (coin even part) × address defect
(dω/dk)² ≤ cos²θ ,  equality only at sin k = 1  cone speed of the sector = cos θ
at k = 0:  U = Exp(θR)                          the uniform sector only turns: rest turn θ
[H, Exp(θR)] = −2 sin θ · K                     curvature of the coin against the cut
```

**T2 (the law of the sector).**

```text
(cone speed)² + (curvature / 2)² = 1 .
```

The even part of the coin turn is the speed; its odd part is half the
curvature and is the rest turn (mass). This is QC3's E² − O² = I read as
propagation.

```text
θ = 0      coin = pure cut reading      curvature 0      speed 1      no rest turn     (flat, light-like)
θ = π/4    R28's coin                   (curv/2)² = ½    speed² = ½   — the value of R38's c_Σ
θ = π/2    coin = pure exchange         curvature 2      speed 0      nothing propagates
```

For R28's own coin (H + K)/√2 the two-event map satisfies
V + V⁻¹ = 2 − ½(2 − S² − S⁻²): a wave law with speed² = ½, certified
without the root. R38's theorem for the two-dimensional source is not
re-derived here; the number is the same.

So: a flat reading propagates at the full unit speed; curvature against
the cut is paid for in speed; what is not propagated is rest turn.

## 2. Response space is velocity space

**T3.** For responses H = m(1 + uK + vS) (RMG1):

```text
collinear:       (1 + u₁K)(1 + u₂K) has odd ratio (u₁ + u₂)/(1 + u₁u₂)
non-collinear:   product = rotation × response,  tan(rotation) = (w₁ × w₂)/(1 + w₁·w₂)
rapidity:        cosh ℓ₃ = [ (1+ρ₁²)(1+ρ₂²) + 4 w₁·w₂ ] / [ (1−ρ₁²)(1−ρ₂²) ]
```

The odd ratio composes as a velocity in units of the sector speed;
positivity u² + v² < 1 is the speed limit; the Klein disc of RMG1 is
velocity space; the return angle of PH1–PH2 is the rotation left by
non-collinear boosts. A fluid's response path is a path of velocities.

## 3. The continuum law and its boost

**T4.** For a unit split generator A (A² = 1) and the circular R,

```text
∂_t ψ = ( c A ∂_x + g R ) ψ      ⇒      ∂_t² ψ = c² ∂_x² ψ − g² ψ ,
```

because A and R anticommute. With B = Exp(ηA/2): B·1·B = cosh η + A sinh η,
B·A·B = A cosh η + sinh η, B·R·B = R. The law keeps its form under the
boost built on **its own** c, with the mass term untouched; a sector of
another speed is not invariant under that boost (exact witness). Each
sector has its own relativity; a universal one needs one speed shared by
all sectors.

## 4. Mass and diffusion: the two parities of one flip

The two components of ψ are the two light-like readings (movers at ±c).
A flip between them at rate g, in its two parities (QC3):

```text
turn    g R          (cA∂_x + gR)² = c²∂_x² − g²        frequency² = c²k² + g²           mass
record  g (S − 1)    (∂_t + g)² = c²∂_x² + g²           ∂_t² + 2g∂_t = c²∂_x²            diffusion D = c²/(2g)
```

**T5.** D · 2g = c². The rate that is the mass in the turn parity is the
rate that sets the diffusion constant in the record parity. This is the
calibration QC3 left open: the turn rate is the mass term of the sector.
With g = mc²/ħ it reads D = ħ/2m (the identification of g is the named
classical input; the identity D·2g = c² is exact).

## 5. Certificate

T1–T2 as Laurent-matrix identities in the shift for five rational coins
(θ = 0 and π/2 included), dispersion and group-speed bound on five
rational address turns each; R28's coin without the root. T3 with exact
products against PH2's polar formula. T4 for two sectors, two rational
boosts, foreign-sector witness. T5 with exact square laws and a bracket
of the slow rate between Dk² and 2Dk². Eight tests; a coin that is not a
turn, an unconditioned shift, a non-boost and a split mass generator are
rejected.

## 6. Claim boundary

```text
SPEED² + (CURVATURE/2)² = 1 FOR THE COIN-TURN WALK; REST TURN θ        PROVED
RESPONSE SPACE = VELOCITY SPACE; RETURN = BOOST ROTATION               PROVED
SECTOR LAW, OWN-BOOST COVARIANCE, FOREIGN SECTOR NOT COVARIANT         PROVED
D · 2g = c² ACROSS THE TWO PARITIES                                    PROVED
THE SPEED OF LIGHT, OR ANY PHYSICAL SPEED, DERIVED                     NO — θ is free; a sector is a choice of coin
WHY ALL PHYSICAL SECTORS WOULD SHARE ONE SPEED                         OPEN
GRAVITY, CURVED SPACETIME                                              NOT TOUCHED
RELEASE OF STORED CURVATURE AFTER A CUT AS A PROPAGATING RIPPLE        NOT BUILT — next stage
```

## 7. Reproduce

```text
python pr1_sector_speed_and_mass.py
python -m unittest test_pr1_exact
```
