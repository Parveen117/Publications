# PR2 — Mass and antimass are the two sheets; the observer who sits with both reads only invariants

Monty Dabas. 7 October 2026. Python 3.12, exact rational arithmetic only.
Continues PR1.

Owner's idea for this stage: in the same accelerated frame, an observer
who sits together with the antimass keeps the curvature — or the
information — invariant.

Sources read before building: **PR1** (sector law, the mass term on R),
**R28** (roles H, K, R = KH), **RMG9-T1** (square law), **QC3-T1**
(E² − O² = I), **QC4-T1** (return = −variance), and from the
information-invariance ledger in Publications **Theorem E** (phase lives
in the cut-odd part, damping in the even part; the conjugate channel is
the opposite sheet of the oriented cut) and **Theorem N** (δ𝓘 = 0 ⟺ χ ≡ 1,
time-reversal invariance of the weight).

## 1. The form and the two sheets

Elements X = xH + yK + tR obey X² = x² + y² − t². Frame changes act by
similarity with boosts Exp(ηA/2) (A a unit split element) and turns
Exp(φR/2).

**T1.** Every frame change preserves X². For timelike X (X² < 0) it also
preserves the sign of t. The timelike elements of a given X² form two
sheets, and no boost history, collinear or not, connects them.

**T2 (the boosted mass is energy and momentum).**

```text
Exp(ηA/2) · gR · Exp(−ηA/2) = g ( cosh η · R + sinh η · A′ ) ,   A′ a unit split element,
energy = g cosh η ,   momentum = g sinh η ,   energy² − momentum² = g² ,   velocity = tanh η .
```

**T3 (antimass).** The reversed turn Exp(−θR) is −gR: the mirror point on
the other sheet, with the same X², the same energy² and momentum². The
second-order law ∂_t² = c²∂_x² − g² of PR1 does not see which sheet.

## 2. The pair observer

Let the record be {Exp(θR), Exp(−θR)} with weights (p, 1−p) — mass and
antimass read together — and carry both through any sequence of frames.

**T4.**

```text
mean    = cos θ + (2p − 1) sin θ · X          X the boosted unit generator
memory  = 4 p (1 − p) sin²θ                   (the variance of the record, QC4)
```

- The even part cos θ and the memory are **scalars**: the same in every
  frame and along every acceleration history.
- The odd part is a vector on the sheet: frame-dependent (it is the
  energy–momentum of T2), with invariant square −((2p−1) sin θ)².
- The weight p — the information about which sheet — is not changed by
  any frame change.
- At p = ½ the odd part is zero in every frame. Nothing the balanced pair
  observer reads depends on the frame.

So the answer to "curvature or information?" is: both are invariant. The
memory (curvature) and the weights (information) are scalars for every p.
What acceleration changes is only the odd part, and the observer who sits
with the antimass in equal weight has no odd part.

**T5 (this is information invariance).** The balanced pair has equal
weight on a record and its reverse: χ ≡ 1. By Theorem N that is δ𝓘 = 0.
By Theorem E the frame-dependent part is exactly the cut-odd (phase)
part, which the balanced record does not carry. A single reading (p = 1)
has χ ≠ 1, zero memory, and the full frame-dependent odd part.

```text
            even part    memory            odd part                   information invariance
p = 1       cos θ        0                 sin θ · X  (moves)         no
p = 3/4     cos θ        (3/4) sin²θ       ½ sin θ · X (moves)        no
p = 1/2     cos θ        sin²θ             0 in every frame           yes
```

The share law of PH3 reappears: memory 4p(1−p), largest for the balanced
pair, zero for either pure sheet.

## 3. Turn against split

**T6.** The generator αA + gR has square α² − g²: for |α| < g it turns
(reduced rate √(g² − α²)), at |α| = g it stops (the dual sector), beyond
it runs away. A split rate equal to the turn rate is the boundary.

If α is read as the boost rate a/2c of a uniformly accelerated frame and
g = mc²/ħ, the boundary is a = 2mc³/ħ. That reading is **not**
established here: in an accelerated frame the split term can be absorbed
by a position-dependent rescaling, and this stage does not decide whether
any of it survives. The value coincides with a maximal acceleration
proposed in the literature (general knowledge, unchecked).

## 4. Certificate

T1 on four elements through five frames (three non-collinear boosts, a
turn, a reverse boost): form, sheet, and the square law. T2–T3 at three
rational rapidities. T4 for p = 1, 3/4, 1/2, 0 through all six frame
stages: even part, odd square, variance, reverse-pair property, and that
the odd part vanishes in every frame exactly at p = ½. T6 in the three
sectors. Seven tests; a non-boost and a one-sided frame change are
rejected.

## 5. Claim boundary

```text
TWO SHEETS; FORM AND SHEET INVARIANT UNDER ALL FRAME CHANGES          PROVED
BOOSTED MASS GENERATOR = (ENERGY, MOMENTUM)                           PROVED
PAIR RECORD: EVEN PART AND MEMORY SCALAR; ODD PART A VECTOR           PROVED
BALANCED PAIR ⇒ χ ≡ 1 ⇒ δ𝓘 = 0, NOTHING FRAME-DEPENDENT               PROVED (with Theorems E, N as cited)
GRAVITY COUPLES ONLY TO THE EVEN PART                                 NOT CLAIMED — the second-order law is sheet-blind; that is all
CRITICAL ACCELERATION 2mc³/ħ AS PHYSICS                               NOT ESTABLISHED
CURVED SPACETIME, EQUIVALENCE PRINCIPLE                               NOT TOUCHED
```

## 6. Reproduce

```text
python pr2_two_sheets_pair_observer.py
python -m unittest test_pr2_exact
```
