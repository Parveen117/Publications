# IN1 — The invariant: observed and lost together keep the uncut, in every frame

Monty Dabas. 7 October 2026. Python 3.12, exact arithmetic over the
Gaussian rationals. Continues FR1.

Owner's statement for this stage: an invariant is an invariant whatever
the frame; only the reading changes; the lost and the observed
information together keep the uncut conserved. The identification left
open in FR1 will be derived.

Sources read before building: **T24** Theorem 6.1 (S = R + D: source =
recognized + memory, for every cut), **GE2-T2** (P + Q = I; memory
through the complement cut), **EMK-1** T2 (the determinant as a sum of
channels), **FR1-T3/T4** (three cuts on the cut-complex carrier, no
fourth), **SY1-S6/S7**, **MS1-T1**, **PR2-T2**.

## 1. The reading tensor

On the cut-complex two-component carrier, with the three cuts
C₁ = K, C₂ = RK, C₃ = ιR, every self-dagger tensor is

```text
ρ = ½ ( n + r₁C₁ + r₂C₂ + r₃C₃ ) ,        n the uncut count ,   r_i the reading of cut i .
```

**T1.** det ρ = ¼ ( n² − r₁² − r₂² − r₃² ).

## 2. A single reading loses nothing

**T2.** For ρ = ψψ†: det ρ = 0, so

```text
n² = r₁² + r₂² + r₃² ,          1 − (r_i/n)² = Σ_{j≠i} (r_j/n)² .
```

What is lost to one cut is exactly what the other two observe. The
"memory" of MS1 was this: lost to one cut, held by the complementary one.

**T3 (the count of cuts is the count needed to lose nothing).** For a
real ψ the third reading is zero and two cuts are complete. For a
cut-complex ψ two cuts are not complete (witness: n = 15, readings
−5, −2, −14); three always are; and there is no fourth (FR1). The
dimension of the frame is the number of readings after which nothing of
a single reading is lost. This is FR1's identification, now a statement
about conservation.

## 3. In every frame

**T4.** Under every frame change ρ → gρg† with det g = 1 — turns, boosts
in any direction, the mixed ones — the readings change and

```text
n² − r₁² − r₂² − r₃²
```

does not. The frame-independent form of "observed and lost keep the
uncut" is a form with one plus and three minus signs. It comes from the
determinant identity of the block; nothing about space or time was put
in.

## 4. What cannot be recovered

**T5.** A record p·ρ₁ + (1−p)·ρ₂ of two single readings has

```text
n² − r·r = 2 p (1−p) ( n₁n₂ − r₁·r₂ ) ≥ 0 ,
```

zero exactly when p is 0 or 1, or the two readings are parallel. So for
every record

```text
(observed by all three cuts)²  +  (what no cut of the frame observes)  =  (uncut)² ,
```

and the second term is the same in every frame. With n read as energy
and r as momentum (PR2-T2) it is the square of the rest mass: a single
reading is light-like; a record of two non-parallel light-like readings
has a rest mass whose square is the share law 2p(1−p) times
(n₁n₂ − r₁·r₂).

```text
recoverable memory     lost to one cut, observed by the others     frame-dependent, sums to the uncut (T2)
unrecoverable memory   observed by no cut of the frame             frame-independent: the invariant (T5)
```

## 5. Reading

- There is one conserved quantity, the uncut. Every cut divides it into
  observed and lost, differently in every frame.
- Information lost to a cut is not gone if a complementary cut holds it.
  Three cuts hold all of a single reading.
- What a record loses, no cut holds — and that amount is what every
  frame agrees on. Rest mass is the information that cannot be read back.
- The form with signature (1, 3) is this conservation law written so
  that no frame is preferred.

## 6. Certificate

T1–T3 on six rational states (two real, four cut-complex): the tensor is
rebuilt from its readings, det = 0, the complement identity for each
cut, completeness with two or three cuts. T4 on three single readings
and one record through five frame changes of determinant one (two
boosts, a turn, a phase turn, a shear with ι). T5 on four pairs and five
shares each, including a parallel pair. Six tests; a missing third cut,
a frame change of determinant two and a wrong law are rejected.

## 7. Claim boundary

```text
det ρ = (n² − r·r)/4 ; SINGLE READING NULL ; LOST TO ONE CUT = OBSERVED BY THE OTHERS    PROVED
THREE CUTS COMPLETE, TWO NOT, FOR CUT-COMPLEX READINGS                                   PROVED
n² − r·r THE SAME IN EVERY FRAME                                                         PROVED
RECORD: n² − r·r = 2p(1−p)(n₁n₂ − r₁·r₂) ≥ 0                                             PROVED
n AS ENERGY, r AS MOMENTUM, THE INVARIANT AS REST MASS²                                  READING (PR2-T2); the algebra does not need it
A LAW OF MOTION; A SCALE; ANY CONSTANT                                                   NO
THAT NATURE'S READINGS ARE CUT-COMPLEX                                                   NOT DERIVED — it is the framework's carrier
```

The correspondence between two-component tensors and a (1,3) form is
known (general knowledge). The stage's content is its origin here — the
certified determinant identity read as conservation of the uncut — and
the split of memory into recoverable and unrecoverable.

## 8. Reproduce

```text
python in1_the_invariant.py
python -m unittest test_in1_exact
```
