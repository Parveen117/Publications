# GE1: Generalized Euler evolution with certified domains

3 October 2026. Runtime: **Python 3.12 only**.

Two unpublished drafts suggested identifying Generalized Euler phase space
with exponential change evolution. Their useful part is made precise here:
which completed carrier an exponential acts on, its generator domain, and
which rational approximations remain stable when the mode cutoff grows.

The RH source already has an infinite Euler-scale flow and energy completion.
GE1 adds scoped domain/approximation proofs and draft corrections; it does
not propose a replacement for that carrier or rename its existing results.

**Result:** real phase evolution has a closed native generator; a reversible
midpoint approximation converges strongly without coupling its step count to
the mode cutoff. The existing native YM heat generator has a closed graph
domain and a positive contractive resolvent approximation with an explicit
core error. A near-resonant phase example has positive decay at every finite
cutoff but no positive infinite-content gap.

Read the [written proofs](THEOREM.md), the [complete draft audit](SOURCE_AUDIT.md)
and the [Yang–Mills application boundary](YM_DOCK.md). These are mathematical
tools and corrections, not a new mass-gap or Riemann Hypothesis proof.

The finite exact certificate checks algebra, stability, adversarial domain
examples, native noncommuting coefficient actions and approximation bounds.
The infinite statements depend on the written proofs and their hypotheses;
they are not mechanically formalized or externally expert-certified.

```bash
python3.12 -B papers/generalized-euler-evolution/certificates/ge1_evolution.py --check
python3.12 -B -m unittest discover -s papers/generalized-euler-evolution/tests -v
```

Default verification is read-only. Only `--write` regenerates GE1 evidence.
Older research and certificate files remain unchanged. Public terminology is
**Generalized Euler** or **Euler** throughout this packet.
