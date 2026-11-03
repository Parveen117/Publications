# Critical native return and a third-probe prediction

Monty Dabas · CR-1–CR-5 · 30 September 2026

[The proof](THEOREM.md) decomposes the exact native return into complementary
aperture channels. Its own algebraic inverse has a simple pole at cut closure,
with normalized coefficient (1+2b)P−. A source or observer that excludes P−
does not see the pole. An inverse on the surviving P+ corner is not a full
inverse at the singular point.

The exact identity (1-x)[1+bz(1+x)]=1-z produces a stable interval for the
normalized pole response and separates proximity error from aperture error.
Native products and the unchanged R2 aperture solver verify these statements.

Using the previously certified two-reading example b=5/7, delta=3/4, a third
relative probe predicts the positive root of 55x²+56x−132=0, without another
fitted parameter. This is a model calculation, not a performed experiment.
Its inverse-pole coefficient is 17/7, not the fine-structure constant.

**Alpha is NOT DERIVED.** The inverse of a native boundary return has not been
identified with a physical photon propagator. Momentum, source, action and
transverse-field matching remain explicit obligations.

## Reproduce

```bash
node papers/native-critical-response/verify.cjs --rkf-root /path/to/Recognition-Kernel-Framework --check
```

Pin: 3cc5a33b05c16d59c90994ddda69dedc0d392424. Default/check is read-only;
--write explicitly creates this packet's source-bound evidence. Canonical
native arithmetic, Pāṇinian proof replay and the original R2 solver are
imported unchanged. Finite evidence, general written proofs and physical
validation are separate. The documented upstream master-pin limitation is
unchanged; this packet does not claim full engine recertification.
