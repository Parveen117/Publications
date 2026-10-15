# HCM-1/2/3 — Source memory, phase records and common-gain recovery

The published thermo-response fixture gives C=diag(5/162,5/72,0).
On an exact three-mode coefficient block its heat generator is tr(C)I-C.
A declared mixed observer has a strictly positive two-time closure defect;
its leading coefficient is 1/2916. One additional retained channel repairs
closure exactly, with memory kernel exp(-t/18)/2916 and tail budget
exp(-T/18)/162 times the past-readout bound.

This uses the existing source generator rather than a newly selected toy
dynamics. The observer is selected explicitly. The result is for a free
single-site coefficient block, not the full interacting instantaneous row.

HCM-2 encodes the retained readings as z=q+i h and derives the exact source
law z'=(-65z+(7+24i)conjugate(z))/1296. A reversible quarter-turn/scale/residual
record carries this law with explicit chart transitions. The residual contains
necessary within-sector information; finite phase labels alone do not close.
The native EMK reflection S=(7K+24RK)/25 recovers both inherited heat rates.
The full historical UGD numeral algebra and its seam ledger are not identified
with this response encoding; their operation and winding boundaries are explicit.

HCM-3 uses the source's coprime gains A=u^9, B=u^4 to recover u=A/B^2 without
roots. Prime valuations give an optional multiplicity audit. The paired-phase
theorem is distinguished from the actual heat response phase. Rational u gives
an exact rational transfer/corrected readout, with explicit clock/noise bounds;
the u=1/2 inverse has worst-case amplification 512.

Replay: **265 exact checks, 3 native symbolic proofs/replays and 28 negative controls**.

- [Source-derived calculation, proofs and boundaries](THEOREM.md)
- [Phase-aware record, EMK transport and UGD interface limits](PHASE_RECORD.md)
- [Coprime gains, prime audits and rational reference correction](PRIME_CLOCK_MAP.md)
- [Source map](SOURCE_MAP.md)
- [Native certificate](CERTIFICATE.json)
- Replay: `python3 papers/native-heat-observer-memory/run.py --check`

Native exact arithmetic and symbolic witnesses are reused. Python 3.12 handles
I/O. General exponential and memory statements have written proofs; no general
proof-assistant or physical four-dimensional Yang–Mills certification is claimed.
