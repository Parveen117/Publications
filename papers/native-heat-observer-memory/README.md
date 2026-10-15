# HCM-1/2 — Source-derived memory and phase-aware response transport

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

Replay: **189 exact checks, 2 native symbolic proofs/replays and 21 negative controls**.

- [Source-derived calculation, proofs and boundaries](THEOREM.md)
- [Phase-aware record, EMK transport and UGD interface limits](PHASE_RECORD.md)
- [Source map](SOURCE_MAP.md)
- [Native certificate](CERTIFICATE.json)
- Replay: `python3 papers/native-heat-observer-memory/run.py --check`

Native exact arithmetic and symbolic witnesses are reused. Python 3.12 handles
I/O. General exponential and memory statements have written proofs; no general
proof-assistant or physical four-dimensional Yang–Mills certification is claimed.
