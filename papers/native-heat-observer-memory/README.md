# HCM-1 — Observer memory in the actual YM-54 free heat protocol

The published thermo-response fixture gives C=diag(5/162,5/72,0).
On an exact three-mode coefficient block its heat generator is tr(C)I-C.
A declared mixed observer has a strictly positive two-time closure defect;
its leading coefficient is 1/2916. One additional retained channel repairs
closure exactly, with memory kernel exp(-t/18)/2916 and tail budget
exp(-T/18)/162 times the past-readout bound.

This uses the existing source generator rather than a newly selected toy
dynamics. The observer is selected explicitly. The result is for a free
single-site coefficient block, not the full interacting instantaneous row.

Replay: **68 exact checks, 1 native symbolic proof/replay and 10 negative controls**.

- [Source-derived calculation, proofs and boundaries](THEOREM.md)
- [Source map](SOURCE_MAP.md)
- [Native certificate](CERTIFICATE.json)
- Replay: `python3 papers/native-heat-observer-memory/run.py --check`

Native exact arithmetic and symbolic witnesses are reused. Python 3.12 handles
I/O. General exponential and memory statements have written proofs; no general
proof-assistant or physical four-dimensional Yang–Mills certification is claimed.
