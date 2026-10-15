# MZC-1/2 — Moving zero, circular operators and the RH sector obligation

The balance set moves with a transported cut: `ker(C0 U^-1)=U ker(C0)`.
The numerical zero remains zero. Affine references, operator similarity,
form congruence and the transported metric are kept distinct.

The application to the existing RH work is precise: moving zero sets and
curvature do not determine the source sign. The next-observation Schur margin
transforms as `sigma' = |a|^2 sigma` under an acquired-frame change preserving
the old subspace. Its sign cannot be repaired by changing coordinates.
This reuses SE49's existing pivot-scaling theorem with an explicit block proof;
a better frame can still make previously inconclusive interval bounds decisive.

MZC-2 reuses SE26's lossless circular operator and SE29's spectral sectors.
RH positivity becomes absence of actual spectral support in the oriented lower
arc. Algebra grades, projector membership and spectral signs are distinguished.
Positive sector compressions require a reducing decomposition to certify the
whole operator; circular composition need not preserve the positive arc.

Replay: **139 exact checks, 3 native symbolic proofs/replays and 14 negative controls**.

- [Written statements, proofs and RH obligation](THEOREM.md)
- [Circular operators and sector-relative positivity](CIRCULAR_SECTORS.md)
- [Pinned source mapping](SOURCE_MAP.md)
- [Exact native controls and symbolic witnesses](CERTIFICATE.json)
- Replay: `python3 papers/native-moving-zero-covariance/run.py --check`

This is a reusable covariance/record-completeness result applied to the RH
adapter. No new zeta values, universal innovation bound or physical quantum
gravity completion is claimed. The original source definitions are retained.
