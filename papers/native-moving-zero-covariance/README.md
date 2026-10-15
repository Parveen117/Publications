# MZC-1 — Moving zero covariance and the RH sign obligation

The balance set moves with a transported cut: `ker(C0 U^-1)=U ker(C0)`.
The numerical zero remains zero. Affine references, operator similarity,
form congruence and the transported metric are kept distinct.

The application to the existing RH work is precise: moving zero sets and
curvature do not determine the source sign. The next-observation Schur margin
transforms as `sigma' = |a|^2 sigma` under an acquired-frame change preserving
the old subspace. Its sign cannot be repaired by changing coordinates.
This reuses SE49's existing pivot-scaling theorem with an explicit block proof;
a better frame can still make previously inconclusive interval bounds decisive.

Replay: **108 exact checks, 3 native symbolic proofs/replays and 10 negative controls**.

- [Written statements, proofs and RH obligation](THEOREM.md)
- [Pinned source mapping](SOURCE_MAP.md)
- [Exact native controls and symbolic witnesses](CERTIFICATE.json)
- Replay: `python3 papers/native-moving-zero-covariance/run.py --check`

This is a reusable covariance/record-completeness result applied to the RH
adapter. No new zeta values, universal innovation bound or physical quantum
gravity completion is claimed. The original source definitions are retained.
