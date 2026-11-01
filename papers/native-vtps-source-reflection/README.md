# VSR-1 — VTPS reflection of the actual five-reading source

VTPS=(V,T,P,S). This packet returns to TLM's actual cubic source, whose
mechanical reading at the cut is 5, and realizes its two signed directions.

- At mu=+5 and mu=-5, the cut reading is 5 and curvature is 125/324.
- The first-layer six-by-six matrices are related by the explicit grading
  D=diag(1,1,-1,-1,-1,-1); all finite depths have the same graded similarity.
- In the actual (V,T) control chart, Maxwell is symmetry of the (P,S)
  response matrix.
- The positive source metric becomes diagonal while its curvature remains
  nonzero. The transported connection, including Q^(-1)dQ, retains it.
- Three off-cut mirror pairs distinguish raw readings from calibrated ones.

[Proofs](THEOREM.md) · [Source map](SOURCE_MAP.md) · [Execution](USAGE.md)

The native adapter and verifier use the unchanged F/Cut and TLM engines.
CERTIFICATE.json, SOURCE_PINS.json, EXPECTED.sha256 and CI_PROVENANCE.json
are generated on GitHub, then imported unchanged. Source-specific reflection
is proved; arbitrary-source symmetry and physical identification are not.
