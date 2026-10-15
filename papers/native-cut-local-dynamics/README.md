# CLD-1 — Dynamics after the cut

Local here means the information sector retained after an observation cut. This packet tracks its evolution as the cut itself changes.

The finite all-state closure criterion is `D=P'+(i/eta)[H,P]=0`. In a moving cut frame the generator is `W*(H-i eta[P',P])W`; omitting its connection can confuse detector motion with source dynamics. Outside-sector return contributes an exact memory term even after an initially selected preparation.

For AGD's infinite geometry, a direct common-core calculation gives the moving-reflection coupling and proves that a rigidly rotating reflection detector is closed for all states only at a stationary integer-pi angle. Starting in cosine mode `n`, angular detector speed `omega` produces initial leakage coefficient `omega^2 n^2`. The previously certified fixed reflection invariant survives unchanged.

Replay result: **167 exact checks, 15 native word certificates/replays and 14 negative controls**.

- [Theorem, proofs and domain boundaries](THEOREM.md)
- [Native certificate](CERTIFICATE.json)
- [Source pins](SOURCE_PINS.json)
- Replay: `python3 papers/native-cut-local-dynamics/run.py --check`

GE3's return memory, NJC's connection and ZCR's fixed-cut estimates are explicitly reused. Finite certificate replay does not formally certify infinite analytic proofs. Spatial gravitational constraints and physical quantum-gravity completion remain open; cut-local dynamics specifies a different, necessary information contract.
