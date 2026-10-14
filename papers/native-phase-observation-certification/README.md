# Native phase observation — executed certification

**PASS — 14 October 2026.** The six newly uploaded source packets have been
rerun, and the Phase Memory v0.7 and Angular Current v0.8 drafts now have a
native executable companion: **60 algebraic certificates, 60 serialized
replays, 373 exact checks and 25 negative controls**. All five separately
uploaded PDFs match their source-packet copies byte for byte.

The certificate covers the upload batch through commit
[`611acc1`](https://github.com/Parveen117/Publications/commit/611acc114db445496de2b1101c22049fc5a95cfc).
Its inputs and engine bytes are fixed by [SOURCE_PINS.json](SOURCE_PINS.json).
The original uploaded documents are preserved. Their statements about
deferred testing describe the original draft stage; this companion records
the subsequent execution and its precise scope.

## Results by source

| Uploaded development | Algebraic certificates / replays | Exact checks | Negative controls | Result |
| --- | ---: | ---: | ---: | --- |
| Clock-Free Native Lattice | 6 / 6 | 257 | 6 | PASS, existing packet rerun |
| Clock-Free Cut Quantisation | 10 / 10 | 354 | 8 | PASS, existing packet rerun |
| Clock-Free Physics Cubits | 13 / 13 | 280 | 9 | PASS, existing packet rerun |
| Clock-Free Scale Compass | 10 / 10 | 610 | 10 | PASS, existing packet rerun |
| Clock-Relation Reconstruction | 7 / 7 | 299 | 15 | PASS, existing packet rerun |
| Native Phase-Word Reconstruction | 13 / 13 | 952 | 13 | PASS, existing packet rerun |
| Phase Memory + Angular Current | 60 / 60 | 373 | 25 | PASS, new scoped certificate |

Counts retain each packet's own convention and can include overlapping
regression coverage. The new 373 checks include presentation audits and
serialized replay checks; the 25 negative controls are counted separately.
The wrapper additionally verifies source-mutation rejection and that the
inherited read-only runs change no extracted bytes.

## What the new run establishes

- The native doubled-phase channel and relative Cayley identities, with the
  anchored integer update exercised across forward, reverse and mixed paths.
  The draft's forward example recovers records `0,0,0,0,1,1,1,1,2`.
- Exact varying-common-gain cancellation, the native operator/current relation,
  and the Cayley differential. The unit-pair tangent rule passes the native
  derivation gate; an incompatible derivative fails.
- The product, dagger and relative-current identities, with 30 exact cases
  using varying gains and a shared rotating frame. Offsets, unequal coordinate
  gains and noncommon frame rotations fail the claimed invariance.
- Endpoint-jet obstruction polynomials through order 6. Matching finite
  endpoint jets can hide an accumulated native phase of pi; the finite data
  do not provide the missing integral.
- Every entry of the displayed frozen two-axis loop coefficient and its
  reversed-order sign. Dependent directions have zero area response.

**Scope:** these are native algebraic certificates and exact finite executions.
General path lifting, integration, periodicity and all-orders claims have
explicit written deductions and dependencies in
[CLAIM_COVERAGE.md](CLAIM_COVERAGE.md). They are not all formalized in the finite
rewriter. The physical source law, an actual angular-current accumulator and
gravitational dynamics remain outside this certificate.

## Reproduce

Use **Python 3.12** and **Node v24.19.0**. No third-party Python or Node packages,
external checkout, network access, or floating-point math are needed for a run.
From the Publications repository root:

```sh
python3.12 papers/native-phase-observation-certification/verify.py --check
```

The default mode is also `--check`. The wrapper validates 18 pinned inputs,
extracts each committed ZIP to a temporary directory, runs its existing
verifier, compares the five PDF copies, and executes the new certificates
using the unchanged native core embedded in the phase-word packet.
It recomputes [CERTIFICATE.json](CERTIFICATE.json) and requires byte-for-byte
agreement, including [EXPECTED.sha256](EXPECTED.sha256). Check mode does not
modify repository files.

Only when intentionally generating a new reviewed certificate:

```sh
python3.12 papers/native-phase-observation-certification/verify.py --write
```

`--write` writes this companion's certificate and checksum. It does not update
source pins, change any draft, replace the native engine, or refresh an
inherited result. A changed pinned source stops either mode.

The pinned native engine is RKF commit
[`e8e3089`](https://github.com/Parveen117/Recognition-Kernel-Framework/commit/e8e3089745bf7dc0532d01c62a5f660f5658c6f0).
[The CI workflow](../../.github/workflows/native-phase-observation-certification.yml)
runs the same read-only command on relevant pushes and pull requests.

## Read the evidence

- [Claim coverage and written proofs](CLAIM_COVERAGE.md): PM1/PM2, AC1–AC4,
  domain restrictions, finite-jet obstruction, and inherited curvature reuse.
- [Native certificate](CERTIFICATE.json): presentations, critical-pair
  witnesses, equation witnesses, finite observations, negative controls and
  all six fresh upstream results with their hashes.
- [Native implementation](verify.cjs) and [orchestration](verify.py).
- Source manuscripts:
  [Phase Memory](../../Native_Phase_Memory_Development_Draft.tex) and
  [Angular Current](../../Native_Angular_Current_Development_Draft.tex).
