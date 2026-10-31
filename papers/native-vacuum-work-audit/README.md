# First vacuum-work feasibility calibration

Author: Monty. 31 October 2026.

**Outcome:** the supplied finite calibration returns zero net energy gain
for each complete ideal drive/recovery cycle. A passive observer-basis change
preserves the energy. The test demonstrates energy bookkeeping on the existing
calculator; it does not demonstrate physical zero-point energy extraction.

One fresh NRQC request, pinned to `fdd39048a6ddaa2dd78aa3197c274d7766f31a0e`,
evaluated ten independent one-qubit cases and 26 queries. **50 checks and six
negative controls pass**, in 0.883 s including verification. No hardware,
Casimir experiment, QET protocol or gravitational field was run.

For the declared Hamiltonian H_E=E_off I+Delta(I-Z)/2, these are the results
in units of the supplied positive level gap Delta:

| Rotation | Drive work | Ideally recoverable work | Exact net gain after return |
|---|---|---|---|
| pi/2 | 1/2 | 1/2 | 0 |
| pi/4 | [0.146446609,0.146446610] | Same enclosure | 0 |
| pi/8 | [0.038060233,0.038060234] | Same enclosure | 0 |

The displayed intervals are outward rounded. The exact zero follows from
W_recover-W_drive=E_initial-E_return and the algebraically exact returned
ground state; it is not a subtraction of independently rounded decimals.
Recovery is an ideal energy decrease available to the prescribed inverse
control, not a physically measured load output or an efficiency measurement.

Offset references E_off/Delta=0,1/2,100 change the absolute readings but not
the cycle work. This fixed-background calculation does not claim that a
vacuum stress term can be removed from gravitational equations. A wrong
return pulse leaves one full gap of excitation, so it fails cycle closure.
Ignoring the drive input gives a false apparent positive yield, as the
retained negative controls show.

## Relation to the research model

The preceding LMD results concern dimensionless probability marginals. A
physical autonomous energy Hamiltonian, actuator, detector/record reset and
load are not yet derived there. The Hamiltonian above is an explicitly added
calibration model. The framework's cut/zero terminology alone supplies no
identification with a physical vacuum's ground-state energy.

[PROTOCOL.md](PROTOCOL.md) gives the basis-transformation identity, the
time-dependent connection term, the ground-state passivity boundary and
primary sources for driven Casimir radiation, quantum energy teleportation
and semiclassical black-hole radiation. These precedents motivate controlled
conversion or recovery tests while identifying their energy sources.

The next physical test should specify a field/oscillator plus controller
and load, track all drive and reset costs, and compare phase-feedback schemes
at the same supplied energy. Ground-state passive changes, drive-off runs,
matched incoherent controls and a full returned-cycle ledger should accompany
any claimed improvement. No vacuum-only net-power result follows from this
calibration or the earlier probability separation.

## Evidence and reproduction

- Protocol commit, published before execution:
  `eb641783c0655ed28e81cdfbbb31314d77f7786a`.
- [Raw exact calculator output](results/calculator.result.json),
  [native verification](results/verification.json),
  [execution/source hashes](results/execution.json), [manifest](results/manifest.json).
- [Cases](cases.json), [input](request.json), [native verifier](verify.cjs).
- Python 3.12.14; Node v24.22.1; one fresh calculator request; no CI dispatch.

From the repository root, use a clean checkout of the pinned calculator and
a new output directory:

```sh
python3.12 papers/native-vacuum-work-audit/run.py \
  --nrqc-root /path/to/nrqc-pinned \
  --out /tmp/vacuum-work-replay \
  --protocol-commit eb641783c0655ed28e81cdfbbb31314d77f7786a
```

The runner checks the calculator commit and clean state, inherited file
hashes, execution-source hashes and calculator Python source hashes. Native
rational intervals perform the scientific verification; Python orchestrates
and hashes. All prior QG packets and the pending joint-law solver are preserved.
