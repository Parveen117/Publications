# Information Invariance — paper under development

`paper/information_invariance_v1_original.tex` — the original draft (kept for provenance).
`paper/information_invariance_v2.tex` — the working paper: proved finite structure separated from conjectured dynamics.
`certificates/l1_certificates.py` — executable certificates for Theorems A and B (exact rational arithmetic).
`LEDGER.md` — what is PROVED / OPEN / CONJECTURE / RETRACTED, with named obligations.

## Certificates

| file | theorems |
|---|---|
| `l1_certificates.py` | A (invariance = flat sector), B (one cycle decides nothing) |
| `d2_phase_origin.py` | C (real flow is Schrödinger iff), D (one-channel obstruction) |
| `d2a_cut_supplies_structure.py` | E (cut supplies conjugate pair; phase is the odd part) |
| `d2b_linearisation_verdict.py` | D2b verdict: damped rotation, not Schrödinger |
| `d3_regulator_scale.py` | F (no Planck length without ħ), G (own regulator is cosmological) |
| `d3prime_emergent_einstein.py` | K (DOF no-go), L (Clausius route imports ħ), M (source not conserved) |
| `d4_decoherence_audit.py` | H (printed rate is not a rate; repair is Diósi–Penrose) |
| `d1_selection_principle.py` | I (invariance selects nothing), J (minimal scale usage selects n=2) |
| `micro_chi_characterisation.py` | N (δ𝓘=0 ⟺ χ≡1), O (curvature sufficient, not necessary) |

Rule for this folder: a sentence may claim a result only if the ledger lists it as PROVED with a reproducing certificate. Everything else is written as a conjecture with its obligation.

```bash
python papers/information-invariance/certificates/l1_certificates.py
python -m pytest papers/information-invariance/tests
```
