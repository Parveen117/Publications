# Information Invariance — paper under development

`paper/information_invariance_v1_original.tex` — the original draft (kept for provenance).
`paper/information_invariance_v2.tex` — the working paper: proved finite structure separated from conjectured dynamics.
`certificates/l1_certificates.py` — executable certificates for Theorems A and B (exact rational arithmetic).
`LEDGER.md` — what is PROVED / OPEN / CONJECTURE / RETRACTED, with named obligations.

Rule for this folder: a sentence may claim a result only if the ledger lists it as PROVED with a reproducing certificate. Everything else is written as a conjecture with its obligation.

```bash
python papers/information-invariance/certificates/l1_certificates.py
python -m pytest papers/information-invariance/tests
```
