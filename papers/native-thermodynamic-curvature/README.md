# Native thermodynamic curvature: NT-1 to NT-8

This native edition connects the existing thermodynamic paper and theorem
register to EMK connection calculus. It constructs the response element and
its connection inside the native algebra, proves the exact Hessian/Riemann
representation, and retains Smriti and observer-return terms.

- [Paper and written proofs](THEOREM.md)
- [Compiled manuscript](pdf/native-thermodynamic-curvature.pdf)
- [LaTeX manuscript](source/main.tex)
- [Original-to-native theorem map](SOURCE_MAP.md)
- [Generated finite certificate](CERTIFICATE.json)
- [Source pins](SOURCE_PINS.json)

The central result is `F_H = -[H^-1 d_i H, H^-1 d_j H]/4`, with the
connection selector and derivative frame explicitly stated. The native open
curvature is `F(A-M) = F(A) - D_A M + M wedge M`. Balanced open curvature
can vanish while the raw curvature remains nonzero. Constant reciprocal
Onsager response is a flat process-form sector; a variable reciprocal matrix
needs the derivative terms too.

This is an additive conversion of the thermo core, not a replacement arXiv
submission. The existing source editions, historical certificates and shared
RKF engine remain preserved. SOURCE_MAP.md specifies exactly which prior
results are represented, extended, or retained with their original hypotheses.

## Reproduce

Use Node 24.19.0 and Python 3.12. No third-party Python package is needed.
Check out the two public source repositories at the commits in SOURCE_PINS.json:

```bash
bash papers/native-thermodynamic-curvature/verify.sh \
  --rkf-root ../Recognition-Kernel-Framework \
  --thermo-root ../thermo-source --check
```

`thermo-source` is Publications at research commit
`f0f1c1650e914b7eaf3bf64ddb7df9f543c92a56`. The publication checkout running
the verifier supplies the local EMK-C1 chapter and this packet. `--check` is
read-only. `--write` explicitly generates a new candidate certificate after
all proof controls pass; it never changes source pins. The CI workflow runs
the same command using fixed commits and checks the native compiler, independent
Christoffel calculation, failure controls and artifact integrity.

The Markdown paper is the authoring surface. To regenerate the LaTeX and PDF,
run `python3.12 papers/native-thermodynamic-curvature/build_paper.py` from the
repository root with Pandoc, pdfLaTeX and Poppler installed. The checked-in PDF
was compiled and visually reviewed locally; the certificate workflow does not
claim a TeX environment or independent proof-assistant verification.

Written proofs establish the stated mathematical implications. Exact finite
checks certify their declared implementations. Physical calibration, independent
review and empirical validation remain separate evidence.
