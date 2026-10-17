# LCT-1 — layer-tower cut transfer

A cut at depth N on the layer tower at a reference anchor transfers the whole
tower exactly:

    f(λ) = R_N(λ)/Q_N(λ) + (λ^{2N}/Q_N(λ)^2) · E_N(λ),   E_N(0) = h_N = det H_{N+1}/det H_N.

- The factor λ^{2N}/Q_N^2 is pre-supplied by the retained readings.
- Everything above the cut is carried by the single bilinear pairing E_N.
- Moving the cut one layer has the pre-supplied factor Q_N/Q_{N+1}, from the
  next pair only.

The sign of E_N is not supplied. For the actual xi source, h_N > 0 for every
N is the existing RH criterion.

- Theorems: [THEOREM.md](THEOREM.md)
- Inputs: [SOURCE_MAP.md](SOURCE_MAP.md)
- Witnesses: [verify.cjs](verify.cjs), 389 exact checks and 65 negative
  controls over Q+iQ. These are regression evidence only.

Usage, with Node 22 and Python 3.12 for I/O:

    python3.12 run.py --check

Branch `layer-cut-transfer-20261017`, based on 167eee5 (PR #125 head). Not a
pull request. No predecessor packet changed.
