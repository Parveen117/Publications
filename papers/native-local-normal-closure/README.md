# LNC-1 — Classical local normal-bracket completion

A selected canonical coframe-plus-scalar action supplies LCF's missing normal geometry transport. The full classical 1+1 secondary constraints, including its massive matter coupling, satisfy:

`{D[f],D[h]}=D[fh'-hf']`, `{D[f],H[N]}=H[fN'-Nf']`, `{H[N],H[M]}=D[NM'-MN']`.

The source `-NJ` in coframe momentum evolution and the local cut variation `-(NJ)'/2` are retained. Closure requires specific derivative coefficients within the declared ansatz but leaves the scalar potential arbitrary. Two exactly constrained initial configurations with different potentials have coframe accelerations `0` and `-1`, exposing that remaining model choice.

Replay: **66 exact checks, 4 native algebraic certificates/replays and 14 negative controls**.

- [Action, bracket derivation and scope](THEOREM.md)
- [Native certificate](CERTIFICATE.json)
- [Source pins](SOURCE_PINS.json)
- Replay: `python3 papers/native-local-normal-closure/run.py --check`

The canonical extension is of the known two-dimensional dilaton-gravity type, explicitly acknowledged. Additional fields and action selection are inputs. Physical reduction, positivity, quantization/anomalies, global topology and 3+1 gravity remain separate obligations.
