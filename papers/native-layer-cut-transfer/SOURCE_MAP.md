# LCT-1 source map

Every retained input, with its exact role. Files are pinned by SHA-256 in
[SOURCE_PINS.json](SOURCE_PINS.json) at base commit
167eee544a2c5cfb5cef5b8032c82ac26a4b2e0d.

| Source | Used for | Not used for |
|---|---|---|
| [ILC-1](../native-infinite-layer-cut/THEOREM.md) (PR #81) | layers c_n, forms H_N, A_N, B_N; positive local model f_N (ILC-T4); RH criterion and identity-theorem completeness (ILC-T2) | any sign of the actual all-depth forms |
| [LCW-1](../native-layer-continuation-window/THEOREM.md) (PR #82) | the next pair (c_{2N}, c_{2N+1}) as the data of the one-step transfer | its window bounds |
| [LAT-1](../native-layer-anchor-transport/THEOREM.md) (PR #83) | anchor hierarchy (LAT-T1) and Moebius transport of the finite model with invariant poles and residues (LAT-T2) | interval certificates |
| [LAA-1](../native-layer-anchor-atlas/THEOREM.md) (PR #84) | separately fitted models at different anchors differ | |
| [LDC-1](../native-layer-depth-compatibility/THEOREM.md) (PR #85) | the monic norm h_N and boundary eigenvalue h_N/h_{N-1} | its record-process limits |
| [LOT-1](../native-layer-observer-transport/THEOREM.md) (PR #86) | continuation position p_N = h_N/W_N | |
| [Finite observer detector](../../RH_FINITE_OBSERVER_DETECTOR.md) (PR #96) | D3: RH iff det H_N > 0 for every N | |
| [CLA-1](../../RH_CLOSURE_AUDIT.md) (PR #112) | source-independent structure gives no sign; finite prefixes are matched by positive and indefinite sources | |
| [RLO-1](../../RH_RELATIONAL_OBSERVER_CONFLUENCE.md) (PR #118) | source contract (1), typed roles | |
| [CCI-1](../../RH_CONJUGATE_CANCELLATION_INDEX.md) (PR #121) | real readings can still carry negative directions | |
| [RCD-1](../../RH_RESIDUAL_COMPLETION.md) (PR #124) | separated-node residual factorization; LCT-T2 is its confluent, node-free form | its filter-positivity range |
| [FWT-1](../../RH_FAITHFUL_WAVE_TRANSPORT.md) (PR #125) | exact recovery of q from complete data | |
| [SPT-1](../../RH_SHUNYA_POTENTIAL.md) (PR #111) | TVSP corners share only the contact one-form | any identification of potentials |
| JED-T5, branch claude/jet-entry-development at 72d6530 | smooth towers can be incomplete (flat deformation) | |
| native engine, papers/native-common-cone-limit/native_engine/native_operator.cjs | exact Q+iQ arithmetic | |

External classical background, credited and not reproved beyond what
THEOREM.md writes out:

- T. J. Stieltjes, *Recherches sur les fractions continues* (1894);
- G. A. Baker Jr. and P. Graves-Morris, *Padé Approximants*, 2nd ed., Cambridge (1996): Padé uniqueness, the Gauss–Christoffel remainder, continued fractions.
