# Thermodynamic response: corrected measurement contract

Repository correction edition of arXiv:2603.20773v3, 2026-09-25.

The source tree in source/ retains the manuscript and corrects the noise-stable
Pluecker test, abstract and conclusion. SOURCE_MANIFEST.json records the exact
original source hashes. The arXiv version remains a historical record.

A bracket reconstructed as u v^T-v u^T satisfies the Pluecker identity even
when the gradients are wrong. It cannot falsify the rank-two hypothesis.
The revised theorem takes six independently observed skew entries, with an
independently justified error budget. Its proof retains the original safe
perturbation bound. The gradient error estimate remains a valid reconstruction
bound and is explicitly separated from the measurement test.

See ../corrections-2026-09-25/ for exact positive and negative controls.
This edition does not assert empirical confirmation, nonlinear quantum gravity,
or wholesale mathematical certification of the manuscript.

Build from source/: latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

[Compiled correction edition](thermodynamic_response_corrected.pdf)
