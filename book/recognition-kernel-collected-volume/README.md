# The Recognition Kernel Framework — Collected Volume

> **Edition notice, 2026-09-25:** This is a historical collected edition.
> Read it with [the correction record](../../papers/corrections-2026-09-25/README.md),
> which supersedes the identified formulas and claim scopes. The PDF and
> generated fragments have not been silently relabelled as a corrected edition.

Single-volume LaTeX compilation of the complete theorem ladders, theory
documents, claim boundaries, and certified papers from three
repositories: `Recognition-Kernel-Framework`, `RH-Framework`, and
`Publications`. 128 source files, nine parts, 600+ pages, nothing
omitted.

- `recognition_kernel_collected_volume.pdf` — the compiled book.
- `book.tex` + `fragments/` — the full LaTeX source.
- `assemble_volume.py` — the assembler that regenerates the volume from
  the three repositories (pandoc + xelatex). Every included source file
  is listed with its SHA-256 in Appendix "Source Manifest"; automated
  byte repairs (control characters that had replaced `\b`, `\f` in a few
  sources) are listed in "Automated Repairs".

The authoritative artefacts remain the repositories themselves; this
volume is a faithful, provenance-tracked rendering.
