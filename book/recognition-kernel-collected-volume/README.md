# The Recognition Kernel Framework — Collected Volume

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
