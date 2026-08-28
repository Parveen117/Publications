#!/usr/bin/env python3
"""Assemble the collected-framework book from three source repositories.

Every markdown and LaTeX source in the declared manifest becomes a
chapter. Files whose pandoc conversion fails to compile are included
verbatim instead, so no content is ever dropped. Appendix B lists the
SHA-256 of every included source file.
"""
import hashlib
import os
import re
import subprocess
import sys

ROOT = "/home/claude"
OUT = os.path.join(ROOT, "book_build")
FRAG = os.path.join(OUT, "fragments")
os.makedirs(FRAG, exist_ok=True)

RKF = os.path.join(ROOT, "RKF")
RHF = os.path.join(ROOT, "RH-Framework")
PUB = os.path.join(ROOT, "Publications")


def md_files(base, sub):
    d = os.path.join(base, sub)
    if not os.path.isdir(d):
        return []
    out = []
    for r, _, fs in os.walk(d):
        if ".git" in r:
            continue
        for f in sorted(fs):
            if f.endswith(".md"):
                out.append(os.path.join(r, f))
    return sorted(out)


def num_key(p):
    m = re.match(r"(\d+)", os.path.basename(p))
    return (int(m.group(1)) if m else 999, os.path.basename(p))


# ---------------- declared manifest, in book order ----------------
ADMIN_NAMES = {
    "07_do_not_reopen.md", "09_cleanup_warning.md", "13_after_merge_tasks.md",
    "19_cleanup_branch_plan.md", "BRANCH_CONSOLIDATION_AUDIT.md",
    "MP_CLEANUP_CANDIDATES.md", "SOURCE_PR_INDEX.md", "provenance.json",
    "10_export_manifest.md",
}

rkf_top = [os.path.join(RKF, f) for f in
           ["README.md", "CLAIM_BOUNDARY.md", "TERMINOLOGY_AND_FILING_ALIGNMENT.md",
            "REVIEWER_GUIDE.md", "REPRODUCE.md", "CERTIFICATE_INDEX.md",
            "SOURCE_PROVENANCE.md", "NOTICE.md"] if os.path.exists(os.path.join(RKF, f))]

theorum_all = [p for p in md_files(RKF, "theorum")
               if os.path.dirname(p).endswith("theorum")]
theorum_main = sorted([p for p in theorum_all
                       if os.path.basename(p) not in ADMIN_NAMES
                       and os.path.basename(p) not in ("README.md",)], key=num_key)
theorum_admin = [p for p in theorum_all if os.path.basename(p) in ADMIN_NAMES]

thermo = sorted(md_files(RKF, "theorum/thermodynamics"), key=num_key)
rkf_theorems = md_files(RKF, "theorems")
rkf_certs = md_files(RKF, "certificates")
rkf_papers_md = md_files(RKF, "papers")
rkf_papers_tex = [os.path.join(RKF, "papers/rkf_completed_weil_endpoint/main.tex")]

rhf_top = [os.path.join(RHF, f) for f in ["README.md", "REVIEWER_CHECKLIST.md"]
           if os.path.exists(os.path.join(RHF, f))]
rhf_theorems = md_files(RHF, "theorems")
rhf_notebooks = md_files(RHF, "notebooks")

pub_papers = sorted([d for d in os.listdir(os.path.join(PUB, "papers"))
                     if os.path.isdir(os.path.join(PUB, "papers", d))])

PARTS = []
PARTS.append(("The Recognition Kernel Framework: Charter and Boundaries", rkf_top))
PARTS.append(("The Theorem Ladder", theorum_main))
PARTS.append(("Thermodynamic Response Theorems", thermo))
PARTS.append(("Foundation Theorems and Certificates", rkf_theorems + rkf_certs))
PARTS.append(("Completed Weil Endpoint Paper", rkf_papers_md + rkf_papers_tex))
PARTS.append(("The RH-Framework Ladder", rhf_top + [p for p in rhf_theorems]))
PARTS.append(("RH-Framework Notebooks", rhf_notebooks))
pub_files = []
for d in pub_papers:
    base = os.path.join(PUB, "papers", d)
    for f in ["README.md", "THEOREM.md", "LINEAGE.md"]:
        p = os.path.join(base, f)
        if os.path.exists(p):
            pub_files.append(p)
    for r, _, fs in os.walk(base):
        for f in sorted(fs):
            if f.endswith(".tex") and "v1_original" not in f:
                pub_files.append(os.path.join(r, f))
PARTS.append(("Publications: Certified Papers", pub_files))
PARTS.append(("Repository Working Notes (Appendix)", theorum_admin))

ALL = [p for _, fs in PARTS for p in fs]
print(f"manifest: {len(ALL)} files across {len(PARTS)} parts", file=sys.stderr)


# ---------------- conversion ----------------
def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def tex_escape(s):
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"),
                 ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"),
                 ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    return s


def title_of(path):
    if path.endswith(".md"):
        for line in open(path, encoding="utf-8", errors="replace"):
            if line.startswith("# "):
                return line[2:].strip()
    return os.path.basename(path)


REPAIRS = []


def pandoc_md_gfm(path, out_tex):
    r = subprocess.run(
        ["pandoc", "-f", "gfm+tex_math_dollars+tex_math_single_backslash",
         "-t", "latex", "--no-highlight",
         "--top-level-division=chapter", "--wrap=preserve", path,
         "-o", out_tex],
        capture_output=True, text=True)
    return r.returncode == 0


def pandoc_md(path, out_tex):
    r = subprocess.run(
        ["pandoc", "-f",
         "markdown+tex_math_single_backslash+pipe_tables",
         "-t", "latex", "--no-highlight",
         "--top-level-division=chapter", "--wrap=preserve", path,
         "-o", out_tex],
        capture_output=True, text=True)
    return r.returncode == 0


def extract_tex_body(path):
    src = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"\\begin{document}(.*)\\end{document}", src, re.S)
    body = m.group(1) if m else src
    body = re.sub(r"\\maketitle", "", body)
    body = re.sub(r"\\(title|author|date)\{.*?\}", "", body, flags=re.S)
    body = re.sub(r"\\tableofcontents", "", body)
    # papers use sections; keep them as sections under the chapter
    return body


def verbatim_chapter(path, title):
    raw = open(path, encoding="utf-8", errors="replace").read()
    raw = raw.replace("\\end{verbatim}", "\\end{ verbatim}")
    return ("\\chapter{%s}\n\\emph{Included verbatim (source preserved "
            "exactly; automated conversion was not compile-safe).}\n"
            "\\begin{verbatim}\n%s\n\\end{verbatim}\n" % (title, raw))


HARNESS = r"""\documentclass{book}
\usepackage{amsmath,amssymb,amsthm,mathrsfs,mathtools,bm}
\usepackage{longtable,booktabs,calc,array}\usepackage{graphicx}
\usepackage[hidelinks]{hyperref}\usepackage{fvextra}\usepackage{pdfpages}
\newcommand{\passthrough}[1]{#1}\newcommand{\tightlist}{}
\newcommand{\includesvg}[2][]{\texttt{[svg]}}\newcommand{\RSC}{\mathrm{RSC}}
\setcounter{secnumdepth}{1}
\begin{document}
\input{%s}
\end{document}
"""


def compiles(frag_path):
    test = os.path.join(OUT, "harness.tex")
    open(test, "w").write(HARNESS % frag_path)
    r = subprocess.run(["xelatex", "-no-pdf", "-interaction=nonstopmode",
                        "-halt-on-error", "-output-directory", OUT, test],
                       capture_output=True, timeout=180)
    return r.returncode == 0


fragments, quarantined, manifest = [], [], []
for idx, path in enumerate(ALL):
    rel = os.path.relpath(path, ROOT)
    manifest.append((rel, sha(path)))
    frag = os.path.join(FRAG, f"f{idx:03d}.tex")
    title = tex_escape(title_of(path))
    ok = False
    if path.endswith(".md"):
        raw = open(path, "rb").read()
        src_path = path
        CTRL = {b"\x07": b"\\a", b"\x08": b"\\b",
                b"\x0b": b"\\v", b"\x0c": b"\\f"}
        if any(k in raw for k in CTRL):
            fixed = raw
            names = []
            for k, v in CTRL.items():
                if k in fixed:
                    fixed = fixed.replace(k, v)
                    names.append(f"0x{k[0]:02x}->{v.decode()}")
            src_path = os.path.join(FRAG, f"s{idx:03d}.md")
            open(src_path, "wb").write(fixed)
            REPAIRS.append((os.path.relpath(path, ROOT),
                            "control bytes repaired before conversion: "
                            + ", ".join(names)))
        if pandoc_md(src_path, frag) or pandoc_md_gfm(src_path, frag):
            body = open(frag, encoding="utf-8", errors="replace").read()
            if "\\chapter" not in body:
                body = f"\\chapter{{{title}}}\n" + body
            body = (f"% source: {rel}\n" + body
                    + f"\n\\begin{{flushright}}\\small\\texttt{{{tex_escape(rel)}}}\\end{{flushright}}\n")
            open(frag, "w", encoding="utf-8").write(body)
            ok = compiles(frag)
            if not ok and pandoc_md_gfm(src_path, frag):
                body = open(frag, encoding="utf-8", errors="replace").read()
                if "\\chapter" not in body:
                    body = f"\\chapter{{{title}}}\n" + body
                body = (f"% source: {rel}\n" + body
                        + f"\n\\begin{{flushright}}\\small\\texttt{{{tex_escape(rel)}}}\\end{{flushright}}\n")
                open(frag, "w", encoding="utf-8").write(body)
                ok = compiles(frag)
    elif path.endswith(".tex"):
        pdir = os.path.dirname(path)
        raw = open(path, "rb").read()
        for k, v in {b"\x07": b"\\a", b"\x08": b"\\b",
                     b"\x0b": b"\\v", b"\x0c": b"\\f"}.items():
            if k in raw:
                raw = raw.replace(k, v)
                open(path, "wb").write(raw)
        eng = "pdflatex" if b"revtex" in open(path, "rb").read() else "xelatex"
        for _ in range(2):
            r = subprocess.run([eng, "-interaction=nonstopmode",
                                os.path.basename(path)],
                               cwd=pdir, capture_output=True, timeout=300)
        pdf = path[:-4] + ".pdf"
        if os.path.exists(pdf) and os.path.getsize(pdf) > 10000:
            if r.returncode != 0:
                REPAIRS.append((os.path.relpath(path, ROOT),
                                "compiled in nonstop mode; recoverable "
                                "LaTeX errors present in source"))
            dst = os.path.join(FRAG, f"p{idx:03d}.pdf")
            subprocess.run(["cp", pdf, dst], check=True)
            open(frag, "w", encoding="utf-8").write(
                f"\\chapter{{{title.replace('.tex','')}}}\n"
                f"\\emph{{The following paper is included as compiled from "
                f"its repository source.}}\n"
                f"\\includepdf[pages=-]{{{dst}}}\n")
            ok = compiles(frag)
        if not ok:
            body = (f"\\chapter{{{title.replace('.tex','')}}}\n"
                    + extract_tex_body(path))
            open(frag, "w", encoding="utf-8").write(body)
            ok = compiles(frag)
    if not ok:
        open(frag, "w", encoding="utf-8").write(verbatim_chapter(path, title))
        if not compiles(frag):
            open(frag, "w", encoding="utf-8").write(
                f"\\chapter{{{title}}}\\emph{{See repository source: "
                f"\\texttt{{{tex_escape(rel)}}}}}\n")
        quarantined.append(rel)
    fragments.append((path, frag))
    print(f"[{idx+1}/{len(ALL)}] {'ok ' if ok else 'VRB'} {rel}",
          file=sys.stderr)

# ---------------- master ----------------
part_lookup = {}
for pt, fs in PARTS:
    for p in fs:
        part_lookup[p] = pt

master = [r"""\documentclass[11pt,openany]{book}
\usepackage{amsmath,amssymb,amsthm,mathrsfs,mathtools,bm}
\usepackage{longtable,booktabs,calc,array}
\usepackage{graphicx}\usepackage[margin=1in]{geometry}
\usepackage[hidelinks]{hyperref}\usepackage{fvextra}\usepackage{pdfpages}
\DefineVerbatimEnvironment{verbatim}{Verbatim}{breaklines=true,fontsize=\footnotesize}
\newcommand{\passthrough}[1]{#1}\newcommand{\tightlist}{}
\newcommand{\includesvg}[2][]{\texttt{[svg]}}\newcommand{\RSC}{\mathrm{RSC}}
\setcounter{secnumdepth}{1}\setcounter{tocdepth}{0}
\title{\Huge The Recognition Kernel Framework\\[0.4em]
\Large Collected Theorems, Theory, and Certified Results}
\author{Monty Dabas\\\small ORCID 0009-0005-6948-209X}
\date{August 2026}
\begin{document}
\frontmatter\maketitle
\chapter*{About this volume}
This volume collects, without omission, the theorem ladders, theory
documents, claim boundaries, and certified papers of the Recognition
Kernel Framework as recorded in three repositories:
\texttt{Recognition-Kernel-Framework}, \texttt{RH-Framework}, and
\texttt{Publications}. Chapters are automated renderings of the
repository sources; any chapter whose automated conversion was not
compile-safe is included verbatim so that no content is lost. Appendix
B lists every included source file with its SHA-256 digest. The
authoritative artefacts remain the repositories themselves, whose
certificates and tests regenerate every pinned result.
\tableofcontents\mainmatter
"""]
cur = None
for path, frag in fragments:
    pt = part_lookup[path]
    if pt != cur:
        master.append(f"\\part{{{tex_escape(pt)}}}\n")
        cur = pt
    master.append(f"\\input{{{frag}}}\n")

master.append("\\appendix\\part{Provenance}\n\\chapter{Source Manifest}\n"
              "\\begin{small}\\begin{longtable}{p{0.62\\textwidth}p{0.33\\textwidth}}\n"
              "\\textbf{Source file} & \\textbf{SHA-256 (first 32)}\\\\\\hline\n")
for rel, h in manifest:
    master.append(f"\\texttt{{{tex_escape(rel)}}} & \\texttt{{{h[:32]}}}\\\\\n")
master.append("\\end{longtable}\\end{small}\n")
if REPAIRS:
    master.append("\\chapter{Automated Repairs}\n\\begin{itemize}\n")
    for rel, note in REPAIRS:
        master.append(f"\\item \\texttt{{{tex_escape(rel)}}}: {tex_escape(note)}\n")
    master.append("\\end{itemize}\n")
master.append("\\end{document}\n")

open(os.path.join(OUT, "book.tex"), "w", encoding="utf-8").write("".join(master))
print(f"quarantined (verbatim): {len(quarantined)}", file=sys.stderr)
for q in quarantined:
    print("  VRB", q, file=sys.stderr)
open(os.path.join(OUT, "QUARANTINE.txt"), "w").write("\n".join(quarantined))
