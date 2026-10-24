"""Python 3.12 I/O only: typeset the manuscript or verify its provenance."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

if sys.version_info[:2] != (3, 12):
    raise SystemExit("Use Python 3.12")

parser = argparse.ArgumentParser(description=__doc__)
mode = parser.add_mutually_exclusive_group(required=True)
mode.add_argument("--build", action="store_true")
mode.add_argument("--check", action="store_true")
args = parser.parse_args()
here = Path(__file__).resolve().parent
repo = here.parents[2]

if args.check:
    manifest = json.loads((here / "MANUSCRIPT_MANIFEST.json").read_text())
    for relative, expected in manifest["sha256"].items():
        actual = hashlib.sha256((repo / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise SystemExit("Manuscript/source hash mismatch: " + relative)
    print(json.dumps({
        "status": "PASS_MANUSCRIPT_PROVENANCE_ONLY",
        "files": len(manifest["sha256"]),
        "evidence_commit": manifest["evidence_commit"],
        "general_proof_assistant_verification": False,
    }))
else:
    engine = shutil.which("pdflatex")
    if not engine:
        raise SystemExit("Install pdfLaTeX and the packages used by paper.tex")
    with tempfile.TemporaryDirectory(prefix="chd-manuscript-") as directory:
        output = Path(directory)
        for _ in range(2):
            run = subprocess.run(
                [engine, "-interaction=nonstopmode", "-halt-on-error",
                 "-output-directory=" + str(output), str(here / "paper.tex")],
                cwd=here, text=True, capture_output=True,
            )
            if run.returncode:
                raise SystemExit(run.stdout[-6000:] + run.stderr[-2000:])
        log = (output / "paper.log").read_text(errors="replace")
        if "undefined references" in log or "Overfull" in log:
            raise SystemExit("Resolve undefined references or overflow before delivery")
        shutil.copyfile(output / "paper.pdf", here / "paper.pdf")
    print("Built " + str(here / "paper.pdf"))
