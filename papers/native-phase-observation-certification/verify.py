#!/usr/bin/env python3
"""Python 3.12 orchestration; all mathematical execution is native Node RKF.

Default --check is read-only. --write explicitly generates only this new packet's
certificate/checksum; it never refreshes source pins or inherited results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PACKETS = (
    ("Clock_Free_Native_Lattice", "CELL_RESULT.json"),
    ("Clock_Free_Cut_Quantisation", "QUANTISATION_RESULT.json"),
    ("Clock_Free_Physics_Cubits", "SCALE_CUBIT_RESULT.json"),
    ("Clock_Free_Scale_Compass", "SCALE_COMPASS_RESULT.json"),
    ("Clock_Relation_Reconstruction", "CLOCK_RELATION_RESULT.json"),
    ("Native_Phase_Word_Reconstruction", "PHASE_WORD_RESULT.json"),
)
CORE = Path("parent/parent/parent/parent/baseline/native/core")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def verify_digest(data: bytes, expected: str, label: str) -> None:
    require(sha256(data) == expected, f"Source pin mismatch: {label}")


def verify_sources(pins: dict) -> None:
    for record in pins["records"]:
        rel = PurePosixPath(record["path"])
        require(not rel.is_absolute() and ".." not in rel.parts, "Invalid pin path")
        source = REPO.joinpath(*rel.parts)
        require(source.is_file() and not source.is_symlink(), f"Missing regular source: {rel}")
        data = source.read_bytes()
        verify_digest(data, record["sha256"], str(rel))
        require(len(data) == record["bytes"], f"Source size mismatch: {rel}")


def extract_packet(archive: Path, target: Path, expected_root: str) -> None:
    with zipfile.ZipFile(archive) as handle:
        seen = set()
        for entry in handle.infolist():
            rel = PurePosixPath(entry.filename)
            require(not rel.is_absolute() and ".." not in rel.parts and "\\" not in entry.filename, "Unsafe archive member")
            require(rel.parts and rel.parts[0] == expected_root, "Unexpected archive root")
            require(entry.filename not in seen, "Duplicate archive member")
            require(not stat.S_ISLNK(entry.external_attr >> 16), "Symlink archive member")
            seen.add(entry.filename)
        handle.extractall(target)


def snapshot(root: Path) -> dict[str, str]:
    return {str(p.relative_to(root)): sha256(p.read_bytes()) for p in sorted(root.rglob("*")) if p.is_file()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Read-only source validation and deterministic replay (default)")
    mode.add_argument("--write", action="store_true", help="Explicitly generate the new frozen certificate and checksum")
    args = parser.parse_args()
    require(sys.version_info[:2] == (3, 12), "Use Python 3.12 only")
    require(not sys.flags.optimize, "Python assertion optimization is forbidden for inherited checks")
    node = subprocess.run(["node", "--version"], text=True, capture_output=True, check=True).stdout.strip()
    require(node == "v24.19.0", "Use pinned Node v24.19.0")
    pins = json.loads((HERE / "SOURCE_PINS.json").read_text())
    verify_sources(pins)
    # A mutation must fail the same validator used on real source inputs.
    first = pins["records"][0]
    try:
        verify_digest((REPO / first["path"]).read_bytes() + b"\nchanged\n", first["sha256"], "intentional negative control")
    except RuntimeError:
        pass
    else:
        raise RuntimeError("Source mutation was incorrectly accepted")

    report = {
        "schema": "native-inherited-packet-replay.v1",
        "source_commit": pins["source_commit"],
        "runtime": {"python": "3.12", "node": node},
        "source_mutation_rejected": True,
        "packets": [],
    }
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.pop("PYTHONOPTIMIZE", None)
    with tempfile.TemporaryDirectory(prefix="native-phase-certification-") as scratch:
        work = Path(scratch)
        for name, result_name in PACKETS:
            archive = REPO / "papers" / (name + "_Source.zip")
            packet = work / name
            extract_packet(archive, work, name)
            require((packet / "verify.py").is_file(), f"Missing packet verifier: {name}")
            before = snapshot(packet)
            run = subprocess.run([sys.executable, "verify.py"], cwd=packet, env=env, text=True, capture_output=True, timeout=180)
            if run.returncode:
                sys.stderr.write(run.stdout + run.stderr)
            require(run.returncode == 0, f"Inherited packet failed: {name}")
            require(snapshot(packet) == before, f"Inherited --check modified packet contents: {name}")
            result_path = packet / result_name
            result = json.loads(result_path.read_text())
            require(result["status"].startswith("PASS"), f"Unexpected recorded status: {name}")
            pdf = REPO / "papers" / (name + ".pdf")
            pdf_match = None
            if pdf.exists():
                pdf_match = pdf.read_bytes() == (packet / pdf.name).read_bytes()
                require(pdf_match, f"Uploaded PDF differs from its source packet: {name}")
            item = {
                "name": name,
                "archive": str(archive.relative_to(REPO)),
                "archive_sha256": sha256(archive.read_bytes()),
                "command": "python3.12 verify.py",
                "exit_code": run.returncode,
                "status": result["status"],
                "summary": result["summary"],
                "result_file": result_name,
                "result_sha256": sha256(result_path.read_bytes()),
                "stdout_sha256": sha256(run.stdout.encode()),
                "stderr_sha256": sha256(run.stderr.encode()),
                "extracted_files_unchanged": True,
                "uploaded_pdf_matches_packet": pdf_match,
            }
            report["packets"].append(item)
            print(f"PASS inherited: {name}", flush=True)
        report_path = work / "UPSTREAM_RUN.json"
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        native_root = work / "Native_Phase_Word_Reconstruction" / CORE
        for filename, expected in pins["native_engine"]["files"].items():
            verify_digest((native_root / filename).read_bytes(), expected, filename)
        subprocess.run([
            "node", str(HERE / "verify.cjs"), "--native-root", str(native_root),
            "--upstream-report", str(report_path), "--write" if args.write else "--check",
        ], cwd=REPO, env=env, check=True, timeout=240)
    verify_sources(pins)
    print("Source pins, five uploaded PDF matches, six read-only packet reruns and new native certification passed.")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, subprocess.SubprocessError, OSError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
