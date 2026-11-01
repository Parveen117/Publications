#!/usr/bin/env python3
"""Python 3.12 orchestration for CID-D2; native RKF performs all mathematics."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
PREVIOUS = REPO / "papers/native-loop-concentration"
PHASE = REPO / "papers/native-phase-observation-certification"
PACKET = "Native_Phase_Word_Reconstruction"
CORE = Path("parent/parent/parent/parent/baseline/native/core")


def require(ok: bool, message: str) -> None:
    if not ok:
        raise RuntimeError(message)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sources(pins: dict) -> None:
    for item in pins["records"]:
        path = REPO / item["path"]
        require(path.resolve().is_relative_to(REPO), "Pin path escapes repository")
        require(path.is_file() and not path.is_symlink(), f"Missing source: {item['path']}")
        data = path.read_bytes()
        require(len(data) == item["bytes"] and sha256(data) == item["sha256"], f"Source pin mismatch: {item['path']}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Read-only deterministic replay (default)")
    mode.add_argument("--write", action="store_true", help="Explicitly generate only CID-D2 certificate/checksum")
    args = parser.parse_args()
    require(sys.version_info[:2] == (3, 12), "Python 3.12 only")
    require(not sys.flags.optimize, "Inherited checks require assertions enabled")
    require(subprocess.check_output(["node", "--version"], text=True).strip() == "v24.19.0", "Node v24.19.0 required")
    pins = json.loads((HERE / "SOURCE_PINS.json").read_text())
    sources(pins)
    sys.dont_write_bytecode = True
    spec = importlib.util.spec_from_file_location("pinned_phase_orchestration", PHASE / "verify.py")
    require(spec is not None and spec.loader is not None, "Missing pinned upstream orchestration")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.pop("PYTHONOPTIMIZE", None)
    inherited = subprocess.run([sys.executable, str(PREVIOUS / "verify.py"), "--check"], cwd=REPO, env=env, text=True, capture_output=True, timeout=480)
    if inherited.returncode:
        sys.stderr.write(inherited.stdout + inherited.stderr)
    require(inherited.returncode == 0, "Prior cut-information, loop-concentration and inherited replay failed")
    previous_data = (PREVIOUS / "CERTIFICATE.json").read_bytes()
    previous_result = json.loads(previous_data)
    upstream = {"status": "PASS_PRIOR_READ_ONLY_REPLAY", "certificate_sha256": sha256(previous_data),
                "summary": previous_result["summary"], "stdout_sha256": sha256(inherited.stdout.encode()),
                "stderr_sha256": sha256(inherited.stderr.encode()), "exit_code": inherited.returncode}
    print("PASS prior cut-information, loop-concentration and complete inherited replay", flush=True)
    with tempfile.TemporaryDirectory(prefix="native-cut-transport-feedback-") as temp:
        work = Path(temp)
        prior.extract_packet(REPO / "papers" / (PACKET + "_Source.zip"), work, PACKET)
        core = work / PACKET / CORE
        for filename, expected in pins["native_engine"]["files"].items():
            require(sha256((core / filename).read_bytes()) == expected, "Native engine pin mismatch: " + filename)
        report = work / "UPSTREAM_RUN.json"
        report.write_text(json.dumps(upstream, indent=2) + "\n")
        subprocess.run(["node", str(HERE / "verify.cjs"), "--native-root", str(core),
                        "--upstream-report", str(report), "--write" if args.write else "--check"],
                       cwd=REPO, env=env, check=True, timeout=300)
    sources(pins)
    print("CID-D2 source integrity, native certification and inherited regression passed.")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
