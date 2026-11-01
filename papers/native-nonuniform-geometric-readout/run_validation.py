"""Python 3.12: process orchestration, source hashes and result I/O only."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError("Use Python 3.12 for this source packet")
    pins = json.loads((HERE / "SOURCE_PINS.json").read_text())
    for relative, pin in pins["sources"].items():
        data = (ROOT / relative).read_bytes()
        assert hashlib.sha256(data).hexdigest() == pin["sha256"], relative
        assert subprocess.check_output(
            ["git", "hash-object", "--stdin"], input=data, cwd=ROOT
        ).decode().strip() == pin["git_blob"], relative
        original = subprocess.check_output(
            ["git", "show", f'{pins["parent"]}:{relative}'], cwd=ROOT
        )
        assert original == data, relative
    for relative, digest in pins["new_executables"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest, relative
    result = subprocess.run(
        ["node", str(HERE / "verify_source.cjs")],
        check=True, text=True, capture_output=True, cwd=ROOT,
    )
    report = json.loads(result.stdout)
    assert report["status"] == "PASS"
    output = HERE / "VALIDATION.json"
    if args.write:
        output.write_text(result.stdout)
    else:
        assert json.loads(output.read_text()) == report, "Saved report differs from replay"
    print(json.dumps({
        "status": report["status"], "source_pins_verified": len(pins["sources"]),
        "exact_checks": report["checks"], "negative_controls": report["negative_controls"],
        "saved_report": "written" if args.write else "matches replay",
        "continuum_POVM_executed": False,
    }))


if __name__ == "__main__":
    main()
