"""AZQ1: GHZ parity-chart audit on Azure Quantum (preregistered).

Run inside Azure Cloud Shell:
  pip install --quiet azure-quantum qiskit qiskit-qir
  python run_azq1.py "<WORKSPACE_RESOURCE_ID>" eastus

Targets (free tiers): quantinuum.sim.h2-1e (noisy emulator),
rigetti.sim.qvm (ideal QVM). Prints the certificate JSON; commit it
back to this folder unchanged.
"""
import hashlib
import json
import os
import sys

from azure.quantum import Workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit

CONN = os.environ.get("AZQ_CONNECTION", "")
if CONN:
    ws = Workspace.from_connection_string(CONN)
else:
    RES_ID = sys.argv[1]
    LOC = sys.argv[2] if len(sys.argv) > 2 else "eastus"
    ws = Workspace(resource_id=RES_ID, location=LOC)
provider = AzureQuantumProvider(ws)
def _bname(b):
    n = getattr(b, "name", None)
    return n if isinstance(n, str) else n()


print("targets:", [_bname(b) for b in provider.backends()])

SHOTS = 512
P0 = {"quantinuum.sim.h2-1e": 0.05, "rigetti.sim.qvm": 0.01}


def ghz(n=3, corrupt=False, xbasis=False):
    qc = QuantumCircuit(n, n)
    qc.h(0)
    qc.cx(0, 1)
    if corrupt:
        qc.x(1)
    qc.cx(1, 2)
    if xbasis:
        for q in range(n):
            qc.h(q)
    qc.measure(range(n), range(n))
    return qc


def forbidden_weight(counts, xbasis):
    tot = sum(counts.values())
    bad = 0
    for bits, c in counts.items():
        b = bits.replace(" ", "")
        if xbasis:
            if b.count("1") % 2 == 1:
                bad += c
        else:
            if b not in ("000", "111"):
                bad += c
    return bad / tot


def run(backend_name):
    be = provider.get_backend(backend_name)
    out = {}
    for label, corrupt in (("clean", False), ("corrupted", True)):
        Ws = {}
        for basis in ("Z", "X"):
            qc = ghz(corrupt=corrupt, xbasis=(basis == "X"))
            job = be.run(qc, shots=SHOTS)
            counts = job.result().get_counts()
            Ws[basis] = forbidden_weight(counts, basis == "X")
        W = (Ws["Z"] + Ws["X"]) / 2
        p0 = P0.get(backend_name, 0.05)
        se = max((p0 * (1 - p0) / (2 * SHOTS)) ** 0.5, 1e-6)
        z = (W - p0) / se
        out[label] = {"W_Z": Ws["Z"], "W_X": Ws["X"], "W": W,
                      "p0": p0, "z": z,
                      "verdict": "NOT_FALSIFIED" if z < 5
                      else "FALSIFIED_MEASUREMENT_CONTRACT"}
        print(f"{backend_name} {label}: W={W:.4f} z={z:+.1f} "
              f"-> {out[label]['verdict']}")
    return out


cert = {"preregistration": "PREREGISTRATION.md", "shots": SHOTS,
        "backends": {}}
for name in ("quantinuum.sim.h2-1e", "rigetti.sim.qvm"):
    try:
        cert["backends"][name] = run(name)
    except Exception as e:
        cert["backends"][name] = {"status": "NOT_EXECUTABLE",
                                  "error": str(e)[:200]}
        print(name, "NOT_EXECUTABLE:", str(e)[:120])

body = json.dumps({k: v for k, v in cert.items() if k != "_pin"},
                  sort_keys=True)
cert["_pin"] = hashlib.sha256(body.encode()).hexdigest()
print(json.dumps(cert, indent=2, sort_keys=True))
open("AZQ1_CERTIFICATE.json", "w").write(
    json.dumps(cert, indent=2, sort_keys=True))
print("\nwritten AZQ1_CERTIFICATE.json — pin:", cert["_pin"])
