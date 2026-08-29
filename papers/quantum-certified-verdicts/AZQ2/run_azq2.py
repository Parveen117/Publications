"""AZQ2: cross-platform concordance certificate (preregistered)."""
import hashlib
import json
import math
import os

from azure.quantum import Workspace
from azure.quantum.qiskit import AzureQuantumProvider
from qiskit import QuantumCircuit

ws = Workspace.from_connection_string(os.environ["AZQ_CONNECTION"])
provider = AzureQuantumProvider(ws)
SHOTS, TAU, ZC = 512, 0.10, 3.0
A, B = "quantinuum.sim.h2-1e", "rigetti.sim.qvm"


def bell(xbasis=False, miscal=False):
    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    if miscal:
        qc.rz(0.60, 0)
    if xbasis:
        qc.h(0)
        qc.h(1)
    qc.measure([0, 1], [0, 1])
    return qc


def parity_E(counts):
    tot = sum(counts.values())
    s = sum(c * (1 if k.replace(" ", "").count("1") % 2 == 0 else -1)
            for k, c in counts.items())
    return s / tot


def run_E(backend_name, xbasis, miscal):
    be = provider.get_backend(backend_name)
    qc = bell(xbasis=xbasis, miscal=miscal)
    counts = be.run(qc, shots=SHOTS).result().get_counts()
    E = parity_E(counts)
    se = math.sqrt(max(1.0 - E * E, 1e-9) / SHOTS)
    return E, se


def concordance(rung, miscal_on_A):
    out = {}
    for obs, xb in (("E_ZZ", False), ("E_XX", True)):
        EA, sA = run_E(A, xb, miscal_on_A)
        EB, sB = run_E(B, xb, False)
        d = EA - EB
        sed = math.sqrt(sA * sA + sB * sB)
        z = (abs(d) - TAU) / sed
        v = "CONCORDANT" if z < ZC else "DISCORDANT"
        out[obs] = {"E_A": EA, "E_B": EB, "d": d, "se_d": sed,
                    "tau": TAU, "z": z, "verdict": v}
        print(f"{rung} {obs}: A={EA:+.3f} B={EB:+.3f} d={d:+.3f} "
              f"z={z:+.1f} -> {v}")
    return out


cert = {"platforms": [A, B], "shots": SHOTS, "tau": TAU,
        "preregistration": "PREREGISTRATION.md",
        "AZQ2A_CLEAN": concordance("AZQ2-A", False),
        "AZQ2B_MISCALIBRATED": concordance("AZQ2-B", True)}
body = json.dumps({k: v for k, v in cert.items() if k != "_pin"},
                  sort_keys=True)
cert["_pin"] = hashlib.sha256(body.encode()).hexdigest()
open("AZQ2_CERTIFICATE.json", "w").write(
    json.dumps(cert, indent=2, sort_keys=True))
print("written AZQ2_CERTIFICATE.json — pin:", cert["_pin"])
