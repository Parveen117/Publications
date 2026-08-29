"""AZQ1-I: GHZ parity-chart audit on real IBM hardware (preregistered)."""
import hashlib
import json
import os

from qiskit import QuantumCircuit, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2

TOKEN = os.environ["IBM_TOKEN"]
svc = None
for ch in ("ibm_quantum_platform", "ibm_cloud", "ibm_quantum"):
    try:
        svc = QiskitRuntimeService(channel=ch, token=TOKEN)
        _ = svc.backends()
        print("channel:", ch)
        break
    except Exception as e:
        print("channel", ch, "failed:", str(e)[:100])
if svc is None:
    raise SystemExit("no channel worked")

be = svc.least_busy(operational=True, simulator=False)
print("backend:", be.name)
SHOTS, P0 = 512, 0.05


def ghz(corrupt=False, xbasis=False):
    qc = QuantumCircuit(3, 3)
    qc.h(0)
    qc.cx(0, 1)
    if corrupt:
        qc.x(1)
    qc.cx(1, 2)
    if xbasis:
        for q in range(3):
            qc.h(q)
    qc.measure(range(3), range(3))
    return qc


def counts_of(res0):
    d = res0.data
    reg = getattr(d, "c", None) or list(vars(d).values())[0]
    return reg.get_counts()


def forbidden(counts, xbasis):
    tot = sum(counts.values())
    bad = sum(c for b, c in counts.items()
              if ((b.replace(" ", "").count("1") % 2 == 1) if xbasis
                  else b.replace(" ", "") not in ("000", "111")))
    return bad / tot


sampler = SamplerV2(mode=be)
cert = {"backend": be.name, "shots": SHOTS,
        "preregistration": "PREREGISTRATION_I.md", "rungs": {}}
for label, corrupt in (("clean", False), ("corrupted", True)):
    Ws = {}
    for basis in ("Z", "X"):
        qc = transpile(ghz(corrupt, basis == "X"), be, optimization_level=1)
        job = sampler.run([qc], shots=SHOTS)
        Ws[basis] = forbidden(counts_of(job.result()[0]), basis == "X")
    W = (Ws["Z"] + Ws["X"]) / 2
    se = max((P0 * (1 - P0) / (2 * SHOTS)) ** 0.5, 1e-6)
    z = (W - P0) / se
    v = "NOT_FALSIFIED" if z < 5 else "FALSIFIED_MEASUREMENT_CONTRACT"
    cert["rungs"][label] = {"W_Z": Ws["Z"], "W_X": Ws["X"], "W": W,
                            "p0": P0, "z": z, "verdict": v}
    print(f"{be.name} {label}: W={W:.4f} z={z:+.1f} -> {v}")

body = json.dumps({k: v for k, v in cert.items() if k != "_pin"},
                  sort_keys=True)
cert["_pin"] = hashlib.sha256(body.encode()).hexdigest()
open("AZQ1_I_CERTIFICATE.json", "w").write(
    json.dumps(cert, indent=2, sort_keys=True))
print("written AZQ1_I_CERTIFICATE.json — pin:", cert["_pin"])
