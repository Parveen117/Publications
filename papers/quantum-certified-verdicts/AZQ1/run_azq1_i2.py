"""AZQ1-I2: GHZ parity audit on real IBM hardware, declared error
suppression (preregistered in PREREGISTRATION_I2.md)."""
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
        print("channel", ch, "failed:", str(e)[:80])
be = svc.least_busy(operational=True, simulator=False)
print("backend:", be.name)
tgt = be.target
SHOTS, P0 = 512, 0.05

# ---- best connected chain a-b-c from live calibration ----
twoq = {}
for gname in ("cz", "ecr", "cx"):
    if gname in tgt.operation_names:
        for qargs, props in tgt[gname].items():
            if props and props.error is not None and len(qargs) == 2:
                twoq[tuple(qargs)] = min(twoq.get(tuple(qargs), 1.0),
                                         props.error)
ro = {q[0]: p.error for q, p in tgt["measure"].items()
      if p and p.error is not None}
adj = {}
for (a, b) in twoq:
    adj.setdefault(a, set()).add(b)
    adj.setdefault(b, set()).add(a)


def pair_err(a, b):
    return twoq.get((a, b), twoq.get((b, a), 1.0))


best, bscore = None, 9e9
for b in adj:
    for a in adj[b]:
        for c in adj[b]:
            if c in (a, b) or a == b:
                continue
            s = (pair_err(a, b) + pair_err(b, c)
                 + ro.get(a, .05) + ro.get(b, .05) + ro.get(c, .05))
            if s < bscore:
                bscore, best = s, (a, b, c)
a, b, c = best
print(f"chain {best} score {bscore:.4f}")


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
    bad = sum(v for k, v in counts.items()
              if ((k.replace(" ", "").count("1") % 2 == 1) if xbasis
                  else k.replace(" ", "") not in ("000", "111")))
    return bad / tot


sampler = SamplerV2(mode=be)
sampler.options.dynamical_decoupling.enable = True
sampler.options.dynamical_decoupling.sequence_type = "XY4"
sampler.options.twirling.enable_gates = True
sampler.options.twirling.enable_measure = True

# descriptive readout-calibration control
cal = {}
for lbl, flip in (("000", False), ("111", True)):
    qc = QuantumCircuit(3, 3)
    if flip:
        qc.x([0, 1, 2])
    qc.measure(range(3), range(3))
    tq = transpile(qc, be, initial_layout=[a, b, c], optimization_level=1)
    cnt = counts_of(sampler.run([tq], shots=256).result()[0])
    good = cnt.get(lbl, 0) / sum(cnt.values())
    cal[lbl] = {"correct_fraction": good}
    print(f"cal {lbl}: correct {good:.3f}")

cert = {"backend": be.name, "chain": list(best), "chain_score": bscore,
        "shots": SHOTS, "preregistration": "PREREGISTRATION_I2.md",
        "readout_calibration": cal, "rungs": {}}
for label, corrupt in (("clean", False), ("corrupted", True)):
    Ws = {}
    for basis in ("Z", "X"):
        qc = transpile(ghz(corrupt, basis == "X"), be,
                       initial_layout=[a, b, c], optimization_level=3)
        Ws[basis] = forbidden(
            counts_of(sampler.run([qc], shots=SHOTS).result()[0]),
            basis == "X")
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
open("AZQ1_I2_CERTIFICATE.json", "w").write(
    json.dumps(cert, indent=2, sort_keys=True))
print("written AZQ1_I2_CERTIFICATE.json — pin:", cert["_pin"])
