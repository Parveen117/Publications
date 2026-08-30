"""AZQ1-I3: GHZ parity audit on real IBM hardware with a
calibration-derived noise budget (preregistered in PREREGISTRATION_I3.md).

Order of operations is the point of this rung:
  1) select chain from live calibration,
  2) derive p0 from that same calibration snapshot,
  3) write + hash AZQ1_I3_PRECOMMIT.json,
  4) only then submit jobs,
  5) pin full per-bitstring histograms in the certificate.
"""
import datetime
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
SHOTS = 8192

# ---- best connected chain a-b-c from live calibration (same rule as I2) ----
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
for (x, y) in twoq:
    adj.setdefault(x, set()).add(y)
    adj.setdefault(y, set()).add(x)


def pair_err(x, y):
    return twoq.get((x, y), twoq.get((y, x), 1.0))


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

# ---- declared noise budget from the SAME calibration snapshot ----
cal_vals = {
    "e2_ab": pair_err(a, b), "e2_bc": pair_err(b, c),
    "ro_a": ro.get(a, .05), "ro_b": ro.get(b, .05), "ro_c": ro.get(c, .05),
}
P0 = max(sum(cal_vals.values()), 0.005)
precommit = {
    "preregistration": "PREREGISTRATION_I3.md",
    "backend": be.name, "chain": list(best),
    "calibration_values": cal_vals, "p0": P0, "shots": SHOTS,
    "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}
pc_body = json.dumps(precommit, sort_keys=True)
precommit["_pin"] = hashlib.sha256(pc_body.encode()).hexdigest()
open("AZQ1_I3_PRECOMMIT.json", "w").write(
    json.dumps(precommit, indent=2, sort_keys=True))
print("PRECOMMIT pinned BEFORE submission:", precommit["_pin"])
print(f"derived p0 = {P0:.5f}")


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

cert = {"backend": be.name, "chain": list(best), "chain_score": bscore,
        "shots": SHOTS, "preregistration": "PREREGISTRATION_I3.md",
        "precommit_pin": precommit["_pin"], "p0": P0,
        "calibration_values": cal_vals,
        "readout_calibration": {}, "histograms": {}, "rungs": {}}

# descriptive readout-calibration control (histograms pinned too)
for lbl, flip in (("000", False), ("111", True)):
    qc = QuantumCircuit(3, 3)
    if flip:
        qc.x([0, 1, 2])
    qc.measure(range(3), range(3))
    tq = transpile(qc, be, initial_layout=[a, b, c], optimization_level=1)
    cnt = counts_of(sampler.run([tq], shots=256).result()[0])
    cert["histograms"][f"cal_{lbl}"] = dict(cnt)
    good = cnt.get(lbl, 0) / sum(cnt.values())
    cert["readout_calibration"][lbl] = {"correct_fraction": good}
    print(f"cal {lbl}: correct {good:.3f}")

for label, corrupt in (("clean", False), ("corrupted", True)):
    Ws = {}
    for basis in ("Z", "X"):
        qc = transpile(ghz(corrupt, basis == "X"), be,
                       initial_layout=[a, b, c], optimization_level=3)
        cnt = counts_of(sampler.run([qc], shots=SHOTS).result()[0])
        cert["histograms"][f"{label}_{basis}"] = dict(cnt)
        Ws[basis] = forbidden(cnt, basis == "X")
    W = (Ws["Z"] + Ws["X"]) / 2
    se = max((P0 * (1 - P0) / (2 * SHOTS)) ** 0.5, 1e-9)
    z = (W - P0) / se
    v = "NOT_FALSIFIED" if z < 5 else "FALSIFIED_MEASUREMENT_CONTRACT"
    cert["rungs"][label] = {"W_Z": Ws["Z"], "W_X": Ws["X"], "W": W,
                            "p0": P0, "z": z, "verdict": v}
    print(f"{be.name} {label}: W={W:.4f} z={z:+.1f} -> {v}")

body = json.dumps({k: v for k, v in cert.items() if k != "_pin"},
                  sort_keys=True)
cert["_pin"] = hashlib.sha256(body.encode()).hexdigest()
open("AZQ1_I3_CERTIFICATE.json", "w").write(
    json.dumps(cert, indent=2, sort_keys=True))
print("written AZQ1_I3_CERTIFICATE.json — pin:", cert["_pin"])
