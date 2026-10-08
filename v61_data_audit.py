# -*- coding: utf-8 -*-
"""V6.1 data audit: re-aggregate every batch from raw CSVs (not from
hand-written reports), write the corrected dB map, and emit the V6.1
consistency checklist.  Output: results_v61/.  Historical results are
never modified."""
import csv
import json
import os
from datetime import datetime

import numpy as np

OUT = "results_v61"
FAIL = []


def rd(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def unit_key(r):
    return (r["seed"], r["cand_idx"], r["rho"], r.get("gamma", "0.5"))


def accept_stats(rows, methods):
    units = {}
    for r in rows:
        units.setdefault(unit_key(r), []).append(r)
    feas = 0
    acc = {m: 0 for m in methods}
    added = None
    n = len(units)
    for _, us in units.items():
        f = all(u["status"] == "clearly_feasible" for u in us)
        feas += f
        ok = {}
        for m in methods:
            ok[m] = all(float(u[m]) > float(u["tau"]) for u in us)
            acc[m] += ok[m]
        if added is None:
            added = 0
        if ok.get("Q_G") and ok.get("Q_old") is False:
            added += 1
    return dict(units=n, feasible=feas, **acc, qg_added=added,
                coverage_qg=(acc["Q_G"] / feas if feas else None),
                coverage_old=(acc["Q_old"] / feas if feas else None))


def main():
    os.makedirs(OUT, exist_ok=True)
    audit = []

    # 1) easy pool (results/: original run, gamma=0.5)
    rows = rd("results/pool_user_details.csv")
    st = accept_stats(rows, ["D", "Q_old", "Q_G", "L_quad"])
    audit.append(dict(batch="easy (seeds 11-55)", source="results/",
                      **st, note="saturated: zero additional accepts"))
    # 2) sensitivity (results_review/stress_details.csv)
    rows = rd("results_review/stress_details.csv")
    st = accept_stats(rows, ["D", "Q_old", "Q_G", "L_quad"])
    audit.append(dict(batch="sensitivity (seeds 11-110)",
                      source="results_review/", **st))
    # 3) validation (results_v41/)
    rows = rd("results_v41/validation_details.csv")
    st = accept_stats(rows, ["D", "Q_old", "Q_G", "L_quad"])
    audit.append(dict(batch="validation (seeds 211-215)",
                      source="results_v41/", **st))
    # 4) heterogeneous-link pool (results_v6/)
    rows = rd("results_v6/fa_realistic_details.csv")
    for r in rows:
        r["Q_G"] = r["Q_G"]
        r["L_quad"] = r["L_quad"]
    st = accept_stats(rows, ["D", "Q_old", "Q_G", "L_quad"])
    audit.append(dict(batch="heterogeneous links (seeds 301-305)",
                      source="results_v6/", **st,
                      note="illustrative pool; dB = 20log10(rho) is a "
                           "deterministic energy-ratio interpretation"))
    # joint accepts (max of Q_G, L_quad) where columns available
    joint_rows = []
    for name, path in [("sensitivity", "results_review/stress_details.csv"),
                       ("validation", "results_v41/validation_details.csv"),
                       ("het", "results_v6/fa_realistic_details.csv")]:
        for r in rd(path):
            r2 = dict(r)
            r2["joint"] = max(float(r["Q_G"]), float(r["L_quad"]))
            joint_rows.append((name, r2))
    for name in ("sensitivity", "validation", "het"):
        rows = [r for n, r in joint_rows if n == name]
        st = accept_stats(rows, ["Q_old", "Q_G", "L_quad", "joint"])
        audit.append(dict(batch=f"joint ({name})", source="(derived)",
                          **{k: v for k, v in st.items()}, joint_only=True))

    fieldnames = sorted({k for a in audit for k in a})
    with open(os.path.join(OUT, "data_audit.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for a in audit:
            w.writerow({k: (f"{a[k]:.4f}" if isinstance(a.get(k), float)
                            else a.get(k, "")) for k in fieldnames})

    # 5) corrected dB map (naming fix only; no recomputation)
    with open(os.path.join(OUT, "fa_error_ratio_db_map.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rho", "max_normalized_error_energy_ratio_db",
                    "interpretation"])
        for rho in (0.06, 0.15, 0.30):
            w.writerow([rho, round(20 * np.log10(rho), 2),
                        "deterministic ball r = rho*||h_hat||: sup ||e||^2/"
                        "||h_hat||^2 = rho^2; dB value 20*log10(rho). NOT a "
                        "statistical NMSE unless an error distribution, "
                        "normalization and coverage probability are "
                        "separately specified"])
        for db in (-20.0, -15.0, -10.0):
            w.writerow([round(10 ** (db / 20.0), 4), db,
                        "reverse mapping dB -> relative radius"])

    # 6) checklist
    checks = []
    easy = audit[0]
    checks.append(("easy pool shows zero additional accepts (no "
                   "'positive in every batch' claim)",
                   easy["qg_added"] == 0 and easy["Q_G"] == easy["Q_old"]))
    sens = [a for a in audit if a["batch"].startswith("sensitivity")][0]
    val = [a for a in audit if a["batch"].startswith("validation")][0]
    checks.append(("sensitivity Q_G added = 10", sens["qg_added"] == 10))
    checks.append(("validation Q_G added = 1", val["qg_added"] == 1))
    for name, cond in checks:
        print(f"[audit:{'PASS' if cond else 'FAIL'}] {name}")
        if not cond:
            FAIL.append(name)
    utils_json = dict(
        written_at=datetime.now().isoformat(timespec="seconds"),
        batches=[{k: v for k, v in a.items()} for a in audit],
        checks=[dict(name=n, ok=bool(c)) for n, c in checks],
        failures=FAIL)
    with open(os.path.join(OUT, "audit_summary.json"), "w",
              encoding="utf-8") as f:
        json.dump(utils_json, f, indent=1, default=str)
    print("\nbatch summary:")
    for a in audit:
        print(" ", {k: a[k] for k in ("batch", "units", "feasible",
                                       "Q_old", "Q_G", "qg_added")
                    if k in a})
    return not FAIL


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
