# -*- coding: utf-8 -*-
"""V6.2 SDP timing re-aggregation from the corrected round only
(results_v41_sdp_corrected/pipeline_timings_v2.csv).

Reports per pipeline: median seconds, ratio of medians vs TR_pure,
median of per-(seed,gamma) pair ratios, pair win counts, spread and calls.
Both ratio statistics are reported separately and never mixed.
Writes results_v62/sdp_timing_summary.csv."""
import csv
import os

import numpy as np

SRC = "results_v41_sdp_corrected/pipeline_timings_v2.csv"
OUT = "results_v62"


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    by = {}
    for r in rows:
        by.setdefault((r["seed"], r["gamma"]), {})[r["pipeline"]] = r
    pipes = sorted({r["pipeline"] for r in rows})
    med = {p: float(np.median([float(r["seconds_med"]) for r in rows
                               if r["pipeline"] == p])) for p in pipes}
    out = []
    for p in pipes:
        base = "TR_pure"
        if p == base:
            continue
        ratios = [float(by[k][base]["seconds_med"])
                  / float(by[k][p]["seconds_med"])
                  for k in by if p in by[k]]
        wins = sum(1 for rr in ratios if rr > 1.0)
        out.append(dict(
            pipeline=p, median_s=med[p],
            ratio_of_medians_vs_TR=med[base] / med[p],
            median_pair_ratio_vs_TR=float(np.median(ratios)),
            min_pair_ratio=float(np.min(ratios)),
            max_pair_ratio=float(np.max(ratios)),
            pair_wins_TR_faster=wins, pairs=len(ratios),
            sdp_calls_total=sum(int(r["sdp_calls"]) for r in rows
                                if r["pipeline"] == p),
            tr_calls_total=sum(int(r["sdp_calls"]) for r in rows
                               if r["pipeline"] == base)))
    out.insert(0, dict(pipeline="TR_pure", median_s=med["TR_pure"],
                       ratio_of_medians_vs_TR=1.0,
                       median_pair_ratio_vs_TR=1.0, min_pair_ratio=1.0,
                       max_pair_ratio=1.0, pair_wins_TR_faster=0,
                       pairs=len(by),
                       sdp_calls_total=sum(int(r["sdp_calls"]) for r in rows
                                           if r["pipeline"] == "TR_pure"),
                       tr_calls_total=0))
    with open(os.path.join(OUT, "sdp_timing_summary.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    for o in out:
        print(f"{o['pipeline']:>12}: median {o['median_s']:.4f} s | "
              f"ratio-of-medians {o['ratio_of_medians_vs_TR']:.3f} | "
              f"median pair ratio {o['median_pair_ratio_vs_TR']:.3f} | "
              f"TR faster in {o['pair_wins_TR_faster']}/{o['pairs']}")
    return True


if __name__ == "__main__":
    main()
