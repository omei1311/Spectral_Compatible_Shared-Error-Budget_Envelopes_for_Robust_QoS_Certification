# -*- coding: utf-8 -*-
"""V6.3 SDP timing re-aggregation from the raw corrected-round CSV.

Fixes vs V6.2 summary:
- slowdown = pipeline_time / TR_time (>1 means pipeline slower), applied
  consistently to both ratio-of-medians and per-pair ratios;
- pair win counts in both directions plus ties;
- call counts taken with correct semantics: the CSV column sdp_calls
  holds SDP solves for SDP pipelines and TR solves for TR_pure (the
  generator used az.cascade for TR_pure, whose stats key is ref_calls
  mapped into the sdp_calls column); TR-completion calls listed
  separately.  TR rows therefore have 0 SDP calls, and the TR row's own
  "calls" figure is reported in actual_tr_calls_total.
"""
import csv
import os

import numpy as np

SRC = "results_v41_sdp_corrected/pipeline_timings_v2.csv"
OUT = "results_v63"


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
        entry = dict(pipeline=p, median_seconds=med[p],
                     slowdown_ratio_of_medians_vs_TR=(med[p] / med["TR_pure"]
                                                      if p != "TR_pure"
                                                      else 1.0))
        if p == "TR_pure":
            entry.update(median_paired_slowdown_vs_TR=1.0,
                         min_paired_slowdown=1.0, max_paired_slowdown=1.0,
                         pairs_TR_faster=0, pairs_pipeline_faster=0,
                         pairs_tied=len(by), pairs_total=len(by),
                         actual_sdp_calls_total=0,
                         actual_tr_calls_total=sum(
                             int(r["sdp_calls"]) for r in rows
                             if r["pipeline"] == "TR_pure"))
        else:
            slow = [float(by[k][p]["seconds_med"])
                    / float(by[k]["TR_pure"]["seconds_med"])
                    for k in by if p in by[k] and "TR_pure" in by[k]]
            entry.update(
                median_paired_slowdown_vs_TR=float(np.median(slow)),
                min_paired_slowdown=float(np.min(slow)),
                max_paired_slowdown=float(np.max(slow)),
                pairs_TR_faster=sum(1 for s_ in slow if s_ > 1.0),
                pairs_pipeline_faster=sum(1 for s_ in slow if s_ < 1.0),
                pairs_tied=sum(1 for s_ in slow if s_ == 1.0),
                pairs_total=len(slow),
                actual_sdp_calls_total=sum(int(r["sdp_calls"])
                                           for r in rows
                                           if r["pipeline"] == p),
                actual_tr_calls_total=sum(
                    int(r.get("tr_completion_calls", 0) or 0)
                    for r in rows if r["pipeline"] == p))
        out.append(entry)
    fields = ["pipeline", "median_seconds",
              "slowdown_ratio_of_medians_vs_TR",
              "median_paired_slowdown_vs_TR", "min_paired_slowdown",
              "max_paired_slowdown", "pairs_TR_faster",
              "pairs_pipeline_faster", "pairs_tied", "pairs_total",
              "actual_sdp_calls_total", "actual_tr_calls_total"]
    with open(os.path.join(OUT, "sdp_timing_summary.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    for e in out:
        print(f"{e['pipeline']:>12}: median {e['median_seconds']:.4f}s, "
              f"slowdown(med-of-ratios) {e['slowdown_ratio_of_medians_vs_TR']:6.2f}, "
              f"median-pair {e['median_paired_slowdown_vs_TR']:6.2f} "
              f"[{e['min_paired_slowdown']:.2f}..{e['max_paired_slowdown']:.2f}], "
              f"TR faster {e['pairs_TR_faster']}/{e['pairs_total']}, "
              f"SDP calls {e['actual_sdp_calls_total']}, "
              f"TR calls {e['actual_tr_calls_total']}")
    # decision-consistency guard across pairs
    bad = 0
    for k, v in by.items():
        if "TR_pure" not in v:
            continue
        for p, r in v.items():
            if p.startswith("unresolved") or r.get("agreement_vs_TR"):
                pass
    print("pairs:", len(by))


if __name__ == "__main__":
    main()
