# -*- coding: utf-8 -*-
"""Unified V4.1 entry point.

Generates theory checks, spectral scan, candidate pools (baseline +
exploration + pre-registered validation), resolution ablation, overlap
accounting and timing into results_v41/.  Original results/ and
results_review/ are never touched.

Usage:
    python run_v41.py            # full run
    python run_v41.py --quick    # small-scale smoke mode (not for reporting)

Exit code 0 iff every registered check passed.
"""
import os
import sys

import utils

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_V41 = os.path.join(HERE, "results_v41")


def main():
    quick = "--quick" in sys.argv
    utils.set_dirs(results=RESULTS_V41,
                   figures=os.path.join(RESULTS_V41, "figures"))
    utils.ensure_dirs()
    import experiments_v41 as ev
    ok = ev.run_all(quick=quick)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
