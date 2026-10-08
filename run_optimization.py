# -*- coding: utf-8 -*-
"""Unified optimization-round entry point.

Writes only to results_v41_optimization/.  Original results/, results_review/
and results_v41/ are never modified.

Usage:
    python run_optimization.py            # full run
    python run_optimization.py --quick    # smoke mode (not for reporting)

Exit code 0 iff every registered check passed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_v41_optimization")


def main():
    quick = "--quick" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    os.chdir(HERE)  # OUT_TAG inside the module is relative
    import optimization_v41 as opt
    ok = opt.run_all(quick=quick)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
