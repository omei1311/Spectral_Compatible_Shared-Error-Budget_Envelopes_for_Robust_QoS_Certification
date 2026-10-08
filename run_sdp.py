# -*- coding: utf-8 -*-
"""S-procedure-SDP baseline entry point.

Writes only to results_v41_sdp/.  Earlier result directories untouched.
Usage:
    python run_sdp.py            # smallcheck + full batch + pipelines
    python run_sdp.py --quick    # pre-fixed smallcheck only
Exit code 0 iff every registered check passed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_v41_sdp")


def main():
    quick = "--quick" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    os.chdir(HERE)
    import utils
    utils.set_dirs(results=OUT)
    import sdp_baseline as sb
    ok = sb.run_all(quick=quick)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
