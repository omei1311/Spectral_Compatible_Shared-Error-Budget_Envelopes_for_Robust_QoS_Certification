# -*- coding: utf-8 -*-
"""Corrected S-procedure-SDP baseline entry point (v2).

Writes only to results_v41_sdp_corrected/.  v1 results in results_v41_sdp/
are retained unchanged as history.
Usage:
    python run_sdp_v2.py            # tests + smallcheck + batch + pipelines
    python run_sdp_v2.py --quick    # tests + smallcheck only
Exit code 0 iff every registered check passed.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results_v41_sdp_corrected")


def main():
    quick = "--quick" in sys.argv
    os.makedirs(OUT, exist_ok=True)
    os.chdir(HERE)
    import utils
    utils.set_dirs(results=OUT)
    import sdp_baseline_v2 as sb
    ok = sb.run_all(quick=quick)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
