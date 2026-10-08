# -*- coding: utf-8 -*-
"""Small I/O helpers (standard library only)."""
import csv
import json
import os

import numpy as np

RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
FIGURES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")


def set_dirs(results=None, figures=None):
    """Redirect output locations (used by the V4.1 entry point so that the
    original results/ and figures/ are never overwritten)."""
    global RESULTS_DIR, FIGURES_DIR
    if results is not None:
        RESULTS_DIR = results
    if figures is not None:
        FIGURES_DIR = figures


def ensure_dirs():
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)


def write_csv(path, fieldnames, rows):
    """rows: iterable of dicts (subset of fieldnames is fine)."""
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (_fmt(r.get(k, ""))) for k in fieldnames})


def _fmt(v):
    if isinstance(v, (bool, np.bool_)):
        return int(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.12g}"
    return v


def write_json(path, obj):
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)

    class Enc(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, (np.floating, np.integer)):
                return o.item()
            if isinstance(o, np.ndarray):
                return o.tolist()
            return super().default(o)

    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, cls=Enc)


def save_fig(fig, name):
    """Save a matplotlib figure as PDF + PNG."""
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fig.savefig(os.path.join(FIGURES_DIR, name + ".pdf"))
    fig.savefig(os.path.join(FIGURES_DIR, name + ".png"), dpi=150)
