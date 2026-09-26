#!/usr/bin/env python
"""Summarise the host CIBERSORT v1.04 runs: per-sample fractions and fit, and primaries v recurrences per cell type.

Groups: primaries IL67B IL68B IL69B; recurrences IL66B NL70B NL71B; controls IL64B N168B N269B.
3 v 3 per cell type: Welch t-test (two-sided) and the exact Mann-Whitney U (two-sided; with 3 v 3 its smallest
possible p is 0.10, so it can never reach 0.05). BH over the cell types of one signature, on the Welch p.

    python summarise_host_deconv.py fractions_v104_host_lm22.tsv fractions_v104_host_zhang2014_brain_mean7_tpm.tsv ... --out host_deconv_summary.tsv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PRIM = ["IL67B", "IL68B", "IL69B"]
REC = ["IL66B", "NL70B", "NL71B"]
CTRL = ["IL64B", "N168B", "N269B"]
FIT = ["P-value", "Correlation", "RMSE"]


def bh(p):
    p = np.asarray(p, float); q = np.full_like(p, np.nan); ok = ~np.isnan(p)
    if ok.sum():
        r = p[ok]; o = np.argsort(r); n = len(r)
        adj = np.minimum.accumulate((r[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        tmp = np.empty(n); tmp[o] = np.minimum(adj, 1); q[ok] = tmp
    return q


def one(path):
    f = pd.read_csv(path, sep="\t", index_col=0)
    name = Path(path).stem.replace("fractions_v104_", "")
    cells = [c for c in f.columns if c not in FIT + ["null_fitted", "null_degenerate"]]
    rows = []
    for c in cells:
        a, b = f.loc[PRIM, c].astype(float), f.loc[REC, c].astype(float)
        if a.isna().any() or b.isna().any():
            continue
        if a.std() == 0 and b.std() == 0:
            tp = np.nan if a.mean() == b.mean() else 0.0
        else:
            tp = stats.ttest_ind(b, a, equal_var=False).pvalue
        up = stats.mannwhitneyu(b, a, alternative="two-sided", method="exact").pvalue
        rows.append({"signature": name, "cell": c, "primary_mean": a.mean(), "recurrence_mean": b.mean(),
                     "difference": b.mean() - a.mean(), "welch_p": tp, "mann_whitney_p": up,
                     "control_mean": f.loc[CTRL, c].astype(float).mean(),
                     "nonzero_tumours": int((f.loc[PRIM + REC, c] > 0).sum())})
    t = pd.DataFrame(rows)
    if len(t):
        t["welch_q_bh"] = bh(t["welch_p"].values)
    return name, f, cells, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fractions", nargs="+")
    ap.add_argument("--out")
    a = ap.parse_args()
    allt = []
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    for p in a.fractions:
        name, f, cells, t = one(p)
        order = PRIM + REC + CTRL
        print(f"\n## {name}  ({len(cells)} cell types)")
        extra = [c for c in ["null_fitted", "null_degenerate"] if c in f.columns]
        print(f.loc[order, FIT + extra].round(3).to_string())
        print("fit, tumours: r %.2f-%.2f, RMSE %.2f-%.2f, P %.3f-%.3f; controls: r %.2f-%.2f, P %.3f-%.3f" % (
            f.loc[PRIM + REC, "Correlation"].min(), f.loc[PRIM + REC, "Correlation"].max(),
            f.loc[PRIM + REC, "RMSE"].min(), f.loc[PRIM + REC, "RMSE"].max(),
            f.loc[PRIM + REC, "P-value"].min(), f.loc[PRIM + REC, "P-value"].max(),
            f.loc[CTRL, "Correlation"].min(), f.loc[CTRL, "Correlation"].max(),
            f.loc[CTRL, "P-value"].min(), f.loc[CTRL, "P-value"].max()))
        keep = [c for c in cells if f.loc[order, c].max() >= 0.02]
        print(f.loc[order, keep].round(3).T.to_string())
        if len(t):
            print(t.drop(columns="signature").round(3).to_string(index=False))
            allt.append(t)
    if a.out and allt:
        pd.concat(allt).to_csv(a.out, sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    sys.exit(main())
