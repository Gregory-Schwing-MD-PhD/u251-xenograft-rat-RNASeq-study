#!/usr/bin/env python
"""ANALYSIS/cibersort/cibersort_neftel.py - glioblastoma cell-state deconvolution of the U251 bulk profiles against the
Neftel 2019 single-cell reference (Smart-seq2, Broad Single Cell Portal SCP393), ahead of an official CIBERSORT run.

Reference cells: the adult malignant cells of the IDH-wt SS2 data (4,916 cells, 20 tumours). Each cell's state is the
largest of the four collapsed meta-module scores in the metadata (MES = max(MESlike1, MESlike2), AC, OPC,
NPC = max(NPClike1, NPClike2)). A cell is "confident" when that score is positive and beats the runner-up by at least
--margin; this cut is ours, not Neftel's hybrid rule, and a variant using every cell is written beside it.
Expression: the SCP393 matrix is log2(TPM/10 + 1); linear TPM = 10 * (2^x - 1). A different base scale would only
rescale the signature, which CIBERSORT standardises away; the per-cell sums are logged as a check.

Signature matrices follow CIBERSORT's published procedure (Newman 2015): per cell type, a Welch t-test on log
expression against all other reference cells, genes with q < 0.3 and higher in the type, ranked by fold change;
the top G per type for G = 50..200, keeping the G that minimises the condition number of the linear-mean matrix.
Only genes present in the bulk mixture are candidates.

Deconvolution here is a re-implementation of CIBERSORT's core (nu-SVR, linear kernel, nu 0.25/0.5/0.75, the fit with
the lowest RMSE; negative weights set to 0 and the rest normalised; signature standardised as one block, each mixture
per sample; no quantile normalisation, as advised for RNA-seq; p from 1000 permuted mixtures drawn from the pooled
mixture values). It exists to check the reference before, and the official CIBERSORT after, the reinstall.

Leave-one-tumour-out check: for each reference tumour, the signature is rebuilt without it, its malignant cells are
averaged into a pseudo-bulk, and the estimated state fractions are compared with the true ones (share of its cells
in each state). Same platform on both sides, so it tests the reference, not the Smart-seq2 to bulk gap.

    python cibersort_neftel.py --neftel <refs/neftel_2019> --tpm <salmon.merged.gene_tpm.tsv> --cohorts <sample_cohorts.csv> --out results
"""
import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.svm import NuSVR

STATES = ["MES", "AC", "OPC", "NPC"]
NONMAL = ["Macrophage", "Oligodendrocyte", "T-cell"]
T0 = time.time()


def say(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def welch(a, b):
    """Welch t-test per row; a, b are genes x cells (log scale). Returns t and two-sided p."""
    n1, n2 = a.shape[1], b.shape[1]
    m1, m2 = a.mean(1, dtype=np.float64), b.mean(1, dtype=np.float64)
    v1, v2 = a.var(1, ddof=1, dtype=np.float64), b.var(1, ddof=1, dtype=np.float64)
    s1, s2 = v1 / n1, v2 / n2
    se = np.sqrt(s1 + s2)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, (m1 - m2) / se, 0.0)
        df = np.where(se > 0, (s1 + s2) ** 2 / (s1 ** 2 / (n1 - 1) + s2 ** 2 / (n2 - 1)), 1.0)
    p = 2 * stats.t.sf(np.abs(t), df)
    return t, np.where(se > 0, p, 1.0)


def build_signature(logm, lin, labels, types, genes, gmin=50, gmax=200, step=10, qcut=0.3):
    """logm, lin: genes x cells (log / linear) restricted to the candidate genes; labels per cell."""
    labels = np.asarray(labels)
    mu = pd.DataFrame({t: lin[:, labels == t].mean(1, dtype=np.float64) for t in types}, index=genes)
    ranked = {}
    for t in types:
        a, b = labels == t, labels != t
        tt, p = welch(logm[:, a], logm[:, b])
        q = bh(p)
        rest = lin[:, b].mean(1, dtype=np.float64)
        fc = (mu[t].values + 1.0) / (rest + 1.0)
        ok = (q < qcut) & (tt > 0) & (fc > 1)
        order = np.argsort(-fc[ok])
        ranked[t] = np.asarray(genes)[ok][order]
    best = None
    for g in range(gmin, gmax + 1, step):
        sel = sorted(set().union(*[set(ranked[t][:g]) for t in types]))
        k = np.linalg.cond(mu.loc[sel].values)
        if best is None or k < best[1]:
            best = (g, k, sel)
    g, k, sel = best
    return mu.loc[sel], {"G": g, "kappa": float(k), "n_genes": len(sel),
                         "candidates_per_type": {t: int(len(ranked[t])) for t in types}}


def core(X, y):
    best = None
    for nu in (0.25, 0.5, 0.75):
        m = NuSVR(kernel="linear", nu=nu, C=1.0, tol=1e-3).fit(X, y)
        w = m.coef_.ravel().copy()
        w[w < 0] = 0
        if w.sum() <= 0:
            continue
        w = w / w.sum()
        k = X @ w
        rmse = float(np.sqrt(np.mean((k - y) ** 2)))
        r = float(np.corrcoef(k, y)[0, 1])
        if best is None or rmse < best[1]:
            best = (w, rmse, r, nu)
    if best is None:
        return np.full(X.shape[1], np.nan), np.nan, np.nan, np.nan
    return best


def deconvolve(sig, mix, perm=1000, seed=1):
    genes = sig.index.intersection(mix.index)
    X = sig.loc[genes].values.astype(float)
    X = (X - X.mean()) / X.std(ddof=1)
    Y = mix.loc[genes].values.astype(float)
    null = np.array([])
    if perm:
        rng = np.random.default_rng(seed)
        pool = Y.ravel()
        rs = []
        for _ in range(perm):
            yr = rng.choice(pool, size=len(genes), replace=False)
            yr = (yr - yr.mean()) / yr.std(ddof=1)
            rs.append(core(X, yr)[2])
        null = np.sort(np.asarray(rs))
    rows = []
    for j, s in enumerate(mix.columns):
        y = Y[:, j]
        y = (y - y.mean()) / y.std(ddof=1)
        w, rmse, r, nu = core(X, y)
        p = float(np.mean(null >= r)) if perm else np.nan
        rows.append({"sample": s, **dict(zip(sig.columns, w)), "P": p, "R": r, "RMSE": rmse, "nu": nu, "n_genes": len(genes)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--neftel", required=True)
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--cohorts", required=True)
    ap.add_argument("--out", default="results")
    ap.add_argument("--margin", type=float, default=0.5)
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--scref-cells", type=int, default=300, help="cells per state in the CIBERSORTx reference")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    log = {}

    # ---- reference cells and states
    meta = pd.read_csv(Path(a.neftel) / "IDHwt.GBM.Metadata.SS2.txt", sep="\t", skiprows=[1], low_memory=False)
    for c in ["MESlike2", "MESlike1", "AClike", "OPClike", "NPClike1", "NPClike2"]:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")
    mal = meta[(meta.CellAssignment == "Malignant") & (meta.GBMType == "Adult")].copy()
    sc = pd.DataFrame({"MES": mal[["MESlike1", "MESlike2"]].max(1), "AC": mal.AClike, "OPC": mal.OPClike,
                       "NPC": mal[["NPClike1", "NPClike2"]].max(1)})
    scored = sc.notna().all(1).values
    n_unscored = int((~scored).sum())
    mal, sc = mal[scored].copy(), sc[scored]
    srt = np.sort(sc.values, 1)
    mal["state"] = sc.idxmax(1).values
    mal["top"] = srt[:, -1]
    mal["margin"] = srt[:, -1] - srt[:, -2]
    mal["confident"] = (mal.top > 0) & (mal.margin >= a.margin)
    for s in STATES:
        mal["score_" + s] = sc[s].values
    mal[["NAME", "Sample", "state", "top", "margin", "confident"] + ["score_" + s for s in STATES]].to_csv(out / "cells_states.tsv", sep="\t", index=False)
    log["cells"] = {"adult_malignant": int(len(mal)), "unscored_dropped": n_unscored, "tumours": int(mal.Sample.nunique()),
                    "state_all": mal.state.value_counts().to_dict(),
                    "state_confident": mal[mal.confident].state.value_counts().to_dict(),
                    "margin": a.margin}
    nm = meta[meta.CellAssignment.isin(NONMAL)]
    nm_ad = nm[nm.GBMType == "Adult"]
    log["cells"]["nonmalignant_adult"] = nm_ad.CellAssignment.value_counts().to_dict()
    log["cells"]["nonmalignant_all"] = nm.CellAssignment.value_counts().to_dict()
    vc = nm_ad.CellAssignment.value_counts()
    nm_use = nm_ad if all(vc.get(t, 0) >= 30 for t in NONMAL) else nm
    log["cells"]["nonmalignant_used"] = "adult" if nm_use is nm_ad else "adult+pediatric"
    say("cells", json.dumps(log["cells"]))

    # ---- expression
    say("reading the SS2 matrix")
    expr = pd.read_csv(Path(a.neftel) / "IDHwtGBM.processed.SS2.logTPM.txt.gz", sep="\t", index_col=0)
    expr = expr[~expr.index.duplicated()]
    keep = list(mal.NAME) + list(nm_use.NAME)
    missing = [c for c in keep if c not in expr.columns]
    assert not missing, f"{len(missing)} metadata cells are not in the matrix"
    expr = expr[keep].astype(np.float32)
    lin_all = (10.0 * (np.exp2(expr.values) - 1.0)).astype(np.float32)
    sums = lin_all[:, : len(mal)].sum(0)
    log["scale_check"] = {"per_cell_sum_linear_tpm_median": float(np.median(sums)), "p5": float(np.percentile(sums, 5)),
                          "p95": float(np.percentile(sums, 95)), "genes": int(expr.shape[0]),
                          "note": "10*(2^x-1); near 1e6 if the matrix is log2(TPM/10+1) over all genes"}
    say("scale", json.dumps(log["scale_check"]))

    # ---- bulk mixture
    tpm = pd.read_csv(a.tpm, sep="\t")
    samples = [c for c in tpm.columns if c not in ("gene_id", "gene_name")]
    mix = tpm.groupby("gene_name")[samples].sum()
    mix.index.name = "GeneSymbol"
    mix.to_csv(out / "mixture_u251_tpm.txt", sep="\t")
    coh = pd.read_csv(a.cohorts).set_index("sample")["cohort"].to_dict() if Path(a.cohorts).exists() else {}
    universe = expr.index.intersection(mix.index)
    gi = expr.index.get_indexer(universe)
    log["genes"] = {"neftel": int(expr.shape[0]), "mixture_symbols": int(mix.shape[0]), "shared": int(len(universe))}
    say("genes", json.dumps(log["genes"]))
    logm = expr.values[gi]
    lin = lin_all[gi]
    ncol = len(mal)
    lab_state = mal.state.values
    conf = mal.confident.values
    tum = mal.Sample.values

    # ---- signatures
    sigs, siglog = {}, {}
    variants = {
        "neftel4_confident": (np.where(conf)[0], lab_state[conf], STATES),
        "neftel4_allcells": (np.arange(ncol), lab_state, STATES),
    }
    idx_nm = np.arange(ncol, ncol + len(nm_use))
    variants["neftel4_confident_plus_nonmalignant"] = (
        np.concatenate([np.where(conf)[0], idx_nm]),
        np.concatenate([lab_state[conf], nm_use.CellAssignment.values]), STATES + NONMAL)
    for name, (cols, labs, types) in variants.items():
        s, info = build_signature(logm[:, cols], lin[:, cols], labs, types, universe)
        sigs[name], siglog[name] = s, info
        s.index.name = "GeneSymbol"
        s.to_csv(out / f"signature_{name}.txt", sep="\t", float_format="%.4f")
        say("signature", name, json.dumps(info))
    log["signatures"] = siglog

    # ---- CIBERSORTx single-cell reference (header = state), sampled across tumours
    rng = np.random.default_rng(7)
    pick = []
    for s in STATES:
        ids = np.where(conf & (lab_state == s))[0]
        pick += list(rng.choice(ids, size=min(a.scref_cells, len(ids)), replace=False))
    ref = pd.DataFrame(lin_all[:, pick], index=expr.index, columns=lab_state[pick])
    ref.index.name = "GeneSymbol"
    with gzip.open(out / "scref_neftel4_cibersortx.txt.gz", "wt") as fh:
        ref.to_csv(fh, sep="\t", float_format="%.3f")
    log["scref"] = {"cells": pd.Series(lab_state[pick]).value_counts().to_dict(), "genes": int(ref.shape[0])}

    # ---- leave-one-tumour-out check (malignant columns only)
    logm_m, lin_m = logm[:, :ncol], lin[:, :ncol]
    rows = []
    for name in ("neftel4_confident", "neftel4_allcells"):
        for t in sorted(set(tum)):
            held = tum == t
            if held.sum() < 30:
                continue
            train = (~held) & (conf if name == "neftel4_confident" else True)
            s, _ = build_signature(logm_m[:, train], lin_m[:, train], lab_state[train], STATES, universe)
            pb = pd.DataFrame({t: lin_m[:, held].mean(1, dtype=np.float64)}, index=universe)
            est = deconvolve(s, pb, perm=0).iloc[0]
            truth = pd.Series(lab_state[held]).value_counts(normalize=True).reindex(STATES).fillna(0)
            for st in STATES:
                rows.append({"variant": name, "tumour": t, "cells": int(held.sum()), "state": st,
                             "true": float(truth[st]), "est": float(est[st]), "R_fit": float(est["R"])})
        say("loto done", name)
    lo = pd.DataFrame(rows)
    lo.to_csv(out / "loto_validation.tsv", sep="\t", index=False)
    summ = {}
    for name, g in lo.groupby("variant"):
        per_t = g.groupby("tumour").apply(lambda d: np.corrcoef(d.true, d.est)[0, 1] if d.true.std() > 0 else np.nan)
        summ[name] = {"tumours": int(g.tumour.nunique()),
                      "pooled_r": float(np.corrcoef(g.true, g.est)[0, 1]),
                      "rmse": float(np.sqrt(np.mean((g.true - g.est) ** 2))),
                      "median_per_tumour_r": float(np.nanmedian(per_t)),
                      "mae_by_state": g.assign(e=(g.est - g.true).abs()).groupby("state").e.mean().round(4).to_dict(),
                      "bias_by_state": g.assign(e=(g.est - g.true)).groupby("state").e.mean().round(4).to_dict()}
    log["loto"] = summ
    say("loto", json.dumps(summ))

    # ---- the bulk samples
    res = []
    for name, s in sigs.items():
        d = deconvolve(s, mix, perm=a.perm)
        d.insert(1, "cohort", d["sample"].map(coh).fillna(""))
        d.insert(0, "signature", name)
        res.append(d)
        say("bulk", name)
    fr = pd.concat(res, ignore_index=True)
    fr.to_csv(out / "fractions_nusvr.tsv", sep="\t", index=False, float_format="%.4f")
    print(fr.to_string(float_format=lambda v: f"{v:.3f}"))
    import sklearn, scipy
    log["versions"] = {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "sklearn": sklearn.__version__}
    json.dump(log, open(out / "build_log.json", "w"), indent=1, default=str)
    say("done")


if __name__ == "__main__":
    main()
