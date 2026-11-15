#!/usr/bin/env python
"""Host (rat) cell-state deconvolution inputs: a mixture and signature matrices in one gene space, human symbols.

Mixture: the rat host TPM (xengsort `host` bin, nf-core/rnaseq salmon), rat genes mapped to human through 1:1
orthologs only (anything 1:many is dropped rather than guessed), summed per human symbol.
Signatures (CIBERSORT format, linear, human symbols):
  LM22 (Newman 2015; human leukocytes) as is. Nude rats (Foxn1-null) have essentially no T cells, so T-cell weight is
    a built-in check on the method.
  Every per-population mean matrix from a mouse reference (Zhang 2014 brain cell types; Bowman 2016 glioma
    microglia vs bone-marrow macrophages; any other --mouse-mean file) is mapped mouse->human by 1:1 orthologs, then
    reduced to a signature the LM22 way: per population, genes >= 5 units and >= 2x every other population, ranked by
    that fold; the top G per population for G = 25..200, keeping the G with the lowest condition number.
  A ready-made mouse signature matrix (--mouse-sig, e.g. seq-ImmuCC) is only mapped to human symbols.

    python host_deconv_prep.py --tpm salmon.merged.gene_tpm.tsv --rat2human rat_human_1to1.tsv --mouse2human mouse_human_1to1.tsv \
        --lm22 lm22_newman2015.txt --mouse-mean zhang2014_mean.txt bowman2016_mean.txt --mouse-sig seqimmucc.txt --out results
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = ["IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B", "IL64B", "N168B", "N269B"]


def one_to_one(path, src_col_hint):
    t = pd.read_csv(path, sep="\t", dtype=str)
    src = next(c for c in t.columns if src_col_hint in c.lower() and "symbol" in c.lower())
    hum = next(c for c in t.columns if "human" in c.lower() and "symbol" in c.lower())
    t = t[[src, hum]].dropna().drop_duplicates()
    t = t[~t[src].duplicated(keep=False) & ~t[hum].duplicated(keep=False)]
    return dict(zip(t[src], t[hum]))


def to_human(df, mapping):
    df = df[df.index.isin(mapping)]
    df.index = df.index.map(mapping)
    return df.groupby(level=0).sum()


def signature(mean, gmin=25, gmax=200, step=25, floor=5.0, fold=2.0):
    ranked = {}
    for c in mean.columns:
        other = mean.drop(columns=c).max(axis=1)
        fc = mean[c] / other.clip(lower=1e-9)
        ok = (mean[c] >= floor) & (fc >= fold)
        ranked[c] = fc[ok].sort_values(ascending=False).index.tolist()
    best = None
    for g in range(gmin, gmax + 1, step):
        genes = sorted(set().union(*[set(v[:g]) for v in ranked.values()]))
        if len(genes) < len(mean.columns) * 5:
            continue
        k = float(np.linalg.cond(mean.loc[genes].values))
        if best is None or k < best[1]:
            best = (g, k, genes)
    g, k, genes = best
    return mean.loc[genes], {"G": g, "kappa": k, "genes": len(genes), "candidates": {c: len(v) for c, v in ranked.items()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--rat2human", required=True)
    ap.add_argument("--mouse2human")
    ap.add_argument("--lm22")
    ap.add_argument("--mouse-mean", nargs="*", default=[])
    ap.add_argument("--mouse-sig", nargs="*", default=[])
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    log = {}
    tpm = pd.read_csv(a.tpm, sep="\t")
    rat = tpm.groupby("gene_name")[SAMPLES].sum()
    r2h = one_to_one(a.rat2human, "rat")
    mix = to_human(rat, r2h)
    mix.index.name = "GeneSymbol"
    mix.to_csv(out / "mixture_host_hs.txt", sep="\t", float_format="%.4f")
    log["mixture"] = {"rat_genes": int(len(rat)), "one_to_one": len(r2h), "human_symbols": int(len(mix)),
                      "share_of_rat_tpm_kept": {s: round(float(mix[s].sum() / rat[s].sum()), 3) for s in SAMPLES}}
    sigs = {}
    if a.lm22:
        sigs["lm22"] = pd.read_csv(a.lm22, sep="\t", index_col=0)
    m2h = one_to_one(a.mouse2human, "mouse") if a.mouse2human else {}
    for f in a.mouse_mean:
        mean = to_human(pd.read_csv(f, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0), m2h)
        s, info = signature(mean)
        sigs[Path(f).stem] = s
        log.setdefault("built", {})[Path(f).stem] = info
    for f in a.mouse_sig:
        sigs[Path(f).stem + "_hs"] = to_human(pd.read_csv(f, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0), m2h)
    jobs = []
    for name, s in sigs.items():
        s.index.name = "GeneSymbol"
        p = out / f"signature_host_{name}.txt"
        s.to_csv(p, sep="\t", float_format="%.4f")
        n_match = int(s.index.isin(mix.index).sum())
        log.setdefault("signatures", {})[name] = {"genes": int(len(s)), "columns": list(s.columns), "genes_in_mixture": n_match}
        jobs.append(f"host_{name}\t{p.resolve()}\t{1000 if len(s) < 1500 else 100}\t{(out / 'mixture_host_hs.txt').resolve()}")
    (out / "v104_host_jobs.txt").write_text("\n".join(jobs) + "\n")
    json.dump(log, open(out / "host_deconv_prep_log.json", "w"), indent=1)
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
