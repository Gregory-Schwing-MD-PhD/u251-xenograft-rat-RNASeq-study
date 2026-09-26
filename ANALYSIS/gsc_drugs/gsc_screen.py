#!/usr/bin/env python
"""Drug discovery aimed at the tumour stem-cell (GSC) program of the U251 xenografts.

GSC gene sets
  varn_stem / varn_prolif_stem: Varn et al. 2022 Cell (GLASS) single-cell signature matrix, Table S3 - genes whose
    highest column is stemcell_tumor (or prolif_stemcell_tumor) and which are at least twice every other column
    (the differentiated tumour state and the nine host cell types). This is the reference that called U251 ~100 %
    proliferating stem-like under CIBERSORT v1.04.
  any extra list in --extra (one gene per line, header 'gene'), e.g. the published GSC programs.
  Every set is restricted to genes expressed in the xenografts (median TPM >= 1 over the six tumours, human reads).

Screens (DSigDB v1.0, the gene sets the published R stage used)
  S1 over-representation: which drug gene sets are enriched in the GSC genes (hypergeometric; universe = expressed
     genes in DSigDB; BH q). Effect = odds ratio.
  S2 recurrence within the stem program: preranked GSEA (gseapy) of DSigDB on the published DESeq2 result, ranked by
     sign(log2FC) x -log10(p) (the table has no Wald statistic) for Recurrent vs Primary, restricted to the GSC genes; sets of 5-500 genes. NES < 0 = the drug's genes fall in
     recurrence within the stem program (the published screen's direction); NES > 0 = they rise.
Then each drug with q < 0.25 in either screen: ChEMBL (max clinical phase, SMILES), ADMET-AI BBB_Martins; ranked
  among clinical-phase compounds (phase 1-4) by effect x BBB probability.

    python gsc_screen.py --dsigdb hs.dsigdb.v1.0.gmt --de therapy_impact.deseq2.results.tsv --tpm salmon.merged.gene_tpm.tsv \
        --varn varn2022_glass_sc12.txt [--extra list1.txt list2.txt] --out results
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

TUM = ["IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B"]


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def read_gmt(path):
    sets = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            sets[f[0]] = {g for g in f[2:] if g}
    return sets


def drug_name(term):
    return term.split("_")[0].strip().lower()


def chembl(name, cache):
    if name in cache:
        return cache[name]
    out = {"chembl_id": "", "max_phase": np.nan, "smiles": ""}
    try:
        r = requests.get("https://www.ebi.ac.uk/chembl/api/data/molecule/search.json",
                         params={"q": name, "limit": 5}, timeout=30).json()
        mols = r.get("molecules") or []
        exact = [m for m in mols if (m.get("pref_name") or "").lower() == name] or mols
        if exact:
            m = exact[0]
            out = {"chembl_id": m.get("molecule_chembl_id", ""),
                   "max_phase": pd.to_numeric(m.get("max_phase"), errors="coerce"),
                   "smiles": ((m.get("molecule_structures") or {}).get("canonical_smiles") or "")}
    except Exception as e:  # network hiccup: leave unresolved
        out["error"] = str(e)
    cache[name] = out
    time.sleep(0.15)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dsigdb", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--varn", required=True)
    ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--out", default="results")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    log = {}

    tpm = pd.read_csv(a.tpm, sep="\t")
    sym = tpm.set_index("gene_id")["gene_name"]
    expr = tpm.groupby("gene_name")[TUM].sum()
    expressed = set(expr.index[expr.median(axis=1) >= 1.0])
    log["expressed_genes"] = len(expressed)

    v = pd.read_csv(a.varn, sep="\t", index_col=0)
    sets = {}
    for col, name in (("stemcell_tumor", "varn_stem"), ("prolif_stemcell_tumor", "varn_prolif_stem")):
        others = v.drop(columns=[col])
        top = (v[col] >= 2 * others.max(axis=1)) & (v[col] > 0)
        sets[name] = set(v.index[top])
    sets["varn_stem_union"] = sets["varn_stem"] | sets["varn_prolif_stem"]
    for f in a.extra:
        g = pd.read_csv(f, sep="\t")
        sets[Path(f).stem] = set(g.iloc[:, 0].astype(str))
    sets = {k: s & expressed for k, s in sets.items()}
    log["gsc_sets"] = {k: len(s) for k, s in sets.items()}
    for k, s in sets.items():
        pd.Series(sorted(s), name="gene").to_csv(out / f"gsc_set_{k}.txt", index=False)

    dsig = read_gmt(a.dsigdb)
    universe = expressed & set().union(*dsig.values())
    log["universe"] = len(universe)

    de = pd.read_csv(a.de, sep="\t")
    idcol = "gene_id" if "gene_id" in de.columns else de.columns[0]
    de["symbol"] = de[idcol].map(sym)
    # the published table carries ashr-shrunken log2FC and p-values but no Wald statistic: rank by signed -log10 p
    de["stat"] = np.sign(de["log2FoldChange"]) * -np.log10(de["pvalue"].clip(lower=1e-300))
    de = de.dropna(subset=["symbol", "stat"]).sort_values("stat", key=np.abs, ascending=False).drop_duplicates("symbol")

    import gseapy as gp
    rows = []
    for k, s in sets.items():
        # S1: over-representation among expressed genes
        su = s & universe
        N, n = len(universe), len(su)
        for term, g in dsig.items():
            gu = g & universe
            K, x = len(gu), len(gu & su)
            if K < 5 or x < 2:
                continue
            p = stats.hypergeom.sf(x - 1, N, K, n)
            orr = (x * (N - K - n + x)) / max(1e-9, (K - x) * (n - x))
            rows.append({"gsc_set": k, "screen": "S1_ORA", "term": term, "drug": drug_name(term), "overlap": x,
                         "set_size": K, "effect": np.log2(max(orr, 1e-9)), "p": p,
                         "genes": ",".join(sorted(gu & su)[:30])})
        # S2: recurrence ranking within the GSC genes
        rnk = de[de.symbol.isin(s)].set_index("symbol")["stat"].sort_values(ascending=False)
        log.setdefault("s2_ranked_genes", {})[k] = int(len(rnk))
        if len(rnk) >= 30:
            res = gp.prerank(rnk=rnk, gene_sets={t: list(g) for t, g in dsig.items()}, min_size=5, max_size=500,
                             permutation_num=1000, seed=1234, threads=4, outdir=None, verbose=False).res2d
            for _, r in res.iterrows():
                rows.append({"gsc_set": k, "screen": "S2_GSEA_recurrence", "term": r["Term"], "drug": drug_name(r["Term"]),
                             "overlap": int(str(r["Tag %"]).split("/")[0]) if "/" in str(r["Tag %"]) else np.nan,
                             "set_size": np.nan, "effect": float(r["NES"]), "p": float(r["NOM p-val"]),
                             "q_gseapy": float(r["FDR q-val"]), "genes": str(r["Lead_genes"])[:300]})
    res = pd.DataFrame(rows)
    res["q"] = np.nan
    for (k, sc), g in res.groupby(["gsc_set", "screen"]):
        res.loc[g.index, "q"] = g["q_gseapy"] if sc.startswith("S2") else bh(g["p"])
    res.to_csv(out / "gsc_screen_all.tsv", sep="\t", index=False)

    hits = res[res.q < 0.25].copy()
    cache_f = out / "chembl_cache.json"
    cache = json.load(open(cache_f)) if cache_f.exists() else {}
    info = {d: chembl(d, cache) for d in sorted(hits.drug.unique())}
    json.dump(cache, open(cache_f, "w"), indent=1, default=str)
    for c in ("chembl_id", "max_phase", "smiles"):
        hits[c] = hits.drug.map(lambda d: info[d].get(c))
    smi = sorted({s for s in hits.smiles if isinstance(s, str) and s})
    if smi:
        from admet_ai import ADMETModel
        pred = pd.DataFrame(ADMETModel().predict(smiles=smi))
        hits["BBB_Martins"] = hits.smiles.map(pred["BBB_Martins"].to_dict())
    hits["phase"] = pd.to_numeric(hits.max_phase, errors="coerce")
    hits["score"] = hits.effect.abs() ** 1.5 * hits.BBB_Martins
    clin = hits[hits.phase.between(1, 4) & hits.BBB_Martins.notna()].sort_values("score", ascending=False)
    hits.to_csv(out / "gsc_screen_hits.tsv", sep="\t", index=False)
    clin.to_csv(out / "gsc_screen_clinical_ranked.tsv", sep="\t", index=False)
    log["hits_q25"] = int(len(hits)); log["clinical_with_bbb"] = int(clin.drug.nunique())
    json.dump(log, open(out / "gsc_screen_log.json", "w"), indent=1)
    pd.set_option("display.width", 200)
    print(json.dumps(log, indent=1))
    for (k, sc), g in clin.groupby(["gsc_set", "screen"]):
        print(f"\n=== {k} / {sc}")
        print(g.drop_duplicates("drug").head(15)[["drug", "phase", "effect", "q", "overlap", "BBB_Martins", "score", "genes"]]
              .to_string(index=False, max_colwidth=60))


if __name__ == "__main__":
    main()
