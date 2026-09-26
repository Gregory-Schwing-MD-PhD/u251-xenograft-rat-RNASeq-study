#!/usr/bin/env python
"""Published vs control vs IL66B-held-out, step by step, written to COMPARISON.md and comparison.json.

control = today's rerun of the published six-tumour inputs through the same R and Python stages; it separates
drift in the software, container or web services (ChEMBL, ADMET-AI) from the effect of holding out IL66B.
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

U = Path(os.environ.get("U251_DIR", "/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study"))
H = U / "ANALYSIS" / "holdout_IL66B"
PUB_RES = U / "ANALYSIS" / "results_therapy_v3"
HO_RES = H / "results_therapy_v3_noIL66B"
PUB = H / "published"                     # the September subtypes/ outputs, copied in by the submit script
LEAD = ["KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
        "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "KEGG_RIBOSOME", "REACTOME_SELENOAMINO_ACID_METABOLISM",
        "REACTOME_CELLULAR_RESPONSE_TO_STARVATION"]
out, md = {}, ["# IL66B held out: the full pipeline against the published run\n",
               "control = the published six-tumour inputs re-run today through the same R and Python stages; holdout = IL66B "
               "removed from the sample sheet and every stage re-run from DESeq2.\n"]


def jacc(a, b):
    a, b = set(a), set(b)
    return len(a & b) / max(1, len(a | b))


def safe(title, fn):
    try:
        fn()
    except Exception as e:                                 # noqa: BLE001
        md.append(f"\n## {title}\n\n**Not compared:** {type(e).__name__}: {e}\n")
        out[title] = {"error": str(e)}


# ---------------------------------------------------------------- 1. differential expression
def de():
    p = pd.read_csv(PUB_RES / "tables/differential/therapy_impact.deseq2.results.tsv", sep="\t")
    h = pd.read_csv(HO_RES / "tables/differential/therapy_impact.deseq2.results.tsv", sep="\t")
    idc = "gene_id" if "gene_id" in p.columns else p.columns[0]
    sym = pd.read_excel(U / "publication_figure" / "Supplementary_Data.xlsx", "S1_DE_all")[["gene_id", "symbol"]].set_index("gene_id")["symbol"]
    sig = lambda d: set(d[(d.padj < 0.05) & (d.log2FoldChange.abs() > 1)][idc])  # noqa: E731
    sp, sh = sig(p), sig(h)
    name = lambda ids: sorted(str(sym.get(i, i)) for i in ids)  # noqa: E731
    m = p.set_index(idc)[["log2FoldChange"]].join(h.set_index(idc)[["log2FoldChange"]], lsuffix="_pub", rsuffix="_ho", how="inner").dropna()
    rho = float(m.corr(method="spearman").iloc[0, 1])
    up = lambda d, s: int((d[d[idc].isin(s)].log2FoldChange > 0).sum())  # noqa: E731
    out["de"] = {"published_n": len(sp), "holdout_n": len(sh), "shared": len(sp & sh), "jaccard": jacc(sp, sh),
                 "published_up": up(p, sp), "holdout_up": up(h, sh), "spearman_lfc_all_genes": rho,
                 "lost": name(sp - sh), "gained": name(sh - sp)}
    d = out["de"]
    md.append("\n## 1. Differential expression (DESeq2, FDR < 0.05, |log2FC| > 1)\n")
    md.append(f"| | published | held out |\n|---|---|---|\n| genes | {d['published_n']} ({d['published_up']} up) | {d['holdout_n']} ({d['holdout_up']} up) |\n"
              f"| shared | {d['shared']} (Jaccard {d['jaccard']:.2f}) | |\n| Spearman of log2FC, all genes | {rho:.3f} | |\n")
    md.append(f"\nLost without IL66B ({len(d['lost'])}): {', '.join(d['lost'][:60])}\n\nGained ({len(d['gained'])}): {', '.join(d['gained'][:60])}\n")


# ---------------------------------------------------------------- 2. Broad GSEA (the pipeline's)
def gsea_pair(res):
    g = res / "report/gsea/therapy_impact/combined_human"
    t = []
    for side, f in (("down", "therapy_impact.combined_human.gsea_report_for_Primary_U2.tsv"), ("up", "therapy_impact.combined_human.gsea_report_for_Recurrent_U2.tsv")):
        d = pd.read_csv(g / f, sep="\t")
        for c in ("NES", "NOM p-val", "FDR q-val"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
        t.append(d.assign(side=side))
    return pd.concat(t).set_index("NAME")


def gsea():
    p, h = gsea_pair(PUB_RES), gsea_pair(HO_RES)
    rows = []
    for s in LEAD:
        rows.append({"set": s, "NES_pub": p.loc[s, "NES"], "q_pub": p.loc[s, "FDR q-val"], "NES_ho": h.loc[s, "NES"] if s in h.index else np.nan,
                     "p_ho": h.loc[s, "NOM p-val"] if s in h.index else np.nan, "q_ho": h.loc[s, "FDR q-val"] if s in h.index else np.nan})
    cnt = lambda d: {"down_q05": int(((d.side == "down") & (d["FDR q-val"] < 0.05)).sum()), "down_q25": int(((d.side == "down") & (d["FDR q-val"] < 0.25)).sum()),  # noqa: E731
                     "up_q25": int(((d.side == "up") & (d["FDR q-val"] < 0.25)).sum())}
    out["gsea"] = {"lead": rows, "published": cnt(p), "holdout": cnt(h)}
    md.append("\n## 2. Gene-set enrichment (the pipeline's Broad GSEA, seed 1234)\n")
    md.append("| set | NES published | q published | NES held out | p held out | q held out |\n|---|---|---|---|---|---|")
    for r in rows:
        md.append(f"| {r['set']} | {r['NES_pub']:.2f} | {r['q_pub']:.3f} | {r['NES_ho']:.2f} | {r['p_ho']:.3f} | {r['q_ho']:.3f} |")
    md.append(f"\nSets at q < 0.05 / q < 0.25 (down): published {out['gsea']['published']['down_q05']} / {out['gsea']['published']['down_q25']}; "
              f"held out {out['gsea']['holdout']['down_q05']} / {out['gsea']['holdout']['down_q25']}. Up side at q < 0.25: {out['gsea']['published']['up_q25']} and {out['gsea']['holdout']['up_q25']}.\n")


# ---------------------------------------------------------------- 3. the DSigDB drug screen (R stage)
def drug_profiles():
    arms = {"published": PUB / "therapy_impact_Drug_Profiles_Comprehensive.csv",
            "control": H / "control/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv",
            "holdout": H / "holdout/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv"}
    d = {k: pd.read_csv(v) for k, v in arms.items()}
    names = {k: set(v.Drug.astype(str).str.lower()) for k, v in d.items()}
    cpx = {k: v[v.Drug.astype(str).str.lower() == "ciclopirox"] for k, v in d.items()}
    out["drug_screen"] = {"n": {k: len(v) for k, v in d.items()},
                          "jaccard_control_vs_published": jacc(names["control"], names["published"]),
                          "jaccard_holdout_vs_published": jacc(names["holdout"], names["published"]),
                          "jaccard_holdout_vs_control": jacc(names["holdout"], names["control"]),
                          "ciclopirox": {k: (None if not len(v) else {"NES": float(v.NES.iloc[0]), "FDR": float(v.FDR.iloc[0]), "rank": int(v.Rank.iloc[0])}) for k, v in cpx.items()},
                          "entered_holdout": sorted(names["holdout"] - names["control"]), "left_holdout": sorted(names["control"] - names["holdout"])}
    s = out["drug_screen"]
    md.append("\n## 3. DSigDB drug screen: the 100 most opposing drug gene sets (R stage)\n")
    md.append(f"Overlap of the 100 compounds (Jaccard): control vs published {s['jaccard_control_vs_published']:.2f}; held out vs published "
              f"{s['jaccard_holdout_vs_published']:.2f}; held out vs control {s['jaccard_holdout_vs_control']:.2f}.\n")
    md.append("| arm | ciclopirox NES | FDR | rank (R integrated score) |\n|---|---|---|---|")
    for k, v in s["ciclopirox"].items():
        md.append(f"| {k} | " + (f"{v['NES']:.3f} | {v['FDR']:.2g} | {v['rank']} |" if v else "absent | | |"))
    md.append(f"\nEntered with IL66B held out ({len(s['entered_holdout'])}): {', '.join(s['entered_holdout'][:40])}\n\nLeft ({len(s['left_holdout'])}): {', '.join(s['left_holdout'][:40])}\n")


# ---------------------------------------------------------------- 4. the final ranking (Python stage)
def final_ranking():
    arms = {"published": PUB / "drug_ranking_final.csv", "control": H / "control/subtypes/drug_ranking_final.csv", "holdout": H / "holdout/subtypes/drug_ranking_final.csv"}
    d = {k: pd.read_csv(v) for k, v in arms.items()}
    for v in d.values():
        v["key"] = v.Drug.astype(str).str.lower()
    top = {k: list(v.sort_values("rank_bbb").key.head(20)) for k, v in d.items()}
    both = {k: set(v[v.both_agree.astype(str) == "True"].key) for k, v in d.items()}
    cpx = {k: v[v.key == "ciclopirox"] for k, v in d.items()}
    out["ranking"] = {"n_clinical": {k: len(v) for k, v in d.items()}, "n_both_bbb": {k: len(s) for k, s in both.items()},
                      "top20": top, "top20_jaccard_holdout_vs_control": jacc(top["holdout"], top["control"]),
                      "top20_jaccard_control_vs_published": jacc(top["control"], top["published"]),
                      "ciclopirox": {k: (None if not len(v) else {"rank_bbb": int(v.rank_bbb.iloc[0]), "rank_nobbb": int(v.rank_nobbb.iloc[0]),
                                                                   "score_bbb": float(v.score_bbb.iloc[0]), "both_agree": str(v.both_agree.iloc[0])}) for k, v in cpx.items()}}
    r = out["ranking"]
    md.append("\n## 4. Final ranking (clinical phase, ADMET-AI, BOILED-Egg)\n")
    md.append("| arm | clinically available | clear both barrier models | ciclopirox rank (weighted / unweighted) | ciclopirox score |\n|---|---|---|---|---|")
    for k in d:
        c = r["ciclopirox"][k]
        md.append(f"| {k} | {r['n_clinical'][k]} | {r['n_both_bbb'][k]} | " + (f"{c['rank_bbb']} / {c['rank_nobbb']} | {c['score_bbb']:.2f} |" if c else "absent | |"))
    md.append(f"\nTop 20 overlap (Jaccard): control vs published {r['top20_jaccard_control_vs_published']:.2f}; held out vs control {r['top20_jaccard_holdout_vs_control']:.2f}.\n")
    md.append("| rank | published | control | held out |\n|---|---|---|---|")
    for i in range(20):
        md.append(f"| {i + 1} | " + " | ".join(top[k][i] if i < len(top[k]) else "" for k in ("published", "control", "holdout")) + " |")


# ---------------------------------------------------------------- 5. subtypes (GSVA)
def subtypes():
    arms = {"published": PUB / "subtype_rerun_results.csv", "control": H / "control/subtypes/subtype_rerun_results.csv", "holdout": H / "holdout/subtypes/subtype_rerun_results.csv"}
    d = {k: pd.read_csv(v).set_index("Signature") for k, v in arms.items()}
    out["subtypes"] = {k: {s: {"change": float(v.loc[s, "Change"]), "p": float(v.loc[s, "p"])} for s in v.index} for k, v in d.items()}
    md.append("\n## 5. Subtype signatures (GSVA; change recurrent minus primary, Welch p)\n")
    md.append("| signature | published | control | held out |\n|---|---|---|---|")
    for s in d["published"].index:
        md.append(f"| {s} | " + " | ".join(f"{d[k].loc[s, 'Change']:+.2f} (p {d[k].loc[s, 'p']:.3f})" if s in d[k].index else "" for k in ("published", "control", "holdout")) + " |")


# ---------------------------------------------------------------- 6. STRING hubs
def hubs():
    arms = {"published": U / "publication_figure" / "Supplementary_Data.xlsx", "control": H / "control/publication_figure/Supplementary_Data.xlsx",
            "holdout": H / "holdout/publication_figure/Supplementary_Data.xlsx"}
    res = {}
    for k, v in arms.items():
        e = pd.read_excel(v, "S6_PPI_edges").drop_duplicates(["gene1", "gene2"])
        deg = pd.concat([e.gene1, e.gene2]).value_counts()
        res[k] = {"nodes": int(deg.size), "edges": int(len(e)), "hubs_deg3": sorted(deg[deg >= 3].index.astype(str))}
    out["hubs"] = res
    md.append("\n## 6. STRING network among the DE genes\n")
    md.append("| arm | nodes | edges | hubs (degree ≥ 3) |\n|---|---|---|---|")
    for k, v in res.items():
        md.append(f"| {k} | {v['nodes']} | {v['edges']} | {', '.join(v['hubs_deg3'])} |")


for title, fn in (("de", de), ("gsea", gsea), ("drug_screen", drug_profiles), ("ranking", final_ranking), ("subtypes", subtypes), ("hubs", hubs)):
    safe(title, fn)
md.append("\nNot rerun: the prior-art audit (a manual literature review of the top 20; compounds entering the held-out top 20 "
          "would need it) and DepMap (a property of the cell line, not of the samples).\n")
json.dump(out, open(H / "comparison.json", "w"), indent=1, default=str)
(H / "COMPARISON.md").write_text("\n".join(md) + "\n")
print("\n".join(md))
