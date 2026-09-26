#!/usr/bin/env python
"""Collect the leave-one-tumour-out GSEA runs (gsea_array.sbatch) into tables and SUMMARY.md.

1. Reproduction: the full six-tumour run with seed 1234 must match the nf-core report (inputs/original) set for set;
   the largest absolute differences in NES, nominal p and FDR q are written to reproduction.json. If they are not zero
   the other runs cannot be read as the same analysis, and SUMMARY.md says so first.
2. For every run: NES, nominal p, FDR q and FWER p of the leading and themed sets (loo_sets.tsv), and the shape of the
   screen (sets at q < 0.05 and q < 0.25 on each side, best q on each side, sets at nominal p < 0.05; loo_screen.tsv).
3. The ten lowest-q down sets of each seed-1234 run (loo_top_down.tsv), to see what leads without each tumour.

    python ANALYSIS/gsea_leave_one_out/collect.py
"""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

U = Path(os.environ.get("U251_DIR", "/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study"))
D = U / "ANALYSIS" / "gsea_leave_one_out"
LEAD = ["KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
        "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "KEGG_RIBOSOME", "REACTOME_SELENOAMINO_ACID_METABOLISM",
        "REACTOME_CELLULAR_RESPONSE_TO_STARVATION"]
THEMED = ["GOBP_RIBOSOME_BIOGENESIS", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_GLYCOLYSIS",
          "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE", "SEMENZA_HIF1_TARGETS",
          "BUFFA_HYPOXIA_METAGENE", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS",
          "REACTOME_CHOLESTEROL_BIOSYNTHESIS"]
NUM = ["SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"]


def read_pair(dn_path, up_path):
    out = []
    for p, side in ((dn_path, "down"), (up_path, "up")):
        t = pd.read_csv(p, sep="\t")
        t = t[[c for c in t.columns if not c.startswith("Unnamed")]]
        for c in NUM:
            t[c] = pd.to_numeric(t[c], errors="coerce")
        out.append(t.assign(side=side))
    return pd.concat(out, ignore_index=True).set_index("NAME")


def run_report(cond, seed):
    runs = sorted((D / "runs" / f"{cond}_s{seed}").glob("therapy_impact.combined_human.Gsea.*"))
    if not runs:
        return None
    r = runs[-1]
    dn = sorted(r.glob("gsea_report_for_Primary_U2_*.tsv"))
    up = sorted(r.glob("gsea_report_for_Recurrent_U2_*.tsv"))
    if not dn or not up:
        return None
    return read_pair(dn[-1], up[-1])


tasks = pd.read_csv(D / "tasks.tsv", sep="\t")
orig = read_pair(D / "inputs" / "original" / "therapy_impact.combined_human.gsea_report_for_Primary_U2.tsv",
                 D / "inputs" / "original" / "therapy_impact.combined_human.gsea_report_for_Recurrent_U2.tsv")

sets_rows, screen_rows, top_rows, missing = [], [], [], []
reports = {}
for _, t in tasks.iterrows():
    rep = run_report(t.condition, t.seed)
    if rep is None:
        missing.append(f"{t.condition}_s{t.seed}")
        continue
    reports[(t.condition, t.seed)] = rep
    dn, up = rep[rep.side == "down"], rep[rep.side == "up"]
    screen_rows.append({"condition": t.condition, "seed": t.seed, "n_genes": t.n_genes, "n_sets": len(rep),
                        "down_q05": int((dn["FDR q-val"] < 0.05).sum()), "down_q25": int((dn["FDR q-val"] < 0.25).sum()),
                        "up_q05": int((up["FDR q-val"] < 0.05).sum()), "up_q25": int((up["FDR q-val"] < 0.25).sum()),
                        "best_down_q": float(dn["FDR q-val"].min()), "best_up_q": float(up["FDR q-val"].min()),
                        "down_p05": int((dn["NOM p-val"] < 0.05).sum()), "up_p05": int((up["NOM p-val"] < 0.05).sum()),
                        "best_down_set": dn["FDR q-val"].idxmin(), "best_up_set": up["FDR q-val"].idxmin()})
    qrank = rep["FDR q-val"].rank(method="min")
    for s in LEAD + THEMED:
        if s in rep.index:
            r = rep.loc[s]
            sets_rows.append({"condition": t.condition, "seed": t.seed, "set": s, "lead": s in LEAD, "side": r.side,
                              "size": r.SIZE, "NES": r.NES, "nom_p": r["NOM p-val"], "fdr_q": r["FDR q-val"],
                              "fwer_p": r["FWER p-val"], "rank_by_q": int(qrank[s])})
    if t.seed == 1234:
        for i, (name, r) in enumerate(dn.sort_values(["FDR q-val", "NES"]).head(10).iterrows(), 1):
            top_rows.append({"condition": t.condition, "rank": i, "set": name, "NES": r.NES, "nom_p": r["NOM p-val"], "fdr_q": r["FDR q-val"]})

sets = pd.DataFrame(sets_rows); screen = pd.DataFrame(screen_rows); top = pd.DataFrame(top_rows)
sets.to_csv(D / "loo_sets.tsv", sep="\t", index=False)
screen.to_csv(D / "loo_screen.tsv", sep="\t", index=False)
top.to_csv(D / "loo_top_down.tsv", sep="\t", index=False)

# ---- reproduction of the published run
repro = {"available": ("full", 1234) in reports}
if repro["available"]:
    f = reports[("full", 1234)]
    j = orig.join(f, lsuffix="_orig", rsuffix="_rerun", how="outer")
    repro.update({"n_sets_orig": int(len(orig)), "n_sets_rerun": int(len(f)), "n_unmatched": int(j["NES_orig"].isna().sum() + j["NES_rerun"].isna().sum())})
    for c, k in (("NES", "NES"), ("NOM p-val", "nom_p"), ("FDR q-val", "fdr_q"), ("FWER p-val", "fwer_p")):
        repro[f"max_abs_diff_{k}"] = float(np.nanmax(np.abs(j[f"{c}_orig"] - j[f"{c}_rerun"])))
    repro["exact"] = repro["n_unmatched"] == 0 and all(repro[f"max_abs_diff_{k}"] < 1e-6 for k in ("NES", "nom_p", "fdr_q", "fwer_p"))
json.dump(repro, open(D / "reproduction.json", "w"), indent=1)

# ---- SUMMARY.md
L = []
L.append("# Leave-one-tumour-out GSEA (U251 LITT recurrence), GSEA 4.3.2, the nf-core command with one tumour removed\n")
if missing:
    L.append(f"**Missing runs ({len(missing)}):** " + ", ".join(missing) + "\n")
if repro.get("available"):
    L.append(f"**Reproduction of the published run (all six, seed 1234):** exact = {repro['exact']}; sets {repro['n_sets_rerun']} vs "
             f"{repro['n_sets_orig']}, unmatched {repro['n_unmatched']}; max |diff| NES {repro['max_abs_diff_NES']:.2g}, "
             f"nominal p {repro['max_abs_diff_nom_p']:.2g}, FDR q {repro['max_abs_diff_fdr_q']:.2g}.\n")
else:
    L.append("**Reproduction run missing: nothing below can be read as the same analysis.**\n")
TI = LEAD[0]
if len(sets):
    L.append("## Translation initiation (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) by condition\n")
    L.append("| condition | genes | NES (seed 1234) | nominal p | FDR q (seed 1234) | FDR q, 5 seeds (min-max) | FWER p | rank by q |")
    L.append("|---|---|---|---|---|---|---|---|")
    for cond in dict.fromkeys(tasks.condition):
        s = sets[(sets.condition == cond) & (sets.set == TI)]
        if not len(s):
            continue
        a = s[s.seed == 1234]
        g = screen[(screen.condition == cond) & (screen.seed == 1234)]
        if not len(a):
            continue
        a = a.iloc[0]
        L.append(f"| {cond} | {int(g.n_genes.iloc[0]) if len(g) else ''} | {a.NES:.2f} | {a.nom_p:.3f} | {a.fdr_q:.3f} | "
                 f"{s.fdr_q.min():.3f}-{s.fdr_q.max():.3f} (n = {len(s)}) | {a.fwer_p:.3f} | {a.rank_by_q} |")
    L.append("\n## The six leading sets, FDR q at seed 1234 (NES in brackets)\n")
    conds = [c for c in dict.fromkeys(tasks.condition) if not c.endswith("_nofilter")]
    L.append("| set | " + " | ".join(conds) + " |")
    L.append("|---|" + "---|" * len(conds))
    for st in LEAD:
        cells = []
        for cond in conds:
            a = sets[(sets.condition == cond) & (sets.seed == 1234) & (sets.set == st)]
            cells.append(f"{a.fdr_q.iloc[0]:.3f} ({a.NES.iloc[0]:.2f})" if len(a) else "")
        L.append(f"| {st} | " + " | ".join(cells) + " |")
    L.append("\n## Shape of the screen (seed 1234; ranges over the five seeds in loo_screen.tsv)\n")
    L.append("| condition | down q<0.05 | down q<0.25 | up q<0.05 | up q<0.25 | best down q (set) | best up q (set) | down p<0.05 | up p<0.05 |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, g in screen[screen.seed == 1234].iterrows():
        L.append(f"| {g.condition} | {g.down_q05} | {g.down_q25} | {g.up_q05} | {g.up_q25} | {g.best_down_q:.3f} ({g.best_down_set}) | "
                 f"{g.best_up_q:.3f} ({g.best_up_set}) | {g.down_p05} | {g.up_p05} |")
(D / "SUMMARY.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
