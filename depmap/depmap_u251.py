# -*- coding: utf-8 -*-
"""Does U-251 MG actually depend on the targets this screen nominates?

DepMap 24Q4 genome-wide CRISPR knockout (Chronos gene effect). U-251 MG is
ACH-000232. Chronos scale: 0 = no effect, -1 = median common-essential gene.
A dependency is conventionally called at < -0.5.

The point of separating "common essential" from selective is that ribosome and
proteasome genes score strongly in every line and are not vulnerabilities of
this tumour in particular.
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
U251 = "ACH-000232"

eff = pd.read_csv(os.path.join(HERE, "CRISPRGeneEffect.csv"), index_col=0)
eff.columns = [c.split(" (")[0] for c in eff.columns]
print("DepMap 24Q4: %d cell lines x %d genes" % eff.shape)
assert U251 in eff.index, "U-251 MG not in the matrix"

ce = pd.read_csv(os.path.join(HERE, "CommonEssentials.csv"))
common = {str(g).split(" (")[0] for g in ce.iloc[:, 0]}
print("common-essential genes: %d\n" % len(common))

u = eff.loc[U251]
med = eff.median(axis=0)
sel = u - med                      # negative => U251 more dependent than typical

# ---------------------------------------------------------------- targets ----
TARGETS = {
    "ciclopirox: ribonucleotide reductase": ["RRM1", "RRM2", "RRM2B"],
    "ciclopirox: eIF5A hypusination":       ["DOHH", "DHPS", "EIF5A"],
    "ciclopirox/DMOG: PHD-HIF axis":        ["EGLN1", "EGLN2", "EGLN3", "HIF1A",
                                             "EPAS1", "ARNT", "VHL"],
    "iron handling / ferroptosis":          ["TFRC", "FTH1", "FTL", "SLC7A11",
                                             "GPX4", "ACSL4", "SLC40A1"],
    "MTAP-PRMT5 synthetic lethality":       ["MTAP", "PRMT5", "MAT2A", "RIOK1",
                                             "WDR77", "CDKN2A"],
    "integrated stress response":           ["EIF2AK4", "ATF4", "EIF2AK1",
                                             "DELE1", "HRI", "DDIT3"],
    "mTORC1 / biosynthesis":                ["MTOR", "RPTOR", "RICTOR", "MYC"],
    "other screen-hit drug targets":        ["KCNH1", "AR", "DRD2", "PIK3CA",
                                             "PIK3CB", "GRK2", "ADRBK1",
                                             "HSPB1", "TSPO"],
}

print("%-46s %8s %8s %8s  %s" % ("gene (pathway)", "U251", "median",
                                 "selective", "call"))
print("-" * 92)
rows = []
for label, genes in TARGETS.items():
    print("\n== %s" % label)
    for g in genes:
        if g not in eff.columns:
            print("   %-43s %s" % (g, "not in DepMap"))
            continue
        gu, gm, gs = u[g], med[g], sel[g]
        tag = []
        if gu < -0.5:
            tag.append("DEPENDENCY")
        if g in common:
            tag.append("common-essential")
        elif gu < -0.5:
            tag.append("selective")
        if gs < -0.25 and gu < -0.3:
            tag.append("U251-enriched")
        print("   %-43s %8.3f %8.3f %8.3f  %s"
              % (g, gu, gm, gs, ", ".join(tag) if tag else "no effect"))
        rows.append({"pathway": label, "gene": g, "u251_gene_effect": gu,
                     "median_across_lines": gm, "selectivity": gs,
                     "dependency": gu < -0.5,
                     "common_essential": g in common})

out = pd.DataFrame(rows)
out.to_csv(os.path.join(HERE, "depmap_u251_targets.csv"), index=False)
print("\nwrote depmap_u251_targets.csv")

# ------------------------------------------------- selective dependencies ----
cand = u[(u < -0.5) & (~u.index.isin(common))].sort_values()
print("\nU251 dependencies that are NOT common-essential: %d" % len(cand))
print("strongest 25:")
for g, v in cand.head(25).items():
    print("   %-14s %7.3f   (median %6.3f, selectivity %6.3f)"
          % (g, v, med[g], sel[g]))
cand.to_frame("u251_gene_effect").assign(
    median_across_lines=med.reindex(cand.index),
    selectivity=sel.reindex(cand.index)).to_csv(
        os.path.join(HERE, "depmap_u251_selective.csv"))
print("\nwrote depmap_u251_selective.csv")
