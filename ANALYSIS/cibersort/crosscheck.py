"""Cross-check the deconvolution against the September GSVA scores of the same Neftel programs (six tumours), and
summarise primary vs recurrent per state. GSVA is an independent method on the same reads; agreement in ranking is
the outside check the fractions need."""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

R = Path(__file__).resolve().parent / "results"
G = pd.read_csv(Path(__file__).resolve().parents[2] / "subtypes" / "subtype_rerun_scores_per_sample.csv", index_col=0)
gs = pd.DataFrame({"MES": G.loc[["Neftel_MES1", "Neftel_MES2"]].mean(), "AC": G.loc["Neftel_AC"],
                   "OPC": G.loc["Neftel_OPC"], "NPC": G.loc[["Neftel_NPC1", "Neftel_NPC2"]].mean()})
tum = list(G.columns)
for f in sorted(R.glob("fractions_nusvr_*.tsv")):
    d = pd.read_csv(f, sep="\t").set_index("sample")
    print("\n===", f.stem.replace("fractions_nusvr_", ""))
    cols = [c for c in ["MES", "AC", "OPC", "NPC", "Macrophage", "Oligodendrocyte", "T-cell"] if c in d.columns]
    print(d[["cohort"] + cols + ["R", "RMSE", "P"]].round(3).to_string())
    t = d.loc[tum]
    prim = t.cohort.str.startswith("Primary")
    print("state   primary mean  recurrent mean  diff   Welch p   Spearman vs GSVA (6 tumours)")
    for s in ["MES", "AC", "OPC", "NPC"]:
        a, b = t.loc[prim, s], t.loc[~prim, s]
        p = stats.ttest_ind(b, a, equal_var=False).pvalue if (a.std() + b.std()) > 0 else float("nan")
        rho = stats.spearmanr(t[s], gs.loc[tum, s]).correlation if t[s].std() > 0 else float("nan")
        print(f"{s:6s}  {a.mean():12.3f}  {b.mean():14.3f}  {b.mean() - a.mean():+.3f}  {p:8.3f}   {rho:+.2f}")
    held = [x for x in tum if x != "IL68B"]
    th = t.loc[held]
    ph = th.cohort.str.startswith("Primary")
    print("without IL68B: " + ", ".join(f"{s} {th.loc[~ph, s].mean() - th.loc[ph, s].mean():+.3f}" for s in ["MES", "AC", "OPC", "NPC"]))
