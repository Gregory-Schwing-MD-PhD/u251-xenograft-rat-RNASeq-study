"""CIBERSORT v1.04 against the Python re-implementation (same signature, same mixture, QN off), then the
primary-vs-recurrent summary and the rank agreement with the September GSVA scores, from the v1.04 numbers.

    python compare_v104.py            (reads results/v104/fractions_v104_*.tsv and results/fractions_nusvr_*.tsv)
"""
from pathlib import Path

import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
COH = pd.read_csv(HERE.parents[0] / "sample_cohorts.csv").set_index("sample")["cohort"].to_dict()
G = pd.read_csv(HERE.parents[1] / "subtypes" / "subtype_rerun_scores_per_sample.csv", index_col=0)
GS = pd.DataFrame({"MES": G.loc[["Neftel_MES1", "Neftel_MES2"]].mean(), "AC": G.loc["Neftel_AC"],
                   "OPC": G.loc["Neftel_OPC"], "NPC": G.loc[["Neftel_NPC1", "Neftel_NPC2"]].mean()})
TUM = list(G.columns)

nus = {}
for f in RES.glob("fractions_nusvr*.tsv"):
    d = pd.read_csv(f, sep="\t")
    for sig, g in d.groupby("signature"):
        nus[sig] = g.set_index("sample")

for f in sorted((RES / "v104").glob("fractions_v104_*.tsv")):
    name = f.stem.replace("fractions_v104_", "")
    v = pd.read_csv(f, sep="\t").set_index("sample")
    types = [c for c in v.columns if c not in ("P-value", "Correlation", "RMSE", "Absolute score (sig.score)")]
    print(f"\n=== {name}  (CIBERSORT v1.04, perm from the job, QN off)")
    out = v[types + ["P-value", "Correlation", "RMSE"]].copy()
    out.insert(0, "cohort", [COH.get(s, "") for s in out.index])
    print(out.round(3).to_string())
    if name in nus:
        n = nus[name].loc[v.index]
        dmax = (v[types] - n[types]).abs().max().max()
        dr = (v["Correlation"] - n["R"]).abs().max()
        drm = (v["RMSE"] - n["RMSE"]).abs().max()
        print(f"against the Python re-implementation: max |fraction difference| {dmax:.4f}; max |R difference| {dr:.4f}; max |RMSE difference| {drm:.4f}")
    t = v.loc[[s for s in TUM if s in v.index]]
    prim = pd.Series([COH[s].startswith("Primary") for s in t.index], index=t.index)
    print("state            primary  recurrent   diff   Welch p   Spearman vs GSVA   diff without IL68B")
    for s in [c for c in ["MES", "AC", "OPC", "NPC"] if c in types]:
        a, b = t.loc[prim, s], t.loc[~prim, s]
        p = stats.ttest_ind(b, a, equal_var=False).pvalue if (a.std() + b.std()) > 0 else float("nan")
        rho = stats.spearmanr(t[s], GS.loc[t.index, s]).correlation if t[s].std() > 0 else float("nan")
        th = t.drop(index="IL68B", errors="ignore")
        ph = prim.drop(index="IL68B", errors="ignore")
        dh = th.loc[~ph, s].mean() - th.loc[ph, s].mean()
        print(f"{s:15s} {a.mean():8.3f} {b.mean():10.3f} {b.mean() - a.mean():+7.3f} {p:9.3f} {rho:+12.2f} {dh:+15.3f}")
    extra = [c for c in types if c not in ("MES", "AC", "OPC", "NPC")]
    if extra:
        print("other columns, mean over the six tumours: " + ", ".join(f"{c} {t[c].mean():.3f}" for c in extra))
