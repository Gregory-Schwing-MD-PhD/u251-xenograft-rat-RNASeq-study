# -*- coding: utf-8 -*-
"""Per-group means/SDs, Cohen's d, Welch-t sample size, and per-sample relation to xengsort graft %.
numpy/scipy only (pandas is broken in this Python)."""
import csv, json, math, sys
from pathlib import Path
import numpy as np
from scipy import stats

R = Path("C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study")
A = R / "ANALYSIS"
PRI = ["IL67B", "IL68B", "IL69B"]; REC = ["IL66B", "NL70B", "NL71B"]; SIX = PRI + REC


def read_tsv(p, sep="\t"):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=sep))


# ---- metrics per sample
metrics = {}
cib = {r["sample"]: r for r in read_tsv(A / "cibersort/results/v104/fractions_v104_neftel4_confident.tsv")}
metrics["CIBERSORT MES fraction (neftel4_confident)"] = {s: float(cib[s]["MES"]) for s in SIX}
metrics["CIBERSORT AC fraction (neftel4_confident)"] = {s: float(cib[s]["AC"]) for s in SIX}
metrics["CIBERSORT fit correlation (neftel4_confident)"] = {s: float(cib[s]["Correlation"]) for s in SIX}
gs = {r["Term"]: r for r in read_tsv(R / "subtypes/subtype_rerun_scores_per_sample.csv", sep=",")}
metrics["GSVA Neftel_AC"] = {s: float(gs["Neftel_AC"][s]) for s in SIX}
metrics["GSVA Neftel_MES1"] = {s: float(gs["Neftel_MES1"][s]) for s in SIX}

# translation initiation: per-sample mean of gene-centred rlog over set members (all 80, and core enrichment only)
sym2id = {}
with open(A / "results_human_final/star_salmon/salmon.merged.gene_counts.tsv", encoding="utf-8") as f:
    rd = csv.reader(f, delimiter="\t"); next(rd)
    for row in rd:
        sym2id.setdefault(row[1], row[0])
rlog = {}
with open(A / "results_therapy_v3/tables/processed_abundance/all.rlog.tsv", encoding="utf-8") as f:
    rd = csv.reader(f, delimiter="\t"); hdr = next(rd)
    idx = [hdr.index(s) for s in SIX]
    for row in rd:
        rlog[row[0]] = np.array([float(row[i]) for i in idx])
gsea_dir = A / "results_therapy_v3/report/gsea/therapy_impact/combined_human"
ti = read_tsv(gsea_dir / "therapy_impact.combined_human.KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION.tsv")
allm = [r["SYMBOL"] for r in ti]; core = [r["SYMBOL"] for r in ti if r["CORE ENRICHMENT"].strip() == "Yes"]


def set_score(symbols):
    ids = [sym2id[s] for s in symbols if s in sym2id and sym2id[s] in rlog]
    M = np.vstack([rlog[i] for i in ids])
    M = M - M.mean(1, keepdims=True)  # centre each gene across the six
    return dict(zip(SIX, M.mean(0))), len(ids)


sc_all, n_all = set_score(allm); sc_core, n_core = set_score(core)
metrics[f"Translation initiation set score, all {n_all} members (mean centred rlog)"] = sc_all
metrics[f"Translation initiation set score, {n_core} core-enrichment genes (mean centred rlog)"] = sc_core

meta = {r["sample"]: r for r in read_tsv(A / "metadata_full.csv", sep=",")}
graft = {s: float(meta[s]["graft_pct"]) for s in SIX}
host = {s: float(meta[s]["host_pct"]) for s in SIX}
both = {s: float(meta[s]["both_pct"]) for s in SIX}
mshare = {s: both[s] / (graft[s] + both[s]) for s in SIX}


# ---- power
def welch_power(delta, s1, s2, n, alpha=0.05):
    v = s1 ** 2 / n + s2 ** 2 / n
    df = v ** 2 / ((s1 ** 2 / n) ** 2 / (n - 1) + (s2 ** 2 / n) ** 2 / (n - 1))
    lam = abs(delta) / math.sqrt(v)
    tc = stats.t.ppf(1 - alpha / 2, df)
    return stats.nct.sf(tc, df, lam) + stats.nct.cdf(-tc, df, lam)


def n_for_power(delta, s1, s2, target=0.8, alpha=0.05, nmax=100000):
    n = 2
    while n <= nmax:
        if welch_power(delta, s1, s2, n, alpha) >= target:
            return n
        n += 1
    return None


def n_normal(d, alpha=0.05, power=0.8):
    z = stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)
    return 2 * z ** 2 / d ** 2


rows = []
for name, v in metrics.items():
    p = np.array([v[s] for s in PRI]); r = np.array([v[s] for s in REC])
    mp, mr = p.mean(), r.mean(); sp, sr = p.std(ddof=1), r.std(ddof=1)
    delta = mr - mp
    spool = math.sqrt((sp ** 2 + sr ** 2) / 2)
    d = delta / spool
    J = 1 - 3 / (4 * (len(p) + len(r) - 2) - 1)  # Hedges small-sample correction, df = 4 -> 0.8
    wt = stats.ttest_ind(r, p, equal_var=False)
    out = dict(metric=name, primary_mean=mp, primary_sd=sp, recurrent_mean=mr, recurrent_sd=sr, diff_rec_minus_pri=delta,
               cohens_d=d, hedges_g=d * J, welch_p_observed=wt.pvalue,
               n_per_group_d=n_for_power(delta, sp, sr), n_per_group_half_d=n_for_power(delta / 2, sp, sr),
               n_per_group_g=n_for_power(delta * J, sp, sr),
               n_normal_approx_d=n_normal(d), n_normal_approx_half_d=n_normal(d / 2))
    # relation to graft %
    x = np.array([graft[s] for s in SIX]); y = np.array([v[s] for s in SIX])
    out["pearson_r_vs_graft"] = stats.pearsonr(x, y)[0]; out["pearson_p_vs_graft"] = stats.pearsonr(x, y)[1]
    out["spearman_rho_vs_graft"] = stats.spearmanr(x, y)[0]; out["spearman_p_vs_graft"] = stats.spearmanr(x, y)[1]
    # within-group: does graft % still track the metric after removing the group mean?
    xc = np.concatenate([x[:3] - x[:3].mean(), x[3:] - x[3:].mean()]); yc = np.concatenate([y[:3] - y[:3].mean(), y[3:] - y[3:].mean()])
    out["within_group_r_vs_graft"] = float(np.corrcoef(xc, yc)[0, 1])
    rows.append(out)

res = {"samples": {s: dict(group="Primary" if s in PRI else "Recurrent", graft_pct=graft[s], host_pct=host[s], both_pct=both[s],
                           both_share_of_human_stream=mshare[s], **{k: v[s] for k, v in metrics.items()}) for s in SIX},
       "graft_by_group": {}, "rows": rows, "n_core": n_core, "n_all": n_all}
gp = np.array([graft[s] for s in PRI]); gr = np.array([graft[s] for s in REC])
res["graft_by_group"] = dict(primary_mean=gp.mean(), primary_sd=gp.std(ddof=1), recurrent_mean=gr.mean(), recurrent_sd=gr.std(ddof=1),
                             welch_p=stats.ttest_ind(gr, gp, equal_var=False).pvalue,
                             mannwhitney_p_exact=stats.mannwhitneyu(gr, gp, alternative="two-sided", method="exact").pvalue)
out_p = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("effects.json")
out_p.write_text(json.dumps(res, indent=1, default=float))
for s in SIX:
    print(s, res["samples"][s]["group"], f"graft {graft[s]:.2f} host {host[s]:.2f} both {both[s]:.2f}")
print("graft by group", {k: round(v, 4) for k, v in res["graft_by_group"].items()})
for o in rows:
    print(f"\n{o['metric']}\n  P {o['primary_mean']:.4f}+/-{o['primary_sd']:.4f}  R {o['recurrent_mean']:.4f}+/-{o['recurrent_sd']:.4f}  "
          f"diff {o['diff_rec_minus_pri']:+.4f}  d {o['cohens_d']:+.3f}  g {o['hedges_g']:+.3f}  welch p {o['welch_p_observed']:.4f}\n"
          f"  n/group: d {o['n_per_group_d']}  g {o['n_per_group_g']}  d/2 {o['n_per_group_half_d']}  (normal approx {o['n_normal_approx_d']:.1f} / {o['n_normal_approx_half_d']:.1f})\n"
          f"  vs graft%: pearson {o['pearson_r_vs_graft']:+.3f} (p {o['pearson_p_vs_graft']:.3f})  spearman {o['spearman_rho_vs_graft']:+.3f} (p {o['spearman_p_vs_graft']:.3f})  within-group r {o['within_group_r_vs_graft']:+.3f}")
