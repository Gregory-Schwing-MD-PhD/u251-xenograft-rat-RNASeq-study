# -*- coding: utf-8 -*-
"""Assemble ANALYSIS/power/power_table.csv from effects.json (part a), power_rnaseq.json + dispersion.json (part b)
and the CIBERSORTx requirements quoted from Newman 2019 / Steen 2020 (part c)."""
import csv, json, math, sys
from pathlib import Path

S = Path(sys.argv[1]); OUT = Path(sys.argv[2])
eff = json.loads((S / "effects.json").read_text())
pw = json.loads((S / "power_rnaseq.json").read_text())
disp = json.loads((S / "dispersion.json").read_text())
cols = ["section", "endpoint", "method", "assumptions", "primary_mean", "primary_sd", "recurrent_mean", "recurrent_sd",
        "diff_rec_minus_pri", "cohens_d", "hedges_g", "welch_p_observed", "n_unit", "n_estimate", "n_conservative", "reference"]
rows = []
f4 = lambda x: f"{x:.4f}"  # noqa: E731
WELCH = ("Welch t, two-sided alpha 0.05, power 0.80, equal n per group; SDs = observed group SDs (unequal allowed); "
         "noncentral t with Welch-Satterthwaite df; d = diff / sqrt((s1^2+s2^2)/2); normal approx n = 2(z.975+z.80)^2/d^2 = 15.70/d^2. "
         "n_estimate at observed d; n_conservative at d/2. Hedges g = 0.8 d (df 4); n at g = ")
REF_WELCH = "Welch 1947 Biometrika 34:28 doi:10.1093/biomet/34.1-2.28; Hedges 1981 J Educ Stat 6:107 doi:10.3102/10769986006002107"
src = {"CIBERSORT": "ANALYSIS/cibersort/results/v104/fractions_v104_neftel4_confident.tsv",
       "GSVA": "subtypes/subtype_rerun_scores_per_sample.csv (GSVA across the six tumours; scores are relative to this sample set)",
       "Translation": "results_therapy_v3 all.rlog.tsv; members from the GSEA report KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION.tsv; per-sample score = mean over genes of rlog centred across the six"}
for r in eff["rows"]:
    if "fit correlation" in r["metric"]:
        continue
    key = "CIBERSORT" if r["metric"].startswith("CIBERSORT") else "GSVA" if r["metric"].startswith("GSVA") else "Translation"
    rows.append(dict(section="a_targeted_endpoint", endpoint=r["metric"], method="Welch t-test sample size",
                     assumptions=WELCH + f"{r['n_per_group_g']}. Source: {src[key]}",
                     primary_mean=f4(r["primary_mean"]), primary_sd=f4(r["primary_sd"]), recurrent_mean=f4(r["recurrent_mean"]),
                     recurrent_sd=f4(r["recurrent_sd"]), diff_rec_minus_pri=f4(r["diff_rec_minus_pri"]), cohens_d=f"{r['cohens_d']:.3f}",
                     hedges_g=f"{r['hedges_g']:.3f}", welch_p_observed=f"{r['welch_p_observed']:.4f}", n_unit="animals per group",
                     n_estimate=r["n_per_group_d"], n_conservative=r["n_per_group_half_d"], reference=REF_WELCH))

a6 = disp["all6"]
bcv_c, bcv_m = a6["common_BCV"], a6["tagwise_BCV_quantiles"]["50%"]
DISPTXT = (f"dispersion from the six tumours' human-stream counts (edgeR 4.8.2 estimateDisp, robust; {a6['n_genes']} genes after filterByExpr): "
           f"common BCV {bcv_c:.3f} (dispersion {a6['common_disp']:.4f}); tagwise BCV median {bcv_m:.3f}, IQR "
           f"{a6['tagwise_BCV_quantiles']['25%']:.3f}-{a6['tagwise_BCV_quantiles']['75%']:.3f}, 90th pct {a6['tagwise_BCV_quantiles']['90%']:.3f}; "
           f"median normalised mean count {a6['mean_norm_count_median']:.0f} at ~25 M human-stream reads/library")
ja = pw["jung_alpha"]
REF_HART = "Hart et al. 2013 J Comput Biol 20:970 doi:10.1089/cmb.2012.0283 (RNASeqPower 1.50.0); Jung 2005 Bioinformatics 21:3097 doi:10.1093/bioinformatics/bti456"
hart = pw["hart"]
for pi0key, lab in [("pi0=0.99", "pi0 0.99 (1% of genes DE)"), ("pi0=0.95", "pi0 0.95 (5% DE)")]:
    tab = {row["depth"]: row for row in hart[pi0key]}
    for depth in ["20", "463"]:
        for cv, cvlab in [("0.2", "our median tagwise BCV 0.20"), ("0.269", "our common BCV 0.27"), ("0.4", "edgeR guide 'human' BCV 0.4"), ("0.1", "edgeR guide 'genetically identical' BCV 0.1")]:
            rows.append(dict(section="b_genome_wide_DE", endpoint="per-gene 1.5-fold change at FDR 0.05, power 0.80",
                             method="Hart 2013 formula n = 2(z_{1-a/2}+z_{.80})^2 (1/depth + CV^2)/(ln 1.5)^2 with Jung 2005 per-test alpha",
                             assumptions=f"{lab}; Jung alpha* = r1 f/(m0 (1-f)) = {ja[pi0key]:.2e} with m = {pw['G']} genes, r1 = 0.8 m1, f = 0.05; depth {depth} reads/gene; CV = {cvlab}",
                             n_unit="animals per group", n_estimate=tab[depth][cv], n_conservative="", reference=REF_HART))
REF_SS = "Bi & Liu 2016 BMC Bioinformatics 17:146 doi:10.1186/s12859-016-0994-9 (ssizeRNA 1.3.3)"
for k, v in pw["ssizeRNA"].items():
    ss = v["ssize"]; n = ss[0][1]
    pw3 = dict((int(a), b) for a, b in v["power"]).get(3); pw6 = dict((int(a), b) for a, b in v["power"]).get(6)
    if k.startswith("vary"):
        p0 = k.split("_")[-1]
        a = f"gene-wise (mean, dispersion) pairs resampled from our {pw['G']} filtered genes; pi0 {p0}; 1.5-fold, half up/half down; FDR 0.05; {DISPTXT}"
    else:
        a = "single mean 463 counts and common dispersion 0.0725 for every gene; pi0 0.95; 1.5-fold; FDR 0.05"
    rows.append(dict(section="b_genome_wide_DE", endpoint="average power 0.80 over DE genes at FDR 0.05, 1.5-fold",
                     method=f"ssizeRNA_{'vary' if k.startswith('vary') else 'single'} (voom-type, FDR-controlled sample size)",
                     assumptions=a, n_unit="animals per group", n_estimate=n, n_conservative="",
                     reference=REF_SS + f"; ssizeRNA power at n=3/group {pw3:.2f}, n=6/group {pw6:.2f}"))
pq = S / "proper_quick.json"
if pq.exists():
    q = json.loads(pq.read_text()); pw["PROPER_summary"] = q["summary"]; pw["PROPER_nsims"] = 5
if "PROPER_summary" in pw:
    REF_P = "Wu, Wang & Wu 2015 Bioinformatics 31:233 doi:10.1093/bioinformatics/btu640 (PROPER 1.42.0; edgeR exact test)"
    summ = pw["PROPER_summary"]
    # summaryPower rows: one per sample size, columns include 'Marginal power'
    def getp(row):
        for kk in row:
            if "arginal power" in kk or "Marginal.power" in kk:
                return row[kk]
    ns, ps = [], []
    for row in summ:
        n1 = row.get("SS1")
        ns.append(n1); ps.append(getp(row))
    first80 = next((n for n, p in zip(ns, ps) if p is not None and p >= 0.8), None)
    curve = "; ".join(f"n={n}: {p:.2f}" for n, p in zip(ns, ps) if p is not None)
    fdrs = "; ".join(f"n={row['SS1']}: {row['Actual FDR']:.3f}" for row in summ)
    rows.append(dict(section="b_genome_wide_DE", endpoint="marginal power over all simulated 1.5-fold DE genes with mean count >10, FDR 0.05",
                     method="PROPER simulation (5 simulations x ~796 DE genes; negative binomial; edgeR exact test; BH FDR)",
                     assumptions=f"p.DE 0.05; lfc +/-log(1.5); baseline means and dispersions = our gene-wise pairs; {DISPTXT}; marginal power by n per group {curve}; actual FDR at nominal 0.05 by n {fdrs}; by-expression power in proper_quick.json",
                     n_unit="animals per group", n_estimate=first80 if first80 else f">{max(ns)}", n_conservative="", reference=REF_P))

# ---- part c: CIBERSORTx
NEW = "Newman et al. 2019 Nat Biotechnol 37:773 doi:10.1038/s41587-019-0114-2 (PMC6610714)"
STE = "Steen et al. 2020 Methods Mol Biol 2117:135 doi:10.1007/978-1-0716-0301-7_7 (PMC7695353), Note 12 and 3.2.1"
sigs = [("Neftel four-state (neftel4_confident; the run in question)", 4), ("Varn 2022 GLASS tumour 3-state", 3),
        ("IvyGAP anatomic 4", 4), ("BayesPrism GBM 6 cell types", 6), ("Neftel 4 + 3 non-malignant", 7),
        ("Mehani 2022 glioma 8", 8), ("Varn 2022 GLASS 12-class", 12), ("LM22", 22)]
for name, c in sigs:
    rows.append(dict(section="c_CIBERSORTx", endpoint=f"group-mode GEP purification, signature {name}, c = {c} cell types",
                     method="Eq. 2 NNLS must be overdetermined (k > c); 'largest gains ... at least 4-5 fold more mixture samples than cell types'",
                     assumptions="k = mixture samples in ONE group-mode run. To contrast primary vs recurrent, the protocol runs group mode on each class separately, so k applies per group.",
                     n_unit="mixture samples per run (per group if run per class)", n_estimate=c + 1, n_conservative=f"{4 * c}-{5 * c}",
                     reference=f"{NEW}; {STE}"))
    rows.append(dict(section="c_CIBERSORTx", endpoint=f"high-resolution (sample-level) purification, signature {name}, c = {c}",
                     method="sliding window w bounded c < w <= k/2; authors set w to 4-5 fold c; docker hires default window = 4 x cell types",
                     assumptions="n_estimate = smallest k allowed (w = c+1 -> k >= 2(c+1)); n_conservative = k at the authors' w = 4c-5c (k >= 8c-10c). Across both groups if run once on the whole set, per group if run per class.",
                     n_unit="mixture samples per run", n_estimate=2 * (c + 1), n_conservative=f"{8 * c}-{10 * c}",
                     reference=f"{NEW} Methods 'Pseudocode for high-resolution expression purification'; window default from the CIBERSORTx hires docker usage text as mirrored at https://open.bioqueue.org/home/knowledge/showKnowledge/sig/cibersortx-hires (third-party mirror, not the official page)"))
rows.append(dict(section="c_CIBERSORTx", endpoint="B-mode or S-mode batch correction (fractions or GEP)", method="ComBat-type correction of mixtures (B-mode) or of the signature (S-mode)",
                 assumptions="B-mode: bulk-sorted or plate-based non-UMI (Smart-seq2) signatures; S-mode: droplet/UMI (10x) signatures; S-mode needs the single-cell reference matrix used to build the signature",
                 n_unit="mixture samples", n_estimate=3, n_conservative=10, reference=f"{NEW} Methods 'B-mode'; {STE} 3.2.1 and Table 2"))
rows.append(dict(section="c_CIBERSORTx", endpoint="general sample-size caveat", method="Discussion",
                 assumptions="'Although further developments are needed to better accommodate smaller sample sizes (e.g., <15) ...'",
                 n_unit="mixture samples", n_estimate=15, n_conservative="", reference=NEW))
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in cols})
print("wrote", OUT, len(rows), "rows")
