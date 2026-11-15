#!/usr/bin/env python3
"""Collect the rnasplice relabelling null: per-task metrics, reproduction of the pipeline by task 1, the rank of the
true split among the 10 balanced splits, and the leave-one-out of the true split.

Writes null/results/{null_metrics.tsv, null_ranks.tsv, null_loo.tsv, null_summary.json, NULL_REPORT.md}.
Exit 0 only if all 16 tasks finished every method and task 1 reproduced every pipeline table it re-derives.
Standard library only.

Pre-declared (README.md, before any result existed) primary count per method:
  DEU    genes with DEXSeq perGeneQValue < 0.05
  DTU    genes with DEXSeq (DTU) perGeneQValue < 0.05
  rMATS  events with FDR < 0.05 and |IncLevelDifference| >= 0.1 in the junction-count (JC) files, summed over
         SE, MXE, A3SS, A5SS and RI. JC, not JCEC: in total RNA the exon/intron-body reads of JCEC include pre-mRNA.
  SUPPA  local events (ioe) with gene-corrected p < 0.05
Secondary counts are reported beside them and are never promoted after the fact.
"""
import csv
import datetime
import json
import math
import os
import sys

M = "/rs/rs_grp_oschome/go2432/u251_multiomics"
D = os.environ.get("RS_D", os.path.join(M, "rnasplice"))   # RS_D only for the local dry run
N = os.path.join(D, "null")
OUT = os.path.join(N, "results")
TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]
TOL = 1e-6
PRIMARY = ["deu_genes_q05", "dtu_genes_q05", "rmats_JC_q05_d01", "suppa_local_p05"]
SECONDARY = ["deu_bins_q05", "dtu_tx_stageR05", "rmats_JCEC_q05_d01", "suppa_local_p05_d01", "suppa_iso_p05"]
METHOD_OF = {"deu": "deu", "dtu": "dtu", "rmats": "rmats", "suppa": "suppa"}
csv.field_size_limit(1 << 30)


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def fnum(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


def env(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1)
            out[k] = v.strip("'")
    return out


def table(path, delim):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def gene_q(path, delim):
    """groupID -> padj (None for NA)"""
    return {r["groupID"]: fnum(r.get("padj")) for r in table(path, delim)}


def exon_q(path):
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out[(r.get("groupID"), r.get("featureID"))] = fnum(r.get("padj"))
    return out


def stager(path):
    return {r["txID"]: (fnum(r.get("gene")), fnum(r.get("transcript"))) for r in table(path, "\t")}


def dpsi(path):
    """SUPPA .dpsi: header has two fields; rows: event, dPSI, p"""
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                out[f[0]] = (fnum(f[1]), fnum(f[2]))
    return out


def mats(path):
    """rMATS MATS file -> {coordinate key: (PValue, FDR, IncLevelDifference)}"""
    out = {}
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        id2 = [i for i, h in enumerate(hdr) if h == "ID"]
        kend = id2[1] if len(id2) > 1 else hdr.index("IJC_SAMPLE_1")
        ip, iq, idf = hdr.index("PValue"), hdr.index("FDR"), hdr.index("IncLevelDifference")
        for f in rd:
            out[tuple(f[1:kend])] = (fnum(f[ip]), fnum(f[iq]), fnum(f[idf]))
    return out


def count(d, pred):
    return sum(1 for v in d.values() if pred(v))


def m_deu(w):
    g = gene_q(os.path.join(w, "perGeneQValue.recurrent_vs_primary.csv"), ",")
    e = exon_q(os.path.join(w, "DEXSeqResults.recurrent_vs_primary.csv"))
    m = {"deu_genes_tested": count(g, lambda v: v is not None),
         "deu_genes_q05": count(g, lambda v: v is not None and v < 0.05),
         "deu_bins_q05": count(e, lambda v: v is not None and v < 0.05)}
    return m, sorted(k for k, v in g.items() if v is not None and v < 0.05), "done"


def m_dtu(w):
    g = gene_q(os.path.join(w, "perGeneQValue.recurrent-primary.tsv"), "\t")
    s = stager(os.path.join(w, "getAdjustedPValues.recurrent-primary.tsv"))
    m = {"dtu_genes_tested": count(g, lambda v: v is not None),
         "dtu_genes_q05": count(g, lambda v: v is not None and v < 0.05),
         "dtu_tx_stageR05": count(s, lambda v: v[1] is not None and v[1] < 0.05)}
    return m, sorted(k for k, v in g.items() if v is not None and v < 0.05), "done"


def m_suppa(w):
    loc = dpsi(os.path.join(w, "recurrent-primary_local_diffsplice.dpsi"))
    iso = dpsi(os.path.join(w, "recurrent-primary_transcript_diffsplice.dpsi"))
    m = {"suppa_local_tested": count(loc, lambda v: v[1] is not None),
         "suppa_local_p05": count(loc, lambda v: v[1] is not None and v[1] < 0.05),
         "suppa_local_p05_d01": count(loc, lambda v: v[1] is not None and v[1] < 0.05 and v[0] is not None and abs(v[0]) >= 0.1),
         "suppa_iso_p05": count(iso, lambda v: v[1] is not None and v[1] < 0.05)}
    return m, sorted(k for k, v in loc.items() if v[1] is not None and v[1] < 0.05), "done"


def m_rmats(w):
    route = open(os.path.join(w, "ROUTE")).read().strip() if os.path.isfile(os.path.join(w, "ROUTE")) else "?"
    m, s_jc = {}, []
    for ct in ("JC", "JCEC"):
        tot, tested = 0, 0
        for t in TYPES:
            d = mats(os.path.join(w, "rmats_out", "%s.MATS.%s.txt" % (t, ct)))
            k = [key for key, v in d.items() if v[1] is not None and v[1] < 0.05 and v[2] is not None and abs(v[2]) >= 0.1]
            m["rmats_%s_%s_q05_d01" % (ct, t)] = len(k)
            tot += len(k)
            tested += count(d, lambda v: v[1] is not None)
            if ct == "JC":
                s_jc += [(t,) + key for key in k]
        m["rmats_%s_q05_d01" % ct] = tot
        m["rmats_%s_tested" % ct] = tested
    return m, sorted("|".join(x) for x in s_jc), "done (%s)" % route


def metrics(run):
    """-> (metrics, significant sets, status) per method. A method that cannot be read is 'missing' or 'error: ...'
    and its counts are absent (never 0), so an incomplete null can never produce a rank."""
    m, sig, st = {}, {}, {}
    for method, fn in (("deu", m_deu), ("dtu", m_dtu), ("suppa", m_suppa), ("rmats", m_rmats)):
        w = os.path.join(run, method)
        if not os.path.isfile(os.path.join(w, "DONE")):
            st[method] = "missing"
            continue
        try:
            mm, ss, tt = fn(w)
        except Exception as e:  # keep going: the report must show which task broke
            st[method] = "error: %s" % e
            continue
        m.update(mm)
        sig[method] = ss
        st[method] = tt
    return m, sig, st


def maxdiff(a, b):
    """a, b: dict key -> number or tuple of numbers. -> (same keys, max abs diff, n NA mismatches)"""
    if set(a) != set(b):
        return False, float("inf"), len(set(a) ^ set(b))
    worst, nam = 0.0, 0
    for k, va in a.items():
        vb = b[k]
        ta = va if isinstance(va, tuple) else (va,)
        tb = vb if isinstance(vb, tuple) else (vb,)
        for x, y in zip(ta, tb):
            if x is None or y is None:
                nam += (x is None) != (y is None)
            else:
                worst = max(worst, abs(x - y))
    return True, worst, nam


def reproduce(E, run):
    """task 1 against the pipeline's own tables"""
    out = {}

    def rec(method, name, a_path, b_path, loader):
        if not (os.path.isfile(a_path) and os.path.isfile(b_path)):
            out.setdefault(method, []).append((name, False, "missing %s" % (a_path if not os.path.isfile(a_path) else b_path)))
            return
        same, worst, nam = maxdiff(loader(a_path), loader(b_path))
        ok = same and worst <= TOL and nam == 0
        out.setdefault(method, []).append((name, ok, "same keys %s, max |diff| %.3g, NA mismatches %d" % (same, worst, nam)))

    ss = os.path.join(D, "results", "star_salmon")
    rec("deu", "perGeneQValue", E["DEU_PIPE_GENE"], os.path.join(run, "deu", "perGeneQValue.recurrent_vs_primary.csv"),
        lambda p: gene_q(p, ","))
    rec("deu", "DEXSeqResults padj", os.path.join(ss, "dexseq_exon", "results", "DEXSeqResults.recurrent_vs_primary.csv"),
        os.path.join(run, "deu", "DEXSeqResults.recurrent_vs_primary.csv"), exon_q)
    rec("dtu", "perGeneQValue", E["DTU_PIPE_GENE"], os.path.join(run, "dtu", "perGeneQValue.recurrent-primary.tsv"),
        lambda p: gene_q(p, "\t"))
    rec("dtu", "stageR", E["STAGER_PIPE"], os.path.join(run, "dtu", "getAdjustedPValues.recurrent-primary.tsv"), stager)
    rec("suppa", "local dpsi", E["SUPPA_PIPE_L"], os.path.join(run, "suppa", "recurrent-primary_local_diffsplice.dpsi"), dpsi)
    rec("suppa", "isoform dpsi", E["SUPPA_PIPE_I"], os.path.join(run, "suppa", "recurrent-primary_transcript_diffsplice.dpsi"), dpsi)
    for t in TYPES:
        for ct in ("JC", "JCEC"):
            f = "%s.MATS.%s.txt" % (t, ct)
            rec("rmats", f, os.path.join(E["RMATS_PIPE_POST"], f), os.path.join(run, "rmats", "rmats_out", f), mats)
    return {k: {"ok": all(x[1] for x in v), "detail": v} for k, v in out.items()}


def rank_of(vals):
    """vals: 10 values, index 0 = true split. rank 1 = most discoveries; ties give a range."""
    if any(v is None for v in vals):
        return None
    v1 = vals[0]
    gt = sum(1 for v in vals[1:] if v > v1)
    eq = sum(1 for v in vals[1:] if v == v1)
    return (1 + gt, 1 + gt + eq)


def main():
    os.makedirs(OUT, exist_ok=True)
    E = env(os.path.join(N, "inputs.env"))
    splits = table(os.path.join(N, "splits.tsv"), "\t")
    res, sigs, status = {}, {}, {}
    for r in splits:
        run = os.path.join(N, "runs", r["label"])
        res[r["label"]], sigs[r["label"]], status[r["label"]] = metrics(run)
    labels = [r["label"] for r in splits]
    split_labels = [r["label"] for r in splits if r["kind"] == "split"]
    loo_labels = [r["label"] for r in splits if r["kind"] == "loo"]
    repro = reproduce(E, os.path.join(N, "runs", "split01"))
    allm = PRIMARY + SECONDARY
    tested = ["deu_genes_tested", "dtu_genes_tested", "rmats_JC_tested", "rmats_JCEC_tested", "suppa_local_tested"]

    with open(os.path.join(OUT, "null_metrics.tsv"), "w") as fh:
        cols = allm + tested + ["rmats_JC_%s_q05_d01" % t for t in TYPES]
        fh.write("label\tkind\trecurrent_group\tprimary_group\t" + "\t".join(cols) + "\tstatus\n")
        for r in splits:
            m = res[r["label"]]
            fh.write("\t".join([r["label"], r["kind"], r["treatment_recurrent"], r["control_primary"]] +
                               [str(m.get(c, "NA")) for c in cols] + ["; ".join("%s %s" % kv for kv in status[r["label"]].items())]) + "\n")

    ranks = {}
    with open(os.path.join(OUT, "null_ranks.tsv"), "w") as fh:
        fh.write("metric\tprimary\tmethod\ttrue_split\tother_nine\trank_range\tp_range\treproduced\n")
        for mname in allm:
            method = mname.split("_")[0]
            vals = [res[l].get(mname) for l in split_labels]
            rk = rank_of(vals)
            rep = repro.get(method, {}).get("ok", False)
            ranks[mname] = {"values": vals, "rank": rk, "p": (rk[0] / 10.0, rk[1] / 10.0) if rk else None,
                            "reproduced": rep, "primary": mname in PRIMARY}
            fh.write("\t".join([mname, "yes" if mname in PRIMARY else "no", method, str(vals[0]),
                                ",".join(str(v) for v in vals[1:]),
                                ("%d-%d" % rk) if rk else "incomplete",
                                ("%.1f-%.1f" % (rk[0] / 10.0, rk[1] / 10.0)) if rk else "NA",
                                "yes" if rep else "NO"]) + "\n")

    loo = {}
    with open(os.path.join(OUT, "null_loo.tsv"), "w") as fh:
        fh.write("dropped\t" + "\t".join("%s\t%s_kept_of_true" % (m, m) for m in PRIMARY) + "\n")
        for l in loo_labels:
            row = [l.replace("loo_", "")]
            loo[l] = {}
            for mname in PRIMARY:
                method = mname.split("_")[0]
                v = res[l].get(mname)
                t = set(sigs["split01"].get(method, []))
                k = len(t & set(sigs[l].get(method, []))) if t else 0
                loo[l][mname] = {"value": v, "true_kept": k, "true_n": len(t)}
                row += [str(v if v is not None else "NA"), "%d/%d" % (k, len(t))]
            fh.write("\t".join(row) + "\n")

    complete = all(s.startswith("done") for l in labels for s in status[l].values())
    reproduced = all(repro.get(k, {}).get("ok", False) for k in ("deu", "dtu", "suppa", "rmats"))
    summary = {"written": now(), "complete": complete, "reproduced": reproduced, "status": status, "metrics": res,
               "ranks": ranks, "loo": loo, "reproduction": repro, "primary": PRIMARY, "secondary": SECONDARY}
    with open(os.path.join(OUT, "null_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=list)

    md = ["# rnasplice relabelling null", "",
          "Written %s by null_collect.py. Recurrent (IL66B, NL70B, NL71B) against primary (IL67B, IL68B, IL69B), graft-only "
          "reads, nf-core/rnasplice 1.0.4. Each method's count for the true labels is ranked among the 10 balanced 3 v 3 "
          "relabellings of the same six tumours (IL67B fixed in one group, so no split appears with its mirror). Rank 1 of "
          "10 is p = 0.10, the smallest this design can give; it is not significance at 0.05. Ties give a rank range, and "
          "the upper end is the one to quote." % now(), "",
          "All 16 tasks complete: **%s**. Task 1 reproduced the pipeline's own tables: **%s**." % (complete, reproduced), ""]
    md += ["## Primary counts (declared before the run)", "",
           "| count | true split | other nine splits | rank of true | p | null reproduces pipeline |", "|---|--:|---|--:|--:|---|"]
    for mname in allm:
        if mname == SECONDARY[0]:
            md += ["", "## Secondary counts", "",
                   "| count | true split | other nine splits | rank of true | p | null reproduces pipeline |", "|---|--:|---|--:|--:|---|"]
        r = ranks[mname]
        rk = r["rank"]
        md.append("| %s | %s | %s | %s | %s | %s |" % (
            mname, r["values"][0], ", ".join(str(v) for v in r["values"][1:]),
            ("%d-%d" % rk if rk[0] != rk[1] else "%d" % rk[0]) if rk else "incomplete",
            ("%.1f-%.1f" % r["p"] if rk[0] != rk[1] else "%.1f" % r["p"][0]) if rk else "NA",
            "yes" if r["reproduced"] else "**NO: not calibrated**"))
    md += ["", "## Leave-one-out of the true split (5 tumours, true labels)", "",
           "Each row drops one tumour and refits. `kept` = how many of the true split's significant genes/events are "
           "still significant. The DTU rows reuse the six-sample DRIMSeq filter (same features), so only the fit changes.", "",
           "| dropped | " + " | ".join("%s (kept of true)" % m for m in PRIMARY) + " |", "|---|" + "---|" * len(PRIMARY)]
    for l in loo_labels:
        md.append("| %s | " % l.replace("loo_", "") + " | ".join(
            "%s (%d/%d)" % (loo[l][m]["value"] if loo[l][m]["value"] is not None else "NA", loo[l][m]["true_kept"], loo[l][m]["true_n"])
            for m in PRIMARY) + " |")
    md += ["", "## Task 1 against the pipeline (the null's own check)", ""]
    for method in ("deu", "dtu", "suppa", "rmats"):
        rp = repro.get(method, {"ok": False, "detail": []})
        md.append("- %s: %s" % (method, "reproduced" if rp["ok"] else "**NOT reproduced**"))
        for name, ok, det in rp["detail"]:
            md.append("  - %s: %s (%s)" % (name, "ok" if ok else "FAIL", det))
    md += ["", "## Task status", ""] + ["- %s: %s" % (l, "; ".join("%s %s" % kv for kv in status[l].items())) for l in labels]
    md += ["", "## What this can and cannot show", "",
           "- A true-split count that ranks 1 of 10 is the strongest statement available: p = 0.10 for that count. "
           "Anything ranked 2 or lower is not distinguishable from an arbitrary regrouping of the same six tumours.",
           "- A rank is about the number of discoveries, not about any one gene or event. A named event is reportable only "
           "if its own statistic also ranks 1 of 10 across the relabellings, and even then only at p = 0.10.",
           "- Duplicates are not removed in rnasplice (65-79 % in these tumours), so rMATS and DEXSeq treat PCR copies as "
           "independent fragments and their nominal FDR is optimistic. Every split carries the same duplicates, so the "
           "rank is the calibration, not the FDR.",
           "- The arms also differ in tumour content and in the IL/NL naming batch; a splicing difference that follows "
           "either would rank the same as one that follows LITT. This null cannot separate them.",
           "- Leave-one-out shows whether a count leans on one tumour (the human 'translation' fall leaned on IL68B)."]
    with open(os.path.join(OUT, "NULL_REPORT.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("complete=%s reproduced=%s" % (complete, reproduced))
    for mname in PRIMARY:
        r = ranks[mname]
        print("%s true=%s rank=%s" % (mname, r["values"][0], r["rank"]))
    return 0 if complete and reproduced else 2


if __name__ == "__main__":
    sys.exit(main())
