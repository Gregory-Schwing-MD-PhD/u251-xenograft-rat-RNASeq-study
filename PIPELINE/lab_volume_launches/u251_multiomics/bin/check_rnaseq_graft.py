#!/usr/bin/env python3
"""Checks for rnaseq_graft (nf-core/rnaseq 3.22.2 on graft-only U251 reads).

Writes $D/CHECKS.md, $D/checks.json and gene-level tables under $D/tables/. Standard library only (runs on any
Python >= 3.6). The numbered checks follow SELECTION.md S3; check 0 (versions) is added.

HARD checks (all must pass for "all_hard_pass", which S4 pass B and invasion_niche gate on, and which alone allows
the work directory to be deleted): 0, 1, 2, 3, 4, 9, 10. REPORTS (numbers and verdicts, never pass/fail): 5, 6, 7, 8.
Check 5 carries one decision fixed before the run: clip_material.

Usage: check_rnaseq_graft.py [--d DIR] [--pub DIR] [--xeng DIR] [--genesets DIR]
"""
import argparse
import csv
import datetime
import glob
import gzip
import json
import math
import os
import re
import shutil
import sys
import zipfile

M = "/rs/rs_grp_oschome/go2432/u251_multiomics"
U = "/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study"

LIBS = ["C2B", "IL64B", "IL66B", "IL67B", "IL68B", "IL69B", "N168B", "N269B", "NL70B", "NL71B"]
# xengsort classify counts (read pairs): host, graft, both. Source: U/ANALYSIS/xengsort_out/<lib>.txt (SCOUT.md s2);
# re-parsed from the logs below and required to agree.
XENG = {
    "C2B": (450649, 84519445, 11022904),
    "IL64B": (90003153, 314037, 4546256),
    "IL66B": (51160559, 25459283, 8298076),
    "IL67B": (36298116, 31886295, 5512307),
    "IL68B": (32216898, 32225887, 5441545),
    "IL69B": (19642823, 49327614, 6013633),
    "N168B": (74411463, 435428, 4086055),
    "N269B": (67885104, 3823923, 4423374),
    "NL70B": (35387952, 31293679, 5345807),
    "NL71B": (52351099, 29210358, 5811628),
}
PRIMARY = ["IL67B", "IL68B", "IL69B"]
RECURRENT = ["IL66B", "NL70B", "NL71B"]
TUMOURS = PRIMARY + RECURRENT
CONTROLS = ["IL64B", "N168B"]
SAMPLES = [l + "_graft" for l in LIBS] + ["IL69B_fullclip"]
BIG = [l for l in LIBS if XENG[l][1] >= 3800000]          # the 8 libraries with >= 3.8 M graft pairs
PRIOR_LEAK = {"primary": (1.9, 5.1), "recurrent": (5.1, 8.0)}  # % of the human file, EVAL_transcripts estimate


def expected_pairs(sample):
    lib, kind = sample.rsplit("_", 1)
    h, g, b = XENG[lib]
    return g + b if kind == "fullclip" else g


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# ------------------------------------------------------------------------------------------------ small utilities
def read_text(path):
    with open(path, "r", errors="replace") as fh:
        return fh.read()


def rank(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx = sum(x) / n
    my = sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def spearman(x, y):
    return pearson(rank(x), rank(y))


def median(v):
    v = sorted(v)
    n = len(v)
    if n == 0:
        return float("nan")
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def fmt(x, nd=4):
    if x is None:
        return "NA"
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return ("%." + str(nd) + "f") % x
    return str(x)


def read_counts(path):
    """salmon.merged.gene_counts.tsv -> (samples, {gene_id: [counts]}, {gene_id: gene_name})"""
    counts, names = {}, {}
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        samples = hdr[2:]
        for row in rd:
            counts[row[0]] = [float(v) for v in row[2:]]
            names[row[0]] = row[1]
    return samples, counts, names


def column(samples, counts, s):
    i = samples.index(s)
    return {g: v[i] for g, v in counts.items()}


def cpm(col):
    tot = sum(col.values())
    return {g: (v * 1e6 / tot if tot > 0 else 0.0) for g, v in col.items()}, tot


def parse_gtf_genes(path):
    """gene_id -> (seqname, gene_type, gene_name) from gene lines; also returns the header lines"""
    genes, header = {}, []
    op = gzip.open if path.endswith(".gz") else open
    rx = re.compile(r'(\S+) "([^"]*)"')
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if len(header) < 20:
                    header.append(line.rstrip("\n"))
                continue
            if "\tgene\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(rx.findall(f[8]))
            genes[a.get("gene_id", "")] = (f[0], a.get("gene_type", a.get("gene_biotype", "")), a.get("gene_name", ""))
    return genes, header


def star_log(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in read_text(path).splitlines():
        if "|" in line:
            k, v = line.split("|", 1)
            out[k.strip()] = v.strip()
    return out


def pct(s):
    try:
        return float(str(s).rstrip("%"))
    except ValueError:
        return float("nan")


def fastqc_zip(path):
    """-> dict with total, minlen, maxlen, filename (from fastqc_data.txt inside the zip)"""
    out = {}
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith("fastqc_data.txt")][0]
        for line in z.read(name).decode("utf-8", "replace").splitlines():
            if line.startswith("Filename\t"):
                out["filename"] = line.split("\t", 1)[1]
            elif line.startswith("Total Sequences\t"):
                out["total"] = int(line.split("\t")[1])
            elif line.startswith("Sequence length\t"):
                v = line.split("\t")[1]
                parts = [int(p) for p in v.split("-")]
                out["minlen"], out["maxlen"] = min(parts), max(parts)
    return out


def trim_report(path):
    out = {"cmd": None, "total": None, "pairs_removed": None, "five_prime_clip": False}
    if not os.path.isfile(path):
        return out
    for line in read_text(path).splitlines():
        if line.startswith("Command line parameters:"):
            out["cmd"] = line.split(":", 1)[1].strip()
        elif line.startswith("Total reads processed:"):
            out["total"] = int(line.split(":", 1)[1].strip().replace(",", ""))
        elif line.startswith("Number of sequence pairs removed because at least one read was shorter"):
            out["pairs_removed"] = int(line.split(":", 1)[1].strip().split()[0])
        elif re.search(r"trimmed by \d+ bp from their 5' end", line):
            out["five_prime_clip"] = True
    return out


def trace_hashes(res):
    """process short name -> list of (tag, hash) of the latest COMPLETED/CACHED attempt, over all traces"""
    out = {}
    for tr in sorted(glob.glob(os.path.join(res, "pipeline_info", "execution_trace_*.txt"))):
        with open(tr) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                if row.get("status") not in ("COMPLETED", "CACHED"):
                    continue
                name = row.get("name", "")
                m = re.match(r"^(\S+?)(?: \((.*)\))?$", name)
                proc = m.group(1).split(":")[-1] if m else name
                tag = (m.group(2) if m else "") or ""
                out.setdefault(proc, {})[tag] = row.get("hash", "")
    return out


def work_dir(work, h):
    if not h or "/" not in h:
        return None
    a, b = h.split("/", 1)
    hits = glob.glob(os.path.join(work, a, b + "*"))
    return hits[0] if len(hits) == 1 else None


def tree_listing(d):
    out = {}
    for root, _, files in os.walk(d, followlinks=True):
        for f in files:
            p = os.path.join(root, f)
            try:
                out[os.path.relpath(p, d)] = os.path.getsize(p)
            except OSError:
                out[os.path.relpath(p, d)] = -1
    return out


def read_list(path):
    if not os.path.isfile(path):
        return []
    return [l.strip() for l in read_text(path).splitlines() if l.strip() and not l.startswith("#")]


def versions_yml(path):
    out, cur = {}, None
    if not os.path.isfile(path):
        return out
    for line in read_text(path).splitlines():
        if not line.strip():
            continue
        if not line.startswith(" "):
            cur = line.strip().rstrip(":").strip("'\"")
            out[cur] = {}
        elif cur is not None and ":" in line:
            k, v = line.strip().split(":", 1)
            out[cur][k.strip().strip("'\"")] = v.strip().strip("'\"")
    return out


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", default=os.path.join(M, "rnaseq_graft"))
    ap.add_argument("--pub", default=os.path.join(U, "ANALYSIS", "results_human_final"))
    ap.add_argument("--xeng", default=os.path.join(U, "ANALYSIS", "xengsort_out"))
    ap.add_argument("--genesets", default=None)
    ap.add_argument("--repair-index", type=int, default=1)
    args = ap.parse_args()
    D = args.d
    RES = os.path.join(D, "results")
    WORK = os.path.join(D, "work")
    PUB = args.pub
    GS = args.genesets or os.path.join(D, "genesets")
    TAB = os.path.join(D, "tables")
    os.makedirs(TAB, exist_ok=True)
    C = {"started": now(), "results": RES, "published": PUB, "checks": {}}
    md = ["# rnaseq_graft checks", "", "Written %s by check_rnaseq_graft.py. Results: `%s`. Published human run: `%s`." % (now(), RES, PUB), ""]

    def record(key, title, passed, lines, data=None, hard=True):
        C["checks"][key] = {"title": title, "hard": hard, "pass": passed, "data": data or {}}
        status = ("PASS" if passed else "FAIL") if hard else "REPORT"
        md.append("## %s. %s: %s" % (key, title, status))
        md.append("")
        md.extend(lines)
        md.append("")

    nfx = os.path.join(D, "NEXTFLOW_EXIT")
    nfx_v = read_text(nfx).strip() if os.path.isfile(nfx) else "missing"
    C["nextflow_exit"] = nfx_v
    md.append("Nextflow exit recorded by the head job: `%s`." % nfx_v)
    md.append("")
    if not os.path.isdir(RES):
        record("run", "results directory", False, ["No results directory: the pipeline did not run (gate or failure)."])
        return finish(C, md, D)

    # xengsort logs re-parsed (outside the pipeline): must equal the table this script carries
    xe_lines, xe_ok = [], True
    for lib in LIBS:
        f = os.path.join(args.xeng, lib + ".txt")
        got = None
        if os.path.isfile(f):
            for line in read_text(f).splitlines():
                p = line.split()
                if len(p) == 6 and p[0] == lib:
                    got = (int(p[1]), int(p[2]), int(p[4]))
                    break
        if got is not None and got != XENG[lib]:
            xe_ok = False
            xe_lines.append("- %s: xengsort log says %s, table says %s" % (lib, got, XENG[lib]))
        elif got is None:
            xe_lines.append("- %s: xengsort log not parsed (%s); the carried table is used" % (lib, f))
    C["xengsort_table_matches_logs"] = xe_ok
    if not xe_ok:
        md.extend(["**xengsort counts disagree with the carried table:**"] + xe_lines + [""])

    # ---------------------------------------------------------------------------------- 0. versions
    vn = versions_yml(os.path.join(RES, "pipeline_info", "nf_core_rnaseq_software_mqc_versions.yml"))
    vp = versions_yml(os.path.join(PUB, "pipeline_info", "nf_core_rnaseq_software_mqc_versions.yml"))
    want = [("Workflow", "nf-core/rnaseq"), ("Workflow", "Nextflow"), ("STAR_ALIGN", "star"),
            ("STAR_GENOMEGENERATE", "star"), ("SALMON_QUANT", "salmon"), ("SALMON_INDEX", "salmon"),
            ("TRIMGALORE", "trimgalore"), ("TRIMGALORE", "cutadapt"), ("FASTQC", "fastqc"),
            ("RSEQC_INFEREXPERIMENT", "rseqc")]
    ok0, l0, d0 = True, ["| process | tool | this run | published |", "|---|---|---|---|"], {}
    for proc, tool in want:
        a = vn.get(proc, {}).get(tool)
        b = vp.get(proc, {}).get(tool)
        d0["%s/%s" % (proc, tool)] = [a, b]
        if a is None or a != b:
            ok0 = False
        l0.append("| %s | %s | %s | %s |" % (proc, tool, a, b))
    pj = sorted(glob.glob(os.path.join(RES, "pipeline_info", "params_*.json")))
    params = {}
    if pj:
        with open(pj[-1]) as fh:
            params = json.load(fh)
    etg = params.get("extra_trimgalore_args", "missing")
    l0.append("")
    l0.append("params: extra_trimgalore_args=%s, save_reference=%s, featurecounts_group_type=%s, min_mapped_reads=%s, gencode=%s"
              % (etg, params.get("save_reference"), params.get("featurecounts_group_type"),
                 params.get("min_mapped_reads"), params.get("gencode")))
    if etg not in (None, "null", ""):
        ok0 = False
    record("0", "Same release and tool versions as the published run", ok0, l0, d0)

    # ---------------------------------------------------------------------------------- 1. clip exactly once
    l1, ok1, d1 = ["| sample | R1 len | R2 len | R1 FastQC n | R2 FastQC n | TrimGalore R1 cmd / R2 cmd 5' clip |", "|---|---|---|---|---|---|"], True, {}
    fq_tot = {}
    for s in SAMPLES:
        r = {}
        for k in (1, 2):
            z = os.path.join(RES, "fastqc", "raw", "%s_raw_%d_fastqc.zip" % (s, k))
            try:
                r[k] = fastqc_zip(z)
            except Exception as e:  # noqa
                r[k] = {"error": str(e)}
        t1 = trim_report(os.path.join(RES, "trimgalore", "%s_trimmed_1.fastq.gz_trimming_report.txt" % s))
        t2 = trim_report(os.path.join(RES, "trimgalore", "%s_trimmed_2.fastq.gz_trimming_report.txt" % s))
        clip_seen = False
        for t in (t1, t2):
            c = t["cmd"] or ""
            if re.search(r"(^|\s)(-U|-u|--cut)\s", c) or "clip_r" in c.lower() or t["five_prime_clip"]:
                clip_seen = True
        good = (r[1].get("maxlen") == 101 and r[2].get("maxlen") == 98 and t1["cmd"] is not None
                and t2["cmd"] is not None and not clip_seen)
        ok1 = ok1 and good
        fq_tot[s] = (r[1].get("total"), r[2].get("total"))
        d1[s] = {"r1": r[1], "r2": r[2], "trim_cmd_r1": t1["cmd"], "trim_cmd_r2": t2["cmd"], "clip_in_trimgalore": clip_seen}
        l1.append("| %s | %s-%s | %s-%s | %s | %s | %s |" % (
            s, r[1].get("minlen"), r[1].get("maxlen"), r[2].get("minlen"), r[2].get("maxlen"),
            r[1].get("total"), r[2].get("total"), "CLIP SEEN" if clip_seen else "none"))
    l1.insert(0, "Rule: raw R1 max length 101 and raw R2 max 98 in all 11 samples (prep clipped R2 by 3 nt; the reads were "
                 "35-101 nt, so R2 should read 32-98), and no TrimGalore report shows a 5' clip. Command line seen: `%s`."
              % (d1.get("C2B_graft", {}).get("trim_cmd_r2")))
    l1.insert(1, "")
    record("1", "Takara R2 clip applied exactly once", ok1, l1, d1)

    # ---------------------------------------------------------------------------------- 2. input counts
    l2, ok2, d2 = ["| sample | expected pairs | FastQC R1/R2 | TrimGalore R1/R2 total | pairs removed (<20 nt) | pairs written | STAR input | ok |",
                   "|---|---|---|---|---|---|---|---|"], True, {}
    for s in SAMPLES:
        e = expected_pairs(s)
        t1 = trim_report(os.path.join(RES, "trimgalore", "%s_trimmed_1.fastq.gz_trimming_report.txt" % s))
        t2 = trim_report(os.path.join(RES, "trimgalore", "%s_trimmed_2.fastq.gz_trimming_report.txt" % s))
        rem = t2["pairs_removed"] if t2["pairs_removed"] is not None else t1["pairs_removed"]
        written = (t1["total"] - rem) if (t1["total"] is not None and rem is not None) else None
        sl = star_log(os.path.join(RES, "star_salmon", "log", "%s.Log.final.out" % s))
        si = int(sl["Number of input reads"]) if "Number of input reads" in sl else None
        good = (t1["total"] == e and t2["total"] == e and fq_tot.get(s, (None, None)) == (e, e)
                and written is not None and si == written)
        ok2 = ok2 and good
        d2[s] = {"expected": e, "fastqc": fq_tot.get(s), "trim_total": [t1["total"], t2["total"]],
                 "pairs_removed": rem, "pairs_written": written, "star_input": si, "ok": good}
        l2.append("| %s | %d | %s/%s | %s/%s | %s | %s | %s | %s |" % (
            s, e, fq_tot.get(s, (None, None))[0], fq_tot.get(s, (None, None))[1], t1["total"], t2["total"], rem,
            written, si, "ok" if good else "MISMATCH"))
    l2.insert(0, "Expected = N_graft from the xengsort logs (N_graft + N_both for IL69B_fullclip), a count made before any of "
                 "this ran. STAR input must equal the pairs TrimGalore kept.")
    l2.insert(1, "")
    record("2", "Input counts carried from xengsort to STAR", ok2, l2, d2)

    # ---------------------------------------------------------------------------------- 3. strand
    gated = ["C2B_graft"] + [t + "_graft" for t in TUMOURS] + ["IL69B_fullclip"]
    l3, ok3, d3 = ["| sample | antisense of assignable | antisense of all | sense | undetermined | salmon format / compatible | gated |",
                   "|---|---|---|---|---|---|---|"], True, {}
    for s in SAMPLES:
        f = os.path.join(RES, "star_salmon", "rseqc", "infer_experiment", "%s.infer_experiment.txt" % s)
        fs = fa = fu = None
        if os.path.isfile(f):
            for line in read_text(f).splitlines():
                if "failed to determine" in line:
                    fu = float(line.split(":")[-1])
                elif '"1++,1--,2+-,2-+"' in line:
                    fs = float(line.split(":")[-1])
                elif '"1+-,1-+,2++,2--"' in line:
                    fa = float(line.split(":")[-1])
        share = fa / (fa + fs) if (fa is not None and fs is not None and fa + fs > 0) else None
        lf = os.path.join(RES, "star_salmon", s, "lib_format_counts.json")
        sal = "NA"
        if os.path.isfile(lf):
            try:
                j = json.load(open(lf))
                sal = "%s / %s" % (j.get("expected_format"), fmt(j.get("compatible_fragment_ratio")))
            except Exception:  # noqa
                pass
        g = s in gated
        good = share is not None and share >= 0.95
        if g:
            ok3 = ok3 and good
        d3[s] = {"antisense_assignable": share, "antisense_all": fa, "sense_all": fs, "undetermined": fu, "salmon": sal}
        l3.append("| %s | %s | %s | %s | %s | %s | %s |" % (s, fmt(share), fmt(fa), fmt(fs), fmt(fu), sal,
                                                            ("yes" if g else "no (reported only)") + ("" if (good or not g) else " FAIL")))
    st = os.path.join(RES, "multiqc", "star_salmon", "multiqc_report_data", "multiqc_fail_strand_check_table.txt")
    if os.path.isfile(st):
        rows = list(csv.DictReader(open(st), delimiter="\t"))
        bad = [r for r in rows if r.get("Status", "").lower() != "pass" and r.get("Sample", "").split("_raw")[0] in gated]
        l3.append("")
        l3.append("Pipeline strand-check table: %d rows; not 'pass' among gated samples: %s" % (
            len(rows), ", ".join("%s(%s, inferred %s)" % (r.get("Sample"), r.get("Status"), r.get("Inferred strandedness")) for r in bad) or "none"))
        if bad:
            ok3 = False
        d3["strand_table"] = rows
    l3.insert(0, "Rule: antisense share of assignable reads >= 0.95 for C2B and the six tumours (published graft+both run: "
                 "97.4-98.4 %). IL64B, N168B and N269B are reported only; a sample under min_mapped_reads (5 %) has no RSeQC "
                 "output but keeps its salmon counts.")
    l3.insert(1, "")
    record("3", "Reverse-stranded (explicit, not inferred)", ok3, l3, d3)

    # ---------------------------------------------------------------------------------- 4. reference identity
    gtfs = sorted(glob.glob(os.path.join(RES, "genome", "*.gtf")))
    l4, ok4, d4 = [], False, {"gtfs": gtfs}
    genes, header, gtf_used = {}, [], None
    for g in gtfs:
        gg, hh = parse_gtf_genes(g)
        v44 = any("version 44 (Ensembl 110)" in h for h in hh)
        l4.append("- `%s`: %d genes; header v44 line %s" % (os.path.basename(g), len(gg), "present" if v44 else "ABSENT"))
        d4[os.path.basename(g)] = {"genes": len(gg), "v44_header": v44}
        if v44:
            ok4 = True
        if len(gg) > len(genes):
            genes, header, gtf_used = gg, hh, g
    if not ok4:
        # the filtered GTF may drop the header; fall back to the published input GTF's header only if a result GTF
        # has the same gene count as it (identity by content, stated as such)
        src = os.path.join(U, "ANALYSIS", "refs", "human", "GRCh38.primary_assembly.annotation.gtf.gz")
        if os.path.isfile(src):
            gs, hs = parse_gtf_genes(src)
            if any("version 44 (Ensembl 110)" in h for h in hs) and gtfs and any(
                    d4[os.path.basename(g)]["genes"] == len(gs) for g in gtfs):
                l4.append("- no result GTF keeps a header; the input GTF (v44 header present) has the same gene count "
                          "(%d), accepted as the same file by content" % len(gs))
                ok4 = True
            if not genes:
                genes, gtf_used = gs, src
    new_path = os.path.join(RES, "star_salmon", "salmon.merged.gene_counts.tsv")
    pub_path = os.path.join(PUB, "star_salmon", "salmon.merged.gene_counts.tsv")
    have_counts = os.path.isfile(new_path) and os.path.isfile(pub_path)
    if have_counts:
        ns, nc, nn = read_counts(new_path)
        ps, pc, pn = read_counts(pub_path)
        same_ids = set(nc) == set(pc)
        l4.append("- gene ids in this run's salmon table %d, in the published table %d, identical sets: %s"
                  % (len(nc), len(pc), same_ids))
        d4["gene_ids_identical_to_published"] = same_ids
        ok4 = ok4 and same_ids
        missing = [s for s in SAMPLES if s not in ns]
        l4.append("- samples in this run's salmon table: %d of 11%s" % (len(SAMPLES) - len(missing),
                                                                           (" (missing: %s)" % ", ".join(missing)) if missing else ""))
        ok4 = ok4 and not missing
    else:
        l4.append("- salmon.merged.gene_counts.tsv missing (this run: %s, published: %s)" % (os.path.isfile(new_path), os.path.isfile(pub_path)))
        ok4 = False
    record("4", "GENCODE v44 reference, same gene set as the published run", ok4, l4, d4)

    if not have_counts:
        for k, t in (("5", "Clip effect"), ("6", "Both-bin effect in IL69B"), ("7", "Rat leak from counts"), ("8", "IL64B and N168B graft reads")):
            record(k, t, False, ["Not computed: count tables missing."], hard=False)
    else:
        gene_info = genes

        def gtype(g):
            return gene_info.get(g, ("", "", ""))[1]

        def gchr(g):
            return gene_info.get(g, ("", "", ""))[0]

        def gname(g):
            return nn.get(g) or pn.get(g) or gene_info.get(g, ("", "", ""))[2]

        rp_rx = re.compile(r"^(RP[LS]\d{1,2}(A|AL|L|L\d)?|RPS4[XY]\d?|RPLP[0-2]|RPSA)$")  # cytosolic RP genes, no pseudogenes/kinases

        def is_rp(g):
            n = gname(g)
            return bool(rp_rx.match(n)) and not n.startswith("RPS6K") and gtype(g) == "protein_coding"

        sets = {}
        for st_ in ("AC", "MES", "NPC", "OPC"):
            sets[st_] = set(read_list(os.path.join(GS, "neftel_%s.txt" % st_)))

        def compare(a_col, b_col, label_a, label_b, fname):
            """a vs b, genes at CPM > 10 in both: log2(CPM_a / CPM_b)"""
            ca, ta = cpm(a_col)
            cb, tb = cpm(b_col)
            keep = [g for g in ca if ca[g] > 10 and cb.get(g, 0) > 10]
            lr = {g: math.log(ca[g] / cb[g], 2) for g in keep}
            rho = spearman([ca[g] for g in keep], [cb[g] for g in keep])
            n = len(keep)
            s01 = sum(1 for g in keep if abs(lr[g]) > 0.1) / float(n) if n else float("nan")
            s02 = sum(1 for g in keep if abs(lr[g]) > 0.2) / float(n) if n else float("nan")
            rp = [g for g in ca if is_rp(g)]
            rp_a = sum(ca[g] for g in rp)
            rp_b = sum(cb.get(g, 0) for g in rp)
            rp_lr = [lr[g] for g in rp if g in lr]
            sc = {}
            for st_, gs in sets.items():
                ids = [g for g in ca if gname(g) in gs]
                if ids:
                    sa = sum(math.log(ca[g] + 1, 2) for g in ids) / len(ids)
                    sb = sum(math.log(cb.get(g, 0) + 1, 2) for g in ids) / len(ids)
                    sc[st_] = {"n": len(ids), label_a: sa, label_b: sb, "diff": sa - sb}
            with open(os.path.join(TAB, fname), "w") as fh:
                fh.write("gene_id\tgene_name\tgene_type\tcount_%s\tcount_%s\tcpm_%s\tcpm_%s\tlog2_ratio\n" % (label_a, label_b, label_a, label_b))
                for g in sorted(keep, key=lambda x: -abs(lr[x])):
                    fh.write("%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.5f\n" % (g, gname(g), gtype(g), a_col.get(g, 0.0), b_col.get(g, 0.0), ca[g], cb[g], lr[g]))
            top = sorted(keep, key=lambda x: -abs(lr[x]))[:15]
            return {"genes_cpm10_both": n, "spearman": rho, "share_abs_log2_gt_0.1": s01, "share_abs_log2_gt_0.2": s02,
                    "median_log2_ratio": median(list(lr.values())), "total_a": ta, "total_b": tb,
                    "rp_genes": len(rp), "rp_cpm_a": rp_a, "rp_cpm_b": rp_b,
                    "rp_block_log2_ratio": math.log(rp_a / rp_b, 2) if rp_a > 0 and rp_b > 0 else float("nan"),
                    "rp_median_gene_log2_ratio": median(rp_lr), "neftel_scores": sc,
                    "top_changed": [(gname(g), round(lr[g], 3), round(ca[g], 1), round(cb[g], 1)) for g in top]}

        def cmp_lines(r, la, lb):
            out = ["- genes at CPM > 10 in both: %d; Spearman %s; median log2(%s/%s) %s" % (
                r["genes_cpm10_both"], fmt(r["spearman"]), la, lb, fmt(r["median_log2_ratio"])),
                "- share with |log2 ratio| > 0.1: %s; > 0.2: %s" % (fmt(100 * r["share_abs_log2_gt_0.1"], 2) + " %",
                                                                   fmt(100 * r["share_abs_log2_gt_0.2"], 2) + " %"),
                "- ribosomal-protein block (%d protein-coding RPL*/RPS* genes): CPM %s -> %s (%s) = log2 %s; median per-gene log2 ratio %s" % (
                    r["rp_genes"], fmt(r["rp_cpm_b"], 0), fmt(r["rp_cpm_a"], 0), lb + " -> " + la,
                    fmt(r["rp_block_log2_ratio"]), fmt(r["rp_median_gene_log2_ratio"]))]
            for st_, v in sorted(r["neftel_scores"].items()):
                out.append("- Neftel %s set (%d genes, mean log2(CPM+1)): %s %s, %s %s, difference %s" % (
                    st_, v["n"], lb, fmt(v[lb], 3), la, fmt(v[la], 3), fmt(v["diff"], 3)))
            out.append("- most changed (gene, log2 ratio, CPM %s, CPM %s): %s" % (la, lb, "; ".join(
                "%s %s (%s/%s)" % t for t in r["top_changed"])))
            return out

        new = lambda s: column(ns, nc, s)  # noqa
        pub = lambda s: column(ps, pc, s)  # noqa

        # ------------------------------------------------------------------------------ 5. clip effect
        r5 = compare(new("IL69B_fullclip"), pub("IL69B"), "IL69B_fullclip", "IL69B_pub", "clip_effect_IL69B.tsv")
        material = r5["share_abs_log2_gt_0.2"] > 0.01
        r5["clip_material"] = material
        C["clip_material"] = material
        l5 = ["Same records as the published IL69B (graft + both), R2 now clipped by 3 nt; everything else identical. "
              "Rule fixed before the run: the clip is **material** if more than 1 % of genes at CPM > 10 in both move by "
              "|log2 ratio| > 0.2 (this folds in run-to-run noise, so it is an upper bound on the clip effect). "
              "Material triggers the host rerun (SELECTION s4).", ""]
        l5 += cmp_lines(r5, "IL69B_fullclip", "IL69B_pub")
        l5 += ["", "**Decision: clip %s** (%s of genes beyond |log2| 0.2; threshold 1 %%)." % (
            "MATERIAL" if material else "not material", fmt(100 * r5["share_abs_log2_gt_0.2"], 2) + " %")]
        record("5", "Clip effect (IL69B_fullclip vs published IL69B)", None, l5, r5, hard=False)

        # ------------------------------------------------------------------------------ 6. both-bin effect
        r6 = compare(new("IL69B_graft"), new("IL69B_fullclip"), "IL69B_graft", "IL69B_fullclip", "both_bin_effect_IL69B.tsv")
        e = XENG["IL69B"]
        lost = 1 - r6["total_a"] / r6["total_b"] if r6["total_b"] else float("nan")
        r6["count_share_removed"] = lost
        r6["pair_share_removed_expected"] = e[2] / float(e[1] + e[2])
        l6 = ["Graft-only against graft + both, both clipped: the effect of dropping xengsort's `both` bin (reads whose "
              "k-mers match rat and human equally).", "",
              "- assigned counts removed: %s %% (the both bin is %s %% of IL69B's human pairs)" % (
                  fmt(100 * lost, 2), fmt(100 * r6["pair_share_removed_expected"], 2))]
        l6 += cmp_lines(r6, "IL69B_graft", "IL69B_fullclip")
        # descriptive, per tumour: ribosomal-protein share of counts, published (graft+both, unclipped) vs graft-only
        l6 += ["", "Ribosomal-protein block share of counts, published (graft+both, unclipped) vs this run (graft-only, "
                   "clipped); descriptive, the clip is folded in (see check 5):", "",
               "| library | arm | RP share published | RP share graft | log2 change |", "|---|---|---|---|---|"]
        rp_tab = {}
        for t in TUMOURS + ["C2B", "N269B"]:
            cp_, _ = cpm(pub(t))
            cg_, _ = cpm(new(t + "_graft"))
            a = sum(v for g, v in cp_.items() if is_rp(g)) / 1e6
            b = sum(v for g, v in cg_.items() if is_rp(g)) / 1e6
            rp_tab[t] = (a, b)
            arm = "primary" if t in PRIMARY else ("recurrent" if t in RECURRENT else "other")
            l6.append("| %s | %s | %s %% | %s %% | %s |" % (t, arm, fmt(100 * a, 2), fmt(100 * b, 2),
                                                           fmt(math.log(b / a, 2) if a > 0 and b > 0 else float("nan"), 3)))
        r6["rp_share_by_library"] = rp_tab
        record("6", "Both-bin effect in IL69B", None, l6, r6, hard=False)

        # ------------------------------------------------------------------------------ 7. rat leak from counts
        c2_new, c2_pub = new("C2B_graft"), pub("C2B")
        rg = {g: c2_new.get(g, 0.0) / c2_pub[g] for g in c2_pub if c2_pub[g] >= 50}
        elig = [g for g, r in rg.items() if r >= 0.2]
        l7 = ["C2B is pure human, so r_g = C2B_graft / published C2B is the share of gene g's human reads that xengsort "
              "keeps as graft (it depends on the read sequence, not on the sample). For a tumour, rat_g = P_g - G_g / r_g, "
              "with P the published count (graft + both, unclipped) and G this run's graft count. The clip enters r_g and "
              "G alike, so it cancels to first order. Genes: r_g >= 0.2 and published C2B >= 50 (%d genes; r_g median %s)."
              % (len(elig), fmt(median([rg[g] for g in elig]), 3)), "",
              "| library | arm | eligible genes' share of published counts | rat share (signed sum) | rat share (genes clipped at 0) | prior range |",
              "|---|---|---|---|---|---|"]
        r7 = {"eligible_genes": len(elig), "by_library": {}}
        for t in TUMOURS + ["N269B", "IL64B", "N168B", "C2B"]:
            P, G = pub(t), new(t + "_graft")
            totP = sum(P.values())
            sP = sum(P.get(g, 0.0) for g in elig)
            rat = [P.get(g, 0.0) - G.get(g, 0.0) / rg[g] for g in elig]
            s_signed = sum(rat) / sP if sP > 0 else float("nan")
            s_clip = sum(max(0.0, x) for x in rat) / sP if sP > 0 else float("nan")
            arm = "primary" if t in PRIMARY else ("recurrent" if t in RECURRENT else "control/other")
            pr = PRIOR_LEAK.get(arm)
            inside = (pr[0] <= 100 * s_signed <= pr[1]) if pr else None
            r7["by_library"][t] = {"arm": arm, "eligible_share": sP / totP if totP else float("nan"),
                                   "rat_share_signed": s_signed, "rat_share_clipped": s_clip, "inside_prior": inside}
            l7.append("| %s | %s | %s %% | %s %% | %s %% | %s |" % (
                t, arm, fmt(100 * sP / totP if totP else float("nan"), 1), fmt(100 * s_signed, 2), fmt(100 * s_clip, 2),
                ("%.1f-%.1f %% (%s)" % (pr[0], pr[1], "inside" if inside else "outside")) if pr else "-"))
        for arm, libs in (("primary", PRIMARY), ("recurrent", RECURRENT)):
            v = [r7["by_library"][t]["rat_share_signed"] for t in libs]
            r7["%s_range" % arm] = [min(v), max(v)]
        l7 += ["", "Primary %s-%s %%, recurrent %s-%s %% (signed). Prior (EVAL_transcripts, xengsort rate model): primary "
                   "1.9-5.1 %%, recurrent 5.1-8.0 %%. C2B reads zero by construction; IL64B and N168B show what the estimator "
                   "gives when almost all of the human file is rat." % (
                       fmt(100 * r7["primary_range"][0], 2), fmt(100 * r7["primary_range"][1], 2),
                       fmt(100 * r7["recurrent_range"][0], 2), fmt(100 * r7["recurrent_range"][1], 2))]
        with open(os.path.join(TAB, "rat_leak_rg.tsv"), "w") as fh:
            fh.write("gene_id\tgene_name\tc2b_published\tc2b_graft\tr_g\n")
            for g in sorted(elig):
                fh.write("%s\t%s\t%.3f\t%.3f\t%.5f\n" % (g, gname(g), c2_pub[g], c2_new.get(g, 0.0), rg[g]))
        record("7", "Rat leak estimated from counts", None, l7, r7, hard=False)

        # ------------------------------------------------------------------------------ 8. IL64B / N168B graft reads
        c2 = new("C2B_graft")
        c2c, c2t = cpm(c2)
        y_ids = [g for g in c2 if gchr(g) == "chrY" and not g.endswith("_PAR_Y")]
        y_cpm_c2b = sum(c2c[g] for g in y_ids)
        u_c2b = pct(star_log(os.path.join(RES, "star_salmon", "log", "C2B_graft.Log.final.out")).get("Uniquely mapped reads %", "nan"))
        hist_rx = re.compile(r"^(H1-\d+|H2A[CJWXZ]\w*|H2B[CUW]\w*|H3C\d+|H3-\d\w*|H4C\d+|H4-16)$")
        rrna_types = {"rRNA", "rRNA_pseudogene", "Mt_rRNA"}
        small_types = {"misc_RNA", "snRNA", "snoRNA", "scaRNA", "vault_RNA", "Y_RNA"}
        def qshare(col):
            tot = sum(col.values())
            mt = sum(v for g, v in col.items() if gname(g).startswith("MT-"))
            rr = sum(v for g, v in col.items() if gtype(g) in rrna_types)
            hi = sum(v for g, v in col.items() if hist_rx.match(gname(g) or ""))
            sm = sum(v for g, v in col.items() if gtype(g) in small_types)
            return tot, mt, rr, hi, sm, ((mt + rr + hi + sm) / tot if tot else float("nan"))

        q_c2b = qshare(c2)[5]
        q_rat = max(0.5, 2 * q_c2b)
        q_hum = max(0.3, 1.5 * q_c2b)
        l8 = ["Rule, fixed before the run (U = STAR unique %%, C2B_graft U = %s %%; Q = share of assigned counts in MT- genes, "
              "rRNA (incl. pseudogenes, Mt_rRNA), replication histones and conserved small RNAs (misc_RNA: 7SL/7SK/Y RNA/vault, "
              "snRNA, snoRNA, scaRNA), C2B_graft Q = %s; Y = reads on non-PAR chrY genes, Y_exp = C2B_graft's chrY CPM (%s) x "
              "the library's counts; rho = Spearman with C2B_graft over genes at C2B CPM > 10):" % (fmt(u_c2b, 2), fmt(q_c2b, 3), fmt(y_cpm_c2b, 1)),
              "- **rat mis-sorted as graft**: (U < 0.5 x U_C2B and too short + too many mismatches > 20 %%) or Q > max(0.5, 2 x Q_C2B) = %s;" % fmt(q_rat, 3),
              "- **non-U251 human**: U >= 0.8 x U_C2B, Q < max(0.3, 1.5 x Q_C2B) = %s, >= 3,000 genes with >= 5 counts, and Y <= max(4, 0.2 x Y_exp);" % fmt(q_hum, 3),
              "- **U251 human (tumour cells in the control)**: U >= 0.8 x U_C2B, Q below the same bound, Y >= 0.5 x Y_exp and rho >= 0.7 "
              "(a fourth outcome, added because U251 is male and the host female, so Y separates it from any other human source);",
              "- otherwise **unclassified** (triggers the competitive alignment, SELECTION s4).", ""]
        r8 = {}
        for lib in CONTROLS:
            s = lib + "_graft"
            col = new(s)
            cc, tot = cpm(col)
            sl = star_log(os.path.join(RES, "star_salmon", "log", "%s.Log.final.out" % s))
            U_ = pct(sl.get("Uniquely mapped reads %", "nan"))
            short = pct(sl.get("% of reads unmapped: too short", "nan"))
            mism = pct(sl.get("% of reads unmapped: too many mismatches", "nan"))
            multi = pct(sl.get("% of reads mapped to multiple loci", "nan"))
            _, mt, rr, hi, sm, Q = qshare(col)
            yobs = sum(col.get(g, 0.0) for g in y_ids)
            yexp = y_cpm_c2b * tot / 1e6
            xist = sum(col[g] for g in col if gname(g) == "XIST")
            det = sum(1 for v in col.values() if v >= 5)
            keep = [g for g in c2c if c2c[g] > 10]
            rho = spearman([cc.get(g, 0.0) for g in keep], [c2c[g] for g in keep]) if tot > 0 else float("nan")
            if (U_ < 0.5 * u_c2b and short + mism > 20) or Q > q_rat:
                verdict = "rat mis-sorted as graft"
            elif U_ >= 0.8 * u_c2b and Q < q_hum and det >= 3000 and yobs <= max(4.0, 0.2 * yexp):
                verdict = "non-U251 human"
            elif U_ >= 0.8 * u_c2b and Q < q_hum and yobs >= 0.5 * yexp and rho >= 0.7:
                verdict = "U251 human (tumour cells in the control)"
            else:
                verdict = "unclassified"
            top = sorted(col, key=lambda g: -col[g])[:50]
            r8[lib] = {"verdict": verdict, "star_unique_pct": U_, "too_short_pct": short, "too_many_mismatches_pct": mism,
                       "multi_pct": multi, "assigned_counts": tot, "genes_ge5": det, "share_MT": mt / tot if tot else None,
                       "share_rRNA": rr / tot if tot else None, "share_histone": hi / tot if tot else None,
                       "share_small_conserved_RNA": sm / tot if tot else None, "Q": Q, "Q_C2B": q_c2b,
                       "chrY_reads": yobs, "chrY_expected_if_U251": yexp, "XIST": xist, "spearman_vs_C2B": rho,
                       "top50": [(gname(g), gtype(g), round(col[g], 1)) for g in top]}
            l8 += ["### %s_graft: %s" % (lib, verdict), "",
                   "- STAR: unique %s %%, multi %s %%, unmapped too short %s %%, too many mismatches %s %%" % (
                       fmt(U_, 2), fmt(multi, 2), fmt(short, 2), fmt(mism, 2)),
                   "- assigned counts %s; genes with >= 5 counts %d; Spearman vs C2B_graft %s" % (fmt(tot, 0), det, fmt(rho, 3)),
                   "- shares of counts: MT- %s %%, rRNA %s %%, histone %s %%, conserved small RNA %s %% (Q %s; C2B_graft %s)" % (
                       fmt(100 * mt / tot if tot else float("nan"), 2), fmt(100 * rr / tot if tot else float("nan"), 2),
                       fmt(100 * hi / tot if tot else float("nan"), 2), fmt(100 * sm / tot if tot else float("nan"), 2),
                       fmt(Q, 3), fmt(q_c2b, 3)),
                   "- chrY (non-PAR) reads %s, expected if these were U251 %s; XIST %s" % (fmt(yobs, 1), fmt(yexp, 1), fmt(xist, 1)),
                   "- top 50 genes: " + "; ".join("%s (%s) %s" % t for t in r8[lib]["top50"]), ""]
            with open(os.path.join(TAB, "%s_graft_top50.tsv" % lib), "w") as fh:
                fh.write("gene_name\tgene_type\tcount\n")
                for t in r8[lib]["top50"]:
                    fh.write("%s\t%s\t%s\n" % t)
        C["control_graft_verdicts"] = {k: v["verdict"] for k, v in r8.items()}
        record("8", "What IL64B's and N168B's graft reads are", None, l8, r8, hard=False)

    # ---------------------------------------------------------------------------------- 9. saved index
    th = trace_hashes(RES)
    idx = os.path.join(RES, "genome", "index")
    star_need = ["Genome", "SA", "SAindex", "genomeParameters.txt", "chrName.txt", "sjdbList.out.tab"]
    sal_need = ["versionInfo.json", "info.json", "pos.bin", "seq.bin", "ctable.bin", "ctg_offsets.bin", "mphf.bin",
                "rank.bin", "refseq.bin", "reflengths.bin", "refAccumLengths.bin", "complete_ref_lens.bin"]
    l9, ok9, d9 = [], True, {}
    for proc, sub, need in (("STAR_GENOMEGENERATE", "star", star_need), ("SALMON_INDEX", "salmon", sal_need)):
        pubd = os.path.join(idx, sub)
        have = tree_listing(pubd) if os.path.isdir(pubd) else {}
        miss = [f for f in need if have.get(f, 0) <= 0]
        wd = None
        for tag, h in th.get(proc, {}).items():
            wd = work_dir(WORK, h)
        wsrc = os.path.join(wd, sub) if wd else None
        wlist = tree_listing(wsrc) if wsrc and os.path.isdir(wsrc) else None
        diff = None
        if wlist is not None:
            diff = sorted(f for f in wlist if have.get(f) != wlist[f])
            if diff and args.repair_index:
                for f in diff:
                    dst = os.path.join(pubd, f)
                    os.makedirs(os.path.dirname(dst), exist_ok=True)
                    shutil.copy2(os.path.join(wsrc, f), dst)
                l9.append("- %s: %d files differed from the task's own output and were copied from `%s`: %s" % (
                    sub, len(diff), wsrc, ", ".join(diff[:20])))
                have = tree_listing(pubd)
                miss = [f for f in need if have.get(f, 0) <= 0]
                diff = sorted(f for f in wlist if have.get(f) != wlist[f])
        good = not miss and (diff is None or not diff)
        ok9 = ok9 and good
        d9[sub] = {"repaired_from_work": bool(wlist is not None and args.repair_index and l9 and l9[-1].startswith("- %s: " % sub) and "were copied" in l9[-1]),
                   "files": len(have), "bytes": sum(v for v in have.values() if v > 0), "missing_required": miss,
                   "work_dir": wsrc, "differs_from_work": diff}
        l9.append("- %s: %d files, %s GB; required missing: %s; identical to the task's output in work: %s" % (
            sub, len(have), fmt(sum(v for v in have.values() if v > 0) / 1e9, 1), ", ".join(miss) or "none",
            "not checkable (work gone)" if diff is None else ("yes" if not diff else "NO (%d differ)" % len(diff))))
    record("9", "Human STAR and salmon indices saved complete", ok9, l9, d9)

    # ---------------------------------------------------------------------------------- 10. no empty outputs
    l10, ok10, d10 = [], True, {}
    for s in SAMPLES:
        q = os.path.join(RES, "star_salmon", s, "quant.genes.sf")
        nz = None
        if os.path.isfile(q):
            with open(q) as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                nz = sum(1 for r in rd if float(r["NumReads"]) > 0)
        lib = s.rsplit("_", 1)[0]
        need = lib in BIG
        good = (nz is not None and nz > 10000) if need else (nz is not None)
        ok10 = ok10 and good
        d10[s] = {"genes_numreads_gt0": nz, "gated": need}
        l10.append("- %s: genes with NumReads > 0 = %s (%s)" % (s, nz, ("needs > 10,000" if need else "reported") + ("" if good else " FAIL")))
    tg = th.get("TRIMGALORE", {})
    empties, checked = [], 0
    for tag, h in tg.items():
        wd = work_dir(WORK, h)
        if not wd:
            continue
        for f in glob.glob(os.path.join(wd, "*_val_[12].fq.gz")):
            checked += 1
            if os.path.getsize(f) <= 20:
                empties.append(f)
    l10.append("- trimmed FASTQs in work checked: %d (expected 22); empty: %s" % (checked, ", ".join(empties) or "none"))
    if checked == 0:
        l10.append("  (work not present; check 2 already requires STAR input = pairs TrimGalore kept, which an empty file fails)")
    if empties or (checked and checked != 22):
        ok10 = False
    d10["trimmed_checked"] = checked
    d10["trimmed_empty"] = empties
    record("10", "No empty outputs", ok10, l10, d10)
    return finish(C, md, D)


def finish(C, md, D):
    hard = [k for k, v in C["checks"].items() if v["hard"]]
    failed = [k for k in hard if not v_pass(C, k)]
    C["all_hard_pass"] = bool(hard) and not failed and len(hard) >= 7
    C["failed_hard"] = failed
    C["finished"] = now()
    head = ["", "**all_hard_pass: %s**%s" % (C["all_hard_pass"], (" (failed: %s)" % ", ".join(failed)) if failed else ""), ""]
    if "clip_material" in C:
        head.append("Clip decision: %s. Control verdicts: %s." % (
            "MATERIAL" if C["clip_material"] else "not material",
            ", ".join("%s %s" % kv for kv in sorted(C.get("control_graft_verdicts", {}).items()))))
        head.append("")
    md[4:4] = head
    with open(os.path.join(D, "checks.json.part"), "w") as fh:
        json.dump(C, fh, indent=1, default=str)
    os.replace(os.path.join(D, "checks.json.part"), os.path.join(D, "checks.json"))
    with open(os.path.join(D, "CHECKS.md.part"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    os.replace(os.path.join(D, "CHECKS.md.part"), os.path.join(D, "CHECKS.md"))
    print("\n".join(md))
    return 0 if C["all_hard_pass"] else 2


def v_pass(C, k):
    return bool(C["checks"][k]["pass"])


if __name__ == "__main__":
    sys.exit(main())
