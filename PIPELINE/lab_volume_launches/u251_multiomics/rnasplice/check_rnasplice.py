#!/usr/bin/env python3
"""Checks for rnasplice (nf-core/rnasplice 1.0.4, Nextflow 23.10.0) on the six graft-only tumours.

Writes $D/CHECKS.md and $D/checks.json. If every HARD check passes and every input the null needs exists, it also
writes $D/null/inputs.env (container images and published files the relabelling null re-uses). Standard library
only (Python >= 3.6); the strand check calls samtools inside the samtools 1.21 image.

HARD checks (all must pass for "hard_pass"; the null launcher refuses to run without it):
  0 run            Nextflow exit 0, Nextflow 23.10.0, pipeline commit = tag 1.0.4, params actually applied
  1 trace          every expected process ran for all six samples; none of the forbidden ones ran
  2 commands       strandedness, read length, unpaired statistics and the DRIMSeq filter in the task scripts
  3 inputs         TrimGalore saw N_graft pairs per library (xengsort counts); STAR saw what TrimGalore wrote
  4 read lengths   FastQC on the input: R1 max 101, R2 max 98 (the 3-nt clip was applied once, upstream)
  5 strand         salmon compatible-fragment ratio and an independent R1-orientation count on named genes
  6 completeness   six samples in every per-sample output; the rMATS prep holds each BAM exactly once
  7 null inputs    every file and container the null re-uses exists
Check 5 has a WARN band (0.90-0.97): a strand error reads near 0.5 (unstranded) or near 0.03 (reversed), so below
0.90 fails; 0.90-0.97 is short of the published RSeQC 97.4-98.4 % and is flagged for a look but does not block.
If the salmon pseudo-alignment branch ran (pseudo_aligner null not honoured), check 1 reports it (not a failure):
STAR+Salmon was declared the analysis of record in README.md before the run, and the null uses only that branch.
"""
import argparse
import csv
import datetime
import glob
import gzip
import json
import os
import re
import subprocess
import sys
import zipfile

M = "/rs/rs_grp_oschome/go2432/u251_multiomics"
C = "/rs/rs_grp_oschome/go2432/singularity_cache"
PIPE_SHA = "1d0494ae3402d1a46e0adadad24f81a0ff855c77"
NF_WANT = "23.10.0"
PRIMARY = ["IL67B", "IL68B", "IL69B"]
RECURRENT = ["IL66B", "NL70B", "NL71B"]
SAMPLES = [l + "_graft" for l in PRIMARY + RECURRENT]
# xengsort graft pairs (U/ANALYSIS/xengsort_out/<lib>.txt; asserted again by prep_reads check b)
N_GRAFT = {"IL66B": 25459283, "IL67B": 31886295, "IL68B": 32225887, "IL69B": 49327614,
           "NL70B": 31293679, "NL71B": 29210358}
# params that must have been applied (read back from results/pipeline_info/params_*.json)
PARAMS_WANT = {
    "source": "fastq", "aligner": "star_salmon", "pseudo_aligner": None, "gencode": True, "igenomes_ignore": True,
    "rmats": True, "rmats_read_len": 100, "rmats_paired_stats": False, "rmats_novel_splice_site": False,
    "dexseq_exon": True, "edger_exon": False, "dexseq_dtu": True, "dtu_txi": "dtuScaledTPM",
    "min_samps_gene_expr": 6, "min_gene_expr": 10, "min_samps_feature_expr": 3, "min_feature_expr": 10,
    "min_samps_feature_prop": 3, "min_feature_prop": 0.1, "suppa": True, "diffsplice_paired": False,
    "diffsplice_method": "empirical", "diffsplice_gene_correction": True, "sashimi_plot": False,
    "save_dexseq_annotation": True, "clip_r2": None,
}
# genes for the independent strand count; strands and exons come from the GTF, never typed in
STRAND_GENES = ["GAPDH", "ACTB", "EEF1A1", "VIM", "FTL", "TPT1", "RPLP0", "EEF2", "PKM", "ENO1", "LDHA", "TUBB",
                "HSP90AB1", "B2M", "MYL6", "CD44", "SPARC", "FN1", "TMSB10", "RPS3"]
SAM_IMG_DEFAULT = os.path.join(C, "depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img")


def now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_text(path):
    with open(path, "r", errors="replace") as fh:
        return fh.read()


def norm_cmd(txt):
    return re.sub(r"\s+", " ", txt.replace("\\\n", " ")).strip()


def fastqc_zip(path):
    out = {}
    with zipfile.ZipFile(path) as z:
        name = [n for n in z.namelist() if n.endswith("fastqc_data.txt")][0]
        for line in z.read(name).decode("utf-8", "replace").splitlines():
            if line.startswith("Total Sequences\t"):
                out["total"] = int(line.split("\t")[1])
            elif line.startswith("Sequence length\t"):
                parts = [int(p) for p in line.split("\t")[1].split("-")]
                out["minlen"], out["maxlen"] = min(parts), max(parts)
    return out


def trim_report(path):
    out = {"cmd": None, "total": None, "pairs_removed": None}
    if not os.path.isfile(path):
        return out
    for line in read_text(path).splitlines():
        if line.startswith("Command line parameters:"):
            out["cmd"] = line.split(":", 1)[1].strip()
        elif line.startswith("Total reads processed:"):
            out["total"] = int(line.split(":", 1)[1].strip().replace(",", ""))
        elif line.startswith("Number of sequence pairs removed because at least one read was shorter"):
            out["pairs_removed"] = int(line.split(":", 1)[1].strip().split()[0])
    return out


def star_log(path):
    out = {}
    if os.path.isfile(path):
        for line in read_text(path).splitlines():
            if "|" in line:
                k, v = line.split("|", 1)
                out[k.strip()] = v.strip()
    return out


def read_traces(res):
    """-> list of rows (dict) of all traces, and {full_name: latest COMPLETED/CACHED row}"""
    rows, ok = [], {}
    for tr in sorted(glob.glob(os.path.join(res, "pipeline_info", "execution_trace_*.txt"))):
        with open(tr) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                row["_trace"] = os.path.basename(tr)
                rows.append(row)
                if row.get("status") in ("COMPLETED", "CACHED"):
                    ok[row.get("name", "")] = row
    return rows, ok


def proc_of(name):
    m = re.match(r"^(\S+?)(?: \((.*)\))?$", name)
    return (m.group(1) if m else name), ((m.group(2) if m else "") or "")


def work_dir(work, h):
    if not h or "/" not in h:
        return None
    a, b = h.split("/", 1)
    hits = glob.glob(os.path.join(work, a, b + "*"))
    return hits[0] if len(hits) == 1 else None


def image_of(wd):
    """the .img the task ran in, from its .command.run"""
    p = os.path.join(wd, ".command.run") if wd else ""
    if not p or not os.path.isfile(p):
        return None
    hits = re.findall(r"(/\S+?\.img)\b", read_text(p))
    return hits[0] if hits else None


def gtf_exons(gtf, names):
    """gene_name -> (chrom, strand, [(start0, end)]) for protein-coding genes on chr1-22,X, from exon lines"""
    want, out = set(names), {}
    rx_n = re.compile(r'gene_name "([^"]+)"')
    rx_t = re.compile(r'gene_type "([^"]+)"')
    op = gzip.open if gtf.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\texon\t" not in line:
                continue
            m = rx_n.search(line)
            if not m or m.group(1) not in want:
                continue
            f = line.split("\t")
            if f[2] != "exon" or "_PAR_Y" in f[8]:
                continue
            t = rx_t.search(f[8])
            if not t or t.group(1) != "protein_coding":
                continue
            g = out.setdefault(m.group(1), [f[0], f[6], []])
            if g[0] != f[0] or g[1] != f[6]:
                g.append("CONFLICT")
                continue
            g[2].append((int(f[3]) - 1, int(f[4])))
    return {k: v for k, v in out.items() if "CONFLICT" not in v}


def sx(img, args, binds):
    cmd = ["singularity", "exec"] + sum([["-B", b] for b in binds], []) + [img] + args
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return p.returncode, p.stdout.strip(), p.stderr.strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", default=os.path.join(M, "rnasplice"))
    ap.add_argument("--gtf", default=os.path.join(M, "refs", "gencode_v44", "GRCh38.primary_assembly.annotation.gtf"))
    ap.add_argument("--sam-img", default=SAM_IMG_DEFAULT)
    args = ap.parse_args()
    D = args.d
    RES = os.path.join(D, "results")
    WORK = os.path.join(D, "work")
    NULL = os.path.join(D, "null")
    GTF = args.gtf
    SAM_IMG = args.sam_img
    SS = os.path.join(RES, "star_salmon")
    Cj = {"started": now(), "results": RES, "checks": {}}
    md = ["# rnasplice checks", "",
          "Written %s by check_rnasplice.py. Results: `%s`. Nominal splicing counts are not reported here: they mean "
          "nothing until ranked against the 10-split null (`null/results/NULL_REPORT.md`)." % (now(), RES), ""]

    def record(key, title, status, lines, data=None):
        """status: PASS, FAIL, WARN (does not block) or INFO"""
        Cj["checks"][key] = {"title": title, "status": status, "data": data or {}}
        md.extend(["## %s. %s: %s" % (key, title, status), ""] + lines + [""])

    # ------------------------------------------------------------------ 0 run
    lines, ok = [], True
    nfx = read_text(os.path.join(D, "NEXTFLOW_EXIT")).strip() if os.path.isfile(os.path.join(D, "NEXTFLOW_EXIT")) else "missing"
    lines.append("- NEXTFLOW_EXIT: `%s` (want 0)" % nfx)
    ok &= nfx == "0"
    sha = read_text(os.path.join(D, "PIPELINE_COMMIT")).strip() if os.path.isfile(os.path.join(D, "PIPELINE_COMMIT")) else "missing"
    lines.append("- pipeline commit: `%s` (want `%s`, GitHub refs/tags/1.0.4)" % (sha, PIPE_SHA))
    ok &= sha == PIPE_SHA
    nflog = os.path.join(D, "launch", ".nextflow.log")
    nfv = None
    if os.path.isfile(nflog):
        with open(nflog, errors="replace") as fh:
            for i, line in enumerate(fh):
                m = re.search(r"version (\d+\.\d+\.\d+)", line)
                if m:
                    nfv = m.group(1)
                    break
                if i > 200:
                    break
    lines.append("- Nextflow version in launch/.nextflow.log: `%s` (want %s)" % (nfv, NF_WANT))
    ok &= nfv == NF_WANT
    pj = sorted(glob.glob(os.path.join(RES, "pipeline_info", "params_*.json")))
    applied = {}
    if pj:
        try:
            applied = json.loads(read_text(pj[-1]))
        except ValueError:
            applied = {}
    bad = []
    for k, v in PARAMS_WANT.items():
        got = applied.get(k, None)
        same = (got == v) or (isinstance(v, float) and isinstance(got, (int, float)) and abs(got - v) < 1e-12)
        if not same:
            bad.append("%s=%r (want %r)" % (k, got, v))
    lines.append("- params actually applied (`%s`): %s" % (os.path.basename(pj[-1]) if pj else "no params_*.json",
                                                          "all %d as intended" % len(PARAMS_WANT) if not bad else "MISMATCH " + "; ".join(bad)))
    ok &= bool(pj) and not bad
    record("0", "run, versions, applied params", "PASS" if ok else "FAIL", lines,
           {"nextflow_exit": nfx, "commit": sha, "nextflow": nfv, "param_mismatch": bad})

    # ------------------------------------------------------------------ 1 trace
    rows, done = read_traces(RES)
    by_proc = {}
    for name, row in done.items():
        p, tag = proc_of(name)
        by_proc.setdefault(p.split(":")[-1], []).append((name, tag, row))
    lines, ok = [], True
    want_per_sample = ["TRIMGALORE", "FASTQC", "STAR_ALIGN", "SAMTOOLS_SORT", "DEXSEQ_COUNT", "SALMON_QUANT_STAR"]
    want_once = ["STAR_GENOMEGENERATE", "DEXSEQ_ANNOTATION", "DEXSEQ_EXON", "CREATE_BAMLIST", "RMATS_PREP", "RMATS_POST",
                 "DRIMSEQ_FILTER", "DEXSEQ_DTU", "STAGER", "GENERATE_EVENTS_IOE", "GENERATE_EVENTS_IOI", "PSIPEREVENT",
                 "PSIPERISOFORM", "DIFFSPLICE_IOE", "DIFFSPLICE_IOI"]
    for p in want_per_sample:
        tags = sorted({t for _, t, _ in by_proc.get(p, [])})
        hit = [s for s in SAMPLES if any(t.startswith(s) for t in tags)]
        good = len(hit) == 6
        ok &= good
        lines.append("- %s: %d task(s) for %d of 6 samples %s" % (p, len(by_proc.get(p, [])), len(hit), "" if good else "<-- FAIL"))
    for p in want_once:
        n = len(by_proc.get(p, []))
        good = n >= 1
        ok &= good
        lines.append("- %s: %d %s" % (p, n, "" if good else "<-- FAIL (missing)"))
    forbidden = [n for n in done if re.search(r"EDGER|SUBREAD|MISO|UCSC_BED|BEDTOOLS", n)]
    salmon_branch = [n for n in done if re.search(r"SALMON_INDEX|SALMON_QUANT_SALMON|DRIMSEQ_DEXSEQ_DTU_SALMON:|SUPPA_SALMON:|TX2GENE_TXIMPORT_SALMON:", n)]
    lines.append("- forbidden processes (edgeR, featureCounts, MISO, bigWig): %s" % (", ".join(forbidden) or "none"))
    ok &= not forbidden
    failed = [r.get("name") for r in rows if r.get("status") == "FAILED"]
    lines.append("- FAILED task attempts across all traces (retried or not): %d%s" % (len(failed), (": " + "; ".join(failed[:10])) if failed else ""))
    if salmon_branch:
        lines.append("- **pseudo_aligner null was NOT honoured**: the salmon pseudo-alignment branch ran (%d tasks). "
                     "Pre-declared rule (README.md, before the run): STAR+Salmon (`results/star_salmon/`) is the analysis "
                     "of record; `results/salmon/` is never used for a claim, and the null uses only the STAR+Salmon files." % len(salmon_branch))
    else:
        lines.append("- salmon pseudo-alignment branch (SALMON_INDEX, SALMON_QUANT_SALMON, *_SALMON DTU/SUPPA): absent, "
                     "so `pseudo_aligner: null` was honoured and there is one result set")
    record("1", "trace: expected processes, nothing forbidden", "PASS" if ok else "FAIL", lines,
           {"salmon_branch_ran": bool(salmon_branch), "failed_attempts": len(failed)})

    # ------------------------------------------------------------------ 2 commands
    lines, ok = [], True
    cmd_of, img = {}, {}

    def cmds(p):
        out = []
        for name, tag, row in by_proc.get(p, []):
            # never read the pseudo-alignment branch (it only exists if pseudo_aligner null was not honoured)
            if re.search(r"(^|:)(DRIMSEQ_DEXSEQ_DTU_SALMON|SUPPA_SALMON|TX2GENE_TXIMPORT_SALMON|SALMON_QUANT_SALMON)(:| |$)", name):
                continue
            wd = work_dir(WORK, row.get("hash"))
            sh = os.path.join(wd, ".command.sh") if wd else None
            txt = norm_cmd(read_text(sh)) if sh and os.path.isfile(sh) else None
            out.append((name, tag, txt))
            if wd and p not in img:
                img[p] = image_of(wd)
        return out

    def need(p, must=(), must_not=(), regex=None):
        nonlocal ok
        got = cmds(p)
        if not got:
            lines.append("- %s: no task script found <-- FAIL" % p)
            ok = False
            return
        for name, tag, txt in got:
            if txt is None:
                lines.append("- %s (%s): .command.sh not found (work dir gone?) <-- FAIL" % (p, tag))
                ok = False
                continue
            miss = [s for s in must if s not in txt]
            extra = [s for s in must_not if re.search(s, txt)]
            rx_ok = True if regex is None else bool(re.search(regex, txt))
            good = not miss and not extra and rx_ok
            ok &= good
            lines.append("- %s%s: %s" % (p, " (%s)" % tag if tag else "",
                                          "ok" if good else "FAIL missing %s; present %s; pattern %s" % (miss, extra, "ok" if rx_ok else "absent")))
            cmd_of.setdefault(p, txt)

    need("RMATS_PREP", must=["--libType fr-firststrand", "--readLength 100", "--variable-read-length", "-t paired", "--task prep"],
         must_not=[r"--novelSS"])
    need("RMATS_POST", must=["--libType fr-firststrand", "--readLength 100", "--variable-read-length", "-t paired", "--task post"],
         must_not=[r"--paired-stats", r"--novelSS"])
    need("DEXSEQ_COUNT", must=["-s reverse", "-p yes", "-r pos", "-a 10", "-f bam"])
    need("DIFFSPLICE_IOE", must=["-m empirical", "-gc"], must_not=[r"(^| )-pa( |$)"])
    need("DIFFSPLICE_IOI", must=["-m empirical", "-gc"], must_not=[r"(^| )-pa( |$)"])
    need("SALMON_QUANT_STAR", must=["--libType=ISR"])
    need("TRIMGALORE", must=["--paired"], must_not=[r"(?i)--clip_r[12]", r"(?i)--three_prime_clip"])
    need("DRIMSEQ_FILTER", regex=r"run_drimseq_filter\.R \S+ \S+ \S+ 6 3 3 10 0\.1 10( |$)")
    need("STAR_ALIGN", must=["--twopassMode Basic", "--quantMode TranscriptomeSAM"])
    for p in ["DEXSEQ_EXON", "DEXSEQ_DTU", "STAGER", "SPLIT_FILES_TPM", "PSIPEREVENT"]:
        cmds(p)
    record("2", "commands: strand, read length, unpaired statistics, DRIMSeq filter", "PASS" if ok else "FAIL", lines)

    # ------------------------------------------------------------------ 3 inputs (xengsort counts -> TrimGalore -> STAR)
    lines, ok = [], True
    for s in SAMPLES:
        lib = s.split("_")[0]
        r1 = trim_report(os.path.join(RES, "trimgalore", s + "_1.fastq.gz_trimming_report.txt"))
        r2 = trim_report(os.path.join(RES, "trimgalore", s + "_2.fastq.gz_trimming_report.txt"))
        sl = star_log(os.path.join(SS, "log", s + ".Log.final.out"))
        star_in = int(sl.get("Number of input reads", "-1") or -1)
        want = N_GRAFT[lib]
        kept = (r2["total"] - r2["pairs_removed"]) if (r2["total"] is not None and r2["pairs_removed"] is not None) else None
        good = r1["total"] == want and r2["total"] == want and kept is not None and star_in == kept
        ok &= good
        lines.append("- %s: TrimGalore R1 %s, R2 %s (xengsort graft %d); pairs kept %s; STAR input %d; unique %s %s"
                     % (s, r1["total"], r2["total"], want, kept, star_in, sl.get("Uniquely mapped reads %", "?"), "" if good else "<-- FAIL"))
    record("3", "inputs: xengsort graft pairs -> TrimGalore -> STAR", "PASS" if ok else "FAIL", lines)

    # ------------------------------------------------------------------ 4 read lengths (FastQC on the input)
    lines, ok = [], True
    for s in SAMPLES:
        lib = s.split("_")[0]
        z1 = glob.glob(os.path.join(RES, "fastqc", s + "_1_fastqc.zip"))
        z2 = glob.glob(os.path.join(RES, "fastqc", s + "_2_fastqc.zip"))
        if not z1 or not z2:
            ok = False
            lines.append("- %s: FastQC zip missing <-- FAIL" % s)
            continue
        f1, f2 = fastqc_zip(z1[0]), fastqc_zip(z2[0])
        good = f1.get("maxlen") == 101 and f2.get("maxlen") == 98 and f1.get("total") == N_GRAFT[lib] == f2.get("total")
        ok &= good
        lines.append("- %s: R1 %s-%s nt, R2 %s-%s nt, %s/%s reads %s" % (s, f1.get("minlen"), f1.get("maxlen"), f2.get("minlen"),
                                                                  f2.get("maxlen"), f1.get("total"), f2.get("total"), "" if good else "<-- FAIL"))
    record("4", "read lengths: R1 max 101, R2 max 98 (clip applied once, upstream)", "PASS" if ok else "FAIL", lines)

    # ------------------------------------------------------------------ 5 strand
    lines, worst = [], 1.0
    sal = {}
    for s in SAMPLES:
        js = glob.glob(os.path.join(SS, "*", s, "lib_format_counts.json"))
        if not js:
            lines.append("- salmon %s: lib_format_counts.json missing" % s)
            worst = 0.0
            continue
        j = json.loads(read_text(js[0]))
        cfr = float(j.get("compatible_fragment_ratio", 0.0))
        isr, isf = float(j.get("ISR", 0) or 0), float(j.get("ISF", 0) or 0)
        share = isr / (isr + isf) if isr + isf > 0 else None
        sal[s] = {"expected_format": j.get("expected_format"), "compatible_fragment_ratio": cfr, "ISR_share": share,
                  "strand_mapping_bias": j.get("strand_mapping_bias"), "file": js[0]}
        worst = min(worst, cfr if j.get("expected_format") == "ISR" else 0.0)
        lines.append("- salmon %s: expected %s, compatible-fragment ratio %.4f, ISR/(ISR+ISF) %s, strand_mapping_bias %s"
                     % (s, j.get("expected_format"), cfr, "%.4f" % share if share is not None else "NA", j.get("strand_mapping_bias")))
    # independent: R1 orientation on the exons of named protein-coding genes (RSeQC infer_experiment logic)
    ind = {}
    ex = gtf_exons(GTF, STRAND_GENES) if os.path.isfile(GTF) else {}
    beds = {"+": os.path.join(D, "tables", "strand_genes_plus.bed"), "-": os.path.join(D, "tables", "strand_genes_minus.bed")}
    os.makedirs(os.path.join(D, "tables"), exist_ok=True)
    for st, path in beds.items():
        with open(path, "w") as fh:
            for g, (chrom, strand, exs) in sorted(ex.items()):
                if strand == st:
                    for a, b in sorted(set(exs)):
                        fh.write("%s\t%d\t%d\t%s\n" % (chrom, a, b, g))
    lines.append("- independent count: %d genes found in the GTF (%s plus, %s minus): %s" % (
        len(ex), sum(1 for v in ex.values() if v[1] == "+"), sum(1 for v in ex.values() if v[1] == "-"), ", ".join(sorted(ex))))
    if not os.path.isfile(SAM_IMG):
        lines.append("- samtools image missing: %s" % SAM_IMG)
        worst = 0.0
    else:
        for s in SAMPLES:
            bam = os.path.join(SS, s + "_sorted.bam")
            if not os.path.isfile(bam):
                lines.append("- %s: genome BAM missing (%s)" % (s, bam))
                worst = 0.0
                continue
            n = {}
            for st in ("+", "-"):
                for orient, flags in (("fwd", ["-f", "64", "-F", "2836"]), ("rev", ["-f", "80", "-F", "2820"])):
                    rc, out, err = sx(SAM_IMG, ["samtools", "view", "-c", "-M", "-q", "30"] + flags + ["-L", beds[st], bam],
                                      [os.path.dirname(M) + ":" + os.path.dirname(M)])
                    n[st + orient] = int(out) if rc == 0 and out.isdigit() else -1
            tot = sum(v for v in n.values() if v >= 0)
            anti = (n["+rev"] + n["-fwd"]) / tot if tot > 0 and min(n.values()) >= 0 else 0.0
            ind[s] = {"counts": n, "antisense_share": anti}
            worst = min(worst, anti)
            lines.append("- %s: R1 antisense to the gene %.4f (plus genes R1 fwd %d / rev %d; minus genes R1 fwd %d / rev %d)"
                         % (s, anti, n["+fwd"], n["+rev"], n["-fwd"], n["-rev"]))
    status = "PASS" if worst >= 0.97 else ("WARN" if worst >= 0.90 else "FAIL")
    lines.append("- worst of all strand measures: %.4f -> %s (PASS >= 0.97, the published RSeQC 97.4-98.4 %%; FAIL < 0.90)" % (worst, status))
    record("5", "strand: salmon ISR and R1 orientation on named genes", status, lines, {"salmon": sal, "independent": ind, "worst": worst})

    # ------------------------------------------------------------------ 6 completeness
    lines, ok = [], True
    cnt_dir = os.path.join(SS, "dexseq_exon", "counts")
    nlines = {}
    for s in SAMPLES:
        p = os.path.join(cnt_dir, s + ".clean.count.txt")
        nlines[s] = sum(1 for _ in open(p)) if os.path.isfile(p) and os.path.getsize(p) > 0 else 0
    good = all(v > 100000 for v in nlines.values()) and len(set(nlines.values())) == 1
    ok &= good
    lines.append("- DEXSeq counts: %s %s" % (", ".join("%s %d rows" % kv for kv in nlines.items()), "" if good else "<-- FAIL (all six, same rows)"))
    bl = {}
    for cond in ("recurrent", "primary"):
        p = os.path.join(SS, "rmats", "bamlist", cond + "_bamlist.txt")
        bl[cond] = read_text(p).strip().split(",") if os.path.isfile(p) else []
    exp_t = sorted(s + "_sorted.bam" for s in SAMPLES if s.split("_")[0] in RECURRENT)
    exp_c = sorted(s + "_sorted.bam" for s in SAMPLES if s.split("_")[0] in PRIMARY)
    good = sorted(bl["recurrent"]) == exp_t and sorted(bl["primary"]) == exp_c
    ok &= good
    lines.append("- rMATS bamlists: recurrent %s; primary %s %s" % (bl["recurrent"], bl["primary"], "" if good else "<-- FAIL"))
    rt = sorted(glob.glob(os.path.join(SS, "rmats", "recurrent-primary", "rmats_temp", "*.rmats")))
    seen = []
    for f in rt:
        with open(f) as fh:
            seen += [b for b in fh.readline().strip().split(",") if b]
    good = sorted(seen) == sorted(exp_t + exp_c)
    ok &= good
    lines.append("- rMATS prep: %d .rmats file(s) listing %s %s" % (len(rt), sorted(seen), "" if good else "<-- FAIL (each BAM exactly once)"))
    post = os.path.join(SS, "rmats", "recurrent-primary", "rmats_post")
    miss = [t + "." + c for t in ("SE", "MXE", "A3SS", "A5SS", "RI") for c in ("MATS.JC.txt", "MATS.JCEC.txt")
            if not os.path.isfile(os.path.join(post, t + "." + c))]
    ok &= not miss
    lines.append("- rMATS post outputs: %s" % ("all 10 MATS files" if not miss else "MISSING " + ", ".join(miss)))
    tpm = os.path.join(SS, "tximport", "suppa_tpm.txt")
    hdr = open(tpm).readline().rstrip("\n").split("\t") if os.path.isfile(tpm) else []
    good = sorted(h for h in hdr if h) == sorted(SAMPLES)
    ok &= good
    lines.append("- SUPPA TPM header: %s %s" % (hdr, "" if good else "<-- FAIL"))
    record("6", "completeness: six samples everywhere", "PASS" if ok else "FAIL", lines, {"dexseq_rows": nlines})

    # ------------------------------------------------------------------ 7 null inputs
    lines, ok = [], True
    files = {
        "DEU_COUNTS": cnt_dir,
        "DEU_GFF": os.path.join(SS, "dexseq_exon", "annotation", "DEXSeq.gff"),
        "DEU_PIPE_GENE": os.path.join(SS, "dexseq_exon", "results", "perGeneQValue.recurrent_vs_primary.csv"),
        "DMDS": os.path.join(SS, "dexseq_dtu", "filter", "drimseq", "dmDSdata.rds"),
        "DTU_PIPE_GENE": os.path.join(SS, "dexseq_dtu", "results", "dexseq", "perGeneQValue.recurrent-primary.tsv"),
        "STAGER_PIPE": os.path.join(SS, "dexseq_dtu", "results", "stager", "getAdjustedPValues.recurrent-primary.tsv"),
        "TPM": tpm,
        "IOE": os.path.join(SS, "suppa", "generate_events", "per_local_event", "events.ioe"),
        "IOI": os.path.join(SS, "suppa", "generate_events", "per_isoform", "events.ioi"),
        "PSI_L": os.path.join(SS, "suppa", "psi_per_local_event", "suppa_local.psi"),
        "PSI_I": os.path.join(SS, "suppa", "psi_per_isoform", "suppa_isoform.psi"),
        "SUPPA_PIPE_L": os.path.join(SS, "suppa", "diffsplice", "per_local_event", "recurrent-primary_local_diffsplice.dpsi"),
        "SUPPA_PIPE_I": os.path.join(SS, "suppa", "diffsplice", "per_isoform", "recurrent-primary_transcript_diffsplice.dpsi"),
        "RMATS_TEMP": os.path.join(SS, "rmats", "recurrent-primary", "rmats_temp"),
        "RMATS_PIPE_POST": post,
        "BAMLIST_T": os.path.join(SS, "rmats", "bamlist", "recurrent_bamlist.txt"),
        "BAMLIST_C": os.path.join(SS, "rmats", "bamlist", "primary_bamlist.txt"),
        "BAM_DIR": SS,
        "GTF": GTF,
        "PB": os.path.join(D, "pipeline_bin"),
    }
    for k, p in files.items():
        good = os.path.exists(p) and (os.path.isdir(p) or os.path.getsize(p) > 0)
        ok &= good
        lines.append("- %s: `%s` %s" % (k, p, "" if good else "<-- MISSING"))
    imgs = {"IMG_DEXSEQ": img.get("DEXSEQ_EXON"), "IMG_DTU": img.get("DEXSEQ_DTU"), "IMG_DRIMSEQ": img.get("DRIMSEQ_FILTER"),
            "IMG_STAGER": img.get("STAGER"), "IMG_RBASE": img.get("SPLIT_FILES_TPM"), "IMG_SUPPA": img.get("DIFFSPLICE_IOE"),
            "IMG_RMATS": img.get("RMATS_POST")}
    for k, p in imgs.items():
        good = bool(p) and os.path.isfile(p)
        ok &= good
        lines.append("- %s: `%s` %s" % (k, p, "" if good else "<-- MISSING"))
    pbok = True
    for line in read_text(os.path.join(NULL, "pipeline_bin.md5")).splitlines():
        h, f = line.split()
        p = os.path.join(files["PB"], f)
        got = subprocess.run(["md5sum", p], stdout=subprocess.PIPE, universal_newlines=True).stdout.split()[:1] if os.path.isfile(p) else []
        good = got == [h]
        pbok &= good
        lines.append("- pipeline_bin/%s md5 %s %s" % (f, got[0] if got else "missing", "= 1.0.4 tarball" if good else "<-- MISMATCH"))
    ok &= pbok
    # the exact option strings the pipeline used, re-used by the null with only inputs and outputs changed
    so = re.search(r"diffSplice (.*?) -i ", cmd_of.get("DIFFSPLICE_IOE", "") or "")
    suppa_opts = so.group(1).strip() if so else ""
    ro = re.search(r"rmats\.py (.*?) 1> ", cmd_of.get("RMATS_POST", "") or "")
    toks = ro.group(1).split() if ro else []
    keep, skip = [], {"--gtf", "--b1", "--b2", "--od", "--tmp", "--nthread", "--tstat", "--task"}
    i = 0
    while i < len(toks):
        if toks[i] in skip:
            i += 2
            continue
        keep.append(toks[i])
        i += 1
    rmats_opts = " ".join(keep)
    good = bool(suppa_opts) and " -pa" not in " " + suppa_opts and "--libType fr-firststrand" in rmats_opts \
        and "--readLength 100" in rmats_opts and "--paired-stats" not in rmats_opts
    ok &= good
    lines.append("- SUPPA diffSplice options (from DIFFSPLICE_IOE): `%s`" % suppa_opts)
    lines.append("- rMATS options (from RMATS_POST, minus inputs/outputs/threads/task): `%s` %s" % (rmats_opts, "" if good else "<-- FAIL"))
    record("7", "null inputs: files, containers, pipeline scripts, option strings", "PASS" if ok else "FAIL", lines)

    hard = [k for k, v in Cj["checks"].items() if v["status"] == "FAIL"]
    warn = [k for k, v in Cj["checks"].items() if v["status"] == "WARN"]
    Cj["hard_pass"] = not hard
    Cj["failed"] = hard
    Cj["warned"] = warn
    Cj["finished"] = now()
    md[3:3] = ["**hard_pass: %s**%s%s" % (not hard, "; failed: " + ", ".join(hard) if hard else "",
                                           "; WARN (not blocking): " + ", ".join(warn) if warn else ""), ""]
    with open(os.path.join(D, "checks.json"), "w") as fh:
        json.dump(Cj, fh, indent=1, sort_keys=True)
    with open(os.path.join(D, "CHECKS.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    envp = os.path.join(NULL, "inputs.env")
    if not hard:
        with open(envp + ".part", "w") as fh:
            fh.write("# written by check_rnasplice.py %s; sourced by null_task.sbatch\n" % now())
            for k, p in list(files.items()) + list(imgs.items()):
                fh.write("%s='%s'\n" % (k, p))
            fh.write("SUPPA_OPTS='%s'\nRMATS_OPTS='%s'\n" % (suppa_opts, rmats_opts))
        os.replace(envp + ".part", envp)
    elif os.path.exists(envp):
        os.remove(envp)
    print("hard_pass=%s failed=%s warned=%s" % (not hard, hard, warn))
    return 0 if not hard else 2


if __name__ == "__main__":
    sys.exit(main())
