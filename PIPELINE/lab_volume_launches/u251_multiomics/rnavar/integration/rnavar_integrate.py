#!/usr/bin/env python3
"""rnavar_integrate.py: variant integration for the U251N LITT rnavar arm.

nf-core/rnavar 1.3.0 ran on the graft-only U251 reads (xengsort graft bin, R2 clipped 3 nt upstream). This script
turns its per-sample calls into two answers: is the model U251 (authentication), and does the clonal make-up move
after LITT (per-site VAF shifts and arm-level allelic imbalance, C2B as the parental baseline, every group statistic
ranked among the 10 balanced relabellings of the six tumours).

Standard library + pysam only (runs in the pinned pysam 0.24.1 container). Subcommands, in the order
integrate.sbatch runs them:
  checks     pipeline checks; the hard ones decide whether anything else runs
  sites      union of PASS sites, VEP annotation, authentication loci from the GENCODE GTF, mpileup regions
  counts     parse the forced bcftools mpileup counts (new graft-only BAMs and old graft+both BAMs)
  select     categories under all 10 splits -> candidate sites for the read-level pass; liftover BED
  readlevel  one BAM: per-read evidence at the candidate sites; alt-read pairs for the rat realignment
  yx         one BAM: fragments on X-degenerate Y genes and their X gametologs
  final      artefact filters, categories, split ranks, arm imbalance, Y/X, authentication, REPORT.md

VAF comes only from the forced mpileup counts, never from HaplotypeCaller's (downsampled) AD.
"""
import argparse
import collections
import glob
import gzip
import itertools
import json
import math
import os
import random
import re
import statistics
import sys

import pysam

# ---------------------------------------------------------------------------------------------------------- design
LIBS = ["C2B", "IL64B", "IL66B", "IL67B", "IL68B", "IL69B", "N168B", "N269B", "NL70B", "NL71B"]
PRIMARY = ["IL67B", "IL68B", "IL69B"]          # not ablated
RECURRENT = ["IL66B", "NL70B", "NL71B"]        # LITT, regrowth
TUMOURS = PRIMARY + RECURRENT
CORE = ["C2B"] + TUMOURS                       # the seven libraries every category is defined on
SMALL = ["IL64B", "N168B"]                     # 0.31 M / 0.44 M graft pairs: descriptive only
LARGE = [x for x in LIBS if x not in SMALL]    # >= 3.8 M graft pairs
RECIPIENTS = ["N269B", "N168B", "IL64B"]       # libraries whose human reads may be carry-over or index hopping
MIN_READS = 20
PASS, FAIL, WARN, INFO = "PASS", "FAIL", "WARN", "INFO"
SEED = 20260927
MAX_CAND = 15000        # read-level pass: at most this many sites (priority order, see select)
MAX_READS_SITE = 3000   # reads examined per site and BAM (pileup max_depth)
MAX_ALT_PAIRS_SITE = 30 # alt-supporting pairs per site and BAM sent to the rat realignment
MAX_ALT_PAIRS_BAM = 60000
GENETIC_CODE = dict(zip((a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"),
                        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))  # NCBI table 1
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
YX_PAIRS = [("RPS4Y1", "RPS4X"), ("DDX3Y", "DDX3X"), ("KDM5D", "KDM5C"), ("UTY", "KDM6A"), ("EIF1AY", "EIF1AX"),
            ("ZFY", "ZFX"), ("USP9Y", "USP9X")]
AUTH = {  # gene -> (codon, expected reference amino acid). Cellosaurus CVCL_0021 (v52).
    "TP53": (273, "R"),   # p.Arg273His, c.818G>A, homozygous
    "IDH1": (132, "R"),   # none reported; engineered R132C/R132H children exist, so R132 must read reference
    "PTEN": (242, "E"),   # p.Glu242ValfsTer15, c.723_724dupTG
}


def splits():
    """The 10 balanced 3-v-3 splits of the six tumours, IL67B always in A so no split appears as its mirror.
    Split 1 is the true labelling (A = primary, B = recurrent)."""
    out = []
    for pair in itertools.combinations(["IL68B", "IL69B", "IL66B", "NL70B", "NL71B"], 2):
        a = ["IL67B"] + list(pair)
        out.append((a, [t for t in TUMOURS if t not in a]))
    assert set(out[0][0]) == set(PRIMARY) and set(out[0][1]) == set(RECURRENT), "split 1 must be the true labels"
    assert len(out) == 10
    return out


# ------------------------------------------------------------------------------------------------------ paths
def p_new_bam(r, lib): return f"{r}/preprocessing/{lib}_graft/{lib}_graft.md.bam"
def p_vcf(r, lib): return f"{r}/variant_calling/{lib}_graft/{lib}_graft.haplotypecaller.filtered.vcf.gz"
def p_vep(r, lib): return f"{r}/annotation/{lib}_graft/{lib}_graft.haplotypecaller.filtered_VEP.ann.vcf.gz"
def p_star_log(r, lib): return f"{r}/reports/star/{lib}_graft/{lib}_graft.Log.final.out"
def p_flagstat(r, lib): return f"{r}/reports/samtools/{lib}_graft/{lib}_graft.md.flagstat"
def p_seq2hla(r, lib, kind="ClassI-class.HLAgenotype4digits"): return f"{r}/seq2hla/{lib}_graft/{lib}_graft-{kind}"
def p_old_bam(u, lib): return f"{u}/ANALYSIS/results_human_final/star_salmon/{lib}.markdup.sorted.bam"


# ------------------------------------------------------------------------------------------------------ helpers
def log(*a):
    print(*a, file=sys.stderr, flush=True)


def write_tsv(path, header, rows):
    tmp = path + ".part"
    op = gzip.open if path.endswith(".gz") else open
    with op(tmp, "wt") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(fmt(x) for x in r) + "\n")
    os.replace(tmp, path)


def fmt(x):
    if x is None:
        return "NA"
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.4g}"
    return str(x)


def read_tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def num(x, cast=int):
    try:
        return cast(x)
    except (TypeError, ValueError):
        return None


def wilson(k, n, z=1.959963984540054):
    """Wilson score 95 % interval for a binomial proportion."""
    if not n:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, c - h), min(1.0, c + h))


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sxx * syy)


def spearman(x, y):
    return pearson(ranks(x), ranks(y))


def fisher_two_sided(a, b, c, d):
    """Two-sided Fisher exact test on [[a, b], [c, d]] (sum of tables no more likely than the observed one)."""
    r1, r2, c1 = a + b, c + d, a + c
    n = r1 + r2
    lf = lambda k: math.lgamma(k + 1)
    def lp(x):
        return lf(r1) + lf(r2) + lf(c1) + lf(n - c1) - lf(n) - lf(x) - lf(r1 - x) - lf(c1 - x) - lf(r2 - c1 + x)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    obs = lp(a)
    tot = 0.0
    for x in range(lo, hi + 1):
        v = lp(x)
        if v <= obs + 1e-7:
            tot += math.exp(v)
    return min(1.0, tot)


def logit_lrt(ks, ns, xs, iters=50):
    """Binomial logistic regression logit(p) = a + b x (IRLS, small ridge on b against separation).
    Returns (b, p) with p from the likelihood-ratio test against the intercept-only model (chi-square, 1 df)."""
    K, N = sum(ks), sum(ns)
    if N == 0 or K == 0 or K == N:
        return (0.0, 1.0)
    def ll(a, b):
        s = 0.0
        for k, n, x in zip(ks, ns, xs):
            eta = a + b * x
            # log p = -log(1 + e^-eta), log(1-p) = -log(1 + e^eta), computed stably
            lp = -math.log1p(math.exp(-eta)) if eta > -30 else eta
            lq = -math.log1p(math.exp(eta)) if eta < 30 else -eta
            s += k * lp + (n - k) * lq
        return s
    p0 = K / N
    a0 = math.log(p0 / (1 - p0))
    a, b = a0, 0.0
    for _ in range(iters):
        g0 = g1 = h00 = h01 = h11 = 0.0
        for k, n, x in zip(ks, ns, xs):
            p = 1 / (1 + math.exp(-max(-30, min(30, a + b * x))))
            w = n * p * (1 - p)
            r = k - n * p
            g0 += r; g1 += r * x
            h00 += w; h01 += w * x; h11 += w * x * x
        h11 += 1e-6; g1 -= 1e-6 * b
        det = h00 * h11 - h01 * h01
        if det <= 0:
            break
        da = (h11 * g0 - h01 * g1) / det; db = (h00 * g1 - h01 * g0) / det
        a += da; b += db
        if abs(da) < 1e-8 and abs(db) < 1e-8:
            break
    stat = max(0.0, 2 * (ll(a, b) - ll(a0, 0.0)))
    return (b, math.erfc(math.sqrt(stat / 2)))


def median(v):
    return statistics.median(v) if v else float("nan")


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def fai_order(fasta):
    order = {}
    with open(fasta + ".fai") as f:
        for i, line in enumerate(f):
            p = line.split("\t")
            order[p[0]] = (i, int(p[1]))
    return order


def normalize(fa, chrom, pos, ref, alt):
    """Minimal, left-aligned representation of one REF/ALT pair (Tan et al. 2015), 1-based pos.
    Returns (chrom, pos, ref, alt) or None for a non-variant or symbolic allele."""
    ref, alt = ref.upper(), alt.upper()
    if not ref or not alt or ref == alt or alt.startswith("<") or alt in ("*", "."):
        return None
    while True:
        changed = False
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
            changed = True
        if (not ref or not alt) and pos > 1:
            pos -= 1
            b = fa.fetch(chrom, pos - 1, pos).upper()
            ref, alt = b + ref, b + alt
            changed = True
        if not changed:
            break
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return (chrom, pos, ref, alt)


def kind_of(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(ref) == 1 and len(alt) > 1 and alt[0] == ref[0]:
        return "ins"
    if len(alt) == 1 and len(ref) > 1 and ref[0] == alt[0]:
        return "del"
    return "complex"


def load_prep(prep_dir, lib):
    path = f"{prep_dir}/{lib}.tsv"
    if not os.path.exists(path):
        return {}
    rows = list(read_tsv(path))
    return rows[0] if rows else {}


def load_xengsort(u, lib):
    """host, graft, ambiguous, both, neither from xengsort's classify log (the block headed 'prefix')."""
    path = f"{u}/ANALYSIS/xengsort_out/{lib}.txt"
    with open(path) as f:
        lines = [x.rstrip("\n") for x in f]
    for i, x in enumerate(lines):
        p = x.split("\t")
        if p[:6] == ["prefix", "host", "graft", "ambiguous", "both", "neither"] and i + 1 < len(lines):
            q = lines[i + 1].split("\t")
            if q[0] == lib:
                return dict(zip(["host", "graft", "ambiguous", "both", "neither"], map(int, q[1:6])))
    raise ValueError(f"no counts block for {lib} in {path}")


# ------------------------------------------------------------------------------------------------------ checks
def parse_star_log(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            if "|" in line:
                k, v = line.split("|", 1)
                out[k.strip()] = v.strip()
    return out


def parse_flagstat(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            m = re.match(r"^(\d+) \+ (\d+) (.+?)(?: \(.*)?$", line.strip())
            if m:
                out[m.group(3)] = int(m.group(1)) + int(m.group(2))
    return out


def csq_fields(vf):
    if "CSQ" not in vf.header.info:
        return None
    d = vf.header.info["CSQ"].description
    return d.split("Format: ", 1)[1].strip().strip('"').split("|") if "Format: " in d else None


def vep_cmdline(vf):
    for r in vf.header.records:
        if r.key == "VEP-command-line":
            return str(r).strip()
    return ""


def csq_af(e):
    v = []
    for k in ("gnomADe_AF", "gnomADg_AF"):
        x = num(e.get(k, ""), float)
        if x is not None:
            v.append(x)
    return max(v) if v else None


def r2_profile(bam_path, fasta, n_per_chrom=10000, chroms=None):
    """R2 reads (primary, MAPQ >= 60, not duplicate): per sequencing cycle 1-10, the fraction soft-clipped at the 5'
    end and the mismatch rate of aligned bases against the reference."""
    fa = pysam.FastaFile(fasta)
    bam = pysam.AlignmentFile(bam_path)
    clip = [0] * 10; mm = [0] * 10; al = [0] * 10; nreads = 0
    names = chroms or [c for c in bam.references if re.match(r"^chr([1-9]|10)$", c)]
    for c in names:
        k = 0
        for r in bam.fetch(c):
            if r.flag & 0xF04 or not r.is_read2 or r.mapping_quality < 60 or r.query_sequence is None:
                continue
            k += 1
            if k > n_per_chrom:
                break
            nreads += 1
            cig = r.cigartuples
            qlen = r.query_length
            c5 = cig[-1][1] if (r.is_reverse and cig[-1][0] == 4) else (cig[0][1] if (not r.is_reverse and cig[0][0] == 4) else 0)
            pairs = dict(r.get_aligned_pairs(matches_only=True))
            need = []
            for cyc in range(10):
                if cyc < c5:
                    clip[cyc] += 1
                    continue
                qi = (qlen - 1 - cyc) if r.is_reverse else cyc
                rp = pairs.get(qi)
                if rp is not None:
                    need.append((cyc, qi, rp))
            if not need:
                continue
            lo = min(x[2] for x in need); hi = max(x[2] for x in need)
            if hi - lo < 2000:
                seg = fa.fetch(c, lo, hi + 1).upper()
                refb = {rp: seg[rp - lo] for _, _, rp in need}
            else:
                refb = {rp: fa.fetch(c, rp, rp + 1).upper() for _, _, rp in need}
            q = r.query_sequence
            for cyc, qi, rp in need:
                al[cyc] += 1
                if q[qi] != refb[rp]:
                    mm[cyc] += 1
    rate = [(clip[i] + mm[i]) / (clip[i] + al[i]) if (clip[i] + al[i]) else float("nan") for i in range(10)]
    return {"reads": nreads, "clip": [clip[i] / nreads if nreads else float("nan") for i in range(10)],
            "mismatch": [mm[i] / al[i] if al[i] else float("nan") for i in range(10)], "rate": rate}


def parse_samtools_mpc(path):
    """MPC (mismatches per cycle, needs -r) over bases per cycle from FFQ+LFQ; the input is R2-only."""
    mpc, bases = {}, collections.Counter()
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] == "MPC":
                mpc[int(p[1])] = sum(int(x) for x in p[3:])
            elif p[0] in ("FFQ", "LFQ"):
                bases[int(p[1])] += sum(int(x) for x in p[2:])
    return [mpc.get(c, 0) / bases[c] if bases.get(c) else float("nan") for c in range(1, 11)]


def cmd_checks(a):
    rows = []
    def row(check, lib, expected, observed, result, outside=""):
        rows.append([check, lib, expected, observed, result, outside])
        log(f"  [{result}] {check} {lib}: {observed}")
    R = a.results
    nfx = open(a.nextflow_exit).read().strip() if os.path.exists(a.nextflow_exit) else "missing"
    row("nextflow_exit", "-", "0", nfx, PASS if nfx == "0" else FAIL)

    # arm labels hard-coded here against the metadata the published DE used
    meta = {}
    if os.path.exists(a.metadata):
        for r in read_tsv_csv(a.metadata):
            meta[r["sample"]] = r["Classification"]
    ok = (sorted(k for k, v in meta.items() if v.startswith("Primary")) == sorted(PRIMARY) and
          sorted(k for k, v in meta.items() if v.startswith("Recurrent")) == sorted(RECURRENT))
    row("arm_labels", "-", f"primary {PRIMARY}, recurrent {RECURRENT}", json.dumps(meta), PASS if ok else FAIL,
        "U/ANALYSIS/metadata_therapy.csv (published DE design)")

    prep = {lib: load_prep(a.prep_checks, lib) for lib in LIBS}
    for lib in LIBS:
        st = prep[lib].get("status", "missing")
        row("prep_reads_status", lib, "OK", st, PASS if st == "OK" else FAIL, "prep_reads checks a-g")

    for lib in LIBS:
        sev = FAIL if lib in LARGE else WARN
        bam = p_new_bam(R, lib)
        ok = os.path.exists(bam) and os.path.getsize(bam) > 0 and os.path.exists(bam + ".bai")
        row("md_bam_present", lib, "md.bam + .bai", bam if ok else "missing", PASS if ok else sev)
        if ok:
            h = pysam.AlignmentFile(bam).header.to_dict()
            sms = sorted({g.get("SM") for g in h.get("RG", [])})
            row("bam_read_group_SM", lib, f"['{lib}_graft']", sms, PASS if sms == [f"{lib}_graft"] else sev,
                "samplesheet sample name")
        slog = parse_star_log(p_star_log(R, lib))
        n_in = num(slog.get("Number of input reads"))
        n_graft = num(prep[lib].get("n_graft"))
        row("star_input_pairs_eq_graft_pairs", lib, n_graft, n_in,
            PASS if (n_in is not None and n_in == n_graft) else sev, "prep_reads graft record count (xengsort graft bin)")
        uniq = num(slog.get("Uniquely mapped reads number")); multi = num(slog.get("Number of reads mapped to multiple loci"))
        fs = parse_flagstat(p_flagstat(R, lib))
        if uniq is not None and multi is not None and fs.get("primary") is not None:
            exp = 2 * (uniq + multi)
            dev = abs(fs["primary"] - exp) / exp if exp else float("inf")
            row("flagstat_primary_vs_2x_star_mapped", lib, exp, f"{fs['primary']} ({dev:.4%})",
                PASS if dev <= 0.01 else sev, "STAR Log.final.out (unique + multi-mapped pairs)")
            row("flagstat_total_vs_2x_star_input", lib, 2 * n_in if n_in else "NA", fs.get("in total"), INFO,
                "the task's literal form; STAR writes no unmapped reads and adds secondary records, so equality is not expected")
        else:
            row("flagstat_primary_vs_2x_star_mapped", lib, "STAR log + flagstat", "missing", sev)
        vcf = p_vcf(R, lib)
        n_rec = sum(1 for _ in pysam.VariantFile(vcf)) if os.path.exists(vcf) else None
        row("vcf_records", lib, "> 0" if lib in LARGE else "descriptive", n_rec,
            PASS if (n_rec or 0) > 0 else (FAIL if lib in LARGE else INFO))
        vep = p_vep(R, lib)
        if os.path.exists(vep):
            vf = pysam.VariantFile(vep)
            cmd = vep_cmdline(vf)
            fields = csq_fields(vf) or []
            n_vep = 0; n_common = 0
            for rec in vf:
                n_vep += 1
                afs = [csq_af(dict(zip(fields, e.split("|")))) for e in (rec.info.get("CSQ") or ())]
                afs = [x for x in afs if x is not None]
                if afs and max(afs) > 0.01:
                    n_common += 1
            row("vep_command_line", lib, "--allele_number, no --filter_common, no --per_gene", cmd[:400],
                PASS if ("allele_number" in cmd and "filter_common" not in cmd and "per_gene" not in cmd) else FAIL)
            row("vep_records_eq_vcf_records", lib, n_rec, n_vep, PASS if n_vep == n_rec else WARN)
            row("vep_rows_gnomAD_AF_gt_0.01", lib, "> 0 (proves --filter_common is gone)", n_common,
                PASS if n_common > 0 else (FAIL if lib in LARGE else INFO), "gnomAD allele frequencies in the VEP 115 cache")
            row("vep_ALLELE_NUM_field", lib, "present", "ALLELE_NUM" in fields, PASS if "ALLELE_NUM" in fields else FAIL)
        else:
            row("vep_vcf_present", lib, vep, "missing", FAIL if lib in LARGE else INFO)
        s4 = p_seq2hla(R, lib)
        row("seq2hla_classI_4digit", lib, "present" if lib in LARGE else "ignored (too few reads)",
            "present" if os.path.exists(s4) else "missing", PASS if os.path.exists(s4) else (WARN if lib in LARGE else INFO))

    # trace: no cache download, one STAR_GENOMEGENERATE, every task finished
    traces = sorted(glob.glob(f"{R}/pipeline_info/execution_trace_*.txt"))
    best = {}          # (process, tag) -> best status
    hashes = collections.defaultdict(set)
    names_hash = {}
    rank = {"COMPLETED": 3, "CACHED": 3, "IGNORED": 1, "FAILED": 0, "ABORTED": 0}
    for t in traces:
        for r in read_tsv(t):
            nm = r.get("name", ""); st = r.get("status", "")
            m = re.match(r"^(\S+)(?: \((.*)\))?$", nm)
            proc = m.group(1).split(":")[-1] if m else nm
            tag = (m.group(2) or "") if m else ""
            key = (proc, tag)
            if rank.get(st, -1) >= rank.get(best.get(key, ""), -1):
                best[key] = st
            if st in ("COMPLETED", "CACHED"):
                hashes[proc].add(r.get("hash"))
                names_hash[(proc, tag)] = r.get("hash")
    row("trace_files", "-", ">= 1", len(traces), PASS if traces else FAIL)
    dl = sorted(p for p in hashes if "DOWNLOAD" in p)
    row("trace_no_cache_download", "-", "no *DOWNLOAD task", dl or "none", PASS if not dl else FAIL)
    ng = len(hashes.get("STAR_GENOMEGENERATE", ()))
    row("trace_one_star_genomegenerate", "-", 1, ng, PASS if ng == 1 else FAIL)
    bad = sorted(f"{p} ({t}): {s}" for (p, t), s in best.items() if s not in ("COMPLETED", "CACHED"))
    row("trace_unfinished_tasks", "-", "none", bad[:40] or "none", PASS if not bad else WARN)
    hc = collections.Counter()
    for (p, t), s in best.items():
        if p == "GATK4_HAPLOTYPECALLER" and s in ("COMPLETED", "CACHED"):
            m = re.match(r"^(\w+?)_graft_", t)
            if m:
                hc[m.group(1)] += 1
    nint = max(hc.values()) if hc else 0
    for lib in LIBS:
        row("haplotypecaller_intervals_complete", lib, nint, hc.get(lib, 0),
            PASS if (hc.get(lib, 0) == nint and nint > 0) else (FAIL if lib in LARGE else WARN))

    # the overrides are in the commands that ran (work directories, when they still exist)
    want = {
        "STAR_ALIGN": (["--sjdbOverhang 100", "--twopassMode Basic", "PL:illumina"], []),
        "STAR_GENOMEGENERATE": (["--sjdbOverhang 100", "--sjdbGTFfile"], []),
        "ENSEMBLVEP_VEP": (["--allele_number", "--everything", "--offline", "--fasta", "--cache_version 115"],
                           ["--filter_common", "--per_gene"]),
        "GATK4_HAPLOTYPECALLER": (["--dont-use-soft-clipped-bases", "--dbsnp"], []),
        "GATK4_BASERECALIBRATOR": (["dbsnp138.primary.vcf.gz", "mills_1000g.primary.vcf.gz"], []),
    }
    for proc, (must, mustnot) in want.items():
        hs = [h for (p, t), h in names_hash.items() if p == proc]
        if not hs:
            row(f"command_sh:{proc}", "-", "task in trace", "none", WARN)
            continue
        wd = glob.glob(f"{a.work}/{hs[0][:2]}/{hs[0][3:]}*/.command.sh")
        if not wd:
            row(f"command_sh:{proc}", "-", "work dir", "not available (work deleted?)", INFO)
            continue
        txt = open(wd[0]).read()
        miss = [x for x in must if x not in txt]; extra = [x for x in mustnot if x in txt]
        row(f"command_sh:{proc}", "-", f"has {must}, lacks {mustnot}", f"missing {miss}, present {extra}",
            PASS if not miss and not extra else FAIL, wd[0])

    # the saved STAR index was built from this FASTA with sjdbOverhang 100
    gp = f"{R}/genome/index/star/genomeParameters.txt"
    if os.path.exists(gp):
        txt = open(gp).read()
        m = re.search(r"^sjdbOverhang\s+(\d+)", txt, re.M)
        row("star_index_sjdbOverhang", "-", "100", m.group(1) if m else "absent", PASS if m and m.group(1) == "100" else FAIL)
        cnl = f"{R}/genome/index/star/chrNameLength.txt"
        if os.path.exists(cnl):
            idx = [tuple(x.split()[:2]) for x in open(cnl) if x.strip()]
            fai = [tuple(x.split("\t")[:2]) for x in open(a.fasta + ".fai") if x.strip()]
            row("star_index_contigs_eq_fasta_fai", "-", f"{len(fai)} contigs", f"{len(idx)} contigs",
                PASS if idx == fai else FAIL, "GENCODE v44 fai (stage_refs)")
    else:
        row("star_index_saved", "-", gp, "missing", WARN)

    # Takara 3-nt R2 artefact: clipped once upstream, so R2 cycles 1-3 must not stand out in the new BAMs.
    # Positive control: the old graft+both BAMs were never clipped.
    prof = {}
    for lib in LARGE:
        for tag, path in ((f"{lib}.new", p_new_bam(R, lib)), (f"{lib}.old", p_old_bam(a.u, lib))):
            if not os.path.exists(path):
                continue
            pr = r2_profile(path, a.fasta, n_per_chrom=a.r2_reads)
            prof[tag] = pr
            e = mean(pr["rate"][0:3]); b = mean(pr["rate"][3:10])
            ratio = e / b if b else float("nan")
            obs = f"reads {pr['reads']}; clip+mismatch rate c1-3 {e:.4f} vs c4-10 {b:.4f} (x{ratio:.2f})"
            if tag.endswith(".new"):
                row("r2_5prime_profile", tag, "c1-3 <= 1.5 x c4-10 + 0.005", obs, PASS if e <= 1.5 * b + 0.005 else FAIL,
                    "old unclipped BAM of the same library (positive control)")
            else:
                row("r2_5prime_profile_positive_control", tag, "elevated (unclipped)", obs, INFO)
            mpc = parse_samtools_mpc(f"{a.mpc_dir}/{tag}.r2.stats")
            if mpc is not None:
                e2 = mean([x for x in mpc[0:3] if not math.isnan(x)]); b2 = mean([x for x in mpc[3:10] if not math.isnan(x)])
                row("samtools_stats_MPC_r2", tag, "c1-3 not elevated vs c4-10" if tag.endswith(".new") else "positive control",
                    f"mismatch/base c1-3 {e2:.4f} vs c4-10 {b2:.4f}",
                    (PASS if e2 <= 1.5 * b2 + 0.005 else WARN) if tag.endswith(".new") else INFO,
                    "samtools stats -r on 300 k R2 records of chr1")
    power = [t for t in prof if t.endswith(".old") and mean(prof[t]["rate"][0:3]) > 1.5 * mean(prof[t]["rate"][3:10])]
    row("r2_profile_has_power", "-", "old BAMs show the artefact", f"{len(power)} of {sum(t.endswith('.old') for t in prof)}",
        INFO if power else WARN, "if no old BAM shows it, the new-BAM PASS says nothing")
    with open(os.path.join(a.out_dir, "r2_profiles.json"), "w") as f:
        json.dump(prof, f, indent=1)

    write_tsv(os.path.join(a.out_dir, "CHECKS.tsv"), ["check", "library", "expected", "observed", "result", "outside_reference"],
              rows)
    nfail = sum(1 for r in rows if r[4] == FAIL)
    log(f"checks: {nfail} FAIL, {sum(1 for r in rows if r[4] == WARN)} WARN")
    return 0 if nfail == 0 else 2


def read_tsv_csv(path):
    with open(path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            if line.strip():
                yield dict(zip(header, line.strip().split(",")))


# ------------------------------------------------------------------------------------------------------ sites
def parse_attrs(s):
    return dict(re.findall(r'(\S+) "([^"]*)"', s))


def cmd_sites(a):
    fa = pysam.FastaFile(a.fasta)
    order = fai_order(a.fasta)
    os.makedirs(a.work, exist_ok=True)

    # ---- one pass over the GTF: gene table, MANE CDS of the authentication genes, exons of the Y/X genes, PAR genes
    yx_names = {n for p in YX_PAIRS for n in p}
    genes, cds, exons, tags = [], collections.defaultdict(list), collections.defaultdict(list), {}
    op = gzip.open if a.gtf.endswith(".gz") else open
    with op(a.gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            feat = p[2]
            if feat == "gene":
                at = parse_attrs(p[8])
                genes.append((at.get("gene_id"), at.get("gene_name"), at.get("gene_type"), p[0], int(p[3]), int(p[4]), p[6]))
            elif feat == "CDS" and any(f'gene_name "{g}";' in p[8] for g in AUTH):
                at = parse_attrs(p[8])
                if 'tag "MANE_Select"' in p[8]:
                    cds[(at["gene_name"], at["transcript_id"])].append((p[0], int(p[3]), int(p[4]), p[6], p[7], at.get("protein_id", "")))
            elif feat == "exon" and any(f'gene_name "{g}";' in p[8] for g in yx_names):
                at = parse_attrs(p[8])
                if at.get("gene_name") in yx_names and not at.get("gene_id", "").endswith("_PAR_Y"):
                    exons[at["gene_name"]].append((p[0], int(p[3]), int(p[4])))
    write_tsv(f"{a.work}/gtf_genes.tsv.gz", ["gene_id", "gene_name", "gene_type", "chrom", "start", "end", "strand"], genes)
    yx_rows = []
    for g, ex in sorted(exons.items()):
        ex.sort()
        merged = []
        for c, s, e in ex:
            if merged and merged[-1][0] == c and s <= merged[-1][2] + 1:
                merged[-1][2] = max(merged[-1][2], e)
            else:
                merged.append([c, s, e])
        yx_rows += [[g] + m for m in merged]
    write_tsv(f"{a.work}/yx_exons.tsv", ["gene_name", "chrom", "start", "end"], yx_rows)

    # ---- authentication loci from the MANE Select CDS, each checked against the known reference residue
    auth = {"loci": {}, "problems": []}
    for gene, (codon, aa) in AUTH.items():
        tx = [k for k in cds if k[0] == gene]
        if len(tx) != 1:
            auth["problems"].append(f"{gene}: {len(tx)} MANE Select transcripts with CDS in the GTF")
            continue
        parts = cds[tx[0]]
        chrom, strand = parts[0][0], parts[0][3]
        segs = sorted((s, e) for _, s, e, _, _, _ in parts)
        first = min(parts, key=lambda x: x[1]) if strand == "+" else max(parts, key=lambda x: x[2])
        if first[4] not in ("0", "."):
            auth["problems"].append(f"{gene}: first CDS phase {first[4]} (incomplete CDS)")
        gpos = []
        if strand == "+":
            for s, e in segs:
                gpos.extend(range(s, e + 1))
        else:
            for s, e in reversed(segs):
                gpos.extend(range(e, s - 1, -1))
        def tbase(i):  # transcript-strand base at CDS position i (1-based)
            b = fa.fetch(chrom, gpos[i - 1] - 1, gpos[i - 1]).upper()
            return b if strand == "+" else b.translate(COMP)
        cod = "".join(tbase(i) for i in range(3 * codon - 2, 3 * codon + 1))
        got = GENETIC_CODE.get(cod, "?")
        rec = {"gene": gene, "transcript": tx[0][1], "protein_id": parts[0][5], "chrom": chrom, "strand": strand,
               "codon": codon, "codon_seq": cod, "aa": got, "aa_expected": aa,
               "codon_gpos": [gpos[i - 1] for i in range(3 * codon - 2, 3 * codon + 1)]}
        if got != aa:
            auth["problems"].append(f"{gene} codon {codon}: reference codon {cod} = {got}, expected {aa} (mapping error)")
        if gene == "TP53":
            if not parts[0][5].startswith("ENSP00000269305"):
                auth["problems"].append(f"TP53 MANE protein {parts[0][5]} is not ENSP00000269305")
            if tbase(818) != "G":
                auth["problems"].append(f"TP53 c.818 is {tbase(818)}, not G")
            g = gpos[817]
            gref = fa.fetch(chrom, g - 1, g).upper()
            galt = "A" if strand == "+" else "T"          # c.818G>A on the transcript strand
            rec["expected_key"] = list(normalize(fa, chrom, g, gref, galt))
        if gene == "PTEN":
            dup = tbase(723) + tbase(724)
            rec["c723_724"] = dup
            if dup != "TG":
                auth["problems"].append(f"PTEN c.723_724 is {dup}, not TG")
            lo, hi = sorted((gpos[722], gpos[723]))
            gseq = fa.fetch(chrom, lo - 1, hi).upper()     # genomic-strand copy of the duplicated segment
            anchor = fa.fetch(chrom, hi - 1, hi).upper()
            rec["expected_key"] = list(normalize(fa, chrom, hi, anchor, anchor + gseq))
            rec["window"] = [lo - 12, hi + 12]
        auth["loci"][gene] = rec
    with open(f"{a.work}/auth.json", "w") as f:
        json.dump(auth, f, indent=1)
    for p in auth["problems"]:
        log("AUTH PROBLEM:", p)

    # ---- union of PASS sites over the ten per-sample calls, VEP annotation tied to its ALT by ALLELE_NUM
    sites = {}
    stats = collections.Counter()
    for lib in LIBS:
        path = p_vep(a.results, lib)
        use_vep = os.path.exists(path)
        if not use_vep:
            path = p_vcf(a.results, lib)
        if not os.path.exists(path):
            log(f"no VCF for {lib}")
            continue
        vf = pysam.VariantFile(path)
        fields = csq_fields(vf) if use_vep else None
        for rec in vf:
            filt = ";".join(rec.filter.keys()) or "."
            smp = rec.samples[0]
            gt = "/".join("." if x is None else str(x) for x in (smp.get("GT") or ()))
            ad = smp.get("AD")
            entries = []
            if fields:
                for e in rec.info.get("CSQ") or ():
                    entries.append(dict(zip(fields, e.split("|"))))
            for i, alt in enumerate(rec.alts or (), start=1):
                key = normalize(fa, rec.chrom, rec.pos, rec.ref, alt)
                if key is None:
                    continue
                s = sites.get(key)
                if s is None:
                    s = sites[key] = {"pass": set(), "filters": {}, "gt": {}, "hc_ad": {}, "csq": None}
                s["filters"][lib] = filt
                s["gt"][lib] = gt
                if ad and len(ad) > i and ad[0] is not None and ad[i] is not None:
                    s["hc_ad"][lib] = f"{ad[0]},{ad[i]}"
                if filt == "PASS":
                    s["pass"].add(lib)
                if entries and s["csq"] is None:
                    mine = [e for e in entries if e.get("ALLELE_NUM") == str(i)]
                    if mine:
                        s["csq"] = mine
                stats["records_alleles"] += 1
    union = {k: v for k, v in sites.items() if v["pass"]}
    log(f"sites: {len(sites)} alleles seen, {len(union)} PASS in >= 1 library")

    def pick(entries):
        if not entries:
            return {}
        for cond in (lambda e: e.get("CANONICAL") == "YES" and e.get("MANE_SELECT"), lambda e: e.get("CANONICAL") == "YES",
                     lambda e: True):
            c = [e for e in entries if cond(e)]
            if c:
                return c[0]
        return {}

    # authentication keys: from VEP (the called allele) and from the GTF (expected), cross-checked in final
    vep_auth = {}
    for key, s in sites.items():
        for e in s["csq"] or ():
            ensp = e.get("ENSP", "").split(".")[0]; hgvsp = e.get("HGVSp", ""); hgvsc = e.get("HGVSc", "")
            sym = e.get("SYMBOL", "")
            if ensp == "ENSP00000269305" and hgvsp.endswith("p.Arg273His"):
                vep_auth.setdefault("TP53_R273H", list(key))
            if sym == "PTEN" and (hgvsp.endswith("p.Glu242ValfsTer15") or re.search(r":c\.723_724dup(TG)?$", hgvsc)):
                vep_auth.setdefault("PTEN_c723_724dupTG", list(key))
            if sym == "IDH1" and re.search(r"p\.Arg132", hgvsp):
                vep_auth.setdefault("IDH1_R132_called", list(key))
    auth["vep_keys"] = vep_auth
    with open(f"{a.work}/auth.json", "w") as f:
        json.dump(auth, f, indent=1)

    extra = {}
    for gene in ("TP53", "PTEN"):
        k = auth["loci"].get(gene, {}).get("expected_key")
        if k:
            extra[tuple(k)] = f"{gene}_expected"
    for lab, k in vep_auth.items():
        extra.setdefault(tuple(k), lab)

    keys = sorted(set(union) | set(extra), key=lambda k: (order.get(k[0], (1e9,))[0], k[1], k[2], k[3]))
    rows = []
    for i, key in enumerate(keys, start=1):
        s = sites.get(key, {"pass": set(), "filters": {}, "gt": {}, "hc_ad": {}, "csq": None})
        e = pick(s["csq"])
        afs = [csq_af(x) for x in (s["csq"] or ())]
        afs = [x for x in afs if x is not None]
        rows.append([f"S{i:07d}", key[0], key[1], key[2], key[3], kind_of(key[2], key[3]), int(key in union),
                     extra.get(key, ""), len(s["pass"]), ",".join(sorted(s["pass"])),
                     ";".join(f"{k}:{v}" for k, v in sorted(s["filters"].items())),
                     ";".join(f"{k}:{v}" for k, v in sorted(s["gt"].items())),
                     ";".join(f"{k}:{v}" for k, v in sorted(s["hc_ad"].items())),
                     max(afs) if afs else None, e.get("SYMBOL", ""), e.get("Gene", ""), e.get("Consequence", ""),
                     e.get("BIOTYPE", ""), e.get("HGVSc", ""), e.get("HGVSp", ""), e.get("Existing_variation", ""),
                     e.get("STRAND", "")])
    write_tsv(f"{a.work}/sites.tsv.gz",
              ["site_id", "chrom", "pos", "ref", "alt", "kind", "in_union", "auth", "n_pass", "pass_libs", "filters", "gt",
               "hc_ad", "gnomad_af", "symbol", "gene_id", "consequence", "biotype", "hgvsc", "hgvsp", "existing",
               "vep_strand"], rows)

    # ---- mpileup regions: SNV = the base; indel = anchor .. anchor + allele length + 10 (aligners place gaps
    # anywhere in a repeat; normalisation brings them back to the same key); authentication codons and PTEN window
    reg = []
    for r in rows:
        c, p, ref, alt = r[1], r[2], r[3], r[4]
        reg.append((c, p, p if r[5] == "snv" else p + max(len(ref), len(alt)) + 10))
    for gene, rec in auth["loci"].items():
        for g in rec.get("codon_gpos", []):
            reg.append((rec["chrom"], g, g))
        if "window" in rec:
            reg.append((rec["chrom"], rec["window"][0], rec["window"][1]))
    reg.sort(key=lambda x: (order.get(x[0], (1e9,))[0], x[1], x[2]))
    merged = []
    for c, s, e in reg:
        if merged and merged[-1][0] == c and s <= merged[-1][2] + 1:
            merged[-1][2] = max(merged[-1][2], e)
        else:
            merged.append([c, s, e])
    with open(f"{a.work}/regions.tsv", "w") as f:
        for c, s, e in merged:
            f.write(f"{c}\t{s}\t{e}\n")
    auth_pos = []
    for gene, rec in auth["loci"].items():
        for j, g in enumerate(rec.get("codon_gpos", [])):
            auth_pos.append([gene, rec["codon"], j + 1, rec["chrom"], g, fa.fetch(rec["chrom"], g - 1, g).upper()])
    write_tsv(f"{a.work}/auth_positions.tsv", ["gene", "codon", "codon_base", "chrom", "pos", "genomic_ref"], auth_pos)
    log(f"regions: {len(merged)} merged intervals for {len(rows)} sites")
    return 0 if not auth["problems"] else 3


# ------------------------------------------------------------------------------------------------------ counts
def tags_all():
    return [f"{lib}.new" for lib in LIBS] + [f"{lib}.old" for lib in LIBS]


def cmd_counts(a):
    fa = pysam.FastaFile(a.fasta)
    sites = list(read_tsv(f"{a.work}/sites.tsv.gz"))
    keyset = {(s["chrom"], int(s["pos"]), s["ref"], s["alt"]): s["site_id"] for s in sites}
    posset = {(s["chrom"], int(s["pos"])) for s in sites}
    auth_pos = list(read_tsv(f"{a.work}/auth_positions.tsv"))
    apos = {(r["chrom"], int(r["pos"])) for r in auth_pos}
    table = {s["site_id"]: {} for s in sites}
    authc = []
    for tag in tags_all():
        path = f"{a.mp_dir}/{tag}.vcf.gz"
        if not os.path.exists(path):
            log(f"no mpileup for {tag}")
            continue
        alle, posdp = {}, {}
        vf = pysam.VariantFile(path)
        for rec in vf:
            smp = rec.samples[0]
            ad = smp.get("AD"); adf = smp.get("ADF"); adr = smp.get("ADR")
            if not ad:
                continue
            ad = [x or 0 for x in ad]; adf = [x or 0 for x in (adf or [0] * len(ad))]; adr = [x or 0 for x in (adr or [0] * len(ad))]
            depth = sum(ad)
            indel = "INDEL" in rec.info
            cp = (rec.chrom, rec.pos)
            if not indel:
                if cp in posset:
                    posdp[cp] = depth
                if cp in apos:
                    authc.append([tag, rec.chrom, rec.pos, depth,
                                  ",".join(f"{al}:{ad[i]}" for i, al in enumerate(rec.alleles) if ad[i])])
            for i, alt in enumerate(rec.alts or (), start=1):
                if alt.startswith("<"):
                    continue
                key = normalize(fa, rec.chrom, rec.pos, rec.ref, alt)
                if key in keyset:
                    alle[key] = (ad[i], depth, adf[i], adr[i])
        for key, sid in keyset.items():
            if key in alle:
                table[sid][tag] = alle[key]
            else:
                table[sid][tag] = (0, posdp.get((key[0], key[1]), 0), 0, 0)
        log(f"counts {tag}: {len(alle)} site alleles observed, {len(posdp)} positions with depth")
    header = ["site_id"]
    for tag in tags_all():
        header += [f"{tag}_alt", f"{tag}_dp", f"{tag}_altF", f"{tag}_altR"]
    rows = []
    for s in sites:
        r = [s["site_id"]]
        for tag in tags_all():
            r += list(table[s["site_id"]].get(tag, (None, None, None, None)))
        rows.append(r)
    write_tsv(f"{a.work}/counts.tsv.gz", header, rows)
    write_tsv(f"{a.work}/auth_position_counts.tsv", ["tag", "chrom", "pos", "depth", "alleles"], authc)
    return 0


def load_counts(work):
    c = {}
    for r in read_tsv(f"{work}/counts.tsv.gz"):
        d = {}
        for tag in tags_all():
            alt = num(r.get(f"{tag}_alt")); dp = num(r.get(f"{tag}_dp"))
            if alt is not None and dp is not None:
                d[tag] = (alt, dp, num(r.get(f"{tag}_altF")) or 0, num(r.get(f"{tag}_altR")) or 0)
        c[r["site_id"]] = d
    return c


# ------------------------------------------------------------------------------------------------------ categories
def state(k, n):
    if n is None or n < MIN_READS:
        return "low"
    v = k / n
    if v >= 0.10 and k >= 3 and wilson(k, n)[0] > 0.05:
        return "present"
    if v <= 0.05:
        return "absent"
    return "grey"


def site_states(cnt, suffix=".new"):
    st, vaf = {}, {}
    for lib in LIBS:
        x = cnt.get(f"{lib}{suffix}")
        if x is None:
            st[lib] = "low"; vaf[lib] = None
        else:
            st[lib] = state(x[0], x[1]); vaf[lib] = x[0] / x[1] if x[1] else None
    return st, vaf


def categories(st, vaf):
    """Categories on the seven core libraries (C2B + six tumours), each needing >= 20 reads in every one used."""
    out = []
    if all(st[x] != "low" for x in CORE):
        if all(vaf[x] >= 0.3 for x in CORE):
            out.append("shared_clonal")
        if st["C2B"] == "present" and all(st[t] == "absent" for t in TUMOURS):
            out.append("c2b_present_lost_in_vivo")
        pres = [t for t in TUMOURS if st[t] == "present"]
        if st["C2B"] == "absent" and len(pres) == 1 and all(st[t] == "absent" for t in TUMOURS if t not in pres):
            out.append("gained_one_animal:" + pres[0])
        absn = [t for t in TUMOURS if st[t] == "absent"]
        if st["C2B"] == "present" and len(absn) == 1 and all(st[t] == "present" for t in TUMOURS if t not in absn):
            out.append("lost_one_animal:" + absn[0])
        if st["C2B"] == "absent" and all(st[t] == "present" for t in RECURRENT) and all(st[t] == "absent" for t in PRIMARY):
            out.append("gained_all_three_recurrences")
        if st["C2B"] == "absent" and all(st[t] == "present" for t in PRIMARY) and all(st[t] == "absent" for t in RECURRENT):
            out.append("gained_all_three_primaries")
    return out


def split_flags(st, vaf, A, B):
    """(gained in all of B, gained in all of A, |mean VAF B - mean VAF A| >= 0.2) for one split."""
    gb = ga = sh = False
    if all(st[x] != "low" for x in CORE):
        gb = st["C2B"] == "absent" and all(st[t] == "present" for t in B) and all(st[t] == "absent" for t in A)
        ga = st["C2B"] == "absent" and all(st[t] == "present" for t in A) and all(st[t] == "absent" for t in B)
    if all(st[t] != "low" for t in TUMOURS):
        sh = abs(mean([vaf[t] for t in B]) - mean([vaf[t] for t in A])) >= 0.2
    return gb, ga, sh


def cmd_select(a):
    sites = {s["site_id"]: s for s in read_tsv(f"{a.work}/sites.tsv.gz")}
    counts = load_counts(a.work)
    sp = splits()
    prio = {}
    shared = []
    for sid, s in sites.items():
        cnt = counts.get(sid, {})
        st, vaf = site_states(cnt)
        reasons = []
        if s["auth"]:
            reasons.append((0, "auth:" + s["auth"]))
        if s["in_union"] == "1":
            cats = categories(st, vaf)
            for k, (A, B) in enumerate(sp, start=1):
                gb, ga, sh = split_flags(st, vaf, A, B)
                if gb or ga:
                    reasons.append((1, f"split{k}_gained"))
                if sh:
                    reasons.append((5, f"split{k}_shift20"))
            for c in cats:
                if c.startswith("c2b_present"):
                    reasons.append((2, c))
                elif c.startswith("gained_one"):
                    reasons.append((3, c))
                elif c.startswith("lost_one"):
                    reasons.append((4, c))
                elif c == "shared_clonal":
                    shared.append(sid)
            if st["N269B"] != "low" and st["IL69B"] != "low" and abs(vaf["N269B"] - vaf["IL69B"]) >= 0.3:
                reasons.append((6, "n269b_vs_il69b"))
        if reasons:
            prio[sid] = (min(r[0] for r in reasons), ";".join(sorted({r[1] for r in reasons})))
    rng = random.Random(SEED)
    for sid in rng.sample(sorted(shared), min(300, len(shared))):
        if sid not in prio:
            prio[sid] = (7, "calibration_shared_clonal")
    ranked = sorted(prio.items(), key=lambda kv: (kv[1][0], kv[0]))
    rows = []
    for i, (sid, (p, why)) in enumerate(ranked):
        s = sites[sid]
        rows.append([sid, s["chrom"], s["pos"], s["ref"], s["alt"], s["kind"], p, why, int(i < MAX_CAND)])
    write_tsv(f"{a.work}/candidates.tsv", ["site_id", "chrom", "pos", "ref", "alt", "kind", "priority", "reasons", "examined"], rows)
    log(f"candidates: {len(rows)} ({sum(r[-1] for r in rows)} examined; cap {MAX_CAND})")
    with open(f"{a.work}/liftover_in.bed", "w") as f:
        for sid, s in sites.items():
            if s["kind"] == "snv":
                f.write(f"{s['chrom']}\t{int(s['pos']) - 1}\t{s['pos']}\t{sid}\t0\t+\n")
    return 0


# ------------------------------------------------------------------------------------------------------ readlevel
def cmd_readlevel(a):
    bam = pysam.AlignmentFile(a.bam)
    names = None
    if a.both_names:
        names = set()
        with gzip.open(a.both_names, "rt") as f:
            for line in f:
                names.add(line.rstrip("\n"))
        log(f"{a.tag}: {len(names)} both-bin read names loaded")
    cands = [c for c in read_tsv(a.candidates) if c["examined"] == "1"]
    out_rows, alt_rows = [], []
    fq1 = open(f"{a.out_prefix}.R1.fq.part", "w"); fq2 = open(f"{a.out_prefix}.R2.fq.part", "w")
    written = set(); n_pairs = 0
    allq = collections.Counter()   # both-bin hit rate over every read examined (name-format check)
    for c in cands:
        chrom, pos, ref, alt, kind = c["chrom"], int(c["pos"]), c["ref"], c["alt"], c["kind"]
        if kind not in ("snv", "ins", "del") or chrom not in bam.references:
            continue
        pos0 = pos - 1
        seen = {}
        trunc = 0
        for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="nofilter", ignore_overlaps=False,
                              ignore_orphans=False, min_base_quality=0, max_depth=MAX_READS_SITE, flag_filter=0):
            if col.reference_pos != pos0:
                continue
            trunc = int(col.nsegments >= MAX_READS_SITE)
            for pr in col.pileups:
                aln = pr.alignment
                if aln.flag & 0xF04 or aln.mapping_quality < 60:
                    continue
                if aln.is_paired and not aln.is_proper_pair:   # bcftools mpileup skips these too (no --count-orphans)
                    continue
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                qp = pr.query_position
                qq = aln.query_qualities
                q = qq[qp] if qq is not None else 40
                if q < 20:
                    continue
                seq = aln.query_sequence
                if kind == "snv":
                    b = seq[qp]
                    allele = "alt" if b == alt else ("ref" if b == ref else "oth")
                elif kind == "ins":
                    ins = alt[1:]
                    if pr.indel == len(ins) and seq[qp + 1:qp + 1 + len(ins)] == ins:
                        allele = "alt"
                    else:
                        allele = "ref" if pr.indel == 0 else "oth"
                else:
                    allele = "alt" if pr.indel == -(len(ref) - 1) else ("ref" if pr.indel == 0 else "oth")
                d5 = (aln.query_length - 1 - qp) if aln.is_reverse else qp
                qn = aln.query_name
                prev = seen.get(qn)
                if prev is None or q > prev[1]:          # overlapping mates: one fragment, the better base
                    seen[qn] = (allele, q, aln.is_read1, d5, aln)
        n = collections.Counter(v[0] for v in seen.values())
        altr = [v for v in seen.values() if v[0] == "alt"]
        refr = [v for v in seen.values() if v[0] == "ref"]
        def cnt(vs, r1, lim):
            return sum(1 for v in vs if v[2] == r1 and v[3] < lim)
        # transcript strand from read orientation (reverse-stranded library: R1 antisense, R2 sense)
        txp = sum(1 for v in seen.values() if (v[2] and v[4].is_reverse) or (not v[2] and not v[4].is_reverse))
        txm = len(seen) - txp
        both_alt = both_ref = None
        if names is not None:
            both_alt = sum(1 for k, v in seen.items() if v[0] == "alt" and k in names)
            both_ref = sum(1 for k, v in seen.items() if v[0] == "ref" and k in names)
            allq["examined"] += len(seen)
            allq["in_both"] += sum(1 for k in seen if k in names)
        out_rows.append([c["site_id"], n["ref"], n["alt"], n["oth"],
                         sum(1 for v in altr if v[2]), sum(1 for v in altr if not v[2]),
                         cnt(altr, True, 6), cnt(altr, False, 3), cnt(altr, False, 6),
                         sum(1 for v in refr if v[2]), sum(1 for v in refr if not v[2]),
                         cnt(refr, True, 6), cnt(refr, False, 3), cnt(refr, False, 6),
                         median([v[3] for v in altr]), median([v[3] for v in refr]), txp, txm, both_alt, both_ref, trunc])
        # alt-supporting pairs for the rat realignment
        k = 0
        for qn, v in seen.items():
            if v[0] != "alt" or k >= MAX_ALT_PAIRS_SITE or n_pairs >= MAX_ALT_PAIRS_BAM:
                continue
            aln = v[4]
            has = aln.has_tag("AS")
            alt_rows.append([c["site_id"], qn, aln.get_tag("AS") if has else None,
                             aln.get_tag("nM") if aln.has_tag("nM") else None, int(aln.is_read1),
                             None if names is None else int(qn in names)])
            k += 1
            if qn in written:
                continue
            try:
                mate = bam.mate(aln)
            except ValueError:
                continue
            recs = {True: None, False: None}
            for r in (aln, mate):
                s = r.query_sequence; qv = r.query_qualities
                if s is None or qv is None:
                    break
                ql = "".join(chr(x + 33) for x in qv)
                if r.is_reverse:
                    s, ql = revcomp(s), ql[::-1]
                recs[r.is_read1] = (s, ql)
            if recs[True] is None or recs[False] is None:
                continue
            nm = f"{a.tag}~{qn}"
            fq1.write(f"@{nm}\n{recs[True][0]}\n+\n{recs[True][1]}\n")
            fq2.write(f"@{nm}\n{recs[False][0]}\n+\n{recs[False][1]}\n")
            written.add(qn); n_pairs += 1
    fq1.close(); fq2.close()
    os.replace(f"{a.out_prefix}.R1.fq.part", f"{a.out_prefix}.R1.fq")
    os.replace(f"{a.out_prefix}.R2.fq.part", f"{a.out_prefix}.R2.fq")
    write_tsv(f"{a.out_prefix}.sites.tsv",
              ["site_id", "n_ref", "n_alt", "n_oth", "alt_r1", "alt_r2", "alt_r1_5p6", "alt_r2_5p3", "alt_r2_5p6",
               "ref_r1", "ref_r2", "ref_r1_5p6", "ref_r2_5p3", "ref_r2_5p6", "alt_d5_median", "ref_d5_median",
               "tx_plus", "tx_minus", "alt_both", "ref_both", "truncated"], out_rows)
    write_tsv(f"{a.out_prefix}.altreads.tsv", ["site_id", "qname", "human_AS", "human_nM", "is_read1", "is_both"], alt_rows)
    namecheck = {"tag": a.tag, "examined_reads": allq["examined"], "in_both_names": allq["in_both"]}
    if names is not None:
        # a random-ish sample of the whole BAM: the first 20,000 primary MAPQ >= 60 reads of chr1-chr5
        tot = hit = 0
        for chrom in [x for x in bam.references if x in ("chr1", "chr2", "chr3", "chr4", "chr5")]:
            k = 0
            for r in bam.fetch(chrom):
                if r.flag & 0xF04 or r.mapping_quality < 60:
                    continue
                k += 1
                if k > 20000:
                    break
                tot += 1; hit += int(r.query_name in names)
        namecheck.update({"sample_reads": tot, "sample_in_both_names": hit})
    with open(f"{a.out_prefix}.namecheck.json", "w") as f:
        json.dump(namecheck, f)
    log(f"readlevel {a.tag}: {len(out_rows)} sites, {n_pairs} alt pairs written")
    return 0


# ------------------------------------------------------------------------------------------------------ yx
def cmd_yx(a):
    bam = pysam.AlignmentFile(a.bam)
    ex = collections.defaultdict(list)
    for r in read_tsv(f"{a.work}/yx_exons.tsv"):
        ex[r["gene_name"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    rows = []
    for g, iv in sorted(ex.items()):
        chrom = iv[0][0]
        lo = min(s for _, s, _ in iv); hi = max(e for _, _, e in iv)
        seen = set()
        if chrom in bam.references:
            for r in bam.fetch(chrom, lo - 1, hi):
                if r.flag & 0xF04 or r.mapping_quality < 60 or not r.is_read1:
                    continue
                if r.query_name in seen:
                    continue
                for bs, be in r.get_blocks():
                    if any(bs < e and be > s - 1 for _, s, e in iv):
                        seen.add(r.query_name)
                        break
        rows.append([a.tag, g, chrom, len(seen)])
    write_tsv(a.out, ["tag", "gene", "chrom", "fragments"], rows)
    return 0


# ------------------------------------------------------------------------------------------------------ final
def rank_true(vals):
    """Rank of split 1 among the 10 (1 = largest); ties give a range. p = worst rank / 10, floor 0.10."""
    t = vals[0]
    if t is None or (isinstance(t, float) and math.isnan(t)):
        return (None, None, None)
    v = [x for x in vals[1:] if x is not None and not (isinstance(x, float) and math.isnan(x))]
    g = sum(1 for x in v if x > t); e = sum(1 for x in v if x == t)
    return (1 + g, 1 + g + e, (1 + g + e) / 10)


def cmd_final(a):
    O = a.out_dir
    os.makedirs(O, exist_ok=True)
    order = fai_order(a.fasta)
    sites = {s["site_id"]: s for s in read_tsv(f"{a.work}/sites.tsv.gz")}
    counts = load_counts(a.work)
    cands = {c["site_id"]: c for c in read_tsv(f"{a.work}/candidates.tsv")}
    auth = json.load(open(f"{a.work}/auth.json"))
    prep = {lib: load_prep(a.prep_checks, lib) for lib in LIBS}
    notes = []            # what the report must say about how each filter was applied
    outside = []          # checks of this computation against something outside it

    # ---- host fraction per library (xengsort), cross-checked with the prep_reads counts
    xs = {}
    for lib in LIBS:
        try:
            x = load_xengsort(a.u, lib)
        except Exception as e:
            notes.append(f"xengsort counts for {lib} unreadable: {e}")
            continue
        tot = sum(x.values())
        xs[lib] = {"host_frac": x["host"] / tot, "graft": x["graft"], "both": x["both"], "total": tot,
                   "both_share_human_file": x["both"] / (x["graft"] + x["both"])}
        same = all(str(x[k]) == prep[lib].get(f"n_{k}") for k in ("host", "graft", "both"))
        outside.append(["xengsort_counts_eq_prep_reads", lib, "host/graft/both equal", same, PASS if same else FAIL])

    # ---- read-level evidence per site (new and old BAMs)
    rl = collections.defaultdict(dict)
    for tag in tags_all():
        p = f"{a.rl_dir}/{tag}.sites.tsv"
        if os.path.exists(p):
            for r in read_tsv(p):
                rl[r["site_id"]][tag] = r
    namechecks = {}
    for p in glob.glob(f"{a.rl_dir}/*.namecheck.json"):
        j = json.load(open(p)); namechecks[j["tag"]] = j
    both_ok = True
    for lib in CORE:
        j = namechecks.get(f"{lib}.old")
        if not j or not j.get("sample_reads"):
            both_ok = False
            outside.append(["both_names_hit_rate", lib, "reads in both-bin names", "not measured", WARN])
            continue
        rate = j["sample_in_both_names"] / j["sample_reads"]
        exp = xs.get(lib, {}).get("both_share_human_file")
        ok = rate > 0 and exp and 0.2 * exp <= rate <= 5 * exp
        both_ok &= bool(ok)
        outside.append(["both_names_hit_rate", lib, f"same order as xengsort both share of the human file ({exp:.3f})" if exp else "NA",
                        f"{rate:.3f} of {j['sample_reads']} reads (chr1-5 sample, old BAM)", PASS if ok else FAIL])
    if not both_ok:
        notes.append("The both-bin read-name lookup failed its hit-rate check for at least one library; the both-bin "
                     "rat-leak and allelic-sorting filters are reported but not applied.")

    # ---- rat realignment: alt-supporting pairs, human AS (BAM) vs rat AS (STAR on mRatBN7.2)
    rat_as = {}
    if a.rat_sam and os.path.exists(a.rat_sam):
        with open(a.rat_sam) as f:
            for line in f:
                if line.startswith("@"):
                    continue
                p = line.split("\t", 12)
                flag = int(p[1]); qn = p[0]
                if flag & 0x100 or flag & 0x800:
                    continue
                if flag & 0x4:
                    rat_as.setdefault(qn, -1)
                    continue
                m = re.search(r"\tAS:i:(-?\d+)", line)
                if m:
                    rat_as[qn] = max(rat_as.get(qn, -1), int(m.group(1)))
    else:
        notes.append("rat realignment output missing: the rat-realignment filter is not applied")
    ratfit = collections.defaultdict(lambda: collections.Counter())
    for tag in tags_all():
        p = f"{a.rl_dir}/{tag}.altreads.tsv"
        if not os.path.exists(p):
            continue
        lib = tag.split(".")[0]
        for r in read_tsv(p):
            h = num(r["human_AS"]); ra = rat_as.get(f"{tag}~{r['qname']}")
            if h is None or ra is None:
                continue
            side = "new" if tag.endswith(".new") else "old"
            if lib in CORE:
                ratfit[r["site_id"]][f"{side}_n"] += 1
                ratfit[r["site_id"]][f"{side}_rat"] += int(ra > h)

    # ---- liftover to rn7 and the rat base (on the human + strand)
    lift = {}
    if a.lifted and os.path.exists(a.lifted):
        rat = pysam.FastaFile(a.rat_fasta)
        ratn = set(rat.references)
        for line in open(a.lifted):
            p = line.rstrip("\n").split("\t")
            c = p[0]; c2 = "MT" if c == "chrM" else (c[3:] if c.startswith("chr") else c)
            if c2 not in ratn:
                continue
            b = rat.fetch(c2, int(p[1]), int(p[2])).upper()
            if p[5] == "-":
                b = b.translate(COMP)
            lift[p[3]] = (b, p[5])
        ident = collections.Counter()
        for sid, (b, strand) in lift.items():
            ident[(strand, "n")] += 1
            ident[(strand, "same")] += int(b == sites[sid]["ref"])
        for strand in "+-":
            n = ident[(strand, "n")]
            outside.append(["liftover_human_ref_eq_rat_base", strand, "high (orthologous, mostly exonic); '-' close to '+'",
                            f"{ident[(strand, 'same')] / n:.3f} of {n}" if n else "none", INFO])
        if ident[("-", "n")] >= 50 and ident[("+", "n")] >= 50:
            fp = ident[("+", "same")] / ident[("+", "n")]; fm = ident[("-", "same")] / ident[("-", "n")]
            if fm < 0.5 * fp:
                notes.append("rat-base filter NOT applied: '-' strand identity is under half the '+' strand identity (strand bug?)")
                lift = {}
    else:
        notes.append("liftover output missing: the rat-base filter is not applied")

    # ---- RNA editing: REDIportal if staged, else the strand rule (and always the read-orientation strand)
    snv_pos = {(s["chrom"], int(s["pos"])): sid for sid, s in sites.items() if s["kind"] == "snv"}
    redi = set()
    use_redi = bool(a.rediportal and os.path.exists(a.rediportal))
    if use_redi:
        with gzip.open(a.rediportal, "rt") as f:
            h = f.readline().rstrip("\n").split("\t")
            ci = {k: i for i, k in enumerate(h)}
            for line in f:
                p = line.split("\t")
                sid = snv_pos.get((p[ci["Region"]], num(p[ci["Position"]])))
                if sid and sites[sid]["ref"] == p[ci["Ref"]] and sites[sid]["alt"] == p[ci["Ed"]]:
                    redi.add(sid)
        edit_rule = f"REDIportal hg38 v3 ({os.path.basename(a.rediportal)})"
    else:
        edit_rule = "strand rule (REDIportal not staged): A>G on + strand genes, T>C on - strand genes"
    gstrand = collections.defaultdict(list)
    for r in read_tsv(f"{a.work}/gtf_genes.tsv.gz"):
        gstrand[(r["chrom"], int(r["start"]) // 100000)].append((int(r["start"]), int(r["end"]), r["strand"]))
    maxlen = max((int(r["end"]) - int(r["start"]) for r in read_tsv(f"{a.work}/gtf_genes.tsv.gz")), default=0)
    nbins = maxlen // 100000 + 2
    def gene_strand(c, p):
        # genes are binned by start; scan back far enough to reach the longest gene in the GTF
        st = set()
        for k in range(nbins):
            st |= {s for (b0, b1, s) in gstrand.get((c, p // 100000 - k), ()) if b0 <= p <= b1}
        return st.pop() if len(st) == 1 else ("both" if st else "none")
    edit_strand = {}
    ag_match = ag_mis = 0
    for sid, s in sites.items():
        if s["kind"] != "snv":
            continue
        gs = gene_strand(s["chrom"], int(s["pos"]))
        ra = (s["ref"], s["alt"])
        e = (gs == "+" and ra == ("A", "G")) or (gs == "-" and ra == ("T", "C"))
        edit_strand[sid] = e
        if s["in_union"] == "1" and gs in "+-" and gs:
            ag_match += int(e)
            ag_mis += int((gs == "+" and ra == ("T", "C")) or (gs == "-" and ra == ("A", "G")))
    outside.append(["editing_strand_excess", "-", "A>G on the transcribed strand exceeds its mirror (A-to-I editing)",
                    f"{ag_match} strand-matched vs {ag_mis} mirrored", INFO])

    # ---- per-site VAF, CI, flags, categories
    tumour_host = [xs.get(t, {}).get("host_frac") for t in TUMOURS]
    rows = []
    flagged = {}
    excluded_main = {}
    sp = splits()
    split_counts = {fs: [[0, 0, 0] for _ in sp] for fs in ("raw", "main", "strict")}
    cat_counts = {fs: collections.Counter() for fs in ("raw", "main", "strict")}
    per_site_cats = {}
    for sid, s in sites.items():
        cnt = counts.get(sid, {})
        st, vaf = site_states(cnt)
        sto, vafo = site_states(cnt, ".old")
        fl = []
        c = cands.get(sid)
        examined = bool(c and c["examined"] == "1")
        r = rl.get(sid, {})
        # (i) rat leak via the both bin (old graft+both BAMs), pooled over the core libraries
        if both_ok and examined:
            ab = sum(num(r[f"{x}.old"]["alt_both"]) or 0 for x in CORE if f"{x}.old" in r)
            an = sum(num(r[f"{x}.old"]["n_alt"]) or 0 for x in CORE if f"{x}.old" in r)
            if an >= 3 and ab / an > 0.5:
                fl.append("both_bin_alt")
        # (i) rat realignment
        rf = ratfit.get(sid)
        if rf and rf["new_n"] >= 3 and rf["new_rat"] / rf["new_n"] > 0.5:
            fl.append("rat_realign")
        # (i) alt equals the rat base
        if sid in lift and s["kind"] == "snv" and lift[sid][0] == s["alt"]:
            fl.append("alt_is_rat_base")
        # (ii) VAF tracks the host fraction: binomial logistic regression of alt/depth on the xengsort host fraction
        # over C2B (host ~0.5 %, anchors the low end) and the six tumours; the rat-leak shape is C2B at or below 5 %
        # and VAF rising with host share. A likelihood-ratio test, not a rank correlation of six noisy VAFs.
        if all(st[x] != "low" for x in CORE) and None not in tumour_host and xs.get("C2B") and vaf["C2B"] <= 0.05:
            b, pv = logit_lrt([cnt[f"{x}.new"][0] for x in CORE], [cnt[f"{x}.new"][1] for x in CORE],
                              [xs["C2B"]["host_frac"]] + tumour_host)
            if b > 0 and pv < 0.001:
                fl.append("tracks_host_fraction")
        # (iii) position in read
        if examined:
            def pool(side, key):
                return sum(num(r[f"{x}.{side}"][key]) or 0 for x in CORE if f"{x}.{side}" in r)
            a1, a15 = pool("new", "alt_r1"), pool("new", "alt_r1_5p6")
            r1, r15 = pool("new", "ref_r1"), pool("new", "ref_r1_5p6")
            if a1 >= 3 and a15 / a1 >= 0.5 and a15 / a1 - (r15 / r1 if r1 else 0) >= 0.3:
                fl.append("r1_5prime_6nt")
            a2, a23 = pool("old", "alt_r2"), pool("old", "alt_r2_5p3")
            r2, r23 = pool("old", "ref_r2"), pool("old", "ref_r2_5p3")
            if a2 >= 3 and a23 / a2 >= 0.5 and a23 / a2 - (r23 / r2 if r2 else 0) >= 0.3:
                fl.append("r2_5prime_3nt_old")
        # (iv) editing
        if (sid in redi) if use_redi else edit_strand.get(sid, False):
            fl.append("editing")
        # (v) allelic sorting: alt share among graft-bin vs both-bin reads (old BAMs, core libraries)
        if both_ok and examined:
            ab = sum(num(r[f"{x}.old"]["alt_both"]) or 0 for x in CORE if f"{x}.old" in r)
            rb = sum(num(r[f"{x}.old"]["ref_both"]) or 0 for x in CORE if f"{x}.old" in r)
            an = sum(num(r[f"{x}.old"]["n_alt"]) or 0 for x in CORE if f"{x}.old" in r)
            rn = sum(num(r[f"{x}.old"]["n_ref"]) or 0 for x in CORE if f"{x}.old" in r)
            ag, rg = an - ab, rn - rb
            if ag + rg >= 10 and ab + rb >= 10:
                vg, vb = ag / (ag + rg), ab / (ab + rb)
                if abs(vg - vb) >= 0.2 and fisher_two_sided(ag, rg, ab, rb) < 0.01:
                    fl.append("allelic_sorting")
        common = (num(s["gnomad_af"], float) or 0) > 0.01
        main_ex = any(f in fl for f in ("both_bin_alt", "rat_realign", "r1_5prime_6nt", "r2_5prime_3nt_old", "editing",
                                        "allelic_sorting")) or ("alt_is_rat_base" in fl and not common)
        strict_ex = main_ex or "alt_is_rat_base" in fl or "tracks_host_fraction" in fl
        flagged[sid] = fl
        excluded_main[sid] = main_ex
        cats = categories(st, vaf) if s["in_union"] == "1" else []
        per_site_cats[sid] = cats
        for fs, ex in (("raw", False), ("main", main_ex), ("strict", strict_ex)):
            if ex or s["in_union"] != "1":
                continue
            for cc in cats:
                cat_counts[fs][cc.split(":")[0]] += 1
                if ":" in cc:
                    cat_counts[fs][cc] += 1
            for k, (A, B) in enumerate(sp):
                gb, ga, sh = split_flags(st, vaf, A, B)
                split_counts[fs][k][0] += gb; split_counts[fs][k][1] += ga; split_counts[fs][k][2] += sh
        row = [sid, s["chrom"], s["pos"], s["ref"], s["alt"], s["kind"], s["in_union"], s["auth"], s["symbol"],
               s["consequence"], s["hgvsp"], s["gnomad_af"], s["n_pass"], s["filters"], ",".join(fl) or "none",
               int(main_ex), int(strict_ex), ",".join(cats) or "none", c["reasons"] if c else "", int(examined)]
        for lib in LIBS:
            x = cnt.get(f"{lib}.new")
            if x:
                lo, hi = wilson(x[0], x[1])
                row += [x[0], x[1], x[0] / x[1] if x[1] else None, lo, hi, st[lib]]
            else:
                row += [None] * 5 + ["low"]
        for lib in LIBS:
            x = cnt.get(f"{lib}.old")
            row += [x[0], x[1], x[0] / x[1] if x and x[1] else None] if x else [None] * 3
        rows.append(row)
    header = ["site_id", "chrom", "pos", "ref", "alt", "kind", "in_union", "auth", "symbol", "consequence", "hgvsp",
              "gnomad_af", "n_pass", "hc_filters", "flags", "excluded_main", "excluded_strict", "categories", "candidate_reasons",
              "read_level_examined"]
    for lib in LIBS:
        header += [f"{lib}_alt", f"{lib}_dp", f"{lib}_vaf", f"{lib}_ci_lo", f"{lib}_ci_hi", f"{lib}_state"]
    for lib in LIBS:
        header += [f"{lib}_old_alt", f"{lib}_old_dp", f"{lib}_old_vaf"]
    write_tsv(f"{O}/sites_vaf.tsv.gz", header, rows)

    # ---- flag and category summaries
    fcount = collections.Counter(f for fl in flagged.values() for f in fl)
    write_tsv(f"{O}/filter_summary.tsv", ["flag", "sites"], sorted(fcount.items()))
    cat_rows = []
    for fs in ("raw", "main", "strict"):
        for k, v in sorted(cat_counts[fs].items()):
            cat_rows.append([fs, k, v])
    write_tsv(f"{O}/categories.tsv", ["filter_set", "category", "sites"], cat_rows)

    # ---- split ranks for site counts
    rank_rows = []
    def add_rank(stat, fs, vals):
        lo, hi, p = rank_true(vals)
        rank_rows.append([stat, fs] + vals + [f"{lo}-{hi}" if lo != hi else lo, max(p, 0.10) if p else None])
    for fs in ("raw", "main", "strict"):
        add_rank("sites_gained_in_all_of_B (B = recurrent in split 1)", fs, [x[0] for x in split_counts[fs]])
        add_rank("sites_gained_in_all_of_A (A = primary in split 1)", fs, [x[1] for x in split_counts[fs]])
        add_rank("sites_gained_in_one_arm_only (A or B)", fs, [x[0] + x[1] for x in split_counts[fs]])
        add_rank("sites_mean_VAF_shift_ge_0.2", fs, [x[2] for x in split_counts[fs]])

    # ---- arm-level allelic imbalance at heterozygous common SNPs
    arms = {}
    cyto_ok = False
    if a.cytoband and os.path.exists(a.cytoband):
        pe, qe = {}, {}
        op = gzip.open if a.cytoband.endswith(".gz") else open
        for line in op(a.cytoband, "rt"):
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            if p[3].startswith("p"):
                pe[p[0]] = max(pe.get(p[0], 0), int(p[2]))
            qe[p[0]] = max(qe.get(p[0], 0), int(p[2]))
        prim = [c for c in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"] if c in order]
        lens_ok = bool(prim) and all(qe.get(c) == order[c][1] for c in prim)
        outside.append(["cytoband_chrom_lengths_eq_gencode_fai", "-", "chr1-22,X,Y equal", lens_ok, PASS if lens_ok else FAIL])
        if lens_ok:
            arms = pe; cyto_ok = True
    if not cyto_ok:
        notes.append("cytoBand unavailable or inconsistent: allelic imbalance is by whole chromosome, not arm")
    # pseudoautosomal regions on chrX, from the X copies of the genes GENCODE also places on chrY as <id>_PAR_Y
    # (PAR2 sits at different coordinates on X and Y, so the Y copies' coordinates cannot be used)
    gl = list(read_tsv(f"{a.work}/gtf_genes.tsv.gz"))
    par_ids = {r["gene_id"][:-len("_PAR_Y")] for r in gl if r["gene_id"].endswith("_PAR_Y")}
    par = [(int(r["start"]), int(r["end"])) for r in gl if r["chrom"] == "chrX" and r["gene_id"] in par_ids]
    del gl
    par_lo = [x for x in par if x[1] < 10_000_000]; par_hi = [x for x in par if x[0] > 100_000_000]
    outside.append(["chrX_PAR_from_GTF", "-", "two blocks (PAR1 near the start, PAR2 near the end of chrX)",
                    f"PAR1 genes {len(par_lo)}, PAR2 genes {len(par_hi)}", INFO])
    par_iv = []
    if par_lo:
        par_iv.append((min(x[0] for x in par_lo), max(x[1] for x in par_lo)))
    if par_hi:
        par_iv.append((min(x[0] for x in par_hi), max(x[1] for x in par_hi)))
    def arm_of(c, p):
        if c == "chrX" and any(lo <= p <= hi for lo, hi in par_iv):
            return None
        if c not in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
            return None
        if cyto_ok and c in arms:
            return c[3:] + ("p" if p <= arms[c] else "q")
        return c[3:]
    snp_ids = [sid for sid, s in sites.items() if s["kind"] == "snv" and (num(s["gnomad_af"], float) or 0) > 0.01
               and "editing" not in flagged.get(sid, []) and s["in_union"] == "1"]
    het_c2b = []
    for sid in snp_ids:
        x = counts[sid].get("C2B.new")
        if x and x[1] >= MIN_READS and 0.1 <= x[0] / x[1] <= 0.9:
            het_c2b.append(sid)
    armvals = collections.defaultdict(lambda: collections.defaultdict(list))
    for sid in het_c2b:
        s = sites[sid]
        arm = arm_of(s["chrom"], int(s["pos"]))
        if not arm:
            continue
        for lib in CORE + ["N269B"]:
            x = counts[sid].get(f"{lib}.new")
            if x and x[1] >= MIN_READS:
                armvals[arm][lib].append(abs(x[0] / x[1] - 0.5))
    arm_rows = []
    arm_med = {}
    for arm in sorted(armvals, key=lambda z: (int(re.sub(r"[pq]$", "", z)) if re.sub(r"[pq]$", "", z).isdigit() else 23, z)):
        for lib in CORE + ["N269B"]:
            v = armvals[arm][lib]
            m = median(v) if len(v) >= 10 else None
            arm_med[(arm, lib)] = m
            c2 = arm_med.get((arm, "C2B"))
            arm_rows.append([arm, lib, len(v), m, (m - c2) if (m is not None and c2 is not None) else None])
    write_tsv(f"{O}/arm_imbalance.tsv", ["arm", "library", "het_common_snps", "median_abs_vaf_minus_half",
                                         "delta_vs_C2B"], arm_rows)
    arm_stats = []
    for A, B in sp:
        tot = 0.0; n10 = 0; used = 0
        for arm in armvals:
            va = [arm_med.get((arm, t)) for t in A]; vb = [arm_med.get((arm, t)) for t in B]
            if None in va or None in vb:
                continue
            d = abs(mean(vb) - mean(va)); tot += d; n10 += int(d >= 0.1); used += 1
        arm_stats.append((tot, n10, used))
    add_rank("arm_imbalance_sum_abs_diff (arms with >= 10 het SNPs in all six)", "common SNPs", [x[0] for x in arm_stats])
    add_rank("arms_with_abs_diff_ge_0.1", "common SNPs", [x[1] for x in arm_stats])
    # het SNPs a tumour has that C2B does not: another way to see non-U251 human reads
    loh_rows = []
    for lib in TUMOURS + ["N269B"]:
        k = n = 0
        for sid in snp_ids:
            xc = counts[sid].get("C2B.new"); xt = counts[sid].get(f"{lib}.new")
            if not xc or not xt or xc[1] < 50 or xt[1] < MIN_READS:
                continue
            vc = xc[0] / xc[1]; vt = xt[0] / xt[1]
            if vc < 0.02 or vc > 0.98:
                n += 1
                k += int(0.2 <= vt <= 0.8)
        loh_rows.append([lib, n, k, k / n if n else None])
    xhet = []
    for lib in CORE + ["N269B"]:
        k = n = 0
        for sid in snp_ids:
            s = sites[sid]
            if s["chrom"] != "chrX" or arm_of("chrX", int(s["pos"])) is None:
                continue
            x = counts[sid].get(f"{lib}.new")
            if x and x[1] >= MIN_READS:
                n += 1; k += int(0.2 <= x[0] / x[1] <= 0.8)
        xhet.append([lib, n, k])
    write_tsv(f"{O}/het_in_tumour_hom_in_c2b.tsv", ["library", "c2b_homozygous_common_snps_ge50_reads",
                                                    "heterozygous_in_library", "fraction"], loh_rows)
    write_tsv(f"{O}/chrX_nonPAR_het_snps.tsv", ["library", "common_snps_ge20_reads", "heterozygous_0.2_0.8"], xhet)

    # ---- Y/X gametolog ratios (male U251 in female rats)
    yx = collections.defaultdict(dict)
    for p in glob.glob(f"{a.yx_dir}/*.yx.tsv"):
        for r in read_tsv(p):
            yx[r["tag"].split(".")[0]][r["gene"]] = int(r["fragments"])
    yx_rows = []
    yx_l2 = {}
    for lib in LIBS:
        g = yx.get(lib)
        if not g:
            continue
        ys = sum(g.get(y, 0) for y, _ in YX_PAIRS); xsum = sum(g.get(x, 0) for _, x in YX_PAIRS)
        l2 = math.log2((ys + 0.5) / (xsum + 0.5))
        yx_l2[lib] = l2
        yx_rows.append([lib, ys, xsum, l2] + [f"{g.get(y, 0)}/{g.get(x, 0)}" for y, x in YX_PAIRS])
    write_tsv(f"{O}/y_over_x.tsv", ["library", "Y_fragments", "X_fragments", "log2_Y_over_X"] + [f"{y}/{x}" for y, x in YX_PAIRS],
              yx_rows)
    if all(t in yx_l2 for t in TUMOURS):
        add_rank("Y_over_X_abs_diff_mean_log2", "BAM counts", [abs(mean([yx_l2[t] for t in B]) - mean([yx_l2[t] for t in A]))
                                                               for A, B in sp])
    inv = sorted(glob.glob(f"{a.invasion_dir}/*y*x*.tsv")) if a.invasion_dir else []
    notes.append(f"invasion_niche Y/X table: {inv[0] if inv else 'not available when this ran; Y/X counted from the rnavar BAMs here'}")

    split_hdr = ["statistic", "filter_set"] + [f"split{k}" for k in range(1, 11)] + ["true_split_rank", "p_floor_0.10"]
    write_tsv(f"{O}/split_ranks.tsv", split_hdr, rank_rows)

    # ---- index-hopping / carry-over probe: private alleles of each donor in the low-human libraries
    idx = {}
    for lib in LIBS:
        p = f"{a.index_dir}/{lib}.txt"
        if os.path.exists(p):
            for line in open(p):
                q = line.rstrip("\n").split("\t")
                if q[0] == "modal_index" and "+" in q[1]:
                    idx[lib] = q[1].split("+")
    hop_rows = []
    for donor in CORE:
        priv = [sid for sid, cc in per_site_cats.items()
                if (f"gained_one_animal:{donor}" in cc or (donor == "C2B" and "c2b_present_lost_in_vivo" in cc))
                and not excluded_main.get(sid)]
        for rec in RECIPIENTS:
            k = n = 0
            for sid in priv:
                x = counts[sid].get(f"{rec}.new")
                if x:
                    k += x[0]; n += x[1]
            shared = []
            if donor in idx and rec in idx:
                if idx[donor][0] == idx[rec][0]:
                    shared.append("i7")
                if idx[donor][1] == idx[rec][1]:
                    shared.append("i5")
            dv = [counts[sid][f"{donor}.new"][0] / counts[sid][f"{donor}.new"][1] for sid in priv
                  if counts[sid].get(f"{donor}.new") and counts[sid][f"{donor}.new"][1]]
            hop_rows.append([donor, rec, len(priv), k, n, k / n if n else None, mean(dv) if dv else None,
                             "+".join(shared) or "none"])
    write_tsv(f"{O}/private_allele_probe.tsv", ["donor", "recipient", "donor_private_sites", "recipient_alt", "recipient_depth",
                                                "recipient_vaf", "donor_mean_vaf", "shared_index"], hop_rows)

    # ---- N269B against IL69B (same rat; descriptive, n = 1)
    n2 = []
    for sid, s in sites.items():
        cnt = counts.get(sid, {})
        x, y = cnt.get("N269B.new"), cnt.get("IL69B.new")
        if s["in_union"] == "1" and x and y and x[1] >= MIN_READS and y[1] >= MIN_READS and abs(x[0] / x[1] - y[0] / y[1]) >= 0.3:
            n2.append([sid, s["chrom"], s["pos"], s["ref"], s["alt"], s["symbol"], x[0], x[1], x[0] / x[1], y[0], y[1], y[0] / y[1],
                       ",".join(flagged.get(sid, [])) or "none"])
    write_tsv(f"{O}/n269b_vs_il69b.tsv", ["site_id", "chrom", "pos", "ref", "alt", "symbol", "N269B_alt", "N269B_dp",
                                         "N269B_vaf", "IL69B_alt", "IL69B_dp", "IL69B_vaf", "flags"], n2)

    # ---- authentication
    auth_rows = []
    def vaf_row(label, sid, rule_lib_min):
        cnt = counts.get(sid, {})
        for lib in LIBS:
            for side in ("new", "old"):
                x = cnt.get(f"{lib}.{side}")
                if not x:
                    continue
                lo, hi = wilson(x[0], x[1])
                auth_rows.append([label, lib, side, x[0], x[1], x[0] / x[1] if x[1] else None, lo, hi, ""])
    sid_by_key = {(s["chrom"], int(s["pos"]), s["ref"], s["alt"]): sid for sid, s in sites.items()}
    verdict = []
    loci = auth.get("loci", {}); vk = auth.get("vep_keys", {})
    for p in auth.get("problems", []):
        verdict.append(["locus_mapping", "FAIL", p])
    # TP53
    tp = loci.get("TP53", {})
    ek = tuple(tp.get("expected_key", [])) or None
    vkey = tuple(vk.get("TP53_R273H", [])) or None
    verdict.append(["TP53_R273H_key", PASS if (vkey and ek == vkey) else (WARN if ek else FAIL),
                    f"GTF-derived {ek}; VEP-called {vkey}"])
    sid = sid_by_key.get(ek) if ek else None
    if not sid:
        verdict.append(["TP53_R273H_homozygous", FAIL, "expected TP53 c.818G>A key not among the counted sites"])
    else:
        vaf_row("TP53 p.Arg273His (c.818G>A)", sid, MIN_READS)
        cnt = counts.get(sid, {})
        bad = []
        for lib in CORE:
            x = cnt.get(f"{lib}.new")
            if not x or x[1] < MIN_READS:
                if lib == "C2B":
                    bad.append(f"C2B: < {MIN_READS} reads")
                continue
            v = x[0] / x[1]
            if v < 0.9:
                bad.append(f"{lib}: VAF {v:.2f}" + (" (0.5-0.8: non-U251 human or rat reads)" if 0.5 <= v <= 0.8 else ""))
        verdict.append(["TP53_R273H_homozygous", PASS if not bad else FAIL,
                        "VAF >= 0.9 in C2B and every tumour with >= 20 reads" + (": " + "; ".join(bad) if bad else "")])
    # PTEN
    pt = loci.get("PTEN", {})
    ek = tuple(pt.get("expected_key", [])) or None
    vkey = tuple(vk.get("PTEN_c723_724dupTG", [])) or None
    verdict.append(["PTEN_c723_724dupTG_key", PASS if (vkey and ek == vkey) else (WARN if ek else FAIL),
                    f"GTF-derived {ek}; VEP-called {vkey}"])
    sid = sid_by_key.get(ek) if ek else None
    if not sid:
        verdict.append(["PTEN_dup_present", FAIL, "expected PTEN c.723_724dupTG key not among the counted sites"])
    else:
        vaf_row("PTEN p.Glu242ValfsTer15 (c.723_724dupTG)", sid, MIN_READS)
        cnt = counts.get(sid, {})
        miss = [lib for lib in CORE if cnt.get(f"{lib}.new") and cnt[f"{lib}.new"][1] >= MIN_READS
                and state(cnt[f"{lib}.new"][0], cnt[f"{lib}.new"][1]) != "present"]
        cov = [lib for lib in CORE if cnt.get(f"{lib}.new") and cnt[f"{lib}.new"][1] >= MIN_READS]
        verdict.append(["PTEN_dup_present", PASS if (cov and not miss and "C2B" in cov) else FAIL,
                        f"present in {len(cov) - len(miss)} of {len(cov)} core libraries with >= 20 reads" +
                        (f"; not present: {miss}" if miss else "")])
    # IDH1 codon 132 reference
    ap = collections.defaultdict(dict)
    for r in read_tsv(f"{a.work}/auth_position_counts.tsv"):
        ap[r["tag"]][(r["chrom"], int(r["pos"]))] = (int(r["depth"]), r["alleles"])
    idh = loci.get("IDH1", {})
    grefs = {int(r["pos"]): r["genomic_ref"] for r in read_tsv(f"{a.work}/auth_positions.tsv") if r["gene"] == "IDH1"}
    bad = []; seen_depth = []
    for lib in CORE:
        for j, g in enumerate(idh.get("codon_gpos", []), start=1):
            d, al = ap.get(f"{lib}.new", {}).get((idh["chrom"], g), (0, ""))
            gref = grefs.get(g)
            refc = sum(int(z.rsplit(":", 1)[1]) for z in al.split(",") if z and z.rsplit(":", 1)[0] == gref)
            auth_rows.append(["IDH1 codon 132 base %d (%s:%d)" % (j, idh["chrom"], g), lib, "new", d - refc, d,
                              (d - refc) / d if d else None, None, None, al])
            seen_depth.append(d)
            if d >= MIN_READS and refc / d < 0.95:
                bad.append(f"{lib} base {j}: ref {refc}/{d}")
    called = vk.get("IDH1_R132_called")
    verdict.append(["IDH1_R132_reference", PASS if (not bad and not called and seen_depth and max(seen_depth) >= MIN_READS) else FAIL,
                    f"ref fraction >= 0.95 at all three codon bases where >= 20 reads; VEP IDH1 p.Arg132 call: {called or 'none'}"
                    + (f"; {bad}" if bad else "")])
    # HLA class I (seq2HLA 4-digit, classical)
    exp = {"A": "A*02:01", "B": "B*18:01", "C": "C*05:01"}
    hla_rows = []
    for lib in LIBS:
        p = p_seq2hla(a.results, lib)
        if not os.path.exists(p):
            hla_rows.append([lib, "", "", "", "", "ignored (too few reads)" if lib in SMALL else "missing"])
            continue
        # seq2HLA 4-digit table, parsed as: locus, allele 1, confidence, allele 2, confidence (tab-separated; an
        # apostrophe marks an allele seq2HLA flags). The raw lines are kept so a reader can check the parse.
        calls, raw = {}, []
        for line in open(p):
            if not line.strip():
                continue
            raw.append(line.rstrip("\n").replace("\t", " "))
            if line.startswith("#"):
                continue
            q = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(q) >= 4:
                calls[q[0]] = (q[1], q[2], q[3], q[4] if len(q) > 4 else "")
        ok = True; desc = {}
        norm = lambda z: z.replace("'", "").replace(":", "").replace("HLA-", "").strip()
        for loc, want in exp.items():
            c = calls.get(loc)
            if not c:
                ok = False; desc[loc] = "no call"; continue
            w = norm(want)
            m1 = w in [norm(x) for x in re.split(r"[,/]", c[0])]; m2 = w in [norm(x) for x in re.split(r"[,/]", c[2])]
            ok &= m1 and m2
            desc[loc] = f"{c[0]} ({c[1]}) / {c[2]} ({c[3]})"
        hla_rows.append([lib, desc.get("A", ""), desc.get("B", ""), desc.get("C", ""), " || ".join(raw),
                         "ignored (too few reads)" if lib in SMALL else ("match" if ok else "MISMATCH")])
    write_tsv(f"{O}/hla_seq2hla.tsv", ["library", "HLA-A", "HLA-B", "HLA-C", "raw_4digit_file",
                                       "vs_A*02:01_B*18:01_C*05:01_homozygous"], hla_rows)
    large_hla = [r for r in hla_rows if r[0] in LARGE]
    nmatch = sum(1 for r in large_hla if r[-1] == "match")
    verdict.append(["HLA_classI_A0201_B1801_C0501_homozygous", PASS if large_hla and nmatch == len(large_hla) else
                    (WARN if nmatch else FAIL), f"{nmatch} of {len(large_hla)} large libraries match (PubMed 25960936)"])
    chrx = {r[0]: r for r in xhet}
    verdict.append(["chrX_nonPAR_het_common_SNPs_C2B (male line)", INFO, f"{chrx.get('C2B', ['', 0, 0])[2]} of "
                    f"{chrx.get('C2B', ['', 0, 0])[1]} common SNPs heterozygous"])
    write_tsv(f"{O}/authentication.tsv", ["locus", "library", "bam", "alt", "depth", "vaf", "ci_lo", "ci_hi", "alleles"], auth_rows)
    write_tsv(f"{O}/authentication_verdict.tsv", ["check", "result", "detail"], verdict)
    write_tsv(f"{O}/outside_checks.tsv", ["check", "scope", "expected", "observed", "result"], outside)

    # ---- report
    rep = []
    rep.append("# rnavar integration: authentication and clonal composition after LITT\n")
    rep.append("Generated by `rnavar_integrate.py final`. Inputs: nf-core/rnavar 1.3.0 on graft-only reads "
               "(`M/rnavar/results`), forced `bcftools mpileup` counts on the new graft-only BAMs and on the published "
               "graft+both BAMs, the rat STAR index and `hg38ToRn7` chain. Tables are in this folder.\n")
    rep.append("## Authentication (Cellosaurus CVCL_0021)\n")
    rep.append("| check | result | detail |\n|---|---|---|")
    for v in verdict:
        rep.append(f"| {v[0]} | {v[1]} | {v[2]} |")
    rep.append("\nThis identifies the U-251MG donor lineage, which includes SNB-19 and U-373MG (ATCC). It does not identify "
               "this lab's stock specifically. A disagreement only with the 1977 serology (A*01,02; B*08,14) is not "
               "misidentification.\n")
    rep.append("## Categories (seven core libraries, >= 20 reads in each, binomial 95 % Wilson intervals in sites_vaf.tsv.gz)\n")
    rep.append("| filter set | category | sites |\n|---|---|--:|")
    for r in cat_rows:
        rep.append(f"| {r[0]} | {r[1]} | {r[2]} |")
    rep.append("\nFilter sets: `raw` = PASS in >= 1 library, no artefact filter; `main` = excludes both-bin rat leak, rat "
               "realignment, read-position, editing, allelic-sorting flags, and alt-equals-rat-base when the allele is not "
               "common in gnomAD; `strict` = main plus every alt-equals-rat-base and every host-fraction-tracking site.\n")
    rep.append("## Ten-split ranks (split 1 = the true labels; p cannot go below 0.10)\n")
    rep.append("| statistic | set | true split | 10 values | rank | p |\n|---|---|--:|---|---|--:|")
    for r in rank_rows:
        rep.append(f"| {r[0]} | {r[1]} | {fmt(r[2])} | {' '.join(fmt(x) for x in r[2:12])} | {r[12]} | {fmt(r[13])} |")
    rep.append("\n## Artefact flags\n")
    rep.append("| flag | sites |\n|---|--:|")
    for k, v in sorted(fcount.items()):
        rep.append(f"| {k} | {v} |")
    rep.append(f"\nEditing rule used: **{edit_rule}**.\n")
    rep.append("`tracks_host_fraction` is confounded with the arm: the recurrent tumours carry more host reads "
               "(xengsort host fraction), so a real recurrence-gained allele can also track host fraction. It is only in "
               "the strict set.\n")
    rep.append("## Checks against something outside this computation\n")
    rep.append("| check | scope | expected | observed | result |\n|---|---|---|---|---|")
    for o in outside:
        rep.append("| " + " | ".join(fmt(x) for x in o) + " |")
    rep.append("\n## Notes\n")
    for n in notes:
        rep.append(f"- {n}")
    rep.append("\n## What this can and cannot show\n")
    rep.append("- It can show whether the tumours, the culture and the contralateral hemisphere carry the U-251MG donor "
               "lineage's TP53, PTEN, IDH1 and HLA class I states, from RNA, at depth.\n"
               "- It can show whether any site-count, arm-imbalance or Y/X difference between the arms stands out among the "
               "10 balanced relabellings. With 3 v 3 the best possible rank gives p = 0.10: nothing here reaches p < 0.05, "
               "by design.\n"
               "- It cannot separate LITT selection from founder bottlenecks: whether the arms and C2B came from the same "
               "culture batch or passage is not recorded.\n"
               "- It cannot see alleles in genes U251 does not express, and RNA allele fractions carry allele-specific "
               "expression and nonsense-mediated decay as well as copy number. Comparing each tumour with C2B at the same "
               "SNP cancels the constant part of that, not all of it.\n"
               "- xengsort's ambiguous and neither bins (0.75-2.49 % of pairs) were never saved; whether alleles are "
               "routed into them differently is not known.\n"
               "- N269B against IL69B is one animal and descriptive. N269B shares an i7 with IL66B and IL69B and an i5 with "
               "C2B and NL71B (prep_reads INDEX_HOPPING.txt); `private_allele_probe.tsv` shows whether private alleles of a "
               "library from another rat reach it.\n")
    with open(f"{O}/REPORT.md.part", "w") as f:
        f.write("\n".join(rep) + "\n")
    os.replace(f"{O}/REPORT.md.part", f"{O}/REPORT.md")
    log("final: REPORT.md written")
    return 0


# ------------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("checks")
    for k in ("--results", "--work", "--prep-checks", "--nextflow-exit", "--fasta", "--u", "--metadata", "--mpc-dir", "--out-dir"):
        s.add_argument(k, required=True)
    s.add_argument("--r2-reads", type=int, default=10000)
    s = sp.add_parser("sites")
    for k in ("--results", "--fasta", "--gtf", "--work"):
        s.add_argument(k, required=True)
    s = sp.add_parser("counts")
    for k in ("--fasta", "--work", "--mp-dir"):
        s.add_argument(k, required=True)
    s = sp.add_parser("select")
    s.add_argument("--work", required=True)
    s = sp.add_parser("readlevel")
    for k in ("--bam", "--tag", "--candidates", "--out-prefix"):
        s.add_argument(k, required=True)
    s.add_argument("--both-names", default="")
    s = sp.add_parser("yx")
    for k in ("--bam", "--tag", "--work", "--out"):
        s.add_argument(k, required=True)
    s = sp.add_parser("final")
    for k in ("--results", "--work", "--fasta", "--u", "--prep-checks", "--rl-dir", "--yx-dir", "--index-dir", "--out-dir",
              "--rat-fasta"):
        s.add_argument(k, required=True)
    for k in ("--rat-sam", "--lifted", "--rediportal", "--cytoband", "--invasion-dir"):
        s.add_argument(k, default="")
    a = ap.parse_args()
    return {"checks": cmd_checks, "sites": cmd_sites, "counts": cmd_counts, "select": cmd_select, "readlevel": cmd_readlevel,
            "yx": cmd_yx, "final": cmd_final}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
