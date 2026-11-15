#!/usr/bin/env python3
"""Checks for stage_refs that are easier to get right in Python than in awk. Standard library only.

Every subcommand prints TSV rows  check <TAB> expected <TAB> observed <TAB> result <TAB> outside_reference
to stdout (result is PASS, FAIL or INFO) and exits 1 if any row is FAIL.

  b64hex    <base64>                          -> hex (GCS x-goog-hash values)
  preflight --base DIR --ver N                rnafusion 4.1.3 exists_not_empty() emulation
  grch38    --fai FAI                         GRCh38 chromosome lengths against the published assembly
  faicmp    --a FAI --b FAI --label L [--primary] [--a-name S --b-name S]
  vcfhdr    --header TXT --fai FAI --label L  ##contig lengths of chr1-22,X,Y,M against the fai
  idxcmp    --a TSV --b TSV --label L         `bcftools index -s` record counts identical
  idxsub    --full TSV --cut TSV --label L    cut file holds every record of chr1-22,X,Y,M in the full file
  chain     --chain GZ --hfai FAI --rat TSV   UCSC hg38ToRn7 chain sizes against GENCODE and mRatBN7.2
  redisample --table GZ --out-regions F --out-expect F --n N
  refsample --out-regions F --out-expect F --label L --total N   (stdin: CHROM POS REF from bcftools query)
  redicmp   --expect F --faidx-out F [--label L --min 0.99]   expected bases against samtools faidx output
"""
import argparse, base64, gzip, os, sys

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
# GRCh38 (GCA_000001405.15) lengths as published by the Genome Reference Consortium / NCBI assembly report
GRCH38 = {"chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555, "chr5": 181538259,
          "chr6": 170805979, "chr7": 159345973, "chr8": 145138636, "chr9": 138394717, "chr10": 133797422,
          "chr11": 135086622, "chr12": 133275309, "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
          "chr16": 90338345, "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
          "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415, "chrM": 16569}

ROWS = []


def row(check, exp, obs, res, ref=""):
    ROWS.append((check, str(exp), str(obs), res, ref))


def done():
    for r in ROWS:
        print("\t".join(r))
    return 1 if any(r[3] == "FAIL" for r in ROWS) else 0


def read_fai(path):
    out = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            out[p[0]] = int(p[1])
    return out


# ------------------------------------------------------------------ rnafusion exists_not_empty emulation
def groovy_files(d):
    """Regular files in the order Groovy's File.eachFileRecurse(FileType.FILES) visits them: File.listFiles() is
    readdir order (os.scandir is readdir order too), directories recursed where they are met (pre-order)."""
    with os.scandir(d) as it:
        entries = [e.path for e in it]
    for p in entries:
        if os.path.isdir(p):          # java.io.File.isDirectory follows symlinks, so does os.path.isdir
            yield from groovy_files(p)
        else:
            yield p


def preflight(a):
    base, v = a.base.rstrip("/"), a.ver
    g = f"{base}/GRCh38/gencode_v{v}"
    dirs = {"starfusion_ref": f"{g}/starfusion/ctat_genome_lib_build_dir",
            "starindex_ref": f"{g}/star",
            "fusioncatcher_ref": f"{g}/fusioncatcher",
            "fusionreport_ref": f"{base}/GRCh38/fusion_report_db"}
    files = {"fasta": f"{g}/gencode/Homo_sapiens_GRCh38_{v}_dna_primary_assembly.fa",
             "fai": f"{g}/gencode/Homo_sapiens_GRCh38_{v}_dna_primary_assembly.fa.fai",
             "gtf": f"{g}/gencode/Homo_sapiens_GRCh38_{v}.gtf",
             "refflat": f"{g}/gencode/Homo_sapiens_GRCh38_{v}.gtf.refflat",
             "rrna_intervals": f"{g}/gencode/Homo_sapiens_GRCh38_{v}.gtf.interval_list",
             "hgnc_ref": f"{base}/GRCh38/hgnc/hgnc_complete_set.txt",
             "hgnc_date": f"{base}/GRCh38/hgnc/HGNC-DB-timestamp.txt",
             "arriba_ref_blacklist": f"{base}/GRCh38/arriba/2.5.0/blacklist_hg38_GRCh38_v2.5.0.tsv.gz",
             "arriba_ref_cytobands": f"{base}/GRCh38/arriba/2.5.0/cytobands_hg38_GRCh38_v2.5.0.tsv",
             "arriba_ref_known_fusions": f"{base}/GRCh38/arriba/2.5.0/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz",
             "arriba_ref_protein_domains": f"{base}/GRCh38/arriba/2.5.0/protein_domains_hg38_GRCh38_v2.5.0.gff3"}
    ref = "rnafusion 4.1.3 subworkflows/local/build_references/main.nf l.238-262 (conf/genomes.config paths)"
    host = os.uname().nodename
    for name, d in dirs.items():
        if not os.path.isdir(d):
            row(f"preflight:{name}", "directory with files", f"missing {d}", "FAIL", ref)
            continue
        fl = list(groovy_files(d))
        if not fl:
            row(f"preflight:{name}", ">=1 file", "no files (the pipeline would stop on a null file)", "FAIL", ref)
            continue
        zero = [p for p in fl if os.path.getsize(p) == 0]
        first, last = fl[0], fl[-1]
        # eachFileRecurse ignores the closure's `return`, so `first_file` ends as the LAST file visited: that is the
        # file the pipeline tests. The first file is checked too (the plan's reading of the code); both must be > 0.
        for tag, p in (("tested_last_visited", last), ("first_visited", first)):
            sz = os.path.getsize(p)
            row(f"preflight:{name}:{tag}", "size > 0", f"{os.path.relpath(p, d)} {sz} bytes (readdir order on {host})",
                "PASS" if sz > 0 else "FAIL", ref)
        row(f"preflight:{name}:zero_byte_files", "info", f"{len(zero)} of {len(fl)} files are 0 bytes"
            + (f" (e.g. {os.path.relpath(zero[0], d)})" if zero else ""), "INFO",
            "any of these landing first/last in readdir order on the head node would trigger a silent rebuild")
    for name, p in files.items():
        sz = os.path.getsize(p) if os.path.isfile(p) else -1
        row(f"preflight:{name}", "size > 0", f"{sz} bytes", "PASS" if sz > 0 else "FAIL", ref)
    gp = f"{g}/star/genomeParameters.txt"
    ov = "missing"
    if os.path.isfile(gp):
        for line in open(gp):
            p = line.split()
            if p and p[0] == "sjdbOverhang":
                ov = p[1]
    row("preflight:star_sjdbOverhang", "99 (rnafusion read_length 100)", ov, "PASS" if ov == "99" else "FAIL",
        "star/genomeParameters.txt")
    return done()


# ------------------------------------------------------------------ genome length checks
def grch38(a):
    fai = read_fai(a.fai)
    bad = [c for c, n in GRCH38.items() if fai.get(c) != n]
    row("grch38_chromosome_lengths", "chr1-22,X,Y,M = GRC GRCh38 lengths", "all 25 equal" if not bad else
        "differ: " + ",".join(f"{c}={fai.get(c)}" for c in bad[:5]), "PASS" if not bad else "FAIL",
        "GRCh38 assembly report (GCA_000001405.15)")
    row("fai_sequence_count", "194 (25 chromosomes + 42 unlocalised + 127 unplaced)", len(fai),
        "PASS" if len(fai) == 194 else "INFO", "GRCh38 primary assembly composition")
    return done()


def faicmp(a):
    missing = [p for p in (a.a, a.b) if not os.path.isfile(p)]
    if missing:
        row(a.label, "both files present", "missing " + " ".join(missing), "FAIL", a.ref)
        return done()
    fa, fb = read_fai(a.a), read_fai(a.b)
    keys = PRIMARY if a.primary else sorted(set(fa) | set(fb))
    bad = [k for k in keys if fa.get(k) != fb.get(k)]
    row(f"{a.label}", f"{len(keys)} sequences, names and lengths equal ({a.a_name} vs {a.b_name})",
        "all equal" if not bad else f"{len(bad)} differ, e.g. " + ",".join(f"{k}:{fa.get(k)}/{fb.get(k)}" for k in bad[:4]),
        "PASS" if not bad else "FAIL", a.ref)
    return done()


def vcfhdr(a):
    fai = read_fai(a.fai)
    hdr = {}
    for line in open(a.header):
        if line.startswith("##contig=<"):
            body = line.strip()[len("##contig=<"):-1]
            kv = dict(x.split("=", 1) for x in body.split(",") if "=" in x)
            if "ID" in kv and "length" in kv:
                hdr[kv["ID"]] = int(kv["length"])
    missing = [c for c in PRIMARY if c not in hdr]
    bad = [c for c in PRIMARY if c in hdr and hdr[c] != fai.get(c)]
    ok = not bad and not missing
    row(a.label, "chr1-22,X,Y,M ##contig lengths = GENCODE v44 fai", "all 25 equal" if ok else
        f"missing={','.join(missing[:5])} differ=" + ",".join(f"{c}:{hdr[c]}/{fai.get(c)}" for c in bad[:5]),
        "PASS" if ok else "FAIL", a.ref)
    return done()


def read_idx(path):
    out = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 3:
            out[p[0]] = int(p[2])
    return out


def idxcmp(a):
    x, y = read_idx(a.a), read_idx(a.b)
    ok = x == y and len(x) > 0
    row(a.label, "per-contig record counts identical before/after", f"{len(x)} contigs, {sum(x.values())} records -> "
        f"{len(y)} contigs, {sum(y.values())} records", "PASS" if ok else "FAIL", "bcftools index -s")
    return done()


def idxsub(a):
    full, cut = read_idx(a.full), read_idx(a.cut)
    want = {c: n for c, n in full.items() if c in PRIMARY}
    ok = want == cut
    dropped = sum(n for c, n in full.items() if c not in PRIMARY)
    row(a.label, "cut = every chr1-22,X,Y,M record of the source", f"source {sum(full.values())} records "
        f"({dropped} on other contigs dropped), cut {sum(cut.values())}", "PASS" if ok else "FAIL", "bcftools index -s")
    return done()


# ------------------------------------------------------------------ liftOver chain
def chain(a):
    hfai = read_fai(a.hfai)
    t, q = {}, {}
    with gzip.open(a.chain, "rt") as f:
        for line in f:
            if line.startswith("chain"):
                p = line.split()
                t[p[2]] = int(p[3])
                q[p[7]] = int(p[8])
    hs = [c for c in PRIMARY if c != "chrM"]
    tb = [c for c in hs if t.get(c) != hfai.get(c)]
    row("chain:hg38_side_sizes", "chr1-22,X,Y tSize = GENCODE v44 fai", "all 24 equal" if not tb else
        ",".join(f"{c}:{t.get(c)}/{hfai.get(c)}" for c in tb[:5]), "PASS" if not tb else "FAIL", "UCSC chain headers")
    rat = {}
    if a.rat and os.path.isfile(a.rat):
        for line in open(a.rat):
            p = line.split()
            if len(p) >= 2:
                rat[p[0]] = int(p[1])
    rc = [str(i) for i in range(1, 21)] + ["X", "Y"]
    if rat:
        rb = [c for c in rc if q.get("chr" + c) != rat.get(c)]
        row("chain:rn7_side_sizes", "chr1-20,X,Y qSize = mRatBN7.2 (Ensembl 110) lengths from the host STAR index",
            "all 22 equal" if not rb else ",".join(f"chr{c}:{q.get('chr' + c)}/{rat.get(c)}" for c in rb[:5]),
            "PASS" if not rb else "FAIL", "H/results_rnaseq/genome/index/star/chrNameLength.txt (host run's genome)")
    else:
        row("chain:rn7_side_sizes", "mRatBN7.2 lengths", "host STAR chrNameLength.txt not readable", "INFO", "")
    row("chain:blocks", "info", f"{len(t)} hg38 sequences, {len(q)} rn7 sequences in chain headers", "INFO", "")
    return done()


# ------------------------------------------------------------------ REDIportal
def redisample(a):
    n_lines = 0
    with gzip.open(a.table, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ci = {k: i for i, k in enumerate(header)}
        need = ["Region", "Position", "Ref", "Ed"]
        if any(k not in ci for k in need):
            row("rediportal:header", "Region,Position,Ref,Ed columns", ",".join(header[:6]), "FAIL", "REDIportal v3")
            return done()
        kept, ag = [], 0
        stride = max(1, a.stride)
        for line in f:
            n_lines += 1
            p = line.split("\t", 6)
            if (p[3], p[4]) in (("A", "G"), ("T", "C")):
                ag += 1
            if (n_lines - 1) % stride == 0 and len(kept) < a.n:
                kept.append((p[1], int(p[2]), p[3]))
    with open(a.out_regions, "w") as r, open(a.out_expect, "w") as e:
        for c, pos, base in kept:
            r.write(f"{c}:{pos}-{pos}\n")
            e.write(f"{c}\t{pos}\t{base}\n")
    row("rediportal:rows", "info", n_lines, "INFO", "REDIportal v3 TABLE1_hg38")
    share = ag / n_lines if n_lines else 0
    row("rediportal:A>G_or_T>C_share", ">= 0.99 (A-to-I editing on + or - strand)", f"{share:.4f}",
        "PASS" if share >= 0.99 else "FAIL", "definition of the database")
    return done()


def refsample(a):
    """Evenly spaced sample of `bcftools query -f '%CHROM\\t%POS\\t%REF\\n'` rows (stdin) whose REF is plain ACGT,
    written as samtools faidx regions plus the expected sequence."""
    step = max(1, a.total // a.n) if a.total else 1
    kept, seen = [], 0
    for line in sys.stdin:
        seen += 1
        if (seen - 1) % step:
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 3 or not p[2] or set(p[2].upper()) - set("ACGT"):
            continue
        kept.append((p[0], int(p[1]), p[2]))
        if len(kept) >= a.n:
            break
    for _ in sys.stdin:   # drain, so the producer does not die on a closed pipe
        pass
    with open(a.out_regions, "w") as r, open(a.out_expect, "w") as e:
        for c, pos, ref in kept:
            r.write(f"{c}:{pos}-{pos + len(ref) - 1}\n")
            e.write(f"{c}\t{pos}\t{ref}\n")
    row(f"{a.label}:sample", "info", f"{len(kept)} records with ACGT REF, every {step}th", "INFO", "")
    return done()


def redicmp(a):
    exp = [l.rstrip("\n").split("\t") for l in open(a.expect)]
    seqs, cur = [], None
    for line in open(a.faidx_out):
        line = line.strip()
        if line.startswith(">"):
            if cur is not None:
                seqs.append(cur)
            cur = ""
        else:
            cur += line
    if cur is not None:
        seqs.append(cur)
    n = min(len(exp), len(seqs))
    match = sum(1 for (c, pos, b), s in zip(exp, seqs) if s.upper() == b.upper())
    frac = match / n if n else 0
    row(a.label, f">= {a.min} of sampled sites", f"{match}/{n} = {frac:.4f}"
        + ("" if len(exp) == len(seqs) else f" (expected {len(exp)} sequences, got {len(seqs)})"),
        "PASS" if n > 0 and frac >= a.min and len(exp) == len(seqs) else "FAIL", a.ref)
    return done()


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("b64hex"); s.add_argument("v")
    s = sp.add_parser("preflight"); s.add_argument("--base", required=True); s.add_argument("--ver", default="46")
    s = sp.add_parser("grch38"); s.add_argument("--fai", required=True)
    s = sp.add_parser("faicmp")
    for k in ("--a", "--b", "--label"):
        s.add_argument(k, required=True)
    s.add_argument("--primary", action="store_true"); s.add_argument("--a-name", default="a")
    s.add_argument("--b-name", default="b"); s.add_argument("--ref", default="")
    s = sp.add_parser("vcfhdr")
    for k in ("--header", "--fai", "--label"):
        s.add_argument(k, required=True)
    s.add_argument("--ref", default="")
    s = sp.add_parser("idxcmp")
    for k in ("--a", "--b", "--label"):
        s.add_argument(k, required=True)
    s = sp.add_parser("idxsub")
    for k in ("--full", "--cut", "--label"):
        s.add_argument(k, required=True)
    s = sp.add_parser("chain")
    for k in ("--chain", "--hfai"):
        s.add_argument(k, required=True)
    s.add_argument("--rat", default="")
    s = sp.add_parser("redisample")
    for k in ("--table", "--out-regions", "--out-expect"):
        s.add_argument(k, required=True)
    s.add_argument("--n", type=int, default=5000); s.add_argument("--stride", type=int, default=3000)
    s = sp.add_parser("refsample")
    for k in ("--out-regions", "--out-expect", "--label"):
        s.add_argument(k, required=True)
    s.add_argument("--n", type=int, default=20000); s.add_argument("--total", type=int, default=0)
    s = sp.add_parser("redicmp")
    for k in ("--expect", "--faidx-out"):
        s.add_argument(k, required=True)
    s.add_argument("--label", default="rediportal:Ref_equals_GENCODE_v44_base")
    s.add_argument("--min", type=float, default=0.99)
    s.add_argument("--ref", default="REDIportal hg38 coordinates vs GENCODE v44 GRCh38 "
                                    "(a build or 0/1-based offset would give ~25 %)")
    a = ap.parse_args()
    if a.cmd == "b64hex":
        print(base64.b64decode(a.v).hex())
        return 0
    return {"preflight": preflight, "grch38": grch38, "faicmp": faicmp, "vcfhdr": vcfhdr, "idxcmp": idxcmp,
            "idxsub": idxsub, "chain": chain, "redisample": redisample, "refsample": refsample,
            "redicmp": redicmp}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
