#!/usr/bin/env python3
"""Audit of nf-core/rnafusion 4.1.3 pass 1 for the U251N LITT xenograft (standard library only).

Stages (each writes <audit>/.done_<stage>; `all` runs the unfinished ones in order):
  collect   parse each caller's OWN output (Arriba fusions.tsv, STAR-Fusion fusion_predictions/abridged,
            FusionCatcher final list, fusion-report lists); FusionInspector is read for information only
  extract   pull every supporting read pair (Arriba read_identifiers, STAR-Fusion JunctionReads/SpanningFrags) out of
            the input FASTQs; test each name against M/reads/both_names/<lib>.names.gz
  rat       realign the supporting pairs to the rat genome (host run's STAR index) with co-linear alignment only
  score     artefact rules (a)-(d), reportability of new fusions, CCLE recall and identity verdict, outside checks,
            whitelist_pass2.tsv, AUDIT_OK

Nothing here reads FusionInspector output as a caller: whitelisting puts a fusion into FusionInspector whether or not a
caller found it.
"""
import argparse, csv, gzip, json, os, re, subprocess, sys, time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

LIBS = ["C2B", "IL64B", "IL66B", "IL67B", "IL68B", "IL69B", "N168B", "N269B", "NL70B", "NL71B"]
CONTROLS = ["IL64B", "N168B"]            # human files 93.5 % and 90.4 % xengsort 'both' (rat); the artefact filter
PRIMARY = ["IL67B", "IL68B", "IL69B"]
RECURRENT = ["IL66B", "NL70B", "NL71B"]
TUMOURS = PRIMARY + RECURRENT
# xengsort bins (pairs), asserted equal by prep_reads (M/reads/checks/<lib>.tsv); human file = graft + both
XENG = {"C2B": (84519445, 11022904), "IL64B": (314037, 4546256), "IL66B": (25459283, 8298076),
        "IL67B": (31886295, 5512307), "IL68B": (32225887, 5441545), "IL69B": (49327614, 6013633),
        "N168B": (435428, 4086055), "N269B": (3823923, 4423374), "NL70B": (31293679, 5345807),
        "NL71B": (29210358, 5811628)}
# CCLE_Fusions_20181130.txt, U251MG rows with an annotated-splice STAR-Fusion call, each private to U251MG among the
# 1,017 lines of the file (checked 2026-09-27). MLF1IP is now CENPU. CCLE coordinates are hg19: match by gene name.
TOP4 = [("DHX38", "AMPD1"), ("CENPU", "ENPP6"), ("PSD3", "RBPMS"), ("MICALL2", "MAD1L1")]
CCLE8 = TOP4 + [("CNPY3", "SPACA1"), ("ZKSCAN7", "KIF15"), ("MRPS18B", "XPO5"), ("PCMTD1", "TCEANC")]
# expected read scale IF expression matched CCLE (poly(A) there, total RNA here, before deduplication): scale only
SCALE = {"DHX38--AMPD1": ("~316", "~112-183"), "CENPU--ENPP6": ("~48", "~17-28"),
         "PSD3--RBPMS": ("~36", "~13-21"), "MICALL2--MAD1L1": ("~20", "~7-11")}
N269B_SCALE = "about 13-27 DHX38--AMPD1 reads (scale only)"
ADJ_BP = 200_000

ROWS = []          # checks: check, expected, observed, result, outside_reference


def row(check, exp, obs, res, ref=""):
    ROWS.append((check, str(exp), str(obs), res, ref))


def sname(lib):
    return lib + "_humclip"


def fkey(g5, g3):
    return f"{g5}--{g3}"


def log(msg):
    print(time.strftime("%Y-%m-%dT%H:%M:%S"), msg, flush=True)


def read_table(path):
    """Tab-separated file whose first line is the header (a leading '#' is dropped). Returns header, list of dicts."""
    with open(path, newline="") as f:
        header = f.readline().rstrip("\n").lstrip("#").split("\t")
        rows = []
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, p + [""] * (len(header) - len(p)))))
    return header, rows


def col(header, prefix):
    for h in header:
        if h.startswith(prefix):
            return h
    return None


def norm_read(n):
    n = n.strip()
    if n.startswith("@"):
        n = n[1:]
    n = n.split()[0] if n else n
    if n.endswith("/1") or n.endswith("/2"):
        n = n[:-2]
    return n


def to_int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def to_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


# ------------------------------------------------------------------ gene symbols
def load_hgnc(path):
    """Current symbols and an unambiguous previous-symbol -> current map (MLF1IP -> CENPU)."""
    cur, prev = set(), defaultdict(set)
    if not os.path.isfile(path):
        return cur, {}
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        r = csv.DictReader(f, delimiter="\t")
        for d in r:
            if d.get("status", "Approved") != "Approved":
                continue
            s = d["symbol"].strip()
            cur.add(s)
            for p in (d.get("prev_symbol") or "").replace('"', "").split("|"):
                p = p.strip()
                if p:
                    prev[p].add(s)
    pm = {p: next(iter(s)) for p, s in prev.items() if len(s) == 1 and p not in cur}
    return cur, pm


class Sym:
    def __init__(self, hgnc_path):
        self.cur, self.prev = load_hgnc(hgnc_path)

    def __call__(self, g):
        g = g.strip()
        if "^" in g:                   # STAR-Fusion LeftGene "SYMBOL^ENSG..."
            g = g.split("^")[0]
        return self.prev.get(g, g) if g not in self.cur else g


def intergenic(g):
    return ("," in g) or ("(" in g) or g in (".", "")


# ------------------------------------------------------------------ stage: collect
def collect(a):
    R, A = a.results, a.audit
    sym = Sym(a.hgnc)
    calls = []                          # one dict per caller x sample x fusion key
    reads = defaultdict(lambda: defaultdict(set))   # lib -> key -> normalised read names (Arriba + STAR-Fusion)
    arriba_idcheck = [0, 0]             # rows whose read_identifiers count equals split1+split2+discordant, rows
    missing = []
    fr_filtered = {}
    fi_info = []
    arriba_discarded = []
    for lib in LIBS:
        s = sname(lib)
        paths = {"arriba": f"{R}/arriba/{s}.arriba.fusions.tsv",
                 "arriba_discarded": f"{R}/arriba/{s}.arriba.fusions.discarded.tsv",
                 "sf_full": f"{R}/starfusion/{s}.starfusion.fusion_predictions.tsv",
                 "sf_abridged": f"{R}/starfusion/{s}.starfusion.abridged.tsv",
                 "fc": f"{R}/fusioncatcher/{s}.fusion-genes.txt",
                 "fr_filtered": f"{R}/fusionreport/{s}/{s}.fusionreport_filtered.tsv",
                 "star_log": f"{R}/star/{s}.Log.final.out"}
        for k, p in paths.items():
            if not os.path.isfile(p):
                missing.append(f"{s}:{k}")
        # ---- Arriba (fusions.tsv only counts as a call; the discarded file is read for the markers)
        if os.path.isfile(paths["arriba"]):
            hdr, rows = read_table(paths["arriba"])
            agg = defaultdict(list)
            for d in rows:
                g5r, g3r = d.get("gene1", ""), d.get("gene2", "")
                g5 = g5r if intergenic(g5r) else sym(g5r)
                g3 = g3r if intergenic(g3r) else sym(g3r)
                fr = to_int(d.get("split_reads1")) + to_int(d.get("split_reads2")) + to_int(d.get("discordant_mates"))
                dup = sum(int(x) for x in re.findall(r"duplicates\((\d+)\)", d.get("filters", "")))
                ids = [norm_read(x) for x in d.get("read_identifiers", "").split(",") if x and x != "."]
                arriba_idcheck[1] += 1
                arriba_idcheck[0] += int(len(set(ids)) == fr)
                splice = "splice-site" in d.get("site1", "") and "splice-site" in d.get("site2", "")
                agg[(g5, g3)].append(dict(frags=fr, dups=dup, splice=splice, conf=d.get("confidence", ""),
                                          type=d.get("type", ""), bp=f"{d.get('breakpoint1')}|{d.get('breakpoint2')}",
                                          site=f"{d.get('site1')}|{d.get('site2')}", ig=intergenic(g5r) or intergenic(g3r)))
                if not intergenic(g5r) and not intergenic(g3r):
                    reads[lib][fkey(g5, g3)].update(ids)
            for (g5, g3), rr in agg.items():
                best = max(rr, key=lambda x: x["frags"])
                conf_rank = {"high": 3, "medium": 2, "low": 1}
                calls.append(dict(lib=lib, caller="arriba", g5=g5, g3=g3, key=fkey(g5, g3), rows=len(rr),
                                  frags=sum(x["frags"] for x in rr),
                                  frags_splice=sum(x["frags"] for x in rr if x["splice"]),
                                  best_splice=best["splice"], dups=sum(x["dups"] for x in rr),
                                  conf=max((x["conf"] for x in rr), key=lambda c: conf_rank.get(c, 0)),
                                  detail=f"type={best['type']};site={best['site']};bp={best['bp']}",
                                  intergenic=any(x["ig"] for x in rr)))
        if os.path.isfile(paths["arriba_discarded"]):
            hdr, rows = read_table(paths["arriba_discarded"])
            for d in rows:
                g5, g3 = sym(d.get("gene1", "")), sym(d.get("gene2", ""))
                if (g5, g3) in CCLE8 or (g3, g5) in CCLE8:
                    arriba_discarded.append(dict(lib=lib, key=fkey(g5, g3),
                                                 frags=to_int(d.get("split_reads1")) + to_int(d.get("split_reads2"))
                                                 + to_int(d.get("discordant_mates")),
                                                 filters=d.get("filters", ""), site=f"{d.get('site1')}|{d.get('site2')}",
                                                 conf=d.get("confidence", "")))
        # ---- STAR-Fusion: calls from the abridged file, read names from the full predictions file
        if os.path.isfile(paths["sf_abridged"]):
            hdr, rows = read_table(paths["sf_abridged"])
            agg = defaultdict(list)
            for d in rows:
                g5, g3 = sym(d.get("LeftGene", "")), sym(d.get("RightGene", ""))
                agg[(g5, g3)].append(dict(frags=to_int(d.get("JunctionReadCount")) + to_int(d.get("SpanningFragCount")),
                                          splice=d.get("SpliceType", "") == "ONLY_REF_SPLICE",
                                          stype=d.get("SpliceType", ""), ffpm=to_float(d.get("FFPM")),
                                          bp=f"{d.get('LeftBreakpoint')}|{d.get('RightBreakpoint')}"))
            for (g5, g3), rr in agg.items():
                best = max(rr, key=lambda x: x["frags"])
                calls.append(dict(lib=lib, caller="starfusion", g5=g5, g3=g3, key=fkey(g5, g3), rows=len(rr),
                                  frags=sum(x["frags"] for x in rr),
                                  frags_splice=sum(x["frags"] for x in rr if x["splice"]),
                                  best_splice=best["splice"], dups="", conf="",
                                  detail=f"SpliceType={best['stype']};FFPM={sum(x['ffpm'] for x in rr):.4g};bp={best['bp']}",
                                  intergenic=False))
        if os.path.isfile(paths["sf_full"]):
            hdr, rows = read_table(paths["sf_full"])
            for d in rows:
                g5, g3 = sym(d.get("LeftGene", "")), sym(d.get("RightGene", ""))
                for c in ("JunctionReads", "SpanningFrags"):
                    reads[lib][fkey(g5, g3)].update(norm_read(x) for x in d.get(c, "").split(",") if x and x != ".")
        # ---- FusionCatcher final list (counts are not deduplicated: information only, never the >= 3 rule)
        if os.path.isfile(paths["fc"]):
            hdr, rows = read_table(paths["fc"])
            c5, c3 = col(hdr, "Gene_1_symbol"), col(hdr, "Gene_2_symbol")
            cp, cu = col(hdr, "Spanning_pairs"), col(hdr, "Spanning_unique_reads")
            cd, ce = col(hdr, "Fusion_description"), col(hdr, "Predicted_effect")
            agg = defaultdict(list)
            for d in rows:
                if not c5 or not c3:
                    break
                g5, g3 = sym(d.get(c5, "")), sym(d.get(c3, ""))
                agg[(g5, g3)].append(dict(frags=to_int(d.get(cp)) + to_int(d.get(cu)), desc=d.get(cd, ""),
                                          eff=d.get(ce, "")))
            for (g5, g3), rr in agg.items():
                best = max(rr, key=lambda x: x["frags"])
                calls.append(dict(lib=lib, caller="fusioncatcher", g5=g5, g3=g3, key=fkey(g5, g3), rows=len(rr),
                                  frags=sum(x["frags"] for x in rr), frags_splice="", best_splice="", dups="",
                                  conf="", detail=f"desc={best['desc']};effect={best['eff']}", intergenic=False))
        # ---- fusion-report filtered list (tools >= tools_cutoff 2): names as the tools gave them
        if os.path.isfile(paths["fr_filtered"]):
            fr_filtered[lib] = [l.strip() for l in open(paths["fr_filtered"]) if l.strip()]
        # ---- FusionInspector (pass 1), information only
        fi = f"{R}/fusioninspector/{s}/{s}.FusionInspector.fusions.abridged.tsv"
        if os.path.isfile(fi):
            hdr, rows = read_table(fi)
            for d in rows:
                fi_info.append(dict(lib=lib, key=d.get("FusionName", ""),
                                    frags=to_int(d.get("JunctionReadCount")) + to_int(d.get("SpanningFragCount")),
                                    ffpm=to_float(d.get("FFPM"))))
    os.makedirs(f"{A}/names", exist_ok=True)
    with open(f"{A}/calls.tsv", "w") as f:
        cols = ["lib", "caller", "g5", "g3", "key", "rows", "frags", "frags_splice", "best_splice", "dups", "conf",
                "intergenic", "detail"]
        f.write("\t".join(cols) + "\n")
        for c in calls:
            f.write("\t".join(str(c[k]) for k in cols) + "\n")
    for lib in LIBS:
        with open(f"{A}/names/{lib}.support.tsv", "w") as f:
            for k, names in reads[lib].items():
                for n in sorted(names):
                    f.write(f"{k}\t{n}\n")
    json.dump(dict(missing=missing, arriba_idcheck=arriba_idcheck, fr_filtered=fr_filtered, fi_info=fi_info,
                   arriba_discarded=arriba_discarded), open(f"{A}/collect.json", "w"), indent=1)
    log(f"collect: {len(calls)} caller calls; missing files: {missing or 'none'}")


# ------------------------------------------------------------------ stage: extract (one process per library)
def _scan(args):
    lib, names_file, r1, r2, both, outdir = args
    need = set()
    for line in open(names_file):
        p = line.rstrip("\n").split("\t")
        if len(p) == 2:
            need.add(p[1].encode())
    res = dict(lib=lib, need=len(need), found_r1=0, found_r2=0, both=0, both_names_total=0, error="")
    if not need:
        open(f"{outdir}/{lib}_R1.fq", "w").close(); open(f"{outdir}/{lib}_R2.fq", "w").close()
        open(f"{outdir}/{lib}.both_hits.txt", "w").close()
        return res
    tag = lib.encode() + b"|"
    for mate, path in (("r1", r1), ("r2", r2)):
        out = open(f"{outdir}/{lib}_{mate.upper()}.fq.part", "wb")
        p = subprocess.Popen(["zcat", path], stdout=subprocess.PIPE, bufsize=1 << 22)
        keep, n, found = 0, 0, 0
        for line in p.stdout:
            if n % 4 == 0:
                nm = line[1:].split(None, 1)[0]
                if nm.endswith(b"/1") or nm.endswith(b"/2"):
                    nm = nm[:-2]
                if nm in need:
                    keep = 4
                    found += 1
                    out.write(b"@" + tag + nm + b"\n")
                    n += 1
                    keep -= 1
                    continue
            if keep:
                out.write(line)
                keep -= 1
            n += 1
        p.wait()
        out.close()
        if p.returncode != 0:
            res["error"] += f"zcat {path} exit {p.returncode}; "
        if n % 4:
            res["error"] += f"{path}: {n} lines is not a multiple of 4; "
        os.replace(f"{outdir}/{lib}_{mate.upper()}.fq.part", f"{outdir}/{lib}_{mate.upper()}.fq")
        res[f"found_{mate}"] = found
    hits = []
    with gzip.open(both, "rb") as f:
        for line in f:
            res["both_names_total"] += 1
            nm = line.strip()
            if nm in need:
                hits.append(nm.decode())
    res["both"] = len(hits)
    with open(f"{outdir}/{lib}.both_hits.txt", "w") as f:
        f.write("\n".join(hits) + ("\n" if hits else ""))
    return res


def extract(a):
    A = a.audit
    os.makedirs(f"{A}/reads", exist_ok=True)
    jobs = []
    for lib in LIBS:
        jobs.append((lib, f"{A}/names/{lib}.support.tsv", f"{a.u}/ANALYSIS/sorted_fastqs/{lib}_human_R1.fq.gz",
                     f"{a.reads}/human_clip/{lib}_human_R2.fq.gz", f"{a.reads}/both_names/{lib}.names.gz",
                     f"{A}/reads"))
    # largest first so the long C2B scan starts at once
    jobs.sort(key=lambda j: -os.path.getsize(j[2]) if os.path.isfile(j[2]) else 0)
    out = []
    with ProcessPoolExecutor(max_workers=a.threads) as ex:
        for r in ex.map(_scan, jobs):
            log(f"extract {r['lib']}: need {r['need']}, R1 {r['found_r1']}, R2 {r['found_r2']}, both {r['both']}"
                f" (of {r['both_names_total']} both names) {r['error']}")
            out.append(r)
    json.dump(out, open(f"{A}/extract.json", "w"), indent=1)


# ------------------------------------------------------------------ stage: rat (co-linear realignment to mRatBN7.2)
def rat(a):
    A = a.audit
    os.makedirs(f"{A}/rat", exist_ok=True)
    ext = {r["lib"]: r for r in json.load(open(f"{A}/extract.json"))}
    # a library whose R1 and R2 extractions differ would give STAR unequal mate files: leave it out (rule c is then
    # 'not assessable' for it, never 0 %)
    use = [lib for lib in LIBS if ext.get(lib, {}).get("found_r1", 0) == ext.get(lib, {}).get("found_r2", -1)
           and not ext.get(lib, {}).get("error")]
    r1s = [f"{A}/reads/{lib}_R1.fq" for lib in use]
    r2s = [f"{A}/reads/{lib}_R2.fq" for lib in use]
    n = sum(os.path.getsize(p) for p in r1s if os.path.isfile(p))
    info = dict(img=a.star_img, index=a.rat_index, ran=False, error="", libraries=use,
                left_out=[lib for lib in LIBS if lib not in use])
    if not a.star_img or not os.path.isfile(a.star_img):
        info["error"] = f"STAR image not found ({a.star_img!r})"
        json.dump(info, open(f"{A}/rat/rat.json", "w"), indent=1)
        return
    if n == 0:
        info["error"] = "no supporting reads extracted"
        json.dump(info, open(f"{A}/rat/rat.json", "w"), indent=1)
        return
    for mate, lst in (("R1", r1s), ("R2", r2s)):
        with open(f"{A}/rat/all_{mate}.fq", "wb") as o:
            for p in lst:
                if os.path.isfile(p):
                    with open(p, "rb") as i:
                        o.write(i.read())
    sx = ["singularity", "exec", "-B", "/rs/rs_grp_oschome/go2432", "-B", f"{a.tmp}:/tmp", a.star_img]
    try:
        info["star_version"] = subprocess.run(sx + ["STAR", "--version"], capture_output=True, text=True).stdout.strip()
    except Exception as e:                      # noqa: BLE001
        info["star_version"] = f"error {e}"
    cmd = sx + ["STAR", "--runThreadN", str(a.threads), "--genomeDir", a.rat_index, "--genomeLoad", "NoSharedMemory",
                "--readFilesIn", f"{A}/rat/all_R1.fq", f"{A}/rat/all_R2.fq", "--outFileNamePrefix", f"{A}/rat/",
                "--outSAMtype", "SAM", "--outSAMunmapped", "Within", "--outSAMattributes", "NH", "HI", "AS", "nM", "NM",
                "--outFilterMultimapNmax", "20", "--outFilterMismatchNoverReadLmax", "0.04",
                "--outFilterScoreMinOverLread", "0.9", "--outFilterMatchNminOverLread", "0.9",
                "--alignIntronMax", "1000000", "--alignMatesGapMax", "1000000", "--chimSegmentMin", "0"]
    info["cmd"] = " ".join(cmd)
    log("rat: " + info["cmd"])
    p = subprocess.run(cmd, capture_output=True, text=True)
    info["ran"] = p.returncode == 0
    info["returncode"] = p.returncode
    info["stderr_tail"] = p.stderr[-2000:]
    if p.returncode != 0:
        info["error"] = f"STAR exit {p.returncode}"
    json.dump(info, open(f"{A}/rat/rat.json", "w"), indent=1)
    for mate in ("R1", "R2"):
        try:
            os.remove(f"{A}/rat/all_{mate}.fq")
        except OSError:
            pass


def rat_contiguous(A):
    """lib -> set of read names whose pair aligned co-linearly to rat (>= 90 % of the pair, <= 4 % mismatches)."""
    out, total = defaultdict(set), defaultdict(set)
    sam = f"{A}/rat/Aligned.out.sam"
    if not os.path.isfile(sam):
        return None, None
    with open(sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            p = line.split("\t", 3)
            q, flag = p[0], int(p[1])
            lib, _, nm = q.partition("|")
            total[lib].add(nm)
            if not flag & 4 and not flag & 0x100 and not flag & 0x800:
                out[lib].add(nm)
    return out, total


# ------------------------------------------------------------------ genome annotation helpers
def gene_spans_gtf(path, key="gene_name", feature="gene"):
    """name -> list of (chrom, start, end, strand) from 'gene' lines (every line when feature is None)."""
    spans = defaultdict(list)
    tmp = {}
    rk = re.compile(key + r' "([^"]+)"')
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if len(p) < 9 or (feature and p[2] != feature):
                continue
            m = rk.search(p[8])
            if not m:
                continue
            g = m.group(1)
            if feature:
                spans[g].append((p[0], int(p[3]), int(p[4]), p[6]))
            else:
                c, s, e, st = tmp.get(g, (p[0], int(p[3]), int(p[4]), p[6]))
                tmp[g] = (c, min(s, int(p[3])), max(e, int(p[4])), st)
    if not feature:
        for g, v in tmp.items():
            spans[g].append(v)
    return spans


def gap(x, y):
    if x[0] != y[0]:
        return None
    return max(0, max(x[1], y[1]) - min(x[2], y[2]))


def adjacency(spans_a, spans_b):
    """smallest same-chromosome gap between two genes, and whether any such pair is on the same strand."""
    best = None
    for x in spans_a:
        for y in spans_b:
            g = gap(x, y)
            if g is None:
                continue
            if best is None or g < best[0]:
                best = (g, x[3] == y[3])
    return best


def load_orthologs(path):
    by_sym, by_id = {}, {}
    with open(path, newline="") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            by_sym.setdefault(d["human_symbol"], d["rat_ensembl_id"])
            by_id.setdefault(d["human_ensembl_id"].split(".")[0], d["rat_ensembl_id"])
    return by_sym, by_id


# ------------------------------------------------------------------ outside checks
def published_table(path, key="Sample"):
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, newline="") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            out[d[key]] = d
    return out


def star_log(path):
    out = {}
    if not os.path.isfile(path):
        return out
    for line in open(path):
        if "|" in line:
            k, v = line.split("|", 1)
            out[k.strip()] = v.strip()
    return out


def picard_dup(path):
    if not os.path.isfile(path):
        return None
    lines = open(path).read().splitlines()
    for i, l in enumerate(lines):
        if l.startswith("LIBRARY") and i + 1 < len(lines):
            h, v = l.split("\t"), lines[i + 1].split("\t")
            d = dict(zip(h, v))
            return to_float(d.get("PERCENT_DUPLICATION"))
    return None


# ------------------------------------------------------------------ stage: score
def score(a):
    R, A = a.results, a.audit
    col_ = json.load(open(f"{A}/collect.json"))
    ext = {r["lib"]: r for r in json.load(open(f"{A}/extract.json"))}
    ratinfo = json.load(open(f"{A}/rat/rat.json")) if os.path.isfile(f"{A}/rat/rat.json") else {"ran": False}
    ratc, rattot = rat_contiguous(A) if ratinfo.get("ran") else (None, None)
    hdr, calls = read_table(f"{A}/calls.tsv")
    for c in calls:
        c["frags"] = to_int(c["frags"])
        c["frags_splice"] = to_int(c["frags_splice"]) if c["frags_splice"] != "" else None
        c["best_splice"] = c["best_splice"] == "True"
        c["dups"] = to_int(c["dups"]) if c["dups"] != "" else 0
    support = defaultdict(lambda: defaultdict(set))
    for lib in LIBS:
        p = f"{A}/names/{lib}.support.tsv"
        if os.path.isfile(p):
            for line in open(p):
                k, n = line.rstrip("\n").split("\t")
                support[lib][k].add(n)
    both_hits = {lib: set(l.strip() for l in open(f"{A}/reads/{lib}.both_hits.txt") if l.strip())
                 if os.path.isfile(f"{A}/reads/{lib}.both_hits.txt") else set() for lib in LIBS}

    # ---- hard checks
    miss = col_["missing"]
    miss_hard = [m for m in miss if not m.endswith(":arriba_discarded")]
    row("outputs_present", "every caller output, fusion-report list and STAR log for 10 libraries",
        "all present" if not miss_hard else "missing: " + ",".join(miss_hard), "PASS" if not miss_hard else "FAIL",
        "rnafusion 4.1.3 publishDir layout (conf/modules.config)")
    star_ok = True
    starinfo = {}
    for lib in LIBS:
        sl = star_log(f"{R}/star/{sname(lib)}.Log.final.out")
        starinfo[lib] = sl
        n_in = to_int(sl.get("Number of input reads"))
        want = sum(XENG[lib])
        ok = n_in == want
        star_ok &= ok
        row(f"star_input_pairs:{lib}", f"{want} (xengsort graft + both = human file pairs)", n_in,
            "PASS" if ok else "FAIL", "xengsort classify logs (U/ANALYSIS/xengsort_out), asserted by prep_reads")
    # name format: supporting names must be found in the FASTQs they came from
    name_ok = {}
    for lib in LIBS:
        e = ext.get(lib, {})
        need, f1, f2 = e.get("need", 0), e.get("found_r1", 0), e.get("found_r2", 0)
        ok = need == 0 or (f1 >= 0.95 * need and f2 == f1 and not e.get("error"))
        name_ok[lib] = ok
        row(f"support_names_found:{lib}", ">= 95 % of supporting read names found in R1, R2 found = R1 found",
            f"need {need}, R1 {f1}, R2 {f2} {e.get('error', '')}", "PASS" if ok else "WARN",
            "the input FASTQ headers (a name-format mismatch would make rules b and c read 0 %)")
        nb = e.get("both_names_total", -1)
        if need:     # the names file is only streamed when there is something to look up
            okb = nb == XENG[lib][1]
            name_ok[lib] = name_ok[lib] and okb
            row(f"both_names_count:{lib}", f"{XENG[lib][1]} (xengsort 'both' pairs)", nb, "PASS" if okb else "WARN",
                "xengsort classify log; a short names file would understate rule b")
    ai = col_["arriba_idcheck"]
    row("arriba_read_identifiers_parse", "read_identifiers count = split_reads1 + split_reads2 + discordant_mates",
        f"{ai[0]} of {ai[1]} rows", "PASS" if ai[1] == 0 or ai[0] >= 0.95 * ai[1] else "WARN",
        "Arriba 2.5.0 output-file definition")
    row("rat_realignment", "STAR realignment of supporting pairs to mRatBN7.2 ran",
        f"ran={ratinfo.get('ran')} {ratinfo.get('star_version', '')} {ratinfo.get('error', '')}"
        f" left out: {','.join(ratinfo.get('left_out', [])) or 'none'}",
        "PASS" if ratinfo.get("ran") else "WARN", "host run's STAR index (H/results_rnaseq/genome/index/star)")
    rat_ok = {}
    for lib in LIBS:
        need = ext.get(lib, {}).get("found_r1", 0)
        seen = len(rattot.get(lib, set())) if rattot is not None else -1
        ok = ratc is not None and lib in ratinfo.get("libraries", []) and seen == need
        rat_ok[lib] = ok
        if need:
            row(f"rat_sam_accounts_for_every_pair:{lib}", f"{need} extracted pairs in the SAM (unmapped kept)",
                seen, "PASS" if ok else "WARN", "STAR --outSAMunmapped Within writes every input pair")
    # published human run (outside this pipeline): STAR unique-mapping % and Picard duplication
    pub_star = published_table(f"{a.u}/ANALYSIS/results_human_final/multiqc/star_salmon/multiqc_report_data/multiqc_star.txt")
    pub_dup = published_table(f"{a.u}/ANALYSIS/results_human_final/multiqc/star_salmon/multiqc_report_data/multiqc_picard_dups.txt")
    for lib in LIBS:
        u_now = to_float(starinfo[lib].get("Uniquely mapped reads %", "").rstrip("%"))
        u_pub = to_float(pub_star.get(lib, {}).get("uniquely_mapped_percent"))
        res = "INFO" if not pub_star.get(lib) else ("PASS" if abs(u_now - u_pub) <= 8 else "WARN")
        row(f"star_unique_pct_vs_published:{lib}", f"within 8 points of the published run ({u_pub:.2f} %)",
            f"{u_now:.2f} %", res, "U/ANALYSIS/results_human_final multiqc_star.txt (GENCODE v44, trimmed, unclipped R2)")
        d_now = picard_dup(f"{R}/picard/{sname(lib)}.bam.metrics")
        d_pub = to_float(pub_dup.get(lib, {}).get("PERCENT_DUPLICATION"))
        res = "INFO" if d_now is None or not pub_dup.get(lib) else ("PASS" if abs(d_now - d_pub) <= 0.08 else "WARN")
        row(f"duplication_vs_published:{lib}", f"within 0.08 of the published Picard rate ({d_pub:.3f})",
            "missing" if d_now is None else f"{d_now:.3f}", res,
            "U/ANALYSIS/results_human_final multiqc_picard_dups.txt")

    # ---- per sample x key: callers and flags
    bykey = defaultdict(lambda: defaultdict(dict))      # lib -> key -> caller -> call
    for c in calls:
        bykey[c["lib"]][c["key"]][c["caller"]] = c
    control_keys = set()
    for lib in CONTROLS:
        for k in bykey[lib]:
            g5, g3 = k.split("--", 1)
            control_keys.add(frozenset((g5, g3)))
    # rat and human gene spans for rule (d)
    log("score: reading gene spans (rat filtered GTF, GENCODE v46 GTF)")
    rat_spans = gene_spans_gtf(a.rat_gtf, key="gene_id", feature=None)
    hum_spans = gene_spans_gtf(a.gtf, key="gene_name", feature="gene")
    orth_sym, orth_id = load_orthologs(a.orthologs)
    ctat_genes = set()
    for line in open(a.ctat_spans):
        for x in line.rstrip("\n").split("\t"):
            ctat_genes.add(x)

    def rule_d(g5, g3):
        r5, r3 = orth_sym.get(g5), orth_sym.get(g3)
        if not r5 or not r3:
            return None, f"no 1:1 rat ortholog for {'/'.join(g for g, r in ((g5, r5), (g3, r3)) if not r)}"
        ad = adjacency(rat_spans.get(r5, []), rat_spans.get(r3, []))
        if ad is None:
            return False, f"rat orthologs {r5},{r3} on different chromosomes"
        flag = ad[0] <= ADJ_BP and ad[1]
        return flag, f"rat {r5},{r3} gap {ad[0]} bp, same strand {ad[1]}"

    def human_adj(g5, g3):
        ad = adjacency(hum_spans.get(g5, []), hum_spans.get(g3, []))
        if ad is None:
            return "different chromosomes or not in v46 GTF"
        return f"gap {ad[0]} bp, same strand {ad[1]}" + (" (read-through distance)" if ad[0] <= ADJ_BP and ad[1] else "")

    scored = []
    dcache = {}
    for lib in LIBS:
        for k, cs in bykey[lib].items():
            g5, g3 = k.split("--", 1)
            ig = any(str(c.get("intergenic")) == "True" for c in cs.values())
            callers = sorted(cs)
            names = support[lib].get(k, set())
            nb = len(names & both_hits[lib])
            bfrac = nb / len(names) if names else None
            if rat_ok[lib] and names:
                rfrac = len(names & ratc.get(lib, set())) / len(names)
            else:
                rfrac = None
            if k not in dcache:
                dcache[k] = (None, "intergenic") if ig else rule_d(g5, g3)
            dflag, dnote = dcache[k]
            marker = (g5, g3) in CCLE8
            recip = (g3, g5) in CCLE8
            a_flag = frozenset((g5, g3)) in control_keys
            b_ass = bool(names) and name_ok[lib]
            c_ass = b_ass and rat_ok[lib]
            b_flag = bool(b_ass and bfrac is not None and bfrac > 0.5)
            c_flag = bool(c_ass and rfrac is not None and rfrac > 0.5)
            dedup_splice = max([cs[c]["frags_splice"] or 0 for c in ("arriba", "starfusion") if c in cs] or [0])
            splice_ok = all(cs[c]["best_splice"] for c in ("arriba", "starfusion") if c in cs) and \
                any(c in cs for c in ("arriba", "starfusion"))
            reasons = []
            if marker or recip:
                reasons.append("CCLE marker" if marker else "reciprocal of a CCLE marker")
            if lib in CONTROLS:
                reasons.append("control library")
            if ig:
                reasons.append("intergenic breakpoint")
            if len(callers) < 2:
                reasons.append(f"{len(callers)} of 3 callers")
            if not any(c in cs for c in ("arriba", "starfusion")):
                reasons.append("no Arriba or STAR-Fusion call: deduplicated count and splice sites unknown")
            else:
                if dedup_splice < 3:
                    reasons.append(f"{dedup_splice} deduplicated fragments at annotated splice sites (< 3)")
                if not splice_ok:
                    reasons.append("an Arriba/STAR-Fusion best breakpoint is not on annotated splice sites")
            if a_flag:
                reasons.append("(a) also called in IL64B or N168B")
            if not b_ass:
                reasons.append("(b) not assessable (no read names or name format mismatch)")
            elif b_flag:
                reasons.append(f"(b) {bfrac:.0%} of supporting reads in xengsort 'both'")
            if not c_ass:
                reasons.append("(c) not assessable")
            elif c_flag:
                reasons.append(f"(c) {rfrac:.0%} of supporting reads align co-linearly to rat")
            if dflag is None and not ig:
                reasons.append(f"(d) not assessable: {dnote}")
            elif dflag:
                reasons.append(f"(d) {dnote}")
            reportable = not reasons
            scored.append(dict(lib=lib, key=k, callers=",".join(callers), n_callers=len(callers),
                               arriba_frags=cs.get("arriba", {}).get("frags", ""),
                               arriba_frags_splice=cs.get("arriba", {}).get("frags_splice", ""),
                               arriba_conf=cs.get("arriba", {}).get("conf", ""),
                               sf_frags=cs.get("starfusion", {}).get("frags", ""),
                               sf_frags_splice=cs.get("starfusion", {}).get("frags_splice", ""),
                               fc_frags_raw=cs.get("fusioncatcher", {}).get("frags", ""),
                               n_support_reads=len(names), both_frac="" if bfrac is None else f"{bfrac:.3f}",
                               rat_frac="" if rfrac is None else f"{rfrac:.3f}",
                               a=int(a_flag), b=int(b_flag), c=int(c_flag), d="" if dflag is None else int(bool(dflag)),
                               d_note=dnote, human_adjacency=human_adj(g5, g3) if not ig else "",
                               marker="CCLE" if marker else ("CCLE_reciprocal" if recip else ""),
                               reportable_new=int(reportable), reasons="; ".join(reasons)))
    cols = list(scored[0].keys()) if scored else ["lib", "key"]
    with open(f"{A}/scored_calls.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for s in sorted(scored, key=lambda x: (x["lib"], -x["n_callers"], x["key"])):
            f.write("\t".join(str(s[c]) for c in cols) + "\n")
    rep = [s for s in scored if s["reportable_new"]]
    with open(f"{A}/reportable_new.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for s in rep:
            f.write("\t".join(str(s[c]) for c in cols) + "\n")

    # ---- CCLE recall, per sample x marker x caller (each caller's own output)
    disc = defaultdict(list)
    for d in col_["arriba_discarded"]:
        disc[(d["lib"], d["key"])].append(d)
    fi = defaultdict(int)
    for d in col_["fi_info"]:
        fi[(d["lib"], d["key"])] += d["frags"]
    rec_rows, present = [], defaultdict(dict)
    for lib in LIBS:
        for g5, g3 in CCLE8:
            k, kr = fkey(g5, g3), fkey(g3, g5)
            cs, csr = bykey[lib].get(k, {}), bykey[lib].get(kr, {})
            called = sorted(cs)
            present[lib][k] = len(called)
            ar = cs.get("arriba", {})
            dd = disc.get((lib, k), [])
            rec_rows.append(dict(lib=lib, marker=k, top4=int((g5, g3) in TOP4), n_callers=len(called),
                                 callers=",".join(called) or "-",
                                 arriba_dedup=ar.get("frags", ""), arriba_dups=ar.get("dups", ""),
                                 arriba_incl_dups=(ar.get("frags", 0) + ar.get("dups", 0)) if ar else "",
                                 arriba_conf=ar.get("conf", ""),
                                 starfusion_dedup=cs.get("starfusion", {}).get("frags", ""),
                                 starfusion_detail=cs.get("starfusion", {}).get("detail", ""),
                                 fusioncatcher_raw=cs.get("fusioncatcher", {}).get("frags", ""),
                                 reciprocal_callers=",".join(sorted(csr)) or "-",
                                 arriba_discarded="; ".join(f"{x['frags']} frags, {x['filters']}" for x in dd) or "-",
                                 fusioninspector_pass1_info=fi.get((lib, k), 0)))
    cols = list(rec_rows[0].keys())
    with open(f"{A}/recall_ccle.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rec_rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- identity verdict
    top4k = [fkey(*p) for p in TOP4]
    c2b_n = sum(1 for k in top4k if present["C2B"].get(k, 0) >= 1)
    lines = ["# Identity verdict (pass 1; each caller's own output; a marker is 'present' when at least one caller "
             "reports it in the CCLE orientation)", ""]
    lines.append(f"- C2B carries {c2b_n} of the top four CCLE U251MG markers: "
                 + ("**at least two, so the pipeline and the line are confirmed**" if c2b_n >= 2 else
                    "**fewer than two: the pipeline and the line are NOT confirmed**"))
    for lib in TUMOURS + ["N269B"] + CONTROLS:
        n = sum(1 for k in top4k if present[lib].get(k, 0) >= 1)
        verdict = ""
        if lib in TUMOURS:
            if n == 0 and c2b_n >= 2:
                verdict = " **MIX-UP until shown otherwise** (C2B has the markers, this tumour has none)"
            elif n < 4:
                verdict = " (a missing marker is not by itself a failure: U251N is not CCLE's U251MG, Torsvik 2014)"
        miss4 = [k for k in top4k if present[lib].get(k, 0) == 0]
        lines.append(f"- {lib}: {n} of 4" + (f"; missing {', '.join(miss4)}" if miss4 else "") + verdict)
    lines.append("")
    lines.append("Per marker, present in C2B and in all six tumours?")
    for k in [fkey(*p) for p in CCLE8]:
        where = [lib for lib in ["C2B"] + TUMOURS if present[lib].get(k, 0) >= 1]
        lines.append(f"- {k}: C2B {'yes' if 'C2B' in where else 'no'}, tumours {sum(1 for t in TUMOURS if t in where)}/6"
                     + (f"; scale if expression matched CCLE: C2B {SCALE[k][0]}, tumours {SCALE[k][1]}" if k in SCALE else ""))
    n2 = [r for r in rec_rows if r["lib"] == "N269B" and r["marker"] == "DHX38--AMPD1"][0]
    lines += ["", f"- N269B DHX38--AMPD1: callers {n2['callers']}; Arriba {n2['arriba_dedup'] or 0} deduplicated "
              f"({n2['arriba_incl_dups'] or 0} with duplicates), STAR-Fusion {n2['starfusion_dedup'] or 0}, FusionCatcher "
              f"{n2['fusioncatcher_raw'] or 0} (raw). Expected {N269B_SCALE}. Present reads show U251-lineage transcripts "
              "in the contralateral library; they do not separate resident cells from index hopping "
              "(M/reads/INDEX_HOPPING.txt): compare the rate in IL64B and N168B below and in the pass-2 evidence."]
    for lib in CONTROLS:
        n = sum(1 for k in top4k if present[lib].get(k, 0) >= 1)
        lines.append(f"- control {lib}: {n} of 4 top markers called"
                     + (" (U251 transcripts in a control: index hopping or residual tumour)" if n else ""))
    open(f"{A}/IDENTITY.md", "w").write("\n".join(lines) + "\n")

    # ---- whitelist for pass 2
    src = defaultdict(set)
    for g5, g3 in CCLE8:
        src[fkey(g5, g3)].add("ccle")
    sym = Sym(a.hgnc)
    for lib, lst in col_["fr_filtered"].items():
        for name in lst:
            parts = name.split("--")
            if len(parts) == 2:
                src[fkey(sym(parts[0]), sym(parts[1]))].add(f"fusionreport_filtered:{lib}")
            else:
                src[name].add(f"fusionreport_filtered:{lib}")
    for c in calls:
        if c["caller"] == "arriba" and c["conf"] == "high":
            src[c["key"]].add(f"arriba_high:{c['lib']}")
    keep, excl = [], []
    for k, s in sorted(src.items()):
        parts = k.split("--")
        why = ""
        if len(parts) != 2 or not parts[0] or not parts[1]:
            why = "not GENE--GENE"
        elif intergenic(parts[0]) or intergenic(parts[1]):
            why = "intergenic breakpoint"
        elif parts[0] == parts[1]:
            why = "same gene"
        elif re.search(r"\s", k):
            why = "whitespace"
        elif parts[0] not in ctat_genes or parts[1] not in ctat_genes:
            why = "gene not in the CTAT library (FusionInspector cannot build it)"
        (excl if why else keep).append((k, ",".join(sorted(s)), why))
    with open(f"{A}/whitelist_pass2.tsv", "w") as f:
        for k, _, _ in keep:
            f.write(k + "\n")
    with open(f"{A}/whitelist_pass2_sources.tsv", "w") as f:
        f.write("fusion\tsources\n")
        for k, s, _ in keep:
            f.write(f"{k}\t{s}\n")
    with open(f"{A}/whitelist_pass2_excluded.tsv", "w") as f:
        f.write("fusion\tsources\treason\n")
        for k, s, w in excl:
            f.write(f"{k}\t{s}\t{w}\n")
    has_ccle = all(fkey(*p) in {k for k, _, _ in keep} for p in CCLE8)
    row("whitelist_pass2", "non-empty and holds the eight CCLE markers", f"{len(keep)} fusions, {len(excl)} excluded,"
        f" CCLE markers all kept {has_ccle}", "PASS" if keep and has_ccle else "FAIL", "")

    # ---- write checks, summary, markers
    with open(f"{A}/checks.tsv", "w") as f:
        f.write("check\texpected\tobserved\tresult\toutside_reference\n")
        for r in ROWS:
            f.write("\t".join(r) + "\n")
    fails = [r for r in ROWS if r[3] == "FAIL"]
    summ = ["# rnafusion pass-1 audit", "", f"Written {time.strftime('%Y-%m-%d %H:%M')} on {os.uname().nodename}.", "",
            open(f"{A}/IDENTITY.md").read(),
            f"## New fusions", "",
            f"{len(rep)} sample-level calls meet every rule (2 of 3 callers; >= 3 deduplicated fragments at annotated "
            "splice sites; every Arriba/STAR-Fusion best breakpoint on annotated splice sites; not in IL64B/N168B; "
            "<= 50 % of supporting reads in 'both'; <= 50 % co-linear in rat; rat orthologs not within 200 kb on one "
            "strand). Table: reportable_new.tsv; every call with its reasons: scored_calls.tsv.", "",
            "Recurrence is not judged here: that needs pass 2's FusionInspector evidence in every sample and the "
            "10-split null (evidence/EVIDENCE.md).", "",
            "## Checks", "", f"{sum(1 for r in ROWS if r[3] == 'PASS')} PASS, "
            f"{sum(1 for r in ROWS if r[3] == 'WARN')} WARN, {len(fails)} FAIL (checks.tsv).", ""]
    summ += [f"- FAIL {r[0]}: expected {r[1]}, observed {r[2]}" for r in fails]
    summ += [f"- WARN {r[0]}: expected {r[1]}, observed {r[2]}" for r in ROWS if r[3] == "WARN"]
    open(f"{A}/AUDIT.md", "w").write("\n".join(summ) + "\n")
    if fails:
        open(f"{A}/AUDIT_FAILED", "w").write("\n".join("\t".join(r) for r in fails) + "\n")
        log(f"score: {len(fails)} hard check(s) failed; AUDIT_OK not written")
        return 1
    open(f"{A}/AUDIT_OK", "w").write(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n")
    log(f"score: AUDIT_OK; {len(keep)} fusions in whitelist_pass2.tsv; {len(rep)} reportable new calls")
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["collect", "extract", "rat", "score", "all"])
    ap.add_argument("--results", required=True)
    ap.add_argument("--audit", required=True)
    ap.add_argument("--reads", required=True, help="M/reads")
    ap.add_argument("--u", required=True, help="grid clone of the u251 repository")
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--gtf", required=True, help="GENCODE v46 GTF under M/refs/rnafusion")
    ap.add_argument("--ctat-spans", required=True)
    ap.add_argument("--rat-gtf", required=True)
    ap.add_argument("--rat-index", required=True)
    ap.add_argument("--orthologs", required=True)
    ap.add_argument("--star-img", default="")
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--xeng-override", default="", help="test fixtures only: JSON {lib: [graft, both]}")
    a = ap.parse_args()
    if a.xeng_override:
        XENG.update({k: tuple(v) for k, v in json.load(open(a.xeng_override)).items()})
    os.makedirs(a.audit, exist_ok=True)
    stages = ["collect", "extract", "rat", "score"] if a.stage == "all" else [a.stage]
    rc = 0
    for st in stages:
        done = f"{a.audit}/.done_{st}"
        if a.stage == "all" and os.path.isfile(done) and st != "score":
            log(f"{st}: done earlier, skipped")
            continue
        log(f"stage {st}")
        r = {"collect": collect, "extract": extract, "rat": rat, "score": score}[st](a)
        rc = r or 0
        open(done, "w").write(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n")
    return rc


if __name__ == "__main__":
    sys.exit(main())
