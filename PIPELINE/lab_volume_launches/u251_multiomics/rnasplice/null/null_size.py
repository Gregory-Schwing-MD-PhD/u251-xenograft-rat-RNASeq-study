#!/usr/bin/env python3
"""Size the null arrays from the pipeline's own trace. Prints: DEU_TIME DEU_MEM FAST_TIME FAST_MEM (slurm syntax).

DEU task = run_dexseq_exon.R once (the pipeline's DEXSEQ_EXON ran the same script on the same counts).
FAST task = run_dexseq_dtu.R + run_stager.R + two suppa diffSplice + rmats post (4 threads, the pipeline used 12),
with room for the --task both fallback (prep at 4 threads). Defaults apply when a trace row is missing."""
import csv
import glob
import os
import re
import sys


def secs(s):
    if not s or s == "-":
        return None
    t = 0.0
    for v, u in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", s):
        t += float(v) * {"ms": 0.001, "d": 86400, "h": 3600, "m": 60, "s": 1}[u]
    return t


def gb(s):
    m = re.match(r"([\d.]+)\s*(B|KB|MB|GB|TB)", s or "")
    if not m:
        return None
    return float(m.group(1)) * {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1.0, "TB": 1e3}[m.group(2)]


def slurm_time(sec):
    sec = int(sec)
    d, r = divmod(sec, 86400)
    h, r = divmod(r, 3600)
    return "%d-%02d:%02d:00" % (d, h, r // 60)


def main(info):
    rt, rss = {}, {}
    for tr in sorted(glob.glob(os.path.join(info, "execution_trace_*.txt"))):
        with open(tr) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                if row.get("status") not in ("COMPLETED", "CACHED"):
                    continue
                name = row.get("name", "")
                if re.search(r"(DRIMSEQ_DEXSEQ_DTU_SALMON|SUPPA_SALMON):", name):
                    continue
                p = name.split(" (")[0].split(":")[-1]
                r, m = secs(row.get("realtime")), gb(row.get("peak_rss"))
                if r is not None:
                    rt[p] = max(rt.get(p, 0), r)
                if m is not None:
                    rss[p] = max(rss.get(p, 0), m)
    h = 3600.0
    deu_t = 2 * rt["DEXSEQ_EXON"] + h if "DEXSEQ_EXON" in rt else 36 * h
    deu_m = 1.5 * rss["DEXSEQ_EXON"] + 4 if "DEXSEQ_EXON" in rss else 48
    parts = [rt.get("DEXSEQ_DTU", 3 * h), rt.get("STAGER", 0.2 * h), rt.get("DIFFSPLICE_IOE", h), rt.get("DIFFSPLICE_IOI", h),
             3 * rt.get("RMATS_POST", 2 * h), 3 * rt.get("RMATS_PREP", 3 * h)]
    fast_t = 1.5 * sum(parts) + h
    fast_m = 1.5 * max([rss.get(k, 0) for k in ("DEXSEQ_DTU", "STAGER", "DIFFSPLICE_IOE", "DIFFSPLICE_IOI", "RMATS_POST", "RMATS_PREP")] + [0]) + 4
    deu_t = min(max(deu_t, 4 * h), 96 * h)
    fast_t = min(max(fast_t, 3 * h), 72 * h)
    deu_m = int(min(max(deu_m, 16), 200))
    fast_m = int(min(max(fast_m, 16), 200))
    print(slurm_time(deu_t), "%dG" % deu_m, slurm_time(fast_t), "%dG" % fast_m)


if __name__ == "__main__":
    main(sys.argv[1])
