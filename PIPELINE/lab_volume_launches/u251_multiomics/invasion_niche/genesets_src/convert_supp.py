#!/usr/bin/env python3
"""Convert the fetched supplementary tables into plain gene-list TSVs, mechanically, for invasion_niche part C.

No gene is typed in: every symbol comes from a cell of the publisher's file. Cells that Excel turned into dates or
numbers are reported and dropped, never repaired by hand.

  python convert_supp.py <dir with the xlsx files>

Writes, in the same directory:
  neftel2019_tableS2_metamodules.tsv    set, rank, gene   (Neftel 2019 Cell, Table S2, meta-module gene lists)
  venkataramani2022_tableS2a_invasivity.tsv  gene, direction, r_S24, r_P3XX, r_T269  (Table S2a, invasivity score genes)
  SOURCES.tsv                            file, url, md5, sheet, rows written, cells dropped
"""
import datetime
import hashlib
import os
import sys

import openpyxl

D = sys.argv[1] if len(sys.argv) > 1 else "."
SRC = {
    "neftel2019_mmc2.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867419306877-mmc2.xlsx",
    "venkataramani2022_mmc2.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S0092867422008479-mmc2.xlsx",
    "darmanis2017_mmc6.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2211124717314626-mmc6.xlsx",
}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


src_rows = []

# ---- Neftel 2019 Table S2
f = os.path.join(D, "neftel2019_mmc2.xlsx")
wb = openpyxl.load_workbook(f, read_only=True, data_only=True)
ws = wb.worksheets[0]
rows = list(ws.iter_rows(values_only=True))
hdr_i = [i for i, r in enumerate(rows) if r and r[0] == "MES2"]
assert len(hdr_i) == 1, "Table S2 header row (MES2, MES1, AC, ...) not found exactly once"
hdr = [c.strip() if isinstance(c, str) else c for c in rows[hdr_i[0]]]
expect = ["MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "G1/S", "G2/M"]
assert hdr[:8] == expect, hdr
out, dropped = [], []
for j, name in enumerate(expect):
    k = 0
    for r in rows[hdr_i[0] + 1:]:
        v = r[j] if j < len(r) else None
        if v is None:
            continue
        if not isinstance(v, str):
            dropped.append("%s:%r" % (name, v))
            continue
        v = v.strip()
        if not v:
            continue
        k += 1
        out.append((name, k, v))
with open(os.path.join(D, "neftel2019_tableS2_metamodules.tsv"), "w", newline="\n") as fh:
    fh.write("set\trank\tgene\n")
    for t in out:
        fh.write("%s\t%d\t%s\n" % t)
src_rows.append(("neftel2019_mmc2.xlsx", SRC["neftel2019_mmc2.xlsx"], md5(f), ws.title, len(out), ";".join(dropped) or "none"))
print("Neftel Table S2:", {n: sum(1 for t in out if t[0] == n) for n in expect}, "dropped:", dropped or "none")

# ---- Venkataramani 2022 Table S2a
f = os.path.join(D, "venkataramani2022_mmc2.xlsx")
wb = openpyxl.load_workbook(f, read_only=True, data_only=True)
ws = [w for w in wb.worksheets if "2a" in w.title]
assert len(ws) == 1, [w.title for w in wb.worksheets]
ws = ws[0]
rows = list(ws.iter_rows(values_only=True))
hdr = list(rows[0])
assert hdr[0] == "Gene" and hdr[-1] == "Direction", hdr
ci = {k: hdr.index(k) for k in ("S24-correlation coefficient estimate", "P3XX-correlation coefficient estimate",
                                "T269-correlation coefficient estimate")}
out, dropped = [], []
for r in rows[1:]:
    g, d = r[0], r[-1]
    if g is None and d is None:
        continue
    if not isinstance(g, str) or d not in ("Correlated", "Anticorrelated"):
        dropped.append("%r/%r" % (g, d))
        continue
    out.append((g.strip(), d, r[ci["S24-correlation coefficient estimate"]], r[ci["P3XX-correlation coefficient estimate"]],
                r[ci["T269-correlation coefficient estimate"]]))
with open(os.path.join(D, "venkataramani2022_tableS2a_invasivity.tsv"), "w", newline="\n") as fh:
    fh.write("gene\tdirection\tr_S24\tr_P3XX\tr_T269\n")
    for t in out:
        fh.write("%s\t%s\t%s\t%s\t%s\n" % t)
src_rows.append(("venkataramani2022_mmc2.xlsx", SRC["venkataramani2022_mmc2.xlsx"], md5(f), ws.title, len(out), ";".join(dropped) or "none"))
print("Venkataramani S2a:", len(out), "genes;", sum(1 for t in out if t[1] == "Correlated"), "correlated,",
      sum(1 for t in out if t[1] == "Anticorrelated"), "anticorrelated; dropped:", dropped or "none")

# ---- Darmanis 2017 Table S6: record what it is (an embedded image), write nothing
f = os.path.join(D, "darmanis2017_mmc6.xlsx")
if os.path.isfile(f):
    wb = openpyxl.load_workbook(f)
    cells = sum(1 for w in wb.worksheets for r in w.iter_rows(values_only=True) for c in r if c is not None)
    src_rows.append(("darmanis2017_mmc6.xlsx", SRC["darmanis2017_mmc6.xlsx"], md5(f), wb.worksheets[0].title, 0,
                     "no text cells (%d); the table is an embedded image, not machine-readable" % cells))

with open(os.path.join(D, "SOURCES.tsv"), "w", newline="\n") as fh:
    fh.write("file\turl\tmd5\tsheet\trows_written\tcells_dropped\n")
    for t in src_rows:
        fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % t)
print("fetched on", datetime.date.today().isoformat())
