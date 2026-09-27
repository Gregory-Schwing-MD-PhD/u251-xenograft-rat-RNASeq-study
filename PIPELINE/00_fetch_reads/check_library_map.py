#!/usr/bin/env python
"""Check that what the archive returned is what this project's analyses expect, by run accession and by name.

The analyses refer to libraries by the laboratory's names (IL67B, NL70B, N269B, C2B ...). The archive returns run
accessions. library_map.csv is the single place those two vocabularies meet, so a rename or a mis-ordered download
cannot quietly shift a sample from one arm to the other - which is the failure this whole chain exists to prevent.

    python check_library_map.py --samplesheet <fetchngs samplesheet.csv> --map library_map.csv
    python check_library_map.py --map library_map.csv          # check the map's own consistency only

Exit 1 on any disagreement: a run in the samplesheet that the map does not name, a library in the map that the
archive did not return, or a duplicate on either side.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ARMS = {"primary", "recurrent", "contralateral", "no_tumour", "culture"}


def read_map(p: Path):
    rows = list(csv.DictReader(open(p, encoding="utf-8-sig", newline="")))
    need = {"library", "arm", "run_accession", "sample_accession", "geo_sample"}
    missing = need - set(rows[0] if rows else {})
    if missing:
        raise SystemExit(f"{p}: missing columns {sorted(missing)}")
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet")
    ap.add_argument("--map", required=True)
    a = ap.parse_args()
    rows = read_map(Path(a.map))
    rc = 0

    pend = [r for r in rows if not r["run_accession"].strip() or r["run_accession"].strip().upper().startswith("PENDING")]
    if pend:
        print(f"MAP INCOMPLETE: {len(pend)} librar(ies) have no run accession yet: "
              + ", ".join(r["library"] for r in pend))
        rc = 1
    for field in ("library", "run_accession"):
        seen = [r[field].strip() for r in rows if r[field].strip()]
        dupes = {x for x in seen if seen.count(x) > 1}
        if dupes:
            print(f"MAP ERROR: duplicate {field}: {sorted(dupes)}")
            rc = 1
    bad_arm = {r["arm"] for r in rows} - ARMS
    if bad_arm:
        print(f"MAP ERROR: unknown arm(s) {sorted(bad_arm)}; expected {sorted(ARMS)}")
        rc = 1

    if a.samplesheet:
        ss = list(csv.DictReader(open(a.samplesheet, encoding="utf-8-sig", newline="")))
        col = next((c for c in ("run_accession", "sample", "fastq_1") if ss and c in ss[0]), None)
        if col is None:
            print(f"SAMPLESHEET: no run/sample column found in {a.samplesheet}")
            return 1
        got = {r[col].split("_")[0].strip() for r in ss}
        want = {r["run_accession"].strip() for r in rows if r["run_accession"].strip()}
        extra, absent = got - want, want - got
        if extra:
            print(f"SAMPLESHEET has runs the map does not name: {sorted(extra)}")
            rc = 1
        if absent:
            print(f"MAP names runs the archive did not return: {sorted(absent)}")
            rc = 1
        if not rc:
            print(f"OK: {len(want)} runs, each named by the map, no extras")
    elif not rc:
        print(f"OK: map is internally consistent ({len(rows)} libraries)")
    return rc


if __name__ == "__main__":
    sys.exit(main())
