#!/usr/bin/env python
"""Build the one table that says which RNA library and which methylation array are the same piece of tissue.

The study deposits two assay types under two different sets of sample names, because two facilities named the same
animals differently: the laboratory's own key and the array core's sample sheet call rats 70 and 71 IL-70/IL-71,
while the sequencing core's FASTQ file names -- and therefore GEO, SRA and every RNA accession -- call them NL70B
and NL71B. The contralateral hemisphere of rat 69 is N269B on the RNA side and "N2" on the array sheet. Neither name
is going to be changed: the RNA names are public and minted in SRA, and the array names are what the chip, the core's
sheet and every analysis in this repository use.

So the correspondence is carried explicitly, here and in the GEO metadata, and never inferred from the names.

    python GEO/make_sample_key.py

Writes GEO/multiomics_sample_key.tsv. It is generated, not maintained: the RNA columns come from
PIPELINE/00_fetch_reads/library_map.csv (checked against SRA and ENA) and the array columns from
GEO/methylation/array_samplesheet.csv (the array core's sheet for plate 1724, chip 205648300021). The two files are
cross-checked against each other and the script refuses to write a key if they disagree.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MAP = REPO / "PIPELINE" / "00_fetch_reads" / "library_map.csv"
SHEET = REPO / "GEO" / "methylation" / "array_samplesheet.csv"
OUT = REPO / "GEO" / "multiomics_sample_key.tsv"

RNA_SERIES = "GSE338105"
ARRAY_SERIES = "PENDING"        # filled in once the methylation Series is issued

COLUMNS = [
    "subject", "arm", "tissue",
    "rna_library", "rna_geo_sample", "rna_sra_experiment", "rna_sra_run", "rna_biosample", "rna_series",
    "array_label", "array_geo_sample", "array_sentrix_id", "array_sentrix_position", "array_beta_column",
    "array_series", "paired", "note",
]


def load() -> tuple[list[dict], dict[str, dict]]:
    rna = list(csv.DictReader(open(MAP, encoding="utf-8-sig", newline="")))
    arrays = {r["sample"]: r for r in csv.DictReader(open(SHEET, encoding="utf-8-sig", newline=""))}
    return rna, arrays


def crosscheck(rna: list[dict], arrays: dict[str, dict]) -> list[str]:
    """Every claim of correspondence must be made by both files, or the key is not written."""
    problems = []
    by_lib = {r["library"]: r for r in rna}

    # 1. every array row names an RNA library that exists
    for s, a in arrays.items():
        if a["rna_library"] not in by_lib:
            problems.append(f"array {s} names RNA library {a['rna_library']}, which is not in {MAP.name}")

    # 2. the two files agree, in both directions, on which array goes with which library
    for s, a in arrays.items():
        lib = a["rna_library"]
        if lib in by_lib and by_lib[lib]["array_label"] != s:
            problems.append(f"array sheet pairs {s} with {lib}, but {MAP.name} pairs {lib} with "
                            f"{by_lib[lib]['array_label'] or '(no array)'}")
    for r in rna:
        lab = r["array_label"]
        if lab and lab not in arrays:
            problems.append(f"{MAP.name} pairs {r['library']} with array {lab}, which is not on the chip sheet")

    # 3. the rat number must agree wherever both files state one
    for s, a in arrays.items():
        lib = a["rna_library"]
        if lib in by_lib and a["rat"] and by_lib[lib]["rat"] and a["rat"] != by_lib[lib]["rat"]:
            problems.append(f"{s}/{lib}: array sheet says rat {a['rat']}, {MAP.name} says rat {by_lib[lib]['rat']}")

    # 4. one chip position per array, and no position used twice
    seen: dict[str, str] = {}
    for s, a in arrays.items():
        pos = f"{a['sentrix_id']}_{a['sentrix_position']}"
        if pos in seen:
            problems.append(f"chip position {pos} is claimed by both {seen[pos]} and {s}")
        seen[pos] = s
    return problems


def main() -> int:
    rna, arrays = load()
    problems = crosscheck(rna, arrays)
    if problems:
        print("REFUSING to write the sample key: the two source files disagree.", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1

    rows = []
    for r in rna:
        lab = r["array_label"]
        a = arrays.get(lab, {})
        paired = "RNA and array from different aliquots of the same piece" if a else "RNA only; no array was run"
        note = r["note"]
        if lab and lab != r["library"]:
            note = f"the same piece is called {r['library']} on the RNA side and {lab} on the array side; {note}"
        rows.append({
            "subject": f"rat {r['rat']}" if r["rat"] else "U251N culture",
            "arm": r["arm"],
            "tissue": r["tissue"],
            "rna_library": r["library"],
            "rna_geo_sample": r["geo_sample"],
            "rna_sra_experiment": r["experiment_accession"],
            "rna_sra_run": r["run_accession"],
            "rna_biosample": r["sample_accession"],
            "rna_series": RNA_SERIES,
            "array_label": lab,
            "array_geo_sample": ARRAY_SERIES if lab else "",
            "array_sentrix_id": a.get("sentrix_id", ""),
            "array_sentrix_position": a.get("sentrix_position", ""),
            "array_beta_column": f"{lab}_EPIC" if lab else "",
            "array_series": ARRAY_SERIES if lab else "",
            "paired": paired,
            "note": note,
        })

    with open(OUT, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    n_paired = sum(1 for r in rows if r["array_label"])
    print(f"wrote {OUT.relative_to(REPO)}")
    print(f"  {len(rows)} RNA libraries, {n_paired} of them with a matched methylation array, "
          f"{len(rows) - n_paired} without")
    print(f"  the three pieces whose two names differ: "
          + "; ".join(f"{r['rna_library']} = {r['array_label']}"
                      for r in rows if r["array_label"] and r["array_label"] != r["rna_library"]))
    print("  cross-check passed: both source files agree on every pairing, rat number and chip position")
    return 0


if __name__ == "__main__":
    sys.exit(main())
