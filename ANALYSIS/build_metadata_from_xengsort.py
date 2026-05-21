#!/usr/bin/env python3
"""Append xengsort contamination fractions to the study metadata.

Intermediate step between xengsort (host/graft read sorting) and the RUVSeq
contamination-adjustment analysis. Reads the per-sample xengsort classify logs
and appends graft_pct / host_pct / both_pct columns to the design metadata --
so those fractions are derived from the pipeline output rather than hand-typed.

xengsort 2.1.0 writes one log per sample to <xengsort-dir>/<sample>.txt. Each
log contains a machine-readable counts block under "## Classification
Statistics":

    ## Classification Statistics

    ```
    prefix  host    graft   ambiguous       both    neither
    IL64B   90003153        314037  175597  4546256 1124564
    ```

This script parses that header + counts row (whitespace/tab separated, column
order taken from the header so it is robust to reordering) and computes each
class percentage as count / total * 100, where total is the sum of all five
classes. That reproduces the percentages xengsort prints in its pipe-table
block (e.g. IL64B graft = 314037 / 96163607 = 0.33%).

Usage:
    python3 build_metadata_from_xengsort.py \
        --base ANALYSIS/metadata_base.csv \
        --xengsort-dir ANALYSIS/xengsort_out \
        --out ANALYSIS/metadata_full.csv
"""

import argparse
import csv
import os
import sys

# Classes xengsort 2.1.0 reports, in the order they normally appear. The parser
# does not rely on this order (it reads the header), but it defines what we
# expect to find and which fractions land in the metadata.
XENGSORT_CLASSES = ("host", "graft", "ambiguous", "both", "neither")
# Columns appended to the metadata, mapped to their xengsort class.
APPENDED = (("graft_pct", "graft"), ("host_pct", "host"), ("both_pct", "both"))

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_xengsort_log(path):
    """Return {class_name: int_count} from a single xengsort classify log.

    Locates the counts block: a header line whose first token is 'prefix' and
    which lists the class names, followed by a data line of integer counts.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()

    for i, line in enumerate(lines):
        tokens = line.split()
        # Header row of the TSV counts block, e.g.
        #   prefix  host  graft  ambiguous  both  neither
        if tokens and tokens[0] == "prefix" and "host" in tokens and "graft" in tokens:
            class_cols = tokens[1:]
            # The next non-empty line is the per-sample counts row.
            for data_line in lines[i + 1:]:
                data_tokens = data_line.split()
                if not data_tokens:
                    continue
                values = data_tokens[1:]  # drop the sample-name prefix
                if len(values) != len(class_cols):
                    raise ValueError(
                        f"{path}: counts row has {len(values)} values but header "
                        f"has {len(class_cols)} classes ({class_cols})")
                try:
                    counts = {cls: int(val) for cls, val in zip(class_cols, values)}
                except ValueError as exc:
                    raise ValueError(f"{path}: non-integer count in row "
                                     f"'{data_line.strip()}'") from exc
                return data_tokens[0], counts
            raise ValueError(f"{path}: found counts header but no data row after it")

    raise ValueError(f"{path}: no '## Classification Statistics' counts block found")


def percentages(counts):
    """Convert a class->count dict into class->percent (of all reads)."""
    total = sum(counts.values())
    if total == 0:
        raise ValueError("total read count is zero")
    return {cls: 100.0 * n / total for cls, n in counts.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.path.join(SCRIPT_DIR, "metadata_base.csv"),
                    help="Design metadata CSV with columns sample,Classification")
    ap.add_argument("--xengsort-dir", default=os.path.join(SCRIPT_DIR, "xengsort_out"),
                    help="Directory of per-sample xengsort logs (<sample>.txt)")
    ap.add_argument("--out", default=os.path.join(SCRIPT_DIR, "metadata_full.csv"),
                    help="Output metadata CSV (base columns + appended fractions)")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Skip (with a warning) samples whose log is missing "
                         "instead of failing")
    args = ap.parse_args()

    if not os.path.isfile(args.base):
        sys.exit(f"ERROR: base metadata not found: {args.base}")
    if not os.path.isdir(args.xengsort_dir):
        sys.exit(f"ERROR: xengsort directory not found: {args.xengsort_dir}")

    with open(args.base, newline="") as fh:
        reader = csv.DictReader(fh)
        base_fields = reader.fieldnames or []
        if "sample" not in base_fields or "Classification" not in base_fields:
            sys.exit("ERROR: base metadata must have 'sample' and "
                     "'Classification' columns")
        rows = list(reader)

    out_fields = list(base_fields) + [c for c, _ in APPENDED
                                      if c not in base_fields]

    print(f"Base metadata : {args.base} ({len(rows)} samples)")
    print(f"xengsort dir  : {args.xengsort_dir}")
    print("-" * 78)
    header = f"{'sample':<8} {'class':<10} " + "  ".join(
        f"{cls:>9}" for cls in XENGSORT_CLASSES) + "   graft%   host%   both%"
    print(header)

    out_rows = []
    failures = []
    for row in rows:
        sample = row["sample"]
        log_path = os.path.join(args.xengsort_dir, f"{sample}.txt")
        if not os.path.isfile(log_path):
            msg = f"{sample}: log not found ({log_path})"
            if args.allow_missing:
                print(f"WARNING: {msg} -- skipped")
                failures.append(sample)
                continue
            sys.exit(f"ERROR: {msg}")

        try:
            prefix, counts = parse_xengsort_log(log_path)
            missing = [c for c in XENGSORT_CLASSES if c not in counts]
            if missing:
                raise ValueError(f"missing classes {missing}")
            pct = percentages(counts)
        except ValueError as exc:
            if args.allow_missing:
                print(f"WARNING: {sample}: parse failed -- {exc} -- skipped")
                failures.append(sample)
                continue
            sys.exit(f"ERROR: {sample}: parse failed -- {exc}")

        if prefix != sample:
            print(f"WARNING: {sample}: log prefix is '{prefix}' (using filename "
                  f"sample name)")

        new_row = dict(row)
        for col, cls in APPENDED:
            new_row[col] = f"{pct[cls]:.2f}"
        out_rows.append(new_row)

        print(f"{sample:<8} {row['Classification']:<10} " +
              "  ".join(f"{counts[c]:>9d}" for c in XENGSORT_CLASSES) +
              f"   {pct['graft']:6.2f}  {pct['host']:6.2f}  {pct['both']:6.2f}")

    if not out_rows:
        sys.exit("ERROR: no samples were successfully processed")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    print("-" * 78)
    print(f"Wrote {args.out} ({len(out_rows)} samples, "
          f"columns: {', '.join(out_fields)})")
    if failures:
        print(f"Samples skipped: {', '.join(failures)}")


if __name__ == "__main__":
    main()
