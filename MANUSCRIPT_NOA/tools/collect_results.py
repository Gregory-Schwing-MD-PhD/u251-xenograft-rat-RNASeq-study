#!/usr/bin/env python
"""Collect the result files the manuscript compiles from into MANUSCRIPT_NOA/results/, and index them by run.

Every file the paper quotes a number from must (a) be a byte copy of a grid output, (b) carry the sha256 of what it
was copied from, and (c) name the job that wrote it. This script makes results/ from the frozen copies that already
sit in SLIDES/figures_cns/v5_data, converts the TSV outputs to JSON so one code path reads everything, and writes
results/RUNS.json: the index from result file -> grid job -> what produced it.

    python tools/collect_results.py            # copy, convert, index
    python tools/collect_results.py --check    # verify results/ still matches its source copies; exit 1 if not
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
SRC = REPO / "SLIDES" / "figures_cns" / "v5_data"
OUT = ROOT / "results"

# file -> the grid job that wrote it, and one line on what it is. A file with no job here is refused: an unsourced
# number is exactly what this mechanism exists to prevent.
RUNS = {
    "composition_main.json": ("40481434 sci_rna_core", "Species-blind virtual bulk, exact mixing identity and Shapley decomposition; the primary endpoint."),
    "composition_S1.json": ("40481475 sci_s1_rnaseq", "Sensitivity S1: the same decomposition on graft-only reads (nf-core/rnaseq 3.22.2 re-quantification)."),
    "composition_S4.json": ("40481616 sci_s4", "Sensitivity S4: the same decomposition over HCOP one-to-one orthologs instead of Compara 110."),
    "c3_ha2_ha3.json": ("40481453/40481454 sci_sbde", "Species-blind differential expression over all ten balanced splits; H-A2 and H-A3."),
    "c4_categories_by_view.tsv": ("40481453 sci_sbde", "The three pre-stated GO categories across four views, Holm across 3 x 4."),
    "c5_l2_t1r.json": ("40481466 sci_rna_r", "GSVA mesenchymal scores with the dilution null, the ESTIMATE lever check, and the DiG robustness scorings."),
    "c6_summary.tsv": ("40481466 sci_rna_r", "The fraction-covariate model: collinearity, variance inflation, and how many tumour-cell changes it recovers."),
    "t1_e1.json": ("40481434 sci_rna_core", "The distal-invasion (DiG) programme, fraction-adjusted, over the twenty labellings; and the abandoned E1 gate."),
    "m4_host.json": ("40481582 sci_host", "Host deconvolution by three published methods with the pre-stated ordering, agreement and athymic controls."),
    "m4_lesion_accounting.tsv": ("40481582 sci_host", "Per-lesion accounting: tumour share plus host cell-type terms, summing to read composition."),
    "composition_samples_main.tsv": ("40481434 sci_rna_core", "Per-library tumour shares on every definition, and human chromosome Y per million human counts."),
    "xengsort_classes.tsv": ("xengsort 2.1.0 (May 2026 sort)", "Read pairs per library assigned graft, host, both, ambiguous or neither."),
    "M2_HM1_summary.tsv": ("40481488 sci_fdna", "Recurrent-minus-primary tumour share on each platform, with exact 90 % intervals and the labelling rank."),
    "fdna_resolved.tsv": ("40481488 sci_fdna", "Every tumour DNA fraction estimate per array, with which one is registered and which gates it passed."),
    "m3_rna_vs_dna.tsv": ("40481434 sci_rna_core", "RNA share against DNA fraction per lesion, and the cell-share range under two DNA indexes."),
    "L1_lump.tsv": ("40481447 sci_dna1", "The LUMP lever check per array: the test of whether the array betas read the human compartment only."),
    "l2_estimate_scores.tsv": ("40481466 sci_rna_r", "ESTIMATE stromal, immune and combined scores on the virtual bulk."),
    "D4_mgmt.json": ("40481447 sci_dna1", "MGMT-STP27 per array, with the flip distance and the observed recurrence shift."),
    "T3_fcpg.json": ("40481447 sci_dna1", "Fluctuating-CpG homogeneity with its two pre-stated gates."),
    "T4_S1.json": ("40481447 sci_dna1", "Convergence of derived arm-level copy-number change, and the smallest shared subclone the rule would detect."),
    "S1_smin.tsv": ("40481447 sci_dna1", "The detection bound: the smallest share of recurrent cells carrying a shared single-copy change that would be flagged."),
    "dilution_tests.json": ("40481489 sci_dilution", "Copy-number amplitude, MGMT promoter mean and fCpG index against the in-silico dilution nulls (UNVALIDATED)."),
    "HC2_dmp_counts_all_splits.tsv": ("40481491 sci_dna_integrate", "Probes at p < 0.001 for all ten balanced splits, per pipeline and design."),
    "HC2_benchlists_criterion.tsv": ("40481545 sci_benchlists", "The fraction-confound criterion under every published cross-species probe-exclusion list."),
    "validation.json": ("40481271 curve", "The human:mouse titration that the in-silico mixing had to reproduce, and its failure against the pre-set bar."),
    "estimates.tsv": ("40481271 curve", "Each array read on the titration curve, raw and mixture-corrected."),
    "criteria.json": ("40480791 spread", "The contralateral case: the pre-registered statistics against the depth-matched dilution null."),
    "c1.json": ("40480666 spread C1", "Arm-level copy number in the contralateral sample against its own core and the culture, with the noise bound."),
    "e3_N2_correlations.tsv": ("40480666 spread C1", "Bin-profile correlation of the contralateral array with the core and the culture."),
    "e2_chrY.tsv": ("40480665 spread E2", "Human chromosome Y counts per million human counts, per library."),
    "t2.json": ("40481604 sci_ai", "Allelic imbalance in species-exact human reads, its gates and its thinning null."),
    "power_closed_only.json": ("40481525 u251_pwr", "The pre-stated targets of the scaled experiment and the animals each needs, closed-form."),
    "c3_splits_summary.tsv": ("40481453 sci_sbde", "Differentially expressed genes per split, with each split's difference in tumour share."),
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tsv_to_json(src: Path, dst: Path) -> dict:
    """A TSV becomes {"_source":…, "_sha256":…, "rows": {first column: {column: value}}, "order": [...]}.

    Values stay strings: the manifest's fmt decides how a number is printed, and a string cannot silently change
    precision on the way in.
    """
    with open(src, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit(f"{src} is empty")
    key = list(rows[0].keys())[0]
    names = [(r[key] or f"row{i}") for i, r in enumerate(rows)]
    if len(set(names)) != len(names):
        # a first column that repeats (S1_smin's dc is -1, -1, 1, 1) would collapse rows on top of each other;
        # key by position instead and say so, so a pointer is never silently pointing at the wrong row
        names = [str(i) for i in range(len(rows))]
        key = f"{key} (repeats; rows keyed by position)"
    doc = {"_source": str(src), "_sha256": sha256(src), "_key_column": key or "(unnamed)",
           "rows": {n: r for n, r in zip(names, rows)},
           "order": names,
           "n_rows": len(rows)}
    dst.write_text(json.dumps(doc, indent=1), encoding="utf-8", newline="\n")
    return doc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    index, missing, changed = {}, [], []
    for name, (job, what) in sorted(RUNS.items()):
        src = SRC / name
        if not src.exists():
            missing.append(name)
            continue
        if name.endswith(".tsv"):
            dst = OUT / (name[:-4] + ".json")
            if a.check and dst.exists():
                old = json.load(open(dst, encoding="utf-8")).get("_sha256")
                if old != sha256(src):
                    changed.append(name)
            else:
                tsv_to_json(src, dst)
        else:
            dst = OUT / name
            if a.check and dst.exists():
                if sha256(dst) != sha256(src):
                    changed.append(name)
            else:
                shutil.copyfile(src, dst)
        index[dst.name] = {"job": job, "what": what, "copied_from": str(src.relative_to(REPO)),
                           "sha256": sha256(dst) if dst.exists() else None}
    (OUT / "RUNS.json").write_text(json.dumps(
        {"about": "Every result file the manuscript compiles from, the grid job that wrote it, and what it is. "
                  "A number in the paper resolves to file#pointer through numbers_manifest.json, and the file "
                  "resolves to a job here.",
         "files": index}, indent=1), encoding="utf-8", newline="\n")
    print(f"results/: {len(index)} files indexed")
    for n in missing:
        print(f"   MISSING source: {n}")
    for n in changed:
        print(f"   CHANGED since it was copied: {n}")
    return 1 if (missing or changed) else 0


if __name__ == "__main__":
    sys.exit(main())
