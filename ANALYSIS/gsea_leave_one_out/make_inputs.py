#!/usr/bin/env python
"""Leave-one-tumour-out inputs for the U251 GSEA (nf-core/differentialabundance GSEA_GSEA step, GSEA 4.3.2).

Why: the translation-initiation, elongation and ribosome sets fall about -0.40 log2 in recurrence, but one primary,
IL68B, is the highest of the six tumours on every translation-initiation and ribosome gene; without it the sets fall
about -0.10 (SLIDES/07_contamination_magnitude.py). The question this answers: does the GSEA result itself (translation
initiation NES -1.99, FDR q = 0.022, the only set with q < 0.05; six sets at q < 0.25, all down) survive when any one
tumour is left out?

Inputs are the pipeline's own GSEA inputs, copied from the nf-core work directory by prep.sbatch: the expression
matrix (therapy_impact.gct = the pipeline's normalised counts of the six tumours, 19,351 genes), the class file, the
gene-set collection (combined_human.gmt) and the chip file. For each left-out tumour the column is dropped and the
pipeline's gene filter is re-applied to the remaining five (>= 10 raw reads in >= 1 sample; with all six this
reproduces the 19,351 genes exactly, asserted below). The five remaining columns keep the six-tumour normalisation:
re-deriving size factors on five moves their relative values by at most 0.3 % (checked on the counts).

Conditions (tasks.tsv, one array task each, 41 in all):
  full                    all six tumours; seed 1234 must reproduce the original report exactly; seeds 1-4 measure the
                          Monte-Carlo spread of the FDR q values (gene-set permutation, 1,000 permutations)
  drop_<tumour>           one tumour left out, pipeline filter re-applied; seeds 1234 and 1-4
  drop_<tumour>_nofilter  one tumour left out, all 19,351 genes kept; seed 1234 (sensitivity to the filter)

Paths can be overridden with U251_DIR, GSEA_WORKDIR and U251_COUNTS (used for a local dry run).

    python ANALYSIS/gsea_leave_one_out/make_inputs.py
"""
import os
from pathlib import Path

import pandas as pd

U = Path(os.environ.get("U251_DIR", "/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study"))
D = U / "ANALYSIS" / "gsea_leave_one_out"
IN = D / "inputs"
COUNTS = Path(os.environ.get("U251_COUNTS", U / "ANALYSIS" / "results_human_final" / "star_salmon" / "salmon.merged.gene_counts.tsv"))
PRI = ["IL67B", "IL68B", "IL69B"]
REC = ["IL66B", "NL70B", "NL71B"]
SEEDS = [1234, 1, 2, 3, 4]

lines = (IN / "therapy_impact.gct").read_text().splitlines()
assert lines[0].startswith("#1.2"), lines[0]
nrow, ncol = (int(x) for x in lines[1].split("\t")[:2])
header = lines[2].split("\t")
samples = header[2:]
assert samples == PRI + REC, samples
rows = [ln.split("\t") for ln in lines[3:] if ln.strip()]
assert len(rows) == nrow == 19351 and ncol == 6, (len(rows), nrow, ncol)
ids = [r[0] for r in rows]
cls = (IN / "therapy_impact.cls").read_text()
assert cls.splitlines()[:3] == ["6 2 1", "#Primary_U2 Recurrent_U2", "Primary_U2 Primary_U2 Primary_U2 Recurrent_U2 Recurrent_U2 Recurrent_U2"], cls

cnt = pd.read_csv(COUNTS, sep="\t").set_index("gene_id")
keep6 = set(cnt.index[(cnt[samples] >= 10).sum(1) >= 1])
assert keep6 == set(ids), "the pipeline's gene filter (>= 10 reads in >= 1 of the six) is not reproduced from these counts"


def write(cond, cols, refilter):
    """A GCT and CLS for the columns given, in the original row order and number formatting."""
    idx = [header.index(c) for c in cols]
    if refilter:
        keep = set(cnt.index[(cnt[cols] >= 10).sum(1) >= 1])
        sel = [r for r in rows if r[0] in keep]
    else:
        sel = rows
    out = IN / cond
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "therapy_impact.gct", "w", newline="\n") as f:
        f.write("#1.2\t\n")
        f.write(f"{len(sel)}\t{len(cols)}\t\n")
        f.write("\t".join(["NAME", "DESCRIPTION"] + cols) + "\n")
        for r in sel:
            f.write("\t".join([r[0], r[1]] + [r[i] for i in idx]) + "\n")
    labels = ["Primary_U2" if c in PRI else "Recurrent_U2" for c in cols]
    assert labels.count("Primary_U2") >= 2 and labels.count("Recurrent_U2") >= 2
    (out / "therapy_impact.cls").write_text(f"{len(cols)} 2 1\n#Primary_U2 Recurrent_U2\n{' '.join(labels)}\n\n", newline="\n")
    return len(sel)


tasks = []
n = write("full", samples, refilter=False)
# the rebuilt full matrix must be the pipeline's matrix, value for value
assert (IN / "full" / "therapy_impact.gct").read_text().splitlines()[2:] == lines[2:], "rebuilt full GCT differs from the pipeline's"
for s in SEEDS:
    tasks.append(("full", s, "full" if s == 1234 else "light", n, ",".join(samples)))
for drop in samples:
    cols = [c for c in samples if c != drop]
    n = write(f"drop_{drop}", cols, refilter=True)
    for s in SEEDS:
        tasks.append((f"drop_{drop}", s, "full" if s == 1234 else "light", n, ",".join(cols)))
for drop in samples:
    cols = [c for c in samples if c != drop]
    n = write(f"drop_{drop}_nofilter", cols, refilter=False)
    tasks.append((f"drop_{drop}_nofilter", 1234, "light", n, ",".join(cols)))
assert len(tasks) == 41, len(tasks)

with open(D / "tasks.tsv", "w", newline="\n") as f:
    f.write("task\tcondition\tseed\toutput\tn_genes\tsamples\n")
    for i, t in enumerate(tasks):
        f.write("\t".join(str(x) for x in (i,) + t) + "\n")
print(f"wrote {len(tasks)} tasks to {D / 'tasks.tsv'}")
for t in tasks:
    if t[1] == 1234:
        print(f"  {t[0]:24s} genes {t[3]:6d}  samples {t[4]}")
