# IL68B held out: the full pipeline, re-run and compared

IL68B, a primary tumour, is the highest of the six on every translation-initiation and ribosome gene; without it the
gene-level fall of those sets is −0.10 instead of −0.40 (`SLIDES/07_contamination_magnitude.py`), although the GSEA
enrichment survives (`ANALYSIS/gsea_leave_one_out`). This folder re-runs every sample-dependent stage with IL68B removed
from the sample sheet and compares each against the published run.

| step | job script | what | notes |
|---|---|---|---|
| 1 | `run_de.sbatch` | nf-core/differentialabundance 1.5.0, `therapy_v3_params.yaml`, counts + gene lengths | the published command with `ANALYSIS/metadata_therapy.csv` (five tumours) and a new outdir |
| 2 | `run_figure.sbatch` | `create_publication_figure_600_dpi.R`: DSigDB drug GSEA, ChEMBL lookups (shared cache), STRING network, Supplementary_Data.xlsx | two arms: **control** (published six-tumour inputs, must reproduce the published drug profiles) and **holdout** |
| 3 | `run_drugs.sbatch` | the September chain from `subtypes/` (ChEMBL SMILES, ADMET-AI, BOILED-Egg + B3DB, final ranking; GSVA subtypes), then `compare.py` | unmodified scripts copied into each arm; environment built on CephFS |

`COMPARISON.md` and `comparison.json` hold the result: differential expression, GSEA leading sets, the 100 opposing
DSigDB compounds, the final ranking and ciclopirox, subtype scores and STRING hubs, each as published / control /
held out. Not rerun: the prior-art audit (manual; new top-20 entrants need it) and DepMap (cell-line property).

Submit from the repository root on the grid:

```bash
H=ANALYSIS/holdout_IL68B; mkdir -p $H/logs
A=$(sbatch --parsable -D $H $H/run_de.sbatch)
B=$(sbatch --parsable -D $H --dependency=afterany:$A $H/run_figure.sbatch)
sbatch -D $H --dependency=afterany:$B $H/run_drugs.sbatch
```

`subtypes_src/` (the `subtypes/` scripts, signature GMTs and the GEO count matrix) and `published/` (the September
`subtypes/` outputs and the manuscript drug-profile table) are copied in before submission.
