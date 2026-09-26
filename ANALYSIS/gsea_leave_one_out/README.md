# Leave-one-tumour-out GSEA

Does the recurrence GSEA result survive when any one of the six tumours is left out? The translation-initiation,
elongation and ribosome sets fall about −0.40 log2 in recurrence, but one primary, IL68B, is the highest of the six
tumours on every translation-initiation and ribosome gene; without it the gene-level fall is about −0.10
(`SLIDES/07_contamination_magnitude.py`). The published GSEA has translation initiation as the only set at FDR
q < 0.05 (NES −1.99, q = 0.022) and six sets at q < 0.25, all down.

The runs repeat the nf-core/differentialabundance `GSEA_GSEA` task of the therapy_v3 run exactly (GSEA 4.3.2, the same
Galaxy-depot image, the same command and the pipeline's own matrix, class file, gene-set collection and chip file,
copied from its work directory) with one tumour removed and the pipeline's gene filter re-applied. The full six-tumour
run with seed 1234 must reproduce the published report set for set (`reproduction.json`).

| file | what |
|---|---|
| `make_inputs.py` | per-condition GCT and CLS files and `tasks.tsv` (41 runs: full and six leave-one-out conditions at seeds 1234 and 1–4; six unfiltered leave-one-out runs at seed 1234) |
| `prep.sbatch` | pulls the image, copies the pipeline inputs, runs a 10-permutation smoke test, then `make_inputs.py` |
| `gsea_array.sbatch` | one GSEA run per task |
| `collect.py`, `collect.sbatch` | `reproduction.json`, `loo_sets.tsv`, `loo_screen.tsv`, `loo_top_down.tsv`, `SUMMARY.md` |

Run on the grid from the repository root:

```bash
D=ANALYSIS/gsea_leave_one_out; mkdir -p $D/logs
P=$(sbatch --parsable -D $D $D/prep.sbatch)
A=$(sbatch --parsable -D $D --dependency=afterany:$P $D/gsea_array.sbatch)
sbatch -D $D --dependency=afterany:$A $D/collect.sbatch
```
