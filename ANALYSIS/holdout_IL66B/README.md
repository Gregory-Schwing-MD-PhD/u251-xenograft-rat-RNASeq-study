# IL66B held out: the full pipeline, re-run and compared

IL66B (a recurrent tumour) is the sample whose removal separates primary from recurrent best on the human reads
(`ANALYSIS/holdout_separation/loo_separation_human.tsv`: 349 genes at padj < 0.05 against 102 with all six;
silhouette 0.43 against 0.22). That choice is post hoc: the table was read before this rerun, so this is a description
of what the data look like without IL66B, not a test. Same stages as `holdout_IL68B` (DE, R drug screen and network,
the September Python chain, GSVA, compare.py); the control arm is the one `holdout_IL68B` produced (identical inputs
and software), copied rather than re-run.
