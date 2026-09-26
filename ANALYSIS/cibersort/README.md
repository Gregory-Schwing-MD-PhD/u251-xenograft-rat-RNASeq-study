# Glioblastoma cell-state deconvolution (CIBERSORT, Neftel 2019 reference)

Deconvolves the human-read bulk profiles (six tumours, the in vitro culture C2B, three control brains) into the four
Neftel 2019 malignant states, MES-like, AC-like, OPC-like and NPC-like.

## Reference

Neftel et al. 2019, *Cell* 178:835, Smart-seq2, Broad Single Cell Portal SCP393, already on the grid at
`ANALYSIS/refs/neftel_2019/` (see its README). Adult IDH-wt malignant cells: 4,916 from 20 tumours.

- **State of a cell.** The metadata gives six meta-module scores, not a state label (`CellAssignment` is
  Malignant / Macrophage / Oligodendrocyte / T-cell only, so the May script that filtered it for "MES-like" etc.
  kept no cells). A cell's state is the largest of MES = max(MESlike1, MESlike2), AC, OPC, NPC = max(NPClike1,
  NPClike2). **Confident** cells: that score is positive and beats the runner-up by at least 0.5. This cut is ours,
  not Neftel's hybrid rule; a variant using every cell is built beside it.
- **Scale.** The matrix is log2(TPM/10 + 1); linear TPM = 10 · (2^x − 1). Checked: every cell then sums to 1,000,000.
- **Signature matrices** (classic CIBERSORT format, GeneSymbol × type, linear TPM), built by CIBERSORT's published
  procedure (Newman 2015): Welch t-test on log expression against all other reference cells, q < 0.3, higher in
  the type, ranked by fold change, top G per type for G = 50–200, keeping the G with the lowest condition number.
  Only genes also present in the bulk matrix are candidates.
  - `signature_neftel4_confident.txt`: the four states from confident cells (main)
  - `signature_neftel4_allcells.txt`: every cell by its top state (sensitivity)
  - `signature_neftel4_confident_plus_nonmalignant.txt`: plus macrophage, oligodendrocyte and T cell. The human
    reads come from the tumour cells only (the host is rat), so these should read near zero; they are a check.
- `scref_neftel4_cibersortx.txt.gz`: a CIBERSORTx single-cell reference (header = state of each cell; 300
  confident cells per state, sampled with a fixed seed), for CIBERSORTx's own signature building.
- `mixture_u251_tpm.txt`: the salmon gene TPM of all ten samples, summed per gene symbol.

## Checks

- **Leave-one-tumour-out** (`loto_validation.tsv`): for each reference tumour, the signature is rebuilt without it,
  its malignant cells are averaged into a pseudo-bulk and the estimated state fractions are compared with the true
  shares. Same platform on both sides: this tests the reference, not the Smart-seq2-to-bulk gap.
- **CIBERSORT-core re-implementation** (`fractions_nusvr.tsv`): nu-SVR, linear kernel, nu 0.25 / 0.5 / 0.75,
  the lowest-RMSE fit, negative weights to zero, signature standardised as one block and each mixture per sample,
  no quantile normalisation (RNA-seq), 1000-permutation p. It is a stand-in until the official CIBERSORT is
  reinstalled, and then an independent check on it: the two should agree closely.

## Running

    sbatch -D ANALYSIS/cibersort ANALYSIS/cibersort/run.sbatch      # ~5 min on reqp

Official CIBERSORT (once reinstalled): signature = `results/signature_neftel4_confident.txt`, mixture =
`results/mixture_u251_tpm.txt`, QN = FALSE, perm = 1000. CIBERSORTx: either the same signature in Impute Cell
Fractions, or `scref_neftel4_cibersortx.txt.gz` in Create Signature Matrix, with S-mode batch correction (full-length
single-cell reference against bulk poly-A RNA-seq).

## Caveats

- U251 is one clonal cell line; the states are programs a cell line can occupy, not the cell types of a patient
  tumour. Fractions say which programs dominate the xenograft, and how LITT recurrence shifts them.
- The control brains carry 0.3–4.9 % human-assigned reads; their fractions are noise and are reported only as such.
- 43 of the 4,916 adult malignant cells have no state scores in the metadata. The first run (job 40467465) left them
  unlabelled but inside their tumour's pseudo-bulk (0.9 % of cells); the script now drops them.
