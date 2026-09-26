# Deconvolution result, 26 Sep 2026 (CIBERSORT-core re-implementation; the official CIBERSORT was not available)

Reference: Neftel 2019 SS2, adult malignant cells, four states. Mixture: salmon TPM of the human reads. Full tables in
`results/crosscheck.txt`; fractions in `results/fractions_nusvr_*.tsv`.

| state | primary mean | recurrent mean | difference | Welch p | rank agreement with September GSVA (Spearman, 6 tumours) |
|---|---|---|---|---|---|
| MES | 0.52 | 0.39 | -0.13 (-0.10 without IL68B) | 0.16 | +0.77 (all-cells signature +0.89) |
| AC | 0.22 | 0.27 | +0.05 | 0.19 | -0.26 (all-cells +0.20) |
| NPC | 0.26 | 0.32 | +0.06 | 0.24 | +0.49 (all-cells +0.60) |
| OPC | 0.00 | 0.02 | +0.02 | 0.49 | not defined (near zero everywhere) |

Main signature (confident cells). The in vitro culture C2B reads MES 0.22, AC 0.43, NPC 0.36; the six xenografts read
MES 0.30-0.56.

What holds and what does not:
- The fit is weak: the correlation between the reconstructed and the observed profile is 0.19-0.26 in every tumour
  (CIBERSORT's own diagnostic; its p-value, 0.012-0.015 here, only rejects "no reference state present", which a
  pure-tumour sample always does). Four patient-derived states explain little of a cell line's profile, and this
  re-implementation has no platform correction (CIBERSORTx B-mode).
- MES is the one state both methods rank the same way (Spearman +0.77 / +0.89 against GSVA), and both put it lower in
  recurrence; the fall survives dropping IL68B. Not significant with three tumours per arm.
- AC does not agree between methods: GSVA had AC falling in recurrence (p 0.007); deconvolution has it slightly up.
  The leave-one-tumour-out check showed this reference confuses AC with MES (AC over-called by 9 points, MES
  under-called by 8), so the AC fraction is not evidence either way.
- Growing in rat brain appears to shift U251 toward MES relative to culture (one culture sample).

Open: the official CIBERSORT / CIBERSORTx run (B-mode) on the same inputs; the non-malignant control signature (still
permuting); `SINGLE_CELL_REFERENCES.md` for GBM-CARE (primary vs recurrent) and GSE134470 (U251 culture vs brain).
