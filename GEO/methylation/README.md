# Depositing the DNA methylation arrays, and how they relate to GSE338105

## How GEO handles a study that used two assay types

GEO keeps one assay type per Series: a Series cannot mix an array platform with a sequencing platform. A study that
did both is held together by a **SuperSeries** — a record whose members are the individual **SubSeries**, one per
assay. That is the standard route for multi-omics deposits (for example GSE209878, RNA-seq + WGBS + ATAC-seq, and
GSE26168, RRBS + RNA-seq), and it is what this study needs.

So:

1. The RNA-seq is already deposited: **GSE338105**, public since 17 July 2026, ten samples GSM9866604–GSM9866613,
   BioProject PRJNA1492510, platforms GPL30173 (human) and GPL37196 (rat).
2. Submit the eight EPIC arrays as **a new Series** using the files this directory prepares.
3. In the submission's cover note, ask GEO to **combine the new Series with GSE338105 under a SuperSeries**. GEO's
   curators create the SuperSeries; it then appears on both records as `!Series_relation = SuperSeries of: …` /
   `SubSeries of: …`. Do not try to add array samples to GSE338105 itself.

The two records are also joined sample by sample, which is the part that matters for reuse: every array row in the
metadata carries `characteristics: matched RNA library` and `characteristics: matched RNA GEO sample`, so a reader
can pair `IL67B_EPIC` with `GSM9866608` without guessing. Note the naming, which is a real trap and is spelled out
in every file here: the arrays called **IL70B, IL71B and N2** are the libraries **NL70B, NL71B and N269B**.

## What to run

```bash
python GEO/methylation/make_metadata.py         # builds metadata_methylation.xlsx from the sample sheet
sbatch -D ~/logs GEO/methylation/stage_submission.sbatch   # stages the IDATs and builds the processed matrices
```

`stage_submission.sbatch` reads the sixteen raw IDATs out of the laboratory archive
(`u251_meth/raw/DNAmethylation-…zip`) rather than out of a Nextflow `work/` tree — the study's only other copy was
in one, and those are swept. It refuses to continue unless all eight positions are present in both channels. It then
builds the processed matrices from those same IDATs with a **pinned** minfi container, so the deposited matrix is
reproducible from the deposited raw files, and writes md5 checksums for everything.

## The checklist

| item | state |
|---|---|
| 16 raw IDATs (8 arrays × Grn/Red), chip 205648300021, R01C01–R08C01 | **exists** in the laboratory archive; staged by the job |
| Sample-to-array mapping | **exists** (`array_samplesheet.csv`; the laboratory's own sheet is `1724 (Raj-8).xlsx` inside the archive) |
| Processed matrix (noob betas) + detection P values | **made by the job** from the IDATs with pinned minfi |
| Metadata workbook (STUDY / SAMPLES / PROTOCOLS) | **made by `make_metadata.py`**, reusing the accepted RNA wording, contributors and contact |
| Extraction, bisulfite, hybridisation protocol detail | **PENDING** — from the laboratory and the core facility |
| Exact EPIC manifest B-version behind GPL21145 | **PENDING** — confirm at submission |
| Upload and the SuperSeries request | **Greg** — GEO submissions are the submitter's to make |

## What is deliberately not submitted

The laboratory's own `beta_values.rda` and the core facility's intensity CSVs are derived files from a different
processing route; the deposit carries matrices this repository can rebuild from the raw IDATs instead. Nothing from
`MBR/` (unpublished manuscript drafts) goes anywhere near GEO.
