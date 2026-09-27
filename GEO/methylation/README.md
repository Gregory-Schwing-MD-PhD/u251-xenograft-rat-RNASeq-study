# Depositing the DNA methylation arrays, and how they relate to GSE338105

## How GEO handles a study that used two assay types

A study that used an array and a sequencer is held together by a **SuperSeries** — a record whose members are the
individual **SubSeries**, one per assay.

Checked, rather than assumed, on 27 September 2026: of the ten GEO Series that carry both a methylation-array assay
type and a sequencing assay type (Entrez `gds`, `"Methylation profiling by array"[DataSet Type] AND "Expression
profiling by high throughput sequencing"[DataSet Type] AND gse[Entry Type]`), **ten out of ten are SuperSeries** —
every one has `!Series_relation = SuperSeries of:` lines in its family SOFT record. No counterexample turned up, so
the practical route is not in doubt, whatever GEO's curation rule is stated to be.

So:

1. The RNA-seq is already deposited: **GSE338105**, public since 17 July 2026, ten samples GSM9866604–GSM9866613,
   BioProject PRJNA1492510, platforms GPL30173 (human) and GPL37196 (rat).
2. Submit the eight EPIC arrays as **a new Series** using the files this directory prepares.
3. In the submission's cover note, ask GEO to **combine the new Series with GSE338105 under a SuperSeries**. GEO's
   curators create the SuperSeries; it then appears on both records as `!Series_relation = SuperSeries of: …` /
   `SubSeries of: …`. Do not try to add array samples to GSE338105 itself.

## The part GEO will not do for you: saying which two samples are the same piece of tissue

A SuperSeries links the two *Series*. It does not link the two *samples*, and neither does anything else at NCBI.
Checked in a real multi-omics SuperSeries (GSE335256): its sequencing SubSeries carries one
`!Sample_relation = BioSample: SAMN...` per sample, and its array SubSeries **carries none at all**. GEO mints a
BioSample for a sequencing sample and not for an array sample, so there is no shared accession for a reuser to join
on. Whatever links the arms has to be written into the metadata by the submitter.

That is made worse here by the names, which we are deliberately **not** changing. Two facilities named the same
animals differently:

| the same piece of tissue | RNA side (GSE338105, SRA, the sequencing core's file names) | array side (this Series, the chip, the core's sheet, every analysis in this repository) |
| --- | --- | --- |
| rat 70, recurrent | `NL70B` | `IL70B` |
| rat 71, recurrent | `NL71B` | `IL71B` |
| rat 69, contralateral hemisphere | `N269B` | `N2` |

The other five (`IL66B`, `IL67B`, `IL68B`, `IL69B`, `C2B`) carry the same name in both arms. Neither vocabulary is
renamed: the RNA names are public and minted in SRA, and the array names are what the chip, the core's sheet and
every analysis in this repository use. `IL` is in fact the laboratory's own animal prefix — its RNA key says IL-70
and IL-71, and `NL` appears only in the sequencing core's file names — so the public record is the one carrying the
odd label, and it is the one that cannot be changed.

So the correspondence is stated four times over, and never left to be inferred from a name:

1. **In each sample's title**, where a reader cannot miss it:
   `U251N recurrent DNA methylation [IL70B = RNA library NL70B]`.
2. **In five characteristics fields per sample** — matched RNA library, GEO sample, SRA experiment, BioSample, and
   the array label on the core's sheet — so the join is machine-readable from the sample record alone.
3. **In the Series summary**, which names all three two-name pieces and says that IL64B and N168B have no array.
4. **In `multiomics_sample_key.tsv`**, deposited as a supplementary file: one row per RNA library with its GSM, SRX,
   SRR, BioSample, and the matched array's label, chip position and column name in the beta matrix. It is generated
   by `python GEO/make_sample_key.py`, which cross-checks the array sheet against the accession map in both
   directions and refuses to write the file if they disagree on any pairing, rat number or chip position.

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
