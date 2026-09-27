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
| Scan protocol | **filled from the data**, not asked for: the IDAT RunInfo block gives scanner N141, iScan Control Software 3.4.8, FPGA 4.0.20, `Extract Algorithm=StandardWithBackground`, decoded 2021-07-11, scanned 2022-04-07 09:24–09:26 |
| EPIC manifest version | **filled from the data**: the core's own `dnamet.R` loads `IlluminaHumanMethylationEPICanno.ilm10b4.hg19`, so EPIC v1.0 **manifest B4, hg19** |
| Per-sample QC | **filled from the core's sheet**: 45 µL DNA per sample, qPCR Ct for HB-313 and HB-365, and 864,815–865,311 probes passing detection at *P* < 0.05 (99.879–99.937 %) |
| Extraction kit, bisulfite kit, core facility name | **PENDING** — the only three fields left; see the email below |
| Upload and the SuperSeries request | **Greg** — GEO submissions are the submitter's to make |

## The three fields that need a human, and the email that gets them

Everything else in `metadata_methylation.xlsx` is filled. These three cannot be recovered from the files: the archive
contains the IDATs, the core's minfi `qcReport` (plots only — no kit or facility text), the sample-array mapping sheet,
the Illumina `.sdf` descriptor, the core's `dnamet.R`, and the processed CSVs. None of them names a kit or the facility.

> Subject: two protocol details for the U251 methylation arrays (plate 1724, "Raj-8")
>
> Hello,
>
> I am depositing the eight Infinium MethylationEPIC arrays from plate 1724 ("Raj-8", chip 205648300021, run April
> 2022) in GEO, and the submission needs three protocol details I cannot recover from the files you sent.
>
> 1. **The DNA extraction kit** used on the eight tissue pieces, and the **DNA mass** corresponding to the 45 µL per
>    sample recorded on the mapping sheet (a concentration is equally good).
> 2. **The bisulfite conversion kit** and the **DNA input mass** taken into conversion.
> 3. **The name of the core facility** that performed the conversion, hybridisation and scanning, as it should be
>    credited in the public record.
>
> For context, so you only have to answer what is missing: I already have the chip and plate identifiers, the well
> positions, the 45 µL volumes, the HB-313 and HB-365 qPCR Ct values, the per-sample detection pass rates, the iScan
> instrument and software versions and the scan date from the IDAT headers, and the EPIC B4 / hg19 annotation from
> `dnamet.R`.
>
> One more, if you happen to know: what the **HB-313 and HB-365 assays** are, so I can describe the pre-conversion
> DNA quality check correctly rather than just quoting Ct values.
>
> Thank you,
> Greg

## What is deliberately not submitted

The laboratory's own `beta_values.rda` and the core facility's intensity CSVs are derived files from a different
processing route; the deposit carries matrices this repository can rebuild from the raw IDATs instead. Nothing from
`MBR/` (unpublished manuscript drafts) goes anywhere near GEO.
