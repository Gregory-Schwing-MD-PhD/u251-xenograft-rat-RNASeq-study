# Cover note for the methylation submission, and the corrections to ask for at the same time

Greg sends this; GEO submissions and corrections are the submitter's to make. Send to **geo@ncbi.nlm.nih.gov** after
the eight arrays are uploaded, quoting the new Series' temporary accession.

Items 2–5 are corrections to records that are **already public**. They are independent of the new submission, but
they go in the same message because they are all one curator's job on one study, and item 2 is the one that makes the
methylation data findable from the RNA side at all.

---

## Why this is an update and not a fresh submission

Greg asked on 2026-09-27 whether, since nothing published cites GSE338105, it would be cleaner to submit everything
again and have the old Series deleted — renaming the samples in the process. Checked against NCBI's own documentation
and records, the answer is no, on five counts.

1. **The premise about citation is correct, and stronger than assumed.** The published Nagaraja et al., J Neurosurg
   145:364–377, 2026 contains no accession of any kind and **no data availability section at all** across its 14
   pages. It could not have cited the Series in any case: GSE338105 was submitted 2026-07-09, **34 days after** the
   article was e-published on 2026-06-05, and the deposit is not the paper's cohort (the paper is 8 rats, 4 versus 4;
   the deposit is 6 tumours, 3 versus 3, plus 4 non-tumour libraries). Nine independent indexes return zero for
   GSE338105, PRJNA1492510 and SRP716631 — Europe PMC full text and its ACCESSION_ID field against a working
   control, PubMed, PMC, Crossref, DataCite, Scholexplorer, BioStudies and Zenodo — both SNO abstracts carry no
   accession, and the paper has no citations.
2. **But deletion is not available for this reason.** GEO documents that submitters cannot delete their own records;
   only GEO staff can, by email, and the page states that "updating records is preferable to deleting records".
   NCBI's documented grounds for *withdrawal* are harm-based — consent, privacy, security, fraud, unauthorised
   submission — and do not include a wish to restructure or rename. NCBI further names "an update was provided as a
   new submission instead of as an update" as a **duplicate-data defect**, resolved by secondary-accession mapping or
   by suppressing the original.
3. **Deleting would not give a clean slate.** The metadata is already mirrored outside NCBI: all ten samples are live
   and PUBLIC in EBI BioSamples, the ten runs are in ENA and DDBJ, GSE338105 is embedded as an external id in the
   mirrored BioProject, OmicsDI has harvested both records, and the runs are served from the AWS Open Data bucket.
   The one withdrawal traceable end to end (CONVERGE, PRJNA289433) shows the shape of the result: six years on, the
   BioProject and SRA study accessions still resolve, carrying "Raw data for this project was withdrawn by request of
   the submitter", while the runs are gone from both NCBI and ENA.
4. **The rename is impossible by any route, so it cannot motivate a resubmission.** The rat identifiers are not in
   any SRA alias or `library_name` — all ten read GSM9866604–GSM9866613 — but they **are** in each Run's original
   FASTQ filenames, and SRA states that once a Run is loaded its files "cannot be replaced nor filenames changed". A
   fresh deposit would re-upload the very same filenames. `IL70B`/`NL70B` survives either way.
5. **Everything actually wanted is a documented in-place edit.** GEO: "You may perform updates and edits at any time
   to any of your submissions", with a self-serve UPDATE button per record and an emailed accession-to-value table for
   non-uniform per-sample edits. The only documented refusals are Platform data tables and deleting published
   accessions. Adding supplementary files and new Samples to an already-public Series is observed in practice
   (GSE206375, public 2022-12-31, then 11 new supplementary files and 22 new Samples in February 2023), and an
   after-the-fact SuperSeries is routine (GSE347440, combined over SubSeries public for 610–638 days).

One constraint to respect when filling the array metadata: GEO requires a Sample title to be unique across every
Sample the submitter has ever deposited, and a Series title across every Series, so the array samples cannot reuse the
RNA titles verbatim. The `_EPIC` suffix already satisfies this.

If a restructure is ever genuinely needed, the route that avoids re-uploading about 100 GB of FASTQ is GEO's
`seq_template_with_sra_accessions.xlsx`, whose SAMPLES section takes `*BioProject`, `*BioSample` and `*SRA Experiment
or Run` per row and builds a new Series over raw data already in SRA. What is **not** documented is whether it may
point at SRA records GEO itself created for an earlier Series — that would have to be asked.

---

## Draft message

> Subject: GSE338105 — companion methylation Series, SuperSeries request, and four corrections
>
> Dear GEO curators,
>
> I have submitted a new Series of eight Illumina Infinium MethylationEPIC arrays (temporary accession
> **\PENDING{temporary accession}**). They are the DNA methylation arm of the study whose RNA-seq is already
> deposited as **GSE338105** (public 17 July 2026; BioProject PRJNA1492510). The DNA and the RNA were taken from
> different aliquots of the same pieces of tissue.
>
> I would be grateful for the following.
>
> **1. Please combine the new Series with GSE338105 under a SuperSeries.** I understand this is done after the fact
> as a matter of routine; GSE347440 was combined over SubSeries that had been public for more than 600 days.
>
> **1b. If it is possible, please attach these eight array samples to the ten BioSamples already registered for the
> RNA-seq arm (SAMN61491943–SAMN61491952), rather than leaving them without one.** The DNA and the RNA came from
> different aliquots of the same eight specimens, and a shared BioSample accession is the only mechanism that would
> let a reuser join the two arms without reading metadata. I appreciate this is unusual for an array submission — in
> a survey of 26,314 array samples deposited since 2010 I found only 17 carrying a BioSample, none after 2014 — but
> NCBI's BioSample FAQ does state that "multiple assay types … may reference the same BioSample if appropriate", and
> GSE18927 contains array samples sharing BioSamples with sequencing samples. If it is not possible, the per-sample
> characteristics fields below carry the pairing instead, and that is sufficient.
>
> **2. Please add a sentence to GSE338105's summary pointing at the methylation Series, and giving the sample
> correspondence.** Without it there is no way to reach the methylation data from the RNA record. Suggested wording:
>
> > DNA methylation arrays on the same tissue are deposited as \PENDING{array Series accession}. Eight of these ten
> > libraries have a matched array. Three of them are named differently in the two Series, because the array core
> > and the sequencing core named the animals differently: RNA libraries NL70B, NL71B and N269B are arrays IL70B,
> > IL71B and N2. The other five (IL66B, IL67B, IL68B, IL69B, C2B) carry the same name in both. IL64B and N168B have
> > no array. The per-sample correspondence, with accessions, is the supplementary file
> > multiomics_sample_key.tsv on the array Series.
>
> **2b. Please also add one characteristics field to the ten existing RNA-seq samples in GSE338105, so that the join
> works from either side.** As deposited, a reader who starts at the array Series can reach the RNA library, but a
> reader who starts at GSE338105 has nothing pointing back. I will send a tab-delimited accession-to-value table in
> the format your update documentation describes:
>
> > `characteristics: matched methylation sample` — for GSM9866608 `IL67B_EPIC`, GSM9866609 `IL68B_EPIC`,
> > GSM9866610 `IL69B_EPIC`, GSM9866611 `IL66B_EPIC`, GSM9866612 `IL70B_EPIC`, GSM9866613 `IL71B_EPIC`,
> > GSM9866607 `N2_EPIC`, GSM9866604 `C2B_EPIC`; and for GSM9866605 and GSM9866606, `none (no array was run)`.
>
> A characteristics field is the right carrier for this: it is the only extensible per-sample field, it is typed in
> MINiML rather than free text, GEOquery pre-parses it into a named column, and it is one of the two keys refine.bio
> harmonises on. Sample titles are not a usable substitute — across four real multi-omic SuperSeries the title
> overlap between the two arms ranged from 0 % to 100 %.
>
> **3. The SRA study descriptor for PRJNA1492510 still says the data are embargoed.** Its STUDY_TITLE reads "GEO
> accession GSE338105 is currently private and is scheduled to be released on Jul 08, 2027", with a matching
> abstract, while every run is public and GEO says Public on 17 July 2026. Anyone who reaches the study through SRA
> rather than GEO will conclude the data are embargoed for another nine months.
>
> I raise this expecting it is not specific to my submission: of 26 recently released GEO Series I checked, 23 carry
> the same kind of stale placeholder in their SRA study descriptor, so this looks like it affects GEO-brokered
> submissions generally rather than mine alone. Please replace the title and abstract with the Series' own if that is
> something you can do per-record; if it needs fixing upstream, please treat this as a report rather than a request.
>
> **4. Two sample descriptions in GSE338105 are wrong, and I would like to correct them.**
>
> - **GSM9866605 (IL64B)** currently carries `treatment: procedural / failed-graft control`, and the overall design
>   says "IL64B failed graft". The laboratory's own key records rat 64 as a primary animal with tumour grown in the
>   brain. What is true of this library is that the piece of tissue sampled held almost no tumour — 0.33 % of its
>   reads were assigned human. Whether the graft failed or the piece was taken beside the tumour **is not recorded**,
>   so please replace both strings with: `no tumour in the sampled piece (0.33 % of reads assigned human); implanted
>   hemisphere of a non-ablated animal`.
> - **GSM9866606 (N168B) and GSM9866607 (N269B)** carry the same `procedural / failed-graft control` string as
>   IL64B, and `cell type: glioblastoma`, although both are contralateral hemispheres. Please set their treatment to
>   `none (U251N implanted in the opposite hemisphere, not ablated)` and their tissue to `contralateral hemisphere`.
>   The three controls are not interchangeable and a reader currently cannot tell them apart. Please also note that
>   N269B is not tumour-free: it is 4.9 % human, and carries Y-linked transcripts of the male U251 line.
>
> **5. Please add `salmon.merged.gene_lengths.tsv` to GSE338105's supplementary files.** The Series currently carries
> `GSE338105_salmon.merged.gene_counts.tsv.gz` only. The differential-expression step needs the matching gene-length
> matrix to give DESeq2 its per-gene effective-length offsets, so as deposited the counts cannot be taken through
> that step without re-quantifying every library. I can upload the file on request.
>
> Thank you,
> \PENDING{signature block}

---

## What to check before sending

| item | state |
|---|---|
| The eight arrays uploaded and the temporary accession in hand | **Greg** |
| `multiomics_sample_key.tsv` among the array Series' supplementary files | built by `stage_submission.sbatch`, staged into `processed/` |
| The stale SRA descriptor still stale | re-check; it may have been refreshed since 27 Sep 2026 |
| `salmon.merged.gene_lengths.tsv` in hand to upload | on the grid under the rnaseq run's `star_salmon/` |
| The wording of item 4 agreed with the laboratory (T.N.N., I.D.) | **ask first** — it corrects a description they supplied |

Item 4 is the one to clear with the laboratory before sending. The other four are ours to assert: they are facts
about the record itself, not about the biology.
