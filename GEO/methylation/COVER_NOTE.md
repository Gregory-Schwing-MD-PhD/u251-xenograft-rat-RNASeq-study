# Cover note for the methylation submission, and the corrections to ask for at the same time

Greg sends this; GEO submissions and corrections are the submitter's to make. Send to **geo@ncbi.nlm.nih.gov** after
the eight arrays are uploaded, quoting the new Series' temporary accession.

Items 2–5 are corrections to records that are **already public**. They are independent of the new submission, but
they go in the same message because they are all one curator's job on one study, and item 2 is the one that makes the
methylation data findable from the RNA side at all.

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
> **1. Please combine the new Series with GSE338105 under a SuperSeries.**
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
> Note that a BioSample accession is minted for a sequencing sample but not for an array sample, so there is no
> shared accession for a reuser to join on; the correspondence exists only in the metadata. That is why I am asking
> for it in the summary as well as in the array samples' own characteristics fields.
>
> **3. The SRA study descriptor for PRJNA1492510 is stale and says the data are embargoed.** Its STUDY_TITLE still
> reads "GEO accession GSE338105 is currently private and is scheduled to be released on Jul 08, 2027", with a
> matching abstract, while every run is public and GEO says Public on 17 July 2026. Anyone who reaches the study
> through SRA rather than GEO will conclude the data are embargoed for another nine months. Please replace the title
> and abstract with the Series' own.
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
