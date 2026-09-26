# Sample key: which library came from which rat

Established 2026-09-26 from the Henry Ford lab's own records. Where this file disagrees with an older label in this
repository, this file is right. The older labels were "Control (failed graft)" for IL64B and "Control (procedural)" /
"procedural controls" / "rat-brain controls" for N168B and N269B, as if those two were separate control animals.
They are not.

## The design in one paragraph

Eight athymic RNU/RNU rats were implanted with U251N cells in the right striatum. Four were not ablated (rats 64,
67, 68, 69: the **primary** arm). Four were treated with MRI-guided LITT and the tumor was allowed to regrow (rats
65, 66, 70, 71: the **recurrent** arm). The design is **unpaired**: every primary tumor and every recurrent tumor
comes from a different animal, and "pre-LITT" in our cohort names means "not ablated", not "sampled before
ablation". At sacrifice the lab also took the **contralateral hemisphere** of two primary-arm rats (68 and 69). One
culture of the same cells was sequenced as well. Ten RNA libraries exist; our differential-expression comparison is
**three primary versus three recurrent**.

## Per library

Human % is the xengsort graft fraction (reads assigned to human only), from `ANALYSIS/metadata_full.csv` and the
xengsort logs. S number is the Illumina sample number in the FASTQ file name.

| Library | S# | Rat | Lab label | Treatment | Tissue | Human % | Role in our analysis |
| ------- | -- | --- | --------- | --------- | ------ | ------: | -------------------- |
| IL67B | S23 | 67 | IL-67, Primary | U251N implanted, not ablated | tumor | 42.52 | primary arm, DE |
| IL68B | S24 | 68 | IL-68, Primary | U251N implanted, not ablated | tumor | 45.38 | primary arm, DE |
| IL69B | S25 | 69 | IL-69, Primary | U251N implanted, not ablated | tumor | 64.40 | primary arm, DE |
| IL64B | S21 | 64 | IL-64, Primary | U251N implanted, not ablated | implanted hemisphere; the sampled tissue held almost no tumor | 0.33 | excluded from DE; control sample |
| IL66B | S22 | 66 | IL-66, Recurrent | U251N implanted, LITT, regrowth | recurrent tumor | 29.44 | recurrent arm, DE |
| NL70B | S26 | 70 | IL-70, Recurrent (IL70B on the DNA sheet) | U251N implanted, LITT, regrowth | recurrent tumor | 42.67 | recurrent arm, DE |
| NL71B | S27 | 71 | IL-71, Recurrent (IL71B on the DNA sheet) | U251N implanted, LITT, regrowth | recurrent tumor | 33.05 | recurrent arm, DE |
| N168B | S30 | 68 (inferred) | N1, Control, "Contralateral" | U251N implanted in the other hemisphere, not ablated | contralateral hemisphere | 0.55 | excluded from DE; control sample; same animal as IL68B |
| N269B | S28 | 69 | N2, Control, "Contralateral" (DNA sheet: "69B control N2") | U251N implanted in the other hemisphere, not ablated | contralateral hemisphere; holds some tumor cells | 4.90 | excluded from DE; control sample; same animal as IL69B |
| C2B | S29 | none | C2, Culture U251, "Isolated from culture" | U251N in vitro | cultured cells | 86.19 | outside every test; exploratory PCA only |

In the lab key but with no library in our data: **IL-65** (Recurrent) and **C1** (Culture).

Counted by animal: eight rats, nine in vivo libraries. Rats 68 and 69 each contribute two libraries (tumor and
contralateral hemisphere).

## Sources

Nothing below is copied into this repository (it is public). The lab holds the originals.

1. **The lab's RNA-seq key** (2023; `metadata.csv` in Raj Nagaraja's RNA-seq analysis folder). Twelve rows, columns
   Sample ID, Classification, Tissue ID: IL-64, IL-67, IL-68, IL-69 "Primary U251, Tumor grown in rat brain";
   IL-65, IL-66, IL-70, IL-71 "Recurrent U251, Post-ablation tumor in rat brain"; N1 and N2 "Control,
   Contralateral"; C1 and C2 "Culture U251, Isolated from culture".
2. **The methylation array core's sample sheet** for the same lab, the same animals and DNA (plate 1724, "Raj-8";
   one BeadChip, 205648300021). Eight positions: IL66B, IL67B, IL68B, IL69B, IL70B, IL71B, "69B control N2", C2B. This
   is the only place where N2 is tied to a rat number (69). It also shows that the lab calls rats 70 and 71 "IL";
   "NL" appears only in the RNA-seq file names.
3. **Nagaraja et al., J Neurosurg 2026** (methods, "RNA-Seq studies"): a separate group of eight rats bearing U251N
   tumors; the tumor was not ablated in four and was ablated and monitored for recurrence in the other four; after
   slicing, "any remaining tumor and contralateral brain tissues were dissected". Its abstract says primary n = 4 and
   recurrent n = 4; its RNA-seq results say four primary and three recurrent were profiled.
4. **The FASTQ read headers** (checked 2026-09-26 on the xengsort-sorted reads, which keep the original read names):
   every library begins `@VH00948:4:AAAMFFWHV:1:...`, i.e. one instrument (VH00948), one run, one flow cell
   (AAAMFFWHV), lane 1, reads of 101 nt. All ten libraries were sequenced together (S21 to S30). A "VH" instrument
   id is a NextSeq 1000/2000 (NextSeq 500 ids begin "NB"), which agrees with the GEO record (NextSeq 2000) and the
   sequencing core, and not with the JNS 2026 methods ("Nextseq500") or a NovaSeq 6000, 150 bp description.
5. **The xengsort classification logs** for the human % column.

## What is not known

- What "IL", "NL" and the "B" suffix stand for. C2B, a culture, also carries "B", so "B" is not an animal code.
- That N168B comes from rat 68 is **inferred** from the naming pattern N\<k\>\<rat\> (N1-68, N2-69). Only N2 = rat
  69 is written down (DNA sheet).
- Why IL64B's sample held almost no tumor: a graft that did not take, or tissue taken beside the tumor. The lab key
  records rat 64 as a primary with "Tumor grown in rat brain". We do not call it a failed graft.
- The interval from LITT to harvest for each recurrent animal. The JNS 2026 paper's imaging cohort was sacrificed 2
  or 4 weeks after LITT; for the RNA-seq rats it says only "monitored for recurrence", and no file we hold records
  the interval per animal.
- Why rat 65 and culture C1 have no RNA library.
- Whether any other contralateral hemispheres (rats 64, 67, 65, 66, 70, 71) were collected.

## Consequences for the analysis

1. **The comparison is 3 v 3, not the paper's 4 v 4.** Primary: IL67B, IL68B, IL69B. Recurrent: IL66B, NL70B, NL71B.
   The fourth primary (IL64B) has no tumor transcriptome to compare (0.33 % human reads) and the fourth recurrent
   (rat 65) has no library. The JNS 2026 abstract's 4 v 4 describes the animals, not the libraries.
2. **The controls are not independent of the primaries.** N168B and N269B are the other hemispheres of rats 68 and
   69, whose tumors are IL68B and IL69B, two of our three primaries. Any test that treats "controls" and "tumors" as
   independent groups (for example tumor versus control brain in the rat host stream) is partly paired and should
   carry the animal. Anything estimated from the controls, such as RUVSeq's factors of unwanted variation, can carry
   signal specific to rats 68 and 69, both in the primary arm; whether this contributes to the RUVSeq factors
   separating the arms has not been tested.
3. **N269B holds tumor cells.** It is 4.9 % human, with human Y-linked transcripts of the male U251 line; rat 69 is
   the most tumor-rich animal (IL69B 64.4 % human). It is not pure rat brain. Contamination estimates that use it as
   pure rat are biased; the sensitivity analyses drop it (`ANALYSIS/graft_relation/SUMMARY.md`,
   `ANALYSIS/human_cohorts/PLAN.md`).
4. **IL64B is excluded from differential expression** and used as a control sample: it came from the implanted
   hemisphere of a non-ablated rat, so it may carry needle-track and graft-site changes the contralateral samples do
   not.
5. **None of the three control samples is a tumor-free animal.** All three come from tumor-bearing, non-ablated rats.
   "Control" in this repository means "a sample with almost no human reads, run through the identical pipeline",
   which is what the contamination checks need; it does not mean a naive or sham-operated rat.
6. **The labels used by code do not change.** The `Classification` value `Control` in `metadata_base.csv`,
   `metadata_full.csv` and `host/metadata_host.csv` stays `Control`. Only the descriptive `cohort` strings in
   `sample_cohorts.csv` changed: IL64B is `Control (no tumour in sample)`, N168B and N269B are
   `Control (contralateral hemisphere)`.

## Dates: what bounds the harvests, and what is not recorded

Searched 2026-09-26: the JNS submission and supplement (every figure page and embedded image, PDF annotations), the
Acta Neurochirurgica 2021 paper, the lab's slide decks (text, notes, embedded media and file properties), the
methylation IDATs and sample sheet, the RNA-seq analysis folder and the sorted reads' headers.

- **No file gives an implantation, ablation or sacrifice date, the LITT-to-harvest interval (2 or 4 weeks), or the
  order of sacrifice for any of rats 64-71.** The RNA-seq rats are not the animals in any MRI or histology figure of
  the paper (a "separate group"; fresh-frozen brains, whereas the figure brains were perfused with tracers and fixed).
- **Upper bound:** the methylation chip (205648300021) holding rats 66, 67, 68, 69, 70, 71, rat 69's contralateral
  hemisphere and C2B was scanned on 7 April 2022 (iScan RunInfo in every IDAT; the core's QC report is dated
  11 April 2022). Those animals were harvested before then. Rats 64 and 65 are not on the chip.
- All ten RNA libraries ran together (NextSeq 2000 VH00948, run 4, flow cell AAAMFFWHV, lane 1); the run date is not in
  the read headers. IL64B already appears in the lab's analysis deck created 2022-06-11.
- **Protocol, not per-animal record** (the paper's imaging cohort; the paper does not say the RNA-seq rats followed
  it): MRI about 2 weeks after implantation, LITT the next day, recurrence visible after about a week, sacrifice 2 or
  4 weeks after LITT; unablated rats survive about 3 weeks after implantation.
- **The lab's MRI files encode animal and date.** A slide in the lab's ablation-video deck embeds
  "IL360705(Ablation) - 20180705_082218_IL36_..." (rat IL36, 5 July 2018): "IL" is the lab's running animal prefix in
  its MRI records, so the MRI study folders for IL64-IL71 would carry the dates.
- **To obtain from the lab:** the MRI study folders for IL64-IL71, the Visualase logs for rats 65, 66, 70 and 71, the
  IACUC #1509 surgery and euthanasia log, and the USC Molecular Genomics Core submission and receipt forms.
