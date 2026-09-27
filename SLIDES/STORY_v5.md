# Abstract 418 talk, v5: what was added to v4, and why

CNS 2026, SSTU02, Monday 2 November 2026, 7:00–7:06 AM. Built 2026-09-27 by `SLIDES/14_build_deck_v5.py` as
`CNS2026_Schwing_Abstract418_CNStemplate_v5.pptx`: **46 slides — v4's 21 unchanged, plus 5 main and 20 backups.**
v4 (`09_build_deck_v4.py`, `STORY_v4.md`) is not edited; v5 opens it, inserts, and saves under a new name.

## Why v5 exists

Between the v4 build and this one, the pre-registered SCIENCE analyses finished: the species-exact composition
decomposition, the methylation side against a measured dilution model, the array-against-RNA cross-validation, the
contralateral case, and the power targets for the scaled experiment. Greg asked for the deck to carry **every
positive result and every pertinent negative**, and to stop being short. It now runs the same 6-minute story with
two new results beats, and everything else is one question away in the backups.

## The running order (new slides in bold)

| # | slide | s |
|---|---|---|
| 1–8 | v4 unchanged: title → margin → design → model → PCA/volcano → GSEA → leave-one-out → "recurrence or less tumor?" | 190 |
| **9** | **How much of the difference is simply less tumor? An exact split** (C = 0.47) | 25 |
| **10** | **Where the mesenchymal shift at recurrence comes from** (host myeloid) | 25 |
| 11–13 | v4 unchanged: subtypes → ciclopirox → why ciclopirox | 85 |
| **14** | **Tumor in the opposite hemisphere — and what it means in the OR** | 35 |
| **15** | **Methylation: nothing survives correction, and what differs is tumor content** | 25 |
| 16–17 | v4 unchanged: what this means → limitations | 55 |
| **18** | **The experiment six animals size** (rats-needed curves) | 30 |
| 19 | v4 acknowledgments | 10 |
| 20–26 | v4 backups B1–B7 | — |
| **27–46** | **new backups B8–B27** | — |

That is about 8 minutes of main slides against a 6-minute slot, so the deck ships with a cut order.

**Cut order (in this order, until it fits):**
1. Slide 15 (methylation) → backup. It is a negative that matters, but B14, B15 and B22 carry it.
2. Slide 10 (mesenchymal) → backup. Slide 9 already makes the composition point.
3. Slide 18 (scaled experiment) → fold its one sentence into slide 17 ("limitations") and move the figure to backup.
4. Then, as in v4: slide 4 to its first and last line; the third line of slides 12 and 13.

**Never cut:** slides 7 and 9 (where the talk earns its trust) and slide 14 (the observation a neurosurgical
audience will remember).

## The two new beats

**Composition, measured rather than estimated (slide 9).** In a human graft in a rat host every read is assignable
by species, so the bulk profile a patient sample would give can be rebuilt and each gene's fold change split
exactly — not by regression — into composition, tumour-cell and host-cell terms. Composition is the largest single
term and still under half (0.47), and only 13 % of the genes a bulk analysis calls differential are
composition-dominant. The pre-registered threshold was 0.50 and it failed; that is reported as the result.

**Tumour in the opposite hemisphere (slide 14).** Rat 69's contralateral hemisphere carries 4.9 % human reads, the
donor genotype at 56 of 59 identity probes, Y dosage 0.88 and the line's copy-number profile, and its transcriptome
is not the core diluted. The design cannot separate cells that crossed in life from core tissue carried across on
the blade — **and the slide says so, because either reading is the finding**: one pass of an instrument left tumour
that sequencing detects and histology would not. The anchor is measured, not asserted: malignant cells were
recovered from gloves or instruments in 14 of 47 canine cancer operations (30 %), 57 % when margins were involved
against 19 % when clear (Orjefelt 2026); iatrogenic seeding after a biopsy or into a distant surgical site is
documented in case reports (Steinmetz 2001; Perrin 1998; Bekar 2010). The scaled design measures the transfer
directly with a cutting-order control (29 pieces bound it at 10 %).

## The twenty new backups

B8 read composition · B9 lesion accounting · B10 fraction on two platforms · B11 the twenty labellings ·
B12 ten splits · B13 four views of the three categories · B14 copy-number amplitude against dilution ·
B15 MGMT · B16 the distal-invasion programme · B17 its floor gate · B18 the detection limit for a shared subclone ·
B19 fluctuating CpGs · B20–B22 array against RNA (dosage, genotype, promoter methylation) · B23 the titration that
failed · B24–B26 the contralateral panels · B27 malignant-state composition.

## Build

```
python SLIDES/13_v5_figures.py      # manuscript panels re-saved at slide type + the R-drawn panels copied
python SLIDES/09_build_deck_v4.py   # only if the v4 deck is missing
python SLIDES/14_build_deck_v5.py
```

`13_v5_figures.py` reads `MBR/nature_draft/figures/` (untracked, as `MBR/ESM_1.xlsx` already is). The builder types
no result number: every value comes from `figures_cns/v5_data/*.json|tsv`, byte copies of the grid outputs, and the
script asserts the load-bearing ones before it draws.

## References added (checked against Europe PMC, 2026-09-27)

21 Shen-Orr 2010 Nat Methods · 22 Spitzer 2025 Nat Genet · 23 Orjefelt 2026 J Small Anim Pract ·
24 Steinmetz 2001 J Neurooncol · 25 Bekar 2010 World Neurosurg · 26 Bady 2016 J Mol Diagn ·
27 Barthel 2019 Nature · 28 Browne 1995 Stat Med. v4's numbering (1–20) is untouched.

## Open, for Greg

- The Shapley attribution follows Shorrocks (J Econ Inequal 2013); that reference is in the notes, not the footer,
  until it is checked.
- Slide 14 says "one pass of an instrument": the dissection record does not state whether a fresh blade was used
  between pieces. If Raj's records say it was, the sentence changes.
- Whether to show slide 14 at all before the manuscript is submitted.
