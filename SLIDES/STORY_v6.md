# Abstract 418 talk, v6: what changed from v5, and how to give it in six minutes

CNS 2026, SSTU02, Monday 2 November 2026, 7:00–7:06 AM. Built by `SLIDES/16_build_deck_v6.py` as
`CNS2026_Schwing_Abstract418_CNStemplate_v6.pptx`: **53 slides — 28 in the talk, 25 backups.** v4's slides and its
builder are untouched; v6 opens the v4 deck, adds, and re-orders.

## Greg's five instructions, and what each one did to the deck

1. **"I definitely want to see slide 14 in the talk."** The contralateral observation is now *three* slides in the
   talk (20, 21, 22), not one: the identity, the dilution null, and what a surgeon can take from it.
2. **"Shift more of the slides from background into the actual talk, in a logical order."** Ten slides that were
   backups in v5 are in the talk in v6: what each lesion is made of, less tumour on both platforms, the composition
   check across ten splits, the four views of the three categories, the invasion programme, MGMT, copy-number
   amplitude, the dilution null, and the sizing curve. The talk runs 4 → 28 slides in a logical arc.
3. **"It's a 6 minute talk."** The deck marks a **6-minute core of 14 slides (303 s)**; the other 14 talk slides are
   the expanded version, in the order they would be given, with a cut list. Nobody has to decide on stage.
4. **"Don't definitively say the instrument contamination is what I'm measuring unless it's definitely true."** It is
   not true, so slide 22 gives the two readings as alternatives this design cannot separate, states what does not
   depend on the choice, and puts the instrument sentence in the conditional. The canine 30 % figure is labelled as
   someone else's measurement.
5. **"Figures too small, too much text, I can't interpret a 4-subplot figure."** Every results slide is now ONE
   panel. `SLIDES/15_slide_panels.py` cuts every multi-panel figure into single panels (the matplotlib ones by each
   axes' own bounding box, the R-drawn ones by their white gutters), and `v6_layout.py` puts each panel wherever it
   is largest — above the text if it is wide, beside it if it is squarish. The builder refuses to finish if a panel
   covers less than 45 % of the content area without being limited by the slide itself, and it reserves text space by
   **wrapped** lines, which is what kept v5's last line clear of the reference footer.

## The 6-minute core (303 s)

| # | slide | s |
|---|---|---|
| 1 | Title | 10 |
| 2 | LITT kills the core; recurrence grows from the margin | 22 |
| 3 | The design: eight rats, RNA and methylation | 18 |
| 4 | A human tumour in a rat brain: every read sorted by species | 20 |
| 6 | Recurrent lesions hold less tumour — on both platforms | 20 |
| 8 | After LITT the regrown tumour turns its ribosomal-protein genes down | 28 |
| 9 | Leave any tumour out: the direction holds, the FDR does not | 25 |
| 10 | Recurrence, or less tumour in the sample? | 25 |
| 11 | How much of the difference is simply less tumour? (C = 0.47) | 28 |
| 17 | Methylation: nothing survives correction | 22 |
| 20 | U251N in the opposite hemisphere | 25 |
| 22 | Two readings, and what a surgeon can take from it | 25 |
| 25 | What this means | 25 |
| 28 | Acknowledgments | 10 |

Fifty-seven seconds of slack against the six minutes: enough for the ablation video on slide 4 and one question.

## The full talk order (28 slides)

1 title · 2 margin · 3 design · 4 model and read sorting · **5 what each lesion is made of** · **6 less tumour on
both platforms** · 7 PCA and volcano · 8 GSEA and leading edge · 9 leave-one-out · 10 the tumour-content confound ·
**11 the exact split** · **12 the composition term behaves** · **13 the three categories in four views** ·
**14 where the mesenchymal shift comes from** · 15 subtype scores · **16 the invasion programme** ·
**17 methylation** · **18 MGMT** · **19 copy-number amplitude** · **20 the opposite hemisphere** ·
**21 not the core diluted** · **22 two readings** · 23 ciclopirox · 24 why ciclopirox · 25 what this means ·
**26 what the next experiment costs** · 27 limitations · 28 acknowledgments. (Bold = new or promoted in v6.)

**Cut order if it runs long:** 12, 19, 5, 18, 21, 13, 16, 26, 14, 7. **Never cut:** 10, 11, 20, 22, 25.

## The backups (25)

Everything not promoted, each now one panel: v4's B1–B7 unchanged, then read composition, the three fraction
estimators, DE counts across the splits, the second mesenchymal signature, the invasion score against tumour
content, the floor gate that stopped the contralateral invasion test, the detection limit for a shared subclone,
fluctuating CpGs, three array-against-RNA panels, the contralateral copy-number dosage check, the titration that
failed, the similarity ranking, the four named arms, malignant-state composition, and two sizing panels.

## Build

```
python SLIDES/13_v5_figures.py      # manuscript panels at slide type
python SLIDES/15_slide_panels.py    # single-panel crops -> figures_cns/v6
python SLIDES/09_build_deck_v4.py   # only if the v4 deck is missing
python SLIDES/16_build_deck_v6.py   # prints each slide's panel area and refuses a panel that is too small
```

Then export and look: PowerPoint COM (`$pres.Export($dir, "PNG", 1920, 1080)`); LibreOffice is not installed on this
machine. Looking at the render is part of the build — it caught the statement slide losing its text to a positional
argument, and the wrapped-line overflow, neither of which the numbers showed.

## Still open

- The panel crop of `fig_a_nullA_p1` clips the last x-axis tick label (0.16 reads as 0.1); harmless on a slide,
  worth re-cropping before the talk.
- Panel letters (a, b, c …) survive on the R-drawn crops; they mean nothing on a slide but are not wrong.
- Whether slide 22 should be given at all before the manuscript is submitted is Greg's call.
