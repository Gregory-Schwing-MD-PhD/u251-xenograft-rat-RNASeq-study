# The contralateral hemisphere as a control: implications and prior art

Deep research, 2026-09-26 (5 search angles, 25 primary sources fetched, 125 claims extracted, 25 verified by three adversarial votes: 18 confirmed, 7 refuted). Question as asked:

> What are the biological and analytical implications of using the CONTRALATERAL (opposite) hemisphere of the same tumour-bearing brain as the control tissue for bulk RNA-seq in an orthotopic glioma xenograft, and what prior art uses or studies this design? Context: human U251N glioblastoma cells implanted unilaterally in athymic RNU/RNU nude rats (outbred); four rats untreated (primary), four ablated with MRI-guided laser interstitial thermal therapy (LITT) and allowed to recur; bulk RNA-seq of the tumour side of each rat plus the contralateral hemisphere of two of the untreated rats; reads species-sorted (human graft vs rat host, xengsort); one contralateral sample carried ~5% human reads including human Y-chromosome reads (tumour cells or carry-over; U251 is male, the rats female), the other ~0.6% (background). Answer with citations (PMID/DOI) and hard numbers: (1) Is the contralateral hemisphere a valid normal-brain control? Evidence that a unilateral glioma or brain lesion changes contralateral gene expression, microglia/astrocyte state, blood-brain barrier, inflammation or metabolism (diaschisis, transcallosal and systemic effects, mass effect), and how large those changes are versus naive or sham animals. (2) How often U251, U87 and patient-derived orthotopic tumours reach the contralateral hemisphere in rodents (corpus callosum migration), how that is detected, and any contralateral effects of LITT or thermal ablation. (3) Prior art that uses contralateral tissue as the within-animal control for transcriptomics or proteomics in glioma, stroke/MCAO, TBI or epilepsy models, and what it reported about adequacy versus sham/naive controls; paired ipsilateral-vs-contralateral designs and the statistics used. (4) What unique biological insight a same-brain control gives (tumour-local vs brain-wide or systemic effects; remote effects of the tumour; pairing that cancels animal-level variation in outbred rats). (5) Analysis consequences: non-independence when the same animals give tumour and control samples (paired or mixed models), and using such samples as negative controls for unwanted-variation removal (RUVSeq RUVs) and for estimating species cross-mapping floors in xenograft sequencing. Prefer primary research and methods papers 2000-2026; say plainly where evidence is thin.

## Summary

The contralateral hemisphere of a unilaterally lesioned or tumour-bearing rodent brain is measurably not a naive brain, so it is a defensible within-animal reference but not a substitute for sham or naive tissue. The only prior art that tested both designs head-to-head for transcriptomics is a rat stroke study: the contralateral subcortex differed from the matching sham hemisphere by 164 DEGs (~6% of the 2,802-DEG ipsilateral response), and swapping sham for contralateral rewrote the lesion-side result in both directions (1,941 shared, 861 sham-only enriched for apoptosis/cell-cycle/neurotransmitter genes, 221 contralateral-only and predominantly immune). In glioma specifically the remote hemisphere carries a tumour-driven myeloid and immunosuppressive program (SB28 mice: contralateral TSPO-PET +18% SUVr, +16% VTr vs sham; CD68, TREM2, TYROBP, GRN, CHI3L1, TGFB1 up on nCounter), so a contralateral-referenced contrast subtracts part of the host inflammatory signature rather than revealing it. Orthotopic glioma cells frequently do reach the opposite hemisphere via the corpus callosum in patient-derived models (6 of 8 HGCC lines), and the confirmed evidence contains no U251-specific exemption and no quantified U251/U87 contralateral burden, so ~5% human reads with human Y reads cannot be dismissed as carry-over on prior-art grounds. Evidence is thin to absent on three asked-about points: contralateral effects of LITT, paired/mixed-model statistics for ipsilateral-versus-contralateral transcriptomics in glioma, and use of such samples for RUVSeq negative controls or species cross-mapping floors.

## Findings

### 1. A unilateral brain lesion leaves the contralateral hemisphere transcriptionally and glially altered relative to sham, at roughly 5-10% of the ipsilateral effect size by DEG count. It is therefore not a naive control, though the perturbation is modest.

Confidence: **high** (vote 3-0 (merged claims 0, 7, 16, 9, 10, 11)).

Rat 90-min transient MCAO, subcortex at 24 h, polyA RNA-seq, n=3/group, >1.5-fold and BH Padj<0.05: contralateral vs matching sham hemisphere (IR-c vs SO-l) = 164 DEGs (96 up, 68 down) against 2,802 DEGs ipsilateral vs sham, i.e. 164/2802 = 5.9%. Breakdown: 114 genes codirectional in both hemispheres, 16 opposite, 34 contralateral-only; the 164 figure is reported in two papers from the same dataset. Independent lesion models agree in direction: 7 d after right motor cortex resection in adult male Sprague-Dawley rats, contralateral motor cortex vs sham gave 70 of 21,106 genes differentially expressed (35 up, 35 down; FDR<0.05, |FC|>=2; 28 rats total, per-group sequencing n not stated), with contralateral GFAP-reactive astrocytes (p=0.001) and microglia (p=0.005) above sham at n=6/group. After focal aspiration lesion of right hindlimb sensorimotor cortex, IBA-1 at 0-500 um was elevated in BOTH hemispheres over time-matched sham at day 3 (group p=0.004, contralateral p=0.0008, no hemisphere or interaction effect), lateralising by day 7 and near baseline by day 28. Every comparator is sham (craniotomy/anaesthesia), not naive, so the gap to a truly unoperated brain is bounded below, not measured.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9266805/ (doi:10.3390/ijms23137308, PMID 35806305)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9834327/ (doi:10.1038/s41598-023-27663-8, PMID 36631528)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12384278/ (doi:10.3390/brainsci15080837, PMID 40867169)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12843123/ (doi:10.3390/life16010142, PMID 41598295)

### 2. Choosing the contralateral hemisphere instead of sham animals systematically rewrites the lesion-side DEG list: it suppresses immune genes elevated bilaterally and adds genes that appear changed only because the reference itself moved.

Confidence: **high** (vote 3-0 (merged claims 1, 2, 17, 3)).

Same rat tMCAO dataset, both controls applied in parallel to the same ipsilateral tissue. Ipsilateral vs contralateral = 2,162 DEGs (1,032 up, 1,130 down); ipsilateral vs sham = 2,802 DEGs (1,390 up, 1,412 down); 1,941 overlap, all codirectional except Parpbp. 861 genes (861/2802 = 31%) were DEGs only against sham, the authors' 'lost results', associated with apoptosis, cell cycle and neurotransmitter response (example Cp). 221 genes (221/2162 = 10%) were DEGs only against contralateral, 'redundant results', predominantly immune-related (example Pla2g3). Authors' conclusion: 'The specific response of the CH transcriptome should be considered when using it as a control in studies of target brain regions in diseases that induce a global bilateral genetic response, such as stroke', with sham samples called 'a rational alternative'. Transfer caveat: the contralateral comparison is within-animal and the sham comparison between-animal, so part of the 2,162-vs-2,802 gap may be variance structure rather than biology; analysis was per-gene t-tests with BH, not a paired or mixed model. One lab, one dataset, n=3/group, one timepoint, stroke not glioma, and the 861/221 split depends on the hard cutoff.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9266805/ (doi:10.3390/ijms23137308, PMID 35806305)

### 3. In glioma specifically, the hemisphere opposite the tumour carries a tumour-driven remote myeloid and immunosuppressive transcriptional program, so a contralateral-referenced bulk RNA-seq contrast subtracts part of the host inflammatory signature rather than revealing it.

Confidence: **high** (vote 3-0 for the expression data (claim 5); 2-1 for the imaging numbers (claim 4)).

SB28 syngeneic glioblastoma in C57BL/6 mice (48 inoculated: 26 tumour, 22 saline sham). Contralateral frontal-pole TSPO-PET in a 2 mm sphere chosen to avoid tumour spill-over: +18% SUVr (P=0.032) and +16% VTr (P=0.0086) vs sham. Two corrections to carry: the test n was 15 vs 14, not 26 vs 22, and both endpoints are pons-normalised, so their agreement is not two independent checks. The orthogonal evidence is expression on snap-frozen contralateral hemispheres (nCounter Immunology panel, n=8 tumour vs 4 sham vs 4 untreated, FDR 0.05): CD68, LAMP1/LAMP2, cathepsins, GRN, AIF1 (IBA1), TREM2, TYROBP, FABP5 and TSPO up, qPCR confirming CD68, SPI1 (PU.1) and HEXB, plus immunosuppression-associated CHI3L1 and TGFB1 raised in distant brain; 40 of the top 50 pathways were cell migration, phagocytosis and immune activation. Sham did not differ from healthy, so the effect is tumour-attributable. Human arm: contralateral TSPO-PET elevated in 60 of 123 sub-regions (FDR-corrected), high contralateral signal an independent survival predictor (HR 2.18, P=0.005; median OS 6.5 vs 13.4 months). Weaker replication in GL261/[18F]GE-180: contralateral uptake 10% higher than sham at day 14 (p=0.02), which those authors themselves caution against over-weighting (pooled sham n=10). Transfer limit: both are syngeneic tumours in immunocompetent mice and the mechanism involves myeloid plus T-cell-dependent immunosuppression, so the magnitude in a T-cell-deficient nude rat carrying a human graft is unmeasured. The 'subtracts part of the signature' step is inference: no located paper analysed contralateral tissue as an RNA-seq control in glioma.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC11474166/ (doi:10.1158/1078-0432.CCR-24-1563, PMID 39150564)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9030822/ (PMID 35453488)

### 4. Two distinct things can ride in a contralateral sample - host-side remote inflammation and physical contamination by graft cells - and they were experimentally separated in SB28, where tumour-cell eGFP appeared only in the tumour hemisphere. But eGFP/qPCR detection is far less sensitive than species-sorted RNA-seq, so that result cannot be used to argue a contralateral xenograft sample is graft-free.

Confidence: **high** (vote 3-0 (claim 6)).

'Relevant eGFP expression could only be detected in the tumor hemisphere', and contralateral hemispheres were validated eGFP-negative by qPCR, which the authors use to conclude the contralateral TSPO signal is myeloid rather than infiltrating tumour. They nonetheless write that contralateral tumour-cell infiltration in patients 'cannot be excluded since autopsy reports and primary human glioblastoma models provided evidence for cross-hemisphere migration of tumor cells along white matter tracts'. The hedge 'relevant' marks a limit-of-detection statement. Mapped onto the U251 design, this supports the conceptual dissection (remote host response is real and does not require graft cells) but gives no quantitative prior for a ~5% human read fraction, which sits far above any eGFP-negative threshold.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC11474166/ (doi:10.1158/1078-0432.CCR-24-1563, PMID 39150564)

### 5. Orthotopic glioma cells reach the contralateral hemisphere often enough that a contralateral sample cannot be assumed graft-free, and the rate is strongly cell-line-intrinsic. The confirmed evidence contains no U251-specific exemption and no quantified U251 or U87 contralateral burden.

Confidence: **medium** (vote 2-1 (claim 14); three competing claims about U251 non-migration and U87/MU20 burdens were refuted 0-3).

Eight patient-derived HGCC glioblastoma stem-cell lines, NOD SCID mice, 5 mice per line, 100,000 cells in 2 uL into striatum at 1.5 mm ML / 0.0 mm AP / 3.0 mm depth, detected on coronal sections with human-specific STEM121 and HuNu: 'Most of them showed extensive invasion, particularly along white matter tracts. For these tumors, we regularly observed migration along corpus callosum or anterior commissure into the contralateral hemisphere.' Six of eight lines invaded white matter and migrated callosally; the exceptions U3013MG and U3054MG stayed locally invasive. Three limits: 'regularly observed' is qualitative with no per-animal or per-line percentage, and the 6/8 mapping is an inference joining the route sentence to the two named exceptions; the paper's own thesis is that route is line-intrinsic and reproducible per line, which argues against generalising; and the host is a NOD SCID mouse with patient-derived stem cells, not serially passaged U251N in an outbred nude rat. Critically, adversarial review REFUTED (0-3) the reassuring claim that serum-cultured U251MG/U373MG form demarcated tumours and do not migrate through the corpus callosum, and also refuted the quantified U87MG (3.9 cells per 5 mg) versus patient-derived MU20 (3898 cells per 5 mg) contralateral burdens and the brain-wide dissemination claim. There is therefore no confirmed prior art either way on U251 crossing the midline: contralateral graft cells are common and line-dependent, and U251N behaviour in rats is unmeasured in this evidence set.

Sources:

- https://www.nature.com/articles/s41598-023-51063-7 (doi:10.1038/s41598-023-51063-7, PMID 38195678, PMC10776844)

### 6. A normal-looking readout from the contralateral side does not establish that it is naive on other axes. Imaging-normal tissue was transcriptionally changed, and blood-flow-normal tissue was metabolically depressed 40-50%.

Confidence: **medium** (vote 3-0 (merged claims 8, 12, 13)).

Stroke: 'We found no pathological changes in the CH of ischemic rats in the MRI data' in the same animals that yielded 164 contralateral DEGs, so MRI-normal is not transcriptionally naive. Glioma metabolism: rats bearing implanted neoplastic glial clones (spherical tumours, mean diameter 6 mm, 2-6 weeks), quantitative double-tracer autoradiography under light barbiturate anaesthesia: 'Flow in the opposite hemisphere was of the same order of magnitude as in normal control rats', yet 'Glucose consumption, in contrast, was distinctly reduced in both hemispheres: in the cortex and putamen, it was 40-50% lower than in normal controls', with conclusion (3) that 'the metabolic rate of glucose is distinctly inhibited in both hemispheres of tumor-bearing animals'. Barbiturate anaesthesia depresses flow and metabolism together, so preserved flow alongside 40-50% lower CMRglc is a genuine uncoupling, not an anaesthesia artifact. Kept at medium because it is a single 1982 study with no replication found, the publisher blocked full text so verification rested on an OpenAlex abstract reconstruction (token-checked against fabrication), the 40-50% figure is reported for both hemispheres jointly rather than contralateral alone, no per-group n is given, and the readout is metabolism not mRNA. Direction is corroborated by human glioma PET (crossed cerebellar diaschisis, ~30% lower rCMRglc contralateral; J Nucl Med 2017;58:768, PMID 27789719).

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9834327/ (doi:10.1038/s41598-023-27663-8, PMID 36631528)
- https://doi.org/10.1038/jcbfm.1982.3 (PMID 7061601)

### 7. There is essentially no published evidence on whether LITT or thermal ablation changes the contralateral hemisphere. The one rat U251 MRI-guided LITT methods paper reports no contralateral, hemisphere, corpus callosum or midline measurement at all, and its authors describe the model as lacking an infiltrative component.

Confidence: **low** (vote 3-0 for the absence (claim 15); a claim that thermal spread was bounded away from the opposite side was refuted 0-3).

Literal counts over the full PMC text (~48.7k characters): 'contralateral' 0, 'corpus callosum' 0, 'hemispher*' 0, 'opposite' 0, 'sham' 0, 'naive' 0, 'non-tumor' 0, 'whole brain' 0; 'midline' occurs twice, only as a stereotactic landmark (incision 2 mm right of midline, burr hole 3.50 mm right of bregma). Reported measurements are tumour and ablation volumes, DCE-MRI Ktrans and plasma volume, CBF, FITC-dextran/Evans blue leakage and H&E/MHC/GFAP histology, all within tumour, ablation zone and peri-ablation rim; the only distant normal tissue sampled was IPSILATERAL parietal cortex, and human-specific MHC staining (the one tool that could have shown graft cells across the midline) was reported only for tumour and peri-ablation fields. Authors' limitation: 'the U251 GBM model employed in this study represented only the space-occupying, contrast-enhancing glioma phenotype ... this model lacks the infiltrative component and genetic variability seen in clinical GBM', adding that 'infiltrating tumor cells may remain even beyond this sublethal temperature zone'. That is a phenotype statement, not a measurement of midline crossing, and the separate claim that temperature 4 mm from the fibre tip stayed sublethal was refuted, so for a LITT-recurrence arm the contralateral hemisphere's status is simply unmeasured in the literature.

Sources:

- https://pubmed.ncbi.nlm.nih.gov/34554269/ (doi:10.1007/s00701-021-05002-y, Acta Neurochir 2021;163(12):3455-3463, PMC8893160)

### 8. The unique value of a same-brain control is that it separates tumour-local from brain-wide change and cancels animal-level variation in an outbred host, but this rationale rests on design reasoning plus one stroke study's side-by-side comparison, not on any glioma transcriptomics paper found.

Confidence: **low** (vote supported indirectly; no confirmed source addresses paired-design transcriptomics in glioma).

The stroke paper is the only located prior art that measured sham and contralateral controls side by side for transcriptomics, and its 861/221 split is the quantified consequence of the choice: a contralateral reference removes bilaterally shared immune response, giving a purer tumour-LOCAL readout, while losing genes whose brain-wide component matters. The glioma paper supplies the complement, that the brain-wide component is real, myeloid and prognostic in patients (HR 2.18). So a contralateral-only design can measure the local contrast but not the remote one, and a sham-only design cannot separate them. The specific pairing argument for outbred RNU/RNU rats (cancelling animal-level genetic and batch variation) is addressed by no confirmed source; with contralateral tissue from only 2 of 8 animals the pairing is in any case unavailable for most of the cohort.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9266805/ (doi:10.3390/ijms23137308, PMID 35806305)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11474166/ (doi:10.1158/1078-0432.CCR-24-1563, PMID 39150564)

### 9. The analysis-side questions asked about - paired or mixed models for non-independent within-animal samples, RUVSeq RUVs negative controls, and species cross-mapping floors in xenograft sequencing - are not covered by any confirmed source in this search. State this as thin evidence rather than filling it in.

Confidence: **low** (vote no claim on these points survived verification).

The only adjacent confirmed observation is methodological: in the stroke study the contralateral comparison is within-animal while the sham comparison is between-animal, and a verifier flagged that the 2,162-versus-2,802 DEG difference may partly reflect that variance structure rather than biology; the authors used per-gene t-tests with Benjamini-Hochberg, not a paired or mixed model. Nothing in the confirmed set addresses RUVSeq RUVs with contralateral samples as negative controls, nor using a contralateral xenograft sample to estimate a species cross-mapping floor. Within the described study the two contralateral samples are the only internal handle: ~0.6% human reads as an empirical background/cross-mapping estimate, and the ~5% sample (human Y reads, female host) as an outlier that cannot be assigned to invasion versus carry-over versus index hopping on this evidence. That is a within-study observation, not a literature-supported method.

Sources:

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9266805/ (doi:10.3390/ijms23137308, PMID 35806305)

## Caveats

Model mismatch is the largest limitation: no confirmed source studied a human line in a T-cell-deficient outbred rat. The quantitative transcriptomic evidence is rat focal stroke (90-min tMCAO), the glioma remote-inflammation evidence is syngeneic tumours in immunocompetent mice (SB28, GL261), the glia evidence is surgical cortical lesions, and the contralateral-invasion evidence is patient-derived stem cells in NOD SCID mice. Effect sizes should not be carried across; direction can be. The stroke numbers (164 DEGs; 2,162 vs 2,802; 1,941/861/221) come from one laboratory and essentially one dataset, n=3 per group, a single 24 h timepoint and hard >1.5-fold cutoffs, so they are imprecise point estimates with no independent replication. Nearly every comparator in this evidence base is SHAM, not naive, so "not naive" strictly means "differs from surgery-alone" - a conservative direction, but never measured against unoperated animals. The contralateral glial response is time-dependent (IBA-1 bilateral at day 3, lateralised by day 7, near baseline by day 28), so timepoint choice changes the answer. Weak spots in individual sources: the cortical-resection RNA-seq never reports its per-group n (only n=3 for Western blot); two histology papers are in moderate-tier MDPI journals with n=3-6 per group and uncorrected t-tests; the 1982 metabolic paper was verified from an abstract reconstruction because the publisher blocked full text, and its tumour clone species is inferred rather than stated. Three tempting reassurances were refuted under adversarial review and must not be used: that U251MG does not migrate through the corpus callosum, the quantified U87 versus patient-derived contralateral cell burdens, and that LITT thermal spread was shown to stay away from the opposite hemisphere. Finally, the step from "the contralateral hemisphere carries a remote myeloid program" to "a contralateral-referenced contrast subtracts part of the tumour-induced signature" is sound reasoning from measured premises but is not itself a measured result: no located paper analysed contralateral tissue as an RNA-seq control in a glioma model, which is the specific gap this study sits in.

## Open questions

- Does U251N in athymic RNU/RNU nude rats cross the corpus callosum, and at what cell burden? No confirmed prior art addresses U251 or rat hosts, and the claim that U251 does not migrate callosally was refuted, so this needs measuring in-house with human-specific detection on serial coronal sections rather than being assumed in either direction.
- Is the ~5% human read fraction with human Y reads in one contralateral sample real graft cells, perfusion/dissection carry-over, or library-level cross-contamination such as index hopping? The other contralateral sample at ~0.6% is the only internal background estimate, no literature gives an expected species cross-mapping floor for xengsort-sorted rat/human brain, and an eight-fold spread across two samples cannot distinguish the mechanisms.
- Does LITT itself perturb the contralateral hemisphere - host inflammation, blood-brain barrier, microglial state - beyond what the untreated tumour already does? The only rat U251 LITT paper reports no contralateral measurement of any kind, so the ablated-and-recurred arm's control tissue carries an unmeasured additional perturbation.
- How large is remote neuroinflammation in a nude rat carrying a human graft, given that the measured SB28/GL261 effect depends on an intact immune compartment? Answering it requires naive or sham rats in the design; without them a contralateral-only reference cannot separate tumour-local from brain-wide change, and with contralateral tissue from 2 of 8 animals the paired design is unavailable for most of the cohort anyway.

## Refuted under adversarial review (do not use)

- The authors say directly that the contralateral hemisphere has often been used as a control and that this is unsafe. They describe the contralateral changes as moderate and centred on inflammation and synaptic remodelling: GSEA found 108 enriched pathways; qRT-PCR confirmed Cacna1c, Slc2a1 and Tgfb1 were up; IL-1b and TGF-b rose while IL-6 did not change (p = 0.37). (https://pmc.ncbi.nlm.nih.gov/articles/PMC12384278/; vote 0-3)
- Astrocytic GFAP rose in the contralateral cortex relative to sham as early as day 3 (p = 0.0443), a smaller rise than on the lesion side (p < 0.0001). By day 28 the contralateral signal nearly matched the ipsilateral: only the group effect (lesion vs sham) stayed significant (p = 0.024), with no hemisphere effect. The contralateral-sham difference therefore grows with time after a unilateral lesion. (https://pmc.ncbi.nlm.nih.gov/articles/PMC12843123/; vote 0-3)
- Conventional serum-cultured human GBM lines behave differently: U251MG (and U373MG) form bulk/demarcated tumours and do NOT migrate through the corpus callosum, so a U251-based orthotopic model is expected to contaminate the contralateral hemisphere far less than patient-derived GSC lines (relevant to whether ~5% human reads contralaterally is true invasion or carry-over). (https://www.nature.com/articles/s41598-023-51063-7; vote 0-3)
- In an orthotopic human GBM xenograft (U87MG in NOD-SCID mice), tumour cells implanted in one hemisphere were already detectable in the contralateral hemisphere 13 days after implantation, so the contralateral hemisphere of a tumour-bearing brain cannot be assumed tumour-cell-free. (https://doi.org/10.3390/cells13020192; vote 0-3)
- Contralateral tumour-cell burden is quantifiable and strongly line-dependent: the expansile U87MG line gave a mean of 3.9 cells per 5 mg of contralateral brain tissue (range 0.8-10.3) whereas the invasive patient-derived MU20 line gave a mean of 3898 cells per 5 mg (range 244-8121), a ~1000-fold difference across cell lines (n=4 mice each). (https://doi.org/10.3390/cells13020192; vote 0-3)
- Graft cells were not confined to the injection site or the contralateral tumour: tumour cells were detected in every sampled tissue of the supratentorial brain, including distal regions, indicating brain-wide (not tumour-local) cell dissemination in both an expansile and an invasive line. (https://doi.org/10.3390/cells13020192; vote 0-3)
- Thermal spread was bounded close to the fibre: with laser at 980 nm, ~1 V for 30-40 s, temperature 4 mm from the fibre tip stayed sublethal (no degrees C reported), arguing against direct thermal injury reaching the opposite hemisphere; the only controls were untreated tumour-bearing rats (n=3), with no naive, sham or contralateral-tissue arm and no follow-up past 24 h. (https://pubmed.ncbi.nlm.nih.gov/34554269/; vote 0-3)

## All sources fetched

- https://pmc.ncbi.nlm.nih.gov/articles/PMC9266805/ (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11474166/ (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC9834327/ (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12384278/ (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12843123/ (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://doi.org/10.1038/jcbfm.1982.3 (primary; Contralateral hemisphere is not naive: transcriptomic and glial changes after a unilateral lesion)
- https://www.nature.com/articles/s41598-023-51063-7 (primary; Orthotopic glioma crossing to the contralateral side, and remote effects of LITT/thermal ablation)
- https://doi.org/10.3390/cells13020192 (primary; Orthotopic glioma crossing to the contralateral side, and remote effects of LITT/thermal ablation)
- https://pubmed.ncbi.nlm.nih.gov/34554269/ (primary; Orthotopic glioma crossing to the contralateral side, and remote effects of LITT/thermal ablation)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC8598735/ (primary; Orthotopic glioma crossing to the contralateral side, and remote effects of LITT/thermal ablation)
- https://doi.org/10.3390/ijms23137308 (primary; Prior art: contralateral tissue as the within-animal control in omics (glioma, MCAO, TBI, epilepsy))
- https://doi.org/10.1038/s41598-023-27663-8 (primary; Prior art: contralateral tissue as the within-animal control in omics (glioma, MCAO, TBI, epilepsy))
- https://doi.org/10.1186/1471-2164-14-282 (primary; Prior art: contralateral tissue as the within-animal control in omics (glioma, MCAO, TBI, epilepsy))
- https://doi.org/10.1158/1078-0432.CCR-24-1563 (primary; Prior art: contralateral tissue as the within-animal control in omics (glioma, MCAO, TBI, epilepsy))
- https://doi.org/10.1186/s40478-021-01232-4 (primary; Prior art: contralateral tissue as the within-animal control in omics (glioma, MCAO, TBI, epilepsy))
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11474166/ (primary; What a same-brain control uniquely shows: tumour-local versus brain-wide and systemic effects)
- https://pubmed.ncbi.nlm.nih.gov/32303613/ (primary; What a same-brain control uniquely shows: tumour-local versus brain-wide and systemic effects)
- https://www.sciencedirect.com/science/article/pii/S1074761326001214 (primary; What a same-brain control uniquely shows: tumour-local versus brain-wide and systemic effects)
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7117758/ (primary; What a same-brain control uniquely shows: tumour-local versus brain-wide and systemic effects)
- https://academic.oup.com/bioinformatics/article/37/2/192/5878955 (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
- https://www.bioconductor.org/packages//release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.html (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
- https://www.biorxiv.org/content/10.1101/2020.05.14.095604v1.full (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12006369/ (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC3667020/ (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2817684/ (primary; Analysis consequences: paired/mixed models, RUVSeq negative controls, xenograft cross-mapping floor)
