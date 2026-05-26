# Characterizing the Transcriptomic Recurrence Signature of Glioblastoma Following Laser Interstitial Thermal Therapy

**Introduction.** Glioblastoma (GBM) recurrence is inevitable after cytoreductive therapy. Laser Interstitial Thermal Therapy (LITT) leaves a sublethal thermal penumbra where resistant cells persist, yet the molecular basis of this adaptation is poorly understood. We hypothesized that post-LITT recurrence carries a distinct transcriptomic signature relative to the pre-LITT primary tumor, and sought druggable targets.

**Methods.** Orthotopic U251 GBM xenografts in RNU rats underwent MRI-guided LITT. Bulk RNA-seq from Primary (pre-LITT, n=3) and Recurrent (post-LITT, n=3) tumors was analyzed after xengsort cross-species read sorting; residual host contamination was tested and found negligible. Global separation was quantified by PERMANOVA (permutational multivariate ANOVA) on variance-stabilized counts. Differentially expressed genes (DESeq2; FDR<0.05 and |log2FC|>1) were complemented by GSEA, a STRING protein–protein interaction (PPI) network of those genes, GSVA subtype scoring, and drug prioritization by an Integrated Score combining therapeutic enrichment (|NES|^1.5) with predicted blood–brain-barrier (BBB) permeability.

**Results.** PERMANOVA attributed 27.8% of transcriptome-wide variance to recurrence status (R²=0.278), confirming a substantial global shift, with significance established at the gene and pathway level. DESeq2 identified 35 differentially expressed genes (23 up, 12 down). GSEA showed dominant downregulation of the translational machinery in recurrence—eukaryotic translation elongation (NES=−3.17, FDR=1.6×10⁻²⁶), ribosome, translation initiation, and the EIF2AK4/GCN2 amino-acid-starvation response—indicating translational shutdown with integrated stress, robust to contamination (which would bias conserved genes oppositely). Subtypes shifted toward Garofano Mitochondrial and away from Neftel Astrocyte-like/NPC-like states. Highest-degree hubs of the PPI network were extracellular-matrix/calcium genes (BGN, COL1A1, IGFBP3, CALB1). Drug prioritization nominated the iron-chelator/HIF agent ciclopirox (Integrated Score 3.40; BBB 1.0; clinical phase 4) and HIF-stabilizer dimethyloxalylglycine (3.04); the PI3K/mTOR inhibitor LY-294002 ranked 7th.

**Conclusions.** Post-LITT recurrence is defined by translational shutdown, an integrated stress/starvation response, and a mitochondrial-metabolic subtype shift—consistent with stress-adapted persister cells. Brain-penetrant HIF/iron-axis agents (ciclopirox, DMOG, deferoxamine) emerge as rational therapeutic candidates.

---

![U251 post-LITT recurrence: 9-panel transcriptomic and pharmacogenomic figure](publication_figure/Publication_Figure_9Panel_VOLCANO_COMPLETE.png)

**Figure 1.** U251 glioblastoma post-LITT recurrence — transcriptomic signature and therapeutic prioritization.
**(A)** Experimental design.
**(B)** Global structure: PCA of Primary vs Recurrent tumors and the differential-expression volcano.
**(C)** GBM molecular-subtype trajectories (Neftel, Garofano, Verhaak) from primary to recurrence; asterisks denote significant shifts (array-weighted limma).
**(D)** Enriched-pathway clustering (GSEA).
**(E)** STRING protein–protein interaction network of differentially expressed genes, hub genes highlighted.
**(F)** Polypharmacology network linking candidate drugs to enriched pathways.
**(G)** Drug–pathway gene-overlap heatmap.
**(H)** Integrated drug scoring (|NES| vs blood–brain-barrier penetration; point size = Integrated Score = |NES|^1.5 × BBB).
**(I)** Top-ranked drug candidates.

<sub>Abstract body ≈ 289 words. Figure: `publication_figure/Publication_Figure_9Panel_VOLCANO_COMPLETE.png` (PDF also in that directory).</sub>
