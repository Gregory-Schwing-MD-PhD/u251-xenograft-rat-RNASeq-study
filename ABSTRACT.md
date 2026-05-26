# Characterizing the Transcriptomic Recurrence Signature of Glioblastoma Following Laser Interstitial Thermal Therapy

**Introduction.** Laser interstitial thermal therapy (LITT) is a minimally invasive treatment for glioblastoma (GBM) cytoreduction via thermal coagulation of the tumor mass. It leaves a sublethal thermal penumbra where resistant tumor cells persist, yet the molecular basis of this adaptation is poorly understood. We hypothesized that post-LITT recurrence carries a distinct transcriptomic signature relative to the pre-LITT primary tumor and sought druggable targets.

**Methods.** Orthotopic U251N GBM xenografts in RNU/RNU rats underwent MRI-guided LITT. Bulk RNA-seq from primary (pre-LITT, n=3) and recurrent (post-LITT, n=3) tumors was analyzed after xengsort cross-species read sorting (k=25). Differential expression used DESeq2 (Wald test, ashr shrinkage, BH/independent filtering; FDR<0.05, |log2FC|>1) via nf-core/differentialabundance, with global structure by PCA/PERMANOVA on variance-stabilized counts, hub genes from a STRING protein–protein interaction network, subtypes by GSVA with limma, and drugs prioritized by an Integrated Score combining therapeutic enrichment (|NES|^1.5) with predicted blood–brain-barrier (BBB) permeability.

**Results.** Primary and recurrent tumors separated by PCA (PC1 41.8%, PC2 26.2%); PERMANOVA attributed 27.8% of variance to recurrence (R²=0.278), with gene- and pathway-level significance. DESeq2 identified 35 genes (23 up, 12 down). GSEA showed dominant downregulation of the translational machinery in recurrence—eukaryotic translation elongation (NES=−3.17, FDR=1.6×10⁻²⁶), ribosome, translation initiation, and the EIF2AK4/GCN2 amino-acid-starvation response—indicating translational shutdown with integrated stress—unlikely a contamination artifact: residual rat-to-human mis-mapping, highest in the rat-richer recurrent tumors, would raise these conserved genes opposite to their observed decrease. Subtypes shifted toward Garofano Mitochondrial, away from Neftel Astrocyte-like/NPC-like. Hubs were ECM/calcium genes (BGN, COL1A1, IGFBP3, CALB1). Drug prioritization nominated the iron-chelator/HIF agent ciclopirox (Integrated Score 3.40, BBB 1.0) and HIF-stabilizer dimethyloxalylglycine (3.04); the PI3K/mTOR inhibitor LY-294002 ranked 7th.

**Conclusions.** Translational shutdown, characterized by an integrated stress/starvation response and a mitochondrial-metabolic subtype shift, seems to be one post-LITT recurrence mechanism, consistent with stress-adapted, persistent cells. HIF/iron-axis agents (ciclopirox, DMOG) and the PI3K/mTOR inhibitor LY-294002 emerge as potential therapeutic candidates.

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

<sub>Abstract body ≈ 299 words. Figure: `publication_figure/Publication_Figure_9Panel_VOLCANO_COMPLETE.png` (PDF also in that directory).</sub>
