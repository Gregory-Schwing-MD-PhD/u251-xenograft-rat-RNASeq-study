# Leave-one-tumour-out GSEA (U251 LITT recurrence), GSEA 4.3.2, the nf-core command with one tumour removed

**Reproduction of the published run (all six, seed 1234):** exact = True; sets 8869 vs 8869, unmatched 0; max |diff| NES 0, nominal p 0, FDR q 0.

## Translation initiation (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) by condition

| condition | genes | NES (seed 1234) | nominal p | FDR q (seed 1234) | FDR q, 5 seeds (min-max) | FWER p | rank by q |
|---|---|---|---|---|---|---|---|
| full | 19351 | -1.99 | 0.000 | 0.022 | 0.022-0.312 (n = 5) | 0.022 | 1 |
| drop_IL67B | 19061 | -1.94 | 0.000 | 0.082 | 0.080-0.590 (n = 5) | 0.236 | 3 |
| drop_IL68B | 18948 | -1.93 | 0.000 | 0.022 | 0.001-0.032 (n = 5) | 0.106 | 6 |
| drop_IL69B | 18845 | -1.95 | 0.000 | 0.597 | 0.106-1.000 (n = 5) | 0.857 | 2041 |
| drop_IL66B | 19176 | -1.95 | 0.000 | 0.360 | 0.129-1.000 (n = 5) | 0.756 | 2 |
| drop_NL70B | 18918 | -1.99 | 0.000 | 0.005 | 0.001-0.015 (n = 5) | 0.013 | 1 |
| drop_NL71B | 18993 | -1.92 | 0.000 | 0.335 | 0.007-0.335 (n = 5) | 0.975 | 9 |
| drop_IL67B_nofilter | 19351 | -1.94 | 0.000 | 0.126 | 0.126-0.126 (n = 1) | 0.123 | 1 |
| drop_IL68B_nofilter | 19351 | -1.97 | 0.000 | 0.002 | 0.002-0.002 (n = 1) | 0.004 | 2 |
| drop_IL69B_nofilter | 19351 | -1.99 | 0.000 | 0.052 | 0.052-0.052 (n = 1) | 0.057 | 1 |
| drop_IL66B_nofilter | 19351 | -2.00 | 0.000 | 0.029 | 0.029-0.029 (n = 1) | 0.029 | 1 |
| drop_NL70B_nofilter | 19351 | -2.01 | 0.000 | 0.001 | 0.001-0.001 (n = 1) | 0.001 | 1 |
| drop_NL71B_nofilter | 19351 | -1.99 | 0.000 | 0.031 | 0.031-0.031 (n = 1) | 0.033 | 1 |

## The six leading sets, FDR q at seed 1234 (NES in brackets)

| set | full | drop_IL67B | drop_IL68B | drop_IL69B | drop_IL66B | drop_NL70B | drop_NL71B |
|---|---|---|---|---|---|---|---|
| KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION | 0.022 (-1.99) | 0.082 (-1.94) | 0.022 (-1.93) | 0.597 (-1.95) | 0.360 (-1.95) | 0.005 (-1.99) | 0.335 (-1.92) |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | 0.158 (-1.96) | 0.017 (-1.97) | 0.010 (-1.96) | 0.563 (-1.94) | 1.000 (-1.86) | 0.059 (-1.94) | 0.070 (-1.98) |
| REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY | 0.146 (-1.95) | 0.225 (-1.90) | 0.040 (-1.91) | 0.863 (-1.92) | 0.822 (-1.91) | 0.031 (-1.96) | 0.059 (-1.97) |
| KEGG_RIBOSOME | 0.204 (-1.94) | 0.304 (-1.89) | 0.011 (-1.95) | 0.086 (-1.98) | 0.533 (-1.93) | 0.033 (-1.95) | 0.081 (-1.97) |
| REACTOME_SELENOAMINO_ACID_METABOLISM | 0.169 (-1.94) | 0.177 (-1.91) | 0.035 (-1.92) | 0.023 (-2.01) | 0.709 (-1.96) | 0.009 (-1.98) | 0.146 (-1.95) |
| REACTOME_CELLULAR_RESPONSE_TO_STARVATION | 0.247 (-1.93) | 0.189 (-1.92) | 0.068 (-1.89) | 1.000 (-1.90) | 1.000 (-1.89) | 0.052 (-1.94) | 0.081 (-1.99) |

## Shape of the screen (seed 1234; ranges over the five seeds in loo_screen.tsv)

| condition | down q<0.05 | down q<0.25 | up q<0.05 | up q<0.25 | best down q (set) | best up q (set) | down p<0.05 | up p<0.05 |
|---|---|---|---|---|---|---|---|---|
| full | 1 | 6 | 0 | 0 | 0.022 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.608 (PEDERSEN_METASTASIS_BY_ERBB2_ISOFORM_1) | 303 | 406 |
| drop_IL67B | 2 | 9 | 0 | 0 | 0.017 (REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION) | 0.601 (GOBP_REGULATION_OF_OSTEOBLAST_PROLIFERATION) | 294 | 385 |
| drop_IL68B | 12 | 29 | 0 | 0 | 0.003 (DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_FIBROBLASTS_STROMAL_CELLS) | 0.765 (GOBP_CALCIUM_ION_TRANSMEMBRANE_IMPORT_INTO_CYTOSOL) | 474 | 199 |
| drop_IL69B | 1 | 2 | 0 | 0 | 0.023 (REACTOME_SELENOAMINO_ACID_METABOLISM) | 0.473 (GOBP_NEGATIVE_REGULATION_OF_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION) | 188 | 576 |
| drop_IL66B | 0 | 0 | 0 | 0 | 0.339 (WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS) | 0.575 (GOBP_SINGLE_STRANDED_VIRAL_RNA_REPLICATION_VIA_DOUBLE_STRANDED_DNA_INTERMEDIATE) | 230 | 426 |
| drop_NL70B | 6 | 13 | 0 | 0 | 0.005 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.708 (REACTOME_VIRAL_MESSENGER_RNA_SYNTHESIS) | 328 | 229 |
| drop_NL71B | 1 | 8 | 0 | 0 | 0.045 (WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS) | 0.595 (SA_CASPASE_CASCADE) | 276 | 397 |
| drop_IL67B_nofilter | 0 | 4 | 0 | 0 | 0.126 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.601 (REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION) | 305 | 382 |
| drop_IL68B_nofilter | 15 | 28 | 0 | 0 | 0.002 (DURANTE_ADULT_OLFACTORY_NEUROEPITHELIUM_FIBROBLASTS_STROMAL_CELLS) | 0.762 (GOBP_SPLICEOSOMAL_SNRNP_ASSEMBLY) | 486 | 201 |
| drop_IL69B_nofilter | 0 | 1 | 0 | 1 | 0.052 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.175 (MITSIADES_RESPONSE_TO_APLIDIN_DN) | 191 | 564 |
| drop_IL66B_nofilter | 1 | 1 | 0 | 0 | 0.029 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.477 (MITSIADES_RESPONSE_TO_APLIDIN_DN) | 233 | 418 |
| drop_NL70B_nofilter | 4 | 13 | 0 | 0 | 0.001 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.717 (GOBP_POSITIVE_REGULATION_OF_CELLULAR_CATABOLIC_PROCESS) | 325 | 232 |
| drop_NL71B_nofilter | 1 | 6 | 0 | 0 | 0.031 (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION) | 0.607 (ZHAN_MULTIPLE_MYELOMA_CD1_DN) | 288 | 405 |
