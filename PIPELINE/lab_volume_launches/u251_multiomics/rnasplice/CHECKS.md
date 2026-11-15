# rnasplice checks

Written 2026-09-27 05:49:24 by check_rnasplice.py. Results: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results`. Nominal splicing counts are not reported here: they mean nothing until ranked against the 10-split null (`null/results/NULL_REPORT.md`).
**hard_pass: False**; failed: 5


## 0. run, versions, applied params: PASS

- NEXTFLOW_EXIT: `0` (want 0)
- pipeline commit: `1d0494ae3402d1a46e0adadad24f81a0ff855c77` (want `1d0494ae3402d1a46e0adadad24f81a0ff855c77`, GitHub refs/tags/1.0.4)
- Nextflow version in launch/.nextflow.log: `23.10.0` (want 23.10.0)
- params actually applied (`params_2026-09-27_05-48-57.json`): all 26 as intended

## 1. trace: expected processes, nothing forbidden: PASS

- TRIMGALORE: 6 task(s) for 6 of 6 samples 
- FASTQC: 6 task(s) for 6 of 6 samples 
- STAR_ALIGN: 6 task(s) for 6 of 6 samples 
- SAMTOOLS_SORT: 6 task(s) for 6 of 6 samples 
- DEXSEQ_COUNT: 6 task(s) for 6 of 6 samples 
- SALMON_QUANT_STAR: 6 task(s) for 6 of 6 samples 
- STAR_GENOMEGENERATE: 1 
- DEXSEQ_ANNOTATION: 1 
- DEXSEQ_EXON: 1 
- CREATE_BAMLIST: 1 
- RMATS_PREP: 1 
- RMATS_POST: 1 
- DRIMSEQ_FILTER: 1 
- DEXSEQ_DTU: 1 
- STAGER: 1 
- GENERATE_EVENTS_IOE: 1 
- GENERATE_EVENTS_IOI: 1 
- PSIPEREVENT: 1 
- PSIPERISOFORM: 1 
- DIFFSPLICE_IOE: 1 
- DIFFSPLICE_IOI: 1 
- forbidden processes (edgeR, featureCounts, MISO, bigWig): none
- FAILED task attempts across all traces (retried or not): 1: NFCORE_RNASPLICE:RNASPLICE:ALIGN_STAR:BAM_SORT_STATS_SAMTOOLS:SAMTOOLS_SORT (IL69B_graft)
- salmon pseudo-alignment branch (SALMON_INDEX, SALMON_QUANT_SALMON, *_SALMON DTU/SUPPA): absent, so `pseudo_aligner: null` was honoured and there is one result set

## 2. commands: strand, read length, unpaired statistics, DRIMSeq filter: PASS

- RMATS_PREP (recurrent-primary): ok
- RMATS_POST (recurrent-primary): ok
- DEXSEQ_COUNT (IL66B_graft): ok
- DEXSEQ_COUNT (NL71B_graft): ok
- DEXSEQ_COUNT (NL70B_graft): ok
- DEXSEQ_COUNT (IL68B_graft): ok
- DEXSEQ_COUNT (IL67B_graft): ok
- DEXSEQ_COUNT (IL69B_graft): ok
- DIFFSPLICE_IOE (recurrent-primary): ok
- DIFFSPLICE_IOI (recurrent-primary): ok
- SALMON_QUANT_STAR (IL66B_graft): ok
- SALMON_QUANT_STAR (NL70B_graft): ok
- SALMON_QUANT_STAR (IL68B_graft): ok
- SALMON_QUANT_STAR (NL71B_graft): ok
- SALMON_QUANT_STAR (IL67B_graft): ok
- SALMON_QUANT_STAR (IL69B_graft): ok
- TRIMGALORE (IL66B_graft): ok
- TRIMGALORE (NL71B_graft): ok
- TRIMGALORE (IL68B_graft): ok
- TRIMGALORE (IL69B_graft): ok
- TRIMGALORE (NL70B_graft): ok
- TRIMGALORE (IL67B_graft): ok
- DRIMSEQ_FILTER (1): ok
- STAR_ALIGN (IL66B_graft): ok
- STAR_ALIGN (NL71B_graft): ok
- STAR_ALIGN (NL70B_graft): ok
- STAR_ALIGN (IL68B_graft): ok
- STAR_ALIGN (IL67B_graft): ok
- STAR_ALIGN (IL69B_graft): ok

## 3. inputs: xengsort graft pairs -> TrimGalore -> STAR: PASS

- IL67B_graft: TrimGalore R1 31886295, R2 31886295 (xengsort graft 31886295); pairs kept 31883963; STAR input 31883963; unique 84.47% 
- IL68B_graft: TrimGalore R1 32225887, R2 32225887 (xengsort graft 32225887); pairs kept 32223637; STAR input 32223637; unique 83.27% 
- IL69B_graft: TrimGalore R1 49327614, R2 49327614 (xengsort graft 49327614); pairs kept 49324067; STAR input 49324067; unique 85.10% 
- IL66B_graft: TrimGalore R1 25459283, R2 25459283 (xengsort graft 25459283); pairs kept 25457488; STAR input 25457488; unique 82.87% 
- NL70B_graft: TrimGalore R1 31293679, R2 31293679 (xengsort graft 31293679); pairs kept 31291440; STAR input 31291440; unique 85.79% 
- NL71B_graft: TrimGalore R1 29210358, R2 29210358 (xengsort graft 29210358); pairs kept 29208313; STAR input 29208313; unique 87.96% 

## 4. read lengths: R1 max 101, R2 max 98 (clip applied once, upstream): PASS

- IL67B_graft: R1 35-101 nt, R2 32-98 nt, 31886295/31886295 reads 
- IL68B_graft: R1 35-101 nt, R2 32-98 nt, 32225887/32225887 reads 
- IL69B_graft: R1 35-101 nt, R2 32-98 nt, 49327614/49327614 reads 
- IL66B_graft: R1 35-101 nt, R2 32-98 nt, 25459283/25459283 reads 
- NL70B_graft: R1 35-101 nt, R2 32-98 nt, 31293679/31293679 reads 
- NL71B_graft: R1 35-101 nt, R2 32-98 nt, 29210358/29210358 reads 

## 5. strand: salmon ISR and R1 orientation on named genes: FAIL

- salmon IL67B_graft: lib_format_counts.json missing
- salmon IL68B_graft: lib_format_counts.json missing
- salmon IL69B_graft: lib_format_counts.json missing
- salmon IL66B_graft: lib_format_counts.json missing
- salmon NL70B_graft: lib_format_counts.json missing
- salmon NL71B_graft: lib_format_counts.json missing
- independent count: 20 genes found in the GTF (11 plus, 9 minus): ACTB, B2M, CD44, EEF1A1, EEF2, ENO1, FN1, FTL, GAPDH, HSP90AB1, LDHA, MYL6, PKM, RPLP0, RPS3, SPARC, TMSB10, TPT1, TUBB, VIM
- IL67B_graft: R1 antisense to the gene 0.9881 (plus genes R1 fwd 5542 / rev 484939; minus genes R1 fwd 937978 / rev 11644)
- IL68B_graft: R1 antisense to the gene 0.9878 (plus genes R1 fwd 6065 / rev 540502; minus genes R1 fwd 1106057 / rev 14277)
- IL69B_graft: R1 antisense to the gene 0.9869 (plus genes R1 fwd 9348 / rev 801474; minus genes R1 fwd 1620553 / rev 22893)
- IL66B_graft: R1 antisense to the gene 0.9890 (plus genes R1 fwd 3621 / rev 353932; minus genes R1 fwd 696463 / rev 8025)
- NL70B_graft: R1 antisense to the gene 0.9880 (plus genes R1 fwd 5210 / rev 477449; minus genes R1 fwd 923761 / rev 11798)
- NL71B_graft: R1 antisense to the gene 0.9882 (plus genes R1 fwd 5356 / rev 437465; minus genes R1 fwd 873618 / rev 10343)
- worst of all strand measures: 0.0000 -> FAIL (PASS >= 0.97, the published RSeQC 97.4-98.4 %; FAIL < 0.90)

## 6. completeness: six samples everywhere: PASS

- DEXSeq counts: IL67B_graft 699611 rows, IL68B_graft 699611 rows, IL69B_graft 699611 rows, IL66B_graft 699611 rows, NL70B_graft 699611 rows, NL71B_graft 699611 rows 
- rMATS bamlists: recurrent ['IL66B_graft_sorted.bam', 'NL71B_graft_sorted.bam', 'NL70B_graft_sorted.bam']; primary ['IL68B_graft_sorted.bam', 'IL67B_graft_sorted.bam', 'IL69B_graft_sorted.bam'] 
- rMATS prep: 6 .rmats file(s) listing ['IL66B_graft_sorted.bam', 'IL67B_graft_sorted.bam', 'IL68B_graft_sorted.bam', 'IL69B_graft_sorted.bam', 'NL70B_graft_sorted.bam', 'NL71B_graft_sorted.bam'] 
- rMATS post outputs: all 10 MATS files
- SUPPA TPM header: ['IL66B_graft', 'IL67B_graft', 'IL68B_graft', 'IL69B_graft', 'NL70B_graft', 'NL71B_graft'] 

## 7. null inputs: files, containers, pipeline scripts, option strings: PASS

- DEU_COUNTS: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_exon/counts` 
- DEU_GFF: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_exon/annotation/DEXSeq.gff` 
- DEU_PIPE_GENE: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_exon/results/perGeneQValue.recurrent_vs_primary.csv` 
- DMDS: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_dtu/filter/drimseq/dmDSdata.rds` 
- DTU_PIPE_GENE: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_dtu/results/dexseq/perGeneQValue.recurrent-primary.tsv` 
- STAGER_PIPE: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/dexseq_dtu/results/stager/getAdjustedPValues.recurrent-primary.tsv` 
- TPM: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/tximport/suppa_tpm.txt` 
- IOE: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/generate_events/per_local_event/events.ioe` 
- IOI: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/generate_events/per_isoform/events.ioi` 
- PSI_L: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/psi_per_local_event/suppa_local.psi` 
- PSI_I: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/psi_per_isoform/suppa_isoform.psi` 
- SUPPA_PIPE_L: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/diffsplice/per_local_event/recurrent-primary_local_diffsplice.dpsi` 
- SUPPA_PIPE_I: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/suppa/diffsplice/per_isoform/recurrent-primary_transcript_diffsplice.dpsi` 
- RMATS_TEMP: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/rmats/recurrent-primary/rmats_temp` 
- RMATS_PIPE_POST: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/rmats/recurrent-primary/rmats_post` 
- BAMLIST_T: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/rmats/bamlist/recurrent_bamlist.txt` 
- BAMLIST_C: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon/rmats/bamlist/primary_bamlist.txt` 
- BAM_DIR: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/results/star_salmon` 
- GTF: `/rs/rs_grp_oschome/go2432/u251_multiomics/refs/gencode_v44/GRCh38.primary_assembly.annotation.gtf` 
- PB: `/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/pipeline_bin` 
- IMG_DEXSEQ: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-bioconductor-dexseq-1.36.0--r40_0.img` 
- IMG_DTU: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-bioconductor-dexseq-1.36.0--r40_0.img` 
- IMG_DRIMSEQ: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-bioconductor-drimseq-1.18.0--r40_0.img` 
- IMG_STAGER: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-bioconductor-stager-1.12.0--r40_0.img` 
- IMG_RBASE: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-r-base-3.4.2.img` 
- IMG_SUPPA: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-suppa-2.3--py36_0.img` 
- IMG_RMATS: `/rs/rs_grp_oschome/go2432/u251_multiomics/singularity/depot.galaxyproject.org-singularity-mulled-v2-8ea76ff0a6a4c7e5c818fd4281abf918f92eeeae-121e48ab4817ec619c157a346458efca1ccf3c0a-0.img` 
- pipeline_bin/run_dexseq_exon.R md5 c5134b464df2577fee97ce693b7e0b83 = 1.0.4 tarball
- pipeline_bin/run_dexseq_dtu.R md5 5220ddf92dc9d46177c99e4b92bea573 = 1.0.4 tarball
- pipeline_bin/run_stager.R md5 34ef168316bb8f7f26d66d85cdf3b3a4 = 1.0.4 tarball
- pipeline_bin/suppa_split_file.R md5 c82be75ec05459de1f586902305ce65b = 1.0.4 tarball
- pipeline_bin/run_drimseq_filter.R md5 c14830c40df37c5597eb92a68d2de708 = 1.0.4 tarball
- SUPPA diffSplice options (from DIFFSPLICE_IOE): `-m empirical -gc -s -c -a 1000 -l 0 -al 0.05 -th 0 -nan 0`
- rMATS options (from RMATS_POST, minus inputs/outputs/threads/task): `-t paired --libType fr-firststrand --readLength 100 --variable-read-length --cstat 0.0001 --allow-clipping` 

