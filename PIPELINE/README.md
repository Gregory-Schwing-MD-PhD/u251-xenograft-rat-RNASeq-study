# The pipeline, from the public accession to the manuscript

Every step here is a published tool driven by a public accession or a committed parameter file. The intent, in
Greg's words (2026-09-27), is that the methods can be *cited* rather than defended: the pipeline should be
predominantly other people's tools, with as little bespoke bioinformatics as the questions allow.

```
  GEO GSE338105 (RNA)                        EPIC arrays (GEO submission in GEO/methylation/)
        |                                              |
   00_fetch_reads      nf-core/fetchngs 1.12.0         |
        |              + library_map.csv               |
        v                                              |
   01_references       GENCODE 44 human, Ensembl 110 rat, xengsort index, gene-set bundle
        |                                              |
        v                                              v
   02_quantify         xengsort 2.1.0 (species split) -> nf-core/rnaseq 3.22.2 (STAR + salmon)
        |              minfi / SeSAMe / conumee / nf-core/methylarray on the array side
        v
   03_analysis         nf-core/differentialabundance 1.5.0, GSEA, GSVA, the pre-registered
        |              decomposition and its nulls (ANALYSIS/science, ANALYSIS/spread, ANALYSIS/power)
        v
   MANUSCRIPT_NOA      each result file -> results/ with its sha256 and job id -> numbers_manifest.json
                       -> \RESULT{} macros -> main.pdf, with a red "Result Pending" for anything unfinished
```

| step | what it is | status |
|---|---|---|
| `00_fetch_reads` | nf-core/fetchngs on `GSE338105`, then the run-accession-to-library check | **committed, runnable** |
| `01_references` | the reference bundle: genomes, annotations, the xengsort index, the gene-set GMTs | **being rebuilt** — the original build script was lost |
| `02_quantify` | the species split and the two quantification runs | exists as launch directories under `ANALYSIS/`; being folded in here |
| `03_analysis` | the pre-registered analyses | `ANALYSIS/science`, `spread`, `xval`, `power`, each with its own pre-registration and hash |
| `MANUSCRIPT_NOA` | the paper, compiled from the result files | **committed, builds** |

## Why this directory exists

An audit on 2026-09-27 found two breaks between the public data and the published numbers:

1. **The reads were not reachable.** The README's route was a private Google Drive folder behind a personal token,
   and the accession appeared nowhere. Fixed in `00_fetch_reads`.
2. **The references were not rebuildable.** `ANALYSIS/create_refs_final.slurm` built the human and rat references and
   the 8,869-set `combined_human.gmt` that is the universe of every GSEA number in the paper — and that script no
   longer exists in the repository or its history. `01_references` reconstructs the recipe from the surviving files
   and records every source URL and checksum; whatever cannot be rebuilt from public sources is archived instead,
   and said so plainly.

Both were breaks in reproducibility, not in the results: the analyses themselves ran from pinned containers and
recorded their inputs' hashes. But a stranger could not have started, which is the standard this directory is held
to.
