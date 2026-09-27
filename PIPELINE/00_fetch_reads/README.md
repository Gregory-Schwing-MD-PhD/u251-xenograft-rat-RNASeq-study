# Step 00 — from the public accession to the reads

**What changed and why.** The repository's previous route to the reads was a private Google Drive folder reached with
a personal OAuth token; the accession `GSE338105` appeared nowhere in the README. Nobody outside this laboratory
could start the pipeline. This step replaces that with a published tool driven by the accession alone.

```bash
sbatch -D ~/logs PIPELINE/00_fetch_reads/run_fetchngs.sbatch
```

| what | how |
|---|---|
| resolve the accession | **nf-core/fetchngs 1.12.0** (pinned revision, not a branch) takes `ids.csv`, which contains one line: `GSE338105` |
| download | `--download_method sratools`, each file checked against the md5 the archive publishes |
| hand off | `--nf_core_pipeline rnaseq` writes a samplesheet already in the format step 02 consumes |
| name the samples | `check_library_map.py` joins the run accessions to this project's library names through `library_map.csv`, and exits 1 on any disagreement |

## `library_map.csv` — the one place the two vocabularies meet

The analyses speak of `IL67B`, `NL70B`, `N269B`, `C2B`; the archive speaks of `SRR…`. Everything downstream is
keyed on the library name, so a rename, a re-ordered download or a swapped pair would move a sample between arms
silently. This file is the join, and `run_fetchngs.sbatch` checks it on every run.

The run accessions were read from this project's own SRA verification run (`u251_science/spread/floor`), which
downloaded each run from SRA and asserted, pair for pair, that the archive's content matched the FASTQ the analyses
had used — so the mapping is not an assumption, it is a checked identity. The GSM and SRX columns are `PENDING`
until they are read back from GEO and ENA and cross-checked.

Two names to keep straight, because both appear in the raw records:

- the libraries `NL70B` and `NL71B` are the arrays `IL70B` and `IL71B`;
- the library `N269B` is the array `N2`.

## What this step does not do

It does not fetch the methylation arrays: those are not in GEO yet (`GEO/methylation/` prepares that submission).
It does not build the references — that is step 01, and it must be run before step 02 can quantify anything.
