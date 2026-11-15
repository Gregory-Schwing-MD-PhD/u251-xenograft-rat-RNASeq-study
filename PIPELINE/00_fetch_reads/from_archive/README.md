# U251: starting the work again from the public archive

Written 2026-11-15 by the line u251-from-geo (part KIT) for the owner's order of that night: the latest U251 work must
be able to go on from GEO (GSE338105) and the public rat sets, not from files held on the lab volume. This folder is
the reproduction kit. Its home on the grid is /rs/rs_grp_oschome/go2432/u251_reproduce/ (copies only: nothing there
was moved, and nothing anywhere was removed by making it).

## What is here

| file | what it is |
|---|---|
| manifest.tsv | every raw read file the work needs, as the archive holds it: study, run, library, layout, file, https address, bytes, md5, read count, GEO sample, first public. 48 files: the ten GSE338105 libraries (PRJNA1492510, 20 files) and the 22 public rat runs of the floor analysis (PRJNA627944 and PRJNA238328, 28 files) |
| build_manifest.py | writes manifest.tsv from ENA's file reports (ena/, read 2026-11-15; --from-ena reads them again) and library_map.csv |
| fetch_from_archive.sh | brings any rows of the manifest back into a folder you name, sized from the archive before any byte lands, all or nothing, every file checked against the archive's md5; on the grid's lab volume it is admitted by the lab grid repository's size rule (grid/size-rule.sh) |
| samplesheet_from_manifest.py | the nf-core sample sheet (sample,fastq_1,fastq_2,strandedness) naming the fetched files, never an old local path |
| library_map.csv | the repository's join of run accessions to library names (PIPELINE/00_fetch_reads/library_map.csv) |
| captured.classes, CAPTURED.sha256.tsv | every code, configuration, sample sheet, parameter file, log and small table copied from the u251_* folders into ../captured/ (same relative paths), with its class and sha256 |
| NEXTFLOW_LAUNCHES.txt, NEXTFLOW_HISTORY.tsv | each Nextflow run's command line, pipeline, revision and status, read from every launch folder's .nextflow.log and .nextflow/history |
| PIPELINE_INFO.list | the pipelines' own records of software versions and parameters (pipeline_info), inside ../captured/ |
| REFERENCE_URLS.txt | every reference address (genome, annotation, index source) named in the captured scripts and configs |
| CONTAINER_IMAGES.tsv | the container images the runs pulled (nf-core images, pullable again by name) |
| R_LIBRARY_VERSIONS.tsv, u251_drug_conda_packages.txt | the R libraries' package versions and the u251_drug conda environment's packages |
| methylarray_3a7ff00_u251patched.diff | the uncommitted patch to nf-core/methylarray at 3a7ff00 that the methylation steps run |
| container/ | the hand-built image's labels, embedded definition, environment, R package list, and bioconductor_pinned.def (a pinned rebuild recipe, not yet built) |
| ../home_clone/ | the grid home clone's uncommitted diff, its modified and untracked files under 5 MB, and the list of larger ones |

## The route, in order

    K=/rs/rs_grp_oschome/go2432/u251_reproduce/kit   # or wherever this folder is
    DEST=/path/of/your/choice

1. Reads. `bash $K/fetch_from_archive.sh $DEST/fastq PRJNA1492510` (the ten libraries, 78.39 GB) and, for the floor
   analysis, `bash $K/fetch_from_archive.sh $DEST/floor 'PRJNA627944|PRJNA238328'` (53.07 GB). On the lab volume run
   it as a Slurm job with SIZE_RULE=<the lab grid repository clone>/grid/size-rule.sh and DATA=/rs/rs_grp_oschome/go2432/data.
   The repository's own route, nf-core/fetchngs 1.12.0 from GSE338105 (PIPELINE/00_fetch_reads/run_fetchngs.sbatch),
   gives the same reads through sratools and is checked by check_library_map.py.
2. Sample sheet. `python3 $K/samplesheet_from_manifest.py $DEST/fastq --out $DEST/samplesheet.csv`: the same shape as
   the repository's ANALYSIS/samplesheet.csv, with the library names, pointing at the fetched files.
3. References. Built from public sources only; the addresses are in REFERENCE_URLS.txt and in the captured params
   (GRCh38 primary assembly with GENCODE, Ensembl mRatBN7.2 for rat, VEP 115, dbSNP 138, REDIportal). The
   repository's PIPELINE/01_references, which would rebuild them in one step, is not written yet
   (PIPELINE/RESTART_FROM_GEO.md, precondition 3).
4. Species split. The repository's main.nf (submit_xengsort.sh, xengsort_params.yaml, image docker://go2432/xengsort,
   recipe DOCKER/Dockerfile.xengsort) with `--input $DEST/samplesheet.csv`.
5. Per-arm reads (graft and human clipped, host): u251_multiomics/bin/prep_reads_run.sbatch and the repository's
   ANALYSIS/host/extract_host.sbatch, both under ../captured/ or the repository, run on the step-4 output.
6. Pipelines. Each launch in NEXTFLOW_LAUNCHES.txt is `nextflow run nf-core/<pipeline> -r <revision> -profile
   singularity -c <captured config> -params-file <captured params>`, with `--input` set to a sheet made by step 2 (or
   by step 5's output for the clipped reads) in place of the old local path. The pipelines themselves are not
   rewritten.

## What this kit does not prove

See ../../RESTART_FROM_GEO.md and records/PROOF.txt: which derived files a kept script is shown to rebuild, which are not yet shown,
and what could not be captured.
