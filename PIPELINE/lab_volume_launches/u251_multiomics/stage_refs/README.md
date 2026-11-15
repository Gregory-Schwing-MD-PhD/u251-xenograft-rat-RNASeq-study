# stage_refs: stage and verify the references the released pipelines read

One Slurm job, no Nextflow. It downloads or unpacks every reference that nf-core/rnavar 1.3.0, nf-core/rnafusion 4.1.3
and the integration steps take as input, and it checks each one against a source outside this build. It exists so that
no pipeline quietly rebuilds a reference or reads the wrong one. Every other item gates on the `.done` markers written
here.

## Question

Are the references on disk the files their publishers released, on the same GRCh38 coordinates as each other, and in
the layout each pipeline expects? And will nf-core/rnafusion use the staged bundle, or will it decide the bundle is
missing and quietly build or download its own? For STAR-Fusion, a rebuild takes about 2 days at 150 GB.

## How to run

```
# grid, from WSL:  scp this folder to /rs/rs_grp_oschome/go2432/u251_multiomics/stage_refs/  (plus bin/crc32c, see below)
sbatch -D /rs/rs_grp_oschome/go2432/u251_multiomics/refs /rs/rs_grp_oschome/go2432/u251_multiomics/stage_refs/run.sbatch
```

- **Resources.** reqp `--qos=requeue`, 4 CPU, 8 GB, 12 h.
- **Log.** `refs/logs/u251_stage_refs_<jobid>.out`.
- **Requeue.** The job is requeue-safe. Each step writes `refs/.done/<step>` only when all of its checks pass. After a
  pre-emption, finished steps are skipped and partial downloads resume (`curl -C -` into `*.part`).
- **Failures.** A failed step does not stop the steps after it. Only `chain` and `known_sites` need `gencode`. The job
  exits 2 if any required step failed.
- **Where results go.**
  - Check rows go to `refs/.checks/<step>.tsv`, gathered into `refs/REFS_CHECKS.tsv` with the columns step, check,
    expected, observed, result, outside_reference.
  - `refs/STATUS.tsv` holds one line per step.
  - Both files are copied to `stage_refs/results/`.

`bin/crc32c` is a static x86-64 build of `crc32c.c` (gcc 13, `-O2 -static`). It is tested against the CRC-32C check
value (`123456789` gives `e3069283`) and against Google's stored hash for the Mills `.tbi`. It is not committed. If
it fails on the node, the job compiles `crc32c.c` with the node's gcc and runs the same test.

## Inputs and outputs (M = /rs/rs_grp_oschome/go2432/u251_multiomics)

| step | source | output under `M/refs/` |
|---|---|---|
| containers | `depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0` (pinned tag, the one nf-core/modules uses); samtools 1.21 image already in the cache | `/rs/rs_grp_oschome/go2432/singularity_cache/depot.galaxyproject.org-singularity-bcftools-1.21--h8b25389_0.img` |
| gencode | `U/ANALYSIS/refs/human/` GENCODE v44 (Ensembl 110) GRCh38 primary assembly `.fa.gz` + `.gtf.gz` (the files the published human run used) | `gencode_v44/GRCh38.primary_assembly.genome.fa{,.fai}`, `.annotation.gtf` |
| chain | UCSC `goldenPath/hg38/liftOver/hg38ToRn7.over.chain.gz` | `chains/hg38ToRn7.over.chain.gz` |
| known_sites | Broad `gcp-public-data--broad-references/hg38/v0/`: `Homo_sapiens_assembly38.dbsnp138.vcf` (10.95 GB, uncompressed), `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` + `.tbi` | `known_sites/dbsnp138.primary.vcf.gz`, `mills_1000g.primary.vcf.gz`, each with `.tbi` |
| vep | Ensembl `release-115/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz` (25.27 GB) + `CHECKSUMS` | `vep_cache/homo_sapiens/115_GRCh38/` (tarball kept in `vep_cache/download/`) |
| rnafusion | `s3://nf-core-awsmegatests/rnafusion/references/`: `GRCh38/gencode_v46/{starfusion/ctat_genome_lib_build_dir,star,fusioncatcher,gencode}`, `GRCh38/arriba/2.5.0`, `GRCh38/hgnc`, `GRCh38/fusion_report_db`. 322 objects, 122.1 GB; not `salmon`, not `gffread` | `rnafusion/GRCh38/...` (the key paths below `references/`), `rnafusion/MANIFEST.tsv`, `rnafusion/REMOVED.tsv` |
| rediportal (optional) | `rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz` (1.51 GB; v3 is the current database version on the download page) | `rediportal/TABLE1_hg38_v3.txt.gz`, or `rediportal/NOT_AVAILABLE` with the reason |

`params.yaml` lists these paths under the parameter names the pipelines use. `samplesheet.csv` lists each reference,
its grid path, the `.done` marker that vouches for it, and the item that reads it.

## Checks, and what each is checked against

Each check is compared with something produced outside this build.

- **gencode.**
  - The md5 of each `.gz` matches the copy of GENCODE's `MD5SUMS` kept with the download. It also matches the
    `MD5SUMS` fetched fresh from ftp.ebi.ac.uk, which confirms the kept copy is GENCODE's.
  - `gzip -d` exits 0, which verifies each member's CRC-32 and length.
  - chr1-22, X, Y and M in the fai have the GRCh38 lengths from the assembly report. Those lengths were checked
    against UCSC `hg38.chrom.sizes` when this job was written.
  - The fai has 194 sequences. This row is INFO only.
  - The GTF header says `version 44 (Ensembl 110)`.
- **chain.**
  - The md5 matches UCSC's `liftOver/md5sum.txt`, fetched at run time.
  - `gzip -t` passes.
  - The hg38 sizes in the chain headers equal the GENCODE v44 fai for chr1-22, X and Y.
  - The rn7 sizes equal the chromosome lengths of the host run's own rat genome, mRatBN7.2 (Ensembl 110), for
    chr1-20, X and Y. The lengths come from `H/results_rnaseq/genome/index/star/chrNameLength.txt`. So rn7 is the
    assembly the host reads were aligned to. A local test on 2026-09-26 found all 24 and all 22 equal.
- **known_sites.**
  - Each downloaded object's CRC-32C matches Google's stored `x-goog-hash`, and so does its MD5 where Google stores
    one. The dbSNP file is a composite object, so it has a CRC-32C only.
  - The Broad header's chr1-22, X, Y and M lengths equal the GENCODE v44 fai. These are two independent providers of
    GRCh38, so a build mismatch shows up here.
  - The cut file keeps every chr1-22, X, Y and M record of the source, per `bcftools index -s`.
  - Record counts per contig are identical before and after the reheader.
  - The final `##contig` lines are exactly the GENCODE fai, in fai order. This is the dictionary GATK compares with
    the reference.
  - `tabix` succeeds on the final file, which proves it is BGZF and sorted.
  - REF alleles agree with the GENCODE genome sequence itself. An evenly spaced sample of up to 20,000 records with a
    plain-ACGT REF (dbSNP: chr21, chr22, chrY, chrM; Mills: every primary chromosome) is compared with
    `samtools faidx` of the same span, and at least 99 % must agree. A build or coordinate mismatch makes most of them
    disagree.
- **vep.**
  - Size is checked against Ensembl's Content-Length.
  - BSD `sum` (1 KiB blocks) must equal Ensembl's `CHECKSUMS` line, which is how Ensembl computes it.
  - `tar` exits 0, and every non-directory member in tar's own listing is on disk.
  - `homo_sapiens/115_GRCh38/info.txt` exists and names GRCh38. That path is rnavar's layout,
    `${vep_cache}/${vep_species}/${vep_cache_version}_${vep_genome}`.
  - The directories 1-22, X, Y and MT exist.
- **rnafusion.**
  - Every listed object is present at its listed size, and its content matches the ETag S3 computed at upload: MD5
    for single-part objects, MD5 of the part MD5s for multipart ones. All common part sizes that give the listed part
    count are tried.
  - The live listing is compared with the listing captured when the pipeline was evaluated
    (`expected_s3_listing.txt`). A change is reported, not blocked.
  - The v46 fai equals the v44 fai from EBI on chr1-22, X, Y and M (FAIL-able); on all 194 sequences this is INFO.
  - The S3 STAR index, the CTAT library's `ref_genome.fa.fai` and the CTAT STAR index were each built on that same
    v46 FASTA.
  - The v46 GTF header says `version 46 (Ensembl 112)`.
  - `star/genomeParameters.txt` has `sjdbOverhang 99`.
- **rnafusion pre-flight.** This emulates `exists_not_empty()`, `subworkflows/local/build_references/main.nf`
  l.238-262. It is described below.
- **rediportal.**
  - Size is checked against the server's Content-Length, and `gzip -t` passes.
  - At least 99 % of rows are A>G or T>C, which is A-to-I editing on either strand.
  - The `Ref` column agrees with the GENCODE v44 base in at least 99 % of 5,000 evenly spaced sites. Positions are
    1-based: chr1:87158 T and chr1:103668 T were confirmed against the UCSC API.

## Settings that are not the plan's defaults, and why

1. **Which file `exists_not_empty()` really tests.** The plan says the function checks the first file it finds. The
   Groovy code, however, is `eachFileRecurse(FILES) { file -> first_file = file; return }`. A `return` inside a
   closure does not stop `eachFileRecurse`, so `first_file` ends up as the **last** file visited in readdir order,
   descending into directories where they appear.
   - The pre-flight therefore records both the last-visited file (the one the pipeline tests) and the first, and
     requires both to be larger than 0 bytes.
   - The walk uses `os.scandir`, which gives the same readdir order as Java's `File.listFiles`.
   - The risk is real. `ctat_genome_lib_build_dir` holds 24 zero-byte files (`__chkpts/*.ok`,
     `ref_genome.fa.star.idx/build.ok`).
   - CephFS readdir order is not alphabetical (`ls -f` on the singularity cache shows hash order). A test on NTFS put
     `cosmic.db` first in `fusion_report_db`.
2. **`fusion_report.log` is set aside, not only `cosmic.db` deleted.** `fusion_report_db` on S3 has two 0-byte files.
   - `cosmic.db` is deleted, as the plan says. `no_cosmic: true` makes FUSIONREPORT_DETECT run with `--no-cosmic`
     (`conf/modules.config` l.131).
   - `fusion_report.log` is a log that FUSIONREPORT_DOWNLOAD moves into the directory, and nothing reads it. It is
     moved to `rnafusion/_set_aside/`, so the directory check no longer depends on readdir order.
   - Both are recorded in `rnafusion/REMOVED.tsv`.
3. **The pre-flight is repeated from a second node.** The node that wrote the files may list a directory in the order
   its cache holds, not the order a fresh client sees. So when `rnafusion` passes, the job submits
   `preflight_recheck.sbatch` with `--exclude=<this node>`. If that node disagrees, it deletes `.done/rnafusion` and
   writes `rnafusion/PREFLIGHT_FAILED_<node>`.
   - **The rnafusion launcher should run the same command on its own head node right before `nextflow run`.** It
     exits 1 when the launch should not go ahead:
     `python3 M/stage_refs/refs_util.py preflight --base M/refs/rnafusion --ver 46`
4. **The REF check samples with `samtools faidx`, not `bcftools norm --check-ref`.** Counting what `norm -c x` drops
   would also count records that norm cannot process for other reasons. A direct comparison of REF with the genome
   sequence has only one way to disagree.
5. **The reheader fallback.** The job uses `bcftools reheader --fai`, as the plan says. If the result still carries
   `##contig` lines for Broad-only contigs (HLA, decoys), it rebuilds the header from the fai with `reheader -h`. The
   final dictionary must be exactly the GENCODE fai either way. `REFS_CHECKS.tsv` records which method was used.
6. **Intermediates are deleted after `known_sites` passes.** That removes the raw 10.95 GB dbSNP VCF, the full bgzipped
   copies and the cut files. The index stats, headers and REF samples stay in `known_sites/work/`. The raw Mills files
   (20 MB) are kept.
7. **The bcftools image is saved under Nextflow's own name for the depot URL**, in the existing CephFS cache. A
   pipeline that uses `bcftools:1.21--h8b25389_0` later finds it there and does not pull it again.
8. **Log path.** The log goes to `M/refs/logs/`, per the item spec, not `M/logs/`. The results are also copied to
   `M/stage_refs/results/`.
9. **REDIportal comes over http.** The https HEAD request for the v3 table did not answer within 30 s. The http URL
   redirects to https and answered with a Content-Length of 1,510,475,996.

## What the result can show, and what it cannot

**It can show:**
- The bytes on disk are the bytes GENCODE, UCSC, the Broad, Ensembl, nf-core's S3 bucket and REDIportal published.
- The human files all sit on one GRCh38 coordinate system, checked by sequence lengths, by REF alleles against the
  sequence, and by the v44/v46 comparison.
- The chain's rat side is the host run's rat genome.
- rnafusion 4.1.3, as written, will take the staged directories, as far as readdir order on the checking nodes shows.

**It cannot show:**
- **It does not know whether the S3 rnafusion bundle is scientifically right.** It is nf-core's own build: CTAT
  library, STAR 2.7.x index, FusionCatcher database. This job checks that it arrived intact and was built on the
  GENCODE v46 FASTA, and nothing more.
- **It does not know whether rnavar's VEP (release 115) and this cache agree beyond the version number.** rnavar
  checks that when it runs.
- **It does not know the order in which rnafusion's head node will list the CTAT directory.** The two-node check
  lowers that risk but does not remove it, hence the launcher gate in point 3 above.
- **Known-sites content** (which variants dbSNP 138 or Mills contains) is taken as published. This job checks only
  that the files agree with the genome.
- **REDIportal** is optional. If it cannot be fetched, the rnavar integration uses the strand rule and says so.

## Files here

- `run.sbatch`: the job.
- `s3_mirror.py`: anonymous S3 ListObjectsV2 plus `curl`, and ETag verification with a requeue-safe cache.
- `refs_util.py`: the checks, using the standard library only.
- `preflight_recheck.sbatch`: the second-node pre-flight.
- `crc32c.c`: the CRC-32C tool.
- `expected_s3_listing.txt`: the S3 listing captured on 2026-09-26.
- `params.yaml`, `samplesheet.csv`, `grid.config`: the reference contract and the container settings.

**Local dry run (2026-09-26).** The dry run was not committed; it lives in the session scratchpad. It ran in WSL with
bcftools/samtools/htslib 1.21 from bioconda and used real downloads for GENCODE v44, the UCSC chain, Broad Mills, the
bcftools image and the small rnafusion prefixes. The VEP tarball, dbSNP and REDIportal were small stand-ins.
