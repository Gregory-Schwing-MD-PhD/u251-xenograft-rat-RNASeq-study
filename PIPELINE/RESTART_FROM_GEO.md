# Wiping the study off the home directory and restarting from the accession

Greg's intent, 2026-09-27: *"wipe this entirely off go2432's home and just restart from GEO download on the ceph
drive — maybe not now but eventually."* Nothing has been deleted. This file records what that costs, what has to be
true before it is safe, and what to do in what order.

## Where the study actually lives today

| location | size | what it is |
| --- | --- | --- |
| `/rs/rs_grp_oschome/go2432/u251_multiomics` | 1.1 TB | Ceph |
| `/rs/rs_grp_oschome/go2432/u251_science` | 529 GB | Ceph |
| `/rs/rs_grp_oschome/go2432/u251_host` | 405 GB | Ceph |
| `/rs/rs_grp_oschome/go2432/u251_meth` | 55 GB | Ceph — includes the laboratory IDAT archive and the R library and sandbox every R step uses |
| `/rs/rs_grp_oschome/go2432/u251_geo_meth` | 179 MB | Ceph — the staged GEO methylation deposit |
| `~/u251-xenograft-murine-RNASeq-study` | **204 GB** | **home** — the git clone; 201 GB of it is `ANALYSIS/` |

So 2.1 TB is already on Ceph and only 204 GB sits in home. Ceph is at 6 of 20 TB (34.7 %); home is at
**4160.69 GB against a 4294.97 GB soft quota (96.9 %)**, hard limit 8589.93 GB.

**Removing U251 from home recovers 204 GB of a 4.16 TB problem** — 96.9 % to about 92.2 %. It is not the fix for the
quota. The four directories that are: `data` 769 GB, `spine-detector` 468 GB, `spinesurg-ct-nnunet` 367 GB,
`FluoroCap` 316 GB. Do the U251 move when the study is *ready* to be reproducible, not because of disk pressure.

Separately, and unrelated to any quota: the **login node's `/tmp` is 100 % full** (25 GB, 28 KB free), filled by other
users' build trees. It breaks shell heredocs and R temporary files. Point `TMPDIR` at `/dev/shm` or a Ceph path for
login-node work; there is nothing of ours to delete there.

## What a wipe would destroy

`git status` in the home clone shows **47 untracked entries**, several of them whole result trees. These exist in no
other place — not in git, not on Ceph, not in GEO:

- `ANALYSIS/gsea_leave_one_out/runs/`, `inputs/`, `smoke/` — the 41 GSEA runs, about 8.7 CPU-hours. The *summaries*
  (`loo_screen.tsv`, `loo_sets.tsv`, `loo_top_down.tsv`, `SUMMARY.md`, `reproduction.json`) **are** committed, so the
  evidence survives a wipe and only the raw run trees would have to be regenerated.
- `ANALYSIS/holdout_IL66B/` and `ANALYSIS/holdout_IL68B/` — `control/`, `holdout/`, `published/`,
  `results_therapy_v3_no*/`, `subtypes_src/`.
- `ANALYSIS/cibersort/results/*` and `signatures_online/` — note the signature files came from the CIBERSORTx website,
  which is a registration, not a download the pipeline can repeat unattended.
- `ANALYSIS/gsc_drugs/results/`, `results_published/`, `refs_online/`.
- `ANALYSIS/holdout_separation/loo_separation_host.tsv`.

Each has to be triaged into one of three buckets before the wipe: **commit** (small, and it is evidence a manuscript
cites), **regenerate** (a chained pipeline step reproduces it), or **move to Ceph** (large and expensive to
regenerate). That triage is the real work and it has not been done.

## Three preconditions, one of them currently false

1. **The methylation arrays must be in GEO.** They are not yet. The deposit is staged and verified at
   `/rs/rs_grp_oschome/go2432/u251_geo_meth` (866,238 probes x 8 arrays, minfi 1.56.0, md5s written) and the metadata
   and cover note are ready, but until Greg uploads it the only copy of the raw IDATs is the laboratory archive at
   `u251_meth/raw/DNAmethylation-20260926T220706Z-1-001.zip` — which is on Ceph, and whose sha256 the staging job
   checks, so it is not at risk from a *home* wipe. This precondition is about the restart, not the wipe.
2. **`salmon.merged.gene_lengths.tsv` must be obtainable from the deposit — and it is not.** GSE338105 carries
   `GSE338105_salmon.merged.gene_counts.tsv.gz` only. The differential-expression step needs the gene-length matrix for
   tximport to give DESeq2 its per-gene effective-length offsets, so **a stranger who downloads GSE338105 today cannot
   reach the DE result**; they must re-quantify every library. Either GEO adds the file (it is item 5 of
   `GEO/methylation/COVER_NOTE.md`) or step 02 must be the documented source of it. **This is the one precondition that
   is false right now, and it is the thing that makes "restart from GEO download" not yet true.**
3. **`PIPELINE/01_references` must exist.** Step 00 (`00_fetch_reads`) fetches the reads from the accession with
   nf-core/fetchngs, and its `library_map.csv` is verified against SRA and ENA. The reference-bundle rebuild is still
   unwritten, so the chain from accession to counts has a hole in it.

## Order of operations, when the time comes

1. Upload the methylation deposit and send the cover note, which also asks GEO to add the gene-length matrix.
2. Write `PIPELINE/01_references`.
3. Triage the 47 untracked entries into commit / regenerate / move-to-Ceph, and record the decision per entry.
4. Re-clone the repository **on Ceph**, not in home, and run the chain from step 00 to the manuscript on the Ceph copy.
5. Only once that run reproduces the committed summaries — `reproduction.json` exact, the numbers manifest resolving
   80/80 — remove the home copy.

Step 5 is the check that makes the wipe safe: the restart has to be demonstrated *before* the original is deleted,
never after. Nothing about this is urgent, because U251 is 5 % of the quota problem.
