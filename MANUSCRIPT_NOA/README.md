# The Neuro-Oncology Advances manuscript, compiled from the pipeline

This directory holds the submission and, more importantly, the chain that makes every number in it auditable. It is
the same mechanism the RSNA AI manuscript uses (`spinesurg-ct-nnunet/SchwingNet/tools/results_macros.py`), carried
over on 2026-09-27 at Greg's instruction: *the tables and results in the paper are placeholders for pending results,
they are stored as JSONs which are imported at compile time, and each JSON maps back to a specific run of a
reproducible pipeline, for full programmatic updating of the paper and reproducibility auditing from GEO download to
manuscript submission.*

## The chain

```
GEO / SRA accession            GSE338105 (RNA), the EPIC arrays (accession pending)
        |
        v   nf-core + published tools, one grid job per step
result file on the grid        e.g. u251_science/rna/composition/composition_main.json   (job 40481434)
        |
        v   tools/collect_results.py      copies it here with its sha256 and the job id
results/<file>.json            + results/RUNS.json   file -> job -> what it is
        |
        v   numbers_manifest.json         key -> file + RFC 6901 pointer + format + job
        |
        v   tools/results_macros.py       at every compile
generated/numbers.tex          \noa@val@KEY, \noa@src@KEY  (or \noa@pend@KEY with the job)
        |
        v   \RESULT{key} in the text
main.pdf                       and generated/numbers_audit.md, the table of every value and its source
```

A key whose file or pointer does not resolve does not silently vanish: it prints in the PDF as a red
**XXX Result Pending-&lt;job&gt;**, so an unfinished result is visible on the page, and `build.sh --check` exits 1.

## Commands

```bash
bash build.sh              # collect, compose the figures, resolve the numbers, compile
bash build.sh --check      # the same, then fail on anything pending, unused, stale or hand-typed
python tools/results_macros.py --check    # just the numbers audit
```

`--check` fails the build when any of these is true:

| check | what it catches |
|---|---|
| a key is pending | a result the paper quotes has not been produced yet |
| a manifest key is unused | a number was removed from the text but left in the manifest |
| a bare result-looking number sits in the Abstract or Results | someone typed a value instead of using `\RESULT` |
| a result file no longer matches the copy it was taken from | the upstream job was re-run and the paper is stale |
| a composed figure is older than one of its panels | a panel was redrawn and the display item was not rebuilt |

## What is here

| path | what it is |
|---|---|
| `main.tex`, `sections/` | the manuscript, in the journal's required order |
| `noaresults.sty` | `\RESULT`, `\RESULTSRC`, `\RESULTFIG`, `\PENDING`; `[review]` prints each number's source beside it |
| `numbers_manifest.json` | every result the paper quotes: key -> file, pointer, format, job, note |
| `results/` | byte copies of the grid outputs, plus `RUNS.json` mapping each to its job |
| `figures/` | the six display items, composed from single panels at 600 dpi, each with a `.json` naming its panels and their hashes |
| `generated/` | written at build time: `numbers.tex`, `numbers_audit.md`, `pending.json` |
| `tools/` | `collect_results.py`, `compose_figures.py`, `results_macros.py` |

## The journal's limits, and where each is enforced

Neuro-Oncology Advances, Basic and Translational Investigation: title <= 160 characters; running title <= 50;
structured abstract <= 250 words; 5 keywords; 2-3 key points totalling <= 260 characters; "Importance of the Study"
<= 150 words; body (Introduction through Discussion) <= 6000 words; **6 display items**; <= 50 references; five
required statements; a 100-word lay summary. `tools/limits.py` \PENDING{to be written}{this pass} will count them at
build time; until then they are counted by hand in the editor's note at the head of `sections/`.

## What is deliberately not here

No IDAT, FASTQ, BAM or lab file. The results files are summary statistics only. The manuscript drafts under `MBR/`
are gitignored; this directory is the submission and is tracked.
