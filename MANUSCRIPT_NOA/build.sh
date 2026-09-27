#!/bin/bash
# Build the Neuro-Oncology Advances manuscript from the pipeline outputs, end to end.
#
#   bash build.sh            # collect -> resolve -> compose -> compile
#   bash build.sh --check    # the same, then refuse (exit 1) on any pending number, unused key, bare number,
#                            # stale figure, or result file that no longer matches the copy it was taken from
#
# The chain this enforces: a grid job writes a result file; tools/collect_results.py copies it into results/ with
# its sha256 and the job id (results/RUNS.json); numbers_manifest.json points a key at a value inside it;
# tools/results_macros.py resolves the key at compile time into generated/numbers.tex; the text says \RESULT{key}.
# Nothing in the manuscript is typed by hand, and every number in the PDF can be walked back to a job.
set -u
cd "$(dirname "$0")"
CHECK=${1:-}
RC=0

echo "== 1. collect the result files the manuscript compiles from"
python tools/collect_results.py ${CHECK:+--check} || RC=1

echo "== 2. compose the display items from their panels"
python tools/compose_figures.py ${CHECK:+--check} || RC=1

echo "== 3. resolve every number in the manifest"
python tools/results_macros.py ${CHECK:+--check} || RC=1

echo "== 4. compile"
TECTONIC="../MBR/.build/tectonic.exe"          # the portable engine MBR/build.sh fetches; MBR/ is gitignored
if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex || RC=1
  latexmk -c >/dev/null 2>&1
elif [ -x "$TECTONIC" ]; then
  echo "   using the portable engine $TECTONIC"
  "$TECTONIC" -X compile --keep-logs --outdir . main.tex || RC=1
else
  echo "   no LaTeX engine found: skipping the PDF (the numbers and figures above are still current)."
  echo "   Get one with:  bash ../MBR/build.sh   (fetches tectonic into MBR/.build/), or compile on the grid."
fi

echo "== done (rc=$RC)"
echo "   generated/numbers_audit.md  every value, its source file, that file's sha256, and the job"
echo "   results/RUNS.json           every result file -> the grid job that wrote it"
exit $RC
