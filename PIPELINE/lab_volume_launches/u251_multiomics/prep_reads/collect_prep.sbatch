#!/bin/bash
# collect_prep: runs afterany on the prep_reads array. Gathers the ten check rows into reads/CHECKS.tsv, writes
# reads/READY only if all ten passed and every output file is present, and reports the index pair of each library
# with the index-hopping verdict for N269B (reads/INDEX_HOPPING.txt). It does not edit RESUME.md.
# Submitted by submit.sh (as $M/bin/collect_prep.sbatch, --dependency=afterany:<array>).
#SBATCH --job-name=u251_prep_collect
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/reads/logs/%x_%j.out
set -eu
export LC_ALL=C
M=${PREP_M:-/rs/rs_grp_oschome/go2432/u251_multiomics}
CONF=$M/prep_reads
R=$M/reads
LIBS=(C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B NL71B)
pv() { awk -v k="$1" '{ sub(/[ \t]*#.*$/, "") } $1 == k ":" { print $2; f = 1; exit } END { if (!f) exit 1 }' "$CONF/params.yaml"; }
RICH=$(pv human_rich_min) || RICH=0.20

echo "[$(date '+%F %T')] collect_prep on $(hostname); prep_reads array ${PREP_ARRAY:-unknown}"
if [ -n "${PREP_ARRAY:-}" ] && command -v sacct > /dev/null; then
  sacct -j "$PREP_ARRAY" --format=JobID%-16,State%-12,ExitCode,Elapsed,MaxRSS,NodeList%-12 2>/dev/null || true
fi

# ---------------------------------------------------------------- CHECKS.tsv and READY
rm -f "$R/READY"
OUT=$R/CHECKS.tsv
: > "$OUT.part"
rows=0; ok=0; problems=""
for S in "${LIBS[@]}"; do
  f=$R/checks/$S.tsv
  if [ -s "$f" ]; then
    if [ ! -s "$OUT.part" ]; then head -n 1 "$f" >> "$OUT.part"; fi
    tail -n +2 "$f" >> "$OUT.part"
    rows=$((rows + 1))
    st=$(awk -F'\t' 'NR == 2 { print $2 }' "$f")
    if [ "$st" = OK ]; then ok=$((ok + 1)); else problems="$problems $S:status=$st"; fi
  else
    problems="$problems $S:no_checks_row"
  fi
  for o in "graft_clip/${S}_graft_R1.fq.gz" "graft_clip/${S}_graft_R2.fq.gz" "human_clip/${S}_human_R2.fq.gz" \
           "both_names/$S.names.gz" "index_seq/$S.txt"; do
    if [ ! -s "$R/$o" ]; then problems="$problems $o:missing"; fi
  done
done
mv -f "$OUT.part" "$OUT"
echo "CHECKS.tsv: $rows rows, $ok OK"
if [ "$rows" -eq 10 ] && [ "$ok" -eq 10 ] && [ -z "$problems" ]; then
  { echo "prep_reads READY $(date '+%F %T')"
    echo "array ${PREP_ARRAY:-unknown}; CHECKS.tsv md5 $(md5sum < "$OUT" | cut -c1-32)"
    echo "graft pairs: $R/graft_clip/<lib>_graft_R{1,2}.fq.gz (R2 clipped 3 nt at 5'; do not clip again)"
    echo "all human R2, clipped: $R/human_clip/<lib>_human_R2.fq.gz"; } > "$R/READY.part"
  mv -f "$R/READY.part" "$R/READY"
  echo "READY written: $R/READY"
else
  echo "NOT READY:$problems"
fi

# ---------------------------------------------------------------- index pairs and the N269B hopping question
TAB=$R/index_pairs.tsv
printf 'lib\tgraft_share\tmodal_index\ti7\ti5\tmodal_share_note\n' > "$TAB.part"
for S in "${LIBS[@]}"; do
  e=$(pv "expected_$S" || echo "0,0,0")
  share=$(awk -F, '{ s = $1 + $2 + $3; printf "%.4f", (s > 0 ? $2 / s : 0) }' <<< "$e")
  idx=$R/index_seq/$S.txt
  mi=NA; note=NA
  if [ -s "$idx" ]; then
    mi=$(awk -F'\t' '$1 == "modal_index" { print $2 }' "$idx"); [ -n "$mi" ] || mi=NA
    note=$(awk -F'\t' '$1 == "modal_index" { print $3 }' "$idx"); [ -n "$note" ] || note=NA
  fi
  case "$mi" in
    *+*) i7=${mi%%+*}; i5=${mi#*+} ;;
    *)   i7=$mi; i5=NA ;;
  esac
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$S" "$share" "$mi" "$i7" "$i5" "$note" >> "$TAB.part"
done
mv -f "$TAB.part" "$TAB"
echo "index pairs (modal index field of every R1 header in the human file):"
column -t -s $'\t' "$TAB" 2>/dev/null || cat "$TAB"

awk -F'\t' -v rich="$RICH" -v now="$(date '+%F %T')" '
function ham(a, b,   i, d) { if (length(a) != length(b)) return -1
                             d = 0; for (i = 1; i <= length(a); i++) if (substr(a, i, 1) != substr(b, i, 1)) d++
                             return d }
function okseq(s) { return s ~ /^[ACGTN]+$/ }
NR == 1 { next }
{ n++; lib[n] = $1; sh[n] = $2 + 0; i7[n] = $4; i5[n] = $5; if ($1 == "N269B") t = n }
END {
  print "N269B index-hopping check (" now ")"
  print "human-rich = graft/(host+graft+both) >= " rich
  if (!t) { print "VERDICT: no N269B row; not assessed."; exit }
  if (!okseq(i7[t])) {
    print "N269B modal index field: " i7[t] (i5[t] != "NA" ? "+" i5[t] : "")
    print "VERDICT: the index sequences are not in these read headers (after xengsort), so index hopping into N269B cannot be assessed from these files. The sequencing core sample sheet would settle it."
    exit }
  printf "N269B: i7 %s, i5 %s, graft share %.4f\n", i7[t], i5[t], sh[t]
  shared_rich = ""; shared_any = ""
  for (r = 1; r <= n; r++) {
    if (r == t) continue
    s7 = (okseq(i7[r]) && i7[r] == i7[t]); s5 = (okseq(i5[t]) && okseq(i5[r]) && i5[r] == i5[t])
    printf "  %-6s graft share %.4f%s  i7 %s (Hamming %d)  i5 %s (Hamming %d)%s\n", lib[r], sh[r], (sh[r] >= rich ? " human-rich" : "           "),
           i7[r], ham(i7[r], i7[t]), i5[r], ham(i5[r], i5[t]), (s7 || s5 ? "  SHARED" : "")
    if (s7) { shared_any = shared_any " " lib[r] "(i7)"; if (sh[r] >= rich) shared_rich = shared_rich " " lib[r] "(i7)" }
    if (s5) { shared_any = shared_any " " lib[r] "(i5)"; if (sh[r] >= rich) shared_rich = shared_rich " " lib[r] "(i5)" }
  }
  if (shared_rich != "")
    print "VERDICT: N269B shares an index with human-rich libraries:" shared_rich ". Index-specific hopping into N269B cannot be excluded; its human reads are not shown to be its own by indexing."
  else if (!okseq(i5[t]))
    print "VERDICT: only one index is in the headers; no human-rich library shares N269B i7, but without i5 the dual-index argument cannot be made."
  else if (shared_any != "")
    print "VERDICT: N269B shares an index only with libraries that are not human-rich:" shared_any ". Hopping from a human-rich library needs both indexes to switch."
  else
    print "VERDICT: N269B i7 and i5 are each unique among the ten libraries. A read from another library is assigned to N269B only if both of its indexes switch. This does not exclude carry-over or cross-contamination before indexing."
}' "$TAB" > "$R/INDEX_HOPPING.txt.part"
mv -f "$R/INDEX_HOPPING.txt.part" "$R/INDEX_HOPPING.txt"
cat "$R/INDEX_HOPPING.txt"
echo "[$(date '+%F %T')] collect_prep done"
