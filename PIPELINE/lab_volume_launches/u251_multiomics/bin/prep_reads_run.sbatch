#!/bin/bash
# prep_reads: the one canonical U251 read set every RNA arm uses. README.md in this folder says why.
#   graft_clip/<lib>_graft_R1.fq.gz   the xengsort `graft` records of the human R1 file, byte for byte
#   graft_clip/<lib>_graft_R2.fq.gz   the same records of R2, with the Takara Pico v2 3-nt 5' artefact removed
#   human_clip/<lib>_human_R2.fq.gz   every human-file R2 record (graft + both), clipped the same way
#   both_names/<lib>.names.gz         read names of the xengsort `both` records (rat-leak audits)
#   index_seq/<lib>.txt               header comment of the first and last R1 record, modal index field, lane tally
#   checks/<lib>.tsv                  one row, written only when every check passes
# The clip is applied here, once. Nothing downstream may clip R2 again (rnaseq: extra_trimgalore_args null).
# Submitted by submit.sh (as $M/bin/prep_reads_run.sbatch, sbatch -D $M/reads).
# Exit: 0 all checks passed; 2 a check failed (outputs renamed *.failed); 3 gate A failed (nothing written);
#       4 setup error (bad params, unreadable input, broken pigz).
#SBATCH --job-name=u251_prep
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --array=1-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/reads/logs/%x_%A_%a.out
set -eu   # no pipefail: `decompress | head` ends on SIGPIPE by design; every output is checked below
export LC_ALL=C

M=${PREP_M:-/rs/rs_grp_oschome/go2432/u251_multiomics}   # PREP_M only for a local dry run
CONF=$M/prep_reads
R=$M/reads
T=${SLURM_ARRAY_TASK_ID:?not an array task}
LIBS=(C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B NL71B)
[ "$T" -ge 1 ] && [ "$T" -le 10 ] || { echo "array task $T outside 1-10"; exit 4; }
S=${LIBS[$((T - 1))]}

GR1=$R/graft_clip/${S}_graft_R1.fq.gz
GR2=$R/graft_clip/${S}_graft_R2.fq.gz
HC2=$R/human_clip/${S}_human_R2.fq.gz
NAMES=$R/both_names/$S.names.gz
IDX=$R/index_seq/$S.txt
CHK=$R/checks/$S.tsv
W=$M/tmp/prep_reads/$S          # never node-local /tmp (swept)

ts()   { date '+%F %T'; }
log()  { echo "[$(ts)] $S $*"; }
die()  { echo "[$(ts)] SETUP ERROR $S: $*"; exit 4; }
quarantine() {   # a failed library must not leave a plausible file under its final name
  local f
  for f in "$GR1" "$GR2" "$HC2" "$NAMES" "$IDX"; do
    if [ -e "$f" ]; then mv -f "$f" "$f.failed"; fi
  done
}
fail() { echo "[$(ts)] CHECK FAILED $S: $*"; quarantine; exit 2; }
pv()   { awk -v k="$1" '{ sub(/[ \t]*#.*$/, "") } $1 == k ":" { print $2; f = 1; exit } END { if (!f) exit 1 }' "$CONF/params.yaml"; }
rd()   { cat "$1" 2>/dev/null || true; }

# ---------------------------------------------------------------- setup
[ -r "$CONF/params.yaml" ] && [ -r "$CONF/samplesheet.csv" ] || die "params.yaml or samplesheet.csv missing in $CONF"
ROW=$(awk -F, -v t="$T" 'NR == t + 1' "$CONF/samplesheet.csv")
IFS=, read -r S_SHEET HR1 HR2 RR1 XLOG <<< "$ROW"
[ "$S_SHEET" = "$S" ] || die "samplesheet row $T is '$S_SHEET', expected $S"

PIGZ=$(pv pigz)            || die "pigz not in params.yaml"
PT=$(pv pigz_threads)      || die "pigz_threads not in params.yaml"
CLIPN=$(pv r2_clip_nt)     || die "r2_clip_nt not in params.yaml"
MAXL=$(pv max_clipped_len) || die "max_clipped_len not in params.yaml"
GATE_N=$(pv gate_records)  || die "gate_records not in params.yaml"
G_C12=$(pv gate_c12_max)   || die "gate_c12_max not in params.yaml"
G_A12=$(pv gate_a12_max)   || die "gate_a12_max not in params.yaml"
G_C59=$(pv gate_c59_min)   || die "gate_c59_min not in params.yaml"
P_C12=$(pv post_c12_min)   || die "post_c12_min not in params.yaml"
EXP=$(pv "expected_$S")    || die "expected_$S not in params.yaml"
[ "$CLIPN" = 3 ] || die "r2_clip_nt is $CLIPN but the clip below is hard-coded to 3 nt (substr(\$0,4))"

# A library whose checks passed in an earlier run is not redone (delete its checks row to force it).
if [ -s "$CHK" ] && [ "$(awk -F'\t' 'NR == 2 { print $2 }' "$CHK")" = OK ] \
   && [ -s "$GR1" ] && [ -s "$GR2" ] && [ -s "$HC2" ] && [ -s "$NAMES" ] && [ -s "$IDX" ]; then
  log "already done: $CHK ($(date -r "$CHK" '+%F %T')); delete it to redo"
  echo "CHECK OK $S (previous run)"
  exit 0
fi

log "task $T on $(hostname), job ${SLURM_ARRAY_JOB_ID:-?}_$T, restart count ${SLURM_RESTART_COUNT:-0}"
log "awk: $(awk --version 2>&1 | head -n 1)"
log "gzip: $(gzip --version 2>&1 | head -n 1); pigz: $("$PIGZ" --version 2>&1 | head -n 1)"
[ -x "$PIGZ" ] || die "pigz not executable: $PIGZ"
[ "$(echo roundtrip | "$PIGZ" -c | "$PIGZ" -dc)" = roundtrip ] || die "pigz round trip failed on $(hostname)"
for f in "$HR1" "$HR2" "$RR1" "$XLOG"; do [ -r "$f" ] || die "cannot read $f"; done
log "inputs: $HR1 ($(stat -c %s "$HR1") B), $HR2 ($(stat -c %s "$HR2") B), $RR1 ($(stat -c %s "$RR1") B)"

mkdir -p "$R/graft_clip" "$R/human_clip" "$R/both_names" "$R/index_seq" "$R/checks" "$R/logs"
rm -f "$CHK" "$GR1.failed" "$GR2.failed" "$HC2.failed" "$NAMES.failed" "$IDX.failed"
rm -rf "$W"; mkdir -p "$W"; export TMPDIR=$W

# xengsort counts, parsed exactly as ANALYSIS/host/extract_host.sbatch does, and asserted against params.yaml
read -r HOST GRAFT BOTH < <(awk -v s="$S" '$1 == s && NF == 6 { print $2, $3, $5; exit }' "$XLOG") \
  || fail "no '$S' row with 6 fields in $XLOG"
[ "$HOST,$GRAFT,$BOTH" = "$EXP" ] || fail "xengsort host,graft,both = $HOST,$GRAFT,$BOTH; params.yaml expects $EXP"
N=$((GRAFT + BOTH))
log "host=$HOST graft=$GRAFT both=$BOTH; the human files should hold $N records each"

# Base composition of R2 positions 1-10. stdout: records, mean C pos 1-2, mean A pos 1-2, mean C pos 5-9.
COMP='function fr(i, b) { return (t[i] > 0) ? c[i, b] / t[i] : 0 }
NR % 4 == 2 { n++; L = length($0)
              for (i = 1; i <= 10 && i <= L; i++) { t[i]++; c[i, substr($0, i, 1)]++ } }
END { print "pos\treads\tA\tC\tG\tT\tN" > tab
      for (i = 1; i <= 10; i++)
        printf "%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", i, t[i], fr(i, "A"), fr(i, "C"), fr(i, "G"), fr(i, "T"), fr(i, "N") > tab
      close(tab)
      c59 = 0; for (i = 5; i <= 9; i++) c59 += fr(i, "C")
      printf "%.0f %.4f %.4f %.4f\n", n, (fr(1, "C") + fr(2, "C")) / 2, (fr(1, "A") + fr(2, "A")) / 2, c59 / 5 }'

# ---------------------------------------------------------------- gate A: is the artefact in these reads?
read -r GN PRE_C12 PRE_A12 PRE_C59 < <("$PIGZ" -dc "$HR2" 2>/dev/null | head -n $((4 * GATE_N)) \
  | awk -v tab="$W/r2_pos1-10_before_clip.tsv" "$COMP") || die "gate A produced nothing"
cat "$W/r2_pos1-10_before_clip.tsv"
log "gate A on $GN records: mean C pos1-2 $PRE_C12 (need < $G_C12), mean A pos1-2 $PRE_A12 (need < $G_A12), mean C pos5-9 $PRE_C59 (need > $G_C59)"
[ "$GN" = "$GATE_N" ] || fail "gate A saw $GN records, wanted $GATE_N"
if ! awk -v c="$PRE_C12" -v a="$PRE_A12" -v c5="$PRE_C59" -v cm="$G_C12" -v am="$G_A12" -v c5m="$G_C59" \
     'BEGIN { exit !(c + 0 < cm + 0 && a + 0 < am + 0 && c5 + 0 > c5m + 0) }'; then
  echo "[$(ts)] GATE A FAILED $S: R2 5' artefact not seen; clipping would cut template. Nothing written."
  exit 3
fi

# ---------------------------------------------------------------- outputs (each written as .part, then mv)
log "1/3 graft R1: first $((4 * GRAFT)) lines of human R1, unchanged"
"$PIGZ" -dc "$HR1" 2>/dev/null | head -n $((4 * GRAFT)) | "$PIGZ" -p "$PT" > "$GR1.part"
mv -f "$GR1.part" "$GR1"

log "2/3 R2: 3 nt off every sequence and quality line; all $N records -> human_clip, first $GRAFT -> graft_clip"
set +e
"$PIGZ" -dc "$HR2" | awk -v ng4=$((4 * GRAFT)) -v side="$W/r2_clip.txt" \
    -v g="'$PIGZ' -p 4 > '$GR2.part'" -v h="'$PIGZ' -p 4 > '$HC2.part'" '
  NR % 4 == 2 { L = length($0); if (mn == "" || L < mn) mn = L; if (L > mx) mx = L; sum += L; if (L < 4) short++ }
  { if (NR % 4 == 2 || NR % 4 == 0) s = substr($0, 4); else s = $0
    print s | h
    if (NR <= ng4) print s | g }
  END { rg = close(g); rh = close(h)
        printf "%.0f %s %s %.0f %d %d %d\n", NR, (mn == "" ? "NA" : mn), mx, sum, short, rg, rh > side }'
ST="${PIPESTATUS[*]}"
set -e
[ "$ST" = "0 0" ] || fail "R2 clip pipeline exit statuses '$ST' (decompress, awk)"
read -r R2_LINES RAW_MIN RAW_MAX RAW_SUM R2_SHORT RG RH <<< "$(rd "$W/r2_clip.txt")"
[ "$RG" = 0 ] && [ "$RH" = 0 ] || fail "R2 compressors exited '$RG' (graft) and '$RH' (human_clip)"
mv -f "$GR2.part" "$GR2"
mv -f "$HC2.part" "$HC2"

log "3/3 R1 scan: names of the both records, index field of every header, lanes, boundary record"
set +e
"$PIGZ" -dc "$HR1" | awk -v ng="$GRAFT" -v side="$W/r1_scan.txt" -v idxf="$IDX.part" \
    -v tally="$W/r1_index_tally.tsv" -v nz="'$PIGZ' -p $PT > '$NAMES.part'" '
  NR % 4 == 1 { k++
    if (substr($0, 1, 1) != "@") badh++
    nm = $1; sub(/^@/, "", nm)
    cm = $0; if (!sub(/^[^ \t]+[ \t]+/, "", cm)) cm = ""
    if (k == 1) { first = cm; firstname = nm }
    last = cm; lastname = nm
    split(cm, w, /[ \t]+/); nf = split(w[1], a, ":"); ix = (nf >= 4 && a[nf] != "") ? a[nf] : "NA"; idx[ix]++
    split(nm, b, ":"); lane[b[4]]++
    if (k > ng) { print nm | nz; if (k == ng + 1) bnd = nm } }
  END { rz = close(nz)
        best = "NA"; bn = 0; nd = 0
        print "index_field\trecords" > tally
        for (x in idx) { nd++; printf "%s\t%.0f\n", x, idx[x] > tally; if (idx[x] > bn) { bn = idx[x]; best = x } }
        close(tally)
        ln = ""; for (x in lane) ln = ln (ln == "" ? "" : ",") "L" x "=" lane[x]
        if (ln == "") ln = "NA"
        print "key\tread_name_or_value\theader_comment_or_note" > idxf
        print "first\t" firstname "\t" first > idxf
        print "last\t" lastname "\t" last > idxf
        printf "modal_index\t%s\t%.0f of %.0f records (%.6f); %d distinct index fields\n", best, bn, k, (k ? bn / k : 0), nd > idxf
        print "lanes\t" ln "\tread count per lane, from the read names" > idxf
        close(idxf)
        printf "%.0f %.0f %d %s %s %d %s %s\n", NR, k, badh, (bnd == "" ? "NA" : bnd), (lastname == "" ? "NA" : lastname), rz, best, ln > side }'
ST="${PIPESTATUS[*]}"
set -e
[ "$ST" = "0 0" ] || fail "R1 scan pipeline exit statuses '$ST' (decompress, awk)"
read -r R1_LINES R1_REC R1_BADH BND R1_LAST RZ MODAL LANES <<< "$(rd "$W/r1_scan.txt")"
[ "$RZ" = 0 ] || fail "names compressor exited '$RZ'"
mv -f "$NAMES.part" "$NAMES"
mv -f "$IDX.part" "$IDX"
cat "$IDX"

# ---------------------------------------------------------------- checks (read back with gzip's zcat, not pigz)
log "checks, part 1: byte identity of the graft files (e)"
( zcat "$HR1" 2>/dev/null | head -n $((4 * GRAFT)) | md5sum | cut -c1-32 > "$W/e_src.md5" ) &
( zcat "$GR1" | md5sum | cut -c1-32 > "$W/e_out.md5" ) &
( zcat "$HC2" 2>/dev/null | head -n $((4 * GRAFT)) | md5sum | cut -c1-32 > "$W/e2_src.md5" ) &
( zcat "$GR2" | md5sum | cut -c1-32 > "$W/e2_out.md5" ) &
wait

# shape of a clipped R2 file: records, last name, min, max, total bases, bad lengths, bad quality lengths, bad lines
SHAPE='NR % 4 == 1 { if (substr($0, 1, 1) != "@") bad++; nm = $1 }
NR % 4 == 2 { L = length($0); if (L < 1 || L > maxl) badlen++; if (mn == "" || L < mn) mn = L; if (L > mx) mx = L; sum += L }
NR % 4 == 3 { if (substr($0, 1, 1) != "+") bad++ }
NR % 4 == 0 { if (length($0) != L) badq++ }
END { sub(/^@/, "", nm)
      printf "%.0f %s %s %s %.0f %d %d %d\n", NR, (nm == "" ? "NA" : nm), (mn == "" ? "NA" : mn), mx, sum, badlen, badq, bad }'
# the same, plus every record name compared with its R1 mate, read in step from the human R1 file
PAIR='NR % 4 == 1 { if (substr($0, 1, 1) != "@") bad++; nm = $1
         got = 0
         for (j = 1; j <= 4; j++) if ((r1cmd | getline x) > 0) { got++; if (j == 1) h1 = x }
         if (got < 4) r1short++
         else { split(h1, q, /[ \t]+/); n1 = q[1]; n2 = nm; sub(/\/[12]$/, "", n1); sub(/\/[12]$/, "", n2); if (n1 != n2) mism++ } }
NR % 4 == 2 { L = length($0); if (L < 1 || L > maxl) badlen++; if (mn == "" || L < mn) mn = L; if (L > mx) mx = L; sum += L }
NR % 4 == 3 { if (substr($0, 1, 1) != "+") bad++ }
NR % 4 == 0 { if (length($0) != L) badq++ }
END { extra = 0; while ((r1cmd | getline x) > 0) extra++; close(r1cmd)
      sub(/^@/, "", nm)
      printf "%.0f %s %s %s %.0f %d %d %d %d %d %d\n", NR, (nm == "" ? "NA" : nm), (mn == "" ? "NA" : mn), mx, sum, badlen, badq, bad, mism, r1short, extra }'

log "checks, part 2: counts, boundary, pairing, clip shape, clip effect (a-d, f, g)"
( zcat "$GR1" | awk 'NR % 4 == 1 { nm = $1; if (substr($0, 1, 1) != "@") bad++ }
    END { sub(/^@/, "", nm); printf "%.0f %s %d\n", NR, (nm == "" ? "NA" : nm), bad }' > "$W/p1.txt" ) &
( zcat "$GR2" | awk -v maxl="$MAXL" "$SHAPE" > "$W/p2.txt" ) &
( zcat "$HC2" | awk -v maxl="$MAXL" -v r1cmd="zcat '$HR1'" "$PAIR" > "$W/p3.txt" ) &
( zcat "$RR1" 2>/dev/null | awk -v L=$((4 * HOST + 1)) 'NR == L { nm = $1; sub(/^@/, "", nm); print nm; exit }' > "$W/c_rat.txt" ) &
( zcat "$NAMES" | awk 'NR == 1 { f = $0 } { l = $0 } END { printf "%.0f %s %s\n", NR, (f == "" ? "NA" : f), (l == "" ? "NA" : l) }' > "$W/names.txt" ) &
( zcat "$HC2" 2>/dev/null | head -n $((4 * GATE_N)) | awk -v tab="$W/r2_pos1-10_after_clip.tsv" "$COMP" > "$W/g.txt" ) &
wait

E_SRC=$(rd "$W/e_src.md5");  E_OUT=$(rd "$W/e_out.md5")
E2_SRC=$(rd "$W/e2_src.md5"); E2_OUT=$(rd "$W/e2_out.md5")
read -r P1_LINES P1_LAST P1_BAD <<< "$(rd "$W/p1.txt")"
read -r P2_LINES P2_LAST P2_MIN P2_MAX P2_SUM P2_BADLEN P2_BADQ P2_BAD <<< "$(rd "$W/p2.txt")"
read -r P3_LINES P3_LAST P3_MIN P3_MAX P3_SUM P3_BADLEN P3_BADQ P3_BAD P3_MISM P3_R1SHORT P3_R1EXTRA <<< "$(rd "$W/p3.txt")"
C_RAT=$(rd "$W/c_rat.txt")
read -r NM_LINES NM_FIRST NM_LAST <<< "$(rd "$W/names.txt")"
read -r GN2 POST_C12 POST_A12 POST_C59 <<< "$(rd "$W/g.txt")"
cat "$W/r2_pos1-10_after_clip.tsv" 2>/dev/null || true

FAILS=""
ck() {   # ck <id> <description> <test command...>
  local id=$1 desc=$2; shift 2
  if "$@" 2>/dev/null; then echo "  ok    $id  $desc"; else echo "  FAIL  $id  $desc"; FAILS="$FAILS $id"; fi
}
isn()  { [[ $1 =~ ^[0-9]+$ ]]; }
eq()   { isn "$1" && isn "$2" && [ "$1" -eq "$2" ]; }
zero_all() { local v; for v in "$@"; do eq "$v" 0 || return 1; done; }
shifted() { isn "$1" && isn "$2" && isn "$3" && [ "$1" -eq $(( $2 - $3 )) ]; }   # $1 = $2 - $3
same() { [ -n "$1" ] && [ "$1" != NA ] && [ "$1" = "$2" ]; }
mate() { local x=${1%/[12]} y=${2%/[12]}; [ -n "$x" ] && [ "$x" != NA ] && [ "$x" = "$y" ]; }
md5eq() { [[ $1 =~ ^[0-9a-f]{32}$ ]] && [ "$1" = "$2" ]; }
gt()   { awk -v a="$1" -v b="$2" 'BEGIN { exit !(a != "" && a != "NA" && a + 0 > b + 0) }'; }

log "check results"
ck a1 "human R1 lines $R1_LINES = 4 x $N"                          eq "$R1_LINES" $((4 * N))
ck a2 "human R1 records $R1_REC = $N"                              eq "$R1_REC" "$N"
ck a3 "human R1 headers all start with @ ($R1_BADH bad)"          eq "$R1_BADH" 0
ck a4 "human R2 lines $R2_LINES = 4 x $N"                          eq "$R2_LINES" $((4 * N))
ck b1 "graft R1 lines $P1_LINES = 4 x $GRAFT, headers not @: $P1_BAD" eq "$P1_LINES" $((4 * GRAFT))
ck b1h "graft R1 headers all start with @"                         eq "$P1_BAD" 0
ck b2 "graft R2 lines $P2_LINES = 4 x $GRAFT"                      eq "$P2_LINES" $((4 * GRAFT))
ck b3 "human_clip R2 lines $P3_LINES = 4 x $N"                     eq "$P3_LINES" $((4 * N))
ck b4 "both names $NM_LINES = $BOTH"                                eq "$NM_LINES" "$BOTH"
ck c1 "boundary: human R1 record $((GRAFT + 1)) '$BND' = rat R1 record $((HOST + 1)) '$C_RAT'" same "$BND" "$C_RAT"
ck c2 "first both name '$NM_FIRST' = boundary record '$BND'"       same "$NM_FIRST" "$BND"
ck c3 "last both name '$NM_LAST' = last human R1 name '$R1_LAST'"  same "$NM_LAST" "$R1_LAST"
ck d1 "last names: graft R1 '$P1_LAST', graft R2 '$P2_LAST'"       mate "$P1_LAST" "$P2_LAST"
ck d2 "last names: human R1 '$R1_LAST', human_clip R2 '$P3_LAST'"  mate "$R1_LAST" "$P3_LAST"
ck d3 "every human_clip R2 name equals its R1 mate ($P3_MISM mismatches)" eq "$P3_MISM" 0
ck d4 "R1 ran out $P3_R1SHORT times, had $P3_R1EXTRA extra lines"  zero_all "$P3_R1SHORT" "$P3_R1EXTRA"
ck e1 "md5 human R1 head $E_SRC = graft R1 $E_OUT"                  md5eq "$E_SRC" "$E_OUT"
ck e2 "md5 human_clip R2 head $E2_SRC = graft R2 $E2_OUT"           md5eq "$E2_SRC" "$E2_OUT"
ck f1 "graft R2: lengths 1-$MAXL (bad $P2_BADLEN), qual = seq length (bad $P2_BADQ), line types (bad $P2_BAD)" \
      zero_all "$P2_BADLEN" "$P2_BADQ" "$P2_BAD"
ck f2 "human_clip R2: lengths 1-$MAXL (bad $P3_BADLEN), qual = seq length (bad $P3_BADQ), line types (bad $P3_BAD)" \
      zero_all "$P3_BADLEN" "$P3_BADQ" "$P3_BAD"
ck f3 "no raw R2 shorter than 4 nt ($R2_SHORT)"                     eq "$R2_SHORT" 0
ck f4 "clipped min $P3_MIN = raw min $RAW_MIN - 3"                  shifted "$P3_MIN" "$RAW_MIN" 3
ck f5 "clipped max $P3_MAX = raw max $RAW_MAX - 3"                  shifted "$P3_MAX" "$RAW_MAX" 3
ck f6 "clipped bases $P3_SUM = raw bases $RAW_SUM - 3 x $N"         shifted "$P3_SUM" "$RAW_SUM" $((3 * N))
ck g1 "clipped R2 mean C pos1-2 $POST_C12 > $P_C12 (was $PRE_C12) on $GN2 records" gt "$POST_C12" "$P_C12"
log "R2 lengths: raw $RAW_MIN-$RAW_MAX; clipped human_clip $P3_MIN-$P3_MAX, graft $P2_MIN-$P2_MAX (expected 32 and 98)"

if [ -n "$FAILS" ]; then
  fail "failed:$FAILS"
fi

# ---------------------------------------------------------------- record
mv -f "$W/r2_pos1-10_before_clip.tsv" "$R/checks/$S.r2_pos1-10_before_clip.tsv"
mv -f "$W/r2_pos1-10_after_clip.tsv" "$R/checks/$S.r2_pos1-10_after_clip.tsv"
mv -f "$W/r1_index_tally.tsv" "$R/index_seq/$S.index_tally.tsv"
H=(lib status n_host n_graft n_both human_R1_records human_R2_records graft_R1_records graft_R2_records
   human_clip_R2_records both_names boundary_name_human boundary_name_rat last_name_graft_R1 last_name_graft_R2
   last_name_human_R1 last_name_human_clip_R2 pair_name_mismatches md5_humanR1_head md5_graft_R1
   md5_humanclipR2_head md5_graft_R2 raw_R2_len_min raw_R2_len_max clip_R2_len_min clip_R2_len_max
   graft_R2_len_min graft_R2_len_max raw_R2_bases clip_R2_bases pre_C12 pre_A12 pre_C59 post_C12 post_A12
   post_C59 modal_index lanes bytes_graft_R1 bytes_graft_R2 bytes_human_clip_R2 bytes_names job)
V=("$S" OK "$HOST" "$GRAFT" "$BOTH" "$R1_REC" $((R2_LINES / 4)) $((P1_LINES / 4)) $((P2_LINES / 4))
   $((P3_LINES / 4)) "$NM_LINES" "$BND" "$C_RAT" "$P1_LAST" "$P2_LAST"
   "$R1_LAST" "$P3_LAST" "$P3_MISM" "$E_SRC" "$E_OUT"
   "$E2_SRC" "$E2_OUT" "$RAW_MIN" "$RAW_MAX" "$P3_MIN" "$P3_MAX"
   "$P2_MIN" "$P2_MAX" "$RAW_SUM" "$P3_SUM" "$PRE_C12" "$PRE_A12" "$PRE_C59" "$POST_C12" "$POST_A12"
   "$POST_C59" "$MODAL" "$LANES" "$(stat -c %s "$GR1")" "$(stat -c %s "$GR2")" "$(stat -c %s "$HC2")" "$(stat -c %s "$NAMES")"
   "${SLURM_ARRAY_JOB_ID:-NA}_$T")
[ "${#H[@]}" -eq "${#V[@]}" ] || fail "internal: ${#H[@]} column names, ${#V[@]} values"
( IFS=$'\t'; echo "${H[*]}"; echo "${V[*]}" ) > "$CHK.part"
mv -f "$CHK.part" "$CHK"
rm -rf "$W"
log "done in $((SECONDS / 60)) min"
echo "CHECK OK $S"
