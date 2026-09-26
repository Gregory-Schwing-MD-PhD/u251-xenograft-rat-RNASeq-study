#!/bin/bash
# Host (rat) reads for the host analysis: xengsort's `host` bin only.
# main.nf wrote <sample>_rat_R{1,2}.fq.gz = cat(host, both), in that order, and the raw FASTQs are no longer on the
# grid, so the host-only reads are exactly the first N_host records of each rat file (N_host from xengsort_out/<s>.txt).
# `both` (conserved k-mers) is left out: in tumour samples much of it is human tumour RNA, in proportion to graft
# fraction, which differs between primary and recurrent; `ambiguous` and `neither` are already absent.
# Checks: the rat file holds host + both records; record N_host+1 of the rat file is the first `both` read, which must
# equal record N_graft+1 of the human file (cat(graft, both)).
#   sbatch -D ANALYSIS/host --array=1-10 ANALYSIS/host/extract_host.sbatch
#SBATCH --job-name=u251_hostx
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/extract_%A_%a.out
set -eu   # no pipefail: zcat | head ends on SIGPIPE by design; the checks below validate the output
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
IN=$U/ANALYSIS/sorted_fastqs
OUT=/rs/rs_grp_oschome/go2432/u251_host/fastq
mkdir -p "$OUT"
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" <<'EOF'
C2B
IL64B
IL66B
IL67B
IL68B
IL69B
N168B
N269B
NL70B
NL71B
EOF
)
read -r HOST GRAFT BOTH < <(awk -v s="$S" '$1==s && NF==6 {print $2, $3, $5; exit}' "$U/ANALYSIS/xengsort_out/$S.txt")
echo "$S host=$HOST graft=$GRAFT both=$BOTH"
PIGZ=/wsu/el7/pre-compiled/python/3.7/bin/pigz
for r in 1 2; do
  zcat "$IN/${S}_rat_R$r.fq.gz" | head -n $((HOST * 4)) | "$PIGZ" -p 6 > "$OUT/${S}_host_R$r.fq.gz.part"
  mv "$OUT/${S}_host_R$r.fq.gz.part" "$OUT/${S}_host_R$r.fq.gz"
done
# checks
n_rat=$(( $(zcat "$IN/${S}_rat_R1.fq.gz" | wc -l) / 4 ))
n_out=$(( $(zcat "$OUT/${S}_host_R1.fq.gz" | wc -l) / 4 ))
first_both_rat=$(zcat "$IN/${S}_rat_R1.fq.gz" | sed -n "$((HOST * 4 + 1)){p;q}" | cut -d' ' -f1)
first_both_hum=$(zcat "$IN/${S}_human_R1.fq.gz" | sed -n "$((GRAFT * 4 + 1)){p;q}" | cut -d' ' -f1)
last_host=$(zcat "$OUT/${S}_host_R1.fq.gz" | tail -n 4 | head -1 | cut -d' ' -f1)
echo "rat records $n_rat (expected host+both $((HOST + BOTH))); host records written $n_out (expected $HOST)"
echo "first both read: rat file $first_both_rat, human file $first_both_hum; last host read $last_host"
[ "$n_rat" -eq $((HOST + BOTH)) ] && [ "$n_out" -eq "$HOST" ] && [ "$first_both_rat" = "$first_both_hum" ] \
  && echo "CHECK OK $S" || { echo "CHECK FAILED $S"; exit 2; }
