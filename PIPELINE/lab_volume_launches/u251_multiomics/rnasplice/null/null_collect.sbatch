#!/bin/bash
# null_collect: afterany on both null arrays. Metrics, task 1 against the pipeline's tables, ranks, leave-one-out,
# NULL_REPORT.md. Deletes $D/work (about 200 GB) only when Nextflow exited 0, every hard check passed, all 16 tasks
# finished every method and task 1 reproduced every pipeline table (CLEAN_WORK=0 keeps it).
#SBATCH --job-name=u251_rs_null_collect
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnasplice
N=$D/null
CLEAN_WORK=${CLEAN_WORK:-1}
echo "null_collect: host $(hostname), job ${SLURM_JOB_ID:-none}, start $(date -Is)"
PY=""
for p in "$HOME/mambaforge/bin/python3" "$HOME/mambaforge/envs/spineps/bin/python" "$(command -v python3 || true)"; do
  if [ -n "$p" ] && [ -x "$p" ]; then PY=$p; break; fi
done
"$PY" "$N/null_collect.py"
rc=$?
echo "null_collect exit $rc (0 = complete and reproduced)"
sed -n '1,40p' "$N/results/NULL_REPORT.md" 2>/dev/null || true
HP=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("hard_pass"))' "$D/checks.json" 2>/dev/null || echo missing)
nfx=$(cat "$D/NEXTFLOW_EXIT" 2>/dev/null || echo missing)
if [ "$rc" -eq 0 ] && [ "$HP" = "True" ] && [ "$nfx" = "0" ] && [ "$CLEAN_WORK" = "1" ] && [ -d "$D/work" ]; then
  sz=$(getfattr --only-values -n ceph.dir.rbytes "$D/work" 2>/dev/null || echo unknown)
  echo "$(date -Is) null complete and reproduced, hard checks passed, Nextflow exit 0: deleting $D/work ($sz bytes)"
  rm -rf "$D/work" && echo "work deleted $(date -Is)" | tee -a "$N/results/NULL_REPORT.md"
else
  echo "work kept (collect rc $rc, hard_pass $HP, nextflow exit $nfx, CLEAN_WORK $CLEAN_WORK)"
fi
exit $rc
