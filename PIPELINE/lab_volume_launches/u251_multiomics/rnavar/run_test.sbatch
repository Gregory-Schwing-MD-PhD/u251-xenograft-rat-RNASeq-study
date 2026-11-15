#!/bin/bash
# rnavar -profile test on the grid, before the real run: Nextflow 25.10.7 from a compute node, the plugins, the
# nf-core/rnavar 1.3.0 pull, slurm submission through grid.config (reqp, requeue, errorStrategy closure), Singularity
# pulls into the shared cache. Separate outdir, work, launch directory and NXF_HOME, so it can never lock or block the
# real run; nothing depends on it. Writes $T/TEST_EXIT.
#SBATCH --job-name=u251_rnavar_test
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnavar/logs/%x_%j.out
set -uo pipefail
umask 0002
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnavar
T=$D/test
echo "rnavar test profile: host $(hostname), job ${SLURM_JOB_ID:-none}, start $(date -Is)"
mkdir -p "$T/launch" "$T/nxf_home"
export NXF_VER=25.10.7
export NXF_BIN=$M/nxf_dist/$NXF_VER/nextflow-$NXF_VER-one.jar
export NXF_HOME=$T/nxf_home
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
export NXF_TEMP=$M/tmp TMPDIR=$M/tmp
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export NXF_OPTS='-Xms1g -Xmx4g'
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH
mkdir -p "$NXF_TEMP" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR" "$(dirname "$NXF_BIN")"
nextflow -version 2>&1 | grep -i version
cd "$T/launch"
nextflow run nf-core/rnavar -r 1.3.0 \
  -profile test,singularity \
  -c ../../grid.config \
  -c ../../test_limits.config \
  --outdir "$T/results" \
  -w "$T/work" \
  -resume \
  -ansi-log false
rc=$?
echo "$rc" > "$T/TEST_EXIT"
echo "test exit $rc at $(date -Is)"
git -C "$NXF_HOME/assets/nf-core/rnavar" rev-parse HEAD 2>/dev/null | sed 's/^/pipeline commit: /'
[ "$rc" = 0 ] && rm -rf "$T/work"     # keep the work directory only when the test failed
exit $rc
