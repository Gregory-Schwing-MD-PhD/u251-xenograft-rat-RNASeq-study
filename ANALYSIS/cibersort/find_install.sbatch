#!/bin/bash
# Where, if anywhere, is a CIBERSORT / CIBERSORTx install on the grid? Searches home and the CephFS volume without a
# time cap (the login-node search was cut off), for files and directories named like CIBERSORT, container images
# and definition files, and R sources that define the CIBERSORT function.
#   sbatch -D ANALYSIS/cibersort ANALYSIS/cibersort/find_install.sbatch
#SBATCH --job-name=u251_find_ciber
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/find_install_%j.out
set -uo pipefail
ROOTS="$HOME /rs/rs_grp_oschome/go2432"
echo "== names containing cibersort, LM22 or fractions (files and directories)"
find $ROOTS -xdev \( -path '*/work' -o -path '*/.nextflow' -o -path '*/node_modules' -o -path '*/.git/objects' \) -prune -o \
  \( -iname '*cibersort*' -o -iname 'LM22*' -o -iname '*cibersortx*' \) -print 2>/dev/null
echo "== container images (.sif .simg) and definition files (.def, Singularity*)"
find $ROOTS -xdev \( -path '*/work' -o -path '*/.nextflow' -o -path '*/node_modules' -o -path '*/.git/objects' \) -prune -o \
  \( -name '*.sif' -o -name '*.simg' -o -name '*.def' -o -name 'Singularity' -o -name 'Singularity.*' \) -print 2>/dev/null
echo "== sandbox containers (directories holding .singularity.d)"
find $ROOTS -xdev -maxdepth 6 -type d -name '.singularity.d' 2>/dev/null
echo "== R or text sources that define CIBERSORT (size < 5 MB)"
find $ROOTS -xdev \( -path '*/work' -o -path '*/.nextflow' -o -path '*/node_modules' -o -path '*/.git/objects' \) -prune -o \
  -type f \( -name '*.R' -o -name '*.r' -o -name '*.txt' -o -name '*.def' -o -name 'Dockerfile*' -o -name '*.dockerfile' \) -size -5M -print 2>/dev/null \
  | xargs -r -d '\n' grep -l -i -E 'CIBERSORT *(<-|=) *function|cibersortx|CoreAlg *(<-|=)' 2>/dev/null
echo "== R libraries with a CIBERSORT package"
find $ROOTS -xdev -maxdepth 8 -type d -iname 'CIBERSORT*' -path '*library*' 2>/dev/null
echo "== containers: R packages inside each .sif found above"
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
for s in $(find $ROOTS -xdev \( -path '*/work' -o -path '*/.nextflow' \) -prune -o \( -name '*.sif' -o -name '*.simg' \) -print 2>/dev/null); do
  hit=$($SING exec "$s" sh -c 'ls -d /usr/local/lib/R/site-library/* /usr/lib/R/site-library/* /usr/local/lib/R/library/* /opt/conda/lib/R/library/* 2>/dev/null | grep -i -E "cibersort|immunedeconv"; ls / /opt 2>/dev/null | grep -i ciber' 2>/dev/null)
  echo "$s: ${hit:-none}"
done
echo "== done"
