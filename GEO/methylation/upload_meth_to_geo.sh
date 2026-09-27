#!/bin/bash
# Upload the eight EPIC arrays to GEO: 16 raw IDATs + 4 processed files, 179 MB, into ONE flat folder.
#
# Modelled on GEO/upload_to_geo.slurm, which is what actually worked for the RNA submission (GSE338105): same host,
# same personalised space, same lftp resume flags. Two deliberate differences:
#
#   1. It runs as a plain shell script, not an sbatch job. The RNA script's own connectivity check found that
#      ftp-private.ncbi.nlm.nih.gov is firewalled from the compute nodes, and its error text says to run from a login
#      node in screen/tmux instead. This is a 179 MB network transfer, not compute, so a login node is the right place.
#      It still works under sbatch if the node can reach the host.
#   2. GEO wants one flat folder per submission, so the idat/ and processed/ subdirectories are flattened here.
#
#   ssh go2432@grid.wayne.edu
#   screen -S geo            # or tmux
#   bash ~/u251-xenograft-murine-RNASeq-study/GEO/methylation/upload_meth_to_geo.sh
#
# Re-running is safe: mput -c / put -c resume partially transferred files.
set -euo pipefail

FTP_HOST="ftp-private.ncbi.nlm.nih.gov"
FTP_USER="geoftp"
REMOTE_SPACE="uploads/go2432@wayne.edu_4TOluI24"   # the same personalised space as the RNA submission
SUBFOLDER="U251_Methylation"                        # select this subfolder at the metadata step on the portal
PASS_FILE="$HOME/.geo_ftp_pass"                     # same file the RNA upload used; never printed, never committed
SRC=/rs/rs_grp_oschome/go2432/u251_geo_meth
export TMPDIR=/rs/rs_grp_oschome/go2432/tmp

[ -f "$PASS_FILE" ] || { echo "ERROR: $PASS_FILE not found (the RNA upload used the same file)." >&2; exit 1; }
FTP_PASS="$(cat "$PASS_FILE")"
[ -d "$SRC/idat" ] && [ -d "$SRC/processed" ] || { echo "ERROR: staged deposit missing at $SRC" >&2; exit 1; }

echo "== 1. the payload, and its checksums, before anything leaves the machine"
n_idat=$(ls -1 "$SRC"/idat/*.idat.gz 2>/dev/null | wc -l)
n_proc=$(ls -1 "$SRC"/processed/* 2>/dev/null | wc -l)
echo "   $n_idat raw IDATs + $n_proc processed files = $((n_idat + n_proc)) (expected 16 + 4 = 20)"
[ "$n_idat" -eq 16 ] || { echo "REFUSING: expected 16 IDATs, found $n_idat" >&2; exit 1; }
[ "$n_proc" -eq 4 ]  || { echo "REFUSING: expected 4 processed files, found $n_proc" >&2; exit 1; }
( cd "$SRC" && md5sum -c md5_checksums.txt ) || { echo "REFUSING: a staged file does not match md5_checksums.txt" >&2; exit 1; }
echo "   all 20 files match their recorded md5"

echo "== 2. the sample key must not still say PENDING once accessions exist"
if grep -q PENDING "$SRC/processed/multiomics_sample_key.tsv"; then
  echo "   NOTE: multiomics_sample_key.tsv still has PENDING in array_geo_sample / array_series."
  echo "   That is expected for the FIRST upload (GEO has not issued them yet). After GEO returns the"
  echo "   accessions: set ARRAY_SERIES in GEO/make_sample_key.py, rerun it, re-stage, and replace this one file."
fi

command -v lftp >/dev/null 2>&1 || module load lftp 2>/dev/null || true
command -v lftp >/dev/null 2>&1 || { echo "ERROR: lftp not found." >&2; exit 1; }

echo "== 3. can this node reach GEO at all (fail fast rather than mid-transfer)"
if ! lftp -u "$FTP_USER,$FTP_PASS" "$FTP_HOST" \
      -e "set dns:order inet; set ftp:ssl-allow no; set net:timeout 25; set net:max-retries 1; cd $REMOTE_SPACE; ls; bye" \
      >/dev/null 2>&1; then
  echo "ERROR: cannot reach $FTP_HOST from $(hostname)." >&2
  echo "       The compute nodes are firewalled from this host. Run this on a login node inside screen/tmux." >&2
  exit 2
fi
echo "   reachable from $(hostname)"

echo "== 4. transfer, flat, into $SUBFOLDER"
lftp -u "$FTP_USER,$FTP_PASS" "$FTP_HOST" <<EOF
set dns:order "inet"
set ftp:ssl-allow no
set net:max-retries 5
set net:timeout 30
set net:reconnect-interval-base 10
cd $REMOTE_SPACE
mkdir -p $SUBFOLDER
cd $SUBFOLDER
lcd $SRC/idat
mput -c *.idat.gz
lcd $SRC/processed
mput -c *
echo "--- remote listing after upload ---"
ls
bye
EOF

echo
echo "Upload finished at $(date)"
echo "Check the listing above shows 20 files. Then on the GEO submission portal:"
echo "  * start a NEW submission of type ARRAY (not high-throughput sequencing - the EPIC BeadChip is an array"
echo "    platform, GPL21145; GEO/seq_template.xlsx is the SEQUENCING template and is the wrong one here)"
echo "  * fill GEO's own array metadata spreadsheet from GEO/methylation/metadata_methylation.xlsx"
echo "  * point the submission at subfolder: $SUBFOLDER"
echo "  * paste the requests from GEO/methylation/COVER_NOTE.md into the submission comment"
