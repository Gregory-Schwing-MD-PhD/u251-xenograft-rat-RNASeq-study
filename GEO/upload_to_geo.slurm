#!/bin/bash
#SBATCH --job-name=geo_upload
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -t 24:00:00
#SBATCH -o geo_upload_%j.out
#SBATCH -e geo_upload_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=go2432@wayne.edu
#
# Uploads the 20 raw FASTQs + processed counts matrix to GEO's FTP server.
# The FTP password is NOT stored here -- it is read from ~/.geo_ftp_pass so
# this script can be committed/shared safely. Create that file once with:
#     printf '%s' 'YOUR_GEO_FTP_PASSWORD' > ~/.geo_ftp_pass && chmod 600 ~/.geo_ftp_pass
#
# Re-running is safe: `mput -c` / `put -c` resume partially transferred files.

set -euo pipefail

# ---------------- configuration ----------------
FTP_HOST="ftp-private.ncbi.nlm.nih.gov"
FTP_USER="geoftp"
REMOTE_SPACE="uploads/go2432@wayne.edu_4TOluI24"   # your personalized GEO upload space
SUBFOLDER="U251_RNAseq"                             # dataset subfolder (select this at metadata step)
PASS_FILE="$HOME/.geo_ftp_pass"

REPO="$HOME/u251-xenograft-murine-RNASeq-study"
FASTQ_DIR="$REPO/DATA/RNASEQ/RAW/downloaded_files"
PROC_FILE="$REPO/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"

# ---------------- pre-flight checks ----------------
if [[ ! -f "$PASS_FILE" ]]; then
  echo "ERROR: password file $PASS_FILE not found." >&2
  echo "Create it with:  printf '%s' 'YOUR_GEO_FTP_PASSWORD' > $PASS_FILE && chmod 600 $PASS_FILE" >&2
  exit 1
fi
FTP_PASS="$(cat "$PASS_FILE")"

module load lftp 2>/dev/null || true
if ! command -v lftp >/dev/null 2>&1; then
  echo "ERROR: lftp not found (try a different module name or install)." >&2
  exit 1
fi

[[ -d "$FASTQ_DIR" ]] || { echo "ERROR: FASTQ dir not found: $FASTQ_DIR" >&2; exit 1; }
[[ -f "$PROC_FILE" ]] || { echo "ERROR: processed file not found: $PROC_FILE" >&2; exit 1; }

n_fastq=$(ls -1 "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
echo "Found $n_fastq FASTQ files + 1 processed file to upload."
echo "Host=$FTP_HOST  Space=$REMOTE_SPACE  Subfolder=$SUBFOLDER"
echo "Job started on host $(hostname) at $(date)"

# ---------------- connectivity check (fails fast if compute node is firewalled) ----------------
if ! lftp -u "$FTP_USER,$FTP_PASS" "$FTP_HOST" -e "set ftp:ssl-allow no; set net:timeout 25; set net:max-retries 1; cd $REMOTE_SPACE; ls; bye" >/dev/null 2>&1; then
  echo "ERROR: could not reach $FTP_HOST from this node." >&2
  echo "Compute nodes are likely firewalled from the internet. Run the upload from a" >&2
  echo "login node inside screen/tmux instead, or use a data-transfer node/partition." >&2
  exit 2
fi
echo "Connectivity OK -- starting transfer at $(date)"

# ---------------- transfer ----------------
lftp -u "$FTP_USER,$FTP_PASS" "$FTP_HOST" <<EOF
set ftp:ssl-allow no
set net:max-retries 5
set net:timeout 30
set net:reconnect-interval-base 10
cd $REMOTE_SPACE
mkdir -f $SUBFOLDER
cd $SUBFOLDER
lcd $FASTQ_DIR
mput -c *.fastq.gz
lcd $(dirname "$PROC_FILE")
put -c $(basename "$PROC_FILE")
echo "--- remote listing after upload ---"
ls
bye
EOF

echo "Upload finished at $(date)"
echo "Verify the remote listing above shows $((n_fastq + 1)) files, then upload the"
echo "metadata .xlsx on the GEO web 'Upload metadata' page (subfolder: $SUBFOLDER)."
