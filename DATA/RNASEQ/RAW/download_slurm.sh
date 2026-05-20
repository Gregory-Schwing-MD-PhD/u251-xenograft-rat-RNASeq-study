#!/bin/bash
#
# Project/Account
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#
# Runtime set to 7 days
#SBATCH --time=7-00:00:00
#
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err

# --- LOAD TOKEN ---
# SLURM jobs run in non-interactive shells and do NOT source ~/.bashrc by default.
# Source the token file explicitly so $GDRIVE_TOKEN is available to download_files.sh.
if [ -f "$HOME/.gdrive_token" ]; then
    source "$HOME/.gdrive_token"
    echo "✅ Loaded GDRIVE_TOKEN from ~/.gdrive_token"
else
    echo "❌ ERROR: ~/.gdrive_token not found."
    echo "   Create it with:"
    echo "     echo 'export GDRIVE_TOKEN=\"ya29.a0Af...\"' > ~/.gdrive_token"
    echo "     chmod 600 ~/.gdrive_token"
    exit 1
fi

if [ -z "$GDRIVE_TOKEN" ]; then
    echo "❌ ERROR: GDRIVE_TOKEN is empty after sourcing ~/.gdrive_token"
    exit 1
fi

# --- SETUP GDOWN ---

# Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "gdown not found. Installing via pip..."
    # Install to local user directory
    pip install --user gdown

    # Add local bin to PATH so the shell can find the 'gdown' command
    export PATH="$HOME/.local/bin:$PATH"
else
    echo "gdown is already installed."
fi

# --- RUN DOWNLOAD SCRIPT ---

# Navigate to the directory where the job was submitted
cd $SLURM_SUBMIT_DIR

if [ -f "./download_files.sh" ]; then
    echo "Found download_files.sh. Making executable and running..."
    chmod +x download_files.sh
    ./download_files.sh
else
    echo "ERROR: download_files.sh not found in $SLURM_SUBMIT_DIR"
    exit 1
fi

echo "Job Complete"
