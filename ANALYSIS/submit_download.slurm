#!/bin/bash
#SBATCH --job-name=download_refs      # Job name
#SBATCH --output=download_refs_%j.log # Standard output and error log
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00               # Time limit hrs:min:sec

# 1. Ensure the script is executable
chmod +x download_phase2_refs.sh

# 2. Run the script
./download_phase2_refs.sh
