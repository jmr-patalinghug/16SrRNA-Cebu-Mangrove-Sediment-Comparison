#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --job-name=Denoise
#SBATCH --output=qiime_dada2_%j.out
#SBATCH --error=qiime_dada2_%j.err

# Notes about this script:
#
# Purpose: Denoise paired-end sequences using DADA2.
# Input: Demultiplexed paired-end sequences (.qza).
# Output: Feature table, representative sequences, and denoising statistics.
# Requirements: QIIME 2 (qiime2-amplicon-2024.10)
#
# The SLURM --output and --error directives above create separate log files containing the standard output and error messages generated during the job.
# These files were used to make troubleshooting easier on the computing cluster.
#
# For this dataset, this step took approximately 9 hours.
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


# Exit the script if a command fails, an undefined variable is used, or a command within a pipeline fails.
set -euo pipefail

# Create a temporary directory specific to this SLURM job.
export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Denoise paired-end sequences using DADA2.
#
# The trim-left values remove the 515F and 806R primer sequences (19 bp forward and 20 bp reverse, respectively).
# Truncation lengths of 150 bp (forward) and 130 bp (reverse) were selected for this dataset based on sequence quality.

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /mnt/scratch/wisnoskilab/jrp889/Output/usc-seq.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 150 \
  --p-trunc-len-r 130 \
  --p-n-threads "$SLURM_CPUS_PER_TASK" \
  --verbose \
  --o-table /mnt/scratch/wisnoskilab/jrp889/Output/table.qza \
  --o-representative-sequences /mnt/scratch/wisnoskilab/jrp889/Output/repseqs.qza \
  --o-denoising-stats /mnt/scratch/wisnoskilab/jrp889/Output/denoising_stats.qza

# The denoising statistics artifact summarizes the number of reads retained or removed during filtering, denoising, merging, and chimera removal.