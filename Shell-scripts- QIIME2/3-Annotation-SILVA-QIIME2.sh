#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --job-name=Annotation-SILVA
#SBATCH --output=qiime_taxonomy_%j.out
#SBATCH --error=qiime_taxonomy_%j.err

# Notes about this script:
#
# Purpose: Assign taxonomic classifications to representative sequences using a SILVA reference classifier.
# Input: Representative sequences generated during the DADA2 denoising step.
# Output: Taxonomic classifications (taxonomy.qza).
# Requirements: QIIME 2 (qiime2-amplicon-2024.10) and a trained SILVA 138 99% Naive Bayes classifier.
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.

set -euo pipefail

export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Assign taxonomy to representative sequences using a SILVA 138
# 99% reference classifier.
qiime feature-classifier classify-sklearn \
  --i-classifier /mnt/scratch/wisnoskilab/dc2484/reference_databases/qiime2_classifiers/silva-138-99-nb-classifier.qza \
  --i-reads /mnt/scratch/wisnoskilab/jrp889/Output/repseqs.qza \
  --o-classification /mnt/scratch/wisnoskilab/jrp889/Output/taxonomy.qza