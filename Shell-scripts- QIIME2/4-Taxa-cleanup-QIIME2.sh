#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task 10
#SBATCH --time=24:00:00
#SBATCH --job-name=Cleaning-unassigned
#SBATCH --output=qiime_cleanup_%j.out
#SBATCH --error=qiime_cleanup_%j.err

# Notes about this script:
#
# Purpose: Retain bacterial ASVs classified to at least the phylum level; remove chloroplast, mitochondrial, and unassigned ASVs; filter the corresponding representative sequences; and generate a taxonomic bar-plot visualization.

# Inputs:
#   - Feature table (table.qza)
#   - Taxonomic classifications (taxonomy.qza)
#   - Representative sequences (repseqs.qza)
#   - Sample metadata (usc-sediment-metadata.tsv)

# Outputs:
#   - Filtered feature table (table.phylumfiltered.qza)
#   - Filtered representative sequences (rep-seqs.phylumfiltered.qza)
#   - Taxonomic bar-plot visualization (taxa-bar-plots-cleaned.qzv)
#
# Requirements: QIIME 2 (qiime2-amplicon-2024.10)
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.

set -euo pipefail
export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Inclusion criterion retains bacterial ASVs with a phylum assignment.
# Excludes Archaea and ASVs classified only at the domain level.

qiime taxa filter-table \
    --i-table /mnt/scratch/wisnoskilab/jrp889/Output/table.qza \
    --i-taxonomy /mnt/scratch/wisnoskilab/jrp889/Output/taxonomy.qza \
    --p-mode contains \
    --p-include 'd__Bacteria;p__' \ #this part also excludes archaea since onldy domain bacteria will be included in the final feature table
    --p-exclude 'c__Chloroplast,f__mitochondria,Unassigned' \
    --o-filtered-table /mnt/scratch/wisnoskilab/jrp889/Output/table.phylumfiltered.qza \
   

qiime taxa filter-seqs \ #cleaning the rep seqs similarly
    --i-sequences /mnt/scratch/wisnoskilab/jrp889/Output/repseqs.qza \
    --i-taxonomy /mnt/scratch/wisnoskilab/jrp889/Output/taxonomy.qza \
    --p-mode contains \
    --p-include 'd__Bacteria;p__' \
    --p-exclude 'c__Chloroplast,f__mitochondria,Unassigned' \
    --o-filtered-sequences /mnt/scratch/wisnoskilab/jrp889/Output/rep-seqs.phylumfiltered.qza \
   

#Taxa barplots, cleaned
qiime taxa barplot \
    --i-table /mnt/scratch/wisnoskilab/jrp889/Output/table.phylumfiltered.qza \
    --i-taxonomy /mnt/scratch/wisnoskilab/jrp889/Output/taxonomy.qza \
    --m-metadata-file /mnt/scratch/wisnoskilab/jrp889/Analysis/usc-sediment-metadata.tsv --o-visualization /mnt/scratch/wisnoskilab/jrp889/Output/taxa-bar-plots-cleaned.qzv 