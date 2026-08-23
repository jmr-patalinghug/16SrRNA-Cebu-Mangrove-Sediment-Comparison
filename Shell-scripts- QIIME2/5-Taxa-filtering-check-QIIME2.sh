#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task 10
#SBATCH --time=24:00:00
#SBATCH --job-name=CleanCheck
#SBATCH --output=qiime_cleancheck_%j.out
#SBATCH --error=qiime_cleancheck_%j.err

# Notes about this script:
#
# Purpose: Compare the unfiltered and taxonomy-filtered feature tables to determine the numbers of ASVs and total reads removed during taxonomic filtering.
#
# Inputs:
#   - Unfiltered feature table (table.qza)
#   - Taxonomy-filtered feature table (table.phylumfiltered.qza)
#
# Outputs:
#   - Exported unfiltered feature table (feature-table.biom)
#   - Exported filtered feature table (feature-table.biom)
#   - BIOM summary of the unfiltered feature table
#     (table_biom_summary.txt)
#   - BIOM summary of the filtered feature table
#     (table_phylumfiltered_biom_summary.txt)
#
# Requirements:
#   - QIIME 2 (qiime2-amplicon-2024.10)
#   - BIOM Format command-line tools
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


set -euo pipefail

export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Export the unfiltered feature table as a BIOM file.
qiime tools export \
  --input-path /mnt/scratch/wisnoskilab/jrp889/Output/table.qza \
  --output-path /mnt/scratch/wisnoskilab/jrp889/Output/exported_table

# Export the taxonomy-filtered feature table as a BIOM file.
qiime tools export \
  --input-path /mnt/scratch/wisnoskilab/jrp889/Output/table.phylumfiltered.qza \
  --output-path /mnt/scratch/wisnoskilab/jrp889/Output/exported_table_phylumfiltered

# Summarize the unfiltered BIOM feature table.
biom summarize-table \
  -i /mnt/scratch/wisnoskilab/jrp889/Output/exported_table/feature-table.biom \
  -o /mnt/scratch/wisnoskilab/jrp889/Output/table_biom_summary.txt

# Summarize the taxonomy-filtered BIOM feature table.
biom summarize-table \
  -i /mnt/scratch/wisnoskilab/jrp889/Output/exported_table_phylumfiltered/feature-table.biom \
  -o /mnt/scratch/wisnoskilab/jrp889/Output/table_phylumfiltered_biom_summary.txt

# Interpretation:
# In each summary file, "Num observations" is the number of ASVs, while
# "Total count" is the total number of reads represented in the feature table.
#
# Subtract the filtered-table values from the corresponding unfiltered-table
# values to calculate the numbers of ASVs and reads removed:
#
#   ASVs removed = unfiltered Num observations - filtered Num observations
#
#   Reads removed = unfiltered Total count - filtered Total count