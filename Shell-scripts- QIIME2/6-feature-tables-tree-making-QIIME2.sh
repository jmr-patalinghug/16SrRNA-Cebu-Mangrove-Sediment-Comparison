#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --job-name=Standard-Feature-Tables
#SBATCH --output=qiime_featuretable_%j.out
#SBATCH --error=qiime_featuretable_%j.err

# Notes about this script:
#
# Purpose: Generate visual summaries of the taxonomy-filtered feature table
#          and representative sequences, and construct a rooted phylogenetic
#          tree for downstream phylogenetic diversity analyses.
#
# Inputs:
#   - Taxonomy-filtered feature table (table.phylumfiltered.qza)
#   - Taxonomy-filtered representative sequences
#     (rep-seqs.phylumfiltered.qza)
#   - Sample metadata (usc-sediment-metadata.tsv)
#
# Outputs:
#   - Feature-table summary visualization
#     (table.phylum.organelles.filtered.qzv)
#   - Representative-sequence visualization
#     (rep-seqs.phylum.organelles.filtered.qzv)
#   - Aligned representative sequences (aligned-rep-seqs-clean.qza)
#   - Masked sequence alignment (masked-aligned-rep-seqs-clean.qza)
#   - Unrooted phylogenetic tree (unrooted-tree-clean.qza)
#   - Rooted phylogenetic tree (rooted-tree-clean.qza)
#
# Requirements: QIIME 2 (qiime2-amplicon-2024.10)
#
# "Clean" in the output filenames indicates that plant organelles
# and ASVs without a phylum-level assignment were removed during the preceding taxonomic filtering step.
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


OutputFolder=/mnt/scratch/wisnoskilab/jrp889/Output
AnalysisFolder=/mnt/scratch/wisnoskilab/jrp889/Analysis

set -euo pipefail

export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Generate a visual and tabular summary of the filtered feature table.
# Sample IDs in the metadata file should correspond to the sample IDs
# in the feature table.
qiime feature-table summarize \
    --i-table "$OutputFolder/table.phylumfiltered.qza" \
    --o-visualization "$OutputFolder/table.phylum.organelles.filtered.qzv" \
    --m-sample-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv"

# Generate a visualization of the filtered representative sequences.
qiime feature-table tabulate-seqs \
    --i-data "$OutputFolder/rep-seqs.phylumfiltered.qza" \
    --o-visualization "$OutputFolder/rep-seqs.phylum.organelles.filtered.qzv"

# Align the representative sequences with MAFFT, mask the alignment, construct an unrooted tree with FastTree, and midpoint-root the tree.
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences "$OutputFolder/rep-seqs.phylumfiltered.qza" \
    --o-alignment "$OutputFolder/aligned-rep-seqs-clean.qza" \
    --o-masked-alignment "$OutputFolder/masked-aligned-rep-seqs-clean.qza" \
    --o-tree "$OutputFolder/unrooted-tree-clean.qza" \
    --o-rooted-tree "$OutputFolder/rooted-tree-clean.qza"

# The rooted tree is required for downstream phylogenetic diversity analyses, including Faith's phylogenetic diversity.