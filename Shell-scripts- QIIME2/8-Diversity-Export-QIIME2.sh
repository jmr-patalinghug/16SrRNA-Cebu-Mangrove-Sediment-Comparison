#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task 10
#SBATCH --time=24:00:00
#SBATCH --job-name=Diversity-and-Export
#SBATCH --output=qiime_diversity_%j.out
#SBATCH --error=qiime_diversity_%j.err

# Notes about this script:
#
# Purpose: Rarefy the taxonomy-filtered feature table to 80 reads per sample;
#          calculate core alpha- and beta-diversity metrics; test differences
#          in diversity among selected metadata groups; and export selected
#          QIIME 2 artifacts for subsequent analyses in R.
#
# Inputs:
#   - Taxonomy-filtered feature table (table.phylumfiltered.qza)
#   - Rooted phylogenetic tree (rooted-tree-clean.qza)
#   - Sample metadata (usc-sediment-metadata.tsv)
#
# Outputs:
#   - Rarefied feature table at a sampling depth of 80
#   - Alpha-diversity vectors:
#       - Faith's phylogenetic diversity
#       - Observed features
#       - Shannon diversity
#       - Pielou's evenness
#   - Beta-diversity distance matrices and PCoA results:
#       - Unweighted UniFrac
#       - Weighted UniFrac
#       - Bray-Curtis
#       - Jaccard
#   - Emperor PCoA visualizations
#   - Alpha-diversity group-significance visualizations
#   - Beta-diversity group-significance visualizations
#   - Exported distance matrices, alpha-diversity vectors, and rarefied
#     feature table for analysis in R
#
# Requirements:
#   - QIIME 2 (qiime2-amplicon-2024.10)
#   - BIOM Format command-line tools
#
# The sampling depth of 80 reads per sample was selected based on the exploratory alpha-rarefaction analyses. 
# This depth retained all 15 samples while providing adequate representation of alpha diversity.
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


OutputFolder=/mnt/scratch/wisnoskilab/jrp889/Output
AnalysisFolder=/mnt/scratch/wisnoskilab/jrp889/Analysis
CoreMetrics80="$OutputFolder/core-metrics-results-80"

set -euo pipefail

export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Calculate phylogenetic and non-phylogenetic diversity metrics after
# rarefying each sample to 80 reads.
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny "$OutputFolder/rooted-tree-clean.qza" \
    --i-table "$OutputFolder/table.phylumfiltered.qza" \
    --p-sampling-depth 80 \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --output-dir "$CoreMetrics80"

# Test differences in Faith's phylogenetic diversity among groups represented
# by categorical variables in the metadata file.
qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CoreMetrics80/faith_pd_vector.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --o-visualization "$CoreMetrics80/faith-pd-group-significance.qzv"

# Test differences in Pielou's evenness among groups represented by
# categorical variables in the metadata file.
qiime diversity alpha-group-significance \
    --i-alpha-diversity "$CoreMetrics80/evenness_vector.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --o-visualization "$CoreMetrics80/evenness-group-significance.qzv"

# Test whether community composition based on unweighted UniFrac distances
# differs by forest type, island, and sampling site.
# The --p-pairwise option also calculates pairwise comparisons among groups.

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column forest-type \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/unweighted-unifrac-forest-type-significance.qzv"

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column island \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/unweighted-unifrac-island-significance.qzv"

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/unweighted_unifrac_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column site \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/unweighted-unifrac-site-significance.qzv"

# Test whether community composition based on Bray-Curtis dissimilarities
# differs by forest type, island, and sampling site.

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/bray_curtis_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column forest-type \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/bray-curtis-forest-type-significance.qzv"

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/bray_curtis_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column island \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/bray-curtis-island-significance.qzv"

qiime diversity beta-group-significance \
    --i-distance-matrix "$CoreMetrics80/bray_curtis_distance_matrix.qza" \
    --m-metadata-file "$AnalysisFolder/usc-sediment-metadata.tsv" \
    --m-metadata-column site \
    --p-pairwise \
    --o-visualization "$CoreMetrics80/bray-curtis-site-significance.qzv"

# Export the unweighted UniFrac distance matrix for analysis in R.
# The exported folder contains a file named distance-matrix.tsv.
qiime tools export \
    --input-path "$CoreMetrics80/unweighted_unifrac_distance_matrix.qza" \
    --output-path "$OutputFolder/exported-unweighted-unifrac-80"

# Export the weighted UniFrac distance matrix for analysis in R.
# The exported folder also contains a file named distance-matrix.tsv.
qiime tools export \
    --input-path "$CoreMetrics80/weighted_unifrac_distance_matrix.qza" \
    --output-path "$OutputFolder/exported-weighted-unifrac-80"

# Rename exported distance-matrix.tsv files after exporting or keep them in
# clearly named folders to distinguish the two distance metrics.

# Export the rarefied feature table. A QIIME 2 FeatureTable[Frequency]
# artifact is exported as a BIOM file.
qiime tools export \
    --input-path "$CoreMetrics80/rarefied_table.qza" \
    --output-path "$CoreMetrics80/exported-rarefied-table-80"

# Convert the exported BIOM feature table to tab-separated format for use in R.
biom convert \
    -i "$CoreMetrics80/exported-rarefied-table-80/feature-table.biom" \
    -o "$CoreMetrics80/exported-rarefied-table-80/rarefied-table-80.tsv" \
    --to-tsv

# Export alpha-diversity vectors for analysis in R.
qiime tools export \
    --input-path "$CoreMetrics80/shannon_vector.qza" \
    --output-path "$CoreMetrics80/exported-shannon-80"

qiime tools export \
    --input-path "$CoreMetrics80/observed_features_vector.qza" \
    --output-path "$CoreMetrics80/exported-observed-features-80"

qiime tools export \
    --input-path "$CoreMetrics80/evenness_vector.qza" \
    --output-path "$CoreMetrics80/exported-evenness-80"

qiime tools export \
    --input-path "$CoreMetrics80/faith_pd_vector.qza" \
    --output-path "$CoreMetrics80/exported-faith-pd-80"