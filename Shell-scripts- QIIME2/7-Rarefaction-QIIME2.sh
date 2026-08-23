#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task 10
#SBATCH --time=24:00:00
#SBATCH --job-name=Rarefaction-plots
#SBATCH --output=qiime_rarefaction_%j.out
#SBATCH --error=qiime_rarefaction_%j.err

# Notes about this script:
#
# Purpose: Generate exploratory alpha-rarefaction curves at several maximum
#          sequencing depths. These visualizations were used to evaluate
#          diversity coverage and sample retention and to select an appropriate
#          sampling depth for downstream diversity analyses.
#
#          The final sampling depth selected for this study was 80 reads per
#          sample. That depth was applied during the subsequent diversity
#          analysis and is not applied directly by this script.
#
# Inputs:
#   - Taxonomy-filtered feature table (table.phylumfiltered.qza)
#   - Rooted phylogenetic tree (rooted-tree-clean.qza)
#   - Sample metadata (usc-sediment-metadata.tsv)
#
# Outputs:
#   - Alpha-rarefaction visualizations at different depths (.qzv)
#
# Requirements: QIIME 2 (qiime2-amplicon-2024.10)
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


OutPutFolder=/mnt/scratch/wisnoskilab/jrp889/Output/
AnalysisFolder=/mnt/scratch/wisnoskilab/jrp889/Analysis/

set -euo pipefail
export TMPDIR="$SLURM_SUBMIT_DIR/qiime_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10


# p-sampling-depth: change per experiment; look across your samples for feature counts or Shannon diversity values -- choose a number as high as possible but lower than the lowest sample
qiime diversity alpha-rarefaction --i-table $OutPutFolder/table.phylumfiltered.qza --i-phylogeny $OutPutFolder/rooted-tree-clean.qza --p-max-depth 50 --m-metadata-file $AnalysisFolder/usc-sediment-metadata.tsv --o-visualization $OutPutFolder/alpha-rarefaction50-clean.qzv
qiime diversity alpha-rarefaction --i-table $OutPutFolder/table.phylumfiltered.qza  --i-phylogeny $OutPutFolder/rooted-tree-clean.qza --p-max-depth 150 --m-metadata-file $AnalysisFolder/usc-sediment-metadata.tsv --o-visualization $OutPutFolder/alpha-rarefaction150-clean.qzv
qiime diversity alpha-rarefaction --i-table $OutPutFolder/table.phylumfiltered.qza  --i-phylogeny $OutPutFolder/rooted-tree-clean.qza --p-max-depth 500 --m-metadata-file $AnalysisFolder/usc-sediment-metadata.tsv  --o-visualization $OutPutFolder/alpha-rarefaction500-clean.qzv
qiime diversity alpha-rarefaction --i-table $OutPutFolder/table.phylumfiltered.qza  --i-phylogeny $OutPutFolder/rooted-tree-clean.qza --p-max-depth 1000 --m-metadata-file $AnalysisFolder/usc-sediment-metadata.tsv --o-visualization $OutPutFolder/alpha-rarefaction1000-clean.qzv
