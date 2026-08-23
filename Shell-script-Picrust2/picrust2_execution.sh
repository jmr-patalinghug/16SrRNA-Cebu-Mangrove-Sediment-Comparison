#!/bin/bash
#SBATCH --job-name=picrust2
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=picrust2_%j.out
#SBATCH --error=picrust2_%j.err
# Notes about this script:
#
# Purpose: Predict the functional potential of the sediment bacterial
#          communities from 16S rRNA gene sequences using the complete PICRUSt2 pipeline.
#
# Inputs:
#   - Taxonomy-filtered representative sequences in FASTA format
#     (dna-sequences.fasta)
#   - Taxonomy-filtered ASV feature table in BIOM format
#     (feature-table.biom)
#
# Outputs:
#   - Predicted gene-family abundances
#   - Predicted metabolic-pathway abundances
#   - Predicted pathway coverage
#   - Sequence-placement and hidden-state prediction results
#   - Quality-control and intermediate PICRUSt2 files
#
# Requirements:
#   - PICRUSt2 conda environment (install properly in cluster)
#   - Representative-sequence FASTA file
#   - ASV abundance table in BIOM format
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


# Activate PICRUSt2 environment
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate picrust2

mkdir -p /mnt/scratch/wisnoskilab/jrp889/picrust2/picrust2_out

picrust2_pipeline.py -s /mnt/scratch/wisnoskilab/jrp889/Output/exported_repseqs-clean/dna-sequences.fasta \
 -i /mnt/scratch/wisnoskilab/jrp889/Output/exported_table_phylumfiltered/feature-table.biom \
 -o /mnt/scratch/wisnoskilab/jrp889/picrust2/picrust2_out \
 -p 1 \
 --verbose

 ## Example provided in picrust2
 #picrust2_pipeline.py -s chemerin_16S_arc/seqs_arc.fna -i chemerin_16S_arc/table.biom -o picrust2_out_pipeline -p 1 --verbose