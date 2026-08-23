#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=20
#SBATCH --job-name=Import_Seq

# Notes about this script:
#
# Purpose: Import raw paired-end sequence files into QIIME 2.
# Input: Raw FASTQ sequence files and a QIIME 2 manifest file.
# Output: QIIME 2 sequence artifact (.qza).
# Requirements: QIIME 2 (qiime2-amplicon-2024.10)
#
# File paths below are specific to the Mississippi State University Department of Biology computing cluster.
# Modify these file paths accordingly when running the script on a different computing system.


# Activate the QIIME 2 conda environment.
source /mnt/scratch/wisnoskilab/shared/bioinformatics/miniconda4/miniconda4/bin/activate
conda activate qiime2-amplicon-2024.10

# Import paired-end FASTQ sequences with Phred33 quality scores.
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path /mnt/scratch/wisnoskilab/jrp889/Analysis/usc-manifest.tsv \
  --output-path /mnt/scratch/wisnoskilab/jrp889/Output/usc-seq.qza \
  --input-format PairedEndFastqManifestPhred33V2

# View the .qzv through QIIME2 view (https://view.qiime2.org/). One can see the quality scores of the reads and some statistics about the reads. 
# Looking at this file is important to determine whether you want to truncate some reads later on for denoising.