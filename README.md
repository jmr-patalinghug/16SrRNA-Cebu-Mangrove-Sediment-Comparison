# 16SrRNA-Cebu-Mangrove-Sediment-Comparison
This repository contains the scripts and supporting files used to analyze sediment bacterial communities from natural and planted mangrove forests in Cebu, Philippines. The study used 16S rRNA gene amplicon sequencing to compare bacterial community composition, diversity, environmental associations, and predicted functional potential among mangrove sites.

This repository supports transparency and reproducibility of the analyses reported in the associated manuscript.

## Study overview

Sediment samples were collected from natural and planted mangrove forests at five sampling sites in Cebu, Philippines:

- Panitugan
- Madridejos
- Santa Fe
- Olango Island Wildlife Sanctuary Site 1 (OIWS-1)
- Olango Island Wildlife Sanctuary Site 2 (OIWS-2)

A total of 15 sediment samples were analyzed.

The main analyses included:

Sequence denoising and ASV inference
Taxonomic classification
Removal of chloroplast, mitochondrial, unassigned, and domain-only features
Alpha- and beta-diversity analyses
Phylogenetic tree construction
Taxonomic composition analysis
Environmental-variable analysis in R
Functional prediction using PICRUSt2

## Folder Descriptions
### Shell-scripts-QIIME2

Contains SLURM-compatible Bash scripts used for sequence processing, taxonomic classification, taxonomic filtering, phylogenetic tree construction, rarefaction analysis, diversity analysis, and data export using QIIME2.

The scripts are numbered according to execution order by the researcher.

### Shell-scripts-Picrust2

Contains the script used to run the complete PICRUSt2 pipeline and predict functional potential from the filtered ASV sequences and abundance table.

### R-scripts

Contains R scripts used for statistical analyses and visualization.

### Input-files

Contains smaller metadata files and exported tables used by the QIIME 2, PICRUSt2, and R workflows.
This repository does not include large raw sequence files or generated QIIME 2 artifacts.

