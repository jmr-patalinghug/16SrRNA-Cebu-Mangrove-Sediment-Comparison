###################################
#Describing taxonomic composition##
###################################



#Load qiime2R package for file conversions without using the cluster
#qiime2R' version '0.99.6'
library(qiime2R)
library(tidyverse)

# Read QIIME2 artifacts
feature_table <- read_qza("table.phylumfiltered.qza")$data
taxonomy_raw <- read_qza("taxonomy.qza")$data

# Convert feature table to data frame
feature_df <- as.data.frame(feature_table) %>%
  rownames_to_column("FeatureID")

# Parse taxonomy into ranks
taxonomy_df <- taxonomy_raw %>%
  parse_taxonomy() %>%
  as.data.frame() %>%
  rownames_to_column("FeatureID")

# Check column names
colnames(taxonomy_df)

#############################
# Checking phylum abundance #
#############################

phylum_abundance <- feature_df %>%
  pivot_longer(
    cols = -FeatureID,
    names_to = "SampleID",
    values_to = "Count"
  ) %>%
  left_join(taxonomy_df, by = "FeatureID") %>%
  filter(!is.na(Phylum)) %>%
  filter(!str_detect(Phylum, regex("chloroplast|mitochondria|unassigned", ignore_case = TRUE))) %>%
  group_by(Phylum) %>%
  summarise(
    TotalReads = sum(Count),
    .groups = "drop"
  ) %>%
  mutate(
    RelativeAbundance = TotalReads / sum(TotalReads) * 100
  ) %>%
  arrange(desc(RelativeAbundance))

phylum_abundance


###################################
# Checking most and least abundant *
####################################

total_phyla <- n_distinct(phylum_abundance$Phylum)

most_abundant <- phylum_abundance %>%
  slice_max(RelativeAbundance, n = 10)

phyla_gt_1 <- phylum_abundance %>%
  filter(RelativeAbundance > 1) %>%
  arrange(desc(RelativeAbundance))

least_abundant <- phylum_abundance %>%
  filter(TotalReads > 0) %>%
  slice_min(RelativeAbundance, n = 10)

total_phyla
most_abundant
least_abundant
phyla_gt_1

