###############################################
#Generating Cleaner Taxa Barplots from QIIME2##
##############################################

#Load libraries
library(qiime2R)
library(tidyverse)
library(ggh4x)
library(colorspace)


# Import QIIME2 files
feature_table <- read_qza("table.phylumfiltered.qza")$data
taxonomy <- read_qza("taxonomy.qza")$data
metadata <- read_tsv("usc-sediment-metadata.tsv")

# Parse taxonomy
tax_split <- parse_taxonomy(taxonomy)

# Convert counts to relative abundance
phylum_table <- feature_table %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  left_join(
    tax_split %>%
      rownames_to_column("Feature.ID"),
    by = "Feature.ID"
  )

# Aggregate at phylum level
phylum_abundance <- phylum_table %>%
  select(-Feature.ID, -Kingdom, -Class, -Order,
         -Family, -Genus, -Species) %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum))

# Convert to long format
phylum_long <- phylum_abundance %>%
  pivot_longer(
    cols = -Phylum,
    names_to = "SampleID",
    values_to = "Abundance"
  )

# Relative abundance
phylum_data <- phylum_long %>%
  group_by(SampleID) %>%
  mutate(
    RelAbundance = Abundance / sum(Abundance) * 100
  ) %>%
  ungroup()

# Join metadata
metadata <- metadata %>%
  rename(SampleID = `sample-id`)

phylum_data <- left_join(
  phylum_data,
  metadata,
  by = "SampleID"
)

# Natural vs Planted grouping
phylum_data$ForestGroup <- ifelse(
  phylum_data$site %in% c(
    "Madredijos",
    "Panitugan",
    "Sta-Fe"
  ),
  "Natural",
  "Planted"
)

# Site order
phylum_data$site <- factor(
  phylum_data$site,
  levels = c(
    "Madredijos",
    "Panitugan",
    "Sta-Fe",
    "OIWS-1",
    "OIWS-2"
  )
)

# Cleaner sample labels
phylum_data$SampleLabel <- stringr::str_remove(
  phylum_data$SampleID,
  "_G_S[0-9]+$"
)

# Re-normalize (prevents rounding issues)
phylum_data <- phylum_data %>%
  group_by(SampleID) %>%
  mutate(
    RelAbundance = RelAbundance / sum(RelAbundance) * 100
  ) %>%
  ungroup()




######################
###### Plotting#######
######################

# Clean spelling and order phyla by total abundance

phylum_data_plot <- phylum_data %>%
  mutate(
    site = recode(site, "Madredijos" = "Madridejos"),
    
    ForestGroup = case_when(
      site %in% c("Madridejos", "Panitugan", "Sta-Fe") ~ "Natural",
      site %in% c("OIWS-1", "OIWS-2") ~ "Planted",
      TRUE ~ ForestGroup
    )
  )

# Order phyla by total relative abundance
phylum_order <- phylum_data_plot %>%
  group_by(Phylum) %>%
  summarise(total_abund = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  pull(Phylum)

phylum_data_plot <- phylum_data_plot %>%
  mutate(
    site = factor(site, levels = c("Madridejos", "Panitugan", "Sta-Fe", "OIWS-1", "OIWS-2")),
    ForestGroup = factor(ForestGroup, levels = c("Natural", "Planted")),
    Phylum = factor(Phylum, levels = phylum_order)
  )

# Colorblind-friendly palette for all phyla
n_phyla <- length(levels(phylum_data_plot$Phylum))

phyla_palette <- qualitative_hcl(
  n = n_phyla,
  palette = "Dark 3"
)

names(phyla_palette) <- levels(phylum_data_plot$Phylum)

taxa_barplot_all <- ggplot(
  phylum_data_plot,
  aes(
    x = SampleLabel,
    y = RelAbundance,
    fill = Phylum
  )
) +
  geom_col(
    width = 0.88,
    color = scales::alpha("grey75", 0.30),
    linewidth = 0.05
  ) +
  ggh4x::facet_nested(
    ~ ForestGroup + site,
    scales = "free_x",
    space = "free_x"
  ) +
  coord_cartesian(
    ylim = c(0, 100)
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = phyla_palette,
    labels = function(x) str_wrap(x, width = 24)
  ) +
  labs(
    x = "Sample ID",
    y = "Relative Abundance (%)",
    fill = "Phylum"
  ) +
  guides(
    fill = guide_legend(
      ncol = 3,
      byrow = TRUE,
      keyheight = unit(0.35, "cm"),
      keywidth = unit(0.35, "cm")
    )
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 9
    ),
    axis.title = element_text(size = 12),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width = unit(0.32, "cm"),
    legend.spacing.y = unit(0.03, "cm"),
    legend.box.spacing = unit(0.25, "cm")
  )

taxa_barplot_all

#######################################
### Plot of taxa (phyla) detected by island ###
#######################################

phylum_data_plot <- phylum_data_plot %>%
  mutate(
    Island = case_when(
      site %in% c("Madridejos", "Panitugan", "Sta-Fe") ~ "Bantayan",
      site %in% c("OIWS-1", "OIWS-2") ~ "Olango",
      TRUE ~ NA_character_
    )
  )


#choose taxonomic level
tax_level <- "Phylum"

#collapse counts by sample and taxon
taxa_data <- phylum_data_plot %>%
  filter(!is.na(Island)) %>%
  mutate(Taxon = .data[[tax_level]]) %>%
  filter(!is.na(Taxon), Taxon != "") %>%
  group_by(SampleID, SampleLabel, Island, site, Taxon) %>%
  summarise(
    Count = sum(Abundance),
    .groups = "drop"
  )

#select taxa that were present in at least two samples from each island

taxa_presence_strict <- taxa_data %>%
  group_by(Taxon, Island) %>%
  summarise(
    n_samples_present = n_distinct(SampleID[Count > 0]),
    total_count = sum(Count),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    Taxon,
    Island = c("Bantayan", "Olango"),
    fill = list(
      n_samples_present = 0,
      total_count = 0
    )
  ) %>%
  group_by(Taxon) %>%
  summarise(
    bantayan_samples = n_samples_present[Island == "Bantayan"],
    olango_samples = n_samples_present[Island == "Olango"],
    bantayan_count = total_count[Island == "Bantayan"],
    olango_count = total_count[Island == "Olango"],
    .groups = "drop"
  ) %>%
  mutate(
    unique_to = case_when(
      bantayan_samples >= 2 & olango_samples == 0 ~ "Bantayan",
      olango_samples >= 2 & bantayan_samples == 0 ~ "Olango",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(unique_to))

#Keep the island-specific taxa only
unique_taxa_data <- taxa_data %>%
  semi_join(taxa_presence_strict, by = "Taxon") %>%
  left_join(
    taxa_presence_strict %>% select(Taxon, unique_to),
    by = "Taxon"
  )

#calculate relative abundance 
sample_totals <- phylum_data_plot %>%
  group_by(SampleID) %>%
  summarise(
    sample_total = sum(Abundance),
    .groups = "drop"
  )

unique_taxa_rel <- unique_taxa_data %>%
  left_join(sample_totals, by = "SampleID") %>%
  mutate(
    RelAbundance = Count / sample_total * 100
  )

#plotting
ggplot(unique_taxa_rel, aes(x = SampleLabel, y = RelAbundance, fill = Taxon)) +
  geom_col(width = 0.8, color = "gray30", linewidth = 0.1) +
  facet_grid(~ Island, scales = "free_x", space = "free_x") +
  labs(
    x = "Sample",
    y = "Relative abundance (%)",
    fill = tax_level,
    title = paste("Strict island-specific", tax_level, "taxa")
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )



################
### By Family ##
################
# Make sure metadata sample ID is renamed only once
if ("sample-id" %in% names(metadata)) {
  metadata <- metadata %>%
    rename(SampleID = `sample-id`)
}

# Rebuild table from original QIIME2 feature table and taxonomy
taxa_table <- feature_table %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  left_join(
    tax_split %>%
      rownames_to_column("Feature.ID"),
    by = "Feature.ID"
  )


#collapse at Family level
tax_level <- "Family"

family_abundance <- taxa_table %>%
  mutate(Taxon = .data[[tax_level]]) %>%
  filter(!is.na(Taxon), Taxon != "") %>%
  select(Feature.ID, Taxon, where(is.numeric)) %>%
  group_by(Taxon) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


#convert to long format
family_long <- family_abundance %>%
  pivot_longer(
    cols = -Taxon,
    names_to = "SampleID",
    values_to = "Abundance"
  )

#Add metadata
family_data <- family_long %>%
  group_by(SampleID) %>%
  mutate(RelAbundance = Abundance / sum(Abundance) * 100) %>%
  ungroup() %>%
  left_join(metadata, by = "SampleID") %>%
  mutate(
    site = as.character(site),
    site = recode(site, "Madredijos" = "Madridejos"),
    site = factor(
      site,
      levels = c("Madridejos", "Panitugan", "Sta-Fe", "OIWS-1", "OIWS-2")
    ),
    Island = case_when(
      site %in% c("Madridejos", "Panitugan", "Sta-Fe") ~ "Bantayan",
      site %in% c("OIWS-1", "OIWS-2") ~ "Olango",
      TRUE ~ NA_character_
    ),
    ForestGroup = case_when(
      site %in% c("Madridejos", "Panitugan", "Sta-Fe") ~ "Natural",
      site %in% c("OIWS-1", "OIWS-2") ~ "Planted",
      TRUE ~ NA_character_
    ),
    SampleLabel = stringr::str_remove(SampleID, "_G_S[0-9]+$")
  )

#run selection for Families present in at least two samples per island
taxa_data <- family_data %>%
  filter(!is.na(Island)) %>%
  group_by(SampleID, SampleLabel, Island, site, Taxon) %>%
  summarise(
    Count = sum(Abundance),
    .groups = "drop"
  )

taxa_presence_strict <- taxa_data %>%
  group_by(Taxon, Island) %>%
  summarise(
    n_samples_present = n_distinct(SampleID[Count > 0]),
    total_count = sum(Count),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    Taxon,
    Island = c("Bantayan", "Olango"),
    fill = list(
      n_samples_present = 0,
      total_count = 0
    )
  ) %>%
  group_by(Taxon) %>%
  summarise(
    bantayan_samples = n_samples_present[Island == "Bantayan"],
    olango_samples = n_samples_present[Island == "Olango"],
    bantayan_count = total_count[Island == "Bantayan"],
    olango_count = total_count[Island == "Olango"],
    .groups = "drop"
  ) %>%
  mutate(
    unique_to = case_when(
      bantayan_samples >= 2 & olango_samples == 0 ~ "Bantayan",
      olango_samples >= 2 & bantayan_samples == 0 ~ "Olango",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(unique_to))


#checking the results
taxa_presence_strict %>%
  count(unique_to)

taxa_presence_strict


#Convert relative abundance

unique_taxa_rel <- taxa_data %>%
  semi_join(taxa_presence_strict, by = "Taxon") %>%
  left_join(
    taxa_presence_strict %>% select(Taxon, unique_to),
    by = "Taxon"
  ) %>%
  left_join(
    family_data %>%
      group_by(SampleID) %>%
      summarise(sample_total = sum(Abundance), .groups = "drop"),
    by = "SampleID"
  ) %>%
  mutate(RelAbundance = Count / sample_total * 100)


#checking
taxa_presence_strict %>%
  arrange(unique_to, desc(bantayan_samples + olango_samples))

library(writexl)

write_xlsx(taxa_presence_strict, "taxa_presence_strict.xlsx")
##########
#plotting#
##########

# Create a more distinct color palette
family_levels <- sort(unique(unique_taxa_rel$Taxon))

family_palette <- colorspace::qualitative_hcl(
  n = length(family_levels),
  palette = "Dark 3"
)

names(family_palette) <- family_levels



sample_order <- c(
  "BN1A", "BN1B", "BN1C",
  "BN2A", "BN2B", "BN2C",
  "BN3A", "BN3B", "BN3C",
  "OL1A", "OL1B", "OL1C",
  "OL2A", "OL2B", "OL2C"
)

plot_scaffold <- family_data %>%
  distinct(SampleLabel, Island) %>%
  mutate(
    SampleLabel = factor(SampleLabel, levels = sample_order),
    unique_to = Island,
    RelAbundance = 0
  )

unique_taxa_rel_plot <- unique_taxa_rel %>%
  filter(Island == unique_to) %>%
  mutate(
    SampleLabel = factor(SampleLabel, levels = sample_order)
  )


#plot
ggplot(unique_taxa_rel_plot, aes(x = SampleLabel, y = RelAbundance, fill = Taxon)) +
  geom_blank(
    data = plot_scaffold,
    aes(x = SampleLabel, y = RelAbundance),
    inherit.aes = FALSE
  ) +
  geom_col(width = 0.8, color = "gray20", linewidth = 0.15) +
  facet_grid(~ unique_to, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = family_palette) +
  labs(
    x = "Sample",
    y = "Relative abundance (%)",
    fill = "Family",
    title = "Bacterial families unique to each island"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(drop = TRUE) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
###################################################################################
##### Comparing sediment bacterial phyla across natural and planted mangroves #####
###################################################################################

# Make sure missing phylum-sample combinations are treated as 0
sample_meta <- phylum_data %>%
  distinct(SampleLabel, ForestGroup, site)

all_phyla <- unique(phylum_data$Phylum)

phylum_complete <- sample_meta %>%
  crossing(Phylum = all_phyla) %>%
  left_join(
    phylum_data %>%
      select(SampleLabel, Phylum, RelAbundance),
    by = c("SampleLabel", "Phylum")
  ) %>%
  mutate(RelAbundance = replace_na(RelAbundance, 0))


taxa_stats <- phylum_complete %>%
  group_by(Phylum) %>%
  summarise(
    mean_abundance = mean(RelAbundance),
    median_natural = median(RelAbundance[ForestGroup == "Natural"]),
    median_planted = median(RelAbundance[ForestGroup == "Planted"]),
    p_value = wilcox.test(RelAbundance ~ ForestGroup, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05
  ) %>%
  arrange(p_adj)

taxa_stats

#save stats
write_xlsx(taxa_stats, "compo_stats.xlsx")

significant_phyla <- taxa_stats %>%
  filter(p_adj < 0.05)

significant_phyla


taxa_stats_gt_1 <- phylum_complete %>%
  group_by(Phylum) %>%
  filter(mean(RelAbundance) > 1) %>%
  summarise(
    mean_abundance = mean(RelAbundance),
    median_natural = median(RelAbundance[ForestGroup == "Natural"]),
    median_planted = median(RelAbundance[ForestGroup == "Planted"]),
    p_value = wilcox.test(RelAbundance ~ ForestGroup, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significant = p_adj < 0.05
  ) %>%
  arrange(p_adj)

taxa_stats_gt_1
