#Picrust2 analysis

library(tidyverse)

#Load the picrust2 output describing pathway abundances 
pathways <- read_tsv("path_abun_unstrat_descrip.tsv")

# Check columns
colnames(pathways)
head(pathways)

#Renaming column names
#Usually first column is pathway ID and one column is description
id_cols <- c(colnames(pathways)[1], "description")
id_cols <- id_cols[id_cols %in% colnames(pathways)]

sample_cols <- setdiff(colnames(pathways), id_cols)

#Checking the top 25 pathways 
top_pathways <- pathways %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample_id",
    values_to = "pathway_abundance"
  ) %>%
  group_by(across(all_of(id_cols))) %>%
  summarize(
    mean_abundance = mean(pathway_abundance, na.rm = TRUE),
    total_abundance = sum(pathway_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) #arranging in descending order 

head(top_pathways, 25)

#save the top 25 pathways
write_csv(top_pathways,"D:/Jenny's Files/Research Related/NRCP Project/Paper 2/PiCrust2 output/top_predicted_pathways.csv")

#save all the pathways as .csv for easier viewing
write_csv(pathways,"D:/Jenny's Files/Research Related/NRCP Project/Paper 2/PiCrust2 output/all_pathways.csv")



########################################
##### Plotting the top 20 pathways ##### 
########################################

top20 <- top_pathways %>%
  slice_head(n = 20)

ggplot(
  top20,
  aes(
    x = reorder(description, mean_abundance),
    y = mean_abundance
  )
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    x = "Predicted MetaCyc pathway",
    y = "Mean predicted pathway abundance",
    title = "Top predicted metabolic pathways"
  )


################################################################
##### Plotting pathways with potential biotech application ##### 
################################################################

#Selecting pathways with potential application to biotech (based on literature search)
#You may replace or add pathway IDs here
target_pathways <- c(
  "AEROBACTINSYN-PWY",  #aerobactin biosynthesis / siderophore
  "PWY-6808",          #dTDP-alpha-D-forosamine biosynthesis
  "PWY-6942",           #dTDP-D-desosamine biosynthesis
  "DTDPRHAMSYN-PWY"    #production of rhamnolipids
)

# Convert from wide to long format and assign island
pathways_long <- pathways %>%
  pivot_longer(
    cols = -c(pathway, description),
    names_to = "sample_id",
    values_to = "abundance"
  ) %>%
  mutate(
    island = case_when(
      str_starts(sample_id, "BN") ~ "Bantayan (Natural)",
      str_starts(sample_id, "OL") ~ "Olango (Planted)",
      TRUE ~ NA_character_
    )
  )

#Convert to relative abundance within each sample
pathways_rel <- pathways_long %>%
  group_by(sample_id) %>%
  mutate(relative_abundance = abundance / sum(abundance, na.rm = TRUE) * 100) %>%
  ungroup()

# Keep only pathways of interest
target_data <- pathways_rel %>%
  filter(pathway %in% target_pathways)

# Summarize by island and pathway
target_summary <- target_data %>%
  group_by(island, pathway, description) %>%
  summarize(
    mean_relative_abundance = mean(relative_abundance, na.rm = TRUE),
    se_relative_abundance = sd(relative_abundance, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

target_summary


#Plotting relative abundance
ggplot(
  target_data,
  aes(
    x = island,
    y = relative_abundance,
    color = island
  )
) +
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8, height = 0) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.4,
    color = "black"
  ) +
  facet_wrap(~ description, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Island",
    y = "Predicted pathway relative abundance (%)",
    title = "Pathways with biotechnological applications by Island",
    color = "Island"
  ) +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )





#############################
target_data <- target_data %>%
  mutate(
    description = str_replace_all(description, "&alpha;", "α"),
    description = str_replace_all(description, "&beta;", "β"),
    island = factor(
      island,
      levels = c("Bantayan (Natural)", "Olango (Planted)")
    )
  )




### plot ###

  ggplot(
    target_data,
    aes(
      x = island,
      y = relative_abundance,
      color = island
    )
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2.5,
    alpha = 0.8
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.18,
    color = "black",
    linewidth = 0.5
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 2.5,
    color = "black"
  ) +
  facet_wrap(~ description, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  labs(
    x = "Island",
    y = "PICRUSt2-predicted pathway relative abundance (%)",
    title = "Predicted biotechnologically relevant pathways by island",
    color = "Island"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1),
    panel.grid.minor = element_blank()
  )
  