#Comparing within siteS/forest-type; Bantayan and Olango at rarefaction 80

#use asv_perm which was previously defined during the RDA analysis preparation
#Convert asv table from counts to relative abundance
asv_rel <- decostand(asv_perm, method = "total")

# Check
rowSums(asv_rel)

#Subset the Bantayan sites or natural forests

natural_meta <- meta_perm[meta_perm$`forest-type` == "natural", ] #metadata
natural_asv <- asv_rel[rownames(natural_meta), ] #asv table

bray_natural <- vegdist(natural_asv, method = "bray") #bray curtis distance 

#PERMANOVA

permanova_natural_site <- adonis2(
  bray_natural ~ site,
  data = natural_meta,
  permutations = 999
)

adonis_natural_site

#PERMDISP
disp_natural_site <- betadisper(
  bray_natural,
  natural_meta$site
)

anova(disp_natural_site)
#not significant; differences not due to dispersion 

#Pairwise PERMANOVA
pairwise_natural_site <- pairwise_permanova(
  distance_matrix = bray_natural,
  group = natural_meta$site,
  permutations = 999,
  p_adjust_method = "BH"
)

pairwise_natural_site[order(pairwise_natural_site$p_adj), ]


pcoa_natural <- cmdscale(bray_natural, eig = TRUE, k = 2)

pcoa_scores <- as.data.frame(pcoa_natural$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$site <- natural_meta$site
pcoa_scores$sample_id <- rownames(natural_meta)

library(ggplot2)

ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = site)) +
  geom_point(size = 3) +
  stat_summary(aes(group = site), fun = mean, geom = "point", size = 5, shape = 4) +
  theme_classic() +
  labs(
    title = "Bray-Curtis PCoA of natural forest sediment microbes",
    x = "PCoA1",
    y = "PCoA2"
  )



###############################################
# within Olango ##############################
##############################################

planted_meta <- meta_perm[meta_perm$`forest-type` == "planted", ]
planted_asv <- asv_rel[rownames(planted_meta), ]

bray_planted <- vegdist(planted_asv, method = "bray")

adonis_planted_site <- adonis2(
  bray_planted ~ site,
  data = planted_meta,
  permutations = 999
)
 
adonis_planted_site #not significant

disp_planted_site <- betadisper(bray_planted, planted_meta$site)
anova(disp_planted_site) #dispersion not significant


#### PCA plots
# Bray-Curtis using all samples
bray_all <- vegdist(asv_rel, method = "bray")

# PCoA
pcoa_all <- wcmdscale(bray_all, eig = TRUE, k = 2, add = "lingoes")

# Scores
pcoa_all_scores <- as.data.frame(pcoa_all$points)
colnames(pcoa_all_scores) <- c("PCoA1", "PCoA2")
pcoa_all_scores$sample_id <- rownames(pcoa_all_scores)

# Add metadata
pcoa_all_scores$site <- meta_perm[rownames(pcoa_all_scores), "site"]
pcoa_all_scores$forest_type <- meta_perm[rownames(pcoa_all_scores), "forest-type"]

# Percent variation explained
eig_pos <- pcoa_all$eig[pcoa_all$eig > 0]
percent_var <- round(100 * eig_pos / sum(eig_pos), 1)

# Plot
pcoa_facet_plot <- ggplot(
  pcoa_all_scores,
  aes(x = PCoA1, y = PCoA2, color = site)
) +
  geom_point(size = 3, alpha = 0.9) +
  stat_summary(
    aes(group = site),
    fun = mean,
    geom = "point",
    size = 5,
    shape = 4,
    stroke = 1.3
  ) +
  facet_wrap(~ forest_type) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Bray-Curtis PCoA by forest type",
    x = paste0("PCoA1 (", percent_var[1], "%)"),
    y = paste0("PCoA2 (", percent_var[2], "%)"),
    color = "Site"
  )

pcoa_facet_plot

########################### UniFrac Distances #######################
######### Within site comparisons ##################################
#Unweighted UniFrac
#Prepare the variables 

meta_uw <- meta[labels(unifrac_dist_uw), , drop = FALSE]

all(labels(unifrac_dist_uw) == rownames(meta_uw))
table(meta_uw$site)
table(meta_uw$`forest-type`)

###########
# Natural/Bantayan sites
###########

natural_samples_uw <- rownames(meta_uw)[meta_uw$`forest-type` == "natural"]

unifrac_uw_natural <- as.dist(
  as.matrix(unifrac_dist_uw)[natural_samples_uw, natural_samples_uw]
)

natural_meta_uw <- meta_uw[natural_samples_uw, , drop = FALSE]

permanova_uw_natural_site <- adonis2(
  unifrac_uw_natural ~ site,
  data = natural_meta_uw,
  permutations = 999
)

permanova_uw_natural_site

# PERMDISP

disp_uw_natural_site <- betadisper(
  unifrac_uw_natural,
  natural_meta_uw$site,
  bias.adjust = TRUE
)

anova(disp_uw_natural_site) #not significant


#pairwise permanova for natural/bantayan
pairwise_uw_natural <- pairwise_permanova(
  distance_matrix = unifrac_uw_natural,
  group = natural_meta_uw$site,
  permutations = 999,
  p_adjust_method = "BH"
)

pairwise_uw_natural[order(pairwise_uw_natural$p_adj), ] #under powered

### Planted sites ########
planted_samples_uw <- rownames(meta_uw)[meta_uw$`forest-type` == "planted"]

unifrac_uw_planted <- as.dist(
  as.matrix(unifrac_dist_uw)[planted_samples_uw, planted_samples_uw]
)

planted_meta_uw <- meta_uw[planted_samples_uw, , drop = FALSE]

permanova_uw_planted_site <- adonis2(
  unifrac_uw_planted ~ site,
  data = planted_meta_uw,
  permutations = 999
)

permanova_uw_planted_site #not significant

#############################
## Weighted ############
# Weighted UniFrac
# Prepare variables as above

meta_w <- meta[labels(unifrac_dist_w), , drop = FALSE]

all(labels(unifrac_dist_w) == rownames(meta_w))
table(meta_w$site)
table(meta_w$`forest-type`)

######
# Natural/Bantayan
#####

natural_samples_w <- rownames(meta_w)[meta_w$`forest-type` == "natural"]

unifrac_w_natural <- as.dist(
  as.matrix(unifrac_dist_w)[natural_samples_w, natural_samples_w]
)

natural_meta_w <- meta_w[natural_samples_w, , drop = FALSE]

permanova_w_natural_site <- adonis2(
  unifrac_w_natural ~ site,
  data = natural_meta_w,
  permutations = 999
)

permanova_w_natural_site

#### PERMDISP ####

disp_w_natural_site <- betadisper(
  unifrac_w_natural,
  natural_meta_w$site,
  bias.adjust = TRUE
)

anova(disp_w_natural_site) #not significant

# Pairwise PERMANOVA
pairwise_w_natural <- pairwise_permanova(
  distance_matrix = unifrac_w_natural,
  group = natural_meta_w$site,
  permutations = 999,
  p_adjust_method = "BH"
)

pairwise_w_natural[order(pairwise_w_natural$p_adj), ]

###########################
## Planted/OIWS sites
########################
planted_samples_w <- rownames(meta_w)[meta_w$`forest-type` == "planted"] #subset planted samples


unifrac_w_planted <- as.dist(
  as.matrix(unifrac_dist_w)[planted_samples_w, planted_samples_w]
)

planted_meta_w <- meta_w[planted_samples_w, , drop = FALSE]

permanova_w_planted_site <- adonis2(
  unifrac_w_planted ~ site,
  data = planted_meta_w,
  permutations = 999
)

permanova_w_planted_site #not significant


