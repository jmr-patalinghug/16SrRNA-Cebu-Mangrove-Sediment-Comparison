library(vegan)
library(glmnet)
library(ggplot2)

######### checking ############

#Loading the rarefied table to double check whether any sites got dropped or not
rarefied_80 <- read.table(
  "rarefied-table-80.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  skip = 1,
  check.names = FALSE,
  comment.char = ""
)
head(rarefied_80, n=3)

#transposing the table with t() function so that ASVS will be columns and sample IDs will be rows
rarefied_80_t <- t(rarefied_80) 
head(rarefied_80_t, n=3)

nrow(rarefied_80_t)
rownames(rarefied_80_t)

#checking the number of reads
#now all at 80; only p sampling depth that does not drop BN2C_G_S159 which affects the statistical analysis 
sort(rowSums(rarefied_80_t)) 



###########################
# Proceeding to analysis####
##########################


#----### UniFrac Distances###----#
########## Weighted ############

#loading the table
unifrac_mat_w <- read.table(
  "distance-matrix-weighted-80.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

unifrac_dist_w <- as.dist(unifrac_mat_w) #converting it to distance measures


# Read metadata
meta <- read.table(
  "usc-sediment-metadata.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

#match unifrac and metadata
meta_unifrac_w <- meta[rownames(unifrac_mat_w), ]

all(rownames(unifrac_mat_w) == rownames(meta_unifrac_w)) #check



######## PERMANOVA ############
permanova_weighted_unifrac <- adonis2(
  unifrac_dist_w ~ site,
  data = meta_unifrac_w,
  permutations = 999
)

permanova_weighted_unifrac 

##### PERMDISP ######
disp_unifrac_w_site <- betadisper(
  unifrac_dist_w,
  meta_unifrac_w$site
)

anova(disp_unifrac_w_site)
plot(disp_unifrac_w_site) #dispersion was not significant which means within site variation did not contribute to the differences across sites

#boxplots to check which groups look different
boxplot(disp_unifrac_w_site) #note: clean this with ggplot


### Pairwise PERMANOVA function ###
pairwise_permanova <- function(distance_matrix, group, permutations = 999, p_adjust_method = "BH") {
  
  group <- as.factor(group)
  group_levels <- levels(group)
  pairs <- combn(group_levels, 2, simplify = FALSE)
  
  results <- lapply(pairs, function(pair) {
    
    idx <- group %in% pair
    
    dist_sub <- as.dist(as.matrix(distance_matrix)[idx, idx])
    group_sub <- droplevels(group[idx])
    
    ad <- adonis2(
      dist_sub ~ group_sub,
      data = data.frame(group_sub = group_sub),
      permutations = permutations
    )
    
    data.frame(
      group1 = pair[1],
      group2 = pair[2],
      n_group1 = sum(group_sub == pair[1]),
      n_group2 = sum(group_sub == pair[2]),
      F = ad$F[1],
      R2 = ad$R2[1],
      p_value = ad$`Pr(>F)`[1]
    )
  })
  
  results <- do.call(rbind, results)
  results$p_adj <- p.adjust(results$p_value, method = p_adjust_method)
  
  results
}

#####

pairwise_Uni_w_site_80 <- pairwise_permanova(
  distance_matrix = unifrac_dist_w,
  group = meta_unifrac_w$site,
  permutations = 999,
  p_adjust_method = "BH"
)

pairwise_Uni_w_site_80[order(pairwise_Uni_w_site_80$p_adj), ]
#underpowered to detect significant differences among individual site pairs after multiple-testing correction, likely due to limited replication with three samples per site.



#################################
########## Unweighted ############
unifrac_mat_uw <- read.table(
  "distance-matrix-unweighted-80.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

unifrac_dist_uw <- as.dist(unifrac_mat_uw)

#match unifrac and metadata
meta_unifrac_uw <- meta[rownames(unifrac_mat_uw), ]

all(rownames(unifrac_mat_uw) == rownames(meta_unifrac_uw)) #check

######## PERMANOVA ############
permanova_unweighted_unifrac <- adonis2(
  unifrac_dist_uw ~ site,
  data = meta_unifrac_uw,
  permutations = 999
)

permanova_unweighted_unifrac 

##### PERMDISP ######
disp_unifrac_uw_site <- betadisper(
  unifrac_dist_uw,
  meta_unifrac_uw$site
)

anova(disp_unifrac_uw_site)
plot(disp_unifrac_uw_site) #dispersion was not significant which means within site variation did not contribute to the differences acrcoss sites

#boxplots to check which groups look different
boxplot(disp_unifrac_uw_site) #note: clean this with ggplot

#####

pairwise_Uni_uw_site_80 <- pairwise_permanova(
  distance_matrix = unifrac_dist_uw,
  group = meta_unifrac_uw$site,
  permutations = 999,
  p_adjust_method = "BH"
)

pairwise_Uni_uw_site_80[order(pairwise_Uni_uw_site_80$p_adj), ]


#####################################
################ BRAY-CURTIS ########

bray_80 <- read.table(
  "distance-matrix-bray-80.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

bray_dist_80 <- as.dist(bray_80)

#match Bray-curtis and metadata
meta_bray_80 <- meta[rownames(bray_80), , drop = FALSE]

all(labels(bray_dist_80) == rownames(meta_bray_80))
table(meta_bray_80$site)

######## PERMANOVA ############
permanova_bray <- adonis2(
  bray_dist_80 ~ site,
  data = meta_bray_80,
  permutations = 999
)

permanova_bray

#### PERMDISP ######
disp_bray_80 <- betadisper(
  bray_dist_80,
  meta_bray_80$site
)

anova(disp_bray_80)

