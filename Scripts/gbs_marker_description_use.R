## ## S2MET Registration
## Marker description and comparisons of the TP and VP using GBS markers
## 
## Author: Jeff Neyhart
## Last Modified: May 23, 2018
## 
## This script will examine the quality of the GBS data by creating a PCA
## and by looking at genomic coverage of markers
## 

# Load libraries and directories`
library(rrBLUP)

repo_dir <- getwd()
source(file.path(repo_dir, "source_use.R"))

## Load the Two-row CAP data
## The file path of this data should be changed for publication
load("path/to/data/2R_CAP_pheno_geno.RData")

# Create a color scheme
colors <- umn_palette(2, 7)[3:7]



## PCA of the relationship matrix of the TP and the VP

A2 <- A.mat(X = s2_imputed_mat, min.MAF = 0, max.missing = 1)
A2 <- A2[c(tp_geno, vp_geno), c(tp_geno, vp_geno)]

# Perform PCA
K2_pca <- svd(A2)

# Get the variance explained by each PC
K2_pca_stdev <- K2_pca$d %>% 
{. / sum(.)} %>% 
{. * 100} %>% 
  round(2) %>% 
  formatC(., digits = 2, format = "f") %>%
  str_c(str_c("PC", seq_along(.)), " (", ., "%)") %>%
  set_names(str_c("PC", seq_along(.)))

# Convert to DF and annotate
K2_pca_df <- K2_pca$u %>% 
  as.data.frame() %>% 
  `row.names<-`(., row.names(A2)) %>%
  `colnames<-`(., str_c("PC", seq(ncol(.)))) %>%
  rownames_to_column("line_name") %>% 
  left_join(., entry_list, by = c("line_name")) %>%
  select(line_name, program_abbr, population = class, parent, pedigree, PC1, PC2, PC3) %>%
  mutate(parent = if_else(parent == "TRUE", "Parent", ""))

# Create a color scheme
colors2 <- c(colors, umn_palette(2, 1))

## Show the BCAP population structure versus S2TP
g_s2_popstr <- K2_pca_df %>% 
  mutate(program = factor(program_abbr, levels = c(sort(setdiff(unique(program_abbr), "UMN")), "UMN"))) %>%
  ggplot(aes(x = PC1, y = PC2, size = parent, shape = program)) + 
  geom_point() +
  xlab(K2_pca_stdev[1]) +
  ylab(K2_pca_stdev[2]) +
  ylim(rev(range(K2_pca_df$PC2))) + 
  scale_size_manual(values = c(0.5,2), name = NULL, guide = FALSE) +
  scale_shape_discrete(name = "Breeding Program", guide = guide_legend(order = 1, ncol = 2)) +
  theme_acs() +
  theme(legend.position = c(0.2, 0.85), legend.key.height = unit(x = 0.75, units = "lines"))

# Save this
ggsave(filename = "s2_geno_comparison.jpg", plot = g_s2_popstr, path = fig_dir,
       height = 3.5, width = 3.5, dpi = 1000)












