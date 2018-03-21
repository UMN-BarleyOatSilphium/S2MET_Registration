## S2MET Registration
## Phenotypic and Genotypic Comparisons of S2TP with BCAP
## 
## Author: Jeff Neyhart
## Last Modified: March 18, 2018
## 
## This script will create some figures to show the relationship of the S2TP
## genotypic data (BOPA only) and phenotypic dataa (only data presumed available
## when creating the S2TP).
## 

# Load libraries and directories`
library(rrBLUP)

repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

## Load the Two-row CAP data
## The file path of this data should be changed for publication
load("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Phenotypic Data/2R_CAP/Base_Example_Phenos/2R_CAP_pheno_geno.RData")


## Phenotypic Comparisons
# Add the annotation of each line to the data.frame
# Add the number of lines in each population to the population annotation
bcap_pheno_ann_toplot <- bcap_pheno_fits %>% 
  left_join(., entry_list_relevant, by = c("line_name" = "name")) %>% 
  select(line_name, population = class, program, parent, line_name, trait, BLUE, BLUP) %>%
  group_by(population) %>% 
  mutate(population_ann = str_c(population, " (n = ", n_distinct(line_name), ")")) %>% 
  ungroup()

# Create a color palette
colors <- umn_palette(2, 4)[2:3]

# Vector to replace trait names
# Use multiple quotes to work with label_parsed
trait_replace <- unique(bcap_pheno_ann_toplot$trait) %>% 
  set_names(x = c("'Alpha Amylase (20°DU)'", "'Beta Glucan (%)'", "'Grain Protein (%)'", 
              "'Grain Yield '(kg~ha^-1)", "'Spot Blotch Severity (%)'"), nm = .)

# Show the BCAP data versus the S2TP data
g_pheno_density <- bcap_pheno_ann_toplot %>% 
  filter(trait != "AlphaAmylase") %>%
  mutate(trait = str_replace_all(trait, trait_replace)) %>%
  # ggplot(aes(x = BLUP, fill = population_ann)) + 
  ggplot(aes(x = BLUE, fill = population_ann)) +
  geom_density(alpha = 0.5) + 
  xlab("Trait Value") +
  ylab("Density") +
  scale_fill_manual(values = colors, name = NULL) + 
  facet_wrap(~ trait, ncol = 2, scales = "free", labeller = label_parsed) + 
  theme_bw() +
  theme(legend.position = c(0.85, 0.32),
        legend.background = element_rect(fill = "grey85"),
        panel.grid = element_blank())
        

# Save this
save_file <- file.path(fig_dir, "bcap_s2tp_pheno_comparison.jpg")
ggsave(filename = save_file, plot = g_pheno_density, height = 5, width = 6, dpi = 1000)


## Compare genotype marker data

# Perform PCA
K_pca <- prcomp(bcap_geno_impute$A)

# Get the variance explained by each PC
K_pca_stdev <- K_pca$sdev %>% 
  {. / sum(.)} %>% 
  {. * 100} %>% 
  round(2) %>% 
  str_c(str_c("PC", seq_along(.)), " (", ., "%)") %>%
  set_names(str_c("PC", seq_along(.)))

# Convert to DF and annotate
K_pca_df <- K_pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name") %>% 
  left_join(., entry_list_relevant, by = c("line_name" = "name")) %>% 
  select(line_name, program, population = class, parent, PC1, PC2, PC3) %>%
  mutate(parent = if_else(parent == "TRUE", "S2TP/Parent", "S2TP")) 
  

# Create a color scheme
colors <- umn_palette(2, 7)[3:7]

## Show the BCAP population structure versus S2TP
g_bcap_popstr <- K_pca_df %>% 
  ggplot(aes(x = PC1, y = PC2, size = parent)) + 
  geom_point(aes(color = program)) +
  geom_point(data = filter(K_pca_df, population == "S2TP"), color = "black") +
  xlab(K_pca_stdev[1]) +
  ylab(K_pca_stdev[2]) +
  scale_size_manual(values = c(1.5,3), name = "Parent") +
  scale_color_manual(values = colors, name = "Breeding\nProgram") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.70))

# Save this
save_file <- file.path(fig_dir, "bcap_s2tp_geno_comparison.jpg")
ggsave(filename = save_file, plot = g_bcap_popstr, height = 6, width = 7.5, dpi = 1000)




