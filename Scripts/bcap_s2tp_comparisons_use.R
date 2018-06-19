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
source(file.path(repo_dir, "source_use.R"))

## Load the Two-row CAP data
## The file path of this data should be changed for publication
load("path/to/data/2R_CAP_pheno_geno.RData")

## Merge entry lists
entry_list_bcap <- entry_list %>%
  filter(name %in% row.names(bcap_geno_impute$A))


## Phenotypic Comparisons
# Add the annotation of each line to the data.frame
# Add the number of lines in each population to the population annotation
bcap_pheno_ann_toplot <- bcap_pheno_fits %>% 
  left_join(., entry_list, by = c("line_name" = "name")) %>% 
  select(line_name, population = class, program, parent, line_name, trait, BLUE, BLUP) %>%
  mutate(population = ifelse(population == "S2TP", "TP", population)) %>%
  group_by(population) %>% 
  mutate(population_ann = str_c(population, " (n = ", n_distinct(line_name), ")")) %>% 
  ungroup()


# Create a color palette
colors <- umn_palette(2, 4)[2:3]

# Vector to replace trait names
# Use multiple quotes to work with label_parsed
trait_replace <- c(GrainYield = "'Grain Yield '(kg~ha^-1)", GrainProtein = "'Grain Protein (%)'",
                   MaltBetaGlucan = "'MaltBetaGlucan (ppm)'", SpotBlotchSeverity = "'Spot Blotch Severity (%)'")

# Show the BCAP data versus the S2TP data
g_pheno_density <- bcap_pheno_ann_toplot %>% 
  mutate(trait = str_replace_all(trait, trait_replace),
         trait = factor(trait, levels = trait_replace)) %>%
  # ggplot(aes(x = BLUP, fill = population_ann)) + 
  ggplot(aes(x = BLUE, fill = population_ann)) +
  geom_density(alpha = 0.5) + 
  xlab("Trait Value") +
  ylab("Density") +
  # scale_fill_manual(values = colors, name = NULL) + 
  scale_fill_grey(name = NULL) + 
  facet_wrap(~ trait, ncol = 2, scales = "free", labeller = label_parsed) + 
  theme_acs() +
  theme(legend.position = c(0.13, 0.89),
        # legend.background = element_rect(fill = "grey85", alpha = 1),
        panel.grid = element_blank())
        

# Save this
save_file <- file.path(fig_dir, "bcap_s2tp_pheno_comparison.jpg")
ggsave(filename = save_file, plot = g_pheno_density, height = 4, width = 5, dpi = 1000)


## Compare genotype marker data

# Perform PCA
K_pca <- svd(bcap_geno_impute$A)

# Get the variance explained by each PC
K_pca_stdev <- K_pca$d %>% 
  {. / sum(.)} %>% 
  {. * 100} %>% 
  round(2) %>% 
  formatC(., digits = 2, format = "f") %>%
  str_c(str_c("PC", seq_along(.)), " (", ., "%)") %>%
  set_names(str_c("PC", seq_along(.)))

# Convert to DF and annotate
K_pca_df <- K_pca$u %>% 
  as.data.frame() %>% 
  `row.names<-`(., row.names(bcap_geno_impute$A)) %>%
  `colnames<-`(., str_c("PC", seq(ncol(.)))) %>%
  rownames_to_column("line_name") %>% 
  left_join(., entry_list_bcap, by = c("line_name" = "name")) %>% 
  select(line_name, program_abbr, population = class, parent, PC1, PC2, PC3) %>%
  mutate(parent = if_else(parent == "TRUE", "TP/Parent", "TP")) 
  

# Create a color scheme
colors <- umn_palette(2, 7)[3:7]

## Show the BCAP population structure versus S2TP
g_bcap_popstr <- K_pca_df %>% 
  mutate(population = ifelse(population == "S2TP", "TP", population)) %>%
  ggplot(aes(x = PC1, y = PC2, size = population, shape = program_abbr)) + 
  geom_point() +
  # geom_point(data = filter(K_pca_df, population == "S2TP"), color = "black") +
  xlab(K_pca_stdev[1]) +
  ylab(K_pca_stdev[2]) +
  scale_size_manual(values = c(1.5,4), name = "Population", guide = guide_legend(order = 2)) +
  # scale_color_manual(values = colors, name = "Breeding\nProgram", guide = guide_legend(order = 1)) +
  scale_shape_discrete(guide = guide_legend(order = 1, title = "Breeding\nProgram",)) +
  theme_acs() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.1, 0.75))

# Save this
ggsave(filename = "bcap_s2tp_geno_comparison.jpg", plot = g_bcap_popstr, path = fig_dir,
       height = 6, width = 7, dpi = 1000)








