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
source(file.path(repo_dir, "source.R"))

## Load the Two-row CAP data
## The file path of this data should be changed for publication
load("C:/Users/Jeff/GoogleDrive/BarleyLab/Breeding/PhenotypicData/2R_CAP/Base_Example_Phenos/2R_CAP_pheno_geno.RData")

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


# ## Plot with lines for pedigree
# # First create a pedigree DF for the VP
# m2_pca_df <- K2_pca_df %>% 
#   filter(population == "S2C1R") %>% 
#   separate(pedigree, into = c("parent1", "parent2"), sep = "/") %>%
#   gather(par_number, parent, parent1, parent2) %>%
#   left_join(., select(K2_pca_df, line_name, PC1:PC3), by = c("parent" = "line_name")) %>% 
#   rename_at(vars(contains(".x")), ~str_replace(string = ., pattern = ".x", replacement = "")) %>%
#   rename_at(vars(contains(".y")), ~str_replace(string = ., pattern = ".y", replacement = "_parent"))
# 
# 
# g_s2_popstr_ped <- K2_pca_df %>% 
#   mutate(program = factor(program, levels = c(sort(setdiff(unique(program), "M2")), "M2"))) %>%
#   ggplot(aes(x = PC1, y = PC2)) + 
#   geom_segment(data = m2_pca_df, aes(x = PC1, xend = PC1_parent, y = PC2, yend = PC2_parent), lwd = 0.1) + 
#   geom_point(aes(color = program, size = parent)) +
#   xlab(K2_pca_stdev[1]) +
#   ylab(K2_pca_stdev[2]) +
#   ylim(rev(range(K2_pca_df$PC2))) + 
#   scale_size_manual(values = c(1,2.5), name = NULL, guide = FALSE) +
#   scale_color_manual(values = colors2, name = "Breeding Program", guide = guide_legend(order = 1, ncol = 2)) +
#   theme_acs() +
#   theme(legend.position = c(0.2, 0.85), legend.key.height = unit(x = 0.75, units = "lines"))
# 
# # Save this
# save_file <- file.path(fig_dir, "s2_geno_comparison_ped.jpg")
# ggsave(filename = save_file, plot = g_s2_popstr_ped, height = 5, width = 5, dpi = 1000)
# 
# 
# ## Breakdown of the VP by family
# (vp_family_count <- entry_list %>% 
#   filter(class == "S2C1R") %>% 
#   group_by(pedigree) %>% 
#   summarize(n = n()) %>% 
#   arrange(desc(n)) %>%
#   separate(pedigree, c("parent1", "parent2"), sep = "/"))
# 
# # Number of parents used
# vp_family_count %>% select(parent1, parent2) %>% unlist() %>% n_distinct()
# 




# ## Calculate and plot the minor allele frequency of markers between the TP and VP
# M_tp <- s2_discrete_mat[tp_geno,]
# M_vp <- s2_discrete_mat[vp_geno,]
# M_list <- list(tp = M_tp, vp = M_vp)
# 
# ## Calculate minor allele frequency for each
# maf_list <- M_list %>% 
#   map(~colMeans(. + 1, na.rm = TRUE) / 2) %>% 
#   map(~pmin(., 1 - .))



## Plot the density of markers across the genome

## Plot the chromosomes
g_snp_density <- barley_lengths %>%
  mutate(chrom = as.factor(chrom)) %>%
  ggplot(aes(y = chrom, yend = chrom, x = 0, xend = length / 1e+6)) +
  geom_segment() +
  geom_point(data = mutate(snp_info, chrom = as.factor(chrom)), aes(x = pos / 1e+6, y = chrom),
             inherit.aes = FALSE) +
  xlab("Position (Mbp)") +
  ylab("Chromosome") +
  # facet_grid(chrom ~ ., switch = "y") +
  theme_acs() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# Save this
ggsave(filename = "s2_gbs_marker_density.jpg", plot = g_snp_density, path = fig_dir,
       height = 3.5, width = 7, dpi = 1000)





## Compare marker density with gene density using the barley gene annotation
## information

## Load the barley annotation
ann_dir <- "C:/Users/Jeff/GoogleDrive/BarleyLab/Projects/Genomics/Annotation"
load(file.path(ann_dir, "barley_genomic_ranges.RData"))


## First plot gene annotation density across the genome
barley_genes <- barley_grange_list$genes %>%
  subset(seqnames != "chrUn")
# Set the lengths of the seqnames to the lengths of the chromosomes
barley_seqlengths <- setNames(object = barley_lengths$length, 
                              nm = str_c("chr", barley_lengths$chrom, "H"))


seqlevels(barley_genes) <- head(seqlevels(barley_genes), -1)
seqlengths(barley_genes) <- barley_seqlengths

# Set a window size and step size
# Here I define both as 10 Mbp
window_size <- step_size <- 10000000

## Split the annotation by chromosome
barley_genes_split <- split(x = barley_genes, f = seqnames(barley_genes)) %>% 
  as.list()

## Split the seqlengths by chromosome
barley_seqlengths_split <- list(names(barley_genes_split), seqlengths(barley_genes)) %>% 
  pmap(~tileGenome(seqlengths = setNames(.y, .x), tilewidth = window_size))

## For each chromosome, subset by overlap with the annotation list
barley_genes_perchrom <- list(barley_seqlengths_split, barley_genes_split) %>%
  pmap_df(~{
    
    # Find the overlaps between the genes and the sequences
    overlaps <- findOverlaps(subject = .y, query = .x) %>%
      as.data.frame()
    
    barley_seqlengths_df <- .x %>% 
      as.data.frame() %>% 
      select(group, seqnames:width)
    
    ## Merge and return
    left_join(overlaps, barley_seqlengths_df, by = c("queryHits" = "group")) %>%
      group_by(queryHits) %>% 
      do({ 
        df <- .
        df %>% 
          distinct(queryHits, seqnames, start, end, width) %>% 
          mutate(n_genes = n_distinct(.y[df$subjectHits, ]$gene_id))
      }) %>%
      ungroup()
  })


## Use the snp info data.frame to create a GRanges object
snp_info_granges <- snp_info %>% 
  mutate(chrom = str_c("chr", chrom, "H")) %>% 
  select(marker = `rs#`, names(.)) %>%
  makeGRangesFromDataFrame(df = ., keep.extra.columns = TRUE, start.field = "pos", end.field = "pos")

# Reset the seqlengths
seqlengths(snp_info_granges) <- barley_seqlengths

## Split the markers by chromosome
snp_info_granges_split <- split(x = snp_info_granges, f = seqnames(snp_info_granges)) %>% 
  as.list()

## For each chromosome, subset by overlap with the marker list
snp_info_per_window <- list(barley_seqlengths_split, snp_info_granges_split) %>%
  pmap_df(~{
    
    # Find the overlaps between the genes and the sequences
    overlaps <- findOverlaps(subject = .y, query = .x) %>%
      as.data.frame()
    
    barley_seqlengths_df <- .x %>% 
      as.data.frame() %>% 
      select(group, seqnames:width)
    
    ## Merge and return
    left_join(overlaps, barley_seqlengths_df, by = c("queryHits" = "group")) %>%
      group_by(queryHits) %>% 
      do({ 
        df <- .
        df %>% 
          distinct(queryHits, seqnames, start, end, width) %>% 
          mutate(n_markers = n_distinct(.y[df$subjectHits, ]$marker))
      }) %>%
      ungroup()
  })

# Combine the marker and gene info
barley_markers_genes <- left_join(barley_genes_perchrom, snp_info_per_window, 
                                  by = c("queryHits", "seqnames", "start", "end", "width")) %>% 
  mutate(n_markers = if_else(is.na(n_markers), 0L, n_markers),
         markers_per_gene = n_markers / n_genes)


## Plot
g_gene_marker_dens <- barley_markers_genes %>% 
  gather(measure, value, n_genes:markers_per_gene) %>%
  filter(measure != "markers_per_gene") %>%
  mutate(pos = (start + end) / 2,
         chrom = parse_number(seqnames),
         measure = str_replace_all(measure, c("n_genes" = "Number of Genes", 
                                              "n_markers" = "Number of Markers"))) %>%
  ggplot(aes(x = pos / 1e6, y = value, col = measure)) + 
  geom_smooth(method = "loess", span = 0.25) +
  # geom_line() +
  xlab("Position (Mbp)") +
  scale_color_discrete(name = NULL) +
  # facet_grid(measure ~ chrom, scales = "free_x", switch = "both", space = "free_x") +
  facet_grid(~ chrom, scales = "free_x", switch = "both", space = "free_x") +
  theme_acs()  %+replace%
  theme_manhattan() +
  theme(axis.title = element_blank(),
        legend.position = "bottom")


# Save this
ggsave(filename = "s2_gbs_marker_density_genes.jpg", plot = g_gene_marker_dens, path = fig_dir,
       height = 3.5, width = 7, dpi = 1000)
















