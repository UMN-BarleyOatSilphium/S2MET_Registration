## ## S2MET Registration
## Prediction accuracy with increasing training population size
## Prediction accuracy across pairs of environments
## 
## Author: Jeff Neyhart
## Last Modified: May 28, 2018
## 
## This script will provide an example of the impact of increasing training population
## size on the accuracy of genomewide predictions. We will use the genotype mean 
## for grain yield in Minnesota and surrounding locations
## 

# Load libraries and directories`
library(rrBLUP)

repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load BLUEs
load(file.path(proj_dir, "Results/all_trait_variance_component_heritability.RData"))
load(file.path(proj_dir, "Results/phenotype_analysis_data.RData"))

## Select the locations of interest
locs <- c("STP", "CRM")
traits <- c("HeadingDate", "GrainProtein")
trait_replace <- set_names(c("Heading Date", "Grain Protein"), traits)



# # Filter the dataset and calculate genotype means
# pheno_tomodel <- S2_MET_BLUEs_use %>% 
#   filter(location %in% locs, trait %in% traits) %>%
#   mutate_at(vars(line_name, environment), as.factor)
# 
# control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
# 
# geno_means <- pheno_tomodel %>%
#   group_by(trait) %>%
#   do({
#     
#     df <- droplevels(.)
#     
#     fit <- lmer(value ~ -1 + line_name + (1|environment) + (1|line_name:environment),
#                 data = df, control = control, weights = df$std.error^2) 
# 
#     # Grab the genotype BLUEs
#     tidy(fit) %>%
#       filter(group == "fixed") %>% 
#       mutate(term = levels(df$line_name)) %>% 
#       select(line_name = term, value = estimate) %>%
#       mutate(line_name = as.factor(line_name))
#     
#   }) %>% ungroup()

geno_means <- trait_main_BLUEs %>% 
  ungroup() %>% 
  rename(line_name = grp)


## Predictions of all traits
all_data_predictions <- geno_means %>%
  group_by(trait) %>%
  do(accuracy = {
    df <- droplevels(.)
    
    train <- model.frame(value ~ line_name, df, subset = line_name %in% tp_geno)
    y <- model.response(train)
    Z <- s2_imputed_mat[train$line_name,,drop = FALSE]
    fit <- mixed.solve(y = y, Z = Z)
    
    val <- subset(df, line_name %in% vp_geno)
    preds <- c(s2_imputed_mat[val$line_name,] %*% fit$u)
    
    cor(val$value, preds)
    
  })
    

# Print
all_data_predictions %>%
  ungroup() %>%
  unnest() %>%
  arrange(accuracy) %>%
  mutate(accuracy = round(accuracy, 2)) %>%
  write_csv(x = ., path = file.path(fig_dir, "all_data_prediction_accuracy.csv"))

geno_means <- filter(geno_means, trait %in% traits)


## Define TP sizes
tp_sizes <- seq(5, 175, by = 5)
samples <- 250

# Filter out the lines without marker data
geno_means_pred <- geno_means %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>% 
  droplevels()


# Set seed for reproducibility
set.seed(1153)
# Iterate over tp_sizes and traits
pred_acc_tp_size <- crossing(trait = unique(geno_means_pred$trait), size = tp_sizes) %>%
  group_by(trait, size) %>%
  do({
    
    tr <- .$trait
    sz <- .$size

    # Subset the means into tp and vp
    tp_means <- geno_means_pred %>%
      filter(trait == tr,
             line_name %in% tp_geno)
    vp_means <- geno_means_pred %>%
      filter(trait == tr,
             line_name %in% vp_geno) %>%
      mutate(line_name = as.character(line_name))

    # Generate random samples of the training set
    tp_samples <- replicate(n = samples, sample_n(tbl = tp_means, size = sz), simplify = FALSE)
  
    pred_out <- tp_samples %>%
      map_dbl(~{
        tp_i <- .
        
        y <- tp_i$value
        Z <- s2_imputed_mat[tp_i$line_name,,drop = FALSE]
        
        fit <- mixed.solve(y = y, Z = Z)
        preds <- c(s2_imputed_mat[vp_geno,] %*% fit$u)
        
        cor(vp_means$value, preds)
        
      })
  
  as.data.frame(cbind(iter = seq(samples), pred_out))

}) %>% ungroup()
  

## Summarize over iterations
pred_acc_tp_size_summ <- pred_acc_tp_size %>% 
  group_by(trait, size) %>% 
  summarize_at(vars(pred_out), funs(mean = mean, sd = sd)) %>%
  ungroup()

## What was the prediction accuracy using the full dataset
pred_acc_tp_size_summ %>% 
  filter(size == max(size))
  
# trait         size  mean            sd
# 1 GrainProtein   175 0.576 0.00000000120
# 2 HeadingDate    175 0.603 0.00000000424


## Plot points
g_pred_acc <- pred_acc_tp_size_summ %>%
  mutate(trait = str_replace_all(trait, trait_replace)) %>%
  ggplot(aes(x = size, y = mean, ymin = mean - sd, ymax = mean + sd)) + 
  geom_point(size = 1) +
  geom_line() +
  geom_errorbar(width = 5) +
  ylab("Prediction accuracy") +
  xlab("Training population size") +
  facet_grid(~ trait) +
  ylim(c(-0.05,1)) +
  theme_acs()

# Save
ggsave(filename = "pred_acc_tp_size.jpg", plot = g_pred_acc, path = fig_dir, height = 3, width = 5, dpi = 1000)



### Prediction accuracy across pairs of environments


# What traits to select
# traits <- c("PlantHeight", "TestWeight", "HeadingDate")
# traits <- c("PlantHeight", "TestWeight", "HeadingDate")


# Filter the dataset for common environments
pheno_tomodel <- S2_MET_BLUEs_use %>% 
  filter(trait %in% traits) %>% 
  group_by(environment) %>%
  filter(n_distinct(trait) == length(traits)) %>% 
  group_by()

# Filter again for environments in which the TP and VP were phenotyped
pheno_tomodel1 <- pheno_tomodel %>% 
  group_by(trait, environment) %>%
  filter(any(line_name %in% tp_geno) & any(line_name %in% vp_geno)) %>%
  ungroup()

# What are the common environments
common_env <- unique(pheno_tomodel1$environment)

# Pairs of environments
pair_env <- crossing(trait = traits, environment1 = common_env, environment2 = common_env)

## Iterate over pairs
pair_env_pred_acc <- pair_env %>%
  group_by(trait, environment1, environment2) %>%
  do(accuracy = {
    
    df <- .
    tp_env <- df$environment1
    vp_env <- df$environment2
    tr <- df$trait
    
    
    # Filter the tp set
    train_data <- pheno_tomodel1 %>%
      filter(environment == tp_env, trait == tr, line_name %in% tp_geno)
    # Filter the vp set
    vp_data <- pheno_tomodel1 %>%
      filter(environment == vp_env, trait == tr, line_name %in% vp_geno)
    
    vp_geno1 <- intersect(vp_geno, vp_data$line_name)
    
    ## Predict and correlate
    mf <- model.frame(value ~ line_name, train_data)
    y <- model.response(mf)
    Z <- s2_imputed_mat[mf$line_name,,drop = FALSE]
    
    fit <- mixed.solve(y = y, Z = Z)
    preds <- c(s2_imputed_mat[vp_geno1,] %*% fit$u)
    
    cor(vp_data$value, preds)
    
  })


## For each trait, determine the order of environments using a heatmap
env_orders <- pair_env_pred_acc %>% 
  unnest() %>%
  spread(environment2, accuracy) %>%
  as.data.frame() %>%
  split(.$trait) %>%
  map(~select(., -trait) %>% remove_rownames() %>% column_to_rownames("environment1")) %>%
  map(~{
    hm <- heatmap(as.matrix(.))
    row.names(.)[hm$rowInd]
  })


## Summarize
pair_env_pred_acc1 <- pair_env_pred_acc %>%
  unnest() %>% 
  mutate(within_env = environment1 == environment2)

pair_env_pred_acc1 %>% 
  group_by(trait, within_env) %>% 
  summarize_at(vars(accuracy), funs(min, max, mean, median))

# trait       within_env     min   max  mean median
# 1 GrainProtein FALSE      -0.361   0.513 0.204  0.225
# 2 GrainProtein TRUE       -0.00846 0.572 0.306  0.305
# 3 HeadingDate  FALSE      -0.0713  0.700 0.416  0.408
# 4 HeadingDate  TRUE        0.0303  0.686 0.448  0.472

# Plot
pair_env_pred_acc1 %>% 
  qplot(x = accuracy, fill = within_env, data = ., facets = ~trait, geom = "density") + 
  theme_acs()

## How many environments could predict a different environment as well as or better than itself?
pair_env_pred_acc1 %>% 
  group_by(trait, environment1) %>% 
  mutate(within_env_acc = accuracy[within_env], 
         better_acc = accuracy > within_env_acc) %>% 
  group_by(trait) %>% 
  summarize(p_better_acc = mean(better_acc))






pair_env_pred_acc1 %>% 
  group_by(trait) %>% 
  summarize_at(vars(accuracy), funs(min, max, mean, median))

# trait          min   max  mean median
# 1 GrainProtein -0.361  0.572 0.214  0.241
# 2 HeadingDate  -0.0713 0.700 0.419  0.411


greys <- grey.colors(n = 3, start = 0, end = 1)

# Plot
env_acc_list <- pair_env_pred_acc %>%
  unnest() %>%
  mutate(trait = str_replace_all(trait, trait_replace)) %>%
  split(.$trait) %>%
  map2(.x = ., .y = env_orders, ~mutate_at(.x, vars(contains("environment")), funs(factor(., levels = .y))) %>%
         ggplot(., aes(x = environment1, y = environment2, fill = accuracy)) + 
         geom_tile() +
         ylab("Validation environment") +
         xlab("Prediction environment") + 
         scale_fill_gradient2(low = greys[3], mid = greys[2], high = greys[1], name = "Prediction\naccuracy") + 
         facet_wrap(~ trait, scales = "free") + 
         theme_acs() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") )

# Combine
g_env_acc <- plot_grid(
  env_acc_list[[1]], 
  env_acc_list[[2]] + theme(axis.title.y = element_blank()),
  # env_acc_list[[3]] + theme(axis.title.y = element_blank()),
  ncol = length(traits), align = "hv")

# Save
ggsave(filename = "env_pair_pred_acc.jpg", plot = g_env_acc, path = fig_dir, height = 4, width = 6, dpi = 1000)


## Combine figures
g_combined <- plot_grid(g_pred_acc, g_env_acc, ncol = 1, labels = LETTERS[1:2], rel_heights = c(0.8, 1))

# Save
ggsave(filename = "prelim_pred_acc_combined.jpg", plot = g_combined, path = fig_dir, height = 6, width = 6, dpi = 1000)
















