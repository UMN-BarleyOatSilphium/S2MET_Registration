## Project: S2MET Registration
## 
## Phenotypic data adjustment and analysis
## Author: Jeff Neyhart
## Last updated: May 29, 2018
## 

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

## Load the optimx package
library(optimx)


## Concatenate the pedigree information
# Load the pedigree data
ped_file <- "C:/Users/Jeff/GoogleDrive/BarleyLab/Breeding/PopulationInformation/Pedigree/UMN_S2_Pedigrees.xlsx"
tp_ped <- read_excel(path = ped_file, sheet = "TP")
vp_ped <- read_excel(path = ped_file, sheet = "S2C1R")

# Combine
s2met_pedigree <- bind_rows(
  select(tp_ped, Line, Pedigree),
  select(vp_ped, Line, Pedigree)
)

# Print to csv
write_csv(x = s2met_pedigree, path = "S2MET_pedigree.csv")



### Adjust Phenotypic Values ###

## Create a data.frame of data to model
# Set a dummy variable for check or test
data_to_model <- S2_MET_tidy_use %>%
  mutate(line_name = as.character(line_name),
         line = ifelse(!line_name %in% checks, line_name, "00check"),
         check = ifelse(line_name %in% checks, line_name, "00line")) %>%
  mutate_at(vars(line, check), as.factor)

# Designate two models, depending on the environment and the experimental design
forms <- formulas(~ value,
                  form1 = ~ -1 + line + check,
                  form2 = ~ -1 + line + check + (1|blk))

# Fit the models
stage_one <- data_to_model %>%
  group_by(trial, trait) %>%
  do({
    
    # Create a separate df object - this greatly improves the speed of the model fitting
    df <- droplevels(.)
    
    # ## Edit the contrasts
    # contrasts(df$line) <- contr.sum(levels(df$line)) %>% `colnames<-`(., head(levels(df$line), -1))
    # contrasts(df$check) <- contr.sum(levels(df$check)) %>% `colnames<-`(., head(levels(df$check), -1))
    # 
    
    # Extract the function and the formula, then fit
    if (is.na(df$blk[1])) {
      f <- forms$form1
      
      fit <- lm(formula = f, data = df)
      
      ## Test
      # fit1 <- lm(value ~ line + check, df)
      
      # Tidy
      fit_tidy <- tidy(fit) %>% 
        mutate(group = "fixed")
      
    } else {
      f <- forms$form2
      f_mean <- str_replace(string = f, pattern = "-1 \\+", "") %>%
        tail(-1) %>% 
        str_c(collapse = " ~")
      
      fit <- lmer(formula = f, data = df)
      
      fit_tidy <- tidy(fit)
      
    }
    
    # Get the levels of the checks
    check_levels <- levels(df$check)
    
    fit_tidy1 <- fit_tidy %>% 
      filter(group == "fixed") %>% 
      mutate(term = if_else(term == "line00check", tail(check_levels, 1), 
                            str_replace_all(term, pattern = "line|check", replacement = "")), 
             estimate = if_else(term %in% head(check_levels, -1), estimate + estimate[1], estimate)) %>% 
      select(line = term, estimate, std.error)
    
    # Find the harmonic mean of the number of replicates
    n_r <- table(df$line_name) %>%
      harm_mean()
    
    # Return the coefficients and the variance of the adjusted means
    data_frame(BLUE = list(fit_tidy1), n_r = n_r) })


# Add metadata to the BLUEs
S2_MET_BLUEs_use <- stage_one %>% 
  ungroup() %>% 
  unnest(BLUE) %>% 
  filter(!line %in% checks) %>% 
  select(trial, trait, line_name = line, estimate, std.error) %>% 
  left_join(., distinct(S2_MET_tidy, trial, environment, location, year), by = "trial") %>%
  select(trial, location, year, environment, trait, line_name, value = estimate, std.error)


## Remove potential outlier environment-trait combinations
S2_MET_BLUEs_use <- S2_MET_BLUEs_use %>% 
  filter(!(trait == "GrainYield" & environment == "HNY16"))





## Variance component table

# Fit mixed models (with random effects for everthing except location and year)
# to determine variance components.
# Use a likelihood ratio test to assess the significance of each of the variance
# comonents
ranef_terms <- c("(1|line_name)", "(1|environment)", "(1|line_name:environment)")
fixef_terms <- c()

## Fit each model to each trait
# We will need to drop models that are not relevant to certain traits (i.e. if 
# those traits were not evaluated in more than one year/location/etc.)


# Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", 
                       calc.derivs = FALSE, optimizer = "optimx",
                       optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))



# Iterate over traits
trait_varcomp <- S2_MET_BLUEs_use %>%
  # Convert to factors
  mutate_at(vars(environment, location, year, line_name), as.factor) %>%
  group_by(trait) %>%
  do({
    
    df <- droplevels(.)
    
    # Print the trait
    print(unique(df$trait))
    
    # Determine the number of locations and years
    # n_loc <- n_distinct(df$location)
    # n_year <- n_distinct(df$year)
    n_env <- n_distinct(df$environment)
    
    # Determine the harmonic mean of the number of locations, years, and reps
    plot_table <- xtabs(~ line_name + environment, df)
    
    ## Harmonic means
    # Locations
    # harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
    #   ifelse(. > 1, 1, .) %>%
    #   rowSums() %>% 
    #   harm_mean()
    # 
    # # Year
    # harm_year <- apply(X = plot_table, MARGIN = c(1,3), sum) %>% 
    #   ifelse(. > 1, 1, .) %>%
    #   rowSums() %>% 
    #   harm_mean()
    # 
    harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      harm_mean()
    
    
    # Remove terms based on the number of locations or years
    # if (all(n_loc == 1, n_year == 1)) {
    #   # Remove any terms with location and year
    #   remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "location|year")]
    #   removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "location|year")]
    #   
    # } else if (n_loc == 1) {
    #   # Remove any term with location
    #   remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "location")]
    #   removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "location")]
    #   
    # } else if (n_year == 1) {
    #   # Remove any term with year
    #   remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "year")]
    #   removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "year")]
    #   
    # } else {
    #   # Remove no terms
    #   remove_fixef_terms <- removed_ranef_terms <- character()
    #   
    # }
    
    if (n_env == 1) {
      # Remove environment
      removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "environment")]
      
    } else {
      # Remove no terms
      remove_fixef_terms <- removed_ranef_terms <- character()
      
    }
    
    
    # Remove any terms
    ranef_terms_use <- setdiff(ranef_terms, removed_ranef_terms)
    fixef_terms_use <- setdiff(fixef_terms, remove_fixef_terms)
    
    # Build the models
    # Iterate over the random effects and drop one
    ranef_terms_drop <- ranef_terms_use %>% 
      set_names(., .) %>% 
      map(~setdiff(ranef_terms_use, .))
    
    # Add the full model
    # Build the sets of random effects for the reduced models
    ranef_model_terms <- c(full = list(ranef_terms_use), ranef_terms_drop)
    # Create the formula strings for these random effect terms
    ranef_model_form <- map(ranef_model_terms, str_c, collapse = " + ")
    
    # Create the formula strings for fixed effects
    fixef_model_form <- str_c(fixef_terms_use, collapse = " + ")
    
    # Combine the model formulas
    forms <- ranef_model_form %>% 
      map(~str_c("value ~ ", str_c(fixef_model_form, ., sep = " + "))) %>% 
      map(as.formula)
    
    # Extract the weights
    wts <- df$std.error^2
    
    # Fit the models
    fits <- fit_with(data = df, .f = lmer, .formulas = forms, control = control,
                     weights = wts)
    
    # Extract the variance components from the full model
    var_comp <- VarCorr(fits$full) %>% 
      as.data.frame()
    
    # Extract the likelihood ratios
    fits_loglik <- map(fits, logLik) %>% map(as.numeric)
    
    # Perform likelihood ratio tests
    fits_lrt <- fits_loglik[-1] %>% 
      map(~(-2) * (. - fits_loglik$full)) %>%
      map_dbl(~pchisq(q = ., df = 1, lower.tail = FALSE) / 2)
    
    # Convert to data.frame
    fits_lrt_df <- fits_lrt %>% 
      set_names(x = ., nm = str_replace_all(string = names(.), pattern = "\\(1\\||\\)", "")) %>% 
      data.frame(grp = names(.), pvalue = ., row.names = NULL, stringsAsFactors = FALSE)
    
    # Combine with the variance components
    var_comp1 <- left_join(x = var_comp, y = fits_lrt_df, by = "grp") %>%
      dplyr::select(term = grp, variance = vcov, pvalue)
    
    # # Calculate heritability
    # h2 <- herit(object = fits$full, n_l = harm_loc, n_y = harm_year, n_r = harm_rep,
    #             exp = "line_name / (line_name + (line_name:year / n_y) + (line_name:location / n_l) + 
    #             (line_name:year:location / (n_l * n_y)) + (Residual / (n_l * n_y * n_r)))")
    
    h2 <- herit(object = fits$full, n_e = harm_env, n_r = harm_rep,
                exp = "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))")
    
    # Return the heritability
    var_comp2 <- var_comp1 %>% mutate(h2 = h2$heritability)
    
    # Return the variance components and the full model fit
    data_frame(var_comp = list(var_comp2), full_fit = fits[1])
    
  })





# ## Try the mohring and piepho 2009 method of heritability
# ## The estimates using this approach are not very different from the ad hoc heritability
# ## method outlined above
# df <- S2_MET_BLUEs_use %>%
#   filter(trait == "AlphaAmylase")
# 
# fit <- lmer(value ~ (1|line_name) + (1|environment) + (1|line_name:environment), 
#             data = df, weights = df$std.error, control = control)
# 
# blups <- ranef(fit, condVar = T)$line_name
# vblup <- 2 * mean(attr(blups, "postVar"))
# 
# varG <- c(VarCorr(fit)[["line_name"]])
# 
# (h2_cullis <- 1 - vblup / (2 * varG))
# 
# # Determine the harmonic mean of the number of locations, years, and reps
# plot_table <- xtabs(~ line_name + environment, df)
# 
# harm_env <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
#   ifelse(. > 1, 1, .) %>%
#   rowSums() %>% 
#   harm_mean()
# 
# # Reps
# harm_rep <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
#   harm_mean()
# 
# h2 <- herit(object = fit, n_e = harm_env, n_r = harm_rep,
#       exp = "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))")
# 
# c(trad = h2$heritability, cullis = h2_cullis)









# Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", 
                       calc.derivs = FALSE, optimizer = "optimx",
                       optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE))

## Fit genotype and environment as fixed effects to measure the range in values
trait_main_BLUEs <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    df <- .
    
    # Model formula
    form <- value ~ line_name + environment + (1|line_name:environment)
    
    # Get the weights
    wts <- df$std.error^2
    
    # Model fitting
    fit <- lmer(formula = form, data = df, control = control, weights = wts)
    
    ## Get the fixed effects
    blues <- fixef(fit) %>% 
      data.frame(grp = names(.), value = ., row.names = NULL, stringsAsFactors = FALSE)
    
    # Scale to intercept and return
    blues %>% 
      mutate(grp = str_replace_all(grp, "line_name|environment", ""), 
             value = ifelse(grp == "(Intercept)", value, value + value[1]), 
             grp = ifelse(grp == "(Intercept)", levels(as.factor(df$line_name))[1], grp))
    
  })


## Calculate the min and max values across environments, the whole S2MET,
## and the two smaller populations
trait_main_BLUEs_grp <- trait_main_BLUEs %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(string = grp, "^[0-9]"), "line_name", "environment"))

# First environments
(trait_main_BLUEs_env <- trait_main_BLUEs_grp %>% 
    filter(type == "environment") %>% 
    group_by(trait) %>% 
    summarize_at(vars(value), funs(min, max, mean)))

(trait_main_BLUEs_gen <- trait_main_BLUEs_grp %>%  
    filter(type == "line_name") %>% 
    mutate(population = ifelse(str_detect(grp, "^0"), "TP", "VP")) %>% 
    group_by(trait, population) %>%  
    summarize_at(vars(value), funs(min, max, mean)))


## Read-in the trait observation info
trait_info <- read_csv(file = file.path(fig_dir, "trait_information_table.csv"))



# Designate significance levels for the variance components
trait_varcomp_sig <- trait_varcomp %>%
  ungroup() %>% 
  unnest(var_comp) %>%
  select(term, trait, variance, pvalue) %>% 
  mutate(
    sig_notation = case_when(
      pvalue <= 0.001 ~ "***",
      pvalue <= 0.01 ~ "**",
      pvalue <= 0.05 ~ "*",
      TRUE ~ ""),
    variance = ifelse(variance < 1, signif(variance, 3), round(variance, 3)), # Choose rounding method based on size
    variance_notation = str_trim(str_c(variance, sig_notation, sep = " "))) %>%
  left_join(., t3_traits, by = c("trait" = "Nickname")) 

# Create the table for printing
trait_varcomp_sig_toprint <- trait_varcomp_sig %>% 
  select(term, Trait, variance_notation) %>% 
  mutate(term = str_to_title(str_replace_all(string = term, pattern = ":", replacement = " %*% ")),
         term = str_replace_all(string = term, "Line_name", "Genotype"),
         Trait = str_to_title(Trait),
         Trait = factor(Trait, levels = trait_info$Trait)) %>%
  spread(term, variance_notation) %>%
  arrange(Trait)


save_file <- file.path(fig_dir, "trait_variance_components.csv")
write_csv(x = trait_varcomp_sig_toprint, path = save_file)








## Summarize the mean, min, max, sd, and heritability of each trait
trait_summary <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>% 
  summarize_at(vars(value), funs(min, max, mean, sd)) %>%
  left_join(., distinct(unnest(trait_varcomp, var_comp), trait, h2))

# Format to print
trait_summary_toprint <- trait_summary %>% 
  left_join(., select(t3_traits, Trait, Nickname), by = c("trait" = "Nickname")) %>% 
  mutate(Trait = str_to_title(Trait)) %>% 
  select(-trait) %>% 
  gather(parameter, value, -Trait) %>% 
  mutate(value = ifelse(value < 0, 0, value),
         value = signif(value, 3)) %>% # Round the values
  spread(parameter, value) %>% 
  mutate(Trait = factor(Trait, levels = trait_info$Trait)) %>% 
  select(Trait, Min = min, Max = max, Mean = mean, S.D. = sd, H = h2) %>%
  arrange(Trait)

save_file <- file.path(fig_dir, "trait_value_summary.csv")
write_csv(x = trait_summary_toprint, path = save_file)








## Calculate genetic correlations among the genotype BLUEs
trait_geno_BLUE <- trait_main_BLUEs_grp %>%
  filter(type == "line_name") %>% 
  select(-type) %>%
  spread(trait, value)




# Create a data.frame for pairwise trait comparisons
trait_phencor <- distinct(trait_main_BLUEs, trait) %>% 
  ungroup() %>%
  crossing(., .) %>%
  mutate(data = pmap(list(trait, trait1), ~{trait_geno_BLUE[,c(.x, .y)]}),
         cor = map(data, ~bootstrap(x = .[,1], y = .[,2], fun = "cor", boot.reps = 1000,
                                    alpha = c(0.001, 0.01, 0.05))))


## Unnest and annotate based on significance (using the confidence interval)
trait_phencor_ann <- trait_phencor %>%
  unnest(cor) %>%
  mutate(significant = !(ci_lower <= 0 & 0 <= ci_upper),
         alpha = rep(c(0.001, 0.01, 0.05), length.out = nrow(.)),
         significant_ann = case_when(significant & alpha == 0.001 ~ "***",
                                     significant & alpha == 0.01 ~ "**",
                                     significant & alpha == 0.05 ~ "*",
                                     TRUE ~ "")) %>%
  group_by(trait, trait1) %>%
  filter(nchar(significant_ann) == max(nchar(significant_ann))) %>%
  # Take the first row of the combinations, if none are signficant
  slice(1) %>%
  ungroup()



# Convert to matrix
trait_phencor_mat <- trait_phencor_ann %>% 
  select(trait, trait1, base) %>% 
  spread(trait1, base) %>% 
  remove_rownames() %>% 
  as.data.frame() %>%
  column_to_rownames("trait") %>% 
  as.matrix()


## Plot the correlation
# Color scheme
color <- umn_palette(2, 4)

# Cluster the correlations using the heatmap function
# Abbreviate the trait names
trait_order <- unique(trait_phencor_ann$trait)[heatmap(trait_phencor_mat)$rowInd] %>%
  abbreviate(minlength = 2)


## Designate an observation as the upper triangle or the lower triangle
trait_phencor_ann_toplot <- trait_phencor_ann %>%
  mutate_at(vars(contains("trait")), ~str_replace_all(., trait_order) %>% 
              factor(., levels = trait_order)) %>%
  arrange(trait, trait1) %>%
  mutate(triangle = apply(X = select(., trait, trait1), MARGIN = 1, FUN = sort) %>% 
           t() %>% duplicated(),
         triangle = ifelse(triangle, "lower", "upper"),
         text = ifelse(significant_ann == "", "", str_c(round(base, 2), significant_ann)),
         text = ifelse(triangle == "lower", text, ""), # Set the lower-triangle text
         base = ifelse(triangle == "lower", NA, base)) # Set the lower-triangle values to 0



greys <- grey.colors(n = 3, start = 0.1, end = 0.9)

# Plot using ggplot
g_trait_phencor <- trait_phencor_ann_toplot %>%
  ggplot(aes(x = trait, y = trait1, fill = base, label = text)) + 
  geom_tile(color = "grey75") + 
  geom_text(size = 3) +
  # scale_fill_gradient2(low = color[3], mid = "white", high = color[1], 
  #                      name = "Phenotypic\nCorrelation\n") +
  scale_fill_gradient2(name = "Phenotypic\nCorrelation\n", low = greys[3], mid = greys[2], high = greys[1],
                       midpoint = 0, na.value = "white") +
  labs(caption = "*Significant at the 0.05 significance level
**Significant at the 0.01 significance level
***Significant at the 0.001 significance level
") + 
  theme_acs() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom")

# Save the graph
ggsave(filename = "trait_genetic_correlations.jpg", plot = g_trait_phencor, path = fig_dir,
       width = 7, height = 7, dpi = 1000)








## Convert the data.frame to a table to print
# First determine the repetitive trait pairs
# Abbreviate the trait names
trait_pairs_remove <- select(trait_phencor_ann, trait, trait1) %>% 
  apply(X = ., MARGIN = 1, FUN = sort) %>%
  t() %>%
  duplicated()

# Create that data.frame
trait_phencor_ann_toprint <- trait_phencor_ann %>% 
  # Need to abbreviate the trait names so the duplicate removal works
  arrange(desc(trait), desc(trait1)) %>%
  mutate(cor_ann = paste(round(cor, 3), significant_ann),
         cor_ann = str_trim(cor_ann)) %>% 
  select(trait, trait1, cor_ann) %>%
  filter(trait_pairs_remove) %>%
  complete(trait, trait1, fill = list(cor_ann = "")) %>%
  spread(trait1, cor_ann)


# Save the correlation matrix
save_file <- file.path(fig_dir, "trait_genetic_correlations.csv")
write.csv(x = trait_phencor_ann_toprint, file = save_file, quote = FALSE, row.names = FALSE)


## Save all analysis data
save_file <- file.path(results_dir, "phenotype_analysis_data.RData")
save("trait_varcomp", "S2_MET_BLUEs_use", "trait_phencor", "trait_main_BLUEs", file = save_file)







## GxE analysis
## 
## Find the top performing lines in each environment for each trait
## 

# First split traits into whether greater or lower values are more favorable
lower_fav_traits <- c("HeadingDate", "PlantHeight", "GrainProtein", "BetaGlucan",
                      "FreeAminoNitrogen", "WortProtein")

# How many of the top n should be filtered?
n_top <- 5

# Change the sign of these traits
enviro_best_lines <- S2_MET_BLUEs_use %>% 
  mutate(value = if_else(trait %in% lower_fav_traits, value * -1, value)) %>% 
  select(environment, trait, line_name, value) %>% 
  group_by(trait, environment) %>% 
  top_n(n = n_top, wt = value)

enviro_best_line <- enviro_best_lines %>% 
  top_n(1)

# Group by trait and line name and find the number of times that each line is
# in the n_top
enviro_best_lines_count <- enviro_best_lines %>% 
  group_by(trait, line_name) %>% 
  summarize(n_best = n())

enviro_best_line_count <- enviro_best_line %>% 
  group_by(trait, line_name) %>% 
  summarize(n_best = n())

# # Look at grain yield
# # Plot a histogram
# g_best_line_gy <- enviro_best_lines_count %>% 
#   filter(trait == "GrainYield") %>% 
#   ggplot(aes(x = n_best)) +
#   geom_histogram(binwidth = 1) + 
#   xlab("Number of Environments in Which the Line Was Superior (Top 5)") +
#   ylab("Number of Lines") +
#   labs(title = "Grain Yield") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.9, margin = margin(t = 10, b = -20)))
# 
# # Save the plot
# save_file <- file.path(fig_dir, "enviro_best_line_grainyield.jpg")
# ggsave(filename = save_file, plot = g_best_line_gy, height = 5, width = 5, dpi = 1000)
# 


## For grain, find the rank of the top 5 genotypes from one environment in all
## other environments.
# First assign rank to all genotypes in all environments
grainyield_env_rank <- S2_MET_BLUEs_use %>%
  filter(trait == "GrainYield") %>% 
  group_by(environment) %>% 
  mutate(rank = rank(-value), n = n(), percentile = rank / n) %>%
  ungroup()

enviro_best_lines_allrank <- enviro_best_lines %>% 
  ungroup() %>%
  filter(trait == "GrainYield") %>%
  left_join(., select(grainyield_env_rank, environment, line_name, value, rank, percentile), by = "line_name")


## Look at some of the lines that were best in many environments
lines_of_interest <- enviro_best_lines_count %>% 
  filter(trait == "GrainYield")  %>% 
  arrange(desc(n_best)) %>% 
  head(3) %>%
  pull(line_name)

## Print the ordered rank of these lines in all environments
enviro_best_lines_allrank %>% 
  filter(line_name %in% lines_of_interest) %>% 
  distinct(line_name, environment.y, value.y, rank, percentile) %>% 
  arrange(line_name, rank) %>% View






