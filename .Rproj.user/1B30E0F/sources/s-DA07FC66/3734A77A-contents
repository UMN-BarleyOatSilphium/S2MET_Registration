## Project: S2MET Registration
## 
## Phenotypic data adjustment and analysis
## Author: Jeff Neyhart
## Last updated: March 15, 2018
## 

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

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
    
    # Extract the function and the formula, then fit
    if (is.na(df$blk[1])) {
      f <- forms$form1
      
      fit <- lm(formula = f, data = df)
      
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



## Variance component table

# Fit mixed models (with random effects for everthing except location and year)
# to determine variance components.
# Use a likelihood ratio test to assess the significance of each of the variance
# comonents
ranef_terms <- c("(1|line_name)", "(1|line_name:year)", "(1|line_name:location)", "(1|line_name:year:location)")
fixef_terms <- c("year", "location")

## Fit each model to each trait
# We will need to drop models that are not relevant to certain traits (i.e. if 
# those traits were not evaluated in more than one year/location/etc.)


# Control for lmer
control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")


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
    n_loc <- n_distinct(df$location)
    n_year <- n_distinct(df$year)
    n_env <- n_distinct(df$environment)
    
    # Determine the harmonic mean of the number of locations, years, and reps
    plot_table <- xtabs(~ line_name + location + year, df)
    
    ## Harmonic means
    # Locations
    harm_loc <- apply(X = plot_table, MARGIN = c(1,2), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Year
    harm_year <- apply(X = plot_table, MARGIN = c(1,3), sum) %>% 
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Reps
    harm_rep <- apply(X = plot_table, MARGIN = c(1,2,3), sum) %>% 
      harm_mean()
    
    
    # Remove terms based on the number of locations or years
    if (all(n_loc == 1, n_year == 1)) {
      # Remove any terms with location and year
      remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "location|year")]
      removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "location|year")]
      
    } else if (n_loc == 1) {
      # Remove any term with location
      remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "location")]
      removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "location")]
      
    } else if (n_year == 1) {
      # Remove any term with year
      remove_fixef_terms <- fixef_terms[str_detect(fixef_terms, "year")]
      removed_ranef_terms <- ranef_terms[str_detect(ranef_terms, "year")]
      
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
    
    # Calculate heritability
    h2 <- herit(object = fits$full, n_l = harm_loc, n_y = harm_year, n_r = harm_rep,
                exp = "line_name / (line_name + (line_name:year / n_y) + (line_name:location / n_l) + 
                (line_name:year:location / (n_l * n_y)) + (Residual / (n_l * n_y * n_r)))")
    
    # Return the heritability
    var_comp2 <- var_comp1 %>% mutate(h2 = h2)
    
    # Return the variance components and the full model fit
    data_frame(var_comp = list(var_comp2), full_fit = fits[1])
    
  })

# Save this data
save_file <- file.path(proj_dir, "Results/all_trait_variance_component_heritability.RData")
save("trait_varcomp", "S2_MET_BLUEs_use", file = save_file)





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
    variance = round(variance, 3),
    variance_notation = str_trim(str_c(variance, sig_notation, sep = " "))) %>%
  left_join(., t3.traits, by = c("trait" = "Nickname")) 

# Create the table
trait_varcomp_sig_toprint <- trait_varcomp_sig %>% 
  select(term, Trait, variance_notation) %>% 
  mutate(term = str_replace_all(string = term, pattern = ":", replacement = " %*% ") %>% str_to_title(),
         term = str_replace_all(string = term, "Line_name", "Genotype"),
         Trait = str_to_title(Trait)) %>%
  spread(Trait, variance_notation) %>%
  rename(Term = term)


save_file <- file.path(fig_dir, "trait_variance_components.csv")
write_csv(x = trait_varcomp_sig_toprint, path = save_file)


## Summarize the mean, min, max, sd, and heritability of each trait
trait_summary <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>% 
  summarize_at(vars(value), funs(min, max, mean, sd)) %>%
  left_join(., distinct(unnest(trait_varcomp, var_comp), trait, h2))

# Format to print
trait_summary_toprint <- trait_summary %>% 
  left_join(., select(t3.traits, Trait, Nickname), by = c("trait" = "Nickname")) %>% 
  mutate(Trait = str_to_title(Trait)) %>% 
  select(-trait) %>% 
  gather(parameter, value, -Trait) %>% 
  mutate(value = round(value, 3)) %>% # Round the values
  spread(Trait, value) %>% 
  rename(Trait = parameter) %>%
  mutate(Trait = factor(Trait, levels = names(trait_summary)[-1])) %>%
  arrange(Trait)

save_file <- file.path(fig_dir, "trait_value_summary.csv")
write_csv(x = trait_summary_toprint, path = save_file)

## Calculate genetic correlations among the genotype BLUPs
trait_geno_BLUP <- trait_varcomp %>% 
  do({.$full_fit[[1]] %>% ranef() %>% as.data.frame() %>% filter(grpvar == "line_name") %>% 
      select(line_name = grp, value = condval) }) %>%
  ungroup()


## Calculate the genetic correlations
trait_gencor <- trait_geno_BLUP %>% 
  spread(trait, value) %>% 
  select(-line_name) %>% 
  cor()

# Convert to df
trait_gencor_df <- trait_gencor %>% 
  as.data.frame() %>% 
  rownames_to_column("trait1") %>% 
  gather(trait2, correlation, -trait1)

## Plot the correlation
# Color scheme
color <- umn_palette(2, 4)

# Cluster the correlations using the heatmap function
# Abbreviate the trait names
trait_order <- row.names(trait_gencor)[heatmap(trait_gencor)$rowInd] %>%
  abbreviate(minlength = 2)

# Plot using ggplot
g_trait_gencor <- trait_gencor_df %>%
  mutate_at(vars(contains("trait")), ~str_replace_all(., trait_order) %>% 
              factor(., levels = trait_order)) %>% 
  ggplot(aes(x = trait1, y = trait2, fill = correlation)) + 
  geom_tile() + 
  scale_fill_gradient2(low = color[3], mid = "white", high = color[1], name = "Genetic\nCorrelation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank())

# Save the correlation matrix
save_file <- file.path(fig_dir, "trait_genetic_correlations.csv")
write.csv(x = trait_gencor, file = save_file, quote = FALSE, row.names = TRUE)

# Save the graph
save_file <- file.path(fig_dir, "trait_genetic_correlations.jpg")
ggsave(filename = save_file, plot = g_trait_gencor, width = 6, height = 5, dpi = 1000)



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

# Look at grain yiedl
# Plot a histogram
g_best_line_gy <- enviro_best_lines_count %>% 
  filter(trait == "GrainYield") %>% 
  ggplot(aes(x = n_best)) +
  geom_histogram(binwidth = 1) + 
  xlab("Number of Environments in Which the Line Was Superior (Top 5)") +
  ylab("Number of Lines") +
  labs(title = "Grain Yield") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.9, margin = margin(t = 10, b = -20)))

# Save the plot
save_file <- file.path(fig_dir, "enviro_best_line_grainyield.jpg")
ggsave(filename = save_file, plot = g_best_line_gy, height = 5, width = 5, dpi = 1000)








