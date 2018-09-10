## S2MET Registration - Source Script
## 
## This script will load necessary directories and files for use in the
## S2MET_Registration project
## 

# Load libraries
library(tidyverse)
library(readxl)
library(knitr)
library(cowplot)
library(lme4)
library(modelr)
library(broom)
library(boot)

# Personal libraries
library(neyhart) # https://github.com/neyhartj/neyhart
library(pbr) # https://github.com/neyhartj/pbr
library(barleypheno) # https://github.com/neyhartj/barleypheno


# Directories
proj_dir <- repo_dir
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
results_dir <- file.path(proj_dir, "Results")
gdrive_dir <- "C:/Users/jln54/GoogleDrive"

## Read in data relevant to the project
# Read in trial information
trial_info <- read_csv(file = file.path(proj_dir, "Data/trial_metadata.csv"))

# Read in the entry list
entry_list <- read_excel(path = file.path(proj_dir, "Data/project_entries.xlsx"))


# Read in the S2MET tidy phenotypic data
load("path/to/phenotypic/data/S2_MET_tidy.RData")
# Read in the marker data
load("path/to/genotype/data/S2_genos_mat.RData")


# Remove some traits
traits_remove <- c("BacterialLeafSteakSeverity", "BarleyColor", "FHBIncidence", 
                   "Lodging", "LodgingDegree", "Nodding", "PowderyMildew", "StrawBreakage",
                   "WortViscosity", "WortClarity", "WortColor", "KernelWeight", "MaturityDate",
                   "ThinGrains")


# Remove the traits from the tidy dataset
S2_MET_tidy_use <- S2_MET_tidy %>%
  filter(!trait %in% traits_remove)

# Assign the checks
checks <- subset(entry_list, class == "Check", line_name, drop = T)


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(class == "S2TP") %>% 
  pull(line_name )

vp <- entry_list %>% 
  filter(class == "S2C1R") %>% 
  pull(line_name )

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))


## Functions
## Bootstrap a correlation coefficient
## alpha can be a vector
boot_cor <- function(x, y, boot.reps = 1000, alpha = 0.05) {
  
  # Error
  boot.reps <- as.integer(boot.reps)
  
  # Prob must be between 0 and 1
  alpha_check <- alpha > 0 | alpha < 1
  
  if (!all(alpha_check))
    stop("'alpha' must be between 0 and 1.")
  
  # Define a function for the correlation
  boot.cor <- function(input.data, i) {
    rep_data <- input.data[i,]
    return(cor(rep_data[,1], rep_data[,2]))
  }
  
  
  # First calculate the base statistic
  base_cor <- suppressWarnings(cor(x, y))
  
  # If the correlation is not NA, proceed
  if (!is.na(base_cor)) {
    
    # Perform the bootstrapping
    boot_results <- boot(data = cbind(x, y), statistic = boot.cor, R = boot.reps)
    
    # Standard error
    se <- sd(boot_results$t)
    # Bias
    bias <- mean(boot_results$t) - base_cor
    
    
    # Confidence interval
    ci_upper <- quantile(boot_results$t, 1 - (alpha / 2))
    ci_lower <- quantile(boot_results$t, (alpha / 2))
    
  } else {
    
    se <- bias <- ci_lower <- ci_upper <- NA
    
  }
  
  # Assemble list and return
  data.frame(cor = base_cor, se = se, bias = bias, alpha = alpha,
             ci_lower = ci_lower, ci_upper = ci_upper, row.names = NULL)
}



