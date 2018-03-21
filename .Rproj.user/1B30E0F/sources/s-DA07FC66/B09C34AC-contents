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

# Personal libraries
library(neyhart)
library(pbr)
library(fbutils)


# Directories
proj_dir <- repo_dir
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")

## Read in data relevant to the project
# Read in trial information
trial_info <- read_csv(file = file.path(proj_dir, "Data/trial_metadata.csv"))

# Read in the entry list
entry_list <- read_excel(path = file.path(proj_dir, "Data/project_entries.xlsx"))


# Read in the S2MET tidy information
load("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Phenotypic Data/Final/Master Phenotypes/S2_MET_tidy.RData")


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

