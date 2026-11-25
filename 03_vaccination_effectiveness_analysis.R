# =============================================================================
# 03_vaccination_effectiveness_analysis.R
# VEBIS training project - influenza A(H1N1)pdm09 (synthetic data)
# =============================================================================

# vaccine effectiveness analysis for influenza A(H1N1)pdm09 using synthetic data
# VE = (1 - aOR) * 100, where aOR is the adjusted odds ratio
#                       of vaccination among cases vs controls
# aOR: cofounder-adjusted odds ratio of vaccination among cases vs controls
#      from logistic regression model.
#      adjuster for age, sex, at least one condition, onset date and country


# 1. Load required packages ----------------------------------------------------
library(dplyr)
library(broom)
library(ggplot2)

# 2. Load data -----------------------------------------------------------------
h1n1_cd1 <- readRDS("data/processed/H1N1_CD1.rds")
h1n1_cd2 <- readRDS("data/processed/H1N1_CD2.rds")


