#===============================================================================
# Reproduction of Figure 3 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 3 of the manuscript by analyzing the *stability*
# and *interpretability* of different tree-based models.
#
# Specifically, it:
# - Loads precomputed experiment results
# - Compares model stability across repetitions and folds
# - Computes distances between predictions and splitting-variable usage
# - Visualizes stability metrics using boxplots
#
# Models compared:
# - Artificial Regression Trees (ARTs)
# - Decision Trees (DTs)
# - CPS-enhanced variants
# - Multiple-tree probability models
#
# Prerequisites:
# - Run `01_prepare_data.R` to generate the NHANES dataset
# - Run `02calculate_results.R` OR use the provided paper results
# - Clone the git repository locally
#
# Output:
# - Stability plots (prediction distance & SV distance)
#   saved to the output directory
#
# Author: Lea Kronziel
#===============================================================================


#------------------------------------------------------------------------------
# Define directories
#------------------------------------------------------------------------------
# Set main directory (root of cloned repository)
main_dir <- getwd()
setwd(main_dir)

# Create processing directory
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")

# Create output directory
dir.create(file.path(main_dir, "output"), showWarnings = FALSE)
out_dir <- file.path(main_dir, "output")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
if (!"pacman" %in% installed.packages()) {
  install.packages("pacman")
}

pacman::p_load(
  ggplot2,
  gridExtra,
  ranger,
  devtools,
  rpart,
  dplyr,
  tidyr,
  reshape2,
  cowplot
)


#------------------------------------------------------------------------------
# Load and prepare results data
#------------------------------------------------------------------------------
# Choose ONE of the following result sources

# Results used in the paper
data_imp      <- readRDS(file.path("data", "registry_art_cps_df_paper.rds"))
data_list_imp <- readRDS(file.path("data", "list_registry_art_cps_paper.rds"))

# Result
