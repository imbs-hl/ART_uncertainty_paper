#================================================================================
# Script to reproduce the results from the paper using batchtools
#================================================================================
# This script sets up and runs a series of experiments using the BatchTools 
# framework to evaluate Artificial Representative Trees (ARTs) and Decision Trees (DTs)
# with conformal predictive systems (CPS) for predicting glycohemoglobin levels 
# and probabilities of prediabetes and diabetes.
#
# Main functionalities:
# 1. Load required R packages and user-defined functions.
# 2. Set paths, registry directories, and data locations.
# 3. Define simulation parameters for cross-validation (CV), random forests (RF), DTs, 
#    and ARTs, including metrics, probability quantiles, min.bucket sizes, 
#    epsilon, and significance levels.
# 4. Create and load BatchTools experiment registries.
# 5. Add problems (data splits for cross-validation) and algorithms (ART, DT, separate ARTs/DTs, RF variable usage).
# 6. Define problem and algorithm designs for parameter sweeps.
# 7. Add experiments and assign jobs to clusters (fast or batch partitions) with appropriate resources.
# 8. Submit, monitor, and resubmit jobs as needed (handling waiting, expired, or failed jobs).
# 9. Collect and save experiment results as data frames and RDS files.
#
# Notes:
# - Single ARTs and DTs are trained for predicting all three outcomes: continuous outcomes (glycohemoglobin) and probabilities of prediabetes and diabetes
# - Separate ARTs or DTs can be used for regression and probability outcomes individually.
# - Performance metrics include MSE, Brier scores, interval coverage, tree depth, number of leaves, and variable usage.
#
# Author: Lea Kronziel
#================================================================================


# Packages
#----------------

if("timbR" %in% installed.packages()){
  warning("Please check, if timbR with at least version 3.3 is installed.")
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}

packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "timbR", "DescTools", "caret", "this.path")
invisible(sapply(packages, library, character.only = TRUE))

# Paths
#----------------
# Define directories
# Please define your main directory here. 
# This should be the directory you cloned the git repository into.
main_dir <- getwd()
setwd(main_dir)

# Create and define registry directory
dir.create(file.path(main_dir, "registries"), showWarnings = FALSE)
reg_dir <- file.path(main_dir, "registries")
# Create and define functions directory
dir.create(file.path(main_dir, "functions"), showWarnings = FALSE)
fun_dir <- file.path(main_dir, "functions")
# Create and define proc directory
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")
# Create and define data directory
dir.create(file.path(main_dir, "data"), showWarnings = FALSE)
data_dir <- file.path(main_dir, "data")

#---------------------------------------
# Load functions needed for the cross validation (cv)
source("functions/get_cv_data.R")

# Functions to generate the ARTs and DTs 
source("functions/get_art.R")
source("functions/get_dt.R")
source("functions/get_seperate_dt.R")
source("functions/get_seperate_art.R")

# Function to get frequency of split variables of random forest (RF)
source("functions/get_rf_variable_usage.R")

#---------------------------------------
# Set method parameters
# Reduced setting selected here, values used in paper are commented behind each parameter

# Number of folds of cv and number of repeated cvs
num_folds <- 4 # 10
repitition <- c(1:2)# c(1:20)

# DT and ART hyperparameters
# ART: Distance metric (ART vs. RF)
metric   <- c("splitting variables") # c("splitting variables", "weighted splitting variables", "prediction")

# ART: Using quantiles instead of all split from RF for ARTs
probs_quantiles <- list(c(0.25,0.5,0.75)) # list(c(0.25,0.5,0.75), c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), NULL) 

# ART: Improvement of prediction quality if similarity to RF stays the same
epsilon <- c(0)

# ART and DT: Stop tree growth so that at least min.bucket observations of train data end in each leaf
min.bucket <- c(150) # c(150, 200, 250)

# ART and DT: significance level for prediction interval of CPS
significance_level <- 0.05

#---------------------------------------
# Create registry
reg_name <- "registry_art_cps" 

reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, reg_name),
  work.dir = proc_dir,
  conf.file = NA, # If you have a batch system, please enter conf file here,
  packages = packages  # Define which packages to use
)

#---------------------------------------
# Add the problem (cv of data set)
addProblem(
  name = "get_cv_data",
  reg = reg, 
  fun = get_cv_data,
  data = 1,
  seed = 12345
)

#---------------------------------------
# Add the algorithms to solve the problem: ARTs and DTs
addAlgorithm(reg = reg, name = "get_art", fun = get_art)
addAlgorithm(reg = reg, name = "get_dt", fun = get_dt)
addAlgorithm(reg = reg, name = "get_seperate_dt", fun = get_seperate_dt)
addAlgorithm(reg = reg, name = "get_seperate_art", fun = get_seperate_art)

# Additionaly calculate frequency of split variables used in RF to compare it with ART and DT
addAlgorithm(reg = reg, name = "get_rf_variable_usage", fun = get_rf_variable_usage)

#---------------------------------------
# Define problem and algorithm designs
prob.designs <- list(
  get_cv_data = expand.grid(
    current_fold = 1:num_folds,
    num_folds = num_folds,
    repitition = repitition,
    stringsAsFactors = FALSE
  )
)

algo.designs <- list(
  get_art = expand.grid(
    metric = metric, 
    stringsAsFactors = FALSE,
    probs_quantiles = probs_quantiles,
    epsilon = epsilon,
    min.bucket = min.bucket,
    significance_level = significance_level
  ),
  get_dt = expand.grid(
    stringsAsFactors = FALSE,
    min.bucket = min.bucket,
    significance_level = significance_level
  ),
  get_seperate_dt = expand.grid(
    stringsAsFactors = FALSE,
    min.bucket = min.bucket
  ),
  get_seperate_art = expand.grid(
    metric = metric, 
    stringsAsFactors = FALSE,
    probs_quantiles = probs_quantiles,
    epsilon = epsilon,
    min.bucket = min.bucket
  ),
  get_rf_variable_usage = data.frame()
)

# Add the experiments
ids = addExperiments(
  reg = reg,
  prob.designs = prob.designs,
  algo.designs = algo.designs,
  repls = 1
)
#---------------------------------------

summarizeExperiments(reg = reg)

# Please change this if you have a batch system. 
submitJobs(ids = ids, reg = reg)

#' With pre selected parameters it will take around 5-10 min to complete.
#' Please note, the run times for the other parameter settings could differ. 
#' Anyway simulating data for the figures in the paper will probably run for several days on you computer. 

# check job status
getStatus()

# Collect and save results ----
results <- reduceResultsList(reg = reg, missing.val = 0)
# Results dataframe
results_df <- lapply(results, function(x){x[[1]]}) %>% bind_rows()
saveRDS(results_df, file = file.path("data", paste0(reg_name, "_df.rds")))
saveRDS(results, file = file.path("data", paste0("list_", reg_name, ".rds")))
