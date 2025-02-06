#' With this script the benchmark evaluation from Kronziel et al. "Increasing the 
#' explainability of artificial representative trees through conformal 
#' prediction to quantify uncertainty" can be reproduced. 
#' Please note, that the simulations in the paper were performed using 
#' batchtools on  a high throughout batch system. This script will implement 
#' the same calculations on your local system, which may lead to a high 
#' computation time. Comments will show you where you can save time or 
#' incorporate your own batch system. 


#---------------------------------------
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

#---------------------------------------
# Load libraries
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}

pacman::p_load(batchtools)
pacman::p_load(ranger)
pacman::p_load(devtools)
pacman::p_load(rpart)
pacman::p_load(dplyr)
pacman::p_load(OpenML)
pacman::p_load(foreign)
pacman::p_load(caret)
pacman::p_load(farff)

if("timbR" %in% installed.packages()){
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "develop")
  library(timbR)
}

# --------------------------------------------------- #
#                  Benchmark Data                     #
# --------------------------------------------------- #
# In this part the data will be simulated and saved 
# for later use

#---------------------------------------
# load functions needed for the simulation

# functions to generate the ART
source("functions/simulate_benchmark_data.R")

# functions to generate the ART
source("functions/calculate_art_uncertainty_benchmark_data.R")

#---------------------------------------
# define parameters

# choose data set from OpenML 
# (https://openml.org/search?type=benchmark&sort=tasks_included&study_type=task&id=269)
task_ids <- c(233215) # for example Mercedes_Benz_Greener_Manufacturing


# parameter of random forest and ART
num.trees <- 500
min_node_size <- c(100)
metric   <- c("weighted splitting variables")
imp.num.var <- c(5) # p ggf später
probs_quantiles <- list(c(0.25,0.5,0.75)) 
epsilon <- c(0.05)
min.bucket <- 100                              # In Publication 25 was also used
significance_level <- seq(0.1, 0.9, 0.2)       # In Publication seq(0.025, 0.975, 0.025) was used

# replicates of simulation (cv is used)
repls <- 1
num_folds <- 7
current_fold <- rep(1:num_folds, length(task_ids))

#---------------------------------------
# Create registry 
reg_name <- "simulate_benchmark"
reg <- batchtools::makeExperimentRegistry(
  file.dir = file.path(reg_dir, reg_name),
  work.dir = main_dir,
  conf.file = NA, # If you have a batch system, please enter conf file here,
  packages = c("ranger", "timbR", "rpart", "dplyr", "OpenML", "foreign", "caret") # Define which packages to use in your simulations
)


#---------------------------------------
# Data Simulation                   



# Add the problem (cross validation for different data sets)
batchtools::addProblem(name = "simulate_benchmark_data",
                       reg = reg, 
                       fun = simulate_benchmark_data,
                       data = 1,
                       seed = 12345)

# Add algorithms to solve the problem 
batchtools::addAlgorithm(reg = reg,
                         name = "artificial_rep_tree",
                         fun = calculate_art_rep_tree
)



# define problem and algorithm designs
prob.designs <- list(
  simulate_benchmark_data = data.frame(sim_id = rep(task_ids, each = num_folds),
                                       current_fold = current_fold,
                                       num_folds = num_folds,
                                       num.trees = num.trees
  )
)


algo.designs <- list(
  artificial_rep_tree = expand.grid(imp.num.var = imp.num.var,
                                    metric = metric, 
                                    stringsAsFactors = FALSE,
                                    probs_quantiles = probs_quantiles,
                                    epsilon = epsilon,
                                    min.bucket = min.bucket, 
                                    significance_level = significance_level)
)

# Add experiments 
ids = batchtools::addExperiments(reg = reg,
                                 prob.designs = prob.designs,
                                 algo.designs = algo.designs,
                                 repls = repls 
)

summarizeExperiments(reg = reg)


# Test jobs before submission
# testJob(id = 1, reg = reg)

# Please change this if you have a batch system. 

submitJobs(ids = ids, reg = reg)

#' With pre selected parameters it will take around 5 min to complete.
#' Please note, the run times for the other settings could differ. 
#' Anyway simulating data for the figures in the paper will probably run for several days on you computer. 

# check job status
getStatus()

# Collect and save results ----
results <- reduceResultsList(reg = reg, missing.val = 0)
# Save results
saveRDS(results, file = file.path(proc_dir, "results_benchmark.Rds"))













