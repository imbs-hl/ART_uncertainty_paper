##' With this script the simulations from Kronziel et al. "Prediction Beyond 
##' Point Estimates: Artificial Trees with Uncertainty" can be reproduced. 
##' Please note, that the simulations in the paper were performed using 
##' batchtools on  a high throughout batch system. This script will implement 
##' the same calculations on your local system, which may lead to a high 
##' computation time. Comments will show you where you can save time or 
##' incorporate your own batch system. 

#---------------------------------------

## Load libraries
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}
library(pacman)
packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "DescTools", "caret", "this.path", "devtools")
p_load(packages, character.only = TRUE)

if("timbR" %in% installed.packages()){
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}

#---------------------------------------
## Define directories
## Please define your main directory here. 
## This should be the directory you cloned the git repository into.
main_dir <- this.dir()
setwd(main_dir)

## Define functions directory
fun_dir <- file.path(main_dir, "functions")
## Create and define registry directory
dir.create(file.path(main_dir, "registries"), showWarnings = FALSE)
reg_dir <- file.path(main_dir, "registries")
## Create and define data directory
dir.create(file.path(main_dir, "data"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "data")


# --------------------------------------------------- #
#                  Data Simulation                    #
# --------------------------------------------------- #
## In this part the data will be simulated and saved 
## for later use

#---------------------------------------
## load functions needed for the simulation

## functions to simulate the data sets
source("functions/simulate_rf_setting_1.R")
source("functions/simulate_rf_setting_2.R")
source("functions/simulate_rf_setting_3.R")
source("functions/simulate_rf_setting_4.R")
source("functions/simulate_rf_setting_5.R")

## functions to build artificial representative trees (ART) and decision trees (DT) 
source("functions/get_art.R")
source("functions/get_dt.R")
source("functions/get_seperate_dt.R")
source("functions/get_seperate_art.R")


#---------------------------------------
## set parameters for simulation

## parameter of data set
n         <- 1000    ## Number of samples in train data set
n_test    <- 100     ## Number of samples in test data set
n_cal     <- 1000    ## Number of samples in calidation data set
p         <- 100     ## Number of variables 
eps       <- 1       ## Simulated noise in data set

## parameter of random forest (RF)
num.trees <- 500     ## Number of trees in random forest
mtry     <- sqrt(p)  ## Mtry for random forest
min_node_size <- 100 ## Minimal node size for random forest

## parameter of ART and DT used in paper are listed as comments below, to save time a reduced selection is used here
metric   <- c("splitting variables") ## Simularity / distance measure for selecting MRT or building ART
probs_quantiles <- list(NULL) ## Use quantiles of split points instead of all split points for continuous variables when creating the ART to save time
min.bucket <- c(150) ## Minimal number of training observations reaching in each leaf
epsilon <- c(0) ## Continue adding more nodes to the ART if the similarity remains the same but the prediction improves by 1 - epsilon

# metric   <- c("weighted splitting variables", "splitting variables", "prediction") ## Simularity / distance measure for selecting MRT or building ART
# probs_quantiles <- list(c(0.25,0.5,0.75), c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9), NULL) ## Use quantiles of split points instead of all split points for continuous variables when creating the ART to save time
# min.bucket <- c(150, 200, 250) ## Minimal number of training observations reaching in each leaf
# epsilon <- c(0) ## Continue adding more nodes to the ART if the similarity remains the same but the prediction improves by 1 - epsilon



## parameter for conformal predictive system (CPS)
significance_level <- 0.05 # Significance level for uncertainty quantification with prediction interval 

## Number of times each experiment is repeated. You can save time here
## (in publication 100 was used, for time reasons 3 is used here)
repls <- 3
#---------------------------------------
## Create registry 
reg_name <- "simulated_results"
reg <- batchtools::makeExperimentRegistry(
  file.dir = file.path(reg_dir, reg_name),
  work.dir = main_dir,
  conf.file = NA, ## If you have a batch system, please enter conf file here,
  packages = c(packages, "timbR") ## Define which packages to use in your simulations
)
# ---------------------------------------
# Data Simulation                   
  
## Add problems ----
## There is a separate function for generating each of the settings from the paper.
## You can save time, excluding settings you are not interested in. 
batchtools::addProblem(name = "simulate_setting_1",
                       reg = reg,
                       fun = simulate_rf_setting_1,
                       data = n,
                       seed = 12345)
# batchtools::addProblem(name = "simulate_setting_2",
#                        reg = reg,
#                        fun = simulate_rf_setting_2,
#                        data = n,
#                        seed = 12345)
# batchtools::addProblem(name = "simulate_setting_3",
#                        reg = reg,
#                        fun = simulate_rf_setting_3,
#                        data = n,
#                        seed = 12345)
# batchtools::addProblem(name = "simulate_setting_4",
#                        reg = reg,
#                        fun = simulate_rf_setting_4,
#                        data = n,
#                        seed = 12345)
# batchtools::addProblem(name = "simulate_setting_5",
#                        reg = reg,
#                        fun = simulate_rf_setting_5,
#                        data = n,
#                        seed = 12345)

## Add algorithms to solve the problem ----
# ART with CPS
batchtools::addAlgorithm(reg = reg,
                         name = "get_art",
                         fun = get_art
)
# DT with CPS
batchtools::addAlgorithm(reg = reg,
                         name = "get_dt",
                         fun = get_dt
)
# One DT for regression, one for probability
batchtools::addAlgorithm(reg = reg,
                         name = "get_seperate_dt",
                         fun = get_seperate_dt
)
# One ART for regression, one for probability
batchtools::addAlgorithm(reg = reg,
                         name = "get_seperate_art",
                         fun = get_seperate_art
)

## define problem and algorithm designs
prob.designs <- list(
  simulate_setting_1 = data.frame(p_eff       = 5, ## Number of true effect variables
                                  beta_eff    = 2, ## Effect size of true effect variables
                                  n_test      = n_test,
                                  n_val       = n_cal,
                                  p           = p,
                                  num.trees   = num.trees,
                                  eps         = eps,
                                  mtry        = mtry,
                                  min_node_size = min_node_size,
                                  stringsAsFactors = FALSE
  )
  # simulate_setting_2 = data.frame(p_eff       = 50,
  #                                 beta_eff    = 0.2,
  #                                 n_test      = n_test,
  #                                 n_val       = n_cal,
  #                                 p           = p,
  #                                 num.trees   = num.trees,
  #                                 eps         = eps,
  #                                 mtry        = mtry,
  #                                 min_node_size = min_node_size,
  #                                 stringsAsFactors = FALSE
  # ),
  # simulate_setting_3 = data.frame(p_eff       = 5,
  #                                 beta_eff    = 2,
  #                                 n_test      = n_test,
  #                                 n_val       = n_cal,
  #                                 p           = p,
  #                                 p_corr      = 5,
  #                                 n_blocks    = 5,
  #                                 cor         = 0.3,
  #                                 num.trees   = num.trees,
  #                                 eps         = eps,
  #                                 mtry        = mtry,
  #                                 min_node_size = min_node_size,
  #                                 stringsAsFactors = FALSE
  # ),
  # simulate_setting_4 = data.frame(p_eff       = 5,
  #                                 beta_eff    = 2,
  #                                 n_test      = n_test,
  #                                 n_val       = n_cal,
  #                                 p           = p,
  #                                 p_int       = 5,
  #                                 beta_int    = 2,
  #                                 num.trees   = num.trees,
  #                                 eps         = eps,
  #                                 mtry        = mtry,
  #                                 min_node_size = min_node_size,
  #                                 stringsAsFactors = FALSE
  # ),
  # simulate_setting_5 = data.frame(p_eff_bin   = 5,
  #                                 p_eff_con   = 5,
  #                                 beta_eff    = 2,
  #                                 n_test      = n_test,
  #                                 n_val       = n_cal,
  #                                 p           = p,
  #                                 num.trees   = num.trees,
  #                                 eps         = eps,
  #                                 mtry        = mtry,
  #                                 min_node_size = min_node_size,
  #                                 stringsAsFactors = FALSE
  # )
)


algo.designs <- list(
  get_art = expand.grid(metric = metric, 
                         stringsAsFactors = FALSE,
                         probs_quantiles = probs_quantiles,
                         epsilon = epsilon,
                         min.bucket = min.bucket,
                         significance_level = significance_level),
  get_dt = expand.grid(stringsAsFactors = FALSE,
                       min.bucket = min.bucket,
                       significance_level = significance_level),
  get_seperate_dt = expand.grid(stringsAsFactors = FALSE,
                                min.bucket = min.bucket),
  get_seperate_art = expand.grid(metric = metric, 
                                stringsAsFactors = FALSE,
                                probs_quantiles = probs_quantiles,
                                epsilon = epsilon,
                                min.bucket = min.bucket)
)

## Add experiments ----
ids = batchtools::addExperiments(reg = reg,
                                 prob.designs = prob.designs,
                                 algo.designs = algo.designs,
                                 repls = repls)
ids[, chunk := 1]

summarizeExperiments(reg = reg)

## Test jobs before submission
# testJob(id = 1, reg = reg)

## Please change hyperparameters of submitJobs() here, if you have a batch system. 
submitJobs(ids = ids, reg = reg)

##' With pre selected parameters it will take around 10 min to complete.
##' Please note, the run times for the other settings could differ. 
##' Anyway simulating data for the figures in the paper will probably run for several days on you computer. 

## check job status
getStatus()

## Collect and save results ----
results <- reduceResultsList(reg = reg)

## Save calculate performance measures and job informations as data.frame
results_df <- lapply(results, function(x){x[[1]]}) %>% bind_rows()
saveRDS(results_df, file = file.path(proc_dir, paste0("results_", reg_name, ".rds")))

## Save regression trees as list
regression_trees <- lapply(results, function(x){
  x$regression_trees
})
saveRDS(regression_trees, file = file.path(proc_dir, paste0("regression_trees_", reg_name, ".rds")))

## Save probability trees as list
probability_trees <- lapply(results, function(x){
  x$probability_trees
})
saveRDS(probability_trees, file = file.path(proc_dir, paste0("probability_trees_", reg_name, ".rds")))

## Save predicted probabilies for CPS methods
pred_prob <- lapply(results, function(x){
  x$prob_df
})
saveRDS(pred_prob, file = file.path(proc_dir, paste0("pred_probabilities_", reg_name, ".rds")))
