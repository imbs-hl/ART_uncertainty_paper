#' With this script the figure 1 from "Uncertainty 
#' quantification enhances the explainability of 
#' artificial representative trees" can be reproduced. 
#' Given a simulated data set. Run simulations.R to get such a data set. 
#' With the standard parameters in simulations.R, only the part for 
#' data scenario 1 is reproduced using less significance levels and repetitions 
#' than in the manusscript, as the runtime without a computing cluster 
#' would otherwise be too high.
#---------------------------------------
# Define directories
# Please define your main directory here. 
# This should be the directory you cloned the git repository into.
main_dir <- getwd()
setwd(main_dir)

# Create and define proc directory
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")
# Create and define output directory
dir.create(file.path(main_dir, "output"), showWarnings = FALSE)
out_dir <- file.path(main_dir, "output")

#---------------------------------------
# Load libraries
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}

pacman::p_load(ggplot2)
pacman::p_load(gridExtra)
pacman::p_load(ranger)
pacman::p_load(devtools)
pacman::p_load(rpart)
pacman::p_load(dplyr)
pacman::p_load(tidyr)


if("timbR" %in% installed.packages()){
  library(timbR)
  warning("Please check, if timbR version 3.1 is installed.")
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}

#---------------------------------------
# Functions to simulate the data sets
source("functions/simulate_rf_setting_1.R")

#---------------------------------------
# Simulate the data from scenario 1
# parameter of data set
n         <- 1000
n_test    <- 1000
n_cal     <- 1000
p         <- 100
eps      <- 1

# Parameter of random forest and ART
num.trees <- 500
mtry     <- sqrt(p)
min_node_size <- c(100)
metric   <- c("weighted splitting variables")
imp.num.var <- c(5) 
epsilon <- c(0.05)
significance_level <- 0.05
min.bucket <- 100

#---------------------------------------
# Build ART
set.seed(1234)

# If you want to use other settings please change also the input parameters
instance <- simulate_rf_setting_1(data = n,
                                  p_eff       = 5,
                                  beta_eff    = 2,
                                  n_test      = n_test,
                                  n_cal       = n_cal, 
                                  p           = p, 
                                  num.trees   = num.trees,
                                  eps         = eps,
                                  mtry        = mtry,
                                  min_node_size = min_node_size
)

## Exctract data from instance 
test_dat       <- instance[[1]]
rf             <- instance[[2]]
params         <- instance[[3]]
cal_dat        <- instance[[4]]
effect_var_ids <- instance[[5]]
noise_var_ids  <- instance[[6]]
train_dat      <- instance[[7]]
dependent_varname <- instance[[8]]

## Generate artificial rep tree
rf_rep <- generate_tree(rf = rf, metric = metric, test_data = test_dat, train_data = train_dat, 
                            importance.mode = TRUE, imp.num.var = imp.num.var, dependent_varname = dependent_varname,
                            probs_quantiles = NULL, epsilon = epsilon,
                            min.node.size = 100, num.splits = 3)

tree_info <- treeInfo(rf_rep)

#---------------------------------------
# Calculate predictions for calibration and test data for uncertainty quantification
y_cal <- cal_dat[,dependent_varname]
y_test <- test_dat[,dependent_varname]
y_cal_pred <- predict(rf_rep, cal_dat)$predictions
y_test_pred <- predict(rf_rep, test_dat)$predictions

#---------------------------------------
# Calibration using Mondrian inductive conformal prediction (ICP)
# Get calibrated prediction for test data
calibrated_mondrian_predictions <- timbR::get_calibrated_prediction_regression_mondrian(y_cal_pred = y_cal_pred, 
                                                                                        y_cal = y_cal, 
                                                                                        y_test_pred = y_test_pred, 
                                                                                        significance_level = significance_level,
                                                                                        tree = rf_rep,
                                                                                        cal_data = cal_dat, 
                                                                                        test_data = test_dat,
                                                                                        dependent_varname = dependent_varname,
                                                                                        show_node_id = TRUE
)

# Round bounds to save space
tree_info_conformal_mondrian_pred <- tree_info %>% left_join(unique(calibrated_mondrian_predictions)) %>% 
  mutate(prediction = round(prediction, 2),
         lower_bound = round(lower_bound, 2),
         upper_bound = round(upper_bound, 2))

# Save plot
plot_tree(tree_info_df = tree_info_conformal_mondrian_pred, train_data_df = train_dat, test_data_df = test_dat, rf_list = rf_rep, tree_number = 1, 
          dependent_var = "y",
          show_sample_size = FALSE, show_prediction_nodes = FALSE, show_uncertainty = TRUE, show_coverage = TRUE, show_intervalwidth = TRUE,
          vert_sep = 5, hor_sep = 5,
          work_dir = out_dir, plot_name = "fig1_tree_scenario1_conformal_mondrian_prediction", colors = NULL)
