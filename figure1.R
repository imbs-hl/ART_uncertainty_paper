#' Reproduction of Figure 1 from:
#' "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#'
#' This script trains an Artificial Regression Tree (ART) based on a random forest
#' and visualizes calibrated predictive intervals for glycohemoglobin.
#'
#' Prerequisites:
#' - Run `01_prepare_data.R` beforehand to generate the preprocessed NHANES dataset.
#' - The git repository must be cloned locally.
#'
#' Output:
#' - A visualization of the ART including uncertainty quantification
#'   is saved to the output directory.
#'
#' Author: Lea Kronziel
#---------------------------------------

#---------------------------------------
# Define directories
#---------------------------------------

# Set main working directory (root of the cloned repository)
main_dir <- getwd()
setwd(main_dir)

# Create directory for intermediate processing results
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")

# Create directory for output (plots, results)
dir.create(file.path(main_dir, "output"), showWarnings = FALSE)
out_dir <- file.path(main_dir, "output")


#---------------------------------------
# Load required libraries
#---------------------------------------

# Install pacman if not available
if (!"pacman" %in% installed.packages()) {
  install.packages("pacman")
}

# Load required packages
pacman::p_load(
  ggplot2,
  gridExtra,
  ranger,
  devtools,
  rpart,
  dplyr,
  tidyr
)

# Install and load timbR (development version if not installed)
if ("timbR" %in% installed.packages()) {
  library(timbR)
  warning("Please verify that timbR version 3.3 is installed.")
} else {
  devtools::install_github("imbs-hl/timbR", ref = "master")
  library(timbR)
}


#---------------------------------------
# Import and prepare data
#---------------------------------------

# Load preprocessed NHANES data
nhannes_data_imp <- read.csv(
  file.path("data", "nhanes_prepandemic_complete.csv")
)

# Remove ID variable and non-used outcome
nhannes_data <- nhannes_data_imp %>%
  select(-SEQN, -prediabetes)


#---------------------------------------
# Split data into training, testing and calibration sets
#---------------------------------------

n <- nrow(nhannes_data)
set.seed(123)

# 50% training data
train_ids <- sample(1:n, n * 0.5, replace = FALSE)

# 25% test data
test_ids <- sample(setdiff(1:n, train_ids), n * 0.25, replace = FALSE)

# Assign datasets
train_data <- nhannes_data[train_ids, ]
test_data  <- nhannes_data[test_ids, ]
cal_data   <- nhannes_data[-c(train_ids, test_ids), ]


#---------------------------------------
# Train Random Forest model
#---------------------------------------

set.seed(123)

# Random forest used as basis for ART generation
rf <- ranger(
  glycohemoglobin ~ .,
  data = train_data,
  min.node.size = nrow(train_data) * 0.1
)


#---------------------------------------
# ART + Conformal Predictive System (Splitting Variables)
#---------------------------------------

set.seed(123)

# Generate Artificial Regression Tree (ART)
# Note: This step may take 1–3 minutes due to dataset size
art_SV <- generate_tree(
  rf,
  metric = "splitting variables",
  train_data,
  test_data = NULL,
  dependent_varname = "glycohemoglobin",
  min.bucket = 150
)

# Build calibrated predictive system using Mondrian CPS
cpd <- get_calibrated_predictive_system_mondrian(
  y_cal_pred  = predict(art_SV, cal_data)$predictions,
  y_cal       = cal_data$glycohemoglobin,
  y_test_pred = predict(art_SV, test_data)$predictions,
  significance_level = 0.05,
  interval_type = "two-tailed",
  direction = NULL,
  show_node_id = TRUE,
  tree = art_SV,
  cal_data = cal_data,
  test_data = test_data,
  dependent_varname = "glycohemoglobin"
) %>%
  mutate(nodeID = leaf - 1)


#---------------------------------------
# Prepare tree information for plotting
#---------------------------------------

tree_info_df <- treeInfo(art_SV) %>%
  left_join(unique(cpd)) %>%
  mutate(
    lower_bound = round(lower_bound, 2),
    upper_bound = round(upper_bound, 2),
    prediction  = round(prediction, 2)
  )


#---------------------------------------
# Plot ART with uncertainty quantification
#---------------------------------------

plot_tree(
  tree_info_df      = tree_info_df,
  train_data_df     = train_data,
  test_data_df      = test_data,
  cal_data_df       = cal_data,
  rf_list           = art_SV,
  tree_number       = 1,
  dependent_var     = "glycohemoglobin",
  threshold         = c(5.7, 6.5),
  significance_level = 0.05,
  interval_type     = "two-tailed",
  direction         = NULL,
  show_sample_size  = FALSE,
  show_prediction_nodes = FALSE,
  show_uncertainty  = TRUE,
  show_coverage     = TRUE,
  show_intervalwidth = TRUE,
  show_cpd          = TRUE,
  cpd_plot_width    = 35,
  show_point_prediction = TRUE,
  show_prediction_interval = TRUE,
  vert_sep          = 15,
  hor_sep           = 25,
  work_dir          = out_dir,
  plot_name         = "fig1_art",
  colors            = NULL
)

#---------------------------------------
# End of script
#---------------------------------------
