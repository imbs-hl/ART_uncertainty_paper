##' With this script the figure 1 from Kronziel et al. "Prediction Beyond Point 
##' Estimates: Artificial Trees with Uncertainty" can be reproduced.
##' Given the prepared nhanes data set 
##' that you get after running 01prepare_data.R.
#---------------------------------------
## Load required libraries

# Install pacman if not available
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}

library(pacman)

# List of required packages
packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "DescTools", "caret", "cowplot","this.path")

# Load all packages
p_load(packages, character.only = TRUE)

# Load timbR (install from GitHub if necessary)
if("timbR" %in% installed.packages()){
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}


#---------------------------------------
## Define directories

# Set main directory (root of the cloned repository)
main_dir <- this.dir()
setwd(main_dir)

# Directory containing processed benchmark results
proc_dir <- file.path(main_dir, "data")

# Create directory for output plots if it does not exist
dir.create(file.path(main_dir, "img"), showWarnings = FALSE)
img_dir <- file.path(main_dir, "img")

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
# Note: This step may take 1–3 minutes due to data set size
art_SV <- generate_tree(
  rf,
  metric = "splitting variables",
  train_data,
  test_data = NULL,
  dependent_varname = "glycohemoglobin",
  min.bucket = 150, 
  num.splits=3
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
  work_dir          = img_dir,
  plot_name         = "fig1_nhanes_application_example_tree",
  colors            = NULL
)

#---------------------------------------
# End of script
#---------------------------------------
