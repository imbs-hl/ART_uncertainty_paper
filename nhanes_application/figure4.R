#===============================================================================
# Reproduction of Figure 4 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 4 of the manuscript by comparing
# variable usage across different tree-based models.
#
# Specifically, it:
# - Loads precomputed results for Artificial Regression Trees (ARTs),
#   Decision Trees (DTs), and Random Forests (RFs)
# - Computes variable usage frequencies for ART + CPS and DT + CPS
# - Extracts and aggregates variable usage from RFs
# - Visualizes variable usage distributions using boxplots
#
# Prerequisites:
# - Run `01_prepare_data.R` to generate the preprocessed NHANES dataset
# - Run `02calculate_results.R` OR use the provided paper results
# - Clone the git repository locally
#
# Output:
# - Figure showing variable usage frequencies across models
#   saved to the output directory
#
# Author: Lea Kronziel
#===============================================================================


#------------------------------------------------------------------------------
# Define directories
#------------------------------------------------------------------------------
# Set the main directory (root of the cloned git repository)
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
  cowplot,
  reshape2
)


#------------------------------------------------------------------------------
# Load results and data
#------------------------------------------------------------------------------
# Precomputed results used in the paper (recommended)
usage_melt <- readRDS(file.path("data", "registry_art_cps_df_paper_fig4.rds"))

# If you want to recompute variable usage from your own experiment runs,
# uncomment and use the following lines instead
# Due to runtime reasons probs_quantiles = c(0.25,0.5,0.75) is used in default setting, you may change this

# data_imp      <- readRDS(file.path("data", "registry_art_cps_df.rds"))
# data_list_imp <- readRDS(file.path("data", "list_registry_art_cps.rds"))
# 
# # Harmonize method names
# data <- data_imp %>%
#   mutate(
#     method = case_when(
#       method == "Regression ART + CPS" ~ "ART + CPS",
#       method == "Regression DT + CPS"  ~ "DT + CPS",
#       method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
#       method == "Regression DT + Probability DTs"  ~ "Mult. DTs",
#       TRUE ~ method
#     )
#   )
# 
# # Load preprocessed NHANES dataset
# nhanes_data_imp <- read.csv(file.path("data", "nhanes_prepandemic_complete.csv"))
# 
# 
# #------------------------------------------------------------------------------
# # Prepare data for prediction-based comparisons
# #------------------------------------------------------------------------------
# # Predictions are computed on the full dataset, as no independent
# # test set remains and all trees were trained on equally sized samples
# 
# nhanes_data <- nhanes_data_imp %>%
#   select(-SEQN, -prediabetes)
# 
# 
# #------------------------------------------------------------------------------
# # Extract trained trees
# #------------------------------------------------------------------------------
# # Collect regression trees for each fold and repetition
# trees <- lapply(data_list_imp, function(x) x$tree)
# 
# 
# #------------------------------------------------------------------------------
# # Variable usage for glycohemoglobin
# #------------------------------------------------------------------------------
# # Step 1: Compute variable usage for ART + CPS and DT + CPS
# #------------------------------------------------------------------------------
# method_values <- c("ART + CPS", "DT + CPS")
# num_features  <- ncol(nhanes_data) - 1
# 
# usage_df <- data.frame(
#   repitition = NA,
#   method     = NA,
#   metric     = NA,
#   min.bucket = NA,
#   type       = NA,
#   t(rep(NA, num_features)),
#   t(rep(NA, num_features))
# )
# 
# for (i in seq_len(max(data$repitition))) {
#   for (method in method_values) {
#     
#     # Select trees corresponding to the current configuration
#     ids <- which(
#       data$method == method &
#         data$repitition == i &
#         (data$metric == "splitting variables" | is.na(data$metric)) &
#         data$min.bucket == 150 &
#         (data$probs_quantiles == "0.25,0.5,0.75" | is.na(data$probs_quantiles))
#     )
#     
#     selected_trees <- trees[ids]
#     
#     # Binary indicator: whether a variable is used at least once
#     usage_binary <- lapply(selected_trees, function(tree) {
#       vars <- treeInfo(tree)$splitvarID
#       indicator <- rep(0, num_features)
#       indicator[vars + 1] <- 1
#       indicator
#     })
#     
#     # Absolute split counts per variable
#     usage_absolute <- lapply(selected_trees, function(tree) {
#       vars <- treeInfo(tree)$splitvarID
#       tabulate(vars + 1, nbins = num_features)
#     })
#     
#     usage_binary   <- bind_cols(usage_binary)
#     usage_absolute <- bind_cols(usage_absolute)
#     
#     usage_mean     <- rowMeans(usage_binary)
#     usage_abs_mean <- rowMeans(usage_absolute)
#     
#     usage_df <- rbind(
#       usage_df,
#       data.frame(
#         repitition = i,
#         method     = method,
#         metric     = "prediction",
#         min.bucket = 150,
#         type       = "glycohemoglobin",
#         t(usage_mean),
#         t(usage_abs_mean)
#       )
#     )
#   }
# }
# 
# 
# #------------------------------------------------------------------------------
# # Reshape ART / DT variable usage for plotting
# #------------------------------------------------------------------------------
# usage_melt_art <- reshape2::melt(
#   usage_df[-1, 1:(num_features + 5)],
#   id.vars = c("repitition", "method", "metric", "min.bucket", "type")
# ) %>%
#   mutate(
#     variable = case_when(
#       variable == "X1" ~ colnames(nhanes_data)[1],
#       variable == "X2" ~ colnames(nhanes_data)[2],
#       variable == "X3" ~ colnames(nhanes_data)[3],
#       variable == "X4" ~ colnames(nhanes_data)[4],
#       variable == "X5" ~ colnames(nhanes_data)[5],
#       variable == "X6" ~ colnames(nhanes_data)[6],
#       variable == "X7" ~ colnames(nhanes_data)[7],
#       variable == "X8" ~ colnames(nhanes_data)[8],
#       variable == "X9" ~ colnames(nhanes_data)[9],
#       TRUE ~ variable
#     )
#   )
# 
# 
# #------------------------------------------------------------------------------
# # Step 2: Variable usage for Random Forests
# #------------------------------------------------------------------------------
# rf_ids <- which(data_imp$method == "RF")
# rf_usage <- data_imp[rf_ids, ]
# 
# rf_usage_avg <- rf_usage[, c(1:2, 7, 26:35)] %>%
#   group_by(method, metric, min.bucket, repitition) %>%
#   dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
# 
# usage_melt_rf <- reshape2::melt(
#   rf_usage_avg,
#   id.vars = c("repitition", "method", "metric", "min.bucket")
# )
# 
# 
# #------------------------------------------------------------------------------
# # Combine ART/DT and RF variable usage
# #------------------------------------------------------------------------------
# usage_melt <- bind_rows(usage_melt_art, usage_melt_rf)


#------------------------------------------------------------------------------
# Step 3: Plot variable usage
#------------------------------------------------------------------------------
plot_variable_usage <- ggplot(
  usage_melt,
  aes(x = method, y = value, col = factor(variable))
) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "",
    y = "Variable usage frequency",
    col = "Variables"
  ) +
  theme(
    text = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.key.size = unit(0.6, "cm")
  )

plot_variable_usage


#------------------------------------------------------------------------------
# Save figure
#------------------------------------------------------------------------------
ggsave(
  plot_variable_usage,
  filename = file.path(out_dir, "fig4_variable_usage.png"),
  width = 15,
  height = 8,
  units = "cm",
  dpi = 200
)
