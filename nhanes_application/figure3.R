#===============================================================================
# Reproduction of Figure 3 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 3 of the manuscript by analyzing the *stability*
# and *interpretability* of different tree-based models.
#
# Specifically, it:
# - Loads precomputed experiment results (distance measures are precomputed
#   since storing all trees exceeds GitHub size limits; please email for full access)
# - Compares model stability across repetitions and cross-validation folds
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
# - Due to runtime constraints, probs_quantiles = c(0.25, 0.5, 0.75) is used by default
# - Clone the git repository locally
#
# Output:
# - Stability plots (prediction distance and SV distance)
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
dist_df <- readRDS(file.path("data", "registry_art_cps_df_paper_fig3.rds"))

# For paper results, precomputed distance measures are used.
# If you want to recreate the plots using your own results, uncomment below.
# Due to runtime constraints, probs_quantiles = c(0.25, 0.5, 0.75) is used by default.

# # Results produced manually via `02calculate_results.R`
# data_imp      <- readRDS(file.path("data", "registry_art_cps_df.rds"))
# data_list_imp <- readRDS(file.path("data", "list_registry_art_cps.rds"))
# 
# 
# # Harmonize method names for consistent labeling across plots
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
# nhanes_imp <- read.csv(file.path("data", "nhanes_prepandemic_complete.csv"))
# 
# #------------------------------------------------------------------------------
# # Load NHANES data used for stability evaluation
# #------------------------------------------------------------------------------
# # Predictions are computed on the full dataset to ensure fairness,
# # since all models were trained on equally sized training samples.
# nhanes_data <- nhanes_imp %>%
#   select(-SEQN, -prediabetes)
# 
# 
# #------------------------------------------------------------------------------
# # Extract trees and prediction objects
# #------------------------------------------------------------------------------
# # List of fitted regression trees
# trees <- lapply(data_list_imp, function(x) x$tree)
# 
# # Probability trees
# trees5.7 <- lapply(data_list_imp, function(x) x$art_prob5.7)
# trees6.5 <- lapply(data_list_imp, function(x) x$art_prob6.5)
# 
# # Node-level probability tables from RF
# pred_prob <- lapply(data_list_imp, function(x) x$prob_df)
# 
# 
# #------------------------------------------------------------------------------
# # Hyperparameter configuration
# #------------------------------------------------------------------------------
# method_values <- c("ART + CPS", "DT + CPS")
# num_features  <- ncol(nhanes_data) - 1
# 
# 
# #------------------------------------------------------------------------------
# # Stability analysis: glycohemoglobin (regression outcome)
# #------------------------------------------------------------------------------
# # Distance between predictions and between splitting-variable usage
# dist_cps <- data.frame(
#   repitition = NA, method = NA, metric = NA, min.bucket = NA,
#   probs_quantiles = NA, dist_pred = NA, dist_sv = NA, type = NA
# )
# 
# for (i in 1:max(data$repitition)) {
#   for (method in method_values) {
#     
#     # Identify models for this configuration
#     ids <- which(
#       data$method == method &
#         data$repitition == i &
#         (data$metric == "splitting variables" | is.na(data$metric)) &
#         data$min.bucket == 150 &
#         (data$probs_quantiles == "0.25,0.5,0.75" | is.na(data$probs_quantiles))
#     )
#     
#     trees_sel <- trees[ids]
#     
#     #--------------------------------------------------------------
#     # Prediction distance (MSE between model predictions)
#     #--------------------------------------------------------------
#     preds <- lapply(trees_sel, function(x) {
#       predict(x, data = nhanes_data, predict.all = TRUE)$predictions
#     }) %>% bind_cols()
#     
#     dist_pred <- mean(
#       as.matrix(dist(t(preds), method = "euclidean"))^2 /
#         nrow(nhanes_data)
#     )
#     
#     #--------------------------------------------------------------
#     # Splitting-variable distance
#     #--------------------------------------------------------------
#     used_vars <- lapply(trees_sel, function(x) {
#       sv <- treeInfo(x)$splitvarID
#       u  <- rep(0, num_features)
#       u[sv + 1] <- 1
#       u
#     }) %>% bind_cols()
#     
#     dist_sv <- mean(
#       as.matrix(dist(t(used_vars), method = "euclidean"))^2 /
#         num_features
#     )
#     
#     # Store results
#     dist_cps <- rbind(
#       dist_cps,
#       data.frame(
#         repitition = i,
#         method = method,
#         metric = "splitting variables",
#         min.bucket = 150,
#         probs_quantiles = "0.25,0.5,0.75",
#         dist_pred = dist_pred,
#         dist_sv = dist_sv,
#         type = "glycohemoglobin"
#       )
#     )
#   }
# }
# 
# 
# #------------------------------------------------------------------------------
# # Probability outcomes: prediabetes and diabetes
# #------------------------------------------------------------------------------
# nhanes_data_5.7 <- nhanes_data %>%
#   mutate(glycohemoglobin = ifelse(glycohemoglobin >= 5.7, 1, 0)) %>%
#   select(-glycohemoglobin)
# 
# nhanes_data_6.5 <- nhanes_data %>%
#   mutate(glycohemoglobin = ifelse(glycohemoglobin >= 6.5, 1, 0)) %>%
#   select(-glycohemoglobin)
# 
# 
# # Initialize result containers
# dist_prediabetes <- data.frame(
#   repitition = NA, method = NA, metric = NA, min.bucket = NA,
#   probs_quantiles = NA, dist_pred = NA, dist_sv = NA, type = NA
# )
# 
# dist_diabetes <- dist_prediabetes
# 
# 
# #------------------------------------------------------------------------------
# # ARTs and DTs with CPS
# #------------------------------------------------------------------------------
# for (i in 1:max(data$repitition)) {
#   for (method in method_values) {
#     
#     # Identify models and corresponding prediction objects
#     ids_arts <- which(
#       data$method == method &
#         data$repitition == i &
#         data$min.bucket == 150 &
#         (data$probs_quantiles == "0.25,0.5,0.75" | is.na(data$probs_quantiles)) &
#         (data$metric == "splitting variables" | is.na(data$metric))
#     )
#     
#     trees_arts     <- trees[ids_arts]
#     pred_prob_arts <- pred_prob[ids_arts]
#     
#     probs_df5.7 <- matrix(NA, nrow = nrow(nhanes_data), ncol = max(data$fold))
#     probs_df6.5 <- matrix(NA, nrow = nrow(nhanes_data), ncol = max(data$fold))
#     
#     for (j in 1:max(data$fold)) {
#       
#       # Determine terminal nodes for each observation
#       nodes_hnannes <- predict(
#         data = nhanes_data,
#         trees_arts[[j]],
#         type = "terminalNodes"
#       )$predictions
#       
#       # Compute predicted probabilities for each model
#       pred_df <- data.frame(
#         nodeID = nodes_hnannes,
#         glycohemoglobin = nhanes_data$glycohemoglobin
#       ) %>%
#         left_join(pred_prob_arts[[j]])
#       
#       probs_df5.7[, j] <- pred_df$prob5.7
#       probs_df6.5[, j] <- pred_df$prob6.5
#     }
#     
#     # Compute prediction distances
#     dist_pred5.7 <- mean(
#       as.matrix(dist(t(probs_df5.7), method = "euclidian"))^2 /
#         nrow(nhanes_data)
#     )
#     dist_pred6.5 <- mean(
#       as.matrix(dist(t(probs_df6.5), method = "euclidian"))^2 /
#         nrow(nhanes_data)
#     )
#     
#     # Compute splitting-variable distance
#     used_variables <- lapply(trees_arts, function(x) {
#       splitting_variables <- treeInfo(x)$splitvarID
#       fu <- rep(0, num_features)
#       fu[splitting_variables + 1] <- 1
#       fu
#     }) %>% bind_cols()
#     
#     dist_sv <- mean(
#       as.matrix(dist(t(as.matrix(used_variables)), method = "euclidian"))^2 /
#         num_features
#     )
#     
#     dist_prediabetes <- rbind(
#       dist_prediabetes,
#       data.frame(
#         repitition = i, method = method,
#         metric = "splitting variables", min.bucket = 150,
#         probs_quantiles = "0.25,0.5,0.75",
#         dist_pred = dist_pred5.7,
#         dist_sv = dist_sv,
#         type = "prediabetes"
#       )
#     )
#     
#     dist_diabetes <- rbind(
#       dist_diabetes,
#       data.frame(
#         repitition = i, method = method,
#         metric = "splitting variables", min.bucket = 150,
#         probs_quantiles = "0.25,0.5,0.75",
#         dist_pred = dist_pred6.5,
#         dist_sv = dist_sv,
#         type = "diabetes"
#       )
#     )
#   }
# }
# 
# 
# #------------------------------------------------------------------------------
# # Multiple ARTs and DTs
# #------------------------------------------------------------------------------
# for (i in 1:max(data$repitition)) {
#   for (method in c("Mult. ARTs", "Mult. DTs")) {
#     
#     # Identify models for this configuration
#     ids_arts <- which(
#       data$method == method &
#         data$repitition == i &
#         data$min.bucket == 150 &
#         (data$probs_quantiles == "0.25,0.5,0.75" | is.na(data$probs_quantiles)) &
#         (data$metric == "splitting variables" | is.na(data$metric))
#     )
#     
#     trees_arts     <- trees[ids_arts]
#     trees_arts5.7  <- trees5.7[ids_arts]
#     trees_arts6.5  <- trees6.5[ids_arts]
#     
#     # Compute predicted probabilities for all models
#     pred5.7 <- lapply(trees_arts5.7, function(x) {
#       predict(x, data = nhanes_data_5.7)$predictions[, 2]
#     }) %>% bind_cols()
#     
#     pred6.5 <- lapply(trees_arts6.5, function(x) {
#       predict(x, data = nhanes_data_6.5)$predictions[, 2]
#     }) %>% bind_cols()
#     
#     # Compute prediction distances
#     dist_pred_prob_arts_5.7 <- mean(
#       as.matrix(dist(t(pred5.7), method = "euclidian"))^2 /
#         nrow(nhanes_data)
#     )
#     dist_pred_prob_arts_6.5 <- mean(
#       as.matrix(dist(t(pred6.5), method = "euclidian"))^2 /
#         nrow(nhanes_data)
#     )
#     
#     # Compute splitting-variable distances
#     used_variables5.7 <- lapply(trees_arts5.7, function(x) {
#       sv <- treeInfo(x)$splitvarID
#       fu <- rep(0, num_features)
#       fu[sv + 1] <- 1
#       fu
#     }) %>% bind_cols()
#     
#     used_variables6.5 <- lapply(trees_arts6.5, function(x) {
#       sv <- treeInfo(x)$splitvarID
#       fu <- rep(0, num_features)
#       fu[sv + 1] <- 1
#       fu
#     }) %>% bind_cols()
#     
#     dist_sv_prob_dts5.7 <- mean(
#       as.matrix(dist(t(as.matrix(used_variables5.7)), method = "euclidian"))^2 /
#         num_features
#     )
#     dist_sv_prob_dts6.5 <- mean(
#       as.matrix(dist(t(as.matrix(used_variables6.5)), method = "euclidian"))^2 /
#         num_features
#     )
#     
#     dist_prediabetes <- rbind(
#       dist_prediabetes,
#       data.frame(
#         repitition = i, method = method,
#         metric = "splitting variables", min.bucket = 150,
#         probs_quantiles = "0.25,0.5,0.75",
#         dist_pred = dist_pred_prob_arts_5.7,
#         dist_sv = dist_sv_prob_dts5.7,
#         type = "prediabetes"
#       )
#     )
#     
#     dist_diabetes <- rbind(
#       dist_diabetes,
#       data.frame(
#         repitition = i, method = method,
#         metric = "splitting variables", min.bucket = 150,
#         probs_quantiles = "0.25,0.5,0.75",
#         dist_pred = dist_pred_prob_arts_6.5,
#         dist_sv = dist_sv_prob_dts6.5,
#         type = "diabetes"
#       )
#     )
#   }
# }
# 
# 
# #------------------------------------------------------------------------------
# # Combine all stability results
# #------------------------------------------------------------------------------
# dist_df <- bind_rows(
#   dist_cps[-1, ],
#   dist_prediabetes[-1, ],
#   dist_diabetes[-1, ]
# )


#------------------------------------------------------------------------------
# Plot stability results
#------------------------------------------------------------------------------
plot_dist_pred <- ggplot(dist_df, aes(x = method, y = dist_pred, col = method)) +
  geom_boxplot() +
  facet_grid(
    . ~ factor(
      type,
      labels = c("glycohemoglobin", "prediabetes", "diabetes"),
      levels = c("glycohemoglobin", "prediabetes", "diabetes")
    ),
    scales = "free_x"
  ) +
  theme_bw() +
  labs(x = "", y = "Prediction distance") +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(size = 18),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255),
    rgb(255, 150, 90, maxColorValue = 255),
    rgb(80, 180, 210, maxColorValue = 255)
  ))

plot_dist_sv <- ggplot(dist_df, aes(x = method, y = dist_sv, col = method)) +
  geom_boxplot() +
  facet_grid(
    . ~ factor(
      type,
      labels = c("glycohemoglobin", "prediabetes", "diabetes"),
      levels = c("glycohemoglobin", "prediabetes", "diabetes")
    ),
    scales = "free_x"
  ) +
  theme_bw() +
  labs(x = "", y = "SV distance") +
  theme(
    text = element_text(size = 15),
    strip.text = element_text(size = 18),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255),
    rgb(255, 150, 90, maxColorValue = 255),
    rgb(80, 180, 210, maxColorValue = 255)
  ))

plot_dist <- plot_grid(plot_dist_sv, plot_dist_pred, labels = "AUTO", nrow = 2)
plot_dist


#------------------------------------------------------------------------------
# Save figure
#------------------------------------------------------------------------------
ggsave(
  plot_dist,
  filename = file.path(out_dir, "fig3_stability.png"),
  width = 33,
  height = 14,
  units = "cm",
  dpi = 200
)
