#===============================================================================
# Reproduction of Figure 3 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 3 of the manuscript by analyzing the *stability*
# and *interpretability* of different tree-based models.
#
# Specifically, it:
# - Loads precomputed experiment results (precalculated disance measures as data 
#   with trees is too big for github, please send an email if you want access to 
#   all trees, train, test, and validation data)
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
dist_df      <- readRDS(file.path("data", "registry_art_cps_df_paper_fig3.rds"))
# For results from paper precalculated disances are used
# If you want to build the plots with your own results uncomment the following code


# Results produced manually via `02calculate_results.R`
# data_imp      <- readRDS(file.path("data", "registry_art_cps_df.rds"))
# data_list_imp <- readRDS(file.path("data", "list_registry_art_cps.rds"))

# 
# # Harmonize method names (consistent labeling across plots)
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
# 
# #------------------------------------------------------------------------------
# # Load NHANES data used for stability evaluation
# #------------------------------------------------------------------------------
# # Predictions are computed on the full dataset to ensure fairness:
# # all models were trained on equally sized training samples
# 
# nhannes_data <- nhannes_data_imp %>%
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
# # Node-level probability tables
# pred_prob <- lapply(data_list_imp, function(x) x$prob_df)
# 
# 
# #------------------------------------------------------------------------------
# # Hyperparameter configuration
# #------------------------------------------------------------------------------
# metric_values       <- "splitting variables"
# min.bucket_values   <- 150
# probs_quantiles     <- ""
# method_values       <- c("ART + CPS", "DT + CPS")
# 
# num_features <- ncol(nhannes_data) - 1
# 
# 
# #------------------------------------------------------------------------------
# # Stability analysis: glycohemoglobin (regression)
# #------------------------------------------------------------------------------
# # Distance between predictions and between splitting-variable usage
# dist_cps <- data.frame(
#   repitition = NA, method = NA, metric = NA, min.bucket = NA,
#   dist_pred = NA, dist_sv = NA, type = NA,
#   t(rep(NA, 9)), t(rep(NA, 9))
# )
# 
# for (i in 1:max(data$repitition)) {
#   for (metric in metric_values) {
#     for (min.bucket in min.bucket_values) {
#       for (method in method_values) {
#         
#         # Identify models for this configuration
#         ids <- which(
#           data$method == method &
#             data$repitition == i &
#             (data$metric == metric | is.na(data$metric)) &
#             data$min.bucket == min.bucket &
#             (data$probs_quantiles == "" | is.na(data$probs_quantiles))
#         )
#         
#         trees_sel <- trees[ids]
#         
#         #--------------------------------------------------------------
#         # Prediction distance (MSE between model predictions)
#         #--------------------------------------------------------------
#         preds <- lapply(trees_sel, function(x) {
#           predict(x, data = nhannes_data, predict.all = TRUE)$predictions
#         }) %>% bind_cols()
#         
#         dist_pred <- mean(
#           as.matrix(dist(t(preds), method = "euclidean"))^2 /
#             nrow(nhannes_data)
#         )
#         
#         #--------------------------------------------------------------
#         # Splitting-variable distance
#         #--------------------------------------------------------------
#         used_vars <- lapply(trees_sel, function(x) {
#           sv <- treeInfo(x)$splitvarID
#           u  <- rep(0, num_features)
#           u[sv + 1] <- 1
#           u
#         }) %>% bind_cols()
#         
#         dist_sv <- mean(
#           as.matrix(dist(t(used_vars), method = "euclidean"))^2 /
#             num_features
#         )
#         
#         # Mean usage per variable
#         used_vars_mean <- rowMeans(used_vars)
#         
#         # Absolute usage counts
#         used_vars_abs <- lapply(trees_sel, function(x) {
#           tabulate(treeInfo(x)$splitvarID + 1, nbins = 9)
#         }) %>% bind_cols()
#         
#         # Store results
#         dist_cps <- rbind(
#           dist_cps,
#           data.frame(
#             repitition = i,
#             method = method,
#             metric = metric,
#             min.bucket = min.bucket,
#             dist_pred = dist_pred,
#             dist_sv = dist_sv,
#             type = "glycohemoglobin",
#             t(used_vars_mean),
#             t(used_vars_abs)
#           )
#         )
#       }
#     }
#   }
# }
# 
# 
# #------------------------------------------------------------------------------
# # Probability outcomes: prediabetes and diabetes
# #------------------------------------------------------------------------------
# nhannes_data_5.7 <- nhannes_data %>%
#   mutate(glycohemoglobin = ifelse(glycohemoglobin >= 5.7, 1, 0)) %>%
#   select(-glycohemoglobin)
# 
# nhannes_data_6.5 <- nhannes_data %>%
#   mutate(glycohemoglobin = ifelse(glycohemoglobin >= 6.5, 1, 0)) %>%
#   select(-glycohemoglobin)
# 
# 
# # Initialize result containers
# dist_prediabetes <- data.frame(repitition = NA, method = NA, metric = NA,
#                                min.bucket = NA, dist_pred = NA,
#                                dist_sv = NA, type = NA)
# 
# dist_diabetes <- dist_prediabetes
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

plot_dist_pred <- ggplot(dist_df, aes(x = method, y = dist_pred, col = method))+#, col = factor(min.bucket)))+
  geom_boxplot() +
  facet_grid(.~factor(type, 
                      labels=c("glycohemoglobin", "prediabetes", "diabetes"), 
                      levels=c("glycohemoglobin", "prediabetes", "diabetes")), scales = "free_x") +
  theme_bw() +
  labs(x = "",
       y = "Prediction distance")+
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        legend.position = "none") +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 


plot_dist_sv <- ggplot(dist_df, aes(x = method, y = dist_sv, col = method))+#, col = factor(min.bucket)))+
  geom_boxplot() +
  facet_grid(.~factor(type, 
                      labels=c("glycohemoglobin", "prediabetes", "diabetes"), 
                      levels=c("glycohemoglobin", "prediabetes", "diabetes")), scales = "free_x") +
  theme_bw() +
  labs(x = "",
       y = "SV distance")+
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        legend.position = "none") +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 

plot_dist <- plot_grid(plot_dist_sv, plot_dist_pred,
                       labels = "AUTO", nrow = 2)

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
