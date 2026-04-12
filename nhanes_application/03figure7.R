##' With this script the figure 7 from Kronziel et al. "Prediction Beyond Point 
##' Estimates: Artificial Trees with Uncertainty" can be reproduced. 
##' Given some results from the nhanes application.
##' Run 02calculate_results.R to get such a data set. 

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

#------------------------------------------------------------------------------
# Load and prepare results data
#------------------------------------------------------------------------------
# Choose ONE of the following result files:

# (1)
# Results produced manually via `02calculate_results.R`
# Due to runtime reasons probs_quantiles = c(0.25,0.5,0.75) is used in default setting, you may change this
data_imp <- readRDS(file.path("data", "results_nhanes_application.rds"))%>%
  filter(probs_quantiles == "0.25,0.5,0.75" | is.na(probs_quantiles))

# (2)
# # Results used in the paper
# data_imp <- read.csv(file.path("data", "results_nhanes_application_results_from_paper.csv")) %>%
#   filter(probs_quantiles == "" | is.na(probs_quantiles))


# Harmonize method names and apply paper-consistent filtering
data <- data_imp %>%
  mutate(
    method = case_when(
      method == "Regression ART + CPS" ~ "ART + CPS",
      method == "Regression DT + CPS"  ~ "DT + CPS",
      method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
      method == "Regression DT + Probability DTs"  ~ "Mult. DTs",
      TRUE ~ method
    )
  ) %>%
  filter(min.bucket == 150) %>%
  filter(metric == "splitting variables" | is.na(metric))


#------------------------------------------------------------------------------
# Predictive accuracy: RMSE (glycohemoglobin)
#------------------------------------------------------------------------------
# RMSE is averaged per repetition
# Regression trees with and without CPS are identical, therefore
# only ART + CPS and DT + CPS are shown (plus RF as reference)

data_rmse <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  mutate(rmse = sqrt(mse_test_dat_tree)) %>%
  filter(method %in% c("ART + CPS", "DT + CPS", "RF")) %>%
  dplyr::summarise(dplyr::across(rmse, mean, na.rm = TRUE), .groups = "drop")

plot_rmse <- ggplot(data_rmse, aes(x = method, y = rmse, col = method)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "RMSE for glycohemoglobin") +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255)
  ))

plot_rmse


#------------------------------------------------------------------------------
# Predictive accuracy: Brier score (prediabetes)
#------------------------------------------------------------------------------
data_brier_57 <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(dplyr::across(brier_score5.7, mean, na.rm = TRUE), .groups = "drop")

plot_brier_prediabetes <- ggplot(
  data_brier_57,
  aes(x = method, y = brier_score5.7, col = method)
) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Brier score for prediabetes") +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255),
    rgb(255, 150, 90, maxColorValue = 255),
    rgb(80, 180, 210, maxColorValue = 255)
  ))

plot_brier_prediabetes


#------------------------------------------------------------------------------
# Predictive accuracy: Brier score (diabetes)
#------------------------------------------------------------------------------
data_brier_65 <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(dplyr::across(brier_score6.5, mean, na.rm = TRUE), .groups = "drop")

plot_brier_diabetes <- ggplot(
  data_brier_65,
  aes(x = method, y = brier_score6.5, col = method)
) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Brier score for diabetes") +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255),
    rgb(255, 150, 90, maxColorValue = 255),
    rgb(80, 180, 210, maxColorValue = 255)
  ))

plot_brier_diabetes


#------------------------------------------------------------------------------
# Uncertainty quality: empirical coverage
#------------------------------------------------------------------------------
data_coverage <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(dplyr::across(coverage, mean, na.rm = TRUE), .groups = "drop")

plot_coverage <- ggplot(
  data_coverage %>% filter(grepl("CPS", method)),
  aes(x = method, y = coverage, col = method)
) +
  geom_boxplot() +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.5) +
  ylim(0.9, 1) +
  theme_bw() +
  labs(x = "", y = "Coverage") +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255)
  ))

plot_coverage


#------------------------------------------------------------------------------
# Uncertainty quality: interval width
#------------------------------------------------------------------------------
data_interval_width <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(dplyr::across(interval_width, mean, na.rm = TRUE), .groups = "drop")

plot_width <- ggplot(
  data_interval_width %>% filter(grepl("CPS", method)),
  aes(x = method, y = interval_width, col = method)
) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "", y = "Interval width") +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_color_manual(values = c(
    rgb(203, 81, 25, maxColorValue = 255),
    rgb(0, 75, 90, maxColorValue = 255)
  ))

plot_width


#------------------------------------------------------------------------------
# Combine and save final figure
#------------------------------------------------------------------------------
prediction_accuracy <- plot_grid(
  plot_rmse,
  plot_brier_prediabetes,
  plot_brier_diabetes,
  plot_width,
  plot_coverage,
  labels = "AUTO",
  rel_widths = c(1, 1.5, 1.5, 1, 1),
  nrow = 1
)

prediction_accuracy

ggsave(
  prediction_accuracy,
  filename = file.path(img_dir, "fig7_nhanes_application_accuracy.png"),
  width = 45,
  height = 9,
  units = "cm",
  dpi = 200
)
