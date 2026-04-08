##' Benchmark evaluation and visualization from:
##' Kronziel et al. (2024): "Predicting Medical Outcomes using Artificial 
##' Representative Trees with Uncertainty Quantification".
##'
##' This script reproduces Figure 4 based on benchmark experiment results.
##' It requires precomputed results from at least one dataset.
##' Run "01run_benchmark.R" beforehand to generate these results.
##'
##' Note:
##' - The original paper uses extensive parameter settings and multiple datasets.
##' - To reduce runtime on local machines, this script typically focuses on a subset
##'   (e.g., dataset "wineW" and reduced parameter configurations).
##'
##' Output:
##' - Boxplots for RMSE, Brier Score, Coverage, and Interval Width
##' - Combined figure summarizing prediction performance


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
              "gridExtra", "DescTools", "caret", "cowplot")

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
## Load and preprocess benchmark results

data <- readRDS(file.path(proc_dir, "results_benchmark_experiments.rds")) %>% 
  # Rename methods for better readability in plots
  mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
                            method == "Regression DT + CPS" ~ "DT + CPS",
                            method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
                            method == "Regression DT + Probability DTs" ~ "Mult. DTs",
                            TRUE ~ method)) %>% 
  filter(min.bucket == 150) %>% 
  filter(metric == "splitting variables" | is.na(metric)) %>% 
  # Filter parameter settings (reduced configuration for runtime reasons)
  # Use quantiles (instead of all split points) for ART + CPS
  filter(probs_quantiles == "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9" & method == "ART + CPS" | 
           probs_quantiles == "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9" & method == "Mult. ARTs" | 
           is.na(probs_quantiles))

# if you want to use the original data from the publication, please download them (see README for details), unpack them, and move them into the data folder: 
# 
# data from publication
#
# data <- readRDS(file.path(proc_dir, "results_benchmark_experiments__benchmark_results_from_paper.rds")) %>% 
#   # Rename methods for better readability in plots
#   mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
#                             method == "Regression DT + CPS" ~ "DT + CPS",
#                             method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
#                             method == "Regression DT + Probability DTs" ~ "Mult. DTs",
#                             TRUE ~ method)) %>% 
#   filter(min.bucket == 150) %>% 
#   filter(metric == "splitting variables" | is.na(metric)) %>% 
#   # Filter parameter settings (reduced configuration for runtime reasons)
#   # Use quantiles (instead of all split points) for ART + CPS
#   filter(probs_quantiles == "" & method == "ART + CPS" | 
#            probs_quantiles == "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9" & method == "Mult. ARTs" | 
#            is.na(probs_quantiles))


#---------------------------------------
## RMSE (Root Mean Squared Error)

data_rmse <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket) %>%
  mutate(rmse = sqrt(mse_test_dat_tree)) 

# Create RMSE boxplot
plot_rmse <- ggplot(data_rmse, aes(x = method, y = rmse, col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ .) +
  theme_bw() +
  labs(x = "", y = "RMSE") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255,150,90, maxColorValue = 255),
                                rgb(80,180,210, maxColorValue = 255)))

plot_rmse


#---------------------------------------
## Brier Score (probability calibration)

data_brier_score <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket)

plot_brier_score <- ggplot(data_brier_score, aes(x = method, y = brier_score0.5, col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ ., nrow = 3) +
  theme_bw() +
  labs(x = "", y = "Brier score") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.text = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255,150,90, maxColorValue = 255),
                                rgb(80,180,210, maxColorValue = 255)))

plot_brier_score


#---------------------------------------
## Coverage (prediction interval validity)

data_coverage <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket) %>% 
  filter(grepl("CPS", method))  # only CPS-based methods

plot_coverage <- ggplot(data_coverage, aes(x = method, y = coverage, col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ ., nrow = 3) +
  theme_bw() +
  labs(x = "", y = "Coverage", col = "") +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.5) +
  ylim(0.75, 1) +
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255)))

plot_coverage


#---------------------------------------
## Interval width (uncertainty size)
data_interval_width <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket) %>% 
  filter(grepl("CPS", method))

# Plot interval widths
plot_width <- ggplot(data_interval_width, aes(x = method, y = interval_width, col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ ., nrow = 3) +
  theme_bw() +
  labs(x = "", y = "Interval width") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.text = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255)))

plot_width

# if you get the Warning message: "Removed xx rows containing non-finite outside the scale range (`stat_boxplot()`). "
# that means, that at least one border of xx predictions intervals is (+/-) infinite


#---------------------------------------
## Combine plots into one figure

prediction_accuracy <- plot_grid(plot_rmse, plot_brier_score, plot_width, plot_coverage,
                                 labels = "AUTO",
                                 rel_widths = c(1, 1, 1, 1),
                                 nrow = 4)

prediction_accuracy

# Save final figure
ggsave(prediction_accuracy,
       filename = file.path(img_dir, "fig4_benchmark_prediction_accuracy.png"),
       width = 50, height = 50, units = "cm", dpi = 200)
