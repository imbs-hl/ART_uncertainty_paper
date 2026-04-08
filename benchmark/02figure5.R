##' Benchmark evaluation and visualization from:
##' Kronziel et al. (2024): "Predicting Medical Outcomes using Artificial 
##' Representative Trees with Uncertainty Quantification".
##'
##' This script reproduces Figure 5 based on benchmark experiment results.
##' It requires precomputed results from at least one dataset.
##' Run "01run_benchmark.R" beforehand to generate these results.
##'
##' Note:
##' - The original paper evaluates multiple datasets and parameter settings.
##' - To reduce runtime on local machines, this script uses a reduced subset
##'   (e.g., fewer parameter configurations).
##'
##' Output:
##' - Boxplots for interpretability measures:
##'   * Maximum tree depth (longest path)
##'   * Number of terminal nodes (leaves)
##' - Combined figure summarizing interpretability


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
  
  # Filter parameter settings (reduced configuration for runtime reasons)
  filter(min.bucket == 150) %>% 
  filter(metric == "splitting variables" | is.na(metric)) %>% 
  
  # Use quantiles (instead of all split points) for ART-based methods
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
## Tree depth (maximum path length)

data_max_depth <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket)

# Boxplot for maximum tree depth
plot_tree_depth <- ggplot(data_max_depth, 
                          aes(x = method, y = as.numeric(max_depth), col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ ., nrow = 3) +
  theme_bw() +
  labs(x = "", y = "Maximum depth (longest path)") +
  theme(text = element_text(size = 12),
        legend.position = "bottom",
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255,150,90, maxColorValue = 255),
                                rgb(80,180,210, maxColorValue = 255)))

plot_tree_depth


#---------------------------------------
## Number of leaves (terminal nodes)

data_num_leaves <- data %>% 
  group_by(method, dataset_name, metric, probs_quantiles, min.bucket)

# Boxplot for number of leaves
plot_num_leaves <- ggplot(data_num_leaves, 
                          aes(x = method, y = as.numeric(num_leaves), col = method)) +
  geom_boxplot() +
  facet_wrap(dataset_name ~ ., nrow = 3) +
  theme_bw() +
  labs(x = "", y = "Number of leaves") +
  theme(text = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255,150,90, maxColorValue = 255),
                                rgb(80,180,210, maxColorValue = 255)))

plot_num_leaves


#---------------------------------------
## Combine interpretability plots

interpretability_tree_depth_num_leaves <- plot_grid(
  plot_num_leaves,
  plot_tree_depth,
  labels = "AUTO",
  rel_widths = c(1, 1),
  nrow = 2
)

interpretability_tree_depth_num_leaves

# Save combined figure
ggsave(interpretability_tree_depth_num_leaves,
       filename = file.path(img_dir, "fig5_benchmark_interpretability_tree_depth_num_leaves.png"),
       width = 25, height = 25, units = "cm", dpi = 200)