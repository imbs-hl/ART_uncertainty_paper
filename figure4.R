#===============================================================================
# Reproduction of Figure 4 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 1 of the manuscript by comparing the
# interpretability of different tree-based models.
#
# Specifically, it:
# - Loads precomputed results for Artificial Regression Trees (ARTs),
#   Decision Trees (DTs), and their CPS-enhanced variants
# - Compares model complexity in terms of
#     (i) maximum tree depth and
#     (ii) number of terminal nodes
# - Visualizes these metrics using boxplots
#
# Prerequisites:
# - Run `01_prepare_data.R` to generate the preprocessed NHANES dataset
# - Run `02calculate_results.R` OR use the provided paper results
# - Clone the git repository locally
#
# Output:
# - Figure showing interpretability metrics (tree depth and number of leaves)
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

# Create and define processing directory
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")

# Create and define output directory
dir.create(file.path(main_dir, "output"), showWarnings = FALSE)
out_dir <- file.path(main_dir, "output")


#------------------------------------------------------------------------------
# Load required libraries
#------------------------------------------------------------------------------
# Install pacman if not available
if (!"pacman" %in% installed.packages()) {
  install.packages("pacman")
}

# Load all required packages
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
# Choose one of the following result files:

# Results used in the paper
data_imp <- readRDS(file.path("data", "registry_art_cps_df_paper.rds"))

# Results produced manually via `02calculate_results.R`
# data_imp <- readRDS(file.path("data", "registry_art_cps_df.rds"))


# Harmonize method names and apply filtering consistent with the paper
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
  filter(metric == "splitting variables" | is.na(metric)) %>%
  filter(probs_quantiles == "" | is.na(probs_quantiles))


#------------------------------------------------------------------------------
# Interpretability metric 1: Maximum tree depth
#------------------------------------------------------------------------------
data_max_depth <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  summarise(
    max_depth = mean(max_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Restrict to CPS-based models
  filter(grepl("CPS", method))

plot_tree_depth <- ggplot(
  data_max_depth,
  aes(x = method, y = as.numeric(max_depth), col = method)
) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "",
    y = "Deepest path"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 15)
  ) +
  scale_color_manual(
    values = c(
      rgb(203, 81, 25, maxColorValue = 255),
      rgb(0, 75, 90, maxColorValue = 255)
    )
  )

plot_tree_depth


#------------------------------------------------------------------------------
# Interpretability metric 2: Number of terminal nodes
#------------------------------------------------------------------------------
data_num_leaves <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  summarise(
    num_leaves = mean(num_leaves, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(grepl("CPS", method))

plot_num_leaves <- ggplot(
  data_num_leaves,
  aes(x = method, y = as.numeric(num_leaves), col = method)
) +
  geom_boxplot() +
  theme_bw() +
  labs(
    x = "",
    y = "Number of leaves"
  ) +
  theme(
    text = element_text(size = 12),
    legend.position = "none",
    strip.text = element_text(size = 15)
  ) +
  scale_color_manual(
    values = c(
      rgb(203, 81, 25, maxColorValue = 255),
      rgb(0, 75, 90, maxColorValue = 255)
    )
  )

plot_num_leaves


#------------------------------------------------------------------------------
# Combine interpretability plots
#------------------------------------------------------------------------------
interpretability_tree_depth_num_leaves <- plot_grid(
  plot_num_leaves,
  plot_tree_depth,
  labels = "AUTO",
  rel_widths = c(1, 1),
  nrow = 1
)

interpretability_tree_depth_num_leaves


#------------------------------------------------------------------------------
# Save figure
#------------------------------------------------------------------------------
ggsave(
  interpretability_tree_depth_num_leaves,
  filename = file.path(out_dir, "fig2_interpretability.png"),
  width = 12,
  height = 6,
  units = "cm",
  dpi = 200
)
