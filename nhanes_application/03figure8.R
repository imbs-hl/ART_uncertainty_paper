##' With this script the figure 8 from Kronziel et al. "Prediction Beyond Point 
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
data_imp <- readRDS(file.path("data", "results_nhanes_application.rds")) %>% 
  filter(probs_quantiles == "0.25,0.5,0.75" | is.na(probs_quantiles))

# (2)
# # Results used in the paper
# data_imp <- read.csv(file.path("data", "results_nhanes_application_results_from_paper.csv")) %>%
#   filter(probs_quantiles == "" | is.na(probs_quantiles))


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
  filter(metric == "splitting variables" | is.na(metric))


#------------------------------------------------------------------------------
# Interpretability metric 1: Maximum tree depth
#------------------------------------------------------------------------------
# Restrict to CPS-based models
data_max_depth <- data %>%
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(
    max_depth = mean(max_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
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
  dplyr::summarise(
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
  filename = file.path(img_dir, "fig8_nhanes_application_interpretability.png"),
  width = 12,
  height = 6,
  units = "cm",
  dpi = 200
)
