#===============================================================================
# Reproduction of Figure 5 from:
# "Predicting Ordinal Outcome via Artificial Trees with Uncertainty Quantification"
#===============================================================================
#
# This script reproduces Figure 5 of the manuscript by comparing the
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


# MSE
#----------------
# pro repitition gemittelt
# Regressionbäume bei + CPs und ohne sind gleich, deswegen reichen 2 Abbildungen, einmal ARTs und einmal DTs
# RF als Vergleich
data_mse <- data %>% 
  # bind_rows(., data.frame(method = "RF", 
  #                         data %>% 
  #                           filter(method == "DT + CPS") %>% 
  #                           select(mse_test_dat_rf, repitition) %>% 
  #                           dplyr::rename(mse_test_dat_tree = mse_test_dat_rf))) %>% 
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  mutate(rmse = sqrt(mse_test_dat_tree)) %>% 
  filter(method == "ART + CPS" | method == "DT + CPS" | method == "RF") %>% 
  dplyr::summarise(across(rmse, mean, na.rm = TRUE), .groups = "drop") # %>% 
# mutate(method = case_when(method == "RF" ~ "RF",
#                           method == "ART + CPS" ~ "ART",
#                           method == "DT + CPS" ~ "DT"))

plot_mse <- ggplot(data_mse, aes(x=method, y= rmse, col = method))+#, col = factor(min.bucket))) +
  geom_boxplot()+
  #facet_wrap(metric~., nrow=1) %>% 
  theme_bw()+
  labs(x = "",
       y = "RMSE for glycohemoglobin")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255))) 

plot_mse


# Brier Score
#----------------
data_brier_score5.7 <- data %>% 
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(across(brier_score5.7, mean, na.rm = TRUE), .groups = "drop")

plot_brier_score_prediabetes <- ggplot(data_brier_score5.7, aes(x=method, y= brier_score5.7, col = method))+#, col = factor(min.bucket))) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "",
       y = "Brier score for prediabetes")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 
plot_brier_score_prediabetes

# ggsave(plot_brier_score_prediabetes, filename = file.path(plot_dir, "diabetes_accuracy_brier_score_prediabetes.png"), width = 15, height = 12, units = "cm", dpi = 200)

data_brier_score6.5 <- data %>% 
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(across(brier_score6.5, mean, na.rm = TRUE), .groups = "drop")

plot_brier_score_diabetes <- ggplot(data_brier_score6.5, aes(x=method, y= brier_score6.5, col = method))+#, col = factor(min.bucket))) +
  geom_boxplot() +
  theme_bw() +
  labs(x = "",
       y = "Brier score on test data",
       col = "min.bucket")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 
plot_brier_score_diabetes



# Coverage
#----------------
data_coverage <- data %>% 
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(across(coverage, mean, na.rm = TRUE), .groups = "drop")

plot_coverage <- ggplot(data_coverage %>% filter(grepl("CPS", method)), aes(x=method, y= coverage, col = method))+#, col = factor(min.bucket))) +
  geom_boxplot() +
  # facet_wrap(metric~., nrow=1) +
  theme_bw() +
  labs(x = "",
       y = "Coverage", 
       col = "min.bucket")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18))+
  geom_hline(yintercept = 0.95,
             linetype = "dashed",
             linewidth = 0.5)+
  ylim(0.9,1) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255))) 
plot_coverage

# interval width
#----------------
data_interval_width <- data %>% 
  group_by(method, repitition, metric, probs_quantiles, min.bucket) %>%
  dplyr::summarise(across(interval_width, mean, na.rm = TRUE), .groups = "drop")

plot_width <- ggplot(data_interval_width %>% filter(grepl("CPS", method)), aes(x=method, y= interval_width, col = method))+#, col = factor(min.bucket))) +
  geom_boxplot() +
  #facet_wrap(metric~., nrow=1) +
  theme_bw() +
  labs(x = "",
       y = "Interval width", 
       col = "min.bucket")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255))) 
plot_width





