##' With this script the figure 10 from Kronziel et al. "Prediction Beyond Point 
##' Estimates: Artificial Trees with Uncertainty" can be reproduced.  
##' Given some results from the nhanes application.
##' Run 02calculate_results.R to get such a data set and the used models. 
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
# Choose ONE of the following result sources
# Results produced manually via `02calculate_results.R`
# Due to runtime reasons probs_quantiles = c(0.25,0.5,0.75) is used in default setting, you may change this
data_imp <- readRDS(file.path(proc_dir, "results_nhanes_application.rds")) %>%
  filter(probs_quantiles == "0.25,0.5,0.75" | is.na(probs_quantiles)) # in Paper no quantiles were used, we used them here due to runtime

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
  ) 


# Lists of fitted regression trees
regression_trees <- readRDS(file.path(proc_dir, "regression_trees_nhanes_application.rds"))

# Load preprocessed NHANES dataset
nhanes_data_imp <- read.csv(file.path("data", "nhanes_prepandemic_complete.csv"))
nhanes_data <- nhanes_data_imp %>%
  select(-SEQN, -prediabetes)
#------------------------------------------------------------------------------
# Variable usage for glycohemoglobin
#------------------------------------------------------------------------------
# Step 1: Compute variable usage for ART + CPS and DT + CPS
#------------------------------------------------------------------------------
method_values <- c("ART + CPS", "DT + CPS")
num_features  <- ncol(nhanes_data) - 1

min.bucket <- 150
probs_quantiles = "0.25,0.5,0.75"

usage_df <- data.frame(
  repitition = NA,
  method     = NA,
  metric     = NA,
  min.bucket = NA,
  type       = NA,
  t(rep(NA, num_features)),
  t(rep(NA, num_features))
)

for (i in seq_len(max(data$repitition))) {
  for (method in method_values) {

    # Select trees corresponding to the current configuration
    ids <- which(
      data$method == method &
        data$repitition == i &
        (data$metric == "splitting variables" | is.na(data$metric)) &
        data$min.bucket == min.bucket &
        (data$probs_quantiles == probs_quantiles | is.na(data$probs_quantiles))
    )

    selected_trees <- regression_trees[ids]

    # Binary indicator: whether a variable is used at least once
    usage_binary <- lapply(selected_trees, function(tree) {
      vars <- treeInfo(tree)$splitvarID
      indicator <- rep(0, num_features)
      indicator[vars + 1] <- 1
      indicator
    })

    # Absolute split counts per variable
    usage_absolute <- lapply(selected_trees, function(tree) {
      vars <- treeInfo(tree)$splitvarID
      tabulate(vars + 1, nbins = num_features)
    })

    usage_binary   <- bind_cols(usage_binary)
    usage_absolute <- bind_cols(usage_absolute)

    usage_mean     <- rowMeans(usage_binary)
    usage_abs_mean <- rowMeans(usage_absolute)

    usage_df <- rbind(
      usage_df,
      data.frame(
        repitition = i,
        method     = method,
        metric     = "prediction",
        min.bucket = 150,
        type       = "glycohemoglobin",
        t(usage_mean),
        t(usage_abs_mean)
      )
    )
  }
}

#------------------------------------------------------------------------------
# Reshape ART / DT variable usage for plotting
#------------------------------------------------------------------------------
usage_melt_art <- reshape2::melt(
  usage_df[-1, 1:(num_features + 5)],
  id.vars = c("repitition", "method", "metric", "min.bucket", "type")
) %>%
  mutate(
    variable = case_when(
      variable == "X1" ~ colnames(nhanes_data)[1],
      variable == "X2" ~ colnames(nhanes_data)[2],
      variable == "X3" ~ colnames(nhanes_data)[3],
      variable == "X4" ~ colnames(nhanes_data)[4],
      variable == "X5" ~ colnames(nhanes_data)[5],
      variable == "X6" ~ colnames(nhanes_data)[6],
      variable == "X7" ~ colnames(nhanes_data)[7],
      variable == "X8" ~ colnames(nhanes_data)[8],
      variable == "X9" ~ colnames(nhanes_data)[9],
      TRUE ~ variable
    )
  )


#------------------------------------------------------------------------------
# Step 2: Variable usage for Random Forests
#------------------------------------------------------------------------------
rf_ids <- which(data_imp$method == "RF")
rf_usage <- data_imp[rf_ids, ]

rf_usage_avg <- rf_usage[, c(1:2, 7, 26:35)] %>%
  group_by(method, metric, min.bucket, repitition) %>%
  dplyr::summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

usage_melt_rf <- reshape2::melt(
  rf_usage_avg,
  id.vars = c("repitition", "method", "metric", "min.bucket")
)


#------------------------------------------------------------------------------
# Combine ART/DT and RF variable usage
#------------------------------------------------------------------------------
usage_melt <- bind_rows(usage_melt_art, usage_melt_rf)


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
  filename = file.path(img_dir, "fig10_nhanes_application_variable_usage.png"),
  width = 15,
  height = 8,
  units = "cm",
  dpi = 200
)
