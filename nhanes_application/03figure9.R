##' With this script the figure 9 from Kronziel et al. "Prediction Beyond Point 
##' Estimates: Artificial Trees with Uncertainty" can be reproduced. 
##' Given some results from the nhanes application.
##' Run 02calculate_results.R to get such a data set and the used models. 
##' You also need the prepared nhanes data set from 01prepare_data.R here. Please run this script first.
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
  warning("Please check, if timbR with at least version 3.3 is installed.")
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
# Choose ONE of the following result sources:

# (1)
# Results produced manually via `02calculate_results.R`
# Due to runtime reasons probs_quantiles = c(0.25,0.5,0.75) is used in default setting, you may change this
data_imp <- readRDS(file.path("data", "results_nhanes_application.rds"))
# Lists of fitted regression trees
regression_trees <- readRDS(file.path(proc_dir, "regression_trees_nhanes_application.rds"))
probability_trees5.7 <- readRDS(file.path(proc_dir, "probability_trees5.7_nhanes_application.rds"))
probability_trees6.5 <- readRDS(file.path(proc_dir, "probability_trees6.5_nhanes_application.rds"))

# Node-level probability tables from RF for methods with CPS
pred_prob <- readRDS(file.path(proc_dir, "pred_probabilities_nhanes_application.rds"))

# (2)
# # Results used in the paper
# data_imp <- read.csv(file.path("data", "results_nhanes_application_results_from_paper.csv"))
# 
# # Lists of fitted regression trees
# regression_trees <- readRDS(file.path(proc_dir, "regression_trees_nhanes_application_results_from_paper.rds"))
# probability_trees5.7 <- readRDS(file.path(proc_dir, "probability_trees5.7_nhanes_application_results_from_paper.rds"))
# probability_trees6.5 <- readRDS(file.path(proc_dir, "probability_trees6.5_nhanes_application_results_from_paper.rds"))
# 
# # Node-level probability tables from RF for methods with CPS
# pred_prob <- readRDS(file.path(proc_dir, "pred_probabilities_nhanes_application_results_from_paper.rds"))




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

#------------------------------------------------------------------------------
# Load NHANES data used for stability evaluation
#------------------------------------------------------------------------------
# Predictions are computed on the full dataset to ensure fairness,
# since all models were trained on equally sized training samples.
nhanes_imp <- read.csv(file.path("data", "nhanes_prepandemic_complete.csv"))
nhanes_data <- nhanes_imp %>%
  select(-SEQN, -prediabetes)

num_features  <- ncol(nhanes_data) - 1
#------------------------------------------------------------------------------
# Stability analysis: glycohemoglobin (regression outcome)
#------------------------------------------------------------------------------
set.seed(123)

# NOTE:
# For stability analysis, the full dataset is used (no train/test split),
# since stability focuses on model variability rather than predictive performance.

#---------------------------------------
## Stability computation

# Maximum number of repetitions used in 02calculate_results.R (reduced for runtime)
num_rep <- max(data$repitition)
# hyperparemeter of ART
probs_quantiles <- "0.25,0.5,0.75" # in Paper no quantiles were used
min.bucket <- 150

# Initialize result data frame
dist_setting <- data.frame(method = NA, metric = NA, 
                           min.bucket = NA, probs_quantiles = NA, rep = NA, 
                           dist_pred = NA, dist_sv = NA, 
                           dist_pred5.7 = NA, dist_sv5.7 = NA, 
                           dist_pred6.5 = NA, dist_sv6.5 = NA)

# Loop over methods and datasets
for(method in c("ART + CPS", "DT + CPS", "Mult. DTs", "Mult. ARTs")){
    
    # Select corresponding results
    ids_dataset_name <- which(data$method == method)
    
    stab_data_setting <- data[ids_dataset_name,]
    trees_setting <- regression_trees[ids_dataset_name]
    probability_trees5.7_setting <- probability_trees5.7[ids_dataset_name]
    probability_trees6.5_setting <- probability_trees6.5[ids_dataset_name]
    pred_prob_setting <- pred_prob[ids_dataset_name]
    
    # Load full dataset
    test_data_setting = nhanes_data
    
    num_features <- ncol(test_data_setting) - 1
    
    # Loop over repetitions
    for(rep in 1:num_rep){
      
      # Select models for current repetition and parameter setting
      ids_setting_i <- which((stab_data_setting$metric == "splitting variables" | is.na(stab_data_setting$metric))  & 
                               (stab_data_setting$min.bucket == min.bucket | is.na(stab_data_setting$min.bucket)) & 
                               (stab_data_setting$probs_quantiles == probs_quantiles | is.na(stab_data_setting$probs_quantiles)) &
                               stab_data_setting$repitition == rep)
      
      trees_setting_i <- trees_setting[ids_setting_i]
      trees5.7_setting_i <- probability_trees5.7_setting[ids_setting_i]
      trees6.5_setting_i <- probability_trees6.5_setting[ids_setting_i]
      
      #---------------------------------------
      # Prediction distance (regression)
      
      models_pred <- lapply(trees_setting_i, function(x){
        predict(x, data = test_data_setting, predict.all = TRUE)$predictions
      })
      
      models_pred <- bind_cols(models_pred)
      
      dist_pred <- colMeans(as.matrix(dist(t(models_pred), 
                                           method = "euclidian"))^2 / nrow(test_data_setting))
      
      #---------------------------------------
      # Prediction distance (probability)
      
      if(grepl("CPS", method)){
        
        pred_prob_i <- pred_prob_setting[ids_setting_i]
        
        probs_df5.7 <- matrix(NA, nrow = nrow(test_data_setting), ncol = length(dist_pred))
        probs_df6.5 <- matrix(NA, nrow = nrow(test_data_setting), ncol = length(dist_pred))
        
        for(j in 1:length(dist_pred)){
          
          nodes <- predict(data = test_data_setting, 
                           trees_setting_i[[j]], 
                           type = "terminalNodes")$predictions
          
          pred_df <- data.frame(nodeID = nodes[,1],
                                y = test_data_setting$glycohemoglobin) %>% 
            left_join(pred_prob_i[[j]], by = "nodeID")
          
          probs_df5.7[, j] <- pred_df$prob5.7
          probs_df6.5[, j] <- pred_df$prob6.5
        }
        
        dist_pred5.7 <- colMeans(as.matrix(dist(t(probs_df5.7), method = "euclidian"))^2 / nrow(test_data_setting))
        dist_pred6.5 <- colMeans(as.matrix(dist(t(probs_df6.5), method = "euclidian"))^2 / nrow(test_data_setting))
        
      } else {
        
        models_pred5.7 <- lapply(trees5.7_setting_i, function(x){
          predict(x, data = test_data_setting, predict.all = TRUE)$predictions[,2,1]
        }) %>% bind_cols()
        models_pred6.5 <- lapply(trees6.5_setting_i, function(x){
          predict(x, data = test_data_setting, predict.all = TRUE)$predictions[,2,1]
        }) %>% bind_cols()
        
        dist_pred5.7 <- colMeans(as.matrix(dist(t(probs_df5.7), method = "euclidian"))^2 / nrow(test_data_setting))
        dist_pred6.5 <- colMeans(as.matrix(dist(t(probs_df6.5), method = "euclidian"))^2 / nrow(test_data_setting))
      }
      
      
      #---------------------------------------
      # Splitting variable distance
      
      # Regression tree
      used_variables <- lapply(trees_setting_i, function(x){
        splitting_variables <- treeInfo(x)$splitvarID
        fu <- rep(0, num_features)
        fu[(splitting_variables + 1)] <- 1
        fu
      })
      
      used_variables <- bind_cols(used_variables)
      
      dist_sv <- colMeans(as.matrix(dist(t(as.matrix(used_variables)), 
                                         method = "euclidian"))^2 / num_features)
      
      if(grepl("CPS", method)){
        dist_sv5.7 <- dist_sv
        dist_sv6.5 <- dist_sv
      }else{
        # Probability trees
        used_variables5.7 <- lapply(trees5.7_setting_i, function(x){
          splitting_variables <- treeInfo(x)$splitvarID
          fu <- rep(0, num_features)
          fu[(splitting_variables + 1)] <- 1
          fu
        })
        used_variables5.7 <- bind_cols(used_variables5.7)
        
        dist_sv5.7 <- colMeans(as.matrix(dist(t(as.matrix(used_variables5.7)), method = "euclidian"))^2 / num_features)
        # Probability trees
        used_variables6.5 <- lapply(trees6.5_setting_i, function(x){
          splitting_variables <- treeInfo(x)$splitvarID
          fu <- rep(0, num_features)
          fu[(splitting_variables + 1)] <- 1
          fu
        })
        used_variables6.5 <- bind_cols(used_variables6.5)
        
        dist_sv6.5 <- colMeans(as.matrix(dist(t(as.matrix(used_variables6.5)), method = "euclidian"))^2 / num_features)
      }
      

      
      #---------------------------------------
      # Store results
      dist_df <- data.frame(method = method,
                            metric = "splitting variables", min.bucket = min.bucket, 
                            probs_quantiles = probs_quantiles, 
                            rep = rep,
                            dist_pred = mean(dist_pred), 
                            dist_sv = mean(dist_sv), 
                            dist_pred5.7 = mean(dist_pred5.7), 
                            dist_sv5.7 = mean(dist_sv5.7), 
                            dist_pred6.5 = mean(dist_pred6.5), 
                            dist_sv6.5 = mean(dist_sv6.5))
      
      dist_setting <- rbind(dist_setting, dist_df)
    
  }
}






#------------------------------------------------------------------------------
# Plot stability results
#------------------------------------------------------------------------------
stability_df <- dist_setting[-1,] %>% 
  select(-metric, -min.bucket, -probs_quantiles, -rep) %>% 
  reshape2::melt() %>% 
  mutate(outcome = case_when(variable == "dist_sv" ~ "glycohemoglobin",
                             variable == "dist_sv5.7" ~ "prediabetes",
                             variable == "dist_sv6.5" ~ "diabetes",
                             variable == "dist_pred" ~ "glycohemoglobin",
                             variable == "dist_pred5.7" ~ "prediabetes",
                             variable == "dist_pred6.5" ~ "diabetes")) %>% 
  mutate(distance = case_when(grepl("_sv", variable)~"SV distance",
                              TRUE ~ "Prediction distance"))

plot_stability_pred <- ggplot(data = stability_df %>% filter(distance == "Prediction distance"),
                         aes(x=method, y = value, col = method)) +
  geom_boxplot()+
  facet_grid(
    . ~ factor(
      outcome,
      labels = c("glycohemoglobin", "prediabetes", "diabetes"),
      levels = c("glycohemoglobin", "prediabetes", "diabetes")
    )) +
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

plot_stability_pred


plot_stability_sv <- ggplot(data = stability_df %>% filter(distance == "SV distance"),
                              aes(x=method, y = value, col = method)) +
  geom_boxplot()+
  facet_grid(
    . ~ factor(
      outcome,
      labels = c("glycohemoglobin", "prediabetes", "diabetes"),
      levels = c("glycohemoglobin", "prediabetes", "diabetes")
    )) +
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

plot_stability_sv




plot_dist <- plot_grid(plot_stability_pred, plot_stability_sv, labels = "AUTO", nrow = 2)
plot_dist


#------------------------------------------------------------------------------
# Save figure
#------------------------------------------------------------------------------
ggsave(
  plot_dist,
  filename = file.path(img_dir, "fig9_nhanes_application_stability.png"),
  width = 33,
  height = 14,
  units = "cm",
  dpi = 200
)
