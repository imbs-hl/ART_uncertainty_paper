##' Benchmark evaluation and visualization from:
##' Kronziel et al. (2024): "Predicting Medical Outcomes using Artificial 
##' Representative Trees with Uncertainty Quantification".
##'
##' This script reproduces Figure 6 based on benchmark experiment results.
##' It requires precomputed results from at least one dataset.
##' Run "01run_benchmark.R" beforehand to generate these results.
##'
##' Note:
##' - The original paper evaluates multiple datasets and parameter settings.
##' - To reduce runtime on local machines, this script uses a reduced subset
##'   (e.g., fewer parameter configurations and a single dataset).
##'
##' Output:
##' - Boxplots for stability measures:
##'   * Prediction distance (regression and probability)
##'   * Splitting variable distance
##' - Combined figure summarizing model stability

# WARNING:
# Computing stability across multiple benchmark datasets or many repetitions
# can be very computationally intensive. The figures in the publication were
# generated on a high-performance computing cluster.


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

# Directory containing external benchmark datasets
data_ext_dir <- file.path(main_dir, "external_data")


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

# Load stored models and prediction probabilities
regression_trees <- readRDS(file.path(proc_dir, "regression_trees_benchmark_experiments.Rds"))
probability_trees <- readRDS(file.path(proc_dir, "probability_trees_benchmark_experiments.Rds"))
pred_prob <- readRDS(file.path(proc_dir, "pred_probabilities_benchmark_experiments.Rds"))

# if you want to use the original data from the publication, please download them (see README for details), unpack them, and move them into the data folder: 
# 
# data from publication
#
# data <- readRDS(file.path(proc_dir, "results_benchmark_experiments_benchmark_results_from_paper.rds")) %>% 
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
# 
# # Load stored models and prediction probabilities
# regression_trees <- readRDS(file.path(proc_dir, "regression_trees_benchmark_experiments_benchmark_results_from_paper.Rds"))
# probability_trees <- readRDS(file.path(proc_dir, "probability_trees_benchmark_experiments_benchmark_results_from_paper.Rds"))
# pred_prob <- readRDS(file.path(proc_dir, "pred_probabilities_benchmark_experiments_benchmark_results_from_paper.Rds"))


#---------------------------------------
## Stability computation setup

set.seed(123)

# NOTE:
# For stability analysis, the full dataset is used (no train/test split),
# since stability focuses on model variability rather than predictive performance.

# Full list of datasets used in the paper (commented out for runtime reasons)
# all_dataset_name <- c("abalone", "airfoil", "deltaA", "deltaE", "friedm", "mortgage", "wineW", "wizmir",
#                      "puma8fh", "puma8fm", "puma8nh", "puma8nm",
#                      "kin8fh", "kin8fm", "kin8nh", "kin8nm",
#                      "bank8fh", "bank8fm", "bank8nh", "bank8nm", "comp")

# Reduced example: only one dataset
all_dataset_name <- "wineW"


#---------------------------------------
## Function to load benchmark datasets

get_bench_data <- function(dataset_name){
  
  # Load dataset depending on its name and rename target variable to "y"
  
  if(dataset_name == "abalone"){
    data <- read.table(file.path(data_ext_dir, "abalone/abalone.data"),
                       sep = ",", header = FALSE) %>% 
      dplyr::rename(y = V9)
    
  } else if(dataset_name == "airfoil"){
    data <- read.table(file.path(data_ext_dir, "airfoil/airfoil_self_noise.dat"),
                       sep = "\t", header = FALSE) %>% 
      dplyr::rename(y = V6)
    
  } else if(dataset_name == "wineW"){
    data <- read.csv(file.path(data_ext_dir, "wineW/winequality-white.csv"),
                     header = TRUE, sep = ";") %>% 
      dplyr::rename(y = quality)
  }
  
  return(data)
}


#---------------------------------------
## Stability computation

# Number of repetitions (reduced for runtime)
num_rep <- 3

# Initialize result data frame
dist_setting <- data.frame(method = NA, dataset_name = NA, metric = NA, 
                           min.bucket = NA, probs_quantiles = NA, 
                           rep = NA, dist_pred = NA, dist_sv = NA, 
                           dist_pred0.5 = NA, dist_sv0.5 = NA)

# Loop over methods and datasets
for(method in unique(data$method)){
  for(i in 1:length(all_dataset_name)){
    
    dataset_name = all_dataset_name[i]
    
    # Select corresponding results
    ids_dataset_name <- which(data$dataset_name == dataset_name & 
                                data$method == method)
    
    stab_data_setting <- data[ids_dataset_name,]
    trees_setting <- regression_trees[ids_dataset_name]
    probability_trees_setting <- probability_trees[ids_dataset_name]
    pred_prob_setting <- pred_prob[ids_dataset_name]
    
    # Load full dataset
    test_data_setting = get_bench_data(dataset_name)
    
    num_features <- ncol(test_data_setting) - 1
    
    # Loop over repetitions
    for(rep in 1:num_rep){
      
      # Select models for current repetition and parameter setting
      ids_setting_i <- which((stab_data_setting$metric == "splitting variables" | is.na(stab_data_setting$metric))  & 
                               (stab_data_setting$min.bucket == 150 | is.na(stab_data_setting$min.bucket)) & 
                               (stab_data_setting$probs_quantiles == "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9" | is.na(stab_data_setting$probs_quantiles)) &
                               stab_data_setting$repitition == rep)
      
      trees_setting_i <- trees_setting[ids_setting_i]
      probability_trees_setting_i <- probability_trees_setting[ids_setting_i]
      
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
        
        probs_df <- matrix(NA, nrow = nrow(test_data_setting), ncol = length(dist_pred))
        
        for(j in 1:length(dist_pred)){
          
          nodes <- predict(data = test_data_setting, 
                           trees_setting_i[[j]], 
                           type = "terminalNodes")$predictions
          
          pred_df <- data.frame(nodeID = nodes,
                                y = test_data_setting$y) %>% 
            left_join(pred_prob_i[[j]], by = "nodeID")
          
          probs_df[, j] <- pred_df$prob0.5
        }
        
        dist_pred0.5 <- colMeans(as.matrix(dist(t(probs_df), 
                                                method = "euclidian"))^2 / nrow(test_data_setting))
        
      } else {
        
        models_pred0.5 <- lapply(probability_trees_setting_i, function(x){
          predict(x, data = test_data_setting, predict.all = TRUE)$predictions[,2,1]
        }) %>% bind_cols()
        
        dist_pred0.5 <- colMeans(as.matrix(dist(t(models_pred0.5), 
                                                method = "euclidian"))^2 / nrow(test_data_setting))
      }
      
      #---------------------------------------
      # Splitting variable distance
      
      used_variables <- lapply(trees_setting_i, function(x){
        splitting_variables <- treeInfo(x)$splitvarID
        fu <- rep(0, num_features)
        fu[(splitting_variables + 1)] <- 1
        fu
      })
      
      used_variables <- bind_cols(used_variables)
      
      dist_sv <- colMeans(as.matrix(dist(t(as.matrix(used_variables)), 
                                         method = "euclidian"))^2 / num_features)
      
      dist_sv0.5 <- dist_sv
      
      #---------------------------------------
      # Store results
      
      dist_df <- data.frame(method = method, dataset_name = dataset_name, 
                            metric = "splitting variables", min.bucket = 150, 
                            probs_quantiles = "0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9", 
                            rep = rep,
                            dist_pred = mean(dist_pred), 
                            dist_sv = mean(dist_sv), 
                            dist_pred0.5 = mean(dist_pred0.5), 
                            dist_sv0.5 = mean(dist_sv0.5))
      
      dist_setting <- rbind(dist_setting, dist_df)
    }
  }
}


#---------------------------------------
## Visualization

data_dist <- dist_setting[-1,] %>% 
  select(method, dataset_name, contains("dist")) %>% 
  reshape2::melt() %>% 
  mutate(variable = case_when(variable == "dist_sv" ~ "SV dist. (regr.)",
                              variable == "dist_sv0.5" ~ "SV dist. (prob.)",
                              variable == "dist_pred" ~ "Prediction dist. (regr.)",
                              variable == "dist_pred0.5" ~ "Prediction dist. (prob.)"))

plot_dist <- ggplot(data_dist, 
                    aes(x = factor(method), y = as.numeric(value), col = method)) +
  geom_boxplot() +
  facet_grid(variable ~ dataset_name, scales = "free", switch = "y") +
  theme_bw() +
  labs(x = "", y = "", col = "") +
  theme(text = element_text(size = 15),
        legend.position = "none",
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255,150,90, maxColorValue = 255),
                                rgb(80,180,210, maxColorValue = 255)))

plot_dist


#---------------------------------------
## Save plot

ggsave(plot_dist, 
       filename = file.path(img_dir, "fig6_benchmark_stability.png"), 
       width = 70, height = 28, units = "cm", dpi = 200)