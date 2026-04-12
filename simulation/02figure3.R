##' With this script the figure 3 from Kronziel et al. "Prediction Beyond Point 
##' Estimates: Artificial Trees with Uncertainty" can be reproduced. 
##' Given a simulated data set. Run 01run_simulations.R to get such a data set. 
##' With the standard parameters in 01run_simulations.R, only the part for 
##' data scenario 1 is reproduced, as the runtime without a computing cluster 
##' would otherwise be too high.


## Load libraries -----------------------------------------------------------
# Install pacman if not already installed (used for convenient package loading)
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}
library(pacman)

# List of required packages for simulation, modeling, and plotting
packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "DescTools", "caret", "stringr", "devtools", "this.path")

# Load (and install if necessary) all required packages
p_load(packages, character.only = TRUE)

# Load timbR package (from GitHub if not installed)
if("timbR" %in% installed.packages()){
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}


#---------------------------------------
## Define directories -------------------------------------------------------
# Define main directory (root of the cloned repository)
main_dir <- this.dir()
setwd(main_dir)

# Directory containing processed data
proc_dir <- file.path(main_dir, "data")

# Create directory for saving plots (if it does not exist)
dir.create(file.path(main_dir, "img"), showWarnings = FALSE)
img_dir <- file.path(main_dir, "img")


#---------------------------------------
## Load and prepare data
# ATTENTION!!!!
warning("Please choose here which data you want to use. The default is the results you got from running 01run_simulations.R")
# if you want to use your simulated results you got from running 01run_simulations.R, use the code at (1) for data preparation
# if you want to use the original data from the publication, please download them (see README for details), unpack them, and move them into the data folder. Use the code at (2) for data preparation

# (1) Results from 01run_simulations.R
# Load simulation results
results <- readRDS(file.path(proc_dir, "results_simulated_results.rds"))  %>%
  # Create readable scenario labels
  mutate(scenario2 = case_when(setting == "Setting 1" ~ "large effects",
                               setting == "Setting 2" ~ "small effects",
                               setting == "Setting 3" ~ "correlations",
                               setting == "Setting 4" ~ "interactions",
                               setting == "Setting 5" ~ "continuous variables"),
         scenario2 = factor(scenario2,
                            levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"),
                            labels = c("1: large effects", "2: small effects", "3: correlations", "4: interactions", "5: continuous \nvariables"))) %>%

  # Harmonize method naming for plotting
  mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
                            method == "Regression DT + CPS" ~ "DT + CPS",
                            method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
                            method == "Regression DT + Probability DTs" ~ "Mult. DTs",
                            TRUE ~ method))

# Load stored trees and predicted probabilities
regression_trees_imp <- readRDS(file.path("data","regression_trees_simulated_results.rds"))
probability_trees_imp <- readRDS(file.path("data","probability_trees_simulated_results.rds"))
pred_prob_imp <- readRDS(file.path("data", "pred_probabilities_simulated_results.rds"))


# 
# # (2) Results from publication
# results <- readRDS(file.path(proc_dir, "results_simulated_results_from_paper.rds")) %>% 
#   # Create readable scenario labels
#   mutate(scenario2 = case_when(setting == "Setting 1" ~ "large effects",
#                                setting == "Setting 2" ~ "small effects",
#                                setting == "Setting 3" ~ "correlations",
#                                setting == "Setting 4" ~ "interactions",
#                                setting == "Setting 5" ~ "continuous variables"),
#          scenario2 = factor(scenario2,
#                             levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"),
#                             labels = c("1: large effects", "2: small effects", "3: correlations", "4: interactions", "5: continuous \nvariables"))) %>%
# 
#   # Harmonize method naming for plotting
#   mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
#                             method == "Regression DT + CPS" ~ "DT + CPS",
#                             method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
#                             method == "Regression DT + Probability DTs" ~ "Mult. DTs",
#                             TRUE ~ method))
# 
# # Load stored trees and predicted probabilities
# regression_trees_imp <- readRDS(file.path("data","regression_trees_simulated_results_from_paper.rds"))
# probability_trees_imp <- readRDS(file.path("data","probability_trees_simulated_results_from_paper.rds"))
# pred_prob_imp <- readRDS(file.path("data", "pred_probabilities_simulated_results_from_paper.rds"))



#----------------
# Stability Analysis --------------------------------------------------------
#----------------

# Load data-generating functions for all simulation settings
source(file.path("functions/simulate_rf_setting_1.R"))
source(file.path("functions/simulate_rf_setting_2.R"))
source(file.path("functions/simulate_rf_setting_3.R"))
source(file.path("functions/simulate_rf_setting_4.R"))
source(file.path("functions/simulate_rf_setting_5.R"))

# Fix random seed for reproducibility
set.seed(123)

# Generate independent test datasets for each setting
test_data_setting1 <- simulate_rf_setting_1(data = 1000,
                                            p_eff       = 5,
                                            beta_eff    = 2,
                                            n_test      = 1000,
                                            n_val       = 1000, 
                                            p           = 100, 
                                            num.trees   = 500,
                                            eps         = 1,
                                            mtry        = 10,
                                            min_node_size = 100)[[2]]
test_data_setting2 <- simulate_rf_setting_2(data = 1000,
                                            p_eff       = 50,
                                            beta_eff    = 0.2,
                                            n_test      = 1000,
                                            n_val       = 1000, 
                                            p           = 100, 
                                            num.trees   = 500,
                                            eps         = 1,
                                            mtry        = 10,
                                            min_node_size = 100)[[2]]
test_data_setting3 <- simulate_rf_setting_3(data = 1000,
                                            p_eff       = 5,
                                            beta_eff    = 2,
                                            n_test      = 1000,
                                            n_val       = 1000,
                                            p           = 100,
                                            p_corr      = 5,
                                            n_blocks    = 5,
                                            cor         = 0.3,
                                            num.trees   = 500,
                                            eps         = 1,
                                            mtry        = 10,
                                            min_node_size = 100)[[2]]

test_data_setting4 <- simulate_rf_setting_4(data = 1000,
                                            p_eff       = 5,
                                            beta_eff    = 2,
                                            n_test      = 1000,
                                            n_val       = 1000,
                                            p           = 100,
                                            p_int       = 5,
                                            beta_int    = 2,
                                            num.trees   = 500,
                                            eps         = 1,
                                            mtry        = 10,
                                            min_node_size = 100)[[2]]

test_data_setting5 <- simulate_rf_setting_5(data = 1000,
                                            p_eff_bin   = 5,
                                            p_eff_con   = 5,
                                            beta_eff    = 2,
                                            n_test      = 1000,
                                            n_val       = 1000,
                                            p           = 100,
                                            num.trees   = 500,
                                            eps         = 1,
                                            mtry        = 10,
                                            min_node_size = 100)[[2]]


# Extract number of settings 
num_settings <- str_extract_all(unique(results$setting), "\\d") %>% as.numeric()

# Sanity checks to ensure correct parameter subset is used
if(all(!grepl("150", results$min.bucket))){
  stop("You need results with min.bucket = 150 for the original figure from the publication.")
}
if(all(!grepl("splitting variables", results$metric))){
  stop("You need results with metric = 'splitting variables' for the original figure from the publication.")
}


# Initialize dataframe to store stability metrics
dist_setting <- data.frame(method = NA, setting = NA, metric = NA, min.bucket = NA, probs_quantiles = NA, 
                           dist_pred = NA, dist_sv = NA, dist_pred0.5 = NA, dist_sv0.5 = NA)

# Loop over methods and simulation settings
for(method in unique(results$method)){
  for(i in num_settings){
    
    # Select results corresponding to current setting and method
    setting = paste0("Setting ", i)
    ids_setting <- which(results$setting== setting & results$method == method)
    
    stab_data_setting <- results[ids_setting,]
    regression_trees <- regression_trees_imp[ids_setting]
    probability_trees <- probability_trees_imp[ids_setting]
    pred_prob_setting <- pred_prob_imp[ids_setting]
    
    # Select corresponding test dataset
    test_data_setting = case_when(i == 1 ~ test_data_setting1,
                                  i == 2 ~ test_data_setting2,
                                  i == 3 ~ test_data_setting3,
                                  i == 4 ~ test_data_setting4,
                                  i == 5 ~ test_data_setting5)
    
    # Fix parameter configuration (as in paper)
    num_features <- ncol(test_data_setting)-1
    metric <- "splitting variables"
    min.bucket <- 150
    quantiles <- ""
    
    # Filter models matching configuration
    ids_setting_i <- which((stab_data_setting$metric==metric | is.na(stab_data_setting$metric))  & 
                             (stab_data_setting$min.bucket==min.bucket | is.na(stab_data_setting$min.bucket)) & 
                             (stab_data_setting$probs_quantiles==quantiles | is.na(stab_data_setting$probs_quantiles)))
    
    num_rep <- length(ids_setting_i)
    
    regression_trees_i <- regression_trees[ids_setting_i]
    probability_trees_setting_i <- probability_trees[ids_setting_i]
    
    
    # --- Prediction distance (continuous outcome) -------------------------
    models_pred_arts_cps <- lapply(regression_trees_i, function(x){
      predict(x, data = test_data_setting, predict.all = TRUE)$predictions
    })
    models_pred_arts_cps <- bind_cols(models_pred_arts_cps)
    
    # Compute pairwise Euclidean distance between model predictions
    dist_pred <- colMeans(as.matrix(dist(t(models_pred_arts_cps), 
                                         method = "euclidian"))^2/nrow(test_data_setting))
    
    
    # --- Prediction distance (binary outcome) -----------------------------
    if(grepl("CPS", method)){
      pred_prob_i <- pred_prob_setting[ids_setting_i]
      
      probs_df0.5 <- matrix(NA, nrow = nrow(test_data_setting), ncol = num_rep)
      for(j in 1:num_rep){
        nodes_test_data <- predict(data = test_data_setting, 
                                   regression_trees_i[[j]], type = "terminalNodes")$predictions
        
        # Map node-wise probabilities to observations
        pred_df <- data.frame(nodeID = nodes_test_data,
                              y = test_data_setting$y) %>% 
          left_join(pred_prob_i[[j]], by = "nodeID")
        
        probs_df0.5[,j] <- pred_df$prob0.5
      }
      
      dist_pred0.5 <- colMeans(as.matrix(dist(t(probs_df0.5), method = "euclidian"))^2/nrow(test_data_setting))
      
    } else {
      # Direct probability predictions from probability trees
      models_pred0.5 <- lapply(probability_trees_setting_i, function(x){
        predict(x, data = test_data_setting)$predictions[,2]
      }) %>% bind_cols()
      
      dist_pred0.5 <- colMeans(as.matrix(dist(t(models_pred0.5), 
                                              method = "euclidian"))^2/nrow(test_data_setting))
    }
    
    
    # --- Variable usage distance -----------------------------------------
    # Compute binary vector indicating whether a variable is used in splits
    
    used_variables <- lapply(regression_trees_i, function(x){
      splitting_variables <- treeInfo(x)$splitvarID
      fu <- rep(0, num_features)
      fu[(splitting_variables + 1)] <- 1
      fu
    })
    used_variables <- bind_cols(used_variables)
    
    dist_sv <- colMeans(as.matrix(dist(t(as.matrix(used_variables)), method = "euclidian"))^2/num_features)
    
    if(grepl("CPS", method)){
      dist_sv0.5 <- dist_sv
    } else {
      used_variables0.5 <- lapply(probability_trees_setting_i, function(x){
        splitting_variables <- treeInfo(x)$splitvarID
        fu <- rep(0, num_features)
        fu[(splitting_variables + 1)] <- 1
        fu
      })
      used_variables0.5 <- bind_cols(used_variables0.5)
      
      dist_sv0.5 <- colMeans(as.matrix(dist(t(as.matrix(used_variables0.5)), method = "euclidian"))^2/num_features)
    }
    
    
    # Store results
    dist_df <- data.frame(method = method, setting = setting, metric = metric, min.bucket = min.bucket, probs_quantiles = quantiles, 
                          dist_pred = dist_pred, dist_sv = dist_sv, dist_pred0.5 = dist_pred0.5, dist_sv0.5 = dist_sv0.5)
    
    colnames(dist_df) <- c("method", "setting", "metric", "min.bucket", "probs_quantiles", "dist_pred", "dist_sv", "dist_pred0.5", "dist_sv0.5")
    
    dist_setting <- rbind(dist_setting, dist_df)
  }
}


#---------------------------------------
## Prepare data for plotting -----------------------------------------------

data_dist <- dist_setting[-1,] %>% 
  dplyr::select(method, setting, contains("dist")) %>% 
  reshape2::melt() %>% 
  
  # Rename metrics for visualization
  mutate(variable = case_when(variable == "dist_sv"~"SV dist. regr.",
                              variable == "dist_sv0.5"~"SV dist. prob.",
                              variable == "dist_pred"~"pred. dist. regr.",
                              variable == "dist_pred0.5"~"pred. dist. prob.",
  )) %>% 
  
  # Add scenario labels
  left_join(results %>% dplyr::select(setting, scenario2) %>% unique())


#---------------------------------------
## Plot Figure 3 (stability comparison) ------------------------------------

plot_dist <- ggplot(data_dist, 
                    aes(x=factor(method), y= as.numeric(value), col = method)) +
  geom_boxplot()+
  facet_grid(variable~scenario2 ,scales="free", switch="y") +
  theme_bw() +
  labs(x = "", y = "", col ="")+
  theme(text = element_text(size = 15), legend.position = "none",
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  # Custom color palette for methods
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 

plot_dist


# Save figure to disk
ggsave(plot_dist, filename = file.path(img_dir, "fig3_simulation_stability.png"),
       width = 36, height = 25, units = "cm", dpi = 200)