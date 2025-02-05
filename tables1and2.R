#' With this script the tables 1 and 2 from Kronziel et al. "Increasing the 
#' explainability of artificial representative trees through conformal 
#' prediction to quantify uncertainty" can be reproduced. 
#' Given a simulated data set. Run simulations.R to get such a data set. 
#' With the standard parameters in simulations.R, only the part for 
#' data scenario 1 is reproduced using less significance levels and repetitions 
#' than in the manusscript, as the runtime without a computing cluster 
#' would otherwise be too high.

#---------------------------------------
# Define directories
# Please define your main directory here. 
# This should be the directory you cloned the git repository into.
main_dir <- getwd()
setwd(main_dir)

# Create and define proc directory
dir.create(file.path(main_dir, "proc"), showWarnings = FALSE)
proc_dir <- file.path(main_dir, "proc")
# Create and define output directory
dir.create(file.path(main_dir, "output"), showWarnings = FALSE)
out_dir <- file.path(main_dir, "output")

#---------------------------------------
# Load libraries
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}

pacman::p_load(ggplot2)
pacman::p_load(gridExtra)
pacman::p_load(ranger)
pacman::p_load(devtools)
pacman::p_load(rpart)
pacman::p_load(dplyr)
pacman::p_load(tidyr)
pacman::p_load(reshape2)

#---------------------------------------
# Load and prepare data
# Data from publication
# min.bucket 25 and 100 were simulated and 100 was used for figures
results <- read.csv2(file.path(proc_dir, "results_simulations.csv"))

# Data produced by simulations.R
# results <- readRDS(file.path(proc_dir, "results.Rds")) %>% bind_rows()

# Change names and calculate confidence level in percent as (1 - significance level)*100 and coverage as 1 - error
results <- results  %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous \nvariables")),
         confidence_level = (1 - significance_level)*100,
         coverage_icp = 100-mean_error_calibration,
         coverage_mondrian_icp = 100-mean_error_mondrian_calibration,
         significance_level = significance_level*100)

#---------------------------------------
# Table 1 (averaged error rate for coverage for mondrian ICP)
min.bucket_mondrian_coverage_tab <- aggregate(mean_error_mondrian_calibration ~ significance_level + scenario + min.bucket_art, 
                                              results, 
                                              function(x) c(mean_error = round(mean(x),2))) %>% 
  pivot_wider(., names_from = c(scenario, min.bucket_art), values_from = c(mean_error_mondrian_calibration)) 

min.bucket_mondrian_coverage_tab <- min.bucket_mondrian_coverage_tab %>% 
  select(colnames(min.bucket_mondrian_coverage_tab)[order(colnames(min.bucket_mondrian_coverage_tab))]) %>% 
  select(significance_level, everything()) %>% 
  arrange(significance_level)

# Please change column names manually afterwards in csv
write.csv2(min.bucket_mondrian_coverage_tab, file.path(out_dir, "tab_simulation_mondrian_coverage.csv"), row.names = F)

#---------------------------------------
# Table 2 (averaged interval width for coverage for mondrian ICP)
min.bucket_mondrian_interval_width_tab <- aggregate(mean_interval_size_mondrian_calibration ~ significance_level + scenario + min.bucket_art, results, function(x) c(mean_error = round(mean(x),2))) %>% 
  pivot_wider(., names_from = c(scenario, min.bucket_art), values_from = c(mean_interval_size_mondrian_calibration))

min.bucket_mondrian_interval_width_tab <- min.bucket_mondrian_interval_width_tab %>% 
  select(colnames(min.bucket_mondrian_interval_width_tab)[order(colnames(min.bucket_mondrian_interval_width_tab))]) %>% 
  select(significance_level, everything()) %>% 
  arrange(significance_level)

# Please change column names manually afterwards in csv
write.csv2(min.bucket_mondrian_interval_width_tab, file.path(out_dir, "tab_simulation_mondrian_interval_width.csv"), row.names = F)


