#' With this script the figure 5 from Kronziel et al. "Uncertainty 
#' quantification enhances the explainability of 
#' artificial representative trees" can be reproduced. 
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
# # Data from publication
# results <- read.csv2(file.path(proc_dir, "results_simulations.csv")) %>% 
#   filter(min.bucket == 100 & metric == "weighted splitting variables") 

results_benchmark <- read.csv2(file.path(proc_dir, "benchmark_data_results.csv")) %>% 
  filter(min.bucket == 100 & metric == "weighted splitting variables") %>%
  mutate(task_name_short = case_when(task_name == "Mercedes_Benz_Greener_Manufacturing" ~ "Mercedes Benz",
                                     task_name == "SAT11-HAND-runtime-regression" ~ "SAT11",
                                     task_name == "Allstate_Claims_Severity" ~ "Allstate",
                                     TRUE ~ task_name)) 

# Data produced by simulations.R
results <- readRDS(file.path(proc_dir, "results.Rds")) %>% 
  bind_rows() %>% 
  filter(min.bucket == 100 & metric == "weighted splitting variables")

results <- results  %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"))) 

# Data produced by benchmark_data.R
results_benchmark <- readRDS(file.path(proc_dir, "results_benchmark.Rds")) %>% 
  bind_rows() %>% 
  filter(min.bucket == 100 & metric == "weighted splitting variables") %>%
  mutate(task_name_short = case_when(task_name == "Mercedes_Benz_Greener_Manufacturing" ~ "Mercedes Benz",
                                     task_name == "SAT11-HAND-runtime-regression" ~ "SAT11",
                                     task_name == "Allstate_Claims_Severity" ~ "Allstate",
                                     TRUE ~ task_name)) 

# Combine simulated data with benchmark data
all_data <- bind_rows(results %>% mutate(task_name_short = scenario),
                      results_benchmark)

#---------------------------------------
# Plot data from figure 5 and save plot
# Plot coverage ICP, Mondrian ICP, CPS, Mondrian CPS
#---------------------------------------
# Transform data in long format
marginal_coverage_df1 <- aggregate(mean_error_icp ~ significance_level + task_name_short + min.bucket + metric, all_data, function(x) c(mean_error = mean(x))) %>% 
  mutate(method = "ICP",
         error = mean_error_icp)
marginal_coverage_df2 <- aggregate(mean_errors_mondrian_icp ~ significance_level + task_name_short + min.bucket + metric, all_data, function(x) c(mean_error = mean(x))) %>% 
  mutate(method = "Mondrian ICP",
         error = mean_errors_mondrian_icp)
marginal_coverage_df3 <- aggregate(mean_error_cps_two_tailed ~ significance_level + task_name_short + min.bucket + metric, all_data, function(x) c(mean_error = mean(x))) %>% 
  mutate(method = "CPS",
         error = mean_error_cps_two_tailed)
marginal_coverage_df4 <- aggregate(mean_error_cps_two_tailed_mondrian ~ significance_level + task_name_short + min.bucket + metric, all_data, function(x) c(mean_error = mean(x))) %>% 
  mutate(method = "Mondrian CPS",
         error = mean_error_cps_two_tailed_mondrian)

marginal_coverage_df <- bind_rows(marginal_coverage_df1, marginal_coverage_df2, marginal_coverage_df3, marginal_coverage_df4)


plot_coverage <- ggplot(marginal_coverage_df, 
                        aes(x = 100-significance_level*100, y = error, col = task_name_short))+
  geom_line()+
  facet_wrap(.~method)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Coverage",
       col = "Data")+
  theme(text = element_text(size = 15), legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18))+
  geom_abline(slope = -1, linetype = "dashed", intercept = 100)

plot_coverage


# save plot
ggsave(plot_coverage, filename = file.path(out_dir, "fig5_coverage.png"), width = 20, height = 12, units = "cm", dpi = 200)

