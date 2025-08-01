#' With this script the figure 7 from "Uncertainty 
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
# Data from publication

# results <- readRDS(file.path(proc_dir, "benchmark_data_results.rds")) %>% 
#   mutate(significance_level2 = (significance_level*100),
#          significance_level3 = paste0(significance_level2, "%"),
#          task_name_short = case_when(task_name == "Mercedes_Benz_Greener_Manufacturing" ~ "Mercedes Benz",
#                                      task_name == "SAT11-HAND-runtime-regression" ~ "SAT11",
#                                      task_name == "Allstate_Claims_Severity" ~ "Allstate",
#                                      TRUE ~ task_name)) 

# Data from benchmark_data.R
results <- readRDS(file.path(proc_dir, "results_benchmark.Rds")) %>%
  bind_rows() %>%
  mutate(task_name_short = case_when(task_name == "Mercedes_Benz_Greener_Manufacturing" ~ "Mercedes Benz",
                                     task_name == "SAT11-HAND-runtime-regression" ~ "SAT11",
                                     task_name == "Allstate_Claims_Severity" ~ "Allstate",
                                     TRUE ~ task_name))


#---------------------------------------
# Plot data from figure 6 and save plot
# Plot interval width for ICP, Mondrian ICP, CPS, Mondrian CPS
#---------------------------------------
# Aggregate and transform data in long format
interval_width_df1 <- aggregate(mean_interval_size_mondrian_icp ~ significance_level + task_name_short + min.bucket + metric, results, function(x) c(mean_error = mean(x)))
interval_width_df2 <- aggregate(mean_interval_size_icp ~ significance_level + task_name_short + min.bucket + metric, results, function(x) c(mean_error = mean(x)))
interval_width_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ significance_level + task_name_short + min.bucket + metric, results, function(x) c(mean_error = mean(x)))
interval_width_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ significance_level + task_name_short + min.bucket + metric, results, function(x) c(mean_error = mean(x)))

plot_data <- left_join(interval_width_df1, interval_width_df2) %>% 
  left_join(interval_width_df3) %>% 
  left_join(interval_width_df4) %>% 
  reshape2::melt(id=c("task_name_short", "significance_level", "min.bucket", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS"))

# Build plot
plot_width = ggplot(plot_data %>% filter(min.bucket == 100 & metric == "weighted splitting variables"), 
                     aes(x=(significance_level*100), y=value, col = variable, group = variable))+
  facet_wrap(task_name_short~., ncol = 5)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "") +
  theme(text = element_text(size = 15), legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))

plot_width

# Insert this if you use more than one benchmark data set
# plot_width2 <- shift_legend2(plot_width)

# Save plot
ggsave(plot_width, filename = file.path(out_dir, "fig7_benchmark_interval_width.png"), width = 30, height = 15, units = "cm", dpi = 200)

