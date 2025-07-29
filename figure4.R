#' With this script the figure 4 from "Uncertainty 
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
pacman::p_load(cowplot)

#---------------------------------------
# Load and prepare data
# # Data from publication
# results <- read.csv2(file.path(proc_dir, "results_simulations.csv"))

# Data produced by simulations.R
results <- readRDS(file.path(proc_dir, "results.Rds")) %>%
   bind_rows()

# Change names 
results <- results  %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")),
         metric = ifelse(metric == "weighted splitting variables", "WSV", metric))

#---------------------------------------
# Plot data from figure 4A 
# Plot interval width for different min.bucket
#---------------------------------------
# Transform data in long format
width_df1 <- aggregate(mean_interval_size_mondrian_icp ~ significance_level + scenario + min.bucket, results %>% filter(metric == "WSV"), function(x) c(mean_error = mean(x)))
width_df2 <- aggregate(mean_interval_size_icp ~ significance_level + scenario + min.bucket, results %>% filter(metric == "WSV"), function(x) c(mean_error = mean(x)))
width_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ significance_level + scenario + min.bucket, results %>% filter(metric == "WSV"), function(x) c(mean_error = mean(x)))
width_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ significance_level + scenario + min.bucket, results %>% filter(metric == "WSV"), function(x) c(mean_error = mean(x)))

plot_data_width <- full_join(width_df1, width_df2) %>% 
  full_join(width_df3) %>% 
  full_join(width_df4) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "min.bucket")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS")) %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))


# Only scenario 1
# if you want to use other scenarios please adjust "filter(scenario == "large effects")" in the line below
plot_minbucket_width_scenario1 <- ggplot(plot_data_width %>% filter(scenario == "large effects"), 
                                         aes(x=(significance_level*100), y=value, col = factor(min.bucket), group = min.bucket))+
  facet_wrap(.~variable, ncol = 2)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "min.bucket") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 16),
        strip.text.x = element_text(size = 16))

#---------------------------------------
# Plot data from figure 4B
# Plot interval width for different metrics
#---------------------------------------
data_width_metric <- results %>% 
  filter(min.bucket == 100) 
width_metric_df1 <- aggregate(mean_interval_size_mondrian_icp ~ significance_level + scenario + metric, data_width_metric, function(x) c(mean_error = mean(x)))
width_metric_df2 <- aggregate(mean_interval_size_icp ~ significance_level + scenario + metric, data_width_metric, function(x) c(mean_error = mean(x)))
width_metric_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ significance_level + scenario + metric, data_width_metric, function(x) c(mean_error = mean(x)))
width_metric_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ significance_level + scenario + metric, data_width_metric, function(x) c(mean_error = mean(x)))

plot_data_metric_width <- full_join(width_metric_df1, width_metric_df2) %>% 
  full_join(width_metric_df3) %>% 
  full_join(width_metric_df4) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS"))

# Only scenario 1
# if you want to use other scenarios please adjust "filter(scenario == "large effects")" in the line below
plot_metric_width_scenario1 <- ggplot(plot_data_metric_width %>% filter(scenario == "large effects"), 
                                      aes(x=(significance_level*100), y=value, col = factor(metric), group = metric))+
  facet_wrap(.~variable, ncol = 2)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "metric") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18))


#---------------------------------------
# Save plot
#---------------------------------------
plots_hyperparameters <- plot_grid(plot_minbucket_width_scenario1, plot_metric_width_scenario1, 
                                   ncol = 1,# align = "hv",
                                   labels = c("A", "B"))

plots_hyperparameters

ggsave(plots_hyperparameters, filename = file.path(out_dir, "fig4_simulation_intervalwidth_hyperparameter_largeeffects.png"), 
       width = 20, height = 20, units = "cm", dpi = 200)