#' With this script the figure 2 from Kronziel et al. "Increasing the 
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
results <- read.csv2(file.path(proc_dir, "results_simulations.csv")) %>% 
  filter(min.bucket_art == 100) 

# Data produced by simulations.R
 results <- readRDS(file.path(proc_dir, "results.Rds")) %>% bind_rows()

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
         coverage_mondrian_icp = 100-mean_error_mondrian_calibration)

#---------------------------------------
# Plot data from figure 2A and save plot

# Aggregate meausures from repetitions using mean
marginal_coverage_mondrian_icp <- aggregate(coverage_mondrian_icp ~ confidence_level + scenario + min.bucket_art, results, function(x) c(mean_error = mean(x)))
marginal_coverage_icp <- aggregate(coverage_icp ~ confidence_level + scenario + min.bucket_art, results, function(x) c(mean_error = mean(x)))

# Transform data in long format
plot_figure2a <- left_join(marginal_coverage_mondrian_icp, marginal_coverage_icp) %>% 
  melt(id=c("scenario", "confidence_level", "min.bucket_art")) %>% 
  mutate(variable = ifelse(variable == "coverage_icp", "ICP", "Mondrian ICP"))


ggplot(plot_figure2a, aes(x = confidence_level, y = value, col = scenario, group = scenario))+
  facet_grid(variable~.)+
  theme_bw()+
  geom_line()+
  geom_abline(slope = 1, linetype = "dashed", intercept = 0)+
  labs(x = "Confidence level (%)",
       y = "Marginal Coverage (%)",
       col = "Simulation scenarios") +
  theme(text = element_text(size = 15), legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18))

# save plot
ggsave(file.path(out_dir, "fig2A_averaged_coverage.png"), units = "cm", width = 20, height = 8)


#---------------------------------------
# Plot data from figure 2B and save plot

# Aggregate meausures from repetitions using mean
interval_width_mondrian_icp <- aggregate(mean_interval_size_mondrian_calibration ~ confidence_level + scenario + min.bucket_art, results, function(x) c(mean_error = mean(x)))
interval_width_icp <- aggregate(mean_interval_size_calibration ~ confidence_level + scenario + min.bucket_art, results, function(x) c(mean_error = mean(x)))

# Transform data in long format
plot_figure2b <- left_join(interval_width_mondrian_icp, interval_width_icp) %>% 
  melt(id=c("scenario", "confidence_level", "min.bucket_art")) %>% 
  mutate(variable = ifelse(variable == "mean_interval_size_calibration", "ICP", "Mondrian ICP"))


ggplot(plot_figure2b, aes(x=confidence_level, y=value, col = scenario, group = scenario))+
  facet_grid(variable~.)+
  theme_bw()+
  geom_line()+
  labs(x = "Confidence level (%)",
       y = "Interval width",
       col = "Simulation scenarios") +
  theme(text = element_text(size = 15), legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18))

# save plot
ggsave(file.path(out_dir, "fig2B_averaged_intervalwidth.png"), units = "cm", width = 20, height = 8)




