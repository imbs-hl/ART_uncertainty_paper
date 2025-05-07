#' With this script the figure 3 from Kronziel et al. "Uncertainty 
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
# min.bucket = 100 and metric = "weighted splitting variables" was used for this figure
results <- read.csv2(file.path(proc_dir, "results_simulations.csv")) %>% 
  filter(min.bucket == 100 & metric == "weighted splitting variables") 

# Data produced by simulations.R
 results <- readRDS(file.path(proc_dir, "results.Rds")) %>% 
   bind_rows() %>% 
   filter(min.bucket == 100 & metric == "weighted splitting variables") 

# Change names and calculate confidence level in percent as (1 - significance level)*100 and coverage as 1 - error
results <- results  %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"))) 

#---------------------------------------
# Plot data from figure 3 and save plot
# Plot interval width ICP, Mondrian ICP, CPS, Mondrian CPS
#---------------------------------------
# Transform data in long format
width_df1 <- aggregate(mean_interval_size_mondrian_icp ~ significance_level + scenario + min.bucket, results, function(x) c(mean_error = mean(x)))
width_df2 <- aggregate(mean_interval_size_icp ~ significance_level + scenario + min.bucket, results, function(x) c(mean_error = mean(x)))
width_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ significance_level + scenario + min.bucket, results, function(x) c(mean_error = mean(x)))
width_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ significance_level + scenario + min.bucket, results, function(x) c(mean_error = mean(x)))

plot_data_width <- left_join(width_df1, width_df2) %>% 
  left_join(width_df3) %>% 
  left_join(width_df4) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "min.bucket")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS")) %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))

plot_width = ggplot(plot_data_width %>% filter(min.bucket == 100), 
                    aes(x=significance_level*100, y=value, col = variable, group = variable))+
  facet_wrap(scenario~.)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "Method")+
  theme(text = element_text(size = 15), #legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18))


# save plot
ggsave(plot_width, filename = file.path(out_dir, "fig3_simulation_intervalwidth.png"), width = 12, height = 12, units = "cm", dpi = 200)

