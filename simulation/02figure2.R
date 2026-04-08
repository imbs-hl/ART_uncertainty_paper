##' With this script the figure 2 from Kronziel et al. "Predicting Medical 
##' Outcomes using Artificial Representative Trees with Uncertainty 
##' Quantification" can be reproduced. Given a simulated data set. 
##' Run 01run_simulations.R to get such a data set. 
##' With the standard parameters in 01run_simulations.R, only the part for 
##' data scenario 1 is reproduced, as the runtime without a computing cluster 
##' would otherwise be too high.



## Load libraries
if (!"pacman" %in% installed.packages()){
  install.packages("pacman")
}
library(pacman)
packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "DescTools", "caret")
p_load(packages, character.only = TRUE)

if("timbR" %in% installed.packages()){
  library(timbR)
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}

#---------------------------------------
## Define directories
## Please define your main directory here. 
## This should be the directory you cloned the git repository into.
main_dir <- this.dir()
setwd(main_dir)

## Define functions directory
proc_dir <- file.path(main_dir, "data")
dir.create(file.path(main_dir, "img"), showWarnings = FALSE)
img_dir <- file.path(main_dir, "img")

#---------------------------------------
## Load and prepare data
#

results <- readRDS(file.path(proc_dir, "results_simulated_results.Rds"))  %>%
  filter(min.bucket == 250) %>%
  filter(metric == "splitting variables" | is.na(metric)) %>%
  filter(probs_quantiles == "" | is.na(probs_quantiles)) %>%
  mutate(scenario2 = case_when(setting == "Setting 1" ~ "large effects",
                               setting == "Setting 2" ~ "small effects",
                               setting == "Setting 3" ~ "correlations",
                               setting == "Setting 4" ~ "interactions",
                               setting == "Setting 5" ~ "continuous variables"),
         scenario2 = factor(scenario2,
                            levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"),
                            labels = c("1: large effects", "2: small effects", "3: correlations", "4: interactions", "5: continuous \nvariables"))) %>%
    mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
                              method == "Regression DT + CPS" ~ "DT + CPS",
                              method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
                              method == "Regression DT + Probability DTs" ~ "Mult. DTs",
                              TRUE ~ method))

# if you want to use the original data from the publication, please download them (see README for details), unoack them, and move them into the data folder: 
# 
# # data from publication
# results <- readRDS(file.path(proc_dir, "results_simulated_results_from_paper.Rds"))  %>% 
#   filter(min.bucket == 150) %>% 
#   filter(metric == "splitting variables" | is.na(metric)) %>% 
#   filter(probs_quantiles == "" | is.na(probs_quantiles)) %>% 
#   mutate(scenario2 = case_when(setting == "Setting 1" ~ "large effects",
#                                setting == "Setting 2" ~ "small effects",
#                                setting == "Setting 3" ~ "correlations",
#                                setting == "Setting 4" ~ "interactions",
#                                setting == "Setting 5" ~ "continuous variables"),
#          scenario2 = factor(scenario2, 
#                             levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
#                             labels = c("1: large effects", "2: small effects", "3: correlations", "4: interactions", "5: continuous \nvariables"))) %>% 
#   mutate(method = case_when(method == "Regression ART + CPS" ~ "ART + CPS",
#                             method == "Regression DT + CPS" ~ "DT + CPS",
#                             method == "Regression ART + Probability ARTs" ~ "Mult. ARTs",
#                             method == "Regression DT + Probability DTs" ~ "Mult. DTs",
#                             TRUE ~ method))

#---------------------------------------
# plot data and save plot

# Transform data to long format
data_accuracy <- bind_rows(results %>% mutate(rmse = sqrt(mse_test_dat_tree)) %>% select(method, scenario2, rmse) %>% dplyr::rename(perf = rmse) %>% mutate(type = "RMSE"),
                           results %>% select(method, scenario2, brier_score0.5) %>% dplyr::rename(perf = brier_score0.5) %>% mutate(type = "Brier score"),
                           results %>% select(method, scenario2, interval_width) %>% dplyr::rename(perf = interval_width) %>% mutate(type = "Interval width") %>% filter(grepl("CPS", method)),
                           results %>% select(method, scenario2, coverage) %>% dplyr::rename(perf = coverage) %>% mutate(type = "Coverage") %>% filter(grepl("CPS", method))
) %>% 
  mutate(type = factor(type, 
                       labels = c("RMSE", "Brier score", "Interval width", "Coverage"), 
                       levels = c("RMSE", "Brier score", "Interval width", "Coverage")))

plot_accuracy <- ggplot(data_accuracy, aes(x=method, y= perf, col = method))+
  geom_boxplot() +
  facet_grid(type~scenario2, scales = "free", switch="y") +
  theme_bw() +
  labs(x = "",
       y = "", 
       col = "")+
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c(rgb(203,81,25, maxColorValue = 255),
                                rgb(0,75,90, maxColorValue = 255),
                                rgb(255, 150, 90, maxColorValue = 255),
                                rgb(80, 180, 210, maxColorValue = 255))) 

plot_accuracy

ggsave(plot_accuracy, filename = file.path(img_dir, "fig2_simulation_prediction_accuracy.png"),
       width = 50, height = 22, units = "cm", dpi = 200)


