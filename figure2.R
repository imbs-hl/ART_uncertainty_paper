#' With this script the figure 2 from "Uncertainty 
#' quantification enhances the explainability of 
#' artificial representative trees" can be reproduced. 
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


if("timbR" %in% installed.packages()){
  library(timbR)
  warning("Please check, if timbR version 3.1 is installed.")
} else {
  devtools::install_github("imbs-hl/timbR", "master")
  library(timbR)
}

#---------------------------------------
# Functions to simulate the data sets
source("functions/simulate_rf_setting_1.R")

#---------------------------------------
# Figure 2 
#---------------------------------------

set.seed(123)
# Simmulate data
instance = simulate_rf_setting_1(data = 1000, n_test = 1000, n_cal = 1000, p = 100, p_eff = 5, beta_eff = 2, eps = 1, num.trees = 500, mtry = 10, min_node_size = 100)

# Exctract data from instance 
test_dat       <- instance[[1]]
rf             <- instance[[2]]
params         <- instance[[3]]
cal_dat        <- instance[[4]]
effect_var_ids <- instance[[5]]
noise_var_ids  <- instance[[6]]
train_dat      <- instance[[7]]
dependent_varname <- instance[[8]]

# Create ART 
rf_rep <- generate_tree(rf = rf, metric = "weighted splitting variables", test_data = test_dat, train_data = train_dat, 
                        importance.mode = TRUE, imp.num.var = 5, dependent_varname = dependent_varname,
                        probs_quantiles = NULL, epsilon = 0.05,
                        min.bucket = 100, num.splits = NULL)

# Calculate predictions for calidation and test data
y_cal <- cal_dat[,dependent_varname]
y_test <- test_dat[,dependent_varname]
y_cal_pred <- predict(rf_rep, cal_dat)$predictions
y_test_pred <- predict(rf_rep, test_dat)$predictions

# Calculate conformal predictive system (CPS)
# Select random index of test observation
idx_test <- sample(c(1:length(y_test)), 1)

cpd <- y_test_pred[idx_test] + sort(y_cal_pred - y_cal)

p <- seq(0, 1 - 1/length(cpd), length.out = length(cpd))  

data_cps <- data.frame(cpd = cpd, p = p)

# For two tailed
lower_index <- tail(which(p <= 0.025), 1)
upper_index <- which(p >= 0.975)[1]

low_percentile <- cpd[lower_index]
high_percentile <- cpd[upper_index]

# Median
mid_index <- which(p >= 0.5)[1]
median <- cpd[mid_index]

# For one tailed
lower_index2 <- tail(which(p <= 0.05), 1)
upper_index2 <- which(p >= 0.95)[1]

low_percentile2 <- cpd[lower_index2]
high_percentile2 <- cpd[upper_index2]

#---------------------------------------
# Plot distribution
#---------------------------------------
plt_two_tailed <- ggplot() +
  theme_bw() +
  geom_line(data = data_cps, aes(x = cpd, y = p))  +
  geom_segment(aes(x = y_test_pred[idx_test], y = 0, yend = 1, color = "prediction", linetype = "prediction")) +
  geom_segment(aes(x = y_test[idx_test], y = 0, yend = 1, color = "y", linetype = "y")) +
  geom_area(data = data_cps %>% filter(cpd >= low_percentile & cpd <= high_percentile),
            aes(x = cpd, y = p), fill = "red", alpha = 0.09) +
  xlab("y") +
  ylab("Q(y)") +
  ylim(0, 1) +
  geom_segment(aes(x = median, y = 0, yend = 1, color = "median", linetype = "median"), show.legend = TRUE) +
  geom_segment(aes(x = min(cpd), xend = 0.4, y = 0.79, color="P(y)≤0.4", linetype = "P(y)≤0.4"), show.legend = TRUE) +
  geom_segment(aes(x = 0.4, y = 0, yend = 0.79, color="P(y)≤0.4", linetype = "P(y)≤0.4"), show.legend = TRUE) +
  geom_segment(aes(x = low_percentile, y = 0, yend = p[lower_index], color="95% prediction interval", linetype="95% prediction interval"), show.legend = TRUE) +
  geom_segment(aes(x = high_percentile, y = 0, yend = p[upper_index], color="95% prediction interval", linetype="95% prediction interval"), show.legend = TRUE) +
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.6, 'cm'))+
  labs(color = "")+
  scale_color_manual(name = "",
                     values=c(rgb(228,32,50, maxColorValue = 255),
                              rgb(236,116,4, maxColorValue = 255),
                              "black", 
                              rgb(149,188,14, maxColorValue = 255),
                              #rgb(0,174,199, maxColorValue = 255),
                              rgb(0,106,163, maxColorValue = 255)))+
  scale_linetype_manual(name = "",values=c(1,1,2,1,1))+
  labs(color = "", linetype = "")

plt_two_tailed

plt_left_tailed <- ggplot() +
  theme_bw() +
  geom_line(data = data_cps, aes(x = cpd, y = p))  +
  geom_segment(aes(x = y_test_pred[idx_test], y = 0, yend = 1, color = "prediction")) +
  geom_segment(aes(x = y_test[idx_test], y = 0, yend = 1, color = "y")) +
  geom_area(data = data_cps %>% filter(cpd <= high_percentile2),
            aes(x = cpd, y = p), fill = "red", alpha = 0.09) +
  xlab("y") +
  ylab("Q(y)") +
  ylim(0, 1) +
  xlim(min(cpd), max(cpd)) +
  geom_segment(aes(x = median, y = 0, yend = 1, color = "median"), show.legend = TRUE) +
  geom_segment(aes(x = high_percentile2, y = 0, yend = p[upper_index2], color="95% prediction interval"), show.legend = TRUE) +
  geom_segment(aes(x = min(cpd), xend = 0.4, y = 0.79, color="P(y)≤0.4"), linetype = "dashed", show.legend = TRUE) +
  geom_segment(aes(x = 0.4, y = 0, yend = 0.79, color="P(y)≤0.4"), linetype = "dashed", show.legend = TRUE) +
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.6, 'cm'))+
  labs(color = "")+
  scale_color_manual(name = "",
                     values=c(rgb(228,32,50, maxColorValue = 255),
                              rgb(236,116,4, maxColorValue = 255),
                              "black", 
                              rgb(149,188,14, maxColorValue = 255),
                              #rgb(0,174,199, maxColorValue = 255),
                              rgb(0,106,163, maxColorValue = 255)))+
  scale_linetype_manual(name = "",values=c(1,1,2,1,1))+
  labs(color = "", linetype = "")

plt_left_tailed

plt_right_tailed <- ggplot() +
  theme_bw() +
  geom_line(data = data_cps, aes(x = cpd, y = p))  +
  geom_segment(aes(x = y_test_pred[idx_test], y = 0, yend = 1,color = "prediction", linetype = "prediction")) +
  geom_segment(aes(x = y_test[idx_test], y = 0, yend = 1, color = "y", linetype = "y")) +
  geom_area(data = data_cps %>% filter(cpd >= low_percentile2),
            aes(x = cpd, y = p), fill = "red", alpha = 0.09) +
  xlab("y") +
  ylab("Q(y)") +
  ylim(0, 1) +
  xlim(min(cpd), max(cpd)) +
  geom_segment(aes(x = median, y = 0, yend = 1, color = "median", linetype = "median"), show.legend = TRUE) +
  geom_segment(aes(x = low_percentile2, y = 0, yend = p[lower_index2], color="95% prediction interval", linetype="95% prediction interval"), show.legend = TRUE) +
  geom_segment(aes(x = min(cpd), xend = 0.4, y = 0.79, color="P(y)≤0.4", linetype="P(y)≤0.4"), show.legend = TRUE) +
  geom_segment(aes(x = 0.4, y = 0, yend = 0.79, color="P(y)≤0.4", linetype="P(y)≤0.4"), show.legend = TRUE) +
  theme(text = element_text(size = 15), legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.6, 'cm'))+
  scale_color_manual(name = "",
                     values=c(rgb(228,32,50, maxColorValue = 255),
                              rgb(236,116,4, maxColorValue = 255),
                              "black", 
                              rgb(149,188,14, maxColorValue = 255),
                              #rgb(0,174,199, maxColorValue = 255),
                              rgb(0,106,163, maxColorValue = 255)))+
  scale_linetype_manual(name = "",values=c(1,1,2,1,1))+
  labs(color = "", linetype = "")

plt_right_tailed


# Combine plots
plots_row <- plot_grid(plt_left_tailed, 
                       plt_two_tailed, 
                       plt_right_tailed + theme(legend.position = "none"), 
                       nrow = 1, align = "hv",
                       labels = c("A", "B", "C"))

# Get legend (copied from https://github.com/wilkelab/cowplot/issues/202)
get_legend <- function(plot) {
  # return all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  # return first non-zero legend if exists, and otherwise first element (which will be a zeroGrob) 
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}
legend_tailed <- get_legend(plt_right_tailed)

# Add legend underneath
final_plot <- plot_grid(plots_row, legend_tailed, ncol = 1, rel_heights = c(1, 0.1))

final_plot

ggsave(final_plot, filename = file.path(out_dir, "fig2_simulation_cps_all_tailed.png"), width = 20, height = 6, units = "cm", dpi = 200)




