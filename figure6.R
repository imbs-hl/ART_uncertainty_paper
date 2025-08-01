#' With this script the figure 6 from "Uncertainty 
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
# Functions for beautiful plots
# Source: https://github.com/wilkelab/cowplot/issues/202
#---------------------------------------
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}
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

#---------------------------------------
# Load and prepare data
# # Data from publication
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
                                     TRUE ~ task_name),
         method = ifelse(method == "artificial tree with uncertainty", "ART", "MRT"))


#---------------------------------------
# Plot data from figure 6 and save plot
# Plot interval width for ICP, Mondrian ICP, CPS, Mondrian CPS
#---------------------------------------

#---------------------------------------
# Figure 3A
# Interval width of ART and MRT for ICP, CPS and Mondrian variants
#---------------------------------------

# Filter and aggregate data for each uncertainty method
data_mrt <- results %>% filter((min.bucket == 100 & method == "ART" & metric == "weighted splitting variables") | 
                                 method == "MRT" & metric == "weighted splitting variables")
mrt_df1 <- aggregate(mean_interval_size_mondrian_icp ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df2 <- aggregate(mean_interval_size_icp ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))

# Combine all aggregated uncertainty quantification to long format
plot_data_mrt <- full_join(mrt_df1, mrt_df2) %>% 
  full_join(mrt_df3) %>% 
  full_join(mrt_df4) %>% 
  reshape2::melt(id=c("task_name_short", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS")) %>% filter(metric == "weighted splitting variables")


# Plot interval width for first data set abalone
plot_mrt1 = ggplot(plot_data_mrt %>% filter(task_name_short == "abalone"), 
                   aes(x=significance_level, y=value, col = method, group = method))+
  facet_wrap(.~variable, ncol = 2)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "")+
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))

plot_mrt1

#---------------------------------------
# Figure 3B
# Jaccard index
#---------------------------------------
# Aggregate data
jaccard_df1 <- aggregate(mean_jaccard_cps ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))
jaccard_df2 <- aggregate(mean_jaccard_icp ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))

# Merge data into long format
plot_data_dist_jaccard <- full_join(jaccard_df1, jaccard_df2) %>% 
  reshape2::melt(id=c("task_name_short", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_jaccard_icp" ~ "ICP",
                              variable == "mean_jaccard_cps" ~ "CPS",
                              TRUE ~ variable)) 

# Build plot
plot_dist_jaccard1 = ggplot(plot_data_dist_jaccard %>% filter(task_name_short == "abalone"), 
                            aes(x=significance_level, y=value, col = method, group = method))+
  facet_wrap(.~variable, ncol = 1)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Jaccard index",
       col = "")+
  theme(text = element_text(size = 15), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))

plot_dist_jaccard1
#---------------------------------------


#---------------------------------------
# Figure 3C
# Earth mover's distance
#---------------------------------------

# Aggregate data
earth_movers_df1 <- aggregate(mean_earth_movers_distance_cps ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))
earth_movers_df2 <- aggregate(mean_earth_movers_distance_icp ~ method + metric + significance_level + task_name_short, data_mrt, function(x) c(mean_error = mean(x)))

# Merge data into long format
plot_data_dist_earth_movers <- full_join(earth_movers_df1, earth_movers_df2) %>% 
  reshape2::melt(id=c("task_name_short", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_earth_movers_distance_icp" ~ "ICP",
                              variable == "mean_earth_movers_distance_cps" ~ "CPS",
                              TRUE ~ variable)) %>% 
  filter(metric == "weighted splitting variables")

# Build plot
plot_dist_earth_movers1 = ggplot(plot_data_dist_earth_movers %>% filter(task_name_short == "abalone"), 
                                 aes(x=significance_level, y=value, col = method, group = method))+
  facet_wrap(.~variable, ncol = 1)+ #, ncol=1
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Earth mover's distance",
       col = "")+
  theme(text = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))+ 
  guides(colour = guide_legend(nrow = 1))

plot_dist_earth_movers1
#---------------------------------------


#---------------------------------------
# Combine plots to one figure
#---------------------------------------
plots_mrts_arts <- plot_grid(plot_mrt1 + theme(legend.position = "none"),
                             plot_dist_jaccard1+ theme(legend.position = "none"),
                             plot_dist_earth_movers1 + theme(legend.position = "none"),
                             ncol = 3, align = "hv",
                             labels = c("A", "B", "C"), rel_widths = c(1.5,1,1))
legend_mrts_arts <- get_legend(plot_dist_earth_movers1)

plots_mrts_arts <- plot_grid(plots_mrts_arts, legend_mrts_arts, nrow = 2, rel_heights = c(1, 0.1))

plots_mrts_arts

# save plot
ggsave(plots_mrts_arts, filename = file.path(out_dir, "fig6_benchmark_abalone.png"),
       width = 25, height = 15, units = "cm", dpi = 200)

