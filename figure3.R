#' With this script the figure 3 from "Uncertainty 
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
pacman::p_load(gtable)
pacman::p_load(lemon)
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
# # # Data from publication
# # # min.bucket = 100 and metric = "weighted splitting variables" was used for this figure
# results <- read.csv2(file.path(proc_dir, "results_simulations.csv")) %>%
#   # Change names and calculate confidence level in percent as (1 - significance level)*100 and coverage as 1 - error
#     mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
#                                 setting == "Setting 2" ~ "small effects",
#                                 setting == "Setting 3" ~ "correlations",
#                                 setting == "Setting 4" ~ "interactions",
#                                 setting == "Setting 5" ~ "continuous variables"),
#            scenario = factor(scenario,
#                              levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"),
#                              labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))

# Data produced by simulations.R
 results <- readRDS(file.path(proc_dir, "results.Rds")) %>%
   bind_rows()  %>%
# Change names and calculate confidence level in percent as (1 - significance level)*100 and coverage as 1 - error
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario,
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"),
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")),
         method = ifelse(method == "artificial tree with uncertainty", "ART", "MRT"))

#---------------------------------------
# Plot data from figure 3 and save plot
# Plot interval width ICP, Mondrian ICP, CPS, Mondrian CPS
#---------------------------------------

#---------------------------------------
# Figure 3A
# Interval width of ART and MRT for ICP, CPS and Mondrian variants
#---------------------------------------

# Filter and aggregate data for each uncertainty method
data_mrt <- results %>% 
  filter((min.bucket == 100 & method == "ART" & metric == "weighted splitting variables") | 
           method == "MRT" & metric == "weighted splitting variables")
mrt_df1 <- aggregate(mean_interval_size_mondrian_icp ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df2 <- aggregate(mean_interval_size_icp ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df3 <- aggregate(mean_interval_width_cps_two_tailed ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))
mrt_df4 <- aggregate(mean_interval_width_cps_two_tailed_mondrian ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))

# Combine all aggregated uncertainty quantification to long format
plot_data_mrt <- full_join(mrt_df1, mrt_df2) %>% 
  full_join(mrt_df3) %>% 
  full_join(mrt_df4) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_interval_size_icp" ~ "ICP",
                              variable == "mean_interval_size_mondrian_icp" ~ "Mondrian ICP",
                              variable == "mean_interval_width_cps_two_tailed" ~ "CPS",
                              TRUE ~ "Mondrian CPS")) %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))%>% 
  filter(scenario == "large effects") %>% filter(metric == "weighted splitting variables")

# Build plot
plot_mrt = ggplot(plot_data_mrt , 
                  aes(x=significance_level, y=value, col = method, group = method))+
  facet_wrap(.~variable, ncol=2)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Interval width",
       col = "")+
  theme(text = element_text(size = 15), #legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))
plot_mrt


#---------------------------------------
# Figure 3B
# Jaccard index
#---------------------------------------
# Aggregate data
jaccard_df1 <- aggregate(mean_jaccard_cps ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))
jaccard_df2 <- aggregate(mean_jaccard_icp ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))

# Merge data into long format
plot_data_dist_jaccard <- full_join(jaccard_df1, jaccard_df2) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_jaccard_icp" ~ "ICP",
                              variable == "mean_jaccard_cps" ~ "CPS",
                              TRUE ~ variable)) %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))

# Build plot
plot_dist_jaccard = ggplot(plot_data_dist_jaccard %>% filter(scenario == "large effects") , 
                           aes(x=significance_level, y=value, col = method, group = method))+
  facet_wrap(.~variable, ncol=1)+
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Jaccard index",
       col = "")+
  theme(text = element_text(size = 15), #legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))
plot_dist_jaccard

#---------------------------------------


#---------------------------------------
# Figure 3C
# Earth mover's distance
#---------------------------------------

# Aggregate data
earth_movers_df1 <- aggregate(mean_earth_movers_distance_cps ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))
earth_movers_df2 <- aggregate(mean_earth_movers_distance_icp ~ method + metric + significance_level + scenario, data_mrt, function(x) c(mean_error = mean(x)))

# Merge data into long format
plot_data_dist_earth_movers <- full_join(earth_movers_df1, earth_movers_df2) %>% 
  reshape2::melt(id=c("scenario", "significance_level", "method", "metric")) %>% 
  mutate(variable = case_when(variable == "mean_earth_movers_distance_icp" ~ "ICP",
                              variable == "mean_earth_movers_distance_cps" ~ "CPS",
                              TRUE ~ variable)) %>% 
  mutate(scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous variables")))

# Build plot
plot_dist_earth_movers = ggplot(plot_data_dist_earth_movers %>% filter(scenario == "large effects"), 
                                aes(x=significance_level, y=value, col = method, group = method))+
  #facet_wrap(.~scenario, ncol=1)+
  facet_wrap(.~variable, ncol=1)+ #, ncol=1
  theme_bw()+
  geom_line()+
  labs(x = "Significance level (%)",
       y = "Earth mover's distance",
       col = "")+
  theme(text = element_text(size = 15), #legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18),
        strip.text.x = element_text(size = 18),
        axis.text = element_text(size = 10))+ 
  guides(colour = guide_legend(nrow = 1))

plot_dist_earth_movers
#---------------------------------------


#---------------------------------------
# Combine plots to one figure
#---------------------------------------
plots_mrts_arts <- plot_grid(plot_mrt + theme(legend.position = "none"),
                             plot_dist_jaccard+ theme(legend.position = "none"),
                             plot_dist_earth_movers + theme(legend.position = "none"),
                             ncol = 3, align = "hv",
                             labels = c("A", "B", "C"), rel_widths = c(1.75,1,1))
legend_mrts_arts <- get_legend(plot_dist_earth_movers)

plots_mrts_arts <- plot_grid(plots_mrts_arts, legend_mrts_arts, nrow = 2, rel_heights = c(1, 0.1))

plots_mrts_arts

ggsave(plots_mrts_arts, filename = file.path(out_dir, "fig3_simulation_mrts_arts_largeeffects.png"),
       width = 25, height = 15, units = "cm", dpi = 200)






