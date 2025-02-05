#---------------------------------------
# packages
packages <- c("batchtools", "checkmate", "data.table", "ggplot2", 
              "ranger", "bindata", "rpart", "plyr", "dplyr", 
              "gridExtra", "timbR", "cowplot", "tidyr")
library(reshape2)
invisible(sapply(packages, library, character.only = TRUE))

#---------------------------------------
# paths
data_dir <- "/Users/kronziel/mywork/diss/ART_uncertainty/R/data"
plot_dir <- "/Users/kronziel/mywork/diss/ART_uncertainty/R/img"
tab_dir <- "/Users/kronziel/mywork/diss/ART_uncertainty/R/tab/"

#---------------------------------------
# load data
all_results <- readRDS(file.path(data_dir, "art_uncertainty7.rds"))

# # TODO nur jobs mit min.node.size = 100 verwenden
# data_pred_coverage <- read.csv2(file.path(data_dir, "data_pred_coverage.csv")) %>% 
#   mutate(uncertainty_level2 = (uncertainty_level*100),
#          uncertainty_level3 = paste0(uncertainty_level2, "%")) %>% 
#   filter(min.node.size_art == 100) %>% 
#   mutate(marginal_coverage = 100-uncertainty_level2,
#          width = upper_bound-lower_bound)

data_simstudy <- lapply(all_results, function(x){x[[1]]}) %>% 
  bind_rows() %>% 
  mutate(uncertainty_level2 = (uncertainty_level*100),
         uncertainty_level3 = paste0(uncertainty_level2, "%")) %>% 
  filter(min.node.size_art == 100) %>% 
  mutate(marginal_coverage = 100-uncertainty_level2) %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                               setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                            levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                            labels = c("large effects", "small effects", "correlations", "interactions", "continuous \nvariables")))


data_simstudy_all <- lapply(all_results, function(x){x[[1]]}) %>% 
  bind_rows() %>% 
  mutate(uncertainty_level2 = (uncertainty_level*100),
         uncertainty_level3 = paste0(uncertainty_level2, "%")) %>% 
  mutate(marginal_coverage = 100-uncertainty_level2) %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large effects", "small effects", "correlations", "interactions", "continuous \nvariables")))


abc = data_simstudy %>% filter(is.na(mean_interval_size_mondrian_calibration))
abcd = data_simstudy_all %>% filter(is.na(mean_interval_size_mondrian_calibration))

#---------------------------------------
# Varianz Überdeckungswahrscheinlichkeit
plot_data <-  data_simstudy %>%
  filter(metric=="weighted splitting variables") %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large \neffects", "small \neffects", "correlations", "interactions", "continuous \nvariables"))) %>%
  select(scenario, contains("mean_error"), marginal_coverage, min.node.size_art) %>%
  melt(id=c("scenario", "marginal_coverage", "min.node.size_art")) %>% 
  mutate(variable = ifelse(variable == "mean_error_calibration", "ICP", "Mondrian ICP"))

plt1 <- ggplot(data = plot_data %>% filter(marginal_coverage == 95), aes(x=scenario, y=(100-value), col = variable)) +
  geom_boxplot() +
  theme_bw()+
  labs(y = "Marginal Coverage (%)")+
  geom_hline(yintercept = 95, linetype = "dashed") +
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18)
  )+
  labs(col = "", x = "Data")

plot_data <-  data_simstudy %>%
  filter(metric=="weighted splitting variables") %>% 
  mutate(scenario = case_when(setting == "Setting 1" ~ "large effects",
                              setting == "Setting 2" ~ "small effects",
                              setting == "Setting 3" ~ "correlations",
                              setting == "Setting 4" ~ "interactions",
                              setting == "Setting 5" ~ "continuous variables"),
         scenario = factor(scenario, 
                           levels = c("large effects", "small effects", "correlations", "interactions", "continuous variables"), 
                           labels = c("large \neffects", "small \neffects", "correlations", "interactions", "continuous \nvariables"))) %>%
  select(scenario, contains("interval"), marginal_coverage, min.node.size_art) %>%
  melt(id=c("scenario", "marginal_coverage", "min.node.size_art")) %>% 
  mutate(variable = ifelse(variable == "mean_interval_size_calibration", "ICP", "Mondrian ICP"))

plt2 <- ggplot(data = plot_data %>% filter(marginal_coverage == 95), aes(x=scenario, y=value, col = variable)) +
  geom_boxplot() +
  theme_bw()+
  labs(y = "Interval width")+
  theme(text = element_text(size = 15), legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18)
  )+
  labs(col = "", x = "Data")

plot_grid(plt1, plt2, labels = "AUTO", rel_widths = c(1, 1.3))

ggsave(filename = file.path(plot_dir, "uncertainty_variance.png"), width = 32, height = 9, units = "cm", dpi = 1200)
ggsave(filename = file.path(plot_dir, "uncertainty_variance.pdf"), width = 32, height = 9, units = "cm", dpi = 1200)

# nur Varianz der Intervalbreite
ggsave(plt2, filename = file.path(plot_dir, "uncertainty_width_variance.pdf"), width = 20, height = 9, units = "cm", dpi = 1200)
#---------------------------------------
# Überdeckungswahrscheinlichkeit & confidence / prediction interval width
marginal_coverage_df1 <- aggregate(mean_error_mondrian_calibration ~ marginal_coverage + scenario + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))
marginal_coverage_df2 <- aggregate(mean_error_calibration ~ marginal_coverage + scenario + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))

plot_data2a <- left_join(marginal_coverage_df1, marginal_coverage_df2) %>% 
  melt(id=c("scenario", "marginal_coverage", "min.node.size_art")) %>% 
  mutate(variable = ifelse(variable == "mean_error_calibration", "ICP", "Mondrian ICP"))

plot_coverage <- ggplot(plot_data2a, aes(x=(100-marginal_coverage), y=value, col = scenario, group = scenario))+
  facet_grid(variable~.)+
  theme_bw()+
  geom_line()+
  geom_abline(slope = 1, linetype = "dashed", intercept = 0)+
  labs(x = "Confidence level (%)",
       y = "Marginal Coverage (%)",
       col = "Data") +
  theme(text = element_text(size = 15), legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'),
        strip.text.y = element_text(size = 18))


interval_width_df1 <- aggregate(mean_interval_size_mondrian_calibration ~ uncertainty_level2 + scenario + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))
interval_width_df2 <- aggregate(mean_interval_size_calibration ~ uncertainty_level2 + scenario + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))

plot_data2b <- left_join(interval_width_df1, interval_width_df2) %>% 
  melt(id=c("scenario", "uncertainty_level2", "min.node.size_art")) %>% 
  mutate(variable = ifelse(variable == "mean_interval_size_calibration", "ICP", "Mondrian ICP"))

plot_intervalwidth <- ggplot(plot_data2b, aes(x=(100-uncertainty_level2), y=value, col = scenario, group = scenario))+
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


plot_grid(plot_coverage, plot_intervalwidth, labels = "AUTO", rel_widths = c(1, 1.6))

ggsave(filename = file.path(plot_dir, "marginal_coverage_intervalwidth.png"), width = 25, height = 10, units = "cm", dpi = 1200)
ggsave(filename = file.path(plot_dir, "marginal_coverage_intervalwidth.pdf"), width = 25, height = 10, units = "cm", dpi = 1200)
ggsave(plot_coverage, filename = file.path(plot_dir, "marginal_coverage_intervalwidth_A.pdf"), width = 10, height = 10, units = "cm", dpi = 1200)
ggsave(plot_intervalwidth, filename = file.path(plot_dir, "marginal_coverage_intervalwidth_B.pdf"), width = 15, height = 10, units = "cm", dpi = 1200)



#--------------------------------------- 
# ggplot(data = data_pred_coverage, aes(x=width, y = error_rate, col = uncertainty_level))+
#   geom_point()+
#   facet_grid(setting~.)
#   
# 
# ggplot(data = data_pred_coverage, aes(x=prediction, y = error_rate, col = uncertainty_level))+
#   geom_point()+
#   facet_grid(setting~.)
# 
# ggplot(data = data_pred_coverage, aes(x=prediction, y = width, col = uncertainty_level))+
#   geom_point()+
#   facet_grid(setting~.)

#---------------------------------------
# # Abbildung Vortrag Floris Ernst
# plt1 <- ggplot(marginal_coverage_df2, aes(x=marginal_coverage, y=mean_error_calibration, col = scenario, group = scenario))+
#   theme_bw()+
#   geom_line()+
#   geom_abline(slope = -1, linetype = "dashed", intercept = 100)+
#   labs(x = "Uncertainty level (%)",
#        y = "Marginal Coverage (%)",
#        col = "Simulation scenarios") +
#   theme(text = element_text(size = 15), legend.position = "none",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.6, 'cm'),
#         strip.text.y = element_text(size = 18))
# 
# plt2 <-   ggplot(interval_width_df2, aes(x=uncertainty_level2, y=mean_interval_size_calibration, col = scenario, group = scenario))+
#   theme_bw()+
# geom_line()+
#   labs(x = "Uncertainty level (%)",
#        y = "Interval width",
#        col = "Simulation scenarios") +
#   theme(text = element_text(size = 15), legend.position = "right",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.6, 'cm'),
#         strip.text.y = element_text(size = 18))
# 
# plot_grid(plt1, plt2, labels = "AUTO", rel_widths = c(1, 1.6))
# 
# ggsave(filename = file.path(plot_dir, "interval_width_conformal_prediction_betreuungsgeespraech.png"), width = 25, height = 10, units = "cm", dpi = 1200)



# ---------------------------------------
min.bucket_mondrian <- aggregate(mean_error_mondrian_calibration ~ marginal_coverage + scenario + min.node.size_art, 
                                 data_simstudy_all, function(x) c(mean_error = mean(x)))

table(min.bucket_mondrian$min.node.size_art)

# min.bucket 5 vs 25 vs 100 als plot
# TODO: nur für mondrian
ggplot(min.bucket_mondrian, aes(x=marginal_coverage, y = mean_error_mondrian_calibration, col = as.factor(min.node.size_art), group = min.node.size_art))+
  theme_bw()+
  geom_line()+
  geom_abline(slope = -1, linetype = "dashed", intercept = 100)+
  facet_grid(scenario~.)

# als Tabelle
min.bucket_mondrian_coverage_tab <- aggregate(mean_error_mondrian_calibration ~ uncertainty_level2 + scenario + min.node.size_art, 
                                              data_simstudy_all, 
                                              function(x) c(mean_error = round(mean(x),2))) %>% 
  pivot_wider(., names_from = c(scenario, min.node.size_art), values_from = c(mean_error_mondrian_calibration)) 

min.bucket_mondrian_coverage_tab <- min.bucket_mondrian_coverage_tab %>% 
  select(colnames(min.bucket_mondrian_coverage_tab)[order(colnames(min.bucket_mondrian_coverage_tab))]) %>% 
  select(uncertainty_level2, everything()) %>% 
  arrange(uncertainty_level2)

# abc=data_simstudy_all %>% filter(is.na(mean_error_mondrian_calibration))

write.csv2(min.bucket_mondrian_coverage_tab, file.path(tab_dir, "tab_simulation_mondrian_coverage.csv"), row.names = F)

min.bucket_mondrian_interval_width_tab <- aggregate(mean_interval_size_mondrian_calibration ~ uncertainty_level2 + scenario + min.node.size_art, data_simstudy_all, function(x) c(mean_error = round(mean(x),2))) %>% 
  pivot_wider(., names_from = c(scenario, min.node.size_art), values_from = c(mean_interval_size_mondrian_calibration))

min.bucket_mondrian_interval_width_tab <- min.bucket_mondrian_interval_width_tab %>% 
  select(colnames(min.bucket_mondrian_interval_width_tab)[order(colnames(min.bucket_mondrian_interval_width_tab))]) %>% 
  select(uncertainty_level2, everything()) %>% 
  arrange(uncertainty_level2)

write.csv2(min.bucket_mondrian_interval_width_tab, file.path(tab_dir, "tab_simulation_mondrian_interval_width.csv"), row.names = F)

# 
# #---------------------------------------
# # # plot MSE
# # ggplot(data_simstudy, aes(x=setting, y = mse_val_dat_rf))+
# #   geom_boxplot()
# 
# #--------------------------------------- 
# 
# 
# 
# 
# ggsave(filename = file.path(plot_dir, "uncertainty_errorrate.png"), width = 20, height = 20, units = "cm", dpi = 1200)
# 
# # three depth
# ggplot(data_simstudy %>% filter(uncertainty_level==0.05), aes(x=setting, y=number_splits, col = factor(min.node.size_art)))+
#   theme_bw()+
#   geom_boxplot()+
#   labs(x = "Data setting",
#        y = "Number of splits in ART",
#        col = "min.bucket in ART") +
#   theme(text = element_text(size = 15), legend.position = "right",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.6, 'cm'),
#         strip.text.y = element_text(size = 18)
#   ) 
# 
# 
# 
# # intereval width
# error_df1 <- aggregate(mean_interval_size_mondrian_calibration ~ uncertainty_level2 + setting + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))
# error_df2 <- aggregate(mean_interval_size_calibration ~ uncertainty_level2 + setting + min.node.size_art, data_simstudy, function(x) c(mean_error = mean(x)))
# 
# plot_data3 <- left_join(error_df1, error_df2) %>% 
#   melt(id=c("setting", "uncertainty_level2", "min.node.size_art")) %>% 
#   mutate(variable = ifelse(variable == "mean_interval_size_calibration", "Conformal regression", "Mondrian conformal regression"))
# 
# # x-achse uncertainty
# # linien mittelwert der settings, im Optimalfall winkelhalbierende
# 
# ggplot(plot_data3, aes(x=uncertainty_level2, y=value, col = setting, group = setting))+
#   facet_grid(variable~.)+
#   theme_bw()+
#   geom_line()+
#   labs(x = "Uncertainty level (%)",
#        y = "Interval width",
#        col = "Data") +
#   theme(text = element_text(size = 15), legend.position = "right",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.6, 'cm'),
#         strip.text.y = element_text(size = 18)
#   ) 
# 
# ggsave(filename = file.path(plot_dir, "uncertainty_intervalwidth.png"), width = 20, height = 20, units = "cm", dpi = 1200)
# 
# 
# plot_data4 <-  data_simstudy %>%
#   filter(metric=="weighted splitting variables") %>%
#   select(setting, contains("mean_interval"), uncertainty_level3, min.node.size_art) %>%
#   melt(id=c("setting", "uncertainty_level3", "min.node.size_art")) %>% 
#   mutate(variable = ifelse(variable == "mean_interval_size_calibration", "Conformal regression", "Mondrian conformal regression"))
# 
# 
# ggplot(data = plot_data4 %>% filter(uncertainty_level3 == "10%"), aes(x=setting, y=value, col = variable)) +
#   geom_boxplot() +
#   theme_bw()+
#   labs(y = "Interval width")+
#   theme(text = element_text(size = 15), legend.position = "right",
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         legend.key.size = unit(0.6, 'cm'),
#         strip.text.y = element_text(size = 18)
#   )+
#   labs(col = "Calibration type", x = "Data")
# 
# ggsave(filename = file.path(plot_dir, "uncertainty_intervalwidth_variance.png"), width = 20, height = 20, units = "cm", dpi = 1200)
# 
# 


#--------------
# tabellen mit einer Nachkommastelle
tab_coverage <- read.csv2(file.path("/Users/kronziel/mywork/diss/ART_uncertainty/R/tab/tab_simulation_mondrian_coverage_paper.csv"),
                          skip = 1)

tab_coverage_rounded <- tab_coverage %>% 
  mutate(across(everything(), ~round(.x, 1))) 

write.csv(tab_coverage_rounded, file.path("/Users/kronziel/mywork/diss/ART_uncertainty/R/tab/tab_simulation_mondrian_coverage_paper_rounded.csv"), row.names = F)

