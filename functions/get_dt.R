get_dt <- function(data, instance, min.bucket = 0, significance_level = 0.1, ...){
  
  # Exctract data from instance 
  #----------------
  train_data     <- instance[[1]]
  test_data      <- instance[[2]]
  cal_data       <- instance[[3]]
  rf             <- instance[[4]]
  information_df <- instance[[5]]
  
  # MSE von RF
  pred_regr_rf <- predict(rf, test_data)$predictions
  mse_regr_rf <- mean((test_data$glycohemoglobin - pred_regr_rf)^2)
  
  # build DT
  #----------------
  start <- proc.time()
  # Regression DT trainieren
  art_regr <- ranger(glycohemoglobin~., data = train_data, 
                     min.bucket = min.bucket, mtry = (ncol(train_data)-1),
                     sample.fraction = 1, replace = FALSE, num.trees = 1)
  

  # Vorhersagen von Regression ART
  pred_regr_art <- predict(art_regr, test_data)$predictions
  
  
  # build cps
  #----------------
  # Probabilities von CPS berechnen
  y_cal <- cal_data$glycohemoglobin
  y_test <- test_data$glycohemoglobin
  y_cal_pred <- predict(art_regr, cal_data)$predictions
  y_test_pred <- predict(art_regr, test_data)$predictions
  
  
  cpd_art <- get_calibrated_predictive_system_mondrian(y_cal_pred = y_cal_pred, y_cal = y_cal, y_test_pred = y_test_pred, 
                                                       significance_level = significance_level,interval_type = "two-tailed", 
                                                       tree = art_regr, 
                                                       cal_data = cal_data, test_data = test_data, show_node_id = TRUE, 
                                                       dependent_varname = "glycohemoglobin") %>% 
    mutate(nodeID = leaf-1)
  
  tree_info_df <- treeInfo(art_regr) %>% 
    left_join(unique(cpd_art)) %>% 
    mutate(lower_bound = round(lower_bound, 2),
           upper_bound = round(upper_bound, 2),
           prediction = round(prediction, 2))
  
  # Vorhersage vs. wahre Probability pro Knoten
  # pro terminalen Knoten die cpd erstellen
  # Vektor mit allen Blättern
  leafs <- tree_info_df %>% filter(terminal) %>% select(nodeID) %>% unlist() %>% unname()
  # Vektor, welche Beobachtung in welchem Blatt landet
  test_data_leafs <- predict(art_regr, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(art_regr, train_data, type = "terminalNodes")$predictions
  cal_data_leafs <- predict(art_regr, cal_data, type = "terminalNodes")$predictions
  
  probs_df <- list()
  
  for(i in 1:length(leafs)){
    # alle Kalibrierungsdaten aus dem Blatt
    cal_dat_leaf <- cal_data[cal_data_leafs==leafs[i],]
    # eine Testbeobachtung aus dem Blatt
    test_dat_leaf <- test_data[test_data_leafs==leafs[i],][1,]
    # cpd
    cpd_dist <- get_predictive_distribution(y_cal = cal_dat_leaf$glycohemoglobin, 
                                            y_cal_pred = predict(art_regr, cal_dat_leaf)$predictions, 
                                            y_test_pred = predict(art_regr, test_dat_leaf)$predictions)[[1]]
    # "p-Wert" berechnen
    p <- seq(0, 1 - 1/length(cpd_dist), length.out = length(cpd_dist))  
    data_cps <- data.frame(cpd = c(cpd_dist,max(cpd_dist+0.000001), max(cal_data$glycohemoglobin)), p = c(p,1,1))
    
    # W' jeweiligen Schwellenwert zu überschreiten
    p_5.7 <- data_cps$p[which(data_cps$cpd > 5.7)[1]]
    p_6.5 <- data_cps$p[which(data_cps$cpd > 6.5)[1]]
    
    probs_df[[i]] <- c(nodeID=leafs[i], prob5.7=1-p_5.7, prob6.5=1-p_6.5)
    
    # Obere Grenze des Prediction Intervals
    # for two tailed
    lower_index <- tail(which(data_cps$p <= 0.025), 1)
    upper_index <- which(data_cps$p >= 0.975)[1]
    
    low_percentile <-  data_cps$cpd[lower_index]
    high_percentile <-  data_cps$cpd[upper_index]
  }
  
  probs_df <- left_join(tree_info_df, bind_rows(probs_df))
  
  # rel. und abs. Häufigkeit der Testdaten pro Blatt über Schwellenwerten zu liegen
  prob_df_i <- data.frame(nodeID = predict(art_regr, test_data, type = "terminalNodes")$prediction,
                          glycohemoglobin = test_data$glycohemoglobin) %>% 
    group_by(nodeID)  %>%
    dplyr::summarise(
      n_total = n(),  # Gesamtanzahl der Beobachtungen pro leaf
      n_over_5_7 = sum(glycohemoglobin > 5.7, na.rm = TRUE),  # Anzahl über 5.7
      proportion_over_5_7_test_data = mean(glycohemoglobin > 5.7, na.rm = TRUE),  # Anteil über 5.7
      n_over_6_5 = sum(glycohemoglobin > 6.5, na.rm = TRUE),  # Anzahl über 6.5
      proportion_over_6_5_test_data = mean(glycohemoglobin > 6.5, na.rm = TRUE)  # Anteil über 6.5
    ) %>% left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound, prob5.7, prob6.5))
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # calculate performance measures
  #----------------
  # MSE tree
  mse_regr_art <- mean((test_data$glycohemoglobin - pred_regr_art)^2)
  
  # Baumtiefe 
  num_splits_regr_art <- nrow(treeInfo(art_regr) %>% filter(!terminal))
  
  # Brier Score pro Testbeobachtung 
  brier_score_df <- data.frame(y = y_test, 
                               nodeID = test_data_leafs) %>% 
    left_join(prob_df_i %>% select(nodeID, prob5.7, prob6.5)) %>% 
    mutate(over5.7 = ifelse(y>= 5.7,1,0),
           over6.5 = ifelse(y>= 6.5,1,0))
  
  brier_score5.7 = BrierScore(brier_score_df$over5.7, brier_score_df$prob5.7)
  brier_score6.5 = BrierScore(brier_score_df$over6.5, brier_score_df$prob6.5)
  
  # Coverage von Testdaten
  coverage_df <- data.frame(y = y_test, 
                            nodeID = test_data_leafs) %>% 
    left_join(probs_df %>% select(nodeID, prediction, lower_bound, upper_bound)) %>% 
    mutate(inside_interval = y >= lower_bound & y <= upper_bound)
  
  coverage <- sum(coverage_df$inside_interval)/nrow(coverage_df)
  
  # interval width
  interval_width_df <- cpd_art %>% 
    select(nodeID, lower_bound, upper_bound) %>% 
    unique()
  
  interval_width <- mean(interval_width_df$upper_bound - interval_width_df$lower_bound)
  
  # längster Pfad
  ti <- treeInfo(art_regr)
  
  # ranger gibt Bäume als Knotenliste mit Child-Knoten an
  # Wir bauen die Tiefe per BFS/DFS nach
  
  # Startknoten (root = NodeID 0)
  depths <- setNames(integer(nrow(ti)), ti$nodeID)
  depths["0"] <- 0
  
  # Für jeden Knoten die Tiefe bestimmen
  for (i in seq_len(nrow(ti))) {
    node <- ti$nodeID[i]
    
    # Wenn Terminal Node -> überspringen
    if (ti$terminal[i]) next
    
    # linken und rechten Kindknoten auslesen
    left  <- ti$leftChild[i]
    right <- ti$rightChild[i]
    
    depths[as.character(left)]  <- depths[as.character(node)] + 1
    depths[as.character(right)] <- depths[as.character(node)] + 1
  }
  
  max_depth <- max(depths)
  
  # Anzahl terminal nodes
  num_leaves <- nrow(treeInfo(art_regr) %>% filter(terminal))
  
  # Sparsity
  num_vars_split <- length(na.omit(unique(treeInfo(art_regr)$splitvarID)))
  
  
  # return model, test_& cal data, performance df
  #----------------
  return(list(info_df = data.frame(method     = "Regression DT + CPS",
                         metric               = NA, 
                         mse_test_dat_tree    = mse_regr_art,
                         mse_test_dat_rf      = mse_regr_rf,
                         probs_quantiles      = NA,
                         epsilon              = NA,
                         min.bucket           = min.bucket,
                         significance_level   = significance_level,
                         r.squared_ranger     = rf$r.squared,
                         num.splits           = paste0(num_splits_regr_art, collapse = ","),
                         runtime              = time,
                         n_test = nrow(test_data),
                         n_train = nrow(train_data),
                         n_cal = nrow(cal_data),
                         brier_score5.7 = brier_score5.7,
                         brier_score6.5 = brier_score6.5, # ab hier
                         interval_width = interval_width,
                         coverage = coverage,
                         max_depth = max_depth,
                         num_leaves = num_leaves,
                         num_vars_split = num_vars_split
  ) %>% cbind(information_df),
  tree = art_regr,
  test_data = test_data,
  cal_data = cal_data,
  prob_df = prob_df_i %>% select(nodeID, prob5.7, prob6.5)
  ))
}