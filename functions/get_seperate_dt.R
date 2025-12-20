get_seperate_dt <- function(data, instance, min.bucket = 0, ...){
  
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
  
  # build DT for regression
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
  
  # Probability DT for prädiabetes trainieren
  train_data_5.7 <- train_data %>% 
    mutate(glycohemoglobin = ifelse(glycohemoglobin>=5.7,1,0))
  
  art_prob5.7 <- ranger(glycohemoglobin~., data = train_data_5.7, 
                     min.bucket = min.bucket, mtry = (ncol(train_data)-1),
                     sample.fraction = 1, replace = FALSE, num.trees = 1,
                     probability = TRUE)
  # Probability DT for Diabetes trainieren
  train_data_6.5 <- train_data %>% 
    mutate(glycohemoglobin = ifelse(glycohemoglobin>=6.5,1,0))
  
  art_prob6.5 <- ranger(glycohemoglobin~., data = train_data_6.5, 
                        min.bucket = min.bucket, mtry = (ncol(train_data)-1),
                        sample.fraction = 1, replace = FALSE, num.trees = 1,
                        probability = TRUE)
  
  end <- proc.time()
  time <- as.numeric((end - start)[1])
  
  # calculate performance measures
  #----------------
  # MSE tree
  mse_regr_art <- mean((test_data$glycohemoglobin - pred_regr_art)^2)
  
  # Baumtiefe 
  num_splits_regr_art <- mean(c(nrow(treeInfo(art_regr) %>% filter(!terminal)),
                                nrow(treeInfo(art_prob5.7) %>% filter(!terminal)),
                                nrow(treeInfo(art_prob6.5) %>% filter(!terminal))
  ))
  
  # Brier Score pro Testbeobachtung 
  # Vektor, welche Beobachtung in welchem Blatt landet
  test_data_leafs <- predict(art_regr, test_data, type = "terminalNodes")$predictions
  train_data_leafs <- predict(art_regr, train_data, type = "terminalNodes")$predictions
  cal_data_leafs <- predict(art_regr, cal_data, type = "terminalNodes")$predictions
  
  brier_score_df <- data.frame(y = y_test, 
                               nodeID = test_data_leafs,
                               prob5.7 = predict(art_prob5.7, test_data)$predictions[,2],
                               prob6.5 = predict(art_prob6.5, test_data)$predictions[,2]
                               )  %>% 
    mutate(over5.7 = ifelse(y>= 5.7,1,0),
           over6.5 = ifelse(y>= 6.5,1,0))
  
  brier_score5.7 = BrierScore(brier_score_df$over5.7, brier_score_df$prob5.7)
  brier_score6.5 = BrierScore(brier_score_df$over6.5, brier_score_df$prob6.5)
  
  
  # längster Pfad
  get_max_depth <- function(rf_object, tree = 1) {
    # Baum-Struktur extrahieren
    ti <- ranger::treeInfo(rf_object, tree = tree)
    
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
    
    # Maximale Tiefe zurückgeben
    max(depths)
  }
  
  max_depth <- mean(max(get_max_depth(art_regr)), max(get_max_depth(art_prob5.7)), max(get_max_depth(art_prob6.5)))
  
  # Anzahl terminal nodes
  num_leaves_regr <- nrow(treeInfo(art_regr) %>% filter(terminal))
  num_leaves_prediabetes <- nrow(treeInfo(art_prob5.7) %>% filter(terminal))
  num_leaves_diabetes <- nrow(treeInfo(art_prob6.5) %>% filter(terminal))
  num_leaves <- mean(num_leaves_regr, num_leaves_prediabetes, num_leaves_diabetes)
  
  # Sparsity
  num_vars_split_regr <- length(na.omit(unique(treeInfo(art_regr)$splitvarID)))
  num_vars_split_prediabetes <- length(na.omit(unique(treeInfo(art_prob5.7)$splitvarID)))
  num_vars_split_diabetes <- length(na.omit(unique(treeInfo(art_prob6.5)$splitvarID)))
  num_vars_split <- mean(num_vars_split_regr, num_vars_split_prediabetes, num_vars_split_diabetes)
  
  
  
  # return model, test_& cal data, performance df
  #----------------
  return(list(info_df = data.frame(method     = "Regression DT + Probability DTs",
                         metric               = NA, 
                         mse_test_dat_tree    = mse_regr_art,
                         mse_test_dat_rf      = mse_regr_rf,
                         probs_quantiles      = NA,
                         epsilon              = NA,
                         min.bucket           = min.bucket,
                         significance_level   = NA,
                         r.squared_ranger     = rf$r.squared,
                         num.splits           = paste0(num_splits_regr_art, collapse = ","),
                         runtime              = time,
                         n_test = nrow(test_data),
                         n_train = nrow(train_data),
                         n_cal = nrow(cal_data),
                         brier_score5.7 = brier_score5.7,
                         brier_score6.5 = brier_score6.5, # ab hier
                         interval_width = NA,
                         coverage = NA,
                         max_depth = max_depth,
                         num_leaves = num_leaves,
                         num_vars_split = num_vars_split
  ) %>% cbind(information_df),
  tree = art_regr,
  art_prob5.7 = art_prob5.7,
  art_prob6.5 = art_prob6.5,
  test_data = test_data,
  cal_data = cal_data
  ))
}